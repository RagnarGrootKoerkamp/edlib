#include "edlib.h"

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <queue>
#include <string>
#include <unistd.h>
#include <vector>

using namespace std;

int readFastaSequences(const char* path, vector<vector<char>>* seqs);

void printAlignment(const char* query, const char* target, const unsigned char* alignment,
                    const int alignmentLength, const int position, const EdlibAlignMode modeCode);

// For debugging
void printSeq(const vector<char>& seq) {
	for(int i = 0; i < static_cast<int>(seq.size()); i++) printf("%d ", seq[i]);
	printf("\n");
}

// Copied from WFA
bool read_input(FILE* input_file, char** line1, char** line2, int* line1_length,
                int* line2_length) {
	// Parameters
	size_t allocated1 = 0, allocated2 = 0;
	// Read queries
	*line1_length = getline(line1, &allocated1, input_file);
	if(*line1_length == -1) return false;
	*line2_length = getline(line2, &allocated2, input_file);
	if(*line2_length == -1) return false;
	return true;
}

int main(int argc, char* const argv[]) {

	//----------------------------- PARSE COMMAND LINE ------------------------//
	// If true, there will be no output.
	bool silent = false;
	// Alignment mode.
	char mode[16] = "NW";
	// How many best sequences (those with smallest score) do we want.
	// If 0, then we want them all.
	int numBestSeqs         = 0;
	bool findAlignment      = false;
	bool findStartLocations = false;
	int option;
	int kArg       = -1;
	int numRepeats = 1;

	// If "STD" or "EXT", cigar string will be printed. if "NICE" nice representation
	// of alignment will be printed.
	char alignmentFormat[16] = "NICE";

	bool invalidOption = false;
	while((option = getopt(argc, argv, "m:n:k:f:r:spl")) >= 0) {
		switch(option) {
		case 'm': strcpy(mode, optarg); break;
		case 'n': numBestSeqs = atoi(optarg); break;
		case 'k': kArg = atoi(optarg); break;
		case 'f': strcpy(alignmentFormat, optarg); break;
		case 's': silent = true; break;
		case 'p': findAlignment = true; break;
		case 'l': findStartLocations = true; break;
		case 'r': numRepeats = atoi(optarg); break;
		default: invalidOption = true;
		}
	}
	if(optind + 1 != argc || invalidOption) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: %s [options...] <queries.fasta> <target.fasta>\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr,
		        "\t-s  If specified, there will be no score or alignment output (silent mode).\n");
		fprintf(stderr, "\t-m HW|NW|SHW  Alignment mode that will be used. [default: NW]\n");
		fprintf(stderr,
		        "\t-n N  Score will be calculated only for N best sequences (best = with smallest "
		        "score)."
		        " If N = 0 then all sequences will be calculated."
		        " Specifying small N can make total calculation much faster. [default: 0]\n");
		fprintf(stderr, "\t-k K  Sequences with score > K will be discarded."
		                " Smaller k, faster calculation. If -1, no sequences will be discarded. "
		                "[default: -1]\n");
		fprintf(stderr, "\t-p  If specified, alignment path will be found and printed. "
		                "This may significantly slow down the calculation.\n");
		fprintf(
		    stderr,
		    "\t-l  If specified, start locations will be found and printed. "
		    "Each start location corresponds to one end location. This may somewhat slow down "
		    "the calculation, but is still faster then finding alignment path and does not consume "
		    "any extra memory.\n");
		fprintf(
		    stderr,
		    "\t-f NICE|CIG_STD|CIG_EXT  Format that will be used to print alignment path,"
		    " can be used only with -p. NICE will give visually attractive format, CIG_STD will "
		    " give standard cigar format and CIG_EXT will give extended cigar format. [default: "
		    "NICE]\n");
		fprintf(stderr, "\t-r N  Core part of calculation will be repeated N times."
		                " This is useful only for performance measurement, when single execution "
		                "is too short to measure."
		                " [default: 1]\n");
		return 1;
	}
	//-------------------------------------------------------------------------//

	if(strcmp(alignmentFormat, "NICE") && strcmp(alignmentFormat, "CIG_STD") &&
	   strcmp(alignmentFormat, "CIG_EXT")) {
		printf("Invalid alignment path format (-f)!\n");
		return 1;
	}

	EdlibAlignMode modeCode;
	if(!strcmp(mode, "SHW"))
		modeCode = EDLIB_MODE_SHW;
	else if(!strcmp(mode, "HW"))
		modeCode = EDLIB_MODE_HW;
	else if(!strcmp(mode, "NW"))
		modeCode = EDLIB_MODE_NW;
	else {
		printf("Invalid mode (-m)!\n");
		return 1;
	}
	printf("Using %s alignment mode.\n", mode);

	EdlibAlignTask alignTask = EDLIB_TASK_DISTANCE;
	if(findStartLocations) alignTask = EDLIB_TASK_LOC;
	if(findAlignment) alignTask = EDLIB_TASK_PATH;

	// int readResult;
	//  Read queries
	char* queriesFilepath = argv[optind];

	FILE* file = fopen(queriesFilepath, "r");

	// ----------------------------- MAIN CALCULATION ----------------------------- //
	printf("\nComparing queries to target...\n");
	int k                    = kArg;
	unsigned char* alignment = NULL;
	int alignmentLength;
	clock_t start = clock();

	char* query;
	int queryLength;
	char* target;
	int targetLength;
	while(read_input(file, &query, &target, &queryLength, &targetLength)) {
		// Calculate score
		EdlibAlignResult result;
		for(int rep = 0; rep < numRepeats;
		    rep++) { // Redundant repetition, for performance measurements.
			result = edlibAlign(query + 1, queryLength - 1, target + 1, targetLength - 1,
			                    edlibNewAlignConfig(k, modeCode, alignTask, NULL, 0));
			edlibFreeAlignResult(result);
		}
		free(query);
		free(target);
	}

	clock_t finish = clock();
	double cpuTime = static_cast<double>(finish - start) / CLOCKS_PER_SEC;
	printf("\nCpu time of searching: %lf\n", cpuTime);
	// ---------------------------------------------------------------------------- //

	return 0;
}

/** Reads sequences from fasta file.
 * @param [in] path Path to fasta file containing sequences.
 * @param [out] seqs Sequences will be stored here, each sequence as vector of letters.
 * @return 0 if all ok, positive number otherwise.
 */
int readFastaSequences(const char* path, vector<vector<char>>* seqs) {
	seqs->clear();

	FILE* file = fopen(path, "r");
	if(file == 0) return 1;

	bool inHeader      = false;
	bool inSequence    = false;
	const int buffSize = 4096;
	char buffer[buffSize];
	while(!feof(file)) {
		int read = fread(buffer, sizeof(char), buffSize, file);
		for(int i = 0; i < read; ++i) {
			char c = buffer[i];
			if(inHeader) { // I do nothing if in header
				if(c == '\n') inHeader = false;
			} else {
				if(c == '>') {
					inHeader   = true;
					inSequence = false;
				} else {
					if(c == '\r' || c == '\n') continue;
					// If starting new sequence, initialize it.
					if(inSequence == false) {
						inSequence = true;
						seqs->push_back(vector<char>());
					}
					seqs->back().push_back(c);
				}
			}
		}
	}

	fclose(file);
	return 0;
}

void printAlignment(const char* query, const char* target, const unsigned char* alignment,
                    const int alignmentLength, const int position, const EdlibAlignMode modeCode) {
	int tIdx = -1;
	int qIdx = -1;
	if(modeCode == EDLIB_MODE_HW) {
		tIdx = position;
		for(int i = 0; i < alignmentLength; i++) {
			if(alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
		}
	}
	for(int start = 0; start < alignmentLength; start += 50) {
		// target
		printf("T: ");
		int startTIdx = -1;
		for(int j = start; j < start + 50 && j < alignmentLength; j++) {
			if(alignment[j] == EDLIB_EDOP_INSERT)
				printf("-");
			else
				printf("%c", target[++tIdx]);
			if(j == start) startTIdx = tIdx;
		}
		printf(" (%d - %d)\n", max(startTIdx, 0), tIdx);

		// match / mismatch
		printf("   ");
		for(int j = start; j < start + 50 && j < alignmentLength; j++) {
			printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
		}
		printf("\n");

		// query
		printf("Q: ");
		int startQIdx = qIdx;
		for(int j = start; j < start + 50 && j < alignmentLength; j++) {
			if(alignment[j] == EDLIB_EDOP_DELETE)
				printf("-");
			else
				printf("%c", query[++qIdx]);
			if(j == start) startQIdx = qIdx;
		}
		printf(" (%d - %d)\n\n", max(startQIdx, 0), qIdx);
	}
}
