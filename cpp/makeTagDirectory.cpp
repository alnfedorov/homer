
// Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include "SeqTag.h"

void printCMD();

int main(int argc, char** argv) {

	std::string directory;
    std::vector<std::string> files;

	int sspeFlag = 0;
	int format = FORMAT_UNKNOWN;
	int mode = MODE_UNIQUE;
    int pairedEndFlag = 0;
	double tbp = 0.0;

    size_t minReadLength = 0;
	size_t maxReadLength = std::numeric_limits<size_t>::max();

    int fragLength = FRAGMENT_LEN_AUTO;

	if (argc < 2) {
		printCMD();
	}

	for (int i=1;i<argc;i++) {
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "\n!!!!!!!!!!!!\n\tNEED to specify directory with first argument!!!\n");
				printCMD();
			}
			directory = argv[i];
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-keep")==0) {
				mode = MODE_KEEPONE;
			} else if (strcmp(argv[i],"-keepOne")==0) {
				mode = MODE_KEEPONE;
			} else if (strcmp(argv[i],"-unique")==0) {
				mode = MODE_UNIQUE;
            } else if (strcmp(argv[i],"-forceBED")==0) {
				mode = MODE_UNIQUE;
			} else if (strcmp(argv[i],"-sspe")==0) {
				sspeFlag = 1;
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%lf", &tbp);
			} else if (strcmp(argv[i],"-len")==0 || strcmp(argv[i],"-fragLength")==0) {
				if (i+1 >= argc ) {
					fprintf(stderr, "Error specifying -len\n");
					printCMD();
				}
				i++;
				if (strcmp(argv[i],"given") == 0) {
					fragLength = FRAGMENT_LEN_GIVEN;
				} else if (strcmp(argv[i],"pe") == 0 || strcmp(argv[i],"PE") == 0) {
					fragLength = FRAGMENT_LEN_PE;
				} else {
					sscanf(argv[i],"%d", &fragLength);
				}
			} else if (strcmp(argv[i],"-format")==0) {
				if (i+1 >= argc ) {
					fprintf(stderr, "Error specifying -format\n");
					printCMD();
				}
				i++;
				if (strcmp(argv[i],"bed") == 0) {
					format = FORMAT_BED;
				} else if (strcmp(argv[i],"sam") == 0) {
					format = FORMAT_SAM;
				} else {
					fprintf(stderr, "Error specifying -format\n");
					printCMD();
				}
			} else {
				fprintf(stderr, "!!! Couldn't recognize: %s !!!\n", argv[i]);
				printCMD();
			}
		} else {
		    files.emplace_back(argv[i]);
			if (pairedEndFlag == 0) {
				int index = 0;
				while (argv[i][index] != '\0') {
					if (argv[i][index] == ',') {
						pairedEndFlag = 1;
						fprintf(stderr, "\tMaking paired end tag directory\n");
						break;
					}
					index++;
				}
			}
			fprintf(stderr,"\tWill parse file: %s\n", files.back().c_str());
		}
	}

	auto tags = new TagLibrary();

	tags->maxReadLength = maxReadLength;
	tags->minReadLength = minReadLength;
    tags->minmapq = 10.0;
    tags->sspeFlag = sspeFlag;
//    tags->setFragLength(fragLength);
    tags->pairedEndFlag = pairedEndFlag;

    tags->parseAlignmentFiles(files,format,mode);

//    if (fragLength == FRAGMENT_LEN_PE) tags->setFragLength(FRAGMENT_LEN_AUTO);
//
//	if (fragLength != FRAGMENT_LEN_AUTO && fragLength != FRAGMENT_LEN_GIVEN) {
//		tags->setFragLength(fragLength);
//	}
	if (true) {
		int max = MAX_TAGS_PER_BP;
		double* d = tags->getTagCountDistribution(NULL, max);
		if (d != NULL) delete []d;
	}

	// now that we've characterized count distribution, force it with tbp if needed
	if (tbp > 0.0) {
		//need to save tbp modification...
		tags->setMaxTBP(tbp);
		fprintf(stderr, "\tRestricting tags per bp...\n");
		tags->setMaxTBP(0);
	}

	if (fragLength != FRAGMENT_LEN_AUTO && fragLength != FRAGMENT_LEN_GIVEN) {
		fprintf(stderr, "\tForcing fragment length = %d\n", fragLength);
//		tags->setFragLength(fragLength);
	}
}


void printCMD() {
	fprintf(stderr, "\n\tUsage: makeTagDirectory <directory> <alignment file 1> [file 2] ... [options]\n");
	fprintf(stderr, "\n\tCreates a platform-independent 'tag directory' for later analysis.\n"); 
	fprintf(stderr, "\tCurrently BED and sam files are accepted. The program will try to\n");
	fprintf(stderr, "\tautomatically detect the alignment format if not specified.  Program will also\n");
	fprintf(stderr, "\tunzip *.gz, *.bz2, and *.zip files and convert *.bam to sam files on the fly\n");
	fprintf(stderr, "\tExisting tag directories can be added or combined to make a new one using -d/-t\n");
	fprintf(stderr, "\tIf more than one format is needed and the program cannot auto-detect it properly,\n");
	fprintf(stderr, "\tmake separate tag directories by running the program separately, then combine them.\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-fragLength <# | given | pe> (Set estimated fragment length or use PE length - given: use read lengths)\n");
	fprintf(stderr, "\t\t\tBy default treats the sample as a single read ChIP-Seq experiment\n");
	fprintf(stderr, "\t\t-format <X> where X can be: (with column specifications underneath)\n");
	fprintf(stderr, "\t\t\tbed - BED format files:\n");
	fprintf(stderr, "\t\t\t\t(1:chr,2:start,3:end,4:+/- or read name,5:# tags,6:+/-)\n");
	fprintf(stderr, "\t\t\tsam - SAM formatted files (use samTools to covert BAMs into SAM if you have BAM)\n");
	fprintf(stderr, "\t\t\t\t-unique (keep if there is a single best alignment based on mapq)\n");
	fprintf(stderr, "\t\t\t\t\t-mapq <#> (Minimum mapq for -unique, default: 10, set negative to use AS:i:/XS:i:)\n");
	fprintf(stderr, "\t\t\t\t-keepOne (keep one of the best alignments even if others exist)\n");
	fprintf(stderr, "\t\t\t\t-sspe (strand specific, paired-end reads[flips strand of 2nd read to match])\n");
	fprintf(stderr, "\t\t\t\t-read1/-read2 (only analyze 1st or 2nd read for PE sequencing)\n");
	fprintf(stderr, "\t\t-t <tag file> [tag file 2] ... (add tag file i.e. *.tags.tsv to new tag directory)\n");
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp, default: no maximum)\n");
	fprintf(stderr, "\t\t-precision <1|2|3> (number of decimal places to use for tag totals, default: 1)\n");
	fprintf(stderr, "\n");

	exit(0);
}
