

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

	char* directory = NULL;
	char* inputDirectory = NULL;

	int peakSize = 200;
	int fragmentLength = TAGADJUST_AUTO;
	int fragmentLengthInput = TAGADJUST_AUTO;
	double maxtbp = -1.0;
	double maxtbpInput = -1.0;
	int centerFlag = 0;
	int nfrFlag = 0;
	int nfrSize = 100;
	int revFlag = 0;
	int normalizeTagThresh = 0;

	long int genomeSize = DEFAULT_GSIZE;
	char* outputfile = NULL;

	
	double normTotal = DEFAULT_NORM_TOTAL;

	char strand = STRAND_BOTH;

	PeakFinder* pf = new PeakFinder();

	if (argc < 2) {
		printCMD();
	}

	char* cmd = new char[10000];
	strcpy(cmd, argv[0]);
	for (int i=1;i<argc;i++) {
		strcat(cmd," ");
		strcat(cmd,argv[i]);
	}
	pf->setCMD(cmd);

	directory = argv[1];
	for (int i=1;i<argc;i++) {
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "!!! First argument needs to be a <tag directory>\n");
				printCMD();
			}
			pf->setDirectory(argv[i]);
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-i")==0) {
				inputDirectory = argv[++i];
			} else if (strcmp(argv[i],"-o")==0) {
				outputfile = argv[++i];
			} else if (strcmp(argv[i],"-fragLength")==0 || strcmp(argv[i],"-len")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragmentLength = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&fragmentLength);
				}
			} else if (strcmp(argv[i],"-inputFragLength")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragmentLengthInput = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&fragmentLengthInput);
				}
			} else if (strcmp(argv[i],"-center")==0) {
				centerFlag = 1;
            } else if (strcmp(argv[i],"-nfr")==0) {
				nfrFlag = 1;
            } else if (strcmp(argv[i],"-region")==0) {
				pf->stitchMode = REGION_MODE_HISTONE;
				pf->regionFlag = 1;
				centerFlag = 0;
            } else if (strcmp(argv[i],"-regionRes")==0) {
				sscanf(argv[++i], "%lf",&(pf->regionSubDivision));
			} else if (strcmp(argv[i],"-gsize")==0) {
				fprintf(stderr, "\tReminder, this recently changed: put actualy genome size, not 2x like before...\n");
				double dsize = 0.0;
				sscanf(argv[++i], "%lf",&dsize);
				genomeSize = (long long int) dsize;
				pf->setGenomeSize(genomeSize);
			} else if (strcmp(argv[i],"-size")==0) {
				sscanf(argv[++i], "%d",&peakSize);
				pf->peakSize = peakSize;
			} else if (strcmp(argv[i],"-norm")==0) {
				sscanf(argv[++i], "%lf",&normTotal);
				pf->normTotal = normTotal;
			} else if (strcmp(argv[i],"-maxBodySize")==0) {
				sscanf(argv[++i], "%d",&(pf->maxBodySize));
			} else if (strcmp(argv[i],"-minDist")==0) {
				int mindist = 0;
				sscanf(argv[++i], "%d",&mindist);
				if (mindist < 1) {
					fprintf(stderr, "!!! -minDist set to a value less than 1 - must be 1 or greater.  Using default...\n");
					mindist = 0;
				}
				pf->minDist = mindist;
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%lf",&maxtbp);
			} else if (strcmp(argv[i],"-inputtbp")==0) {
				sscanf(argv[++i], "%lf",&maxtbpInput);
			} else if (strcmp(argv[i],"-style")==0) {
				i++;
				if (strcmp(argv[i],"histone")==0) {
					pf->style = PEAK_STYLE_HISTONE;
					pf->stitchMode = REGION_MODE_HISTONE;
					pf->regionFlag = 1;
					pf->peakSize = 500;
					pf->minDist = 1000;
					pf->poisson = 0.001;
					pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
				} else {
					fprintf(stderr, "Didn't recognize -style %s !!!\n", argv[i]);
					printCMD();
				}
			} else if (strcmp(argv[i],"-strand")==0) {
				i++;
				if (strcmp(argv[i],"both")==0) {
					strand = STRAND_BOTH;
				} else if (strcmp(argv[i],"separate")==0) {
					strand = STRAND_SEPARATE;
				}
				pf->strand = strand;
			} else if (strcmp(argv[i],"-fdr")==0) {
				double f = 0.001;;
				sscanf(argv[++i], "%lf",&f);
				pf->fdr = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_FDR;
			} else if (strcmp(argv[i],"-poisson")==0) {
				double f = 0.001;;
				sscanf(argv[++i], "%lf",&f);
				pf->poisson = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
			} else if (strcmp(argv[i],"-minTagThreshold")==0) {
				double f = 0.0;
				sscanf(argv[++i], "%lf",&f);
				pf->minTagThresh = f;
			} else if (strcmp(argv[i],"-ntagThreshold")==0) {
				double f = 0.0;
				sscanf(argv[++i], "%lf",&f);
				pf->tagThresh = f;
				pf->minTagThresh = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_THRESH;
				normalizeTagThresh = 1;
			} else if (strcmp(argv[i],"-P")==0) {
				double f = 1.0;
				sscanf(argv[++i], "%lf",&f);
				pf->poissonInput = f;
			} else if (strcmp(argv[i],"-F")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->inputFold = f;
			} else if (strcmp(argv[i],"-L")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->localFold = f;
			} else if (strcmp(argv[i],"-LP")==0) {
				double f = 1.0;
				sscanf(argv[++i], "%lf",&f);
				pf->poissonLocal = f;
			} else if (strcmp(argv[i],"-C")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->clonalFold = f;
			} else if (strcmp(argv[i],"-inputSize")==0) {
				int s = 0;
				sscanf(argv[++i], "%d",&s);
				pf->inputSize = s;
			} else if (strcmp(argv[i],"-localSize")==0) {
				int s = 0;
				sscanf(argv[++i], "%d",&s);
				pf->localSize = s;
			} else {
				printCMD();
			}
		}
	}

	if (maxtbp >= -0.1 || maxtbpInput >= -0.1) {
		pf->setMaxTBP(maxtbp,maxtbpInput);
	}

	if (outputfile != NULL) {
		if (strcmp(outputfile,"auto") == 0) {
			outputfile = new char[10000];
			if (pf->style == PEAK_STYLE_HISTONE) {
				sprintf(outputfile,"%s/regions.txt",directory);
			} else {
				sprintf(outputfile,"%s/peaks.txt",directory);
			}
		}
		pf->setOutputFile(outputfile);
	}

//
	TagLibrary* tags = new TagLibrary(directory);
//	tags->readTagDirectory();
	tags->revStrand = revFlag;
	if (normalizeTagThresh) {
		pf->tagThresh *= tags->totalTags/pf->normTotal;
		pf->minTagThresh = pf->tagThresh;
	}

	tags->setSingleRead(1);
	int tagAdjust;
	if (fragmentLength == TAGADJUST_AUTO) {
		tagAdjust=TAGADJUST_AUTO;
	} else {
		tagAdjust = (int) floor(((double)fragmentLength)/2.0);
		tags->fragmentLengthEstimate = fragmentLength;
	}
	tags->setTagAdjust(tagAdjust);

	TagLibrary* inputTags = NULL;
	if (inputDirectory != NULL) {
		inputTags = new TagLibrary(inputDirectory);
//		inputTags->readTagDirectory();
		inputTags->revStrand = revFlag;
		inputTags->setSingleRead(1);
		int tagAdjustInput;
		if (fragmentLengthInput == TAGADJUST_AUTO) {
			tagAdjustInput=TAGADJUST_AUTO;
		} else {
			tagAdjustInput = (int) floor(((double)fragmentLengthInput)/2.0);
			inputTags->fragmentLengthEstimate = fragmentLengthInput;
		}
		inputTags->setTagAdjust(tagAdjustInput);
	}

	pf->setTagLibraries(tags,inputTags);

	fprintf(stderr, "\tFragment Length = %d\n", tags->fragmentLengthEstimate);

	FILE* fp = stdout;
	if (pf->outputFileName != NULL) {
		fp = fopen(pf->outputFileName, "w");
		if (fp == NULL) {
			fprintf(stderr, "Could not open %s for writing!!!\n", pf->outputFileName);
			exit(1);
		}
	}

	PeakLibrary* peaks = NULL;
	peaks = pf->findPeaks();
    if (centerFlag) {
        peaks->centerPeaks(tags,pf->peakSize,pf->strand);
    } else if (nfrFlag) {
        peaks->centerNFR(tags,pf->peakSize,pf->strand,nfrSize);
    }

	peaks->print(fp);

	if (pf->outputFileName != NULL) {
		fclose(fp);
	}
}


void printCMD() {
	fprintf(stderr, "\n\tUsage: findPeaks <tag directory> [options]\n");
	fprintf(stderr, "\n\tFinds peaks in the provided tag directory.  By default, peak list printed to stdout\n"); 
	fprintf(stderr, "\n\tGeneral analysis options:\n");
	fprintf(stderr, "\t\t-o <filename|auto> (file name for to output peaks, default: stdout)\n"); 
	fprintf(stderr, "\t\t\t\"-o auto\" will send output to \"<tag directory>/peaks.txt\", \".../regions.txt\",\n");
	fprintf(stderr, "\t\t\tor \".../transcripts.txt\" depending on the \"-style\" option\n");
	fprintf(stderr, "\t\t-style <option> (Specialized options for specific analysis strategies)\n");
	fprintf(stderr, "\t\t\thistone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)\n");

	fprintf(stderr, "\n\tchipseq/histone options:\n");
	fprintf(stderr, "\t\t-i <input tag directory> (Experiment to use as IgG/Input/Control)\n"); 
	fprintf(stderr, "\t\t-size <#> (Peak size, default: auto)\n"); 
	fprintf(stderr, "\t\t-minDist <#> (minimum distance between peaks, default: peak size x2)\n"); 
	fprintf(stderr, "\t\t-gsize <#> (Set effective mappable genome size, default: 2e9)\n"); 
	fprintf(stderr, "\t\t-fragLength <#|auto> (Approximate fragment length, default: auto)\n"); 
	fprintf(stderr, "\t\t-inputFragLength <#|auto> (Approximate fragment length of input tags, default: auto)\n"); 
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp to count, 0 = no limit, default: auto)\n");
	fprintf(stderr, "\t\t-inputtbp <#> (Maximum tags per bp to count in input, 0 = no limit, default: auto)\n");
	fprintf(stderr, "\t\t-strand <both|separate> (find peaks using tags on both strands or separate, default:both)\n");
	fprintf(stderr, "\t\t-norm # (Tag count to normalize to, default 10000000)\n");

	fprintf(stderr, "\t\t-region (extends start/stop coordinates to cover full region considered \"enriched\")\n");
	fprintf(stderr, "\t\t\t-regionRes <#> (number of fractions peaks are divided in when extending 'regions', def: 4)\n");
	fprintf(stderr, "\t\t-center (Centers peaks on maximum tag overlap and calculates focus ratios)\n");
	fprintf(stderr, "\t\t-nfr (Centers peaks on most likely nucleosome free region [works best with mnase data])\n");
	fprintf(stderr, "\t\t\t(-center and -nfr can be performed later with \"getPeakTags\"\n");

	fprintf(stderr, "\n\tPeak Filtering options: (set -F/-L/-C to 0 to skip)\n");
	fprintf(stderr, "\t\t-F <#> (fold enrichment over input tag count, default: 4.0)\n");
	fprintf(stderr, "\t\t  -P <#> (poisson p-value threshold relative to input tag count, default: 0.0001)\n");
	fprintf(stderr, "\t\t-L <#> (fold enrichment over local tag count, default: 4.0)\n");
	fprintf(stderr, "\t\t  -LP <#> (poisson p-value threshold relative to local tag count, default: 0.0001)\n");
	fprintf(stderr, "\t\t-C <#> (fold enrichment limit of expected unique tag positions, default: 2.0)\n");
	fprintf(stderr, "\t\t-localSize <#> (region to check for local tag enrichment, default: 10000)\n");
	fprintf(stderr, "\t\t-inputSize <#> (Size of region to search for control tags, default: 2x peak size)\n"); 
	fprintf(stderr, "\t\t-fdr <#> (False discovery rate, default = 0.001)\n"); 
	fprintf(stderr, "\t\t-poisson <#> (Set poisson p-value cutoff, default: uses fdr)\n"); 
	fprintf(stderr, "\t\t-tagThreshold <#> (Set # of tags to define a peak, default: 25)\n"); 
	fprintf(stderr, "\t\t-ntagThreshold <#> (Set # of normalized tags to define a peak, by default uses 1e7 for norm)\n"); 
	fprintf(stderr, "\t\t-minTagThreshold <#> (Absolute minimum tags per peak, default: expected tags per peak)\n"); 

	exit(0);
}
