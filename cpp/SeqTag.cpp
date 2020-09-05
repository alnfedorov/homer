#include "SeqTag.h"

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

#define BUFFER 100000

void split(char *string, char **cols, int &numCols, char delim) {
    cols[0] = string;
    numCols = 1;
    char delim2 = 0;
    if (delim == WHITE_SPACE) {
        delim = '\t';
        delim2 = 32;
    }
    int len = strlen(string);
    for (int i = 0; i < len; i++) {
        if (string[i] == delim || string[i] == delim2) {
            string[i] = '\0';
            cols[numCols] = &(string[i + 1]);
            numCols++;
        } else if (string[i] == '\n') {
            string[i] = '\0';
        } else if (string[i] == '\r') {
            string[i] = '\0';
        }
    }
}


// class PeakFinder ---------------------------------------------------------------------

PeakFinder::PeakFinder() {
    name = NULL;
    directory = NULL;
    outputFileName = NULL;
    excludePeaksFile = NULL;
    peakSize = 0;
    localSize = 10000;
    inputSize = 0;
    tagThresh = 0;
    minTagThresh = -1.0;
    minDist = 0;
    totalTags = 0.0;
    tagsUsedForClustering = 0.0;
    maxtbp = 0.0;
    maxtbpInput = 0.0;
    strand = BOTH_STRANDS;
    regionSubDivision = PEAKFRACTION_REGION_MINDIST;
    gsize = DEFAULT_GSIZE;
    //gsize *= 2;
    fdr = 0.001;
    fdrThresh = 0.0;
    normTotal = DEFAULT_NORM_TOTAL;
    poisson = 0.0;
    poissonInput = 0.0001;
    poissonLocal = 0.0001;
    poissonThresh = 0.0;
    tbpAuto = 1;
    tbpThreshold = 0.01;
    tbp = 0.0;
    tbpInput = 0.0;
    diffMode = DIFFPEAK_MODE_DIFF;
    filterMode = PEAKFINDER_FILTER_MODE_FDR;
    fdrTable = NULL;
    poissonTable = NULL;
    tagsInPeaks = 0.0;
    fdrSize = PEAKFINDER_FDRSIZE;
    extraheader = NULL;
    stitchMode = REGION_MODE_HISTONE;
    maxBodySize = GROSEQ_MAXBODYSIZE;

    style = PEAK_STYLE_HISTONE;

    inputFold = 4.0;
    localFold = 4.0;
    clonalFold = 2.0;
    numPeaks = 0;

    cmd = NULL;

    tags = NULL;
    input = NULL;
}

PeakFinder::~PeakFinder() {
    if (name != NULL) delete[]name;
    if (extraheader != NULL) delete[]extraheader;
    if (directory != NULL) delete[]directory;
}

void PeakFinder::addHeader(char *str) {
    if (str == NULL) return;
    char *header = extraheader;
    int len = 0;
    if (header != NULL) len = strlen(header) + 1;
    extraheader = new char[len + strlen(str) + 2];
    extraheader[0] = '\0';
    if (header != NULL) {
        strcpy(extraheader, header);
        delete[]header;
    }
    strcat(extraheader, "\t");
    strcat(extraheader, str);
}

void PeakFinder::setTagLibraries(TagLibrary *exp, TagLibrary *i) {
    tags = exp;
    input = i;
}

void PeakFinder::setMaxTBP(double maxT, double maxI) {
    tbpAuto = 0;
    maxtbp = maxT;
    maxtbpInput = maxI;
    if (maxtbp < 0) maxtbp = maxtbpInput;
    if (maxtbpInput < 0) maxtbpInput = maxtbp;
}

void PeakFinder::setGenomeSize(long long int genomeSize) {
    gsize = genomeSize;
}

void PeakFinder::setCMD(char *str) {
    if (str == NULL) return;
    cmd = new char[strlen(str) + 1];
    strcpy(cmd, str);
}

void PeakFinder::setOutputFile(char *fname) {
    if (fname == NULL) return;
    outputFileName = new char[strlen(fname) + 1];
    strcpy(outputFileName, fname);
}

void PeakFinder::setDirectory(char *newname) {
    if (newname == NULL) return;
    directory = new char[strlen(newname) + 1];
    strcpy(directory, newname);
}

void PeakFinder::determineMaxTBP() {
    maxtbp = 1.0;
    maxtbpInput = 1.0;
    fprintf(stderr, "\tTotal Tags = %.1lf\n", totalTags);
    fprintf(stderr, "\tTags per bp = %.6lf\n", tbp);
    if (tbp > tbpThreshold) maxtbp = floor(tbp / tbpThreshold);
    if (tbpInput > tbpThreshold) maxtbpInput = floor(tbp / tbpThreshold);
    fprintf(stderr, "\tMax tags per bp set automatically to %.1f\n", maxtbp);
}

void PeakFinder::checkParameters() {
    if (tags != NULL) {

        if (tags->gsizeEstimate > 10 && tags->gsizeEstimate < gsize &&
            gsize == ((long long int) DEFAULT_GSIZE)) {
            gsize = tags->gsizeEstimate;
            fprintf(stderr, "\t!!! Estimated genome size (from tag directory) is smaller than default\n");
            fprintf(stderr, "\t    genome size.  Using estimate (%lld) [to change specify -gsize]\n", gsize);
            //gsize *= 2;
        }
        tbp = tags->totalTags / (double) gsize;
        //fprintf(stderr, "%lf\t%lf\t%lf\n", tbp,tags->totalTags,(double)gsize);
        if (strand != STRAND_BOTH) tbp /= 2.0;
        tags->tbp = tbp;
        totalTags = tags->totalTags;
    }
    if (input != NULL) {
        tbpInput = input->totalTags / (double) gsize;
        if (strand != STRAND_BOTH) tbpInput /= 2.0;
        input->tbp = tbpInput;
        //tagAdjustInput = input->fragmentLengthEstimate;
    }
    if (peakSize == 0) {
        peakSize = 200;
    }
    if (minDist == 0) {
        minDist = (int) floor(peakSize * FINDPEAKS_MINDIST_DEFAULT);
    }
    if (tbpAuto) determineMaxTBP();
    tags->setMaxTBP(maxtbp);
    if (input != NULL) {
        tags->setMaxTBP(maxtbp);
    }

    if (strand == STRAND_SEPARATE) {
        fprintf(stderr, "\tFinding tags on separate strands: doubling effective genome size\n");
    }

}

PeakLibrary *PeakFinder::findPeaks() {
    checkParameters();

    int curMinDist = minDist;
    if (regionFlag == 0) {
        fprintf(stderr, "\tFinding peaks of size %d, no closer than %d\n", peakSize, minDist);
    } else {
        fprintf(stderr, "\tInitially finding peaks of size %d bp for stitching into regions", peakSize);
        fprintf(stderr, " with %d bp stitching size)\n", minDist);
        curMinDist = (int) (-1 * (peakSize / regionSubDivision));
    }

    if (minTagThresh < 0) {
        minTagThresh = tbp * peakSize - 1.0;
    }

    PeakLibrary *filteredPeaks = NULL;
    int halfPeakSize = (peakSize + 1) / 2;

    PeakLibrary *putativePeaks = NULL;
    putativePeaks = tags->findPutativePeaks(peakSize, curMinDist, strand, minTagThresh);

    tagsUsedForClustering = tags->getAdjustedTagTotal();

    fprintf(stderr, "\t\tTags Used for cluster (less clonal tags) = %.1lf / %.1lf\n",
            tagsUsedForClustering, totalTags);

    // remove peaks not meeting fdr levels or poisson p-value or absolute tag thresh
    approxFdrTable();
    addHeader((char *) "findPeaks Score");
    filteredPeaks = filterPeaks(putativePeaks);
    delete putativePeaks;
    putativePeaks = NULL;

    filteredPeaks->setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfPeakSize, halfPeakSize);

    if (input != NULL && inputFold > 0) {
        tags->setMaxTBP(0);
        input->setMaxTBP(0);
        if (tags->totalTags > input->totalTags) {
            addHeader((char *) "Total Tags (normalized to Control Experiment)");
            addHeader((char *) "Control Tags");
        } else {
            addHeader((char *) "Total Tags");
            addHeader((char *) "Control Tags (normalized to IP Experiment)");
        }
        addHeader((char *) "Fold Change vs Control");
        addHeader((char *) "p-value vs Control");
        int halfInputSize = halfPeakSize;
        if (inputSize == 0) {
            halfInputSize = halfPeakSize * 2;
        } else {
            halfInputSize = inputSize / 2;
        }
        int strOutputFlag = 1;
        if (regionFlag) strOutputFlag = 0;
        PeakLibrary *inputFiltered = filteredPeaks->getDifferentialPeaks(tags, input,
                                                                         inputFold, poissonInput, diffMode,
                                                                         -1 * halfInputSize, halfInputSize, strand,
                                                                         strOutputFlag);
        delete filteredPeaks;
        filteredPeaks = inputFiltered;
    }

    if (regionFlag) {
        fprintf(stderr, "\tStitching together putative peaks into regions\n");
        PeakLibrary *regions = filteredPeaks->stitchRegions(minDist, stitchMode);
        delete filteredPeaks;
        filteredPeaks = regions;

        if (input != NULL && inputFold > 0) {
            fprintf(stderr, "\tChecking regions against input...\n");
            tags->setMaxTBP(0);
            input->setMaxTBP(0);
            filteredPeaks->setPeakTagSizeFixed(0, 0);
            int strOutputFlag = 1;
            PeakLibrary *inputFiltered = filteredPeaks->getDifferentialPeaks(tags, input,
                                                                             inputFold, poissonInput, diffMode,
                                                                             0, 0, strand, strOutputFlag);
            delete filteredPeaks;
            filteredPeaks = inputFiltered;
        }

    } else {
        if (localFold > 0) {
            tags->setMaxTBP(0);
            addHeader((char *) "Fold Change vs Local");
            addHeader((char *) "p-value vs Local");
            PeakLibrary *localFiltered = filteredPeaks->filterLocalPeaks(tags,
                                                                         peakSize, localSize, localFold,
                                                                         poissonLocal, diffMode, strand);
            delete filteredPeaks;
            filteredPeaks = localFiltered;
        }
    }

    if (clonalFold > 0) {
        tags->setMaxTBP(0);
        addHeader((char *) "Clonal Fold Change");
        PeakLibrary *clonalFiltered = filteredPeaks->filterClonalPeaks(tags,
                                                                       peakSize, clonalFold, diffMode, strand);
        delete filteredPeaks;
        filteredPeaks = clonalFiltered;
    }

    PeakLibrary *excludePeaks = NULL;
    if (excludePeaksFile != NULL) {
        fprintf(stderr, "\tUsing provided excluded peaks to modify super enhancers\n");
        excludePeaks = new PeakLibrary(excludePeaksFile, PEAK_READ_MODE_NORMAL);
    }

    //normalize tag counts
    tags->setMaxTBP(0);

    filteredPeaks->setPeakTagSizeFixed(0, 0);
    //filteredPeaks->setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);

    filteredPeaks->setDefaultPeakOrder();
    filteredPeaks->tagsInPeaks = 0.0;

    //int tagsIndex = filteredPeaks->addTagLibrary(tags);
    //Doubletable* expTags = filteredPeaks->countPeakTags(tagsIndex,0,0,strand,COUNT_MODE_TOTAL);
    Doubletable *expTags = filteredPeaks->countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, excludePeaks);

    for (int i = 0; i < filteredPeaks->numPeaks; i++) {
        filteredPeaks->peakOrder[i]->v = expTags->search(filteredPeaks->peakOrder[i]->name);
        filteredPeaks->tagsInPeaks += filteredPeaks->peakOrder[i]->v;
    }
    filteredPeaks->normalizePeakScore(normTotal / tags->totalTags);

    tagsInPeaks = filteredPeaks->tagsInPeaks;
    delete expTags;

    numPeaks = filteredPeaks->numPeaks;
    fprintf(stderr, "\tTotal Peaks identified = %d\n", numPeaks);
    filteredPeaks->setDefaultPeakOrder();

    return filteredPeaks;
}

void PeakLibrary::normalizePeakScore(float normFactor) {
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->v *= normFactor;
        delete[](keys[i]);
    }
    delete[]keys;
}

PeakLibrary *PeakFinder::filterPeaks(PeakLibrary *putativePeaks) {

    double *peakHeights = new double[fdrSize];
    for (int i = 0; i < fdrSize; i++) {
        peakHeights[i] = 0;
    }

    char **keys = putativePeaks->peaks->keys();
    for (int i = 0; i < putativePeaks->peaks->total; i++) {
        Peak *p = (Peak *) putativePeaks->peaks->search(keys[i]);
        int index = (int) p->v;
        if (index >= fdrSize) {
            index = fdrSize - 1;
        }
        peakHeights[index] += 1.0;
    }

    fprintf(stderr, "\t\tThreshold\tPeak Count\tExpected Peak Count\tFDR\tPoisson\n");
    for (int i = fdrSize - 2; i >= 0; i--) {
        peakHeights[i] += peakHeights[i + 1];
        double cfdr = 0.0;
        if (peakHeights[i] > 0) {
            cfdr = fdrTable[i] / peakHeights[i];
        }
        if (cfdr < fdr) {
            fdrThresh = (float) i;
        }
        if (cfdr > fdr / 10000.0) {// && i <= 50) {
            if (cfdr > 1.0) cfdr = 1.0;
            fprintf(stderr, "\t\t%d\t%.3lf\t%.3lf\t%.2le\t%.2le\n", i, peakHeights[i], fdrTable[i], cfdr,
                    poissonTable[i]);
        }
    }

    if (filterMode == PEAKFINDER_FILTER_MODE_FDR) {
        poisson = poissonTable[(int) fdrThresh];
        fprintf(stderr, "\t%.2lf%% FDR Threshold set at %.1f (poisson pvalue ~ %.2le)\n",
                100.0 * fdr, fdrThresh, poisson);
    }

    auto *goodPeaks = new PeakLibrary();
    int numGood = 0;
    char *strScore = new char[10000];
    tagsInPeaks = 0.0;

    double threshold2Use = fdrThresh;
    if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
        threshold2Use = poissonThresh;
    } else if (filterMode == PEAKFINDER_FILTER_MODE_THRESH) {
        threshold2Use = tagThresh;
    }

    for (int i = 0; i < putativePeaks->peaks->total; i++) {
        Peak *p = (Peak *) putativePeaks->peaks->search(keys[i]);
        if (p->v >= threshold2Use - 000000.1) {
            tagsInPeaks += (double) p->v;
            sprintf(strScore, "%f", p->v);
            p->addData(strScore);
            goodPeaks->addPeak(p);
            numGood++;
        }
        delete[](keys[i]);
    }
    delete[]keys;
    delete[]peakHeights;
    delete[]strScore;
    fprintf(stderr, "\t%d peaks passed threshold\n", numGood);
    goodPeaks->sortChr();

    return goodPeaks;
}

void PeakFinder::approxFdrTable() {

    tbp = tagsUsedForClustering / ((double) gsize);

    //tpp = tags per peak
    tpp = tbp * ((double) peakSize);
    if (strand != STRAND_BOTH) {
        tpp = tbp * ((double) peakSize / 2.0);
    }
    fprintf(stderr, "\tExpected tags per peak = %lf (tbp = %lf)\n", tpp, tbp);

    //this is the toughest part - and is a little empirical based on randomized tag positions
    double numTests = 2.0 * gsize / ((double) peakSize);
    if (regionFlag) {
        numTests = 2.0 * gsize / ((double) (peakSize / regionSubDivision / 2.0));
    } else {
        numTests = 2.0 * gsize / ((double) (minDist / 2.0));
    }
    if (strand == STRAND_SEPARATE) {
        numTests *= 2;
        fprintf(stderr, "\tFinding tags on separate strands: doubling effective genome size\n");
    }


    fdrSize = PEAKFINDER_FDRSIZE;
    fdrTable = new double[fdrSize];
    poissonTable = new double[fdrSize];


    double cum = 0.0;
    int minFDRIndex = -1;

    for (int i = PEAKFINDER_FDRSIZE - 1; i >= 0; i--) {
        double lp = logPoisson(i, tpp);
        lp = exp(lp);
        cum += lp;
        if (numTests * cum > fdrThresh) {
            if (minFDRIndex < 0) {
                fdrSize = i;
                minFDRIndex = i;
            }
            fdrTable[i] = numTests * cum;
            poissonTable[i] = cum;
        }
        if (cum < poisson) {
            poissonThresh = (float) i;
        }
    }
    if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
        fprintf(stderr, "\tPoisson Threshold set at %d tags\n", (int) poissonThresh);
    }
    return;
    //exit(0);

    cum = 0.0;

    for (int i = 0; i < fdrSize; i++) {
        double lp = logPoisson(i, tpp);
        lp = exp(lp);
        //fprintf(stderr, "\t%d\t%.9lf\t%.9lf\t%.3f\n",i, lp, 1-cum,numTests*(1-cum));
        double poissonCum = 1 - cum;
        fdrTable[i] = numTests * poissonCum;
        poissonTable[i] = poissonCum;
        fprintf(stderr, "\t%d\t%lf\t%le\t%le\t%le\n", i, fdrTable[i], numTests, poissonCum, lp);
        if (poissonThresh < 0.01 && poissonCum < poisson) {
            poissonThresh = (float) i;
            if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
                fprintf(stderr, "\tPoisson Threshold set at %d tags\n", i);
            }
        }
        cum += lp;
    }

}


// class PeakLibrary --------------------------------------------------------------------

PeakLibrary::PeakLibrary() {
    initialize(PEAK_LIBRARY_DEFAULT_SIZE);
}

PeakLibrary::PeakLibrary(int expectedNumPeaks) {
    initialize(expectedNumPeaks);
}

PeakLibrary::PeakLibrary(char *fname, int mode) {
    if (fname == NULL) {
        initialize(PEAK_LIBRARY_DEFAULT_SIZE);
        return;
    }
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "!!! Could not open peak file: %s !!!\n", fname);
        initialize(PEAK_LIBRARY_DEFAULT_SIZE);
        return;
    }
    int numLines = 0;
    char *buf = new char[BUFFER];
    while (fgets(buf, BUFFER, fp) != NULL) {
        numLines++;
    }
    fclose(fp);
    delete[]buf;
    initialize(numLines + 1000);
    readPeakFile(fname, mode);
}

void PeakLibrary::initialize(int expectedNumPeaks) {
    chrs = new Hashtable(10000);
    peaks = new Hashtable(expectedNumPeaks);
    duplicates = new Inttable((int) (expectedNumPeaks / 5));
    exps = NULL;
    name = NULL;
    genome = NULL;
    numPeaks = 0;
    numExps = 0;
    avgPeakSize = 0.0;
    tagsInPeaks = 0.0;
    fixedFlag = 0;
    peakOrder = NULL;
    duplicateWarningFlag = 1;
}

PeakLibrary::~PeakLibrary() {
    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrPeaks *ct = (ChrPeaks *) chrs->search(keys[i]);
        delete ct;
        delete[](keys[i]);
    }
    delete[]keys;
    delete chrs;
    keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        delete p;
        delete[](keys[i]);
    }
    delete[]keys;
    delete peaks;
    delete duplicates;
    if (name != NULL) delete[]name;
    if (genome != NULL) delete[]genome;
    if (exps != NULL) {
        delete[]exps;
    }
    if (peakOrder != NULL) delete[]peakOrder;

}

void PeakLibrary::readPeakFile(char *filename, int mode) {

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open peak file (%s)\n", filename);
        return;
    }

    char *buf = new char[BUFFER];
    char *countingName = new char[BUFFER];
    char *posname = new char[BUFFER];
    char **line = new char *[BUFFER];
    char *chr = NULL;
    char *name = NULL;
    char *ann = NULL;
    int numCols = 0;


    avgPeakSize = 0.0;
    int NN = 0;
    int start = 0;
    int end = 0;
    int midpoint = NULL_REF;
    float tagCount = 0;
    float focusRatio = 0;
    int dir = 0;
    //unsigned int currentPriority = 0;
    unsigned int priority = 0;
    unsigned int intergenicPriority = 1000000000;
    char defaultPrefix[8] = "default";
    char defaultName[50] = "";
    int defaultCount = 0;
    int mappability = -1;
    int forceFormat = FORMAT_UNKNOWN;
    if (mode == PEAK_READ_MODE_2DBED) {
        forceFormat = FORMAT_BED;
        //mode=PEAK_READ_MODE_NORMAL;
    }

    while (fgets(buf, BUFFER, fp) != NULL) {
        split(buf, line, numCols, '\t');
        if (numCols < 3) continue;
        if (line[0][0] == '#') continue;

        //check if there is a header (i.e. 3rd column should be integer
        if (line[2][0] != '-' && (line[2][0] < 48 || line[2][0] > 57)) {
            continue;
        }

        int format = forceFormat;

        name = NULL;
        chr = NULL;
        start = -1;
        end = -1;
        dir = -1;


        int col1notNum = 0;
        int col2notNum = 0;
        int col3notNum = 0;
        col1notNum = checkInt(line[1]);
        col2notNum = checkInt(line[2]);
        if (col2notNum) {
            continue;
        }
        int col3strand = -1;
        int col4strand = -1;
        int col5strand = -1;

        if (numCols == 3) {
            if (col1notNum) {
                continue;
            }
            format = FORMAT_BED;
        } else {
            col3notNum = checkInt(line[3]);
            if (col1notNum && col3notNum) {
                continue;
            }
            if (numCols == 4) {
                if (col3notNum) {
                    format = FORMAT_BED;
                    col3strand = checkStrand(line[3]);
                    if (col3strand == -1) {
                        name = line[3];
                    } else {
                        dir = col3strand;
                    }
                } else {
                    format = FORMAT_PEAK;
                }
            } else {
                col4strand = checkStrand(line[4]);
                col3strand = checkStrand(line[3]);
                if (col4strand != -1 && col3notNum == 0) {
                    format = FORMAT_PEAK;
                    dir = col4strand;
                } else if (col3notNum) {
                    format = FORMAT_BED;
                    if (col3strand != -1) {
                        dir = col3strand;
                    } else {
                        name = line[3];
                    }
                    if (numCols > 5) {
                        col5strand = checkStrand(line[5]);
                        if (col5strand != -1) {
                            dir = col5strand;
                        }
                    }
                } else {
                    continue;
                }
            }
        }

        if (format == FORMAT_BED) {
            chr = line[0];
            sscanf(line[1], "%d", &start);
            sscanf(line[2], "%d", &end);
            start++;
        } else if (format == FORMAT_PEAK) {
            name = line[0];
            chr = line[1];
            sscanf(line[2], "%d", &start);
            sscanf(line[3], "%d", &end);
        } else {
            continue;
        }

        if (dir == -1) {
            dir = 0;
        }
        if (name == NULL) {
            defaultCount++;
            sprintf(defaultName, "%s-%d", defaultPrefix, defaultCount);
            name = defaultName;
        }

        //initialize optional items
        midpoint = NULL_REF;
        tagCount = 0.0;
        focusRatio = 0.0;
        ann = NULL;

        avgPeakSize += (end - start);
        NN++;

        if (mode == PEAK_READ_MODE_NORMAL && format == FORMAT_PEAK) {
            if (numCols > 5) sscanf(line[5], "%f", &tagCount);
            if (numCols > 6) sscanf(line[6], "%f", &focusRatio);
        } else if (mode == PEAK_READ_MODE_ANNOTATION && format == FORMAT_PEAK) {
            if (numCols > 5) {
                ann = line[5];
                if (strncmp(ann, "I", 1) == 0) {
                    //currentPriority = intergenicPriority++;
                    intergenicPriority++;
                } else {
                    //currentPriority = priority++;
                    priority++;
                }
            } else {
                //currentPriority = priority++;
                priority++;
            }
        } else if (mode >= PEAK_READ_MODE_COUNTING) {
            sprintf(countingName, "%d-%d", mode, NN);
            name = countingName;
        }
        if (mode == PEAK_READ_MODE_2DBED) {
            //if (start % 100 == 1) start--;
            sprintf(posname, "%s:%d-%d", chr, start, end);
            name = posname;
            if (numCols > 7) sscanf(line[7], "%f", &tagCount);
            if (numCols > 8) sscanf(line[8], "%f", &focusRatio);
        }

        addPeak(name, chr, start, end, midpoint, (char) dir, tagCount, focusRatio, ann, mappability, priority);
        //if (numPeaks % 1000000==0) fprintf(stderr, "\t\tloaded %d peaks\n", numPeaks);
    }
    fclose(fp);

    avgPeakSize /= (double) NN;

    delete[]buf;
    delete[]line;
    delete[]countingName;
    delete[]posname;
    sortChr();
}

void PeakLibrary::sortChr() {

    numPeaks = 0;

    char **chrKeys = chrs->keys();
    qsort(chrKeys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrPeaks *ct = (ChrPeaks *) chrs->search(chrKeys[i]);
        ct->sort();
        numPeaks += ct->numPeaks;
        delete[](chrKeys[i]);
    }
    delete[]chrKeys;
}

Peak *PeakLibrary::addPeak(char *name, char *chr, int start, int end, int midpoint, char dir,
                           float value, float ratio, char *extraData, int mappability, unsigned int priority) {
    Peak *oldPeak = NULL;
    Peak *p = NULL;
    char *tmpName = NULL;
    if (name == NULL) {
        int L = strlen(chr) + 26 + 1 + 6;
        tmpName = new char[L];
        sprintf(tmpName, "%s:%d-%d", chr, start, end);
        name = tmpName;
    }
    if (name != NULL) oldPeak = (Peak *) peaks->search(name);
    static int warningIssued = 0;
    if (oldPeak != NULL) {
        int dupCount = duplicates->search(name);
        if (dupCount == EMPTY_INT) {
            dupCount = 1;
        }
        dupCount++;
        //fprintf(stderr, "dupCount=%d \t%s\n",dupCount,name);
        duplicates->insert(dupCount, name);

        char *newname = new char[strlen(name) + 20];
        sprintf(newname, "%s--%d", name, dupCount);
        if (priority < 1) {
            if (warningIssued == 0 && duplicateWarningFlag) {
                fprintf(stderr, "\tDuplicate peak name (%s) - this could potentially cause problems\n", name);
                fprintf(stderr, "\t\tSometimes unavoidable for BED/2DBED formats\n");
                fprintf(stderr, "\t\tNew name for this peak is %s\n", newname);
            } else if (warningIssued % 1000 == 0 && duplicateWarningFlag) {
                fprintf(stderr, "\t\tWarning over %d peaks with duplicate names\n", warningIssued);
            }
            warningIssued++;
        }
        oldPeak = (Peak *) peaks->search(newname);
        if (oldPeak != NULL && duplicateWarningFlag) {
            fprintf(stderr, "There's a problem!!! %s - %d %s\n", newname, dupCount, name);
        }
        p = new Peak(newname, name, chr, start, end, midpoint, dir, value, ratio, extraData, mappability, priority);
        delete[]newname;
    } else {
        p = new Peak(name, NULL, chr, start, end, midpoint, dir, value, ratio, extraData, mappability, priority);
    }
    peaks->insert(p, p->name);
    ChrPeaks *cp = (ChrPeaks *) chrs->search(chr);
    if (cp == NULL) {
        cp = new ChrPeaks();
        chrs->insert(cp, chr);
    }
    cp->addPeak(p);
    numPeaks++;
    if (tmpName != NULL) delete[]tmpName;
    return p;
}

void PeakLibrary::addPeak(Peak *p) {
    char *n = p->name;
    if (p->ogname != NULL) {
        n = p->ogname;
    }
    addPeak(n, p->chr, p->start, p->end, p->refPos, p->strand, p->v, p->focusRatio, p->data,
            p->uniqMap, p->priority);
}

void PeakLibrary::print(FILE *fp) {
    char **keys = peaks->keys();
    for (int i = 0; i < numPeaks; i++) {
        //fprintf(stderr, "\t%d\t%s\n", i, keys[i]);
    }
    sortKeys(keys);
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->print(fp);
        delete[](keys[i]);
    }
    delete[]keys;
}

void PeakLibrary::sortKeys(char **keys) {
    std::sort(keys, keys + numPeaks, [&](const char *a, const char *b) {
        Peak *p1 = (Peak *) peaks->search(*((char **) a));
        Peak *p2 = (Peak *) peaks->search(*((char **) b));
        float v1 = p1->v;
        float v2 = p2->v;
        if (v1 < v2) return 1;
        if (v2 < v1) return -1;
        char *c1 = p1->chr;
        char *c2 = p2->chr;
        int c = chrcmp((void *) &c1, (void *) &c2);
        if (c != 0) return c;
        int pos1 = p1->start;
        int pos2 = p2->start;
        if (pos1 > pos2) return 1;
        if (pos1 < pos2) return -1;
        return 0;
    });
}

Doubletable *
PeakLibrary::countPeakTagsLowMemory(TagLibrary *tags, char direction, int mode, PeakLibrary *excludePeaks) {

    Doubletable *results = new Doubletable();
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        results->insert(0.0, keys[i]);
        delete[](keys[i]);
    }
    delete[]keys;

    keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (tags->singleFile) tags->readSingleTagFile();

    for (int i = 0; i < chrs->total; i++) {
        ChrPeaks *cp = (ChrPeaks *) chrs->search(keys[i]);
        ChrTags *ct = (ChrTags *) tags->chrs->search(keys[i]);
        ChrPeaks *cte = NULL;
        if (excludePeaks != NULL) cte = (ChrPeaks *) excludePeaks->chrs->search(keys[i]);
        delete[](keys[i]);
        if (ct == NULL || cp == NULL) {
            continue;
        }
        cp->countPeakTagsLowMemory(results, ct, direction, mode, cte);
    }
    delete[]keys;

    return results;
}

int PeakLibrary::addTagLibrary(TagLibrary *t) {
    TagLibrary **newexps = new TagLibrary *[numExps + 1];
    if (exps != NULL) {
        for (int i = 0; i < numExps; i++) {
            newexps[i] = exps[i];
        }
        delete[]exps;
    }
    newexps[numExps] = t;
    int rv = numExps;
    numExps++;
    exps = newexps;

    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->addExp();  // this forces a blank space to be added (linkedlist)
        delete[](keys[i]);
    }
    delete[]keys;
    keys = NULL;

    keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (t->singleFile) t->readSingleTagFile();
    for (int i = 0; i < chrs->total; i++) {
        ChrPeaks *cp = (ChrPeaks *) chrs->search(keys[i]);
        ChrTags *ct = (ChrTags *) t->chrs->search(keys[i]);
        delete[](keys[i]);
        if (ct == NULL || cp == NULL) {
            continue;
        }
        cp->addTagLibrary(ct, numExps - 1);
    }
    delete[]keys;

    sortPeakTags(numExps - 1); //here we optimize the linkedlist for access
    return rv;
}

void PeakLibrary::sortPeakTags(int expIndex) {
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->sortTags(expIndex);
        delete[](keys[i]);
    }
    delete[]keys;
}

void PeakLibrary::setPeakTagSizeRefPos(int newOffset, int startOffset, int endOffset) {
//fprintf(stderr, "%d\t%d\t%d\n", newOffset, startOffset, endOffset);
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->setPeakTagSizeRefPos(newOffset, startOffset, endOffset);
        delete[](keys[i]);
    }
    delete[]keys;
    sortChr();
}

void PeakLibrary::setPeakTagSizeFixed(int startOffset, int endOffset) {
    fixedFlag = 1;
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->setPeakTagSizeFixed(startOffset, endOffset);
        delete[](keys[i]);
    }
    delete[]keys;
    sortChr();
}

void PeakLibrary::setDefaultPeakOrder() {
    if (peakOrder != NULL) {
        delete[]peakOrder;
    }
    if (numPeaks < 1) return;
    peakOrder = new Peak *[numPeaks];
    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        peakOrder[i] = (Peak *) peaks->search(keys[i]);
        delete[](keys[i]);
    }
    delete[]keys;

    qsort(peakOrder, numPeaks, sizeof(Peak *), &cmpPeaks);
    for (int i = 0; i < numPeaks; i++) {
    }

}

PeakLibrary *PeakLibrary::filterClonalPeaks(TagLibrary *tags, int peakSize,
                                            double threshold, int mode, char strand) {

    char *outputstr = new char[10000];
    int halfPeakSize = (peakSize) / 2;

    float currentMaxTBP = tags->maxtbp;
    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfPeakSize, halfPeakSize);

    tags->setMaxTBP(0);
    //int tagsIndex = addTagLibrary(tags);
    Doubletable *expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    //Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);

    tags->setMaxTBP(1);
    //int posIndex = addTagLibrary(tags);
    Doubletable *posTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);

    tags->setMaxTBP(currentMaxTBP);

    double avgTagsPerPosition = tags->averageTagsPerPosition;
    if (avgTagsPerPosition < 1) avgTagsPerPosition = 1.0;

    //Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);
    //Doubletable* posTags = countPeakTags(posIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);

    int expectedSize = peakSize * 10;
    double *expected = new double[expectedSize];
    for (int i = 0; i < expectedSize; i++) {
        if (i == 0) {
            expected[i] = 0;
        } else if (i == 1) {
            expected[i] = 1;
        } else {
            expected[i] = expected[i - 1] + ((double) peakSize - expected[i - 1]) / ((double) peakSize);
        }
    }

    PeakLibrary *goodPeaks = new PeakLibrary();
    char **keys = expTags->keys();
    int numGood = 0;
    int totalChecked = 0;
    goodPeaks->tagsInPeaks = 0.0;
    for (int i = 0; i < peaks->total; i++) {
        double pt = expTags->search(keys[i]);
        double ct = posTags->search(keys[i]);
        if (pt < EMPTY_DOUBLE_CHECK || ct < EMPTY_DOUBLE_CHECK) {
            delete[](keys[i]);
            continue;
        }
        totalChecked++;

        int index = (int) (pt / avgTagsPerPosition);
        if (index > expectedSize - 1) index = expectedSize - 1;
        double et = expected[index];
        double fold = et / ct;
//fprintf(stderr, "%s\t%lf\t%lf\t%d\t%lf\n", keys[i], pt, ct, index,et);
        if (ct > 0 && fold < threshold) {
            goodPeaks->tagsInPeaks += pt;
            numGood++;
            Peak *p = (Peak *) peaks->search(keys[i]);
            sprintf(outputstr, "%.2lf", fold);
            p->addData(outputstr);
            goodPeaks->addPeak(p);
        }
        delete[](keys[i]);

    }
    delete[]keys;
    delete[]outputstr;
    delete posTags;
    delete expTags;
    if (totalChecked == 0) {
        fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
    } else {
        double ratio = ((double) numGood) / ((double) totalChecked);
        fprintf(stderr, "\tClonal filtering: %d of %d (%.2lf%% passed)\n", numGood, totalChecked, ratio * 100.0);
    }
    goodPeaks->sortChr();
    return goodPeaks;
}

PeakLibrary *PeakLibrary::filterLocalPeaks(TagLibrary *tags, int peakSize, int localSize,
                                           double foldThresh, double poissonThresh, int mode, char strand) {

    char *outputstr = new char[10000];

    int halfPeakSize = (peakSize) / 2;
    int halfLocalSize = (localSize) / 2;

    //int tagsIndex = addTagLibrary(tags);

    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfPeakSize, halfPeakSize);
    Doubletable *expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfLocalSize, halfLocalSize);
    Doubletable *localTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);

    //Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);
    //Doubletable* localTags = countPeakTags(tagsIndex,-1*halfLocalSize,halfLocalSize,strand,COUNT_MODE_TOTAL);

    double peakLength = (double) peakSize;
    double localLength = (double) (localSize - peakSize);
    double localAdjustFactor = peakLength / localLength;

    double tagPseudoCount = 0.0;
    if (tags->tbp > 0.0) tagPseudoCount = 0.5;
    //if (tags->tbp > 0.0) tagPseudoCount = localLength*tags->tbp;

    PeakLibrary *goodPeaks = new PeakLibrary();
    char **keys = peaks->keys();
    int numGood = 0;
    goodPeaks->tagsInPeaks = 0.0;
    int totalChecked = 0;
    for (int i = 0; i < peaks->total; i++) {
        double pt = expTags->search(keys[i]);
        double lt = localTags->search(keys[i]);
        if (pt < EMPTY_DOUBLE_CHECK || lt < EMPTY_DOUBLE_CHECK) {
            delete[](keys[i]);
            continue;
        }
        totalChecked++;
        lt -= pt;
        if (lt < tagPseudoCount) lt = tagPseudoCount;

        double ltadjusted = lt * localAdjustFactor;
        double fold = pt / ltadjusted;
        int foldGood = 0;

        if (foldThresh > 0.000001) {
            if (fold > foldThresh) {
                foldGood = 1;
            }
        } else {
            foldGood = 1;
        }

        int ptInt = (int) pt;
        //since the input is our "expected" meausrement, that becomes lambda
        //double cumPvalue = cumulativePoisson(ptInt, ltadjusted);
        //double pp = 1-cumPvalue;
        double pp = exp(ilogCumulativePoisson(ptInt, ltadjusted));
        //fprintf(stderr, "%lf\t%lf\n", pp, threshold);
        int poissonGood = 0;
        if (poissonThresh < 1.0) {
            if (pp < poissonThresh) {
                poissonGood = 1;
            }
        } else {
            poissonGood = 1;
        }

        if (foldGood == 1 && poissonGood == 1) {
            goodPeaks->tagsInPeaks += pt;
            numGood++;
            Peak *p = (Peak *) peaks->search(keys[i]);
            sprintf(outputstr, "%.2lf\t%.2le", fold, pp);
            p->addData(outputstr);
            goodPeaks->addPeak(p);
        }
        delete[](keys[i]);
    }
    delete[]keys;
    delete[]outputstr;
    delete localTags;
    delete expTags;
    if (totalChecked == 0) {
        fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
    } else {
        double ratio = ((double) numGood) / ((double) totalChecked);
        fprintf(stderr, "\tLocal Background Filtering: %d of %d (%.2lf%% passed)\n", numGood, totalChecked,
                ratio * 100.0);
    }
    goodPeaks->sortChr();
    return goodPeaks;
}


PeakLibrary *PeakLibrary::getDifferentialPeaks(TagLibrary *tags, TagLibrary *input,
                                               double foldThresh, double poissonThresh, int mode, int start, int end,
                                               char strand, int strFlag) {


    char *outputstr = new char[10000];
    //int tagsIndex = addTagLibrary(tags);
    //int inputIndex = addTagLibrary(input);
    //Doubletable* expTags = countPeakTags(tagsIndex,start,end,strand,COUNT_MODE_TOTAL);
    //Doubletable* inputTags = countPeakTags(inputIndex,start,end,strand,COUNT_MODE_TOTAL);
    Doubletable *expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    Doubletable *inputTags = countPeakTagsLowMemory(input, strand, COUNT_MODE_TOTAL, NULL);

    //double peakSize = ((double)(end-start));
    if (fixedFlag) {
        //peakSize = avgPeakSize;
    }
    //fprintf(stderr, "peakSize = %lf\n", peakSize);

    double inputPseudoCount = 0.0;
    double tagPseudoCount = 0.0;
    if (input->tbp > 0.0) inputPseudoCount = 0.5;
    if (tags->tbp > 0.0) tagPseudoCount = 0.5;
    //if (input->tbp > 0.0) inputPseudoCount = ((double)(peakSize))*input->tbp;
    //if (tags->tbp > 0.0) tagPseudoCount = ((double)(peakSize))*tags->tbp;
    //fprintf(stderr, "\tpeudo= %lf\n",  tagPseudoCount);

    double minTotal = tags->totalTags;
    if (input->totalTags < tags->totalTags) minTotal = input->totalTags;

    double tagsRatio = minTotal / tags->totalTags;
    double inputRatio = minTotal / input->totalTags;

    int numGood = 0;
    int totalChecked = 0;

    PeakLibrary *goodPeaks = new PeakLibrary();
    goodPeaks->tagsInPeaks = 0.0;
    char **keys = expTags->keys();
    for (int i = 0; i < expTags->total; i++) {
        double tp = expTags->search(keys[i]);
        double ip = inputTags->search(keys[i]);
        if (tp < EMPTY_DOUBLE_CHECK || ip < EMPTY_DOUBLE_CHECK) {
            delete[](keys[i]);
            continue;
        }

        totalChecked++;
        if (ip < inputPseudoCount) ip = inputPseudoCount;
        if (tp < tagPseudoCount) tp = tagPseudoCount;
        double ipn = ip * inputRatio;
        double tpn = tp * tagsRatio;
        double foldchange = tpn;
        if (ipn > 0) {
            foldchange = tpn / ipn;
        }
        int foldGood = 0;
        if (foldThresh > 0.000001) {
            if (mode == DIFFPEAK_MODE_DIFF) {
                if (foldchange > foldThresh) {
                    foldGood = 1;
                }
            } else if (mode == DIFFPEAK_MODE_SAME) {
                if (foldchange < foldThresh && foldchange > 1.0 / foldThresh) {
                    foldGood = 1;
                }
            } else if (mode == DIFFPEAK_MODE_REV) {
                if (foldchange < 1.0 / foldThresh) {
                    foldGood = 1;
                }
            }
        } else {
            foldGood = 1;
        }

        int tpInt = (int) tpn;
        //since the input is our "expected" meausrement, that becomes lambda
        //double cumPvalue = cumulativePoisson(tpInt, ipn);
        //double pp = 1-cumPvalue;
        double pp = exp(ilogCumulativePoisson(tpInt, ipn));
        //fprintf(stderr, "%.1lf\t%.1lf\t%lf\n", tpn, ipn, pp);
        int poissonGood = 0;
        if (poissonThresh < 1.0) {
            if (mode == DIFFPEAK_MODE_DIFF) {
                if (pp < poissonThresh && foldchange > 1) {
                    poissonGood = 1;
                }
            } else if (mode == DIFFPEAK_MODE_SAME) {
                int ipInt = (int) ipn;
                //double cumPvalueR = cumulativePoisson(ipInt, tpn);
                //double ppR = 1-cumPvalue;
                double ppR = exp(ilogCumulativePoisson(ipInt, tpn));
                if (ppR > poissonThresh && pp > poissonThresh) {
                    poissonGood = 1;
                }
            } else if (mode == DIFFPEAK_MODE_REV) {
                int ipInt = (int) ipn;
                //double cumPvalueR = cumulativePoisson(ipInt, tpn);
                //double ppR = 1-cumPvalue;
                double ppR = exp(ilogCumulativePoisson(ipInt, tpn));
                if (ppR < poissonThresh && foldchange < 1) {
                    poissonGood = 1;
                    pp = ppR;
                }
            }
        } else {
            poissonGood = 1;
        }

        if (poissonGood == 1 && foldGood == 1) {
            goodPeaks->tagsInPeaks += tp;
            Peak *p = (Peak *) peaks->search(keys[i]);
            if (strFlag) {
                sprintf(outputstr, "%.1lf\t%.1lf\t%.2lf\t%.2le", tpn, ipn, foldchange, pp);
                p->addData(outputstr);
            }
            goodPeaks->addPeak(p);
            numGood++;
        }
        delete[](keys[i]);
    }
    delete[]keys;
    delete[]outputstr;
    delete expTags;
    delete inputTags;

    if (totalChecked == 0) {
        fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
    } else {
        double ratio = ((double) numGood) / ((double) totalChecked);
        fprintf(stderr, "\tDifferential Peaks: %d of %d (%.2lf%% passed)\n", numGood, totalChecked, ratio * 100.0);
    }
    goodPeaks->sortChr();
    //fprintf(stderr, "totalTagsInPeaks=%lf\n", goodPeaks->tagsInPeaks);
    return goodPeaks;

}

PeakLibrary *PeakLibrary::stitchRegions(int maxDistance, int mode) const {
    PeakLibrary *regions = new PeakLibrary();
    //regions->tagsInPeaks = 0.0;
    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrPeaks *cp = (ChrPeaks *) chrs->search(keys[i]);
        cp->stitchRegions(regions, maxDistance, mode);
        delete[](keys[i]);
    }
    delete[]keys;
    regions->sortChr();
    return regions;
}

void PeakLibrary::centerPeaks(TagLibrary *tags, int peakSize, char strand) {

    fprintf(stderr, "\tCentering peaks of size %d using a fragment length of %d\n",
            peakSize, tags->fragmentLengthEstimate);
    int fragLength = tags->fragmentLengthEstimate;
    int halfPeakSize = (int) (peakSize / 2.0);

    tags->setTagAdjust(0);
    if (peakSize == 0) {
        setPeakTagSizeFixed(-halfPeakSize - fragLength, halfPeakSize + fragLength);
    } else {
        setPeakTagSizeRefPos(NULL_OFFSET, -halfPeakSize - fragLength, halfPeakSize + fragLength);
    }
    int expIndex = addTagLibrary(tags);

    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->centerPeak(expIndex, fragLength, strand);
        delete[](keys[i]);
    }
    delete[]keys;
}

void PeakLibrary::centerNFR(TagLibrary *tags, int peakSize, char strand, int nfrSize) {

    fprintf(stderr, "\tCentering peaks on Nucleosome Free Region of size %d using a fragment length of %d\n",
            nfrSize, tags->fragmentLengthEstimate);
    int fragLength = tags->fragmentLengthEstimate;
    int halfPeakSize = (int) ((nfrSize + peakSize) / 2.0);
    int nucSize = 150;

    tags->setTagAdjust(0);
    if (peakSize == 0) {
        setPeakTagSizeFixed(-1 * halfPeakSize - fragLength - nucSize, halfPeakSize + fragLength + nucSize);
    } else {
        setPeakTagSizeRefPos(NULL_OFFSET, -halfPeakSize - fragLength - nucSize, halfPeakSize + fragLength + nucSize);
    }
    int expIndex = addTagLibrary(tags);

    char **keys = peaks->keys();
    for (int i = 0; i < peaks->total; i++) {
        Peak *p = (Peak *) peaks->search(keys[i]);
        p->centerNFR(expIndex, fragLength, strand, nfrSize, nucSize);
        delete[](keys[i]);
    }
    delete[]keys;
}

// class ChrPeaks -----------------------------------------------

ChrPeaks::ChrPeaks() {
    //peaks = new Peak*[PEAK_INC];
    peaks = NULL;
    peakList = NULL;
    numPeaks = 0;
}

ChrPeaks::~ChrPeaks() {
    if (peaks != NULL) delete[]peaks;
    if (peakList != NULL) delete peakList;
}

void ChrPeaks::addPeak(Peak *p) {
    if (peakList == NULL) {
        peakList = new LinkedList();
    }
    peakList->add(p);
}

void ChrPeaks::sort() {

    int numListPeaks = 0;
    Peak **array = NULL;
    if (peakList != NULL) {
        array = (Peak **) peakList->toArray(numListPeaks);
    }

    int newNumPeaks = numPeaks + numListPeaks;
    Peak **newpeaks = new Peak *[newNumPeaks];
    for (int i = 0; i < numPeaks; i++) {
        newpeaks[i] = peaks[i];
    }
    for (int i = 0; i < numListPeaks; i++) {
        newpeaks[i + numPeaks] = array[i];
    }
    delete[]peaks;
    delete[]array;
    delete peakList;
    peakList = NULL;

    peaks = newpeaks;
    numPeaks = newNumPeaks;
    if (numPeaks < 2) return;
    qsort(peaks, numPeaks, sizeof(Peak *), &cmpPeaks);
}

//global clusterF**k variable
static Inttable *chrIndex = NULL;

int chrcmp(const void *chr1, const void *chr2) {
    if (chr1 == chr2) return 0;
    if (chr1 == NULL) return -1;
    if (chr2 == NULL) return 1;
    char *c1 = *((char **) chr1);
    char *c2 = *((char **) chr2);
    int sc = strcmp(c1, c2);
    if (sc == 0) return 0;

    if (chrIndex == NULL) {
        chrIndex = new Inttable(10000);
        char *tmp = new char[1000];
        int index = 1;
        for (int i = 0; i <= 100; i++) {
            sprintf(tmp, "chr%d", i);
            chrIndex->insert(index++, tmp);
            sprintf(tmp, "chr%d_random", i);
            chrIndex->insert(index++, tmp);
            sprintf(tmp, "chr%dL", i);
            chrIndex->insert(index++, tmp);
            sprintf(tmp, "chr%dL_random", i);
            chrIndex->insert(index++, tmp);
            sprintf(tmp, "chr%dR", i);
            chrIndex->insert(index++, tmp);
            sprintf(tmp, "chr%dR_random", i);
            chrIndex->insert(index++, tmp);
        }
        chrIndex->insert(index++, (char *) "chrI");
        chrIndex->insert(index++, (char *) "chrII");
        chrIndex->insert(index++, (char *) "chrIII");
        chrIndex->insert(index++, (char *) "chrIV");
        chrIndex->insert(index++, (char *) "chrV");
        chrIndex->insert(index++, (char *) "chrVI");
        chrIndex->insert(index++, (char *) "chrVII");
        chrIndex->insert(index++, (char *) "chrVIII");
        chrIndex->insert(index++, (char *) "chrIX");
        chrIndex->insert(index++, (char *) "chrU");
        chrIndex->insert(index++, (char *) "chrU_random");
        chrIndex->insert(index++, (char *) "chrX");
        chrIndex->insert(index++, (char *) "chrX_random");
        chrIndex->insert(index++, (char *) "chrXI");
        chrIndex->insert(index++, (char *) "chrXII");
        chrIndex->insert(index++, (char *) "chrXIII");
        chrIndex->insert(index++, (char *) "chrXIV");
        chrIndex->insert(index++, (char *) "chrXV");
        chrIndex->insert(index++, (char *) "chrXVI");
        chrIndex->insert(index++, (char *) "chrXVII");
        chrIndex->insert(index++, (char *) "chrXVIII");
        chrIndex->insert(index++, (char *) "chrXIX");
        chrIndex->insert(index++, (char *) "chrY");
        chrIndex->insert(index++, (char *) "chrY_random");
        chrIndex->insert(index++, (char *) "chrZ");
        chrIndex->insert(index++, (char *) "chrZ_random");
        chrIndex->insert(index++, (char *) "chrM");
        chrIndex->insert(index++, (char *) "chrM_random");
        chrIndex->insert(index++, (char *) "chrMT");
        chrIndex->insert(index++, (char *) "chrMT_random");
        chrIndex->insert(index++, (char *) "chrUn");
        chrIndex->insert(index++, (char *) "chrUn_random");
        chrIndex->insert(index++, (char *) "genome");
        chrIndex->insert(index++, (char *) "null");
        delete[]tmp;
    }
    int i1 = chrIndex->search(c1);
    int i2 = chrIndex->search(c2);
    //if (i1 == EMPTY_INT) fprintf(stderr, "Couldn't find anything for %s\n", c1);
    //if (i2 == EMPTY_INT) fprintf(stderr, "Couldn't find anything for %s\n", c2);
    if (i1 != EMPTY_INT && i2 != EMPTY_INT) {
        if (i1 < i2) return -1;
        if (i1 > i2) return 1;
        return 0;
    } else if (i1 != EMPTY_INT) {
        return -1;
    } else if (i2 != EMPTY_INT) {
        return 1;
    }
    return sc;
}

int cmpPeaks(const void *a, const void *b) {
    char *ac = (*(Peak **) a)->chr;
    char *bc = (*(Peak **) b)->chr;
    int cc = chrcmp(&ac, &bc);
    if (cc != 0) return cc;
    int ap = (*(Peak **) a)->tagStart;
    int bp = (*(Peak **) b)->tagStart;
    if (ap < bp) return -1;
    if (ap > bp) return 1;
    int ad = (*(Peak **) a)->tagEnd;
    int bd = (*(Peak **) b)->tagEnd;
    if (ad < bd) return -1;
    if (ad > bd) return 1;
    char al = (*(Peak **) a)->strand;
    char bl = (*(Peak **) b)->strand;
    if (al < bl) return -1;
    if (al > bl) return 1;
    return 0;
}

void ChrPeaks::countPeakTagsLowMemory(Doubletable *results, ChrTags *ct, char direction, int mode,
                                      ChrPeaks *excludePeaks) {

    ct->loadTags();

    int peakIndex = 0;
    //establish when we can start forgetting about peaks
    // peaks are already sorted by their starting positions;
    int *finishedLength = new int[numPeaks];
    int *excludeFinishedLength = NULL;
    int excludePeakIndex = 0;
    if (excludePeaks != NULL) {
        excludeFinishedLength = new int[excludePeaks->numPeaks];
    }
    double *peakTotals = new double[numPeaks];
    double *peakPositions = new double[numPeaks];
    for (int i = 0; i < numPeaks; i++) {
        peakTotals[i] = 0.0;
        peakPositions[i] = 0.0;
        if (i == 0) {
            finishedLength[i] = peaks[i]->tagEnd;
        } else {
            if (finishedLength[i - 1] > peaks[i]->tagEnd) {
                finishedLength[i] = finishedLength[i - 1];
            } else {
                finishedLength[i] = peaks[i]->tagEnd;
            }
        }
    }

    if (excludePeaks != NULL) {
        for (int i = 0; i < excludePeaks->numPeaks; i++) {
            if (i == 0) {
                excludeFinishedLength[i] = excludePeaks->peaks[i]->tagEnd;
            } else {
                if (excludeFinishedLength[i - 1] > excludePeaks->peaks[i]->tagEnd) {
                    excludeFinishedLength[i] = excludeFinishedLength[i - 1];
                } else {
                    excludeFinishedLength[i] = excludePeaks->peaks[i]->tagEnd;
                }
            }
        }
    }

    for (int i = 0; i < ct->totalPositions; i++) {
        int p = ct->tags[i].p;
        int d = ct->tags[i].d;
        float v = ct->tags[i].v;

        //if excludePeaks are used, don't bother counting read if it overlaps with the exclude peaks
        if (excludePeaks != NULL && excludePeakIndex < excludePeaks->numPeaks) {
            if (p < excludePeaks->peaks[excludePeakIndex]->tagStart) {
            } else {
                int stop = 0;
                while (p > excludeFinishedLength[excludePeakIndex]) {
                    excludePeakIndex++;
                    if (excludePeakIndex >= excludePeaks->numPeaks) {
                        stop = 1;
                        break;
                    }
                }
                if (stop == 0) {
                    int bad = 0;
                    for (int j = excludePeakIndex; j < excludePeaks->numPeaks; j++) {
                        if (direction == STRAND_BOTH) {
                        } else if (direction == STRAND_SEPARATE || direction == POSITIVE_STRAND) {
                            if (excludePeaks->peaks[j]->strand != d) {
                                continue;
                            }
                        } else if (direction == NEGATIVE_STRAND) {
                            if (excludePeaks->peaks[j]->strand == d) {
                                continue;
                            }
                        }
                        if (p >= excludePeaks->peaks[j]->tagStart && p <= excludePeaks->peaks[j]->tagEnd) {
                            bad = 1;
                            break;
                        }
                        if (p <= excludePeaks->peaks[j]->tagStart) {
                            break;
                        }
                    }
                    //means we're overlapping an excludePeak so continue the loop
                    if (bad) continue;
                }
            }
        }


        if (p < peaks[peakIndex]->tagStart) continue;
        int stop = 0;
        while (p > finishedLength[peakIndex]) {
            peakIndex++;
            if (peakIndex >= numPeaks) {
                stop = 1;
                break;
            }
        }
        if (stop) break;
        for (int j = peakIndex; j < numPeaks; j++) {
            if (direction == STRAND_BOTH) {
            } else if (direction == STRAND_SEPARATE || direction == POSITIVE_STRAND) {
                if (peaks[j]->strand != d) {
                    continue;
                }
            } else if (direction == NEGATIVE_STRAND) {
                if (peaks[j]->strand == d) {
                    continue;
                }
            }
            if (p >= peaks[j]->tagStart && p <= peaks[j]->tagEnd) {
                peakTotals[j] += v;
                peakPositions[j] += 1.0;
            }
            if (p <= peaks[j]->tagStart) {
                break;
            }
        }
    }

    for (int i = 0; i < numPeaks; i++) {

        double total = peakTotals[i];
        double totalPosUsed = peakPositions[i];
        double totalBP = (double) peaks[i]->end - peaks[i]->start + 1;

        if (mode == COUNT_MODE_TBP) {
            if (totalBP > 0) total /= totalBP;
        }
        if (mode == COUNT_MODE_RATIO) {
            if (totalPosUsed > 0.0) {
                total /= totalPosUsed;
            } else {
                total = EMPTY_DOUBLE;
            }
        }
        results->insert(total, peaks[i]->name);
    }
    delete[]finishedLength;
    delete[]peakTotals;
    delete[]peakPositions;
    if (excludeFinishedLength != NULL) delete[]excludeFinishedLength;
    ct->freeTags();

}

void ChrPeaks::addTagLibrary(ChrTags *ct, int expIndex) {
    ct->loadTags();

    int peakIndex = 0;
    //int tagIndex = 0;
    //establish when we can start forgetting about peaks
    // peaks are already sorted by their starting positions;
    int *finishedLength = new int[numPeaks];
    for (int i = 0; i < numPeaks; i++) {
        if (i == 0) {
            finishedLength[i] = peaks[i]->tagEnd;
        } else {
            if (finishedLength[i - 1] > peaks[i]->tagEnd) {
                finishedLength[i] = finishedLength[i - 1];
            } else {
                finishedLength[i] = peaks[i]->tagEnd;
            }
        }
    }

    //fprintf(stderr, "TotalPeaks=%d, totalTags=%d\n",numPeaks,ct->totalPositions);

    for (int i = 0; i < ct->totalPositions; i++) {
        int p = ct->tags[i].p;

        if (p < peaks[peakIndex]->tagStart) continue;
        int stop = 0;
        while (p > finishedLength[peakIndex]) {
            peakIndex++;
            if (peakIndex >= numPeaks) {
                stop = 1;
                break;
            }
        }
        if (stop) break;
        for (int j = peakIndex; j < numPeaks; j++) {
            if (p >= peaks[j]->tagStart && p <= peaks[j]->tagEnd) {
                peaks[j]->addTag(&(ct->tags[i]), expIndex);
            }
            if (p <= peaks[j]->tagStart) {
                break;
            }
        }
    }
    delete[]finishedLength;

    ct->freeTags();
}

void ChrPeaks::stitchRegions(PeakLibrary *regions, int maxDistance, int mode) {

    double *totals = new double[numPeaks];
    int *totalDist = new int[numPeaks];
    int *uniqmap = new int[numPeaks];
    for (int i = 0; i < numPeaks; i++) totals[i] = 0.0;
    for (int i = 0; i < numPeaks; i++) totalDist[i] = 0;
    for (int i = 0; i < numPeaks; i++) uniqmap[i] = 0;
    int totalIndex = 0;
    int maxNumEstDist = 500;

    for (int k = 0; k < 2; k++) {

        for (int index = 0; index < numPeaks; index++) {
            int i = index;
            int peakLimit = numPeaks;
            if (k == 1) {
                i = numPeaks - index - 1;
                peakLimit = i + i + 1;
            }

            if (peaks[i]->strand != k) continue;
            int curStart = peaks[i]->start;
            int curEnd = peaks[i]->end;
            //char curStrand = peaks[i]->strand;
            float v = peaks[i]->v;
            int lastAdded = i;
            totalIndex = 0;
            totals[totalIndex] = v;
            uniqmap[totalIndex] = peaks[i]->uniqMap;
            if (uniqmap[totalIndex] < 0) uniqmap[totalIndex] = peaks[i]->end - peaks[i]->start;
            if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;
            totalDist[totalIndex] = peaks[i]->refPos;

            for (int indexj = i + 1; indexj < peakLimit; indexj++) {
                int j = indexj;
                if (k == 1) {
                    j = i - (indexj - i);
                }
                if (peaks[j]->strand != k) continue;
                if (peaks[j]->refPos - peaks[lastAdded]->refPos < maxDistance
                    && peaks[j]->refPos - peaks[lastAdded]->refPos > -1 * maxDistance) {

                    int add = 0;
                    if (mode == REGION_MODE_GROSEQ) {
                        double score = 0.0;
                        double NN = 0.0;
                        //int diff = 0;
                        int dindex = totalIndex;
                        while (dindex >= 0 && abs(totalDist[totalIndex] - totalDist[dindex])
                                              < maxNumEstDist) {
                            //score+=totals[dindex];
                            double curScore = totals[dindex] / ((double) uniqmap[dindex]);
//							if (curScore > score) {
//								score = curScore;
//							}
                            score += curScore;
                            NN += 1.0;
                            //diff = totalDist[totalIndex]-totalDist[dindex];
                            dindex--;
                        }
                        if (NN > 0.0) {
                            //score /= NN;
                            score /= NN;
                        } else {
                            score = 1;
                        }

                        double pmap = (double) peaks[j]->uniqMap;
                        if (pmap < 0) pmap = (double) (peaks[j]->end - peaks[j]->start);
                        if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;

                        double fold = (peaks[j]->v / pmap) / (score);

                        //fprintf(stderr, "%d\t%lf\t%lf\t%lf\t%lf\n", peaks[j]->refPos,fold,peaks[j]->v,pmap,score);
                        if (fold > REGION_GROSEQ_FOLDCHANGE) {
                            break;
                        } else {
                            add = 1;
                        }
                    } else if (mode == REGION_MODE_HISTONE) {
                        add = 1;
                    }

                    if (add) {
                        totalIndex++;
                        if (k == 1) {
                            curStart = peaks[j]->start;
                        } else {
                            curEnd = peaks[j]->end;
                        }
                        totalDist[totalIndex] = peaks[j]->refPos;
                        totals[totalIndex] = peaks[j]->v;
                        uniqmap[totalIndex] = peaks[j]->uniqMap;
                        if (uniqmap[totalIndex] < 0) uniqmap[totalIndex] = peaks[j]->end - peaks[j]->start;
                        if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;
                        v += peaks[j]->v;
                        lastAdded = j;
                    }
                } else {
                    break;
                }
            }
            peaks[i]->start = curStart;
            peaks[i]->end = curEnd;
            //fprintf(stderr, "done = %d\t%d\n", curStart, curEnd);
            peaks[i]->focusRatio = curEnd - curStart;
            peaks[i]->v = v;
            regions->addPeak(peaks[i]);
            index = lastAdded;
            if (k == 1) {
                index = numPeaks - lastAdded - 1;
            }
        }
    }
    delete[]totals;
    delete[]totalDist;
}

// class Peaks -----------------------------------------------

Peak::Peak() {
    name = NULL;
    chr = NULL;
    refPos = 0;
    start = 0;
    end = 0;
    tagStart = 0;
    tagEnd = 0;
    priority = 0;
    strand = STRAND_POSITIVE;
    seq = NULL;
    uniqMap = 0;
    v = 0.0;
    focusRatio = 0.0;
    data = NULL;
    ogname = NULL;
    exps = NULL;
    numTags = NULL;
    maxTags = NULL;
    numExps = 0;
}

Peak::Peak(char *newname, char *originalName, char *newchr, int newstart, int newend, int newRef, char dir,
           float value, float ratio, char *otherdata, int mappability, unsigned int newpriority) {
    if (newname == NULL) {
        int L = strlen(newchr) + 13 + 1 + 6;
        name = new char[L];
        sprintf(name, "%s-%d-%d", newchr, newstart, dir);
        //fprintf(stderr, "name=%s \t %d\n",name,L);
    } else {
        name = new char[strlen(newname) + 1];
        strcpy(name, newname);
    }
    ogname = NULL;
    if (originalName != NULL) {
        ogname = new char[strlen(originalName) + 1];
        strcpy(ogname, originalName);
    }
    chr = new char[strlen(newchr) + 1];
    strcpy(chr, newchr);
    start = newstart;
    end = newend;
    tagStart = newstart;
    tagEnd = newend;
    refPos = newRef;
    strand = dir;
    v = value;
    seq = NULL;
    priority = newpriority;
    focusRatio = ratio;
    uniqMap = mappability;
    exps = NULL;
    numTags = NULL;
    maxTags = NULL;
    numExps = 0;
    if (otherdata != NULL) {
        data = new char[strlen(otherdata) + 1];
        strcpy(data, otherdata);
    } else {
        data = NULL;
    }
    if (refPos == NULL_REF) {
        refPos = (start + end) / 2;
    }
}

Peak::~Peak() {
    if (name != NULL) delete[]name;
    if (ogname != NULL) delete[]ogname;
    if (chr != NULL) delete[]chr;
    if (data != NULL) delete[]data;
    if (exps != NULL) {
        for (int i = 0; i < numExps; i++) {
            delete exps[i];
        }
        delete[]exps;
    }
    if (numTags != NULL) delete[]numTags;
    if (maxTags != NULL) delete[]maxTags;
    if (seq != NULL) delete[]seq;
}

void Peak::print(FILE *fp) {
    char dir = '+';
    if (strand == 1) dir = '-';
    if (v < 1.0) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.3f\t%.3f", name, chr, start, end, dir, v, focusRatio);
    } else if (v < 10.0) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.2f\t%.3f", name, chr, start, end, dir, v, focusRatio);
    } else {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.1f\t%.3f", name, chr, start, end, dir, v, focusRatio);
    }
    //fprintf(fp, "\t%d", uniqMap);
    if (data != NULL) {
        fprintf(fp, "\t%s", data);
    }
    fprintf(fp, "\n");
}

void Peak::addExp() {
    Tag **newexps = new Tag *[numExps + 1];
    int *newnumTags = new int[numExps + 1];
    if (exps != NULL) {
        for (int i = 0; i < numExps; i++) {
            newexps[i] = exps[i];
            newnumTags[i] = numTags[i];
        }
        delete[]exps;
        delete[]numTags;
    }

    //we're going to hide this in the variable for now...
    LinkedList *linkedlist = new LinkedList();
    newexps[numExps] = (Tag *) linkedlist;
    newnumTags[numExps] = 0;

    numTags = newnumTags;
    exps = newexps;

    numExps++;
}

void Peak::addData(char *str) {
    if (str == NULL) return;
    char *olddata = data;
    int len = 0;
    if (olddata != NULL) len = strlen(olddata) + 1;
    data = new char[len + strlen(str) + 2];
    data[0] = '\0';
    if (olddata != NULL) {
        strcpy(data, olddata);
        strcat(data, "\t");
        delete[]olddata;
    }
    strcat(data, str);
}

void Peak::addTag(Tag *t, int expIndex) {

    Tag *nt = new Tag();
    nt->copy(t);
    if (strand == 0) {
        nt->p = nt->p - refPos;
    } else {
        nt->p = refPos - nt->p;
        if (nt->d == 0) {
            nt->d = 1;
        } else {
            nt->d = 0;
        }
    }
    LinkedList *linkedlist = (LinkedList *) exps[expIndex];
    linkedlist->add(nt);
    //fprintf(stderr, "+ %.1f ",nt->v);
    //numTags[expIndex]++;
}

void Peak::setPeakTagSizeFixed(int startOffset, int endOffset) {
    setOffset(0);
    if (strand == 0) {
        tagStart = start + startOffset;
        tagEnd = end + endOffset;
        refPos = start;
    } else {
        tagStart = start - endOffset;
        tagEnd = end - startOffset;
        refPos = end;
    }
}

void Peak::setPeakTagSizeRefPos(int newOffset, int startOffset, int endOffset) {
    setOffset(newOffset);
    if (strand == 0) {
        tagStart = refPos + startOffset;
        tagEnd = refPos + endOffset;
    } else {
        tagStart = refPos - endOffset;
        tagEnd = refPos - startOffset;
    }
}

void Peak::setOffset(int newoffset) {
    if (newoffset == NULL_OFFSET) {
        //refPos = (int)floor(((float)(start+end))/2.0);
        refPos = (start + end) / 2;
        return;
    }
    if (strand == 0) {
        refPos = start - newoffset;
    } else {
        refPos = end + newoffset;
    }
}

void Peak::sortTags(int expIndex) {

    LinkedList *list = (LinkedList *) exps[expIndex];
    Tag **array = (Tag **) list->toArray(numTags[expIndex]);
    Tag *tagset = new Tag[numTags[expIndex]];
    delete list;
    exps[expIndex] = tagset;
    for (int i = 0; i < numTags[expIndex]; i++) {
        exps[expIndex][i].copy(array[i]);
        delete array[i];
    }
    delete[]array;

    for (int i = 0; i < numExps; i++) {
        if (i == expIndex || expIndex == ALL_PEAK_EXPS) {
            qsort(exps[i], numTags[i], sizeof(Tag), &cmpTags);
        }
    }
}

void Peak::centerPeak(int expIndex, int fragLength, char strandInfo) {

    int numCovTags = 0;
    Tag *coverageTags = getCoverageTags(expIndex, fragLength, numCovTags, strandInfo);
    if (numCovTags < 1) {
        focusRatio = 0.0;
        return;
    }


    float centerSum = 0;
    float centerN = 0;

    float max = 0;
    for (int i = 0; i < numCovTags; i++) {
        //if (coverageTags[i].v >= max) {
        if (coverageTags[i].v > max) {
            int diff = 1;
            int p = coverageTags[i].p;
            if (i < numCovTags - 1) {
                diff = coverageTags[i + 1].p - coverageTags[i].p;
                p += diff / 2;
            }
            if (max < coverageTags[i].v) {
                max = coverageTags[i].v;
                centerSum = diff * p;
                centerN = diff;
            } else {
                centerSum += diff * p;
                centerN += diff;
            }
        }
    }
    delete[]coverageTags;

    int centerP = refPos;
    int offset = 0;
    if (centerN == 0) {
        //fprintf(stderr, "Prolem\n");
    } else {
        centerP = (int) (((double) centerSum) / ((double) centerN));
        offset = centerP;
        if (strand == STRAND_POSITIVE) {
            refPos += offset;
        } else {
            refPos -= offset;
        }
        int halfPeakSize = (int) ((end - start) / 2.0);
        start = refPos - halfPeakSize;
        end = refPos + halfPeakSize;
    }

    double totalTags = 1.0;
    double goodTags = 0.0;
    for (int i = 0; i < numTags[expIndex]; i++) {
        if (exps[expIndex][i].d == 0) {
            if (exps[expIndex][i].p - PEAKRATIO_TOLERANCE_BP < offset) {
                goodTags += exps[expIndex][i].v;
            }
        } else {
            if (exps[expIndex][i].p + PEAKRATIO_TOLERANCE_BP > offset) {
                goodTags += exps[expIndex][i].v;
            }
        }
        totalTags += exps[expIndex][i].v;
    }

    if (totalTags > 0) {
        focusRatio = goodTags / totalTags;
    } else {
        focusRatio = 0.0;
    }

}

void Peak::centerNFR(int expIndex, int fragLength, char strandInfo, int nfrSize, int nucSize) {
    int numCovTags = 0;
    Tag *coverageTags = getCoverageTags(expIndex, fragLength, numCovTags, strandInfo);

    if (numCovTags < 1) {
        focusRatio = 0.0;
        return;
    }

    int nfrSizeHalf = nfrSize / 2;
    int covOffset = -1 * (nucSize + nfrSizeHalf) + coverageTags[0].p;
    int covEnd = coverageTags[numCovTags - 1].p - covOffset + (nucSize + nfrSizeHalf);

    double *cov = new double[covEnd + 1];
    int index = 0;
    for (int i = 0; i <= covEnd; i++) {
        int p = i + covOffset;
        if (index > numCovTags - 1) {
            cov[i] = 0;
            continue;
        }
        if (p < coverageTags[0].p) {
            cov[i] = 0.0;
        } else if (p < coverageTags[index + 1].p) {
            cov[i] = coverageTags[index].v;
        } else {
            index++;
            if (index > numCovTags - 1) {
                cov[i] = 0.0;
            } else {
                cov[i] = coverageTags[index].v;
            }
        }
    }


    //initialize Scores
    double nuc1Score = 0;
    double nuc2Score = 0;
    double nfrScore = 0;
    for (int i = 0; i < nucSize; i++) nuc1Score += cov[i];
    for (int i = nucSize; i < nucSize + nfrSize; i++) nfrScore += cov[i];
    for (int i = nucSize + nfrSize; i < 2 * nucSize + nfrSize; i++) nuc2Score += cov[i];

    int nuc1edge1Offset = -nucSize - nfrSizeHalf;
    int nuc1edge2Offset = -nfrSizeHalf;
    int nuc2edge1Offset = nfrSizeHalf;
    int nuc2edge2Offset = nucSize + nfrSizeHalf;

    double centerSum = 0;
    double centerN = 0;
    double max = -1e10;
    for (int i = nucSize + nfrSizeHalf; i < covEnd - (nucSize + nfrSizeHalf); i++) {
        nuc1Score -= cov[i + nuc1edge1Offset];
        nuc1Score += cov[i + nuc1edge2Offset];
        nuc2Score -= cov[i + nuc2edge1Offset];
        nuc2Score += cov[i + nuc2edge2Offset];
        nfrScore -= cov[i + nuc1edge2Offset];
        nfrScore += cov[i + nuc2edge1Offset];
        double ss = (nuc1Score + nuc2Score) / (2 * (double) nucSize) - nfrScore / ((double) nfrSize);

        if (ss >= max) {
            int p = i + covOffset;
            if (max < ss) {
                max = ss;
                centerSum = (double) p;
                centerN = 1.0;
            } else {
                centerSum += (double) p;
                centerN += 1.0;
            }
        }
//		fprintf(stderr, "%d\t%lf\t%lf\n", i+covOffset, cov[i],ss);
    }
    focusRatio = max;
    delete[]coverageTags;
    delete[]cov;

    int centerP = refPos;
    int offset = 0;
    if (centerN < 0.5) {
        //fprintf(stderr, "Prolem\n");
    } else {
        centerP = (int) (((double) centerSum) / ((double) centerN));
        offset = centerP;
        if (strand == STRAND_POSITIVE) {
            refPos += offset;
        } else {
            refPos -= offset;
        }
        int halfPeakSize = (int) ((end - start) / 2.0);
        start = refPos - halfPeakSize;
        end = refPos + halfPeakSize;
    }
//print(stderr);
//exit(0);
}

Tag *Peak::getCoverageTags(int expIndex, int fragLength, int &coveragePositions, char strandInfo) {

    coveragePositions = numTags[expIndex] * 2;
    Tag *coverageTags = new Tag[coveragePositions];

    int cIndex = 0;
    for (int i = 0; i < numTags[expIndex]; i++) {
        Tag *t = &(exps[expIndex][i]);
        if (strandInfo == STRAND_SEPARATE && t->d == 1) {
            //fprintf(stderr, "skipping...\n");
            continue;
        }
        if (strandInfo == POSITIVE_STRAND && t->d == 1) {
            //fprintf(stderr, "skipping...\n");
            continue;
        }
        if (strandInfo == NEGATIVE_STRAND && t->d == 0) {
            //fprintf(stderr, "skipping...\n");
            continue;
        }
        coverageTags[cIndex].v = t->v;
        coverageTags[cIndex].p = t->p;
        coverageTags[cIndex].d = 0;
        coverageTags[cIndex].len = t->len;
        cIndex++;
        coverageTags[cIndex].d = 0;
        coverageTags[cIndex].len = t->len;
        if (t->d == 0) {
            coverageTags[cIndex].p = t->p + fragLength;
            coverageTags[cIndex].v = -1 * t->v;
        } else {
            coverageTags[cIndex].p = t->p - fragLength;
            coverageTags[cIndex - 1].v = -1 * t->v;
            coverageTags[cIndex].v = t->v;
        }
        //fprintf(stderr, "%d\n", coverageTags[cIndex].p+refPos);
        cIndex++;
    }
    coveragePositions = cIndex;
    if (coveragePositions < 1) {
        delete[]coverageTags;
        return NULL;
    }

    qsort(coverageTags, coveragePositions, sizeof(Tag), &cmpTags);

    int last = 0;
    for (int i = 1; i < coveragePositions; i++) {
        if (coverageTags[i].p == coverageTags[last].p) {
            coverageTags[last].v += coverageTags[i].v;
        } else {
            last++;
            if (last == i) continue;
            coverageTags[last].copy(&(coverageTags[i]));
        }
    }
    coveragePositions = last + 1;

    float total = 0.0;
    for (int i = 0; i < coveragePositions; i++) {
        total += coverageTags[i].v;
        coverageTags[i].v = total;
    }

    return coverageTags;
}

//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################

TagLibrary::TagLibrary(char *dir) {
    directory = NULL;
    if (dir != NULL) {
        directory = new char[strlen(dir) + 1];
        strcpy(directory, dir);
    }
    chrs = new Hashtable(100000);
    chrNames = NULL;
    name = NULL;
    totalTags = 0;
    parseAlignmentCpGMinValue = 0.0;
    mCflag = 0;
    sspeFlag = 0;
    totalPositions = 0;
    averageTagsPerPosition = 0.0;
    averageTagLength = 0.0;
    minmapq = 10.0;
    peReadFlag = PE_READ_FILTER_KEEPBOTH;
    tbp = 0.0;
    maxtbp = 0.0;
    mintbp = 0.0;
    singleFile = 0;
    singleFilename = NULL;
    pairedEndFlag = 0;
    localInteractionFraction = -1.0;
    interChrInteractionFraction = -1.0;
    medianTagsPerPosition = -1;

    restrictionSite = NULL;
    tagFile = NULL;
    fragmentLengthSetting = FRAGMENT_LEN_AUTO;
    fragmentLengthEstimate = NULL_INT_VALUE;
    gsizeEstimate = 0;
    revStrand = 0;
    totalTagsMultiMappers = 0.0;
    totalPosMultiMappers = 0;
    manualTagTotal = TOTAL_READS_TAGDIR_DEFAULT;
}

TagLibrary::~TagLibrary() {
    if (chrs != NULL) {
        char **keys = chrs->keys();
        qsort(keys, chrs->total, sizeof(char *), &chrcmp);
        for (int i = 0; i < chrs->total; i++) {
            ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
            delete ct;
            delete[](keys[i]);
        }
        delete[]keys;
        delete chrs;
    }
    if (chrNames != NULL) delete chrNames;
    if (name != NULL) delete[]name;
}

void TagLibrary::setName(char *newname) {
    if (newname == NULL) return;
    if (name != NULL) delete[]name;
    name = new char[strlen(newname) + 1];
    strcpy(name, newname);
}

void TagLibrary::setSingleFile(int singleFileFlag) {
    singleFile = singleFileFlag;
    singleFilename = new char[10000];
    sprintf(singleFilename, "%s/genome.tags.tsv", directory);
}

void TagLibrary::setFragLength(int fraglength) {
    fragmentLengthSetting = fraglength;
    fragmentLengthEstimate = fraglength;
}

void TagLibrary::makeDirectory() {
    fprintf(stderr, "\tCreating directory: %s and removing existing *.tags.tsv\n", directory);
    char *filename = new char[10000];
    strcpy(filename, "mkdir -p \"");
    strcat(filename, directory);
    strcat(filename, "\"");
    (void) system(filename);
    strcpy(filename, "rm -f \"");
    strcat(filename, directory);
    strcat(filename, "\"/*.tags.tsv");
    (void) system(filename);
    delete[]filename;
}

void TagLibrary::setSingleRead(int flag) {
    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->forceSingleReadFlag = flag;
        delete[](keys[i]);
    }
    delete[]keys;
}

char *unzipFileIfNeeded(char *file, int &zipFlag, int &curFormat) {
    zipFlag = 0;
    if (file == NULL) {
        return NULL;
    }
    int len = strlen(file);
    char *newName = new char[len + 1];
    char *command = new char[100 + 2 * len + 1];
    strcpy(newName, file);
    if (len > 4) {
        if (file[len - 3] == '.' && file[len - 2] == 'g' && file[len - 1] == 'z') {
            //gzipped file
            fprintf(stderr, "\tTreating %s as a GNU zip file\n", file);
            newName[len - 3] = '\0';
            zipFlag = ZIPPED_FLAG_GZ;
            sprintf(command, "gunzip -c \"%s\" > \"%s\"", file, newName);
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'b' && file[len - 2] == 'z' && file[len - 1] == '2') {
            //bz2 zipped file
            fprintf(stderr, "\tTreating %s as a bz2 zip file\n", file);
            newName[len - 4] = '\0';
            zipFlag = ZIPPED_FLAG_BZ2;
            sprintf(command, "bunzip2 -c \"%s\" > \"%s\"", file, newName);
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'z' && file[len - 2] == 'i' && file[len - 1] == 'p') {
            //zip file
            fprintf(stderr, "\tTreating %s as a zip file\n", file);
            newName[len - 4] = '\0';
            zipFlag = ZIPPED_FLAG_ZIP;
            sprintf(command, "unzip -p \"%s\" > \"%s\"", file, newName);
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'b' && file[len - 2] == 'a' && file[len - 1] == 'm') {
            //zip file
            fprintf(stderr, "\tTreating %s as a bam file\n", file);
            newName[len - 4] = '\0';
            zipFlag = ZIPPED_FLAG_BAM;
            curFormat = FORMAT_SAM;
            sprintf(command, "samtools view -h \"%s\" > \"%s\"", file, newName);
            (void) system(command);
        } else {
            //not recognized, or not zipped...
        }

    }
    delete[]command;
    return newName;
}

void rezipFileIfNeeded(char *file, int zipFlag) {
    int len = strlen(file);
    char *command = new char[100 + 2 * len + 1];
    if (zipFlag) {
        sprintf(command, "rm \"%s\"", file);
        //fprintf(stderr, "Ending: %s\n", command);
        (void) system(command);
    }
    delete[]file;
    delete[]command;
}

void TagLibrary::parseAlignmentFiles(char **files, int numFiles, int format, int mode,
                                     char **tagDirs, int numDirs, char **tagFiles, int numTagFiles) {

    makeDirectory();


    for (int i = 0; i < numFiles; i++) {
        fprintf(stderr, "\n");
        if (pairedEndFlag) {
            fprintf(stderr, "\tReading paired end alignment files %s\n", files[i]);
            readPEAlignment(files[i], format, mode);
        } else {
            int zippedFlag = 0;
            int currentFormat = format;
            char *workingFilename = unzipFileIfNeeded(files[i], zippedFlag, currentFormat);
            fprintf(stderr, "\tReading alignment file %s\n", files[i]);
            readAlignment(workingFilename, currentFormat, mode, TAGLIBRARY_SINGLE_FLAG);
            rezipFileIfNeeded(workingFilename, zippedFlag);
        }
    }
    for (int i = 0; i < numDirs; i++) {
        fprintf(stderr, "\tAdding tag directory %s\n", tagDirs[i]);
        addTagDirectory(tagDirs[i]);
    }
    for (int i = 0; i < numTagFiles; i++) {
        fprintf(stderr, "\tAdding tag file %s\n", tagFiles[i]);
        readNewTagFile(tagFiles[i]);
    }
    if (totalTags < 1.0) {
        fprintf(stderr, "\t!!! Something is wrong - no reads were added to tag directory !!!\n");
        fprintf(stderr, "\t!!! Check your input files or the makeTagDirectory command options... !!!\n");
        exit(0);
    }
    //fprintf(stderr, "Troubleshooting...\n");
    //int ok = system("sleep 4");
    //ok = 1;
    optimizeTagFiles();
    if (manualTagTotal < TOTAL_READS_TAGDIR_ALL + 0.5) {
        totalTags = totalTagsMultiMappers;
        fprintf(stderr, "\t! Setting tag total to be equal to the total(including multimappers)\n");
    } else if (manualTagTotal > 0) {
        totalTags = manualTagTotal;
        fprintf(stderr, "\t! Setting tag total to be equal to %.2lf\n", totalTags);
    }
    fprintf(stderr, "\tTotal Tags = %.1f\n", totalTags);
    fprintf(stderr, "\tTotal Positions = %lld\n", totalPositions);
    printTagInfo();
}

void TagLibrary::addTagDirectory(char *tagDir) {

    DIR *dir = opendir(tagDir);
    char *file = NULL;
    struct dirent *entry = NULL;
    char *fullfile = new char[10000];
    while ((entry = readdir(dir))) {
        file = entry->d_name;
        int len = strlen(file);
        if (len < 9) continue;
        if (strcmp(".tags.tsv", &(file[len - 9])) == 0) {
            strcpy(fullfile, tagDir);
            strcat(fullfile, "/");
            strcat(fullfile, file);

            readNewTagFile(fullfile);

        }
    }
    closedir(dir);
    delete[]fullfile;


}

void TagLibrary::readPEAlignment(char *file, int format, int mode) {
    char **files = new char *[100];
    int numFiles = 0;
    split(file, files, numFiles, ',');
    if (numFiles < 2) {
        fprintf(stderr, "!!! Missing comma between file names for Paired End data !!!\n");
        fprintf(stderr, "!!! %s !!!\n", file);
        exit(0);
    }

    int zippedFlag = 0;
    int currentFormat = format;
    char *workingFilename = unzipFileIfNeeded(files[0], zippedFlag, currentFormat);
    //fprintf(stderr, "\tReading alignment file %s\n", files[i]);
    readAlignment(workingFilename, currentFormat, mode, TAGLIBRARY_PE_READ1_FLAG);
    rezipFileIfNeeded(workingFilename, zippedFlag);
    zippedFlag = 0;
    currentFormat = format;
    workingFilename = unzipFileIfNeeded(files[1], zippedFlag, currentFormat);
    //fprintf(stderr, "\tReading alignment file %s\n", files[i]);
    readAlignment(workingFilename, currentFormat, mode, TAGLIBRARY_PE_READ2_FLAG);
    rezipFileIfNeeded(workingFilename, zippedFlag);

    delete[]files;

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (chrNames == NULL) {
        chrNames = new Hashtable(1000);
    }
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        if (ct->tagfpR1 != NULL) fclose(ct->tagfpR1);
        if (ct->tagfpR2 != NULL) fclose(ct->tagfpR2);
        ct->tagfpR1 = NULL;
        ct->tagfpR2 = NULL;
        if (NULL == chrNames->search(keys[i])) {
            chrNames->insert(keys[i], keys[i]);
        }
    }

    //FILE* badfp = fopen("nomatch.tsv","w");

    //joinPEreads
    char *fname = new char[10000];
    char *buf = new char[BUFFER];
    char **line = new char *[BUFFER];
    int numCols = 0;
    fprintf(stderr, "\t\tMatching paired reads...\n");
    for (int i = 0; i < chrs->total; i++) {
        fprintf(stderr, "\t\t\t%s\n", keys[i]);
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        sprintf(fname, "%s.R1", ct->tagFile);
        FILE *fp = fopen(fname, "r");
        if (fp == NULL) {
            //fprintf(stderr, "Could not open %s\n", fname);
            continue;
        }


        Hashtable *read1tags = new Hashtable(10000000);
        //Inttable* found = new Inttable(1000000);
        while (fgets(buf, BUFFER, fp) != NULL) {
            split(buf, line, numCols, '\t');
            int p = 0;
            int len = 0;
            int d = 0;
            float v = 0.0;
            char *c = (char *) chrNames->search(line[1]);
            //fprintf(stderr, "|%s|\t|%s|\n", c,line[1]);
            sscanf(line[2], "%d", &p);
            sscanf(line[3], "%d", &d);
            sscanf(line[4], "%f", &v);
            sscanf(line[5], "%d", &len);
            PETag *petag = new PETag(line[0], c, p, (char) d, v, len);
            read1tags->insert(petag, line[0]);
            //found->insert(0,line[0]);
        }
        fclose(fp);

        for (int j = 0; j < chrs->total; j++) {
            ChrTags *ct2 = (ChrTags *) chrs->search(keys[j]);
            sprintf(fname, "%s.R2", ct2->tagFile);
            fp = fopen(fname, "r");
            if (fp == NULL) {
                //fprintf(stderr, "Could not open %s for matching to %s\n", fname,keys[i]);
                continue;
            }

            while (fgets(buf, BUFFER, fp) != NULL) {
                split(buf, line, numCols, '\t');
                PETag *petag = (PETag *) read1tags->search(line[0]);
                if (petag != NULL) {
                    int p = 0;
                    int d = 0;
                    float v = 0.0;
                    int len = 0;
                    char *c = (char *) chrNames->search(line[1]);
                    sscanf(line[2], "%d", &p);
                    sscanf(line[3], "%d", &d);
                    sscanf(line[4], "%f", &v);
                    sscanf(line[5], "%d", &len);
                    petag->p2 = p;
                    petag->d2 = (char) d;
                    petag->len2 = len;
                    petag->chr2 = c;
                    ct->printAlignedPETag(petag, 0);
                    ct2->printAlignedPETag(petag, 1);
                    //int fnum = found->search(line[0]);
                    //found->insert(fnum+1,line[0]);
                } else {
                }
            }
            fclose(fp);
        }
        char **hkeys = read1tags->keys();
        for (int j = 0; j < read1tags->total; j++) {
            PETag *petag = (PETag *) read1tags->search(hkeys[j]);
            //int fnum = found->search(hkeys[j]);
            //if (fnum == 0) petag->print(badfp);
            //if (fnum > 1) {
            //	fprintf(stderr, "Multiple matches: %d %s ",fnum , hkeys[j]);
            //	petag->print(stderr);
            //}
            delete petag;
            delete[](hkeys[j]);
        }
        delete[]hkeys;
        delete read1tags;
        //delete found;
    }

    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        sprintf(fname, "rm -f \"%s.R1\"", ct->tagFile);
        (void) system(fname);
        sprintf(fname, "rm -f \"%s.R2\"", ct->tagFile);
        (void) system(fname);
    }
    //joinPEreads
    delete[]fname;
    delete[]buf;
    delete[]line;
    delete[]keys;


}

void TagLibrary::readAlignment(char *file, int format, int mode, int PEflag) {

    FILE *fp = fopen(file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s\n", file);
        return;
    }

    char *buf = new char[BUFFER];
    char *name = new char[BUFFER];
    char *lastname = new char[BUFFER];
    char *chr = new char[BUFFER];
    char *lastchr = new char[BUFFER];
    char **line = new char *[BUFFER];
    char **line2 = new char *[BUFFER];
    int numCols = 0;
    char cigarCodes[1000];
    int cigarLens[1000];
    int numCodes;
    name[0] = '\0';
    lastname[0] = '\0';
    chr[0] = '\0';
    lastchr[0] = '\0';


    PETag *petag = new PETag();
    int pos = 0;
    int dir = 0;
    float numTags = 0;
    int len = -1;

    double avgReadsPerPos = 0;
    int totalReads = 0;

    while (fgets(buf, BUFFER, fp) != NULL) {
        //sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
        split(buf, line, numCols, '\t');

        if (numCols < 3) continue;


        if (format == FORMAT_UNKNOWN) {
            if (numCols > 2 && format == FORMAT_UNKNOWN) {
                if (line[0][0] == 'c' && line[0][1] == 'h' && line[0][2] == 'r') {
                    fprintf(stderr, "\tGuessing that your alignment file is BED format\n");
                    format = FORMAT_BED;
                }
            }
            if (numCols > 20 && format == FORMAT_UNKNOWN) {
                //fprintf(stderr, "\tGuessing that your alignment file is eland_export format\n");
                //format = FORMAT_ELANDEXPORT;
                fprintf(stderr, "\tGuessing that your alignment file is SAM format (lots of columns though...)\n");
                format = FORMAT_SAM;
            }
            if (line[0][0] == '@' && format == FORMAT_UNKNOWN) {
                if (strcmp(line[0], "@SQ") == 0 || strcmp(line[0], "@HD") == 0 || strcmp(line[0], "@RG") == 0) {
                    fprintf(stderr, "\tGuessing that your alignment file is SAM format\n");
                    format = FORMAT_SAM;
                }
            }
            if (format == FORMAT_UNKNOWN) {
                format = FORMAT_SAM;
                fprintf(stderr, "\tGuessing (by default) that your alignment file is sam format\n");
            }
        }

        if (format == FORMAT_SAM) {
            if (line[0][0] == '@') continue;
            if (line[2][0] == '*') continue;
            int samFlag = 0;
            sscanf(line[1], "%d", &samFlag);
            if (samFlag & 0x4) {
                //unmapped flag
                continue;
            }
            strcpy(name, line[0]);
            strcpy(chr, line[2]);
            int initLen = 0;
            len = strlen(line[9]);

            if (fragmentLengthSetting == FRAGMENT_LEN_PE && numCols > 8) {
                int pelen = 0;
                sscanf(line[8], "%d", &pelen);
                if (pelen < 0) pelen *= -1;
                if (pelen > 0) {
                    len = pelen;
                }
            }

            if (len < minReadLength || len > maxReadLength) continue;
            sscanf(line[3], "%d", &pos);
            if (pos == 0) continue;
            if (line[5][0] == '*') continue;

            dir = STRAND_POSITIVE;
            int peFlag = 0;
            if (samFlag & 0x10) {
                dir = STRAND_NEGATIVE;
            }

            if (sspeFlag && (samFlag & 0x80)) {
                if (dir == STRAND_POSITIVE) {
                    dir = STRAND_NEGATIVE;
                } else {
                    dir = STRAND_POSITIVE;
                }
            }

            if (samFlag & 0x1) {
                peFlag = 1;
            }

            if (peReadFlag == PE_READ_FILTER_KEEPFIRST && (samFlag & 0x80)) {
                continue;
            } else if (peReadFlag == PE_READ_FILTER_KEEPSECOND && !(samFlag & 0x80)) {
                continue;
            }

            if (samFlag & 0x100) {
            } else {
                if (peFlag) totalTagsMultiMappers += 0.5;
                else totalTagsMultiMappers += 1.0;
            }

            if (mode == MODE_UNIQUE) {
                if (samFlag & 0x100) continue;
                if (minmapq < 0) {
                    int good = 1;
                    double alnScore = -1e8;
                    double secScore = -1e9;
                    for (int i = 11; i < numCols; i++) {
                        if (strlen(line[i]) < 6) continue;
                        if (strncmp(line[i], "X0:i:", 5) == 0) {
                            if (strcmp(&(line[i][5]), "1") != 0) {
                                good = 0;
                            }
                        } else if (strncmp(line[i], "AS:i:", 5) == 0) {
                            sscanf(&(line[i][5]), "%lf", &alnScore);
                        } else if (strncmp(line[i], "XS:i:", 5) == 0) {
                            sscanf(&(line[i][5]), "%lf", &secScore);
                        }
                    }
                    if (alnScore > secScore) {
                    } else {
                        good = 0;
                    }
                    if (good == 0) continue;
                } else {
                    double mapqScore = 0.0;
                    sscanf(line[4], "%lf", &mapqScore);
                    if (mapqScore < minmapq) continue;
                }
            } else if (mode == MODE_KEEPONE) {
                if (samFlag & 0x100) continue;
            }
            int startPos = getRightCoordFromCIGAR(line[5], dir, cigarCodes, cigarLens, numCodes, initLen);
            if (startPos == CIGAR_ERROR) {
                fprintf(stderr, "\tError in read: %s\n", name);
                startPos = 0;
            }
            pos += startPos;
            if (fragmentLengthSetting != FRAGMENT_LEN_PE) {
                len = initLen;
            }
            float v = 1.0;
            if (peFlag) v = 0.5;
            addAlignedTag(name, chr, pos, (char) dir, len, v, PEflag);
        }
        else if (format == FORMAT_BED) {
            strcpy(chr, line[0]);
            name[0] = '\0';

            dir = STRAND_POSITIVE;
            numTags = 1.0;
            if (numCols > 3) {
                if (line[3][0] == '+' || line[3][0] == '-') {
                    if (line[3][0] == '-') dir = STRAND_NEGATIVE;
                } else if (numCols > 5) {
                    strcpy(name, line[3]);
                    if (line[5][0] == '-') {
                        dir = STRAND_NEGATIVE;
                    }
                }
            }
            if (numCols > 4 && mode == MODE_BED_FORCE5TH) {
                sscanf(line[4], "%f", &numTags);
            }
            int start = 0;
            int end = 0;
            sscanf(line[1], "%d", &start);
            sscanf(line[2], "%d", &end);
            if (dir == 0) {
                pos = start + 1; //0 base index
            } else {
                pos = end;
            }
            len = end - start;
            if (len < minReadLength || len > maxReadLength) continue;
            addAlignedTag(name, chr, pos, (char) dir, len, numTags, PEflag);
            avgReadsPerPos += numTags;
            totalReads++;
        }
    }
    fclose(fp);

    if (format == FORMAT_BED) {
        if (totalReads > 1) {
            avgReadsPerPos /= (double) totalReads;
            if (avgReadsPerPos > 10.0) {
                fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                fprintf(stderr, "Average reads per BED file: %.1lf\n", avgReadsPerPos);
                fprintf(stderr, "Good chance that the 5th column of your BED file has a weird value in it!\n");
                fprintf(stderr, "By default, this is read as the number of reads for that position\n");
                fprintf(stderr, "To count each entry as only one read (ignore the 5th column) use -forceBED\n");
                fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            }
        }
    }


    delete petag;
    delete[]buf;
    delete[]chr;
    delete[]lastchr;
    delete[]name;
    delete[]lastname;
    delete[]line;
    delete[]line2;
}

int
TagLibrary::getRightCoordFromCIGAR(char *str, int dir, char *cigarCodes, int *cigarLens, int &numCodes, int &initLen) {
    numCodes = 0;
    int i = 0;
    char *start = str;
    int lastStartI = 0;
    int totalLen = 0;

    int firstStart = -1;
    int firstEnd = -1;
    int lastStart = -1;
    int lastEnd = -1;
    int firstActive = 0;
    int lastActive = 0;
    //fprintf(stderr, "%s\n",str);
    while (str[i] != '\0') {
        if (str[i] < 48 || str[i] > 57) {
            if (i - lastStartI == 0) {
                fprintf(stderr, "Parsing problem with CIGAR string...\n");
                return CIGAR_ERROR;
            }
            cigarCodes[numCodes] = str[i];
            str[i] = '\0';
            sscanf(start, "%d", &(cigarLens[numCodes]));
            //fprintf(stderr, " %c|%s|%d ", cigarCodes[numCodes],start,cigarLens[numCodes]);
            if (cigarCodes[numCodes] == 'M' || cigarCodes[numCodes] == '=') {
                if (firstStart < 0) {
                    firstStart = totalLen;
                    firstActive = 1;
                }
                if (firstActive) {
                    firstEnd = totalLen + cigarLens[numCodes];
                }
                if (lastActive == 0) {
                    lastStart = totalLen;
                    lastActive = 1;
                }
                if (lastActive) {
                    lastEnd = totalLen + cigarLens[numCodes];
                }
                totalLen += cigarLens[numCodes];
            } else if (cigarCodes[numCodes] == 'D' || cigarCodes[numCodes] == 'N' || cigarCodes[numCodes] == 'X'
                       || cigarCodes[numCodes] == 'H' || cigarCodes[numCodes] == 'P') {
                totalLen += cigarLens[numCodes];
                firstActive = 0;
                lastActive = 0;
            } else if (cigarCodes[numCodes] == 'I' || cigarCodes[numCodes] == 'S') {
                totalLen -= cigarLens[numCodes];
            }
            start = &(str[i + 1]);
            lastStartI = i + 1;
            numCodes++;
        }
        i++;
    }
    int rv = totalLen;
    if (dir == STRAND_POSITIVE) {
        rv = firstStart;
        initLen = firstEnd - firstStart;
    } else {
        rv = lastEnd - 1;
        initLen = lastEnd - lastStart;
    }
    //fprintf(stderr, "\t%d\t%d\t%d\n", rv, initLen,dir);
    return rv;
}


void TagLibrary::readSingleTagFile() {

    if (chrs != NULL && chrs->total > 0) {
        char **keys = chrs->keys();
        qsort(keys, chrs->total, sizeof(char *), &chrcmp);
        for (int i = 0; i < chrs->total; i++) {
            ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
            ct->freeTags();
            delete[](keys[i]);
        }
        delete[]keys;
    }

    if (tagFile != NULL) {
        fclose(tagFile);
    }
    if (singleFilename == NULL) {
        fprintf(stderr, "Single tag file is not active\n");
        exit(0);
    }
    tagFile = fopen(singleFilename, "r");
    if (tagFile == NULL) {
        fprintf(stderr, "!!! Could not open %s for reading !!!\n", singleFilename);
        exit(0);
    }

    char *buf = new char[BUFFER];
    char *chr = NULL;
    char *name = new char[BUFFER];
    char **line = new char *[BUFFER];
    int numCols = 0;
    int pos = 0;
    float tagCount = 0;
    int dir = 0;
    int len = -1;

    int optimizedFlag = 1;
    int numChangedChr = 0;
    char *lastChr = new char[10000];
    int lastPos = 0;


    while (fgets(buf, BUFFER, tagFile) != NULL) {
        //sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
        split(buf, line, numCols, '\t');
        if (numCols < 5) continue;

        chr = line[1];
        sscanf(line[2], "%d", &pos);
        sscanf(line[3], "%d", &dir);
        sscanf(line[4], "%f", &tagCount);
        len = -1;
        if (numCols > 5) {
            sscanf(line[5], "%d", &len);
        }


        ChrTags *ct = (ChrTags *) chrs->search(chr);
        if (ct == NULL) {
            ct = new ChrTags(chr);
            ct->singleFile = 1;
            ct->mCflag = mCflag;
            ct->pairedEndFlag = pairedEndFlag;
            ct->tagfp = tagFile;
            chrs->insert(ct, chr);
        }
        if (pos > ct->appearentSize) {
            ct->appearentSize = pos;
        }
        ct->addTag(pos, (char) dir, len, tagCount);
        if (strcmp(chr, lastChr) != 0) {
            numChangedChr++;
            strcpy(lastChr, chr);
        } else {
            if (lastPos > pos) optimizedFlag = 0;
        }
        lastPos = pos;

    }
    fclose(tagFile);
    tagFile = NULL;

    if (numChangedChr > chrs->total) {
        optimizedFlag = 0;
    }
    if (optimizedFlag == 0) {
        fprintf(stderr, "\tOptimizing single genome.tags.tsv file...\n");
        tagFile = fopen(singleFilename, "w");
    }


    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);

    gsizeEstimate = 0;
    totalTags = 0.0;
    totalPositions = 0;
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->optimizedFlag = optimizedFlag;
        ct->optimizeTags();
        if (optimizedFlag == 0) {
            ct->tagfp = tagFile;
            ct->print();
        }
        ct->loaded = 1;
        totalTags += ct->totalTags;
        totalPositions += ct->totalPositions;
        gsizeEstimate += ct->appearentSize;
        //fprintf(stderr, "\t%s\t%lld\n", keys[i],ct->appearentSize);
    }
    if (optimizedFlag == 0) {
        fprintf(stderr, "\tEstimated genome size = %lld\n", gsizeEstimate);
        fclose(tagFile);
        tagFile = NULL;
    }

    delete[]lastChr;
    delete[]buf;
    delete[]name;
    delete[]line;
}


void TagLibrary::readNewTagFile(char *filename) {


    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s\n", filename);
        return;
    }

    char *buf = new char[BUFFER];
    char *chr = NULL;
    char *name = NULL;
    char **line = new char *[BUFFER];
    int numCols = 0;
    int pos = 0;
    int pos2 = 0;
    float tagCount = 0;
    int dir = 0;
    int dir2 = 0;
    int len = -1;
    int len2 = -1;

    PETag *petag = new PETag();

    //pairedEndFlag = -1;

    while (fgets(buf, BUFFER, fp) != NULL) {
        //sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
        split(buf, line, numCols, '\t');
        if (numCols < 5) continue;

        name = line[0];
        chr = line[1];
        sscanf(line[2], "%d", &pos);
        sscanf(line[3], "%d", &dir);
        sscanf(line[4], "%f", &tagCount);
        len = -1;
        if (numCols > 5) {
            sscanf(line[5], "%d", &len);
        }
        if (numCols < 9) {
            if (len < minReadLength || len > maxReadLength) continue;
            addAlignedTag(line[0], line[1], pos, (char) dir, len, tagCount, TAGLIBRARY_SINGLE_FLAG);
        } else {
            petag->name = name;
            petag->chr1 = chr;
            petag->p1 = pos;
            petag->d1 = dir;
            petag->v = tagCount;
            petag->len1 = len;
            sscanf(line[7], "%d", &pos2);
            sscanf(line[8], "%d", &dir2);
            if (numCols > 9) {
                sscanf(line[9], "%d", &len2);
            }
            petag->chr2 = line[6];
            petag->p2 = pos2;
            petag->d2 = dir2;
            petag->len2 = len2;
            if (len2 < minReadLength || len2 > maxReadLength) continue;
            addAlignedPETag(petag);
            pairedEndFlag = 1;
        }
    }
    fclose(fp);

    delete petag;
    delete[]buf;
    delete[]line;
}


void TagLibrary::optimizeTagFiles() {

    fprintf(stderr, "\n\tOptimizing tag files...\n");
    totalTags = 0;
    totalPositions = 0;

    setMaxTBP(maxtbp);
    setMinTBP(mintbp);

    if (singleFile) {
        readSingleTagFile();
        return;
    }

    char **chrKeys = chrs->keys();
    qsort(chrKeys, chrs->total, sizeof(char *), &chrcmp);
    gsizeEstimate = 0;
    for (int i = 0; i < chrs->total; i++) {
        //fprintf(stderr, "Optimizing %s\n", chrKeys[i]);
        ChrTags *ct = (ChrTags *) chrs->search(chrKeys[i]);
        ct->optimizeTagFile();
        totalTags += ct->totalTags;
        totalPositions += ct->totalPositions;
        gsizeEstimate += ct->appearentSize;
        //fprintf(stderr, "\t\t%s: %lld\n", chrKeys[i], ct->appearentSize);
        delete[](chrKeys[i]);
    }
    fprintf(stderr, "\tEstimated genome size = %lld\n", gsizeEstimate);

    if (gsizeEstimate > 1) {
        tbp = totalTags / gsizeEstimate;
        fprintf(stderr, "\tEstimated average read density = %.6lf per bp\n", tbp);
    }

    delete[]chrKeys;

}

void TagLibrary::addAlignedPETag(PETag *petag) {
    char *chr = petag->chr1;
    ChrTags *chrtags = (ChrTags *) chrs->search(chr);
    if (chrtags == NULL) {
        chrtags = new ChrTags(chr);
        chrtags->pairedEndFlag = 1;
        chrs->insert(chrtags, chr);

        if (singleFile == 0) {
            char *file = getTagFileName(chr);
            chrtags->setTagFile(file);
            chrtags->openTagFile((char *) "w");
            delete[]file;
        } else {
            if (tagFile == NULL) {
                singleFilename = new char[10000];
                sprintf(singleFilename, "%s/genome.tags.tsv", directory);
                tagFile = fopen(singleFilename, "w");
            }
            chrtags->tagfp = tagFile;
            chrtags->singleFile = 1;
        }
    }
    chrtags->printAlignedPETag(petag, 0);
    totalTags += petag->v;
}

void TagLibrary::addAlignedTag(char *name, char *chr, int pos, char dir, int length,
                               float value, int PEflag) {

    static int warning = 0;

    if (pos < 1) pos = 1;
    ChrTags *chrtags = (ChrTags *) chrs->search(chr);
    if (chrtags == NULL) {
        chrtags = new ChrTags(chr);
        chrtags->pairedEndFlag = pairedEndFlag;
        chrtags->mCflag = mCflag;
        chrs->insert(chrtags, chr);

        if (singleFile == 0) {
            if (chrs->total > 1000 && warning == 0) {
                warning = 1;
                fprintf(stderr,
                        "!!! Warning: more than 1000 chromosomes detected.  This can casue file I/O problems\n");
                fprintf(stderr, "!!! If the command fails, consider reruning makeTagDirectories with \"-single\"\n");
                fprintf(stderr, "!!! so that it uses a single tag file instead of one per chromosome\n");
            }

            char *file = getTagFileName(chr);
//fprintf(stderr, "||%s||\n", file);
            chrtags->setTagFile(file);
            chrtags->openTagFile((char *) "w");
            delete[]file;
        } else {

            if (tagFile == NULL) {
                char *file = new char[10000];
                sprintf(file, "%s/genome.tags.tsv", directory);
                tagFile = fopen(file, "w");
                delete[]file;
            }
            chrtags->tagfp = tagFile;
            chrtags->singleFile = 1;
        }
    }
    if (pairedEndFlag && chrtags->tagfpR1 == NULL) {
        char *readFile = new char[10000];
        sprintf(readFile, "%s.R1", chrtags->tagFile);
        chrtags->tagfpR1 = fopen(readFile, "w");
        sprintf(readFile, "%s.R2", chrtags->tagFile);
        chrtags->tagfpR2 = fopen(readFile, "w");
        delete[]readFile;
    }
    chrtags->printAlignedTag(name, pos, dir, length, value, PEflag);
    totalTags += value;
}

char *TagLibrary::getTagFileName(char *chrname) {
    char *file = new char[10000];
    sprintf(file, "%s/%s.tags.tsv", directory, chrname);
    return file;
}

char *TagLibrary::getDirFileName(char *filename) {
    char *file = new char[10000];
    strcpy(file, directory);
    strcat(file, "/");
    strcat(file, filename);
    return file;
}


void TagLibrary::printTagInfo() {
    printTagInfo(NULL);
}

void TagLibrary::printTagInfo(FILE *fp) {

    char *filename = new char[10000];
    strcpy(filename, directory);
    strcat(filename, "/tagInfo.txt");

    FILE *info = NULL;
    if (fp == NULL) {
        info = fopen(filename, "w");
    } else {
        info = fp;
    }
    if (info == NULL) {
        fprintf(stderr, "Could not open file for output in directory %s\n", directory);
        return;
    }

    if (name != NULL) {
        fprintf(info, "name=%s\tUnique Positions\tTotal Tags\n", name);
    } else {
        fprintf(info, "name\tUnique Positions\tTotal Tags\n");
    }

    fprintf(info, "genome\t%lld\t%.1lf\n", totalPositions, totalTags);

    if (pairedEndFlag) {
        fprintf(info, "pairedEnd=true\n");
    }
    if (localInteractionFraction > -0.5) {
        fprintf(info, "localInteractionFraction=%lf\n", localInteractionFraction);
    }
    if (interChrInteractionFraction > -0.5) {
        fprintf(info, "interChrInteractionFraction=%lf\n", interChrInteractionFraction);
    }
    if (mCflag) {
        fprintf(info, "mC=true\n");
    }
    if (restrictionSite != NULL) {
        fprintf(info, "restrictionSite=%s\n", restrictionSite);
    }
    fprintf(info, "fragmentLengthEstimate=%d\t\t\n", fragmentLengthEstimate);
    fprintf(info, "tagsPerBP=%lf\t\t\n", tbp);
    fprintf(info, "averageTagsPerPosition=%.3lf\t\t\n", averageTagsPerPosition);
    fprintf(info, "medianTagsPerPosition=%d\t\t\n", medianTagsPerPosition);
    fprintf(info, "averageTagLength=%.3lf\t\t\n", averageTagLength);
    fprintf(info, "gsizeEstimate=%lld\t\t\n", gsizeEstimate);
    if (singleFile) {
        fprintf(info, "singleTagFile=true\t\t\n");
    }

    char **chr = chrs->keys();
    qsort(chr, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(chr[i]);
        fprintf(info, "%s\t%d\t%.1f\n", chr[i], ct->totalPositions, ct->totalTags);
        delete[](chr[i]);
    }

    if (fp == NULL) fclose(info);
    delete[]chr;
    delete[]filename;

}


double *TagLibrary::getPETagDistribution(int windowSize, int largeWindowSize, int largeResolution,
                                         char *outputPrefix, int &distLength) {
    FILE *fp = NULL;
    FILE *localfp = NULL;

    if (outputPrefix != NULL) {
        char *fname = new char[100000];
        sprintf(fname, "%s.FreqDistribution_%d.txt", outputPrefix, largeResolution);
        char *file = getDirFileName(fname);
        fp = fopen(file, "w");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open %s for writing\n", file);
            return NULL;
        }
        delete[]file;

        sprintf(fname, "%s.LocalDistribution.txt", outputPrefix);
        file = getDirFileName(fname);
        localfp = fopen(file, "w");
        if (localfp == NULL) {
            fprintf(stderr, "Cannot open %s for writing\n", file);
            return NULL;
        }
        delete[]fname;
    }

    double *sameStrand = new double[windowSize];
    double *diffStrand = new double[windowSize];
    double *smoothed = new double[windowSize];
    for (int i = 0; i < windowSize; i++) {
        sameStrand[i] = 0.0;
        diffStrand[i] = 0.0;
        smoothed[i] = 0.0;
    }


    int largeLength = (largeWindowSize / largeResolution) + 2;
    double *largeWindow = new double[largeLength];
    double *largeWindowN = new double[largeLength];
    for (int i = 0; i < largeLength; i++) {
        largeWindow[i] = 0.0;
        largeWindowN[i] = 0.0;
    }
    distLength = largeLength;


    //int genomeIndexSize = gsizeEstimate/largeResolution+1;

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (singleFile) readSingleTagFile();
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        int chrLength = ct->getPETagDistribution(sameStrand, diffStrand, windowSize,
                                                 largeWindow, largeResolution, largeLength);
        delete[](keys[i]);
        int indexLen = chrLength / largeResolution + 1;
        for (int j = 0; j < indexLen; j++) {
            int z = j;
            if (j > largeLength - 2) z = largeLength - 2;
            largeWindowN[z] += (double) (indexLen - j);
        }
        largeWindowN[largeLength - 1] += (double) indexLen;
    }
    delete[]keys;

    double sum = 0.0;
    for (int i = 0; i < largeLength; i++) {
        if (largeWindowN[i] > 0.0) {
            largeWindow[i] /= largeWindowN[i];
        }
        sum += largeWindow[i];
    }
    for (int i = 0; i < largeLength; i++) {
        largeWindow[i] /= sum;
    }

    interChrInteractionFraction = largeWindow[largeLength - 1];
    localInteractionFraction = largeWindow[0];
    fprintf(stderr, "\tLocal interaction fraction (< 1kb): %.2lf%%\n", localInteractionFraction * 100.0);
    fprintf(stderr, "\tInterchromosomal interaction fraction: %.2lf%%\n", interChrInteractionFraction * 100.0);
    //int minFragLen = AUTOCORRELATION_OFFSETMIN;
    //if (averageTagLength > 0.0) minFragLen = (int) averageTagLength+3;

    fragmentLengthEstimate = 0;

    if (localfp != NULL) {
        int halfWindow = windowSize / 2;
        fprintf(localfp, "Local Distance in bp between PE tags\tSame Strand\tOpposite Strands\n");
        for (int i = 0; i < windowSize; i++) {
            int offset = i - halfWindow;
            fprintf(localfp, "%d\t%lf\t%lf\n", offset, sameStrand[i], diffStrand[i]);
        }
        fclose(localfp);
    }
    if (fp != NULL) {
        fprintf(fp, "Distance between PE tags\tFraction of total PE tags");
        //fprintf(fp, "(Interchromosomal:%le)\n", largeWindow[largeLength-1]/totalTags);
        fprintf(fp, "(Interchromosomal:%le)\n", largeWindow[largeLength - 1]);
        for (int i = 0; i < largeLength - 1; i++) {
            if (i == largeLength - 2) fprintf(fp, "More than ");
            //fprintf(fp, "%d\t%le\n", i*largeResolution,largeWindow[i]/totalTags);
            fprintf(fp, "%d\t%le\n", i * largeResolution, largeWindow[i]);
        }
        fclose(fp);
    }

    delete[]sameStrand;
    delete[]diffStrand;
    delete[]smoothed;

    return largeWindow;
}

void TagLibrary::setMaxTBP(float max) {

    maxtbp = max;

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->setMaxTBP(max);
        delete[](keys[i]);
    }
    if (keys != NULL) delete[]keys;
}

void TagLibrary::setMinTBP(float min) {

    mintbp = min;

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->setMinTBP(min);
        delete[](keys[i]);
    }
    if (keys != NULL) delete[]keys;
}

void TagLibrary::setTagAdjust(int dist) {
    if (dist == TAGADJUST_AUTO) {
        if (0) { //fragmentLengthEstimate <= AUTOCORRELATION_OFFSETMIN) {
            dist = TAGADJUST_DEFAULT;
            fragmentLengthEstimate = dist * 2;
        } else {
            dist = (fragmentLengthEstimate) / 2;
        }
    }
    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->setTagAdjust(dist);
        ct->revStrand = revStrand;
        delete[](keys[i]);
    }
    delete[]keys;
}

PeakLibrary *TagLibrary::findPutativePeaks(int peakSize, int minDist, char strand, float minCount) {

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    PeakLibrary *putativePeaks = new PeakLibrary(10000000);
    if (singleFile) readSingleTagFile();
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->loadTags();
        ct->findPutativePeaks(putativePeaks, peakSize, minDist, strand, minCount);
        ct->freeTags();
        delete[](keys[i]);
    }
    delete[]keys;

    putativePeaks->sortChr();
    return putativePeaks;
}

//GRO-Seq transcript identification


double *TagLibrary::getTagLengthDistribution(FILE *nfp, int &max) {
    FILE *fp = nfp;
    if (fp == NULL) {
        char *file = NULL;
        if (mCflag) {
            file = getDirFileName((char *) "mCreadCoverageDistribution.txt");
        } else {
            file = getDirFileName((char *) "tagLengthDistribution.txt");
        }

        fp = fopen(file, "w");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open %s for writing\n", file);
            return NULL;
        }
    }

    if (max == 0) {
        max = MAX_READ_LENGTH;
    }
    double *dist = new double[max];
    for (int i = 0; i < max; i++) {
        dist[i] = 0;
    }


    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (singleFile) readSingleTagFile();
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->getTagLengthDistribution(dist, max);
        delete[](keys[i]);
    }
    delete[]keys;


    float totalDist = 0;
    averageTagLength = 0.0;
    int maxValue = 0;
    for (int i = 0; i < max; i++) {
        dist[i] /= (double) totalTags;
        if (dist[i] > 0) maxValue = i;
        averageTagLength += ((double) i) * dist[i];
        totalDist += dist[i];
    }

    if (mCflag) {
        fprintf(stderr, "\tAverage read depth (from read depth of data exceeding %.1lf threshold) = %.1lf\n",
                parseAlignmentCpGMinValue, averageTagLength);
        fprintf(fp, "Read depth (Average read depth = %lf)", averageTagLength);
        fprintf(fp, "\tFraction of Positions\n");
    } else {
        fprintf(stderr, "\tAverage tag length = %.1lf\n", averageTagLength);
        fprintf(fp, "Tag Length (Average tag length = %lf)", averageTagLength);
        fprintf(fp, "\tFraction of Tags\n");
    }
    for (int i = 0; i < max; i++) {
        fprintf(fp, "%d\t%lf\n", i, dist[i]);
        if (i >= maxValue) break;
    }
    if (nfp == NULL) {
        fclose(fp);
    }
    return dist;

}


double *TagLibrary::getTagCountDistribution(FILE *nfp, int &max) {
    FILE *fp = nfp;
    if (fp == NULL) {
        char *file = NULL;
        if (mCflag) {
            file = getDirFileName((char *) "mCratioDistribution.txt");
        } else {
            file = getDirFileName((char *) "tagCountDistribution.txt");
        }
        fp = fopen(file, "w");
        if (fp == NULL) {
            fprintf(stderr, "Cannot open %s for writing\n", file);
            return NULL;
        }
    }

    if (max == 0) {
        max = MAX_TAGS_PER_BP;
    }
    double *dist = new double[max];
    for (int i = 0; i < max; i++) {
        dist[i] = 0;
    }

    int scaleFactor = 1;
    if (mCflag) {
        scaleFactor = 20;
    }

    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (singleFile) readSingleTagFile();
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        ct->getTagCountDistribution(dist, max, scaleFactor);
        delete[](keys[i]);
    }
    delete[]keys;


    float totalDist = 0;
    averageTagsPerPosition = 0;
    int maxValue = 0;
    for (int i = 0; i < max; i++) {
        dist[i] /= (double) totalPositions;
        if (dist[i] > 0) maxValue = i;
        averageTagsPerPosition += ((double) i) / ((double) scaleFactor) * dist[i];
        if (totalDist < 0.5 && totalDist + dist[i] > 0.5) {
            medianTagsPerPosition = i;
        }
        totalDist += dist[i];
    }

    if (mCflag) {
        double median = ((double) medianTagsPerPosition) / ((double) scaleFactor);
        fprintf(stderr, "\tMedian mC/C = %.2lf\n", median);
        fprintf(stderr, "\tAverage mC/C = %.3lf\n", averageTagsPerPosition);
        fprintf(fp, "Tags per tag position (Median = %.2lf, tags per genomic bp = %.3lf)",
                median, tbp);
        fprintf(fp, "\tFraction of Positions\n");
        for (int i = 0; i < max; i++) {
            fprintf(fp, "%.2lf\t%lf\n", ((double) i) / ((double) scaleFactor), dist[i]);
            if (i >= maxValue) break;
        }
    } else {

        int idealTbp = (int) ceil(tbp * 2 + 0.001);

        fprintf(stderr, "\tMedian tags per position = %d (ideal: %d)\n", medianTagsPerPosition, idealTbp);
        fprintf(stderr, "\tAverage tags per position = %.3lf\n", averageTagsPerPosition);
        if (medianTagsPerPosition > idealTbp * 3) {
            fprintf(stderr, "\t\t!! Might have some clonal amplification in this sample if sonication was used\n");
            fprintf(stderr,
                    "\t\tIf this is ChIP-Seq using sonicated fragments, consider adding the option \"-tbp %d\"\n",
                    idealTbp);
            fprintf(stderr, "\t\tIgnore if analyzing RNA, MNase, etc. data\n");
        } else if (medianTagsPerPosition > idealTbp) {
            fprintf(stderr, "\t\t!! Might have some clonal amplification in this sample if sonication was used\n");
            fprintf(stderr, "\t\tIf using BED files, try using -forceBED to ignore the 5th column\n");
            fprintf(stderr, "\t\tIgnore if analyzing RNA, MNase, etc. data\n");
        }
        fprintf(fp, "Tags per tag position (Median = %d, tags per genomic bp = %.3lf)",
                medianTagsPerPosition, tbp);
        fprintf(fp, "\tFraction of Positions\n");
        for (int i = 0; i < max; i++) {
            fprintf(fp, "%d\t%lf\n", i, dist[i]);
            if (i >= maxValue) break;
        }
    }

    if (nfp == NULL) {
        fclose(fp);
    }

    return dist;
}

double TagLibrary::getAdjustedTagTotal() {
    double total = 0;
    char **keys = chrs->keys();
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        total += ct->totalTags;
        delete[](keys[i]);
    }
    delete[]keys;
    return total;
}


void TagLibrary::readAndSave() {
    char **keys = chrs->keys();
    totalTags = 0.0;
    totalPositions = 0;
    qsort(keys, chrs->total, sizeof(char *), &chrcmp);
    if (singleFile) {
        readSingleTagFile();
        tagFile = fopen(singleFilename, "w");
        //fprintf(stderr, "\t\t%s\t%d\n",singleFilename,tagFile);
    }
    for (int i = 0; i < chrs->total; i++) {
        ChrTags *ct = (ChrTags *) chrs->search(keys[i]);
        if (singleFile) ct->tagfp = tagFile;
        ct->readAndSave();
        totalTags += ct->totalTags;
        totalPositions += ct->totalPositions;
        delete[](keys[i]);
    }
    delete[]keys;
    if (singleFile) {
        fclose(tagFile);
        tagFile = NULL;
    }
    printTagInfo();
}



// class ChrTags ------------------------------------------------------------


ChrTags::ChrTags(char *newchr) {
    chr = NULL;
    if (newchr != NULL) {
        chr = new char[strlen(newchr) + 1];
        strcpy(chr, newchr);
    }
    linkTag = NULL;
    firstTag = NULL;
    tags = NULL;
    petags = NULL;
    mCflag = 0;
    totalTags = 0.0;
    totalPositions = 0;
    adjusted = 0;
    optimizeOverride = 0;
    maxtbp = 0.0;
    mintbp = 0.0;
    tagAdjust = 0;
    numLinks = 0;
    loaded = 0;
    tagFile = NULL;
    tagfp = NULL;
    tagfpR1 = NULL;
    tagfpR2 = NULL;
    dontSAVE = 0;
    singleFile = 0;
    gcFreq = NULL;
    seq = NULL;
    chrNames = new Hashtable(1000);
    firstPETag = NULL;
    linkPETag = NULL;
    pairedEndFlag = 0;
    forceSingleReadFlag = 0;
    appearentSize = 0;
    revStrand = 0;
}

ChrTags::~ChrTags() {
    if (tags != NULL) delete[]tags;
    if (chr != NULL) delete[]chr;
    if (tagFile != NULL) delete[]tagFile;
    if (tagfp != NULL) fclose(tagfp);
    if (gcFreq != NULL) delete[]gcFreq;
    if (seq != NULL) delete[]seq;
    linkTag = NULL;
    firstTag = NULL;
    firstPETag = NULL;
    linkPETag = NULL;
    tags = NULL;
    if (chrNames != NULL) {
        char **keys = chrNames->keys();
        for (int i = 0; i < chrNames->total; i++) {
            char *cname = (char *) chrNames->search(keys[i]);
            if (cname != NULL) delete[]cname;
            delete[](keys[i]);
        }
        delete[]keys;
        delete chrNames;
    }
}

void ChrTags::freeTags() {
    if (tags != NULL) delete[]tags;
    if (petags != NULL) delete[]petags;
    tags = NULL;
    petags = NULL;
    loaded = 0;
    adjusted = 0;
}

void ChrTags::setTagFile(char *file) {
    if (file == NULL) return;
    tagFile = new char[strlen(file) + 1];
    strcpy(tagFile, file);
}

void ChrTags::openTagFile(char *mode) {
    if (tagFile == NULL) return;
    closeTagFile();
    tagfp = fopen(tagFile, mode);
    if (tagfp == NULL) {
        fprintf(stderr, "!!!!! Could not open file %s for printing tags!!!!!\n", tagFile);
        fprintf(stderr, "!!!!! Is this a valid file name?  May need to:\n");
        fprintf(stderr, "\tTry a different name for the tag directory (is the same name as an existing file?)\n");
        fprintf(stderr, "\t-or- you may need to rename your chromosomes if they have weird characters!\n");
        fprintf(stderr, "\t-or- try using the \"-single\" flag.\n");
        exit(0);
    }
}

void ChrTags::closeTagFile() {
    if (tagfp != NULL) fclose(tagfp);
    if (tagfpR1 != NULL) fclose(tagfpR1);
    if (tagfpR2 != NULL) fclose(tagfpR2);
    tagfp = NULL;
    tagfpR1 = NULL;
    tagfpR2 = NULL;
}

void ChrTags::printAlignedTag(char *rname, int pos, char dir, int length, float value, int PEflag) {
    FILE *fp = tagfp;
    if (PEflag == TAGLIBRARY_PE_READ1_FLAG) {
        fp = tagfpR1;
    } else if (PEflag == TAGLIBRARY_PE_READ2_FLAG) {
        fp = tagfpR2;
    }
    if (Tag::precision == 1) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%.1f\t%d\n", rname, chr, pos, dir, value, length);
    } else if (Tag::precision == 2) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%.2f\t%d\n", rname, chr, pos, dir, value, length);
    } else {
        fprintf(fp, "%s\t%s\t%d\t%d\t%.3f\t%d\n", rname, chr, pos, dir, value, length);
    }
}

void ChrTags::printAlignedPETag(PETag *petag, int revFlag) {
    petag->print(tagfp, revFlag);
}


void ChrTags::loadTags() {
    if (loaded && adjusted) return;
    if (pairedEndFlag && forceSingleReadFlag == 1) {
        optimizeOverride = 1;
        //fprintf(stderr, "Reading single...\n");
    }
    if (singleFile == 0 && loaded == 0) {
        readTagFile();
    }
    adjustTags();
    if (pairedEndFlag && forceSingleReadFlag == 1) optimizeOverride = 0;
}

void ChrTags::readTagFile() {
    openTagFile((char *) "r");
    if (tagfp == NULL) {
        fprintf(stderr, "Could not open %s\n", tagFile);
        return;
    }
    char *buf = new char[BUFFER];
    char *name = new char[BUFFER];
    char **line = new char *[BUFFER];
    int numCols = 0;
    int pos = 0;
    float tagCount = 0;
    int dir = 0;
    int len = -1;
    char *chr1 = NULL;
    char *chr2 = NULL;
    int pos2 = 0;
    int dir2 = 0;
    int len2 = -1;
    appearentSize = 0;

    optimizedFlag = 1;
    int lastPos = -1000000000;

    while (fgets(buf, BUFFER, tagfp) != NULL) {
        //sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
        //
        //fprintf(stderr, "buf=%s\n", buf);
        split(buf, line, numCols, '\t');
        if (numCols < 5) continue;

        sscanf(line[2], "%d", &pos);
        sscanf(line[3], "%d", &dir);
        sscanf(line[4], "%f", &tagCount);

        //if (tagCount < 0.00001) continue;
        len = -1;
        if (numCols > 5) {
            sscanf(line[5], "%d", &len);
        }
        if (pos > appearentSize) appearentSize = (long long int) pos;

        if (pairedEndFlag && forceSingleReadFlag == 0) {
            if (numCols < 10) {
                fprintf(stderr, "!!! Expecting paired end tags (%s) !!!\n", tagFile);
                exit(0);
            }
            chr1 = (char *) chrNames->search(line[1]);
            if (chr1 == NULL) {
                char *cc = new char[strlen(line[1]) + 1];
                strcpy(cc, line[1]);
                chrNames->insert(cc, line[1]);
                chr1 = cc;
            }
            chr2 = (char *) chrNames->search(line[6]);
            if (chr2 == NULL) {
                char *cc = new char[strlen(line[6]) + 1];
                strcpy(cc, line[6]);
                chrNames->insert(cc, line[6]);
                chr2 = cc;
            }
            sscanf(line[7], "%d", &pos2);
            sscanf(line[8], "%d", &dir2);
            sscanf(line[9], "%d", &len2);
            addPETag(chr1, pos, dir, len, chr2, pos2, dir2, len2, tagCount);

        } else {
            if (revStrand) {
                if (dir == STRAND_POSITIVE) {
                    pos += len;
                    dir = STRAND_NEGATIVE;
                } else {
                    pos -= len;
                    dir = STRAND_POSITIVE;
                }
                if (pos < 0) pos = 0;
                if (pos > appearentSize) appearentSize = (long long int) pos;
            }
            addTag(pos, (char) dir, len, tagCount);
        }

        if (pos < lastPos) optimizedFlag = 0;
        lastPos = pos;
    }


    //fprintf(stderr, "appearentSize = %lld\n", appearentSize);

    delete[]buf;
    delete[]name;
    delete[]line;

    closeTagFile();

    optimizeTags();
    loaded = 1;
}

void ChrTags::addTag(int pos, char dir, int length, float value) {
    if (linkTag == NULL) {
        firstTag = new LinkedTag(pos, dir, length, value, NULL);
        linkTag = firstTag;
    } else {
        LinkedTag *next = new LinkedTag(pos, dir, length, value, NULL);
        linkTag->tag = next;
        linkTag = next;
    }
    numLinks++;
}

void ChrTags::addPETag(char *c1, int p1, char d1, int len1, char *c2, int p2,
                       char d2, int len2, float v) {
    if (linkPETag == NULL) {
        firstPETag = new LinkedPETag(c1, p1, d1, len1, c2, p2, d2, len2, v, NULL);
        linkPETag = firstPETag;
    } else {
        LinkedPETag *next = new LinkedPETag(c1, p1, d1, len1, c2, p2, d2, len2, v, NULL);
        linkPETag->tag = next;
        linkPETag = next;
    }
    numLinks++;
}

void ChrTags::print() {
    if (dontSAVE) {
        fprintf(stderr, "\t!!!! Tags were adjusted - shouldn't save them!!! Error in code\n");
        return;
    }
    if (singleFile) {
        if (tagfp == NULL) {
            //	tagfp = fopen();
        }
        for (int i = 0; i < totalPositions; i++) {
            tags[i].print(tagfp, chr);
        }
    } else {
        openTagFile((char *) "w");
        if (pairedEndFlag) {
            for (int i = 0; i < totalPositions; i++) {
                petags[i].print(tagfp);
            }
        } else {
            for (int i = 0; i < totalPositions; i++) {
                tags[i].print(tagfp, chr);
            }
        }
    }
}

void ChrTags::optimizeTagFile() {
    optimizeOverride = 1;
    readTagFile();
    adjustTags();
    dontSAVE = 0;
    print();
    freeTags();
    optimizeOverride = 0;
}

void ChrTags::optimizeTags() {

    if (pairedEndFlag && forceSingleReadFlag == 0) {
        optimizePETags();
        return;
    }

    if (numLinks == 0 || linkTag == NULL) return;

    totalPositions = 0;
    totalTags = 0;
    if (tags != NULL) {
        delete[]tags;
    }
    tags = new Tag[numLinks];

    LinkedTag *current = firstTag;
    LinkedTag *next = NULL;
    while (current != NULL) {
        tags[totalPositions].copy(current);
        totalTags += tags[totalPositions].v;
        next = current->tag;
        delete current;
        totalPositions++;
        current = next;
    }
    linkTag = NULL;
    firstTag = NULL;
    numLinks = 0;

    //fprintf(stderr, "%s  - optimizedFlag = %d\n", chr, optimizedFlag);
    if (optimizedFlag && !optimizeOverride) return;

    qsort(tags, totalPositions, sizeof(Tag), &cmpTags);
    int last = 0;
    totalTags = 0;
    if (totalPositions > 0) {
        totalTags = tags[0].v;
    } else {
        fprintf(stderr, "No tags in tag file!!!\n");
        return;
    }
    for (int i = 1; i < totalPositions; i++) {
        totalTags += tags[i].v;
        int j = i - 1;
        if (mCflag) {
            if (tags[i].p == tags[j].p && tags[i].d == tags[j].d) {
                double sum = tags[last].v * tags[last].len + tags[i].v * tags[i].len;
                double d = (double) (tags[last].len + tags[i].len);
                if (d > 0.0) sum /= d;
                tags[last].v = sum;
                tags[last].len = (int) d;
            } else {
                last++;
                if (last == i) continue;
                tags[last].copy(&(tags[i]));
            }
        } else {
            if (tags[i].p == tags[j].p && tags[i].d == tags[j].d && tags[i].len == tags[j].len) {
                tags[last].v += tags[i].v;
            } else {
                last++;
                if (last == i) continue;
                tags[last].copy(&(tags[i]));
            }
        }
    }
    totalPositions = last + 1;

    if (maxtbp > 0.0) {
        totalTags = 0.0;
        for (int i = 0; i < totalPositions; i++) {
            if (tags[i].v > maxtbp) tags[i].v = maxtbp;
            totalTags += tags[i].v;
        }
    }
    if (mintbp > 0.0) {
        totalTags = 0.0;
        for (int i = 0; i < totalPositions; i++) {
            if (tags[i].v < mintbp) tags[i].v = mintbp;
            totalTags += tags[i].v;
        }
    }
}

void ChrTags::optimizePETags() {

    if (numLinks == 0 || linkPETag == NULL) return;

    totalPositions = 0;
    totalTags = 0;
    if (petags != NULL) {
        delete[]petags;
    }
    petags = new PETag[numLinks];

    LinkedPETag *current = firstPETag;
    LinkedPETag *next = NULL;
    while (current != NULL) {
        petags[totalPositions].copy(current);
        totalTags += petags[totalPositions].v;
        next = current->tag;
        delete current;
        totalPositions++;
        current = next;
    }
    linkPETag = NULL;
    firstPETag = NULL;
    numLinks = 0;

    //fprintf(stderr, "%s  - optimizedFlag = %d\n", chr, optimizedFlag);
    if (optimizedFlag && !optimizeOverride) return;
    //fprintf(stderr, "Optimizing\n");
    qsort(petags, totalPositions, sizeof(PETag), &cmpPETags);

    int last = 0;
    totalTags = 0;
    if (totalPositions > 0) {
        totalTags = petags[0].v;
    }
    for (int i = 1; i < totalPositions; i++) {
        totalTags += petags[i].v;
        int j = i - 1;
        if (j >= 0 && petags[i].p1 == petags[j].p1 && petags[i].d1 == petags[j].d1 &&
            strcmp(petags[i].chr2, petags[j].chr2) == 0 && petags[i].p2 == petags[j].p2 &&
            petags[i].d2 == petags[j].d2) {
            petags[last].v += petags[i].v;
        } else {
            last++;
            if (last == i) continue;
            petags[last].copy(&(petags[i]));
        }
    }
    totalPositions = last + 1;

    if (maxtbp > 0.0) {
        totalTags = 0.0;
        for (int i = 0; i < totalPositions; i++) {
            if (petags[i].v > maxtbp) petags[i].v = maxtbp;
            totalTags += petags[i].v;
        }
    }
    if (mintbp > 0.0) {
        totalTags = 0.0;
        for (int i = 0; i < totalPositions; i++) {
            if (petags[i].v < mintbp) petags[i].v = mintbp;
            totalTags += petags[i].v;
        }
    }
}

int cmpTags(const void *a, const void *b) {
    int ap = ((Tag *) a)->p;
    int bp = ((Tag *) b)->p;
    if (ap < bp) return -1;
    if (ap > bp) return 1;
    char ad = ((Tag *) a)->d;
    char bd = ((Tag *) b)->d;
    if (ad < bd) return -1;
    if (ad > bd) return 1;
    char al = ((Tag *) a)->len;
    char bl = ((Tag *) b)->len;
    if (al < bl) return -1;
    if (al > bl) return 1;
    return 0;
}

int cmpPETags(const void *a, const void *b) {
    int ap = ((PETag *) a)->p1;
    int bp = ((PETag *) b)->p1;
    if (ap < bp) return -1;
    if (ap > bp) return 1;
    char ad = ((PETag *) a)->d1;
    char bd = ((PETag *) b)->d1;
    if (ad < bd) return -1;
    if (ad > bd) return 1;
    char *ac = ((PETag *) a)->chr2;
    char *bc = ((PETag *) b)->chr2;
    int x = strcmp(ac, bc);
    if (x < 0) return -1;
    if (x > 0) return 1;
    ap = ((PETag *) a)->p2;
    bp = ((PETag *) b)->p2;
    if (ap < bp) return -1;
    if (ap > bp) return 1;
    ap = ((PETag *) a)->d2;
    bp = ((PETag *) b)->d2;
    if (ap < bp) return -1;
    if (ap > bp) return 1;
    return 0;
}

int ChrTags::getPETagDistribution(double *sameStrand, double *diffStrand, int windowSize,
                                  double *largeWindow, int resolution, int largeLength) {

    loadTags();

    int rvLen = 0;
    if (totalPositions > 0) {
        rvLen = petags[totalPositions - 1].p1 - petags[0].p1;
    }

    int halfWindow = (windowSize) / 2;
    for (int i = 0; i < totalPositions; i++) {

        float v = petags[i].v;

        if (strcmp(petags[i].chr1, petags[i].chr2) != 0) {
            largeWindow[largeLength - 1] += v;
        } else {
            int diff = petags[i].p2 - petags[i].p1;
            int bin = (int) (fabs(((double) diff) / ((double) resolution)) + 0.5);
            if (bin > largeLength - 2) {
                largeWindow[largeLength - 2] += v;
            } else {
                largeWindow[bin] += v;
            }

            //if (diff > halfWindow || diff < -1*halfWindow) continue;
            if (halfWindow - diff < 0 || halfWindow - diff >= windowSize) continue;
            if (halfWindow + diff < 0 || halfWindow + diff >= windowSize) continue;

            if (petags[i].d1 == petags[i].d2) {
                if (petags[i].d1 == 0) {
                    sameStrand[halfWindow + diff] += v;
                } else {
                    sameStrand[halfWindow - diff] += v;
                }
            } else {
                if (petags[i].d1 == 0) {
                    diffStrand[halfWindow + diff] += v;
                } else {
                    diffStrand[halfWindow - diff] += v;
                }
            }
        }
    }
    freeTags();
    forceSingleReadFlag = 0;
    return rvLen;
}


void ChrTags::setMaxTBP(float max) {
    maxtbp = max;
}

void ChrTags::setMinTBP(float min) {
    mintbp = min;
}

void ChrTags::setTagAdjust(int dist) {
    tagAdjust = dist;
}

void ChrTags::adjustTags() {
    if (totalPositions < 1) return;
    if (adjusted) return;
//fprintf(stderr, "maxtbp=%f\tadj=%d\n",maxtbp,tagAdjust);
    if (pairedEndFlag && forceSingleReadFlag == 0) {
        if (maxtbp > FLOAT_ZERO || mintbp > FLOAT_ZERO) {
            totalTags = 0.0;
            for (int i = 0; i < totalPositions; i++) {
                if (petags[i].v > maxtbp) {
                    petags[i].v = maxtbp;
                }
                if (petags[i].v < mintbp) {
                    petags[i].v = mintbp;
                }
                totalTags += petags[i].v;
            }
        }
        return;
    } else {
        int minPosition = tags[0].p;
        int maxPosition = tags[totalPositions - 1].p;
        if (maxtbp > FLOAT_ZERO || mintbp > FLOAT_ZERO) {
            totalTags = 0;
            for (int i = 0; i < totalPositions; i++) {
                if (tags[i].v > maxtbp) {
                    tags[i].v = maxtbp;
                }
                if (tags[i].v < mintbp) {
                    tags[i].v = mintbp;
                }
                totalTags += tags[i].v;
            }
        }
        dontSAVE = 1;
        if (tagAdjust != 0) {
            for (int i = 0; i < totalPositions; i++) {
                if (tags[i].d == 0) {
                    tags[i].p += tagAdjust;
                } else {
                    tags[i].p -= tagAdjust;
                }
                if (tags[i].p > maxPosition) {
                    tags[i].p = maxPosition;
                }
                if (tags[i].p < minPosition) {
                    tags[i].p = minPosition;
                }
            }
            qsort(tags, totalPositions, sizeof(Tag), &cmpTags);
        }
    }
    adjusted = 1;
}


// negative value for minDist is used to indicate region based peak finding
void ChrTags::findPutativePeaks(PeakLibrary *putativePeaks, int peakSize, int minDist, char strand, double minCount) {

    doubleIndex *tagCount = new doubleIndex[totalPositions];
    int *centers = new int[totalPositions];
    char *mask = new char[totalPositions];
    char *pname = new char[100];
    int pid = 1;

    int regionFlag = 0;
    if (minDist < 0) {
        minDist *= -1;
        regionFlag = 1;
    }

    int lastPosIndex = -1;
    int lastNegIndex = -1;
    int lastPosFinish = -1;
    int lastNegFinish = -1;
    int lastPosPosition = -1;
    int maxEND = 0;
    if (totalPositions > 0) {
        maxEND = tags[totalPositions - 1].p;
    }

    fprintf(stderr, "\t\tFinding peaks on %s (minCount=%.1lf), total tags positions = %d\n",
            chr, minCount, totalPositions);

    for (int i = 0; i < totalPositions; i++) {
        float v = tags[i].v;
        char d = tags[i].d;
        tagCount[i].index = i;
        int p = tags[i].p;

        if (strand == STRAND_BOTH) {
            int start = i + 1;
            if (lastPosFinish >= i) {
                if (lastPosPosition < p) {
                    tagCount[i].v = tagCount[i - 1].v - tags[i - 1].v;
                    tagCount[i].vp = tagCount[i - 1].vp - tags[i - 1].v * (double) ((tags[i - 1].p) - p);
                    start = lastPosIndex + 1;
                } else {
                    //some tags could be at the same location...
                    tagCount[i].v = tagCount[i - 1].v;
                    tagCount[i].vp = tagCount[i - 1].vp;
                    continue;
                }
            } else {
                tagCount[i].v = v;
                //tagCount[i].vp = v*(double)(p-p); i.e. 0
                tagCount[i].vp = 0;
            }
            for (int j = start; j < totalPositions; j++) {
                int diff = tags[j].p - p;
                if (diff > peakSize) break;
                tagCount[i].v += tags[j].v;
                tagCount[i].vp += tags[j].v * (double) ((tags[j].p) - p);
                lastPosFinish = i;
            }
            lastPosPosition = p;
        } else if (strand == STRAND_SEPARATE) {
            int start = i + 1;
            if (d == 0) {
                if (lastPosFinish >= i) {
                    if (tags[lastPosFinish].p == p) {
                        tagCount[i].v = tagCount[lastPosIndex].v;
                        tagCount[i].vp = tagCount[lastPosIndex].vp;
                        lastPosIndex = i;
                        continue;
                    } else {
                        tagCount[i].v = tagCount[lastPosIndex].v - tags[lastPosIndex].v;
                        tagCount[i].vp = tagCount[lastPosIndex].vp
                                         - tags[lastPosIndex].v * (double) (tags[lastPosIndex].p - p);
                        start = lastPosFinish + 1;
                    }
                } else {
                    tagCount[i].v = v;
                    tagCount[i].vp = 0;
                    //tagCount[i].vp = v*(double)(p-p);
                }
                lastPosIndex = i;
            } else {
                if (lastNegFinish >= i) {
                    if (tags[lastNegFinish].p == p) {
                        tagCount[i].v = tagCount[lastNegIndex].v;
                        tagCount[i].vp = tagCount[lastNegIndex].vp;
                        lastNegIndex = i;
                        continue;
                    } else {
                        tagCount[i].v = tagCount[lastNegIndex].v - tags[lastNegIndex].v;
                        tagCount[i].vp = tagCount[lastNegIndex].vp
                                         - tags[lastNegIndex].v * (double) (tags[lastNegIndex].p - p);
                        start = lastNegFinish + 1;
                    }
                } else {
                    tagCount[i].v = v;
                    tagCount[i].vp = 0;
                    //tagCount[i].vp = v*(double)(p-p); i.e. 0
                }
                lastNegIndex = i;
            }
            for (int j = start; j < totalPositions; j++) {
                int diff = tags[j].p - p;
                if (diff > peakSize) break;
                if (d != tags[j].d) continue;
                tagCount[i].v += tags[j].v;
                tagCount[i].vp += tags[j].v * (double) (tags[j].p - p);
                if (d == 0) {
                    lastPosFinish = i;
                } else {
                    lastNegFinish = i;
                }
            }
        }
    }

    for (int i = 0; i < totalPositions; i++) {
        mask[i] = 0;
        if (tagCount[i].v < 0.000001) tagCount[i].v = 1.0;
        tagCount[i].position = tags[i].p + (int) floor(tagCount[i].vp / tagCount[i].v);
        if (tagCount[i].position < 0) {
            fprintf(stderr, "%d\t%d\t%d\t%lf\t%lf\n", i, tags[i].p, tagCount[i].position, tagCount[i].vp,
                    tagCount[i].v);
        }
        centers[i] = tagCount[i].position;
    }

    qsort(tagCount, totalPositions, sizeof(doubleIndex), &cmpDoubleIndex);


    int d = minDist;
    if (d < 1) d = 1;
    int expectedNumPeaks = appearentSize / d * 2;
    if (expectedNumPeaks < 10) expectedNumPeaks = 10;

    int halfSize = (peakSize) / 2;

    for (int i = 0; i < totalPositions; i++) {

        unsigned int index = tagCount[i].index;
        double value = tagCount[i].v;
        int center = tagCount[i].position;

        if (value < minCount) break;

        char d = tags[index].d;
//fprintf(stderr, "%d\n", d);
        if (strand == STRAND_BOTH) d = 0;

        if (!mask[index]) {
            int start = center - halfSize;
            int end = start + peakSize;
            if (start < 1) {
                start = 1;
            }
            if (end > maxEND) {
                end = maxEND;
            }

            sprintf(pname, "%s-%d", chr, pid++);
            //sprintf(pname, "%d",rand());
            //fprintf(stderr, "%s\t%s\t%d\t%d\t%d\t%d\t%f\n",pname,chr,start,end,center,d,value);
            putativePeaks->addPeak(pname, chr, start, end, center, d, value, 0, NULL, -1, 0);
        }
        if (regionFlag == 0 || !mask[index]) {
            for (int j = index; j < totalPositions; j++) {
                if (centers[j] - center < minDist) {
                    if (strand == STRAND_BOTH || d == tags[j].d) {
                        mask[j] = 1;
                    }
                } else {
                    break;
                }
            }
            for (int j = index - 1; j >= 0; j--) {
                if (center - centers[j] < minDist) {
                    if (strand == STRAND_BOTH || d == tags[j].d) {
                        mask[j] = 1;
                    }
                } else {
                    break;
                }
            }
        }
    }

    delete[]tagCount;
    delete[]centers;
    delete[]mask;
}

int cmpDoubleIndex(const void *a, const void *b) {
    double av = ((doubleIndex *) a)->v;
    double bv = ((doubleIndex *) b)->v;
    if (av < bv) return 1;
    if (av > bv) return -1;
    return 0;
}

void ChrTags::getTagLengthDistribution(double *dist, int max) {
    forceSingleReadFlag = 1;
    loadTags();
    for (int i = 0; i < totalPositions; i++) {
        int len = tags[i].len;
        if (len < 0) len = 0;
        if (len >= max) {
            len = max - 1;
        }
        dist[len] += tags[i].v;
    }
    freeTags();
    forceSingleReadFlag = 0;
}

void ChrTags::getTagCountDistribution(double *dist, int max, int scaleFactor) {
    //forceSingleReadFlag = 1;
    loadTags();
    for (int i = 0; i < totalPositions; i++) {
        int v = 0;
        if (pairedEndFlag) {
            v = (int) floor((double) petags[i].v * scaleFactor);
        } else {
            v = (int) floor((double) tags[i].v * scaleFactor);
        }
        if (v >= max) {
            v = max - 1;
        }
        dist[v]++;
    }
    freeTags();
    //forceSingleReadFlag = 0;
}


void ChrTags::readAndSave() {
    loadTags();
    dontSAVE = 0;
    print();
    freeTags();
}


// class UniqMapChrs -------------------------------------------------------------------


// class Tag ---------------------------------------------------------------

int Tag::precision = TAG_VALUE_RESOLUTION;
int PETag::precision = TAG_VALUE_RESOLUTION;

Tag::Tag() {
    p = 0;
    len = -1;
    v = 0.0;
    d = 0;
}

PETag::PETag() {
    init();
}

PETag::PETag(char *n, char *c, int np, char nd, float nv, int nlen) {
    init();
    name = n;
    chr1 = c;
    p1 = np;
    d1 = nd;
    v = nv;
    len1 = nlen;
}

void PETag::init() {
    name = NULL;
    chr1 = NULL;
    p1 = -1;
    d1 = 0;
    len1 = -1;
    chr2 = NULL;
    p2 = -1;
    d2 = 0;
    len2 = -1;
    v = 0.0;
}

PETag::~PETag() {
}

void PETag::print(FILE *fp) {
    print(fp, 0);
}

void PETag::copy(PETag *src) {
    chr1 = src->chr1;
    p1 = src->p1;
    d1 = src->d1;
    len1 = src->len1;
    chr2 = src->chr2;
    p2 = src->p2;
    d2 = src->d2;
    len2 = src->len2;
    v = src->v;
}

void PETag::print(FILE *fp, int revFlag) {
    /*if (percision == 1) {
		if (v <  0.05) return;
	} else if (percision == 2) {
		if (v <  0.005) return;
	} else if (percision == 3) {
		if (v <  0.005) return;
	}*/
    if (name != NULL) fprintf(fp, "%s", name);
    if (revFlag) {
        fprintf(fp, "\t%s\t%d\t%d", chr2, p2, d2);
    } else {
        fprintf(fp, "\t%s\t%d\t%d", chr1, p1, d1);
    }
    if (precision == 1) {
        fprintf(fp, "\t%.1f", v);
    } else if (precision == 2) {
        fprintf(fp, "\t%.2f", v);
    } else {
        fprintf(fp, "\t%.3f", v);
    }
    if (revFlag) {
        fprintf(fp, "\t%d\t%s\t%d\t%d\t%d\n", len2, chr1, p1, d1, len1);
    } else {
        fprintf(fp, "\t%d\t%s\t%d\t%d\t%d\n", len1, chr2, p2, d2, len2);
    }
}

void Tag::copy(Tag *src) {
    p = src->p;
    len = src->len;
    v = src->v;
    d = src->d;
}

void Tag::print(FILE *fp, char *chr) {
    if (precision == 1) {
        fprintf(fp, "\t%s\t%d\t%d\t%.1f\t%d\n", chr, p, d, v, len);
    } else if (precision == 2) {
        //if (v <  0.005) return;
        fprintf(fp, "\t%s\t%d\t%d\t%.2f\t%d\n", chr, p, d, v, len);
    } else {
        //if (v <  0.05) return;
        fprintf(fp, "\t%s\t%d\t%d\t%.3f\t%d\n", chr, p, d, v, len);
    }
}

LinkedTag::LinkedTag(int pos, char dir, int length, float value, LinkedTag *link) {
    p = pos;
    v = value;
    d = dir;
    len = length;
    tag = link;
}

LinkedPETag::LinkedPETag(char *nc1, int np1, char nd1, int nlen1, char *nc2, int np2, char nd2,
                         int nlen2, float nv, LinkedPETag *link) {
    chr1 = nc1;
    p1 = np1;
    d1 = nd1;
    len1 = nlen1;
    chr2 = nc2;
    p2 = np2;
    d2 = nd2;
    len2 = nlen2;
    v = nv;
    tag = link;
}


int checkInt(char *str) {
    int bad = 0;
    int len = strlen(str);
    for (int i = 0; i < len; i++) {
        if (str[i] < 45) bad = 1;
        if (str[i] == 46 || str[i] == 47) bad = 1;
        if (str[i] > 57) bad = 1;
    }
    return bad;
}

int checkStrand(char *str) {
    int value = -1;
    if (strcmp(str, "0") == 0) {
        value = 0;
    } else if (strcmp(str, "1") == 0) {
        value = 1;
    } else if (strcmp(str, "+") == 0) {
        value = 0;
    } else if (strcmp(str, ".") == 0) {
        value = 0;
    } else if (strcmp(str, "-") == 0) {
        value = 1;
    } else {
        value = -1;
    }
    return value;
}
