#include <iostream>
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

std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

// class PeakFinder ---------------------------------------------------------------------

PeakFinder::PeakFinder() {
    outputFileName = "";
    peakSize = 0;
    localSize = 10000;
    inputSize = 0;
    tagThresh = 0;
    minTagThresh = -1.0;
    minDist = 0;
    totalTags = 0.0;
    tagsUsedForClustering = 0.0;
    maxtbp = 0.0;
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

    style = PEAK_STYLE_HISTONE;

    inputFold = 4.0;
    localFold = 4.0;
    clonalFold = 2.0;
    numPeaks = 0;

    tags = NULL;
    input = NULL;
}

PeakFinder::~PeakFinder() {
    delete[]extraheader;
    delete[]fdrTable;
    delete[]poissonTable;
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

void PeakFinder::setOutputFile(std::string fname) {
    outputFileName = std::move(fname);
}

void PeakFinder::determineMaxTBP() {
    maxtbp = 1.0;
    fprintf(stderr, "\tTotal Tags = %.1lf\n", totalTags);
    fprintf(stderr, "\tTags per bp = %.6lf\n", tbp);
    if (tbp > tbpThreshold) maxtbp = floor(tbp / tbpThreshold);
    if (tbpInput > tbpThreshold) floor(tbp / tbpThreshold);
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
        //tagAdjustInput = input->fragmentLength;
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
    tagsUsedForClustering = tags->getAdjustedTagTotal();

    {
        PeakLibrary *putativePeaks = tags->findPutativePeaks(peakSize, curMinDist, strand, minTagThresh);

        fprintf(stderr, "\t\tTags Used for cluster (less clonal tags) = %.1lf / %.1lf\n",
                tagsUsedForClustering, totalTags);

        // remove peaks not meeting fdr levels or poisson p-value or absolute tag thresh
        approxFdrTable();
        addHeader((char *) "findPeaks Score");
        filteredPeaks = filterPeaks(putativePeaks);
        delete putativePeaks;
    }

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

    //normalize tag counts
    tags->setMaxTBP(0);

    filteredPeaks->setPeakTagSizeFixed(0, 0);
    //filteredPeaks->setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);

    filteredPeaks->setDefaultPeakOrder();
    filteredPeaks->tagsInPeaks = 0.0;

    //int tagsIndex = filteredPeaks->addTagLibrary(tags);
    //Doubletable* expTags = filteredPeaks->countPeakTags(tagsIndex,0,0,strand,COUNT_MODE_TOTAL);
    auto expTags = filteredPeaks->countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, excludePeaks);

    for (int i = 0; i < filteredPeaks->numPeaks; i++) {
        filteredPeaks->peakOrder[i]->v = expTags[filteredPeaks->peakOrder[i]->name];
        filteredPeaks->tagsInPeaks += filteredPeaks->peakOrder[i]->v;
    }
    filteredPeaks->normalizePeakScore(normTotal / tags->totalTags);

    tagsInPeaks = filteredPeaks->tagsInPeaks;

    numPeaks = filteredPeaks->numPeaks;
    fprintf(stderr, "\tTotal Peaks identified = %d\n", numPeaks);
    filteredPeaks->setDefaultPeakOrder();

    return filteredPeaks;
}

void PeakLibrary::normalizePeakScore(float normFactor) {
    for (auto& it: peaks) {
        it.second->v *= normFactor;
    }
}

PeakLibrary *PeakFinder::filterPeaks(PeakLibrary *putativePeaks) {
    double *peakHeights = new double[fdrSize];
    for (int i = 0; i < fdrSize; i++) {
        peakHeights[i] = 0;
    }

    for (auto& it: putativePeaks->peaks) {
        auto p = it.second;
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
    tagsInPeaks = 0.0;

    double threshold2Use = fdrThresh;
    if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
        threshold2Use = poissonThresh;
    } else if (filterMode == PEAKFINDER_FILTER_MODE_THRESH) {
        threshold2Use = tagThresh;
    }

    for (auto& it: putativePeaks->peaks) {
        auto p = it.second;
        if (p->v >= threshold2Use - 000000.1) {
            tagsInPeaks += (double) p->v;

            static auto fmt = "%f";
            int sz = std::snprintf(nullptr, 0, fmt, p->v);
            auto buf = new char[sz+1];
            std::snprintf(buf, sz+1, fmt, p->v);
            p->addData(buf);

            goodPeaks->addPeak(p);
            numGood++;
        }
    }
    delete[]peakHeights;
    fprintf(stderr, "\t%d peaks passed threshold\n", numGood);
    goodPeaks->sortChr();

    return goodPeaks;
}

void PeakFinder::print(FILE *fp) {
    if (fp == NULL) return;

    fprintf(fp, "# HOMER Peaks\n");
    fprintf(fp, "# Peak finding parameters:\n");
    fprintf(fp, "#\n");

    fprintf(fp, "# total peaks = %d\n", numPeaks);
    fprintf(fp, "# peak size = %d\n", peakSize);
    if (strand == STRAND_BOTH) {
        fprintf(fp, "# peaks found using tags on both strands\n");
    } else {
        fprintf(fp, "# peaks found separately for each strand\n");
    }
    fprintf(fp, "# minimum distance between peaks = %d\n", minDist);
    //fprintf(fp, "# genome = %s\n",genome);
    fprintf(fp, "# genome size = %.0lld\n", gsize);
    fprintf(fp, "# Total tags = %.1lf\n", totalTags);
    fprintf(fp, "# Total tags in peaks = %.1lf\n", tagsInPeaks);
    double ratio = tagsInPeaks / totalTags * 100.0;
    fprintf(fp, "# Approximate IP efficiency = %.2lf%%\n", ratio);
    fprintf(fp, "# tags per bp = %.6f\n", tbp);
    fprintf(fp, "# expected tags per peak = %.3f\n", tpp);
    fprintf(fp, "# maximum tags considered per bp = %.1f\n", maxtbp);
    fprintf(fp, "# effective number of tags used for normalization = %.1f\n", normTotal);
    if (regionFlag) {
        fprintf(fp, "# Individual peaks have been stitched together into variable length regions\n");
    }
    if (filterMode == PEAKFINDER_FILTER_MODE_FDR) {
        fprintf(fp, "# FDR rate threshold = %.9lf\n", fdr);
        fprintf(fp, "# FDR effective poisson threshold = %le\n", poisson);
        fprintf(fp, "# FDR tag threshold = %.1f\n", fdrThresh);
    } else if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
        fprintf(fp, "# Poisson p-value threshold = %le\n", poisson);
        fprintf(fp, "# Poisson tag threshold = %.1f\n", poissonThresh);
    } else if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
        fprintf(fp, "# Manual tag threshold = %.1f\n", tagThresh);
    }
    if (input != NULL && (inputFold > 0 || poissonInput < 1.0)) {
        fprintf(fp, "#\n");
        if (inputFold > 0) fprintf(fp, "# Fold over input required = %.2lf\n", inputFold);
        if (poissonInput < 1.0) fprintf(fp, "# Poisson p-value over input required = %.2le\n", poissonInput);
    }
    if (localFold > 0 || poissonLocal < 1.0) {
        fprintf(fp, "#\n");
        fprintf(fp, "# size of region used for local filtering = %d\n", localSize);
        if (localFold > 0) fprintf(fp, "# Fold over local region required = %.2lf\n", localFold);
        if (poissonLocal < 1.0) fprintf(fp, "# Poisson p-value over local region required = %.2le\n", poissonLocal);
    }
    if (clonalFold > 0) {
        fprintf(fp, "#\n");
        fprintf(fp, "# Maximum fold under expected unique positions for tags = %.2lf\n", clonalFold);
    }

//    fprintf(fp, "#\n# cmd = %s\n", cmd);
    fprintf(fp, "#\n# Column Headers:\n");
    fprintf(fp, "#PeakID\tchr\tstart\tend\tstrand");
    fprintf(fp, "\tNormalized Tag Count");

    if (regionFlag) {
        fprintf(fp, "\tregion size");
    } else {
        fprintf(fp, "\tNot used");
    }
    if (extraheader != NULL) {
        fprintf(fp, "%s", extraheader);
    }
    fprintf(fp, "\n");
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
       //exit(0);

//    cum = 0.0;
//
//    for (int i = 0; i < fdrSize; i++) {
//        double lp = logPoisson(i, tpp);
//        lp = exp(lp);
//        //fprintf(stderr, "\t%d\t%.9lf\t%.9lf\t%.3f\n",i, lp, 1-cum,numTests*(1-cum));
//        double poissonCum = 1 - cum;
//        fdrTable[i] = numTests * poissonCum;
//        poissonTable[i] = poissonCum;
//        fprintf(stderr, "\t%d\t%lf\t%le\t%le\t%le\n", i, fdrTable[i], numTests, poissonCum, lp);
//        if (poissonThresh < 0.01 && poissonCum < poisson) {
//            poissonThresh = (float) i;
//            if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
//                fprintf(stderr, "\tPoisson Threshold set at %d tags\n", i);
//            }
//        }
//        cum += lp;
//    }

}


// class PeakLibrary --------------------------------------------------------------------

PeakLibrary::PeakLibrary() {
    initialize(PEAK_LIBRARY_DEFAULT_SIZE);
}

PeakLibrary::PeakLibrary(int expectedNumPeaks) {
    initialize(expectedNumPeaks);
}

void PeakLibrary::initialize(int expectedNumPeaks) {
    duplicates.reserve(expectedNumPeaks / 5);
    numPeaks = 0;
    tagsInPeaks = 0.0;
    fixedFlag = 0;
    duplicateWarningFlag = 1;
}

PeakLibrary::~PeakLibrary() {
    for(auto& it: chrs)
        if (it.second != NULL)
            delete it.second;
    for(auto& it: peaks)
        if (it.second != NULL)
            delete it.second;
}

void PeakLibrary::sortChr() {
    numPeaks = 0;
    for (auto& it: chrs)
    {
        auto ct = it.second;
        ct->sort();
        numPeaks += ct->peaks.size();
    }
}

Peak *PeakLibrary::addPeak(const std::string &name, const std::string &chr, int start, int end, int midpoint, char dir,
                           float value, float ratio, char *extraData, int mappability, unsigned int priority) {
    Peak *oldPeak = NULL;
    Peak *p = NULL;
    if (!name.empty()) {
        auto it = peaks.find(name);
        if (it != peaks.end())
            oldPeak = it->second;
    }

    static int warningIssued = 0;
    if (oldPeak != NULL) {
        auto it = duplicates.find(name);
        int dupCount;
        if (it == duplicates.end())
            dupCount = 1;
        else
            dupCount = it->second;
        dupCount++;
        duplicates[name] = dupCount;

        auto newname = name + "--" + std::to_string(dupCount);
        if (priority < 1) {
            if (warningIssued == 0 && duplicateWarningFlag) {
                fprintf(stderr, "\tDuplicate peak name (%s) - this could potentially cause problems\n", name.c_str());
                fprintf(stderr, "\t\tSometimes unavoidable for BED/2DBED formats\n");
                fprintf(stderr, "\t\tNew name for this peak is %s\n", newname.c_str());
            } else if (warningIssued % 1000 == 0 && duplicateWarningFlag) {
                fprintf(stderr, "\t\tWarning over %d peaks with duplicate names\n", warningIssued);
            }
            warningIssued++;
        }
        auto it2 = peaks.find(newname);
        if (it2 != peaks.end() && duplicateWarningFlag) {
            fprintf(stderr, "There's a problem!!! %s - %d %s\n", newname.c_str(), dupCount, name.c_str());
        }
        p = new Peak(newname, name, chr, start, end, midpoint, dir, value, ratio, extraData, mappability, priority);
    } else {
        p = new Peak(name, "", chr, start, end, midpoint, dir, value, ratio, extraData, mappability, priority);
    }
    peaks[p->name] = p;
    auto it = chrs.find(chr);
    if (it == chrs.end())
        chrs[chr] = new ChrPeaks();
    chrs[chr]->addPeak(p);
    numPeaks++;
    return p;
}

void PeakLibrary::addPeak(Peak *p) {
    auto& n = p->name;
    if (!p->ogname.empty()) {
        n = p->ogname;
    }
    addPeak(n, p->chr, p->start, p->end, p->refPos, p->strand, p->v, p->focusRatio, p->data,
            p->uniqMap, p->priority);
}

void PeakLibrary::print(FILE *fp) {
//    auto keys = std::vector<std::pair<std::string, Peak*>>(peaks.begin(), peaks.end());
//    std::sort(keys.begin(), keys.end(), [&](const auto& a, const auto& b) {mus
//        Peak *p1 = a.second;
//        Peak *p2 = b.second;
//        char *c1 = p1->chr;
//        char *c2 = p2->chr;
//        int c = chrcmp((void *) &c1, (void *) &c2);
//        if (c < 0) return true;
//        if (p1->start < p2->start) return true;
//        if (p1->start > p2->start) return false;
//        if (p1->v < p2->v) return true;
//        if (p1->v > p2->v) return false;
//        return false;
//    });

    for (auto& it: peaks) {
        it.second->print(fp);
    }
}

std::unordered_map<std::string, double>
PeakLibrary::countPeakTagsLowMemory(TagLibrary *tags, char direction, int mode, PeakLibrary *excludePeaks) {
    std::unordered_map<std::string, double> results;
    for (auto& it: peaks)
        results[it.first] = 0.0;

    for (auto& it: chrs) {
        ChrPeaks *cp = it.second;
        ChrTags *ct = tags->chrs[it.first];

        ChrPeaks *cte = NULL;
        if (excludePeaks != NULL) {
            auto it2 = excludePeaks->chrs.find(it.first);
            if (it2 != excludePeaks->chrs.end()) {
                cte = it2->second;
            }
        }
        if (ct == NULL || cp == NULL) {
            continue;
        }
        cp->countPeakTagsLowMemory(results, ct, direction, mode, cte);
    }
    return results;
}

void PeakLibrary::setPeakTagSizeRefPos(int newOffset, int startOffset, int endOffset) {
//fprintf(stderr, "%d\t%d\t%d\n", newOffset, startOffset, endOffset);
    for (auto& it: peaks)
        it.second->setPeakTagSizeRefPos(newOffset, startOffset, endOffset);
    sortChr();
}

void PeakLibrary::setPeakTagSizeFixed(int startOffset, int endOffset) {
    fixedFlag = 1;
    for (auto& it: peaks)
        it.second->setPeakTagSizeFixed(startOffset, endOffset);
    sortChr();
}

void PeakLibrary::setDefaultPeakOrder() {
    if (numPeaks < 1) return;
    peakOrder = std::vector<Peak*>();
    peakOrder.resize(numPeaks);

    size_t i = 0;
    for (auto& it: peaks) {
        peakOrder[i] = it.second;
        i += 1;
    }
}

PeakLibrary *PeakLibrary::filterClonalPeaks(TagLibrary *tags, int peakSize,
                                            double threshold, int mode, char strand) {
    int halfPeakSize = (peakSize) / 2;

    float currentMaxTBP = tags->maxtbp;
    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfPeakSize, halfPeakSize);

    tags->setMaxTBP(0);
    //int tagsIndex = addTagLibrary(tags);
    auto expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    //Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);

    tags->setMaxTBP(1);
    //int posIndex = addTagLibrary(tags);
    auto posTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);

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

    int numGood = 0;
    int totalChecked = 0;
    goodPeaks->tagsInPeaks = 0.0;
    for (auto& it: peaks) {
        double pt = expTags[it.first];
        double ct = posTags[it.first];
        if (pt < EMPTY_DOUBLE_CHECK || ct < EMPTY_DOUBLE_CHECK) {
            continue;
        }
        totalChecked++;

        int index = (int) (pt / avgTagsPerPosition);
        if (index > expectedSize - 1) index = expectedSize - 1;
        double et = expected[index];
        double fold = et / ct;
        if (ct > 0 && fold < threshold) {
            goodPeaks->tagsInPeaks += pt;
            numGood++;
            Peak *p = it.second;

            static auto fmt = "%.2lf";
            int sz = std::snprintf(nullptr, 0, fmt, fold);
            auto buf = new char[sz+1];
            std::snprintf(buf, sz+1, fmt, fold);
            p->addData(buf);

            goodPeaks->addPeak(p);
        }
    }
    delete[] expected;

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
    int halfPeakSize = (peakSize) / 2;
    int halfLocalSize = (localSize) / 2;

    //int tagsIndex = addTagLibrary(tags);

    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfPeakSize, halfPeakSize);
    auto expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    setPeakTagSizeRefPos(NULL_OFFSET, -1 * halfLocalSize, halfLocalSize);
    auto localTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);

    //Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);
    //Doubletable* localTags = countPeakTags(tagsIndex,-1*halfLocalSize,halfLocalSize,strand,COUNT_MODE_TOTAL);

    double peakLength = (double) peakSize;
    double localLength = (double) (localSize - peakSize);
    double localAdjustFactor = peakLength / localLength;

    double tagPseudoCount = 0.0;
    if (tags->tbp > 0.0) tagPseudoCount = 0.5;
    //if (tags->tbp > 0.0) tagPseudoCount = localLength*tags->tbp;

    PeakLibrary *goodPeaks = new PeakLibrary();
    int numGood = 0;
    goodPeaks->tagsInPeaks = 0.0;
    int totalChecked = 0;
    for (auto& it: peaks) {
        double pt = expTags[it.first];
        double lt = localTags[it.first];
        if (pt < EMPTY_DOUBLE_CHECK || lt < EMPTY_DOUBLE_CHECK) {
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
            Peak *p = it.second;

            static auto fmt = "%.2lf\t%.2le";
            int sz = std::snprintf(nullptr, 0, fmt, fold, pp);
            auto buf = new char[sz+1];
            std::snprintf(buf, sz+1, fmt, fold, pp);
            p->addData(buf);

            goodPeaks->addPeak(p);
        }
    }
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
    auto expTags = countPeakTagsLowMemory(tags, strand, COUNT_MODE_TOTAL, NULL);
    auto inputTags = countPeakTagsLowMemory(input, strand, COUNT_MODE_TOTAL, NULL);

    double inputPseudoCount = 0.0;
    double tagPseudoCount = 0.0;
    if (input->tbp > 0.0) inputPseudoCount = 0.5;
    if (tags->tbp > 0.0) tagPseudoCount = 0.5;

    double minTotal = tags->totalTags;
    if (input->totalTags < tags->totalTags) minTotal = input->totalTags;

    double tagsRatio = minTotal / tags->totalTags;
    double inputRatio = minTotal / input->totalTags;

    int numGood = 0;
    int totalChecked = 0;

    auto *goodPeaks = new PeakLibrary();
    goodPeaks->tagsInPeaks = 0.0;
    for (auto& it: expTags) {
        double tp = it.second;
        double ip = inputTags[it.first];
        if (tp < EMPTY_DOUBLE_CHECK || ip < EMPTY_DOUBLE_CHECK) {
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
        double pp = exp(ilogCumulativePoisson(tpInt, ipn));
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
            Peak *p = peaks[it.first];
            if (strFlag) {
                static auto fmt = "%.1lf\t%.1lf\t%.2lf\t%.2le";
                int sz = std::snprintf(nullptr, 0, fmt, tpn, ipn, foldchange, pp);
                auto buf = new char[sz+1];
                std::snprintf(buf, sz+1, fmt, tpn, ipn, foldchange, pp);
                p->addData(buf);
            }
            goodPeaks->addPeak(p);
            numGood++;
        }
    }

    if (totalChecked == 0) {
        fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
    } else {
        double ratio = ((double) numGood) / ((double) totalChecked);
        fprintf(stderr, "\tDifferential Peaks: %d of %d (%.2lf%% passed)\n", numGood, totalChecked, ratio * 100.0);
    }
    goodPeaks->sortChr();
    return goodPeaks;
}

PeakLibrary *PeakLibrary::stitchRegions(int maxDistance, int mode) const {
    auto *regions = new PeakLibrary();
    //regions->tagsInPeaks = 0.0;
    for (auto& it: chrs)
        it.second->stitchRegions(regions, maxDistance, mode);
    regions->sortChr();
    return regions;
}

// class ChrPeaks -----------------------------------------------

ChrTags::~ChrTags() {
    if (callback) {
        callback();
    }
    for (auto& it: _tags_for_maxtbp) {
        if (it.second == NULL || (it.first == 0 && it.second == _zeromaxtbp.data()))
            continue;
        delete[] it.second;
    }
}

ChrTags::ChrTags(std::string newchr, std::function<void()> callback) {
    chr = newchr;
    totalTags = 0.0;
    maxtbp = 0.0;
    _zeromaxtbp = std::vector<Tag>();
    _zeromaxtbp.reserve(500000);
    pairedEndFlag = 0;
    forceSingleReadFlag = false;
    appearentSize = 0;
}

void ChrPeaks::addPeak(Peak *p) {
    peaks.push_back(p);
}

void ChrPeaks::sort() {
    std::sort(peaks.begin(), peaks.end(), [](const Peak* a, const Peak* b) {
        auto cc = chrcmp(a->chr, b->chr);
        if (cc < 0) return true;
        if (cc > 0) return false;
        if (a->tagStart < b->tagStart) return true;
        if (a->tagStart > b->tagStart) return false;
        if (a->tagEnd < b->tagEnd) return true;
        if (a->tagEnd > b->tagEnd) return false;
        if (a->strand < b->strand) return true;
        if (a->strand > b->strand) return false;
        return false;
    });
//    qsort(peaks, peaks.size(), sizeof(Peak *), &cmpPeaks);
}

int chrcmp(const std::string &chr1, const std::string &chr2) {
    static std::unordered_map<std::string, int> *chrIndex;

    if (chr1 == chr2) return 0;
    if (chr1.empty()) return -1;
    if (chr2.empty()) return 1;

    if (chrIndex == NULL) {
        chrIndex = new std::unordered_map<std::string, int>(10000);
        int index = 1;
        std::string tmp;
        for (int i = 0; i <= 100; i++) {
            tmp = "chr"+std::to_string(i);
            (*chrIndex)[tmp] = index++;
            tmp = "chr"+std::to_string(i)+"_random";
            (*chrIndex)[tmp] = index++;
            tmp = "chr"+std::to_string(i)+"L";
            (*chrIndex)[tmp] = index++;
            tmp = "chr"+std::to_string(i)+"L_random";
            (*chrIndex)[tmp] = index++;
            tmp = "chr"+std::to_string(i)+"R";
            (*chrIndex)[tmp] = index++;
            tmp = "chr"+std::to_string(i)+"R_random";
            (*chrIndex)[tmp] = index++;
        }
        (*chrIndex)["chrI"] = index++;
        (*chrIndex)["chrI"] = index++;
        (*chrIndex)["chrII"] = index++;
        (*chrIndex)["chrIII"] = index++;
        (*chrIndex)["chrIV"] = index++;
        (*chrIndex)["chrV"] = index++;
        (*chrIndex)["chrVI"] = index++;
        (*chrIndex)["chrVII"] = index++;
        (*chrIndex)["chrVIII"] = index++;
        (*chrIndex)["chrIX"] = index++;
        (*chrIndex)["chrU"] = index++;
        (*chrIndex)["chrU_random"] = index++;
        (*chrIndex)["chrX"] = index++;
        (*chrIndex)["chrX_random"] = index++;
        (*chrIndex)["chrXI"] = index++;
        (*chrIndex)["chrXII"] = index++;
        (*chrIndex)["chrXIII"] = index++;
        (*chrIndex)["chrXIV"] = index++;
        (*chrIndex)["chrXV"] = index++;
        (*chrIndex)["chrXVI"] = index++;
        (*chrIndex)["chrXVII"] = index++;
        (*chrIndex)["chrXVIII"] = index++;
        (*chrIndex)["chrXIX"] = index++;
        (*chrIndex)["chrY"] = index++;
        (*chrIndex)["chrY_random"] = index++;
        (*chrIndex)["chrZ"] = index++;
        (*chrIndex)["chrZ_random"] = index++;
        (*chrIndex)["chrM"] = index++;
        (*chrIndex)["chrM_random"] = index++;
        (*chrIndex)["chrMT"] = index++;
        (*chrIndex)["chrMT_random"] = index++;
        (*chrIndex)["chrUn"] = index++;
        (*chrIndex)["chrUn_random"] = index++;
        (*chrIndex)["genome"] = index++;
        (*chrIndex)["null"] = index++;
    }
    auto it1 = chrIndex->find(chr1);
    auto it2 = chrIndex->find(chr2);
    auto end = chrIndex->end();

    if (it1 != end && it2 != end) {
        if (it1->second < it2->second) return -1;
        if (it1->second > it2->second) return 1;
        return 0;
    } else if (it1 != end) {
        return -1;
    } else if (it2 != end) {
        return 1;
    }
    return strcmp(chr1.c_str(), chr2.c_str());
}

void ChrPeaks::countPeakTagsLowMemory(std::unordered_map<std::string, double> &results, ChrTags *ct, char direction, int mode,
                                      ChrPeaks *excludePeaks) {
    ct->loadTags();

    int peakIndex = 0;
    //establish when we can start forgetting about peaks
    // peaks are already sorted by their starting positions;
    auto numPeaks = peaks.size();

    int *finishedLength = new int[numPeaks];
    int *excludeFinishedLength = NULL;

    int excludePeakIndex = 0;
    if (excludePeaks != NULL) {
        excludeFinishedLength = new int[excludePeaks->peaks.size()];
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
        for (int i = 0; i < excludePeaks->peaks.size(); i++) {
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

    auto tags = ct->tags();
    auto size = ct->size();
    for (auto i = tags; i < tags + size; ++i) {
        int p = i->p;
        int d = i->d;
        float v = i->v;

        //if excludePeaks are used, don't bother counting read if it overlaps with the exclude peaks
        if (excludePeaks != NULL && excludePeakIndex < excludePeaks->peaks.size()) {
            if (p < excludePeaks->peaks[excludePeakIndex]->tagStart) {
            } else {
                int stop = 0;
                while (p > excludeFinishedLength[excludePeakIndex]) {
                    excludePeakIndex++;
                    if (excludePeakIndex >= excludePeaks->peaks.size()) {
                        stop = 1;
                        break;
                    }
                }
                if (stop == 0) {
                    int bad = 0;
                    for (int j = excludePeakIndex; j < excludePeaks->peaks.size(); j++) {
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
        results[peaks[i]->name] = total;
    }
    delete[]finishedLength;
    delete[]peakTotals;
    delete[]peakPositions;
    delete[]excludeFinishedLength;
}

void ChrPeaks::stitchRegions(PeakLibrary *regions, int maxDistance, int mode) {
    auto numPeaks = peaks.size();

    auto *totals = new double[numPeaks];
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

                        auto pmap = (double) peaks[j]->uniqMap;
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
    delete[]uniqmap;
}

// class Peaks -----------------------------------------------

Peak::Peak() {
    refPos = 0;
    start = 0;
    end = 0;
    tagStart = 0;
    tagEnd = 0;
    priority = 0;
    strand = STRAND_POSITIVE;
    uniqMap = 0;
    v = 0.0;
    focusRatio = 0.0;
    data = NULL;
}

Peak::Peak(const std::string &newname, const std::string &originalName, const std::string &newchr, int newstart, int newend, int newRef, char dir,
           float value, float ratio, char *otherdata, int mappability, unsigned int newpriority) {
    if (newname.empty()) {
        int L = newchr.size() + 13 + 1 + 6;
        name = newchr + "-" + std::to_string(newstart) + "-" + dir;
    } else {
        name = newname;
    }

    if (!originalName.empty()) {
        ogname = originalName;
    }
    chr = newchr;
    start = newstart;
    end = newend;
    tagStart = newstart;
    tagEnd = newend;
    refPos = newRef;
    strand = dir;
    v = value;
    priority = newpriority;
    focusRatio = ratio;
    uniqMap = mappability;
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
    if (data != NULL) delete[]data;
}


void Peak::print(FILE *fp) {
    char dir = '+';
    if (strand == 1) dir = '-';
    if (v < 1.0) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.3f\t%.3f", name.c_str(), chr.c_str(), start, end, dir, v, focusRatio);
    } else if (v < 10.0) {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.2f\t%.3f", name.c_str(), chr.c_str(), start, end, dir, v, focusRatio);
    } else {
        fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.1f\t%.3f", name.c_str(), chr.c_str(), start, end, dir, v, focusRatio);
    }
    //fprintf(fp, "\t%d", uniqMap);
    if (data != NULL) {
        fprintf(fp, "\t%s", data);
    }
    fprintf(fp, "\n");
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

//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################

TagLibrary::~TagLibrary() {
    for (auto &it: chrs) {
        delete it.second;
        it.second = NULL;
    }
}

TagLibrary::TagLibrary() {
    totalTags = 0;
    sspeFlag = 0;
    totalPositions = 0;
    averageTagsPerPosition = 0.0;
    fragmentLengthEstimate = 0;
    minmapq = 10.0;
    tbp = 0.0;
    maxtbp = 0.0;
    pairedEndFlag = 0;
    averageTagLength = 0.0;
    gsizeEstimate = 0;
}

void TagLibrary::setSingleRead(int flag) {
    for (auto& t: chrs) {
        t.second->forceSingleReadFlag = flag;
    }
}

std::string random_string( size_t length )
{
    auto randchar = []() -> char
    {
        const char charset[] =
                "0123456789"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}

char *unzipFileIfNeeded(const std::string &file, int &zipFlag, int &curFormat) {
    zipFlag = 0;
    if (file.empty()) {
        return NULL;
    }
    int len = file.size();
    std::string newName = file;
    char *command = new char[100 + 2 * len + 1];
    if (len > 4) {
        if (file[len - 3] == '.' && file[len - 2] == 'g' && file[len - 1] == 'z') {
            //gzipped file
            fprintf(stderr, "\tTreating %s as a GNU zip file\n", file.c_str());
            newName.resize(len - 3);
            newName += random_string(10);
            zipFlag = ZIPPED_FLAG_GZ;
            sprintf(command, "gunzip -c \"%s\" > \"%s\"", file.c_str(), newName.c_str());
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'b' && file[len - 2] == 'z' && file[len - 1] == '2') {
            //bz2 zipped file
            fprintf(stderr, "\tTreating %s as a bz2 zip file\n", file.c_str());
            newName.resize(len - 4);
            newName += random_string(10);
            zipFlag = ZIPPED_FLAG_BZ2;
            sprintf(command, "bunzip2 -c \"%s\" > \"%s\"", file.c_str(), newName.c_str());
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'z' && file[len - 2] == 'i' && file[len - 1] == 'p') {
            //zip file
            fprintf(stderr, "\tTreating %s as a zip file\n", file.c_str());
            newName.resize(len - 4);
            newName += random_string(10);
            zipFlag = ZIPPED_FLAG_ZIP;
            sprintf(command, "unzip -p \"%s\" > \"%s\"", file.c_str(), newName.c_str());
            (void) system(command);
        } else if (file[len - 4] == '.' && file[len - 3] == 'b' && file[len - 2] == 'a' && file[len - 1] == 'm') {
            //zip file
            fprintf(stderr, "\tTreating %s as a bam file\n", file.c_str());
            newName.resize(len - 4);
            newName += random_string(10);
            zipFlag = ZIPPED_FLAG_BAM;
            curFormat = FORMAT_SAM;
            sprintf(command, "samtools view -h \"%s\" > \"%s\"", file.c_str(), newName.c_str());
            (void) system(command);
        } else {
            //not recognized, or not zipped...
        }
    }
    delete[]command;
    char *newNam_c_str = new char[newName.size() + 10];
    strcpy(newNam_c_str, newName.c_str());
    return newNam_c_str;
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

void TagLibrary::parseAlignmentFiles(const std::vector<std::string>& files, const int& format, const int& mode) {
//    makeDirectory();
    srand(time(NULL));
    for (auto& file: files) {
        fprintf(stderr, "\n");
        int zippedFlag = 0;
        int currentFormat = format;
        char *workingFilename = unzipFileIfNeeded(file.c_str(), zippedFlag, currentFormat);
        fprintf(stderr, "\tReading alignment file %s\n", workingFilename);
        readAlignment(workingFilename, currentFormat, mode);
        rezipFileIfNeeded(workingFilename, zippedFlag);
    }

    if (totalTags < 1.0) {
        fprintf(stderr, "\t!!! Something is wrong - no reads were added to tag directory !!!\n");
        fprintf(stderr, "\t!!! Check your input files or the makeTagDirectory command options... !!!\n");
        exit(0);
    }
//    printTagInfo();
}

void TagLibrary::readAlignment(char *file, int format, int mode) {

    FILE *fp = fopen(file, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s\n", file);
        return;
    }

    char *buf = new char[BUFFER];

    char cigarCodes[1000];
    int cigarLens[1000];
    int numCodes;

    float numTags = 0;

    double avgReadsPerPos = 0;
    int totalReads = 0;

    while (fgets(buf, BUFFER, fp) != NULL) {
        //sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
        auto line = split(buf, '\t');
        if (line.size() < 3) continue;

        if (format == FORMAT_UNKNOWN) {
            if (line.size() > 2 && format == FORMAT_UNKNOWN) {
                if (line[0][0] == 'c' && line[0][1] == 'h' && line[0][2] == 'r') {
                    fprintf(stderr, "\tGuessing that your alignment file is BED format\n");
                    format = FORMAT_BED;
                }
            }
            if (line.size() > 20 && format == FORMAT_UNKNOWN) {
                fprintf(stderr, "\tGuessing that your alignment file is SAM format (lots of columns though...)\n");
                format = FORMAT_SAM;
            }
            if (line[0][0] == '@' && format == FORMAT_UNKNOWN) {
                if (strcmp(line[0].c_str(), "@SQ") == 0 || strcmp(line[0].c_str(), "@HD") == 0 || strcmp(line[0].c_str(), "@RG") == 0) {
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
            int samFlag = std::stoi(line[1]);
            if (samFlag & 0x4) {
                //unmapped flag
                continue;
            }

            auto& name = line[0];
            auto& chr = line[2];
            int initLen = 0;
            auto len = line[9].size();

            bool peFlag = false;
            if (samFlag & 0x1) {
                peFlag = true;
            }

//            if (peFlag && line.size() > 8) {
//                len = std::abs(std::stoi(line[8]));
//            }

            if (len < minReadLength || len > maxReadLength)
                continue;
            auto pos = std::stoi(line[3]);
            if (pos == 0) continue;
            if (line[5][0] == '*') continue;

            auto dir = STRAND_POSITIVE;
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

            if (mode == MODE_UNIQUE) {
                if (samFlag & 0x100) continue;
                if (minmapq < 0) {
                    int good = 1;
                    double alnScore = -1e8;
                    double secScore = -1e9;
                    for (size_t i = 11; i < line.size(); i++) {
                        if (line[i].size() < 6) continue;
                        if (strncmp(line[i].c_str(), "X0:i:", 5) == 0) {
                            if (strcmp(&(line[i][5]), "1") != 0) {
                                good = 0;
                            }
                        } else if (strncmp(line[i].c_str(), "AS:i:", 5) == 0) {
                            sscanf(&(line[i][5]), "%lf", &alnScore);
                        } else if (strncmp(line[i].c_str(), "XS:i:", 5) == 0) {
                            sscanf(&(line[i][5]), "%lf", &secScore);
                        }
                    }
                    if (alnScore > secScore) {
                    } else {
                        good = 0;
                    }
                    if (good == 0) continue;
                } else {
                    double mapqScore = std::stod(line[4]);
                    if (mapqScore < minmapq) continue;
                }
            } else if (mode == MODE_KEEPONE) {
                if (samFlag & 0x100) continue;
            }

            int startPos = getRightCoordFromCIGAR(line[5], dir, cigarCodes, cigarLens, numCodes, initLen);
            if (startPos == CIGAR_ERROR) {
                fprintf(stderr, "\tError in read: %s\n", name.c_str());
                startPos = 0;
            }

            pos += startPos;
            len = initLen;
//            if (!peFlag) {
//                len = initLen;
//            }
            float v = 1.0;
            if (peFlag) v = 0.5;
            addAlignedTag(name, chr, pos, (char) dir, len, v);
        }
        else if (format == FORMAT_BED) {
            auto& chr = line[0];
            std::string name;

            auto dir = STRAND_POSITIVE;
            numTags = 1.0;
            if (line.size() > 3) {
                if (line[3][0] == '+' || line[3][0] == '-') {
                    if (line[3][0] == '-') dir = STRAND_NEGATIVE;
                } else if (line.size() > 5) {
                    name = line[3];
                    if (line[5][0] == '-') {
                        dir = STRAND_NEGATIVE;
                    }
                }
            }
            if (line.size() > 4) {
                numTags = std::stof(line[4]);
            }
            int start = std::stoi(line[1]);
            int end = std::stoi(line[2]);
            int pos;
            if (dir == 0) {
                pos = start + 1; //0 base index
            } else {
                pos = end;
            }
            size_t len = abs(end - start);
            if (len < minReadLength || len > maxReadLength) continue;
            addAlignedTag(name, chr, pos, (char) dir, len, numTags);
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

    delete[]buf;
}

int
TagLibrary::getRightCoordFromCIGAR(std::string& str, int dir, char *cigarCodes, int *cigarLens, int &numCodes, int &initLen) {
    numCodes = 0;
    int i = 0;
    const char *start = str.c_str();
    int lastStartI = 0;
    int totalLen = 0;

    int firstStart = -1;
    int firstEnd = -1;
    int lastStart = -1;
    int lastEnd = -1;
    int firstActive = 0;
    int lastActive = 0;

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


void TagLibrary::optimizeTags() {

    fprintf(stderr, "\n\tOptimizing tag files...\n");
    totalTags = 0;
    totalPositions = 0;

    int tagAdjust = fragmentLengthEstimate / 2;

    gsizeEstimate = 0;
    for (auto& it: chrs) {
        ChrTags *ct = it.second;

        tagAdjust = fragmentLengthEstimate / 2;

        ct->optimizeTags(tagAdjust);
        totalTags += ct->totalTags;
        totalPositions += ct->size();
        gsizeEstimate += ct->appearentSize;
    }
    fprintf(stderr, "\tEstimated genome size = %lld\n", gsizeEstimate);
    fprintf(stderr, "\tTotal Tags = %.1f\n", totalTags);
    fprintf(stderr, "\tTotal Positions = %lld\n", totalPositions);

    setMaxTBP(maxtbp);

    if (gsizeEstimate > 1) {
        tbp = totalTags / gsizeEstimate;
        fprintf(stderr, "\tEstimated average read density = %.6lf per bp\n", tbp);
    }
}

void
TagLibrary::addAlignedTag(const std::string &name, const std::string &chr, int pos, char dir, int length, float value) {

    if (pos < 1) pos = 1;
    auto it = chrs.find(chr);
    ChrTags *chrtags;
    if (it == chrs.end()) {
        chrs[chr] = new ChrTags(chr);
        chrtags = chrs[chr];
        chrtags->pairedEndFlag = pairedEndFlag;
    }
    else {
        chrtags = it->second;
    }
    chrtags->addAlignedTag(name, pos, dir, length, value);
    totalTags += value;
}


void TagLibrary::setMaxTBP(float max) {
    maxtbp = max;
    for (auto& it: chrs) {
        it.second->setMaxTBP(max);
    }
}

PeakLibrary *TagLibrary::findPutativePeaks(int peakSize, int minDist, char strand, float minCount) {
    PeakLibrary *putativePeaks = new PeakLibrary(10000000);
    for (auto& it: chrs) {
        it.second->loadTags();
        it.second->findPutativePeaks(putativePeaks, peakSize, minDist, strand, minCount);
        it.second->freeTags();
    }
    putativePeaks->sortChr();
    return putativePeaks;
}


void ChrTags::getTagLengthDistribution(double* dist,int max) {
    auto tags = this->tags();
    auto size = this->size();
    for (auto it = tags; it < tags + size; ++it) {
        int len = it->len;
        if (len < 0) len = 0;
        if (len >= max) {
            len = max-1;
        }
        dist[len]+=it->v;
    }
}


double* TagLibrary::getTagLengthDistribution(FILE* nfp, int &max) {
    if (max == 0) {
        max = 32000;
    }

    double* dist = new double[max];
    for (int i=0;i<max;i++) {
        dist[i] = 0;
    }

    for (auto& it: chrs)
        it.second->getTagLengthDistribution(dist, max);

    float totalDist = 0;
    averageTagLength = 0.0;
    for (int i=0;i<max;i++) {
        dist[i] /= (double)totalTags;
        averageTagLength += ((double)i)*dist[i];
        totalDist += dist[i];
    }

    fprintf(stderr, "\tAverage tag length = %.1lf\n", averageTagLength);
    return dist;
}




double *TagLibrary::getTagCountDistribution(FILE *nfp, int &max) {
    if (max == 0) {
        max = MAX_TAGS_PER_BP;
    }
    double *dist = new double[max];
    for (int i = 0; i < max; i++) {
        dist[i] = 0;
    }

    int scaleFactor = 1;

    for (auto& it: chrs) {
        it.second->getTagCountDistribution(dist, max, scaleFactor);
    }

    float totalDist = 0;
    averageTagsPerPosition = 0;
    for (int i = 0; i < max; i++) {
        dist[i] /= (double) totalPositions;
        averageTagsPerPosition += ((double) i) / ((double) scaleFactor) * dist[i];
        totalDist += dist[i];
    }

    fprintf(stderr, "\tAverage tags per position = %.3lf\n", averageTagsPerPosition);
    return dist;
}

double TagLibrary::getAdjustedTagTotal() {
    double total = 0;
    for (auto& it: chrs) {
        total += it.second->totalTags;
    }
    return total;
}

void
ChrTags::autoCorrelateTags(double *sameStrand, double *diffStrand, int windowSize, double maxTags, double &totalCount,
                           double *sameStrandN, double *diffStrandN) {
    auto tags = this->tags();
    auto size = this->size();
    int halfWindow = (windowSize)/2;

    for (int i=0;i<size;i++) {
        double refV = (double)tags[i].v;
        for (int j=i+1;j<size;j++) {
            int diff = tags[j].p - tags[i].p;
            if (diff > halfWindow) break;
            double v = 0.0; // (double)(tags[i].v*tags[i].v);
            v = 1.0;
            if (tags[i].d == tags[j].d) {
                if (tags[i].d == 0) {
                    sameStrand[halfWindow+diff] += v;
                    sameStrandN[halfWindow+diff] += 1.0;
                } else {
                    sameStrand[halfWindow-diff] += v;
                    sameStrandN[halfWindow-diff] += 1.0;
                }
            } else {
                if (tags[i].d == 0) {
                    diffStrand[halfWindow+diff] += v;
                    diffStrandN[halfWindow+diff] += 1.0;
                } else {
                    diffStrand[halfWindow-diff] += v;
                    diffStrandN[halfWindow-diff] += 1.0;
                }
            }
            totalCount += v;
        }
        if (totalCount > maxTags) break;
    }
}

void TagLibrary::autoCorrelateTags(int windowSize, double maxTags) {
    const int avgWindow = 7;
    double totalCount = 0;
    auto* sameStrand = new double[windowSize+1];
    auto* diffStrand = new double[windowSize+1];
    auto* sameStrandN = new double[windowSize+1];
    auto* diffStrandN = new double[windowSize+1];
    auto* smoothed = new double[windowSize+1];

    for (int i=0;i<windowSize+1;i++) {
        sameStrand[i] = 0.0;
        diffStrand[i] = 0.0;
        sameStrandN[i] = 0.0;
        diffStrandN[i] = 0.0;
        smoothed[i] = 0.0;
    }

    for (auto& it: chrs) {
        it.second->autoCorrelateTags(sameStrand, diffStrand,windowSize,maxTags,totalCount,
                                      sameStrandN, diffStrandN);
    }

    int halfWindow = windowSize/2;
    double backEst = 0;
    double maximum = 0;
    int lastEstimate = 0;

    int minFragLen = 40;
    if (averageTagLength > 0.0)
        minFragLen = (int) averageTagLength+3;
    int bottomFlag = 0;

    fragmentLengthEstimate = 0;
    peakSizeEstimate = 0;

    const int autocorr_backoffset = -2;
    const int autocorr_halfwindow = 7;

    for (int i=0;i<windowSize;i++) {
        int offset = i-halfWindow;
        double avg = 0;
        double n= 0;
        for (int j=-avgWindow;j<avgWindow;j++) {
            if (j+i < 0) continue;
            if (j+i >= windowSize) break;
            avg += diffStrand[j+i];
            n+=1.0;
        }
        if (n>0.5) avg /= n;
        smoothed[i] = avg;

        if (offset == autocorr_backoffset-autocorr_halfwindow) {
            backEst = smoothed[i];
        }
        if (offset > minFragLen) {
            if (smoothed[i] > maximum && bottomFlag) {
                fragmentLengthEstimate = offset;
                maximum = smoothed[i];
            } else if (bottomFlag == 0) {
                if (smoothed[i] > smoothed[i-1]) bottomFlag =1;
            }
            if (smoothed[i] < backEst) {
                if (lastEstimate != fragmentLengthEstimate) {
                    peakSizeEstimate = offset-fragmentLengthEstimate;
                    lastEstimate = fragmentLengthEstimate;
                }
            }
        }
    }

    fprintf(stderr, "\tFragment Length Estimate: %d\n", fragmentLengthEstimate);

    if (fragmentLengthEstimate <= minFragLen) {
        fprintf(stderr, "\t\t!!! No reliable estimate for fragment length\n");
        fprintf(stderr, "\t\t!!! PLEASE SET MANUALLY (currently set to 150, edit tagInfo.txt to change)\n");
        fragmentLengthEstimate = 150;
    }
    fprintf(stderr, "\tPeak Width Estimate: %d\n", peakSizeEstimate);
    if (peakSizeEstimate < fragmentLengthEstimate) {
        peakSizeEstimate = fragmentLengthEstimate;
        fprintf(stderr, "\t\t!!! No reliable estimate for peak size\n");
        fprintf(stderr, "\t\tSetting Peak width estimate to be equal to fragment length estimate\n");
    }


    double posBackSignal = 0.0;
    double negBackSignal = 0.0;
    double posBackN=0.0;
    double negBackN=0.0;
    double posSignal = 0.0;
    double negSignal = 0.0;
    double posN=0.0;
    double negN=0.0;

    for (int i=0;i<windowSize;i++) {
        int offset = i-halfWindow;
        if (offset < fragmentLengthEstimate*-1 || offset > fragmentLengthEstimate*2) {
            posBackSignal += sameStrand[i];
            posBackN += 1.0;
            negBackSignal += diffStrand[i];
            negBackN += 1.0;
        } else {
            posSignal += sameStrand[i];
            posN += 1.0;
            negSignal += diffStrand[i];
            negN += 1.0;
        }
    }

    if (posN > 0.0) posSignal /= posN;
    if (negN > 0.0) negSignal /= negN;
    if (posBackN > 0.0) posBackSignal /= posBackN;
    if (negBackN > 0.0) negBackSignal /= negBackN;
    double posEnrichment = 1.0;
    double negEnrichment = 1.0;
    double posVsNegEnrichment = 1.0;
    if (posBackSignal > 0.0) posEnrichment = posSignal/posBackSignal;
    if (negBackSignal > 0.0) negEnrichment = negSignal/negBackSignal;
    if (negSignal > 0.0) posVsNegEnrichment = posSignal/negSignal;

    fprintf(stderr, "\tAutocorrelation quality control metrics:\n");
    fprintf(stderr, "\t\tSame strand fold enrichment: %.1lf\n", posEnrichment);
    fprintf(stderr, "\t\tDiff strand fold enrichment: %.1lf\n", negEnrichment);
    fprintf(stderr, "\t\tSame / Diff fold enrichment: %.1lf\n", posVsNegEnrichment);
    fprintf(stderr, "\n");

    delete []sameStrand;
    delete []diffStrand;
    delete []sameStrandN;
    delete []diffStrandN;
    delete []smoothed;

}

// class ChrTags ------------------------------------------------------------


ChrTags::ChrTags(std::string newchr) {
    chr = newchr;
    totalTags = 0.0;
    maxtbp = 0.0;
    _zeromaxtbp = std::vector<Tag>();
    _zeromaxtbp.reserve(500000);
    pairedEndFlag = 0;
    forceSingleReadFlag = false;
    appearentSize = 0;
}

void ChrTags::addAlignedTag(const std::string &name, int pos, char dir, int length, float value) {
//    if (_tags_for_maxtbp.size() == 1)
    auto tag = Tag();
    tag.p = pos;
    tag.len = length;
    tag.v = value;
    tag.d = dir;
    _zeromaxtbp.push_back(tag);
//    _zeromaxtbp.emplace_back(pos, length, value, dir);
//    else
//        tags().emplace_back(pos, length, value, dir);
    if (pos > appearentSize) appearentSize = (long long int) pos;
}

void ChrTags::loadTags() {
//    if (loaded && adjusted) return;
//    if (pairedEndFlag && forceSingleReadFlag == 1) {
//        optimizeOverride = 1;
//        fprintf(stderr, "Reading single...\n");
//    }
//    adjustTags();
//    if (pairedEndFlag && forceSingleReadFlag == 1) optimizeOverride = 0;
}

void ChrTags::clipTags() {
    if (int(maxtbp) == 0)
        return;

    auto tags = this->tags();
    auto size = this->size();
    totalTags = 0.0;
    for (auto tag = tags; tag < tags + size; ++tag) {
        if (tag->v > maxtbp) tag->v = maxtbp;
        totalTags += tag->v;
    }
    _total_tags_for_maxtbp[maxtbp] = totalTags;
}

void ChrTags::optimizeTags(int tagAdjust) {
//    sortMergeTags();
    if (int(maxtbp) != 0)
        throw std::runtime_error("It shouldn't have happened.");
    sortMergeTags();
    adjustTags(tagAdjust);
    sortMergeTags();
    clipTags();
}

void ChrTags::sortMergeTags() {
    Tag *tags;
    size_t size;

    auto it = _tags_for_maxtbp.find(0);
    if (it == _tags_for_maxtbp.end()){
        tags = _zeromaxtbp.data();
        size = _zeromaxtbp.size();
    }
    else {
        tags = it->second;
        size = _tags_size_for_maxtbp[0];
    }
    //fprintf(stderr, "%s  - optimizedFlag = %d\n", chr, optimizedFlag);

    std::sort(tags, tags + size, &cmpTags);

    // merge tags if needed
    size_t last = 0;
    totalTags = 0;
    if (size > 0) {
        totalTags = tags[0].v;
    } else {
        fprintf(stderr, "No tags in tag file!!!\n");
        return;
    }

    for (size_t i = 1; i < size; i++) {
        totalTags += tags[i].v;
        auto j = i - 1;
//        if (tags[i].p == tags[j].p && tags[i].d == tags[j].d) {
        if (tags[i].p == tags[j].p && tags[i].d == tags[j].d && tags[i].len == tags[j].len) {
            tags[last].v += tags[i].v;
        } else {
            last++;
            if (last == i) continue;
            tags[last].copy(&tags[i]);
        }
    }
    if (it == _tags_for_maxtbp.end()) {
        _zeromaxtbp.resize(last + 1);
        _tags_size_for_maxtbp[0] = _zeromaxtbp.size();
        _tags_for_maxtbp[0] = _zeromaxtbp.data();
    }
    else {
        _tags_size_for_maxtbp[0] = last + 1;
    }

    _total_tags_for_maxtbp[0] = totalTags;
}

bool cmpTags(const Tag& a, const Tag& b) {
    if (a.p < b.p) return true;
    if (a.p > b.p) return false;
    if (a.d < b.d) return true;
    if (a.d > b.d) return false;
    if (a.len < b.len) return true;
    if (a.len > b.len) return false;
    return false;
}


void ChrTags::setMaxTBP(float max) {
    maxtbp = max;
    auto it = _tags_for_maxtbp.find(maxtbp);
    if (it == _tags_for_maxtbp.end()) {
        auto size = _tags_size_for_maxtbp[0];
        auto tags = _tags_for_maxtbp[0];

        Tag *arr = new Tag[size];
        std::copy(tags, tags + size, arr);
        _tags_for_maxtbp[maxtbp] = arr;
        _tags_size_for_maxtbp[maxtbp] = size;
        clipTags();
    }
    totalTags = _total_tags_for_maxtbp[maxtbp];
}

void ChrTags::adjustTags(int tagAdjust) {
    auto tags = this->tags();
    auto size = this->size();
    if (size == 0) return;

    if (pairedEndFlag && !forceSingleReadFlag)
        return;

    int minPosition = tags[0].p;
    int maxPosition = tags[size - 1].p;

//    if (pairedEndFlag) {
//        for (auto & tag : tags) {
////            tag.p += tag.len / 2;
//            if (tag.d == 0) {
//                tag.p += tag.len / 2;
//            } else {
//                tag.p -= tag.len / 2;
//            }
//            if (tag.p > maxPosition) {
//                tag.p = maxPosition;
//            }
//            if (tag.p < minPosition) {
//                tag.p = minPosition;
//            }
//        }
//    }
//    else
    if (tagAdjust != 0) {
        for (auto tag = tags; tag < tags + size; ++tag) {
            if (tag->d == 0) {
                tag->p += tagAdjust;
            } else {
                tag->p -= tagAdjust;
            }
            if (tag->p > maxPosition) {
                tag->p = maxPosition;
            }
            if (tag->p < minPosition) {
                tag->p = minPosition;
            }
        }
    }
}


// negative value for minDist is used to indicate region based peak finding
void ChrTags::findPutativePeaks(PeakLibrary *putativePeaks, int peakSize, int minDist, char strand, double minCount) {
    auto tags = this->tags();
    auto size = this->size();
    doubleIndex *tagCount = new doubleIndex[size];

    int regionFlag = 0;
    if (minDist < 0) {
        minDist *= -1;
        regionFlag = 1;
    }

    int lastPosIndex = -1;
    int lastPosFinish = -1;
    int lastPosPosition = -1;
    int maxEND = 0;

    if (size > 0) {
        maxEND = tags[size - 1].p;
    }

    fprintf(stderr, "\t\tFinding peaks on %s (minCount=%.1lf), total tags positions = %zu\n",
            chr.c_str(), minCount, size);

    for (int i = 0; i < size; i++) {
        float v = tags[i].v;
        char d = tags[i].d;
        tagCount[i].index = i;
        int p = tags[i].p;

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
        for (int j = start; j < size; j++) {
            int diff = tags[j].p - p;
            if (diff > peakSize) break;
            tagCount[i].v += tags[j].v;
            tagCount[i].vp += tags[j].v * (double) ((tags[j].p) - p);
            lastPosFinish = i;
        }
        lastPosPosition = p;
    }

    char *mask = new char[size];
    int *centers = new int[size];
    int pid = 1;

    for (int i = 0; i < size; i++) {
        mask[i] = 0;
        if (tagCount[i].v < 0.000001) tagCount[i].v = 1.0;
        tagCount[i].position = tags[i].p + (int) floor(tagCount[i].vp / tagCount[i].v);
        centers[i] = tagCount[i].position;
    }

//    std::sort(tagCount, tagCount + tags.size(), [](const doubleIndex& a, const doubleIndex& b) {
//        return a.v > b.v;
//    });
    qsort(tagCount, size, sizeof(doubleIndex), &cmpDoubleIndex);

    int d = minDist;
    if (d < 1) d = 1;
    int expectedNumPeaks = appearentSize / d * 2;
    if (expectedNumPeaks < 10) expectedNumPeaks = 10;

    int halfSize = (peakSize) / 2;

    for (int i = 0; i < size; i++) {
        const unsigned int& index = tagCount[i].index;
        const double& value = tagCount[i].v;
        const int& center = tagCount[i].position;

        if (value < minCount) break;

        const char d = strand == STRAND_BOTH ? 0 : tags[index].d;
        if (!mask[index]) {
            int start = center - halfSize;
            int end = start + peakSize;
            if (start < 1) {
                start = 1;
            }
            if (end > maxEND) {
                end = maxEND;
            }

//            sprintf(pname, "%s-%d", chr.c_str(), pid++);
            //sprintf(pname, "%d",rand());
            //fprintf(stderr, "%s\t%s\t%d\t%d\t%d\t%d\t%f\n",pname,chr,start,end,center,d,value);
            putativePeaks->addPeak(chr+"-"+std::to_string(pid), chr, start, end, center, d, value, 0, NULL, -1, 0);
            pid += 1;
        }
        if (regionFlag == 0 || !mask[index]) {
            for (int j = index; j < size; j++) {
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

void ChrTags::getTagCountDistribution(double *dist, int max, int scaleFactor) {
    //forceSingleReadFlag = 1;
    loadTags();
    auto tags = this->tags();
    auto size = this->size();
    for (size_t i = 0; i < size; i++) {
        int v = 0;
        v = (int) floor((double) tags[i].v * scaleFactor);
        if (v >= max) {
            v = max - 1;
        }
        dist[v]++;
    }
    freeTags();
    //forceSingleReadFlag = 0;
}

void ChrTags::freeTags() {

}

// class Tag ---------------------------------------------------------------

//Tag::Tag() {
//    p = 0;
//    len = -1;
//    v = 0.0;
//    d = 0;
//}

void Tag::copy(Tag *src) {
    p = src->p;
    len = src->len;
    v = src->v;
    d = src->d;
}

