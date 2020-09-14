#pragma once
#include "SeqTag.h"

TagLibrary* load(const std::vector<std::string>& files) {
    auto library = new TagLibrary();
    library->maxReadLength = std::numeric_limits<size_t>::max();
    library->minReadLength = 0;
    library->minmapq = 10.0;
    library->sspeFlag = false;
    library->pairedEndFlag = false;

    int format = FORMAT_UNKNOWN;
    int mode = MODE_UNIQUE;
    library->parseAlignmentFiles(files, format, mode);

    library->setSingleRead(1);
    library->optimizeTags();

    // it is needed to feel the average tags something field.
    auto max = 32000;
    auto d = library->getTagLengthDistribution(NULL, max);
    if (d != NULL) delete []d;

    max = MAX_TAGS_PER_BP;
    d = library->getTagCountDistribution(NULL, max);
    if (d != NULL) delete []d;

    constexpr int autoCorrRange = 4000;
    constexpr double autoCorrMaxTags = 1e9;
    library->autoCorrelateTags(autoCorrRange, autoCorrMaxTags);
    library->optimizeTags();
    return library;
}


void findPeaks(PeakFinder& pf) {
    FILE* fp = stdout;
    fp = fopen(pf.outputFileName.c_str(), "w");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s for writing!!!\n", pf.outputFileName.c_str());
        exit(1);
    }
    auto peaks = pf.findPeaks();
    pf.print(fp);
    peaks->print(fp);
    delete peaks;
    fclose(fp);
}
