#pragma once

#include <cstdio>

#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <ctime>
#include <climits>
#include <algorithm>
#include <pthread.h>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <functional>
#include <sstream>

#include <list>
#include <cstring>

#include "statistics.h"

#define MAX_TAGS_PER_BP 32000
#define PEAK_LIBRARY_DEFAULT_SIZE 1000000
#define EMPTY_DOUBLE_CHECK -1.20e100
#define EMPTY_INT -1080706050
#define EMPTY_DOUBLE -1.23e100

#define FORMAT_UNKNOWN 0
#define FORMAT_BED 3
#define FORMAT_SAM 7

#define MODE_UNIQUE 0
#define MODE_KEEPONE 1

#define NULL_REF -123456789
#define NULL_OFFSET -123456789
#define POSITIVE_STRAND 0
#define NEGATIVE_STRAND 1
#define BOTH_STRANDS 2

#define CIGAR_ERROR -12345678

#define COUNT_MODE_TOTAL 0
#define COUNT_MODE_TBP 1
#define COUNT_MODE_RATIO 2

#define FINDPEAKS_MINDIST_DEFAULT 2.0

#define STRAND_POSITIVE 0
#define STRAND_NEGATIVE 1
#define STRAND_BOTH 2
#define STRAND_SEPARATE 3

#define DEFAULT_NORM_TOTAL 1e7

#define PEAKFINDER_FDRSIZE 100000
#define PEAKFINDER_FILTER_MODE_FDR 1
#define PEAKFINDER_FILTER_MODE_POISSON 2
#define PEAKFINDER_FILTER_MODE_THRESH 3

#define DIFFPEAK_MODE_DIFF 1
#define DIFFPEAK_MODE_SAME 2
#define DIFFPEAK_MODE_REV 3

#define PEAKFRACTION_REGION_MINDIST 4.0
#define REGION_MODE_HISTONE 0
#define REGION_MODE_GROSEQ 1
#define REGION_GROSEQ_FOLDCHANGE 3.0

#define PEAK_STYLE_HISTONE 1

#define DEFAULT_GSIZE 2000000000

#define ZIPPED_FLAG_GZ 1
#define ZIPPED_FLAG_BZ2 2
#define ZIPPED_FLAG_ZIP 3
#define ZIPPED_FLAG_BAM 4

class TagLibrary;

class ChrTags;

struct Tag;

class Peak;

class ChrPeaks;

class PeakLibrary;


struct PeakFinder {
    std::string outputFileName;
    int peakSize;
    int localSize;
    int inputSize;
    float tagThresh;
    float minTagThresh;
    int minDist;
    double tagsUsedForClustering;
    double totalTags;
    double normTotal;
    float maxtbp;
    float tbpThreshold;
    float tbp;
    float tpp;
    int stitchMode;
    float tbpInput;
    char strand;
    long long int gsize;
    double fdr;
    float fdrThresh;
    double poisson;
    double poissonInput;
    double poissonLocal;
    float poissonThresh;
    int filterMode;
    int diffMode;
    int tbpAuto;
    int numPeaks;
    double regionSubDivision;
    TagLibrary *tags;
    TagLibrary *input;
    double *fdrTable;
    int fdrSize;
    double *poissonTable;
    double tagsInPeaks;
    int regionFlag;
    char *extraheader;

    int style;

    double inputFold;
    double localFold;
    double clonalFold;

    PeakFinder();

    ~PeakFinder();

    void addHeader(char *);

    void print(FILE *fp);

    void setOutputFile(std::string name);

    void setTagLibraries(TagLibrary *exp, TagLibrary *input);

    void determineMaxTBP();

    void approxFdrTable();

    void checkParameters();

    PeakLibrary *findPeaks();

    PeakLibrary *filterPeaks(PeakLibrary *);
};

class TagLibrary {
public:
    std::map<std::string, ChrTags *> chrs;

    double totalTags;
    long long int totalPositions;

    double tbp;
    float maxtbp;

    double minmapq;
    size_t minReadLength;
    size_t maxReadLength;

    int fragmentLengthEstimate;
    int peakSizeEstimate;

    double averageTagsPerPosition;
    long long int gsizeEstimate;
    double averageTagLength;

    bool sspeFlag;
    bool pairedEndFlag;

    TagLibrary();

    ~TagLibrary();

    void autoCorrelateTags(int windowSize, double maxTags);

    void parseAlignmentFiles(const std::vector<std::string> &files, const int &format, const int &mode);

    void setMaxTBP(float maxTagsPerBp);

    void setSingleRead(int singleReadFlag);

    double *getTagCountDistribution(FILE *nfp, int &max);

    double *getTagLengthDistribution(FILE *nfp, int &max);

    double getAdjustedTagTotal();

    //internal functions
    void readAlignment(char *file, int format, int mode);

    int
    getRightCoordFromCIGAR(std::string &str, int dir, char *cigarCodes, int *cigarLens, int &numCodes, int &initLen);

    void addAlignedTag(const std::string &name, const std::string &chr, int pos, char dir, int length, float value);

    void optimizeTags();

    PeakLibrary *findPutativePeaks(int peakSize, int minDist, char strand, float minCount);
};


class ChrTags {
public:
    std::map<int, Tag *> _tags_for_maxtbp;
    std::map<int, size_t> _tags_size_for_maxtbp;
    std::map<int, double> _total_tags_for_maxtbp;
    std::vector<Tag> _zeromaxtbp;
    std::function<void()> callback;

    double totalTags;
    float maxtbp;
    long long int appearentSize;
    std::string chr;
    int pairedEndFlag;
    bool forceSingleReadFlag;

    explicit ChrTags(std::string newchr);

    ChrTags(std::string newchr, std::function<void()> callback);

    ~ChrTags();

    Tag *tags() {
        return _tags_for_maxtbp[maxtbp];
    };

    size_t size() {
        return _tags_size_for_maxtbp[maxtbp];
    };

    void addAlignedTag(const std::string &name, int pos, char dir, int length, float value);


    void optimizeTags(int tagAdjust);

    void setMaxTBP(float maxTagsPerBp);

    void findPutativePeaks(PeakLibrary *putativePeaks, int peakSize, int minDist, char strand, double minCount);

    void autoCorrelateTags(double *sameStrand, double *diffStrand, int windowSize,
                           double maxTags, double &totalCount, double *sameStrandN, double *diffStrandN);

    void getTagCountDistribution(double *d, int max, int scaleFactor);

    void getTagLengthDistribution(double *dist, int max);

    void loadTags();

    void clipTags();

    void adjustTags(int tagAdjust);

    void sortMergeTags();

    void freeTags();
};

bool cmpTags(const Tag &, const Tag &);

struct Tag {
    int p; // position
    int len; // position
    float v; // value = number of tags
    char d; // direction 0=+,1=-
    Tag() = default;

    void copy(Tag *src);
};

class PeakLibrary {
public:
    std::unordered_map<std::string, ChrPeaks *> chrs; // holds pointers to ChrPeak objects
    std::unordered_map<std::string, Peak *> peaks; // holds pointers to Peak objects
    int numPeaks;
    double tagsInPeaks;
    int fixedFlag;
    int duplicateWarningFlag;
    std::vector<Peak *> peakOrder;
    std::unordered_map<std::string, int> duplicates;

    PeakLibrary *getDifferentialPeaks(TagLibrary *tags, TagLibrary *input,
                                      double foldThreshold, double poissonThresh, int mode, int start, int end,
                                      char strand, int strFlag);
    PeakLibrary *filterClonalPeaks(TagLibrary *tags, int peakSize,
                                   double threshold, int mode, char strand);

    PeakLibrary *filterLocalPeaks(TagLibrary *tags, int peakSize, int localSize,
                                  double threshold, double poissonThresh, int mode, char strand);
    PeakLibrary *stitchRegions(int maxDistance, int mode) const;

    Peak *
    addPeak(const std::string &name, const std::string &chr, int start, int end, int midpoint, char dir, float value,
            float ratio, char *extraData, int mappability, unsigned int priority);

    PeakLibrary();

    ~PeakLibrary();

    void print(FILE *);

    explicit PeakLibrary(int expectedNumberOfPeaks); // default is 100000

    void initialize(int expectedNumPeaks);

    void setDefaultPeakOrder();

    std::unordered_map<std::string, double>
    countPeakTagsLowMemory(TagLibrary *tags, char direction, int mode, PeakLibrary *excludePeaks);

    void setPeakTagSizeRefPos(int offset, int startOffset, int endOffset);

    void setPeakTagSizeFixed(int startOffset, int endOffset);

    void normalizePeakScore(float normFactor);

    void addPeak(Peak *p);

    void sortChr();

};

class ChrPeaks {
public:
    std::vector<Peak *> peaks;

    void addPeak(Peak *p);

    void sort();

    void stitchRegions(PeakLibrary *regions, int maxDistance, int mode);

    void countPeakTagsLowMemory(std::unordered_map<std::string, double> &results, ChrTags *ct, char direction, int mode,
                                ChrPeaks *excludePeaks);
};

class Peak {
public:
    std::string name;
    std::string chr;
    int refPos;
    int start;
    int end;
    int tagStart;
    int tagEnd;
    unsigned int priority;
    char strand;
    int uniqMap;
    float v;
    float focusRatio;
    char *data;
    std::string ogname;

    Peak();

    Peak(const std::string &newname, const std::string &originalName, const std::string &newchr,
         int newstart, int newend, int newRef, char dir, float value,
         float ratio, char *otherdata, int mappability, unsigned int newpriority);

    ~Peak();

    void print(FILE *fp);

    void setPeakTagSizeRefPos(int newOffset, int startOffset, int endOffset);

    void setPeakTagSizeFixed(int startOffset, int endOffset);

    void setOffset(int offset);

    void addData(char *str);
};


// support classes

class doubleIndex {
public:
    double v;
    double vp;
    unsigned int index;
    int position;
};

int cmpDoubleIndex(const void *a, const void *b);

int cmpPeaks(const void *, const void *);

char *unzipFileIfNeeded(const std::string &file, int &zipFlag, int &curFormat);

void rezipFileIfNeeded(char *file, int zipFlag);

int chrcmp(const std::string& chr1, const std::string& chr2);
