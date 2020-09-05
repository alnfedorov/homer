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

#include "Hashtable.h"
#include "statistics.h"

#ifndef SEQTAG_H
#define SEQTAG_H

#define WHITE_SPACE 7

#define MAX_READ_LENGTH 32000
#define MAX_TAGS_PER_BP 32000
#define PEAK_LIBRARY_DEFAULT_SIZE 1000000

#define PEAK_READ_MODE_NORMAL 0
#define PEAK_READ_MODE_ANNOTATION 1
#define PEAK_READ_MODE_COUNTING 2
#define PEAK_READ_MODE_2DBED 3

#define PE_READ_FILTER_KEEPBOTH 0
#define PE_READ_FILTER_KEEPFIRST 1
#define PE_READ_FILTER_KEEPSECOND 2

#define FORMAT_UNKNOWN 0
#define FORMAT_BED 3
#define FORMAT_SAM 7 
#define FORMAT_PEAK 8

#define MODE_UNIQUE 0
#define MODE_KEEPONE 1
#define MODE_KEEPALL 2
#define MODE_BED_FORCE5TH 4

#define NULL_REF -123456789
#define NULL_OFFSET -123456789
#define ALL_PEAK_EXPS -1
#define POSITIVE_STRAND 0 
#define NEGATIVE_STRAND 1 
#define BOTH_STRANDS 2
#define NULL_INT_VALUE -123456789

#define TAGADJUST_DEFAULT 75
#define TAGADJUST_AUTO -123456789
#define FRAGMENT_LEN_AUTO -123456789
#define FRAGMENT_LEN_GIVEN -12345
#define FRAGMENT_LEN_PE -123456

#define CIGAR_ERROR -12345678

#define TOTAL_READS_TAGDIR_DEFAULT -1.0
#define TOTAL_READS_TAGDIR_ALL -2.0

#define COUNT_MODE_TOTAL 0
#define COUNT_MODE_TBP 1
#define COUNT_MODE_RATIO 2

#define FINDPEAKS_MINDIST_DEFAULT 2.0

#define TAGLIBRARY_SINGLE_FLAG 0
#define TAGLIBRARY_PE_READ1_FLAG 1
#define TAGLIBRARY_PE_READ2_FLAG 2

#define FLOAT_ZERO 1e-30

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

#define TAG_VALUE_RESOLUTION 1

#define PEAKRATIO_TOLERANCE_BP 5

#define PEAKFRACTION_REGION_MINDIST 4.0
#define REGION_MODE_HISTONE 0
#define REGION_MODE_GROSEQ 1
#define REGION_GROSEQ_FOLDCHANGE 3.0
#define GROSEQ_MAXBODYSIZE 10000

#define PEAK_STYLE_HISTONE 1

#define DEFAULT_GSIZE 2000000000

#define ZIPPED_FLAG_GZ 1
#define ZIPPED_FLAG_BZ2 2
#define ZIPPED_FLAG_ZIP 3
#define ZIPPED_FLAG_BAM 4

class TagLibrary;
class ChrTags;
class LinkedTag;
class LinkedPETag;
class Tag;
class PETag;
class Peak;
class ChrPeaks;
class PeakLibrary;


class PeakFinder {
public:
	char* name;
	char* directory;
    char* outputFileName;
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
	float maxtbpInput;
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
	char* excludePeaksFile;
	TagLibrary* tags;
	TagLibrary* input;
	double * fdrTable;
	int fdrSize;
	double * poissonTable;
	double tagsInPeaks;
    int regionFlag;
    char* extraheader;
    int maxBodySize;

    int style;

    double inputFold;
	double localFold;
	double clonalFold;
    char* cmd;
		
	PeakFinder();
	~PeakFinder();
	PeakLibrary* findPeaks();

    void addHeader(char*);

    void setCMD(char* name);
	void setDirectory(char* name);
	void setOutputFile(char* name);
	void setTagLibraries(TagLibrary* exp, TagLibrary* input);
	void setMaxTBP(double,double);
	void setGenomeSize(long long int);
	void determineMaxTBP();
    void approxFdrTable();
	void checkParameters();

    PeakLibrary* filterPeaks(PeakLibrary*);
};

class TagLibrary {
public:
	Hashtable* chrs;
    double totalTags;
	long long int totalPositions;
    char* name;
    char* directory;
    int medianTagsPerPosition;
	int fragmentLengthEstimate;
	int fragmentLengthSetting;
	long long int gsizeEstimate;
    double averageTagsPerPosition;
	double averageTagLength;
	double tbp;
	double parseAlignmentCpGMinValue;
    float maxtbp;
	float mintbp;
	double minmapq;
	int revStrand;
	int singleFile;
	int minReadLength;
	int maxReadLength;
	int peReadFlag;
    double manualTagTotal;
	double totalTagsMultiMappers;
	long long int totalPosMultiMappers;

    char* singleFilename;
	FILE* tagFile;
    char* restrictionSite;
    int sspeFlag;

    int pairedEndFlag;
	int mCflag;
    double localInteractionFraction;
	double interChrInteractionFraction;
	Hashtable* chrNames;
	
	TagLibrary(char* directory);
	~TagLibrary();

    void readAndSave();
	void parseAlignmentFiles(char** files, int numFiles, int format, int mode,
								char** tagDirs, int numDirs,char** tagFiles, int numTagFiles);
	void setMaxTBP(float maxTagsPerBp);
	void setMinTBP(float minTagsPerBp);
	void setSingleFile(int flag);
	void setTagAdjust(int centerDistance);
	void setFragLength(int fragLength);
	void setSingleRead(int singleReadFlag);
	void addAlignedPETag(PETag* petag);

    double* getTagCountDistribution(FILE* nfp, int &max);
	double* getTagLengthDistribution(FILE* nfp, int &max);

    double getAdjustedTagTotal();


    void printTagInfo();
	void printTagInfo(FILE*);

    void setName(char*);

    //internal functions
	void readAlignment(char* files, int format, int mode, int PEflag);
	void readPEAlignment(char* files, int format, int mode);
	void addTagDirectory(char* tagDir);
	void readSingleTagFile();
	void readNewTagFile(char* filename);

    void makeDirectory();

    int getRightCoordFromCIGAR(char* str, int dir, char* cigarCodes, int* cigarLens, int &numCodes, int &initLen);

    void addAlignedTag(char* name, char* chr,int pos,char dir,int length, float value,int PEflag);
	char* getTagFileName(char* chr);
	char* getDirFileName(char* filename);
	void optimizeTagFiles();

    PeakLibrary* findPutativePeaks(int peakSize, int minDist, char strand, float minCount);

    double* getPETagDistribution(int windowSize, int largeWindowSize, int largeResolution,
					char* outputPrefix, int &arrayLength);
};


class ChrTags {
public:
	Tag* tags;
	PETag* petags;
	int totalPositions;
	double totalTags;
    float* gcFreq;
	int mCflag;
	int singleFile;
    int loaded;
	int adjusted;
	int optimizedFlag;
	int optimizeOverride;
	float maxtbp;
	float mintbp;
	int tagAdjust;
	int dontSAVE;
	long long int appearentSize;
	int revStrand;
	char* chr;
	char* tagFile;
	FILE* tagfp;
	FILE* tagfpR1;
	FILE* tagfpR2;
	LinkedTag* firstTag;
	LinkedTag* linkTag;
	LinkedPETag* firstPETag;
	LinkedPETag* linkPETag;
	int numLinks;
	int pairedEndFlag;
	int forceSingleReadFlag;
    Hashtable* chrNames;

	char* seq;


    ChrTags(char* chr);
	~ChrTags();

	void adjustTags();
	void readAndSave();

    void readTagFile();

    void addTag(int pos, char dir,int length, float value);
	void addPETag(char* c1, int p1, char d1,int len1, char* c2, int p2, char d2, int len2,float value);
	void printAlignedTag(char* name, int pos, char dir,int length, float value,int PEflag);
	void printAlignedPETag(PETag* petag, int revFlag);
	void optimizeTags();
	void optimizePETags();
	void optimizeTagFile();
	void setMaxTBP(float maxTagsPerBp);
	void setMinTBP(float minTagsPerBp);
	void setTagAdjust(int centerDistance);
	void setTagFile(char* file);
	void openTagFile(char* mode);
	void closeTagFile();
	void print();

    void findPutativePeaks(PeakLibrary* putativePeaks, int peakSize, int minDist, char strand, double minCount);

    void getTagCountDistribution(double* d, int max,int scaleFactor);
	void getTagLengthDistribution(double* d, int max);

    int getPETagDistribution(double* sameStrand, double* diffStrand, int windowSize,
                            double * largeWindow,int resolution, int largeLength);

    void loadTags();
	void freeTags();

};
int cmpTags(const void*, const void*);
int cmpPETags(const void*, const void*);

class Tag {
public:
	int p; // position
	int len; // position
	float v; // value = number of tags
	char d; // direction 0=+,1=-
	Tag();
	void copy(Tag* src);
	void print(FILE* fp,char* chr);

	static int precision;
};

class LinkedTag : public Tag {
public:
	LinkedTag* tag;
	LinkedTag(int,char,int,float,LinkedTag*);
	~LinkedTag() = default;
};

class PETag {
public: 
	char* name;
	char* chr1;
	int p1; 
	char d1;
	int len1;
	char* chr2;
	int p2; 
	char d2;
	int len2;
	float v;
	PETag();
	PETag(char* name,char* chr, int p, char d, float v,int len);
	~PETag();
	void init();
	void copy(PETag* src);
	void print(FILE*);

    void print(FILE*,int revFlag);

	static int precision;
};


class LinkedPETag : public PETag {
public:
	LinkedPETag* tag;
	LinkedPETag(char* nc1, int np1, char nd1,int nlen1,char* nc2, int np2, char nd2,int,float,LinkedPETag*);
	~LinkedPETag() = default;
};

class PeakLibrary {
public:
	Hashtable* chrs; // holds pointers to ChrPeak objects
	Hashtable* peaks; // holds pointers to Peak objects
	char* name;
	char* genome;
	int numPeaks;
	double tagsInPeaks;
	double avgPeakSize;
	int fixedFlag;
	int duplicateWarningFlag;
	Peak** peakOrder;

    Inttable* duplicates;

    TagLibrary** exps;
	int numExps;


	PeakLibrary();
	PeakLibrary(char* file, int mode);
	PeakLibrary(int expectedNumberOfPeaks); // default is 100000
	void initialize(int expectedNumPeaks);
	~PeakLibrary();
	void print(FILE*);

    void readPeakFile(char* filename,int mode);
	void setDefaultPeakOrder();

    //pass NULL to trim peak positions using the current sizes
	int addTagLibrary(TagLibrary* t);

    PeakLibrary* stitchRegions(int maxDistance,int mode) const;
	PeakLibrary* getDifferentialPeaks(TagLibrary* tags, TagLibrary* input,
                        double foldThreshold, double poissonThresh, int mode, int start, int end, char strand,int strFlag);
	PeakLibrary* filterLocalPeaks(TagLibrary* tags, int peakSize, int localSize,
                        double threshold, double poissonThresh, int mode, char strand);
	PeakLibrary* filterClonalPeaks(TagLibrary* tags, int peakSize,
                        double threshold, int mode, char strand);

    void centerPeaks(TagLibrary* tags, int peakSize,char strand);

    void centerNFR(TagLibrary* tags, int peakSize,char strand,int nfrSize);

    Doubletable* countPeakTagsLowMemory(TagLibrary* tags, char direction,int mode,PeakLibrary* exclude);

    void setPeakTagSizeRefPos(int offset, int startOffset, int endOffset);
	void setPeakTagSizeFixed(int startOffset, int endOffset);

    void sortPeakTags(int expIndex);
	Peak* addPeak(char* name, char* chr,int start, int end, int midpoint, char dir, float value, 
						float ratio, char* extradata, int mappability, unsigned int priority);

    void normalizePeakScore(float normFactor);
	void addPeak(Peak* p);

    void sortChr();

    void sortKeys(char**);

};

class ChrPeaks {
public:	
	Peak** peaks;
	int numPeaks;
    LinkedList* peakList;

    ChrPeaks();
	~ChrPeaks();
	void addPeak(Peak* p);
	void sort();

    void addTagLibrary(ChrTags* t,int libraryIndex);

    void stitchRegions(PeakLibrary* regions, int maxDistance, int mode);

    void countPeakTagsLowMemory(Doubletable* results, ChrTags* ct,char direction, int mode,ChrPeaks* exclude);

};

class Peak {
public:
	char* name;
	char* chr;
	int refPos;
	int start;
	int end;
	int tagStart;
	int tagEnd;
    unsigned int priority;
	char strand;
	char* seq;
	int uniqMap;
	float v;
	float focusRatio;
	char* data;
	char* ogname;


    Tag** exps;
	int* numTags;
	int* maxTags;
	int numExps;

	Peak();

    Peak(char* name,char* originalName, char* chr, int start, int end,int midpoint, char dir, float value,
				float focusRatio, char* otherdata,int mappability,unsigned int priority);
	~Peak();

    Tag* getCoverageTags(int expIndex, int fragLength, int& coveragePositions,char strand);
	void centerPeak(int expIndex,int fragLength,char strand);

    void centerNFR(int expIndex,int fragLength,char strand,int nfrSize,int nucSize);

    void print(FILE* fp);

    void addExp();
	void addTag(Tag* t, int expIndex);
	void sortTags(int expIndex);

    void setPeakTagSizeRefPos(int newOffset,int startOffset,int endOffset);
	void setPeakTagSizeFixed(int startOffset,int endOffset);
	void setOffset(int offset);
	void addData(char* str);

};


// support classes

class doubleIndex {
public:
	double v;
	double vp;
	unsigned int index;
	int position;
};

int cmpDoubleIndex(const void* a, const void* b);
int cmpPeaks(const void*, const void*);

void split(char* string, char** cols, int &numCols, char delim);
int checkInt(char* str);
int checkStrand(char* str);

char* unzipFileIfNeeded(char* file, int &zipFlag, int &format);
void rezipFileIfNeeded(char* file, int zipFlag);

int chrcmp(const void* chr1, const void* chr2);


#endif
