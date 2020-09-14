#include "callPeaks.h"

#include "iostream"
int main() {
    std::cout << std::boolalpha;
    std::cout << std::is_pod<Tag>::value << std::endl;
//    std::vector <std::string> tfiles = {
//            "/home/aleksander/data/bioinformatics/encode/H3K4me3/ENCSR361FWQ/original/filtered-bam/ENCFF719TTA-h3k4me3.bam"};
    std::vector<std::string> tfiles = {"/home/aleksander/projects/homer/data/raw/se/treatment.bam"};
    auto treatment = load(tfiles);

//    std::vector <std::string> cfiles = {
//            "/home/aleksander/data/bioinformatics/encode/H3K4me3/ENCSR361FWQ/original/filtered-bam/ENCFF552BJI-control.bam"};
    std::vector<std::string> cfiles = {"/home/aleksander/projects/homer/data/raw/se/control.bam"};
    auto control = load(cfiles);

    std::string saveto = "/home/aleksander/data/bioinformatics/callPeaks.bed";
    PeakFinder pf;
    pf.stitchMode = REGION_MODE_HISTONE;
    pf.style = PEAK_STYLE_HISTONE;
    pf.regionFlag = 1;
    pf.peakSize = 500;
    pf.minDist = 1000;
    pf.poisson = 0.001;
    pf.filterMode = PEAKFINDER_FILTER_MODE_POISSON;
    pf.setTagLibraries(treatment, control);
    pf.setOutputFile(saveto);
    findPeaks(pf);
    delete treatment;
    delete control;
}
