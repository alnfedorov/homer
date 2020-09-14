#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "callPeaks.h"
#include <iostream>

namespace py = pybind11;

PYBIND11_MODULE(pyhomer, m) {
    m.doc() = "Python bindings for homer findPeaks algorithm.";

    py::class_<TagLibrary>(m, "TagLibrary")
            .def("__eq__", [](const TagLibrary &t1, const TagLibrary &t2) {
                if (t1.pairedEndFlag != t2.pairedEndFlag || t1.sspeFlag != t2.sspeFlag || t1.maxtbp != t2.maxtbp ||
                    t1.fragmentLengthEstimate != t2.fragmentLengthEstimate ||
                    t1.averageTagLength != t2.averageTagLength ||
                    t1.averageTagsPerPosition != t2.averageTagsPerPosition || t1.totalTags != t2.totalTags ||
                    t1.totalPositions != t2.totalPositions || t1.gsizeEstimate != t2.gsizeEstimate ||
                    t1.tbp != t2.tbp || t1.peakSizeEstimate != t2.peakSizeEstimate || t1.chrs.size() != t2.chrs.size())
                    return false;

                for (const auto &it1: t1.chrs) {
                    auto it2 = t2.chrs.find(it1.first);
                    if (it2 == t2.chrs.end())
                        return false;

                    auto &chr1 = *it1.second;
                    auto &chr2 = *it2->second;
                    if (chr1.forceSingleReadFlag != chr2.forceSingleReadFlag ||
                        chr1.pairedEndFlag != chr2.pairedEndFlag ||
                        chr1.chr != chr2.chr || chr1.appearentSize != chr2.appearentSize ||
                        chr1._total_tags_for_maxtbp != chr2._total_tags_for_maxtbp ||
                        chr1._tags_size_for_maxtbp != chr2._tags_size_for_maxtbp)
                        return false;

                    auto eq = [](const Tag &t1, const Tag &t2) {
                        return t1.v == t2.v && t1.len == t2.len && t1.d == t2.d && t1.p == t2.p;
                    };

                    for (auto &it1: chr1._tags_for_maxtbp) {
                        auto it2 = chr2._tags_for_maxtbp.find(it1.first);
                        if (it2 == chr2._tags_for_maxtbp.end())
                            return false;

                        auto size = chr1._tags_size_for_maxtbp[it1.first];
                        for (size_t i = 0; i < size; ++i)
                            if (!eq(it1.second[i], it2->second[i]))
                                return false;
                    }
                }
                return true;
            })
            .def(py::pickle(
                    [](const TagLibrary &p) { // __getstate__
                        static auto pickle_chrtags = [](ChrTags *tags) {
                            auto result = py::list();
                            result.append(tags->chr);
                            result.append(tags->pairedEndFlag);
                            result.append(tags->appearentSize);
                            result.append(tags->maxtbp);
                            result.append(tags->forceSingleReadFlag);
                            result.append(tags->totalTags);

                            for (const auto &it: tags->_tags_for_maxtbp) {
                                auto ptr = it.second;
                                auto size = tags->_tags_size_for_maxtbp[it.first];
                                auto total = tags->_total_tags_for_maxtbp[it.first];
                                auto capsule = py::capsule((void *) ptr, [](void *vec) {});
                                auto array = py::array(size * sizeof(Tag), (uint8_t *) ptr, std::move(capsule));
                                result.append(py::make_tuple(it.first, size, total, array));
                            }
                            return result;
                        };

                        auto l = py::list();
                        l.append(p.totalTags);
                        l.append(p.totalPositions);

                        l.append(p.tbp);
                        l.append(p.maxtbp);

                        l.append(p.minmapq);
                        l.append(p.minReadLength);
                        l.append(p.maxReadLength);

                        l.append(p.fragmentLengthEstimate);
                        l.append(p.peakSizeEstimate);

                        l.append(p.averageTagLength);
                        l.append(p.averageTagsPerPosition);
                        l.append(p.gsizeEstimate);

                        l.append(p.sspeFlag);
                        l.append(p.pairedEndFlag);

                        for (const auto &it: p.chrs) {
                            auto tags = pickle_chrtags(it.second);
                            l.append(py::make_tuple(it.first, tags));
                        }
                        return l;
                    },
                    [](const py::list &l) { // __setstate__
                        static auto unpickle_chrtags = [](const py::list &chrtags) {
                            auto ptr = new ChrTags(chrtags[0].cast<std::string>());
                            ptr->pairedEndFlag = chrtags[1].cast<bool>();
                            ptr->appearentSize = chrtags[2].cast<long long>();
                            ptr->maxtbp = chrtags[3].cast<float>();
                            ptr->forceSingleReadFlag = chrtags[4].cast<bool>();
                            ptr->totalTags = chrtags[5].cast<double>();

                            std::vector<py::object> arrays;
                            for (size_t i = 6; i < chrtags.size(); ++i) {
                                auto tuple = chrtags[i].cast<py::tuple>();
                                auto maxtbp = tuple[0].cast<int>();
                                auto size = tuple[1].cast<size_t>();
                                auto total = tuple[2].cast<double>();
                                auto array = tuple[3].cast<py::array>();

                                ptr->_total_tags_for_maxtbp[maxtbp] = total;
                                ptr->_tags_size_for_maxtbp[maxtbp] = size;
                                ptr->_tags_for_maxtbp[maxtbp] = (Tag *) array.mutable_data();
                                arrays.push_back(array.base());
                            }
                            ptr->callback = [arrays, ptr]() {
                                // they will be dropped by arrays destructor
                                for (auto& it: ptr->_tags_for_maxtbp)
                                    it.second = NULL;
                            };
                            return ptr;
                        };

                        auto result = new TagLibrary();
                        result->totalTags = l[0].cast<double>();
                        result->totalPositions = l[1].cast<long long int>();
                        result->tbp = l[2].cast<double>();
                        result->maxtbp = l[3].cast<float>();
                        result->minmapq = l[4].cast<double>();
                        result->minReadLength = l[5].cast<size_t>();
                        result->maxReadLength = l[6].cast<size_t>();
                        result->fragmentLengthEstimate = l[7].cast<int>();
                        result->peakSizeEstimate = l[8].cast<int>();
                        result->averageTagLength = l[9].cast<double>();
                        result->averageTagsPerPosition = l[10].cast<double>();
                        result->gsizeEstimate = l[11].cast<long long>();
                        result->sspeFlag = l[12].cast<bool>();
                        result->pairedEndFlag = l[13].cast<bool>();

                        for (size_t i = 14; i < l.size(); ++i) {
                            auto tuple = l[i].cast<py::tuple>();
                            auto chr = tuple[0].cast<std::string>();
                            auto ptr = unpickle_chrtags(tuple[1].cast<py::list>());
                            result->chrs[chr] = ptr;
                        }
                        return result;
                    }
            ));

    py::class_<PeakFinder>(m, "PeakFinder")
            .def(py::init())
            .def_property("stitchMode",
                          [](const PeakFinder &finder) { return "histone"; },
                          [](PeakFinder &finder, const std::string &mode) {
                              if (mode == "histone") finder.stitchMode = REGION_MODE_HISTONE;
                          })
            .def_property("style",
                          [](const PeakFinder &finder) { return "histone"; },
                          [](PeakFinder &finder, const std::string &mode) {
                              if (mode == "histone") finder.style = PEAK_STYLE_HISTONE;
                          })
            .def_property("extendRegions",
                          [](const PeakFinder &finder) { return finder.regionFlag; },
                          [](PeakFinder &finder, const bool &flag) { finder.regionFlag = flag; })
            .def_property("peakSize",
                          [](const PeakFinder &finder) { return finder.peakSize; },
                          [](PeakFinder &finder, const int &size) { finder.peakSize = size; })
            .def_property("minDist",
                          [](const PeakFinder &finder) { return finder.minDist; },
                          [](PeakFinder &finder, const int &size) { finder.minDist = size; })
            .def_property("localSize",
                          [](const PeakFinder &finder) { return finder.localSize; },
                          [](PeakFinder &finder, const int &size) { finder.localSize = size; })
            .def_property("inputSize",
                          [](const PeakFinder &finder) { return finder.inputSize; },
                          [](PeakFinder &finder, const int &size) { finder.inputSize = size; })
            .def_property("fdr",
                          [](const PeakFinder &finder) { return finder.fdr; },
                          [](PeakFinder &finder, const double &fdr) {
                              finder.fdr = fdr;
                              finder.filterMode = PEAKFINDER_FILTER_MODE_FDR;
                          })
            .def_property("poisson",
                          [](const PeakFinder &finder) { return finder.poisson; },
                          [](PeakFinder &finder, const double &poisson) {
                              finder.poisson = poisson;
                              finder.filterMode = PEAKFINDER_FILTER_MODE_POISSON;
                          })
            .def_property("expectedFold",
                          [](const PeakFinder &finder) { return finder.clonalFold; },
                          [](PeakFinder &finder, const double &fold) { finder.clonalFold = fold; })
            .def_property("poissonLocal",
                          [](const PeakFinder &finder) { return finder.poissonLocal; },
                          [](PeakFinder &finder, const double &cutoff) { finder.poissonLocal = cutoff; })
            .def_property("localFold",
                          [](const PeakFinder &finder) { return finder.localFold; },
                          [](PeakFinder &finder, const double &cutoff) { finder.localFold = cutoff; })
            .def_property("inputFold",
                          [](const PeakFinder &finder) { return finder.inputFold; },
                          [](PeakFinder &finder, const double &cutoff) { finder.inputFold = cutoff; })
            .def_property("poissonInput",
                          [](const PeakFinder &finder) { return finder.poissonInput; },
                          [](PeakFinder &finder, const double &cutoff) { finder.poissonInput = cutoff; })
            .def_property("saveto",
                          [](const PeakFinder &finder) { return finder.outputFileName; },
                          [](PeakFinder &finder, const std::string &savetp) { finder.setOutputFile(savetp); })
            .def("setTagLibraries",
                 [](PeakFinder &finder, TagLibrary *treatment, TagLibrary *control) {
                     finder.setTagLibraries(treatment, control);
                 })
            .def("findPeaks", [](PeakFinder &finder) { findPeaks(finder); });
    m.def("load", &load);

}