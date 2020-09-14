This repository contains a stripped version of the homer v4.10.1 with python bindings for findPeaks, makeTagDirectory routines.

Please, refer to the official homer software [page](http://homer.ucsd.edu/homer/) for the fresh homer sources.  
**Don't use this repository unless you know what you are doing.**  

Fork changelog:
1. All .pl files are dropped.
2. CMakeLists to build python bindings.
3. All non relevant to peak calling functionality is dropped.
4. Refactoring.

The final algorithm closely resembles the original one (IOU for called peaks is about 0.95-1.0).
At the moment, paired-end reads are treated as a single-end. Homer can parse paired-end reads, but cant shift them properly later(same as the original algorithm).
