cmake_minimum_required(VERSION 3.11)

project(
        homer-findPeaks
        DESCRIPTION "The peak caller for ChIP-Seq data"
        LANGUAGES CXX
)

set(LIBRARIES cpp/SeqTag.cpp cpp/statistics.cpp cpp/Hashtable.cpp)

# findPeaks
set(SOURCES cpp/findPeaks.cpp ${LIBRARIES})
add_executable(findPeaks ${SOURCES})
target_compile_options(findPeaks PRIVATE -Wall -Wextra -Wno-unused-parameter)
target_compile_options(findPeaks PRIVATE -std=c++17 -flto -O3)
target_link_options(findPeaks PRIVATE -flto)

# makeTagDirectory
set(SOURCES cpp/makeTagDirectory.cpp ${LIBRARIES})
add_executable(makeTagDirectory ${SOURCES})
target_compile_options(makeTagDirectory PRIVATE -Wall -Wextra -Wno-unused-parameter)
target_compile_options(makeTagDirectory PRIVATE -std=c++17 -flto -O3)
target_link_options(makeTagDirectory PRIVATE -flto)


# python bindings - pybind
add_subdirectory(third_party)

set(SOURCES pybindings/pyhomer.cpp ${LIBRARIES})
#
pybind11_add_module(pyhomer MODULE THIN_LTO ${SOURCES})
target_include_directories(pyhomer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp)

#target_compile_options(pyhomer PRIVATE -Wall -Wextra -Wno-unused-parameter)
#target_compile_options(pyhomer PRIVATE -std=c++17 -flto -O3)
