## How to Contribute Code to MaCh3
The gitflow development policy should be followed, with the master and develop branches only being committed to via pull request.

New features should be developed on branches in this repository with the branch naming: `feature_FEATURENAME`

Please see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more details

## Doxygen
TODO

## Formatting
To ensure a unified style in MaCh3 software you can use a clang-format file which has instructions about formatting code.
```
clang-format --assume-filename=/path/to/your/.clang-format=${MaCh3_ROOT}/../.clang-format blarb.cpp
```
Please see [here](https://clang.llvm.org/docs/ClangFormat.html) and [here](https://root.cern/contribute/coding_conventions/) for more details.

## Logger
MaCh3 is using spdlog logger see [here](https://github.com/gabime/spdlog/tree/master). It is advised to use it. For example rather than
```
std::cout<< " Error: Something is wrong" <<std::endl;
```
it is advised to used logger
```
SPDLOG_ERROR("Something is wrong");
```
