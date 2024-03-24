## How to Contribute Code to MaCh3
The gitflow development policy should be followed, with the master and develop branches only being commited to via pull request.

New features should be developed on branches in this repository with the branch naming: `feature_FEATURENAME`

Please see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more details

## Doxygen
TODO

## Formating
To ensure a unified style in MaCh3 software you can use a clang-format file which has instructions about formatting code.
```
clang-format --assume-filename=/path/to/your/.clang-format=${MaCh3_ROOT}/../.clang-format blarb.cpp
```
Please see [here](https://clang.llvm.org/docs/ClangFormat.html) and [here](https://root.cern/contribute/coding_conventions/) for more details.
