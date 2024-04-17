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
MaCh3 is using spdlog logger see [here](https://github.com/gabime/spdlog/tree/master). And wraps it around MaCh3 names [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/MaCh3Logger.h)
It is advised to use it. For example rather than
```
std::cout<< " Error: Something is wrong" <<std::endl;
```
it is advised to used logger
```
MACH3LOG_ERROR("Something is wrong");
```
To pass argument to logger use following syntax:
```
MACH3LOG_INFO("I like {}, do you like {}", FirsString, SecondString);
```
Logger by default will print whole float. Normally to show only several significant figures you would use `std::precision(2)`. To obtain same thing with logger please use `{:.2f}` like this:
```
MACH3LOG_INFO("Here is full LLH but I only show 2 significant figures {:.2f}", LLH);
```
