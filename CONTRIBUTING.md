# How to Contribute Code to MaCh3
The gitflow development policy should be followed, with the master and develop branches only being committed to via pull request.

New features should be developed on branches in this repository with the branch naming: `feature_FEATURENAME`

Please see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more details

## Doxygen
When making comments try following Doxygen type of comments

```cpp
// I like comments
void Foo(){}
```
try
```cpp
/// @brief I like comments
void Foo(){}
```
After making release or tag please
```bash
cd Doc
doxygen Doxyfile
```
to produce both html and latex file. To produce book like pdf file:
```bash
cd latex
pdflatex refman.tex
```

## Logger
MaCh3 is using spdlog logger see [here](https://github.com/gabime/spdlog/tree/master). And wraps it around MaCh3 names [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/MaCh3Logger.h)
It is advised to use it. For example rather than
```cpp
std::cout<< " Error: Something is wrong" <<std::endl;
```
it is advised to used logger
```cpp
MACH3LOG_ERROR("Something is wrong");
```
To pass argument to logger use following syntax:
```cpp
MACH3LOG_INFO("I like {}, do you like {}", FirsString, SecondString);
```
Logger by default will print whole float. Normally to show only several significant figures you would use `std::precision(2)`. To obtain same thing with logger please use `{:.2f}` like this:
```cpp
MACH3LOG_INFO("Here is full LLH but I only show 2 significant figures {:.2f}", LLH);
```

## Config Syntax
MaCh3 currently uses yaml as config handler. To help unify syntax over the code there are several YamlHelper function available [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/YamlHelper.h). Most important is `GetFromManager`. For code below which checks if config entry exist and if doesn't set some default value

```cpp
bool AsimovFit = false;

if(config[AsimovFit])
{
  AsimovFit = config[AsimovFit].as<bool>;
}
```
This can be replaced with:
```cpp
bool AsimovFit = GetFromManager<bool>(config[AsimovFit], false);
```

## double vs float?
Some fits require a lot of RAM. The easiest and fastest solution to reduce RAM
is to use `float` instead of `double`.

MaCh3 has a custom type defined as `_float_`, which is usually a `double`
unless the `_LOW_MEMORY_STRUCTS_` directive is defined at the compilation
level. When defined, `_float_` will be an actual `float`.

By using `_float_`, one can flexibly change between these types. During
development, it is advised to use these data types unless specific data
types are necessary due to desired precision, code safety, etc.

## Error handling
MaCh3 uses custom error handling implemented [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/MaCh3Exception.h)
Instead of throw
```cpp
throw;
```
it is advised to use:

```cpp
throw MaCh3Exception(__FILE__ , __LINE__ , "Some message");
```
or
```cpp
throw MaCh3Exception(__FILE__ , __LINE__ );
```
This way we can ensure error messages are unified and user always get hints where in the code problem occurred. If current MaCh3Exception is not sufficient consider implementing new or expanding current exceptions in MaCh3Exception.h.

## Formatting
To ensure a unified style in MaCh3 software you can use a clang-format file which has instructions about formatting code.
```bash
clang-format --assume-filename=/path/to/your/.clang-format=${MaCh3_ROOT}/../.clang-format blarb.cpp
```
Please see [here](https://clang.llvm.org/docs/ClangFormat.html) and [here](https://root.cern/contribute/coding_conventions/) for more details.
