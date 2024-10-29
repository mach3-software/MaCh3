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
/// @details anything longer than ~80 should go in a details
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
If you want to mimic std::setw<<10
```cpp
MACH3LOG_INFO("Some break {:<10}", blarb);
```
Laslty if you want combine precision and std::setw like format
```cpp
MACH3LOG_INFO("Some break {:<10.2f}", blarb);
```

## Config Syntax
MaCh3 currently uses yaml as config handler. To help unify syntax over the code there are several YamlHelper function available [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/YamlHelper.h). Most important is `GetFromManager`. For code below which checks if config entry exist and if doesn't set some default value

```cpp
bool AsimovFit = false;

if(config[AsimovFit])
{
  AsimovFit = config[AsimovFit].as<bool>();
}
```
This can be replaced with:
```cpp
bool AsimovFit = GetFromManager<bool>(config[AsimovFit], false);
```

## double vs float?
Some fits require a lot of RAM. The easiest and fastest solution to reduce RAM
is to use `float` instead of `double`.

MaCh3 has a custom type defined as `M3::float_t`, which is usually a `double`
unless the `_LOW_MEMORY_STRUCTS_` directive is defined at the compilation
level. When defined, `M3::float_t` will be an actual `float`.

By using `M3::float_t`, one can flexibly change between these types. During
development, it is advised to use these data types unless specific data
types are necessary due to desired precision, code safety, etc.

## Error handling
MaCh3 uses custom error handling implemented [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/MaCh3Exception.h)

Never ever ever bare throw. Always throw an exception, preferably one that subclasses one defined by the standard library in `<stdexcept>`.
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
This way we can ensure error messages are unified and user always get hints where in the code problem occurred. If current `MaCh3Exception` is not sufficient consider implementing new or expanding current exceptions in MaCh3Exception.h.

## Compiler warning levels getting you down?

If you are trying to compile some new development and it is failing because of some innocuous warning that has been elevated to an error by the compiler flags, please don't just turn off the flags. A much better approach is to disable the diagnostic locally. This makes it easier to keep most of the code stringently checked, while giving you, the developer, the ability to stay in flow. It also allows for later 'fixing' of these warnings, if they need to be fixed, to be done systematically by greping for the relevant directives.

The way to turn off diagnostics is, as below:

```c++
#pragma GCC diagnostic ignored "-Wfloat-conversion"
```

N.B. that clang also understands these directives, so don't panic that they have GCC in the directive.

This will disable that diagnostic for the rest of the compilation unit (usually a .cc file). Note that this means if you include these in headerfiles, they will disable diagnostics more widely, please try and disable the diagnostics over as little code as possible.

### An example

We got this compiler error:

```
```

for this code:

```c++
```

We might update it to the below to disable the diagnostic just for the relevant line.

## Formatting
To ensure a unified style in MaCh3 software you can use a clang-format file which has instructions about formatting code.
```bash
clang-format --assume-filename=/path/to/your/.clang-format=${MaCh3_ROOT}/../.clang-format blarb.cpp
```
Please see [here](https://clang.llvm.org/docs/ClangFormat.html) and [here](https://root.cern/contribute/coding_conventions/) for more details.
