# How to Contribute Code to MaCh3
The gitflow development policy should be followed, with the master and develop branches only being committed to via pull request.

New features should be developed on branches in this repository with the branch naming: `feature_FEATURENAME`

Please see [here](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more details

## Versioning of this package

This package will follow this numbering convention for releases
- Major `X.0.0` for paradigm-shifting breaking changes
- Minor `0.X.0` for more minor breaking changes (or lots of non-breaking changes + a long time since last release)
- Patch `0.0.X` for non-breaking changes

## PR titles
When creating a Pull Request (PR), make sure the title starts with one of the following prefixes to indicate the type of PR:

- **`feat:`**: For new features
- **`bugfix:`** For fixing bugs
- **`tidy:`** For minor changes like refactoring or updating dependencies
- **`breaking:`** For changes that break backward compatibility

## Continuous Integration

MaCh3 undergoes rigorous checks via multiple continuous integration (CI) processes. You can find the full list of CI checks [here](https://github.com/mach3-software/MaCh3/tree/develop/.github/workflows#readme).

CI heavily depends on the setup in the [MaCh3 Tutorial](https://github.com/mach3-software/MaCh3Tutorial).

### PR Merge Requirements
- **CI Success**: All CI checks must pass before a PR can be merged.
- **Backward Compatibility**: If a PR introduces breaking changes, it will fail the CI tests. In this case, the following steps are required:
    1. **Update MaCh3 Tutorial**: A separate PR must be made to update the [MaCh3 Tutorial](https://github.com/mach3-software/MaCh3Tutorial).
    2. **Merge Tutorial PR**: The tutorial update PR must be merged first.
    3. **CI Restart**: After merging the tutorial update, restart the CI tests until they pass.

This process ensures that:
- Breaking changes are properly registered.
- The [MaCh3 Tutorial](https://github.com/mach3-software/MaCh3Tutorial) remains up-to-date.

### Monitoring Jobs
Some bots are tasked with weekly or long-term monitoring tasks.
For example:
- [Valgrind outputs](https://github.com/mach3-software/MaCh3/tree/gh-profiling) are updated on a weekly basis to track memory-related issues and visualize the CMake dependency graph.
- [Benchmarking results](https://mach3-software.github.io/MaCh3Tutorial/) are stored and published periodically to monitor performance metrics over time.

These outputs are primarily intended for expert developers. However, they **should still be regularly monitored** to identify:
- Persistent memory leaks
- Performance regressions
- Issues in the CMake linking flow

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
You can
```cpp
/// @param blarb this is parameter needed for function
/// @return returning some fancy number
int Fool(int blarb){}
```
Doxygen documentation will be updated after each commit to develop.

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
MaCh3 currently uses yaml as config handler. To help unify syntax over the code there are several YamlHelper function available [here](https://github.com/mach3-software/MaCh3/blob/develop/manager/YamlHelper.h).

We suggest instead of using default yaml syntax
```cpp
bool AsimovFit = config[AsimovFit].as<bool>();
```
is to use `Get`:
```cpp
bool AsimovFit = Get<bool>(config[AsimovFit], __FILE__ , __LINE__);
```
`Get` has additional protection, checks if node exist, checks conversion and throw proper errors. This allow to troubleshoot code much easier.

Another useful one is `GetFromManager`. For code below which checks if config entry exist and if doesn't set some default value
```cpp
bool AsimovFit = false;

if(config[AsimovFit])
{
  AsimovFit = config[AsimovFit].as<bool>();
}
```
This can be replaced with:
```cpp
bool AsimovFit = GetFromManager<bool>(config[AsimovFit], false, __FILE__, __LINE__);
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

If you are trying to compile some new development and it is failing because of some innocuous warning that has been elevated to an error by the compiler flags, please don't just turn off the flags. A much better approach is to disable the diagnostic locally. This makes it easier to keep most of the code stringently checked, while giving you, the developer, the ability to stay in flow.
It also allows for later 'fixing' of these warnings, if they need to be fixed, to be done systematically by greping for the relevant directives.

The way to turn off diagnostics is, as below:

```c++
#pragma GCC diagnostic ignored "-Wfloat-conversion"
```

N.B. that clang also understands these directives, so don't panic that they have `GCC` in them.

This will disable that diagnostic for the rest of the compilation unit (usually a .cc file). Note that this means if you include these in headerfiles, they will disable diagnostics more widely, please try and disable the diagnostics over as little code as possible.

If a specific error is really getting you down and its showing up everywhere, the serious option is to disable it repo-wide by modifying the `MaCh3Warnings` interface target, defined in the top-level project [CMakeLists.txt](../CMakeLists.txt) like so:

```cmake
target_compile_options(MaCh3Warnings INTERFACE
    # ...
    -Wno-conversion
    # ...
)
```

Please attempt more localised options before reaching for this, but sometimes this represents the best way to proceed with development without 'fixing' innocuous warnings.

The really serious option is to configure with: `-DMaCh3_WERROR_ENABLED=OFF`, which will disable the `-Werror` flag.

### An example

We got this compiler error:

```shell
/root/software/MaCh3/MaCh3_splitpr/splines/splineFDBase.cpp: In member function ‘virtual void splineFDBase::CalcSplineWeights()’:
/root/software/MaCh3/MaCh3_splitpr/splines/splineFDBase.cpp:349:35: error: useless cast to type ‘double’ [-Werror=useless-cast]
  349 |     weightvec_Monolith[iSpline] = double(weight);
      |                                   ^~~~~~~~~~~~~~
cc1plus: all warnings being treated as errors
```

for this code:

```c++
    weightvec_Monolith[iSpline] = double(weight);
```

The compiler is right, that this is a useless cast, but `weight` can sometimes be a float, in which case we would get a conversion warning/error, so it seems like a no-win situation. We can 'save' the current diagnostics with `#pragma GCC diagnostic push`, disable the relevant one as above, and then revert to the saved set of diagnostics with `#pragma GCC diagnostic pop`.
Putting it all together might look like:

```c++
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
    weightvec_Monolith[iSpline] = double(weight);
#pragma GCC diagnostic pop
```

This allows us to disable the diagnostic just for the relevant line.

## C++ Standard
As of May 2025, the minimum supported C++ standard for MaCh3 is **C++14**, although the codebase is compatible with newer standards such as **C++17** and beyond.
These are often referred to collectively as *modern C++*.
MaCh3 aims to use modern C++ features to encourage safer, cleaner, and more maintainable code.

### Pointers
Instead of using raw pointer
```c++
    Class* example = new Class;

    // ... some code ...
    delete example;
```
Use *std::unique_ptr* to manage memory automatically:
```c++
    auto example = std::make_unique<Class>();
```
This eliminates the need to manually delete the object and helps prevent memory leaks.

### Dynamic Arrays
Instead of manually managing dynamic arrays:

```c++
    double* example = new double[N];

    // ... some code ...

    delete[] example;
```
Use a *std::vector* to automatically manage the memory:
```c++
    std::vector<double> example(N);
```
### Use `nullptr` Instead of `NULL` or `0`

Modern C++ introduces `nullptr` as a type-safe null pointer constant.
Avoid legacy `NULL` or `0` which can lead to ambiguous or unsafe behavior.
```c++
    int* ptr = nullptr;

    if (ptr == nullptr) {
        std::cout << "Pointer is null." << std::endl;
    }
```
### Casting

Use `static_cast<T>(value)` instead of `(T)value` for type conversions.
It makes the intent clear and avoids unsafe conversions allowed by C-style casts.
```c++
    double x = 3.14;
    int y = static_cast<int>(x);
```

## Formatting
To ensure a unified style in MaCh3 software you can use a clang-format file which has instructions about formatting code.
```bash
clang-format --assume-filename=/path/to/your/.clang-format=${MaCh3_ROOT}/../.clang-format blarb.cpp
```
Please see [here](https://clang.llvm.org/docs/ClangFormat.html) and [here](https://root.cern/contribute/coding_conventions/) for more details.

## Code-style
We want to try to make code more readable in general so please try to stick to using name style and conventions given in the google style guide [here](https://google.github.io/styleguide/cppguide.html). A few examples are given below:

- When naming function please capitalise the first letter and use CamelCase.
- For member variables try to use snake case
- Don't unnecessarily abbreviate variable names for the sake of saving a few characters
```c++
void GetMyFaviouriteVariable();
int my_name_member_variable;
int parameter_counter; //Not par_count
```

Please note that a lot of this has not been followed in MaCh3 but we are trying to improve the coding style so please bare this in mind when contributing!
