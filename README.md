# MaCh3 <img src="Doc/mach3logo.png" alt="MaCh3" align="center" width="100"/>

The Markov Chain 3 flavour is a framework born in 2013 as a Bayesian MCMC
fitter for [T2K](https://t2k-experiment.org/pl/) oscillation analysis. It has
now been used for multiple T2K Oscillation analyses at both the Near and Far
detectors throughout the years. The framework is also utilized by the
[DUNE](https://www.dunescience.org/) and [HK](https://www-sk.icrr.u-tokyo.ac.jp/en/hk/)
oscillation analysis groups. Additionally, it supports joint fits between T2K
and NOvA, as well as T2K and SK's atmospheric data.

The framework has also evolved to allow non-MCMC modules to interrogate the
likelihoods implemented.

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://github.com/mach3-software/MaCh3/blob/develop/LICENSE.txt)
[![DOI](https://zenodo.org/badge/331049416.svg)](https://doi.org/10.5281/zenodo.7608367)
[![Release](https://img.shields.io/github/release/mach3-software/MaCh3.svg)](https://github.com/mach3-software/MaCh3/releases/latest)
[![Container Image](https://img.shields.io/badge/Container-Image-brightgreen)](https://github.com/mach3-software/MaCh3/pkgs/container/mach3)
[![Code - Documented](https://img.shields.io/badge/Code-Documented-2ea44f)](https://github.com/mach3-software/MaCh3/wiki)
[![Code - Doxygen](https://img.shields.io/badge/Code-Doxygen-2ea44f)](https://mach3-software.github.io/MaCh3/index.html)
[![Build CI](https://github.com/mach3-software/MaCh3/actions/workflows/CIBuild.yml/badge.svg)](https://github.com/mach3-software/MaCh3/actions/workflows/CIBuild.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/mach3-software/mach3/badge/develop)](https://www.codefactor.io/repository/github/mach3-software/mach3/overview/develop)

## Famous Plots
Example of plots made using MaCh3 apparent in scientific publications, for more see [here](https://github.com/mach3-software/MaCh3/wiki/14.-MaCh3-in-the-Field)
<img src="Doc/Plots/delta.png" alt="MaCh3" align="left" width="200"/>
<img src="Doc/Plots/Jarlskog.png" alt="MaCh3" align="center" width="200"/>

## Cite
When using MaCh3 you must cite our doi from Zenodo. The bibtex file can be found by exporting the citation from this link: [on Zenodo](https://zenodo.org/records/7608367) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7608367.svg)](https://doi.org/10.5281/zenodo.7608367).

## Help and Guidelines 📄
- [Tutorial](https://github.com/mach3-software/MaCh3Tutorial)
- [Indico](https://indico.global/category/1289/) If you need a password, please reach out to MaCh3-leadership for access.
- [How to contribute](https://github.com/mach3-software/MaCh3/blob/develop/.github/CONTRIBUTING.md)
- [Wiki](https://github.com/mach3-software/MaCh3/wiki)
- [Slack](https://t2k-experiment.slack.com/archives/C06EM0C6D7W/p1705599931356889)
- [Discussions](https://github.com/mach3-software/MaCh3/discussions)
- [Benchmark](https://mach3-software.github.io/MaCh3Tutorial/)

## How to Compile
MaCh3 follows standard cmake pattern. By default you should get most optimal build setting although below we list many configurable options:
```bash
mkdir build;
cd build;
cmake ../
make -jN [Where N is number of threads]
make install
```

Don't forget to:
```bash
source bin/setup.MaCh3.sh
```

If the build fails, check your [system/environment requirements](#system-requirements).

## Building against MaCh3
To include MaCh3 in your cmake project you can use following syntax
```cmake
CPMFindPackage(
  NAME MaCh3
  GIT_TAG "blarb"
  GITHUB_REPOSITORY mach3-software/MaCh3
)
```
Where "blarb" is the MaCh3 version. You can find a list of releases [here](https://github.com/mach3-software/MaCh3/wiki/0.1.-History)  
If you compiled MaCh3 and sourced it you can simply call
```cmake
find_package(MaCh3)
```

Once you found MaCh3 you might want to link your library against MaCh3. You can do this as follows:
```cmake
target_link_libraries(blarb MaCh3::All)
```

Some functionalities rely on setting `Env{MACH3}` which should point to path experiment specific MaCh3. This way MaCh3 can easily find `Env{MACH3}/inputs/SomeInput.root` for example.

## Python 🐍

MaCh3 has an optional python interface (pyMaCh3) which provides much of the same functionality as the c++ interface (see [here](https://mach3-software.github.io/MaCh3/pyMaCh3/mainpage.html) for documentation).

You can tell the build system to set up the pyMaCh3 interface by specifying

```bash
cmake ../ -DMaCh3_PYTHON_ENABLED=ON
make && make install
```

when building

### Building with Pip

Additionally, you can build just the Python module by doing:

```bash
pip install -t <install location> .
```
The (optional) -t option specifies an install location which can be useful if you are on a computing cluster and don't have write access to the default install location. If you specify a non-standard location you will need to add it to your `PYTHONPATH` as above so that python can find the module.

## Multithreading
MaCh3 quite heavily relies on Multithreading, it is turned on by default. If for debugging purposes you would like to turn it off please use
```bash
cmake ../ -DMaCh3_MULTITHREAD_ENABLED=OFF
```

## CUDA
If the system has access to GPU, MaCh3 will enable GPU functionality automatically. If you would like to CPU only despite having access to [CUDA](https://developer.nvidia.com/cuda-toolkit)
```bash
mkdir build; cd build;
cmake ../ -DMaCh3_GPU_ENABLED=OFF
```
MaCh3 supports quite a high range of CUDA architectures if something doesn't work on your GPU let us know. MaCh3 supports only NVIDIA GPUs.

## Oscillator
MaCh3 has access to several neutrino oscillation calculators via [NuOscillator](https://github.com/dbarrow257/NuOscillator/tree/main) framework.

Following neutrino oscillation calculators are available:
|Oscillator        | Hardware   | Source     | Reference  |
|------------------|------------|------------|------------|
| CUDAProb3Linear  | CPU/GPU    | Beam       |            |
| CUDAProb3        | CPU/GPU    | Atm        | [Ref](https://doi.org/10.1016/j.cpc.2018.07.022)        |
| ProbGPULinear    | GPU        | Beam       | [Ref](http://dx.doi.org/10.3204/DESY-PROC-2014-05/23)   |
| Prob3++Linear    | CPU        | Beam       |            |
| NuFastLinear     | CPU        | Beam       | [Ref](https://doi.org/10.48550/arXiv.2405.02400)        |
| OscProb          | CPU        | Beam/Atm   | [Ref](https://doi.org/10.5281/zenodo.6347002)           |
| NuSQUIDSLinear   | CPU        | Beam       | [Ref](https://doi.org/10.1016/j.cpc.2022.108346)        |
| GLoBESLinear     | CPU        | Beam       | [Ref](https://doi.org/10.1016/j.cpc.2005.01.003)        |

If nothing is specified in cmake build then NuFastLinear_ENABLED will be used. To control which oscillation calculators you want to use here is syntax:

```bash
cmake ../ -DCUDAProb3Linear_ENABLED=ON -DCUDAProb3_ENABLED=ON -DProbGPULinear_ENABLED=ON -DProb3ppLinear_ENABLED=ON -DNuFastLinear_ENABLED=ON -DOscProb_ENABLED=ON
```
You can specify more than one engine, and MaCh3 support different engines for different samples.
For example, you can use NuFast for beam samples and CUDAProb3 for atmospheric samples.

## Fitting algorithms
The following fitting algorithms are available:

| Algorithm    | Reference        |Need Ext Lib  |
|--------------|------------------|--------------|
| MR2T2        | [Ref](https://doi.org/10.1063/1.1699114)          | No       |
| DelayedMR2T2 | [Ref](https://doi.org/10.1007/s11222-006-9438-0)  | No       |
| MINUIT2      | [Ref](https://cds.cern.ch/record/2296388/)        | Yes      |
| PSO          | [Ref](https://doi.org/10.1162/EVCO_r_00180)       | No       |

## Debug
Several debugging options are available which are heavy for RAM and performance and, therefore not used by default. To enable it:
```bash
cmake ../ -DMaCh3_DEBUG_ENABLED=<ON,OFF>
```
There are several debug modes, to enable more detailed but very heavy specific debug levels. Level 1 is the default debug activated by the above.

```bash
cmake ../ -DMaCh3_DEBUG_ENABLED=<ON,OFF> -DDEBUG_LEVEL=<1,2,3>
```
## Logger
Several logging options are available. By default, MaCh3 uses `INFO` level logging.
You can enable more verbose logging by specifying the desired level during configuration:
```bash
cmake ../ -DLOG_LEVEL="TRACE"
```
You can find more [here](https://github.com/mach3-software/MaCh3/blob/develop/cmake/Modules/Logger.cmake).

## Other CMake Options

| Option                               | Meaning                                                                         |
| ------                               | -------                                                                         |
| `MaCh3_NATIVE_ENABLED`               | Enables native CPU optimizations for improved performance. Not recommended on clusters with multiple CPU configurations due to potential compatibility issues.   |
| `MaCh3_NuOsc_GPU_ENABLED`            | By default MaCh3 will use NuOscillator with GPU if MaCh3 is compiled with GPU, this flag allows disabling GPU for NuOscillator even if MaCh3 has GPU enabled     |
| `MaCh3_LOW_MEMORY_STRUCTS_ENABLED`   | This will use float/short int for many structures |


## System Requirements
MaCh3 requires a C++ compiler (e.g. [gcc](https://gcc.gnu.org)), [CMake](https://cmake.org), and [ROOT](https://root.cern/). Based on several tests, recommended versions are:
```bash
  GCC:   >= 8.5   [lower versions may work]
  C++:   >= 14
  CMake: >= 3.14
  ROOT:  >= 6.20
```

If you don't already have the correct dependencies, you can install them in a variety of ways, e.g. building from source, via a package manager. If you have access to `cvmfs`, the easiest solution is probably to source an [LCG release](https://lcgdocs.web.cern.ch/lcgdocs/lcgreleases/introduction/). For example, on an x86_64 Alma9 machine:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/setup.sh
```
Alternatively, you can create a conda environment which *should* provide a ready-to-go environment. If you have micromamba, you can do:
```bash
micromamba env create -n MaCh3 -c conda-forge root cmake
```
and if you don't have micromamba installed, you can follow the instructions [here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

Before you build MaCh3 following the instructions at the [start](#building-against-mach3) of this README, make sure to activate the environment:
```bash
micromamba activate MaCh3
```
and do this everytime you begin from a fresh terminal and want to use MaCh3.

> [!TIP]
> For MacOS users: the conda environment approach has been validated to work for MacOS. Other avenues are likely to be more difficult.

Other external libraries, which currently include:
1. [yaml-cpp](https://github.com/jbeder/yaml-cpp)
2. [spdlog](https://github.com/gabime/spdlog)

are being handled through [CPM](https://github.com/cpm-cmake/CPM.cmake).

Based on several test here are recommended version:
```bash
  GCC:   >= 8.5   [lower versions may work]
  C++:   >= 14
  CMake: >= 3.14
  ROOT:  >= 6.18
```
### Supported operational systems
| Name        | Status |
|-------------|--------|
| Alma9       | ✅     |
| Rocky9      | ✅     |
| Ubi9        | ✅     |
| Ubuntu22.04 | ✅     |
| Ubuntu25.04 | ✅     |
| Fedora32    | ✅     |
| CentOS7     | ❔     |
| MacOS       | ❔     |
| Windows     | ❌     |

✅ - Part of CI/CD <br>
❔ - Not part of CI/CD but used by some users/developers so it might work <br>
❌ - Not supported and no plans right now <br>

### Supported compilers
| Name        | Status |
|-------------|--------|
| GNU         | ✅     |
| CLANG       | ❔     |
| INTELLLVM   | ❔     |
| MSVC        | ❌     |

✅ - Fully working with every feature fully tested by CI/CD <br>
❔ - Not every feature may work, only compilation being tested by CI/CD <br>
❌ - Not supported and no plans right now <br>

## Plotting and Diagnostic 📊
Example of chain diagnostic utils can be found [here](https://github.com/mach3-software/MaCh3/tree/develop/Diagnostics) with example of config.
The MaCh3 core plotting library code can be found [here](https://github.com/mach3-software/MaCh3/tree/develop/plotting) along with example config files and some apps for making standard plots.

## How To Use
This is an example how your executable can look like using MaCh3:
```cpp
  //Manager is responsible for reading from config
  std::unique_ptr<manager> fitMan = MaCh3ManagerFactory(argc, argv);

  std::vector<SampleHandlerBase*> sample; //vector storing information about sample for different detector
  std::vector<ParameterHandlerBase*> Cov; // vector with systematic implementation
  MakeMaCh3Instance(fitMan.get(), sample, Cov); //Factory like function which initialises everything

  // FitterBase class, can be replaced with other fitting method
  std::unique_ptr<FitterBase> MarkovChain = MaCh3FitterFactory(FitManager.get());

  //Adding samples and covariances to the Fitter class could be in the factory
  for(unsigned int i = 0; sample.size(); i++)
    MarkovChain->AddSampleHandler(sample[i]);
  for(unsigned int i = 0; Cov.size(); i++)
    MarkovChain->AddSystObj(Cov[i]);

  MarkovChain->RunLLHScan(); // can run LLH scan
  MarkovChain->RunMCMC(); //or run actual fit
```
For more see [here](https://github.com/mach3-software/MaCh3Tutorial/blob/main/Tutorial/MCMCTutorial.cpp)
