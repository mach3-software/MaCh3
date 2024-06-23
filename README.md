# MaCh3 <img src="Doc/mach3logo.png" alt="MaCh3" align="center" width="100"/>
The Markov Chain 3 flavour is a framework born in 2013 as a Bayesian MCMC fitter for [T2K](https://t2k-experiment.org/pl/) oscillation analysis. It has now been used for multiple T2K Oscillation analyses both at the Near and Far detectors throughout the years and is also used by the DUNE and HK oscillation analysis groups as well as for joint fits between T2K and NOvA and T2K and SK's atmospheric data.

The framework has also evolved to allow non MCMC modules to interrogate the likelihoods implemented.

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://github.com/mach3-software/MaCh3/blob/develop/LICENSE.txt)
[![Release - v1.0.0](https://img.shields.io/badge/Release-v1.0.0-2ea44f)](https://github.com/mach3-software/MaCh3/releases)
[![Container Image](https://img.shields.io/badge/Container-Image-brightgreen)](https://github.com/mach3-software/MaCh3/pkgs/container/mach3)
[![Code - Documented](https://img.shields.io/badge/Code-Documented-2ea44f)](https://github.com/mach3-software/MaCh3/wiki)
[![Code - Doxygen](https://img.shields.io/badge/Code-Doxygen-2ea44f)](https://mach3-software.github.io/MaCh3/index.html)
[![Build Status](https://github.com/mach3-software/MaCh3/workflows/Docker%20CI%20Alma9/badge.svg)](https://github.com/mach3-software/MaCh3/actions?query=workflow%3A%22Docker+CI+Alma9%22)

## Famous Plots
Example of plots made using MaCh3 apparent in scientific publications, for more see [here](https://github.com/mach3-software/MaCh3/wiki/14.-MaCh3-in-the-Field)
<img src="Doc/Plots/delta.png" alt="MaCh3" align="left" width="200"/>
<img src="Doc/Plots/Jarlskog.png" alt="MaCh3" align="center" width="200"/>

## Cite
When citing MaCh3, please use [on Zenodo](https://zenodo.org/records/10949376) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10949376.svg)](https://doi.org/10.5281/zenodo.10949376).

# How to Compile
MaCh3 follows standard cmake pattern. By default you should get most optimal build setting although below we list many configurable options:
```
mkdir build;
cd build;
cmake ../
make -jN [set number of threads]
make install
```

Don't forget to:
```
source bin/setup.MaCh3.sh
```
## Building against MaCh3
If you compiled MaCh3 and sourced it you can simply call
```
find_package(MaCh3)
```
and cmake will find it. Alternatively, you can use CPM, for example:
```
CPMFindPackage(
  NAME MaCh3
  GIT_TAG "blarb"
  GITHUB_REPOSITORY mach3-software/MaCh3
)
```
Once you found MaCh3 you might want to link your library against MaCh3. You can do this as follows:
```
target_link_libraries(blarb MaCh3::All)
```

Some functionalities rely on setting `Env{MACH3}` which should point to path experiment specific MaCh3. This way MaCh3 can easily find `Env{MACH3}/inputs/SomeInput.root` for example.

## Multithreading
MaCh3 quite heavily relies on Multithreading, it is turned on by default. If for debugging purposes you would like to turn it off please use

```
cmake ../ [-DMaCh3_MULTITHREAD_ENABLED=<OFF>]
```

## CUDA
If the system has access to GPU, MaCh3 will enable GPU functionality automatically. If you would like to CPU only despite having access to CUDA
```
mkdir build; cd build;
cmake ../ [-USE_CPU=ON]
```
MaCh3 supports quite a high range of CUDA architectures if something doesn't work on your GPU let us know. MaCh3 supports only NVIDIA GPUs.


## Oscillator
MaCh3 uses several neutrino oscillation calculators. By default, CUDAProb3 is used. If you would like to use Prob3++

```
cmake ../ [-DUSE_PROB3=<ON,OFF>]
```
Following neutrino oscillation calculators are available:

|Oscillator  | Hardware   | Source     |
|------------|------------|------------|
| CUDAProb3  | CPU/GPU    | Beam/Atm   |
| Prob3++    | CPU        | Beam       |
| probGPU    | GPU        | Beam       |

## Fitting algorithms
The following fitting algorithms are available:

| Algorithm  | Need Ext Library |
|------------|------------------|
| MR2T2      | No               |
| MINUIT2    | Yes              |
| PSO        | No               |


## Debug
Several debugging options are available which are heavy for RAM and performance and, therefore not used by default. To enable it:
```
cmake ../ [-DMaCh3_DEBUG_ENABLED=<ON,OFF>]
```
There are several debug modes, to enable more detailed but very heavy specific debug levels. Level 1 is the default debug activated by the above.

```
cmake ../ [-DMaCh3_DEBUG_ENABLED=<ON,OFF>] [-DDEBUG_LEVEL=<1,2,3>]
```
## System Requirements
Most of external libraries are being handled through CPM. The only external library that is not being handled through CPM and is required is [ROOT](https://root.cern/). Currently used external dependencies include:

1. [yaml-cpp](https://github.com/jbeder/yaml-cpp)
2. [spdlog](https://github.com/gabime/spdlog)

Based on several test here are recommended version:
```
  GCC: ...
  CMake: >= 3.14
  ROOT: >= 6.18
```

# How To Use
This is an example how your executable can look like using MaCh3:
```
  manager *fitMan = nullptr; //Manager is responsible for reading from config

  std::vector<samplePDFBase*> sample; //vector storing information about sample for different detector
  std::vector<covarianceBase*> Cov; // vector with systematic implementation
  mcmc *markovChain = nullptr; // MCMC class, can be replaced with other fitting method
  MakeMaCh3Instance(fitMan, sample, Cov, markovChain); //Factory like function which initialises everything

  //Adding samples and covariances to the Fitter class could be in the factory
  for(unsigned int i = 0; sample.size(); i++)
    markovChain->addSamplePDF(sample[i]);
  for(unsigned int i = 0; Cov.size(); i++)
    markovChain->addSystObj(Cov[i]);

  markovChain->RunLLHScan(); // can run LLH scan
  markovChain->runMCMC(); //or run actual fit
```

## Help and Guidelines
- [How to contribute](https://github.com/mach3-software/MaCh3/blob/develop/CONTRIBUTING.md)
- [Wiki](https://github.com/mach3-software/MaCh3/wiki)
- [Mailing lists](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=MACH3)
- [Slack](https://t2k-experiment.slack.com/archives/C06EM0C6D7W/p1705599931356889)


### Plotting and Diagnostic
Example of chain diagnostic utils can be found [here](https://github.com/mach3-software/MaCh3/tree/develop/Diagnostics) with example of config.
