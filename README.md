# MaCh3 <img src="Doc/mach3logo.png" alt="MaCh3" align="center" width="100"/>
Markov Chain 3 flavour is frameworks which was born in 2013 as [T2K](https://t2k-experiment.org/pl/) Bayesian MCMC fitter for oscillation analysis. It has been used for multiple T2K Oscillation analysis both at Near and Far detectors throughout years.

TODO: Add more history: maybe stuff about T2K+SK and T2K+NOvA and mention HK and DUNE

It has been used to  Since then framework evolved and has non MCMC modules.

## Cite
When citing MaCh3, please use [on Zenodo](https://zenodo.org/records/7608419#.Y-BgaC8RrpA).

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
MaCh3 supports quite a high range of CUDA architectures if something doesn't work on your GPU let us know.


## Oscillator
MaCh3 uses several neutrino oscillation calculators. By default, CUDAProb3 is used. If you would like to use Prob3++

```
cmake ../ [-DUSE_PROB3=<ON,OFF>]
```
Following neutrino oscillation calculators are available:
<ol>
<li> CUDAProb3 [**CPU/GPU**][**Beam/Atm**]  </li>
<li> Prob3++ [**CPU**][**Beam**] </li>
<li> probGPU [**GPU**][**Beam**] </li>
</ol>

## Fitting algorithms
The following fitting algorithms are available:
<ol>
<li> MR2T2  </li>
<li> MINUIT2  </li>
<li> PSO  </li>
</ol>

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
Most of external libraries are being handled through CPM. The only external library that is not being handled through CPM and is required is [ROOT](https://root.cern/).
```
  GCC: ...
  CMake: ...
  ROOT: ...
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

## Plotting and Diagnostic
Example of chain diagnostic utils can be found [here](https://github.com/mach3-software/MaCh3/tree/develop/Diagnostics).

<p align="left">
  <img src="Doc/delta.png" alt="delta" width="200"/>
</p>

TODO this should be expanded


## Help and Guidelines
- [How to contribute](https://github.com/mach3-software/MaCh3/blob/develop/CONTRIBUTING.md)
- [Wiki](https://github.com/mach3-software/MaCh3/wiki)
- [Mailing lists](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=MACH3)


