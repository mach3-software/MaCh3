<p align="center">
  <img src="Doc/mach3logo.png" alt="Mach3 Logo" width="100"/>
</p>

# MaCh3
MaCh3 is ...

# How to Compile

```
mkdir build; cd build;
cmake ../
```

Don't forget to
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

Some fucnitonalities relies on setting `Env{MACH3}` which should point to path experiment specyfic MaCh3. This way MaCh3 can easily find `Env{MACH3}/inputs/SomeInput.root` for example.

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


## Oscialtor
MaCh3 uses several neutrino oscillation calculators. By default, CUDAProb3 is used. If you would like to use Prob3++

```
cmake ../ [-DUSE_PROB3=<ON,OFF>]
```
Following neutrino oscillation calculators are available:
<ol>
<li> CUDAProb3 [CPU/GPU][Beam/Atm]  </li>
<li> Prob3++ [CPU][Beam] </li>
<li> probGPU [GPU][Beam] </li>
</ol>

## Fitting algorithms
Following fitting algorithms are available:
<ol>
<li> MR2T2  </li>
<li> MINUIT2  </li>
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
```
  GCC: ...
  CMake: ...
  ROOT: ...
```

# How To Use
This is example how your exectuable can look like using MaCh3:
```
  manager *fitMan = nullptr; //Manager is responsible for reading from config

  std::vector<samplePDFBase*> sample; //vector storing infomation about sample for different detector
  std::vector<covarianceBase*> Cov; // vector with systematic implementation
  mcmc *markovChain = nullptr; // MCMC class, can be repalced with other fitting method
  MakeMaCh3Instance(fitMan, sample, Cov, markovChain); //Factory like function whihc initialsies everything

  //Adding samles and covariances to Fitter class, could be in factory
  for(unsigned int i = 0; sample.size(); i++)
    markovChain->addSamplePDF(sample[i]);
  for(unsigned int i = 0; Cov.size(); i++)
    markovChain->addSystObj(Cov[i]);

  markovChain->RunLLHScan(); // can run LLH scan
  markovChain->runMCMC(); //or run actual fit
```


## Help and Guidelines
- [How to contribute](https://github.com/mach3-software/MaCh3/develop/CONTRIBUTING.md)
- [Wiki](https://github.com/mach3-software/MaCh3/wiki)
