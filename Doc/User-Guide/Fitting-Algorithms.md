# Fitting Algorithms {#fitting-algorithms}
There are a few fitting algorithms available. 
## MCMC
Standard Metropolis–Rosenbluth–Rosenbluth–Teller–Teller
## Minuit
Uses ROOT-based library, for more see: https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html

# Choice
To use **MCMC** you need:
```
mcmc* Fitter = new mcmc(Manager);
```
If you are interested in **MINUIT**
```
MinuitFit* Fitter = new MinuitFit(Manager);
```
Then simply
```
  for(unsigned int i = 0; sample.size(); i++)
    Fitter ->addSamplePDF(sample[i]);
  for(unsigned int i = 0; Cov.size(); i++)
    Fitter ->addSystObj(Cov[i]);

  Fitter ->RunLLHScan(); // can run LLH scan
  Fitter ->runMCMC(); //or run actual fit
```
## Factory
There is implemented Factory method which allow to select algorithm based on config setting
```
General:
  FittingAlgorithm: "MCMC"
```
or
```
General:
  FittingAlgorithm: "PSO"
```
[MaCh3Factory.cpp](https://github.com/mach3-software/MaCh3/blob/develop/mcmc/MaCh3Factory.cpp)