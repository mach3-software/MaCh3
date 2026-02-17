# Manager and Config Handling {#manager-and-config-handling}

Manager class is meant to handle or manage configs. MaCh3 uses yaml configs. There are several necessary entires your config should have

## General
Here you specify name of the output file and whether you want to run data fit, Asimov fit or fake data fit.
```
General:
  OutputFile: "blarb.root"
  Asimov: True
  RealData: False
  FakeData: False
```
## Systematics
Here you load yaml configs describing each systematic type. What is important you load vector configs which are later merged into one Covariance Class which describes cross-section, flux etc.
```
  Systematics:
    XsecCovFile: ["inputs/SystematicsOA2024XsecSplineParams_noParName.yaml", "inputs/SystematicsOA2024XsecNormParams_noParName.yaml", "inputs/SystematicsOA2024XsecFunctionalParams_noParName.yaml", "inputs/SystematicsOA2024FluxParams_noParName.yaml"]
    XsecCovName: "xsec_cov"
    OscCovFile: "inputs/osc_covariance_2021_PDG2021_v1.root"
    OscCovName: "osc_cov"
```

## MCMC
Number of steps for MCMC algorithm, and step scale. To learn more about step scale please visit chapter about diagnostic and step size tunning.
```
  MCMC:
    NSteps: 2000000
    #KS: Sets how frequent (number of steps) you want to autosave your chain, if you autosave rarely chain will be sliglhy faster, if you job wall is small you might want more freqeunt autosave to reduce number of lost steps
    #AutoSave: 500
    AutoSave: 10000
    UseReactorConstraint: No
    #Burn in for posterior predictive code
    BurnInSteps: 200000
    XsecStepScale: 0.01
    NdDetStepScale: 0.035
```

## Likelihood type
Likelihood type which will be used in the fit.
```
LikelihoodOptions:
  #False means you calculate W2 histogram only during the first reweight, advisable for Barlow-Beesto
  UpdateW2: false
  TestStatistic: "Barlow-Beeston"
  #TestStatistic: "Poisson"
  #TestStatistic: "DembinskiAbdelmottele"
  #TestStatistic: "IceCube"
```

## Overriding
There are two ways how to easily override configs:
```
./bin/MCMCTutorial TutorialConfigs/FitterConfig.yaml General:MCMC:NSteps:100000
```
This is useful for overriding single arguments. If you are interested in overriding multiple arguments, you can override the default config with an additional configuration.

```
./bin/MCMCTutorial TutorialConfigs/FitterConfig.yaml --override TutorialConfigs/Override.yaml
```
