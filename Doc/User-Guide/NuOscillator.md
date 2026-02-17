# NuOscillator {#nuoscillator}

MaCh3 is using [NuOscillator framework](https://github.com/dbarrow257/NuOscillator) to perform oscillation probabilities calculations. 
Almost everything is being handled via NuOscillator config example can be found [here](https://github.com/dbarrow257/NuOscillator/tree/main/NuOscillatorConfigs). 

## Engines 
Multiple engines are available, ranging from beam, and atmospheric to BSM.
To see available engines we recommend [this](https://github.com/dbarrow257/NuOscillator/tree/main?tab=readme-ov-file#implemented-engines)

## MaCh3 Specyfic Settings
NuOscConfigFile contains information like binned vs. unbinned, engine, etc., which is handled by NuOscillator. EqualBinningPerOscChannel is a setting used for binned/subsampling oscillations only. This setting tells whether each OSC channel should be using its own NuOscillator object and thus different binning or all OSC channels should be using the same binning and single object (making fit overall much faster).

```
NuOsc:
  NuOscConfigFile: "TutorialConfigs/NuOscillator/NuFASTLinear.yaml"
  EqualBinningPerOscChannel: true
```


