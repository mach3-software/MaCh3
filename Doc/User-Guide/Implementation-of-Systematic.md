# Implementation of Systematics {#implementation-of-systematics}

Most of the implementation is being handled via config which is being passed to the covariance object. It contains information like the name of the parameter, prior value prior error and type of systematic. 

# Special Proposals

## CircularBounds
```
- Systematic:
    Names:
      FancyName: delta_cp
    SpecialProposal:
        CircularBounds: [-3.141592, 3.141592]
```
## Flipping
```
- Systematic:
    Names:
      FancyName: sin2th_23

    # flip octant around point of maximal disappearance (0.5112)
    # this ensures we move to a parameter value which has the same oscillation probability
    SpecialProposal:
        FlipParameter: 0.5112
```
# Types
There are four types of systematic within MaCh3 framework:
<ol>
<li> Splines  </li>
<li> Normalisation </li>
<li> Function </li>
<li> Oscillation </li>
</ol>


## Splines  
Splines are being handled via a dedicated spline class. The most important thing is to ensure the correct name of the spline which will be loaded via input file.
```
- Systematic:
    SampleNames: ["ND"]
    Error: 0.06
    FlatPrior: false
    IsFixed: false
    Mode: []
    Names:
      FancyName: MaQE
    ParameterBounds: [0.0, 9999]
    ParameterGroup: Xsec
    ParameterValues:
      Generated: 1.21
      PreFitValue: 1.03
    SplineInformation:
      Mode: [0]
      SplineName: MaCCQE
      InterpolationType: TSpline3
    StepScale:
      MCMC: 0.4
    Type: Spline
```
You can select multiple interpolation types:
![image](https://github.com/user-attachments/assets/04dc2a6e-d599-4eae-a299-4cf04c3017c0)

In addition, you can cap spline knots. If for example there was a mistake in production and some knots have weights at 1000, and you don't want to remake spines:
```
    SplineInformation:
      SplineKnotUpBound: 10
      SplineKnotLowBound: 0
```
## Normalisation 
Implementation of these relies on _KinematicCuts_. Based on which event is affected by a given norm param or not.

```
- Systematic:
    DetID: 985
    Error: 0.11
    FlatPrior: false
    IsFixed: false
    Mode: [ 0 ]
    Names:
      FancyName: Q2_norm_5
      ParameterName: Q2_norm_5
    ParameterBounds:
      - 0
      - 999.
    ParameterGroup: Xsec
    TargetNuclei: [12, 16]
    KinematicCuts:
      - TrueQ2:
        - 0.25
        - 0.5
    ParameterValues:
      Generated: 1.
      PreFitValue: 1.
    StepScale:
      MCMC: 0.2
    Type: Norm
```

### KinematicCuts
TODO!!!
  
## Functional
These need to be implemented in experiment specific SamplePDF via _CalcXsecWeightFunc_.
```
- Systematic:
    DetID: 1
    Error: 6
    FlatPrior: false
    IsFixed: false
    Names:
      FancyName: EB_dial_C_nu
      ParameterName: EB_dial_C_nu
    ParameterBounds:
      - -10.
      - 15.
    ParameterGroup: Xsec
    ParameterValues:
      Generated: 0
      PreFitValue: 2.
    Correlations:
      - EB_dial_C_nubar: 0.778
      - EB_dial_O_nu: 0.870
      - EB_dial_O_nubar: 0.653
    StepScale:
      MCMC: 0.2
    Type: Functional
```
## Oscillation
Oscillation parameters don't have any special fields compared with other types. However, order and number are crucial here. This is what is being passed to the NuOscillator framework. Therefore you must be careful to pass them properly.


# Parameter Tunes
TODO!!!