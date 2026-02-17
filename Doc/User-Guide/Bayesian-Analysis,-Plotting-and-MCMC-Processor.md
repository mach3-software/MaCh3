# Bayesian Analysis {#bayesian-analysis}
MCMC Processor is a class responsible for processing MCMC and producing validation plots.

## Posteriors
Marginalised posterior is the standard output of MCMC. PDF is mean, Gauss indicates Gaussian fit to posterior while HPD (Highest Posterior Density Point). In most cases parameter posteriors are Gaussian, then all 3 values would give the same result. However, the strength of MCMC is that there is no assumption about gaussiantiy, and it can handle non-Gaussian parameters.

<img width="350" src="https://github.com/user-attachments/assets/e77f7287-c71e-49c6-9f61-161fd54e4229">

This plot will all be produced when running ProcessMCMC.

## Parameter Plot
This plot is summarising 1D posteriors in more compact way. It is being printed by default.

<img width="350" src="https://github.com/user-attachments/assets/37c15953-1926-44a4-9a87-d232ba2aaf06">

Since these plots are often targeted for publication, there is a tool called **GetPosftitParams** which make this plot with ability group parameter in whatever fashion you like. You can read more [here](https://github.com/mach3-software/MaCh3/wiki/15.-Plotting)

## Credible Intervals
This plot helps to tell which values are excluded based on X-credible intervals

<img width="350" src="https://github.com/user-attachments/assets/2199f9ff-1901-4f7f-a5e8-43433049027e">

To produce, make sure the following field is true.
```
ProcessMCMC:
  MakeCredibleIntervals: true
```

# 2D Posteriors
It is possible to produce 2D posteriors. They are very useful to identify if parameters are correlated or not. In this example, there are strong correlations.

<img width="350" src="https://github.com/user-attachments/assets/91bb1c28-e3a5-4418-a077-ad9b968a29b5">

This can take some time, though. There are two ways: faster (using multithreading) but requiring lots of RAM, or slower but without RAM requirements. Once you obtain 2D posteriors, you can produce multiple additional plots.

To enable, this option must be on.
```
ProcessMCMC:
  PlotCorr: true
```

### Threshold
Not every 2D plot will be made. MaCh3 uses configurable threshold to print only more interesting plot. Threshold applies to correlation factor.
```
ProcessMCMC:
  Post2DPlotThreshold: 0.2
```

## Credible Region

<img width="350" src="https://github.com/user-attachments/assets/fccae6a7-d18b-455f-9ed7-a413cc6eb746">

```
ProcessMCMC:
  PlotCorr: true
  MakeCredibleRegions: true
```

## Correlation Matrix
Correlation matrix etc are calculated based on 2D Posteriors

<img width="350" src="https://github.com/user-attachments/assets/71d4eacd-10f8-4738-b32c-4de6a6d6a3a3">

<img width="350" src="https://github.com/user-attachments/assets/a77f7e24-e381-46f5-b924-f6e034676b74">

## Triangle plot
```
  TrianglePlot:
    - ["Test", ["Norm_Param_0", "Norm_Param_1", "Spline_Param_0", "Spline_Param_1"]]
```
You can specify as many parameters as you like. But also as many combinations as you like
```
  TrianglePlot:
    - ["Test", ["Norm_Param_0", "Norm_Param_1", "Spline_Param_0", "Spline_Param_1"]]
    - ["Test2", ["Norm_Param_0", "Norm_Param_1"]]
```

<img width="350" src="https://github.com/user-attachments/assets/06322f8e-8a4f-4ce4-81bb-e8865e1e47e9">

## Violin plot

<img width="350" src="https://github.com/user-attachments/assets/168c224d-11bc-45c9-a3a6-416d878cfbe9">

```
ProcessMCMC:
  PlotCorr: true
  MakeViolin: true
```

## Bayes factor and Savage-Dickey
It is possible to obtain the Bayes factor for different hypothesis
```
  BayesFactor:
    # Goes as follows: ParamName Name[Model 1, Model 2], Model1[lower, upper ], Model2[lower, upper ]
    - ["sin2th_23", ["LO", "UO"], [0, 0.5], [0.5, 1]]
```

<img width="350" src="https://github.com/user-attachments/assets/258d88c5-5dc0-4e69-a1ce-282918570737">

or calculate savage Dickey, which is Bayes factor for point-like hypothesis
```
  SavageDickey:
    - ["Alpha_q3", 0.0001, [0, 1]]
```
<img width="350" src="https://github.com/user-attachments/assets/c021450b-7634-4e18-b9d6-fd2695a6a3e4">

## Parameter Evolution
```
  ParameterEvolution:
    - ["Norm_Param_0", 20]
```
Select parameter name and how many frames you want. The more, the longer it takes, so be careful

<img width="350" src="https://github.com/user-attachments/assets/84fce0eb-450f-437b-904f-9123300c3cf8">

## Bipolar plot

<img width="350" src="https://github.com/user-attachments/assets/a1c9b296-e05e-4d3e-8923-f39e2704f759">

## Ridgeline plot

<img width="350" src="https://github.com/user-attachments/assets/cc0a15f5-2a46-46f7-9f78-5d115d56d9c9" />


## Reweighting
TODO!!!
