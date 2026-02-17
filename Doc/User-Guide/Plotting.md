# Plotting {#plotting}

The goal with MaCh3s plotting library is to be as flexible as possible and to abstract away all of the most annoying parts about plot making like reading in the data from fitter output files and keeping track of parameter names across different fitters, while allowing the user as much freedom as possible when it comes plot style and formatting. 

The plotting Library consists primarily of 3 main manager classes: 
- [PlottingManager](https://mach3-software.github.io/MaCh3/classMaCh3Plotting_1_1PlottingManager.html) - This controls the high level stuff. It deals with application level options, command line options, and holds references to the other managers. This is the main class you should interact with on the user level.
- [InputManager](https://mach3-software.github.io/MaCh3/classMaCh3Plotting_1_1InputManager.html) - This deals with the details of how input files (the outputs of the fitters) should be read and stored. This should be accessed via `<Plotting Manager Instance>.input()`.
- [StyleManager](https://mach3-software.github.io/MaCh3/classMaCh3Plotting_1_1StyleManager.html) - This provides some helpful utilities for creating plots. This should be accessed via `<Plotting Manager Instance>.style()`.

# Config Files

In order to achieve a high level of flexibility, MaCh3s plotting library makes heavy use of config files so that the behavior of the code can be adapted easily to individual experiments without having to totally reinvent the wheel and rewrite large amounts of boilerplate code. In line with the rest of MaCh3, the plotting library uses [YAML](https://yaml.org/spec/1.2.2/) to specify its config files. Each of the manager classes above has a corresponding config file to control its functionality. The format of these config files are detailed below.

## Plotting Config

This is the highest level config file and corresponds to the behavior of the [Plotting Manager](#plottingmanager). It should be used to control high level variables like where to find other config files and various options for the plotting [applications](#standard-apps). This config acts as a sort of fixed reference point for the plotting library and so it's location is semi-hardcoded using the $MACH3 environment variable. In your experiment specific MaCh3 repository which is built against this MaCh3 core repository, you should define a plotting directory which contains the config file "PlottingConfig.yaml":

```
{MACH3}/plotting/PlottingConfig.yaml
```

You can find an example of such a config in this repository [here](https://github.com/mach3-software/MaCh3/plotting/PlottingConfig.yaml)

The most important options are the ones under `ManagerOptions`:

```
ManagerOptions:
    translationConfig: "" 
    styleConfig: ""
```

 These tell the plotting manager where to find the other config files that are needed. `translationConfig` tells the manager where to find the [translation config](#translation-config) and `styleConfig` where to find the [style config](#style-config). These should be paths relative to the `MACH3` environment variable e.g. if your MaCh3 directory contained your configs in a subdirectory called `Cool-Configs` then it would look something like:

```
${MaCh3}/..
${MaCh3}/Cool-Configs/Cool-Translation-Config.yaml
${MaCh3}/Cool-Configs/Cool-Style-Config.yaml
${MaCh3}/..
```
and your config options would look like:
 
```
ManagerOptions:
    translationConfig: "Cool-Configs/Cool-Translation-Config.yaml"
    styleConfig: "Cool-Configs/Cool-Style-Config.yaml"
```

If the options are left blank as they are above then when running the plotting code they will default to

```
translationConfig: "${MACH3}/plotting/universalTranslator.yaml"
styleConfig: "${MACH3}/plotting/StyleConfig.yaml"
```

Additionally, this config file contains options specific to each of the [applications](#apps) described below. They are described in more detail in that section but the general pattern they follow is 
```
{ApplicationName}:
    {Option1}: <value>
    {Option2}: <value>
    ...
```


## Translation Config

This config defines the behaviour of the InputManager. It is used to define the output format of files from different fitters, and the parameters and samples that exist within the analysis, and what they are called in different fitters. This means that you can easily plot and compare files from different fitters without having to worry about the nitty gritty details of how to read from each specific file type and can focus on making beautiful plots! The components which are defined in this config are described below and you can find an example of such a config [here](https://github.com/mach3-software/MaCh3/plotting/UniversalTranslator.yaml).

### Fitter Specification

In this section of the config file we specify what the output of each fitter we are interested in looks like so that the plotting code knows where to look for different outputs. We start this section in typical yaml style with

```
FitterSpec:
```

Now we declare the names of fitters the plotting code should be aware of

```
    fitters: ["FITTER_1", "FITTER_2"]
```

Now for each fitter specified we give some details. For each fitter declared in `fitters: [blabla]` you need to now include a yaml node like

```
    FITTER_1:
        ....

    FITTER_2:
        ...
```

Under each of these headers you must now give details. 

#### LLH Scan 
Let's start with LLH scan information. 

You must specify the `LLHObjectType`. This is the type of root object that the plotter should expect to find in the output file. i.e. how this fitter stores their LLH scans. This can at present be either TH1D or TGraph. For MaCh3, TH1D is typically used so we would have
```
        LLHObjectType: "TH1D"
```

Now we specify where to look for the scans for each type of LLH scan (sample, penalty, and total). Each of these gets its own header and the locations are specified using the location strings described in more detain [here](#location-strings). This may look something like 

```
        ## tell the code where to find likelihood scans for sample, penalty and total likelihood
        sample_LLH:
            location: ["likelihood_scans/Sample_LLH/{PARAMETER}_sam"]
            
        penalty_LLH:
            location: ["likelihood_scans/Penalty_LLH:{PARAMETER}_pen"]
            
        total_LLH:
            location: 
            [
                "likelihood_scans/Total_LLH_DIR_1:{PARAMETER}_tot",
                "likelihood_scans/Total_LLH_DIR_2:{PARAMETER}_tot",
            ]
```

We often also like to make likelihood scans broken down by sample. Here we can also tell the code where to find these per-sample likelihood scans. This might look something like
```
        bySample_LLH:
            location: 
            [
                "{SAMPLE}_LLH:{PARAMETER}{SAMPLE}",
            ]
```

#### MCMC

We can now also specify where to look for things relating to the MCMC posteriors.

under the `1dPosteriors` heading we can specify the location to look for one dimensional posterior objects that have been produced from a raw chain using the MCMCProcessor. For example

```
        1dPosteriors:
            location: [
                "Post_TH1Ds:{PARAMETER}",
            ]
```

We can also specify where to look for the TTree containing the raw MCMC steps e.g.

```
        MCMCsteps:
            location: [
                "posteriors",
            ]
```

Note that the way that these are found is different to most of the other objects. This will use the MCMCProcessor to find the names of the branches in the posterior tree. You as a user do not need to worry about this detail but you should be aware that this use of the MCMCProcessor means that the `MCMCsteps` option is only usable for MaCh3 based fitters and you will need to have access to the yaml configs that were used initially to produce the chain.

### Parameter Specification

Another important part of the translation config is to define the parameters that are modelled in the experiment. This section of the yaml is imaginatively marked with the `Parameters` header:

```
Parameters:
```

The first step is to define the master list of all of the parameters that the plotting library should be aware of, this is defined using the `Parameters` variable. An example of this could look something like

```
    Parameters:
    [  
        "Norm_Param_0",
        "Norm_Param_1",
        "Norm_Param_2",

        "Spline_Param_0",
        "Spline_Param_1",

        "Func_Param_0",
        "Func_Param_1",
        "Func_Param_2",
        "Func_Param_3",

        "sin2th_12",
        "sin2th_23",
        "sin2th_13",

        "delm2_12",
        "delm2_23",

        "delta_cp"
    ]
```

Note that the names specified here are internal to the plotting library and do not need to correspond to the ones used in other parts of MaCh3 (although if they do it makes things a bit easier). These are simply labels used to uniquely identify the parameters across all files and fitters, as well as within the config files.

#### Parameter Specific Options

We can then also specify options specific to each parameter. The options exist under headings which should match the labels specified in the list above e.g.

```
    Norm_Param_0:

        ## We can specify "tags" to apply to a parameter
        ## This allows you to later on in your plotting scripts get, for example, all parameters which have the tag "pion"
        ## and plot them together, or all parameters with both the tags "pion" and "FSI"
        tags: ["xsec", "FSI", "pion"]

        ## We can also specify options for this parameter that are specific to each fitter.
        FITTER_2:
            ## In particular we can specify specific names for this parameter in the different fitters
            ## for different types of object. 
            LLH: "XSEC-PARAMETER-1" 
            PostFit: "PARAMETER-1"

            ## We can also provide a post fit error location for this parameter which will 
            ## override the more general locations specified earlier for this fitter. 
            ## This can be useful if e.g. the parameter appears in more than one of the possible post fit error
            ## TH1Ds. The code doesn't like this and so the override will break the degeneracy.
            postFitLoc: ["errors/Cross-Section Systematics Group 2/values/postFitErrors"]

```

If no fitter specific names are specified for a parameter then the name will default to the label in the list above. Hence why it is useful to use the regular MaCh3 names for those.

### Sample Specification

We must also do the same for the samples that exist in the fitter. The format for doing this is essentially the same as for the parameters as described above. An example of this might look like

```
Samples:

    ## Here we define the master list of all the available samples
    Samples: 
    [
        "SAMPLE_1",
        "SAMPLE_2",
        "SAMPLE_3",
    ]

    ## now we can specify particular overrides for each parameter for different fitters
    ## For example, lets assume that FITTER_1 uses the same names that are laid out in the list above
    ## but FITTER_2 uses a slightly different naming convention
    SAMPLE_1: 
        FITTER_1:

        FITTER_2:
            ## this is the name of the sample when dealing with LLH scans
            ## currently this is the only one that matters, however in future it might be useful
            ## to also be able to specify this for different file types, hence why LLH is specified
            ## and not just something general
            LLH: "Sample-1"

        ## We can specify "tags" to apply to a sample
        ## This allows you to later on in your plotting scripts get, for example, all the samples which, for example,
        ## are related to some particular sub-detector by adding a tag for that detector
        ## Or to samples which correspond to a particular tagged particle
        tags: ["subDet1", "particleTag"]
```

### Location Strings

The translation config makes use of a custom location specifier format which we will describe here. 

Generally speaking multiple locations can be specified for each object type by using a list in the config file like

```
location: [
    "locationString1",
    "locationString2",
    "locationString3",
    ...
```

Each specified location will be checked when looking for objects, which can be very useful e.g. if a fitter saves different parameter types in different locations.

There are two types of location strings that can be specified. The first is the simplest case where you simply specify the exact location of an object with something like

```
location: ["path\to\objects\{PARAMETER}"]
```

(See [Special Tokens](#special-tokens) below on how to use tokens like `{PARAMETER}`)
This will look for objects at that specific location with no flexibility. This has the advantage of being slightly faster as only one location needs to be checked per object. But sometimes this is not flexible enough and a slightly broader search is needed. This brings us to the second option in which the directory and object names can be specified separately by breaking up the string using `:` as a delimiter like

```
location: ["path\to\objects:{PARAMETER}"]
```

In this case the TDirectory `path\to\objects` will be loaded. Then the parameter {PARAMETER} will be searched for inside this directory by trying to match the expanded {PARAMETER} to the end of the object names in the directory. This can be useful if for example some unknown string like a parameter ID gets prepended to the names of objects when saving, so that the TDirectory looks something like:

```
path\to\objects
    flux1_parameter_1
    flux2_parameter_2
    xsec1_parameter_3
    xsec2_parameter_4
```

Having to keep track of the fluxA, xsecB labels and updating the fitter specific parameter names manually would be incredibly tedious and not very sustainable as they are liable to change as more parameters of different types are added. Instead by specifying the location `path\to\objects:{PARAMETER}`, these parameters would be found correctly without having to worry about these additional labels.

(NB: This approach will very likely change to become more flexible in the future. As it stands it is only able to deal with cases where the end of the parameter name matches, this isn't very flexible. It would likely be a lot better to use regex or something like that to be able to specify more general cases with wildcards. regex is likely a good approach since root has built in functionality for matching these to object names)


#### Special Tokens

There are a number of special tokens that can be used when specifying a location. These tokens, which take the form `{SOME_LABEL}` will be replaced in the code when searching for a particular object in the input file. Currently the following tokens exist:

- `{PARAMETER}` - This will be replaced with the fitter specific name of the parameter currently being looked for. e.g. if you specify the location string "path_to\{PARAMETER}_LLH" for an LLH scan then when the code is looking for LLH scans for parameter_1 the location will be expanded to "path_to\parameter_1_LLH"

- `{SAMPLE}` - This will be replaced by the fitter specific name of the sample currently being looked for. e.g. if you specify the location string "path_to\{SAMPLE}\{PARAMETER}_LLH" for the sample specific LLH scans then when looking for the scan of parameter_2 for sample sample_1 the location will be expanded to "path_to\sample_1\parameter_2_LLH"

- `{PARAMETER2}` - This can be used for cases where an object is labelled by 2 parameters, e.g. 2D LLH scans, or 2D posterior distributions. If you specified the location string "path_to\{PARAMETER}_{PARAMETER2}_2D_object" then when looking for a 2D object for parameter_1 and parameter_45 then it would be expanded to "path_to\parameter_1_parameter_45_2D_object". (Note that currently no 2d objects are actually looked for but this could be useful in the future) 


## Style Config

This config corresponds to the behavior of the [Style Manager](#stylemanager). Here you can specify options relating to the style of your plots. 

### Colour Palettes

You can specify colour palettes in the same style as root as follows:

```
ColorPallettes: # <- Tells yaml that the following are colour palettes

  RedWhiteBlue: # <- The name of this palette
  [
      ## blue to white to red gradient, usually used for cov matrices and sigma vatiations
      [255.0], ## Number of colours
      [ 0.00, 0.25, 0.50, 0.75, 1.00 ], ## stops
      [ 0.00, 0.25, 1.00, 1.00, 0.50 ], ## Reds
      [ 0.00, 0.25, 1.00, 0.25, 0.00 ], ## greens
      [ 0.50, 1.00, 1.00, 0.25, 0.00 ], ## blues
  ]

  AnotherPalette:
  ...
```
Which will give you a red-white-blue palette. The palettes you specify here can then be used in your plotting scripts by referencing the name set here using, for example, `<plotting manager instance>.style()->SetPalette("RedWhiteBlue")`

### TH1 Styles

You can define styles that can be applied to TH1s as follows:

```
TH1Styles:
## define custom styles for TH1 type hists

  redHatchedError:     ## <- name of the style
    MarkerColor: 632   ## marker colour 
    MarkerStyle: 7     ## marker style
    FillColor: 632     ## fill colour
    FillStyle: 3003    ## fill style
    LineColor: 632     ## line Colour
    LineStyle: 1       ## line style
```
Which can then be applied to a TH1 using `<Plotting Manager Instance>.style()->setTH1Style(/*TH1*/ <histogram>, /*std::string*/ <Style Name>)` which will apply the style defined under <Style Name> to the TH1 <histogram>.

### Pretty Names

You can define fancy names to use in place of the parameter and sample tags defined in the [translation config](#translation-config). You can even use latex here for added fanciness. This should look something like:

```
    ## first nice names for parameters
    parameters:

        XSEC_PAR_1: "{#phi}ancy xsec parameter: 1"
        XSEC_PAR_2: "{#phi}ancy xsec parameter: 2"
        XSEC_PAR_3: "{#phi}ancy xsec parameter: 3"

        DETSYST_PAR_1: "{#Delta}etector systematic 1"
        DETSYST_PAR_2: "{#Delta}etector systematic 2"
        
        FLUX_PAR_1: "{#Phi}_1"
        FLUX_PAR_2: "{#Phi}_2"
        
        OSC_PAR_1: "{#Delta}m^2_{23}"
        OSC_PAR_2: "{#delta}_{CP}"
        OSC_PAR_3: "{#theta}_{14}"

    ## now same for samples
    samples:
        
        SAMPLE_1: "sample 1"
        SAMPLE_2: "sample 2"
        SAMPLE_3: "sample 3"

```

# Usage

In this section we will cover the general usage of the plotting library. The [Standard Apps](#standard-apps) sections covers the usage of the apps that are used to make the "standard" MaCh3 plots that often appear in tech notes and papers. If you are not interested in doing anything fancy and simply want to reproduce some of these standard plots for comparisons and validations against a previous result, then you can just use these. 

If you are interested in doing more advanced things and need a more custom solution, then see the [Custom Plotting Scripts](#custom-plotting-scripts) section for information on using the plotting utilities in your own scripts.

## Command Line Interface

Through the use of the PlottingManager class, we are able to define a common command line interface across all plotting apps (which you can also use in your own custom plotting apps, and even python scripts!). The general usage pattern looks like 

```
<Plotting App> [Optional parameters] <Input File 1> <Input File 2> ...
```

where the optional parameters are detailed in [Available Options](#available-options).

In the case of a python script this simply becomes

```
python <Plotting Script>.py [Optional parameters] <Input File 1> <Input File 2>
```

### Available Options

There are a number of pre-defined general options which can be specified when plotting. These are described below:

```
    -r overlay ratio plots

    -s also make plots split by sample contribution, to use this option, the LLH scan must have been run with the option `LLH_SCAN_BY_SAMPLE = true` in the config file

    -g draw a grid on the plots

    -l a string specifying the labels for each scan to plot in the legent. this should be a string of labels separated by semi colons, e.g.: -`l "label1;label2;label3"`

    -o the name of the output pdf

    -d a string specifying additional drawing options to pass to the histogram draw calls, e.g. `-d "C"` will plot smooth curves through the histogram bins. See https://root.cern/doc/master/classTHistPainter.html#HP01a for possible options.
```

Note that not all of these options are valid for all plotting purposes. For example, the `-s` option only makes sense for plotting objects which can be split by sample, this does not apply for example to plotting post fit errors. Currently specifying an invalid option will simply not have any effect on the produced plots however in future it would likely be useful to produce some sort of error or warning.

## Standard Apps

### PlotLLH

Plots log likelihood scans which have been made by scanning each parameter through its prior range and calculating the likelihood. can take any number of LLH scans as input, will use the first one as a baseline when making e.g. ratio plots.


### GetPostFitParamPlots

This will plot output from ProcessMCMC for nice plots which can go to TN. Bits are hardcoded to make plots nicer users should be careful when using the non-conventional xsec model. If you used `ProcessMCMC` with `PlotDet` you will also get an overlay of detector parameters (ND or ND+FD depending on chain type). If Violin plot was produced in `ProcessMCMC` you will get fancy plots here as well.

## Custom Plotting Scripts

If you want more flexibility with your plotting, you can use the manager classes in your own plotting scripts so that you have all the benefits of MaCh3s plotting library (Generalised input reading, unified command line interface, configurability, a number of general plotting utility functions) but with total freedom of how to actually make plots. In general this is as easy as including the [PlottingManager](https://mach3-software.github.io/MaCh3/classMaCh3Plotting_1_1PlottingManager.html) in your script, initialising it with the command line inputs, and then you're good to go! This section will give some more details on how to do this, along with some examples using fitter outputs that you can generate by following the [MaCh3 Tutorial](https://github.com/mach3-software/MaCh3Tutorial). The tutorial also covers usage of the plotting library in a more hands on way, so it is highly recommended to follow that in addition to reading the wiki.

### c++

### Python
