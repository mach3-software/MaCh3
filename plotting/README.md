## Plotting

## Config Files

The goal with MaCh3s plotting library is to be as flexible as possible and to abstract away all of the most annoying parts about plot making like reading in the data from fitter output files and keeping track of parameter names across different fitters, while allowing the user as much freedom as possible when it comes plot style and formatting . In order to do this, MaCh3s plotting library makes heavy use of config files so that the behaviour of the code can be adapted easily to individual experiments without having to totally reinvent the wheel and rewrite large ammounts of boilerplate code. In line with the rest of MaCh3, the plotting library uses [YAML](https://yaml.org/spec/1.2.2/) to specify its config files. There are a few different config files which control different parts of the plotting library:

### Plotting Config

This is the highest level config file used to control high level variables like where to find other config files and various options for the plotting [applications](#apps). This config acts as a sort of fixed reference point for the plotting library and so it's location is semi-hardcoded using the $MACH3 environment variable. In your experiment specific MaCh3 repository which is built against this MaCh3 core repository, you should define a plotting directory which contains the config file "PlottingConfig.yaml":

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

Additionally, this config file contains options specific to each of the (applications)[#apps] described below. They are described in more detail in that section but the general parrern they follow is 
```
{ApplicationName}:
    {Option1}: <value>
    {Option2}: <value>
    ...
```

### Translation Config

### Style Config

## Fitter Specification

## Apps

**GetPostfitParamPlots** - This will plot output from ProcessMCMC for nice plots which can go to TN. Bits are hardcoded to make plots nicer users should be careful when using the non-conventional xsec model. If you used `ProcessMCMC` with `PlotDet` you will also get an overlay of detector parameters (ND or ND+FD depending on chain type). If Violin plot was produced in `ProcessMCMC` you will get fancy plots here as well.


**PlotLLH** - Plot LLH scans, flexible and configurable in command line. can take any number of LLH scans as input, will use the first one as a baseline when making e.g. ratio plots. The first file must be a MaCh3 scan.
options:

    -r overlay ratio plots

    -s also make plots split by sample contribution, to use this option, the LLH scan must have been run with the option `LLH_SCAN_BY_SAMPLE = true` in the config file

    -g draw a grid on the plots

    -l a string specifying the labels for each scan to plot in the legent. this should be a string of labels separated by semi colons, e.g.: -`l "label1;label2;label3"`

    -o the name of the output pdf

    -d a string specifying additional drawing options to pass to the histogram draw calls, e.g. `-d "C"` will plot smooth curves through the histogram bins. See https://root.cern/doc/master/classTHistPainter.html#HP01a for possible options.
