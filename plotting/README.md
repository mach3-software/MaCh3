## Plotting

This directory contains the core MaCh3 plotting tools. `plottingUtils` contains all the code defining the plotting library, the majority of which basically consists of a highly flexible root file reader which can be configured by different experiments to read from their MaCh3 output files and put everything in a nice common format that can then be accessed in plotting applications and allows for easy iteration over parameters and samples. It is also possible to configure the plotting library to be able to read files from other fitters, which is especially useful for cross-fitter validations (a task that tends to take up a lot of time during an analsis). The aim here is to avoid having to have vast numbers of un-maintained plotting scripts for dealing with every possible combination of fitter and analysis iteration and instead maintain this core code and allow experiments to simply maintain their own config files.

## Structure

The plotting library consists of 3 main manager classes: `PlottingManager`, `inputManager` and `styleManager` each of which is responsible for different aspects of plot making.

### PlottingManager

This is intended to be the main class that the user interacts with. It is responsible for parsing user specified options, keeping track of which config files to use, and managing instances of the other manager classes.

### inputManager

This is responsible for keeping track of the files whose data is to be plotted and reading in that data from them. It will attempt to automatically detect which fitter the file came from based on the config-defined file structures. It will read in the data accordingly, storing in a common format which it provides functions to easily obtain.

### styleManager 

This is responsible for the styling of the plots produces. It contains the functionality to apply user specified styles, colour palettes and fancy names for parameters and samples.

## Config Files

The goal with MaCh3s plotting library is to be as flexible as possible and to abstract away all of the most annoying parts about plot making like reading in the data from fitter output files and keeping track of parameter names across different fitters, while allowing the user as much freedom as possible when it comes plot style and formatting . In order to do this, MaCh3s plotting library makes heavy use of config files so that the behaviour of the code can be adapted easily to individual experiments without having to totally reinvent the wheel and rewrite large ammounts of boilerplate code. In line with the rest of MaCh3, the plotting library uses [YAML](https://yaml.org/spec/1.2.2/) to specify its config files. There are a few different config files which control different parts of the plotting library. Some example yaml configs found in this directory.

## Apps

There are a few pre-written applications to plot some standard things that we typically look at (there will be some more to follow but they need to be "modernised" to use plotting library)

**GetPostfitParamPlots** - This will plot output from ProcessMCMC for nice plots which can go to TN. Bits are hardcoded to make plots nicer users should be careful when using the non-conventional xsec model. If you used `ProcessMCMC` with `PlotDet` you will also get an overlay of detector parameters (ND or ND+FD depending on chain type). If Violin plot was produced in `ProcessMCMC` you will get fancy plots here as well.


**PlotLLH** - Plot LLH scans, flexible and configurable in command line. can take any number of LLH scans as input, will use the first one as a baseline when making e.g. ratio plots. The first file must be a MaCh3 scan.
options:

    -r overlay ratio plots

    -s also make plots split by sample contribution, to use this option, the LLH scan must have been run with the option `LLH_SCAN_BY_SAMPLE = true` in the config file

    -g draw a grid on the plots

    -l a string specifying the labels for each scan to plot in the legent. this should be a string of labels separated by semi colons, e.g.: -`l "label1;label2;label3"`

    -o the name of the output pdf

    -d a string specifying additional drawing options to pass to the histogram draw calls, e.g. `-d "C"` will plot smooth curves through the histogram bins. See https://root.cern/doc/master/classTHistPainter.html#HP01a for possible options.
