## Plotting

## Config Files

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
