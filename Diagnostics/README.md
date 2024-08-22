### Plotting and Diagnostic
**ProcessMCMC** - The main app you want to use for analysing the ND280 chain. It prints posterior distribution after burn-in the cut. Moreover, you can compare two/three different chains. There are a few options you can modify easily inside the app like selection, burn-in cut, and whether to plot xse+flux or only flux. Other functionality
<ol>
<li> Produce a covariance matrix with multithreading (which requires lots of RAM due to caching) </li>
<li> Violin plots </li>
<li> Credible intervals and regions </li>
<li> Calculate Bayes factor and give significance based on Jeffreys scale </li>
<li> Produce triangle plots </li>
<li> Study covariance matrix stability </li>
</ol>

**GetPostfitParamPlots** - This will plot output from ProcessMCMC for nice plots which can go to TN. Bits are hardcoded to make plots nicer users should be careful when using the non-conventional xsec model. If you used `ProcessMCMC` with `PlotDet` you will also get an overlay of detector parameters (ND or ND+FD depending on chain type). If Violin plot was produced in `ProcessMCMC` you will get fancy plots here as well.

**GetPenaltyTerm** - Since xsec and flux and ND spline systematic are treated as the same systematic object we cannot just take log_xsec, hence we need this script, use `GetFluxPenaltyTerm MCMCChain.root config`. Parameters of relevance are loaded via config, thus you can study any combination you want. Time needed increases with number of sets :(

**DiagMCMC** - Perform MCMC diagnostic like autocorrelation or trace plots.

**RHat** - Performs RHat diagnostic to study if all used chains converged to the same stationary distribution.
```
./RHat Ntoys MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]
```

**PlotLLH** - Plot LLH scans, flexible and configurable in command line. can take any number of LLH scans as input, will use the first one as a baseline when making e.g. ratio plots. The first file must be a MaCh3 scan.
options:

    -r overlay ratio plots

    -s also make plots split by sample contribution, to use this option, the LLH scan must have been run with the option `LLH_SCAN_BY_SAMPLE = true` in the config file

    -g draw a grid on the plots

    -l a string specifying the labels for each scan to plot in the legent. this should be a string of labels separated by semi colons, e.g.: -`l "label1;label2;label3"`

    -o the name of the output pdf

    -d a string specifying additional drawing options to pass to the histogram draw calls, e.g. `-d "C"` will plot smooth curves through the histogram bins. See https://root.cern/doc/master/classTHistPainter.html#HP01a for possible options.


**CombineMaCh3Chains** - will combine chains files produced by **MCMC**, enforcing the condition that all the files to combine were made using the exact same software versions and config files
```
CombineMaCh3Chains [-h] [-c [0-9]] [-f] [-o <output file>] file1.root [file2.root, file3.root ...]
```
*fileX.root* are the individual spline files to combine, can specify any number, need at least one

-c target compression level for the combined file, default is 1, in line with hadd

-f force overwrite of the combined file if it exists already

-h print usage message and exit

*Output file* (optional) name of the combined file. If not specified, will just use *file1.root*, the first in the list of files, same as *hadd*.
