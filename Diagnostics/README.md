# Diagnostic
**ProcessMCMC** - The main app you want to use for analysing the ND280 chain. It prints posterior distribution after burn-in the cut. Moreover, you can compare two/three different chains. There are a few options you can modify easily inside the app like selection, burn-in cut, and whether to plot xsec+flux or only flux. Other functionality
<ol>
<li> Produce a covariance matrix with multithreading (which requires lots of RAM due to caching) </li>
<li> Violin plots </li>
<li> Credible intervals and regions </li>
<li> Calculate Bayes factor and give significance based on Jeffreys scale </li>
<li> Produce triangle plots </li>
<li> Study covariance matrix stability </li>
</ol>

You can find more options [here](https://github.com/mach3-software/MaCh3/wiki/09.-Bayesian-Analysis,-Plotting-and-MCMC-Processor).
**GetPenaltyTerm** - Since xsec and flux and ND spline systematic are treated as the same systematic object we cannot just take log_xsec, hence we need this script, use `GetFluxPenaltyTerm MCMCChain.root config`. Parameters of relevance are loaded via config, thus you can study any combination you want. Time needed increases with number of sets :(

**DiagMCMC** - Perform MCMC diagnostic like autocorrelation or trace plots.

**PlotMCMCDiag** - Plot output of **DiagMCMC**, can compare multiple files if more than one argument is used
```bash
./PlotMCMCDiag DiagMCMC_file_1.root Label_1, DiagMCMC_file_2.root Label_2 ...
```
Can process up to 4 files

**RHat** - Performs RHat diagnostic to study if all used chains converged to the same stationary distribution.
```bash
./RHat Ntoys MCMCchain_1.root MCMCchain_2.root MCMCchain_3.root ... [how many you like]
```

**CombineMaCh3Chains** - will combine chains files produced by **MCMC**, enforcing the condition that all the files to combine were made using the exact same software versions and config files
```bash
CombineMaCh3Chains [-h] [-c [0-9]] [-f] [-o <output file>] file1.root [file2.root, file3.root ...]
```
*fileX.root* are the individual spline files to combine, can specify any number, need at least one

-c target compression level for the combined file, default is 1, in line with hadd

-f force overwrite of the combined file if it exists already

-h print usage message and exit

*Output file* (optional) name of the combined file. If not specified, will just use *file1.root*, the first in the list of files, same as *hadd*.


**SmearChain** - Allows you to smear contour. For example after performing sets of study one finds out that used sets of uncertainty doesn't fully cover analysis need. Then one can smear additionally contour.
