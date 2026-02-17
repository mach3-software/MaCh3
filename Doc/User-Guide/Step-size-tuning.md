# Step-Size Tuning {#step-size-tuning}

# What is and why you need to know what step size tunning <br />
It is good practice to run an MCMC diagnostic before looking at posterior distributions, although in most cases you will stumble on the need to diagnose after checking posteriors. If your posterior looks like this either you do have not enough steps or wrongly tune the chain. <br />
<img src="https://user-images.githubusercontent.com/45295406/213859228-8dbf752d-8c3c-432a-a3d4-2fd33d0f1863.png" width="400" height="400"/>

Before discussing step size tunning first we need to understand how a step is proposed.
* Gaussian throw - random number drawn from Gauss distribution starting at the previous step multiplied with a spread equal to parameter error.
* Correlated throw - the proposed step for two correlated parameters should be more likely to change direction in the same way, hence we include correlation using the Cholesky matrix.
* Individual step scale - user selected value for each parameter by a default =1, by step size tunning in most cases we mean modifying this value.
* Global step scale - user selected value same cor all parameters in a given covariance class. For example, the value is the same for all xsec-based parameters.

## MCMC diagnostic <br />
There are several plots worth studying for MCMC diagnostic you can find the executables which produce them in the `Diagnostics` folder, while the text is based on **[1]**.
### Autocorealtions <br />
we can study chain autocorrelation, which tells us how much particular steps are correlated. To test it, we introduce a variable _Lag(n)= corr(Xi, Xi−n)_ which tells us how much correlated are steps that are _n_ steps apart, where _i_ is the maximal considered distanced, here _i = 25000_. Fig. shows autocorrelations for studied chains since we want our steps to be more random, and less correlated to quickly converge. The rule of thumb is for autocorrelation to drop below 0.2 for _Lag(n = 10000)_. This isn’t a strict criterion so if sometimes autocorrelations drop slightly slower than the blue line in our Figure it’s not a problem. <br />

Example of well-tuned step scale (colours represent chains which have different starting positions) <br />
<img src="https://user-images.githubusercontent.com/45295406/213926048-896ff3e2-741c-483b-89e6-2c7d025c984e.png" width="400" height="400"/>

If your autocorrelation looks like this though you really should increase step size. One exception would be parameter which has no effect. Imagine you run ND only fit but the parameter affects FD only. Then it is expected autocorrelation to look badly. <br />
<img src="https://user-images.githubusercontent.com/45295406/213938757-ede5079b-dd6f-4beb-88a7-cd2525ef5805.png" width="400" height="400"/>

### Trace <br />
Fig. shows the trace which is the value of a chosen parameter at each step. It can be seen that at first, the chains have different traces but after a thousand steps, they start to stabilise and oscillates around a very similar value, indicating that the chain converged and a stationary state was achieved. <br />

<img src="https://user-images.githubusercontent.com/45295406/213926271-c8ae49c6-4ee3-45d5-a377-79e959ffdc4a.png" width="400" height="400"/>

### AcceptanceProbability Batch <br />
Fig. shows a mean value of acceptance probability (_A(θ',θ)_) in an interval of
5k steps (batched mean). This quantity is quite high at the beginning indicating the chain didn’t
converge. When the chain gets close to the stationary state it starts to stabilise. Orange 
stabilised fastest while blue and green are slowly catching up, but the red didn’t converge yet. <br />

<img src="https://user-images.githubusercontent.com/45295406/213926134-2cecc97b-4e03-4b47-ad58-4cda0d2090ed.png" width="400" height="400"/>

### R Hat <br />
Usually, we run several chains which are later combined. There is a danger that not all chains will converge then using them will bias results. R hat is meant to estimate whether chains converged successfully or not. According to Gelman, you should calculate R hat for at least 4 chains and if R hat > 1.1 then it might indicate wrong chains convergence. Below you can find an example of chains which wrongly converged and one which successfully 

Chains covnerged to different values
<img src="https://user-images.githubusercontent.com/45295406/219701548-da171c94-9fa6-4857-8e22-372e278b91dc.png" width="400" height="400"/>

Successfuly converged chains
<img src="https://user-images.githubusercontent.com/45295406/219701496-abd8070c-616a-4637-ac64-4e1af14dd459.png" width="400" height="400"/>

### Geweke <br />
Geweke Diagnostic helps to define what burn it should it be. You should select burn-in around 15% in this case as this is where distirbution stabilizes.
<img src="https://github.com/mach3-software/MaCh3/assets/45295406/61f6f86b-abda-489d-afbd-f758b7ecef27.png" width="400" height="400"/>

## Global Step Scale <br />

According to **[2]** (Section 3.2.1), the global step size should be equal to 2.38^2/N_params. Keep in mind this refers to the global step scale for a given covariance object like xsec covariance or detector covariance.

## Manual Tunning <br />
This procedure is very tedious and requires intuition of how a given parameter behaves. It is a bit of dark magic however skilled users should be able to tune it relatively fast compared with the non-skilled user. The process is as follows, you run the chain, run the `diagnostic` executable look at plots adjust the step scale then run again fit and the process repeats. Each time you should look at autocorrelations, traces, etc. (see discussion above). Another important trick is not to run full fit. Instead of running 10M chain, you might run 200k. Number of steps depends on number of parameters and _Lag(n = ?)_ you are interested in. <br />

There are a few things you should be aware of when tunning:
* Parameters with a broad range may have a higher step scale, while those with a narrow should have smaller ones to reduce the probability of going out of bounds.
* Highly correlated should have a similiar step-scale, for edge cases like ~100% step scale should be identical!
* Autocorelations should drop below 0.2 for _Lag(n = 10000)_. If it drops immediately then the step scale is too big.
* Study trace, if is converging, exploring phase space fast enough. Exploring too fast is wrong.
* Study acceptance probability. If every step is accepted then the scale is too small, while if barely any step is getting accepted you might consider decreasing the step scale.
* Doing LLH scan and assigning a step scale based on such a result is also a good idea.

The last point is that data fit may require a different tuning than the Asimov fit. Still, if you tune for Asimov it should be easy to re-tune for a data fit. 

### Manually Tuning Individual Step Scales in MaCh3
Currently MaCh3 configures systematics in two ways. The first is simply through YAML configs which are used for cross-section systematics. Changing the step scale of parameter in a YAML config is simply a case of modifying the "StepScale" option for the parameter in the YAML file.

Some systematics, for example oscillation parameters, still use an `XML->ROOT` pipeline. Here the systematics are initially defined in an XML file first, for example in DUNE-MaCh3 they are located [here](https://github.com/DUNE/MaCh3_DUNE/blob/main/utils/oscMatrixMaker/osc_covariance_DUNE_PDG2021_v1.xml). Individual parameter step scales can then be modified through the `stepscale` option in the XML (although this may vary between systematics as it's a legacy system!). Finally, these can then be converted into a ROOT file with a suitable script, in the DUNE example this is located [here](https://github.com/DUNE/MaCh3_DUNE/blob/main/utils/oscMatrixMaker/makeOscMatrix.py).
<br />

## Adaptive MCMC <br />
Hopefully by this point you've realised that step size tuning is
1. Hard
2. Tedious

Thankfully we can automate it! It turns out that, provided you Markov chain satisfies the Markov chain central limit theorem [2] it's optimal to propose steps from the posterior covariance matrix (multiplied by a scale factor). 

The config then has the following options which let you tune this
```Yaml
AdaptionOptions:
  Settings:
    StartThrow:      # [int] At which step do we start throwing from our new matrix?
    StartUpdate:     # [int] At which step do we start adding steps to the covariance matrix?
    EndUpdate:       # [int] At which step do we stop adding steps + fix the covariance matrix?
    UpdateStep:      # [int] How often do we want to update i.e. change matrix we're throwing from?
    SaveNIterations: # [int] How often do we want to save this to file
    OutputFileName:  # [str] Name of the file containing the matrix

  Covariance:
    xsec_cov: # Name of covariance matrix
      DoAdaption:   # [bool] Adapt this parameter
      MatrixBlocks: # [list[list[int]]] Split the matrix into blocks. Form of this is [[block start, block end]]. Note you can nest this like [[start, end, start, end]] to define a single block over multiple bits of the matrix

      # Robbins-Monro Settings
      UseRobbinsMonro: # [bool] Do we you want to use Robbins-Monro to adapt matrix-scale 
      TargetAcceptance: 0.234 # Target acceptance rate
      TotalRobbinsMonroRestarts: 5 # How many times do we restart the Robbins-Monro count
      AcceptanceRateBatchSize: 5000 # After how many steps we restart the Robbins-Monro count

      # External stuff 
      UseExternalMatrix:      # [bool] Use a matrix from a file 
      ExternalMatrixFileName: # [str] Name of file containing matrix
      ExternalMatrixName:     # [str] Name of matrix in file
      ExternalMeansName:      # [str] Name of vector containing means, useful if you want to continue adapting from a previous chain
```
Adaptive step size tuning is still a little bit fiddly but works well if you have "Gaussian-ish" parameters [for example cross-section!]. Generally I'd recommend updating every few 1000 steps and stop updating after around 1,000,000 steps. <br />

For more info on adaptive MCMC please see [section 16](https://github.com/mach3-software/MaCh3/wiki/16.-Adaptive-MCMC) in the wiki!

# References <br />
[1] Kamil Skwarczynski PhD thesis <br />
[2] https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/s13634-020-00675-6.pdf  <br />

If you have complaints blame: Kamil Skwarczynski