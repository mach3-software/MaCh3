# Gradient based sampling and the fast reweight {#GradientSampling}

This page covers four things added together:
- the likelihood gradient;
- the Hamiltonian samplers that use it;
- the reweight backends that make each reweight cheap;
- how an experiment switches each of them on.

Everything is opt-in. An experiment that changes nothing runs exactly as before: the same event loop, the same samplers, the same results.

## Quick start

1. Add the two hooks to your sample handler, if they hold for it ([Experiment hooks](#gs_hooks)). They switch on the fast CPU backend.
2. Run `CheckReweight` on your fit config ([Checking](#gs_checks)). It compares the fast backend with the plain event loop, and the gradient with finite differences.
3. Set `FittingAlgorithm: "Auto"`. The sampler then chooses everything else:
   - parameters with a gradient go to NUTS, and the rest to a random walk;
   - the NUTS step size and the random walk width are tuned during warmup;
   - so is the number of random walk sweeps per step.

   All of these are fixed when warmup ends.

## The likelihood gradient {#gs_gradient}

`SampleHandlerInterface` gains these methods. Each has a default, so existing handlers compile unchanged.

| Method | What it does |
|---|---|
| `SetComputeGradients(bool)` | turns the gradient on or off. Off by default, so gradient-free fitters pay nothing |
| `GetGradLikelihood()` | gives \f$\partial(-\log L)/\partial x\f$ at the last `Reweight()`, for every parameter of `GetGradParHandler()` |
| `GetGradientCoverage()` | says which parameters have a gradient at all |
| `SetGradientsPaused(bool)` | leaves the gradient out of the next reweights. The hybrid sampler uses it during random-walk sweeps |

`SampleHandlerBase` implements them for the Poisson test statistic. The other statistics also depend on the sum of squared weights, whose derivative isn't tracked, so asking for a gradient with them is an error.

Which parameters have a gradient:

| Parameter type | Gradient |
|---|---|
| Normalisation | always (the weight is the parameter, so the derivative is 1) |
| Binned spline | yes. The derivative of the same cubic segment is evaluated with the weight, when a gradient is asked for. In a GPU build the splines, and so their derivatives, are evaluated on the GPU in single precision |
| Oscillation | yes, for NuOscillator's `CUDAProb3Linear`, binned or unbinned ([Oscillation](#gs_osc)) |
| Functional (shift) | no: shifts move events between bins, so a binned likelihood is piecewise constant in them |

A weight's derivative is found from the weight's address, so nothing is stored per event weight. See `WeightDerivativeSource`.

## Samplers {#gs_samplers}

Set them with `General.FittingAlgorithm`:

| Name | Class | What it is |
|---|---|---|
| `HMC` | `HamiltonianMC` | fixed length leapfrog trajectories |
| `NUTS` | `NoUTurnMC` | the No-U-Turn sampler, with dual averaging of the step size |
| `Hybrid` | `HybridMCMC` | NUTS on the parameters with a gradient, a correlated random walk on the rest |
| `Auto` | `HybridMCMC` | Hybrid with every choice made automatically, and plain NUTS when every parameter has a gradient |

Settings, in `General.MCMC` (Auto needs none of them):

| Key | Default | Meaning |
|---|---|---|
| `MaxTreeDepth` | 5 (Auto: 10) | most NUTS doublings per step |
| `DualAveraging` | true | tune the step size during warmup |
| `TargetAcceptance` | 0.8 | NUTS mean acceptance target |
| `WarmupSteps` | NSteps / 2 | steps of tuning before everything is frozen |
| `DeltaT` | 0.1 | starting step size, in units of the metric |
| `AdaptMetric` | false | rebuild the metric when adaption replaces the throw matrix |
| `MetricShrinkage` | 0 | shrink an adapted metric towards its diagonal, from 0 to 1. The prior covariance is never shrunk: its correlations are exact |
| `MetricFromPrior` | false | keep the metric at each handler's first throw matrix |

And in `General.MCMC.Hybrid`:

| Key | Default | Meaning |
|---|---|---|
| `Scheme` | composition | `composition`: sweeps then one NUTS step, every step. `mixture`: one kernel per step at random |
| `RandomWalkTypes` | automatic | parameter types for the random walk. Left out, every parameter without a gradient goes there |
| `RandomWalkNames` | none | names (wildcards allowed) forced into the random walk |
| `RandomWalkSweeps` | 1 (Auto: `auto`) | sweeps per step, or `auto` to choose them during warmup |
| `RandomWalkTimeShare` | 0.4 | with `auto` sweeps, the share of each step's time spent on the random walk |
| `RandomWalkScale` | 1 | starting multiple of the optimal \f$2.38/\sqrt{d}\f$ width |
| `RandomWalkTuneScale` | true | tune the width onto `RandomWalkTargetAcceptance` (0.234) during warmup |
| `RandomWalkCheckEvery` | 1000 (Auto: warmup / 10) | steps per tuning batch |

Everything tuned is frozen at the end of warmup, so the production chain is plain, unadapted MCMC.

## Reweight backends {#gs_backends}

Each sample config can carry a `Reweight` section:

```yaml
Reweight:
  Backend: auto        # auto | legacy | cpu | cuda | metal
  SelfCheckPoints: 3   # points compared with the plain event loop at start-up, 0 to skip
  FastAccess: true     # read shifted cuts and bins through pointers
  LinearShifts: false  # work out linear shifts as arithmetic (see below)
```

- **`legacy`**: the plain event loop, unchanged. It is also the reference the others are checked against.
- **`cpu`** (`CpuReweight`): built from a `ReweightModel`.
  - Every distinct varying weight is stored once, events that share them are grouped, and events that can never be selected are dropped.
  - When the shifts haven't moved since the last full pass (every leapfrog), the groups are summed bin by bin, with no work per event.
  - When they have moved (every random-walk sweep), the experiment's own shift, cut and bin code runs, but each weight costs one multiply.
- **`cuda`** (`CudaReweight`): the cpu backend with its cached pass and gradient on an NVIDIA GPU, in double precision.
  - Full passes, where the shifts moved, stay on the CPU, since they run the experiment's own code.
  - Each pass reads the weight pool into pinned memory, copies it over once and copies the bins (or the gradient) back. The sums run in a fixed order, so the results match the cpu backend to rounding and repeat exactly.
  - It is built with `-DMaCh3_GPU_ENABLED=ON`. The kernels are in `Samples/Reweight/GPU/CudaReweightDevice.cu`, which includes no MaCh3 headers.
- **`metal`** (`MetalReweight`): the cpu backend with its cached pass and gradient on an Apple GPU.
  - It works in single precision, for testing on Macs only.
  - It is built with `-DMaCh3_METAL_ENABLED=ON` and is never picked by `auto`.
- **`auto`**: if the handler declares its additional weights constant, `cuda` in a CUDA build with a GPU present, otherwise `cpu`; `legacy` if the handler doesn't declare it. If the start-up self-check finds a difference, it falls back to `legacy` with a warning.

Every backend gives the plain loop's sums to rounding (single precision: about \f$10^{-6}\f$), the gradient included. Fixed-seed chains are reproducible bit for bit on one thread. With several threads the event sums are not added in a fixed order, so two runs differ in the last digits, and a NUTS chain amplifies that within a few steps into a different (equally valid) chain.

### Linear shifts

With `LinearShifts: true`, shifts that are linear in their parameter (energy scales, resolutions, CVN shifts) are measured once into a table of coefficients. A sweep then works out each event's shifted variables as arithmetic. It doesn't call the shift code or write to the event.

- The model is checked per event at random points; an event it doesn't hold for (a clamp, say) keeps going through the shift code.
- **Compatibility:** with linear shifts on, the events' own shifted variables can be stale after a reweight. Core's plotting and integral functions refresh them first. Experiment code that reads shifted variables itself after a reweight has to call `SampleHandlerBase::SyncShiftedKinematics()`.

## Experiment hooks {#gs_hooks}

Both have defaults that keep the old behaviour:

```cpp
// Every weight AddAdditionalWeightPointers adds is constant for the whole fit, and CalcWeightFunc
// does nothing. Switches on the cpu backend under Backend: auto
bool AdditionalWeightsAreConstant() const override { return true; }

// Which kinematic variables a shift can change. Events failing a cut on any other variable can
// never be selected, and are dropped. Default: every variable counts as shifted
bool IsKinematicParameterShifted(const int KinematicVariable) const override {
  return KinematicVariable == kRecoNeutrinoEnergy || KinematicVariable == kCVNNumu || KinematicVariable == kCVNNue;
}
```

What the fast paths assume, and what checks them:

| Assumption | Checked by |
|---|---|
| A shift depends only on its own parameter values and the nominal event. Otherwise call `InvalidateReweightCache()` whenever what it depends on changes | the start-up self-check and `CheckReweight` |
| Cut variables can be read through `GetPointerToKinematicParameter` | every event, at start-up; the backend falls back if not |
| Binning is uniform and one-dimensional for the flat bin lookup | every event, at start-up; the backend falls back if not |
| Shifts registered as plain functions of one parameter, for the flat shift table and linear shifts | at registration; otherwise the generic path |

`Functional::event_shifts` keeps its reading interface (`size()`, `operator[]`, iteration), but stores each distinct list of shifts once. Code that pushed into it directly should call `Add()` instead.

## Oscillation {#gs_osc}

With NuOscillator's `CUDAProb3Linear`, `OscillationGradient` works out the oscillation weights again with the same formulae carried through dual numbers. That gives each weight's exact derivative with respect to every free oscillation parameter. It is used only if it reproduces NuOscillator's weights; otherwise the oscillation parameters have no gradient, and the hybrid sampler puts them in the random walk.

- **Binned mode** (recommended for gradient fits): one weight per energy bin and channel, so a few thousand evaluations however many events there are.
- **Bin width matters.** Evaluating a bin at its centre is an approximation. DUNE's `OscBinning.root` (10 MeV bins up to 10 GeV) differs from unbinned by up to 2% per bin. A binning refined to steps of 0.05 rad in phase and 0.3% in \f$\Delta E/E\f$ differs by under 0.05% (`Scripts/RefineOscBinning.C` in MaCh3_DUNE). Check this for your own binning before relying on it.
- The parameters need `NuOscName` entries naming CUDAProb3Linear's eight: `sin2_th12`, `sin2_th23`, `sin2_th13`, `dm2_12`, `dm2_23`, `delta_cp`, `path_length`, `matter_density`.

## Where the code lives {#gs_layout}

The sample handlers experiments build on (`SampleHandlerInterface`, `SampleHandlerBase` and the rest of develop's `Samples/`) are where they always were. Everything added for gradient sampling is in subfolders of `Samples/` by topic, apart from `Samples/SampleHandlerReweight.cpp`, which implements `SampleHandlerBase`'s side of it.

| Area | Files | What they hold |
|---|---|---|
| Samplers | `Fitters/Algorithms/Hamiltonian/`: `HamiltonianMC`, `NoUTurnMC`, `HamiltonianMetric`, `DualAveraging.h` | leapfrog with reflection at bounds, NUTS tree building, the metric, step size tuning |
| Hybrid sampler | `Fitters/Algorithms/Hybrid/`: `HybridMCMC.h`, `HybridMCMC.cpp` | set-up and the two schemes (composition, mixture) |
| | `HybridSettings` | reading and checking `General.MCMC.Hybrid` |
| | `HybridBlocks` | which parameters the random walk owns and which NUTS owns |
| | `HybridRandomWalk.cpp`, `RandomWalkTuner.h` | the random walk kernel, its width and the sweep count, and their tuning in warmup |
| | `HybridDelayedAcceptance.cpp` | the random walk screened on a surrogate (`RandomWalkSurrogateStride`) |
| | `HybridFlipMoves.cpp` | flip moves between separated modes |
| Gradient | `Samples/Gradient/LikelihoodGradient`, `Splines/WeightDerivativeSource.h`, `Parameters/SparsePrior.h`, `Samples/Gradient/GradientCheck` | the per-event gradient pass, where each weight's derivative lives, the prior and its gradient, the finite difference check |
| Reweight backends | `Samples/Reweight/`: `ReweightBackend.h`, `ReweightModel`, `EventGroups`, `CpuReweight`; `Samples/Reweight/GPU/`: `CudaReweight*`, `MetalReweight*` | the backend interface, the shared description of the events, and each backend |
| | `Samples/Reweight/ReweightOptions` | reading the `Reweight` section and building the backend it asks for |
| | `Samples/SampleHandlerReweight.cpp` | `SampleHandlerBase`'s side: choosing and checking the backend, the gradient, the surrogate |
| | `Samples/Reweight/FastEventAccess.h`, `SweepTables.h` | shifted cuts and bins read through pointers and flat lookups |
| | `Samples/Reweight/LinearShifts.h`, `LinearShiftBuilder` | the linear shift table, and measuring and checking it |
| | `Samples/Reweight/WeightParameterSnapshot.h` | skipping the oscillation and spline weights when only shifts moved |
| Oscillation | `Samples/Oscillation/`: `OscillationProbability.h`, `OscillationGradient`, `OscillationGradientDevice*`, `NuFastProbability.h` | the oscillation formulae in dual numbers, matched against NuOscillator |
| Flips | `Parameters/FlipMap` | a flip group's map, its self-inverse check and Jacobian, for `ParameterHandlerBase::ProposeFlips()` |

## Checking {#gs_checks}

- **`SampleHandlerBase::CheckReweightBackend`**: compares the backend in use with the plain event loop at random points. It cycles through moving only the shifts, only the rest, and both, and compares the MC prediction, W2, the likelihood and, with gradients on, the gradient.
- **`CheckLikelihoodGradient`** (`Samples/Gradient/GradientCheck.h`): compares the gradient with central finite differences. It picks parameters from each type in turn, and checks the prior gradient for every parameter.
  - Its default step, \f$10^{-4}\f$ of each prior error, suits double precision weights. With single precision spline weights (any GPU build), the finite differences are rounding noise at that step: spline parameters disagree by up to 20% on every backend, the plain loop included. A step of \f$10^{-2}\f$ brings that down to about \f$10^{-4}\f$, with a few small gradients near \f$2\times10^{-3}\f$. `CheckReweight` takes the step as `GradientRelativeStep`.
- **`CheckReweight`** in MaCh3_DUNE runs both and exits non-zero on a failure. Run it on any new configuration, and after any change to the reweight.
