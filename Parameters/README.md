# Parameters

## ParameterHandlerBase Design Decisions

Elementary objects needed for sampling:
- ParameterList
  - Read from YAML. Knows everything about each systematic parameter: name, error, possibly-correlated prior, restrictions on parameter applicability etc...

- StepProposer
  - Knows nothing about the meaning of parameters.
  - Generates parameter values that are steps away from some current point according to a symmetric proposal function. In this case specified fully by a multivariate gaussian.
  - Keeps track of parameter step scales and possibly updates the proposal matrix as steps are accepted.
  - Possibly steps in a rotated parameter space.
  - Possibly has 'special proposals': flips and circular
  - Needs to know how to translate from its parameter basis back to the systematic parameter basis, that is the job of ParameterList

Open Questions:
- Should the likelihood be calculated in the systematic (full) parameter space or the proposal (possibly truncated) parameter space. NPar might differ significantly?
