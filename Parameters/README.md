# Parameters

## ParameterHandlerBase Design Decisions

Elementary objects needed for sampling:
- Parameter list
  - Read from YAML. Knows everything about each systematic parameter: name, error, possibly-correlated prior, restrictions on parameter applicability

- StepProposer
  - Knows nothing about the meaning of parameters
  - Generates parameter values that are steps away from some current point according to a proposal matrix.
  - Keeps track of parameter step scales and possibly updates the proposal matrix as steps are accepted.
  - Possibly steps in a rotated parameter space.
  - Possibly has 'special proposals': flips and circular
