#include "Parameters/ParameterHandlerBase.h"
#include <regex>

ParameterHandlerBase::ParameterHandlerBase()
    : pca.enabled{false},
special_proposal.enabled{false}, rng.gaus(0), settings.use_adaptive{false},
settings.PrintLength{35} {}

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string const &name,
                                           std::string const &inputFile)
    : ParameterHandlerBase() {
  // ********************************************
  MACH3LOG_DEBUG("Constructing instance of ParameterHandler, named: {}", name);

  SetName(name);
  config.inputFile = file;
}

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string name, std::string file,
                                           double threshold, int FirstPCA,
                                           int LastPCA)
    : ParameterHandlerBase(name, file) {
  // ********************************************

  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  TFile infile(file.c_str(), "READ");
  if (infile.IsZombie()) {
    MACH3LOG_ERROR("Could not open input covariance ROOT file {} !!!", file);
    MACH3LOG_ERROR("Was about to retrieve matrix with name {}", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto *CovMat = infile.Get<TMatrixDSym>(name.c_str());

  if (!CovMat) {
    MACH3LOG_ERROR("Could not find covariance matrix name {} in file {}", name,
                   file);
    MACH3LOG_ERROR("Are you really sure {} exists in the file?", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  params.covariance = M3::MakePosDef(M3::ROOTToEigen(*CovMat));

  for (int iThread = 0; iThread < M3::GetNThreads(); iThread++) {
    random_number.emplace_back(0);
  }

  MACH3LOG_INFO("Created covariance matrix named: {}", GetName());
  MACH3LOG_INFO("from file: {}", file);
}

void ParameterHandlerBase::AddParameters(std::vector<ParamInfo> const &params) {

  for (auto const &p : params) {
    params.name.push_back(p.name);
    params.fancy_name.push_back(p.fancy_name);
    params.samples.push_back(p.affected_samples);
  }

  size_t new_block_start = params.prefit.size();

  params.prefit.conservativeResize(params.name.size());
  params.error.conservativeResize(params.name.size());
  params.lowbound.conservativeResize(params.name.size());
  params.upbound.conservativeResize(params.name.size());
  params.flatprior.conservativeResize(params.name.size());
  params.fixed.conservativeResize(params.name.size());
  steps.scale.conservativeResize(params.name.size());

  params.covariance.conservativeResize(params.name.size(), params.name.size());

  for (size_t i = 0; i < params.size(); ++i) {
    params.prefit[new_block_start + i] = params[i].prefit;
    params.error[new_block_start + i] = params[i].error;
    params.lowbound[new_block_start + i] = params[i].bounds[0];
    params.upbound[new_block_start + i] = params[i].bounds[1];
    params.flatprior[new_block_start + i] = params[i].flatprior;
    params.fixed[new_block_start + i] = false;
    steps.scale[new_block_start + i] = params[i].stepscale;
    params.covariance(new_block_start + i, new_block_start + i) =
        params[i].error * params[i].error;
  }

  params.inv_covariance = MatrixXd(0);
}

void ParameterHandlerBase::SetParameterCorrelation(int pidi, int pidj,
                                                   double corr) {
  if (pidi == pidj) {
    MACH3LOG_ERROR(
        "AddParameterCorrelation cannot be used to set covariance "
        "matrix diagonal elements: ({0},{0}) attempted to be set to {1}",
        pidi, corr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  params.covariance(pidi, pidj) =
      corr *
      std::sqrt(params.covariance(pidi, pidi) * params.covariance(pidj, pidj));
  params.covariance(pidj, pidi) = params.covariance(pidi, pidj);

  params.inv_covariance = MatrixXd(0);
}

void ParameterHandlerBase::SetParameterAllCorrelations(
    int paramid, std::map<std::string, double> const &correlations) {

  params.covariance.row(paramid) =
      Eigen::VectorXd::Zeros(params.covariance.rows());
  params.covariance.col(paramid) =
      Eigen::VectorXd::Zeros(params.covariance.rows());

  params.covariance(paramid, paramid) =
      params.error[paramid] * params.error[paramid];

  for (auto const &[other_name, corr] : correlations) {
    SetParameterCorrelation(paramid, GetParIndex(other_name), corr);
  }
}

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(
    const std::vector<std::string> &YAMLFiles, std::string name,
    double threshold, int FirstPCA, int LastPCA)
    : ParameterHandlerBase(name, YAMLFiles[0]) {
  // ********************************************

  MACH3LOG_INFO("Constructing instance of ParameterHandler using: ");
  for (auto const &yf : YAMLFile) {
    MACH3LOG_INFO("  {}", yf);
  }
  MACH3LOG_INFO("as an input.");

  settings.pca = true;
  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("Principal component analysis but given the threshold for "
                  "the principal components to be less than 0, or greater than "
                  "(or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: {}", threshold);
    MACH3LOG_INFO("Disabling PCA...");
    settings.pca = false;
  }

  config.YAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for (auto const &yfn : YAMLFiles) {
    for (const auto &[idx, node] : M3::OpenConfig(yf)["Systematics"]) {
      config.YAMLDoc["Systematics"].push_back(node);
    }
  }

  std::map<int, std::map<std::string, double>> parameter_correlations;

  std::vector<ParamInfo> param_infos;

  // ETA - read in the systematics.
  for (auto const &[idx, node] : config.YAMLDoc["Systematics"]) {
    size_t param_id = param_infos.size();

    auto const &pardef = node["Systematic"];

    auto fancy_name = GetFromManager<std::string>(pardef["Names"]["FancyName"],
                                                  __FILE__, __LINE__);
    auto prefit = GetFromManager<double>(
        pardef["ParameterValues"]["PreFitValue"], __FILE__, __LINE__);

    auto error = GetFromManager<double>(pardef["Error"], __FILE__, __LINE__);

    if (error <= 0) {
      MACH3LOG_ERROR("Error for param {}({}) is negative and equal to {}",
                     fancy_name, param_id, error);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    auto stepscale =
        GetFromManager<double>(pardef["StepScale"]["MCMC"], __FILE__, __LINE__);

    auto bounds = GetBounds(pardef["ParameterBounds"]);

    auto flatprior =
        GetFromManager<bool>(pardef["FlatPrior"], false, __FILE__, __LINE__);

    auto samplenames = GetFromManager<std::vector<std::string>>(
        pardef["SampleNames"], {}, __FILE__, __LINE__);

    // Allow to fix param, this setting should be used only for params which are
    // permanently fixed like baseline, please use global config for fixing
    // param more flexibly
    auto fixed =
        GetFromManager<bool>(pardef["FixParam"], false, __FILE__, __LINE__);

    if (pardef["SpecialProposal"]) {
      special_proposals.push_back(param_id);
      EnableSpecialProposal(pardef["SpecialProposal"]);
    }

    if (pardef["Correlations"]) {
      for (auto const &[_, it] : pardef["Correlations"]) {
        for (auto const &[key, corr] : it) {
          parameter_correlations[param_id][key.as<std::string>()] =
              corr.as<double>();
        }
      }
    }

    param_infos.emplace_back({fancy_name,
                              fancy_name,
                              prefit,
                              error,
                              stepscale,
                              {bounds[0], bounds[1]},
                              flatprior,
                              samplenames});
  }

  AddParameters(param_infos);

  for (auto const &[paramid, correlations] : parameter_correlations) {
    SetParameterCorrelations(paramid, correlations);
  }

  params.covariance = M3::MakePosDef(params.covariance);

  Tunes = ParameterTunes(config.YAMLDoc["Systematics"]);

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for (const auto &file : YAMLFile) {
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total",
                params.prefit.size());
  MACH3LOG_INFO("----------------");
}

// ********************************************
void ParameterHandlerBase::ConstructPCA(const double eigen_threshold,
                                        int FirstPCAdpar, int LastPCAdpar) {
  // ********************************************
  if (AdaptiveHandler) {
    MACH3LOG_ERROR("Adaption has been enabled and now trying to enable PCA. "
                   "Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Check whether first and last pcadpar are set and if not just PCA everything
  if (FirstPCAdpar == -999 || LastPCAdpar == -999) {
    if (FirstPCAdpar == -999 && LastPCAdpar == -999) {
      FirstPCAdpar = 0;
      LastPCAdpar = covMatrix->GetNrows() - 1;
    } else {
      MACH3LOG_ERROR("You must either leave FirstPCAdpar and LastPCAdpar at "
                     "-999 or set them both to something");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  PCAObj.ConstructPCA(params.covariance, FirstPCAdpar, LastPCAdpar,
                      eigen_threshold);
}

// ********************************************
void ParameterHandlerBase::EnableSpecialProposal(const YAML::Node &param,
                                                 const int Index) {
  // ********************************************
  special_proposal.enabled = true;

  bool CircEnabled = false;
  bool FlipEnabled = false;

  if (param["CircularBounds"]) {
    CircEnabled = true;
  }

  if (param["FlipParameter"]) {
    FlipEnabled = true;
  }

  if (!CircEnabled && !FlipEnabled) {
    MACH3LOG_ERROR("None of Special Proposal were enabled even though param "
                   "{}, has SpecialProposal entry in Yaml",
                   GetParFancyName(Index));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (CircEnabled) {
    special_proposal.CircularBoundsIndex.push_back(Index);
    special_proposal.CircularBoundsValues.push_back(
        Get<std::pair<double, double>>(param["CircularBounds"], __FILE__,
                                       __LINE__));
    MACH3LOG_INFO(
        "Enabling CircularBounds for parameter {} with range [{}, {}]",
        GetParFancyName(Index),
        special_proposal.CircularBoundsValues.back().first,
        special_proposal.CircularBoundsValues.back().second);
    // KS: Make sure circular bounds are within physical bounds. If we are
    // outside of physics bound MCMC will never explore such phase space region
    if (special_proposal.CircularBoundsValues.back().first <
            _fLowBound.at(Index) ||
        special_proposal.CircularBoundsValues.back().second >
            _fUpBound.at(Index)) {
      MACH3LOG_ERROR("Circular bounds [{}, {}] for parameter {} exceed "
                     "physical bounds [{}, {}]",
                     special_proposal.CircularBoundsValues.back().first,
                     special_proposal.CircularBoundsValues.back().second,
                     GetParFancyName(Index), _fLowBound.at(Index),
                     _fUpBound.at(Index));
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  if (FlipEnabled) {
    special_proposal.FlipParameterIndex.push_back(Index);
    special_proposal.FlipParameterPoint.push_back(
        Get<double>(param["FlipParameter"], __FILE__, __LINE__));
    MACH3LOG_INFO("Enabling Flipping for parameter {} with value {}",
                  GetParFancyName(Index),
                  special_proposal.FlipParameterPoint.back());
  }

  if (CircEnabled && FlipEnabled) {

    const double fp = special_proposal.FlipParameterPoint.back();
    const double low = special_proposal.CircularBoundsValues.back().first;
    const double high = special_proposal.CircularBoundsValues.back().second;

    if (fp < low || fp > high) {
      MACH3LOG_ERROR("FlipParameter value {} for parameter {} is outside the "
                     "CircularBounds [{}, {}]",
                     fp, GetParFancyName(Index), low, high);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Sanity check: ensure flipping any x in [low, high] keeps the result in
    // [low, high]
    const double flipped_low = 2 * fp - low;
    const double flipped_high = 2 * fp - high;
    const double min_flip = std::min(flipped_low, flipped_high);
    const double max_flip = std::max(flipped_low, flipped_high);

    if (min_flip < low || max_flip > high) {
      MACH3LOG_ERROR("Flipping about point {} for parameter {} would leave "
                     "circular bounds [{}, {}]",
                     special_proposal.FlipParameterPoint.back(),
                     GetParFancyName(Index), low, high);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
}

// ********************************************
// Set the covariance matrix for this class
void ParameterHandlerBase::SetCovMatrix(TMatrixDSym *cov) {
  // ********************************************
  if (cov == nullptr) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to {}",
                   __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (cov->GetNrows() != params.prefit.size()) {
    MACH3LOG_ERROR(
        "Passed covariance matrix size {0}x{0}, but have {1} parameters.",
        ov->GetNrows(), params.prefit.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  params.covariance = M3::ROOTToEigen(*cov);
  params.inv_covarance = params.covariance.invert();

  SetThrowMatrix(cov);
}

// ********************************************
// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void ParameterHandlerBase::SetPar(int i, double val) {
  // ********************************************
  MACH3LOG_INFO("Overriding {}: ", GetParName(i));
  MACH3LOG_INFO(
      "steps.proposed ({}), steps.current ({}), params.prefit ({}) to ({})",
      steps.proposed[i], steps.current[i], params.prefit[i], val);

  steps.proposed[i] = val;
  steps.current[i] = val;
  params.prefit[i] = val;
}

// *************************************
// Throw the parameters according to the covariance matrix
// This shouldn't be used in MCMC code ase it can break Detailed Balance;
void ParameterHandlerBase::ThrowParameters() {
  // *************************************

  // First draw a new random_vector
  Randomize();

  Eigen::VectorXd deltas = steps.l_proposal * throws.random_vector;

  for (int i = 0; i < throws.values.size(); ++i) {

    int tries = 0;
    while (((params.prefit[i] + deltas[i]) < params.lowbound) ||
           ((params.prefit[i] + deltas[i]) >= params.upbound)) {
      deltas[i] = steps.l_proposal.col(i) * rng.gaus(rng.e1);
      tries++;
      if (tries > 10000) {
        // KS: Since we are multithreading there is danger that those messages
        // will be all over the place, small price to pay for faster code
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", tries,
                      i);
        MACH3LOG_WARN("Matrix: {}", settings.name);
        MACH3LOG_WARN("Param: {}", params.name[i]);
        MACH3LOG_WARN("Setting _fPropVal:  {} to {}", steps.proposed[i],
                      params.prefit[i]);
        MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
        deltas[i] = 0;
      }
    }

    steps.proposed = params.prefit + deltas;
    steps.current = steps.proposed;
  }
}

// *************************************
// Throw each parameter within their 1 sigma range
// Used to start the chain in different states
void ParameterHandlerBase::RandomConfiguration() {
  // *************************************
  // Have the 1 sigma for each parameter in each covariance class, sweet!
  // Don't want to change the prior array because that's what determines our
  // likelihood Want to change the _fPropVal, _fCurrVal, _fPreFitValue
  // _fPreFitValue and the others will already be set
  for (int i = 0; i < _fNumPar; ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (IsParameterFixed(i))
      continue;
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    const double paramrange = _fUpBound[i] - _fLowBound[i];
    const double sigma = sqrt((*covMatrix)(i, i));
    double throwrange = sigma;
    if (paramrange < sigma)
      throwrange = paramrange;

    _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1) * throwrange;
    // Try again if we the initial parameter proposal falls outside of the range
    // of the parameter
    int throws = 0;
    while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
      if (throws > 1000) {
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws,
                      i);
        MACH3LOG_WARN("Matrix: {}", matrixName);
        MACH3LOG_WARN("Param: {}", _fNames[i]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      _fPropVal[i] =
          _fPreFitValue[i] + random_number[0]->Gaus(0, 1) * throwrange;
      throws++;
    }
    MACH3LOG_INFO("Setting current step in {} param {} = {} from {}",
                  matrixName, i, _fPropVal[i], _fCurrVal[i]);
    _fCurrVal[i] = _fPropVal[i];
  }
  if (pca)
    PCAObj.TransferToPCA();
}

// *************************************
// Set a single parameter
void ParameterHandlerBase::SetSingleParameter(const int parNo,
                                              const double parVal) {
  // *************************************
  steps.current[parNo] = parVal;
  steps.proposed[parNo] = parVal;
  MACH3LOG_DEBUG("Setting {} (parameter {}) to {})", GetParName(parNo), parNo,
                 parVal);
  if (pca)
    PCAObj.TransferToPCA();
}

// ********************************************
void ParameterHandlerBase::SetParCurrProp(const int parNo,
                                          const double parVal) {
  // ********************************************
  SetSingleParameter(parNo, parVal);
}

// ************************************************
// Propose a step for the set of systematics parameters this covariance class
// holds
void ParameterHandlerBase::ProposeStep() {
  // ************************************************
  // Make the random numbers for the step proposal
  Randomize();
  CorrelateSteps();

  // KS: According to Dr Wallace we update using previous not proposed step
  // this way we do special proposal after adaptive after.
  // This way we can shortcut and skip rest of proposal
  if (!doSpecialStepProposal)
    return;

  SpecialStepProposal();
}

// ************************************************
void ParameterHandlerBase::SpecialStepProposal() {
  // ************************************************
  /// @warning KS: Following Asher comment we do "Step->Circular Bounds->Flip"

  // HW It should now automatically set dcp to be with [-pi, pi]
  for (size_t i = 0; i < special_proposal.CircularBoundsIndex.size(); ++i) {
    const int index = special_proposal.CircularBoundsIndex[i];
    if (!IsParameterFixed(index))
      CircularParBounds(index, special_proposal.CircularBoundsValues[i].first,
                        special_proposal.CircularBoundsValues[i].second);
  }

  // Okay now we've done the standard steps, we can add in our nice flips
  // hierarchy flip first
  for (size_t i = 0; i < special_proposal.FlipParameterIndex.size(); ++i) {
    const int index = special_proposal.FlipParameterIndex[i];
    if (!IsParameterFixed(index))
      FlipParameterValue(special_proposal.FlipParameterIndex[i],
                         special_proposal.FlipParameterPoint[i]);
  }
}

// ************************************************
// "Randomize" the parameters in the covariance class for the proposed step
// Used the proposal kernel and the current parameter value to set proposed step
// Also get a new random number for the randParams
void ParameterHandlerBase::Randomize() _noexcept_ {
  // ************************************************
  throws.random_vector = Eigen::VectorXd::NullaryExpr(
      params.prefit.size(), [&](int) { return rng.gaus(rng.e1); });
}

// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its
// current value + some correlated throw
void ParameterHandlerBase::CorrelateSteps() _noexcept_ {
  // ************************************************

  if (!pca) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i) {
      if (!IsParameterFixed(i) > 0.) {
        _fPropVal[i] = _fCurrVal[i] +
                       corr_throw[i] * _fGlobalStepScale * _fIndivStepScale[i];
      }
    }
    // If doing PCA throw uncorrelated in PCA basis (orthogonal basis by
    // definition)
  } else {
    PCAObj.CorrelateSteps(_fIndivStepScale, _fGlobalStepScale, randParams,
                          corr_throw);
  }
}
// ********************************************
// Update so that current step becomes the previously proposed step
void ParameterHandlerBase::AcceptStep() _noexcept_ {
  // ********************************************
  steps.current = steps.proposed;

  if (AdaptiveHandler) {
    AdaptiveHandler.IncrementAcceptedSteps();
  }
}

// *************************************
// HW: This method is a tad hacky but modular arithmetic gives me a headache.
void ParameterHandlerBase::CircularParBounds(const int index,
                                             const double LowBound,
                                             const double UpBound) {
  // *************************************
  if (_fPropVal[index] > UpBound) {
    _fPropVal[index] =
        LowBound + std::fmod(_fPropVal[index] - UpBound, UpBound - LowBound);
  } else if (_fPropVal[index] < LowBound) {
    _fPropVal[index] =
        UpBound - std::fmod(LowBound - _fPropVal[index], UpBound - LowBound);
  }
}

// *************************************
void ParameterHandlerBase::FlipParameterValue(const int index,
                                              const double FlipPoint) {
  // *************************************
  if (random_number[0]->Uniform() < 0.5) {
    _fPropVal[index] = 2 * FlipPoint - _fPropVal[index];
  }
}

// ********************************************
// Throw the proposed parameter by mag sigma
// Should really just have the user specify this throw by having argument double
void ParameterHandlerBase::ThrowParProp(const double mag) {
  // ********************************************
  Randomize();
  if (!pca) {
    // Make the correlated throw
    M3::MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams,
                          _fNumPar);
    // Number of sigmas we throw
    for (int i = 0; i < _fNumPar; i++) {
      if (!IsParameterFixed(i) > 0.)
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i] * mag;
    }
  } else {
    PCAObj.ThrowParProp(mag, randParams);
  }
}
// ********************************************
// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void ParameterHandlerBase::ThrowParCurr(const double mag) {
  // ********************************************
  Randomize();
  if (!pca) {
    // Get the correlated throw vector
    M3::MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams,
                          _fNumPar);
    // The number of sigmas to throw
    // Should probably have this as a default parameter input to the function
    // instead
    for (int i = 0; i < _fNumPar; i++) {
      if (!IsParameterFixed(i) > 0.) {
        _fCurrVal[i] = corr_throw[i] * mag;
      }
    }
  } else {
    PCAObj.ThrowParCurr(mag, randParams);
  }
}
// ********************************************
// Function to print the prior values
void ParameterHandlerBase::PrintNominal() const {
  // ********************************************
  MACH3LOG_INFO("Prior values for {} ParameterHandler:", GetName());
  for (int i = 0; i < params.name.size(); i++) {
    MACH3LOG_INFO("    {}   {} ", GetParFancyName(i), GetParInit(i));
  }
}

// ********************************************
// Function to print the prior, current and proposed values
void ParameterHandlerBase::PrintNominalCurrProp() const {
  // ********************************************
  MACH3LOG_INFO("Printing parameters for {}", GetName());
  MACH3LOG_INFO("{:<30} {:<10} {:<10} {:<10}", "Name", "Prior", "Current",
                "Proposed");
  for (int i = 0; i < params.prefit.size(); ++i) {
    MACH3LOG_INFO("{:<30} {:<10.2f} {:<10.2f} {:<10.2f}", GetParFancyName(i),
                  params.prefit[i], params.current[i], params.proposed[i]);
  }
}

// ********************************************
// Get the likelihood in the case where we want to include priors on the
// parameters _fFlatPrior stores if we want to evaluate the likelihood for the
// given parameter
//                    true = don't evaluate likelihood (so run without a prior)
//                    false = evaluate likelihood (so run with a prior)

double ParameterHandlerBase::CalcLikelihood() const _noexcept_ {
  // ********************************************
  // filter out parameters with flat priors from the pull
  auto diff =
      params.flatprior.select(Eigen::VectorXd::Zeros(params.prefit.size()),
                              params.proposed - params.prefit);
  return diff.T * params.inv_covariance * diff;
}

// ********************************************
int ParameterHandlerBase::CheckBounds() const _noexcept_ {
  // ********************************************
  return ((params.proposed.array() < params.lowbound.array()) ||
          (params.proposed.array() >= params.upbound.array()))
      .count();
}

// ********************************************
double ParameterHandlerBase::GetLikelihood() {
  // ********************************************
  // Default behaviour is to reject negative values + do std llh calculation
  const int NOutside = CheckBounds();

  if (NOutside > 0)
    return NOutside * M3::_LARGE_LOGL_;

  return CalcLikelihood();
}

// ********************************************
// Sets the proposed parameters to the prior values
void ParameterHandlerBase::SetParameters(const std::vector<double> &pars) {
  // ********************************************
  // If empty, set the proposed to prior
  if (pars.empty()) {
    steps.proposed = params.prefit;
    // If not empty, set the parameters to the specified
  } else {
    if (pars.size() != size_t(_fNumPar)) {
      MACH3LOG_ERROR("Parameter arrays of incompatible size! Not changing "
                     "parameters! {} has size {} but was expecting {}",
                     matrixName, pars.size(), _fNumPar);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    for (int i = 0; i < pars.size(); i++) {
      // Make sure that you are actually passing a number to set the parameter
      // to
      if (std::isnan(pars[i])) {
        MACH3LOG_ERROR("Trying to set parameter value to a nan for parameter "
                       "{} in matrix {}. This will not go well!",
                       GetParName(i), matrixName);
        throw MaCh3Exception(__FILE__, __LINE__);
      } else {
        steps.proposed[i] = pars[i];
      }
    }
  }
}

// ********************************************
void ParameterHandlerBase::SetBranches(TTree &tree, bool SaveProposal) {
  // ********************************************
  // loop over parameters and set a branch
  for (int i = 0; i < params.name.size(); ++i) {
    tree.Branch(params.name[i].c_str(), params.current.data() + i,
                Form("%s/D", params.name[i].c_str()));
  }
  // When running PCA, also save PCA parameters
  if (pca) {
    PCAObj.SetBranches(tree, SaveProposal, _fNames);
  }
  if (SaveProposal) {
    // loop over parameters and set a branch
    for (int i = 0; i < params.name.size(); ++i) {
      tree.Branch(Form("%s_Prop", params.name[i].c_str()),
                  params.proposed.data() + i,
                  Form("%s_Prop/D", params.name[i].c_str()));
    }
  }
  if (use_adaptive && AdaptiveHandler.GetUseRobbinsMonro()) {
    tree.Branch(Form("GlobalStepScale_%s", GetName().c_str()),
                &_fGlobalStepScale,
                Form("GlobalStepScale_%s/D", GetName().c_str()));
  }
}

// ********************************************
void ParameterHandlerBase::SetStepScale(const double scale,
                                        const bool verbose) {
  // ********************************************
  if (scale <= 0) {
    MACH3LOG_ERROR("Invalid step scale {}", scale);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (verbose) {
    MACH3LOG_INFO("{} setStepScale() = {}", GetName(), scale);
    const double SuggestedScale = 2.38 / std::sqrt(_fNumPar);
    if ((std::fabs(scale - SuggestedScale) / SuggestedScale) > 1) {
      MACH3LOG_WARN(
          "Defined Global StepScale is {}, while suggested suggested {}", scale,
          SuggestedScale);
    }
  }
  steps.global_scale = scale;
}

// ********************************************
int ParameterHandlerBase::GetParIndex(const std::string &name) const {
  // ********************************************
  auto idx = std::find(params.name.begin(), params.name.end(), name);
  if (idx == params.name.end()) {
    MACH3LOG_ERROR("No parameter with the name {} exists.", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return std::distance(params.name.begin(), idx);
}

// ********************************************
void ParameterHandlerBase::SetFixAllParameters() {
  // ********************************************
  params.fixed = Eigen::VectorXi::Ones(params.prefit.size());
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const int i) {
  // ********************************************
  params.fixed[i] = 1;
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const std::string &name) {
  // ********************************************
  SetFixParameter(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::SetFreeAllParameters() {
  // ********************************************
  params.fixed = Eigen::VectorXi::Zeros(params.prefit.size());
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const int i) {
  // ********************************************
  params.fixed[i] = 0;
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const std::string &name) {
  // ********************************************
  SetFreeParameter(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::ToggleFixAllParameters() {
  // ********************************************
  throw std::runtime_error("Don't Toggle parameters");
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const int i) {
  // ********************************************
  throw std::runtime_error("Don't Toggle parameters");
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const std::string &name) {
  // ********************************************
  throw std::runtime_error("Don't Toggle parameters");
}

// ********************************************
bool ParameterHandlerBase::IsParameterFixed(const std::string &name) const {
  // ********************************************
  return IsParameterFixed(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::SetFlatPrior(const int i, const bool eL) {
  // ********************************************
  if (i > params.flatprior.size()) {
    MACH3LOG_INFO(
        "Can't {} for Cov={}/Param={} because size of Covariance = {}",
        __func__, GetName(), i, params.flatprior.size());
    MACH3LOG_ERROR("Fix this in your config file please!");
    throw MaCh3Exception(__FILE__, __LINE__);
  } else {
    if (eL) {
      MACH3LOG_INFO("Setting {} (parameter {}) to flat prior", GetParName(i),
                    i);
    } else {
      // HW :: This is useful
      MACH3LOG_INFO("Setting {} (parameter {}) to non-flat prior",
                    GetParName(i), i);
    }
    params.flatprior[i] = eL;
  }
}

// ********************************************
void ParameterHandlerBase::SetIndivStepScale(
    const std::vector<double> &stepscale) {
  // ********************************************
  if (int(stepscale.size()) != params.stepscale.size()) {
    MACH3LOG_WARN(
        "Stepscale vector not equal to number of parameters. Quitting..");
    MACH3LOG_WARN("Size of argument vector: {}", stepscale.size());
    MACH3LOG_WARN("Expected size: {}", params.stepscale.size());
    return;
  }

  params.stepscale = M3::StdVectorToEigen(stepscale);
  PrintIndivStepScale();
}

// ********************************************
void ParameterHandlerBase::PrintIndivStepScale() const {
  // ********************************************
  MACH3LOG_INFO("============================================================");
  MACH3LOG_INFO("{:<{}} | {:<11}", "Parameter:", PrintLength, "Step scale:");
  for (int iParam = 0; iParam < params.fancyname.size(); iParam++) {
    MACH3LOG_INFO("{:<{}} | {:<11}", params.fancyname[iParam].c_str(),
                  PrintLength, params.stepscale[iParam]);
  }
  MACH3LOG_INFO("============================================================");
}

// ********************************************
void ParameterHandlerBase::ResetIndivStepScale() {
  // ********************************************
  steps.global_scale = 1.0;
  SetIndivStepScale(std::vector<double>(1, params.stepscale.size()));
}

// ********************************************
// HW: Code for throwing from separate throw matrix, needs to be set after init
// to ensure pos-def
void ParameterHandlerBase::SetThrowMatrix(TMatrixDSym *cov) {
  // ********************************************

  if (cov == nullptr) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to {}",
                   __func__);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (params.l_proposal.rows() != cov->GetNrows()) {
    MACH3LOG_ERROR("Matrix given for throw Matrix is not the same size as the "
                   "proposal matrix stored in object!");
    MACH3LOG_ERROR("Stored proposal matrix size: {}", params.l_proposal.rows());
    MACH3LOG_ERROR("Given matrix size: {}", cov->GetNrows());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (use_adaptive && AdaptiveHandler.AdaptionUpdate()) {
    M3::MakeClosestPosDef(cov);
    params.l_proposal = Eigen::LLT<MatrixXd>(M3::ROOTToEigen(cov)).L();
  } else {
    params.l_proposal =
        Eigen::LLT<MatrixXd>(M3::MakePosDef(M3::ROOTToEigen(*cov))).L();
  }
}

// ********************************************
void ParameterHandlerBase::UpdateThrowMatrix(TMatrixDSym *cov) {
  // ********************************************
  SetThrowMatrix(cov);
}

// ********************************************
// HW : Here be adaption
void ParameterHandlerBase::InitialiseAdaption(const YAML::Node &adapt_manager) {
  // ********************************************
  if (PCAObj) {
    MACH3LOG_ERROR("PCA has been enabled and now trying to enable Adaption. "
                   "Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (AdaptiveHandler) {
    MACH3LOG_ERROR("Adaptive Handler has already been initialise can't do it "
                   "again so skipping.");
    return;
  }

  // Now we read the general settings [these SHOULD be common across all
  // matrices!]
  bool success = AdaptiveHandler.InitFromConfig(adapt_manager, matrixName,
                                                &_fCurrVal, &_fError);
  if (!success) {
    return;
  }

  AdaptiveHandler.Print();

  // Next let"s check for external matrices
  // We"re going to grab this info from the YAML manager
  if (!GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"]
                                         [matrixName]["UseExternalMatrix"],
                            false, __FILE__, __LINE__)) {
    MACH3LOG_WARN(
        "Not using external matrix for {}, initialising adaption from scratch",
        matrixName);
    // If we don't have a covariance matrix to start from for adaptive tune we
    // need to make one!
    use_adaptive = true;
    AdaptiveHandler.CheckMatrixValidityForAdaption(GetCovMatrix());
    AdaptiveHandler.CreateNewAdaptiveCovariance();
    return;
  }

  // Finally, we accept that we want to read the matrix from a file!
  auto external_file_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][matrixName]
                   ["ExternalMatrixFileName"],
      "", __FILE__, __LINE__);
  auto external_matrix_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][matrixName]
                   ["ExternalMatrixName"],
      "", __FILE__, __LINE__);
  auto external_mean_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][matrixName]
                   ["ExternalMeansName"],
      "", __FILE__, __LINE__);

  AdaptiveHandler.SetThrowMatrixFromFile(external_file_name,
                                         external_matrix_name,
                                         external_mean_name, use_adaptive);
  SetThrowMatrix(AdaptiveHandler.GetAdaptiveCovariance());

  ResetIndivStepScale();

  MACH3LOG_INFO("Successfully Set External Throw Matrix Stored in {}",
                external_file_name);
}

// ********************************************
// Truely adaptive MCMC!
void ParameterHandlerBase::UpdateAdaptiveCovariance() {
  // ********************************************
  // Updates adaptive matrix
  // First we update the total means

  // Skip this if we're at a large number of steps
  if (AdaptiveHandler.SkipAdaption()) {
    AdaptiveHandler.IncrementNSteps();
    return;
  }

  /// Need to adjust the scale every step
  if (AdaptiveHandler.GetUseRobbinsMonro()) {
    bool verbose = false;
#ifdef DEBUG
    verbose = true;
#endif
    AdaptiveHandler.UpdateRobbinsMonroScale();
    SetStepScale(AdaptiveHandler.GetAdaptionScale(), verbose);
  }

  // Call main adaption function
  AdaptiveHandler.UpdateAdaptiveCovariance();

  // Set scales to 1 * optimal scale
  if (AdaptiveHandler.IndivStepScaleAdapt()) {
    ResetIndivStepScale();
    SetStepScale(AdaptiveHandler.GetAdaptionScale());
  }

  if (AdaptiveHandler.UpdateMatrixAdapt()) {
    UpdateThrowMatrix(
        AdaptiveHandler.GetAdaptiveCovariance()); // Now we update and continue!
    // Also Save the adaptive to file
    AdaptiveHandler.SaveAdaptiveToFile(AdaptiveHandler.GetOutFileName(),
                                       GetName());
  }

  AdaptiveHandler.IncrementNSteps();
}

// ********************************************
// KS: After step scale, prefit etc. value were modified save this modified
// config.
void ParameterHandlerBase::SaveUpdatedMatrixConfig(
    std::string const &filename) {
  // ********************************************
  if (!config.YAMLDoc) {
    MACH3LOG_CRITICAL("Yaml node hasn't been initialised for matrix {}, "
                      "something is not right",
                      matrixName);
    MACH3LOG_CRITICAL("I am not throwing error but should be investigated");
    return;
  }

  YAML::Node copyNode = config.YAMLDoc;
  int i = 0;

  for (YAML::Node param : copyNode["Systematics"]) {
    // KS: Feel free to update it, if you need updated prefit value etc
    param["Systematic"]["StepScale"]["MCMC"] =
        MaCh3Utils::FormatDouble(params.stepscale[i], 4);
    i++;
  }
  // Save the modified node to a file
  std::ofstream fout(filename);
  fout << copyNode;
  fout.close();
}

// ********************************************
bool ParameterHandlerBase::AppliesToSample(
    const int SystIndex, const std::string &SampleName) const {
  // ********************************************
  // Empty means apply to all
  if (!params.samples[SystIndex].size()) {
    return true;
  }

  // Make a copy and to lower case to not be case sensitive
  std::string SampleNameCopy = SampleName;
  std::transform(SampleNameCopy.begin(), SampleNameCopy.end(),
                 SampleNameCopy.begin(), ::tolower);

  // Check for unsupported wildcards in SampleNameCopy
  if (SampleNameCopy.find('*') != std::string::npos) {
    MACH3LOG_ERROR("Wildcards ('*') are not supported in sample name: '{}'",
                   SampleName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  bool Applies = false;

  for (size_t i = 0; i < _fSampleNames[SystIndex].size(); i++) {
    // Convert to low case to not be case sensitive
    std::string pattern = _fSampleNames[SystIndex][i];
    std::transform(pattern.begin(), pattern.end(), pattern.begin(), ::tolower);

    // Replace '*' in the pattern with '.*' for regex matching
    std::string regexPattern =
        "^" + std::regex_replace(pattern, std::regex("\\*"), ".*") + "$";
    try {
      std::regex regex(regexPattern);
      if (std::regex_match(SampleNameCopy, regex)) {
        Applies = true;
        break;
      }
    } catch (const std::regex_error &e) {
      // Handle regex error (for invalid patterns)
      MACH3LOG_ERROR("Regex error: {}", e.what());
    }
  }
  return Applies;
}

// ********************************************
// Set proposed parameter values vector to be base on tune values
void ParameterHandlerBase::SetTune(const std::string &TuneName) {
  // ********************************************
  SetParameters(Tunes->GetTune(TuneName));
}

// *************************************
/// @brief Matches branches in a TTree to parameters in a systematic handler.
///
/// @param PosteriorFile Pointer to the ROOT TTree from MaCh3 fit.
/// @param Systematic Pointer to the systematic parameter handler.
/// @param[out] BranchValues Vector to store the values of the branches (resized
/// inside).
/// @param[out] BranchNames Vector to store the names of the branches (resized
/// inside).
///
/// @throws MaCh3Exception if any parameter branch is uninitialized.
void ParameterHandlerBase::MatchMaCh3OutputBranches(
    TTree *PosteriorFile, std::vector<double> &BranchValues,
    std::vector<std::string> &BranchNames) {
  // *************************************
  BranchValues.resize(GetNumParams());
  BranchNames.resize(GetNumParams());

  for (int i = 0; i < GetNumParams(); ++i) {
    BranchNames[i] = GetParName(i);
    if (!PosteriorFile->GetBranch(BranchNames[i].c_str())) {
      MACH3LOG_ERROR("Branch '{}' does not exist in the TTree!",
                     BranchNames[i]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    PosteriorFile->SetBranchStatus(BranchNames[i].c_str(), true);
    PosteriorFile->SetBranchAddress(BranchNames[i].c_str(), &BranchValues[i]);
  }
}
