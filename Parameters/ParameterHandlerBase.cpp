#include "Parameters/ParameterHandlerBase.h"
#include <regex>

ParameterHandlerBase::ParameterHandlerBase() {
  pca.enabled = false;
  special_proposal.enabled = false;
  rng.gaus = std::normal_distribution<double>(0);
  rng.unif = std::uniform_real_distribution<double>(0, 1);
  settings.use_adaptive = false;
  settings.PrintLength = 35;
}

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string const &name,
                                           std::string const &file)
    : ParameterHandlerBase() {
  // ********************************************

  SetName(name);
  config.inputFiles = {
      file,
  };

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

  params.covariance = M3::MakeMatrixPosDef(M3::ROOTToEigen(*CovMat));

  MACH3LOG_INFO("Created covariance matrix named: {}", GetName());
  MACH3LOG_INFO("from file: {}", file);
}

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(
    const std::vector<std::string> &YAMLFiles, std::string const &name)
    : ParameterHandlerBase() {
  // ********************************************

  SetName(name);
  config.inputFiles = YAMLFiles;

  MACH3LOG_INFO("Constructing instance of ParameterHandler using: ");
  for (auto const &yf : config.inputFiles) {
    MACH3LOG_INFO("  {}", yf);
  }
  MACH3LOG_INFO("as an input.");

  config.YAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for (auto const &yfn : config.inputFiles) {
    for (const auto &node : M3OpenConfig(yfn)["Systematics"]) {
      config.YAMLDoc["Systematics"].push_back(node);
    }
  }

  std::map<int, std::map<std::string, double>> parameter_correlations;

  std::vector<ParamInfo> param_infos;

  // ETA - read in the systematics.
  for (auto const &node : config.YAMLDoc["Systematics"]) {
    int param_id = int(param_infos.size());

    auto const &pardef = node["Systematic"];

    auto fancy_name =
        Get<std::string>(pardef["Names"]["FancyName"], __FILE__, __LINE__);
    auto prefit = Get<double>(pardef["ParameterValues"]["PreFitValue"],
                              __FILE__, __LINE__);

    auto error = Get<double>(pardef["Error"], __FILE__, __LINE__);

    if (error <= 0) {
      MACH3LOG_ERROR("Error for param {}({}) is negative and equal to {}",
                     fancy_name, param_id, error);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    auto stepscale =
        Get<double>(pardef["StepScale"]["MCMC"], __FILE__, __LINE__);

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
      EnableSpecialProposal(pardef["SpecialProposal"], param_id);
    }

    if (pardef["Correlations"]) {
      for (auto const &corrn : pardef["Correlations"]) {
        for (auto const &corr : corrn) {
          parameter_correlations[param_id][corr.first.as<std::string>()] =
              corr.second.as<double>();
        }
      }
    }

    param_infos.emplace_back(ParamInfo{fancy_name,
                                       fancy_name,
                                       prefit,
                                       error,
                                       stepscale,
                                       {bounds[0], bounds[1]},
                                       flatprior,
                                       !fixed,
                                       samplenames});
  }

  AddParameters(param_infos);

  for (auto const &[paramid, correlations] : parameter_correlations) {
    SetParameterAllCorrelations(paramid, correlations);
  }

  params.covariance = M3::MakeMatrixPosDef(params.covariance);

  Tunes = ParameterTunes(config.YAMLDoc["Systematics"]);

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for (const auto &file : config.inputFiles) {
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total",
                params.prefit.size());
  MACH3LOG_INFO("----------------");
}

// ********************************************
void ParameterHandlerBase::ConstructPCA(const double threshold, int first,
                                        int last) {
  // ********************************************
  if (settings.use_adaptive) {
    MACH3LOG_ERROR("Adaption has been enabled and now trying to enable PCA. "
                   "Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  block_indices[0] = first;
  block_indices[1] = last;

  int nphysics_parameters = (last - first) + 1;

  double block_trace =
      params.covariance.block(first, first, last - first + 1, last - first + 1)
          .trace();

  // a covariance is real symmetric, so self adjoint
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EigenSolver(
      params.covariance.block(first, first, nphysics_parameters,
                              nphysics_parameters));
  Eigen::VectorXd eigen_val = EigenSolver.eigenvalues();
  Eigen::MatrixXd eigen_vect = EigenSolver.eigenvectors();
  std::vector<std::pair<int, double>> evals;
  for (int i = 0; i < eigen_val.size(); ++i) {
    evals.emplace_back(i, eigen_val[i]);
  }
  std::sort(
      evals.begin(), evals.end();
      [](std::pair<int, double> const &l, std::pair<int, double> const &r) {
        return l.second < r.second;
      });

  double evsum = 0;
  int northo_parameters = 0;

  for (auto &[i, ev] : evals) {
    if (std::fabs(ev) > threshold) {
      evsum += std::fabs(ev);
      northo_parameters++;
    } else {
      break;
    }
  }

  //rows of this matrix correspond to eigenvectors of the input matrix
  //  we can go from parameters defined in the orthogonal basis back to the
  // 'physics' parameters by RowVectOfOrthoParamVals * pca.ortho_to_physics
  pca.ortho_to_physics =
      Eigen::MatrixXd::Zero(evals.size(), nphysics_parameters);

  for (size_t i = 0; i < northo_parameters; ++i) {
    for (size_t j = 0; j < nphysics_parameters; ++j) {
      pca.ortho_to_physics.row(i)[j] =
          eigen_vect.col(evals[i].first)[evals[j].first] *
          std::sqrt(evals[i].second);
    }
  }

  MACH3LOG_INFO("Threshold of {} on eigen values, kept {}/{} for a total "
                "variance of {}/{}",
                threshold, nkept, evals.size(), evsum, block_trace);
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
    auto bounds = Get<std::pair<double, double>>(param["CircularBounds"],
                                                 __FILE__, __LINE__);
    // KS: Make sure circular bounds are within physical bounds. If we are
    // outside of physics bound MCMC will never explore such phase space region
    if (bounds.first < params.lowbound[Index] ||
        bounds.second > params.upbound[Index]) {
      MACH3LOG_ERROR("Circular bounds [{}, {}] for parameter {} exceed "
                     "physical bounds [{}, {}]",
                     bounds.first, bounds.second, GetParFancyName(Index),
                     params.lowbound[Index], params.upbound[Index]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    special_proposal.circ_bounds.emplace_back(
        std::make_tuple(Index, bounds.first, bounds.second));

    MACH3LOG_INFO(
        "Enabling CircularBounds for parameter {} with range [{}, {}]",
        GetParFancyName(Index), bounds.first, bounds.second);
  }

  if (FlipEnabled) {
    special_proposal.flips.emplace_back(std::make_pair(
        Index, Get<double>(param["FlipParameter"], __FILE__, __LINE__)));
    MACH3LOG_INFO("Enabling Flipping for parameter {} with value {}",
                  GetParFancyName(Index), special_proposal.flips.back().second);
  }

  if (CircEnabled && FlipEnabled) {

    const double fp = special_proposal.flips.back().second;
    const double low = std::get<1>(special_proposal.circ_bounds.back());
    const double high = std::get<2>(special_proposal.circ_bounds.back());

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
                     special_proposal.flips.back().second,
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
        cov->GetNrows(), params.prefit.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  params.covariance = M3::ROOTToEigen(*cov);
  params.inv_covariance = params.covariance.inverse();

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

  throws.values = steps.l_proposal * throws.random_vector;

  steps.proposed = params.prefit + throws.values;

  for (int i = 0; i < throws.values.size(); ++i) {

    int tries = 0;
    while ((steps.proposed[i] < params.lowbound[i]) ||
           (steps.proposed[i] >= params.upbound[i])) {
      steps.proposed[i] =
          params.prefit[i] + (steps.l_proposal.col(i) * rng.gaus(rng.e1))[i];
      tries++;
      if (tries > 10000) {
        // KS: Since we are multithreading there is danger that those messages
        // will be all over the place, small price to pay for faster code
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", tries,
                      i);
        MACH3LOG_WARN("Matrix: {}", settings.name);
        MACH3LOG_WARN("Param: {}", params.name[i]);
        MACH3LOG_WARN("Setting steps.proposed:  {} to {}", steps.proposed[i],
                      params.prefit[i]);
        MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
        steps.proposed[i] = params.prefit[i];
      }
    }

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
  // likelihood Want to change the steps.proposed, steps.current, params.prefit
  // params.prefit and the others will already be set

  for (int i = 0; i < params.prefit.size(); ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (IsParameterFixed(i)) {
      continue;
    }
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    double throwrange = std::min(params.upbound[i] - params.lowbound[i],
                                 std::sqrt(params.covariance(i, i)));

    steps.proposed[i] = params.prefit[i] + rng.gaus(rng.e1) * throwrange;
    // Try again if we the initial parameter proposal falls outside of the range
    // of the parameter
    int nthrows = 0;
    while ((steps.proposed[i] > params.upbound[i]) ||
           (steps.proposed[i] < params.lowbound[i])) {
      if (nthrows > 1000) {
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed",
                      nthrows, i);
        MACH3LOG_WARN("Matrix: {}", settings.name);
        MACH3LOG_WARN("Param: {}", params.name[i]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      steps.proposed[i] = params.prefit[i] + rng.gaus(rng.e1) * throwrange;
      nthrows++;
    }
    MACH3LOG_INFO("Setting current step in {} param {} = {} from {}",
                  settings.name, i, steps.proposed[i], steps.current[i]);
    steps.current[i] = steps.proposed[i];
  }
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
  if (!special_proposal.enabled) {
    return;
  }

  SpecialStepProposal();
}

// ************************************************
void ParameterHandlerBase::SpecialStepProposal() {
  // ************************************************
  /// @warning KS: Following Asher comment we do "Step->Circular Bounds->Flip"

  // HW It should now automatically set dcp to be with [-pi, pi]
  for (auto const &[idx, low, up] : special_proposal.circ_bounds) {
    if (IsParameterFixed(idx)) {
      continue;
    }
    if (steps.proposed[idx] > up) {
      steps.proposed[idx] = low + std::fmod(steps.proposed[idx] - up, up - low);
    } else if (steps.proposed[idx] < low) {
      steps.proposed[idx] = up - std::fmod(low - steps.proposed[idx], up - low);
    }
  }

  // Okay now we've done the standard steps, we can add in our nice flips
  // hierarchy flip first
  for (auto const &[idx, pivot] : special_proposal.flips) {
    if (IsParameterFixed(idx)) {
      continue;
    }
    if (rng.unif(rng.e1) < 0.5) {
      steps.proposed[idx] = 2 * pivot - steps.proposed[idx];
    }
  }
}

// ************************************************
// "Randomize" the parameters in the covariance class for the proposed step
// Used the proposal kernel and the current parameter value to set proposed step
// Also get a new random number for the randParams
void ParameterHandlerBase::Randomize() _noexcept_ {
  // ************************************************
  throws.random_vector =
      Eigen::VectorXd::NullaryExpr(params.prefit.size(), [&](int i) {
        return params.isfree[i] ? rng.gaus(rng.e1) : 0;
      });
}

// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its
// current value + some correlated throw
void ParameterHandlerBase::CorrelateSteps() _noexcept_ {
  // ************************************************

  Randomize();

  // array makes the multiplaction by steps.scale component-wise rather than
  // vector-ie
  steps.proposed =
      steps.current + ((steps.l_proposal * throws.random_vector).array() *
                       steps.scale.array() * steps.global_scale)
                          .matrix();
}
// ********************************************
// Update so that current step becomes the previously proposed step
void ParameterHandlerBase::AcceptStep() _noexcept_ {
  // ********************************************
  steps.current = steps.proposed;

  if (settings.use_adaptive) {
    AdaptiveHandler.IncrementAcceptedSteps();
  }
}

// ********************************************
// Function to print the prior values
void ParameterHandlerBase::PrintNominal() const {
  // ********************************************
  MACH3LOG_INFO("Prior values for {} ParameterHandler:", GetName());
  for (int i = 0; i < int(params.name.size()); i++) {
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
                  params.prefit[i], steps.current[i], steps.proposed[i]);
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
      params.flatprior.select(Eigen::VectorXd::Zero(params.prefit.size()),
                              steps.proposed - params.prefit);
  return diff.transpose() * params.inv_covariance * diff;
}

// ********************************************
int ParameterHandlerBase::CheckBounds() const _noexcept_ {
  // ********************************************
  return int(((steps.proposed.array() < params.lowbound.array()) ||
              (steps.proposed.array() >= params.upbound.array()))
                 .count());
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
    if (pars.size() != size_t(steps.proposed.size())) {
      MACH3LOG_ERROR("Parameter arrays of incompatible size! Not changing "
                     "parameters! {} has size {} but was expecting {}",
                     settings.name, pars.size(), steps.proposed.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    for (int i = 0; i < int(pars.size()); i++) {
      // Make sure that you are actually passing a number to set the parameter
      // to
      if (std::isnan(pars[i])) {
        MACH3LOG_ERROR("Trying to set parameter value to a nan for parameter "
                       "{} in matrix {}. This will not go well!",
                       GetParName(i), settings.name);
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
  for (int i = 0; i < int(params.name.size()); ++i) {
    tree.Branch(params.name[i].c_str(), steps.current.data() + i,
                Form("%s/D", params.name[i].c_str()));
  }
  // // When running PCA, also save PCA parameters
  // if (pca.enabled) {
  //   PCAObj.SetBranches(tree, SaveProposal, _fNames);
  // }
  if (SaveProposal) {
    // loop over parameters and set a branch
    for (int i = 0; i < int(params.name.size()); ++i) {
      tree.Branch(Form("%s_Prop", params.name[i].c_str()),
                  steps.proposed.data() + i,
                  Form("%s_Prop/D", params.name[i].c_str()));
    }
  }
  if (settings.use_adaptive && AdaptiveHandler.GetUseRobbinsMonro()) {
    tree.Branch(Form("GlobalStepScale_%s", GetName().c_str()),
                &steps.global_scale,
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
    const double SuggestedScale = 2.38 / std::sqrt(steps.scale.size());
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
  return int(std::distance(params.name.begin(), idx));
}

// ********************************************
void ParameterHandlerBase::SetFixAllParameters() {
  // ********************************************
  params.isfree = Eigen::VectorXi::Zero(params.prefit.size());
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const int i) {
  // ********************************************
  params.isfree[i] = 0;
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const std::string &name) {
  // ********************************************
  SetFixParameter(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::SetFreeAllParameters() {
  // ********************************************
  params.isfree = Eigen::VectorXi::Ones(params.prefit.size());
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const int i) {
  // ********************************************
  params.isfree[i] = 1;
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
void ParameterHandlerBase::ToggleFixParameter(const int) {
  // ********************************************
  throw std::runtime_error("Don't Toggle parameters");
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const std::string &) {
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
  if (int(steps.scale.size()) != steps.scale.size()) {
    MACH3LOG_WARN(
        "Stepscale vector not equal to number of parameters. Quitting..");
    MACH3LOG_WARN("Size of argument vector: {}", steps.scale.size());
    MACH3LOG_WARN("Expected size: {}", steps.scale.size());
    return;
  }

  steps.scale = M3::StdVectorToEigen(stepscale);
  PrintIndivStepScale();
}

// ********************************************
void ParameterHandlerBase::PrintIndivStepScale() const {
  // ********************************************
  MACH3LOG_INFO("============================================================");
  MACH3LOG_INFO("{:<{}} | {:<11}", "Parameter:", settings.PrintLength,
                "Step scale:");
  for (int iParam = 0; iParam < int(params.fancy_name.size()); iParam++) {
    MACH3LOG_INFO("{:<{}} | {:<11}", params.fancy_name[iParam].c_str(),
                  settings.PrintLength, steps.scale[iParam]);
  }
  MACH3LOG_INFO("============================================================");
}

// ********************************************
void ParameterHandlerBase::ResetIndivStepScale() {
  // ********************************************
  steps.global_scale = 1.0;
  SetIndivStepScale(std::vector<double>(steps.scale.size(), 1.0));
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

  if (steps.l_proposal.rows() != cov->GetNrows()) {
    MACH3LOG_ERROR("Matrix given for throw Matrix is not the same size as the "
                   "proposal matrix stored in object!");
    MACH3LOG_ERROR("Stored proposal matrix size: {}", steps.l_proposal.rows());
    MACH3LOG_ERROR("Given matrix size: {}", cov->GetNrows());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (settings.use_adaptive && AdaptiveHandler.AdaptionUpdate()) {
    M3::MakeMatrixClosestPosDef(cov);
    Eigen::LLT<Eigen::MatrixXd> lltproposal(M3::ROOTToEigen(*cov));
    if (lltproposal.info() != Eigen::Success) {
      MACH3LOG_ERROR("Failed to LLT decompose proposal matrix");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    steps.l_proposal = lltproposal.matrixL();
  } else {
    Eigen::LLT<Eigen::MatrixXd> lltproposal(
        M3::MakeMatrixPosDef(M3::ROOTToEigen(*cov)));
    if (lltproposal.info() != Eigen::Success) {
      MACH3LOG_ERROR("Failed to LLT decompose proposal matrix");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    steps.l_proposal = lltproposal.matrixL();
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
  if (pca.enabled) {
    MACH3LOG_ERROR("PCA has been enabled and now trying to enable Adaption. "
                   "Right now both configuration don't work with each other");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (settings.use_adaptive) {
    MACH3LOG_ERROR("Adaptive Handler has already been initialise can't do it "
                   "again so skipping.");
    return;
  }

  // Now we read the general settings [these SHOULD be common across all
  // matrices!]
  // bool success = AdaptiveHandler.InitFromConfig(adapt_manager, settings.name,
  //                                               &steps.current,
  //                                               &params.error);
  // if (!success) {
  //   return;
  // }

  AdaptiveHandler.Print();

  // Next let"s check for external matrices
  // We"re going to grab this info from the YAML manager
  if (!GetFromManager<bool>(adapt_manager["AdaptionOptions"]["Covariance"]
                                         [settings.name]["UseExternalMatrix"],
                            false, __FILE__, __LINE__)) {
    MACH3LOG_WARN(
        "Not using external matrix for {}, initialising adaption from scratch",
        settings.name);
    // If we don't have a covariance matrix to start from for adaptive tune we
    // need to make one!
    settings.use_adaptive = true;
    AdaptiveHandler.CheckMatrixValidityForAdaption(GetCovMatrix());
    AdaptiveHandler.CreateNewAdaptiveCovariance();
    return;
  }

  // Finally, we accept that we want to read the matrix from a file!
  auto external_file_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][settings.name]
                   ["ExternalMatrixFileName"],
      "", __FILE__, __LINE__);
  auto external_matrix_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][settings.name]
                   ["Externalsettings.name"],
      "", __FILE__, __LINE__);
  auto external_mean_name = GetFromManager<std::string>(
      adapt_manager["AdaptionOptions"]["Covariance"][settings.name]
                   ["ExternalMeansName"],
      "", __FILE__, __LINE__);

  AdaptiveHandler.SetThrowMatrixFromFile(
      external_file_name, external_matrix_name, external_mean_name,
      settings.use_adaptive);
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
                      settings.name);
    MACH3LOG_CRITICAL("I am not throwing error but should be investigated");
    return;
  }

  YAML::Node copyNode = config.YAMLDoc;
  int i = 0;

  for (YAML::Node param : copyNode["Systematics"]) {
    // KS: Feel free to update it, if you need updated prefit value etc
    param["Systematic"]["StepScale"]["MCMC"] =
        MaCh3Utils::FormatDouble(steps.scale[i], 4);
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

  for (size_t i = 0; i < params.samples[SystIndex].size(); i++) {
    // Convert to low case to not be case sensitive
    std::string pattern = params.samples[SystIndex][i];
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
  SetParameters(Tunes.GetTune(TuneName));
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
