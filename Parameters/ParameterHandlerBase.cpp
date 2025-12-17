#include "Parameters/ParameterHandlerBase.h"

#include <regex>

// ********************************************
ParameterHandlerBase::ParameterHandlerBase(std::string name, std::string file,
                                           double, int, int)
    : inputFile(file) {
  // ********************************************
  MACH3LOG_DEBUG("Constructing instance of ParameterHandler");

  parlist = ParameterList::MakeFromTMatrix(name, file);
  current = parlist.params.prefit;
  proposed = parlist.params.prefit;

  covMatrix = new TMatrixDSym(int(parlist.params.covariance.rows()),
                              int(parlist.params.covariance.cols()));
  M3::EigenToROOT(parlist.params.covariance, *covMatrix);
  invCovMatrix = new TMatrixDSym(int(parlist.params.covariance.rows()),
                                 int(parlist.params.covariance.cols()));
  M3::EigenToROOT(parlist.params.inv_covariance, *invCovMatrix);

  proposer = parlist.MakeProposer();

  throwMatrix = new TMatrixDSym(int(proposer.params.proposal.rows()),
                                int(proposer.params.proposal.cols()));
  M3::EigenToROOT(proposer.params.proposal, *throwMatrix);
}
// ********************************************
ParameterHandlerBase::ParameterHandlerBase(
    const std::vector<std::string> &YAMLFile, std::string name,
    double threshold, int FirstPCA, int LastPCA)
    : matrixName(name), inputFile(YAMLFile[0].c_str()) {
  // ********************************************
  MACH3LOG_INFO("Constructing instance of ParameterHandler using");

  for (unsigned int i = 0; i < YAMLFile.size(); i++) {
    MACH3LOG_INFO("{}", YAMLFile[i]);
  }
  MACH3LOG_INFO("as an input");

  bool pca = true;
  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("Principal component analysis but given the threshold for "
                  "the principal components to be less than 0, or greater than "
                  "(or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }

  parlist = ParameterList::MakeFromYAML(YAMLFile);
  current = parlist.params.prefit;
  proposed = parlist.params.prefit;

  covMatrix = new TMatrixDSym(int(parlist.params.covariance.rows()),
                              int(parlist.params.covariance.cols()));
  M3::EigenToROOT(parlist.params.covariance, *covMatrix);
  invCovMatrix = new TMatrixDSym(int(parlist.params.covariance.rows()),
                                 int(parlist.params.covariance.cols()));
  M3::EigenToROOT(parlist.params.inv_covariance, *invCovMatrix);

  if (pca) {
    ConstructPCA(threshold, FirstPCA, LastPCA);
  }

  proposer = parlist.MakeProposer();

  throwMatrix = new TMatrixDSym(int(proposer.params.proposal.rows()),
                                int(proposer.params.proposal.cols()));
  M3::EigenToROOT(proposer.params.proposal, *throwMatrix);
}

// ********************************************
void ParameterHandlerBase::ConstructPCA(const double threshold,
                                        int FirstPCAdpar, int LastPCAdpar) {
  // ********************************************
  parlist.ConstructTruncatedPCA(threshold, FirstPCAdpar, LastPCAdpar);
}

// ********************************************
// Set the covariance matrix for this class
void ParameterHandlerBase::SetCovMatrix(TMatrixDSym *) {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::SetCovMatrix Removed");
}

// ********************************************
// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void ParameterHandlerBase::SetPar(int i, double val) {
  // ********************************************

  current[i] = val;
  proposed[i] = val;

  MACH3LOG_INFO("Overriding {}: ", GetParName(i));
  MACH3LOG_INFO("_fPropVal ({}), _fCurrVal ({}), _fPreFitValue ({}) to ({})",
                proposed[i], current[i], parlist.params.prefit[i], val);

  parlist.RotateSystematicParameterValuesToPCBasis(current,
                                                   proposer.params.current);
  parlist.RotateSystematicParameterValuesToPCBasis(proposed,
                                                   proposer.params.proposed);
}

// ********************************************
std::vector<double> ParameterHandlerBase::GetProposed() const {
  // ********************************************
  return M3::EigenToStdVector(proposer.params.proposed);
}

// *************************************
// Throw the parameters according to the covariance matrix
// This shouldn't be used in MCMC code ase it can break Detailed Balance;
void ParameterHandlerBase::ThrowParameters() {
  // *************************************

  int tries = 0;
  proposer.Propose();
  while (CheckBounds()) {
    proposer.Propose();

    tries++;
    if (tries > 10000) {
      // KS: Since we are multithreading there is danger that those messages
      // will be all over the place, small price to pay for faster code
      MACH3LOG_ERROR(
          "Tried {} times to throw parameters within bounds but failed.",
          tries);
      MACH3LOG_ERROR("Matrix: {}", GetName());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  proposer.Accept();

  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.current,
                                                   current);
  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.proposed,
                                                   proposed);
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

  for (int i = 0; i < parlist.NumSystematicBasisParameters(); ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (IsParameterFixed(i)) {
      continue;
    }
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    const double paramrange =
        parlist.params.upbound[i] - parlist.params.lowbound[i];
    const double sigma = parlist.params.error[i];
    double throwrange = sigma;
    if (paramrange < sigma) {
      throwrange = paramrange;
    }

    current[i] = parlist.params.prefit[i] +
                 proposer.rng.gaus(proposer.rng.e1) * throwrange;
    // Try again if we the initial parameter proposal falls outside of the range
    // of the parameter
    int throws = 0;
    while (current[i] > parlist.params.upbound[i] ||
           current[i] < parlist.params.lowbound[i]) {
      if (throws > 1000) {
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws,
                      i);
        MACH3LOG_WARN("Matrix: {}", GetName());
        MACH3LOG_WARN("Param: {}", GetParName(i));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      current[i] = parlist.params.prefit[i] +
                   proposer.rng.gaus(proposer.rng.e1) * throwrange;
      throws++;
    }
    MACH3LOG_INFO("Setting current step in {} param {}, {} = {} from {}",
                  GetName(), i, GetParName(i), current[i], current[i]);
  }
  parlist.RotateSystematicParameterValuesToPCBasis(current,
                                                   proposer.params.proposed);
  proposer.Accept();
  proposed = current;
}

// *************************************
// Set a single parameter
void ParameterHandlerBase::SetSingleParameter(const int parNo,
                                              const double parVal) {
  // *************************************
  SetPar(parNo, parVal);
}

// ********************************************
void ParameterHandlerBase::SetParCurrProp(const int parNo,
                                          const double parVal) {
  // ********************************************
  SetPar(parNo, parVal);
}

void ParameterHandlerBase::SetParProp(const int i, const double val) {
  proposed[i] = val;
  parlist.RotateSystematicParameterValuesToPCBasis(proposed,
                                                   proposer.params.proposed);
}

void ParameterHandlerBase::SetRandomThrow(const int, const double) {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::SetRandomThrow Removed");
}
double ParameterHandlerBase::GetRandomThrow(const int) const {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::GetRandomThrow Removed");
}

// ************************************************
// Propose a step for the set of systematics parameters this covariance class
// holds
void ParameterHandlerBase::ProposeStep() {
  // ************************************************
  proposer.Propose();
  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.proposed,
                                                   proposed);
}

// ********************************************
// Update so that current step becomes the previously proposed step
void ParameterHandlerBase::AcceptStep() _noexcept_ {
  // ********************************************
  proposer.Accept();
  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.current,
                                                   current);
}

// ********************************************
// Throw the proposed parameter by mag sigma
// Should really just have the user specify this throw by having argument double
void ParameterHandlerBase::ThrowParProp(const double mag) {
  // ********************************************
  Eigen::VectorXd new_proposed =
      ((proposer.Propose() - proposer.params.current) * mag) +
      proposer.params.current;
  // probably fine but avoid ailiasing for now
  proposer.params.proposed = new_proposed;
  parlist.RotatePCParameterValuesToSystematicBasis(new_proposed, proposed);
}
// ********************************************
// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void ParameterHandlerBase::ThrowParCurr(const double mag) {
  // ********************************************

  Eigen::VectorXd curr_proposed = proposer.params.proposed;
  Eigen::VectorXd new_curr =
      ((proposer.Propose() - proposer.params.current) * mag);
  proposer.params.current = new_curr;
  // probably fine but avoid ailiasing for now
  proposer.params.proposed = curr_proposed;

  parlist.RotatePCParameterValuesToSystematicBasis(new_curr, current);
}
// ********************************************
// Function to print the prior values
void ParameterHandlerBase::PrintNominal() const {
  // ********************************************
  MACH3LOG_INFO("Prior values for {} ParameterHandler:", GetName());
  for (int i = 0; i < parlist.NumSystematicBasisParameters(); i++) {
    MACH3LOG_INFO("    {}   {} ", GetParFancyName(i), GetParInit(i));
  }
}

// ********************************************
// Function to print the prior, current and proposed values
void ParameterHandlerBase::PrintNominalCurrProp() const {
  // ********************************************

  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.current,
                                                   current);

  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.proposed,
                                                   proposed);

  MACH3LOG_INFO("Printing parameters for {}", GetName());
  MACH3LOG_INFO("{:<30} {:<10} {:<10} {:<10}", "Name", "Prior", "Current",
                "Proposed");
  for (int i = 0; i < parlist.NumSystematicBasisParameters(); ++i) {
    MACH3LOG_INFO("{:<30} {:<10.2f} {:<10.2f} {:<10.2f}", GetParFancyName(i),
                  parlist.params.prefit[i], proposed[i], current[i]);
  }
}

// ********************************************
// Get the likelihood in the case where we want to include priors on the
// parameters _fFlatPrior stores if we want to evaluate the likelihood for the
// given parameter
//                    true = don't evaluate likelihood (so run without a prior)
//                    false = evaluate likelihood (so run with a prior)
double ParameterHandlerBase::CalcLikelihood() _noexcept_ {
  // ********************************************
  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.proposed,
                                                   proposed);

  return parlist.Chi2(proposed);
}

// ********************************************
int ParameterHandlerBase::CheckBounds() const _noexcept_ {
  // ********************************************

  parlist.RotatePCParameterValuesToSystematicBasis(proposer.params.proposed,
                                                   proposed);

  return int(((proposed.array() < parlist.params.lowbound.array()) ||
              (proposed.array() >= parlist.params.upbound.array()))
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

double ParameterHandlerBase::GetInvCovMatrix(const int, const int) const {

  throw std::runtime_error("ParameterHandlerBase::GetInvCovMatrix Removed");
}

double ParameterHandlerBase::GetCorrThrows(const int) const {

  throw std::runtime_error("ParameterHandlerBase::GetCorrThrows Removed");
}

bool ParameterHandlerBase::GetFlatPrior(const int i) const {
  int propidx = parlist.SystematicParameterIndexToPCIndex(i);
  if (propidx != ParameterList::ParameterInPCABlock) {
    return parlist.params.flatprior[i];
  } else {
    MACH3LOG_ERROR("You are trying Fix parameter {}, {}, which has been PCA'd.",
                   i, GetParName(i));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ********************************************
// Sets the proposed parameters to the prior values
void ParameterHandlerBase::SetParameters(const std::vector<double> &pars) {
  // ********************************************
  // If empty, set the proposed to prior
  if (pars.empty()) {
    // For xsec this means setting to the prior (because prior is the prior)
    parlist.RotateSystematicParameterValuesToPCBasis(parlist.params.prefit,
                                                     proposer.params.proposed);
    // If not empty, set the parameters to the specified
  } else {
    if (int(pars.size()) != parlist.NumSystematicBasisParameters()) {
      MACH3LOG_ERROR("Parameter arrays of incompatible size! Not changing "
                     "parameters! {} has size {} but was expecting {}",
                     GetName(), pars.size(),
                     parlist.NumSystematicBasisParameters());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    parlist.RotateSystematicParameterValuesToPCBasis(M3::StdVectorToEigen(pars),
                                                     proposer.params.proposed);
  }
}

// ********************************************
void ParameterHandlerBase::SetBranches(TTree &tree, bool SaveProposal) {
  // ********************************************
  // loop over parameters and set a branch
  for (int i = 0; i < parlist.NumSystematicBasisParameters(); ++i) {
    tree.Branch(GetParName(i).c_str(), &current.data()[i],
                Form("%s/D", GetParName(i).c_str()));
  }
  if (SaveProposal) {
    // loop over parameters and set a branch
    for (int i = 0; i < parlist.NumSystematicBasisParameters(); ++i) {
      tree.Branch(Form("%s_Prop", GetParName(i).c_str()), &proposed.data()[i],
                  Form("%s_Prop/D", GetParName(i).c_str()));
    }
  }
}

// ********************************************
void ParameterHandlerBase::SetStepScale(const double scale,
                                        const bool verbose) {
  // ********************************************
  if (scale <= 0) {
    MACH3LOG_ERROR(
        "You are trying so set StepScale to 0 or negative this will not work");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (verbose) {
    MACH3LOG_INFO("{} setStepScale() = {}", GetName(), scale);
    const double SuggestedScale =
        2.38 / std::sqrt(parlist.NumPCBasisParameters());
    if (std::fabs(scale - SuggestedScale) / SuggestedScale > 1) {
      MACH3LOG_WARN(
          "Defined Global StepScale is {}, while suggested suggested {}", scale,
          SuggestedScale);
    }
  }
  proposer.params.global_scale = scale;
}

// ********************************************
int ParameterHandlerBase::GetParIndex(const std::string &name) const {
  // ********************************************
  return parlist.FindParameter(name);
}

// ********************************************
void ParameterHandlerBase::SetFixAllParameters() {
  // ********************************************
  proposer.params.isfree = Eigen::ArrayXi::Ones(proposer.NumParameters());
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const int i) {
  // ********************************************
  int propidx = parlist.SystematicParameterIndexToPCIndex(i);
  if (propidx != ParameterList::ParameterInPCABlock) {
    proposer.params.isfree[propidx] = false;
  } else {
    MACH3LOG_ERROR("You are trying Fix parameter {}, {}, which has been PCA'd.",
                   i, GetParName(i));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ********************************************
void ParameterHandlerBase::SetFixParameter(const std::string &name) {
  // ********************************************
  SetFixParameter(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::SetFreeAllParameters() {
  // ********************************************
  proposer.params.isfree = Eigen::ArrayXi::Ones(proposer.NumParameters());
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const int i) {
  // ********************************************
  int propidx = parlist.SystematicParameterIndexToPCIndex(i);
  if (propidx != ParameterList::ParameterInPCABlock) {
    proposer.params.isfree[propidx] = true;
  } else {
    MACH3LOG_ERROR(
        "You are trying Free parameter {}, {}, which has been PCA'd.", i,
        GetParName(i));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// ********************************************
void ParameterHandlerBase::SetFreeParameter(const std::string &name) {
  // ********************************************
  SetFreeParameter(GetParIndex(name));
}

// ********************************************
void ParameterHandlerBase::ToggleFixAllParameters() {
  // ********************************************
  throw std::runtime_error(
      "ParameterHandlerBase::ToggleFixAllParameters Removed");
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const int) {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::ToggleFixParameter Removed");
}

// ********************************************
void ParameterHandlerBase::ToggleFixParameter(const std::string &) {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::ToggleFixParameter Removed");
}

bool ParameterHandlerBase::IsParameterFixed(const int) const {
  throw std::runtime_error("ParameterHandlerBase::IsParameterFixed Removed");
}

// ********************************************
bool ParameterHandlerBase::IsParameterFixed(const std::string &) const {
  // ********************************************
  throw std::runtime_error("ParameterHandlerBase::IsParameterFixed Removed");
}

// ********************************************
void ParameterHandlerBase::SetFlatPrior(const int i, const bool eL) {
  // ********************************************

  int propidx = parlist.SystematicParameterIndexToPCIndex(i);
  if (propidx != ParameterList::ParameterInPCABlock) {
    parlist.params.flatprior[i] = eL;
  } else {
    MACH3LOG_ERROR("You are trying set a flat prior for parameter {}, {}, "
                   "which has been PCA'd.",
                   i, GetParName(i));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

void ParameterHandlerBase::SetIndivStepScale(const int ParameterIndex,
                                             const double StepScale) {
  int idx = parlist.SystematicParameterIndexToPCIndex(ParameterIndex);
  if (idx == ParameterList::ParameterInPCABlock) {
    MACH3LOG_ERROR("You are trying to set the step scale of parameter {}, {}, "
                   "which has been PCA'd.",
                   ParameterIndex, GetParName(ParameterIndex));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  proposer.params.scale[idx] = StepScale;
}

double ParameterHandlerBase::GetIndivStepScale(const int ParameterIndex) const {
  int idx = parlist.SystematicParameterIndexToPCIndex(ParameterIndex);
  if (idx == ParameterList::ParameterInPCABlock) {
    MACH3LOG_ERROR("You are trying to get the step scale of parameter {}, {}, "
                   "which has been PCA'd.",
                   ParameterIndex, GetParName(ParameterIndex));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return proposer.params.scale[idx];
}

// ********************************************
void ParameterHandlerBase::SetIndivStepScale(
    const std::vector<double> &stepscale) {
  // ********************************************
  if (int(stepscale.size()) != parlist.NumSystematicBasisParameters()) {
    MACH3LOG_WARN(
        "Stepscale vector not equal to number of parameters. Quitting..");
    MACH3LOG_WARN("Size of argument vector: {}", stepscale.size());
    MACH3LOG_WARN("Expected size: {}", parlist.NumSystematicBasisParameters());
    return;
  }

  for (int i = 0; i < parlist.NumSystematicBasisParameters(); ++i) {

    int idx = parlist.SystematicParameterIndexToPCIndex(i);
    if (idx == ParameterList::ParameterInPCABlock) {
      continue;
    }

    proposer.params.scale[idx] = stepscale[i];
  }

  PrintIndivStepScale();
}

// ********************************************
void ParameterHandlerBase::PrintIndivStepScale() const {
  // ********************************************
  MACH3LOG_INFO("============================================================");
  MACH3LOG_INFO("{:<{}} | {:<11}", "Parameter:", PrintLength, "Step scale:");
  for (int iParam = 0; iParam < parlist.NumSystematicBasisParameters();
       iParam++) {
    int idx = parlist.SystematicParameterIndexToPCIndex(iParam);
    if (idx == ParameterList::ParameterInPCABlock) {
      MACH3LOG_INFO("{:<{}} | {:<11}", GetParName(idx), PrintLength, "PCA'd");
    } else {
      MACH3LOG_INFO("{:<{}} | {:<11}", GetParName(idx), PrintLength,
                    proposer.params.scale[idx]);
    }
  }
  MACH3LOG_INFO("============================================================");
}

// ********************************************
void ParameterHandlerBase::ResetIndivStepScale() {
  // ********************************************
  proposer.params.global_scale = 1;
  proposer.params.scale = Eigen::ArrayXd::Ones(proposer.NumParameters());
}

// ********************************************
// HW: Code for throwing from separate throw matrix, needs to be set after init
// to ensure pos-def
void ParameterHandlerBase::SetThrowMatrix(TMatrixDSym *cov) {
  // ********************************************

  M3::MakeMatrixPosDef(cov);
  Eigen::MatrixXd covariance = M3::ROOTToEigen(*cov);

  if (parlist.pca.enabled) {
    Eigen::MatrixXd covariance_pc = Eigen::MatrixXd::Zero(
        proposer.NumParameters(), proposer.NumParameters());

    covariance_pc.topLeftCorner(parlist.pca.first_index,
                                parlist.pca.first_index) =
        parlist.params.covariance.topLeftCorner(parlist.pca.first_index,
                                                parlist.pca.first_index);

    covariance_pc.block(parlist.pca.first_index, parlist.pca.first_index,
                        parlist.pca.npc_parameters(),
                        parlist.pca.npc_parameters()) =
        Eigen::MatrixXd::Identity(parlist.pca.npc_parameters(),
                                  parlist.pca.npc_parameters());

    covariance_pc.bottomRightCorner(parlist.pca.ntail, parlist.pca.ntail) =
        parlist.params.covariance.bottomRightCorner(parlist.pca.ntail,
                                                    parlist.pca.ntail);

    proposer.SetProposalMatrix(covariance_pc);
  } else {
    proposer.SetProposalMatrix(covariance);
  }
}

void ParameterHandlerBase::UpdateThrowMatrix(TMatrixDSym *) {
  throw std::runtime_error("ParameterHandlerBase::UpdateThrowMatrix Removed");
}

double ParameterHandlerBase::GetThrowMatrix(const int, const int) const {
  throw std::runtime_error("ParameterHandlerBase::GetThrowMatrix(i,j) Removed");
}

TH2D *ParameterHandlerBase::GetCorrelationMatrix() {
  return M3::GetCorrelationMatrix(GetName(), parlist.params.name,
                                  parlist.params.covariance);
}

const double *ParameterHandlerBase::RetPointer(const int) {
  throw std::runtime_error("ParameterHandlerBase::RetPointer Removed");
}

// ********************************************
// HW : Here be adaption
void ParameterHandlerBase::InitialiseAdaption(const YAML::Node &) {
  // ********************************************
  throw std::runtime_error(
      "ParameterHandlerBase::InitialiseAdaption Not Yet Implemented");
}

void ParameterHandlerBase::SaveAdaptiveToFile(const std::string &,
                                              const std::string &) {
  throw std::runtime_error(
      "ParameterHandlerBase::SaveAdaptiveToFile Not Yet Implemented");
}

// ********************************************
// Truely adaptive MCMC!
void ParameterHandlerBase::UpdateAdaptiveCovariance() {
  // ********************************************
  // Updates adaptive matrix
  throw std::runtime_error(
      "ParameterHandlerBase::UpdateAdaptiveCovariance Not Yet Implemented");
}

// ********************************************
// KS: After step scale, prefit etc. value were modified save this modified
// config.
void ParameterHandlerBase::SaveUpdatedMatrixConfig() {
  // ********************************************
  throw std::runtime_error(
      "ParameterHandlerBase::SaveUpdatedMatrixConfig Not Yet Implemented");
}

void ParameterHandlerBase::SetNumberOfSteps(const int) {
  throw std::runtime_error(
      "ParameterHandlerBase::SetNumberOfSteps Not Yet Implemented");
}

// ********************************************
bool ParameterHandlerBase::AppliesToSample(
    const int SystIndex, const std::string &SampleName) const {
  // ********************************************
  // Empty means apply to all
  if (parlist.params.samples[SystIndex].size() == 0)
    return true;

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

  for (size_t i = 0; i < parlist.params.samples[SystIndex].size(); i++) {
    // Convert to low case to not be case sensitive
    std::string pattern = parlist.params.samples[SystIndex][i];
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
void ParameterHandlerBase::SetTune(const std::string &) {
  // ********************************************
  throw std::runtime_error(
      "ParameterHandlerBase::SaveUpdatedMatrixConfig Not Yet Implemented");
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
