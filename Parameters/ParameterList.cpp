#include "Parameters/Parameterlist.h"

#include "Parameters/PCA.h"
#include "Parameters/ParameterHandlerUtils.h"

void ParameterList::InsertParameters(int insert_before,
                                     std::vector<ParamInfo> const &new_params) {

  for (auto it = new_params.rbegin(); it < new_params.rend(); ++it) {
    params.name.insert(params.name.begin() + insert_before, it->name);
    params.fancy_name.insert(params.fancy_name.begin() + insert_before,
                             it->fancy_name);
    params.samples.insert(params.samples.begin() + insert_before,
                          it->affected_samples);
    params.flip_pivot.insert(params.flip_pivot.begin() + insert_before,
                             it->flip_pivot);
    params.circ_bounds.insert(params.circ_bounds.begin() + insert_before,
                              it->circ_bounds);
  }

  auto N = params.prefit.size();

  params.prefit.conservativeResize(params.name.size());
  params.error.conservativeResize(params.name.size());
  params.lowbound.conservativeResize(params.name.size());
  params.upbound.conservativeResize(params.name.size());
  params.stepscale.conservativeResize(params.name.size());
  params.flatprior.conservativeResize(params.name.size());
  params.isfree.conservativeResize(params.name.size());

  params.prefit.tail(N - insert_before) =
      params.prefit.segment(insert_before, N - insert_before).eval();
  params.error.tail(N - insert_before) =
      params.error.segment(insert_before, N - insert_before).eval();
  params.lowbound.tail(N - insert_before) =
      params.lowbound.segment(insert_before, N - insert_before).eval();
  params.upbound.tail(N - insert_before) =
      params.upbound.segment(insert_before, N - insert_before).eval();
  params.stepscale.tail(N - insert_before) =
      params.stepscale.segment(insert_before, N - insert_before).eval();
  params.flatprior.tail(N - insert_before) =
      params.flatprior.segment(insert_before, N - insert_before).eval();
  params.isfree.tail(N - insert_before) =
      params.isfree.segment(insert_before, N - insert_before).eval();

  Eigen::MatrixXd ncovariance =
      Eigen::MatrixXd::Zero(params.name.size(), params.name.size());
  ncovariance.topLeftCorner(N, N) = params.covariance;

  ncovariance.bottomRightCorner(N - insert_before, N - insert_before) =
      ncovariance.block(insert_before, insert_before, N - insert_before,
                        N - insert_before);
  params.covariance = ncovariance;

  for (size_t i = 0; i < new_params.size(); ++i) {
    params.prefit[insert_before + i] = new_params[i].prefit;
    params.error[insert_before + i] = new_params[i].error;
    params.lowbound[insert_before + i] = new_params[i].bounds[0];
    params.upbound[insert_before + i] = new_params[i].bounds[1];
    params.stepscale[insert_before + i] = new_params[i].stepscale;
    params.flatprior[insert_before + i] = new_params[i].flatprior;
    params.isfree[insert_before + i] = new_params[i].isfree;

    params.covariance.row(insert_before + i).setZero();
    params.covariance.col(insert_before + i).setZero();

    params.covariance(insert_before + i, insert_before + i) =
        new_params[i].error * new_params[i].error;
  }

  params.inv_covariance = Eigen::MatrixXd();
}

void ParameterList::SetParameterCorrelation(int i, int j, double corr) {
  if (i == j) {
    MACH3LOG_ERROR(
        "AddParameterCorrelation cannot be used to set covariance "
        "matrix diagonal elements: ({0},{0}) attempted to be set to {1}",
        i, corr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  params.covariance(i, j) =
      corr * std::sqrt(params.covariance(i, i) * params.covariance(j, j));
  params.covariance(j, i) = params.covariance(i, j);

  params.inv_covariance = Eigen::MatrixXd();
}

void ParameterList::SetParameterAllCorrelations(
    int paramid, std::map<std::string, double> const &correlations) {

  params.covariance.row(paramid).setZero();
  params.covariance.col(paramid).setZero();

  params.covariance(paramid, paramid) =
      params.error[paramid] * params.error[paramid];

  for (auto const &[other_name, corr] : correlations) {
    SetParameterCorrelation(paramid, FindParameter(other_name), corr);
  }
}

int ParameterList::FindParameter(std::string const &name) const {
  auto pit = std::find(params.name.begin(), params.name.end(), name);
  if (pit == params.name.end()) {
    MACH3LOG_ERROR("ParameterList manages no parameter named {}.", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return int(std::distance(params.name.begin(), pit));
}

int ParameterList::FindParameterByFancyName(
    std::string const &fancy_name) const {
  auto pit =
      std::find(params.fancy_name.begin(), params.fancy_name.end(), fancy_name);
  if (pit == params.fancy_name.end()) {
    MACH3LOG_ERROR("ParameterList manages no parameter named {}.", fancy_name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return int(std::distance(params.fancy_name.begin(), pit));
}

double ParameterList::Chi2(Eigen::ArrayXd const &systematic_vals) {

  if (!params.inv_covariance.size()) {
    params.inv_covariance = params.covariance.inverse();
  }

  auto diff =
      params.flatprior.select(Eigen::VectorXd::Zero(params.prefit.size()),
                              systematic_vals - params.prefit);
  return diff.transpose() * params.inv_covariance * diff;
}

StepProposer ParameterList::MakePCAProposer(double threshold, int first,
                                            int last) const {
  StepProposer proposer;

  if ((first >= 0) && (first < params.prefit.size())) {

    if (first >= last) {
      MACH3LOG_ERROR("On constructing PCA'd StepProposer, given first index: "
                     "{} and incompatible last index: {}.",
                     first, last);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    M3::EnsureNoOutOfBlockCorrelations(params.covariance, {first, last});

    proposer.systematic_basis.pca.first_index = first;
    proposer.systematic_basis.pca.last_index = last;
    proposer.systematic_basis.pca.ntail =
        int(params.prefit.size() - (last + 1));
    proposer.systematic_basis.pca.offset = params.prefit;

    int blocksize = (last + 1) - first;

    std::tie(proposer.systematic_basis.pca.pc_to_syst_rotation,
             proposer.systematic_basis.pca.syst_to_pc_rotation) =
        CalculateTruncatedPCARotation(
            params.covariance.block(proposer.systematic_basis.pca.first_index,
                                    proposer.systematic_basis.pca.first_index,
                                    blocksize, blocksize),
            threshold);

    MACH3LOG_INFO(
        "ParameterList::MakePCAProposer: threshold = {}, "
        "nsystparams = {} truncated to {} PC params.",
        threshold,
        proposer.systematic_basis.pca.nrotated_systematic_parameters(),
        proposer.systematic_basis.pca.npc_parameters());
  } else {
    proposer.systematic_basis.pca.first_index = NumParameters();
    proposer.systematic_basis.pca.last_index = NumParameters();
    proposer.systematic_basis.pca.ntail = 0;
  }
  // if we are using a PCA rotation then the proposer should only be aware of
  // the PC parameter basis.

  // copy the unrotated parameters and covariance blocks
  int nproposer_parameters = proposer.NumProposalBasisParameters();

  proposer.SetSystematicParameterValues(params.prefit);

  Eigen::MatrixXd covariance_pc =
      Eigen::MatrixXd::Identity(nproposer_parameters, nproposer_parameters);

  covariance_pc.topLeftCorner(proposer.systematic_basis.pca.first_index,
                              proposer.systematic_basis.pca.first_index) =
      params.covariance.topLeftCorner(
          proposer.systematic_basis.pca.first_index,
          proposer.systematic_basis.pca.first_index);

  covariance_pc.bottomRightCorner(proposer.systematic_basis.pca.ntail,
                                  proposer.systematic_basis.pca.ntail) =
      params.covariance.bottomRightCorner(proposer.systematic_basis.pca.ntail,
                                          proposer.systematic_basis.pca.ntail);

  proposer.SetProposalMatrix(covariance_pc);

  // PCA parameters get a step scale of 1
  proposer.proposal_basis.scale = Eigen::ArrayXd::Ones(nproposer_parameters);
  proposer.proposal_basis.isfree = Eigen::ArrayXi::Ones(nproposer_parameters);

  for (int i = 0; i < NumParameters(); ++i) {

    int idx = proposer.GetProposalParameterIndexFromSystematicIndex(i);

    if (params.flip_pivot[i].first) {
      if (idx == StepProposer::ParameterInPCABlock) {
        MACH3LOG_ERROR("Parameter, {}, index: {} is in the PCA block but "
                       "contained a special flip proposal definition",
                       params.name[i], i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      proposer.special_proposal.flips.push_back(
          {idx, params.flip_pivot[i].second});
    }

    if (std::get<0>(params.circ_bounds[i])) {
      if (idx == StepProposer::ParameterInPCABlock) {
        MACH3LOG_ERROR(
            "Parameter, {}, index: {} is in the PCA block but "
            "contained a special circular bounds proposal definition",
            params.name[i], i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      proposer.special_proposal.circ_bounds.push_back(
          {idx, std::get<1>(params.circ_bounds[i]),
           std::get<2>(params.circ_bounds[i])});
    }

    if (idx == StepProposer::ParameterInPCABlock) {
      continue;
    }

    proposer.proposal_basis.scale[idx] = params.stepscale[i];
    proposer.proposal_basis.isfree[idx] = params.isfree[i];
  }

  proposer.proposal_basis.global_scale = 1;

  return proposer;
}

StepProposer ParameterList::MakeProposer() const {
  return MakePCAProposer(0, std::numeric_limits<int>::max(), 0);
}

std::string ParameterList::SystematicParameterToString(int i) const {
  return fmt::format(
      R"(name: {}, fancy_name: {}
    prefit: {:.3g}, error: {:.3g}, bounds: [ {:.3g}, {:.3g} ],
      stepscale: {:.3g}
      flatprior: {}, isfree: {} has_flip: {}, has_circ_bounds: {})",
      params.name[i], params.fancy_name[i], params.prefit[i], params.error[i],
      params.lowbound[i], params.upbound[i], params.stepscale[i],
      params.flatprior[i], params.isfree[i], params.flip_pivot[i].first,
      std::get<0>(params.circ_bounds[i]));
}

ParameterList ParameterList::MakeFromYAML(YAML::Node const &config) {

  ParameterList parlist;

  parlist.config.YAMLDoc = config;

  std::map<int, std::map<std::string, double>> parameter_correlations;
  std::vector<ParameterList::ParamInfo> param_infos;

  // ETA - read in the systematics.
  for (auto const &node : parlist.config.YAMLDoc["Systematics"]) {
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

    std::tuple<bool, double, double> param_circ_bounds{false, 0, 0};
    std::pair<bool, double> param_flip_pivot{false, 0};
    if (pardef["SpecialProposal"]) {

      bool CircEnabled = bool(pardef["SpecialProposal"]["CircularBounds"]);
      bool FlipEnabled = bool(pardef["SpecialProposal"]["FlipParameter"]);

      if (!CircEnabled && !FlipEnabled) {
        MACH3LOG_ERROR(
            "None of Special Proposal were enabled even though param "
            "{}, has SpecialProposal entry in Yaml",
            fancy_name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      if (CircEnabled) {
        auto circular_bounds = Get<std::array<double, 2>>(
            pardef["SpecialProposal"]["CircularBounds"], __FILE__, __LINE__);
        // KS: Make sure circular circular_bounds are within physical
        // circular_bounds. If we are outside of physics bound MCMC will never
        // explore such phase space region
        if ((circular_bounds[0] < bounds[0]) ||
            (circular_bounds[1] >= bounds[1])) {
          MACH3LOG_ERROR(
              "Circular circular_bounds [{}, {}] for parameter {} exceed "
              "physical circular_bounds [{}, {}]",
              circular_bounds[0], circular_bounds[1], fancy_name, bounds[0],
              bounds[1]);
          throw MaCh3Exception(__FILE__, __LINE__);
        }

        param_circ_bounds =
            std::make_tuple(true, circular_bounds[0], circular_bounds[1]);

        MACH3LOG_INFO(
            "Enabling CircularBounds for parameter {} with range [{}, {}]",
            fancy_name, circular_bounds[0], circular_bounds[1]);
      }

      if (FlipEnabled) {
        param_flip_pivot = std::make_pair(
            true, Get<double>(pardef["SpecialProposal"]["FlipParameter"],
                              __FILE__, __LINE__));

        MACH3LOG_INFO("Enabling Flipping for parameter {} with value {}",
                      fancy_name, param_flip_pivot.second);
      }

      if (CircEnabled && FlipEnabled) {

        const double fp = param_flip_pivot.second;
        const double low = std::get<1>(param_circ_bounds);
        const double high = std::get<2>(param_circ_bounds);

        if (fp < low || fp > high) {
          MACH3LOG_ERROR(
              "FlipParameter value {} for parameter {} is outside the "
              "CircularBounds [{}, {}]",
              fp, fancy_name, low, high);
          throw MaCh3Exception(__FILE__, __LINE__);
        }

        // Sanity check: ensure flipping any x in [low, high] keeps the result
        // in [low, high]
        const double flipped_low = 2 * fp - low;
        const double flipped_high = 2 * fp - high;
        const double min_flip = std::min(flipped_low, flipped_high);
        const double max_flip = std::max(flipped_low, flipped_high);

        if (min_flip < low || max_flip > high) {
          MACH3LOG_ERROR("Flipping about point {} for parameter {} would leave "
                         "circular bounds [{}, {}]",
                         param_flip_pivot.second, fancy_name, low, high);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    } // end special proposal handling

    if (pardef["Correlations"]) {
      for (auto const &corrn : pardef["Correlations"]) {
        for (auto const &corr : corrn) {
          parameter_correlations[param_id][corr.first.as<std::string>()] =
              corr.second.as<double>();
        }
      }
    }

    param_infos.emplace_back(ParameterList::ParamInfo{fancy_name,
                                                      fancy_name,
                                                      prefit,
                                                      error,
                                                      stepscale,
                                                      {bounds[0], bounds[1]},
                                                      flatprior,
                                                      !fixed,
                                                      samplenames,
                                                      param_flip_pivot,
                                                      param_circ_bounds});
  }

  parlist.AddParameters(param_infos);

  for (auto const &[paramid, correlations] : parameter_correlations) {

    // check that we don't have incompatible correlations defined the other way
    for (auto const &[oparam, corr] : correlations) {
      int oparamid = parlist.FindParameter(oparam);
      auto const &pname = parlist.params.name[paramid];
      auto const &opname = parlist.params.name[oparamid];
      if (!parameter_correlations.count(oparamid) ||
          !parameter_correlations.at(oparamid).count(pname)) {
        MACH3LOG_ERROR("Encountered inconsistent correlations defined between "
                       "parameters {0} and {1}. {0} reports a correlation of "
                       "{2} with {1} that is not reciprocated.",
                       pname, opname, corr);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      if (parameter_correlations.count(oparamid) &&
          parameter_correlations.at(oparamid).count(pname) &&
          std::fabs(parameter_correlations.at(oparamid).at(pname) - corr) >
              1E-8) {
        MACH3LOG_ERROR("Encountered inconsistent correlations defined between "
                       "parameters {} and {}: {} and {}.",
                       pname, opname, corr,
                       parameter_correlations.at(oparamid).at(pname));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }

    parlist.SetParameterAllCorrelations(paramid, correlations);
  }

  parlist.params.covariance = M3::MakeMatrixPosDef(parlist.params.covariance);

  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total",
                parlist.params.prefit.size());
  MACH3LOG_INFO("----------------");

  return parlist;
}

ParameterList
ParameterList::MakeFromYAML(const std::vector<std::string> &YAMLFiles) {

  MACH3LOG_INFO("Constructing instance of ParameterHandler using: ");
  for (auto const &yf : YAMLFiles) {
    MACH3LOG_INFO("  {}", yf);
  }
  MACH3LOG_INFO("as an input.");

  YAML::Node YAMLDoc;
  YAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for (auto const &yfn : YAMLFiles) {
    for (const auto &node : M3OpenConfig(yfn)["Systematics"]) {
      YAMLDoc["Systematics"].push_back(node);
    }
  }

  auto parlist = MakeFromYAML(YAMLDoc);
  parlist.config.inputFiles = YAMLFiles;

  MACH3LOG_INFO("Created parameter list from files: ");
  for (const auto &file : parlist.config.inputFiles) {
    MACH3LOG_INFO("{} ", file);
  }

  return parlist;
}

ParameterList ParameterList::MakeFromTMatrix(std::string const &name,
                                             std::string const &file) {

  ParameterList parlist;

  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  std::unique_ptr<TFile> infile{TFile::Open(file.c_str(), "READ")};
  if (infile->IsZombie()) {
    MACH3LOG_ERROR("Could not open input covariance ROOT file {} !!!", file);
    MACH3LOG_ERROR("Was about to retrieve matrix with name {}", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TMatrixDSym *CovMat = static_cast<TMatrixDSym *>(infile->Get(name.c_str()));

  if (!CovMat) {
    MACH3LOG_ERROR("Could not find covariance matrix name {} in file {}", name,
                   file);
    MACH3LOG_ERROR("Are you really sure {} exists in the file?", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::vector<ParameterList::ParamInfo> param_infos;
  for (int i = 0; i < CovMat->GetNrows(); ++i) {
    param_infos.emplace_back(ParameterList::ParamInfo{
        std::string(CovMat->GetName()) + "_" + std::to_string(i),
        "",
        0,
        std::sqrt((*CovMat)(i, i)),
        1,
        {-std::numeric_limits<double>::max(),
         std::numeric_limits<double>::max()},
        false,
        true,
        {},
        {false, 0},
        {false, 0, 0}});
  }
  parlist.AddParameters(param_infos);

  parlist.params.covariance = M3::MakeMatrixPosDef(M3::ROOTToEigen(*CovMat));

  return parlist;
}
