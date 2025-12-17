#include "Parameters/Parameterlist.h"

#include "Parameters/PCA.h"
#include "Parameters/ParameterHandlerUtils.h"

void ParameterList::AddParameters(std::vector<ParamInfo> const &new_params) {

  for (auto const &p : new_params) {
    params.name.push_back(p.name);
    params.fancy_name.push_back(p.fancy_name);
    params.samples.push_back(p.affected_samples);
    params.flip_pivot.push_back(p.flip_pivot);
    params.circ_bounds.push_back(p.circ_bounds);
  }

  size_t new_block_start = params.prefit.size();

  params.prefit.conservativeResize(params.name.size());
  params.error.conservativeResize(params.name.size());
  params.lowbound.conservativeResize(params.name.size());
  params.upbound.conservativeResize(params.name.size());
  params.stepscale.conservativeResize(params.name.size());
  params.flatprior.conservativeResize(params.name.size());
  params.isfree.conservativeResize(params.name.size());

  params.covariance.conservativeResize(params.name.size(), params.name.size());

  for (size_t i = 0; i < new_params.size(); ++i) {
    params.prefit[new_block_start + i] = new_params[i].prefit;
    params.error[new_block_start + i] = new_params[i].error;
    params.lowbound[new_block_start + i] = new_params[i].bounds[0];
    params.upbound[new_block_start + i] = new_params[i].bounds[1];
    params.stepscale[new_block_start + i] = new_params[i].stepscale;
    params.flatprior[new_block_start + i] = new_params[i].flatprior;
    params.isfree[new_block_start + i] = new_params[i].isfree;
    params.covariance(new_block_start + i, new_block_start + i) =
        new_params[i].error * new_params[i].error;
  }

  params.inv_covariance = Eigen::MatrixXd();
}

void ParameterList::SetParameterCorrelation(int pidi, int pidj, double corr) {
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

  params.inv_covariance = Eigen::MatrixXd();
}

void ParameterList::SetParameterAllCorrelations(
    int paramid, std::map<std::string, double> const &correlations) {

  params.covariance.row(paramid) =
      Eigen::VectorXd::Zero(params.covariance.rows());
  params.covariance.col(paramid) =
      Eigen::VectorXd::Zero(params.covariance.rows());

  params.covariance(paramid, paramid) =
      params.error[paramid] * params.error[paramid];

  for (auto const &[other_name, corr] : correlations) {
    SetParameterCorrelation(paramid, FindParameter(other_name), corr);
  }
}

int ParameterList::FindParameter(std::string const &name) {
  auto pit = std::find(params.name.begin(), params.name.end(), name);
  if (pit == params.name.end()) {
    MACH3LOG_ERROR("ParameterList manages no parameter named {}.", name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return int(std::distance(params.name.begin(), pit));
}

void ParameterList::ConstructTruncatedPCA(double threshold, int first,
                                          int last) {

  M3::EnsureNoOutOfBlockCorrelations(params.covariance, {first, last});

  pca.enabled = true;
  pca.first_index = first;
  pca.last_index = last;
  pca.ntail = int(params.prefit.size() - (last + 1));

  int blocksize = (last + 1) - first;

  pca.pc_to_syst_rotation = CalculateTruncatedPCARotation(
      params.covariance.block(pca.first_index, pca.first_index, blocksize,
                              blocksize),
      threshold);

  MACH3LOG_INFO("ParameterList::ConstructTruncatedPCA: threshold = {}, "
                "nsystparams = {} truncated to {} PC params.",
                threshold, pca.nrotated_syst_parameters(),
                pca.npc_parameters());
}

void ParameterList::RotatePCParameterValuesToSystematicBasis(
    Eigen::ArrayXd const &pc_vals, Eigen::ArrayXd &systematic_vals) {

  systematic_vals.head(pca.first_index) = pc_vals.head(pca.first_index);

  // 1) rotate the principle component parameter values into 'deltas' in the
  //  systematic basis into the by using the PCA matrix.
  // 2) add on the prefit values from to these deltas to return to the
  // systematic basis.
  systematic_vals.segment(pca.first_index, pca.nrotated_syst_parameters()) =
      params.prefit.segment(pca.first_index, pca.nrotated_syst_parameters()) +
      (pc_vals.segment(pca.first_index, pca.npc_parameters()).matrix() *
       pca.pc_to_syst_rotation)
          .array();

  systematic_vals.tail(pca.ntail) = pc_vals.tail(pca.ntail);
}

void ParameterList::RotateSystematicParameterValuesToPCBasis(
    Eigen::ArrayXd const &systematic_vals, Eigen::ArrayXd &pc_vals) {

  pc_vals.head(pca.first_index) = systematic_vals.head(pca.first_index);

  // 1) subtract the prefit values from the current parameter values so that the
  //   PC basis has a CV of 0
  // 2) rotate the resulting parameter 'deltas' in the systematic basis into the
  //  PC basis by using the transpose of the PCA matrix
  pc_vals.segment(pca.first_index, pca.npc_parameters()) =
      ((systematic_vals.segment(pca.first_index,
                                pca.nrotated_syst_parameters()) -
        params.prefit.segment(pca.first_index, pca.nrotated_syst_parameters()))
           .matrix() *
       pca.pc_to_syst_rotation.transpose())
          .array();

  pc_vals.tail(pca.ntail) = systematic_vals.tail(pca.ntail);
}

int ParameterList::SystematicParameterIndexToPCIndex(int i) {
  if ((!pca.enabled) || (i < pca.first_index)) {
    return i;
  } else if ((i >= pca.first_index) && (i <= pca.last_index)) {
    return ParameterInPCABlock;
  } else {
    return i - pca.nrotated_syst_parameters() + pca.npc_parameters();
  }
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

StepProposer ParameterList::MakeProposer() {
  StepProposer proposer;

  if (!pca.enabled) {
    proposer.SetProposalMatrix(params.covariance);
    proposer.SetParameterValues(params.prefit);
  } else {
    // if we are using a PCA rotation then the proposer should only be aware of
    // the PC parameter basis.

    // copy the unrotated parameters and covariance blocks
    int nproposer_parameters =
        pca.first_index + pca.npc_parameters() + pca.ntail;

    Eigen::ArrayXd prefit_pc = Eigen::ArrayXd::Zero(nproposer_parameters);

    RotateSystematicParameterValuesToPCBasis(params.prefit, prefit_pc);

    proposer.SetParameterValues(prefit_pc);

    Eigen::MatrixXd covariance_pc =
        Eigen::MatrixXd::Zero(nproposer_parameters, nproposer_parameters);

    covariance_pc.topLeftCorner(pca.first_index, pca.first_index) =
        params.covariance.topLeftCorner(pca.first_index, pca.first_index);

    covariance_pc.block(pca.first_index, pca.first_index, pca.npc_parameters(),
                        pca.npc_parameters()) =
        Eigen::MatrixXd::Identity(pca.npc_parameters(), pca.npc_parameters());

    covariance_pc.bottomRightCorner(pca.ntail, pca.ntail) =
        params.covariance.bottomRightCorner(pca.ntail, pca.ntail);

    proposer.SetProposalMatrix(covariance_pc);
  }

  for (int i = 0; i < params.prefit.size(); ++i) {

    int idx = SystematicParameterIndexToPCIndex(i);

    if (params.flip_pivot[i].first) {
      if (idx == ParameterInPCABlock) {
        MACH3LOG_ERROR("Parameter, {}, index: {} is in the PCA block but "
                       "contained a special flip proposal definition",
                       params.name[i], i);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      proposer.special_proposal.flips.push_back(
          {idx, params.flip_pivot[i].second});
    }

    if (std::get<0>(params.circ_bounds[i])) {
      if (idx == ParameterInPCABlock) {
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
  }

  return proposer;
}

std::string ParameterList::SystematicParameterToString(int i) {
  return fmt::format(
      R"(name: {}, fancy_name: {}
  prefit: {:.3g}, error: {:.3g}, bounds: [ {:.3g}, {:.3g} ], stepscale: {:.3g}
  flatprior: {}, isfree: {}
  is_pca: {}, has_flip: {}, has_circ_bounds: {})",
      params.name[i], params.fancy_name[i], params.prefit[i], params.error[i],
      params.lowbound[i], params.upbound[i], params.stepscale[i],
      params.flatprior[i], params.isfree[i],
      SystematicParameterIndexToPCIndex(i) == ParameterInPCABlock,
      params.flip_pivot[i].first, std::get<0>(params.circ_bounds[i]));
}

ParameterList MakeFromYAML(YAML::Node const &config) {

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

ParameterList MakeFromYAML(const std::vector<std::string> &YAMLFiles) {

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

ParameterList MakeFromTMatrix(std::string const &name,
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
