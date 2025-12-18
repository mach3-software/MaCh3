#pragma once

#include "Parameters/StepProposer.h"

#include "Manager/Manager.h"

#include "TMatrixDSym.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Walloca"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include "Eigen/Dense"
#pragma GCC diagnostic pop

#include <array>
#include <map>
#include <random>
#include <string>
#include <vector>

struct ParameterList {

  struct {
    std::vector<std::string> name, fancy_name;
    Eigen::ArrayXd prefit, error, lowbound, upbound, stepscale;
    Eigen::ArrayXi flatprior, isfree;
    std::vector<std::vector<std::string>> samples;

    Eigen::MatrixXd covariance;
    Eigen::MatrixXd inv_covariance;

    // enabled, flip pivot value
    std::vector<std::pair<bool, double>> flip_pivot;
    // enabled, low bound, up bound
    std::vector<std::tuple<bool, double, double>> circ_bounds;
  } params;

  ParameterList() {};

  struct ParamInfo {
    std::string name, fancy_name;
    double prefit, error, stepscale;
    std::array<double, 2> bounds;
    bool flatprior, isfree;
    std::vector<std::string> affected_samples;
    std::pair<bool, double> flip_pivot;
    std::tuple<bool, double, double> circ_bounds;
  };

  void InsertParameters(int insert_at,
                        std::vector<ParamInfo> const &new_params);
  void InsertParameter(int insert_at, ParamInfo const &new_param) {
    InsertParameters(insert_at, {
                                    new_param,
                                });
  }

  void AddParameters(std::vector<ParamInfo> const &new_params) {
    InsertParameters(NumSystematicBasisParameters(), new_params);
  }
  void AddParameter(ParamInfo const &new_param) {
    InsertParameter(NumSystematicBasisParameters(), new_param);
  }

  std::string SystematicParameterToString(int i) const;

  void SetParameterCorrelation(int pidi, int pidj, double corr);
  void SetParameterAllCorrelations(
      int paramid, std::map<std::string, double> const &correlations);

  int FindParameter(std::string const &name) const;
  int FindParameterByFancyName(std::string const &fancy_name) const;
  int NumSystematicBasisParameters() const { return int(params.prefit.size()); }
  int NumPCBasisParameters() const {
    return NumSystematicBasisParameters() - pca.nrotated_syst_parameters() +
           pca.npc_parameters();
  }

  struct {
    int first_index, last_index, ntail;
    Eigen::MatrixXd pc_to_syst_rotation;
    Eigen::MatrixXd syst_to_pc_rotation;

    int nrotated_syst_parameters() const {
      return int(pc_to_syst_rotation.rows());
    }
    int npc_parameters() const { return int(pc_to_syst_rotation.cols()); }
  } pca;

  void ConstructTruncatedPCA(double threshold, int first, int last);

  void RotatePCParameterValuesToSystematicBasis(
      Eigen::ArrayXd const &pc_vals, Eigen::ArrayXd &systematic_vals) const;

  void RotateSystematicParameterValuesToPCBasis(
      Eigen::ArrayXd const &systematic_vals, Eigen::ArrayXd &pc_vals) const;

  constexpr static int ParameterInPCABlock = std::numeric_limits<int>::max();
  int SystematicParameterIndexToPCIndex(int i) const;

  double Chi2(Eigen::ArrayXd const &systematic_vals);

  StepProposer MakeProposer() const;

  struct {
    /// The input root file we read in
    std::vector<std::string> inputFiles;
    /// Stores config describing systematics
    YAML::Node YAMLDoc;
  } config;

  static ParameterList MakeFromYAML(YAML::Node const &config);
  static ParameterList MakeFromYAML(const std::vector<std::string> &YAMLFiles);
  static ParameterList MakeFromTMatrix(std::string const &name,
                                       std::string const &file);
};
