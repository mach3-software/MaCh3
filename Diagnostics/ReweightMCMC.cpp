//MaCh3 includes
#include "Manager/Manager.h"
#include "Fitters/MCMCProcessor.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TGraph.h"
_MaCh3_Safe_Include_End_ //}

// C++ includes
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <map>

/// @file ReweightMCMC.cpp
/// @brief This executable allow to reweight MCMC Chain, such technique is used to study impact of different priors without rerunning MCMC
/// @ingroup MaCh3DiagnosticProcessing
///
/// @author David Riley
/// @author Evan Goodman

/// Structure to hold reweight configuration
struct ReweightConfig {
  std::string key;       ///< The YAML key for this reweight
  std::string name;
  std::string type;  ///< "Gaussian", "TGraph2D"
  int dimension;     ///< 1 or 2
  std::vector<std::string> paramNames; ///< Parameter names.
  std::vector<std::vector<double>> newPriorValues; ///< new [mean, sigma] pairs
  std::vector<std::vector<double>> oldPriorValues; ///< new [mean, sigma] pairs
  std::vector<bool> flatPrior;

  std::string weightBranchName; ///< Output weight branch name.
  bool enabled;

  // For TGraph 1D or 2D
  std::string fileName;         ///< ROOT file containing graph data.
  std::string graphName;        ///< Graph name in the ROOT file.

  // For TGraph1D
  std::unique_ptr<TGraph> graph_1D; ///< 1D interpolation graph.

  // For TGraph2D
  std::string hierarchyType; ///< "NO", "IO", or "auto"
  std::unique_ptr<TGraph2D> graph_NO; ///< Normal Ordering graph.
  std::unique_ptr<TGraph2D> graph_IO; ///< Inverted Ordering graph.
};

/// @brief Function to interpolate 2D graph for Normal Ordering
double Graph_interpolateNO(TGraph2D* graph, double theta13, double dm32)
{
  if (!graph) {
    MACH3LOG_ERROR("Graph pointer is null");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  double xmax = graph->GetXmax();
  double xmin = graph->GetXmin();
  double ymin = graph->GetYmin();
  double ymax = graph->GetYmax();

  double chiSquared, prior;

  if (theta13 < xmax && theta13 > xmin && dm32 < ymax && dm32 > ymin) {
    chiSquared = graph->Interpolate(theta13, dm32);
    prior = std::exp(-0.5 * chiSquared);
  } else {
    prior = 0.0;
  }

  return prior;
}

/// @brief Function to interpolate 2D graph for Inverted Ordering  
double Graph_interpolateIO(TGraph2D* graph, double theta13, double dm32)
{
  if (!graph) {
    MACH3LOG_ERROR("Graph pointer is null");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  double xmax = graph->GetXmax();
  double xmin = graph->GetXmin();
  double ymax = graph->GetYmax();
  double ymin = graph->GetYmin();

  // The dm32 value is positive for in the TGraph2D so we should compare the abs value of the -delM32 values to get the chisq
  double mod_dm32 = std::abs(dm32);
  double chiSquared, prior;

  if (theta13 < xmax && theta13 > xmin && mod_dm32 < ymax && mod_dm32 > ymin) {
    chiSquared = graph->Interpolate(theta13, mod_dm32);
    prior = std::exp(-0.5 * chiSquared);
  } else {
    prior = 0.0;
  }

  return prior;
}

/// @brief Function to interpolate 1D graph
double Graph_interpolate1D(TGraph* graph, double theta13)
{
  /// @todo double check implementation of TGraph interpolation for 1D
  if (!graph) {
    MACH3LOG_ERROR("Graph pointer is null");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  double xmax = -999999999;
  double xmin =  999999999;

  for (int i = 0; i < graph->GetN(); i++) {
    double x = graph->GetX()[i];
    if (x > xmax) xmax = x;
    if (x < xmin) xmin = x;
  }

  double chiSquared, prior;

  if (theta13 < xmax && theta13 > xmin) {
    chiSquared = graph->Eval(theta13);
    prior = std::exp(-0.5 * chiSquared);
  } else {
    prior = 0.0;
  }

  return prior;
}

/// @brief Get parameter information from MCMCProcessor
bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName,
                      double& mean, double& sigma, bool& isFlat)
{
  // Try to find the parameter index
  int paramIndex = processor->GetParamIndexFromName(paramName);

  if (paramIndex == M3::_BAD_INT_) { // This indicate parameter not found
    return false;
  }

  // Get parameter information
  TString title;
  processor->GetNthParameter(paramIndex, mean, sigma, title);
  isFlat = processor->GetParamFlat(paramIndex);

  return true;
}

/// @brief Load reweighting setting like 1D or 2D from YAML config
void LoadReweightingSettings(std::vector<ReweightConfig>& reweightConfigs, const YAML::Node& reweight_settings) {
  // iterate through the keys in the reweighting yaml creating and storing the ReweightConfig as we go
  for (const auto& reweight : reweight_settings) {
    const std::string& reweightKey = reweight.first.as<std::string>();
    const YAML::Node& reweightConfigNode = reweight.second;

    // Check if this particular reweight is enabled !!! Currently only support one reweight at a time so this defaults to enabled
    if (!GetFromManager<bool>(reweightConfigNode["Enabled"], true, __FILE__ , __LINE__)) {
      MACH3LOG_INFO("Skipping disabled reweight: {}", reweightKey);
      continue;
    }

    ReweightConfig reweightConfig;
    reweightConfig.key = reweightKey;
    reweightConfig.name = GetFromManager<std::string>(reweightConfigNode["ReweightName"], reweightKey, __FILE__ , __LINE__);
    reweightConfig.type = GetFromManager<std::string>(reweightConfigNode["ReweightType"], "Gaussian", __FILE__ , __LINE__);
    reweightConfig.dimension = GetFromManager<int>(reweightConfigNode["ReweightDim"], 1);
    // reweightConfig.weightBranchName = "Weight_" + reweightKey; // for now all weights will be stored in branches called Weight until multi weight support
    reweightConfig.weightBranchName = "Weight";
    reweightConfig.enabled = true;

    auto paramNames = GetFromManager<std::vector<std::string>>(reweightConfigNode["ReweightVar"], {}, __FILE__ , __LINE__);
    reweightConfig.paramNames = paramNames;
    reweightConfig.oldPriorValues.resize(paramNames.size());
    reweightConfig.flatPrior.resize(paramNames.size());

    // Handle different reweight types as they fill different members
    if (reweightConfig.dimension == 1) {
      if (reweightConfig.type == "Gaussian") {
        // For Gaussian reweights, we need the parameter name(s) and prior values (mean, sigma pairs)
        // Get prior values - handle both single [mean, sigma] pair and list of pairs for safety
        auto priorNode = reweightConfigNode["ReweightPrior"];
        std::vector<std::vector<double>> allPriorValues;

        if (priorNode.IsSequence() && priorNode.size() > 0) {
          // Check if first element is a number (single [mean, sigma] pair) or sequence (list of pairs)
          if (priorNode[0].IsScalar()) {
            // Single [mean, sigma] pair - convert to list format
            auto newPriorValues = GetFromManager<std::vector<double>>(priorNode, {}, __FILE__ , __LINE__);
            if (newPriorValues.size() == 2) {
              allPriorValues.push_back(newPriorValues);
            }
          } else {
            // List of [mean, sigma] pairs
            for (const auto& priorPair : priorNode) {
              auto newPriorValues = GetFromManager<std::vector<double>>(priorPair, {}, __FILE__ , __LINE__);
              if (newPriorValues.size() == 2) {
                allPriorValues.push_back(newPriorValues);
              }
            }
          }
        }

        reweightConfig.newPriorValues = allPriorValues;

        if (paramNames.empty() || allPriorValues.empty() || paramNames.size() != allPriorValues.size()) {
          MACH3LOG_ERROR("Invalid Gaussian reweight configuration for {}: {} parameters, {} prior pairs",
                          reweightKey, paramNames.size(), allPriorValues.size());
          continue;
        }
      } else if (reweightConfig.type == "TGraph") {
        // For TGraph reweights, we need the parameter name and the TGraph file and name
        auto fileName = GetFromManager<std::string>(reweightConfigNode["ReweightPrior"]["file"], "", __FILE__ , __LINE__);
        auto graphName = GetFromManager<std::string>(reweightConfigNode["ReweightPrior"]["graph_name"], "", __FILE__ , __LINE__);
        reweightConfig.fileName = fileName;
        reweightConfig.graphName = graphName;

        if (paramNames.empty() || paramNames.size() != 1 || fileName.empty() || graphName.empty()) {
          MACH3LOG_ERROR("Invalid TGraph reweight configuration for {}", reweightKey);
          continue;
        }

        // Load the 1D graph
        MACH3LOG_INFO("Loading 1D constraint from file: {} (graph: {})", reweightConfig.fileName, reweightConfig.graphName);
        auto constraintFile = std::unique_ptr<TFile>(TFile::Open(reweightConfig.fileName.c_str(), "READ"));
        if (!constraintFile || constraintFile->IsZombie()) {
          MACH3LOG_ERROR("Failed to open constraint file: {}", reweightConfig.fileName);
          continue;
        }

        std::unique_ptr<TGraph> graph(constraintFile->Get<TGraph>(reweightConfig.graphName.c_str()));
        if (graph) {
          // Create a completely independent copy
          auto cloned_graph = static_cast<TGraph*>(graph->Clone());
          cloned_graph->SetBit(kCanDelete, true); // Allow ROOT to delete it when we're done
          reweightConfig.graph_1D = std::unique_ptr<TGraph>(cloned_graph);
          MACH3LOG_INFO("Loaded 1D graph: {}", reweightConfig.graphName);
        } else {
          MACH3LOG_ERROR("Failed to load graph: {}", reweightConfig.graphName);
          continue;
        }
      } else {
        MACH3LOG_ERROR("Unknown 1D reweight type: {} for {}", reweightConfig.type, reweightKey);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } else if (reweightConfig.dimension == 2) {
        // 2D reweights need 2 parameter names
        if (paramNames.size() != 2) {
          MACH3LOG_ERROR("2D reweighting requires exactly 2 parameter names for {}", reweightKey);
          continue;
        }

        if (reweightConfig.type == "TGraph2D") {
          auto priorConfig = reweightConfigNode["ReweightPrior"];
          reweightConfig.fileName = GetFromManager<std::string>(priorConfig["file"], "", __FILE__ , __LINE__);
          reweightConfig.graphName = GetFromManager<std::string>(priorConfig["graph_name"], "", __FILE__ , __LINE__);
          reweightConfig.hierarchyType = GetFromManager<std::string>(priorConfig["hierarchy"], "auto", __FILE__ , __LINE__);

          if (reweightConfig.fileName.empty() || reweightConfig.graphName.empty()) {
            MACH3LOG_ERROR("Invalid TGraph2D configuration for {}", reweightKey);
            continue;
          }

          // Load the 2D graphs
          MACH3LOG_INFO("Loading 2D constraint from file: {} (graph: {})", reweightConfig.fileName, reweightConfig.graphName);
          auto constraintFile = std::unique_ptr<TFile>(TFile::Open(reweightConfig.fileName.c_str(), "READ"));
          if (!constraintFile || constraintFile->IsZombie()) {
            MACH3LOG_ERROR("Failed to open constraint file: {}", reweightConfig.fileName);
            continue;
          }

          // Load both NO and IO graphs if hierarchy is auto
          if (reweightConfig.hierarchyType == "auto" || reweightConfig.hierarchyType == "NO") {
            std::string graphName_NO = reweightConfig.graphName + "_NO";
            MACH3LOG_INFO("Loading NO graph: {}", graphName_NO);

            std::unique_ptr<TGraph2D> graph_NO(constraintFile->Get<TGraph2D>(graphName_NO.c_str()));
            if (graph_NO) {
              // Create a completely independent copy
              auto cloned_graph = static_cast<TGraph2D*>(graph_NO->Clone());
              cloned_graph->SetDirectory(nullptr); // Detach from file
              cloned_graph->SetBit(kCanDelete, true); // Allow ROOT to delete it when we're done
              reweightConfig.graph_NO = std::unique_ptr<TGraph2D>(cloned_graph);
              MACH3LOG_INFO("Loaded NO graph: {}", graphName_NO);
            } else {
              MACH3LOG_ERROR("Failed to load NO graph: {}", graphName_NO);
            }
          }

          if (reweightConfig.hierarchyType == "auto" || reweightConfig.hierarchyType == "IO") {
            std::string graphName_IO = reweightConfig.graphName + "_IO";
            MACH3LOG_INFO("Loading IO graph: {}", graphName_IO);
            std::unique_ptr<TGraph2D> graph_IO(constraintFile->Get<TGraph2D>(graphName_IO.c_str()));
            if (graph_IO) {
              // Create a completely independent copy
              auto cloned_graph = static_cast<TGraph2D*>(graph_IO->Clone());
              cloned_graph->SetDirectory(nullptr); // Detach from file
              cloned_graph->SetBit(kCanDelete, true); // Allow ROOT to delete it when we're done
              reweightConfig.graph_IO = std::unique_ptr<TGraph2D>(cloned_graph);
              MACH3LOG_INFO("Loaded IO graph: {}", graphName_IO);
            } else {
              MACH3LOG_ERROR("Failed to load IO graph: {}", graphName_IO);
            }
          }

          constraintFile->Close();
        } else {
          MACH3LOG_ERROR("Unknown 2D reweight type: {} for {}", reweightConfig.type, reweightKey);
          continue;
        }
    } else {
      MACH3LOG_ERROR("Unsupported reweight dimension: {} for {}", reweightConfig.dimension, reweightKey);
      continue;
    } // end check over dimensions

    reweightConfigs.push_back(std::move(reweightConfig));
    MACH3LOG_INFO("Added reweight configuration: {} ({}D, type: {})", reweightConfigs.back().name, reweightConfigs.back().dimension, reweightConfigs.back().type);
  }

  if (reweightConfigs.empty()) {
    MACH3LOG_ERROR("No valid reweight configurations found in config file");
    throw MaCh3Exception(__FILE__, __LINE__);
  } else if (reweightConfigs.size() > 1) {    // check number of ReweightConfigs, currently maximum supported is 1 due to structure of ProcessMCMC
    MACH3LOG_ERROR("Currently only one reweight configuration is supported at a time, found {}", reweight_settings.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

/// @brief Calculate 1D weight
[[nodiscard]] double Get1DWeight(const ReweightConfig& rwConfig,
                                 const std::map<std::string, double>& paramValues) {
  double weight = 1.;
  if(rwConfig.type == "Gaussian") {
    auto& paramNames = rwConfig.paramNames;
    //KS: Calculate reweight weight. Weights are multiplicative so we can do several reweights at once.
    /// @warning Big limitation is that code only works for uncorrelated parameters :(
    for (unsigned int j = 0; j < paramNames.size(); ++j)
    {
      auto name = rwConfig.paramNames.at(j);
      // Extract means and sigmas from the prior pairs
      auto newPriorPair = rwConfig.newPriorValues[j];
      double NewCentral = newPriorPair[0]; // mean
      double NewError  = newPriorPair[1];   // sigma

      double new_chi = (paramValues.at(name) - NewCentral)/NewError;
      double new_prior = std::exp(-0.5 * new_chi * new_chi);

      double old_chi = -1;
      double old_prior = -1;
      if(rwConfig.flatPrior[j]) {
        old_prior = 1.0;
      } else {
        auto oldPriorPair = rwConfig.oldPriorValues[j];
        double OldCentral = newPriorPair[0]; // mean
        double OldError  = newPriorPair[1];   // sigma

        old_chi = (paramValues.at(name) - OldCentral)/OldError;
        old_prior = std::exp(-0.5 * old_chi * old_chi);
      }
      weight *= new_prior/old_prior;
    }
  } else if (rwConfig.type == "TGraph") {
    double paramValue = paramValues.at(rwConfig.paramNames.at(0));
    weight = Graph_interpolate1D(rwConfig.graph_1D.get(), paramValue);
  } else {
    MACH3LOG_ERROR("Unsupported 1D reweight type: {} for {}", rwConfig.type, rwConfig.key);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return weight;
}

/// @brief Calculate 2D weight
[[nodiscard]] double Get2DWeight(const ReweightConfig& rwConfig,
                                 const std::map<std::string, double>& paramValues) {
  double weight = 1.;
  if (rwConfig.type == "TGraph2D") {
    double dm32 = paramValues.at(rwConfig.paramNames.at(0));
    double theta13 = paramValues.at(rwConfig.paramNames.at(1));
    if (dm32 > 0) {
      // Normal Ordering
      if (rwConfig.graph_NO) {
        weight = Graph_interpolateNO(rwConfig.graph_NO.get(), theta13, dm32);
      } else {
        MACH3LOG_ERROR("NO graph not available for {}", rwConfig.key);
        weight = 0.0;
      }
    } else {
      // Inverted Ordering
      if (rwConfig.graph_IO) {
        weight = Graph_interpolateIO(rwConfig.graph_IO.get(), theta13, dm32);
      } else {
        MACH3LOG_ERROR("IO graph not available for {}", rwConfig.key);
        weight = 0.0;
      }
    }
  }
  return weight;
}

/// @brief Main executable responsible for reweighting MCMC chains
/// @param inputFile MCMC Chain file path
/// @param configFile Config file with reweighting settings
/// @author David Riley
/// @author Evan Goodman
/// @todo add a generic 2D reweight that is not dm32 and theta13 specific DWR
void ReweightMCMC(const std::string& configFile, const std::string& inputFile)
{
  MACH3LOG_INFO("File for reweighting: {} with config {}", inputFile, configFile);
  // Load configuration
  YAML::Node reweight_yaml = M3OpenConfig(configFile);
  YAML::Node reweight_settings = reweight_yaml["ReweightMCMC"];

  // Parse all reweight configurations first
  std::vector<ReweightConfig> reweightConfigs;

  LoadReweightingSettings(reweightConfigs, reweight_settings);

  // Create MCMCProcessor to get parameter information
  auto processor = std::make_unique<MCMCProcessor>(inputFile);
  processor->Initialise();

  // Validate that all required parameters exist in the chain
  /// @todo Get list only of unique parameters, this is repeating unnecessarily when adding more than 1 weight DWR
  for (auto& rwConfig : reweightConfigs) {
    for (size_t i = 0; i < rwConfig.paramNames.size(); ++i) {
      const auto& paramName = rwConfig.paramNames[i];
      int paramIndex = processor->GetParamIndexFromName(paramName);
      if (paramIndex == M3::_BAD_INT_) {
        MACH3LOG_ERROR("Parameter {} not found in MCMC chain", paramName);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MACH3LOG_INFO("Parameter {} found in chain", paramName);

      double mean, sigma;
      bool isFlat;
      GetParameterInfo(processor.get(), paramName, mean, sigma, isFlat);

      rwConfig.oldPriorValues[i] = {mean, sigma};
      rwConfig.flatPrior[i] = isFlat;
    }
  }

  /// @todo Finish Asimov shifting implementation, for now just warn that Asimovs are not being properly handled
  // Get the settings for the MCMC
  auto TempFile = std::unique_ptr<TFile>(TFile::Open(inputFile.c_str(), "READ"));
  if (!TempFile || TempFile->IsZombie()) {
    MACH3LOG_ERROR("Cannot open MCMC file: {}", inputFile);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  std::unique_ptr<TMacro> Config(TempFile->Get<TMacro>("MaCh3_Config"));
  if (!Config) {
    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", inputFile.c_str());
    TempFile->ls();
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  MACH3LOG_INFO("Loading YAML config from MCMC chain");
  YAML::Node Settings = TMacroToYAML(*Config);
  bool asimovfit = GetFromManager<bool>(Settings["General"]["Asimov"], false, __FILE__ , __LINE__);
  if (asimovfit) {
    MACH3LOG_WARN("MCMC chain was produced from an Asimov fit");
    MACH3LOG_WARN("ReweightMCMC does not currently handle Asimov shifting, results may be incorrect!");
  } else {
    MACH3LOG_INFO("Not an Asimov fit, proceeding with reweighting");
  }

  // Open input file and get tree
  auto inFile = std::unique_ptr<TFile>(TFile::Open(inputFile.c_str(), "READ"));
  if (!inFile || inFile->IsZombie()) {
    MACH3LOG_ERROR("Cannot open input file: {}", inputFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::unique_ptr<TTree> inTree(inFile->Get<TTree>("posteriors"));
  if (!inTree) {
    MACH3LOG_ERROR("Cannot find 'posteriors' tree in input file");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Create output file
  std::string configString = configFile.substr(configFile.find_last_of('/') + 1, configFile.find_last_of('.') - configFile.find_last_of('/') - 1);
  std::string outputFile = inputFile.substr(0, inputFile.find_last_of('.')) + "_reweighted_" + configString + ".root";
  auto outFile = std::unique_ptr<TFile>(TFile::Open(outputFile.c_str(), "RECREATE"));
  if (!outFile || outFile->IsZombie()) {
    MACH3LOG_ERROR("Cannot create output file: {}", outputFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  MACH3LOG_INFO("Output file will be: {}", outputFile);

  // Copy all the remaining objects into the out file (i.e. all but posteriors tree)
  TIter next(inFile->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next())) {
    inFile->cd();
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      // It's a folder, create and copy its contents
      TDirectory* srcDir = static_cast<TDirectory*>(obj.get());
      TDirectory* destDir = outFile->mkdir(srcDir->GetName());
      TIter nextSubKey(srcDir->GetListOfKeys());
      while (TKey* subKey = dynamic_cast<TKey*>(nextSubKey())) {
        srcDir->cd();
        std::unique_ptr<TObject> subObj(subKey->ReadObj());
        destDir->cd();
        subObj->Write();
      }
    } else if (std::string(key->GetName()) != "posteriors") {
      // Regular object, skip "posteriors" tree
      outFile->cd();
      obj->Write();
    }
  }

  // Clone the tree structure
  outFile->cd();
  std::unique_ptr<TTree> outTree(inTree->CloneTree(0));

  // Set up parameter reading
  std::map<std::string, double> paramValues;
  for (const auto& rwConfig : reweightConfigs) {
    for (const auto& paramName : rwConfig.paramNames) {
      if (paramValues.find(paramName) == paramValues.end()) {
        paramValues[paramName] = 0.0;
        // KS: Params other than Osc have branch like Param_, so we need to match it
        auto idx = processor->GetParamIndexFromName(paramName);
        auto BranchName = processor->GetBranchNames()[idx];
        inTree->SetBranchAddress(BranchName, &paramValues[paramName]);
      }
    }
  }

  // Add weight branches
  std::map<std::string, double> weights;
  std::map<std::string, TBranch*> weightBranches;

  for (const auto& rwConfig : reweightConfigs) {
    weights[rwConfig.weightBranchName] = 1.0;
    weightBranches[rwConfig.weightBranchName] = outTree->Branch(
        rwConfig.weightBranchName.c_str(),
        &weights[rwConfig.weightBranchName],
        (rwConfig.weightBranchName + "/D").c_str());
    MACH3LOG_INFO("Added weight branch: {}", rwConfig.weightBranchName);
  }

  // For 2D reweight and non-gaussian (ie TGraph) 1D reweight we need to do it ourselves
  // Process all entries
  Long64_t nEntries = inTree->GetEntries();
  MACH3LOG_INFO("Processing {} entries", nEntries);

  /// @todo add tracking for how many events are outside the graph ranges for diagnostics DWR
  for (Long64_t i = 0; i < nEntries; ++i) {
    if(i % (nEntries/20) == 0) M3::Utils::PrintProgressBar(i, nEntries);

    inTree->GetEntry(i);

    // Calculate weights for all configurations
    for (const auto& rwConfig : reweightConfigs) {
      double weight = 1.0;
      if (rwConfig.dimension == 1) {
        weight = Get1DWeight(rwConfig, paramValues);
      } else if (rwConfig.dimension == 2) {
        weight = Get2DWeight(rwConfig, paramValues);
      }
      weights[rwConfig.weightBranchName] = weight;
    }
    // Fill the output tree
    outTree->Fill();
  } // end loop over entries

  // Write and close
  outFile->cd();
  outTree->Write();

  // once we have finished the reweight save its configuration (reweightConfigNode) to the root file as a macro
  TMacro reweightMacro;
  reweightMacro.SetName("Reweight_Config");
  reweightMacro.SetTitle("ReweightMCMC configuration");
  std::stringstream ss;
  ss << reweight_settings;
  reweightMacro.AddLine(ss.str().c_str());
  reweightMacro.Write();

  MACH3LOG_INFO("Reweighting completed successfully!");
  MACH3LOG_INFO("Final reweighted file is: {}", outputFile);
}

/// @brief Main function
int main(int argc, char *argv[])
{
  SetMaCh3LoggerFormat();

  if (argc != 3) {
    MACH3LOG_ERROR("How to use: {} <config.yaml> <input_file.root>", argv[0]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::string configFile = argv[1];
  std::string inputFile = argv[2];

  ReweightMCMC(configFile, inputFile);

  return 0;
}
