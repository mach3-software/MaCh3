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

// C++ includes
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <map>
#include <chrono>
#include <omp.h>
_MaCh3_Safe_Include_End_ //}

/// @file ReweightMCMC.cpp
/// @brief This executable allow you to reweight MCMC Chain, such technique is used to study impact of different priors without rerunning MCMC
/// It can currently handle two types of reweights:
/// - 1D Gaussian prior reweighting on any parameter(s) in the chain, where the weight is calculated as the ratio of the new Gaussian prior to the original flat prior (which is just the value of the Gaussian at that point)
/// - 2D reweighting using a TGraph2D, currently set up to handle the specific case of applying a 2D reactor constraint in the (theta13, dm32) plane, where the weight is calculated as exp(-0.5 * chi^2) with chi^2 obtained from interpolating the TGraph2D at the point corresponding to the current values of theta13 and dm32 in the chain. The TGraph2D can contain separate surfaces for Normal and Inverted ordering and will automatically apply the correct one based on the sign of dm32 in each point in the chain.
/// @author David Riley
/// @author Evan Goodman

/// Structure to hold reweight configuration
struct ReweightConfig {
    std::string key;       // The YAML key for this reweight
    std::string name;
    std::string type;  // "Gaussian", "TGraph", "TGraph1D", "TGraph2D"
    int dimension;     // 1 or 2
    std::vector<std::string> paramNames;
    std::vector<std::vector<double>> priorValues; // Changed to handle multiple [mean, sigma] pairs
    std::string weightBranchName;
    bool enabled;
   
    // For TGraph 1D or 2D
    std::string fileName;
    std::string graphName;

    // For TGraph1D
    std::unique_ptr<TGraph> graph_1D;
    double graph1D_xmin = 0.0;
    double graph1D_xmax = 0.0;
    std::unique_ptr<TGraph> graph_NO_1D;
    std::unique_ptr<TGraph> graph_IO_1D;
    double graphNO_1D_xmin = 0.0;
    double graphNO_1D_xmax = 0.0;
    double graphIO_1D_xmin = 0.0;
    double graphIO_1D_xmax = 0.0;

    // For TGraph2D
    std::string hierarchyType; // "NO", "IO", or "auto"
    std::unique_ptr<TGraph2D> graph_NO;
    std::unique_ptr<TGraph2D> graph_IO;
    double graphNO_xmin = 0.0;
    double graphNO_xmax = 0.0;
    double graphNO_ymin = 0.0;
    double graphNO_ymax = 0.0;
    double graphIO_xmin = 0.0;
    double graphIO_xmax = 0.0;
    double graphIO_ymin = 0.0;
    double graphIO_ymax = 0.0;
};

/// @brief Main executable responsible for reweighting MCMC chains
/// @param configFile Config file with reweighting settings. Example in MaCh3Tutorial/TutorialConfigs/ReweightingConfigs/
/// @param inputFile MCMC Chain file path
/// @author David Riley
/// @author Evan Goodman

void ReweightMCMC(const std::string& configFile, const std::string& inputFile);
void ReweightMCMC(const std::string& configFile, const std::string& inputFile, bool reducedChain); // overload to handle temporary fix for reduced chains
/// @todo add a generic 2D reweight that is not dm32 and theta13 specific DWR

namespace {
bool gVerboseLogging = false;
}

/// @brief Function to interpolate 2D graph for Normal Ordering
double Graph_interpolateNO(TGraph2D* graph, double theta13, double dm32, double xmin, double xmax, double ymin, double ymax);

/// @brief Function to interpolate 2D graph for Inverted Ordering  
double Graph_interpolateIO(TGraph2D* graph, double theta13, double dm32, double xmin, double xmax, double ymin, double ymax);

/// @brief Function to interpolate 1D graph
double Graph_interpolate1D(TGraph* graph, double theta13, double xmin, double xmax);

/// @brief Get parameter information from MCMCProcessor
bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName, 
                     double& mean, double& sigma);

/// @brief Load reweighting setting like 1D or 2D from YAML config
void LoadReweightingSettings(std::vector<ReweightConfig>& reweightConfigs, const YAML::Node& reweight_settings);

/// @brief Main function
int main(int argc, char *argv[]) {
    SetMaCh3LoggerFormat();
    
    if (argc != 3 && argc != 4) {
        MACH3LOG_ERROR("How to use: {} <config.yaml> <input_file.root> [verbose]", argv[0]);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    std::string configFile = argv[1]; 
    std::string inputFile = argv[2];

    if (argc == 4) {
        std::string mode = argv[3];
        if (mode == "verbose") {
            gVerboseLogging = true;
            MACH3LOG_INFO("Verbose progress logging enabled (prints every 100000 entries)");
        } else {
            MACH3LOG_ERROR("Unknown optional argument: {}. Only 'verbose' is supported.", mode);
            MACH3LOG_ERROR("How to use: {} <config.yaml> <input_file.root> [verbose]", argv[0]);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
    }
    
    // Check input for posteriors vs osc_posteriors, fail gracefully if osc_posteriors since ProcessMCMC doesn't support reduced chain naming conventions
    auto tempFile = std::unique_ptr<TFile>(TFile::Open(inputFile.c_str(), "READ"));
    // if (tempFile->Get<TTree>("osc_posteriors")) {
    //     MACH3LOG_ERROR("Sorry, it seems like your posteriors have been reduced (TTree is named osc_posteriors) please use unreduced posterior files as MCMCProcessor can't handle the reduced naming conventions");
    //     throw MaCh3Exception(__FILE__, __LINE__);
    // } else if (!tempFile->Get<TTree>("posteriors")) {
    //     MACH3LOG_ERROR("Cannot find 'posteriors' tree in input file, found the following trees instead:");
    //     tempFile->ls();
    //     throw MaCh3Exception(__FILE__, __LINE__);
    // }

    // experimentally, accept osc_posteriors and try work around this for the purpose of reweighting OAR11B
    // short term fix, maybe there is a more elegant way to handle this in the long term, but for now just check if we have osc_posteriors and if so pass a flag to ReweightMCMC to handle this case in an overload
    if (tempFile->Get<TTree>("posteriors")) {
        ReweightMCMC(configFile, inputFile);
    } else if (tempFile->Get<TTree>("osc_posteriors")) {
        ReweightMCMC(configFile, inputFile, true);
    } else {
        MACH3LOG_ERROR("Cannot find 'posteriors' or 'osc_posteriors' tree in input file, found the following trees instead:");
        tempFile->ls();
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    return 0;
}

void LoadReweightingSettings(std::vector<ReweightConfig>& reweightConfigs, const YAML::Node& reweight_settings) {
    // Experimental mode: set true only if you intentionally want multi-reweight support.
    constexpr bool kEnableExperimentalMultiReweight = true;
    if (kEnableExperimentalMultiReweight) {
        MACH3LOG_WARN("EXPERIMENTAL multi-reweight mode is enabled. This path is not validated and will break your other plotting scripts");
    }

    // iterate through the keys in the reweighting yaml creating and storing the ReweightConfig as we go
    for (const auto& reweight : reweight_settings) {
        const std::string& reweightKey = reweight.first.as<std::string>();
        const YAML::Node& reweightConfigNode = reweight.second;

        // Check if this particular reweight is enabled !!! Currently only support one reweight at a time so this defaults to enabled
        if (!GetFromManager<bool>(reweightConfigNode["Enabled"], true)) {
            MACH3LOG_INFO("Skipping disabled reweight: {}", reweightKey);
            continue;
        }

        ReweightConfig reweightConfig;
        reweightConfig.key = reweightKey;
        reweightConfig.name = GetFromManager<std::string>(reweightConfigNode["ReweightName"], reweightKey);
        reweightConfig.type = GetFromManager<std::string>(reweightConfigNode["ReweightType"], "Gaussian");
        reweightConfig.dimension = GetFromManager<int>(reweightConfigNode["ReweightDim"], 1);
        reweightConfig.weightBranchName = kEnableExperimentalMultiReweight ? ("Weight_" + reweightKey) : "Weight";
        reweightConfig.enabled = true;

        // Handle different reweight types as they fill different members
        if (reweightConfig.dimension == 1) {
            if (reweightConfig.type == "Gaussian") {
                // For Gaussian reweights, we need the parameter name(s) and prior values (mean, sigma pairs)
                auto paramNames = GetFromManager<std::vector<std::string>>(reweightConfigNode["ReweightVar"], {});

                // Get prior values - handle both single [mean, sigma] pair and list of pairs for safety
                auto priorNode = reweightConfigNode["ReweightPrior"];
                std::vector<std::vector<double>> allPriorValues;

                if (priorNode.IsSequence() && priorNode.size() > 0) {
                    // Check if first element is a number (single [mean, sigma] pair) or sequence (list of pairs)
                    if (priorNode[0].IsScalar()) {
                        // Single [mean, sigma] pair - convert to list format
                        auto priorValues = GetFromManager<std::vector<double>>(priorNode, {});
                        if (priorValues.size() == 2) {
                            allPriorValues.push_back(priorValues);
                        }
                    } else {
                        // List of [mean, sigma] pairs
                        for (const auto& priorPair : priorNode) {
                            auto priorValues = GetFromManager<std::vector<double>>(priorPair, {});
                            if (priorValues.size() == 2) {
                                allPriorValues.push_back(priorValues);
                            }
                        }
                    }
                }

                reweightConfig.paramNames = paramNames;
                reweightConfig.priorValues = allPriorValues;

                if (paramNames.empty() || allPriorValues.empty() || paramNames.size() != allPriorValues.size()) {
                    MACH3LOG_ERROR("Invalid Gaussian reweight configuration for {}: {} parameters, {} prior pairs",
                                   reweightKey, paramNames.size(), allPriorValues.size());
                    continue;
                }
            } else if (reweightConfig.type == "TGraph" || reweightConfig.type == "TGraph1D") {
                // For TGraph reweights, we need the parameter name and the TGraph file and name
                auto paramNames = GetFromManager<std::vector<std::string>>(reweightConfigNode["ReweightVar"], {});
                auto priorNode = reweightConfigNode["ReweightPrior"];
                std::string fileName = GetFromManager<std::string>(priorNode["file"], "");
                std::string graphName = GetFromManager<std::string>(priorNode["graph_name"], "");
                std::string graphNameNO = GetFromManager<std::string>(priorNode["graph_name_no"], "");
                std::string graphNameIO = GetFromManager<std::string>(priorNode["graph_name_io"], "");
                reweightConfig.paramNames = paramNames;
                reweightConfig.fileName = fileName;
                reweightConfig.graphName = graphName;

                if (paramNames.empty() || paramNames.size() != 1 || fileName.empty()) {
                    MACH3LOG_ERROR("Invalid TGraph reweight configuration for {}", reweightKey);
                    continue;
                }

                // Load the 1D graph
                MACH3LOG_INFO("Loading 1D constraint from file: {}", reweightConfig.fileName);
                auto constraintFile = std::unique_ptr<TFile>(TFile::Open(reweightConfig.fileName.c_str(), "READ"));
                if (!constraintFile || constraintFile->IsZombie()) {
                    MACH3LOG_ERROR("Failed to open constraint file: {}", reweightConfig.fileName);
                    continue;
                }

                const bool useHierarchyGraphs = !graphNameNO.empty() && !graphNameIO.empty();
                if (useHierarchyGraphs) {
                    std::unique_ptr<TGraph> graphNO(constraintFile->Get<TGraph>(graphNameNO.c_str()));
                    std::unique_ptr<TGraph> graphIO(constraintFile->Get<TGraph>(graphNameIO.c_str()));

                    if (!graphNO || !graphIO) {
                        MACH3LOG_ERROR("Failed to load 1D NO/IO graphs: NO='{}', IO='{}'", graphNameNO, graphNameIO);
                        continue;
                    }

                    auto clonedGraphNO = static_cast<TGraph*>(graphNO->Clone());
                    clonedGraphNO->SetBit(kCanDelete, true);
                    reweightConfig.graph_NO_1D = std::unique_ptr<TGraph>(clonedGraphNO);
                    reweightConfig.graphNO_1D_xmin = 999999999;
                    reweightConfig.graphNO_1D_xmax = -999999999;
                    for (int i = 0; i < reweightConfig.graph_NO_1D->GetN(); ++i) {
                        const double x = reweightConfig.graph_NO_1D->GetX()[i];
                        if (x < reweightConfig.graphNO_1D_xmin) reweightConfig.graphNO_1D_xmin = x;
                        if (x > reweightConfig.graphNO_1D_xmax) reweightConfig.graphNO_1D_xmax = x;
                    }

                    auto clonedGraphIO = static_cast<TGraph*>(graphIO->Clone());
                    clonedGraphIO->SetBit(kCanDelete, true);
                    reweightConfig.graph_IO_1D = std::unique_ptr<TGraph>(clonedGraphIO);
                    reweightConfig.graphIO_1D_xmin = 999999999;
                    reweightConfig.graphIO_1D_xmax = -999999999;
                    for (int i = 0; i < reweightConfig.graph_IO_1D->GetN(); ++i) {
                        const double x = reweightConfig.graph_IO_1D->GetX()[i];
                        if (x < reweightConfig.graphIO_1D_xmin) reweightConfig.graphIO_1D_xmin = x;
                        if (x > reweightConfig.graphIO_1D_xmax) reweightConfig.graphIO_1D_xmax = x;
                    }

                    MACH3LOG_INFO("Loaded 1D NO graph: {}", graphNameNO);
                    MACH3LOG_INFO("Loaded 1D IO graph: {}", graphNameIO);
                } else {
                    if (graphName.empty()) {
                        MACH3LOG_ERROR("Invalid TGraph reweight configuration for {}: missing graph_name or graph_name_no/graph_name_io", reweightKey);
                        continue;
                    }

                    std::unique_ptr<TGraph> graph(constraintFile->Get<TGraph>(reweightConfig.graphName.c_str()));
                    if (graph) {
                        auto clonedGraph = static_cast<TGraph*>(graph->Clone());
                        clonedGraph->SetBit(kCanDelete, true);
                        reweightConfig.graph_1D = std::unique_ptr<TGraph>(clonedGraph);
                        reweightConfig.graph1D_xmin = 999999999;
                        reweightConfig.graph1D_xmax = -999999999;
                        for (int i = 0; i < reweightConfig.graph_1D->GetN(); ++i) {
                            const double x = reweightConfig.graph_1D->GetX()[i];
                            if (x < reweightConfig.graph1D_xmin) reweightConfig.graph1D_xmin = x;
                            if (x > reweightConfig.graph1D_xmax) reweightConfig.graph1D_xmax = x;
                        }
                        MACH3LOG_INFO("Loaded 1D graph: {}", reweightConfig.graphName);
                    } else {
                        MACH3LOG_ERROR("Failed to load graph: {}", reweightConfig.graphName);
                        continue;
                    }
                }
            } else {
                MACH3LOG_ERROR("Unknown 1D reweight type: {} for {}", reweightConfig.type, reweightKey);
                throw MaCh3Exception(__FILE__, __LINE__);
            }

        } else if (reweightConfig.dimension == 2) {
            auto paramNames = GetFromManager<std::vector<std::string>>(reweightConfigNode["ReweightVar"], {});

            // 2D reweights need 2 parameter names
            if (paramNames.size() != 2) {
                MACH3LOG_ERROR("2D reweighting requires exactly 2 parameter names for {}", reweightKey);
                continue;
            }

            reweightConfig.paramNames = paramNames;

            if (reweightConfig.type == "TGraph2D") {
                auto priorConfig = reweightConfigNode["ReweightPrior"];
                reweightConfig.fileName = GetFromManager<std::string>(priorConfig["file"], "");
                reweightConfig.graphName = GetFromManager<std::string>(priorConfig["graph_name"], "");
                reweightConfig.hierarchyType = GetFromManager<std::string>(priorConfig["hierarchy"], "auto");

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
                        reweightConfig.graphNO_xmin = reweightConfig.graph_NO->GetXmin();
                        reweightConfig.graphNO_xmax = reweightConfig.graph_NO->GetXmax();
                        reweightConfig.graphNO_ymin = reweightConfig.graph_NO->GetYmin();
                        reweightConfig.graphNO_ymax = reweightConfig.graph_NO->GetYmax();
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
                        reweightConfig.graphIO_xmin = reweightConfig.graph_IO->GetXmin();
                        reweightConfig.graphIO_xmax = reweightConfig.graph_IO->GetXmax();
                        reweightConfig.graphIO_ymin = reweightConfig.graph_IO->GetYmin();
                        reweightConfig.graphIO_ymax = reweightConfig.graph_IO->GetYmax();
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
        }

        reweightConfigs.push_back(std::move(reweightConfig));
        MACH3LOG_INFO("Added reweight configuration: {} ({}D, type: {})", reweightConfigs.back().name, reweightConfigs.back().dimension, reweightConfigs.back().type);
    }

    if (reweightConfigs.empty()) {
        MACH3LOG_ERROR("No valid reweight configurations found in config file");
        throw MaCh3Exception(__FILE__, __LINE__);
    } else if (!kEnableExperimentalMultiReweight && reweightConfigs.size() > 1) {
        MACH3LOG_ERROR("Currently only one reweight configuration is supported at a time, found {}", reweightConfigs.size());
        MACH3LOG_ERROR("To try multi-reweight anyway, set kEnableExperimentalMultiReweight=true in ReweightMCMC.cpp and recompile MaCh3.");
        throw MaCh3Exception(__FILE__, __LINE__);
    } else if (kEnableExperimentalMultiReweight && reweightConfigs.size() > 1) {
        MACH3LOG_WARN("Proceeding with {} reweights in EXPERIMENTAL mode.", reweightConfigs.size());
        MACH3LOG_WARN("Check your yaml to save as combined or separate weight branches");
    }
}

void ReweightMCMC(const std::string& configFile, const std::string& inputFile){
    MACH3LOG_INFO("File for reweighting: {} with config {}", inputFile, configFile);
    // Load configuration
    YAML::Node reweight_yaml = M3OpenConfig(configFile);
    YAML::Node reweight_settings = reweight_yaml["ReweightMCMC"]["WeightConfigs"];
    YAML::Node general_settings = reweight_yaml["ReweightMCMC"]["Settings"];

   
    std::string yaml_reweight_dump = YAML::Dump(reweight_settings); // print the config to the log for transparency
    std::cout << "Loaded reweighting configuration:\n" << yaml_reweight_dump << std::endl;
    std::string yaml_general_dump = YAML::Dump(general_settings); // print the config to the log for transparency
    std::cout << "Loaded general settings configuration:\n" << yaml_general_dump << std::endl;

    
    bool combineWeights = GetFromManager<bool>(general_settings["CombineReweights"], false);
    if (combineWeights) {
        MACH3LOG_WARN("CombineWeights is enabled, final weight will be product of all individual weights. This is EXPERIMENTAL and may not be validated, use with caution.");
    } else {
        MACH3LOG_WARN("CombineWeights is disabled, all weights will appear in the output file and this may break your plotting/processing scripts!");
    }

    // Parse all reweight configurations first
    std::vector<ReweightConfig> reweightConfigs;
   
    LoadReweightingSettings(reweightConfigs, reweight_settings);

    // Create MCMCProcessor to get parameter information 
    auto processor = std::make_unique<MCMCProcessor>(inputFile);
    processor->Initialise();
    
    // Validate that all required parameters exist in the chain 
    /// @todo Get list only of unique parameters, this is repeating unnecessarily when adding more than 1 weight DWR
    for (const auto& rwConfig : reweightConfigs) {
        for (const auto& paramName : rwConfig.paramNames) {
            int paramIndex = processor->GetParamIndexFromName(paramName);
            if (paramIndex == M3::_BAD_INT_) {
                MACH3LOG_ERROR("Parameter {} not found in MCMC chain", paramName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            MACH3LOG_INFO("Parameter {} found in chain", paramName);
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
    bool asimovfit = GetFromManager<bool>(Settings["General"]["Asimov"], false);
    if (asimovfit) {
        MACH3LOG_WARN("MCMC chain was produced from an Asimov fit");
        MACH3LOG_WARN("ReweightMCMC does not currently handle Asimov shifting, results may be incorrect depending on what you want from me!");
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
                inTree->SetBranchAddress(paramName.c_str(), &paramValues[paramName]);
            }
        }
    }

    // Add weight branches
    std::map<std::string, double> weights;
    std::map<std::string, TBranch*> weightBranches;
    if (combineWeights){
        weights["CombinedWeight"] = 1.0;
        weightBranches["CombinedWeight"] = outTree->Branch("Weight", &weights["CombinedWeight"], "Weight/D");
        MACH3LOG_INFO("Added combined weight branch: Weight");
    } else {
        for (const auto& rwConfig : reweightConfigs) {
            weights[rwConfig.weightBranchName] = 1.0;
            weightBranches[rwConfig.weightBranchName] = outTree->Branch(
                rwConfig.weightBranchName.c_str(), 
                &weights[rwConfig.weightBranchName], 
                (rwConfig.weightBranchName + "/D").c_str()
            );
            MACH3LOG_INFO("Added weight branch: {}", rwConfig.weightBranchName);
        }
    }

    bool processMCMCreweighted = false;

    // Offload only the single 1D Gaussian case to MCMCProcessor.
    // Any mixed/multiple reweight setup is handled locally so all branches stay in one tree.
    const bool useProcessMCMCForGaussian = (reweightConfigs.size() == 1 && reweightConfigs[0].dimension == 1 && reweightConfigs[0].type == "Gaussian");

    if (useProcessMCMCForGaussian) {
        const auto& rwConfig = reweightConfigs[0];
        const std::vector<std::string>& paramNames = rwConfig.paramNames;
        std::vector<double> priorCentral;
        std::vector<double> priorSigma;

        for (const auto& priorPair : rwConfig.priorValues) {
            priorCentral.push_back(priorPair[0]);
            priorSigma.push_back(priorPair[1]);
        }

        processor->ReweightPrior(paramNames, priorCentral, priorSigma);
        MACH3LOG_INFO("Applied Gaussian reweighting with MCMCProcessor for {} parameters", paramNames.size());
        for (size_t i = 0; i < paramNames.size(); ++i) {
            MACH3LOG_INFO("  {}: mean={}, sigma={}", paramNames[i], priorCentral[i], priorSigma[i]);
        }
        processMCMCreweighted = true;
    }

    // Process all entries for local reweighting path.
    // Process all entries
    Long64_t nEntries = inTree->GetEntries();
    const Long64_t progressStep = (nEntries >= 20) ? (nEntries / 20) : 1;
    MACH3LOG_INFO("Processing {} entries", nEntries);
    

    /// @todo add tracking for how many events are outside the graph ranges for diagnostics DWR
    
    if (processMCMCreweighted) {
        MACH3LOG_INFO("MCMCProcessor has reweighted, skipping duplicate reweighting");
    } else {
        for (Long64_t i = 0; i < nEntries; ++i) {
            if (gVerboseLogging && nEntries > 100000 && (i % 100000 == 0)) {
                MACH3LOG_INFO("Reweight progress (posteriors): entry {}/{}", i, nEntries);
            }
            if (i % progressStep == 0) MaCh3Utils::PrintProgressBar(i, nEntries);
            inTree->GetEntry(i);

            if (combineWeights) {
                weights["CombinedWeight"] = 1.0; // reset combined weight for this entry
            }

            // Calculate weights for all configurations
            for (const auto& rwConfig : reweightConfigs) {
                double weight = 1.0;

                if (rwConfig.dimension == 1 && rwConfig.type == "Gaussian") {
                    for (size_t j = 0; j < rwConfig.paramNames.size(); ++j) {
                        const std::string& paramName = rwConfig.paramNames[j];
                        const double paramValue = paramValues[paramName];
                        const double newCentral = rwConfig.priorValues[j][0];
                        const double newError = rwConfig.priorValues[j][1];
                        if (newError <= 0.0) {
                            MACH3LOG_ERROR("Invalid Gaussian sigma={} for parameter {}", newError, paramName);
                            throw MaCh3Exception(__FILE__, __LINE__);
                        }

                        const double newChi = (paramValue - newCentral) / newError;
                        const double newPrior = std::exp(-0.5 * newChi * newChi);
                        const double oldPrior = 1.0;
                        weight *= (newPrior / oldPrior);
                    }
                } else if (rwConfig.dimension == 1 && rwConfig.type != "Gaussian") {
                    if (rwConfig.type == "TGraph" || rwConfig.type == "TGraph1D") {
                        double paramValue = paramValues[rwConfig.paramNames[0]];
                        if (rwConfig.graph_NO_1D && rwConfig.graph_IO_1D) {
                            if (paramValue >= 0.0) {
                                weight = Graph_interpolate1D(rwConfig.graph_NO_1D.get(), paramValue, rwConfig.graphNO_1D_xmin, rwConfig.graphNO_1D_xmax);
                            } else {
                                weight = Graph_interpolate1D(rwConfig.graph_IO_1D.get(), paramValue, rwConfig.graphIO_1D_xmin, rwConfig.graphIO_1D_xmax);
                            }
                        } else {
                            weight = Graph_interpolate1D(rwConfig.graph_1D.get(), paramValue, rwConfig.graph1D_xmin, rwConfig.graph1D_xmax);
                        }
                    } else {
                        MACH3LOG_ERROR("Unsupported 1D reweight type: {} for {}", rwConfig.type, rwConfig.key);
                    }
                } else if (rwConfig.dimension == 2) {
                    if (rwConfig.type == "TGraph2D") {
                        double dm32 = paramValues[rwConfig.paramNames[0]];
                        double theta13 = paramValues[rwConfig.paramNames[1]];
                        if (dm32 > 0) {
                            // Normal Ordering
                            if (rwConfig.graph_NO) {
                                weight = Graph_interpolateNO(rwConfig.graph_NO.get(), theta13, dm32, rwConfig.graphNO_xmin, rwConfig.graphNO_xmax, rwConfig.graphNO_ymin, rwConfig.graphNO_ymax);
                            } else {
                                MACH3LOG_ERROR("NO graph not available for {}", rwConfig.key);
                                weight = 0.0;
                            }
                        } else {
                            // Inverted Ordering
                            if (rwConfig.graph_IO) {
                                weight = Graph_interpolateIO(rwConfig.graph_IO.get(), theta13, dm32, rwConfig.graphIO_xmin, rwConfig.graphIO_xmax, rwConfig.graphIO_ymin, rwConfig.graphIO_ymax);
                            } else {
                                MACH3LOG_ERROR("IO graph not available for {}", rwConfig.key);
                                weight = 0.0;
                            }
                        }
                    }
                }

                if (combineWeights) {
                    weights["CombinedWeight"] *= weight;
                } else {
                    weights[rwConfig.weightBranchName] = weight;
                }
            }

            outTree->Fill();
        }
    }
    
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

    // finally give final diagnostics to reader, and handle file cleanup if MCMCProcessor already did the reweighting
    if (processMCMCreweighted) {
        MACH3LOG_INFO("MCMCProcessor reweighting applied, Final reweighted file is: {}_reweighted.root", inputFile.substr(0, inputFile.find_last_of('.')));
        // delete the file we just created since MCMCProcessor already did the reweighting and saved it to a file
        outTree.reset();
        outFile->Close();
        outFile.reset();

        if (std::remove(outputFile.c_str()) != 0) {
            MACH3LOG_ERROR("Error deleting temporary file: {}", outputFile);
        } else {
            MACH3LOG_INFO("Deleted temporary file: {}", outputFile);
        }
    } else {
        MACH3LOG_INFO("Reweighting completed successfully!");
        MACH3LOG_INFO("Final reweighted file is: {}", outputFile);
    }
}

void ReweightMCMC(const std::string& configFile, const std::string& inputFile, bool reducedChain){
    MACH3LOG_INFO("File for reweighting: {} with config {}", inputFile, configFile);
    if (!reducedChain) {
        ReweightMCMC(configFile, inputFile);
        return;
    }
    
    // Load configuration
    YAML::Node reweight_yaml = M3OpenConfig(configFile);
    YAML::Node reweight_settings = reweight_yaml["ReweightMCMC"]["WeightConfigs"];
    YAML::Node general_settings = reweight_yaml["ReweightMCMC"]["Settings"];

   
    std::string yaml_reweight_dump = YAML::Dump(reweight_settings); // print the config to the log for transparency
    std::cout << "Loaded reweighting configuration:\n" << yaml_reweight_dump << std::endl;
    std::string yaml_general_dump = YAML::Dump(general_settings); // print the config to the log for transparency
    std::cout << "Loaded general settings configuration:\n" << yaml_general_dump << std::endl;

    bool combineWeights = GetFromManager<bool>(general_settings["CombineReweights"], false);
    if (combineWeights) {
        MACH3LOG_WARN("CombineWeights is enabled, final weight will be product of all individual weights. This is EXPERIMENTAL and may not be validated, use with caution.");
    } else {
        MACH3LOG_WARN("CombineWeights is disabled, all weights will appear in the output file and this may break your plotting/processing scripts!");
    }

    // Parse all reweight configurations first
    std::vector<ReweightConfig> reweightConfigs;
   
    LoadReweightingSettings(reweightConfigs, reweight_settings);

    // Create MCMCProcessor to get parameter information 
    // auto processor = std::make_unique<MCMCProcessor>(inputFile);
    // processor->Initialise();

    // since we don't have access to MCMCProcessor for reduced chains, we will check the tree directly for params
    // Open input file and get tree
    auto inFile = std::unique_ptr<TFile>(TFile::Open(inputFile.c_str(), "READ"));
    if (!inFile || inFile->IsZombie()) {
        MACH3LOG_ERROR("Cannot open input file: {}", inputFile);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    std::unique_ptr<TTree> inTree(inFile->Get<TTree>("osc_posteriors"));
    if (!inTree) {
        MACH3LOG_ERROR("Cannot find 'osc_posteriors' tree in input file");
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // need to map MaCh3 style parameter names to reduced chain style names
    std::map<std::string, std::string> paramMapping;
    paramMapping["sin2th_23"] = "theta23";
    paramMapping["sin2th_13"] = "theta13";
    paramMapping["sin2th_12"] = "theta12";
    paramMapping["delm2_23"] = "dm23";
    paramMapping["delm2_12"] = "dm12";
    paramMapping["delta_cp"] = "dcp";

    // Validate that all required parameters exist in the chain 
    for (const auto& rwConfig : reweightConfigs) {
        for (const auto& paramName : rwConfig.paramNames) {
            const auto mapIt = paramMapping.find(paramName);
            if (mapIt == paramMapping.end()) {
                MACH3LOG_ERROR("No reduced-chain mapping found for parameter {}", paramName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            if (inTree->GetBranch(mapIt->second.c_str()) == nullptr) {
                MACH3LOG_ERROR("Parameter {} not found in MCMC chain", paramName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            MACH3LOG_INFO("Parameter {} found in chain", paramName);
        }
    }
    ///////// Reduced chains don't contain the fit configurations so we cant check if theyre asimov, just trust the user
    /// @todo Finish Asimov shifting implementation, for now just warn that Asimovs are not being properly handled
    // Get the settings for the MCMC
    //auto TempFile = std::unique_ptr<TFile>(TFile::Open(inputFile.c_str(), "READ"));
    //if (!TempFile || TempFile->IsZombie()) {
    //    MACH3LOG_ERROR("Cannot open MCMC file: {}", inputFile);
    //    throw MaCh3Exception(__FILE__ , __LINE__ );
    //}
    //std::unique_ptr<TMacro> Config(TempFile->Get<TMacro>("MaCh3_Config"));
    //if (!Config) {
    //    MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", inputFile.c_str());
    //    TempFile->ls();
    //    throw MaCh3Exception(__FILE__ , __LINE__ );
    //}
    //MACH3LOG_INFO("Loading YAML config from MCMC chain");
    //YAML::Node Settings = TMacroToYAML(*Config);
    //bool asimovfit = GetFromManager<bool>(Settings["General"]["Asimov"], false);
    //if (asimovfit) {
    //    MACH3LOG_WARN("MCMC chain was produced from an Asimov fit");
    //    MACH3LOG_WARN("ReweightMCMC does not currently handle Asimov shifting, results may be incorrect!");
    //} else {
    //    MACH3LOG_INFO("Not an Asimov fit, proceeding with reweighting");
    //}
    
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
        } else if (std::string(key->GetName()) != "osc_posteriors") {
            // Regular object, skip "osc_posteriors" tree
            outFile->cd();
            obj->Write();
        }
    }

    // Clone the tree structure
    outFile->cd();
    std::unique_ptr<TTree> outTree(inTree->CloneTree(0));
    
    // Set up parameter reading using the mapped parameter names for reduced chains
    std::map<std::string, double> paramValues;
    for (const auto& rwConfig : reweightConfigs) {
        for (const auto& paramName : rwConfig.paramNames) {
            const auto mapIt = paramMapping.find(paramName);
            if (mapIt == paramMapping.end()) {
                MACH3LOG_ERROR("No reduced-chain mapping found for parameter {}", paramName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            const std::string& reducedParamName = mapIt->second;
            if (paramValues.find(reducedParamName) == paramValues.end()) {
                paramValues[reducedParamName] = 0.0;
                inTree->SetBranchAddress(reducedParamName.c_str(), &paramValues[reducedParamName]);
            }
        }
    }

    // Add weight branches
    std::map<std::string, double> weights;
    std::map<std::string, TBranch*> weightBranches;
    if (combineWeights){
        weights["CombinedWeight"] = 1.0;
        weightBranches["CombinedWeight"] = outTree->Branch("Weight", &weights["CombinedWeight"], "Weight/D");
        MACH3LOG_INFO("Added combined weight branch: Weight");
    } else {
        for (const auto& rwConfig : reweightConfigs) {
            weights[rwConfig.weightBranchName] = 1.0;
            weightBranches[rwConfig.weightBranchName] = outTree->Branch(
                rwConfig.weightBranchName.c_str(), 
                &weights[rwConfig.weightBranchName], 
                (rwConfig.weightBranchName + "/D").c_str()
            );
            MACH3LOG_INFO("Added weight branch: {}", rwConfig.weightBranchName);
        }
    }
    
    // bool processMCMCreweighted=false;

    // // If a given reweight is 1D Gaussian we can just let MCMCProcessor method do the reweight
    // for (const auto& rwConfig : reweightConfigs){
    //     if (rwConfig.dimension == 1 && rwConfig.type == "Gaussian"){
    //         // Extract the parameter names and convert priorValues to the format processor needs
    //         const std::vector<std::string>& paramNames = rwConfig.paramNames;
    //         std::vector<double> priorCentral;
    //         std::vector<double> priorSigma;
            
    //         // Extract means and sigmas from the prior pairs
    //         for (const auto& priorPair : rwConfig.priorValues) {
    //             priorCentral.push_back(priorPair[0]); // mean
    //             priorSigma.push_back(priorPair[1]);   // sigma
    //         }
            
    //         processor->ReweightPrior(paramNames, priorCentral, priorSigma);
    //         MACH3LOG_INFO("Applied Gaussian reweighting for {} parameters", paramNames.size());
    //         for (size_t i = 0; i < paramNames.size(); ++i) {
    //             MACH3LOG_INFO("  {}: mean={}, sigma={}", paramNames[i], priorCentral[i], priorSigma[i]);
    //         }
    //         processMCMCreweighted=true;
    //     }
    // }

    // ProcessMCMC cannot handle the 1D rewighting on a reduced chain so we need to do it ourselves even for the 1D gaussian case.
    bool processMCMCreweighted=false;

    const auto totalTimingStart = std::chrono::steady_clock::now(); //REMOVE

    // For 2D reweight and non-gaussian (ie TGraph) 1D reweight we need to do it ourselves.
    // Stage 1: cache all required parameter values outside multithread block.
    // this need to be done in a chunked way, or youll run out of memory
    Long64_t nEntries = inTree->GetEntries();
    MACH3LOG_INFO("Processing {} entries", nEntries);
    const Long64_t progressStep = (nEntries >= 20) ? (nEntries / 20) : 1;
    const Long64_t cacheChunkSize = 10000000;
    MACH3LOG_INFO("Reduced-chain chunk size set to {} entries", cacheChunkSize);
    MACH3LOG_INFO("If your Jobs are \"Killed\" try decreasing this size!");

    // Validate Gaussian priors once before entering the threaded section.
    for (const auto& rwConfig : reweightConfigs) {
        if (rwConfig.dimension == 1 && rwConfig.type == "Gaussian") {
            for (size_t j = 0; j < rwConfig.priorValues.size(); ++j) {
                const double sigma = rwConfig.priorValues[j][1];
                if (sigma <= 0.0) {
                    MACH3LOG_ERROR("Invalid Gaussian sigma={} for parameter {}", sigma, rwConfig.paramNames[j]);
                    throw MaCh3Exception(__FILE__, __LINE__);
                }
            }
        }
    }

    // Warm up ROOT interpolation internals before entering parallel region.
    for (const auto& rwConfig : reweightConfigs) {
        if (rwConfig.dimension == 1 && (rwConfig.type == "TGraph" || rwConfig.type == "TGraph1D")) {
            if (rwConfig.graph_NO_1D && rwConfig.graph_IO_1D) {
                const double noMid = 0.5 * (rwConfig.graphNO_1D_xmin + rwConfig.graphNO_1D_xmax);
                const double ioMid = 0.5 * (rwConfig.graphIO_1D_xmin + rwConfig.graphIO_1D_xmax);
                (void)Graph_interpolate1D(rwConfig.graph_NO_1D.get(), noMid, rwConfig.graphNO_1D_xmin, rwConfig.graphNO_1D_xmax);
                (void)Graph_interpolate1D(rwConfig.graph_IO_1D.get(), ioMid, rwConfig.graphIO_1D_xmin, rwConfig.graphIO_1D_xmax);
            } else if (rwConfig.graph_1D && rwConfig.graph_1D->GetN() > 0) {
                (void)Graph_interpolate1D(rwConfig.graph_1D.get(), rwConfig.graph_1D->GetX()[0], rwConfig.graph1D_xmin, rwConfig.graph1D_xmax);
            }
        } else if (rwConfig.dimension == 2 && rwConfig.type == "TGraph2D") {
            if (rwConfig.graph_NO) {
                const double xmid = 0.5 * (rwConfig.graphNO_xmin + rwConfig.graphNO_xmax);
                const double ymid = 0.5 * (rwConfig.graphNO_ymin + rwConfig.graphNO_ymax);
                (void)Graph_interpolateNO(rwConfig.graph_NO.get(), xmid, ymid, rwConfig.graphNO_xmin, rwConfig.graphNO_xmax, rwConfig.graphNO_ymin, rwConfig.graphNO_ymax);
            }
            if (rwConfig.graph_IO) {
                const double xmid = 0.5 * (rwConfig.graphIO_xmin + rwConfig.graphIO_xmax);
                const double ymid = 0.5 * (rwConfig.graphIO_ymin + rwConfig.graphIO_ymax);
                (void)Graph_interpolateIO(rwConfig.graph_IO.get(), xmid, -ymid, rwConfig.graphIO_xmin, rwConfig.graphIO_xmax, rwConfig.graphIO_ymin, rwConfig.graphIO_ymax);
            }
        }
    }

    if (processMCMCreweighted) {
        MACH3LOG_INFO("MCMCProcessor has reweighted, skipping duplicate reweighting");
    } else {
        /// @todo add tracking for how many events are outside the graph ranges for diagnostics DWR
        for (Long64_t chunkStart = 0; chunkStart < nEntries; chunkStart += cacheChunkSize) {
            const auto chunkTimingStart = std::chrono::steady_clock::now(); //REMOVE
            const Long64_t chunkEnd = (chunkStart + cacheChunkSize < nEntries) ? (chunkStart + cacheChunkSize) : nEntries;
            const size_t chunkN = static_cast<size_t>(chunkEnd - chunkStart);

            std::map<std::string, std::vector<double>> cachedParamValues;
            for (const auto& kv : paramValues) {
                cachedParamValues.emplace(kv.first, std::vector<double>(chunkN, 0.0));
            }

            // Stage 1: cache this chunk serially (ROOT tree reads are not thread-safe).
            for (Long64_t i = chunkStart; i < chunkEnd; ++i) {
                if (gVerboseLogging && (i % 100000 == 0)) {
                    MACH3LOG_INFO("Caching progress (osc_posteriors): entry {}/{}", i, nEntries);
                }
                inTree->GetEntry(i);
                const size_t localIdx = static_cast<size_t>(i - chunkStart);
                for (auto& kv : cachedParamValues) {
                    kv.second[localIdx] = paramValues.at(kv.first);
                }
            }

            const auto cacheTimingEnd = std::chrono::steady_clock::now(); //REMOVE

            // Stage 2: compute weights for this chunk in parallel.
            std::vector<std::vector<double>> cachedWeights(reweightConfigs.size(), std::vector<double>(chunkN, 1.0));
            #pragma omp parallel for schedule(static)
            for (Long64_t local = 0; local < static_cast<Long64_t>(chunkN); ++local) {
                const Long64_t i = chunkStart + local;
                if (gVerboseLogging && (i % 100000 == 0)) {
                    #pragma omp critical
                    MACH3LOG_INFO("Threaded reweight progress (osc_posteriors): entry {}/{}", i, nEntries);
                }
                for (size_t cfgIdx = 0; cfgIdx < reweightConfigs.size(); ++cfgIdx) {
                    const auto& rwConfig = reweightConfigs[cfgIdx];
                    double weight = 1.0;

                    if (rwConfig.dimension == 1 && rwConfig.type == "Gaussian") {
                        // Match MCMCProcessor::ReweightPrior: product of new/old priors for each requested parameter.
                        for (size_t j = 0; j < rwConfig.paramNames.size(); ++j) {
                            const std::string& mach3ParamName = rwConfig.paramNames[j];
                            const auto mapIt = paramMapping.find(mach3ParamName);
                            if (mapIt == paramMapping.end()) {
                                weight = 0.0;
                                break;
                            }

                            const std::string& reducedParamName = mapIt->second;
                            const double paramValue = cachedParamValues.at(reducedParamName)[static_cast<size_t>(local)];
                            const double newCentral = rwConfig.priorValues[j][0];
                            const double newError = rwConfig.priorValues[j][1];
                            if (newError <= 0.0) {
                                weight = 0.0;
                                break;
                            }

                            const double newChi = (paramValue - newCentral) / newError;
                            const double newPrior = std::exp(-0.5 * newChi * newChi);
                            const double oldPrior = 1.0; // reduced-chain fallback assumes flat old prior.
                            weight *= (newPrior / oldPrior);
                        }
                    } else if (rwConfig.dimension == 1 && rwConfig.type != "Gaussian") {
                        if (rwConfig.type == "TGraph" || rwConfig.type == "TGraph1D") {
                            const auto mapIt = paramMapping.find(rwConfig.paramNames[0]);
                            if (mapIt != paramMapping.end()) {
                                const double paramValue = cachedParamValues.at(mapIt->second)[static_cast<size_t>(local)];
                                if (rwConfig.graph_NO_1D && rwConfig.graph_IO_1D) {
                                    if (paramValue >= 0.0) {
                                        weight = Graph_interpolate1D(rwConfig.graph_NO_1D.get(), paramValue, rwConfig.graphNO_1D_xmin, rwConfig.graphNO_1D_xmax);
                                    } else {
                                        weight = Graph_interpolate1D(rwConfig.graph_IO_1D.get(), paramValue, rwConfig.graphIO_1D_xmin, rwConfig.graphIO_1D_xmax);
                                    }
                                } else {
                                    weight = Graph_interpolate1D(rwConfig.graph_1D.get(), paramValue, rwConfig.graph1D_xmin, rwConfig.graph1D_xmax);
                                }
                            } else {
                                weight = 0.0;
                            }
                        } else {
                            weight = 0.0;
                        }
                    } else if (rwConfig.dimension == 2) {
                        if (rwConfig.type == "TGraph2D") {
                            const auto dmMapIt = paramMapping.find(rwConfig.paramNames[0]);
                            const auto thMapIt = paramMapping.find(rwConfig.paramNames[1]);
                            if (dmMapIt == paramMapping.end() || thMapIt == paramMapping.end()) {
                                weight = 0.0;
                            } else {
                                const double dm32 = cachedParamValues.at(dmMapIt->second)[static_cast<size_t>(local)];
                                const double theta13 = cachedParamValues.at(thMapIt->second)[static_cast<size_t>(local)];
                                if (dm32 > 0) {
                                    // Normal Ordering
                                    if (rwConfig.graph_NO) {
                                        weight = Graph_interpolateNO(rwConfig.graph_NO.get(), theta13, dm32, rwConfig.graphNO_xmin, rwConfig.graphNO_xmax, rwConfig.graphNO_ymin, rwConfig.graphNO_ymax);
                                    } else {
                                        weight = 0.0;
                                    }
                                } else {
                                    // Inverted Ordering
                                    if (rwConfig.graph_IO) {
                                        weight = Graph_interpolateIO(rwConfig.graph_IO.get(), theta13, dm32, rwConfig.graphIO_xmin, rwConfig.graphIO_xmax, rwConfig.graphIO_ymin, rwConfig.graphIO_ymax);
                                    } else {
                                        weight = 0.0;
                                    }
                                }
                            }
                        }
                    }

                    cachedWeights[cfgIdx][static_cast<size_t>(local)] = weight;
                }
            }

            const auto computeTimingEnd = std::chrono::steady_clock::now(); //REMOVE

            // Stage 3: fill this chunk serially.
            for (Long64_t i = chunkStart; i < chunkEnd; ++i) {
                if (gVerboseLogging && (i % 100000 == 0)) {
                    MACH3LOG_INFO("Fill progress (osc_posteriors): entry {}/{}", i, nEntries);
                }
                if (i % progressStep == 0) MaCh3Utils::PrintProgressBar(i, nEntries);
                inTree->GetEntry(i);
                const size_t localIdx = static_cast<size_t>(i - chunkStart);
                if (combineWeights) {
                    weights["CombinedWeight"] = 1.0;
                    for (size_t cfgIdx = 0; cfgIdx < reweightConfigs.size(); ++cfgIdx) {
                        weights["CombinedWeight"] *= cachedWeights[cfgIdx][localIdx];
                    }
                } else {
                    for (size_t cfgIdx = 0; cfgIdx < reweightConfigs.size(); ++cfgIdx) {
                        const auto& rwConfig = reweightConfigs[cfgIdx];
                        weights[rwConfig.weightBranchName] = cachedWeights[cfgIdx][localIdx];
                    }
                }
                outTree->Fill();
            }

            const auto fillTimingEnd = std::chrono::steady_clock::now(); //REMOVE
            const auto cacheMs = std::chrono::duration_cast<std::chrono::milliseconds>(cacheTimingEnd - chunkTimingStart).count();
            const auto computeMs = std::chrono::duration_cast<std::chrono::milliseconds>(computeTimingEnd - cacheTimingEnd).count();
            const auto fillMs = std::chrono::duration_cast<std::chrono::milliseconds>(fillTimingEnd - computeTimingEnd).count();
            MACH3LOG_INFO("Chunk {}-{} timing: cache={} ms, compute={} ms, fill={} ms", chunkStart, chunkEnd, cacheMs, computeMs, fillMs); //REMOVE
        }
    }
    MaCh3Utils::PrintProgressBar(nEntries, nEntries);

    const auto totalTimingEnd = std::chrono::steady_clock::now(); //REMOVE
    const auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(totalTimingEnd - totalTimingStart).count(); //REMOVE
    MACH3LOG_INFO("Reduced-chain reweight total wall time: {} ms", totalMs); //REMOVE
    
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

    // finally give final diagnostics to reader, and handle file cleanup if MCMCProcessor already did the reweighting 
    if (processMCMCreweighted){
        MACH3LOG_INFO("MCMCProcessor reweighting applied, Final reweighted file is: {}_reweighted.root", inputFile.substr(0, inputFile.find_last_of('.')));
        // delete the file we just created since MCMCProcessor already did the reweighting and saved it to a file
        outTree.reset(); // Release TTree before closing TFile
        outFile->Close();
        outFile.reset(); 
        
        if (std::remove(outputFile.c_str()) != 0) {
            MACH3LOG_ERROR("Error deleting temporary file: {}", outputFile);
        } else {
            MACH3LOG_INFO("Deleted temporary file: {}", outputFile);
        }
    } else {
        MACH3LOG_INFO("Reweighting completed successfully!");
        MACH3LOG_INFO("Final reweighted file is: {}", outputFile);
    }
}

double Graph_interpolateNO(TGraph2D* graph, double theta13, double dm32, double xmin, double xmax, double ymin, double ymax)
{
    if (!graph) {
        MACH3LOG_ERROR("Graph pointer is null");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    double chiSquared, prior; 
 
    if (theta13 < xmax && theta13 > xmin && dm32 < ymax && dm32 > ymin) {
        chiSquared = graph->Interpolate(theta13, dm32);
        prior = std::exp(-0.5 * chiSquared);
    } else {
        prior = 0.0;
    }
    
    return prior;
}

double Graph_interpolateIO(TGraph2D* graph, double theta13, double dm32, double xmin, double xmax, double ymin, double ymax)
{
    if (!graph) {
        MACH3LOG_ERROR("Graph pointer is null");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
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

double Graph_interpolate1D(TGraph* graph, double theta13, double xmin, double xmax)
{
    /// @todo double check implementation of TGraph interpolation for 1D
    if (!graph) {
        MACH3LOG_ERROR("Graph pointer is null");
        throw MaCh3Exception(__FILE__, __LINE__);
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

bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName, 
                     double& mean, double& sigma)
{
    // Try to find the parameter index
    int paramIndex = processor->GetParamIndexFromName(paramName);
    
    if (paramIndex == M3::_BAD_INT_) { // This indicate parameter not found
        return false;
    }
    
    // Get parameter information
    TString title;
    processor->GetNthParameter(paramIndex, mean, sigma, title);
    
    return true;
}
