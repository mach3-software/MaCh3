//MaCh3 includes
#include "manager/manager.h"
#include "mcmc/MCMCProcessor.h"

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

// Structure to hold reweight configuration
struct ReweightConfig {
    std::string key;       // The YAML key for this reweight
    std::string name;
    std::string type;  // "Gaussian", "TGraph2D"
    int dimension;     // 1 or 2
    std::vector<std::string> paramNames;
    std::vector<double> priorValues;
    std::string weightBranchName;
    bool enabled;
    
    // For 2D TGraph2D reweighting
    std::string fileName;
    std::string graphName;
    std::string hierarchyType; // "NO", "IO", or "auto"
    std::unique_ptr<TGraph2D> graph_NO;
    std::unique_ptr<TGraph2D> graph_IO;
};

/// @brief Main executable responsible for reweighting MCMC chains
/// @param inputFile MCMC Chain file path
/// @param configFile Config file with reweighting settings
/// @author David Riley
/// @author Evan Goodman

void ReweightMCMC(const std::string& inputFile, const std::string& configFile);

// TODO: Maybe add a generic 2D reweight that is not dm32 and theta13 specific? DWR

/// @brief Function to interpolate 2D graph for Normal Ordering
double Graph_interpolateNO(TGraph2D* graph, double theta13, double dm32);

/// @brief Function to interpolate 2D graph for Inverted Ordering  
double Graph_interpolateIO(TGraph2D* graph, double theta13, double dm32);

/// @brief Function to interpolate 1D graph
double Graph_interpolate1D(TGraph* graph, double theta13);

/// @brief Calculate Gaussian weight for a given parameter value
double CalculateGaussianWeight(double value, double newMean, double newSigma);

/// @brief Get parameter information from MCMCProcessor
bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName, 
                     double& mean, double& sigma);

/// @brief Main function
int main(int argc, char *argv[]) 
{
    SetMaCh3LoggerFormat();
    
    if (argc != 3) {
        MACH3LOG_ERROR("How to use: {} <input_file.root> <config.yaml>", argv[0]);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    std::string inputFile = argv[1];
    std::string configFile = argv[2];
    
    ReweightMCMC(inputFile, configFile);
    
    return 0;
}

void ReweightMCMC(const std::string& inputFile, const std::string& configFile)
{
    MACH3LOG_INFO("File for reweighting: {} with config {}", inputFile, configFile);

    // Load configuration
    YAML::Node reweight_yaml = M3OpenConfig(configFile);
    YAML::Node reweight_settings = reweight_yaml["ReweightMCMC"];

    // Parse all reweight configurations first
    std::vector<ReweightConfig> reweightConfigs;
   
    // iterate through the keys in the reweighting yaml creating and storing the ReweightConfig as we go
    for (const auto& reweight : reweight_settings) {
        const std::string& reweightKey = reweight.first.as<std::string>();
        const YAML::Node& reweightConfigNode = reweight.second;
        
        // Check if this particular reweight is enabled
        if (!GetFromManager<bool>(reweightConfigNode["Enabled"], false)) {
            MACH3LOG_INFO("Skipping disabled reweight: {}", reweightKey);
            continue;
        }
        
        ReweightConfig reweightConfig;
        reweightConfig.key = reweightKey;
        reweightConfig.name = GetFromManager<std::string>(reweightConfigNode["ReweightName"], reweightKey);
        reweightConfig.type = GetFromManager<std::string>(reweightConfigNode["ReweightType"], "Gaussian");
        reweightConfig.dimension = GetFromManager<int>(reweightConfigNode["ReweightDim"], 1);
        reweightConfig.weightBranchName = "Weight_" + reweightKey;
        reweightConfig.enabled = true;
        
        // Handle different reweight types as they fill different members
        if (reweightConfig.dimension == 1) {

            std::string paramName = GetFromManager<std::string>(reweightConfigNode["ReweightVar"], "");
            auto priorValues = GetFromManager<std::vector<double>>(reweightConfigNode["ReweightPrior"], {});
            
            if (paramName.empty() || priorValues.size() < 2) {
                MACH3LOG_ERROR("Invalid 1D reweight configuration for {}", reweightKey);
                continue;
            }
            
            reweightConfig.paramNames = {paramName};
            reweightConfig.priorValues = priorValues;
            
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
                    auto graph_NO = constraintFile->Get<TGraph2D>(graphName_NO.c_str());
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
                    auto graph_IO = constraintFile->Get<TGraph2D>(graphName_IO.c_str());
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
        }
        
        reweightConfigs.push_back(std::move(reweightConfig));
        MACH3LOG_INFO("Added reweight configuration: {} ({}D, type: {})", reweightConfigs.back().name, reweightConfigs.back().dimension, reweightConfigs.back().type);
    }
    
    if (reweightConfigs.empty()) {
        MACH3LOG_WARN("No valid reweight configurations found");
        return;
    }
    
    // Create MCMCProcessor to get parameter information 
    auto processor = std::make_unique<MCMCProcessor>(inputFile);
    processor->Initialise();
    
    // Validate that all required parameters exist in the chain TODO: Get list only of UNIQUE parameters, this is repeating unnecessarily DWR
    for (const auto& rwConfig : reweightConfigs) {
        for (const auto& paramName : rwConfig.paramNames) {
            int paramIndex = processor->GetParamIndexFromName(paramName);
            if (paramIndex == -1) {
                MACH3LOG_ERROR("Parameter {} not found in MCMC chain", paramName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            MACH3LOG_INFO("Parameter {} found in chain", paramName);
        }
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
    std::string outputFile = inputFile.substr(0, inputFile.find_last_of('.')) + "_reweighted.root";
    auto outFile = std::unique_ptr<TFile>(TFile::Open(outputFile.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        MACH3LOG_ERROR("Cannot create output file: {}", outputFile);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    MACH3LOG_INFO("Output file will be: {}", outputFile);
    
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
    
    // We should set up the branches for reweighting first and then iterate over the entries
    // Adding one branch at a time is inefficient as you have to copy and replace the same file over and over
    // We will work out all the weights we need to add first and add them event by event

    // Add weight branches
    std::map<std::string, double> weights;
    std::map<std::string, TBranch*> weightBranches;
    
    for (const auto& rwConfig : reweightConfigs) {
        weights[rwConfig.weightBranchName] = 1.0;
        weightBranches[rwConfig.weightBranchName] = outTree->Branch(
            rwConfig.weightBranchName.c_str(), 
            &weights[rwConfig.weightBranchName], 
            (rwConfig.weightBranchName + "/D").c_str()
        );
        MACH3LOG_INFO("Added weight branch: {}", rwConfig.weightBranchName);
    }
    
    // Process all entries
    Long64_t nEntries = inTree->GetEntries();
    MACH3LOG_INFO("Processing {} entries", nEntries);
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % (nEntries/10) == 0) {
            MACH3LOG_INFO("Processing entry {}/{}", i, nEntries);
        }
        
        inTree->GetEntry(i);
        
        // Calculate weights for all configurations
        for (const auto& rwConfig : reweightConfigs) {
            double weight = 1.0;
            
            if (rwConfig.dimension == 1) {
                const std::string& paramName = rwConfig.paramNames[0];
                double paramValue = paramValues[paramName];
                
                if (rwConfig.type == "Gaussian") {
                    // TODO : just use MCMCProcessors gaussian reweight DWR
                    weight = CalculateGaussianWeight(paramValue, rwConfig.priorValues[0], rwConfig.priorValues[1]);
                }
                
            } else if (rwConfig.dimension == 2) {
                if (rwConfig.type == "TGraph2D") {
                    double dm32 = paramValues[rwConfig.paramNames[0]];
                    double theta13 = paramValues[rwConfig.paramNames[1]];
                    std::cout << "Parameters for reweight " << rwConfig.key << ": " << rwConfig.paramNames[0] << "=" << dm32 << ", " << rwConfig.paramNames[1] << "=" << theta13 << std::endl;
                    if (dm32 > 0) {
                        // Normal Ordering
                        if (rwConfig.graph_NO) {
                            weight = Graph_interpolateNO(rwConfig.graph_NO.get(), theta13, dm32);
                            std::cout << "rwConfig.graph_NO min x: " << rwConfig.graph_NO->GetXmin() << ", max x: " << rwConfig.graph_NO->GetXmax() << std::endl;
                            std::cout << "rwConfig.graph_NO min y: " << rwConfig.graph_NO->GetYmin() << ", max y: " << rwConfig.graph_NO->GetYmax() << std::endl;
                            std::cout << "NO weight for theta13=" << theta13 << ", dm32=" << dm32 << " is " << weight << std::endl;
                        } else {
                            MACH3LOG_ERROR("NO graph not available for {}", rwConfig.key);
                            weight = 0.0;
                        }
                    } else {
                        // Inverted Ordering
                        if (rwConfig.graph_IO) {
                            weight = Graph_interpolateIO(rwConfig.graph_IO.get(), theta13, dm32);
                            std::cout << "IO weight for theta13=" << theta13 << ", dm32=" << dm32 << " is " << weight << std::endl;
                        } else {
                            MACH3LOG_ERROR("IO graph not available for {}", rwConfig.key);
                            weight = 0.0;
                        }
                    }
                }
            }
            
            weights[rwConfig.weightBranchName] = weight;
        }
        
        // Fill the output tree
        outTree->Fill();
    }
    
    // Write and close
    outFile->cd();
    outTree->Write();
    
    MACH3LOG_INFO("Reweighting completed successfully!");
    MACH3LOG_INFO("Final reweighted file is: {}", outputFile);
}

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

double Graph_interpolate1D(TGraph* graph, double theta13)
{
    // TODO implement TGraph interpolation for 1D
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

double CalculateGaussianWeight(double value, double newMean, double newSigma)
{
    // Calculate Gaussian prior weight
    // Weight = P_new(x) = exp(-(x-mu)^2/(2*sigma^2)) / sqrt(2*pi*sigma^2)
    // Since we're often going from flat priors, just calculate the new Gaussian probability
    // TODO handle non flat priors
    
    if (newSigma <= 0) {
        MACH3LOG_ERROR("Invalid sigma value for Gaussian weight calculation");
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    double chi2 = std::pow((value - newMean) / newSigma, 2);
    double weight = std::exp(-0.5 * chi2);
    
    return weight;
}

bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName, 
                     double& mean, double& sigma)
{
    // Try to find the parameter index
    int paramIndex = processor->GetParamIndexFromName(paramName);
    
    if (paramIndex == -1) { // Assuming -1 or similar indicates not found
        return false;
    }
    
    // Get parameter information
    TString title;
    processor->GetNthParameter(paramIndex, mean, sigma, title);
    
    return true;
}
