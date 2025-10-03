//MaCh3 includes
#include "Manager/Manager.h"
#include "Fitters/MCMCProcessor.h"

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
   
    // For TGraph 1D or 2D
    std::string fileName;
    std::string graphName;

    // For TGraph1D
    std::unique_ptr<TGraph> graph_1D;

    // For TGraph2D
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

/// @brief Get parameter information from MCMCProcessor
bool GetParameterInfo(MCMCProcessor* processor, const std::string& paramName, 
                     double& mean, double& sigma);

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

void ReweightMCMC(const std::string& configFile, const std::string& inputFile)
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
        
        // Check if this particular reweight is enabled !!! Curently only support one reweight at a time so this defaults to enabled
        if (!GetFromManager<bool>(reweightConfigNode["Enabled"], true)) {
           MACH3LOG_INFO("Skipping disabled reweight: {}", reweightKey);
           continue;
        }
        
        ReweightConfig reweightConfig;
        reweightConfig.key = reweightKey;
        reweightConfig.name = GetFromManager<std::string>(reweightConfigNode["ReweightName"], reweightKey);
        reweightConfig.type = GetFromManager<std::string>(reweightConfigNode["ReweightType"], "Gaussian");
        reweightConfig.dimension = GetFromManager<int>(reweightConfigNode["ReweightDim"], 1);
        // reweightConfig.weightBranchName = "Weight_" + reweightKey; // for now all weights will be stored in branches called Weight until multi weight support
        reweightConfig.weightBranchName = "Weight";
        reweightConfig.enabled = true;
        
        // Handle different reweight types as they fill different members
        if (reweightConfig.dimension == 1) {
            if (reweightConfig.type == "Gaussian") {
                // For Gaussian reweights, we need the parameter name and prior values (mean, sigma)
                std::string paramName = GetFromManager<std::string>(reweightConfigNode["ReweightVar"], "");
                auto priorValues = GetFromManager<std::vector<double>>(reweightConfigNode["ReweightPrior"], {});
                
                reweightConfig.paramNames = {paramName};
                reweightConfig.priorValues = priorValues;

                if (paramName.empty() || priorValues.size() != 2) {
                    MACH3LOG_ERROR("Invalid Gaussian reweight configuration for {}", reweightKey);
                    continue;
                }
            } else if (reweightConfig.type == "TGraph") {
                // For TGraph reweights, we need the parameter name and the TGraph file and name
                std::string paramName = GetFromManager<std::string>(reweightConfigNode["ReweightVar"], "");
                std::string fileName = GetFromManager<std::string>(reweightConfigNode["ReweightPrior"]["file"], "");
                std::string graphName = GetFromManager<std::string>(reweightConfigNode["ReweightPrior"]["graph_name"], ""); 
                reweightConfig.paramNames = {paramName};
                reweightConfig.fileName = fileName;
                reweightConfig.graphName = graphName;

                if (paramName.empty() || fileName.empty() || graphName.empty()) {
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
                auto graph = constraintFile->Get<TGraph>(reweightConfig.graphName.c_str());
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
        MACH3LOG_ERROR("No valid reweight configurations found in config file");
        throw MaCh3Exception(__FILE__, __LINE__);
    } else if (reweightConfigs.size() > 1) {    // check number of ReweightConfigs, currently maximum supported is 1 due to structure of ProcessMCMC
        MACH3LOG_ERROR("Currently only one reweight configuration is supported at a time, found {}", reweight_settings.size());
        throw MaCh3Exception(__FILE__, __LINE__);
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
   
    // TODO Finish Asimov shifting implementation, for now just warn that Asimovs are not being properly handled
    // Get the settings for the MCMC
    TFile *TempFile = TFile::Open(inputFile.c_str(), "READ");
    if (TempFile == nullptr || TempFile->IsZombie()) {
        MACH3LOG_ERROR("Cannot open MCMC file: {}", inputFile);
        throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    TMacro *Config = TempFile->Get<TMacro>("MaCh3_Config");
    if (Config == nullptr) {
        MACH3LOG_ERROR("Didn't find MaCh3_Config tree in MCMC file! {}", inputFile.c_str());
        TempFile->ls();
        throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    MACH3LOG_INFO("Loading YAML config from MCMC chain");
    YAML::Node Settings = TMacroToYAML(*Config);
    bool asimovfit = GetFromManager<bool>(Settings["General"]["Asimov"], false);
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
    std::string outputFile = inputFile.substr(0, inputFile.find_last_of('.')) + "_reweighted.root";
    auto outFile = std::unique_ptr<TFile>(TFile::Open(outputFile.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        MACH3LOG_ERROR("Cannot create output file: {}", outputFile);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    MACH3LOG_INFO("Output file will be: {}", outputFile);
    
    // Copy all the remaining objects into the out file (ie all but posteriors tree)

    TIter next(inFile->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next())) {
        inFile->cd();
        TObject* obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            // It's a folder, create and copy its contents
            TDirectory* srcDir = static_cast<TDirectory*>(obj);
            TDirectory* destDir = outFile->mkdir(srcDir->GetName());
            TIter nextSubKey(srcDir->GetListOfKeys());
            while (TKey* subKey = dynamic_cast<TKey*>(nextSubKey())) {
                srcDir->cd();
                TObject* subObj = subKey->ReadObj();
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
    
    for (const auto& rwConfig : reweightConfigs) {
        weights[rwConfig.weightBranchName] = 1.0;
        weightBranches[rwConfig.weightBranchName] = outTree->Branch(
            rwConfig.weightBranchName.c_str(), 
            &weights[rwConfig.weightBranchName], 
            (rwConfig.weightBranchName + "/D").c_str()
        );
        MACH3LOG_INFO("Added weight branch: {}", rwConfig.weightBranchName);
    }
   
    // If a given reweight is 1D Gaussian we can just let MCMCProcessor method do the reweight
    for (const auto& rwConfig : reweightConfigs){
        if (rwConfig.dimension == 1 && rwConfig.type == "Gaussian"){
            // case the rwConfig parameters to the specific format processor needs
            const std::vector<std::string>& paramName = {rwConfig.paramNames[0]};
            const std::vector<double>& priorCentral = {rwConfig.priorValues[0]};
            const std::vector<double>& priorSigma = {rwConfig.priorValues[1]};
            processor->ReweightPrior(paramName, priorCentral, priorSigma);
            MACH3LOG_INFO("Applied Gaussian reweighting for {} with mean={} and sigma={}", paramName[0], priorCentral[0], priorSigma[0]);
        }
    }
    // For 2D reweight and non-gaussian (ie TGraph) 1D reweight we need to do it ourselves
    // Process all entries
    Long64_t nEntries = inTree->GetEntries();
    MACH3LOG_INFO("Processing {} entries", nEntries);
    
    // TODO: add tracking for how many events are outside the graph ranges for diagnostics DWR

    for (Long64_t i = 0; i < nEntries; ++i) {
        if(i % 10000 == 0) MaCh3Utils::PrintProgressBar(i, nEntries);
       
        inTree->GetEntry(i);
        
        // Calculate weights for all configurations
        for (const auto& rwConfig : reweightConfigs) {
            double weight = 1.0;
            
            if (rwConfig.dimension == 1 && rwConfig.type != "Gaussian") {
                if (rwConfig.type == "TGraph") {
                        double paramValue = paramValues[rwConfig.paramNames[0]];
                        weight = Graph_interpolate1D(nullptr, paramValue); // TODO replace nullptr with actual TGraph pointer when implemented
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
            }
            
            weights[rwConfig.weightBranchName] = weight;
        }
        
        // Fill the output tree
        outTree->Fill();
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
