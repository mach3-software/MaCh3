/// @file CombineModule.cpp
/// @brief Implementation of the CombineModule class
///
/// @author Clarence Wret
/// @author Kamil Skwarczynski
#include "Fitters/Processing/MCMCProcessor.h"
#include "Manager/Manager.h"
#include "CLI/Modules/CombineModule.hpp"

namespace M3{
    
    CombineModule::~CombineModule() = default;
    
    MaCh3ArgumentParser* CombineModule::get_parser(){
        m_parser = std::make_unique<MaCh3ArgumentParser>("combine", "1.0", argparse::default_arguments::help);
        m_parser->add_description("Tool for combining multiple chains or ROOT files into a single output.");
        m_parser->add_argument("-o", "--output")
        .help("Output ROOT file.")
        .store_into(OutFileName)
        .metavar("OUTPUT_FILE");
        m_parser->add_argument("-c", "--compression-level")
        .help("Compression level for the output ROOT file.")
        .store_into(targetCompression)
        .metavar("COMPRESSION_LEVEL")
        .default_value(1);
        m_parser->add_argument("-f", "--force")
        .help("Force overwrite of the output file if it exists.")
        .store_into(forceOverwrite)
        .flag();
        m_parser->add_argument("-m", "--merge-mode")
        .help("Merge files in-spite of differences")
        .store_into(forceMerge)
        .flag();
        m_parser->add_argument("input")
        .help("Input ROOT files.")
        .metavar("INPUT")
        .store_into(inpFileList)
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
        return m_parser.get();
    }
   
    /// @brief Main function for combining chains or ROOT files
    int CombineModule::Run() {
        SetMaCh3LoggerFormat();
        M3::Utils::MaCh3Welcome();
        if(OutFileName == ""){
            MACH3LOG_INFO("Using first file in list as output: ", inpFileList[0].c_str());
            OutFileName = inpFileList[0];
            inpFileList.erase(inpFileList.begin());
        }

        if(forceOverwrite){
            MACH3LOG_INFO("Will overwrite {} if it exists already", OutFileName.c_str());
        }
        MACH3LOG_INFO("Combining a total of {} files into {}", inpFileList.size(), OutFileName.c_str());
        CombineChain();
        return 0;
    }
    
    /// @brief KS: This allow us to skip output name etc in config. We expect Output name will be different but this doesn't invalidate chain merging
    bool CombineModule::ShouldSkipLine(const std::string& line, const std::vector<std::string>& SkipVector) {
        // Otherwise, check if the line contains any word from SkipVector
        for (const auto& word : SkipVector) {
            MACH3LOG_TRACE("{} : {}",line, word);
            if (line.find(word) != std::string::npos) {
                MACH3LOG_TRACE("Found matching word, therefore Skipping");
                return true;
            }
        }
        return false;
    }
    
    /// @brief make sure two configs are identical but skip specified fields. For example when comparing two chains nsteps or output name might be different and this is still fine to merge
    /// @param File1 Config from chain1
    /// @param File2 Config from chain2
    /// @param SkipVector Fields in yaml file to skip
    bool CombineModule::CompareTwoConfigs(const std::string& File1, const std::string& File2, const std::vector<std::string>& SkipVector) {
        std::istringstream file1(File1);
        std::istringstream file2(File2);
        
        std::string line1, line2;
        int lineNumber = 1;
        bool areEqual = true;
        
        while (std::getline(file1, line1) && std::getline(file2, line2)) {
            if (ShouldSkipLine(line1, SkipVector) || ShouldSkipLine(line2, SkipVector)) {
                ++lineNumber;
                continue;
            }
            if (line1 != line2) {
                areEqual = false;
                MACH3LOG_WARN("Difference found on line {}:", lineNumber);
                MACH3LOG_WARN("Config1: {}", line1);
                MACH3LOG_WARN("Config2: {}", line2);
            }
            ++lineNumber;
        }
        // Check if one file has extra lines
        while (std::getline(file1, line1)) {
            MACH3LOG_WARN("Extra line in {} on line {}: {}", File1, lineNumber, line1);
            ++lineNumber;
        }
        while (std::getline(file2, line2)) {
            MACH3LOG_WARN("Extra line in {} on line {}: {}", File2, lineNumber, line2);
            ++lineNumber;
        }
        return areEqual;
    }
    
    /// @brief EM: Will compare the version header contained in the two provided files and shout if they don't match
    bool CombineModule::checkSoftwareVersions(TFile *file, TFile *prevFile, const std::string& ConfigName, const std::vector<std::string>& SkipVector)
    {
        bool weirdFile = false;
        
        TMacro *versionHeader = file->Get<TMacro>(ConfigName.c_str());
        TMacro *prevVersionHeader = prevFile->Get<TMacro>(ConfigName.c_str());
        
        // EM: compare the digest of the version header file in this file, with the previous one
        if(!CompareTwoConfigs(TMacroToString(*versionHeader), TMacroToString(*prevVersionHeader), SkipVector)){
            MACH3LOG_ERROR("Looks like the {} embedded config for file {} is different to the previous ones", ConfigName, file->GetName());
            MACH3LOG_ERROR("This strongly suggests that this file was made with different software versions than the previous ones");
            weirdFile = true;
        }
        
        return weirdFile;
    }
    
    /// @brief When we merge two chains they have TDirectory ROOT didn't provide method for this so here we have this bad boy
    void CombineModule::CopyDir(TDirectory *source) {
        //copy all objects and subdirs of directory source as a subdir of the current directory
        source->ls();
        TDirectory *savdir = gDirectory;
        TDirectory *adir = savdir->Get<TDirectory>(source->GetName());
        // if directory doesn't exist make it
        if (!adir) {
            adir = savdir->mkdir(source->GetName());
        }
        adir->cd();
        //loop on all entries of this directory
        TKey *key;
        TIter nextkey(source->GetListOfKeys());
        while ((key = static_cast<TKey*>(nextkey()))) {
            const char *classname = key->GetClassName();
            TClass *cl = gROOT->GetClass(classname);
            if (!cl) continue;
            if (cl->InheritsFrom("TDirectory")) {
                source->cd(key->GetName());
                TDirectory *subdir = gDirectory;
                adir->cd();
                CopyDir(subdir);
                adir->cd();
            } else if (cl->InheritsFrom("TTree")) {
                TTree *T = source->Get<TTree>(key->GetName());
                adir->cd();
                TTree *newT = T->CloneTree();
                newT->Write();
            } else {
                source->cd();
                TObject *obj = key->ReadObj();
                adir->cd();
                obj->Write();
                delete obj;
            }
        }
        adir->SaveSelf(kTRUE);
        savdir->cd();
    }
    
    /// @brief Compare two histograms if they are identical
    /// @todo add checks for stuff like bin content etc
    bool CombineModule::CompareHistograms(const TH1* h1, const TH1* h2, const std::string& histName, const std::string& folderName)
    {
        if (!h1 || !h2) {
            MACH3LOG_ERROR("Null pointer passed to CompareHistograms for '{}'", histName);
            return false;
        }
        
        const double int1 = h1->Integral();
        const double int2 = h2->Integral();
        if (std::abs(int1 - int2) > 1e-6) {
            MACH3LOG_ERROR("Histogram '{}' in folder '{}' has different integrals: current = {}, previous = {}",
                histName, folderName, int1, int2);
                return false;
            }
            return true;
        }
        
        /// @brief Loop through TH1 and TMacro objects in FolderName in 'file' and compare with those in 'prevFile'
        bool CombineModule::CheckFolder(TFile* file, TFile* prevFile, const std::string& FolderName, const std::vector<std::string>& SkipVector)
        {
            bool mismatch = false;
            TDirectory* dir = file->GetDirectory(FolderName.c_str());
            TDirectory* prevDir = prevFile->GetDirectory(FolderName.c_str());
            
            if (!dir || !prevDir) {
                MACH3LOG_ERROR("Could not find folder '{}' in one or both files", FolderName);
                return true;
            }
            
            TIter nextKey(dir->GetListOfKeys());
            TKey* key;
            
            while ((key = static_cast<TKey*>(nextKey()))) {
                const std::string objName = key->GetName();
                TObject* obj = key->ReadObj();
                if (!obj) continue;
                
                // Handle TH1 comparison
                if (obj->InheritsFrom("TH1")) {
                    TH1* hist = static_cast<TH1*>(obj);
                    TH1* prevHist = dynamic_cast<TH1*>(prevDir->Get(objName.c_str()));
                    if (!prevHist) {
                        MACH3LOG_ERROR("Missing histogram '{}' in previous file (folder '{}')", objName, FolderName);
                        mismatch = true;
                        continue;
                    }
                    if (!CompareHistograms(hist, prevHist, objName, FolderName)) {
                        mismatch = true;
                    }
                }
                // Handle TMacro comparison
                else if (obj->InheritsFrom("TMacro")) {
                    TMacro* macro = static_cast<TMacro*>(obj);
                    TMacro* prevMacro = dynamic_cast<TMacro*>(prevDir->Get(objName.c_str()));
                    if (!prevMacro) {
                        MACH3LOG_ERROR("Missing TMacro '{}' in previous file (folder '{}')", objName, FolderName);
                        mismatch = true;
                        continue;
                    }
                    if (!CompareTwoConfigs(TMacroToString(*macro), TMacroToString(*prevMacro), SkipVector)) {
                        mismatch = true;
                    }
                }
            }
            return mismatch;
        }
        
        /// @brief custom function for merging TTree, should be similar to what HADD is using
        /// @warning KS: for some reason if "fast" is enable then I cannot open in ROOT5, no one should use R5 at this point..
        void CombineModule::FastMergeTTrees(const std::vector<std::string>& files, const std::string& outFile, const std::string& TTreeName) {
            TChain chain(TTreeName.c_str());
            for (const auto& f : files) chain.Add(f.c_str());
            
            TFile* outF = TFile::Open(outFile.c_str(), "UPDATE");
            
            TTree* newTree = chain.CloneTree(-1, "fast");
            newTree->SetName(TTreeName.c_str());
            outF->cd();
            newTree->Write("", TObject::kOverwrite);
            outF->Close();
            delete outF;
        }
        
        void CombineModule::CombineChain()
        {
            std::string outFileOption;
            if(forceOverwrite) outFileOption = "RECREATE";
            else outFileOption = "CREATE";
            
            TFile *prevFile = nullptr;
            
            // EM: loop through all the files in the provided list, compare the embedded version and config files
            //     If they match, we add the file to the list of files to be merged.
            //     If not, we throw an error and provide a (hopefully) helpful message telling the user why the files couldn't be merged.
            for(uint fileId = 0; fileId < inpFileList.size(); fileId++)
            {
                std::string fileName = inpFileList[fileId];
                TFile *file = new TFile(fileName.c_str());
                
                if(file->Get<TTree>("posteriors")->GetEntries() == 0){
                    MACH3LOG_WARN("Hmmm, file {} Doesn't seem to have any entries", fileName.c_str());
                    MACH3LOG_WARN("That's weird but I guess there's no rule that says a file can't be empty");
                    MACH3LOG_WARN("I'll skip it but maybe double check that this doesn't indicate some deeper problem");
                    continue;
                }
                
                // EM: need to set this in the initial case
                if(prevFile == nullptr) {
                    prevFile = file;
                }
                
                MACH3LOG_DEBUG("############ File {} #############", fileId);
                
                bool weirdFile = false;
                if(checkSoftwareVersions(file, prevFile, "MaCh3Engine/version_header")) weirdFile = true;
                if(checkSoftwareVersions(file, prevFile, "MaCh3_Config", {"OutputFile:", "NSteps:"})) weirdFile = true;
                if(CheckFolder(file, prevFile, "SampleFolder")) weirdFile = true;
                if(CheckFolder(file, prevFile, "CovarianceFolder")) weirdFile = true;
                
                if(weirdFile && !forceMerge){
                    MACH3LOG_ERROR("");
                    MACH3LOG_ERROR("=====================================================================================");
                    MACH3LOG_ERROR("This is not a great idea and could lead to weird outputs and cause some big headaches");
                    MACH3LOG_ERROR("further down the road. But if you reeeeally wanna do it and you know what you're");
                    MACH3LOG_ERROR("doing you can come here and remove the 'throw'");
                    MACH3LOG_ERROR("Or use -m option");
                    MACH3LOG_ERROR("{}:{}", __FILE__, __LINE__ + 2);
                    MACH3LOG_ERROR("=====================================================================================");
                    throw MaCh3Exception(__FILE__ , __LINE__ );
                }
                
                if(prevFile != file) {
                    prevFile->Close();
                    delete prevFile;
                }
                
                // EM: set these for the next iteration
                prevFile = file;
            }
            
            if (!forceOverwrite && access(OutFileName.c_str(), F_OK) != -1) {
                MACH3LOG_ERROR("Output file '{}' already exists. Use -f to force overwrite.", OutFileName);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
            //KS: Create new file
            TFile* outputFile = M3::Open(OutFileName, "recreate", __FILE__, __LINE__);
            outputFile->Close();
            delete outputFile;
            
            TStopwatch clock;
            clock.Start();
            
            MACH3LOG_INFO("Starting merging");
            FastMergeTTrees(inpFileList, OutFileName, "posteriors");
            FastMergeTTrees(inpFileList, OutFileName, "Settings");
            
            clock.Stop();
            MACH3LOG_INFO("Merging of took {:.2f}s to finish", clock.RealTime());
            
            //KS: Sadly we need to open file to save TDirectories to not have weird copy of several obejcts there...
            outputFile = M3::Open(OutFileName, "UPDATE", __FILE__, __LINE__);
            outputFile->cd();
            
            // EM: Write out the version and config files to the combined file
            std::vector<std::string> configNames = {"MaCh3_Config", "Reweight_Config", "Smearing_Config"};
            for (std::size_t i = 0; i < configNames.size(); ++i) {
                const std::string& name = configNames[i];
                TMacro* macro = prevFile->Get<TMacro>(name.c_str());
                if (macro != nullptr) {
                    macro->Write();
                    delete macro;
                }
            }
            
            // Get the source directory
            TDirectory *MaCh3EngineDir = prevFile->Get<TDirectory>("MaCh3Engine");
            TDirectory *CovarianceFolderDir = prevFile->Get<TDirectory>("CovarianceFolder");
            TDirectory *SampleFolderDir = prevFile->Get<TDirectory>("SampleFolder");
            
            CopyDir(MaCh3EngineDir);
            CopyDir(CovarianceFolderDir);
            CopyDir(SampleFolderDir);
            
            outputFile->Close();
            delete outputFile;
            
            delete prevFile;
            MACH3LOG_INFO("Done!");
        }
    }