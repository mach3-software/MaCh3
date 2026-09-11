/// @file CombineModule.hpp
/// @brief Allows you to combine multiple chains or ROOT files into a single output.

#pragma once
#include "CLI/API/plugin.hpp"
// C++ includes
#include <unistd.h>

// MaCh3 includes
#include "Manager/Manager.h"
#include "Samples/SampleStructs.h"
#include "Samples/HistogramUtils.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TList.h"
#include "TFile.h"
#include "TMacro.h"
#include "TTree.h"
#include "TMD5.h"
#include "TFileMerger.h"
#include "TKey.h"
#include "TROOT.h"
_MaCh3_Safe_Include_End_ //}

namespace M3{

  /// @class CombineModule
  /// @brief Allows you to combine multiple chains or ROOT files into a single output.
  class CombineModule: public IModuleBase{

    public:
      /// @brief Destructor
      virtual ~CombineModule();

      /// @brief Get the argument parser for this module
      /// @return Pointer to the configured MaCh3ArgumentParser
      MaCh3ArgumentParser* get_parser() override;

      /// @brief Execute the application
      /// @return Exit code (0 on success)
      int Run() override;

      /// @brief KS: This allow us to skip output name etc in config. We expect Output name will be different but this doesn't invalidate chain merging
      bool ShouldSkipLine(const std::string& line, const std::vector<std::string>& SkipVector);

      /// @brief make sure two configs are identical but skip specified fields. For example when comparing two chains nsteps or output name might be different and this is still fine to merge
      /// @param File1 Config from chain1
      /// @param File2 Config from chain2
      /// @param SkipVector Fields in yaml file to skip
      bool CompareTwoConfigs(const std::string& File1, const std::string& File2, const std::vector<std::string>& SkipVector);

      /// @brief EM: Will compare the version header contained in the two provided files and shout if they don't match
      bool checkSoftwareVersions(TFile *file, TFile *prevFile, const std::string& ConfigName, const std::vector<std::string>& SkipVector = {});

      /// @brief When we merge two chains they have TDirectory ROOT didn't provide method for this so here we have this bad boy
      void CopyDir(TDirectory *source);

      /// @brief Compare two histograms if they are identical
      /// @todo add checks for stuff like bin content etc
      bool CompareHistograms(const TH1* h1, const TH1* h2, const std::string& histName, const std::string& folderName);

      /// @brief Loop through TH1 and TMacro objects in FolderName in 'file' and compare with those in 'prevFile'
      bool CheckFolder(TFile* file, TFile* prevFile, const std::string& FolderName, const std::vector<std::string>& SkipVector = {});

      /// @brief custom function for merging TTree, should be similar to what HADD is using
      /// @warning KS: for some reason if "fast" is enable then I cannot open in ROOT5, no one should use R5 at this point..
      void FastMergeTTrees(const std::vector<std::string>& files, const std::string& outFile, const std::string& TTreeName);

      void CombineChain();
          
    private:
      std::string OutFileName = "";
      int targetCompression = 1;
      std::vector<std::string> inpFileList;
      bool forceOverwrite = false;
      bool forceMerge = false;

  };
}
