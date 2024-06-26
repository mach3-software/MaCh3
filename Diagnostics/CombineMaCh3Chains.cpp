#include "TList.h"
#include "TFile.h"
#include "TMacro.h"
#include "TTree.h"
#include "TMD5.h"
#include "TFileMerger.h"
#include "TKey.h"
#include "TROOT.h"

#include "manager/manager.h"

std::string OutFileName = "";
int targetCompression = 1;
std::vector<std::string> inpFileList;
bool forceOverwrite = false;

/// EM: Will compare the version header contained in the two provided files and shout if they don't match
bool checkSoftwareVersions(TFile *file, TFile *prevFile, std::string ConfigName)
{
  bool weirdFile = false;

  TMacro *versionHeader = (TMacro*)file->Get(ConfigName.c_str());
  TMacro *prevVersionHeader = (TMacro*)prevFile->Get(ConfigName.c_str());

  // EM: compare the digest of the version header file in this file, with the previous one
  if((versionHeader == NULL) && (prevVersionHeader == NULL)){
    MACH3LOG_WARN("files don't contain an embedded version header, indicating they were made before this feature was added");
    MACH3LOG_WARN("  I can still combine them but I can't guarantee they were made with matching software versions");
    MACH3LOG_WARN("  This is ok but you'll just have to be extra careful and make sure you check this yourself!");
  }
  else if((versionHeader != NULL) && (prevVersionHeader == NULL)){
    MACH3LOG_ERROR("looks like file {} has a version header file but previous ones do not", file->GetName());
    MACH3LOG_ERROR("This is odd and suggests this file was made with an MaCh3 version");
    MACH3LOG_ERROR("from after this feature was added, whereas other files were made with an older MaCh3");
    weirdFile = true;
  }
  else if((versionHeader == NULL) && (prevVersionHeader != NULL)){
    MACH3LOG_ERROR("looks like file {} doesn't have a version header file but previous ones do", file->GetName());
    MACH3LOG_ERROR("This is odd and suggests this file was made with an MaCh3 version");
    MACH3LOG_ERROR("from before this feature was added, whereas other files were made with newer MaCh3");
    weirdFile = true;
  }
  else{
    MACH3LOG_DEBUG("  Prev header digest: {} :: current: {}", (prevVersionHeader->Checksum())->AsString(), (versionHeader->Checksum())->AsString());

    if(std::strcmp((versionHeader->Checksum())->AsString(), (prevVersionHeader->Checksum())->AsString()) != 0){
      MACH3LOG_ERROR("Looks like the version header for file {} is different to the previous ones", file->GetName());
      MACH3LOG_ERROR("This strongly suggests that this file was made with different software versions than the previous ones");
      weirdFile = true;
    }
  }
  return weirdFile;
}

void CopyDir(TDirectory *source) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  source->ls();
  TDirectory *savdir = gDirectory;
  TDirectory *adir = (TDirectory*)savdir->Get(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
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
      TTree *T = (TTree*)source->Get(key->GetName());
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

void CombineChain()
{
  TFileMerger *fileMerger = new TFileMerger();

  // EM: If we ever add new trees to the chain files they will need to be added here too
  fileMerger->AddObjectNames("posteriors");
  fileMerger->AddObjectNames("Settings");

  MACH3LOG_INFO("These objects will be merged: {}", fileMerger->GetObjectNames());

  std::string outFileOption;
  if(forceOverwrite) outFileOption = "RECREATE";
  else outFileOption = "CREATE";

  // EM: Attempt to open the output file
  bool openedFile = fileMerger->OutputFile(OutFileName.c_str(), outFileOption.c_str(), targetCompression);
  if (!openedFile){
      MACH3LOG_ERROR("Failed to create output file.");
      throw MaCh3Exception(__FILE__ , __LINE__ );
  }

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


    if(checkSoftwareVersions(file, prevFile, "MaCh3Engine/version_header"))
        weirdFile = true;

    if(checkSoftwareVersions(file, prevFile, "MaCh3_Config"))
      weirdFile = true;

    if(weirdFile){
      MACH3LOG_ERROR("");
      MACH3LOG_ERROR("=====================================================================================");
      MACH3LOG_ERROR("This is not a great idea and could lead to weird outputs and cause some big headaches");
      MACH3LOG_ERROR("further down the road. But if you reeeeally wanna do it and you know what you're");
      MACH3LOG_ERROR("doing you can come here and remove the 'throw'");
      MACH3LOG_ERROR("{}:{}", __FILE__, __LINE__ + 2);
      MACH3LOG_ERROR("=====================================================================================");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    // EM: file seems good, we'll add the trees to the lists
    fileMerger->AddFile(file);

    // EM: set these for the next iteration
    prevFile = file;
  }

  TFile *outputFile = fileMerger->GetOutputFile();
  outputFile->cd();

  // EM: write out the version and config files to the combined file
  TMacro *MaCh3_Config = (TMacro*)prevFile->Get("MaCh3_Config");

  if(MaCh3_Config != NULL) MaCh3_Config->Write();
  delete MaCh3_Config;

  // EM: now let's combine all the trees and write to the output file
  bool mergeSuccess = fileMerger->PartialMerge(TFileMerger::kRegular | TFileMerger::kAll | TFileMerger::kOnlyListed);
  if(mergeSuccess){
    MACH3LOG_INFO("Files merged successfully");
  }
  else{
    MACH3LOG_ERROR("Failed to merge files");
  }
  delete fileMerger;

  //KS: Sadly we need to open file to save TDirectories to not have weird copy of several obejcts there...
  outputFile = new TFile(OutFileName.c_str(), "UPDATE");

  // Get the source directory
  TDirectory *MaCh3EngineDir = (TDirectory*)prevFile->Get("MaCh3Engine");
  TDirectory *CovarianceFolderDir = (TDirectory*)prevFile->Get("CovarianceFolder");

  outputFile->cd();
  CopyDir(MaCh3EngineDir);
  CopyDir(CovarianceFolderDir);

  delete prevFile;

  MACH3LOG_INFO("Done!");
}
    
void usage(){
  MACH3LOG_INFO("Combine MaCh3 Chains files, very similar to hadd, but will compare embedded version info in the files to avoid accidentally combining files made with different software versions. Also avoids having a hige dump of separate version files in the output that happens with hadd.");
  MACH3LOG_INFO("Cmd line syntax should be:");
  MACH3LOG_INFO("combineND280Splines [-h] [-c [0-9]] [-f] [-o <output file>] input1.root [input2.root, input3.root ...]");
  MACH3LOG_INFO("inputX.root    : names of individual spline files to combine, can specify any number, need at least one");
  MACH3LOG_INFO("output file    : name of combined spline file. optional: if not specified, the app will just use the first input file as the output, the same as hadd'");
  MACH3LOG_INFO("-c             : target compression level for the combined file, default is 1, in line with hadd");
  MACH3LOG_INFO("-f             : force overwrite the output file if it exists already");
  MACH3LOG_INFO("-h             : print this message and exit");
}

void ParseArg(int argc, char *argv[]){
  if(argc < 2){
    MACH3LOG_ERROR("Too few arguments!!");
    MACH3LOG_ERROR("USAGE:");
    usage();
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  int c;
  for(;;) {
    c = getopt(argc, argv, "o:c:hf");
    if (c == -1){ // loop over the remaining arguments
      while (optind < argc){
        // any non option input is assumed to be a root file
        std::string fName = std::string(argv[optind]);
        MACH3LOG_DEBUG("adding {} to file list", fName.c_str());
        inpFileList.push_back(fName);
        optind ++;
      }
      break;
    }
    else{
      switch (c) {
        case 'o': {
          OutFileName = optarg;
          break;
        }
        case 'f': {
          forceOverwrite = true;
          break;
        }
        case 'c': {
          targetCompression = atoi(optarg);
          break;
        }
        case 'h': {
          usage();
          exit(0);
        }
        default: {
          MACH3LOG_ERROR("Un recognised option");
          usage();
          exit(1);
        }
      }
    }
  }

  if(OutFileName == ""){
    MACH3LOG_INFO("Using first file in list as output: ", inpFileList[0].c_str());
    OutFileName = inpFileList[0];
    inpFileList.erase(inpFileList.begin());
  }

  if(forceOverwrite){
    MACH3LOG_INFO("Will overwrite {} if it exists already", OutFileName.c_str());
  }

  MACH3LOG_INFO("Combining a total of {} files into {}", inpFileList.size(), OutFileName.c_str());
}

int main(int argc, char *argv[])
{
  SetMaCh3LoggerFormat();
  ParseArg(argc, argv);
  CombineChain();
  return 1;
}
