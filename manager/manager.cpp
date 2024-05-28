#include "TFile.h"

#include "manager/manager.h"

// *************************
manager::manager(std::string const &filename)
    : config(YAML::LoadFile(filename)) {
// *************************

  FileName = filename;
  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();

  MACH3LOG_INFO("Setting config to be: {}", filename);

  MACH3LOG_INFO("Config is now: ");
  MaCh3Utils::PrintConfig(config);

  if (config["LikelihoodOptions"])
  {
    std::string likelihood = GetFromManager<std::string>(config["LikelihoodOptions"]["TestStatistic"], "Barlow-Beeston");
    if (likelihood == "Barlow-Beeston")                 mc_stat_llh = TestStatistic(kBarlowBeeston);
    else if (likelihood == "IceCube")                   mc_stat_llh = TestStatistic(kIceCube);
    else if (likelihood == "Poisson")                   mc_stat_llh = TestStatistic(kPoisson);
    else if (likelihood == "Pearson")                   mc_stat_llh = TestStatistic(kPearson);
    else if (likelihood == "Dembinski-Abdelmotteleb")   mc_stat_llh = TestStatistic(kDembinskiAbdelmottele);
    else {
      MACH3LOG_ERROR("Wrong form of test-statistic specified!");
      MACH3LOG_ERROR("You gave {} and I only support:", likelihood);
      for(int i = 0; i < kNTestStatistics; i++)
      {
        MACH3LOG_ERROR("{}", TestStatistic_ToString(TestStatistic(i)));
      }
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  } else {
    mc_stat_llh = kPoisson;
  }

  Modes = nullptr;
  std::string ModeInput = GetFromManager<std::string>(config["General"]["MaCh3Modes"], "null");
  if(ModeInput != "null") Modes = new MaCh3Modes(ModeInput);
}

// *************************
// Empty destructor, for now...
manager::~manager() {
// *************************

  if(!Modes) delete Modes;

}

// *************************
// Save all the settings of the class to an output file
// Reflection in C++ is a bit of a pain :(
// Outputfile is the TFile pointer we write to
void manager::SaveSettings(TFile* const OutputFile) {
// *************************

  std::string OutputFilename = std::string(OutputFile->GetName());
  OutputFile->cd();

  // EM: embed the config used for this app
  TMacro MaCh3Config("MaCh3_Config", "MaCh3_Config");
  MaCh3Config.ReadFile(FileName.c_str());
  MaCh3Config.Write();

  // The Branch!
  TTree *SaveBranch = new TTree("Settings", "Settings");

  // Fill the doubles
  SaveBranch->Branch("Output", &OutputFilename);
  SaveBranch->Branch("TestStatistic", &mc_stat_llh);

  // Get settings defined by pre-processor directives, e.g. CPU MP and GPU
  #ifdef MULTITHREAD
  bool cpu_mp_on = true;
  int n_cpus = omp_get_max_threads();
  #else
  bool cpu_mp_on = false;
  int n_cpus = 1;
  #endif

  #ifdef CUDA
  bool gpu_on = true;
  #else
  bool gpu_on = false;
  #endif

  SaveBranch->Branch("GPU",   &gpu_on);
  SaveBranch->Branch("CPUMP", &cpu_mp_on);
  SaveBranch->Branch("nCPUs", &n_cpus);

  SaveBranch->Fill();
  SaveBranch->Write();

  delete SaveBranch;
}

// *************************
void manager::Print() {
// *************************

  MACH3LOG_INFO("---------------------------------");
  std::cout << config << "\n";
  MACH3LOG_INFO("---------------------------------");
}
