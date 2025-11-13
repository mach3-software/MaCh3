#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT include
#include "TFile.h"
_MaCh3_Safe_Include_End_ //}

// *************************
manager::manager(std::string const &filename)
: config(M3OpenConfig(filename)) {
// *************************
  FileName = filename;

  Initialise();
}

// *************************
manager::manager(const YAML::Node ConfigNode) {
// *************************
  config = ConfigNode;
  FileName = "unknown";

  Initialise();
}

// *************************
void manager::Initialise() {
// *************************


  #ifdef MPIENABLED
  InitialiseMPI();
  #endif

    MACH3LOG_INFO("Setting config to be: {}", FileName);

  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();

  MACH3LOG_INFO("Config is now: ");
  MaCh3Utils::PrintConfig(config);
}

#ifdef MPIENABLED
void manager::InitialiseMPI(){
  int ini = 0;
  MPI_Initialized(&ini);
  if(!ini){
    MACH3LOG_WARN("MPI not initialised, not running as an MPI process");
  }


  // Print out MPI information
  int total_processes = 0;
  SetMaCh3MPILogging();
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
  MACH3LOG_INFO("Starting MaCh3 with MPI support, total processes: {}", total_processes);

  // Register a safe finalizer to be run at program exit. This avoids
  // calling MPI_Finalize() from the manager destructor where static
  // destruction order can cause "MPI_FINALIZE CANNOT..." errors.
  std::atexit([]()
  {
    int initialised = 0;
    int finalised = 0;
    MPI_Initialized(&initialised);
    MPI_Finalized(&finalised);
    if (initialised && !finalised) {
      MPI_Finalize();
    } });
    

  }
  #endif

// *************************
// Empty destructor, for now...
manager::~manager() {
// *************************
}

// *************************
// Save all the settings of the class to an output file
// Reflection in C++ is a bit of a pain :(
// Outputfile is the TFile pointer we write to
void manager::SaveSettings(TFile* const OutputFile) const {
// *************************
  std::string OutputFilename = std::string(OutputFile->GetName());
  OutputFile->cd();

  // EM: embed the config used for this app
  TMacro ConfigSave = YAMLtoTMacro(config, "MaCh3_Config");
  ConfigSave.Write();

  // The Branch!
  TTree *SaveBranch = new TTree("Settings", "Settings");

  // Fill the doubles
  SaveBranch->Branch("Output", &OutputFilename);

  // Get settings defined by pre-processor directives, e.g. CPU MP and GPU
  #ifdef MULTITHREAD
  bool cpu_mp_on = true;
  int n_cpus = omp_get_max_threads();
  #else
  bool cpu_mp_on = false;
  int n_cpus = 1;
  #endif

  #ifdef MaCh3_CUDA
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
void manager::Print() const {
// *************************
  MACH3LOG_INFO("---------------------------------");
  MaCh3Utils::PrintConfig(config);
  MACH3LOG_INFO("---------------------------------");
}

// *************************
int manager::GetMCStatLLH() const {
// *************************
  int mc_stat_llh = kNTestStatistics;
  if (config["LikelihoodOptions"])
  {
    auto likelihood = GetFromManager<std::string>(config["LikelihoodOptions"]["TestStatistic"], "Barlow-Beeston", __FILE__ , __LINE__);

    for(int i = 0; i < kNTestStatistics; i++) {
      if(likelihood == TestStatistic_ToString(TestStatistic(i))) {
        mc_stat_llh = TestStatistic(i);
        break;
      }
    }
    if(mc_stat_llh == kNTestStatistics)
    {
      MACH3LOG_ERROR("Wrong form of test-statistic specified!");
      MACH3LOG_ERROR("You gave {} and I only support:", likelihood);
      for(int i = 0; i < kNTestStatistics; i++)
      {
        MACH3LOG_ERROR("{}", TestStatistic_ToString(TestStatistic(i)));
      }
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  } else {
    MACH3LOG_WARN("Didn't find a TestStatistic specified");
    MACH3LOG_WARN("Defaulting to using a {} likelihood", TestStatistic_ToString(kPoisson));
    mc_stat_llh = kPoisson;
  }
  return mc_stat_llh;
}


