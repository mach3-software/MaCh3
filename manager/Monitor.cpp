#include "manager/Monitor.h"

namespace MaCh3Utils {

// *************************
void MaCh3Welcome() {
// *************************

  // KS: Just make sure we only call it once
  static bool MaCh3WelcomeInitialised = false;

  if(MaCh3WelcomeInitialised) return;

  std::string MaCh3_VERSION = GetMaCh3Version();

  MACH3LOG_INFO("##################################");
  MACH3LOG_INFO("Welcome to:  ");
  MACH3LOG_INFO("  __  __        _____ _     ____  ");
  MACH3LOG_INFO(" |  \\/  |      / ____| |   |___ \\ ");
  MACH3LOG_INFO(" | \\  / | __ _| |    | |__   __) |");
  MACH3LOG_INFO(" | |\\/| |/ _` | |    | '_ \\ |__ < ");
  MACH3LOG_INFO(" | |  | | (_| | |____| | | |___) |");
  MACH3LOG_INFO(" |_|  |_|\\__,_|\\_____|_| |_|____/ ");
  MACH3LOG_INFO("Version: {}", MaCh3_VERSION);
  MACH3LOG_INFO("##################################");

  GetCPUInfo();

  GetGPUInfo();

  #ifdef DEBUG
  GetOSInfo();
  GetDiskUsage();
  #endif

  MaCh3WelcomeInitialised = true;
}

// ************************
// KS: Get version of MaCh3
std::string GetMaCh3Version() {
// ************************

  //KS: Find MaCh3 version based on header file. There could be better way to just include version.h but as long as we don't have to hardcode version I am content
  std::string MaCh3_VERSION = "";

  std::string file = std::string(std::getenv("MaCh3_ROOT")) + "/version.h";
  // Open the version.h file
  std::ifstream versionFile(file);

  // Check if the file is opened successfully
  if (!versionFile.is_open()) {
    MACH3LOG_ERROR("Error: Couldn't open version.h {}", file);
    throw;
  }

  std::string line;
  const std::string searchKey = "MaCh3_VERSION=";

  // Read each line from the file
  while (std::getline(versionFile, line)) {
    // Search for the line containing MaCh3_VERSION
    auto pos = line.find(searchKey);
    if (pos != std::string::npos) {
      // Extract the version string
      MaCh3_VERSION = line.substr(pos + searchKey.length());
      MaCh3_VERSION.erase(0, MaCh3_VERSION.find_first_not_of("\"")); // Remove leading quote
      MaCh3_VERSION.erase(MaCh3_VERSION.find_last_not_of("\";") + 1); // Remove trailing quote and semicolon
      break; // Stop searching once found
    }
  }
  // Close the file
  versionFile.close();

  return MaCh3_VERSION;
}


// ************************
// KS: Find out more about operational system
void GetOSInfo() {
// ************************

  MACH3LOG_INFO("Operating System Information:");

  // Distribution and version
  MACH3LOG_INFO("Distribution: {}", TerminalToString("lsb_release -d | awk -F':' '{print $2}'"));
  MACH3LOG_INFO("Kernel Version: {}", TerminalToString("uname -r"));
}


// ************************
//KS: Simple function retrieving CPU info
void GetCPUInfo(){
// ************************

  //KS: Use -m 1 to limit to only one grep because in one computing node there is a lot of CPU which are the same
  MACH3LOG_INFO("Using following CPU:");

  MACH3LOG_INFO("{}", TerminalToString("cat /proc/cpuinfo | grep -m 1 name"));
  MACH3LOG_INFO("{}", TerminalToString("cat /proc/cpuinfo | grep -m 1 MHz"));
  //KS: Below code is convoluted because I mostly work on English based Linux but sometimes on Polish based Linux, this ensures it works on both. We can add support for other languages if needed
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -i Archit"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -i 'Cache L1d'"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -i 'Cache L1i'"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -i 'Cache L2'"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -i 'Cache L3'"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -m 1 -E 'Thread.* per core:|Wątków na rdzeń:'"));
  MACH3LOG_INFO("{}", TerminalToString("lscpu | grep -m 1 -E '^CPU(:|\\(s\\)):?\\s+[0-9]+'"));

  //KS: /proc/cpuinfo and lscpu holds much more info I have limited it but one can expand it if needed
}


// ************************
//KS: Simple function retrieving GPU info
void GetGPUInfo(){
// ************************

#ifdef CUDA
  MACH3LOG_INFO("Using following GPU:");
  // Print GPU name
  MACH3LOG_INFO("GPU Name: {}", TerminalToString("nvidia-smi --query-gpu=name --format=csv,noheader"));
  // Print number of GPUs
  MACH3LOG_INFO("Number of GPUs: {}", TerminalToString("nvidia-smi --query-gpu=count --format=csv,noheader"));
  // Print total VRAM
  MACH3LOG_INFO("Total VRAM: {} MB", TerminalToString("nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits"));
  // Print Driver Version
  MACH3LOG_INFO("Driver Version: {}", TerminalToString("nvidia-smi --query-gpu=driver_version --format=csv,noheader"));
#endif
  return;
}

// ************************
// KS: Find out about Disk usage
void GetDiskUsage() {
// ************************
  MACH3LOG_INFO("Disk Usage:");

  // Get disk usage
  MACH3LOG_INFO("{}", TerminalToString("df -h --total | grep total"));
}


// ************************
// KS: Convoluted code to grab output from terminal to string
std::string TerminalToString(const char* cmd) {
// ************************

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  // Remove trailing newline characters
  result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
  return result;
}

// ************************
//KS: Simple to retrieve speed of get entry inspired by
void EstimateDataTransferRate(TChain* chain, const int entry){
// ************************

  TStopwatch timer;

  timer.Start();
  Int_t bytesProcessed{ chain->GetEntry(entry) };

  timer.Stop();

  Double_t timeInSeconds = timer.RealTime();
  Double_t dataRateMBps = (double(bytesProcessed) / (1024.0 * 1024.0)) / timeInSeconds;

  MACH3LOG_INFO("Data transfer: {} B, rate: {:.2f} MB/s", bytesProcessed, dataRateMBps);
}



// ************************
//KS: Simply print progress bar
void PrintProgressBar(const int Done, const int All){
// ************************

  double progress = double((double(Done)/double(All)));
  const int barWidth = 20;
  std::ostringstream progressBar;

  progressBar << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      progressBar << "=";
    else if (i == pos)
      progressBar << ">";
    else
      progressBar << " ";
  }

  progressBar << "] " << std::setw(3) << Done <<"/"<< All<<" ("<<static_cast<int>(progress * 100.0)<<"%)\r";
  MACH3LOG_INFO("{}", progressBar.str());
}

// ***************************************************************************
//CW: Get memory, which is probably silly
int getValue(std::string Type){ //Note: this value is in KB!
// ***************************************************************************
  std::ifstream file("/proc/self/status");
  int result = -1;
  std::string line;

  if (Type == "VmSize")
  {
    while (std::getline(file, line))
    {
      if (line.compare(0, 7, "VmSize:") == 0)
      {
        result = parseLine(line.substr(7));
        break;
      }
    }
  }
  else if (Type == "VmRSS")
  {
    while (std::getline(file, line))
    {
      if (line.compare(0, 6, "VmRSS:") == 0)
      {
        result = parseLine(line.substr(6));
        break;
      }
    }
  }
  else if (Type == "MemTotal")
  {
    std::ifstream meminfo("/proc/meminfo");
    while (std::getline(meminfo, line))
    {
      if (line.find("MemTotal:") != std::string::npos) {
        result = parseLine(line.substr(9));
        break;
      }
    }
  }
  else
  {
    std::cerr << "Not supported getValue: " << Type << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  return result;
}

// ***************************************************************************
//CW: Get memory, which is probably silly
int parseLine(const std::string& line){
// ***************************************************************************
  std::istringstream iss(line);
  int value;
  iss >> value;
  return value;
}


// ***************************************************************************
//KS: Print Yaml config using logger
void PrintConfig(const YAML::Node& node){
// ***************************************************************************
  std::stringstream ss;
  ss << node;
  std::string yamlString = ss.str();

  std::istringstream iss(yamlString);
  std::string line;
  while (std::getline(iss, line)) {
    MACH3LOG_INFO("{}", line);
  }
}



}
