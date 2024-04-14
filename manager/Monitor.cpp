#include "manager/Monitor.h"

namespace MaCh3Utils {

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

}
