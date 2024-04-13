#include "manager/Monitor.h"



namespace MaCh3Utils {


// ************************
//KS: Simple function retrieving CPU info
void GetCPUInfo(){
// ************************
  //KS: Use -m 1 to limit to only one grep becasue in one computing node there is a lot of CPU which are the same
  std::cout << "================" << std::endl;
  std::cout << "Using following CPU: " << std::endl;
  system("cat /proc/cpuinfo | grep -m 1 name");
  system("cat /proc/cpuinfo | grep -m 1 MHz");
  //KS: Below code is convoluted because I moslty work on Enlgish based linux but sometiems on Polish based Linux, this ensures it works on both. We can add support for othre langauges if needed
  std::system("lscpu | grep -i Archit");
  std::system("lscpu | grep -i cache");
  std::system("lscpu | grep -m 1 -E 'Thread.* per core:|Wątków na rdzeń:'");
  std::system("lscpu | grep -m 1 -E '^CPU(:|\\(s\\)):?\\s+[0-9]+'");
  std::cout << "================" << std::endl;

  //KS: /proc/cpuinfo and lscpu holds much more info I have limited it but one can expand it if needed
}


// ************************
//KS: Simple to rtrieve speed of get entry inspired by
void EstimateDataTransferRate(TChain* chain, const int entry){
// ************************

  TStopwatch timer;

  timer.Start();
  Int_t bytesProcessed{ chain->GetEntry(entry) };

  timer.Stop();

  Double_t timeInSeconds = timer.RealTime();
  Double_t dataRateMBps = (double(bytesProcessed) / (1024.0 * 1024.0)) / timeInSeconds;

  std::cout << "Data transfer: " << bytesProcessed << " B, rate: " << dataRateMBps << " MB/s"<<std::endl;
}



// ************************
//KS: Simply print progress bar
void PrintProgressBar(const double progress){
// ************************
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

  progressBar << "] " << std::setw(3) << static_cast<int>(progress * 100.0) << "%\r";

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
