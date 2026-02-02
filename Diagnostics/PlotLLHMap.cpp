// MaCh3 includes
#include "Manager/Manager.h"

// ROOT includes
#include <ROOT/RDataFrame.hxx>
#include "TCanvas.h"
#include "TChain.h"


std::vector<std::string> GetParams(std::vector<std::string> &PoIs, ROOT::RDataFrame Map)
{
  std::vector<std::string> PtPs;

  if(PoIs.empty())
    MACH3LOG_INFO("No parameters requested, processing all parameters found in the LLHMap!");

  for(auto p : Map.GetColumnNames())
  {
    if(std::string(p).find("_LLH") != std::string::npos)
      continue;

    PtPs.push_back(std::string(p));
  }

  for(auto pit = PoIs.begin(); pit != PoIs.end();)
  {
    if(std::find(PtPs.begin(), PtPs.end(), *pit) != PtPs.end())
    {
      ++pit;
    } else {
      MACH3LOG_WARN("Parameter {} not found in the LLHMap!", *pit);
      pit = PoIs.erase(pit);
    }
  }

  for(auto pit = PtPs.begin(); pit != PtPs.end();)
  {
    if(std::find(PoIs.begin(), PoIs.end(), *pit) != PoIs.end() || PoIs.empty())
      ++pit;
    else
      pit = PtPs.erase(pit);
  }

  for(auto p : PtPs)
    MACH3LOG_INFO("Param: {}",p);

  return PtPs;
}

// *******************
int main(int argc, char *argv[]) {
// *******************
  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();

  if(argc < 3)
  {
    MACH3LOG_ERROR("No arguments! Usage: {} <config.yaml> <file1.root> <file2.root> ...", argv[0]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }


  YAML::Node Settings = M3OpenConfig(std::string(argv[1]));
  auto OutFileName = GetFromManager<std::string>(Settings["General"]["OutputFile"], "LLHMap.root");

  auto OutFile = new TFile(OutFileName.c_str(),"UPDATE");
  OutFile->cd();
  TDirectory* DirProfile1D = OutFile->mkdir("Profiled1D_LLH", "profile1D", true);
  TDirectory* DirProfile2D = OutFile->mkdir("Profiled2D_LLH", "profile2D", true);

  // Get the list of input files
  std::vector<std::string> inpFileList;
  for(int i = 2; i < argc; ++i)
  {
    inpFileList.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file(s): {}", inpFileList.back().c_str());
  }

  ROOT::RDataFrame LLHMap("llhmap", inpFileList);

  auto ParamsOfInterest = GetFromManager<std::vector<std::string>>(Settings["LLHScan"]["LLHParameters"],{});
  std::vector<std::string> ParamsToProfile = GetParams(ParamsOfInterest, LLHMap);

  MACH3LOG_INFO("!!Assuming a uniform binning of the LLHScan!!");

  MACH3LOG_INFO("... Starting generating 1D profiled log likelihoods ...");

  DirProfile1D->cd();
  for(auto p : ParamsToProfile)
  {
    // Find the binning

    // Expected number of bins
    int nExpBins = GetFromManager<int>(Settings["LLHScan"]["LLHScanPoints"], 20, __FILE__, __LINE__);
    if(CheckNodeExists(Settings,"LLHScan","ScanPoints"))
      nExpBins = GetFromManager<int>(Settings["LLHScan"]["ScanPoints"][p], nExpBins, __FILE__, __LINE__);

    // Calculate the number of bins from the LLHMap
    std::set<double> cBins;
    auto checkBin = [&cBins](double var) {
      if(!(std::find(cBins.begin(), cBins.end(), var) != cBins.end())) cBins.emplace(var);
    };

    LLHMap.Foreach(checkBin, {p.c_str()});

    if(nExpBins != int(cBins.size())) {
      MACH3LOG_WARN("The config expects different number of {} bins for parameter {} than {} values included in LLHMap!", nExpBins, p, cBins.size());
      nExpBins = int(cBins.size());
    }

    MACH3LOG_INFO("There are {} values (bins) for parameter {} inside LLHMap.", nExpBins, p);


    double minx = LLHMap.Min(p.c_str()).GetValue();
    double maxx = LLHMap.Max(p.c_str()).GetValue();
    double binw = std::abs(maxx-minx)/double(nExpBins-1);

    minx = minx-.5*binw;
    maxx = maxx+.5*binw;

    MACH3LOG_INFO("Generating a Profile1DLLH histogram for {} of {} bins from {} to {} (bin center at {} and {})", p, nExpBins, minx, maxx, minx+.5*binw, maxx-.5*binw);

    std::string hTitle = p+" profiled LLH";
    std::string hName = p+"_LLHProf1D";
    TH1D* hprof1d = new TH1D(hName.c_str(), hTitle.c_str(), nExpBins, minx, maxx);

    for(int bidx = 1; bidx < hprof1d->GetNbinsX() + 1; ++bidx)
    {
      const int count = int(double(hprof1d->GetNbinsX())/double(5));
      if (bidx % count == 0)
        MaCh3Utils::PrintProgressBar(bidx, hprof1d->GetNbinsX());

      auto b_lo = hprof1d->GetXaxis()->GetBinLowEdge(bidx);
      auto b_hi = b_lo + hprof1d->GetXaxis()->GetBinWidth(bidx);

      //MACH3LOG_INFO("Reading bin from {} to {}", b_lo, b_hi);

      double llhmin = LLHMap.Filter(p+">"+std::to_string(b_lo)+"&& "+p+"<"+std::to_string(b_hi)).Min("Total_LLH").GetValue();

      //MACH3LOG_INFO("llhmin is {}", llhmin);

      if(llhmin > M3::_LARGE_LOGL_) llhmin = M3::_LARGE_LOGL_;
      hprof1d->SetBinContent(bidx, llhmin);

      //MACH3LOG_INFO("Bin number {} has {}", bidx, hprof1d->GetBinContent(bidx));
    }

    hprof1d->Write(hName.c_str(), TObject::kOverwrite);
  }

  MACH3LOG_INFO("... Starting generating 2D profiled log likelihoods ...");
  DirProfile2D->cd();

  OutFile->Close();

  return .0;
}
