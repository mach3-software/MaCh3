// MaCh3 includes
#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include <ROOT/RDataFrame.hxx>
_MaCh3_Safe_Include_End_ //}

/// @file PlotLLHMap.cpp
/// @brief Processing n-dimensional LLHMap outputs generating 1D and 2D profiled likelihoods as defined in config

/// @author Tomas Nosek

// *******************
// Process the parameters as given by the config file and present in the input LLHMap
std::vector<std::string> GetParams(std::vector<std::string> &PoIs, ROOT::RDataFrame Map)
// *******************
{
  std::vector<std::string> PtPs;

  if(PoIs.empty())
    MACH3LOG_INFO("No parameters requested, processing all parameters found in the LLHMap!");

  // First, filter out all non-parametric (LLH) columns of LLHMap
  for(auto p : Map.GetColumnNames())
  {
    if(std::string(p).find("_LLH") != std::string::npos)
      continue;

    PtPs.push_back(std::string(p));
  }

  // Check if the parameters user wants to plot, actually live in the LLHMap
  // Remove from the list if not
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

  // Check what LLHMap parameters the user wants to plot
  // Remove from the list if not
  for(auto pit = PtPs.begin(); pit != PtPs.end();)
  {
    if(std::find(PoIs.begin(), PoIs.end(), *pit) != PoIs.end() || PoIs.empty())
      ++pit;
    else
      pit = PtPs.erase(pit);
  }

  // Return the vector of parameters we are about to plot
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

  // Open the settings and output file
  YAML::Node Settings = M3OpenConfig(std::string(argv[1]));
  auto OutFileName = GetFromManager<std::string>(Settings["General"]["OutputFile"], "LLHMap.root");

  auto OutFile = new TFile(OutFileName.c_str(),"UPDATE");
  OutFile->cd();

  // Prepare dirs for profiled LLHs
  TDirectory* DirProfile1D = OutFile->mkdir("Profiled1D_LLH", "profile1D", true);
  TDirectory* DirProfile2D = OutFile->mkdir("Profiled2D_LLH", "profile2D", true);

  // Get the list of input files
  std::vector<std::string> inpFileList;
  for(int i = 2; i < argc; ++i)
  {
    inpFileList.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file(s): {}", inpFileList.back().c_str());
  }

  // Read the llhmap trees from the input files
  ROOT::RDataFrame LLHMap("llhmap", inpFileList);

  // Process what parameters to plot
  auto ParamsOfInterest = GetFromManager<std::vector<std::string>>(Settings["LLHScan"]["LLHParameters"],{});
  std::vector<std::string> ParamsToProfile = GetParams(ParamsOfInterest, LLHMap);

  // This now only works for uniform LLHMaps
  // TODO: Think how to allow for non-uniform LLHMaps
  MACH3LOG_INFO("!!Assuming a uniform binning of the LLHScan!!");

  MACH3LOG_INFO("... Starting generating 1D profiled log likelihoods ...");

  DirProfile1D->cd();
  for(auto p : ParamsToProfile)
  {
    // Find the binning for each histogram

    // Expected number of bins based on the config
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

    // Preparing a histogram. Again, this now works only with uniform steps in LLHMap
    double minx = LLHMap.Min(p.c_str()).GetValue();
    double maxx = LLHMap.Max(p.c_str()).GetValue();
    double binw = std::abs(maxx-minx)/double(nExpBins-1);

    minx = minx-.5*binw;
    maxx = maxx+.5*binw;

    MACH3LOG_INFO("Generating a Profile1DLLH histogram for {} of {} bins from {} to {} (bin center at {} and {})", p, nExpBins, minx, maxx, minx+.5*binw, maxx-.5*binw);

    std::string hTitle = p+" profiled LLH";
    std::string hName = p+"_LLHProf1D";
    TH1D* hprof1d = new TH1D(hName.c_str(), hTitle.c_str(), nExpBins, minx, maxx);

    // Now loop over the bins. For each bin, find the minimal Total_LLH over the rest of the parameter space.
    // TODO: Switch between different LLHs (sample, xsec, etc.)
    for(int bidx = 1; bidx < hprof1d->GetNbinsX() + 1; ++bidx)
    {
      const int count = int(double(hprof1d->GetNbinsX())/double(5));
      if (bidx % count == 0)
        MaCh3Utils::PrintProgressBar(bidx, hprof1d->GetNbinsX());

      auto b_lo = hprof1d->GetXaxis()->GetBinLowEdge(bidx);
      auto b_hi = b_lo + hprof1d->GetXaxis()->GetBinWidth(bidx);

      double llhmin = LLHMap.Filter(p+">"+std::to_string(b_lo)+"&&"+p+"<"+std::to_string(b_hi)).Min("Total_LLH").GetValue();

      if(llhmin > M3::_LARGE_LOGL_) llhmin = M3::_LARGE_LOGL_;
      hprof1d->SetBinContent(bidx, llhmin);
    }

    // Save the 1D profiled histogram
    hprof1d->Write(hName.c_str(), TObject::kOverwrite);
  }

  MACH3LOG_INFO("... Starting generating 2D profiled log likelihoods ...");

  DirProfile2D->cd();
  std::vector<std::string> Strings2D;

  // Similar as with 1D profiles, but looping over parameters twice
  for(auto p1 : ParamsToProfile)
  {
    for(auto p2 : ParamsToProfile)
    {
      // Skip whenever p1 == p2 or already profiled in reversed order (p1-p2 or p2-p1)
      if(p1 == p2)
        continue;
      if(std::find(Strings2D.begin(), Strings2D.end(), p2+"_"+p1) != Strings2D.end())
        continue;

      std::string h1Name = p1+"_LLHProf1D";
      std::string h2Name = p2+"_LLHProf1D";

      // Get the binning info from 1D histograms
      auto h1 = DirProfile1D->Get<TH1D>(h1Name.c_str());
      auto h2 = DirProfile1D->Get<TH1D>(h2Name.c_str());

      // Auxiliary
      int nx = h1->GetNbinsX();
      double minx = h1->GetBinLowEdge(1);
      double maxx = h1->GetBinLowEdge(nx)+h1->GetBinWidth(nx);

      int ny = h2->GetNbinsX();
      double miny = h2->GetBinLowEdge(1);
      double maxy = h2->GetBinLowEdge(ny)+h2->GetBinWidth(ny);

      std::string hTitle = p1+" vs. "+p2+" profiled LLH";
      std::string hName = p1+"_"+p2+"_LLHProf2D";
      TH2D* hprof2d = new TH2D(hName.c_str(), hTitle.c_str(), nx, minx, maxx, ny, miny, maxy);

      MACH3LOG_INFO("The 2D histogram has {} x-bins from {} to {} and {} y-bins from {} to {}", nx, minx, maxx, ny, miny, maxy);

      for(int bidx = 1; bidx < hprof2d->GetNbinsX() + 1; ++bidx)
      {
        for(int bidy = 1; bidy < hprof2d->GetNbinsY() + 1; ++bidy)
        {
          const int count = int(double(hprof2d->GetNbinsX()*hprof2d->GetNbinsY())/double(5));
          if ( ((bidx-1)*hprof2d->GetNbinsY() + bidy) % count == 0)
            MaCh3Utils::PrintProgressBar((bidx-1)*hprof2d->GetNbinsY() + bidy, hprof2d->GetNbinsX()*hprof2d->GetNbinsY());

          auto bx_lo = hprof2d->GetXaxis()->GetBinLowEdge(bidx);
          auto bx_hi = bx_lo + hprof2d->GetXaxis()->GetBinWidth(bidx);

          auto by_lo = hprof2d->GetYaxis()->GetBinLowEdge(bidy);
          auto by_hi = by_lo + hprof2d->GetYaxis()->GetBinWidth(bidy);

          double llhmin = LLHMap.Filter(p1+">"+std::to_string(bx_lo)+"&&"+p1+"<"+std::to_string(bx_hi)+"&&"+p2+">"+std::to_string(by_lo)+"&&"+p2+"<"+std::to_string(by_hi)).Min("Total_LLH").GetValue();

          if(llhmin > M3::_LARGE_LOGL_) llhmin = M3::_LARGE_LOGL_;
          hprof2d->SetBinContent(bidx, bidy, llhmin);
        }
      }

      hprof2d->Write(hName.c_str(), TObject::kOverwrite);

      Strings2D.push_back(p1+"_"+p2);
    }
  }

  OutFile->Close();

  return 0;
}
