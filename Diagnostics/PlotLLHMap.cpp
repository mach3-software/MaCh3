// MaCh3 includes
#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include <ROOT/RDataFrame.hxx>
_MaCh3_Safe_Include_End_ //}

/// @file PlotLLHMap.cpp
/// @brief Processing n-dimensional LLHMap outputs generating 1D and 2D profiled likelihoods as defined in config
/// @ingroup MaCh3DiagnosticProcessing
///
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
    if(p.find("_LLH") != std::string::npos)
      continue;

    PtPs.push_back(p);
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
// Extract binning
std::pair<bool,std::vector<double>> ExtractBinning(std::string param, YAML::Node Settings, ROOT::RDataFrame Map)
// *******************
{
  // Expected number of bins based on the config
  unsigned int nExpBins = GetFromManager<unsigned int>(Settings["LLHScan"]["LLHScanPoints"], 20, __FILE__, __LINE__);
  if(CheckNodeExists(Settings,"LLHScan","ScanPoints"))
    nExpBins = GetFromManager<unsigned int>(Settings["LLHScan"]["ScanPoints"][param], nExpBins, __FILE__, __LINE__);


  // Preparing a histogram. Again, this now works only with uniform steps in LLHMap
  auto values = Map.Take<double>(param.c_str());

  // remove duplicates
  std::set<double> unique(values->begin(), values->end());

  MACH3LOG_INFO("There are {} values (bins) for parameter {} inside LLHMap.", unique.size(), param);

  if(nExpBins != static_cast<unsigned int>(unique.size())) {
    MACH3LOG_WARN("The config expects different number of {} bins for parameter {} than {} values included in LLHMap!", nExpBins, param, unique.size());
    nExpBins = static_cast<unsigned int>(unique.size());
  }

  // Extract minimum and mximum
  double minx = *unique.begin();
  double maxx = *unique.rbegin();

  // What is the minimal difference between two unique values to construct a bin width
  double minDiff = maxx - minx;

  for (auto it = std::next(unique.begin()); it != unique.end(); ++it) {
    auto prev = std::prev(it);
    minDiff = std::min(minDiff, *it - *prev);
  }

  if (unique.size() > 1) {
    unsigned int bins = static_cast<unsigned int>( 1 + std::abs(maxx-minx) / minDiff );
    if (nExpBins != bins) {
      MACH3LOG_WARN("The scan is not uniform! Instead, we will use bins based on unique scan values and minimal {:.3e} distance between them!", minDiff);
      nExpBins = 0;
    }
  }

  double half_w = std::max(1e-12*minx, .5*minDiff);
  minx = minx-half_w;
  maxx = maxx+half_w;

  std::vector<double> binning = {minx, maxx};


  if (nExpBins == 0) {
    for (auto it = std::next(unique.begin()); std::next(it) != unique.end(); ++it)
    {
      // first, insert next edge
      double next_edge = *it + half_w;
      binning.insert(binning.end()-1, next_edge);

      // now check, if it aligns with the next unique value
      auto next = std::next(it);
      if ( std::abs((next_edge - *next + half_w) / (next_edge + *next - half_w )) < 1e-12 )
        continue;
      else
        binning.insert(binning.end()-1, *next-half_w);
    }
  } else {
    for (auto it = unique.begin(); std::next(it) != unique.end(); ++it)
      binning.insert(binning.end()-1, *it + half_w);
  }

  return std::make_pair(nExpBins>0, binning);
}

// *******************
int main(int argc, char *argv[]) {
// *******************
  SetMaCh3LoggerFormat();
  M3::Utils::MaCh3Welcome();

  if(argc < 3)
  {
    MACH3LOG_ERROR("No arguments! Usage: {} <config.yaml> <file1.root> <file2.root> ...", argv[0]);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  // Open the settings and output file
  YAML::Node Settings = M3OpenConfig(std::string(argv[1]));
  auto OutFileName = GetFromManager<std::string>(Settings["General"]["OutputFile"], "LLHMap.root");
  auto Plot2D = GetFromManager<bool>(Settings["LLHScan"]["Plot2D"], false, __FILE__, __LINE__);

  auto OutFile = new TFile(OutFileName.c_str(),"UPDATE");
  OutFile->cd();

  // Prepare dirs for profiled LLHs
  TDirectory* DirProfile1D = OutFile->mkdir("Profiled1D_LLH", "profile1D", true);
  TDirectory* DirProfile2D;
  if (Plot2D) DirProfile2D = OutFile->mkdir("Profiled2D_LLH", "profile2D", true);

  TDirectory* DirMarginal1D = OutFile->mkdir("Marginal1D_L", "marginal1D", true);
  TDirectory* DirMarginal2D;
  if (Plot2D) DirMarginal2D = OutFile->mkdir("Marginal2D_L", "marginal2D", true);


  // Get the list of input files
  std::vector<std::string> inpFileList;
  for(int i = 2; i < argc; ++i)
  {
    inpFileList.push_back(std::string(argv[i]));
    MACH3LOG_INFO("Adding file(s): {}", inpFileList.back().c_str());
  }

  // Read the llhmap trees from the input files
  auto Map = ROOT::RDataFrame("llhmap", inpFileList);

  // Now prepare the L column
  // TODO: Let
  auto LLHMap = Map.Define("L", "exp(-0.5*Total_LLH)");

  // Process what parameters to plot
  auto ParamsOfInterest = GetFromManager<std::vector<std::string>>(Settings["LLHScan"]["LLHParameters"],{});
  std::vector<std::string> ParamsToProfile = GetParams(ParamsOfInterest, Map);

  MACH3LOG_INFO("... Starting generating 1D profiled and marginalized likelihoods ...");
  MACH3LOG_WARN("!!! LLHMap numerical marginalization assumes uncorrelated priors !!!");

  std::map<std::string, TH1D*> hProfiles1d;
  std::map<std::string, TH1D*> hMarginals1d;
  for(auto p = ParamsToProfile.begin(); p != ParamsToProfile.end(); ++p)
  {
    // Find the binning for the histograms.
    // ExtractBinning gives a pair of bool, bin_edges.
    // The bool denotes whether binning is uniform or not.
    std::pair<bool, std::vector<double>> binning = ExtractBinning(*p, Settings, Map);

    std::string hProfTitle = *p+" profiled -2LogL";
    std::string hProfName = *p+"_LLHProf1D";
    std::string hMargTitle = *p+" marginalized L";
    std::string hMargName = *p+"_LMarg1D";

    TH1D* hprof1d = new TH1D(hProfName.c_str(), hProfTitle.c_str(), int(binning.second.size()-1), binning.second.data());
    TH1D* hmarg1d = new TH1D(hMargName.c_str(), hMargTitle.c_str(), int(binning.second.size()-1), binning.second.data());
    if (binning.first)
    {
      MACH3LOG_INFO("Initializing 1D profiled -2LogL and marginalized L histograms for {} of {} bins from {:.3e} to {:.3e} (bin center at {:.3e} and {:.3e})", *p, hprof1d->GetNbinsX(), hprof1d->GetXaxis()->GetXmin(), hprof1d->GetXaxis()->GetXmax(), hprof1d->GetBinCenter(1), hprof1d->GetBinCenter(hprof1d->GetNbinsX()));
    } else {
      std::string binnies;
      // TN: Waiting for C++ 20 std::format() function
      for (auto binnie : binning.second)
      {
        std::ostringstream os;
        os << std::scientific << std::setprecision(4) << binnie;
        binnies += " "+os.str();

      }

      MACH3LOG_INFO("Initializing 1D profiled -2LogL and marginalized L histograms of {} bins with edges{}.", binning.second.size(), binnies);
    }

    hProfiles1d[*p]  = hprof1d;
    hMarginals1d[*p] = hmarg1d;
  }

  for(auto p = ParamsToProfile.begin(); p != ParamsToProfile.end(); ++p)
  {

    // Now loop over the bins. For each bin, find the minimal Total_LLH over the rest of the parameter space.
    // TODO: Switch between different LLHs (sample, xsec, etc.)
    const int count = hProfiles1d[*p]->GetNbinsX() > 5 ? int(double(hProfiles1d[*p]->GetNbinsX())/double(5)) : 1;

    // Profile over the rest of the parameters
    MACH3LOG_INFO("Profiling 1D -2LogL and numerically profiling 1D L for parameter {}!", *p);
    for(int bidx = 1; bidx < hProfiles1d[*p]->GetNbinsX() + 1; ++bidx)
    {
      if (bidx % count == 0)
        M3::Utils::PrintProgressBar(bidx, hProfiles1d[*p]->GetNbinsX());

      auto b_lo = hProfiles1d[*p]->GetXaxis()->GetBinLowEdge(bidx);
      auto b_hi = b_lo + hProfiles1d[*p]->GetXaxis()->GetBinWidth(bidx);

      double llhmin = LLHMap.Filter(*p+">"+std::to_string(b_lo)+"&&"+*p+"<"+std::to_string(b_hi)).Min("Total_LLH").GetValue();

      // Rather store a non-sensensical value if out of bounds
      if(llhmin >= M3::_LARGE_LOGL_) llhmin = -12345;
      hProfiles1d[*p]->SetBinContent(bidx, llhmin);

      auto L = LLHMap.Filter(*p+">"+std::to_string(b_lo)+"&&"+*p+"<"+std::to_string(b_hi)).Sum("L");

      hMarginals1d[*p]->SetBinContent(bidx, *L);
    }

    // Save the 1D histograms
    DirProfile1D->cd();
    hProfiles1d[*p]->Write(hProfiles1d[*p]->GetName(), TObject::kOverwrite);

    DirMarginal1D->cd();
    hMarginals1d[*p]->Scale(1./hMarginals1d[*p]->Integral());
    hMarginals1d[*p]->Write(hMarginals1d[*p]->GetName(), TObject::kOverwrite);
  }


  if (Plot2D)
  {
    MACH3LOG_INFO("... Starting generating 2D profiled and marginalized likelihoods ...");
    MACH3LOG_WARN("!!! THIS MIGHT TAKE QUITE A WHILE !!!");
    MACH3LOG_WARN("!!! LLHMap numerical marginalization assumes uncorrelated priors !!!");

    std::vector<std::string> Keys2D;
    std::vector<std::string> ParamsFiltered;
    std::map<std::string, TH2D*> hProfiles2d;
    std::map<std::string, TH2D*> hMarginals2d;

    for(auto p : ParamsToProfile)
    {
      auto h = DirProfile1D->Get<TH1D>((p+"_LLHProf1D").c_str());
      if(h->GetNbinsX()<2)
      {
        MACH3LOG_WARN("There is less than 2 bins for {}, 2D is equivalent to 1D! Removing from 2D plots ...", p);
      } else {
        ParamsFiltered.push_back(p);
      }
    }

    // Similar as with 1D profiles, but looping over parameters twice
    for(auto p1 = ParamsFiltered.begin(); p1 != ParamsFiltered.end(); ++p1)
    {
      for(auto p2 = std::next(p1); p2 != ParamsFiltered.end(); ++p2)
      {
        // Skip whenever p1 == p2 or already profiled in reversed order (p1-p2 or p2-p1)
        // if(p1 == p2)
        //   continue;
        // if(std::find(Strings2D.begin(), Strings2D.end(), p2+"_"+p1) != Strings2D.end())
        //   continue;

        // std::string hProf1Name = p1+"_LLHProf1D";
        // std::string hProf2Name = p2+"_LLHProf1D";

        // Get the binning info from 1D histograms

        auto h1 = DirProfile1D->Get<TH1D>((*p1+"_LLHProf1D").c_str());
        auto h2 = DirProfile1D->Get<TH1D>((*p2+"_LLHProf1D").c_str());


        std::string key = *p1+"_"+*p2;
        Keys2D.push_back(key);

        MACH3LOG_INFO("Initializing 2D profiled -2LogL and marginalized L histograms for {} vs. {} based on previously generated 1D histograms.", *p1, *p2);

        std::string hProfTitle = *p1+" vs. "+*p2+" profiled -2LogL";
        std::string hProfName = key+"_LLHProf2D";
        TH2D* hprof2d = new TH2D(
            hProfName.c_str(), hProfTitle.c_str(),
            h1->GetXaxis()->GetNbins(), h1->GetXaxis()->GetXbins()->GetArray(),
            h2->GetXaxis()->GetNbins(), h2->GetXaxis()->GetXbins()->GetArray()
        );

        hProfiles2d[key] = hprof2d;

        std::string hMargTitle = *p1+" vs. "+*p2+" marginalized L";
        std::string hMargName = key+"_LMarg1D";
        TH2D* hmarg2d = new TH2D(
            hMargName.c_str(), hMargTitle.c_str(),
            h1->GetXaxis()->GetNbins(), h1->GetXaxis()->GetXbins()->GetArray(),
            h2->GetXaxis()->GetNbins(), h2->GetXaxis()->GetXbins()->GetArray()
        );

        hMarginals2d[key] = hmarg2d;
      }
    }

    for(auto p1 = ParamsFiltered.begin(); p1 != ParamsFiltered.end(); ++p1)
    {
      for(auto p2 = std::next(p1); p2 != ParamsFiltered.end(); ++p2)
      {
        std::string key = *p1+"_"+*p2;

        MACH3LOG_INFO("Numerically profiling 2D -2LogL and marginalizing 2D L for parameters {} and {}!",*p1, *p2);
        const Long64_t nBinsX = static_cast<Long64_t>(hProfiles2d[key]->GetNbinsX());
        const Long64_t nBinsY = static_cast<Long64_t>(hProfiles2d[key]->GetNbinsY());

        const Long64_t TotalBins = nBinsX * nBinsY;
        const int count = TotalBins > 5 ? int(double(TotalBins)/double(5)) : 1;

        for(int bidx = 1; bidx < nBinsX + 1; ++bidx)
        {
          for(int bidy = 1; bidy < nBinsY + 1; ++bidy)
          {
            if ( ((bidx-1)*hProfiles2d[key]->GetNbinsY() + bidy) % count == 0)
               M3::Utils::PrintProgressBar((bidx-1)*hProfiles2d[key]->GetNbinsY() + bidy, TotalBins);

            auto bx_lo = hProfiles2d[key]->GetXaxis()->GetBinLowEdge(bidx);
            auto bx_hi = bx_lo + hProfiles2d[key]->GetXaxis()->GetBinWidth(bidx);

            auto by_lo = hProfiles2d[key]->GetYaxis()->GetBinLowEdge(bidy);
            auto by_hi = by_lo + hProfiles2d[key]->GetYaxis()->GetBinWidth(bidy);

            double llhmin = LLHMap.Filter(*p1+">"+std::to_string(bx_lo)+"&&"+*p1+"<"+std::to_string(bx_hi)+"&&"+*p2+">"+std::to_string(by_lo)+"&&"+*p2+"<"+std::to_string(by_hi)).Min("Total_LLH").GetValue();

            // Rather store a non-sensensical value if out of bounds
            // TODO: Think how to do this smarter and faster
            if(llhmin >= M3::_LARGE_LOGL_) llhmin = -12345;
            hProfiles2d[key]->SetBinContent(bidx, bidy, llhmin);

            auto L = LLHMap.Filter(*p1+">"+std::to_string(bx_lo)+"&&"+*p1+"<"+std::to_string(bx_hi)+"&&"+*p2+">"+std::to_string(by_lo)+"&&"+*p2+"<"+std::to_string(by_hi)).Sum("L");

            hMarginals2d[key]->SetBinContent(bidx,bidy, *L);
          } // end bins y loop
        } // end bins x loop

        // Save the 2D histograms
        DirProfile2D->cd();
        hProfiles2d[key]->Write(hProfiles2d[key]->GetName(), TObject::kOverwrite);

        DirMarginal2D->cd();
        hMarginals2d[key]->Scale(1./hMarginals2d[key]->Integral());
        hMarginals2d[key]->Write(hMarginals2d[key]->GetName(), TObject::kOverwrite);
      } // end p2 loop
    } // end p1 loop

  } // end 2D plot

  OutFile->Close();

  return 0;
}
