//MaCh3 Includes
#include "plottingUtils/plottingUtils.h"
#include "plottingUtils/plottingManager.h"

#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PlotSigmaVariation.cpp
/// @todo Integrate within StylePlotting to get fancy labels etc
/// @todo Add support for 2D and by mode etc.
/// @author Kamil Skwarczynski

std::vector<std::string> DialNameVector;
std::vector<std::string> SampleNameVector;
std::vector<double> sigmaArray;

int PriorKnot = M3::_BAD_INT_;

constexpr const int NVars = 5;
constexpr Color_t Colours[NVars] = {kRed, kGreen+1, kBlack, kBlue+1, kOrange+1};
constexpr ELineStyle Style[NVars] = {kDotted, kDashed, kSolid, kDashDotted, kDashDotted};

/// @brief Scan inputs to figure out dial name and used sample names
void ScanInput(std::vector<std::string>& DialNameVecr,
               std::vector<std::string>& SampleNameVec,
               std::vector<double>& SigmaValues,
               const std::string& filename)
{
  MACH3LOG_DEBUG("Starting {}", __func__);
  TFile *infile = TFile::Open(filename.c_str());
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  //Get all entries in input file
  TIter next(SigmaDir->GetListOfKeys());
  TKey *key = nullptr;

  // Loop through all entries
  while ((key = static_cast<TKey*>(next()))) {
    // get directory names, ignore flux
    auto classname = std::string(key->GetClassName());
    auto dirname = std::string(key->GetName());

    if (classname != "TDirectoryFile") continue;
    dirname = std::string(key->GetName());

    SigmaDir->cd(dirname.c_str());
    TIter nextsub(gDirectory->GetListOfKeys());
    TKey *subkey = nullptr;

    DialNameVecr.push_back(dirname);
    MACH3LOG_DEBUG("Entering Dial {}", dirname);

    if(SampleNameVec.size() != 0) continue;
    //loop over items in directory, hard code which th2poly we want
    while ((subkey = static_cast<TKey*>(nextsub())))
    {
      auto subdirname = std::string(subkey->GetName());
      SampleNameVec.push_back(subdirname);
      MACH3LOG_DEBUG("Entering Sample {}", subdirname);

      if(SigmaValues.size() != 0) continue;
      SigmaDir->cd((dirname + "/" +  subdirname).c_str());

      TKey *subsubkey = nullptr;
      TIter nextsubsub(gDirectory->GetListOfKeys());
      while ((subsubkey = static_cast<TKey*>(nextsubsub())))
      {
        auto subsubdirname = std::string(subsubkey->GetTitle());
        MACH3LOG_DEBUG("Entering Hist {}", subsubdirname);
        std::string histname = subsubdirname;

        auto StripPrefix = [](std::string& str, const std::string& prefix) {
          if (str.find(prefix + "_") == 0) {
            str.erase(0, prefix.length() + 1);
            if (str.find(prefix) == 0) {
              MACH3LOG_ERROR("Failed to strip prefix '{}' from string '{}'", prefix, str);
              throw MaCh3Exception(__FILE__, __LINE__);
            }
          } else {
            MACH3LOG_ERROR("String '{}' does not start with expected prefix '{}'", str, prefix);
            throw MaCh3Exception(__FILE__, __LINE__);
          }
        };
        // Remove sample and dial name to avoid potential issues
        StripPrefix(histname, subdirname);
        StripPrefix(histname, dirname);

        MACH3LOG_DEBUG("Name afters striping {}", histname);

        double sigma = 0.0;
        // Find the "_sig_" part in the name
        size_t sig_pos = histname.find("_sig_");
        // Extract the part after "_sig_"
        std::string sigma_part = histname.substr(sig_pos + 5);

        // Check if it starts with 'n' (negative) or 'p' (positive) or 'nom' (0 or prior)
        if (histname.find("nom_") != std::string::npos) {
          sigma = 0.0;
          PriorKnot = static_cast<int>(SigmaValues.size());
          MACH3LOG_DEBUG("Found prior knot {}", PriorKnot);
        } else if (sigma_part.size() > 0 && sigma_part[0] == 'n') {
          sigma = -std::stod(sigma_part.substr(1));
        } else if (sigma_part.size() > 0 && sigma_part[0] == 'p') {
          sigma = std::stod(sigma_part.substr(1));
        } else {
          MACH3LOG_ERROR("Weirdly formatted string {}", sigma_part);
          throw MaCh3Exception(__FILE__ , __LINE__ );
        }
        MACH3LOG_DEBUG("Adding sigma {}", sigma);

        SigmaValues.push_back(sigma);
      }
    }
  }

  if(PriorKnot == M3::_BAD_INT_){
    MACH3LOG_ERROR("Didn't find prior knot, something is not right...");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  if(sigmaArray.size() < NVars){
    MACH3LOG_ERROR("Found sigma {}, while I have some hardcoding for {}",sigmaArray.size(), NVars);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  infile->Close();
  delete infile;
}

/// @brief Check whether to skip directory or not based on defined strings
bool SkipDirectory(const std::vector<std::string>& ExcludeString, const std::vector<std::string>& IncludeString, const std::string& dirname)
{
  bool Skip = false;
  for(unsigned int i = 0; i < ExcludeString.size(); i++)
  {
    if (dirname.find(ExcludeString[i]) != std::string::npos){ Skip = true; break; }
  }
  for(unsigned int i = 0; i < IncludeString.size(); i++)
  {
    if (!(dirname.find(IncludeString[i]) != std::string::npos)){ Skip = true; break; }
  }
  return Skip;
}

std::vector<double> GetDialValues(const std::vector<std::unique_ptr<TH1D>>& Poly) {
  std::vector<double> values;
  for (const auto& hist : Poly) {
      std::string title = hist->GetTitle();
      auto pos = title.rfind("_val_");
      if (pos != std::string::npos) {
        std::string val_str = title.substr(pos + 5); // skip "_val_"
        double val = std::stod(val_str);
        values.push_back(val);
        MACH3LOG_DEBUG("Extracted dial value {} from title '{}'", val, title);
      } else {
        MACH3LOG_DEBUG("Failed to extract dial value from title '{}'", title);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
  }
  return values;
}

void PlotRatio(const std::vector<std::unique_ptr<TH1D>>& Poly,
               const std::unique_ptr<TCanvas>& canv,
               const std::string& Title,
               const std::string& outfilename)
{
  canv->Clear();

  TPad* pad1 = new TPad("pad1","pad1",0.,0.25,1.,1.0);
  pad1->AppendPad();
  TPad* pad2 = new TPad("pad2","pad2",0.,0.,1.,0.25);
  pad2->AppendPad();

  pad1->SetLeftMargin(canv->GetLeftMargin());
  pad1->SetRightMargin(canv->GetRightMargin());
  pad1->SetTopMargin(canv->GetTopMargin());
  pad1->SetBottomMargin(0);

  pad2->SetLeftMargin(canv->GetLeftMargin());
  pad2->SetRightMargin(canv->GetRightMargin());
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.30);

  pad1->SetGrid();
  pad2->SetGrid();

  pad1->cd();

  auto DialValues = GetDialValues(Poly);
  double max = 0;
  for(int ik = 0; ik < static_cast<int>(sigmaArray.size()); ++ik)
  {
    Poly[ik]->SetLineWidth(2.);
    Poly[ik]->SetLineColor(Colours[ik]);
    Poly[ik]->SetLineStyle(Style[ik]);

    max = std::max(max, Poly[ik]->GetMaximum());
  }
  Poly[0]->SetTitle(Title.c_str());
  Poly[0]->SetMaximum(max*1.2);
  Poly[0]->Draw("HIST");
  for(int ik = 1; ik < static_cast<int>(sigmaArray.size()); ++ik)
  {
    Poly[ik]->Draw("HIST SAME");
  }

  std::vector<double> Integral(sigmaArray.size());
  for(int ik = 0; ik < static_cast<int>(sigmaArray.size()); ++ik)
    Integral[ik] = Poly[ik]->Integral();

  auto leg = std::make_unique<TLegend>(0.55, 0.55, 0.8, 0.88);
  leg->SetTextSize(0.04);
  for (int j = 0; j < static_cast<int>(sigmaArray.size()); j++)
  {
    if(j == PriorKnot) {
      leg->AddEntry(Poly[j].get(), Form("Prior (%.2f), #int=%.2f", DialValues[j], Integral[j]), "l");
    } else {
      leg->AddEntry(Poly[j].get(), Form("%.0f#sigma (%.2f), #int=%.2f", sigmaArray[j], DialValues[j], Integral[j]), "l");
    }
  }
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");

  pad2->cd();

  auto line = std::make_unique<TLine>(Poly[0]->GetXaxis()->GetBinLowEdge(Poly[0]->GetXaxis()->GetFirst()), 1.0, Poly[0]->GetXaxis()->GetBinUpEdge(Poly[0]->GetXaxis()->GetLast()), 1.0);
  std::vector<std::unique_ptr<TH1D>> Ratio(sigmaArray.size()-1);

  size_t ratio_index = 0; // Track position in Ratio vector
  for (int i = 0; i < static_cast<int>(Poly.size()); ++i) {
    if (i == PriorKnot) continue; // Skip PriorKnot
    Ratio[ratio_index] = M3::Clone(Poly[i].get());
    Ratio[ratio_index]->Divide(Poly[PriorKnot].get());
    ratio_index++;
  }

  Ratio[0]->GetYaxis()->SetTitle("Ratio to Prior");
  Ratio[0]->SetBit(TH1D::kNoTitle);
  Ratio[0]->GetXaxis()->SetTitleSize(0.12);
  Ratio[0]->GetYaxis()->SetTitleOffset(0.4);
  Ratio[0]->GetYaxis()->SetTitleSize(0.10);

  Ratio[0]->GetXaxis()->SetLabelSize(0.10);
  Ratio[0]->GetYaxis()->SetLabelSize(0.10);

  Ratio[0]->SetBit(TH1D::kNoTitle);

  double maxz = -999;
  double minz = +999;
  for (int j = 0; j < static_cast<int>(sigmaArray.size())-1; j++)
  {
    for (int i = 1; i < Ratio[0]->GetXaxis()->GetNbins(); i++)
    {
      maxz = std::max(maxz, Ratio[j]->GetBinContent(i));
      minz = std::min(minz, Ratio[j]->GetBinContent(i));
    }
  }
  maxz = maxz*1.001;
  minz = minz*1.001;

  if (std::fabs(1 - maxz) > std::fabs(1-minz))
    Ratio[0]->GetYaxis()->SetRangeUser(1-std::fabs(1-maxz),1+std::fabs(1-maxz));
  else
    Ratio[0]->GetYaxis()->SetRangeUser(1-std::fabs(1-minz),1+std::fabs(1-minz));

  Ratio[0]->Draw("HIST");

  for(int ik = 1; ik < static_cast<int>(sigmaArray.size())-1; ++ik)
  {
    Ratio[ik]->Draw("HIST SAME");
  }

  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  canv->Print((outfilename).c_str());

  delete pad1;
  delete pad2;
}

void CompareSigVar1D(const std::string& filename, const YAML::Node& Settings)
{
  //Get input file, make canvas and output file
  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1080, 1080);
  TFile *infile = TFile::Open(filename.c_str());
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  std::string outfilename = filename.substr(0, filename.find(".root"));
  outfilename = outfilename + "_RatioPlots1d.pdf";
  gErrorIgnoreLevel = kWarning;
  gStyle->SetDrawBorder(0);
  gStyle->SetTitleBorderSize(2);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  canvas->SetGrid();
  canvas->SetTopMargin(0.10);
  canvas->SetBottomMargin(0.08);
  canvas->SetRightMargin(0.05);
  canvas->SetLeftMargin(0.12);
  canvas->Print((outfilename+"[").c_str());

  auto IncludeString = GetFromManager<std::vector<std::string>>(Settings["IncludeString"], {});
  auto ExcludeString = GetFromManager<std::vector<std::string>>(Settings["ExcludeString"], {});
  TDirectory *dir = nullptr;

  for(size_t id = 0; id < DialNameVector.size(); id++)
  {
    for(size_t is = 0; is < SampleNameVector.size(); is++)
    {
      if(SkipDirectory(ExcludeString, IncludeString, (DialNameVector[id] + "/" + SampleNameVector[is]).c_str())) continue;
      MACH3LOG_INFO("Entering {}/{}", DialNameVector[id], SampleNameVector[is]);
      SigmaDir->cd((DialNameVector[id] + "/" +  SampleNameVector[is]).c_str());

      //set dir to current directory
      dir = gDirectory;

      //make -3,-1,0,1,3 polys
      std::vector<std::unique_ptr<TH1D>> Projection;
      TIter nextsub(dir->GetListOfKeys());
      TKey *subsubkey = nullptr;

      //loop over items in directory, hard code which th2poly we want
      while ((subsubkey = static_cast<TKey*>(nextsub())))
      {
        auto name = std::string(subsubkey->GetName());
        auto classname = std::string(subsubkey->GetClassName());

        if (classname == "TH1D")
        {
          name = DialNameVector[id] + "/" + SampleNameVector[is] + "/" + name;
          Projection.emplace_back(M3::Clone(SigmaDir->Get<TH1D>(name.c_str())));
        }
      }
      std::string Title = DialNameVector[id] + " " + SampleNameVector[is];
      PlotRatio(Projection, canvas, Title, outfilename);
      gDirectory->cd("..");
    }
  }

  canvas->Print((outfilename+"]").c_str());
  infile->Close();
  delete infile;
}

int main(int argc, char **argv)
{
  SetMaCh3LoggerFormat();
  if (argc != 3)
  {
    MACH3LOG_ERROR("Need two inputs: output of sigma var and config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  std::string filename = argv[1];
  std::string ConfigName = argv[2];
  MACH3LOG_INFO("Running {} with {} {}", argv[0], filename, ConfigName);
  // Load the YAML file
  YAML::Node Config = M3OpenConfig(ConfigName);

  // Access the "MatrixPlotter" section
  YAML::Node settings = Config["PlotSigmaVariation"];

  ScanInput(DialNameVector, SampleNameVector, sigmaArray, filename);

  CompareSigVar1D(filename, settings);

  return 0;
}
