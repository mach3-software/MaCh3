//MaCh3 Includes
#include "PlottingUtils/PlottingUtils.h"
#include "PlottingUtils/PlottingManager.h"

#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PlotSigmaVariation.cpp
/// @author Kamil Skwarczynski

std::vector<std::string> DialNameVector;
std::vector<std::string> SampleNameVector;
std::vector<int> SampleMaxDim;
std::vector<double> sigmaArray;

int PriorKnot = M3::_BAD_INT_;

constexpr const int NVars = 5;
constexpr Color_t Colours[NVars] = {kRed, kGreen+1, kBlack, kBlue+1, kOrange+1};
constexpr ELineStyle Style[NVars] = {kDotted, kDashed, kSolid, kDashDotted, kDashDotted};

/// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
MaCh3Plotting::PlottingManager* PlotMan;

/// @brief Histograms have name like ND_CC0pi_1DProj0_Norm_Param_0_sig_n3.00_val_0.25. This code is trying to extract sigma names
void FindKnot(std::vector<double>& SigmaValues,
              const std::string& dirname,
              const std::string& subdirname,
              const std::string& ProjName,
              std::string histname) {
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
  StripPrefix(histname, ProjName);
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


/// @brief Scan inputs to figure out dial name and used sample names
void ScanInput(std::vector<std::string>& DialNameVecr,
               std::vector<std::string>& SampleNameVec,
               std::vector<int>& SampleDimVec,
               std::vector<double>& SigmaValues,
               const std::string& filename)
{
  MACH3LOG_DEBUG("Starting {}", __func__);
  TFile *infile = M3::Open(filename, "OPEN", __FILE__, __LINE__);
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  //Get all entries in input file
  TIter next(SigmaDir->GetListOfKeys());
  TKey *key = nullptr;

  // Loop over directory with dial names
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
    //loop over directories with sample names
    while ((subkey = static_cast<TKey*>(nextsub())))
    {
      auto subdirname = std::string(subkey->GetName());
      SampleNameVec.push_back(subdirname);
      SampleDimVec.push_back(0);
      MACH3LOG_DEBUG("Entering Sample {}", subdirname);
      SigmaDir->cd((dirname + "/" +  subdirname).c_str());

      TKey *subsubkey = nullptr;
      TIter nextsubsub(gDirectory->GetListOfKeys());

      // Check if we already filled sigma vector
      bool FillSigma = false;
      if(SigmaValues.size() == 0) FillSigma = true;
      // loop over histograms
      while ((subsubkey = static_cast<TKey*>(nextsubsub())))
      {
        auto subsubdirname = std::string(subsubkey->GetTitle());
        MACH3LOG_DEBUG("Entering Hist {}", subsubdirname);
        std::string histname = subsubdirname;

        classname = std::string(subsubkey->GetClassName());

        if (classname != "TH1D") continue;
        // Find if there is more dimensions
        size_t proj_pos = histname.find("_1DProj");
        int proj_number = -1;

        if (proj_pos != std::string::npos) {
          size_t number_start = proj_pos + 7; // skip "_1DProj"
          size_t number_end = histname.find_first_not_of("0123456789", number_start);
          proj_number = std::stoi(histname.substr(number_start, number_end - number_start));
        }
        SampleDimVec.back() = std::max(proj_number, SampleDimVec.back());
        MACH3LOG_DEBUG("Found dimension {} with dimension", proj_number);

        // KS: Extract knot position from hist only we haven't done this before and for projection X
        // Sigma are same for projection Y and Z and beyond
        if(FillSigma && proj_number == 0) FindKnot(SigmaValues, dirname, subdirname, "1DProj" + std::to_string(proj_number), histname);
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

  if(SampleDimVec.size() != SampleDimVec.size()) {
    MACH3LOG_ERROR("Sample name vec ({}) and sample dimension vec ({}) have different sizes, something is not right");
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

/// @brief Extracts dial value for from histogram title
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


void MakeRatio(const std::vector<std::unique_ptr<TH1D>>& Poly,
               std::vector<std::unique_ptr<TH1D>>& Ratio) {
  size_t ratio_index = 0; // Track position in Ratio vector
  for (int i = 0; i < static_cast<int>(Poly.size()); ++i) {
    if (i == PriorKnot) continue; // Skip PriorKnot
    Ratio[ratio_index] = M3::Clone(Poly[i].get());
    Ratio[ratio_index]->Divide(Poly[PriorKnot].get());
    ratio_index++;
  }

  Ratio[0]->GetYaxis()->SetTitle("Ratio to Prior");
  Ratio[0]->SetBit(TH1D::kNoTitle);
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
}


void PlotRatio(const std::vector<std::unique_ptr<TH1D>>& Poly,
               const std::unique_ptr<TCanvas>& canv,
               const std::string& Title,
               const std::string& outfilename)
{
  canv->Clear();
  gStyle->SetDrawBorder(0);
  gStyle->SetTitleBorderSize(2);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  canv->SetGrid();
  canv->SetTopMargin(0.10);
  canv->SetBottomMargin(0.08);
  canv->SetRightMargin(0.05);
  canv->SetLeftMargin(0.12);

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

  TLine line(Poly[0]->GetXaxis()->GetBinLowEdge(Poly[0]->GetXaxis()->GetFirst()),
             1.0, Poly[0]->GetXaxis()->GetBinUpEdge(Poly[0]->GetXaxis()->GetLast()), 1.0);
  std::vector<std::unique_ptr<TH1D>> Ratio(sigmaArray.size()-1);

  MakeRatio(Poly, Ratio);

  auto PrettyX = PlotMan->style().prettifyKinematicName(Ratio[0]->GetXaxis()->GetTitle());
  Ratio[0]->GetXaxis()->SetTitle(PrettyX.c_str());
  Ratio[0]->GetXaxis()->SetTitleSize(0.12);
  Ratio[0]->GetYaxis()->SetTitleOffset(0.4);
  Ratio[0]->GetYaxis()->SetTitleSize(0.10);

  Ratio[0]->GetXaxis()->SetLabelSize(0.10);
  Ratio[0]->GetYaxis()->SetLabelSize(0.10);
  Ratio[0]->Draw("HIST");
  for(int ik = 1; ik < static_cast<int>(sigmaArray.size())-1; ++ik)
  {
    Ratio[ik]->Draw("HIST SAME");
  }

  line.SetLineWidth(2);
  line.SetLineColor(kBlack);
  line.Draw("SAME");

  canv->Print((outfilename).c_str());

  delete pad1;
  delete pad2;
}

void CompareSigVar1D(const std::string& filename, const YAML::Node& Settings)
{
  //Get input file, make canvas and output file
  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1080, 1080);
  TFile *infile = M3::Open(filename, "OPEN", __FILE__, __LINE__);
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  std::string outfilename = filename.substr(0, filename.find(".root"));
  outfilename = outfilename + "_RatioPlots1d.pdf";
  gErrorIgnoreLevel = kWarning;
  canvas->Print((outfilename+"[").c_str());

  auto IncludeString = GetFromManager<std::vector<std::string>>(Settings["IncludeString"], {});
  auto ExcludeString = GetFromManager<std::vector<std::string>>(Settings["ExcludeString"], {});
  TDirectory *dir = nullptr;

  for(size_t id = 0; id < DialNameVector.size(); id++)
  {
    for(size_t is = 0; is < SampleNameVector.size(); is++)
    {
      if(SkipDirectory(ExcludeString, IncludeString, (DialNameVector[id] + "/" + SampleNameVector[is]).c_str())) continue;
      MACH3LOG_INFO("{} Entering {}/{}", __func__, DialNameVector[id], SampleNameVector[is]);
      SigmaDir->cd((DialNameVector[id] + "/" +  SampleNameVector[is]).c_str());

      //set dir to current directory
      dir = gDirectory;

      // Loop over dimensions
      for(int iDim = 0; iDim <= SampleMaxDim[is]; iDim++)
      {
        MACH3LOG_DEBUG("Starting loop over dimension {}", iDim);
        //make -3,-1,0,1,3 polys
        std::vector<std::unique_ptr<TH1D>> Projection;
        TIter nextsub(dir->GetListOfKeys());
        TKey *subsubkey = nullptr;

        //loop over items in directory, hard code which th2poly we want
        while ((subsubkey = static_cast<TKey*>(nextsub())))
        {
          auto name = std::string(subsubkey->GetName());
          auto classname = std::string(subsubkey->GetClassName());
          // Looking
          const std::string ProjectionName = "_1DProj" + std::to_string(iDim);
          const bool IsProjection = (name.find(ProjectionName) != std::string::npos);
          if (classname == "TH1D" && IsProjection)
          {
            name = DialNameVector[id] + "/" + SampleNameVector[is] + "/" + name;
            Projection.emplace_back(M3::Clone(SigmaDir->Get<TH1D>(name.c_str())));
            MACH3LOG_DEBUG("Adding hist {}", name);
          }
        }
        std::string Title = PlotMan->style().prettifyParamName(DialNameVector[id]) + " "
                            + PlotMan->style().prettifySampleName(SampleNameVector[is]);
        PlotRatio(Projection, canvas, Title, outfilename);
        gDirectory->cd("..");
      }
    }
  }

  canvas->Print((outfilename+"]").c_str());
  infile->Close();
  delete infile;
}

void PlotRatio2D(const std::vector<std::unique_ptr<TH2>>& Poly,
                 const std::unique_ptr<TCanvas>& canv,
                 const std::string& Title,
                 const std::string& outfilename)
{
  canv->Clear();
  gStyle->SetDrawBorder(0);
  gStyle->SetTitleBorderSize(2);
  gStyle->SetOptStat(0); //Set 0 to disable statistic box
  canv->SetGrid();
  canv->SetTopMargin(0.10);
  canv->SetBottomMargin(0.10);
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.20);

  constexpr int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  for (int i = 0; i < static_cast<int>(Poly.size()); ++i) {
    if (i == PriorKnot) continue; // Skip PriorKnot
    std::unique_ptr<TH2> Ratio = M3::Clone(Poly[i].get());
    Ratio->Divide(Poly[PriorKnot].get());
    Ratio->SetTitle((Title + " " + std::to_string(static_cast<int>(sigmaArray[i])) + "sigma").c_str());

    const double maxz = Ratio->GetMaximum();
    const double minz = Ratio->GetMinimum();
    if (std::fabs(1-maxz) > std::fabs(1-minz))
      Ratio->GetZaxis()->SetRangeUser(1-std::fabs(1-maxz),1+std::fabs(1-maxz));
    else
      Ratio->GetZaxis()->SetRangeUser(1-std::fabs(1-minz),1+std::fabs(1-minz));
    Ratio->GetXaxis()->SetTitleOffset(1.1);
    Ratio->GetYaxis()->SetTitleOffset(1.1);
    Ratio->GetZaxis()->SetTitleOffset(1.5);
    Ratio->GetZaxis()->SetTitle("Ratio to Prior");

    auto PrettyX = PlotMan->style().prettifyKinematicName(Ratio->GetXaxis()->GetTitle());
    Ratio->GetXaxis()->SetTitle(PrettyX.c_str());
    auto PrettyY = PlotMan->style().prettifyKinematicName(Ratio->GetYaxis()->GetTitle());
    Ratio->GetYaxis()->SetTitle(PrettyY.c_str());

    Ratio->Draw("COLZ");
    canv->Print((outfilename).c_str());
  }
}

void CompareSigVar2D(const std::string& filename, const YAML::Node& Settings)
{
  //Get input file, make canvas and output file
  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1080, 1080);
  TFile *infile = M3::Open(filename, "OPEN", __FILE__, __LINE__);
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  std::string outfilename = filename.substr(0, filename.find(".root"));
  outfilename = outfilename + "_RatioPlots2d.pdf";
  gErrorIgnoreLevel = kWarning;
  canvas->Print((outfilename+"[").c_str());

  auto IncludeString = GetFromManager<std::vector<std::string>>(Settings["IncludeString"], {});
  auto ExcludeString = GetFromManager<std::vector<std::string>>(Settings["ExcludeString"], {});
  TDirectory *dir = nullptr;

  for(size_t id = 0; id < DialNameVector.size(); id++)
  {
    for(size_t is = 0; is < SampleNameVector.size(); is++)
    {
      if(SkipDirectory(ExcludeString, IncludeString, (DialNameVector[id] + "/" + SampleNameVector[is]).c_str())) continue;
      MACH3LOG_INFO("{} Entering {}/{}", __func__, DialNameVector[id], SampleNameVector[is]);
      SigmaDir->cd((DialNameVector[id] + "/" +  SampleNameVector[is]).c_str());

      //set dir to current directory
      dir = gDirectory;

      int nDim = SampleMaxDim[is];

      TIter nextsub(dir->GetListOfKeys());
      TKey *subsubkey = nullptr;

      // Loop over all unique dimension pairs
      for (int iDim1 = 0; iDim1 <= nDim; ++iDim1) {
        for (int iDim2 = iDim1 + 1; iDim2 <= nDim; ++iDim2) {

          // Reset iterator for each dimension pair
          nextsub.Reset();
          //make -3,-1,0,1,3 polys
          std::vector<std::unique_ptr<TH2>> Projection;

          //loop over items in directory, hard code which th2poly we want
          while ((subsubkey = static_cast<TKey*>(nextsub())))
          {
            auto name = std::string(subsubkey->GetName());
            auto classname = std::string(subsubkey->GetClassName());
            // Looking
            const std::string ProjectionName = "_2DProj_" + std::to_string(iDim1) + "_vs_" + std::to_string(iDim2);
            const bool IsProjection = (name.find(ProjectionName) != std::string::npos);
            if ((classname == "TH2D" ||  classname == "TH2Poly")&& IsProjection)
            {
              name = DialNameVector[id] + "/" + SampleNameVector[is] + "/" + name;
              Projection.emplace_back(M3::Clone(SigmaDir->Get<TH2>(name.c_str())));
              MACH3LOG_DEBUG("Adding hist {}", name);
            }
          }
          std::string Title = PlotMan->style().prettifyParamName(DialNameVector[id]) + " "
                              + PlotMan->style().prettifySampleName(SampleNameVector[is]);
          if(Projection.size() == sigmaArray.size()) PlotRatio2D(Projection, canvas, Title, outfilename);
        }
      }
      gDirectory->cd("..");
    }
  }
  canvas->Print((outfilename+"]").c_str());
  infile->Close();
  delete infile;
}


void PlotEventRate(const std::vector<std::vector<std::unique_ptr<TH1D>>>& Poly,
                   const std::unique_ptr<TCanvas>& canv,
                   const std::string& Title,
                   const std::string& outfilename)
{
  std::vector<std::unique_ptr<TH1D>> EvenRates(sigmaArray.size());
  for(int ih = 0; ih < static_cast<int>(sigmaArray.size()); ih++)
  {
    EvenRates[ih] = std::make_unique<TH1D>(Title.c_str(), Title.c_str(), SampleNameVector.size(), 0, SampleNameVector.size());
    EvenRates[ih]->SetDirectory(nullptr);
    EvenRates[ih]->GetYaxis()->SetTitleOffset(1.4);
    EvenRates[ih]->GetYaxis()->SetTitle("Events");
    EvenRates[ih]->SetLineWidth(2.);
    EvenRates[ih]->SetLineStyle(Style[ih]);
    EvenRates[ih]->SetLineColor(Colours[ih]);

    for(size_t iSample = 0; iSample < SampleNameVector.size(); iSample ++) {
      EvenRates[ih]->SetBinContent(iSample+1, Poly[iSample][ih]->Integral());

      std::string SamName = PlotMan->style().prettifySampleName(SampleNameVector[iSample]);
      EvenRates[ih]->GetXaxis()->SetBinLabel(iSample+1, SamName.c_str());
    }
  }

  TPad* pad1 = new TPad("pad1","pad1",0.,0.50,1.,1.0);
  pad1->AppendPad();
  TPad* pad2 = new TPad("pad2","pad2",0.,0.,1.,0.50);
  pad2->AppendPad();

  pad1->SetLeftMargin(canv->GetLeftMargin());
  pad1->SetRightMargin(canv->GetRightMargin());
  pad1->SetTopMargin(canv->GetTopMargin());
  pad1->SetBottomMargin(0);

  pad2->SetLeftMargin(canv->GetLeftMargin());
  pad2->SetRightMargin(canv->GetRightMargin());
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.60);

  pad1->SetGrid();
  pad2->SetGrid();

  pad1->cd();
  double max = 0;
  for(int ik = 0; ik < static_cast<int>(sigmaArray.size()); ++ik) {
    max = std::max(max, EvenRates[ik]->GetMaximum());
  }
  EvenRates[0]->SetTitle(Title.c_str());
  EvenRates[0]->SetMaximum(max*1.2);
  EvenRates[0]->Draw("HIST");
  for(int ik = 1; ik < static_cast<int>(sigmaArray.size()); ++ik)
  {
    EvenRates[ik]->Draw("HIST SAME");
  }

  auto DialValues = GetDialValues(Poly[0]);
  auto leg = std::make_unique<TLegend>(0.55, 0.55, 0.8, 0.88);
  leg->SetTextSize(0.04);
  for (int j = 0; j < static_cast<int>(sigmaArray.size()); j++)
  {
    if(j == PriorKnot) {
      leg->AddEntry(EvenRates[j].get(), Form("Prior (%.2f)", DialValues[j]), "l");
    } else {
      leg->AddEntry(EvenRates[j].get(), Form("%.0f#sigma (%.2f)", sigmaArray[j], DialValues[j]), "l");
    }
  }
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");

  pad2->cd();
  TLine line(EvenRates[0]->GetXaxis()->GetBinLowEdge(EvenRates[0]->GetXaxis()->GetFirst()),
             1.0, EvenRates[0]->GetXaxis()->GetBinUpEdge(EvenRates[0]->GetXaxis()->GetLast()), 1.0);
  std::vector<std::unique_ptr<TH1D>> Ratio(sigmaArray.size()-1);
  MakeRatio(EvenRates, Ratio);
  Ratio[0]->GetXaxis()->SetTitleSize(0.08);
  Ratio[0]->GetYaxis()->SetTitleOffset(0.4);
  Ratio[0]->GetYaxis()->SetTitleSize(0.06);

  Ratio[0]->GetXaxis()->SetLabelSize(0.08);
  Ratio[0]->GetYaxis()->SetLabelSize(0.04);
  Ratio[0]->GetXaxis()->LabelsOption("v");
  Ratio[0]->Draw("HIST");

  for(int ik = 1; ik < static_cast<int>(sigmaArray.size())-1; ++ik)
  {
    Ratio[ik]->Draw("HIST SAME");
  }

  line.SetLineWidth(2);
  line.SetLineColor(kBlack);
  line.Draw("SAME");

  canv->Print((outfilename).c_str());

  delete pad1;
  delete pad2;
}

void MakeEventRatePlot(const std::string& filename, const YAML::Node& Settings)
{
  (void) Settings;
  //Get input file, make canvas and output file
  auto canvas = std::make_unique<TCanvas>("canv", "canv", 1080, 1080);
  TFile *infile = M3::Open(filename, "OPEN", __FILE__, __LINE__);
  TDirectoryFile *SigmaDir = infile->Get<TDirectoryFile>("SigmaVar");

  std::string outfilename = filename.substr(0, filename.find(".root"));
  outfilename = outfilename + "_EventRate.pdf";
  gErrorIgnoreLevel = kWarning;
  canvas->Print((outfilename+"[").c_str());

  TDirectory *dir = nullptr;
  for(size_t id = 0; id < DialNameVector.size(); id++)
  {
    std::vector<std::vector<std::unique_ptr<TH1D>>> Projection(SampleNameVector.size());
    for(size_t is = 0; is < SampleNameVector.size(); is++)
    {
      MACH3LOG_INFO("{} Entering {}/{}", __func__, DialNameVector[id], SampleNameVector[is]);
      SigmaDir->cd((DialNameVector[id] + "/" +  SampleNameVector[is]).c_str());

      //set dir to current directory
      dir = gDirectory;

      //make -3,-1,0,1,3 polys
      TIter nextsub(dir->GetListOfKeys());
      TKey *subsubkey = nullptr;

      //loop over items in directory, hard code which th2poly we want
      while ((subsubkey = static_cast<TKey*>(nextsub())))
      {
        auto name = std::string(subsubkey->GetName());
        auto classname = std::string(subsubkey->GetClassName());
        // Looking
        const std::string ProjectionName = "_1DProj0";
        const bool IsProjection = (name.find(ProjectionName) != std::string::npos);
        if (classname == "TH1D" && IsProjection)
        {
          name = DialNameVector[id] + "/" + SampleNameVector[is] + "/" + name;
          Projection[is].emplace_back(M3::Clone(SigmaDir->Get<TH1D>(name.c_str())));
          MACH3LOG_DEBUG("Adding hist {}", name);
        }
      }
      gDirectory->cd("..");
    }
    std::string Title = PlotMan->style().prettifyParamName(DialNameVector[id]);
    PlotEventRate(Projection, canvas, Title, outfilename);
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

  ScanInput(DialNameVector, SampleNameVector, SampleMaxDim, sigmaArray, filename);

  PlotMan = new MaCh3Plotting::PlottingManager();
  PlotMan->initialise();

  CompareSigVar1D(filename, settings);
  CompareSigVar2D(filename, settings);
  MakeEventRatePlot(filename, settings);

  if(PlotMan) delete PlotMan;
  return 0;
}
