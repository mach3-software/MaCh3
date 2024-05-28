//MaCh3 includes
#include "mcmc/MCMCProcessor.h"
#include "manager/manager.h"

inline void ProcessMCMC(std::string inputFile);
inline void MultipleProcessMCMC();
inline void CalcBayesFactor(MCMCProcessor* Processor);
inline void CalcSavageDickey(MCMCProcessor* Processor);
inline void CalcBipolarPlot(MCMCProcessor* Processor);
inline void GetTrianglePlot(MCMCProcessor* Processor);
inline void DiagnoseCovarianceMatrix(MCMCProcessor* Processor, std::string inputFile);
inline void ReweightPrior(MCMCProcessor* Processor);
inline TH2D* TMatrixIntoTH2D(TMatrixDSym* Matrix, std::string title);
inline void KolmogorovSmirnovTest(MCMCProcessor** Processor, TCanvas* Posterior, TString canvasname);

int nFiles;
std::vector <std::string> FileNames;
std::vector <std::string> TitleNames;
std::string config;

int main(int argc, char *argv[]) 
{
    nFiles = 0;
    if (argc != 3 && argc !=6 && argc != 8)
    {
      std::cerr << "How to use: "<< argv[0] <<"<Config> <MCMM_ND_Output.root>" << std::endl;
      exit(-1);
    }
  
    if (argc == 3)
    {
      std::cout << "Producing single fit output" << std::endl;
      config = argv[1];
      std::string filename = argv[2];
      ProcessMCMC(filename);
    } 
    // If we want to compare two or more fits (e.g. binning changes or introducing new params/priors)
    else if (argc == 6 || argc == 8)
    {
      std::cout << "Producing two fit comparison" << std::endl;
      config = argv[1];

      FileNames.push_back(argv[2]);
      TitleNames.push_back(argv[3]);
          
      FileNames.push_back(argv[4]);
      TitleNames.push_back(argv[5]);
      //KS: If there is third file add it
      if(argc == 8)
      {
        FileNames.push_back(argv[6]);
        TitleNames.push_back(argv[7]);
      }

      MultipleProcessMCMC();
    }
    
  return 0;
}

void ProcessMCMC(std::string inputFile)
{
  std::cout << "File for study: " << inputFile << " with config  "<<config<<std::endl;
    
  // Make the processor
  MCMCProcessor* Processor = new MCMCProcessor(inputFile, false);

  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  bool PlotCorr = GetFromManager<bool>(Settings["PlotCorr"], false);

  Processor->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["ExcludedTypes"], {""}));
  Processor->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["ExcludedNames"], {""}));
  //Apply additional cuts to 1D posterior
  Processor->SetPosterior1DCut(GetFromManager<std::string>(Settings["Posterior1DCut"], ""));

  if(PlotCorr) Processor->SetOutputSuffix("_drawCorr");
  //KS:Turn off plotting detector and some other setting, should be via some config
  Processor->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["PlotRelativeToPrior"], false));
  Processor->SetPrintToPDF(GetFromManager<bool>(Settings["PrintToPDF"], true));

  //KS: Whether you want prior error bands for parameters with flat prior or not
  Processor->SetPlotErrorForFlatPrior(GetFromManager<bool>(Settings["PlotErrorForFlatPrior"], true));
  Processor->SetFancyNames(GetFromManager<bool>(Settings["FancyNames"], true));
  Processor->SetPlotBinValue(GetFromManager<bool>(Settings["PlotBinValue"], false));
  //KS: Plot only 2D posteriors with correlations greater than 0.2
  Processor->SetPost2DPlotThreshold(GetFromManager<double>(Settings["Post2DPlotThreshold"], 0.2));

  Processor->Initialise();
  
  // Make the postfit
  Processor->MakePostfit();
  Processor->DrawPostfit();
  //KS: Should set via config whether you want below or not
  if(GetFromManager<bool>(Settings["MakeCredibleIntervals"], true)) {
    Processor->MakeCredibleIntervals(GetFromManager<std::vector<double>>(Settings["CredibleIntervals"], {0.99, 0.90, 0.68}),
                                     GetFromManager<std::vector<short int>>(Settings["CredibleIntervalsColours"], {436, 430, 422}),
                                     GetFromManager<bool>(Settings["CredibleInSigmas"], false));
  }
  if(GetFromManager<bool>(Settings["CalcBayesFactor"], true))  CalcBayesFactor(Processor);
  if(GetFromManager<bool>(Settings["CalcSavageDickey"], true)) CalcSavageDickey(Processor);
  if(GetFromManager<bool>(Settings["CalcBipolarPlot"], false)) CalcBipolarPlot(Processor);

  if(PlotCorr)
  {
    Processor->SetSmoothing(GetFromManager<bool>(Settings["Smoothing"], true));
    // Make the covariance matrix
    //We have different treatment for multithread
//#ifdef MULTITHREAD
    Processor->CacheSteps();
    //KS: Since we cached let's make fancy violins :)
    if(GetFromManager<bool>(Settings["MakeViolin"], true)) Processor->MakeViolin();
    Processor->MakeCovariance_MP();
//#else
    //Processor->MakeCovariance();
//#endif
    Processor->DrawCovariance();

    auto const &MakeSubOptimality = Settings["MakeSubOptimality"];
    if(MakeSubOptimality[0].as<bool>()) Processor->MakeSubOptimality(MakeSubOptimality[1].as<int>());

    if(GetFromManager<bool>(Settings["MakeCredibleRegions"], false)) {
      Processor->MakeCredibleRegions(GetFromManager<std::vector<double>>(Settings["CredibleRegions"], {0.99, 0.90, 0.68}),
                                     GetFromManager<std::vector<short int>>(Settings["CredibleRegionStyle"], {2, 1, 3}),
                                     GetFromManager<std::vector<short int>>(Settings["CredibleRegionColor"], {413, 406, 416}),
                                     GetFromManager<bool>(Settings["CredibleInSigmas"], false)
                                     );
    }
    if(GetFromManager<bool>(Settings["GetTrianglePlot"], true)) GetTrianglePlot(Processor);

    //KS: When creating covariance matrix longest time is spend on caching every step, since we already cached we can run some fancy covariance stability diagnostic
    if(GetFromManager<bool>(Settings["DiagnoseCovarianceMatrix"], false)) DiagnoseCovarianceMatrix(Processor, inputFile);
  }
  if(GetFromManager<bool>(Settings["ReweightPrior"], false)) ReweightPrior(Processor);

  delete Processor;
}

void MultipleProcessMCMC()
{
  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  const Color_t PosteriorColor[] = {kBlue-1, kRed, kGreen+2};
  //const Style_t PosteriorStyle[] = {kSolid, kDashed, kDotted};
  nFiles = FileNames.size();
  MCMCProcessor** Processor; 
  Processor = new MCMCProcessor*[nFiles];
  for (int ik = 0; ik < nFiles;  ik++)
  {
    std::cout << "File for study:       " << FileNames[ik] << std::endl;
    // Make the processor
    Processor[ik] = new MCMCProcessor(FileNames[ik], false);
    Processor[ik]->SetOutputSuffix(("_" + std::to_string(ik)).c_str());

    Processor[ik]->SetExcludedTypes(GetFromManager<std::vector<std::string>>(Settings["ExcludedTypes"], {""}));
    Processor[ik]->SetExcludedNames(GetFromManager<std::vector<std::string>>(Settings["ExcludedNames"], {""}));

    //Apply additional cuts to 1D posterior
    Processor[ik]->SetPosterior1DCut(GetFromManager<std::string>(Settings["Posterior1DCut"], ""));

    Processor[ik]->SetPlotRelativeToPrior(GetFromManager<bool>(Settings["PlotRelativeToPrior"], false));
    Processor[ik]->SetFancyNames(GetFromManager<bool>(Settings["FancyNames"], true));
    Processor[ik]->Initialise();
  }
  //KS: Multithreading here is very tempting but there are some issues with root that need to be resovled :(
  for (int ik = 0; ik < nFiles;  ik++)
  {
    // Make the postfit
    Processor[ik]->MakePostfit();
    Processor[ik]->DrawPostfit();
  }

  // Open a TCanvas to write the posterior onto
  TCanvas* Posterior = new TCanvas("PosteriorMulti", "PosteriorMulti", 0, 0, 1024, 1024);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Posterior->SetGrid();
  Posterior->SetBottomMargin(0.1);
  Posterior->SetTopMargin(0.05);
  Posterior->SetRightMargin(0.03);
  Posterior->SetLeftMargin(0.10);

  FileNames[0] = FileNames[0].substr(0, FileNames[0].find(".root")-1);
  TString canvasname = FileNames[0];
  for (int ik = 1; ik < nFiles;  ik++)
  {
    while (FileNames[ik].find("/") != std::string::npos)
    {
      FileNames[ik] = FileNames[ik].substr(FileNames[ik].find("/")+1, FileNames[ik].find(".root")-1);
    }
    canvasname = canvasname + "_"+FileNames[ik];
  }

  canvasname = canvasname +".pdf[";

  Posterior->Print(canvasname);
  // Once the pdf file is open no longer need to bracket
  canvasname.ReplaceAll("[","");

  for(int i = 0; i < Processor[0]->GetNParams(); ++i) 
  {
    // This holds the posterior density
    TH1D **hpost = new TH1D*[nFiles];
    TLine **hpd = new TLine*[nFiles];
    hpost[0] = (TH1D *) (Processor[0]->GetHpost(i))->Clone();

    bool Skip = false;
    for (int ik = 1 ; ik < nFiles;  ik++)
    {
      // KS: If somehow this chain doesn't given params we skip it
      const int Index = Processor[ik]->GetParamIndexFromName(hpost[0]->GetTitle());
      if(Index == _UNDEF_)
      {
        Skip = true;
        break;
      }
      hpost[ik] = (TH1D *)(Processor[ik]->GetHpost(Index))->Clone();
    }

    // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if(hpost[0]->GetMaximum() == hpost[0]->Integral()*1.5 || Skip)
    {
      for (int ik = 0; ik < nFiles;  ik++)
        delete hpost[ik];

      delete[] hpost;
      delete[] hpd;
      continue;
    }
    for (int ik = 0; ik < nFiles;  ik++)
    {
      RemoveFitter(hpost[ik], "Gauss");

      // Set some nice colours
      hpost[ik]->SetLineColor(PosteriorColor[ik]);
      //hpost[ik]->SetLineStyle(PosteriorStyle[ik]);
      hpost[ik]->SetLineWidth(2);

      // Area normalise the distributions
      hpost[ik]->Scale(1./hpost[ik]->Integral(), "width");
    }
    TString Title;
    double Prior = 1.0;
    double PriorError = 1.0;

    Processor[0]->GetNthParameter(i, Prior, PriorError, Title);

    // Now make the TLine for the Asimov
    TLine *Asimov = new TLine(Prior,  hpost[0]->GetMinimum(), Prior,  hpost[0]->GetMaximum());
    Asimov->SetLineColor(kRed-3);
    Asimov->SetLineWidth(2);
    Asimov->SetLineStyle(kDashed);

    // Make a nice little TLegend
    TLegend *leg = new TLegend(0.12, 0.7, 0.6, 0.97);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    TString asimovLeg = Form("#splitline{Prior}{x = %.2f , #sigma = %.2f}", Prior, PriorError);
    leg->AddEntry(Asimov, asimovLeg, "l");

    for (int ik = 0; ik < nFiles;  ik++)
    {
      TString rebinLeg = Form("#splitline{%s}{#mu = %.2f, #sigma = %.2f}", TitleNames[ik].c_str(), hpost[ik]->GetMean(), hpost[ik]->GetRMS());
      leg->AddEntry(hpost[ik],  rebinLeg, "l");

      hpd[ik] = new TLine(hpost[ik]->GetBinCenter(hpost[ik]->GetMaximumBin()), hpost[ik]->GetMinimum(), hpost[ik]->GetBinCenter(hpost[ik]->GetMaximumBin()), hpost[ik]->GetMaximum());
      hpd[ik]->SetLineColor(hpost[ik]->GetLineColor());
      hpd[ik]->SetLineWidth(2);
      hpd[ik]->SetLineStyle(kSolid);
    }

    // Find the maximum value to nicley resize hist
    double maximum = 0;
    for (int ik = 0; ik < nFiles;  ik++) maximum = std::max(maximum, hpost[ik]->GetMaximum());
    for (int ik = 0; ik < nFiles;  ik++) hpost[ik]->SetMaximum(1.3*maximum);

    hpost[0]->Draw("hist");
    for (int ik = 1; ik < nFiles;  ik++) hpost[ik]->Draw("hist same");
    Asimov->Draw("same");
    for (int ik = 0; ik < nFiles;  ik++) hpd[ik]->Draw("same");
    leg->Draw("same");
    Posterior->cd();
    Posterior->Print(canvasname);

    delete Asimov;
    delete leg;
    for (int ik = 0; ik < nFiles;  ik++)
    {
      delete hpost[ik];
      delete hpd[ik];
    }
    delete[] hpost;
    delete[] hpd;
  }//End loop over paramters
    
  // Finally draw the parameter plot onto the PDF
  // Close the .pdf file with all the posteriors
  Posterior->cd();
  Posterior->Clear();

  if(GetFromManager<bool>(Settings["PerformKStest"], true)) KolmogorovSmirnovTest(Processor, Posterior, canvasname);
  
  // Close the pdf file
  std::cout << "Closing pdf " << canvasname << std::endl;
  canvasname+="]";
  Posterior->Print(canvasname);
  
  delete Posterior;
  for (int ik = 0; ik < nFiles;  ik++) delete Processor[ik];
  delete[] Processor;
}

// KS: Calculate Bayes factor for a given hiphothesis, most onformative are those related to osc params. However, it make realtive easy interpreation for switch dials
void CalcBayesFactor(MCMCProcessor* Processor)
{
  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  std::vector<std::string> ParNames;
  std::vector<std::vector<double>> Model1Bounds;
  std::vector<std::vector<double>> Model2Bounds;
  std::vector<std::vector<std::string>> ModelNames;
  for (const auto& dg : Settings["BayesFactor"])
  {
    ParNames.push_back(dg[0].as<std::string>());
    ModelNames.push_back(dg[1].as<std::vector<std::string>>());
    Model1Bounds.push_back(dg[2].as<std::vector<double>>());
    Model2Bounds.push_back(dg[3].as<std::vector<double>>());
  }

  Processor->GetBayesFactor(ParNames, Model1Bounds, Model2Bounds, ModelNames);
  return;
}

void CalcSavageDickey(MCMCProcessor* Processor)
{
  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  std::vector<std::string> ParNames;
  std::vector<double> EvaluationPoint;
  std::vector<std::vector<double>> Bounds;
  
  for (const auto& d : Settings["SavageDickey"])
  {
    ParNames.push_back(d[0].as<std::string>());
    EvaluationPoint.push_back(d[1].as<double>());
    Bounds.push_back(d[2].as<std::vector<double>>());
  }

  Processor->GetSavageDickey(ParNames, EvaluationPoint, Bounds);
  return;
}


void CalcBipolarPlot(MCMCProcessor* Processor)
{
  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  std::vector<std::string> ParNames;

  for (const auto& d : Settings["BipolarPlot"])
  {
    ParNames.push_back(d[0].as<std::string>());
  }
  Processor->GetPolarPlot(ParNames);
  return;
}


void GetTrianglePlot(MCMCProcessor* Processor){
  std::cout<<std::endl;

  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  for (const auto& dg : Settings["TrianglePlot"])
  {
    std::string ParName = dg[0].as<std::string>();

    std::vector<std::string> NameVec = dg[1].as<std::vector<std::string>>();
    Processor->MakeTrianglePlot(NameVec,
                                GetFromManager<std::vector<double>>(Settings["CredibleIntervals"], {0.99, 0.90, 0.68}),
                                GetFromManager<std::vector<short int>>(Settings["CredibleIntervalsColours"], {436, 430, 422}),
                                GetFromManager<std::vector<double>>(Settings["CredibleRegions"], {0.99, 0.90, 0.68}),
                                GetFromManager<std::vector<short int>>(Settings["CredibleRegionStyle"], {2, 1, 3}),
                                GetFromManager<std::vector<short int>>(Settings["CredibleRegionColor"], {413, 406, 416}),
                                GetFromManager<bool>(Settings["CredibleInSigmas"], false));
  }
}

//KS: You validate stability of posterior covariance matrix, you set burn calc cov matrix increase burn calc again and compare. By performing such operation several hundred times we can check when matrix becomes stable
void DiagnoseCovarianceMatrix(MCMCProcessor* Processor, std::string inputFile)
{
  //Turn of plots from Processor
  Processor->SetPrintToPDF(false);
  // Open a TCanvas to write the posterior onto
  TCanvas* Canvas = new TCanvas("Canvas", "Canvas", 0, 0, 1024, 1024);
  Canvas->SetGrid();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Canvas->SetTickx();
  Canvas->SetTicky();
  Canvas->SetBottomMargin(0.1);
  Canvas->SetTopMargin(0.05);
  Canvas->SetRightMargin(0.15);
  Canvas->SetLeftMargin(0.10);
  
  //KS: Fancy colots
  const int NRGBs = 10;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.35, 0.50, 0.60, 0.65, 0.75, 0.90, 1.00 };
  Double_t red[NRGBs]   = { 0.50, 1.00, 1.00, 0.25, 0.00, 0.10, 0.50, 1.00, 0.75, 0.55 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00, 0.60, 0.90, 1.00, 0.75, 0.75 };
  Double_t blue[NRGBs]  = { 0.00, 0.25, 1.00, 1.00, 0.50, 0.60, 0.90, 1.00, 0.05, 0.05 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
  
  std::string OutName = inputFile;
  OutName = OutName.substr(0, OutName.find(".root"));
  Canvas->Print(Form("Correlation_%s.pdf[", OutName.c_str()), "pdf");
  Canvas->Print(Form("Covariance_%s.pdf[", OutName.c_str()), "pdf");
  
  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  const int entries = Processor->GetnSteps();
  const int NIntervals = GetFromManager<int>(Settings["NIntervals"], 5);
  const int IntervalsSize = entries/NIntervals;
  //We start with burn from 0 (no burn in at all)
  int BurnIn = 0;
  std::cout<<"Diagnosing matrices with entries="<< entries<<", NIntervals="<<NIntervals<<" and IntervalsSize="<<IntervalsSize<<std::endl;

  TMatrixDSym *Covariance = nullptr;
  TMatrixDSym *Correlation = nullptr;
  
  TH2D *CovariancePreviousHist = nullptr;
  TH2D *CorrelationPreviousHist = nullptr;

  TH2D *CovarianceHist = nullptr;
  TH2D *CorrelationHist = nullptr;

  //KS: Get first covariances, we need two for comparison...
  Processor->SetStepCut(BurnIn);
  Processor->GetCovariance(Covariance, Correlation);
  
  CovariancePreviousHist = TMatrixIntoTH2D(Covariance, "Covariance"); 
  CorrelationPreviousHist = TMatrixIntoTH2D(Correlation, "Correlation");
      
  delete Covariance;
  Covariance = nullptr;
  delete Correlation;
  Correlation = nullptr;
  
  //KS: Loop over alls dsired cuts
  for(int k = 1; k < NIntervals; ++k)
  {
    BurnIn = k*IntervalsSize;
    Processor->SetStepCut(BurnIn);
    Processor->GetCovariance(Covariance, Correlation);
    Processor->ResetHistograms();
    
    CovarianceHist = TMatrixIntoTH2D(Covariance, "Covariance"); 
    CorrelationHist = TMatrixIntoTH2D(Correlation, "Correlation");
            
    TH2D *CovarianceDiff = (TH2D*)CovarianceHist->Clone("Covariance_Ratio");
    TH2D *CorrelationDiff = (TH2D*)CorrelationHist->Clone("Correlation_Ratio");
    
    //KS: Bit messy but quite often covariance is 0 is divided by 0 is problemiatic so
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int j = 1; j < CovarianceDiff->GetXaxis()->GetNbins()+1; ++j)
    {
      for (int i = 1; i < CovarianceDiff->GetYaxis()->GetNbins()+1; ++i)
      {
        if( std::fabs (CovarianceDiff->GetBinContent(j, i)) < 1.e-5 && std::fabs (CovariancePreviousHist->GetBinContent(j, i)) < 1.e-5)
        {
          CovarianceDiff->SetBinContent(j, i, _UNDEF_);
          CovariancePreviousHist->SetBinContent(j, i, _UNDEF_);
        }
        if( std::fabs (CorrelationDiff->GetBinContent(j, i)) < 1.e-5 && std::fabs (CorrelationPreviousHist->GetBinContent(j, i)) < 1.e-5)
        {
          CorrelationDiff->SetBinContent(j, i, _UNDEF_);
          CorrelationPreviousHist->SetBinContent(j, i, _UNDEF_);
        }
      }
    }
    //Divide matrices
    CovarianceDiff->Divide(CovariancePreviousHist);
    CorrelationDiff->Divide(CorrelationPreviousHist);
    
    //Now it is time for fancy names etc.
    for (int j = 0; j < CovarianceDiff->GetXaxis()->GetNbins(); ++j)
    {
      TString Title = "";
      double Prior = 1.0;
      double PriorError = 1.0;
  
      Processor->GetNthParameter(j, Prior, PriorError, Title);
      
      CovarianceDiff->GetXaxis()->SetBinLabel(j+1, Title);
      CovarianceDiff->GetYaxis()->SetBinLabel(j+1, Title);
      CorrelationDiff->GetXaxis()->SetBinLabel(j+1, Title);
      CorrelationDiff->GetYaxis()->SetBinLabel(j+1, Title);
    }
    CovarianceDiff->GetXaxis()->SetLabelSize(0.015);
    CovarianceDiff->GetYaxis()->SetLabelSize(0.015);
    CorrelationDiff->GetXaxis()->SetLabelSize(0.015);
    CorrelationDiff->GetYaxis()->SetLabelSize(0.015);
    
    std::stringstream ss;
    ss << "BCut_";
    ss << BurnIn;
    ss << "/";
    ss << "BCut_";
    ss << (k-1)*IntervalsSize;
    std::string str = ss.str();
    
    TString Title = "Cov " + str;
    CovarianceDiff->GetZaxis()->SetTitle( Title );
    Title = "Corr " + str;
    CorrelationDiff->GetZaxis()->SetTitle(Title);
    
    CovarianceDiff->SetMinimum(-2);
    CovarianceDiff->SetMaximum(2);
    CorrelationDiff->SetMinimum(-2);
    CorrelationDiff->SetMaximum(2);

    Canvas->cd();
    CovarianceDiff->Draw("colz");
    Canvas->Print(Form("Covariance_%s.pdf", OutName.c_str()), "pdf");

    CorrelationDiff->Draw("colz");
    Canvas->Print(Form("Correlation_%s.pdf", OutName.c_str()), "pdf");
        
    //KS: Current hist become previous as we need it for further comparison
    delete CovariancePreviousHist;
    CovariancePreviousHist = (TH2D*)CovarianceHist->Clone();
    delete CorrelationPreviousHist;
    CorrelationPreviousHist = (TH2D*)CorrelationHist->Clone();;
    
    delete CovarianceHist;
    CovarianceHist = nullptr;
    delete CorrelationHist;
    CorrelationHist = nullptr;
    
    delete CovarianceDiff;
    delete CorrelationDiff;
    delete Covariance;
    Covariance = nullptr;
    delete Correlation;
    Correlation = nullptr;
  }
  Canvas->cd();
  Canvas->Print(Form("Covariance_%s.pdf]", OutName.c_str()), "pdf");
  Canvas->Print(Form("Correlation_%s.pdf]", OutName.c_str()), "pdf");
  
  Processor->SetPrintToPDF(true);
  if(Covariance != nullptr)              delete Covariance;
  if(Correlation != nullptr)             delete Correlation;
  if(CovariancePreviousHist != nullptr)  delete CovariancePreviousHist;
  if(CorrelationPreviousHist != nullptr) delete CorrelationPreviousHist;
  if(CovarianceHist != nullptr)          delete CovarianceHist;
  if(CorrelationHist != nullptr)         delete CorrelationHist;
  delete Canvas;
}

void ReweightPrior(MCMCProcessor* Processor)
{
  std::cout<<std::endl;

  YAML::Node card_yaml = YAML::LoadFile(config.c_str());
  YAML::Node Settings = card_yaml["ProcessMCMC"];

  const auto& Prior = Settings["PriorReweighting"];

  std::vector<std::string> Names = Prior[0].as<std::vector<std::string>>();
  std::vector<double> NewCentral = Prior[1].as<std::vector<double>>();
  std::vector<double> NewError = Prior[2].as<std::vector<double>>();

  Processor->ReweightPrior(Names, NewCentral, NewError);
}

//KS: Convert TMatrix to TH2D, mostly useful for making fancy plots
TH2D* TMatrixIntoTH2D(TMatrixDSym* Matrix, std::string title)       
{
  TH2D* hMatrix = new TH2D(title.c_str(), title.c_str(), Matrix->GetNrows(), 0.0, Matrix->GetNrows(), Matrix->GetNcols(), 0.0, Matrix->GetNcols());
  for(int i = 0; i < Matrix->GetNrows(); i++)
  {
    for(int j = 0; j < Matrix->GetNcols(); j++)
    {
      //KS: +1 becasue there is offset in histogram realtive to TMatrix
      hMatrix->SetBinContent(i+1,j+1, (*Matrix)(i,j));
    }
  }
  
  return hMatrix;
}

//KS: Perform KS test to check if two posteriors for the same parameter came from the same distribution
void KolmogorovSmirnovTest(MCMCProcessor** Processor, TCanvas* Posterior, TString canvasname)
{
  const Color_t CumulativeColor[] = {kBlue-1, kRed, kGreen+2};
  const Style_t CumulativeStyle[] = {kSolid, kDashed, kDotted};

  for(int i = 0; i < Processor[0]->GetNParams(); ++i) 
  {
    // This holds the posterior density
    TH1D **hpost = new TH1D*[nFiles];
    TH1D **CumulativeDistribution = new TH1D*[nFiles];
       
    TString Title;
    double Prior = 1.0;
    double PriorError = 1.0;

    Processor[0]->GetNthParameter(i, Prior, PriorError, Title);
    bool Skip = false;
    for (int ik = 0 ; ik < nFiles;  ik++)
    {
      int Index = 0;
      if(ik == 0 ) Index = i;
      else
      {
        // KS: If somehow this chain doesn't given params we skip it
        Index = Processor[ik]->GetParamIndexFromName(hpost[0]->GetTitle());
        if(Index == _UNDEF_)
        {
          Skip = true;
          break;
        }
      }
      hpost[ik] = (TH1D*) (Processor[ik]->GetHpost(Index))->Clone();
      CumulativeDistribution[ik] = (TH1D*) (Processor[ik]->GetHpost(Index))->Clone();
      CumulativeDistribution[ik]->Fill(0., 0.);
      CumulativeDistribution[ik]->Reset();
      CumulativeDistribution[ik]->SetMaximum(1.);
      TString TempTittle = Title+" Kolmogorov Smirnov";
      CumulativeDistribution[ik]->SetTitle(TempTittle);
      
      TempTittle = Title+" Value";
      CumulativeDistribution[ik]->GetXaxis()->SetTitle(TempTittle);
      CumulativeDistribution[ik]->GetYaxis()->SetTitle("Cumulative Probability");
      
      CumulativeDistribution[ik]->SetLineWidth(2);
      CumulativeDistribution[ik]->SetLineColor(CumulativeColor[ik]);
      CumulativeDistribution[ik]->SetLineStyle(CumulativeStyle[ik]);
    }

    // Don't plot if this is a fixed histogram (i.e. the peak is the whole integral)
    if(hpost[0]->GetMaximum() == hpost[0]->Integral()*1.5 || Skip)
    {
      for (int ik = 0; ik < nFiles;  ik++)
      {
        delete hpost[ik];
        delete CumulativeDistribution[ik];
      }
      delete[] hpost;
      delete[] CumulativeDistribution;
      continue;
    }
    
    for (int ik = 0 ; ik < nFiles;  ik++)
    {
      const int NumberOfBins = hpost[ik]->GetXaxis()->GetNbins();
      double Cumulative = 0;
      const double Integral = hpost[ik]->Integral();
      for (int j = 1; j < NumberOfBins+1; ++j)
      {
        Cumulative += hpost[ik]->GetBinContent(j)/Integral;
        
        CumulativeDistribution[ik]->SetBinContent(j, Cumulative);
      }
      //KS: Set overflow to 1 just in case
      CumulativeDistribution[ik]->SetBinContent(NumberOfBins+1, 1.);
    }
    
    int* TestStatBin = new int[nFiles];
    double* TestStatD = new double[nFiles];
    TLine **LineD = new TLine*[nFiles];

    for (int ik = 0 ; ik < nFiles;  ik++) { TestStatBin[ik] = 0; TestStatD[ik] = -999;}

    //Find KS statistic
    for (int ik = 1 ; ik < nFiles;  ik++)
    {
      const int NumberOfBins = CumulativeDistribution[0]->GetXaxis()->GetNbins();
      for (int j = 1; j < NumberOfBins+1; ++j)
      {
        double BinValue = CumulativeDistribution[0]->GetBinCenter(j);
        int BinNumber = CumulativeDistribution[ik]->FindBin(BinValue);
        //KS: Calculate D statistic for this bin, only save it if it's bigger than previosly found value
        double TempDstat = std::fabs(CumulativeDistribution[0]->GetBinContent(j) - CumulativeDistribution[ik]->GetBinContent(BinNumber));
        if(TempDstat > TestStatD[ik])
        {
          TestStatD[ik] = TempDstat;
          TestStatBin[ik] = j;
        }
      }
    }

    for (int ik = 0 ; ik < nFiles;  ik++)
    {
      LineD[ik] =  new TLine(CumulativeDistribution[0]->GetBinCenter(TestStatBin[ik]), 0, CumulativeDistribution[0]->GetBinCenter(TestStatBin[ik]), CumulativeDistribution[0]->GetBinContent(TestStatBin[ik]));
      LineD[ik]->SetLineColor(CumulativeColor[ik]);
      LineD[ik]->SetLineWidth(2.0);
    }
    CumulativeDistribution[0]->Draw();
    for (int ik = 0 ; ik < nFiles;  ik++)
      CumulativeDistribution[ik]->Draw("SAME");
    
    TLegend *leg = new TLegend(0.15, 0.7, 0.5, 0.90);
    leg->SetTextSize(0.04);
    for (int ik = 0; ik < nFiles;  ik++)
      leg->AddEntry(CumulativeDistribution[ik], TitleNames[ik].c_str(), "l");
    for (int ik = 1; ik < nFiles;  ik++)
      leg->AddEntry(LineD[ik], Form("#Delta D = %.4f", TestStatD[ik]), "l");
    
    leg->SetLineColor(0);
    leg->SetLineStyle(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->Draw("SAME");

    for (int ik = 1; ik < nFiles;  ik++)
      LineD[ik]->Draw("sam");
      
    Posterior->cd();
    Posterior->Print(canvasname);
    
    delete leg;
    for (int ik = 0; ik < nFiles;  ik++)
    {
        delete hpost[ik];
        delete CumulativeDistribution[ik];
        delete LineD[ik];
    }
    delete[] hpost;
    delete[] CumulativeDistribution;
    delete[] LineD;
    delete[] TestStatBin;
    delete[] TestStatD;
  } //End loop over parameter
}
