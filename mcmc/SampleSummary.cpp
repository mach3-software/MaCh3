#include "mcmc/SampleSummary.h"

// *******************
// The constructor
SampleSummary::SampleSummary(const int n_Samples, const std::string &Filename, samplePDFBase* const sample, const int nSteps) {
// *******************

  std::cout << "Making sample summary class..." << std::endl;
  #ifdef MULTITHREAD
  std::cout << "with OpenMP and " << omp_get_max_threads() << " threads" << std::endl;
  #endif
  
  StandardFluctuation = true;
  
  if(StandardFluctuation) std::cout << "Using standard method of statistical fluctuation" << std::endl;
  else std::cout << "Using alternative method of statistical fluctuation, which is much slower" << std::endl;
  
  //KS: If true it will print posterior predictive for every beta parameter it is quite useful but make root big number of plots
  DoBetaParam = true;
  if(DoBetaParam) std::cout << "I will calculate #beta parameters from Barlow-Beeston formalism" << std::endl;

  //If true code will normalise each histogram, this way you can calculate shape only error. etc. pvalue will be completely wrong unfortunately
  doShapeOnly = false;

  nChainSteps = nSteps;
  //KS: nChainSteps == 0 means we run PriorPredcitive
  if(nChainSteps == 0) isPriorPredictive = true;
  else isPriorPredictive = false;
  
  OutputName = Filename;
  nSamples = n_Samples;
  SamplePDF = sample;

  //Get mach3 modes from manager
  Modes = SamplePDF->GetMaCh3Modes();

  nThrows = 0;
  first_pass = true;
  Outputfile = nullptr;
  OutputTree = nullptr;
  Dir = nullptr;

  rnd = new TRandom3();

  DataHist = new TH2Poly*[nSamples];
  DataHist_ProjectX = new TH1D*[nSamples];
  DataHist_ProjectY = new TH1D*[nSamples];
  NominalHist = new TH2Poly*[nSamples];
  PosteriorHist = new TH1D**[nSamples];
  W2NomHist = new TH2Poly*[nSamples];
  w2Hist = new TH1D**[nSamples];

  ViolinHists_ProjectX = new TH2D*[nSamples];
  ViolinHists_ProjectY = new TH2D*[nSamples];
    
  if(DoBetaParam) BetaHist = new TH1D**[nSamples];
  else BetaHist = nullptr;

  maxBins = new int[nSamples];
  
  lnLHist_Mean = new TH2Poly*[nSamples];
  lnLHist_Mode = new TH2Poly*[nSamples];
  lnLHist_Mean_ProjectX = new TH1D*[nSamples];
  MeanHist = new TH2Poly*[nSamples];
  if(DoBetaParam) MeanHistCorrected = new TH2Poly*[nSamples];
  else MeanHistCorrected = NULL;
  ModeHist = new TH2Poly*[nSamples];
  W2MeanHist = new TH2Poly*[nSamples];
  W2ModeHist = new TH2Poly*[nSamples];
  lnLHist_Mean1D = new TH1D*[nSamples];
  lnLHist_Mode1D = new TH1D*[nSamples];
  lnLHist_Sample_DrawData = new TH1D*[nSamples];
  lnLHist_Sample_DrawflucDraw = new TH1D*[nSamples];
  lnLHist_Sample_PredflucDraw = new TH1D*[nSamples];
    
  //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automatically adjust histogram range
  // https://root.cern.ch/doc/master/classTH1.html#auto-bin
  lnLHist = new TH1D("lnLHist_predpredfluc", "lnLHist_predpredfluc", 100, 1, -1);
  lnLHist->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Pred)");
  lnLHist->GetYaxis()->SetTitle("Counts");

  lnLHist_drawdata = new TH1D("lnLHist_drawdata", "lnLHist_drawdata", 100, 1, -1);
  lnLHist_drawdata->GetXaxis()->SetTitle("-2LLH (Data, Draw)");
  lnLHist_drawdata->GetYaxis()->SetTitle("Counts");

  lnLHist_drawfluc = new TH1D("lnLHist_drawpredfluc", "lnLHist_drawpredfluc", 100, 1, -1);
  lnLHist_drawfluc->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Draw)");
  lnLHist_drawfluc->GetYaxis()->SetTitle("Counts");

  lnLHist_drawflucdraw = new TH1D("lnLHist_drawflucdraw", "lnLHist_drawflucdraw", 100, 1, -1);
  lnLHist_drawflucdraw->GetXaxis()->SetTitle("-2LLH (Draw Fluc, Draw)");
  lnLHist_drawflucdraw->GetYaxis()->SetTitle("Counts");

  lnLDrawHist = new TH2D("lnLDrawHist", "lnLDrawHist", 50, 1, -1, 50, 1, -1);
  lnLDrawHist->GetXaxis()->SetTitle("-2LLH_{Pred Fluc, Draw}");
  lnLDrawHist->GetYaxis()->SetTitle("-2LLH_{Data, Draw}");

  lnLFlucHist = new TH2D("lnLFlucHist", "lnLFlucHist", 50, 1, -1, 50, 1, -1);
  lnLFlucHist->GetXaxis()->SetTitle("-2LLH_{Draw Fluc, Draw}");
  lnLFlucHist->GetYaxis()->SetTitle("-2LLH_{Data, Draw}");

  lnLDrawHistRate = new TH2D("lnLDrawHistRate", "lnLDrawHistRate", 50, 1, -1, 50, 1, -1);
  lnLDrawHistRate->GetXaxis()->SetTitle("-2LLH_{Pred Fluc, Draw}");
  lnLDrawHistRate->GetYaxis()->SetTitle("-2LLH_{Data, Draw}");

  //KS: This is silly as it assumes all samples uses same kinematics
  lnLFlucHist_ProjectX = new TH2D("lnLFlucHist_ProjectX", "lnLFlucHist_ProjectX", 50, 1, -1, 50, 1, -1);
  lnLFlucHist_ProjectX->GetXaxis()->SetTitle(("-2LLH_{Draw Fluc, Draw} for " + SamplePDF->GetKinVarLabel(0, 0)).c_str());
  lnLFlucHist_ProjectX->GetYaxis()->SetTitle(("-2LLH_{Data, Draw} for " + SamplePDF->GetKinVarLabel(0, 0)).c_str());
  
  // Holds the hist of random number draws, only works for posterior predictive
  if(!isPriorPredictive)
  {
    RandomHist = new TH1D("RandomHist", "RandomHist", 100, 0, nChainSteps);
    RandomHist->GetXaxis()->SetTitle("Step");
    const double binwidth = nChainSteps/RandomHist->GetNbinsX();
    std::stringstream ss;
    ss << "Draws/" << binwidth;
    RandomHist->GetYaxis()->SetTitle(ss.str().c_str());
    RandomHist->SetLineWidth(2);
  }
  else RandomHist = NULL;

  for (_int_ i = 0; i < nSamples; ++i)
  {
    PosteriorHist[i] = NULL;
    w2Hist[i] = NULL;
      
    if(DoBetaParam) BetaHist[i] = nullptr;
    DataHist[i] = NULL;
    DataHist_ProjectX[i] = NULL;
    DataHist_ProjectY[i] = NULL;
    NominalHist[i] = NULL;
    
    ViolinHists_ProjectX[i] = NULL;
    ViolinHists_ProjectY[i] = NULL;

    MeanHist[i] = NULL;
    if(DoBetaParam) MeanHistCorrected[i] = NULL;
    W2MeanHist[i] = NULL;
    W2ModeHist[i] = NULL;
    lnLHist_Mean[i] = NULL;
    lnLHist_Mode[i] = NULL;
    lnLHist_Mean_ProjectX[i] = NULL;
    lnLHist_Mean1D[i] = NULL;
    lnLHist_Mode1D[i] = NULL;
    lnLHist_Sample_DrawData[i] = NULL;
    lnLHist_Sample_DrawflucDraw[i] = NULL;
    lnLHist_Sample_PredflucDraw[i] = NULL;
  }//end loop over samples
  
  llh_data_draw = NULL;
  llh_drawfluc_draw = NULL;
  llh_predfluc_draw = NULL;
  llh_rate_data_draw = NULL;
  llh_rate_predfluc_draw = NULL;

  llh_data_drawfluc = NULL;
  llh_data_predfluc = NULL;
  llh_draw_pred = NULL;
  llh_drawfluc_pred = NULL;

  llh_predfluc_pred = NULL;
  llh_drawfluc_predfluc = NULL;
  llh_datafluc_draw = NULL;

  DoByModePlots = false;
  MeanHist_ByMode = NULL;
  PosteriorHist_ByMode = NULL;

  nModelParams = 0;

  Debug = 0;
}

// *******************
//  Destructor
SampleSummary::~SampleSummary() {
// *******************
  Outputfile->cd();

  //ROOT is wierd and once you write TFile claim ownership of histograms. Best is to first delete histograms and then close file
  Outputfile->Close();
  delete Outputfile;

  if(rnd != nullptr) delete rnd;
  
  if(lnLHist != NULL) delete lnLHist;
  if(lnLHist_drawdata != NULL) delete lnLHist_drawdata;
  if(lnLHist_drawfluc != NULL) delete lnLHist_drawfluc;
  if(lnLHist_drawflucdraw != NULL) delete lnLHist_drawflucdraw;
  if(lnLDrawHist != NULL) delete lnLDrawHist;
  if(lnLFlucHist != NULL) delete lnLFlucHist;
  if(lnLDrawHistRate != NULL) delete lnLDrawHistRate;
  if(RandomHist != NULL)  delete RandomHist;

  if(lnLFlucHist_ProjectX != NULL) delete lnLFlucHist_ProjectX;
  if(DoByModePlots)
  {
    for (_int_ i = 0; i < nSamples; ++i)
    {
      if(DataHist[i] == NULL) continue;
      for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
      {
        for (int k = 1; k <= maxBins[i]; ++k)
        {
          if(PosteriorHist_ByMode[i][j][k] != NULL) delete PosteriorHist_ByMode[i][j][k];
        }
        delete[] PosteriorHist_ByMode[i][j];
        if(MeanHist_ByMode[i][j] != NULL) delete MeanHist_ByMode[i][j];
      }
      delete[] PosteriorHist_ByMode[i];
      delete[] MeanHist_ByMode[i];
    }
    delete[] PosteriorHist_ByMode;
    delete[] MeanHist_ByMode;
  }

  for (_int_ i = 0; i < nSamples; ++i)
  {
    if(DataHist[i] == NULL) continue;
    if(DataHist[i] != NULL) delete DataHist[i];
    if(NominalHist[i] != NULL) delete NominalHist[i];
    if(MeanHist[i] != NULL) delete MeanHist[i];
    if(DoBetaParam && MeanHistCorrected[i] != NULL) delete MeanHistCorrected[i];
    if(W2MeanHist[i] != NULL) delete W2MeanHist[i];
    if(W2ModeHist[i] != NULL) delete W2ModeHist[i];
    
    if(ViolinHists_ProjectX[i] != NULL) delete ViolinHists_ProjectX[i];
    if(ViolinHists_ProjectY[i] != NULL) delete ViolinHists_ProjectY[i];
    
    if(lnLHist_Mean[i] != NULL) delete lnLHist_Mean[i];
    if(lnLHist_Mode[i] != NULL) delete lnLHist_Mode[i];
    if(lnLHist_Mean_ProjectX[i] != NULL) delete lnLHist_Mean_ProjectX[i];
    if(lnLHist_Mean1D[i] != NULL) delete lnLHist_Mean1D[i];
    if(lnLHist_Mode1D[i] != NULL) delete lnLHist_Mode1D[i];
    if(lnLHist_Sample_DrawData[i] != NULL) delete lnLHist_Sample_DrawData[i];
    if(lnLHist_Sample_DrawflucDraw[i] != NULL) delete lnLHist_Sample_DrawflucDraw[i];
    if(lnLHist_Sample_PredflucDraw[i] != NULL) delete lnLHist_Sample_PredflucDraw[i];

    for (int j = 1; j <= maxBins[i]; ++j)
    {
      if(PosteriorHist[i][j] != NULL) delete PosteriorHist[i][j];
      if(w2Hist[i][j] != NULL) delete w2Hist[i][j];
      if(DoBetaParam && BetaHist[i][j] != NULL) delete BetaHist[i][j];
    }
    delete[] PosteriorHist[i];
    delete[] w2Hist[i];
    if(DoBetaParam) delete[] BetaHist[i];
  }

  delete[] DataHist;
  delete[] NominalHist;
  delete[] MeanHist;
  if(DoBetaParam) delete[] MeanHistCorrected;
  delete[] W2MeanHist;
  delete[] W2ModeHist;

  delete[] ViolinHists_ProjectX;
  delete[] ViolinHists_ProjectY;
      
  delete[] lnLHist_Mean;
  delete[] lnLHist_Mode;
  delete[] lnLHist_Mean_ProjectX;
  delete[] lnLHist_Mean1D;
  delete[] lnLHist_Mode1D;
  delete[] lnLHist_Sample_DrawData;
  delete[] lnLHist_Sample_DrawflucDraw;
  delete[] lnLHist_Sample_PredflucDraw;
  
  delete[] maxBins;
      
  delete[] PosteriorHist;
  delete[] w2Hist;
  if(DoBetaParam) delete[] BetaHist;

  delete[] llh_data_draw;
  delete[] llh_data_drawfluc;
  delete[] llh_data_predfluc;
  delete[] llh_rate_data_draw;
  delete[] llh_rate_predfluc_draw;
  delete[] llh_draw_pred;
  delete[] llh_drawfluc_pred;
  delete[] llh_drawfluc_predfluc;
  delete[] llh_drawfluc_draw;
  delete[] llh_predfluc_pred;
  delete[] llh_predfluc_draw;
  delete[] llh_datafluc_draw;
  
  delete[] llh_data_draw_ProjectX;
  delete[] llh_drawfluc_draw_ProjectX;
}

// *******************
// Check size of sample against size of vectors
bool SampleSummary::CheckSamples(int Length) {
// *******************
  bool ok = (nSamples == Length);
  if (!ok) {
    std::cerr << "Size of SampleVector input != number of defined samples" << std::endl;
    std::cout << "Size of SampleVector: " << Length << std::endl;
    std::cout << "Size of defined samples: " << nSamples << std::endl;
    std::cerr << "Something has gone wrong with making the Samples" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  return ok;
}

// *******************
// Add a data histogram to the list (will have N_samples of these)
// Since the data doesn't change with varying the MC
void SampleSummary::AddData(std::vector<TH2Poly*> &Data) {
// *******************
  const int Length = Data.size();
  // Check length of samples are OK
  if (!CheckSamples(Length)) throw;
  for (int i = 0; i < Length; ++i) {
    if (Data[i] == NULL) {
      DataHist[i] = NULL;
      DataHist_ProjectX[i] = NULL;
      DataHist_ProjectY[i] = NULL;
      maxBins[i] = 0;
    } else {
      std::string classname = std::string(DataHist[i]->Class_Name());
      if(classname == "TH2Poly")
      {
        DataHist[i] = (TH2Poly*)(Data[i]->Clone());
        if(doShapeOnly) NormaliseTH2Poly(DataHist[i]);
        DataHist_ProjectX[i] = ProjectPoly(DataHist[i], true, i);
        DataHist_ProjectY[i] = ProjectPoly(DataHist[i], false, i);
        maxBins[i] = DataHist[i]->GetNumberOfBins();
      } else {
        std::cerr<<"Somehow sample "<<SamplePDF->GetSampleName(i)<<"doesn't use TH2Poly"<<std::endl;
        std::cerr<<"Right now I only support TH2Poly but I am ambitious piece of code and surely will have more support in the future"<<std::endl;
        throw;
      }
    }
  }
}

// *******************
// Add the nominal histograms to the list (will have N_samples of these)
void SampleSummary::AddNominal(std::vector<TH2Poly*> &Nominal, std::vector<TH2Poly*> &NomW2) {
// *******************

  const int Length = Nominal.size();
  if (!CheckSamples(Length)) throw;
  
  //KS: ROOT is super annoying and you cannot use clone with openMP, hence we have another loop below
  for (int i = 0; i < Length; ++i) 
  {
    if (Nominal[i] == NULL) {
      PosteriorHist[i] = NULL;
      NominalHist[i] = NULL;
      W2NomHist[i] = NULL;
      lnLHist_Mean[i] = NULL;
      lnLHist_Mode[i] = NULL;
      lnLHist_Mean_ProjectX[i] = NULL;
      MeanHist[i] = NULL;
      if(DoBetaParam) MeanHistCorrected[i] = NULL;
      ModeHist[i] = NULL;
      W2MeanHist[i] = NULL;
      W2ModeHist[i] = NULL;
      lnLHist_Sample_DrawData[i] = NULL;
      lnLHist_Sample_DrawflucDraw[i] = NULL;
      lnLHist_Sample_PredflucDraw[i] = NULL;
    // If not NULL it indicates the selection was turned on, so initialise the privates
    } else {
      NominalHist[i] = (TH2Poly*)(Nominal[i]->Clone());
      if(doShapeOnly) NormaliseTH2Poly(NominalHist[i]);
      W2NomHist[i] = (TH2Poly*)(NomW2[i]->Clone());
      
      lnLHist_Mean[i] = (TH2Poly*)(NominalHist[i]->Clone());
      lnLHist_Mean[i]->SetDirectory(0);
      lnLHist_Mode[i] = (TH2Poly*)(NominalHist[i]->Clone());
      lnLHist_Mode[i]->SetDirectory(0);
      lnLHist_Mean_ProjectX[i] = (TH1D*)(DataHist_ProjectX[i]->Clone());
      MeanHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      if(DoBetaParam) MeanHistCorrected[i] = (TH2Poly*)(NominalHist[i]->Clone());
      ModeHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      W2MeanHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      W2ModeHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
    }
  }
  
  // Loop over the length of nominal and set the initial distributions up
  //KS: Don't multithread, mostly due to fact that we initialise histograms
  for (int i = 0; i < Length; ++i) {
  // If NULL it indicates the selection was turned off, so initialise all the hists to NULL
    if (Nominal[i] != NULL) 
    {
      std::string name = std::string(NominalHist[i]->GetName());
      name = name.substr(0, name.find("_nom"));

      PosteriorHist[i] = new TH1D*[maxBins[i]+1];
      w2Hist[i] = new TH1D*[maxBins[i]+1];

      if(DoBetaParam) BetaHist[i] = new TH1D*[maxBins[i]+1];

      for (int j = 0; j <= maxBins[i]; ++j)
      {
        PosteriorHist[i][j] = NULL;
        w2Hist[i][j] = NULL;

        if(DoBetaParam) BetaHist[i][j] = nullptr;
      }
      lnLHist_Mean[i]->SetNameTitle((name+"_MeanlnL").c_str(), (name+"_MeanlnL").c_str());
      lnLHist_Mean[i]->Reset("");
      lnLHist_Mean[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");
      
      lnLHist_Mode[i]->SetNameTitle((name+"_ModelnL").c_str(), (name+"_ModelnL").c_str());
      lnLHist_Mode[i]->Reset("");
      lnLHist_Mode[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");

      lnLHist_Mean_ProjectX[i]->SetNameTitle((name+"_MeanlnL_ProjectX").c_str(), (name+"_MeanlnL_ProjectX").c_str()); 
      lnLHist_Mean_ProjectX[i]->Reset("");
      lnLHist_Mean_ProjectX[i]->GetYaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");

      MeanHist[i]->SetNameTitle((name+"_mean").c_str(), (name+"_mean").c_str());
      MeanHist[i]->Reset("");
      MeanHist[i]->GetZaxis()->SetTitle("Mean");
      
      if(DoBetaParam)
      {
        MeanHistCorrected[i]->SetNameTitle((name+"_mean_corrected").c_str(), (name+"_mean_corrected").c_str());
        MeanHistCorrected[i]->Reset("");
        MeanHistCorrected[i]->GetZaxis()->SetTitle("Mean");
      }
      std::vector<double> xbins;
      std::vector<double> ybins;

      SamplePDF->SetupBinning(i, xbins, ybins);
      
      //KS: Y axis is number of events to get estimate of maximal number we use integral
      const int MaxBinning = (doShapeOnly) ? 1 : NoOverflowIntegral(NominalHist[i])/4;
      ViolinHists_ProjectX[i] = new TH2D((name+"_Violin_ProjectX").c_str(), (name+"_Violin_ProjectX").c_str(), xbins.size()-1, &xbins[0] , 400, 0, MaxBinning);
      ViolinHists_ProjectX[i]->GetYaxis()->SetTitle("Events");
      ViolinHists_ProjectX[i]->GetXaxis()->SetTitle(std::string(NominalHist[i]->GetXaxis()->GetTitle()).c_str() );
      ViolinHists_ProjectX[i]->SetDirectory(0);

      ViolinHists_ProjectY[i] = new TH2D((name+"_Violin_ProjectY").c_str(), (name+"_Violin_ProjectY").c_str(), ybins.size()-1, &ybins[0] , 400, 0, MaxBinning);
      ViolinHists_ProjectY[i]->GetYaxis()->SetTitle("Events");
      ViolinHists_ProjectY[i]->GetXaxis()->SetTitle(std::string(NominalHist[i]->GetYaxis()->GetTitle()).c_str());
      ViolinHists_ProjectY[i]->SetDirectory(0);

      ModeHist[i]->SetNameTitle((name+"_mode").c_str(), (name+"_mode").c_str());
      ModeHist[i]->Reset("");
      ModeHist[i]->GetZaxis()->SetTitle("Mode");

      W2MeanHist[i]->SetNameTitle((name+"_w2mean").c_str(), (name+"_w2mean").c_str());
      W2MeanHist[i]->Reset("");
      W2MeanHist[i]->GetZaxis()->SetTitle("W2Mean");

      W2ModeHist[i]->SetNameTitle((name+"_w2mode").c_str(), (name+"_w2mode").c_str());
      W2ModeHist[i]->Reset("");
      W2ModeHist[i]->GetZaxis()->SetTitle("W2Mode");

      // Declare the lnL histograms
      lnLHist_Mean1D[i] = new TH1D((name+"_MeanlnL1D").c_str(),(name+"_MeanlnL1D").c_str(), 50, 1, -1);
      lnLHist_Mean1D[i]->GetXaxis()->SetTitle("-2LLH (Data, Pred)");
      lnLHist_Mean1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Mode1D[i] = new TH1D((name+"_ModelnL1D").c_str(),(name+"_ModelnL1D").c_str(), 50, 1, -1);
      lnLHist_Mode1D[i]->GetXaxis()->SetTitle("-2LLH (Data, Pred)");
      lnLHist_Mode1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Sample_DrawData[i] = new TH1D((name+"_lnLdrawdata").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_DrawData[i]->GetXaxis()->SetTitle("-2LLH (Data, Draw)");
      lnLHist_Sample_DrawData[i]->GetYaxis()->SetTitle("Counts");
      
      lnLHist_Sample_DrawflucDraw[i] = new TH1D((name+"_lnLdrawfluc").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_DrawflucDraw[i]->GetXaxis()->SetTitle("-2LLH (Draw Fluc, Draw)");
      lnLHist_Sample_DrawflucDraw[i]->GetYaxis()->SetTitle("Counts");
      
      lnLHist_Sample_PredflucDraw[i] = new TH1D((name+"_lnLpredfluc").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_PredflucDraw[i]->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Draw)");
      lnLHist_Sample_PredflucDraw[i]->GetYaxis()->SetTitle("Counts");
    }
  }
  //KS: Separate loop for thread safe reasons
  for (int i = 0; i < Length; ++i)
  {
    //KS: We copy histograms so delete original
    delete Nominal[i];
    delete NomW2[i];
  }
}

// *******************
// Add a throw from the MCMC to the posterior predictive
// The input here is nSamples long
void SampleSummary::AddThrow(std::vector<TH2Poly*> &SampleVector, std::vector<TH2Poly*> &W2Vec, const double LLHPenalty, const double Weight, const int DrawNumber) {
// *******************

  nThrows++;
  //KS: Only make sense for PosteriorPredictive
  if( !isPriorPredictive )RandomHist->Fill(DrawNumber);

  const int size = SampleVector.size();
  if (!CheckSamples(size)) throw;

  // Push back the throw
  MCVector.push_back(SampleVector);
  LLHPenaltyVector.push_back(LLHPenalty);
  WeightVector.push_back(Weight);
  W2MCVector.push_back(W2Vec);

  // Initialise the posterior hist
  if (first_pass)
  {
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      const int nXBins = 500;
      //Initialise TH1D which corresponds to each bin in the sample's th2poly
      std::string name = std::string(SampleVector[SampleNum]->GetName());
      for (int i = 1; i <= maxBins[SampleNum]; ++i)
      {
        //Get PolyBin
        TH2PolyBin* bin = (TH2PolyBin*) SampleVector[SampleNum]->GetBins()->At(i-1);

        // Just make a little fancy name
        std::stringstream ss2;
        ss2 << name << "_";
        ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
        ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";

        PosteriorHist[SampleNum][i] = new TH1D(ss2.str().c_str(), ss2.str().c_str(),nXBins, 1, -1);
        PosteriorHist[SampleNum][i]->SetDirectory(0);
        w2Hist[SampleNum][i] = new TH1D(("w2_"+ss2.str()).c_str(), ("w2_"+ss2.str()).c_str(),nXBins, 1, -1);
        w2Hist[SampleNum][i]->SetDirectory(0);
        if(DoBetaParam)
        {
          std::string betaName = "#beta_param_";
          BetaHist[SampleNum][i] = new TH1D((betaName+ss2.str()).c_str(), (betaName+ss2.str()).c_str(), 70, 1, -1);
          BetaHist[SampleNum][i]->GetXaxis()->SetTitle("#beta parameter value");
          BetaHist[SampleNum][i]->GetYaxis()->SetTitle("Counts");
        }
      } //end loop over bins
    }//end loop over samples
  }
  first_pass = false;

  // Loop over the sameples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
  {
    if (SampleVector[SampleNum] == NULL) continue;
    if(doShapeOnly) NormaliseTH2Poly(SampleVector[SampleNum]);
    // Loop over the distribution and fill the prior/posterior predictive
    for (int i = 1; i <= maxBins[SampleNum]; ++i) {
      const double Content = SampleVector[SampleNum]->GetBinContent(i);
      PosteriorHist[SampleNum][i]->Fill(Content, Weight);
      const double w2 = W2Vec[SampleNum]->GetBinContent(i);
      w2Hist[SampleNum][i]->Fill(w2, Weight);
      if(DoBetaParam)
      {
        const double data = DataHist[SampleNum]->GetBinContent(i);
        const double BetaParam = GetBetaParameter(data, Content, w2, likelihood);
        BetaHist[SampleNum][i]->Fill(BetaParam, Weight);
      }
    } // end bin loop
  } // end samples loop
} // end AddThrow

// *******************
// Add a throw from the MCMC to the posterior predictive
// The input here is has dimension [nsample][nMaCh3Modes]
void SampleSummary::AddThrowByMode(std::vector<std::vector<TH2Poly*>> &SampleVector_ByMode) {
// *******************

  MCVectorByMode.push_back(SampleVector_ByMode);

  //KS: This means this is first time
  if(!DoByModePlots)
  {
    std::cout<< "Turning reaction breadkwon mode, brum brum"<<std::endl;
    PosteriorHist_ByMode = new TH1D***[nSamples];
    MeanHist_ByMode = new TH2Poly**[nSamples];

    for (_int_ SampleNum = 0;  SampleNum < nSamples; SampleNum++)
    {
      if (DataHist[SampleNum] == NULL) continue;

      PosteriorHist_ByMode[SampleNum] = new TH1D**[Modes->GetNModes()+1];
      MeanHist_ByMode[SampleNum] = new TH2Poly*[Modes->GetNModes()+1];
      for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
      {
        PosteriorHist_ByMode[SampleNum][j] = new TH1D*[maxBins[SampleNum]+1];
        const int nXBins = 500;

        std::string name = std::string(NominalHist[SampleNum]->GetName());
        name = name.substr(0, name.find("_nom"));
        name = name + "_"+Modes->GetMaCh3ModeName(j);

        for (int i = 1; i <= maxBins[SampleNum]; i++)
        {
          //Get PolyBin
          TH2PolyBin* bin = (TH2PolyBin*)NominalHist[SampleNum]->GetBins()->At(i-1);

          // Just make a little fancy name
          std::stringstream ss2;
          ss2 << name << "_";
          ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
          ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";

          //Initialise TH1D which corresponds to each bin in the sample's th2poly
          PosteriorHist_ByMode[SampleNum][j][i] = new TH1D((name+ss2.str()).c_str(),(name+ss2.str()).c_str(),nXBins, 1, -1);
        }
        MeanHist_ByMode[SampleNum][j] = (TH2Poly*)(NominalHist[SampleNum]->Clone());
        MeanHist_ByMode[SampleNum][j]->SetNameTitle((name+"_mean").c_str(), (name+"_mean").c_str());
        MeanHist_ByMode[SampleNum][j]->Reset("");
        MeanHist_ByMode[SampleNum][j]->GetZaxis()->SetTitle("Mean");
      }
    }
  }
  DoByModePlots = true;    
  // Loop over the sameples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (_int_ SampleNum = 0;  SampleNum < nSamples; SampleNum++)
  {
    if (DataHist[SampleNum] == NULL) continue;
    
    for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
    {
      if(doShapeOnly) NormaliseTH2Poly(SampleVector_ByMode[SampleNum][j]);
      // Loop over the distribution and fill the prior/posterior predictive
      for (int i = 1; i <= maxBins[SampleNum]; ++i)
      {
        const double Content = SampleVector_ByMode[SampleNum][j]->GetBinContent(i);
        const int Entries = PosteriorHist_ByMode[SampleNum][j][i]->GetEntries();
        PosteriorHist_ByMode[SampleNum][j][i]->Fill(Content, WeightVector[Entries]);
      }
    }
  }
} // end AddThrowByMode

// **********************
void SampleSummary::PrepareOutput() {
// **********************

  // Make the output file (MakePosterioPredictive call writes to this)
  std::string TempString = OutputName;
  TempString.replace(TempString.find(".root"), 5, std::string("_procsW2.root"));
  Outputfile = new TFile(TempString.c_str(), "RECREATE");

  // The array of doubles we write to the TTree
  // Data vs Draw
  llh_data_draw = new double[nSamples];
  // Fluctuated Draw vs Draw
  llh_drawfluc_draw = new double[nSamples];
  // Fluctuated Predicitve vs Draw
  llh_predfluc_draw = new double[nSamples];

  // Data vs Draw using Rate
  llh_rate_data_draw = new double[nSamples];
  // Data vs Fluctuated Predictive using Rate
  llh_rate_predfluc_draw = new double[nSamples];

  // Data vs Fluctuated Draw
  llh_data_drawfluc = new double[nSamples];
  // Data vs Fluctuated Predictive
  llh_data_predfluc = new double[nSamples];
  // Draw vs Predictive
  llh_draw_pred = new double[nSamples];
  // Fluctuated Draw vs Predictive
  llh_drawfluc_pred = new double[nSamples];
  // Fluctuated Draw vs Fluctuated Predictive
  llh_drawfluc_predfluc = new double[nSamples];

  // Fluctuated Predictive vs Predictive
  llh_predfluc_pred = new double[nSamples];
  // Fluctuated Data vs Draw
  llh_datafluc_draw = new double[nSamples];
  
  // Data vs Draw for 1D projection
  llh_data_draw_ProjectX = new double[nSamples];
  llh_drawfluc_draw_ProjectX = new double[nSamples];
    
  // The output tree we're going to write to
  OutputTree = new TTree("LLH_draws", "LLH_draws");
  SampleNames.resize(nSamples);
  // Loop over the samples and set the addresses of the variables to write to file
  for (_int_ i = 0; i < nSamples; ++i)
  {
    // Get the name
    std::string SampleName = SamplePDF->GetSampleName(i);
    // Strip out spaces
    while (SampleName.find(" ") != std::string::npos) {
      SampleName.replace(SampleName.find(" "), 1, std::string("_"));
    }
    SampleNames[i] = SampleName;
    //CW: Also strip out - signs because it messes up TBranches
    while (SampleName.find("-") != std::string::npos) {
      SampleName.replace(SampleName.find("-"), 1, std::string("_"));
    }
//  All LLH below are used for actual p-value calculations
    OutputTree->Branch((SampleName+"_data_draw").c_str(),     &llh_data_draw[i]);
    OutputTree->Branch((SampleName+"_drawfluc_draw").c_str(), &llh_drawfluc_draw[i]);
    OutputTree->Branch((SampleName+"_predfluc_draw").c_str(), &llh_predfluc_draw[i]);

//  All LLH below are used for actual p-value calculations however using rate only
    OutputTree->Branch((SampleName+"_rate_data_draw").c_str(), &llh_rate_data_draw[i]);
    OutputTree->Branch((SampleName+"_rate_predfluc_draw").c_str(), &llh_rate_predfluc_draw[i]);

//  All LLH below are for validation reason but not used for final P-Value
    OutputTree->Branch((SampleName+"_data_drawfluc").c_str(), &llh_data_drawfluc[i]);
    OutputTree->Branch((SampleName+"_data_predfluc").c_str(), &llh_data_predfluc[i]);
    OutputTree->Branch((SampleName+"_draw_pred").c_str(),     &llh_draw_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_pred").c_str(), &llh_drawfluc_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_predfluc").c_str(), &llh_drawfluc_predfluc[i]);
    OutputTree->Branch((SampleName+"_predfluc_pred").c_str(), &llh_predfluc_pred[i]);
    OutputTree->Branch((SampleName+"_datafluc_draw").c_str(), &llh_datafluc_draw[i]);
 
//  All LLH below are used for calcauting P-Value but using 1D projections
    OutputTree->Branch((SampleName+"_data_draw_ProjectX").c_str(), &llh_data_draw_ProjectX[i]);
    OutputTree->Branch((SampleName+"_drawfluc_draw_ProjectX").c_str(), &llh_drawfluc_draw_ProjectX[i]);
  }
//All LLH below are used for actual p-value calculations
  OutputTree->Branch("LLH_Penalty",         &llh_penalty);
  OutputTree->Branch("Total_LLH_Data_Draw", &total_llh_data_draw);
  OutputTree->Branch("Total_LLH_DrawFluc_Draw", &total_llh_drawfluc_draw);
  OutputTree->Branch("Total_LLH_PredFluc_Draw", &total_llh_predfluc_draw);
  
//  All LLH below are used for actual p-value calculations however using rate only
  OutputTree->Branch("Total_LLH_Rate_PredFluc_Draw", &total_llh_rate_predfluc_draw);

//All LLH below are for validation reason but not used for final P-Value
  OutputTree->Branch("Total_LLH_Data_DrawFluc", &total_llh_data_drawfluc);
  OutputTree->Branch("Total_LLH_Data_PredFluc", &total_llh_data_predfluc);
  OutputTree->Branch("Total_LLH_Draw_Pred",     &total_llh_draw_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_Pred", &total_llh_drawfluc_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_PredFluc", &total_llh_drawfluc_predfluc);
  OutputTree->Branch("Total_LLH_PredFluc_Pred", &total_llh_predfluc_pred);
  OutputTree->Branch("Total_LLH_DataFluc_Draw", &total_llh_datafluc_draw);
  
//All LLH below are used for calcauting P-Value but 1D projections
  OutputTree->Branch("total_llh_data_draw_ProjectX", &total_llh_data_draw_ProjectX);
  OutputTree->Branch("total_llh_drawfluc_draw_ProjectX", &total_llh_drawfluc_draw_ProjectX);

  Outputfile->cd();
  Dir = new TDirectory*[nSamples];
  for (_int_ i = 0; i < nSamples; ++i)
  {
    // Make a new directory
    Dir[i] = Outputfile->mkdir((SampleNames[i]).c_str());
  }
}

// *******************
// Write the contents to the file
void SampleSummary::Write() {
// *******************

  // Prepare the output tree
  PrepareOutput();

  std::cout << "Summarising " << nThrows << " throws..." << std::endl;
  // After all the throws are added finalise the sample
  TStopwatch timer;
  timer.Start();
  MakePredictive();
  timer.Stop();
  std::cout << "Made Prior/Posterior Predictive, it took " << timer.RealTime() << "s" <<", now writing..." << std::endl;

  // Study Bayesian Information Criterion
  StudyBIC();

  // Study Deviance Information Criterion
  StudyDIC();

  OutputTree->Write();

  // Make the various distributions
  lnLHist->Write();
  lnLHist_drawfluc->Write();
  lnLHist_drawflucdraw->Write();
  lnLHist_drawdata->Write();
  lnLDrawHist->Write();
  lnLFlucHist->Write();
  lnLDrawHistRate->Write();
  //KS: Only avaible for Posterior Predictive
  if(!isPriorPredictive) RandomHist->Write();

  lnLFlucHist_ProjectX->Write();
  
  // Loop over each sample and write to file
  //KS: Multithreading is tempting here but we also write to ROOT file, separating all LLH and poly projections from write could work well
  for (_int_ i = 0; i < nSamples; ++i)
  {
     // Skip the null histograms
    if (DataHist[i] == NULL || NoOverflowIntegral(DataHist[i]) == 0) continue;
    Dir[i]->cd();

    // Make the data/MC ratio histogram
    TH2Poly *RatioHistMean = RatioPolys(DataHist[i], MeanHist[i]);
    RatioHistMean->GetZaxis()->SetTitle("Data/Mean");
    TH2Poly *RatioHistMode = RatioPolys(DataHist[i], ModeHist[i]);
    RatioHistMode->GetZaxis()->SetTitle("Data/Mode");
    TH2Poly *RatioHistNom = RatioPolys(DataHist[i], NominalHist[i]);
    RatioHistNom->GetZaxis()->SetTitle("Data/Nom");

    // And the normalised data histogram
    TH2Poly *DataNormHist = NormalisePoly(DataHist[i]);
    // Last true refers to if project along x or y
    TH2Poly *MeanNormHist = NormalisePoly(MeanHist[i]);
    TH2Poly *ModeNormHist = NormalisePoly(ModeHist[i]);
    TH1D *MeanProjectX = ProjectPoly(MeanHist[i], true, i, true);
    TH1D *MeanProjectY = ProjectPoly(MeanHist[i], false, i, true);
    TH1D *ModeProjectX = ProjectPoly(ModeHist[i], true, i, true);
    TH1D *ModeProjectY = ProjectPoly(ModeHist[i], false, i, true);

    TH1D *MeanHistCorrectedProjectX = NULL;
    if(DoBetaParam) MeanHistCorrectedProjectX = ProjectPoly(MeanHistCorrected[i], true, i, true);
    TH1D *MeanHistCorrectedProjectY = NULL;
    if(DoBetaParam) MeanHistCorrectedProjectY = ProjectPoly(MeanHistCorrected[i], false, i, true);

    TH1D *W2MeanProjectX = ProjectPoly(W2MeanHist[i], true, i);
    TH1D *W2MeanProjectY = ProjectPoly(W2MeanHist[i], false, i);
    TH1D *W2ModeProjectX = ProjectPoly(W2ModeHist[i], true, i);
    TH1D *W2ModeProjectY = ProjectPoly(W2ModeHist[i], false, i);

    TH2Poly *NomNormHist = NormalisePoly(NominalHist[i]);
    TH1D *NomProjectX = ProjectPoly(NominalHist[i], true, i);
    TH1D *NomProjectY = ProjectPoly(NominalHist[i], false, i);

    TH1D *W2NomProjectX = ProjectPoly(W2NomHist[i], true, i);
    TH1D *W2NomProjectY = ProjectPoly(W2NomHist[i], false, i);

    // Same for the TH2Ds
    CalcLLH(DataHist[i], NominalHist[i], W2NomHist[i]);
    CalcLLH(DataHist[i], MeanHist[i], W2MeanHist[i]);
    CalcLLH(DataHist[i], ModeHist[i], W2ModeHist[i]);

    // Calculate the log likelihood for the 1D dists
    // Sets the title of the second TH1D to the -2LLH
    CalcLLH(DataHist_ProjectX[i], NomProjectX, W2NomProjectX);
    CalcLLH(DataHist_ProjectX[i], MeanProjectX, W2MeanProjectX);
    CalcLLH(DataHist_ProjectX[i], ModeProjectX, W2ModeProjectX);
    CalcLLH(DataHist_ProjectY[i], NomProjectY, W2NomProjectY);
    CalcLLH(DataHist_ProjectY[i], MeanProjectY, W2MeanProjectY);
    CalcLLH(DataHist_ProjectY[i], ModeProjectY, W2ModeProjectY);

    std::string SampleName = SampleNames[i];
    // Also strip out - signs because it messes up TBranches
    while (SampleName.find("-") != std::string::npos) {
      SampleName.replace(SampleName.find("-"), 1, std::string("_"));
    } 
    OutputTree->Draw((SampleName+"_data_draw:"+SampleName+"_drawfluc_draw>>htemp").c_str());
    TH2D *TempHistogram = (TH2D*)((gDirectory->Get("htemp"))->Clone());
    TempHistogram->GetXaxis()->SetTitle("-2LLH(Draw Fluc, Draw)");
    TempHistogram->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram->SetNameTitle((SampleNames[i]+"_drawfluc_draw").c_str(), (SampleNames[i]+"_drawfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram);
    TempHistogram->Write();
    delete TempHistogram;

    // Also write the 2D histograms for the p-value
    OutputTree->Draw((SampleName+"_data_draw:"+SampleName+"_predfluc_draw>>htemp2").c_str());
    TH2D *TempHistogram2 = (TH2D*)((gDirectory->Get("htemp2"))->Clone());
    TempHistogram2->GetXaxis()->SetTitle("-2LLH(Pred Fluc, Draw)");
    TempHistogram2->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram2->SetNameTitle((SampleNames[i]+"_predfluc_draw").c_str(), (SampleNames[i]+"_predfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram2);
    TempHistogram2->Write();
    delete TempHistogram2;
   
    // finally p-value for 1D projection
    OutputTree->Draw((SampleName+"_rate_data_draw:"+SampleName+"_rate_predfluc_draw>>htemp3").c_str());
    TH2D *TempHistogram3 = (TH2D*)((gDirectory->Get("htemp3"))->Clone());
    TempHistogram3->GetXaxis()->SetTitle("-2LLH(Pred Fluc, Draw)");
    TempHistogram3->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram3->SetNameTitle((SampleNames[i]+"_rate_predfluc_draw").c_str(), (SampleNames[i]+"_rate_predfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram3);
    TempHistogram3->Write();
    delete TempHistogram3;

    // finally p-value for 1D projection
    OutputTree->Draw((SampleName+"_data_draw_ProjectX:"+SampleName+"_drawfluc_draw_ProjectX>>htemp4").c_str());
    TH2D *TempHistogram4 = (TH2D*)((gDirectory->Get("htemp4"))->Clone());
    TempHistogram4->GetXaxis()->SetTitle(("-2LLH_{Draw Fluc, Draw} for " + SamplePDF->GetKinVarLabel(i, 0)).c_str());
    TempHistogram4->GetYaxis()->SetTitle(("-2LLH_{Data, Draw} for " + SamplePDF->GetKinVarLabel(i, 0)).c_str());
    TempHistogram4->SetNameTitle((SampleNames[i]+"_drawfluc_draw_ProjectX").c_str(), (SampleNames[i]+"_drawfluc_draw_ProjectX").c_str());
    MakeCutLLH2D(TempHistogram4);
    TempHistogram4->Write();
    delete TempHistogram4;
    
    // Write the Histograms to each folder
    DataHist[i]->Write();
    NominalHist[i]->Write();
    MeanHist[i]->Write();
    ModeHist[i]->Write();
    RatioHistMean->Write();
    RatioHistMode->Write();
    RatioHistNom->Write();
    if(DoBetaParam) MeanHistCorrected[i]->Write();

    W2NomHist[i]->Write();
    W2MeanHist[i]->Write();
    W2ModeHist[i]->Write();

    DataNormHist->Write();
    NomNormHist->Write();
    MeanNormHist->Write();
    ModeNormHist->Write();

    DataHist_ProjectX[i]->Write();
    NomProjectX->Write();
    MeanProjectX->Write();
    ModeProjectX->Write();
    if(DoBetaParam) MeanHistCorrectedProjectX->Write();
    ViolinHists_ProjectX[i]->Write();
    
    DataHist_ProjectY[i]->Write();
    NomProjectY->Write();
    MeanProjectY->Write();
    ModeProjectY->Write();
    if(DoBetaParam) MeanHistCorrectedProjectY->Write();
    ViolinHists_ProjectY[i]->Write();

    W2NomProjectX->Write();
    W2MeanProjectX->Write();
    W2ModeProjectX->Write();

    W2NomProjectY->Write();
    W2MeanProjectY->Write();
    W2ModeProjectY->Write();
    
    //KS: This will dump lots of hists, use it only for debugging
    if(Debug > 0)
    {
      for (int b = 1; b <= maxBins[i]; ++b)
      {
        PosteriorHist[i][b]->Write();
        //This isn't useful check only in desperation
        if(Debug > 1) w2Hist[i][b]->Write();
      }
    }
    lnLHist_Mean[i]->Write();
    lnLHist_Mode[i]->Write();

    lnLHist_Mean_ProjectX[i]->Write();

    lnLHist_Mean1D[i]->Write();
    lnLHist_Mode1D[i]->Write();

    MakeCutLLH1D(lnLHist_Sample_DrawData[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_DrawData[i]->Write();
    MakeCutLLH1D(lnLHist_Sample_DrawflucDraw[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_DrawflucDraw[i]->Write();
    MakeCutLLH1D(lnLHist_Sample_PredflucDraw[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_PredflucDraw[i]->Write();
    
    if(DoByModePlots)
    {
      for (_int_ j = 0; j < Modes->GetNModes()+1; ++j)
      {
        MeanHist_ByMode[i][j]->Write();
        TH1D *MeanProjectX_ByMode = ProjectPoly(MeanHist_ByMode[i][j], true, i, true);
        TH1D *MeanProjectY_ByMode = ProjectPoly(MeanHist_ByMode[i][j], false, i, true);
        MeanProjectX_ByMode->Write();
        MeanProjectY_ByMode->Write();
        //KS: This will dump lots of hists, use it only for debugging
        if(Debug > 0)
        {
          for (int b = 1; b <= maxBins[i]; ++b)
          {
            PosteriorHist_ByMode[i][j][b]->Write();
          }
        }
        delete MeanProjectX_ByMode;
        delete MeanProjectY_ByMode;
      } // End loop over bins
    }
    
    // Delete temporary objects
    delete RatioHistMean;
    delete RatioHistMode;
    delete RatioHistNom;

    delete DataNormHist;
    delete MeanNormHist;
    delete ModeNormHist;
    delete NomNormHist;

    delete DataHist_ProjectX[i];
    delete MeanProjectX;
    delete ModeProjectX;
    if(DoBetaParam) delete MeanHistCorrectedProjectX;
    delete NomProjectX;

    delete DataHist_ProjectY[i];
    delete MeanProjectY;
    delete ModeProjectY;
    if(DoBetaParam) delete MeanHistCorrectedProjectY;
    delete NomProjectY;

    delete W2NomProjectX;
    delete W2MeanProjectX;
    delete W2ModeProjectX;
  
    delete W2NomProjectY;
    delete W2MeanProjectY;
    delete W2ModeProjectY;
    
    std::cout << std::endl;
  } //end loop over samples
  delete[] DataHist_ProjectX;
  delete[] DataHist_ProjectY;

  if(DoBetaParam) PlotBetaParameters();

  StudyKinematicCorrelations();
  std::cout << "Wrote to " << Outputfile->GetName() << std::endl;
}

// *******************
// Make the posterior predictive distributions: fit Poissons etc
void SampleSummary::MakePredictive() {
// *******************

  // First make the projection on the z axis of the TH3D* for every pmu cosmu bin
  double llh_total_temp = 0.0;

  // Loop over the samples
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:llh_total_temp)
  #endif
  for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
  {
    // Skip disabled samples
    if (DataHist[SampleNum] == NULL || NoOverflowIntegral(DataHist[SampleNum]) == 0) continue;

    // Count the -2LLH for each histogram
    double negLogL_Mean = 0.0;
    double negLogL_Mode = 0.0;

    // Loop over each pmu cosmu bin
    for (int j = 1; j < maxBins[SampleNum]+1; ++j)
    {
      TH1D *Projection = (TH1D*) PosteriorHist[SampleNum][j];
      TH1D *W2Projection = (TH1D*) w2Hist[SampleNum][j];

      // Data content for the j,kth bin
      const double nData = DataHist[SampleNum]->GetBinContent(j);
      
      // Get the mean for this projection for all the samples
      // This is the mean prediction for this given j,k bin
      const double nMean = Projection->GetMean();
      const double nMeanError = Projection->GetRMS();
      const double nMode = Projection->GetBinCenter(Projection->GetMaximumBin());
      const double nModeError = GetModeError(Projection);
      
      const double nW2Mean = W2Projection->GetMean();
      const double nW2Mode = W2Projection->GetBinCenter(W2Projection->GetMaximumBin());
      
      double TempLLH_Mean = 0.0;
      double TempLLH_Mode = 0.0;
      
      //KS:Get LLH contribution getTestStatLLH can calculate Barlow Beeston/IceCube or Poisson
      TempLLH_Mean = SamplePDF->getTestStatLLH(nData, nMean, nW2Mean);
      TempLLH_Mode = SamplePDF->getTestStatLLH(nData, nMode, nW2Mode);

      // Increment -2LLH
      //KS: do times 2 because banff reports chi2
      negLogL_Mean += 2*TempLLH_Mean;
      negLogL_Mode += 2*TempLLH_Mode;
      
      // Set the content and error to the mean in the bin
      MeanHist[SampleNum]->SetBinContent(j, MeanHist[SampleNum]->GetBinContent(j)+nMean);
      // KS: This -1 is only needed for root older than 6.18 for more see https://t2k-experiment.slack.com/archives/CU9CBG6NS/p1714551365661589
      MeanHist[SampleNum]->SetBinError(j, nMeanError);

      if(DoBetaParam)
      {
        TH1D *BetaTemp = (TH1D*)BetaHist[SampleNum][j];
        const double nBetaMean = BetaTemp->GetMean();
        const double nBetaMeanError = BetaTemp->GetRMS();
        //KS: Here we modify predictions by beta parameter from Barlow-Beeston
        MeanHistCorrected[SampleNum]->SetBinContent(j, MeanHistCorrected[SampleNum]->GetBinContent(j)+nMean*nBetaMean);
        //KS: Use error propagation to calcuate error
        const double ErrorTemp = std::sqrt( (nBetaMean*nMeanError) * (nBetaMean*nMeanError) + (nMean*nBetaMeanError) * (nMean*nBetaMeanError));
        // KS: This -1 is only needed for root older than 6.18 for more see https://t2k-experiment.slack.com/archives/CU9CBG6NS/p1714551365661589
        MeanHistCorrected[SampleNum]->SetBinError(j, ErrorTemp);
      }
      // Set the content to the mode in the bin
      ModeHist[SampleNum]->SetBinContent(j, ModeHist[SampleNum]->GetBinContent(j)+nMode);
      // KS: This -1 is only needed for root older than 6.18 for more see https://t2k-experiment.slack.com/archives/CU9CBG6NS/p1714551365661589
      ModeHist[SampleNum]->SetBinError(j, nModeError);
      // Set the content to the mean in the bin
      W2MeanHist[SampleNum]->SetBinContent(j, W2MeanHist[SampleNum]->GetBinContent(j)+nW2Mean);
      // Set the content to the mode in the bin
      W2ModeHist[SampleNum]->SetBinContent(j, W2ModeHist[SampleNum]->GetBinContent(j)+nW2Mode);
      
      // Set the mean and average LLH for this given bin
      // Can use these hists to see where the largest -2LLH hists come from
      lnLHist_Mean[SampleNum]->SetBinContent(j, 2.0*TempLLH_Mean);
      lnLHist_Mode[SampleNum]->SetBinContent(j, 2.0*TempLLH_Mode);
      
      lnLHist_Mean1D[SampleNum]->Fill(2.0*TempLLH_Mean);
      lnLHist_Mode1D[SampleNum]->Fill(2.0*TempLLH_Mode);
    } // End loop over bins
    if(DoByModePlots)
    {
      for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
      {
        // Loop over each pmu cosmu bin
        for (int i = 1; i < maxBins[SampleNum]+1; ++i)
        {
          // Make the posterior/prior predictive projection on z
          // The z axis of Predictive is the bin content
          // Essentially zooming in on one bin and looking at the mean and mode of that bin
          TH1D *Projection = (TH1D*)PosteriorHist_ByMode[SampleNum][j][i];

          // Get the mean for this projection for all the samples
          const double nMean = Projection->GetMean();
          const double nMeanError = Projection->GetRMS();

          // Set the content and error to the mean in the bin
          MeanHist_ByMode[SampleNum][j]->SetBinContent(i, MeanHist_ByMode[SampleNum][j]->GetBinContent(i)+nMean);
          // KS: This -1 is only needed for root older than 6.18 for more see https://t2k-experiment.slack.com/archives/CU9CBG6NS/p1714551365661589
          MeanHist_ByMode[SampleNum][j]->SetBinError(i, nMeanError);
        } // End loop over bins
      }
    }
    llh_total_temp += negLogL_Mean;
  } // End loop over samples

  // This is not multithreaded as due to ProjectPoly it is not safe
  for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
  {
    // Skip disabled samples
    if (DataHist[SampleNum] == NULL || NoOverflowIntegral(DataHist[SampleNum]) == 0) continue;

    //KS:: Might consider caching it as we use it once agian much later
    TH1D *MeanProjectX = ProjectPoly(MeanHist[SampleNum], true, SampleNum, true);
    TH1D *W2MeanProjectX = ProjectPoly(W2MeanHist[SampleNum], true, SampleNum);
    // Loop over each pmu bin for 1D projection
    for (int j = 1; j <= lnLHist_Mean_ProjectX[SampleNum]->GetXaxis()->GetNbins(); ++j)
    {
      // Data content for the j,kth bin
      const double nData = DataHist_ProjectX[SampleNum]->GetBinContent(j);
      const double nMean = MeanProjectX->GetBinContent(j);
      const double nW2Mean = W2MeanProjectX->GetBinContent(j);

      double TempLLH_Mean = 0.0;
      TempLLH_Mean = SamplePDF->getTestStatLLH(nData, nMean, nW2Mean);

      //KS: do times 2 because banff reports chi2
      lnLHist_Mean_ProjectX[SampleNum]->SetBinContent(j, 2.0*TempLLH_Mean);
    }// End loop over  bins

    delete MeanProjectX;
    delete W2MeanProjectX;
  } // End loop over samples

  llh_total = llh_total_temp;
  // Now we have our posterior predictive histogram and it's LLH
  std::cout << "Prior/Posterior predictive LLH mean (sample only) = " << llh_total << std::endl;
  std::stringstream ss;
  ss << llh_total;
  lnLHist->SetTitle((std::string(lnLHist->GetTitle())+"_"+ss.str()).c_str());

  // Now make the fluctuated hists of the MeanHist and ModeHist
  MakeChi2Hists();

  // Get the 1D LLH dists
  MakeCutLLH();

} // End MakePredictive() function

// *******************
// Make the fluctuated histograms (2D and 1D) for the chi2s
// Essentially taking the MCMC draws and calculating their LLH to the Posterior predictive distribution
// And additionally taking the data histogram and calculating the LLH to the predictive distribution
// Additionally we calculate the chi2 of the draws (fluctuated) of  the MC with the prior/posterior predictive and plot it vs the chi2 from the draws of MCMC and the data
void SampleSummary::MakeChi2Hists() {
// *******************

  std::cout << "Making the chi2 histograms..." << std::endl;
  // Have this to signify if we're doing the first pass
  first_pass = true;

  double AveragePenalty = 0;

  // Vectors to hold exact LLH
  std::vector<double> LLH_PredFluc_V(nThrows);
  std::vector<double> LLH_DataDraw_V(nThrows);
  std::vector<double> LLH_DrawFlucDraw_V(nThrows);
    
  // Loop over the draws
  // Should look into multi-threading this. Would require temporary THxx structures from arrays
  //KS: Update above would be ideal but currently we loop over samples (see loop below) which isn't as efficient as loop over throws but much much easier to implement
  for (unsigned int i = 0; i < nThrows; ++i)
  {
    if (i % (nThrows/10) == 0) {
      std::cout << "   On throw " << i << "/" << nThrows << " (" << int(double(i)*100.0/double(nThrows)) << "%)"<< std::endl;
    }

    // Set the total LLH to zero to initialise
    double total_llh_data_draw_temp = 0.0;
    double total_llh_drawfluc_draw_temp = 0.0;
    double total_llh_predfluc_draw_temp = 0.0;

    double total_llh_rate_data_draw_temp = 0.0;
    double total_llh_rate_predfluc_draw_temp = 0.0;

    double total_llh_data_drawfluc_temp = 0.0;
    double total_llh_data_predfluc_temp = 0.0;
    double total_llh_draw_pred_temp = 0.0;
    double total_llh_drawfluc_pred_temp = 0.0;
    double total_llh_drawfluc_predfluc_temp = 0.0;
    double total_llh_predfluc_pred_temp = 0.0;
    double total_llh_datafluc_draw_temp = 0.0;

    double total_llh_data_draw_ProjectX_temp = 0.0;
    double total_llh_drawfluc_draw_ProjectX_temp = 0.0;

    // Save the double that gets written to file
    llh_penalty = LLHPenaltyVector[i];
    AveragePenalty += llh_penalty;

    // Make the Poisson fluctuated hist
    TH2Poly **FluctHist = new TH2Poly*[nSamples];
    // Also Poisson fluctuate the drawn MCMC hist
    TH2Poly **FluctDrawHist = new TH2Poly*[nSamples];
    // Finally Poisson fluctuate the data histogram
    TH2Poly **DataFlucHist = new TH2Poly*[nSamples];

    // Finally Poisson fluctuate the data histogram
    TH1D **FluctDrawHistProjectX = new TH1D*[nSamples];
    TH1D **DrawHistProjectX = new TH1D*[nSamples];
    TH1D **DrawHistProjectY = new TH1D*[nSamples];
    TH1D **DrawW2HistProjectX = new TH1D*[nSamples];

    //KS: We have to clone histograms here to avoid cloning in MP loop, we have to make sure binning matches, content doesn't have to
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      FluctHist[SampleNum] = (TH2Poly*)(MeanHist[SampleNum]->Clone());
      FluctDrawHist[SampleNum] = (TH2Poly*)(MeanHist[SampleNum]->Clone());
      DataFlucHist[SampleNum] = (TH2Poly*)(MeanHist[SampleNum]->Clone());

      FluctDrawHistProjectX[SampleNum] = (TH1D*)(DataHist_ProjectX[SampleNum]->Clone());

      // Get the ith draw for the jth sample
      TH2Poly *DrawHist = (TH2Poly*)(MCVector[i][SampleNum]);
      TH2Poly *DrawW2Hist = (TH2Poly*)(W2MCVector[i][SampleNum]);

      //ProjectPoly calls new TH1D under the hood, never define new ROOT object under MP...
      DrawHistProjectX[SampleNum] = ProjectPoly(DrawHist, true, SampleNum);
      DrawW2HistProjectX[SampleNum] = ProjectPoly(DrawW2Hist, true, SampleNum);
      DrawHistProjectY[SampleNum] = ProjectPoly(DrawHist, false, SampleNum);
    }
    #ifdef MULTITHREAD 
    //KS: might be most obscure OMP reduction I have ever made...
    #pragma omp parallel for reduction(+:total_llh_data_draw_temp, total_llh_drawfluc_draw_temp, total_llh_predfluc_draw_temp, total_llh_rate_data_draw_temp, total_llh_rate_predfluc_draw_temp, total_llh_data_drawfluc_temp, total_llh_data_predfluc_temp, total_llh_draw_pred_temp, total_llh_drawfluc_pred_temp, total_llh_drawfluc_predfluc_temp, total_llh_predfluc_pred_temp, total_llh_datafluc_draw_temp, total_llh_data_draw_ProjectX_temp, total_llh_drawfluc_draw_ProjectX_temp)
    #endif
    // Loop over the samples
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      // Get the ith draw for the jth sample
      TH2Poly *DrawHist = (TH2Poly*)(MCVector[i][SampleNum]);
      TH2Poly *DrawW2Hist = (TH2Poly*)(W2MCVector[i][SampleNum]);
      // Skip empty samples
      if (DrawHist == NULL) continue;

      // Add LLH penalties from the systematics to the LLH that use the drawn histogram
      // Data vs Draw
      llh_data_draw[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Draw
      llh_drawfluc_draw[SampleNum] = llh_penalty;
      // Fluctuated Predicitve vs Draw
      llh_predfluc_draw[SampleNum] = llh_penalty;
      
      // Data vs Draw using rate
      llh_rate_data_draw[SampleNum] = llh_penalty;
      // Fluctuated Predicitve vs Draw using rate
      llh_rate_predfluc_draw[SampleNum] = llh_penalty;

       // Data vs Fluctuated Draw
      llh_data_drawfluc[SampleNum] = llh_penalty;
      // Draw vs Predictive
      llh_draw_pred[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Predictive
      llh_drawfluc_pred[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Fluctuated Predictive
      llh_drawfluc_predfluc[SampleNum] = llh_penalty;
      // Fluctuated Data vs Draw
      llh_datafluc_draw[SampleNum] = llh_penalty;
      
      //Some LLH for 1D projections
      llh_data_draw_ProjectX[SampleNum] = llh_penalty;
      llh_drawfluc_draw_ProjectX[SampleNum] = llh_penalty;
      
      //Other get 0 penalty term
      // Fluctuated Predictive vs Predictive
      llh_predfluc_pred[SampleNum] = 0.0;
      // Data vs Fluctuated Predictive
      llh_data_predfluc[SampleNum] = 0.0;

      // Make the Poisson fluctuated hist
      MakeFluctuatedHistogram(FluctHist[SampleNum], MeanHist[SampleNum]);
      // Also Poisson fluctuate the drawn MCMC hist
      MakeFluctuatedHistogram(FluctDrawHist[SampleNum], DrawHist);
      // Finally Poisson fluctuate the data histogram
      MakeFluctuatedHistogram(DataFlucHist[SampleNum], DataHist[SampleNum]);

      // Likelihood between the drawn histogram and the data
      const double DataDrawLLH = GetLLH(DataHist[SampleNum], DrawHist, DrawW2Hist);
      llh_data_draw[SampleNum] += DataDrawLLH;
      total_llh_data_draw_temp += DataDrawLLH;
      
      // Likelihood between drawn histogram and fluctuated drawn histogram
      const double DrawFlucDrawLLH = GetLLH(FluctDrawHist[SampleNum], DrawHist, DrawW2Hist);
      llh_drawfluc_draw[SampleNum] += DrawFlucDrawLLH;
      total_llh_drawfluc_draw_temp += DrawFlucDrawLLH;
      
      // Likelihood between drawn histogram and fluctuated posterior predictive distribution
      const double PredFlucDrawLLH = GetLLH(FluctHist[SampleNum], DrawHist, DrawW2Hist);
      llh_predfluc_draw[SampleNum] += PredFlucDrawLLH;
      total_llh_predfluc_draw_temp += PredFlucDrawLLH;

//Rate Based p-value
      // Likelihood between the drawn histogram and the data
      const double RateDataDrawLLH = SamplePDF->getTestStatLLH(NoOverflowIntegral(DataHist[SampleNum]), NoOverflowIntegral(DrawHist), NoOverflowIntegral(DrawW2Hist));
      llh_rate_data_draw[SampleNum] += RateDataDrawLLH;
      total_llh_rate_data_draw_temp += RateDataDrawLLH;

      // Likelihood between drawn histogram and fluctuated posterior predictive distribution using rate
      const double RatePredFlucDrawLLH = SamplePDF->getTestStatLLH(NoOverflowIntegral(FluctHist[SampleNum]), NoOverflowIntegral(DrawHist), NoOverflowIntegral(DrawW2Hist));
      llh_rate_predfluc_draw[SampleNum] += RatePredFlucDrawLLH;
      total_llh_rate_predfluc_draw_temp += RatePredFlucDrawLLH;

//    All LLH below are for validation reason but not used for final P-Value
      // Likelihood between the fluctuated drawn histogram and the data
      const double DataDrawFlucLLH = GetLLH(DataHist[SampleNum], FluctDrawHist[SampleNum], DrawW2Hist);
      llh_data_drawfluc[SampleNum] += DataDrawFlucLLH;
      total_llh_data_drawfluc_temp += DataDrawFlucLLH;

      // Likelihood between the drawn histogram and the data
      const double DataPredFlucLLH = GetLLH(DataHist[SampleNum], FluctHist[SampleNum], W2MeanHist[SampleNum]);
      llh_data_predfluc[SampleNum] += DataPredFlucLLH;
      total_llh_data_predfluc_temp += DataPredFlucLLH;

      // Likelihood between the drawn hist and the Posterior Predictive
      const double DrawPredLLH = GetLLH(DrawHist, MeanHist[SampleNum], W2MeanHist[SampleNum]);
      llh_draw_pred[SampleNum] += DrawPredLLH;
      total_llh_draw_pred_temp += DrawPredLLH;

      // Likelihood between fluctuated drawn and predictive
      const double DrawFlucPredLLH = GetLLH(FluctDrawHist[SampleNum], MeanHist[SampleNum], W2MeanHist[SampleNum]);
      llh_drawfluc_pred[SampleNum]  += DrawFlucPredLLH;
      total_llh_drawfluc_pred_temp  += DrawFlucPredLLH;

      // Likelihood between drawn histogram and fluctuated drawn histogram
      const double DrawFlucPredFlucLLH = GetLLH(FluctDrawHist[SampleNum], FluctHist[SampleNum], W2MeanHist[SampleNum]);
      llh_drawfluc_predfluc[SampleNum] += DrawFlucPredFlucLLH;
      total_llh_drawfluc_predfluc_temp += DrawFlucPredFlucLLH;

      // Likelihood between the fluctuated drawn histogram and the posterior predictive
      const double PredFlucPredLLH = GetLLH(FluctHist[SampleNum], MeanHist[SampleNum], W2MeanHist[SampleNum]);
      llh_predfluc_pred[SampleNum] += PredFlucPredLLH;
      total_llh_predfluc_pred_temp += PredFlucPredLLH;

      // Likelihood between fluctuated data histogram and drawn histogram 
      const double DataFlucDrawLLH = GetLLH(DataFlucHist[SampleNum], DrawHist, DrawW2Hist);
      llh_datafluc_draw[SampleNum] += DataFlucDrawLLH;
      total_llh_datafluc_draw_temp += DataFlucDrawLLH;

      lnLHist_Sample_DrawData[SampleNum]->Fill(DataDrawLLH);
      lnLHist_Sample_DrawflucDraw[SampleNum]->Fill(DrawFlucDrawLLH);
      lnLHist_Sample_PredflucDraw[SampleNum]->Fill(PredFlucDrawLLH);
      
//    At the end we leave LLH for 1D projections
      MakeFluctuatedHistogram(FluctDrawHistProjectX[SampleNum], DrawHistProjectX[SampleNum]);
      
      // Likelihood between the drawn histogram and the data for muon momentum
      const double DataDrawLLH_ProjectX = GetLLH(DataHist_ProjectX[SampleNum], DrawHistProjectX[SampleNum], DrawW2HistProjectX[SampleNum]);
      llh_data_draw_ProjectX[SampleNum] += DataDrawLLH_ProjectX;
      total_llh_data_draw_ProjectX_temp += DataDrawLLH_ProjectX;
      
      const double DrawFlucDrawLLH_ProjectX = GetLLH(FluctDrawHistProjectX[SampleNum], DrawHistProjectX[SampleNum], DrawW2HistProjectX[SampleNum]);
      llh_drawfluc_draw_ProjectX[SampleNum] += DrawFlucDrawLLH_ProjectX;
      total_llh_drawfluc_draw_ProjectX_temp += DrawFlucDrawLLH_ProjectX;
      
      //KS: This might seem complicated but we make X and Y projection for each sample. Then we add this to the Violin plot making nice Gaussian in each kineatmical bin of x and y axis
      FastViolinFill(ViolinHists_ProjectX[SampleNum], DrawHistProjectX[SampleNum]);
      FastViolinFill(ViolinHists_ProjectY[SampleNum], DrawHistProjectY[SampleNum]);
    } // End loop over samples (still looping throws)

    // Delete the temporary histograms
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      delete FluctHist[SampleNum];
      delete FluctDrawHist[SampleNum];
      delete DataFlucHist[SampleNum];
      delete FluctDrawHistProjectX[SampleNum];
      delete DrawHistProjectX[SampleNum];
      delete DrawHistProjectY[SampleNum];
      delete DrawW2HistProjectX[SampleNum];
    }
    delete[] FluctHist;
    delete[] FluctDrawHist;
    delete[] DataFlucHist;
    delete[] FluctDrawHistProjectX;
    delete[] DrawHistProjectX;
    delete[] DrawHistProjectY;
    delete[] DrawW2HistProjectX;

    total_llh_data_draw = total_llh_data_draw_temp;
    total_llh_drawfluc_draw = total_llh_drawfluc_draw_temp;
    total_llh_predfluc_draw = total_llh_predfluc_draw_temp;

    total_llh_rate_data_draw = total_llh_rate_data_draw_temp;
    total_llh_rate_predfluc_draw = total_llh_rate_predfluc_draw_temp;

    total_llh_data_drawfluc = total_llh_data_drawfluc_temp;
    total_llh_data_predfluc = total_llh_data_predfluc_temp;
    total_llh_draw_pred = total_llh_draw_pred_temp;
    total_llh_drawfluc_pred = total_llh_drawfluc_pred_temp;
    total_llh_drawfluc_predfluc = total_llh_drawfluc_predfluc_temp;
    total_llh_predfluc_pred = total_llh_predfluc_pred_temp;
    total_llh_datafluc_draw = total_llh_datafluc_draw_temp;

    total_llh_data_draw_ProjectX = total_llh_data_draw_ProjectX_temp;
    total_llh_drawfluc_draw_ProjectX = total_llh_drawfluc_draw_ProjectX_temp;

    // Add LLH penalties from the systematics to the LLH that use the drawn histogram
    total_llh_data_draw     += llh_penalty;
    total_llh_drawfluc_draw += llh_penalty;
    total_llh_predfluc_draw += llh_penalty;
    //Rate based
    total_llh_rate_data_draw += llh_penalty;
    total_llh_rate_predfluc_draw += llh_penalty;

    total_llh_data_drawfluc += llh_penalty;
    total_llh_draw_pred     += llh_penalty;
    total_llh_drawfluc_pred += llh_penalty;
    total_llh_drawfluc_predfluc += llh_penalty;

    total_llh_data_draw_ProjectX += llh_penalty;
    total_llh_drawfluc_draw_ProjectX += llh_penalty;

    lnLHist->Fill(total_llh_predfluc_pred);
    lnLHist_drawdata->Fill(total_llh_data_draw);
    lnLHist_drawfluc->Fill(total_llh_predfluc_draw);
    lnLHist_drawflucdraw->Fill(total_llh_drawfluc_draw);

    lnLDrawHist->Fill(total_llh_predfluc_draw, total_llh_data_draw);
    lnLFlucHist->Fill(total_llh_drawfluc_draw, total_llh_data_draw);

    lnLDrawHistRate->Fill(total_llh_rate_predfluc_draw, total_llh_rate_data_draw);

    lnLFlucHist_ProjectX->Fill(total_llh_drawfluc_draw_ProjectX, total_llh_data_draw_ProjectX);
    
    // Also save to arrays to make sure we have the utmost super accuracy
    LLH_PredFluc_V[i] = total_llh_predfluc_draw;
    LLH_DataDraw_V[i] = total_llh_data_draw;
    LLH_DrawFlucDraw_V[i] = total_llh_drawfluc_draw;

    // Write to the output tree
    OutputTree->Fill();
  } // End loop over throws

  AveragePenalty = AveragePenalty/double(nThrows);
  std::cout << "Average LLH penalty over toys is " << AveragePenalty << std::endl;

  // Calculate exact p-value instead of binned
  unsigned int Accept_PredFluc = 0;
  unsigned int Accept_DrawFluc = 0;
  for (unsigned int i = 0; i < nThrows; ++i)
  {
    if (LLH_DataDraw_V[i] > LLH_DrawFlucDraw_V[i]) Accept_DrawFluc++;
    if (LLH_DataDraw_V[i] > LLH_PredFluc_V[i]) Accept_PredFluc++;
  }
  const double pvalue_DrawFluc = double(Accept_DrawFluc)/double(nThrows);
  const double pvalue_PredFluc = double(Accept_PredFluc)/double(nThrows);

  std::cout << "Calculated exact p-value using Fluctuation of Draw: "       << pvalue_DrawFluc << std::endl;
  std::cout << "Calculated exact p-value using Fluctuation of Prediction: " << pvalue_PredFluc << std::endl;
}

// *******************
// Make the cut LLH histogram
void SampleSummary::MakeCutLLH() {
// *******************
  Outputfile->cd();
  MakeCutLLH1D(lnLHist);
  MakeCutLLH1D(lnLHist_drawfluc);
  MakeCutLLH1D(lnLHist_drawdata);
  MakeCutLLH1D(lnLHist_drawflucdraw);
  
  MakeCutLLH2D(lnLDrawHist);
  MakeCutLLH2D(lnLFlucHist);
  MakeCutLLH2D(lnLDrawHistRate);
  MakeCutLLH2D(lnLFlucHist_ProjectX);
}

// ****************
// Make the 1D cut distribution and give the 1D p-value
void SampleSummary::MakeCutLLH1D(TH1D *Histogram, double llh_ref) {
// ****************
  const double TotalIntegral = Histogram->Integral();
  double Above = 0.0;
  // Get the LLH reference from total llh or some reference histogram
  double llh_reference = 0.0;
  if (llh_ref >= 0) {
    llh_reference = llh_ref;
  } else {
    llh_reference = llh_total;
  }
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    const double xvalue = Histogram->GetBinCenter(i+1);
    if (xvalue >= llh_reference) {
      Above += Histogram->GetBinContent(i+1);
    }
  }
  const double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Write a TCanvas and make a line and a filled histogram
  TLine *TempLine = new TLine(llh_reference , Histogram->GetMinimum(), llh_reference, Histogram->GetMaximum());
  TempLine->SetLineColor(kBlack);
  TempLine->SetLineWidth(2);

  // Make the fill histogram
  TH1D *TempHistogram = (TH1D*)(Histogram->Clone());
  TempHistogram->SetFillStyle(1001);
  TempHistogram->SetFillColor(kRed);
  for (int i = 0; i < TempHistogram->GetNbinsX(); ++i) 
  {
    if (TempHistogram->GetBinCenter(i+1) < llh_reference) 
    {
      TempHistogram->SetBinContent(i+1, 0.0);
    }
  }

  TLegend *Legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);
  Legend->AddEntry(TempLine, Form("Reference LLH, %.0f, p-value=%.2f", llh_reference, pvalue), "l");
  Legend->AddEntry(Histogram, Form("LLH, #mu=%.1f#pm%.1f", Histogram->GetMean(), Histogram->GetRMS()), "l");
  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  Histogram->Draw();
  TempHistogram->Draw("same");
  TempLine->Draw("same");
  Legend->Draw("same");

  TempCanvas->Write();

  delete TempLine;
  delete TempHistogram;
  delete TempCanvas;
  delete Legend;
}


// ****************
// Make the 2D cut distribution and give the 2D p-value
void SampleSummary::MakeCutLLH2D(TH2D *Histogram) {
// ****************

  const double TotalIntegral = Histogram->Integral();
  // Count how many fills are above y=x axis
  // This is the 2D p-value
  double Above = 0.0;
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    const double xvalue = Histogram->GetXaxis()->GetBinCenter(i+1);
    for (int j = 0; j < Histogram->GetYaxis()->GetNbins(); ++j) {
      const double yvalue = Histogram->GetYaxis()->GetBinCenter(j+1);
      // We're only interested in being _ABOVE_ the y=x axis
      if (xvalue >= yvalue) {
        Above += Histogram->GetBinContent(i+1, j+1);
      }
    }
  }
  const double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Now add the TLine going diagonally
  double minimum = Histogram->GetXaxis()->GetBinLowEdge(1);
  if (Histogram->GetYaxis()->GetBinLowEdge(1) > minimum) {
    minimum = Histogram->GetYaxis()->GetBinLowEdge(1);
  }
  double maximum = Histogram->GetXaxis()->GetBinLowEdge(Histogram->GetXaxis()->GetNbins());
  if (Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins()) < maximum) {
    maximum = Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins());
    //KS: Extend by bin width to perfectly fit canvas
    maximum += Histogram->GetYaxis()->GetBinWidth(Histogram->GetYaxis()->GetNbins());
  }
  else maximum += Histogram->GetXaxis()->GetBinWidth(Histogram->GetXaxis()->GetNbins());
  TLine *TempLine = new TLine(minimum, minimum, maximum, maximum);
  TempLine->SetLineColor(kRed);
  TempLine->SetLineWidth(2);

  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  TempCanvas->cd();
  gStyle->SetPalette(81);
  Histogram->SetMinimum(-0.01);
  Histogram->Draw("colz");
  TempLine->Draw("same");

  TempCanvas->Write();
  delete TempLine;
  delete TempCanvas;
}

// ****************
// Make the 1D Event Rate Hist
void SampleSummary::MakeCutEventRate(TH1D *Histogram, const double DataRate) {
// ****************
    // For the event rate histogram add a TLine to the data rate
    TLine *TempLine = new TLine(DataRate, Histogram->GetMinimum(), DataRate, Histogram->GetMaximum());
    TempLine->SetLineColor(kRed);
    TempLine->SetLineWidth(2);
    // Also fit a Gaussian because why not?
    TF1 *Fitter = new TF1("Fit", "gaus", Histogram->GetBinLowEdge(1), Histogram->GetBinLowEdge(Histogram->GetNbinsX()+1));
    Histogram->Fit(Fitter, "RQ");
    Fitter->SetLineColor(kRed-5);
    // Calculate a p-value
    double Above = 0.0;
    for (int z = 0; z < Histogram->GetNbinsX(); ++z) {
      const double xvalue = Histogram->GetBinCenter(z+1);
      if (xvalue >= DataRate) {
        Above += Histogram->GetBinContent(z+1);
      }
    }
    const double pvalue = Above/Histogram->Integral();
    TLegend *Legend = new TLegend(0.4, 0.75, 0.98, 0.90);
    Legend->SetFillColor(0);
    Legend->SetFillStyle(0);
    Legend->SetLineWidth(0);
    Legend->SetLineColor(0);
    Legend->AddEntry(TempLine, Form("Data, %.0f, p-value=%.2f", DataRate, pvalue), "l");
    Legend->AddEntry(Histogram, Form("MC, #mu=%.1f#pm%.1f", Histogram->GetMean(), Histogram->GetRMS()), "l");
    Legend->AddEntry(Fitter, Form("Gauss, #mu=%.1f#pm%.1f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
    std::string TempTitle = std::string(Histogram->GetName());
    TempTitle += "_canv";
    TCanvas *TempCanvas = new TCanvas(TempTitle.c_str(), TempTitle.c_str(), 1024, 1024);
    TempCanvas->SetGridx();
    TempCanvas->SetGridy();
    TempCanvas->SetRightMargin(0.03);
    TempCanvas->SetBottomMargin(0.08);
    TempCanvas->SetLeftMargin(0.10);
    TempCanvas->SetTopMargin(0.06);
    TempCanvas->cd();
    Histogram->Draw();
    TempLine->Draw("same");
    Fitter->Draw("same");
    Legend->Draw("same");
    TempCanvas->Write();
    Histogram->Write();
        
    delete TempLine;
    delete TempCanvas;
    delete Fitter;
    delete Legend;
}

// ****************
// Make a ratio histogram
template<class HistType>
HistType* SampleSummary::RatioHists(HistType *NumHist, HistType *DenomHist) {
// ****************

  HistType *NumCopy = (HistType*)(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());
  NumCopy->Divide(DenomHist);

  return NumCopy;
}

// ****************
// Make a ratio th2poly
TH2Poly* SampleSummary::RatioPolys(TH2Poly *NumHist, TH2Poly *DenomHist) {
// ****************

  TH2Poly *NumCopy = (TH2Poly*)(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());

  for(int i = 1; i < NumCopy->GetNumberOfBins()+1; ++i)
  {
    NumCopy->SetBinContent(i,NumHist->GetBinContent(i)/DenomHist->GetBinContent(i));
  }

  return NumCopy;
}

// ****************
// Normalise a TH2Poly
void SampleSummary::NormaliseTH2Poly(TH2Poly* Histogram){
// ****************
  const double Integral = NoOverflowIntegral(Histogram);

  for(int j = 1; j < Histogram->GetNumberOfBins()+1; j++)
  {
    Histogram->SetBinContent(j, Histogram->GetBinContent(j)/Integral);
  }
}

// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH1D * const & DatHist, TH1D * const & MCHist, TH1D * const & W2Hist) {
// ****************
  const double llh = GetLLH(DatHist, MCHist, W2Hist);
  std::stringstream ss;
  ss << "_2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(55) << std::left << MCHist->GetName() << std::setw(10) << DatHist->Integral() << std::setw(10) << MCHist->Integral() << std::setw(10) << llh << std::endl;
}

// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH2Poly * const & DatHist, TH2Poly * const & MCHist, TH2Poly * const & W2Hist) {
// ****************
  const double llh = GetLLH(DatHist, MCHist, W2Hist);
  std::stringstream ss;
  ss << "_2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(55) << std::left << MCHist->GetName() << std::setw(10) << NoOverflowIntegral(DatHist) << std::setw(10) << NoOverflowIntegral(MCHist) << std::setw(10) << llh << std::endl;
}

// ****************
double SampleSummary::GetLLH(TH2Poly * const & DatHist, TH2Poly * const & MCHist, TH2Poly * const & W2Hist) {
// ****************
  double llh = 0.0;
  for (int i = 1; i < DatHist->GetNumberOfBins()+1; ++i)
  {
    const double data = DatHist->GetBinContent(i);
    const double mc = MCHist->GetBinContent(i);
    const double w2 = W2Hist->GetBinContent(i);
    llh += SamplePDF->getTestStatLLH(data, mc, w2);
  }
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// ****************
double SampleSummary::GetLLH(TH1D * const & DatHist, TH1D * const & MCHist, TH1D * const & W2Hist) {
// ****************
  double llh = 0.0;
  for (int i = 1; i <= DatHist->GetXaxis()->GetNbins(); ++i)
  {
    const double data = DatHist->GetBinContent(i);
    const double mc = MCHist->GetBinContent(i);
    const double w2 = W2Hist->GetBinContent(i);
    llh += SamplePDF->getTestStatLLH(data, mc, w2);
  }
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// ****************
void SampleSummary::PlotBetaParameters() {
// ****************

  // Make a new directory
  TDirectory *BetaDir = Outputfile->mkdir("BetaParameters");
  BetaDir->cd();

  std::cout<<"Writing Beta parameters"<<std::endl;
  TDirectory **DirBeta = new TDirectory*[nSamples];
  for (_int_ i = 0; i < nSamples; ++i)
  {
    // Make a new directory
    DirBeta[i] = BetaDir->mkdir((SampleNames[i]).c_str());
    DirBeta[i]->cd();

    // Loop over each pmu cosmu bin
    for (int j = 1; j < maxBins[i]+1; ++j)
    {
      const double data = DataHist[i]->GetBinContent(j);
      const double mc = NominalHist[i]->GetBinContent(j);
      const double w2 = W2NomHist[i]->GetBinContent(j);

      const double BetaPrior = GetBetaParameter(data, mc, w2, likelihood);

      TLine *TempLine = new TLine(BetaPrior, BetaHist[i][j]->GetMinimum(), BetaPrior, BetaHist[i][j]->GetMaximum());
      TempLine->SetLineColor(kRed);
      TempLine->SetLineWidth(2);

      // Also fit a Gaussian because why not?
      TF1 *Fitter = new TF1("Fit", "gaus", BetaHist[i][j]->GetBinLowEdge(1), BetaHist[i][j]->GetBinLowEdge(BetaHist[i][j]->GetNbinsX()+1));
      BetaHist[i][j]->Fit(Fitter, "RQ");
      Fitter->SetLineColor(kRed-5);

      TLegend *Legend = new TLegend(0.4, 0.75, 0.98, 0.90);
      Legend->SetFillColor(0);
      Legend->SetFillStyle(0);
      Legend->SetLineWidth(0);
      Legend->SetLineColor(0);
      Legend->AddEntry(TempLine, Form("Prior #mu=%.4f, N_{data}=%.0f", BetaPrior, data), "l");
      Legend->AddEntry(BetaHist[i][j], Form("Post, #mu=%.4f#pm%.4f", BetaHist[i][j]->GetMean(), BetaHist[i][j]->GetRMS()), "l");
      Legend->AddEntry(Fitter, Form("Gauss, #mu=%.4f#pm%.4f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
      std::string TempTitle = std::string(BetaHist[i][j]->GetName());

      TempTitle += "_canv";
      TCanvas *TempCanvas = new TCanvas(TempTitle.c_str(), TempTitle.c_str(), 1024, 1024);
      TempCanvas->SetGridx();
      TempCanvas->SetGridy();
      TempCanvas->SetRightMargin(0.03);
      TempCanvas->SetBottomMargin(0.08);
      TempCanvas->SetLeftMargin(0.10);
      TempCanvas->SetTopMargin(0.06);
      TempCanvas->cd();
      BetaHist[i][j]->Draw();
      TempLine->Draw("same");
      Fitter->Draw("same");
      Legend->Draw("same");
      TempCanvas->Write();
      BetaHist[i][j]->Write();

      delete TempLine;
      delete TempCanvas;
      delete Fitter;
      delete Legend;
    }
    DirBeta[i]->Write();
    delete DirBeta[i];
  }
  delete[] DirBeta;
  BetaDir->Write();
  delete BetaDir;

  Outputfile->cd();
}

// ****************
// Make a projection
void SampleSummary::StudyKinematicCorrelations() {
// ****************
  std::cout<<" Calculating Correlations"<<std::endl;

  TStopwatch timer;
  timer.Start();

    // Data vs Draw for 1D projection
  double* NEvents_Sample = new double[nSamples];
  double event_rate = 0.;

  // The output tree we're going to write to
  TTree* Event_Rate_Tree = new TTree("Event_Rate_draws", "Event_Rate_draws");
  Event_Rate_Tree->Branch("Event_Rate", &event_rate);
  // Loop over the samples and set the addresses of the variables to write to file
  for (_int_ i = 0; i < nSamples; ++i)
  {
    // Get the name
    std::string SampleName = SampleNames[i];
    //CW: Also strip out - signs because it messes up TBranches
    while (SampleName.find("-") != std::string::npos) {
      SampleName.replace(SampleName.find("-"), 1, std::string("_"));
    }
    Event_Rate_Tree->Branch((SampleName+"_Event_Rate").c_str(), &NEvents_Sample[i]);
  }

  // Holds the total event rate
  TH1D *EventHist = new TH1D("EventHist", "Total Event Rate", 100, 1, -1);
  EventHist->GetXaxis()->SetTitle("Total event rate");
  EventHist->GetYaxis()->SetTitle("Counts");
  EventHist->SetLineWidth(2);

  // Holds the event rate for the distribution
  TH1D **SumHist = new TH1D*[nSamples];

  for (int i = 0; i < nSamples; ++i)
  {
    std::string name = std::string(NominalHist[i]->GetName());
    name = name.substr(0, name.find("_nom"));

    SumHist[i] = new TH1D((name+"_sum").c_str(),(name+"_sum").c_str(), 100, 1, -1);
    SumHist[i]->GetXaxis()->SetTitle("N_{events}");
    SumHist[i]->GetYaxis()->SetTitle("Counts");
    double Integral = NoOverflowIntegral(DataHist[i]);
    std::stringstream ss;
    ss << Integral;
    SumHist[i]->SetTitle((std::string(SumHist[i]->GetTitle())+"_"+ss.str()).c_str());
  }

  for (unsigned int it = 0; it < nThrows; ++it)
  {
    double event_rate_temp = 0.;
    // Loop over the samples
    #ifdef MULTITHREAD
    #pragma omp parallel for reduction(+:event_rate_temp)
    #endif
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      NEvents_Sample[SampleNum] = NoOverflowIntegral(MCVector[it][SampleNum]);
      // Fill the sum histogram with the integral of the sampled distribution
      SumHist[SampleNum]->Fill(NEvents_Sample[SampleNum], WeightVector[it]);

      event_rate_temp += NEvents_Sample[SampleNum];
    } // end samples loop
    event_rate = event_rate_temp;
    EventHist->Fill(event_rate);
    Event_Rate_Tree->Fill();
  } //end loops over throws
  Event_Rate_Tree->Write();
  delete Event_Rate_Tree;

  double DataRate = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:DataRate)
  #endif
  for (_int_ i = 0; i < nSamples; ++i)
  {
    DataRate += NoOverflowIntegral(DataHist[i]);
  }
  MakeCutEventRate(EventHist, DataRate);
  delete EventHist;

  for (_int_ SampleNum = 0; SampleNum < nSamples; ++SampleNum)
  {
    Dir[SampleNum]->cd();
    //Make fancy event rate histogram
    MakeCutEventRate(SumHist[SampleNum], NoOverflowIntegral(DataHist[SampleNum]));
  }

  // Make a new direcotry
  TDirectory *CorrDir = Outputfile->mkdir("Correlations");
  CorrDir->cd();

  TMatrixDSym* SampleCorrelation = new TMatrixDSym(nSamples);

  TH2D*** SamCorr = new TH2D**[nSamples]();

  for (int i = 0; i < nSamples; ++i)
  {
    SamCorr[i] = new TH2D*[nSamples]();

    (*SampleCorrelation)(i,i) = 1.0;

    const double Min_i = SumHist[i]->GetXaxis()->GetBinLowEdge(1);
    const double Max_i = SumHist[i]->GetXaxis()->GetBinUpEdge(SumHist[i]->GetNbinsX()+1);
    for (int j = 0; j < nSamples; ++j)
    {
      const double Min_j = SumHist[j]->GetXaxis()->GetBinLowEdge(1);
      const double Max_j = SumHist[j]->GetXaxis()->GetBinUpEdge(SumHist[j]->GetNbinsX()+1);

      // TH2D to hold the Correlation
      SamCorr[i][j] = new TH2D(Form("SamCorr_%i_%i",i,j), Form("SamCorr_%i_%i",i,j), 70, Min_i, Max_i, 70, Min_j, Max_j);
      SamCorr[i][j]->SetDirectory(0);
      SamCorr[i][j]->SetMinimum(0);
      SamCorr[i][j]->GetXaxis()->SetTitle(SampleNames[i].c_str());
      SamCorr[i][j]->GetYaxis()->SetTitle(SampleNames[j].c_str());
      SamCorr[i][j]->GetZaxis()->SetTitle("Events");
    }
  }

  // Now we are sure we have the diagonal elements, let's make the off-diagonals
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < nSamples; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;

      for (unsigned int it = 0; it < nThrows; ++it)
      {
        SamCorr[i][j]->Fill(NoOverflowIntegral(MCVector[it][i]), NoOverflowIntegral(MCVector[it][j]));
      }
      SamCorr[i][j]->Smooth();

      // Get the Covariance for these two parameters
      (*SampleCorrelation)(i,j) = SamCorr[i][j]->GetCorrelationFactor();
      (*SampleCorrelation)(j,i) = (*SampleCorrelation)(i,j);

    }// End j loop
  }// End i loop

  TH2D* hSamCorr = new TH2D("Sample Correlation", "Sample Correlation", nSamples, 0, nSamples, nSamples, 0, nSamples);
  hSamCorr->SetDirectory(0);
  hSamCorr->GetZaxis()->SetTitle("Correlation");
  hSamCorr->SetMinimum(-1);
  hSamCorr->SetMaximum(1);
  hSamCorr->GetXaxis()->SetLabelSize(0.015);
  hSamCorr->GetYaxis()->SetLabelSize(0.015);

  // Loop over the Covariance matrix entries
  for (int i = 0; i < nSamples; ++i)
  {
    hSamCorr->GetXaxis()->SetBinLabel(i+1, SampleNames[i].c_str());

    for (int j = 0; j < nSamples; ++j)
    {
      hSamCorr->GetYaxis()->SetBinLabel(j+1, SampleNames[j].c_str());
      // The value of the Covariance
      const double corr = (*SampleCorrelation)(i,j);
      hSamCorr->SetBinContent(i+1, j+1, corr);
    }
  }

  hSamCorr->Draw("colz");
  hSamCorr->Write("Sample_Corr");
  delete hSamCorr;

  SampleCorrelation->Write("Sample_Correlation");
  delete SampleCorrelation;

  for (int i = 0; i < nSamples; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      // Skip the diagonal elements which we've already done above
      if (j == i) continue;
      SamCorr[i][j]->Write();
    }// End j loop
  }// End i loop

  for (int i = 0; i < nSamples; ++i)
  {
    for (int j = 0; j < nSamples; ++j)
    {
      delete SamCorr[i][j];
    }
    delete[] SamCorr[i];
  }
  delete[] SamCorr;

  //KS: This can take ages so better turn it off by default
  bool DoPerKinemBin = false;
  if(DoPerKinemBin)
  {
    //KS: Now the same but for kinematic bin of each sample
    for (int SampleNum = 0; SampleNum < nSamples; ++SampleNum)
    {
      TMatrixDSym* KinCorrelation = new TMatrixDSym(maxBins[SampleNum]);

      TH2D*** KinCorr = new TH2D**[maxBins[SampleNum]]();
      for (int i = 0; i < maxBins[SampleNum]; ++i)
      {
        KinCorr[i] = new TH2D*[maxBins[SampleNum]]();
        (*KinCorrelation)(i,i) = 1.0;

        const double Min_i = PosteriorHist[SampleNum][i+1]->GetXaxis()->GetBinLowEdge(1);
        const double Max_i = PosteriorHist[SampleNum][i+1]->GetXaxis()->GetBinUpEdge(PosteriorHist[SampleNum][i+1]->GetNbinsX()+1);

        //Get PolyBin
        TH2PolyBin* bin = (TH2PolyBin*)NominalHist[SampleNum]->GetBins()->At(i);
        // Just make a little fancy name
        std::stringstream ss2;
        ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
        ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";

        for (int j = 0; j < maxBins[SampleNum]; ++j)
        {
          const double Min_j = PosteriorHist[SampleNum][j+1]->GetXaxis()->GetBinLowEdge(1);
          const double Max_j = PosteriorHist[SampleNum][j+1]->GetXaxis()->GetBinUpEdge(PosteriorHist[SampleNum][j+1]->GetNbinsX()+1);

          // TH2D to hold the Correlation
          KinCorr[i][j] = new TH2D(Form("Kin_%i_%i_%i",SampleNum,i,j), Form("Kin_%i_%i_%i",SampleNum,i,j), 70, Min_i, Max_i, 70, Min_j, Max_j);
          KinCorr[i][j]->SetDirectory(0);
          KinCorr[i][j]->SetMinimum(0);

          KinCorr[i][j]->GetXaxis()->SetTitle(ss2.str().c_str());

          bin = (TH2PolyBin*)NominalHist[SampleNum]->GetBins()->At(j);
          // Just make a little fancy name
          std::stringstream ss3;
          ss3 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
          ss3 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";
          KinCorr[i][j]->GetYaxis()->SetTitle(ss3.str().c_str());

          KinCorr[i][j]->GetZaxis()->SetTitle("Events");
        }
      }
      // Now we are sure we have the diagonal elements, let's make the off-diagonals
      #ifdef MULTITHREAD
      #pragma omp parallel for
      #endif
      for (int i = 0; i < maxBins[SampleNum]; ++i)
      {
        for (int j = 0; j <= i; ++j)
        {
          // Skip the diagonal elements which we've already done above
          if (j == i) continue;

          for (unsigned int it = 0; it < nThrows; ++it)
          {
            KinCorr[i][j]->Fill(MCVector[it][SampleNum]->GetBinContent(i+1), MCVector[it][SampleNum]->GetBinContent(j+1));
          }
          KinCorr[i][j]->Smooth();

          // Get the Covariance for these two parameters
          (*KinCorrelation)(i,j) = KinCorr[i][j]->GetCorrelationFactor();
          (*KinCorrelation)(j,i) = (*KinCorrelation)(i,j);

        }// End j loop
      }// End i loop


      TH2D* hKinCorr = new TH2D(SampleNames[SampleNum].c_str(), SampleNames[SampleNum].c_str(), maxBins[SampleNum], 0, maxBins[SampleNum], maxBins[SampleNum], 0, maxBins[SampleNum]);
      hKinCorr->SetDirectory(0);
      hKinCorr->GetZaxis()->SetTitle("Correlation");
      hKinCorr->SetMinimum(-1);
      hKinCorr->SetMaximum(1);
      hKinCorr->GetXaxis()->SetLabelSize(0.015);
      hKinCorr->GetYaxis()->SetLabelSize(0.015);

      // Loop over the Covariance matrix entries
      for (int i = 0; i < maxBins[SampleNum]; ++i)
      {
        //Get PolyBin
        TH2PolyBin* bin = (TH2PolyBin*)NominalHist[SampleNum]->GetBins()->At(i);
        // Just make a little fancy name
        std::stringstream ss2;
        ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
        ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";
        hKinCorr->GetXaxis()->SetBinLabel(i+1, ss2.str().c_str());

        for (int j = 0; j < maxBins[SampleNum]; ++j)
        {
          bin = (TH2PolyBin*)NominalHist[SampleNum]->GetBins()->At(j);
          // Just make a little fancy name
          std::stringstream ss3;
          ss3 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
          ss3 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";
          KinCorr[i][j]->GetYaxis()->SetTitle(ss3.str().c_str());

          hKinCorr->GetYaxis()->SetBinLabel(j+1, ss3.str().c_str());
          // The value of the Covariance
          const double corr = (*KinCorrelation)(i,j);
          hKinCorr->SetBinContent(i+1, j+1, corr);
        }
      }

      hKinCorr->Draw("colz");
      hKinCorr->Write((SampleNames[SampleNum] + "_Corr").c_str());
      delete hKinCorr;

      KinCorrelation->Write((SampleNames[SampleNum] + "_Correlation").c_str());
      delete KinCorrelation;

      /*
      for (int i = 0; i < maxBins[SampleNum]; ++i)
      {
        for (int j = 0; j <= i; ++j)
        {
          // Skip the diagonal elements which we've already done above
          if (j == i) continue;
          KinCorr[i][j]->Write();
        }// End j loop
      }// End i loop
      */
      for (int i = 0; i < maxBins[SampleNum]; ++i)
      {
        for (int j = 0; j < maxBins[SampleNum]; ++j)
        {
          delete KinCorr[i][j];
        }
        delete[] KinCorr[i];
      }
      delete[] KinCorr;
    }//end loop over samples
  }//end if DoPerKinemBin
  else
  {
    std::cout<<" Not calculating correlations per each kinematic bin"<<std::endl;
  }

  if(DoByModePlots)
  {
    // Holds the total event rate by mode
    TH1D ** EventHist_ByMode = new TH1D*[Modes->GetNModes()+1];

    for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
    {
      std::string ModeName = Modes->GetMaCh3ModeName(j);
      EventHist_ByMode[j] = new TH1D(Form("EventHist_%s", ModeName.c_str()), Form("Total Event Rate %s", ModeName.c_str()), 100, 1, -1);
      EventHist_ByMode[j]->GetXaxis()->SetTitle("Total event rate");
      EventHist_ByMode[j]->GetYaxis()->SetTitle("Counts");
      EventHist_ByMode[j]->SetLineWidth(2);
    }

    //KS: Here we calculate total event rates for each mode, maybe not most efficient but can be improved in the future
    for (unsigned int it = 0; it < nThrows; ++it)
    {
      for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
      {
        double event_rate_temp = 0.;
        #ifdef MULTITHREAD
        #pragma omp parallel for reduction(+:event_rate_temp)
        #endif
        for (_int_ SampleNum = 0;  SampleNum < nSamples; SampleNum++)
        {
          event_rate_temp += NoOverflowIntegral(MCVectorByMode[it][SampleNum][j]);
        }
        EventHist_ByMode[j]->Fill(event_rate_temp);
      }
    }

    for (_int_ i = 0; i < Modes->GetNModes()+1; ++i)
    {
      MakeCutEventRate(EventHist_ByMode[i], DataRate);
    }

    TMatrixDSym* ModeCorrelation = new TMatrixDSym(Modes->GetNModes()+1);

    TH2D*** ModeCorr = new TH2D**[Modes->GetNModes()+1]();
    for (int i = 0; i < Modes->GetNModes()+1; ++i)
    {
      ModeCorr[i] = new TH2D*[Modes->GetNModes()+1]();

      (*ModeCorrelation)(i,i) = 1.0;

      const double Min_i = EventHist_ByMode[i]->GetXaxis()->GetBinLowEdge(1);
      const double Max_i = EventHist_ByMode[i]->GetXaxis()->GetBinUpEdge(EventHist_ByMode[i]->GetNbinsX()+1);
      for (int j = 0; j < Modes->GetNModes()+1; ++j)
      {
        const double Min_j = EventHist_ByMode[j]->GetXaxis()->GetBinLowEdge(1);
        const double Max_j = EventHist_ByMode[j]->GetXaxis()->GetBinUpEdge(EventHist_ByMode[j]->GetNbinsX()+1);

        // TH2D to hold the Correlation
        ModeCorr[i][j] = new TH2D(Form("ModeCorr_%i_%i",i,j), Form("ModeCorr_%i_%i",i,j), 70, Min_i, Max_i, 70, Min_j, Max_j);
        ModeCorr[i][j]->SetDirectory(0);
        ModeCorr[i][j]->SetMinimum(0);
        ModeCorr[i][j]->GetXaxis()->SetTitle(Modes->GetMaCh3ModeName(i).c_str());
        ModeCorr[i][j]->GetYaxis()->SetTitle(Modes->GetMaCh3ModeName(j).c_str());
        ModeCorr[i][j]->GetZaxis()->SetTitle("Events");
      }
    }

    // Now we are sure we have the diagonal elements, let's make the off-diagonals
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (int i = 0; i < Modes->GetNModes()+1; ++i)
    {
      for (int j = 0; j <= i; ++j)
      {
        // Skip the diagonal elements which we've already done above
        if (j == i) continue;

        for (unsigned int it = 0; it < nThrows; ++it)
        {
          double Integral_X = 0.;
          double Integral_Y = 0.;
          for (int SampleNum = 0; SampleNum < nSamples; ++SampleNum)
          {
            Integral_X += NoOverflowIntegral(MCVectorByMode[it][SampleNum][i]);
            Integral_Y += NoOverflowIntegral(MCVectorByMode[it][SampleNum][j]);
          }
          ModeCorr[i][j]->Fill(Integral_X, Integral_Y);
        }
        ModeCorr[i][j]->Smooth();

        // Get the Covariance for these two parameters
        (*ModeCorrelation)(i,j) = ModeCorr[i][j]->GetCorrelationFactor();
        (*ModeCorrelation)(j,i) = (*ModeCorrelation)(i,j);

      }// End j loop
    }// End i loop

    TH2D* hModeCorr = new TH2D("Mode Correlation", "Mode Correlation", Modes->GetNModes()+1, 0, Modes->GetNModes()+1, Modes->GetNModes()+1, 0, Modes->GetNModes()+1);
    hModeCorr->SetDirectory(0);
    hModeCorr->GetZaxis()->SetTitle("Correlation");
    hModeCorr->SetMinimum(-1);
    hModeCorr->SetMaximum(1);
    hModeCorr->GetXaxis()->SetLabelSize(0.015);
    hModeCorr->GetYaxis()->SetLabelSize(0.015);

    // Loop over the Covariance matrix entries
    for (int i = 0; i < Modes->GetNModes()+1; ++i)
    {
      hModeCorr->GetXaxis()->SetBinLabel(i+1, Modes->GetMaCh3ModeName(i).c_str());

      for (int j = 0; j < Modes->GetNModes()+1; ++j)
      {
        hModeCorr->GetYaxis()->SetBinLabel(j+1, Modes->GetMaCh3ModeName(j).c_str());
        // The value of the Covariance
        const double corr = (*ModeCorrelation)(i,j);
        hModeCorr->SetBinContent(i+1, j+1, corr);
      }
    }
    hModeCorr->Draw("colz");
    hModeCorr->Write("Mode_Corr");
    delete hModeCorr;

    for (int i = 0; i < Modes->GetNModes()+1; ++i)
    {
      for (int j = 0; j <= i; ++j)
      {
        // Skip the diagonal elements which we've already done above
        if (j == i) continue;
        ModeCorr[i][j]->Write();
      }// End j loop
    }// End i loop

    for (int i = 0; i < Modes->GetNModes()+1; ++i)
    {
      for (int j = 0; j < Modes->GetNModes()+1; ++j)
      {
        delete ModeCorr[i][j];
      }
      delete[] ModeCorr[i];
    }
    delete[] ModeCorr;
    ModeCorrelation->Write("Mode_Correlation");
    delete ModeCorrelation;

    for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
    {
      delete EventHist_ByMode[j];
    }
    delete[] EventHist_ByMode;
  }

  for (_int_ i = 0; i < nSamples; ++i)
  {
    delete SumHist[i];
  }
  delete[] SumHist;

  CorrDir->Close();
  delete CorrDir;

  timer.Stop();
  std::cout << "Calculating correlations took " << timer.RealTime() << "s" << std::endl;

  Outputfile->cd();
}


// ****************
// Make a projection
TH1D* SampleSummary::ProjectHist(TH2D* Histogram, bool ProjectX) {
// ****************

  TH1D* Projection = NULL;
  std::string name;
  if (ProjectX) {
    name = std::string(Histogram->GetName()) + "_x";
    Projection = Histogram->ProjectionX(name.c_str(), 1, Histogram->GetYaxis()->GetNbins(), "e");
  } else {
    name = std::string(Histogram->GetName()) + "_y";
    Projection = Histogram->ProjectionY(name.c_str(), 1, Histogram->GetXaxis()->GetNbins(), "e");
  }

  return Projection;
}

// ****************
// Make a projection
TH1D* SampleSummary::ProjectPoly(TH2Poly* Histogram, const bool ProjectX, const _int_ selection, const bool MakeErrorHist) {
// ****************

  std::vector<double> xbins;
  std::vector<double> ybins;

  SamplePDF->SetupBinning(selection, xbins, ybins);
  TH1D* Projection = NULL;
  std::string name;
  if (ProjectX) {
    name = std::string(Histogram->GetName()) + "_x";
    Projection = PolyProjectionX(Histogram, name.c_str(), xbins, MakeErrorHist);
  } 
  else {
    name = std::string(Histogram->GetName()) + "_y";
    Projection = PolyProjectionY(Histogram, name.c_str(), ybins, MakeErrorHist);
  }

  return Projection;
}


// ****************
//KS: We have two methods how to apply statistical fluctuation standard is faster hence is default
void SampleSummary::MakeFluctuatedHistogram(TH1D *FluctHist, TH1D* PolyHist){
// ****************
    if(StandardFluctuation) MakeFluctuatedHistogramStandard(FluctHist, PolyHist);
    else MakeFluctuatedHistogramAlternative(FluctHist, PolyHist);
}

// ****************
//KS: We have two methods how to apply statistical fluctuation standard is faster hence is default
void SampleSummary::MakeFluctuatedHistogram(TH2Poly *FluctHist, TH2Poly* PolyHist){
// ****************
    if(StandardFluctuation) MakeFluctuatedHistogramStandard(FluctHist, PolyHist);
    else MakeFluctuatedHistogramAlternative(FluctHist, PolyHist);
}

// ****************
// Make Poisson Fluctuation of TH1D hist
void SampleSummary::MakeFluctuatedHistogramStandard(TH1D *FluctHist, TH1D* PolyHist){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0);

  for (int i = 1; i <= PolyHist->GetXaxis()->GetNbins(); ++i)  
  {
    // Get the posterior predictive bin content
    const double MeanContent = PolyHist->GetBinContent(i);
    // Get a Poisson fluctuation of the content
    const double Random = rnd->PoissonD(MeanContent);
    // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
    FluctHist->SetBinContent(i,Random);
  }
}

// ****************
// Make Poisson Fluctuation of TH2Poly hist
void SampleSummary::MakeFluctuatedHistogramStandard(TH2Poly *FluctHist, TH2Poly* PolyHist){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0, 0.0);

  for (int i = 1; i < FluctHist->GetNumberOfBins()+1; ++i) 
  {
    // Get the posterior predictive bin content
    const double MeanContent = PolyHist->GetBinContent(i);
    // Get a Poisson fluctuation of the content
    const double Random = rnd->PoissonD(MeanContent);
    // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
    FluctHist->SetBinContent(i,Random);
  }
}

// ****************
// Make Poisson Fluctuation of TH1D hist
void SampleSummary::MakeFluctuatedHistogramAlternative(TH1D* FluctHist, TH1D* PolyHist){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0);

  const double evrate = PolyHist->Integral();
  const int num = rnd->PoissonD(evrate);
  int count = 0;
  while(count < num)
  {
    const double candidate = PolyHist->GetRandom();
    FluctHist->Fill(candidate);                 
    count++;
  }
}

// ****************
// Make Poisson fluctuation of TH2Poly hist
void SampleSummary::MakeFluctuatedHistogramAlternative(TH2Poly *FluctHist, TH2Poly* PolyHist){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0, 0.0);

  const double evrate = NoOverflowIntegral(PolyHist);
  const int num = rnd->PoissonD(evrate);
  int count = 0;
  while(count < num)
  {
    const int iBin = GetRandomPoly2(PolyHist);
    FluctHist->SetBinContent(iBin, FluctHist->GetBinContent(iBin) + 1);                
    count++;
  }
}

// ****************
//KS: ROOT developers were too lazy do develop getRanom2 for TH2Poly, this implementation is based on:
// https://root.cern.ch/doc/master/classTH2.html#a883f419e1f6899f9c4255b458d2afe2e
int SampleSummary::GetRandomPoly2(const TH2Poly* PolyHist){
// ****************
  const int nbins = PolyHist->GetNumberOfBins();
  const double r1 = rnd->Rndm();
  
  double* fIntegral = new double[nbins+2];
  fIntegral[0] = 0.0;

  //KS: This is custom version of ComputeIntegral, once again ROOT was lazy :(
  for (int i = 1; i < nbins+1; ++i) 
  {
    fIntegral[i] = 0.0;
    const double content = PolyHist->GetBinContent(i);
    fIntegral[i] += fIntegral[i - 1] + content;
  }
  for (Int_t bin = 1; bin < nbins+1; ++bin)  fIntegral[bin] /= fIntegral[nbins];
  fIntegral[nbins+1] = PolyHist->GetEntries();
  
  //KS: We just return one rather then X and Y, this way we can use SetBinContent rather than Fill, which is faster 
  int iBin = TMath::BinarySearch(nbins, fIntegral, (Double_t) r1);
  //KS: Have to increment because TH2Poly has stupid offset arghh
  iBin += 1;
  
  delete[] fIntegral;
  return iBin;
}

// ****************
// Get the mode error from a TH1D
double SampleSummary::GetModeError(TH1D* hpost){
// ****************
  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();

  // The total integral of the posterior
  const double Integral = hpost->Integral();

  double sum = 0.0;
  int LowBin = MaxBin;
  int HighBin = MaxBin;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum/Integral < 0.6827 && (LowBin >= 0 || HighBin < hpost->GetNbinsX()+1) )
  {
    LowCon = 0.0;
    HighCon = 0.0;
    //KS:: Move further only if you haven't reached histogram end
    if(LowBin > 1)
    {
      LowBin--;
      LowCon = hpost->GetBinContent(LowBin);
    }
    if(HighBin < hpost->GetNbinsX())
    {
      HighBin++;
      HighCon = hpost->GetBinContent(HighBin);
    }

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/Integral > 0.6827 && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/Integral > 0.6827 && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }
  }

  double Mode_Error = 0.0;
  if (LowCon > HighCon) {
    Mode_Error = std::fabs(hpost->GetBinCenter(LowBin)-hpost->GetBinCenter(MaxBin));
  } else {
    Mode_Error = std::fabs(hpost->GetBinCenter(HighBin)-hpost->GetBinCenter(MaxBin));
  }

  return Mode_Error;
}

// ****************
void SampleSummary::StudyBIC(){
// ****************
  //make fancy event rate histogram
  double DataRate = 0.0;
  double BinsRate = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:DataRate, BinsRate)
  #endif
  for (_int_ i = 0; i < nSamples; ++i)
  {
    if (DataHist[i] == NULL) continue;
    DataRate += NoOverflowIntegral(DataHist[i]);
    BinsRate += maxBins[i];
  }

  const double EventRateBIC = GetBIC(llh_total, DataRate, nModelParams);
  const double BinBasedBIC = GetBIC(llh_total, BinsRate, nModelParams);
  std::cout << "******************************" << std::endl;
  std::cout << "Calculated Bayesian Information Criterion using global number of events: " << EventRateBIC << std::endl;
  std::cout << "Calculated Bayesian Information Criterion using global number of bins: " << BinBasedBIC << std::endl;
  std::cout << "Additional info: nModelParams "<<nModelParams<<" DataRate: "<<DataRate<<" BinsRate: "<<BinsRate<<std::endl;
  std::cout << "******************************" << std::endl;
}

// ****************
// Get the Deviance Information Criterion (DIC), for more see https://doi.org/10.1111/1467-9868.00353 or https://search.r-project.org/CRAN/refmans/BRugs/html/dic.html
void SampleSummary::StudyDIC() {
// ****************

  //The posterior mean of the deviance
  double Dbar = 0.;

  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:Dbar)
  #endif
  for (unsigned int i = 0; i < nThrows; ++i)
  {
    double LLH_temp = 0.;
    for (_int_ SampleNum = 0;  SampleNum < nSamples; ++SampleNum)
    {
      LLH_temp += GetLLH(DataHist[SampleNum], MCVector[i][SampleNum], W2MCVector[i][SampleNum]);
    }
    Dbar += LLH_temp;
  }
  Dbar = Dbar / nThrows;

  // A point estimate of the deviance
  const double Dhat = llh_total;

  //Effective number of parameters
  const double p_D = std::fabs(Dbar - Dhat);

  //Actual test stat
  const double DIC_stat = Dhat + 2 * p_D;
  std::cout<<"Effective number of parameters following DIC formalism is equal to "<< p_D <<std::endl;
  std::cout<<"DIC test statistic = "<< DIC_stat <<std::endl;
  std::cout << "******************************" << std::endl;

  return;
}


// ****************
//Fast and thread safe fill of violin histogram, it assumes both histograms have the same binning
void SampleSummary::FastViolinFill(TH2D* violin, TH1D* hist_1d){
// ****************
  for (int x = 0; x < violin->GetXaxis()->GetNbins(); ++x) 
  {
    const double y = violin->GetYaxis()->FindBin(hist_1d->GetBinContent(x+1));
    violin->SetBinContent(x+1, y,  violin->GetBinContent(x+1, y)+1);
  }
  return;
}
