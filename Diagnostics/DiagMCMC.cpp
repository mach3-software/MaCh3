// Classes etc declared in DiagMCMC
#include "DiagMCMC.h"


// The main program
int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << "./DiagMCMC MaCh3_output_file.root" << std::endl;
    exit(-1);
  }

  // TODO:
  // Limit the number of steps that go in
  // Print recommended burn-in?
  // Add red to bad plots?
  // Read covariance matrix from MCMC output instead

  std::string filename = argv[1];

  // Construct the diagnosis handler
  MCMCDiag_Handler handler(filename);

  // Set the number of entries
  //handler.SetNEntries(1E5);

  // Diagnose the MCMC
  handler.DiagMCMC();

  return 0;
}


// **************************
// Class implementations
// The constructor
MCMCDiag_Handler::MCMCDiag_Handler(std::string FileName) {
  // **************************

  OutputFileName = FileName;
  while (OutputFileName.find(".root") != std::string::npos) {
    OutputFileName = OutputFileName.substr(0, OutputFileName.find(".root"));
  }
  OutputFileName += "_MCMC_diag.root";


  InputFile = new TFile(FileName.c_str(), "OPEN");
  Chain = new TChain("posteriors", "");
  Chain->Add(FileName.c_str());

  // Get the settings for the MCMC
  TTree *Settings = (TTree*)(InputFile->Get("Settings"));
  if (Settings == NULL) {
    std::cerr << "Didn't find Settings tree in MCMC file " << InputFile << std::endl;
    std::cerr << "Will try lowercase" << std::endl;
    InputFile->ls();
    Settings = (TTree*)(InputFile->Get("settings"));
    if (Settings == NULL) throw;
  }

  // Get the xsec Covariance matrix
  std::string *XSecInput = 0;
  if (Settings->SetBranchAddress("XsecCov", &XSecInput) < 0) {
    Settings->Print();
    std::cerr << "Couldn't find XsecCov branch in output" << std::endl;
    throw;
  }

  // Write the XSecCov and TempFile
  Settings->GetEntry(0);

  // Delete the TTrees and the input file handle since we've now got the settings we need
  delete Settings;

  XsecCov = std::string("../")+(*XSecInput);

  std::cout << "--------------------" << std::endl;
  std::cout << "File for study:       " << FileName << std::endl;
  std::cout << "Starting MCMC diagnostics checker" << std::endl;
  std::cout << "I will write to " << OutputFileName << std::endl;
  std::cout << "--------------------" << std::endl;

  // Initialise the parameters
  nParams   = 0;
  nXsec     = 0;
  nFlux     = 0;
  nSamples  = 0;
  nSysts    = 0;
  nEntries  = -1;

}



// **************************
// Diagnose the MCMC
void MCMCDiag_Handler::DiagMCMC() {
  // **************************

  // Check the input file for reasonable branches
  CountBranches();

  // Get names, central values, priors from the input covariance matrix
  ReadXsecCov();

  // Read the input tree and copy entries from disk into memory
  ReadInputTree();

  // Draw the simple trace matrices
  ParamTraces();

  // Get the batched means
  BatchedMeans();

  // Draw the auto-correlations
  AutoCorrelation();

  // Draw acceptance Probability
  AcceptanceProbabilities();
  
  // Draw the logL
  // DrawL
  // Draw 

  // Close the output file
  OutputFile->Close();
}

// **************************
// Check the input branches
void MCMCDiag_Handler::CountBranches() {
  // **************************

  std::cout << "Counting branches..." << std::endl;

  // Get the list of branches
  TObjArray* brlis = (TObjArray*)Chain->GetListOfBranches();

  // Get the number of branches
  nBranches = brlis->GetEntries();

  // Loop over the number of branches
  // Find the name and how many of each systematic we have
  for (int i = 0; i < nBranches; i++) {

    // Get the TBranch and its name
    TBranch *br = (TBranch*)brlis->At(i);
    std::string BranchName = std::string(br->GetName());

    // Count how many entries we have in branch
    // Can also be specified by user: only set nEntries if user hasn't specified it
    if (nEntries == -1) {
      nEntries = br->GetEntries();
    }

    // Count the number of cross-section we have
    if (BeginString(BranchName,"xsec_")) {
      ParamName_v.push_back(BranchName);
      nParams++;
      nXsec++;
    }

    // Why not include the flux parameters too...
    else if (BeginString(BranchName,"b_")) {
      ParamName_v.push_back(BranchName);
      nParams++;
      nFlux++;
    }

    else if (BeginString(BranchName,"LogL_sample_")) {
      SampleName_v.push_back(BranchName);
      nSamples++;
    }

    else if (BeginString(BranchName,"LogL_systematic_")) {
      SystName_v.push_back(BranchName);
      nSysts++;
    }
  }

  std::cout << "# of branches:                " << nBranches << std::endl;
  std::cout << "# of systematics of interest: " << nParams << std::endl;

  // Initialise the ParamValues
  ParamValues = new double*[nEntries]();
  SampleValues = new double*[nEntries]();
  SystValues = new double*[nEntries]();
  AccProbValues = new double[nEntries]();
  for (int i = 0; i < nEntries; ++i) {
    ParamValues[i] = new double[nParams]();
    SampleValues[i] = new double[nSamples]();
    SystValues[i] = new double[nSysts]();
    for (int j = 0; j < nParams; ++j) {
      ParamValues[i][j] = -999.99;
    }
    for (int j = 0; j < nSamples; ++j) {
      SampleValues[i][j] = -999.99;
    }
    for (int j = 0; j < nSysts; ++j) {
      SystValues[i][j] = -999.99;
    }
    AccProbValues[i] = -999.99;
  }

  // Initialise the sums
  ParamSums = new double[nParams]();
  for (int i = 0; i < nParams; ++i) {
    ParamSums[i] = 0.0;
  }
}

// *****************
// Compare strings
bool MCMCDiag_Handler::BeginString(std::string NameString, std::string CompareString) {
// *****************
  bool StartsWith = false;
  if (NameString.compare(0, CompareString.length(), CompareString) == 0) {
    StartsWith = true;
  }
  return StartsWith;
}

// *****************
// We need the different steps for each parameter
// So let's get all the parameter variations for one given parameter
void MCMCDiag_Handler::ReadInputTree() {
  // *****************

  std::cout << "Reading input tree..." << std::endl;
  TStopwatch clock;
  clock.Start();

  // Set all the branches to off
  Chain->SetBranchStatus("*", false);
  // Turn on the branches which we want for parameters
  for (int i = 0; i < nParams; ++i) {
    Chain->SetBranchStatus(ParamName_v[i].c_str(), true);
  }

  // Turn on the branches which we want for LogL sample
  for (int i = 0; i < nSamples; ++i) {
    Chain->SetBranchStatus(SampleName_v[i].c_str(), true);
  }

  // Turn on the branches which we want for LogL systs
  for (int i = 0; i < nSysts; ++i) {
    Chain->SetBranchStatus(SystName_v[i].c_str(), true);
  }

  // Turn on the branches which we want for acc prob
  Chain->SetBranchStatus("accProb", true);
  
  // 10 entries output
  int countwidth = nEntries/10;

  // Can also do the batched means here to minimize excessive loops
  nBatches = 20;
  // The length of each batch
  int BatchLength = nEntries/nBatches+1;
  BatchedAverages = new double*[nBatches]();
  AccProbBatchedAverages = new double[nBatches]();
  for (int i = 0; i < nBatches; ++i) {
    BatchedAverages[i] = new double[nParams];
    AccProbBatchedAverages[i] = 0;
    for (int j = 0; j < nParams; ++j) {
      BatchedAverages[i][j] = 0.0;
    }
  }

  // Loop over the entries
  for (int i = 0; i < nEntries; ++i) {

    if (i % countwidth == 0) {
      std::cout << i << "/" << nEntries << " (" << double(i)/double(nEntries)*100. << "%)" << std::endl;
    }

    // Set the branch addresses for params
    for (int j = 0; j < nParams; ++j) {
      Chain->SetBranchAddress(ParamName_v[j].c_str(), &ParamValues[i][j]);
    }

    // Set the branch addresses for samples
    for (int j = 0; j < nSamples; ++j) {
      Chain->SetBranchAddress(SampleName_v[j].c_str(), &SampleValues[i][j]);
    }

    // Set the branch addresses for systematics
    for (int j = 0; j < nSysts; ++j) {
      Chain->SetBranchAddress(SystName_v[j].c_str(), &SystValues[i][j]);
    }
      
    // Set the branch addresses for Acceptance Probability
    Chain->SetBranchAddress("accProb", &AccProbValues[i]);

    
    // Fill up the ParamValues array
    Chain->GetEntry(i);

    // Find which batch the event belongs in
    int BatchNumber = -1;
    // I'm so lazy! But it's OK, the major overhead here is GetEntry: saved by ROOT!
    for (int j = 0; j < nBatches; ++j) {
      if (i < (j+1)*BatchLength) {
        BatchNumber = j;
        break;
      }
    }

    // Fill up the sum for each j param
    for (int j = 0; j < nParams; ++j) {
      ParamSums[j] += ParamValues[i][j];
      BatchedAverages[BatchNumber][j] += ParamValues[i][j];
    }
    
      //KS: Could easyli add this to above loop but I accProb is different beast so better keep it like this
      AccProbBatchedAverages[BatchNumber] += AccProbValues[i];
  }

  clock.Stop();

  std::cout << "Took " << clock.RealTime() << "s to finish for " << nEntries << " events" << std::endl;

  // Make the sums into average
  for (int i = 0; i < nParams; ++i) {
    ParamSums[i] /= nEntries;
    for (int j = 0; j < nBatches; ++j) {
      // Divide by the total number of events in the batch
      BatchedAverages[j][i] /= BatchLength;
      if(i==0) AccProbBatchedAverages[j] /= BatchLength; //KS: we have only one accProb, keep it like this for now
    }
  }

  delete Chain;

  // And make our sweet output file
  OutputFile = new TFile(OutputFileName.c_str(), "RECREATE");
}


// *****************
// Draw trace plots of the parameters
// i.e. parameter vs step
void MCMCDiag_Handler::ParamTraces() {
  // *****************

  std::cout << "Making trace plots..." << std::endl;

  // Make the TH1Ds
  TraceParamPlots = new TH1D*[nParams];
  TraceSamplePlots = new TH1D*[nSamples];
  TraceSystsPlots = new TH1D*[nSysts];

  // Set the titles and limits for TH2Ds
  for (int j = 0; j < nParams; ++j) {

    std::string HistName = XsecNames[j] +"_"+ ParamName_v[j]+"_Trace";
    //if (BeginString(ParamName_v[j],"xsec_") && j >= nFlux) {
    //  HistName = XsecNames[j-nFlux];
    //}

    TraceParamPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceParamPlots[j]->GetXaxis()->SetTitle("Step");
    TraceParamPlots[j]->GetYaxis()->SetTitle("Parameter Variation");
  }

  for (int j = 0; j < nSamples; ++j) {
    std::string HistName = SampleName_v[j];
    TraceSamplePlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceSamplePlots[j]->GetXaxis()->SetTitle("Step");
    TraceSamplePlots[j]->GetYaxis()->SetTitle("Sample -logL");
  }

  for (int j = 0; j < nSysts; ++j) {
    std::string HistName = SystName_v[j];
    TraceSystsPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nEntries, 0, nEntries);
    TraceSystsPlots[j]->GetXaxis()->SetTitle("Step");
    TraceSystsPlots[j]->GetYaxis()->SetTitle("Systematic -logL");
  }

  // Have now made the empty TH2Ds, now for writing content to them!

  // Loop over the number of parameters to draw their traces
  // Each histogram
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
#pragma omp parallel for
#endif
  for (int i = 0; i < nEntries; ++i) {
    // Set bin content for the ith bin to the parameter values
    for (int j = 0; j < nParams; ++j) {
      TraceParamPlots[j]->SetBinContent(i, ParamValues[i][j]);
    }

    for (int j = 0; j < nSamples; ++j) {
      TraceSamplePlots[j]->SetBinContent(i, SampleValues[i][j]);
    }

    for (int j = 0; j < nSysts; ++j) {
      TraceSystsPlots[j]->SetBinContent(i, SystValues[i][j]);
    }
  }


  // Write the output and delete the TH2Ds
  TDirectory *TraceDir = OutputFile->mkdir("Trace");
  TraceDir->cd();
  for (int j = 0; j < nParams; ++j) {
    // Fit a linear function to the traces
    TF1 *Fitter = new TF1("Fitter","[0]", int(nEntries/2), nEntries);
    Fitter->SetLineColor(kRed);
    TraceParamPlots[j]->Fit("Fitter","Rq");
    TraceParamPlots[j]->Write();
    delete Fitter;
    delete TraceParamPlots[j];
  }
  delete[] TraceParamPlots;

  TDirectory *LLDir = OutputFile->mkdir("LogL");
  LLDir->cd();
  for (int j = 0; j < nSamples; ++j) {
    TraceSamplePlots[j]->Write();
    delete TraceSamplePlots[j];
    delete SampleValues[j];
  }
  delete[] TraceSamplePlots;
  delete[] SampleValues;

  for (int j = 0; j < nSysts; ++j) {
    TraceSystsPlots[j]->Write();
    delete TraceSystsPlots[j];
    delete SystValues[j];
  }
  delete[] TraceSystsPlots;
  delete[] SystValues;
}


// *********************************
void MCMCDiag_Handler::AutoCorrelation() {
  // *********************************

  std::cout << "Making auto-correlations..." << std::endl;
  const int nLags = 25000;

  // The sum of (Y-Ymean)^2 over all steps for each parameter
  double **DenomSum = new double*[nParams]();
  double **NumeratorSum = new double*[nParams]();
  for (int i = 0; i < nParams; ++i) {
    DenomSum[i] = new double[nLags];
    NumeratorSum[i] = new double[nLags];
  }
  LagKPlots = new TH1D*[nParams];

  // Loop over the parameters of interacts
  for (int j = 0; j < nParams; ++j) {

    // Loop over each lag
    for (int k = 0; k < nLags; ++k) {
      NumeratorSum[j][k] = 0.0;
      DenomSum[j][k] = 0.0;
    }

    // Make TH1Ds for each parameter which hold the lag
    std::string HistName = XsecNames[j] +"_"+ ParamName_v[j]+"_Lag";
    LagKPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nLags, 0.0, nLags);
    LagKPlots[j]->GetXaxis()->SetTitle("Lag");
    LagKPlots[j]->GetYaxis()->SetTitle("Auto-correlation function");
  }

  // Loop over the lags
  // Each lag is indepdent so might as well multi-thread them!
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
#pragma omp parallel for
#endif
  for (int k = 0; k < nLags; ++k) {

    // Loop over the number of entries
    for (int i = 0; i < nEntries; ++i) {

      // Loop over the number of parameters
      for (int j = 0; j < nParams; ++j) {

        double Diff = ParamValues[i][j]-ParamSums[j];

        // Only sum the numerator up to i = N-k
        if (i < nEntries-k) {
          double LagTerm = ParamValues[i+k][j]-ParamSums[j];
          double Product = Diff*LagTerm;
          NumeratorSum[j][k] += Product;
        }

        // Square the difference to form the denominator
        double Denom = Diff*Diff;
        DenomSum[j][k] += Denom;
      }
    }
  }

  OutputFile->cd();
  TDirectory *AutoCorrDir = OutputFile->mkdir("Auto_corr");
  // Now fill the LagK auto-correlation plots
  for (int j = 0; j < nParams; ++j) {
    for (int k = 0; k < nLags; ++k) {
      LagKPlots[j]->SetBinContent(k, NumeratorSum[j][k]/DenomSum[j][k]);
    }
    AutoCorrDir->cd();
    LagKPlots[j]->Write();
    delete LagKPlots[j];
  }
  delete[] LagKPlots;

  for (int i = 0; i < nEntries; ++i) {
    delete ParamValues[i];
  }
  delete[] ParamValues;

  delete ParamSums;

}


// **************************
// Batched means, literally read from an array and chuck into TH1D
void MCMCDiag_Handler::BatchedMeans() {
  // **************************

  BatchedParamPlots = new TH1D*[nParams];
  for (int j = 0; j < nParams; ++j) {
      std::string HistName = XsecNames[j] +"_"+ ParamName_v[j]+"_batch";
    BatchedParamPlots[j] = new TH1D(HistName.c_str(), HistName.c_str(), nBatches, 0, nBatches);
  }

  for (int i = 0; i < nBatches; ++i) {
    for (int j = 0; j < nParams; ++j) {
      BatchedParamPlots[j]->SetBinContent(i+1, BatchedAverages[i][j]);
      int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
      int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
      std::stringstream ss;
      ss << BatchRangeLow << " - " << BatchRangeHigh;
      BatchedParamPlots[j]->GetXaxis()->SetBinLabel(i+1, ss.str().c_str());
    }
  }

  TDirectory *BatchDir = OutputFile->mkdir("Batched_means");
  BatchDir->cd();
  for (int j = 0; j < nParams; ++j) {
    TF1 *Fitter = new TF1("Fitter","[0]", 0, nBatches);
    Fitter->SetLineColor(kRed);
    BatchedParamPlots[j]->Fit("Fitter","Rq");
    BatchedParamPlots[j]->Write();
    delete Fitter;
    delete BatchedParamPlots[j];
  }
  delete[] BatchedParamPlots;

  for (int i = 0; i < nBatches; ++i) {
    delete BatchedAverages[i];
  }

  delete[] BatchedAverages;

}

// **************************
// Acceptance Probability
void MCMCDiag_Handler::AcceptanceProbabilities() {
  // **************************

    std::cout << "Making AccProb plots..." << std::endl;

    // Set the titles and limits for TH2Ds

    AcceptanceProbPlot = new TH1D("AcceptanceProbability", "Acceptance Probability", nEntries, 0, nEntries);
    AcceptanceProbPlot->GetXaxis()->SetTitle("Step");
    AcceptanceProbPlot->GetYaxis()->SetTitle("Acceptance Probability");

    BatchedAcceptanceProblot = new TH1D("AcceptanceProbability_Batch", "AcceptanceProbability_Batch", nBatches, 0, nBatches);
    BatchedAcceptanceProblot->GetYaxis()->SetTitle("Acceptance Probability");
    
  for (int i = 0; i < nBatches; ++i) {
      BatchedAcceptanceProblot->SetBinContent(i+1, AccProbBatchedAverages[i]);
      int BatchRangeLow = double(i)*double(nEntries)/double(nBatches);
      int BatchRangeHigh = double(i+1)*double(nEntries)/double(nBatches);
      std::stringstream ss;
      ss << BatchRangeLow << " - " << BatchRangeHigh;
      BatchedAcceptanceProblot->GetXaxis()->SetBinLabel(i+1, ss.str().c_str());
  }
  
  
#ifdef MULTITHREAD
  std::cout << "Using multi-threading..." << std::endl;
#pragma omp parallel for
#endif
  for (int i = 0; i < nEntries; ++i) {
    // Set bin content for the ith bin to the parameter values
      AcceptanceProbPlot->SetBinContent(i, AccProbValues[i]);
  }
    
  TDirectory *probDir = OutputFile->mkdir("AccProb");
  probDir->cd();
  
  AcceptanceProbPlot->Write();
  BatchedAcceptanceProblot->Write();
  
  delete AcceptanceProbPlot;  
  delete BatchedAcceptanceProblot; 
}

// **************************
// Read the input covariance matrix entries
// Get stuff like parameter input errors, names, and so on
void MCMCDiag_Handler::ReadXsecCov() {
  // **************************

  CheckXsecCov();

  // Open the input covariance file to get the pre-fit error
  TFile *covFile = new TFile(XsecCov.c_str(), "OPEN");
  covFile->cd();

  // Get the array of cross-section parameter names
  TObjArray *xsec_param_names = (TObjArray*)(covFile->Get("xsec_param_names"));

  TMatrixDSym *covMatrix = (TMatrixDSym*)(covFile->Get("xsec_cov"));

  // The stuff for the xsec limits
  TVectorD* prior   = (TVectorD*)(covFile->Get("xsec_param_prior"));
  TVectorD* nominal = (TVectorD*)(covFile->Get("xsec_param_nom"));
  TVectorD* lower   = (TVectorD*)(covFile->Get("xsec_param_lb"));
  TVectorD* upper   = (TVectorD*)(covFile->Get("xsec_param_ub"));
  TMatrixT<double> *id = (TMatrixD*)(covFile->Get("xsec_param_id"));

  int nXsecPars = xsec_param_names->GetEntries();
  if (nXsecPars != nXsec) {
    std::cerr << "Number of cross-section parameters in found MCMC input file disagrees with number of branches I found in the MCMC output file" << std::endl;
    std::cerr << "Points to something wrong with the name of the input cross-section file or when writing the MCMC to file" << std::endl;
    std::cerr << "Name of file is " << XsecCov << std::endl;
    std::cout << std::endl;
    xsec_param_names->Print();
    std::cout << std::endl;
    for (int i = 0; i < nXsec; ++i) {
      int i_upd = i;
      if (nFlux > 0) {
        i_upd = nFlux + i;
      }
      std::cout << i_upd << " " << ParamName_v[i_upd] << std::endl;
    }
    throw;
  }

  // Loop over the xsec params
  for (int i = 0; i < nXsecPars; ++i) {

    std::string TempString = std::string(((TObjString*)xsec_param_names->At(i))->GetString());
    XsecNames.push_back(TempString);
    // Also get the priors etc from the input file

    // Central value and error
    double central = 0.0;
    double error = 0.0;
    double down_error = 0.0;
    double up_error = 0.0;

    // How many 1 sigma we want to move away from central
    double sigmas = 5.0;

    if ( ((*nominal)(i)) > 0 && ((*id)(i,0)) != -2) {
      central = ( (*prior)(i)) / ((*nominal)(i));
      error = sigmas * sqrt(((*covMatrix)(i,i))) / ((*nominal)(i));
    } else { 
      central = ((*prior)(i));
      error = sigmas * sqrt(((*covMatrix)(i,i)));
    }

    // We might be passed the valid range of ieter
    // Do a check to see if this is true
    if (central - error <  ((*lower)(i))) {
      down_error = (*lower)(i);
    } else {
      down_error = central - error;
    }

    if (central + error > ((*upper)(i))) {
      up_error = (*upper)(i);
    } else {
      up_error = central + error;
    }

    // Normalisation ieter are unlikely of going above 2 sigma
    if ((*id)(i,0) == -1) {
      up_error = central + 2*sqrt((*covMatrix)(i, i));
    }

    XsecLimLow.push_back(down_error);
    XsecLimHigh.push_back(up_error);
    XsecCentral.push_back(central);
  }

  covFile->Close();
  delete prior;
  delete nominal;
  delete lower;
  delete upper;
  delete id;

} // End function

// **************************
// Function to check the cross-section covariance matrix
void MCMCDiag_Handler::CheckXsecCov() {
  // **************************
  if (XsecCov.empty()) {
    std::cerr << "Have not got a name for the cross-section covariance matrix" << std::endl;
    throw;
  }
}
