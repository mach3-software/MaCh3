//MaCh3 includes
#include "mcmc/StatisticalUtils.h"

// **************************
std::string GetJeffreysScale(const double BayesFactor){
// **************************
  std::string JeffreysScale = "";
  //KS: Get fancy Jeffreys Scale as I am to lazy to look into table every time
  if(BayesFactor < 0)        JeffreysScale = "Negative";
  else if( 5 > BayesFactor)  JeffreysScale = "Barely worth mentioning";
  else if( 10 > BayesFactor) JeffreysScale = "Substantial";
  else if( 15 > BayesFactor) JeffreysScale = "Strong";
  else if( 20 > BayesFactor) JeffreysScale = "Very strong";
  else JeffreysScale = "Decisive";

  return JeffreysScale;
}

// **************************
std::string GetDunneKaboth(const double BayesFactor){
// **************************
  std::string DunneKaboth = "";
  //KS: Get fancy DunneKaboth Scale as I am to lazy to look into table every time

  if(2.125 > BayesFactor)         DunneKaboth = "< 1 #sigma";
  else if( 20.74 > BayesFactor)   DunneKaboth = "> 1 #sigma";
  else if( 369.4 > BayesFactor)   DunneKaboth = "> 2 #sigma";
  else if( 15800 > BayesFactor)   DunneKaboth = "> 3 #sigma";
  else if( 1745000 > BayesFactor) DunneKaboth = "> 4 #sigma";
  else DunneKaboth = "> 5 #sigma";

  return DunneKaboth;
}

// *********************
double GetSigmaValue(const int sigma) {
// *********************
  double width = 0;
  switch (std::abs(sigma))
  {
    case 1:
      width = 0.682689492137;
      break;
    case 2:
      width = 0.954499736104;
      break;
    case 3:
      width = 0.997300203937;
      break;
    case 4:
      width = 0.999936657516;
      break;
    case 5:
      width = 0.999999426697;
      break;
    case 6:
      width = 0.999999998027;
      break;
    default:
      MACH3LOG_ERROR("{}  is unsupported value of sigma", sigma);
      throw MaCh3Exception(__FILE__ , __LINE__ );
      break;
  }
  return width;
}

// ****************
double GetBIC(const double llh, const int data, const int nPars){
// ****************
  if(nPars == 0)
  {
    MACH3LOG_ERROR("You haven't passed number of model parameters as it is still zero");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  const double BIC = double(nPars * logl(data) + llh);

  return BIC;
}

// ****************
double GetNeffective(const int N1, const int N2) {
// ****************
  const double Nominator = (N1+N2);
  const double Denominator = (N1*N2);
  const double N_e = Nominator/Denominator;
  return N_e;
}

// ****************
void CheckBonferoniCorrectedpValue(const std::vector<std::string>& SampleNameVec,
                                          const std::vector<double>& PValVec,
                                          const double Threshold) {
// ****************
  MACH3LOG_INFO("");
  if(SampleNameVec.size() != PValVec.size())
  {
    MACH3LOG_ERROR("Size of vectors do not match");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  const size_t NumberOfStatisticalTests = SampleNameVec.size();
  //KS: 0.05 or 5% is value used by T2K.
  const double StatisticalSignificanceDown = Threshold / double(NumberOfStatisticalTests);
  const double StatisticalSignificanceUp = 1 - StatisticalSignificanceDown;
  MACH3LOG_INFO("Bonferroni-corrected statistical significance level: {:.2f}", StatisticalSignificanceDown);

  int Counter = 0;
  for(unsigned int i = 0; i < SampleNameVec.size(); i++)
  {
    if(  (PValVec[i] < 0.5 && PValVec[i] < StatisticalSignificanceDown) ) {
      MACH3LOG_INFO("Sample {} indicates disagreement between the model predictions and the data", SampleNameVec[i]);
      MACH3LOG_INFO("Bonferroni-corrected statistical significance level: {:.2f} p-value: {:.2f}", StatisticalSignificanceDown, PValVec[i]);
      Counter++;
    } else if( (PValVec[i] > 0.5 && PValVec[i] > StatisticalSignificanceUp) ) {
      MACH3LOG_INFO("Sample {} indicates disagreement between the model predictions and the data", SampleNameVec[i]);
      MACH3LOG_INFO("Bonferroni-corrected statistical significance level: {:.2f} p-value: {:.2f}", StatisticalSignificanceUp, PValVec[i]);
      Counter++;
    }
  }
  if(Counter == 0) {
    MACH3LOG_INFO("Every sample passed Bonferroni-corrected statistical significance level test");
  } else {
    MACH3LOG_WARN("{} samples didn't pass Bonferroni-corrected statistical significance level test", Counter);
  }
  MACH3LOG_INFO("");
}

// ****************
double GetAndersonDarlingTestStat(const double CumulativeData, const double CumulativeMC, const double CumulativeJoint) {
// ****************
  double ADstat = std::fabs(CumulativeData - CumulativeMC)/ std::sqrt(CumulativeJoint*(1 - CumulativeJoint));

  if( std::isinf(ADstat) || std::isnan(ADstat)) return 0;
  return ADstat;
}

// ****************
int GetNumberOfRuns(const std::vector<int>& GroupClasifier) {
// ****************
  int NumberOfRuns = 0;
  int PreviousGroup = -999;

  //KS: If group changed increment run
  for (unsigned int i = 0; i < GroupClasifier.size(); i++)
  {
    if(GroupClasifier[i] != PreviousGroup)
      NumberOfRuns++;
    PreviousGroup = GroupClasifier[i];
  }

  return NumberOfRuns;
}

// ****************
double GetBetaParameter(const double data, const double mc, const double w2, TestStatistic TestStat) {
// ****************
  double Beta = 0.0;

  if (TestStat == kDembinskiAbdelmottele) {
    //the so-called effective count
    const double k = mc*mc / w2;
    //Calculate beta which is scaling factor between true and generated MC
    Beta = (data + k) / (mc + k);
  }
  //KS: Below is technically only true for Cowan's BB, which will not be true for Poisson or IceCube, because why not...
  else {
    // CW: Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
    const double fractional = std::sqrt(w2)/mc;
    // CW: -b/2a in quadratic equation
    const double temp = mc*fractional*fractional-1;
    // CW: b^2 - 4ac in quadratic equation
    const double temp2 = temp*temp + 4*data*fractional*fractional;
    if (temp2 < 0) {
      MACH3LOG_ERROR("Negative square root in Barlow Beeston coefficient calculation!");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    // CW: Solve for the positive beta
    Beta = (-1*temp+std::sqrt(temp2))/2.;
  }
  return Beta;
}


// *********************
double GetSubOptimality(const std::vector<double>& EigenValues, const int TotalTarameters) {
// *********************
  double sum_eigenvalues_squared_inv = 0.0;
  double sum_eigenvalues_inv = 0.0;
  for (unsigned int j = 0; j < EigenValues.size(); j++)
  {
    //KS: IF Eigen values are super small skip them
    //if(EigenValues[j] < 0.0000001) continue;
    sum_eigenvalues_squared_inv += std::pow(EigenValues[j], -2);
    sum_eigenvalues_inv += 1.0 / EigenValues[j];
  }
  const double SubOptimality = TotalTarameters * sum_eigenvalues_squared_inv / std::pow(sum_eigenvalues_inv, 2);
  return SubOptimality;
}


// **************************
void GetArithmetic(TH1D * const hist, double& Mean, double& Error) {
// **************************
  Mean = hist->GetMean();
  Error = hist->GetRMS();
}

// **************************
void GetGaussian(TH1D*& hist, TF1* gauss, double& Mean, double& Error) {
// **************************
  const double meanval = hist->GetMean();
  const double err = hist->GetRMS();
  const double peakval = hist->GetBinCenter(hist->GetMaximumBin());

  // Set the range for the Gaussian fit
  gauss->SetRange(meanval - 1.5*err , meanval + 1.5*err);
  // Set the starting parameters close to RMS and peaks of the histograms
  gauss->SetParameters(hist->GetMaximum()*err*std::sqrt(2*3.14), peakval, err);

  // Perform the fit
  hist->Fit(gauss->GetName(),"Rq");
  hist->SetStats(0);

  Mean = gauss->GetParameter(1);
  Error = gauss->GetParameter(2);
}

// ***************
void GetHPD(TH1D* const hist, double& Mean, double& Error, double& Error_p, double& Error_m, const double coverage) {
// ****************
  // Get the bin which has the largest posterior density
  const int MaxBin = hist->GetMaximumBin();
  // And it's value
  const double peakval = hist->GetBinCenter(MaxBin);

  // The total integral of the posterior
  const long double Integral = hist->Integral();
  //KS: and integral of left handed and right handed parts
  const long double LowIntegral = hist->Integral(1, MaxBin-1) + hist->GetBinContent(MaxBin)/2.0;
  const long double HighIntegral = hist->Integral(MaxBin+1, hist->GetNbinsX()) + hist->GetBinContent(MaxBin)/2.0;

  // Keep count of how much area we're covering
  //KS: Take only half content of HPD bin as one half goes for right handed error and the other for left handed error
  long double sum = hist->GetBinContent(MaxBin)/2.0;

  // Counter for current bin
  int CurrBin = MaxBin;
  while (sum/HighIntegral < coverage && CurrBin < hist->GetNbinsX()) {
    CurrBin++;
    sum += hist->GetBinContent(CurrBin);
  }
  const double sigma_p = std::fabs(hist->GetBinCenter(MaxBin)-hist->GetXaxis()->GetBinUpEdge(CurrBin));
  // Reset the sum
  //KS: Take only half content of HPD bin as one half goes for right handed error and the other for left handed error
  sum = hist->GetBinContent(MaxBin)/2.0;

  // Reset the bin counter
  CurrBin = MaxBin;
  // Counter for current bin
  while (sum/LowIntegral < coverage && CurrBin > 1) {
    CurrBin--;
    sum += hist->GetBinContent(CurrBin);
  }
  const double sigma_m = std::fabs(hist->GetBinCenter(CurrBin)-hist->GetBinLowEdge(MaxBin));

  // Now do the double sided HPD
  //KS: Start sum from the HPD
  sum = hist->GetBinContent(MaxBin);
  int LowBin = MaxBin;
  int HighBin = MaxBin;
  long double LowCon = 0.0;
  long double HighCon = 0.0;

  while (sum/Integral < coverage && (LowBin > 0 || HighBin < hist->GetNbinsX()+1))
  {
    LowCon = 0.0;
    HighCon = 0.0;
    //KS:: Move further only if you haven't reached histogram end
    if(LowBin > 1)
    {
      LowBin--;
      LowCon = hist->GetBinContent(LowBin);
    }
    if(HighBin < hist->GetNbinsX())
    {
      HighBin++;
      HighCon = hist->GetBinContent(HighBin);
    }

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/Integral > coverage && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/Integral > coverage && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }
  }

  double sigma_hpd = 0.0;
  if (LowCon > HighCon) {
    sigma_hpd = std::fabs(hist->GetBinLowEdge(LowBin)-hist->GetBinCenter(MaxBin));
  } else {
    sigma_hpd = std::fabs(hist->GetXaxis()->GetBinUpEdge(HighBin)-hist->GetBinCenter(MaxBin));
  }

  Mean = peakval;
  Error = sigma_hpd;
  Error_p = sigma_p;
  Error_m = sigma_m;
}

// ***************
void GetCredibleInterval(TH1D* const hist, TH1D* hpost_copy, const double coverage) {
// ***************
  if(coverage > 1)
  {
    MACH3LOG_ERROR("Specified Credible Interval is greater that 1 and equal to {} Should be between 0 and 1", coverage);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  //KS: Reset first copy of histogram
  hpost_copy->Reset("");
  hpost_copy->Fill(0.0, 0.0);

  //KS: Temporary structure to be thread save
  std::vector<double> hist_copy(hist->GetXaxis()->GetNbins()+1);
  std::vector<bool> hist_copy_fill(hist->GetXaxis()->GetNbins()+1);
  for (int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
  {
    hist_copy[i] = hist->GetBinContent(i);
    hist_copy_fill[i] = false;
  }

  /// Loop over histogram bins with highest number of entries until covered 90 or 68.3%
  const long double Integral = hist->Integral();
  long double sum = 0;

  while ((sum / Integral) < coverage)
  {
    /// Get bin of highest content and save the number of entries reached so far
    int max_entry_bin = 0;
    double max_entries = 0.;
    for (int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
    {
      if (hist_copy[i] > max_entries)
      {
        max_entries = hist_copy[i];
        max_entry_bin = i;
      }
    }
    /// Replace bin value by -1 so it is not looped over as being maximum bin again
    hist_copy[max_entry_bin] = -1.;
    hist_copy_fill[max_entry_bin] = true;

    sum += max_entries;
  }
  //KS: Now fill our copy only for bins which got included in coverage region
  for(int i = 0; i <= hist->GetXaxis()->GetNbins(); ++i)
  {
    if(hist_copy_fill[i]) hpost_copy->SetBinContent(i, hist->GetBinContent(i));
  }
}

// ***************
void GetCredibleIntervalSig(TH1D* const hist, TH1D* hpost_copy, const bool CredibleInSigmas, const double coverage) {
// ***************
  //KS: Slightly different approach depending if intervals are in percentage or sigmas
  if(CredibleInSigmas) {
    //KS: Convert sigmas into percentage
    const double CredReg = GetSigmaValue(int(std::round(coverage)));
    GetCredibleInterval(hist, hpost_copy, CredReg);
  } else {
    GetCredibleInterval(hist, hpost_copy, coverage);
  }
}

// ***************
void GetCredibleRegion(TH2D* const hist2D, const double coverage) {
// ***************
  if(coverage > 1)
  {
    MACH3LOG_ERROR("Specified Credible Region is greater than 1 and equal to {:.2f} Should be between 0 and 1 {}", coverage);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }

  //KS: Temporary structure to be thread save
  std::vector<std::vector<double>> hist_copy(hist2D->GetXaxis()->GetNbins()+1,
                                             std::vector<double>(hist2D->GetYaxis()->GetNbins()+1));
  for (int i = 0; i <= hist2D->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j <= hist2D->GetYaxis()->GetNbins(); ++j) {
      hist_copy[i][j] = hist2D->GetBinContent(i, j);
    }
  }

  /// Loop over histogram bins with highest number of entries until covered 90 or 68.3%
  const long double Integral = hist2D->Integral();
  long double sum = 0;

  //We need to as ROOT requires array to set to contour
  double Contour[1];
  while ((sum / Integral) < coverage)
  {
    /// Get bin of highest content and save the number of entries reached so far
    int max_entry_bin_x = 0;
    int max_entry_bin_y = 0;
    double max_entries = 0.;
    for (int i = 0; i <= hist2D->GetXaxis()->GetNbins(); ++i)
    {
      for (int j = 0; j <= hist2D->GetYaxis()->GetNbins(); ++j)
      {
        if (hist_copy[i][j] > max_entries)
        {
          max_entries = hist_copy[i][j];
          max_entry_bin_x = i;
          max_entry_bin_y = j;
        }
      }
    }
    /// Replace bin value by -1 so it is not looped over as being maximum bin again
    hist_copy[max_entry_bin_x][max_entry_bin_y] = -1.;

    sum += max_entries;
    Contour[0] = max_entries;
  }
  hist2D->SetContour(1, Contour);
}

// ***************
void GetCredibleRegionSig(TH2D* const hist2D, const bool CredibleInSigmas, const double coverage) {
// ***************
  if(CredibleInSigmas) {
    //KS: Convert sigmas into percentage
    const double CredReg = GetSigmaValue(int(std::round(coverage)));
    GetCredibleRegion(hist2D, CredReg);
  } else {
    GetCredibleRegion(hist2D, coverage);
  }
}

// *********************
double GetIQR(TH1D *Hist) {
// *********************
  if(Hist->Integral() == 0) return 0.0;

  constexpr double quartiles_x[3] = {0.25, 0.5, 0.75};
  double quartiles[3];

  Hist->GetQuantiles(3, quartiles, quartiles_x);

  return quartiles[2] - quartiles[0];
}

// ********************
double ComputeKLDivergence(TH2Poly* DataPoly, TH2Poly* PolyMC) {
// *********************
  double klDivergence = 0.0;
  double DataIntegral = NoOverflowIntegral(DataPoly);
  double MCIntegral = NoOverflowIntegral(PolyMC);
  for (int i = 1; i < DataPoly->GetNumberOfBins()+1; ++i)
  {
    if (DataPoly->GetBinContent(i) > 0 && PolyMC->GetBinContent(i) > 0) {
      klDivergence += DataPoly->GetBinContent(i) / DataIntegral *
      std::log((DataPoly->GetBinContent(i) / DataIntegral) / ( PolyMC->GetBinContent(i) / MCIntegral));
    }
  }
  return klDivergence;
}
// ********************
double FisherCombinedPValue(const std::vector<double>& pvalues) {
// ********************
  double testStatistic = 0;
  for(size_t i = 0; i < pvalues.size(); i++)
  {
    const double pval = std::max(0.00001, pvalues[i]);
    testStatistic += -2.0 * std::log(pval);
  }
  // Degrees of freedom is twice the number of p-values
  int degreesOfFreedom = int(2 * pvalues.size());
  double pValue = TMath::Prob(testStatistic, degreesOfFreedom);

  return pValue;
}

// ********************
void ThinningMCMC(const std::string& FilePath, const int ThinningCut) {
// ********************
  // Define the path for the temporary thinned file
  std::string TempFilePath = "Thinned_" + FilePath;
  int ret = system(("cp " + FilePath + " " + TempFilePath).c_str());
  if (ret != 0) {
    MACH3LOG_WARN("Error: system call to copy file failed with code {}", ret);
  }

  TFile *inFile = TFile::Open(TempFilePath.c_str(), "UPDATE");
  if (!inFile || inFile->IsZombie()) {
    MACH3LOG_ERROR("Error opening file: {}", TempFilePath);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TTree *inTree = inFile->Get<TTree>("posteriors");
  if (!inTree) {
    MACH3LOG_ERROR("Error: TTree 'posteriors' not found in file.");
    inFile->ls();
    inFile->Close();
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Clone the structure without data
  TTree *outTree = inTree->CloneTree(0);

  // Loop over entries and apply thinning
  Long64_t nEntries = inTree->GetEntries();
  double retainedPercentage = (double(nEntries) / ThinningCut) / double(nEntries) * 100;
  MACH3LOG_INFO("Thinning will retain {:.2f}% of chains", retainedPercentage);
  for (Long64_t i = 0; i < nEntries; i++) {
    if (i % (nEntries/10) == 0) {
      MaCh3Utils::PrintProgressBar(i, nEntries);
    }
    if (i % ThinningCut == 0) {
      inTree->GetEntry(i);
      outTree->Fill();
    }
  }
  inFile->WriteTObject(outTree, "posteriors", "kOverwrite");
  inFile->Close();
  delete inFile;

  MACH3LOG_INFO("Thinned TTree saved and overwrote original in: {}", TempFilePath);
}

// ********************
double GetZScore(const double value, const double mean, const double stddev) {
// ********************
  return (value - mean) / stddev;
}

// ********************
double GetPValueFromZScore(const double zScore) {
// ********************
  return 0.5 * std::erfc(-zScore / std::sqrt(2));
}

// ****************
// Get the mode error from a TH1D
double GetModeError(TH1D* hpost) {
// ****************
  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();

  // The total integral of the posterior
  const double Integral = hpost->Integral();

  int LowBin = MaxBin;
  int HighBin = MaxBin;
  double sum = hpost->GetBinContent(MaxBin);;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum/Integral < 0.6827 && (LowBin > 0 || HighBin < hpost->GetNbinsX()+1) )
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
