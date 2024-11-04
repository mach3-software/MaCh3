#pragma once

// C++ includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

//MaCh3 includes
#include "samplePDF/Structs.h"
#include "manager/manager.h"

/// @file StatisticalUtils.h
/// @brief Utility functions for statistical interpretations in MaCh3

// **************************
/// @brief  KS: Following H. Jeffreys \cite jeffreys1998theory
/// @param BayesFactor Obtained value of Bayes factor
inline std::string GetJeffreysScale(const double BayesFactor){
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
/// @brief  KS: Based on Table 1 in https://www.t2k.org/docs/technotes/435
/// @param BayesFactor Obtained value of Bayes factor
inline std::string GetDunneKaboth(const double BayesFactor){
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
/// @brief KS: Convert sigma from normal distribution into percentage
inline double GetSigmaValue(const int sigma) {
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
/// @brief Get the Bayesian Information Criterion (BIC) or Schwarz information criterion (also SIC, SBC, SBIC)
inline double GetBIC(const double llh, const int data, const int nPars){
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
/// @brief KS: See 14.3.10 in Numerical Recipes in C
/// \cite press1992numerical
inline double GetNeffective(const int N1, const int N2) {
// ****************
  const double Nominator = (N1+N2);
  const double Denominator = (N1*N2);
  const double N_e = Nominator/Denominator;
  return N_e;
}

// ****************
/// @brief KS: For more see https://www.t2k.org/docs/technotes/429/TN429_v8#page=63
/// @param SampleNameVec vector of sample names
/// @param PValVec pvalue for each sample
/// @param Threshold pvalue accepted threshold, usually 5%
inline void CheckBonferoniCorrectedpValue(const std::vector<std::string>& SampleNameVec,
                                          const std::vector<double>& PValVec,
                                          const double Threshold = 0.05) {
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
    if(  (PValVec[i] < 0.5 && PValVec[i] < StatisticalSignificanceDown) )
    {
      MACH3LOG_INFO("Sample {} indicates disagreement between the model predictions and the data", SampleNameVec[i]);
      MACH3LOG_INFO("Bonferroni-corrected statistical significance level: {:.2f} p-value: {:.2f}", StatisticalSignificanceDown, PValVec[i]);
      Counter++;
    }
    else if( (PValVec[i] > 0.5 && PValVec[i] > StatisticalSignificanceUp) )
    {
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
/// @param CumulativeData Value of CDF for data
/// @param CumulativeMC Value of CDF for MC
/// @param CumulativeJoint Value of CDF for joint data and MC distribution
inline double GetAndersonDarlingTestStat(const double CumulativeData, const double CumulativeMC, const double CumulativeJoint) {
// ****************
  double ADstat = std::fabs(CumulativeData - CumulativeMC)/ std::sqrt(CumulativeJoint*(1 - CumulativeJoint));

  if( std::isinf(ADstat) || std::isnan(ADstat)) return 0;
  return ADstat;
}

// ****************
/// @brief KS: https://esjeevanand.uccollege.edu.in/wp-content/uploads/sites/114/2020/08/NON-PARAMTERIC-TEST-6.pdf
inline int GetNumberOfRuns(const std::vector<int>& GroupClasifier) {
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
/// @brief KS: Calculate Beta parameter which will be different based on specified test statistic
/// @param data Number of data events in a bin
/// @param mc Number of MC events in a bin
/// @param w2 Value of weight squared in a bin
/// @param TestStat Test statistic based on which we calculate beta
inline double GetBetaParameter(const double data, const double mc, const double w2, TestStatistic TestStat) {
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
/// @brief Based on \cite roberts2009adaptive
/// @param EigenValues Eigen values of covariance matrix
inline double GetSubOptimality(const std::vector<double>& EigenValues, const int TotalTarameters) {
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
/// @brief CW: Get Arithmetic mean from posterior
/// @param hist histograms from which we extract arithmetic mean
/// @param Mean Arithmetic Mean value
/// @param Error Arithmetic Error value
inline void GetArithmetic(TH1D * const hist, double& Mean, double& Error) {
// **************************
  Mean = hist->GetMean();
  Error = hist->GetRMS();
}

// **************************
/// @brief CW: Fit Gaussian to posterior
/// @param hist histograms to which we fit gaussian
/// @param gauss tf1 with gaussian, we pass pointer to make things faster
/// @param Mean Gaussian Mean value
/// @param Error Gaussian Error value
inline void GetGaussian(TH1D*& hist, TF1* gauss, double& Mean, double& Error) {
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
/// @brief Get Highest Posterior Density (HPD)
/// @param hist histograms from which we HPD
/// @param Mean HPD Mean value
/// @param Error HPD Error value
/// @param Error_p HPD Negative (left hand side) Error value
/// @param Error_m HPD Positive (right hand side) Error value
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
inline void GetHPD(TH1D* const hist, double& Mean, double& Error, double& Error_p, double& Error_m, const double coverage = 0.6827) {
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
/// @brief KS: Get 1D histogram within credible interval, hpost_copy has to have the same binning, I don't do Copy() as this will lead to problems if this is used under multithreading
/// @param hist histograms based on which we calculate credible interval
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
inline void GetCredibleInterval(TH1D* const hist, TH1D* hpost_copy, const double coverage = 0.6827) {
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
/// @brief KS: Set 2D contour within some coverage
/// @param hist2D histograms based on which we calculate credible regions
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
inline void GetCredibleRegion(TH2D* const hist2D, const double coverage = 0.6827) {
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

// *********************
/// @brief Interquartile Range (IQR)
/// @param hist histograms from which we IQR
inline double GetIQR(TH1D *Hist) {
// *********************
  if(Hist->Integral() == 0) return 0.0;

  constexpr double quartiles_x[3] = {0.25, 0.5, 0.75};
  double quartiles[3];

  Hist->GetQuantiles(3, quartiles, quartiles_x);

  return quartiles[2] - quartiles[0];
}

// ********************
/// @brief Compute the Kullback-Leibler divergence between two TH2Poly histograms.
///
/// @param DataPoly Pointer to the data histogram (TH2Poly).
/// @param PolyMC Pointer to the Monte Carlo histogram (TH2Poly).
/// @return The Kullback-Leibler divergence value. Returns 0 if the data or MC integral is zero.
inline double ComputeKLDivergence(TH2Poly* DataPoly, TH2Poly* PolyMC) {
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
/// @brief KS: Combine p-values using Fisher's method.
///
/// @param pvalues A vector of individual p-values to combine.
/// @return The combined p-value, representing the overall significance.
inline double FisherCombinedPValue(const std::vector<double>& pvalues) {
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
