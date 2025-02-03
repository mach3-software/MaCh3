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
#include "samplePDF/HistogramUtils.h"

/// @file StatisticalUtils.h
/// @brief Utility functions for statistical interpretations in MaCh3
/// @author Kamil Skwarczynski

/// @brief  KS: Following H. Jeffreys \cite jeffreys1998theory
/// @param BayesFactor Obtained value of Bayes factor
std::string GetJeffreysScale(const double BayesFactor);

/// @brief  KS: Based on Table 1 in https://www.t2k.org/docs/technotes/435
/// @param BayesFactor Obtained value of Bayes factor
std::string GetDunneKaboth(const double BayesFactor);

/// @brief KS: Convert sigma from normal distribution into percentage
double GetSigmaValue(const int sigma);

/// @brief Get the Bayesian Information Criterion (BIC) or Schwarz information criterion (also SIC, SBC, SBIC)
double GetBIC(const double llh, const int data, const int nPars);

/// @brief KS: See 14.3.10 in Numerical Recipes in C
/// \cite press1992numerical
double GetNeffective(const int N1, const int N2);

/// @brief KS: For more see https://www.t2k.org/docs/technotes/429/TN429_v8#page=63
/// @param SampleNameVec vector of sample names
/// @param PValVec pvalue for each sample
/// @param Threshold pvalue accepted threshold, usually 5%
void CheckBonferoniCorrectedpValue(const std::vector<std::string>& SampleNameVec,
                                          const std::vector<double>& PValVec,
                                          const double Threshold = 0.05);

/// @param CumulativeData Value of CDF for data
/// @param CumulativeMC Value of CDF for MC
/// @param CumulativeJoint Value of CDF for joint data and MC distribution
double GetAndersonDarlingTestStat(const double CumulativeData, const double CumulativeMC, const double CumulativeJoint);

/// @brief KS: https://esjeevanand.uccollege.edu.in/wp-content/uploads/sites/114/2020/08/NON-PARAMTERIC-TEST-6.pdf
int GetNumberOfRuns(const std::vector<int>& GroupClasifier);

/// @brief KS: Calculate Beta parameter which will be different based on specified test statistic
/// @param data Number of data events in a bin
/// @param mc Number of MC events in a bin
/// @param w2 Value of weight squared in a bin
/// @param TestStat Test statistic based on which we calculate beta
double GetBetaParameter(const double data, const double mc, const double w2, TestStatistic TestStat);

/// @brief Based on \cite roberts2009adaptive
/// @param EigenValues Eigen values of covariance matrix
double GetSubOptimality(const std::vector<double>& EigenValues, const int TotalTarameters);

/// @brief CW: Get Arithmetic mean from posterior
/// @param hist histograms from which we extract arithmetic mean
/// @param Mean Arithmetic Mean value
/// @param Error Arithmetic Error value
void GetArithmetic(TH1D * const hist, double& Mean, double& Error);

/// @brief CW: Fit Gaussian to posterior
/// @param hist histograms to which we fit gaussian
/// @param gauss tf1 with gaussian, we pass pointer to make things faster
/// @param Mean Gaussian Mean value
/// @param Error Gaussian Error value
void GetGaussian(TH1D*& hist, TF1* gauss, double& Mean, double& Error);

/// @brief Get Highest Posterior Density (HPD)
/// @param hist histograms from which we HPD
/// @param Mean HPD Mean value
/// @param Error HPD Error value
/// @param Error_p HPD Negative (left hand side) Error value
/// @param Error_m HPD Positive (right hand side) Error value
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
void GetHPD(TH1D* const hist, double& Mean, double& Error, double& Error_p, double& Error_m, const double coverage = 0.6827);

/// @brief KS: Get 1D histogram within credible interval, hpost_copy has to have the same binning, I don't do Copy() as this will lead to problems if this is used under multithreading
/// @param hist histograms based on which we calculate credible interval
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
void GetCredibleInterval(TH1D* const hist, TH1D* hpost_copy, const double coverage = 0.6827);

/// @brief KS: Get 1D histogram within credible interval, hpost_copy has to have the same binning, I don't do Copy() as this will lead to problems if this is used under multithreading
/// @param hist histograms based on which we calculate credible interval
/// @param CredibleInSigmas Whether interval is in sigmas or percentage
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
void GetCredibleIntervalSig(TH1D* const hist, TH1D* hpost_copy, const bool CredibleInSigmas, const double coverage = 0.6827);

/// @brief KS: Set 2D contour within some coverage
/// @param hist2D histograms based on which we calculate credible regions
/// @param CredibleInSigmas Whether interval is in sigmas or percentage
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
void GetCredibleRegion(TH2D* const hist2D, const double coverage = 0.6827);

/// @brief KS: Set 2D contour within some coverage
/// @param hist2D histograms based on which we calculate credible regions
/// @param coverage What is defined coverage, by default 0.6827 (1 sigma)
void GetCredibleRegionSig(TH2D* const hist2D, const bool CredibleInSigmas, const double coverage = 0.6827);

/// @brief Interquartile Range (IQR)
/// @param hist histograms from which we IQR
double GetIQR(TH1D *Hist);

/// @brief Compute the Kullback-Leibler divergence between two TH2Poly histograms.
///
/// @param DataPoly Pointer to the data histogram (TH2Poly).
/// @param PolyMC Pointer to the Monte Carlo histogram (TH2Poly).
/// @return The Kullback-Leibler divergence value. Returns 0 if the data or MC integral is zero.
double ComputeKLDivergence(TH2Poly* DataPoly, TH2Poly* PolyMC);

/// @brief KS: Combine p-values using Fisher's method.
///
/// @param pvalues A vector of individual p-values to combine.
/// @return The combined p-value, representing the overall significance.
double FisherCombinedPValue(const std::vector<double>& pvalues);

/// @brief Thin MCMC Chain, to save space and maintain low autocorrelations.
///
/// @param FilePath Path to MCMC chain you want to thin
/// @param ThinningCut every which entry you want to thin
/// @cite 2011ThinningMCMC
/// @warning Thinning is done over entry not steps, it may now work very well for merged chains
void ThinningMCMC(const std::string& FilePath, const int ThinningCut);

/// @brief Compute the Z-score for a given value.
///
/// The Z-score indicates how many standard deviations a value is from the mean.
/// A positive Z-score means the value is above the mean, while a negative Z-score
/// means it is below the mean.
///
/// @param value The data point for which to compute the Z-score.
/// @param mean The mean of the data set.
/// @param stddev The standard deviation of the data set. Must be non-zero.
/// @return The Z-score of the given value.
/// @warning Ensure that stddev is not zero to avoid division by zero.
double GetZScore(const double value, const double mean, const double stddev);

/// @brief Compute the P-value from a given Z-score.
///
/// The P-value represents the probability of observing a value as extreme as
/// the given Z-score under the null hypothesis of a standard normal distribution.
///
/// @param zScore The Z-score for which to compute the P-value.
/// @return The P-value corresponding to the given Z-score.
double GetPValueFromZScore(const double zScore);

/// @brief Get the mode error from a TH1D
/// @param hpost hist from which we extract mode error
double GetModeError(TH1D* hpost);
