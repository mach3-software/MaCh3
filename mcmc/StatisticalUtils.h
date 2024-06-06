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

// **************************
/// @brief  KS: Following H. Jeffreys. The theory of probability. UOP Oxford, 1998. DOI: 10.2307/3619118.
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
inline std::string GetDunneKaboth(const double BayesFactor){
// **************************
  std::string DunneKaboth = "";
  //KS: Get fancy DunneKaboth Scale as I am to lazy to look into table every time

  if(2.125 > BayesFactor)    DunneKaboth = "< 1 #sigma";
  else if( 20.74 > BayesFactor)  DunneKaboth = "> 1 #sigma";
  else if( 369.4 > BayesFactor) DunneKaboth = "> 2 #sigma";
  else if( 15800 > BayesFactor) DunneKaboth = "> 3 #sigma";
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
      throw;
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
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  const double BIC = nPars * logl(data) + llh;

  return BIC;
}

// ****************
/// @brief KS: See 14.3.10 in Numerical Recipes in C
inline double GetNeffective(const int N1, const int N2) {
// ****************

  const double Nominator = (N1+N2);
  const double Denominator = (N1*N2);
  const double N_e = Nominator/Denominator;
  return N_e;
}

// ****************
/// @brief KS: For more see https://www.t2k.org/docs/technotes/429/TN429_v8#page=63
inline void CheckBonferoniCorrectedpValue(std::vector<std::string> SampleNameVec, std::vector<double> PValVec, double Threshold = 0.05) {
// ****************

  std::cout<<""<<std::endl;

  if(SampleNameVec.size() != PValVec.size())
  {
    MACH3LOG_ERROR("Size of vectors do not match");
    throw;
  }

  const int NumberOfStatisticalTests = SampleNameVec.size();
  //KS: 0.05 or 5% is value used by T2K.
  const double StatisticalSignificanceDown = Threshold / NumberOfStatisticalTests;
  const double StatisticalSignificanceUp = 1 - StatisticalSignificanceDown;
  std::cout<<"Bonferroni-corrected statistical significance level "<<StatisticalSignificanceDown<<std::endl;
  int Counter = 0;
  for(unsigned int i = 0; i < SampleNameVec.size(); i++)
  {
    if(  (PValVec[i] < 0.5 && PValVec[i] < StatisticalSignificanceDown) )
    {
      std::cout<<"Sample "<<SampleNameVec[i]<<" indicates of any disagreement between the model predictions and the data"<<std::endl;
      std::cout<<"Bonferroni-corrected statistical significance level "<<StatisticalSignificanceDown<<" p-value "<<PValVec[i]<<std::endl;
      Counter++;
    }
    else if( (PValVec[i] > 0.5 && PValVec[i] > StatisticalSignificanceUp) )
    {
      std::cout<<"Sample "<<SampleNameVec[i]<<" indicates of any disagreement between the model predictions and the data"<<std::endl;
      std::cout<<"Bonferroni-corrected statistical significance level "<<StatisticalSignificanceUp<<" p-value "<<PValVec[i]<<std::endl;
      Counter++;
    }
  }

  if(Counter == 0)
  {
    std::cout<<" Every sample passed Bonferroni-corrected statistical significance level test"<<std::endl;
  }
  else
  {
    std::cout<<Counter<<" samples didn't pass Bonferroni-corrected statistical significance level test"<<std::endl;
  }
  std::cout<<""<<std::endl;
}


// ****************
inline double GetAndersonDarlingTestStat(double CumulativeData, double CumulativeMC, double CumulativeJoint) {
// ****************

  double ADstat = std::fabs(CumulativeData - CumulativeMC)/ std::sqrt(CumulativeJoint*(1 - CumulativeJoint));

  if( std::isinf(ADstat) || std::isnan(ADstat)) return 0;
  return ADstat;
}


// ****************
/// @brief KS: https://esjeevanand.uccollege.edu.in/wp-content/uploads/sites/114/2020/08/NON-PARAMTERIC-TEST-6.pdf
inline int GetNumberOfRuns(std::vector<int> GroupClasifier) {
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
      std::cerr << "Negative square root in Barlow Beeston coefficient calculation!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
    // CW: Solve for the positive beta
    Beta = (-1*temp+std::sqrt(temp2))/2.;
  }
  return Beta;
}


// *********************
/// @brief Based on http://probability.ca/jeff/ftpdir/adaptex.pdf
inline double GetSubOptimality(std::vector<double> EigenValues, const int TotalTarameters) {
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

// *********************
// Interquartile Range (IQR)
inline double GetIQR(TH1D *Hist){
// *********************

  if(Hist->Integral() == 0) return 0.0;

  double quartiles_x[3] = {0.25, 0.5, 0.75};
  double quartiles[3];

  Hist->GetQuantiles(3, quartiles, quartiles_x);

  return quartiles[2] - quartiles[0];
}
