#include "SampleHandlerBase.h"

// ***************************************************************************
SampleHandlerBase::SampleHandlerBase() {
// ***************************************************************************
  nEvents = 0;
  nSamples = 0;
}

// ***************************************************************************
SampleHandlerBase::~SampleHandlerBase() {
// ***************************************************************************

}

// ***************************************************************************
// Poisson likelihood calc for data and MC event rates
double SampleHandlerBase::GetTestStatLLH(const double data, const double mc) const {
// ***************************************************************************
  // Return MC if there are no data, returns 0 for data == 0 && mc == 0  
  if(data == 0) return mc;

  // If there are some data, but the prediction falls below the MC bound => return Poisson LogL for the low MC bound
  if(mc < M3::_LOW_MC_BOUND_) return (M3::_LOW_MC_BOUND_ - data + data * std::log(data/M3::_LOW_MC_BOUND_));

  // Otherwise, just return usual Poisson LogL using Stirling's approximation
  // http://hyperphysics.phy-astr.gsu.edu/hbase/math/stirling.html
  return (mc - data + data * std::log(data/mc));
}

// *************************
// data is data, mc is mc, w2 is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
double SampleHandlerBase::GetTestStatLLH(const double data, const double mc, const double w2) const {
// *************************
  switch (fTestStatistic)
  {
    //CW: Not full Barlow-Beeston or what is referred to as "light": we're not introducing any more parameters
    // Assume the MC has a Gaussian distribution around generated
    // As in https://arxiv.org/abs/1103.0354 eq 10, 11
    //CW: Calculate the Barlow-Beeston likelihood contribution from MC statistics
    // Assumes the beta scaling parameters are Gaussian distributed
    // Follows arXiv:1103.0354 section 5 and equation 8, 9, 10, 11 on page 4/5
    // Essentially solves equation 11
    case (kBarlowBeeston):
    {
      // The MC used in the likelihood calculation is allowed to be changed by Barlow Beeston beta parameters
      double newmc = mc;

      // If MC falls below the low MC bound, use low MC bound for newmc
      if (mc < M3::_LOW_MC_BOUND_) newmc = M3::_LOW_MC_BOUND_;
        
      // Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
      const double fractional = std::sqrt(w2)/newmc;
      // fractional^2 to avoid doing same operation several times
      const double fractional2 = fractional*fractional;
      // b in quadratic equation
      const double temp = newmc*fractional2-1;
      // b^2 - 4ac in quadratic equation
      const double temp2 = temp*temp + 4*data*fractional2;
      if (temp2 < 0) {
        MACH3LOG_ERROR("Negative square root in Barlow Beeston coefficient calculation!");
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
      // Solve for the positive beta
      const double beta = (-1*temp+sqrt(temp2))/2.;

      // With data, test-stat shall return LogL for newMC*beta which has the low MC bound
      newmc *= beta;
      double stat = newmc - data + data * std::log(data/newmc);
      // If there are no data, test-stat shall return only MC*beta
      if (data == 0) stat = mc*beta;

      // Now, MC stat penalty
      // The penalty from MC statistics using Conways approach (https://cds.cern.ch/record/1333496?)
      double penalty = 0;
      if (fractional > 0) penalty = (beta-1)*(beta-1)/(2*fractional*fractional);

      // Returns test-stat plus the MC stat penalty 
      return stat+penalty;
    }
    break;
    //KS: Alternative calculation of Barlow-Beeston following Hans Dembinski and Ahmed Abdelmottele arXiv:2206.12346v2
    case (kDembinskiAbdelmotteleb):
    {
      //KS: code follows authors implementation from:
      //https://github.com/scikit-hep/iminuit/blob/059d06b00cae097ebf340b218b4eb57357111df8/src/iminuit/cost.py#L274-L300

      //the so-called effective count
      const double k = mc*mc / w2;
      //Calculate beta which is scaling factor between true and generated MC
      const double beta = (data + k) / (mc + k);

      newmc = mc*beta;
      // And penalise the movement in beta relative the mc uncertainty
      const double penalty = k*beta-k+k*std::log(k/(k*beta));

      // Calculate the new Poisson likelihood
      // For Barlow-Beeston newmc is modified, so can only calculate Poisson likelihood after Barlow-Beeston
      // For the Poisson likelihood, this is just the usual calculation
      // For IceCube likelihood, we calculate it later
      double stat = 0;
      // All likelihood calculations may use the bare Poisson likelihood, so calculate here
      if (data == 0) stat = newmc;
      else if (newmc > 0) stat = newmc-data+data*std::log(data/newmc);

      // Return the statistical contribution and penalty
      return stat+penalty;
    }
    break;
    //CW: Also try the IceCube likelihood
    // It does not modify the MC content
    // https://arxiv.org/abs/1901.04645
    // Argüelles, C.A., Schneider, A. & Yuan, T. J. High Energ. Phys. (2019) 2019: 30. https://doi.org/10.1007/JHEP06(2019)030
    // We essentially construct eq 3.16 and take the logarithm
    // in eq 3.16, mu is MC, sigma2 is w2, k is data
    case (kIceCube):
    {
      double stat = 0.0;

      // If there for some reason is 0 mc uncertainty, return the Poisson LLH
      if (w2 == 0)
      {
        // Calculate the new Poisson likelihood
        if (data == 0) stat = newmc;
        else if (newmc > 0) stat = newmc-data+data*std::log(data/newmc);

        return stat;
      }
      // Auxiliary variables
      const long double b = mc/w2;
      const long double a = mc*b+1;
      const long double k = data;
      // Use C99's implementation of log of gamma function to not be C++11 dependent
      stat = double(-1*(a * logl(b) + lgammal(k+a) - lgammal(k+1) - ((k+a)*log1pl(b)) - lgammal(a)));

      // Return the statistical contribution and penalty
      return stat;
    }
    break;
    //KS: Pearson works on assumption that event distribution in each bin is described by a Gaussian which in our case is not fulfilled for all bins, hence use it at your own risk
    case (kPearson):
    {
      //KS: 2 is because this function returns -LLH not -2LLH
      const double stat = (data-mc)*(data-mc)/(2*mc);

      // Return the statistical
      return stat;
    }
    break;
    case (kPoisson):
    {
      //Just call getTestStatLLH which doesn't take in weights
      //and is a Poisson likelihood comparison.
      return GetTestStatLLH(data, mc);
    }
    break;
    case TestStatistic::kNTestStatistics:
      MACH3LOG_ERROR("kNTestStatistics is not a valid TestStatistic!");
      throw MaCh3Exception(__FILE__, __LINE__);
    default:
    MACH3LOG_ERROR("Couldn't find TestStatistic {} exiting!", static_cast<int>(fTestStatistic));
    throw MaCh3Exception(__FILE__ , __LINE__ );
  } // end switch
}

// ***************************************************************************
// CW: Silence cout and cerr. Last is risky but psyche persists on spamming both
void SampleHandlerBase::QuietPlease() {
// ***************************************************************************
  #if DEBUG > 0
  return;
  #else
  buf = std::cout.rdbuf();
  errbuf = std::cerr.rdbuf();
  std::cout.rdbuf( nullptr );
  std::cerr.rdbuf( nullptr );
  #endif
}

// ***************************************************************************
// CW: Reset cout and cerr
void SampleHandlerBase::NowTalk() {
// ***************************************************************************
  #if DEBUG > 0
  return;
  #else
  std::cout.rdbuf(buf);
  std::cerr.rdbuf(errbuf);
  #endif
}
