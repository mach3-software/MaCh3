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

      // If there are no data, test-stat shall return only MC*beta
      double stat = mc*beta;
       // With data, test-stat shall return LogL for newMC*beta which includes the low MC bound
      if (data > 0) {
        newmc *= beta;
        double stat = newmc - data + data * std::log(data/newmc);
      }

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
      
      // If w2 == 0 for any reason, return Poisson LogL
      if (w2 == 0) return getTestStatLLH(data,mc);
      
      // The MC can be changed
      double newmc = mc;
      
      // If MC falls below the low MC bound, use low MC bound for newmc
      if (mc < M3::_LOW_MC_BOUND_) newmc = M3::_LOW_MC_BOUND_;
      
      //the so-called effective count
      const double k = newmc*newmc / w2;
      //Calculate beta which is scaling factor between true and generated MC
      const double beta = (data + k) / (newmc + k);
      
      newmc *= beta;
      // And penalise the movement in beta relative the mc uncertainty
      const double penalty = k * beta - k + k * std::log( k / ( k * beta ) );

      // If there are no data, this shall return newmc
      double stat = newmc;
      // All likelihood calculations may use the bare Poisson likelihood, so calculate here
      // Only if there are some data
      if(data > 0) stat = newmc - data + data * std::log( data / newmc );

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
      // IceCube low MC bound is implemented to return Poisson(data, _LOW_MC_BOUND_)
      // up until the IceCube(data, mc) test-statistic is less than Poisson(data, _LOW_MC_BOUND_)

      // If there is 0 MC uncertainty, or the MC is less than low MC bound while having some data
      // => Return Poisson(data, mc)
      if ( w2 == 0 || ( mc < M3::_LOW_MC_BOUND_ && data > 0 ) ) return getTestStatLLH(data,mc);
      
      // Auxiliary variables
      const long double b = mc/w2;
      const long double a = mc*b+1;
      
      // Use C99's implementation of log of gamma function to not be C++11 dependent
      const double stat = double(-1*(a * logl(b) + lgammal(data+a) - lgammal(data+1) - ((data+a)*log1pl(b)) - lgammal(a)));

      // Check whether the stat is more than Poisson-like bound for low mc (mc < data)
      const double poisson = getTestStatLLH_test(data, M3::_LOW_MC_BOUND_);
      if(mc < data && stat > poisson) return poisson;

      // Otherwise, return IceCube test stat
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
