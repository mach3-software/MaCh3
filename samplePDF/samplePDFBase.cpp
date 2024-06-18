#include "samplePDFBase.h"

samplePDFBase::samplePDFBase(double pot) 
{
  nDims = 0;
  rnd = new TRandom3(0);
  MCthrow = false;
  dathist = NULL;
  dathist2d = NULL;
  dataSample = NULL;
  dataSample2D = NULL;
}

samplePDFBase::~samplePDFBase()
{
  if(dathist != NULL) delete dathist;
  if(dathist2d != NULL) delete dathist2d;
  if(dataSample != NULL) delete dataSample;
  if(dataSample2D != NULL) delete dataSample2D;  
  delete rnd;
}

void samplePDFBase::init(double pot)
{
}

void samplePDFBase::init(double pot, std::string mc_version)
{
    
    
  //TODO KS: Need to set test stat from config file
  // Set the test-statistic
  //SetTestStatistic(static_cast<TestStatistic>(FitManager->GetMCStatLLH()));
}

void samplePDFBase::addData(std::vector<double> &data)
{
  if(nDims != 0 && nDims != 1)
  {
    std::cerr<<"You have initialised this sample with "<<nDims<<" dimensions already and now trying to set dimentison to 1"<<std::endl;
    std::cerr<<"This will not work, you can find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  nDims = 1;
  dataSample = new std::vector<double>(data);
  if(dathist == NULL)
  {
      std::cerr<<"dathist not initialised"<<std::endl;
      std::cerr<<"Find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  dathist->Reset();                                                         
  for (int i = 0; i < int(dataSample->size()); i++)                         
    dathist->Fill(dataSample->at(i));
}

void samplePDFBase::addData(std::vector< std::vector <double> > &data)
{
  if(nDims != 0 && nDims != 2)
  {
    std::cerr<<"You have initialised this sample with "<<nDims<<" dimensions already and now trying to set dimentison to 2"<<std::endl;
    std::cerr<<"This will not work, you can find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  nDims = 2;  
  dataSample2D = new std::vector< std::vector <double> >(data);
  if(dathist2d == NULL)
  {
      std::cerr<<"dathist2d not initialised"<<std::endl;
      std::cerr <<"Find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  dathist2d->Reset();                                                       
  for (int i = 0; i < int(dataSample2D->size()); i++)                       
    dathist2d->Fill(dataSample2D->at(0)[i],dataSample2D->at(1)[i]); 
}

void samplePDFBase::addData(TH1D* binneddata)
{
  if(nDims != 0 && nDims != 1)
  {
    std::cerr<<"You have initialised this sample with "<<nDims<<" dimensions already and now trying to set dimentison to 1"<<std::endl;
    std::cerr<<"This will not work, you can find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  nDims = 1;
  std::cout << "adding 1D data histogram : " << binneddata -> GetName() << " with " << binneddata->Integral() << " events." << std::endl;
  //KS: If exist delete to avoid memory leak
  if(dathist != NULL) delete dathist;
  dathist = binneddata;
}

void samplePDFBase::addData(TH2D* binneddata)
{
  if(nDims != 0 && nDims != 2)
  {
    std::cerr<<"You have initialised this sample with "<<nDims<<" dimensions already and now trying to set dimentison to 2"<<std::endl;
    std::cerr<<"This will not work, you can find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  nDims = 2;
  std::cout << "adding 2D data histogram : " << binneddata -> GetName() << " with " << binneddata->Integral() << " events." << std::endl;
  //KS: If exist delete to avoid memory leak
  if(dathist2d != NULL) delete dathist;
  dathist2d = binneddata;
}

std::vector<double> samplePDFBase::generate()
{
  std::vector<double> data;
  TH1D *pdf = (TH1D*)get1DHist();
  double evrate = getEventRate();
  int num = rnd->Poisson(evrate);
  std::cout << std::endl << "sampling " << num << " events from " << evrate << std::endl;

  // rejection sampling
  //double M = 6; // *** DO NOT HARDCODE THIS, WILL ALTER RESULTS WHEN POT IS CHANGED ***
  int count = 0;

  dathist->Reset();

  while(count < num)
  {
    /*double candidate = gRandom->Uniform(upp);
      double accProb = pdf->GetBinContent(pdf->FindBin(candidate)) / M;
      double rand = gRandom->Uniform(1);
      if (accProb >= rand)
      {
      std::cout << candidate << " " << std::flush;
      data.push_back(candidate);
      dathist->Fill(candidate);
      count++;
      }*/
    double candidate = pdf->GetRandom();
    std::cout << candidate << " " << std::flush;                                                                                                                       
    data.push_back(candidate);                                                                                                                                         
    dathist->Fill(candidate);                                                                                                                                          
    count++;
  }

  std::cout << "sampling complete" << std::endl;
  return data;
}

std::vector< std::vector <double> > samplePDFBase::generate2D(TH2D* pdf)
{
  std::vector< std::vector <double> > data;
  if(!pdf) pdf = (TH2D*)get2DHist();

  if(MCthrow)
  {
    for(int i=1; i<=pdf->GetNbinsX(); i++)
    {
      for(int j=1; j<=pdf->GetNbinsY(); j++)
      {
        pdf->SetBinContent(i,j,rnd->Gaus(pdf->GetBinContent(i,j),pdf->GetBinError(i,j)));
      }
    }
  }

  double evrate = pdf->Integral();
  int num = rnd->Poisson(evrate);
  std::cout << "sampling " << num << " events from " << evrate << std::endl;

  std::vector<double> var1;
  std::vector<double> var2;
  double x,y;

  dathist2d->Reset();

  for(int i=0; i < num; i++)
  {
    pdf->GetRandom2(x,y);
    var1.push_back(x);
    var2.push_back(y);
    dathist2d->Fill(x, y);
  }
  data.push_back(var1);
  data.push_back(var2);

  std::cout << "sampling complete " << data[0].size() << std::endl;
  return data;
}

TH1D* samplePDFBase::get1DHist()
{
  fill1DHist();
  return _hPDF1D;
}
TH2D* samplePDFBase::get2DHist()
{
  fill2DHist();
  return _hPDF2D;
}

double samplePDFBase::getEventRate()
{
  //if (_hErec == NULL) 
  fill1DHist();
  return _hPDF1D->Integral();
}

// ***************************************************************************
//KS: So far only Poisson LLH, in future Barlow-Beeston and IceCube
double samplePDFBase::getTestStatLLH(double data, double mc) {
// ***************************************************************************
    double negLogL = 0;
    if(mc == 0) mc = 1E-8;
    if(mc > 0 && data > 0)
    {
        //http://hyperphysics.phy-astr.gsu.edu/hbase/math/stirling.html
        negLogL += (mc - data + data * TMath::Log(data/mc));
    }
    else if(mc > 0 && data == 0) negLogL += mc;
    
    return negLogL; 
}


// *************************
// data is data, mc is mc, w2 is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
double samplePDFBase::getTestStatLLH(const double data, const double mc, const double w2) {
// *************************

  // Need some MC
  if (mc == 0) return 0.0;

  // The MC used in the likeliihood calculation
  // Is allowed to be changed by Barlow Beeston beta parameters
  double newmc = mc;

  switch (fTestStatistic)
  {
    //CW: Not full Barlow-Beeston or what is referred to as "light": we're not introducing any more parameters
    // Assume the MC has a Gaussian distribution around generated
    // As in https://arxiv.org/abs/1103.0354 eq 10, 11
    //CW: Calculate the Barlow-Beeston likelhood contribution from MC statistics
    // Assumes the beta scaling parameters are Gaussian distributed
    // Follows arXiv:1103.0354 section 5 and equation 8, 9, 10, 11 on page 4/5
    // Essentially solves equation 11
    case (kBarlowBeeston):
    {
      // The penalty from MC statistics using Conways approach (https://cds.cern.ch/record/1333496?)
      double penalty = 0;
      // Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
      const double fractional = std::sqrt(w2)/mc;
      // b in quadratic equation
      const double temp = mc*fractional*fractional-1;
      // b^2 - 4ac in quadratic equation
      const double temp2 = temp*temp + 4*data*fractional*fractional;
      if (temp2 < 0) {
        std::cerr << "Negative square root in Barlow Beeston coefficient calculation!" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
      // Solve for the positive beta
      const double beta = (-1*temp+sqrt(temp2))/2.;
      newmc = mc*beta;
      // And penalise the movement in beta relative the mc uncertainty
      if (fractional > 0) penalty = (beta-1)*(beta-1)/(2*fractional*fractional);
      else penalty = 0;

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
    //KS: Alterantive calcaution of Barlow-Beeston following Hans Dembinski and Ahmed Abdelmottele arXiv:2206.12346v2
    case (kDembinskiAbdelmottele):
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
      // Auxillary variables
      const long double b = mc/w2;
      const long double a = mc*b+1;
      const long double k = data;
      // Use C99's implementation of log of gamma function to not be C++11 dependent
      stat = -1*(a * logl(b) + lgammal(k+a) - lgammal(k+(long double)1) - ((k+a)*log1pl(b)) - lgammal(a));

      // Return the statistical contribution and penalty
      return stat;
    }
    break;
    //KS: Pearson works on assumption that event distribution in each bin is described by a Gaussian which in our case is not fulfilled for all bins, hence use it at your own risk
    case (kPearson):
    {
      //KS: 2 is beacuese this function returns -LLH not -2LLH
      const double stat = (data-mc)*(data-mc)/(2*mc);

      // Return the statistical
      return stat;
    }
    break;
    case (kPoisson):
    {
      double stat = 0.0;
      // All likelihood calculations may use the bare Poisson likelihood, so calculate here
      if (data == 0) stat = newmc;
      else if (newmc > 0) stat = newmc-data+data*std::log(data/newmc);

      // Return the statistical contribution and penalty
      return stat;
    }
    break;

    default:
    std::cerr << "Couldn't find TestStatistic " << fTestStatistic << " exiting!" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  } // end switch
}

/*
// **************************************************
// Helper function to set LLH type used in the fit
void samplePDFBase::SetTestStatistic(TestStatistic test_stat) {
// **************************************************
  fTestStatistic = test_stat;

  std::string name = TestStatistic_ToString((TestStatistic)test_stat);
  std::cout << "Using "<< name <<" likelihood in ND280" << std::endl;
  //if(UpdateW2) std::cout << "With updating W2" << std::endl;
  //else  std::cout << "Without updating W2" << std::endl;
}
*/

void samplePDFBase::set1DBinning(int nbins, double* boundaries)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins,boundaries);
  dathist->SetBins(nbins,boundaries);
}

void samplePDFBase::set1DBinning(int nbins, double low, double high)
{
  _hPDF1D->Reset();
  _hPDF1D->SetBins(nbins,low,high);
  dathist->SetBins(nbins,low,high);
}
void samplePDFBase::set2DBinning(int nbins1, double* boundaries1, int nbins2, double* boundaries2)
{
  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins1,boundaries1,nbins2,boundaries2);
  dathist2d->SetBins(nbins1,boundaries1,nbins2,boundaries2);
}

void samplePDFBase::set2DBinning(int nbins1, double low1, double high1, int nbins2, double low2, double high2)
{
  _hPDF2D->Reset();
  _hPDF2D->SetBins(nbins1,low1,high1,nbins2,low2,high2);
  dathist2d->SetBins(nbins1,low1,high1,nbins2,low2,high2);
}


// ***************************************************************************
//KS: Sample getter
std::string samplePDFBase::GetSampleName(int Sample) {
// ***************************************************************************

  if(Sample > nSamples)
  {
   std::cerr<<" You are asking for sample "<< Sample <<" I only have "<< nSamples<<std::endl;
   throw;
  }

  return SampleName[Sample];
}
// ***************************************************************************
void samplePDFBase::GetSampleNames(std::vector<std::string> &sampleNameVect) {
// ***************************************************************************
  if(sampleNameVect.size() !=0)
    sampleNameVect.clear() ;

  for(int i = 0; nSamples; i++)
  {
    sampleNameVect.push_back(GetSampleName(i));
  }
}

// ***************************************************************************
void samplePDFBase::GetModeName(std::vector<std::string> &modeNameVect) {
// ***************************************************************************

  if(modeNameVect.size() !=0)
    modeNameVect.clear() ;

  for(int i = 0; Modes->GetNModes()+1; i++)
  {
    modeNameVect.push_back(Modes->GetMaCh3ModeName(i));
  }

}

// ***************************************************************************
// CW: Silence cout and cerr. Last is risky but psyche persists on spamming both
void samplePDFBase::QuietPlease() {
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
void samplePDFBase::NowTalk() {
// ***************************************************************************
  #if DEBUG > 0
  return;
  #else
  std::cout.rdbuf(buf);
  std::cerr.rdbuf(errbuf);
  #endif
}
