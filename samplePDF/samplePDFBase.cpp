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

void samplePDFBase::addData(std::vector< vector <double> > &data)
{
  if(nDims != 0 && nDims != 2)
  {
    std::cerr<<"You have initialised this sample with "<<nDims<<" dimensions already and now trying to set dimentison to 2"<<std::endl;
    std::cerr<<"This will not work, you can find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  nDims = 2;  
  dataSample2D = new std::vector< vector <double> >(data);
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

std::vector< vector <double> > samplePDFBase::generate2D(TH2D* pdf)
{
  std::vector< vector <double> > data;
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

double samplePDFBase::GetLikelihood()
{
  if(nDims == 0)    
  {
    std::cerr << "data sample is empty!" << std::endl;
    return -1;
  }
  double negLogL = 0;

  // for now, bin up dataset and do a binned fit
  // mc
  if(nDims == 1)
  {
    TH1D* pdf = (TH1D*)get1DHist();
    // data
    if (!dathist)
      std::cerr << "***Data histogram empty!***" << std::endl;

    /*      if(dataSample)
            {
            dathist->Reset();

            for (int i = 0; i < int(dataSample->size()); i++)
            dathist->Fill(dataSample->at(i));
            }*/

    // get likelihood
    //      std::cout << pdf->Integral() << " " << dathist->Integral() << std::endl;
    for (int i = 1; i <= pdf->GetNbinsX(); i++)
    {
      double mc = pdf->GetBinContent(i);
      double dat = dathist->GetBinContent(i);
      negLogL += getTestStatLLH(dat, mc);
    }
  }
  if(nDims == 2)
  {
    TH2D* pdf = (TH2D*)get2DHist();
    // data
    if(dataSample2D && !dathist2d)
      std::cerr << "***data histogram empty!***" << std::endl;
    /*if(dataSample2D)
      {
      dathist2d->Reset();
      for (int i = 0; i < int(dataSample2D->size()); i++)
      dathist2d->Fill(dataSample2D->at(0)[i],dataSample2D->at(1)[i]);
      }*/

    // get likelihood
    for (int i = 1; i <= pdf->GetNbinsX(); i++)
    {
      for(int j = 1; j <= pdf->GetNbinsY(); j++)
      {
        double dat = dathist2d->GetBinContent(i,j); 
        double mc = pdf->GetBinContent(i,j);
        negLogL += getTestStatLLH(dat, mc);
      }
    }
  }
  return negLogL;
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

// this function will compute the likelihood against a nominal dataset (histogram)
/*double samplePDFBase::GetLikelihoodNominal()
  {
  TH1D* pdf = get1DHist();
  double mc = pdf->GetBinContent(i,j);
  double dat = //dathist->GetBinContent(i,j);
  if(mc > 0 && dat > 0)
  {
  negLogL += (mc - dat + dat * TMath::Log(dat/mc));// + 1 / (12 * dat) + 0.5 *  TMath::Log(2*TMath::Pi() * dat));
  }
  else if(mc > 0 && dat == 0)
  negLogL += mc;

  }*/

double samplePDFBase::GetLikelihood_kernel(std::vector<double> &dataSet)
{
  // this doesnt work
  /*  std::cout << "kernel estimation likelihood" << std::endl;
      double sig = 0.5;
      double sum = 0;
      double norm = 1 / (sig * TMath::Sqrt(2*TMath::Pi()));

      for(int d = 0; d < dataSet.size(); d++)
      {
      for(int i = 0; i < skmcSamples.size(); i ++)
      for(int j = 0; j < skmcSamples[i].nEvents; j++)
      {
  //sum += TMath::Power(dataSet[d] - skmcSamples[i].rw_erec[j]);
  sum += skmcSamples[i].pot_s * skmcSamples[i].norm_s * skmcSamples[i].osc_w[j] * skmcSamples[i].flux_w[j] * skmcSamples[i].skdet_w[j] * skmcSamples[i].energyscale_w[j] * norm * TMath::Exp(-0.5 * TMath::Power((dataSet[d] - skmcSamples[i].rw_erec[j])/sig, 2));
  }
  }
  //sum /= dataSet[d];
  //  sum /= sig * sig;
  // sum *= -1 * dataSet[d] * skmcSamples[i].pot_s * skmcSamples[i].norm_s * skmcSamples[i].osc_w[j] * skmcSamples[i].flux_w[j] * skmcSamples[i].skdet_w[j] * skmcSamples[i].energyscale_w[j];
  // sum += -1 * dataSet[d] * TMath::Log(sig);

  std::cout << "finished." << std::endl;
  return -1 * TMath::Log(sum); */
  return 0;
}


// *************************
// Calculate the Barlow-Beeston likelhood contribution from MC statistics
// Assumes the beta scaling parameters are Gaussian distributed
// Follows arXiv:1103.0354 section 5 and equation 8, 9, 10, 11 on page 4/5
// Essentially solves equation 11
// data is data, mc is mc, w2 is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
double samplePDFBase::getTestStatLLH(double data, double mc, double w2) {
// *************************

  // Need some MC
  if (mc == 0) return 0.0;

  // The MC used in the likeliihood calculation
  // Is allowed to be changed by Barlow Beeston beta parameters
  double newmc = mc;

  // Not full Barlow-Beeston or what is referred to as "light": we're not introducing any more parameters
  // Assume the MC has a Gaussian distribution around generated
  // As in https://arxiv.org/abs/1103.0354 eq 10, 11

  // The penalty from MC statistics using Barlow-Beeston
  double penalty = 0;
  if (fTestStatistic == kBarlowBeeston) {
    // Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
    const double fractional = sqrt(w2)/mc;
    // -b/2a in quadratic equation
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
  }

  // Calculate the new Poisson likelihood
  // For Barlow-Beeston newmc is modified, so can only calculate Poisson likelihood after Barlow-Beeston
  // For the Poisson likelihood, this is just the usual calculation
  // For IceCube likelihood, we calculate it later
  double stat = 0;
  // All likelihood calculations may use the bare Poisson likelihood, so calculate here
  if (data == 0) stat = newmc;
  else if (newmc > 0) stat = newmc-data+data*TMath::Log(data/newmc);

  // Also try the IceCube likelihood
  // It does not modify the MC content
  // https://arxiv.org/abs/1901.04645
  // Argüelles, C.A., Schneider, A. & Yuan, T. J. High Energ. Phys. (2019) 2019: 30. https://doi.org/10.1007/JHEP06(2019)030
  // We essentially construct eq 3.16 and take the logarithm
  // in eq 3.16, mu is MC, sigma2 is w2, k is data
  if (fTestStatistic == kIceCube) {
    // If there for some reason is 0 mc uncertainty, return the Poisson LLH
    if (w2 == 0) return stat;

    // Reset the penalties if there is mc uncertainty
     stat = 0.0;
     penalty = 0.0;
     // Auxillary variables
     const long double b = mc/w2;
     const long double a = mc*b+1;
     const long double k = data;
     // Use C99's implementation of log of gamma function to not be C++11 dependent
     stat = -1*(a * logl(b) + lgammal(k+a) - lgammal(k+(long double)1) - ((k+a)*log1pl(b)) - lgammal(a));
   }

   //KS: Pearson works on assumption that event distribution in each bin is described by a Gaussian which in our case is not fulfilled for all bins, hence use it at your own risk
   if (fTestStatistic == kPearson)
   {
      stat = 0;
      stat = (data-mc)*(data-mc)/mc;
   }

   // Return the statistical contribution and penalty
   return stat+penalty;
}

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

  for(int i = 0; ModeStruct->GetNModes()+1; i++)
  {
    modeNameVect.push_back(ModeStruct->Mode_ToString(i));
  }

}
