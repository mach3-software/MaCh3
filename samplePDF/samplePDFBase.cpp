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
}

void samplePDFBase::addData(std::vector<double> &data)
{
  nDims = 1;
  dataSample = new std::vector<double>(data);
  if(dathist == NULL)
  {
      std::cerr<<"dathist not initialised"<<std::endl;
      std::cerr <<"Find me here "<< __FILE__ << ":" << __LINE__ << std::endl;
      throw;
  }
  dathist->Reset();                                                         
  for (int i = 0; i < int(dataSample->size()); i++)                         
    dathist->Fill(dataSample->at(i));
}

void samplePDFBase::addData(std::vector< vector <double> > &data)
{
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
  nDims = 1;
  std::cout << "adding 1D data histogram : " << binneddata -> GetName() << " with " << binneddata->Integral() << " events." << std::endl;
  //KS: If exist delete to avoid memory leak
  if(dathist != NULL) delete dathist;
  dathist = binneddata;
}

void samplePDFBase::addData(TH2D* binneddata)
{
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
  //double upp = up_bnd; // upper energy bound for sampling
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

double samplePDFBase::getLikelihood()
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
double getTestStatLLH(double data, double mc) {
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
/*double samplePDFBase::getLikelihoodNominal()
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

double samplePDFBase::getLikelihood_kernel(std::vector<double> &dataSet)
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


void samplePDFBase::getSampleName(std::vector<std::string> &sampleNameVect, bool latex) 
{
  if(sampleNameVect.size() !=0)
    sampleNameVect.clear() ;

  if(!latex) {
    sampleNameVect.push_back("beam nu_mu") ;
    sampleNameVect.push_back("beam nu_e") ;
    sampleNameVect.push_back("beam antinu_mu") ;
    sampleNameVect.push_back("beam antinu_e") ;
    sampleNameVect.push_back("oscillated nu_e") ;
    sampleNameVect.push_back("oscillated antinu_e") ;
  }
  else{
    sampleNameVect.push_back("beam $\\nu_{\\mu}$") ;
    sampleNameVect.push_back("beam $\\nu_{e}$") ;
    sampleNameVect.push_back("beam $\\bar{\\nu_{\\mu}}$") ;
    sampleNameVect.push_back("beam $\\bar{\\nu_{e}}$") ;
    sampleNameVect.push_back("oscillated $\\nu_{e}$") ;
    sampleNameVect.push_back("oscillated $\\bar{\\nu_{e}}$") ;
  }
}


void samplePDFBase::getModeName(std::vector<std::string> &modeNameVect, bool latex) 
{
  if(modeNameVect.size() !=0)
    modeNameVect.clear() ;

  if(!latex) {
    modeNameVect.push_back("CCQE") ;
    modeNameVect.push_back("CC-1pi") ;
    modeNameVect.push_back("CC-coherent") ;
    modeNameVect.push_back("CC-Npi") ;
    modeNameVect.push_back("CC-other") ;
    modeNameVect.push_back("NC pi0") ;
    modeNameVect.push_back("NC-1pi +/-") ;
    modeNameVect.push_back("NC-coherent") ;
    modeNameVect.push_back("NC-other") ;
    modeNameVect.push_back("2p2h") ;
    modeNameVect.push_back("NC-1gamma") ;
  }
  else{
    modeNameVect.push_back("CCQE");
    modeNameVect.push_back("CC 1$\\pi$");
    modeNameVect.push_back("CC coherent");
    modeNameVect.push_back("CC n$\\pi$");
    modeNameVect.push_back("CC other");
    modeNameVect.push_back("NC $\\pi^{0}$");
    modeNameVect.push_back("NC $\\pi^{+/-}$");
    modeNameVect.push_back("NC coherent");
    modeNameVect.push_back("NC other");
    modeNameVect.push_back("2p-2h");
    modeNameVect.push_back("NC 1$\\gamma$");
  }
}
