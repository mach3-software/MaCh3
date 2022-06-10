#ifndef _samplePDFBase_h_
#define _samplePDFBase_h_

#include <iostream>
#include <vector>

#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "samplePDFInterface.h"
#include "splines/splineBase.h"
#include "Structs.h"

class samplePDFBase : public samplePDFInterface {
 public:
  samplePDFBase(){};
  samplePDFBase(double pot);

  virtual ~samplePDFBase();

  void getModeName(std::vector<std::string> &modeNameVect, bool latex=false) ;
  void getSampleName(std::vector<std::string> &sampleNameVect, bool latex=false) ;


  TH1D* get1DHist();                                               
  TH2D* get2DHist();
  TH1D* get1DDataHist(){return dathist;}
  TH2D* get2DDataHist(){return dathist2d;}
  void set1DBinning(int nbins, double* boundaries);
  void set1DBinning(int nbins, double low, double high);
  void set2DBinning(int nbins1, double* boundaries1, int nbins2, double* boundaries2);
  void set2DBinning(int nbins1, double low1, double high1, int nbins2, double low2, double high2);
  double getEventRate();
  void setMCthrow(bool mc){MCthrow=mc;}
      
  // generate fake dataset based on rejection sampling    
  vector< vector <double> > generate2D(TH2D* pdf=0);
  vector<double> generate();
  virtual double getLikelihood();
  virtual std::vector<double>* getDataSample()
    {return dataSample;};
  // nominal spectrum things
  //  double getLikelihoodNominal(); // computes the likelihood against a nominal spectra
  /*  TH1D *generateNominal1D();
  TH2D *generateNominal2D();
  TH1D *nominalSpectrum1D; 
  TH2D *nominalSpectrum2D;*/

  void addData(std::vector<double> &dat);
  void addData(std::vector< vector <double> > &dat);
  void addData(TH1D* binneddata);
  void addData(TH2D* binneddata);

  void addXsecSplines(splineBase* splines){xsecsplines=splines;}
  //virtual void whatAmI(){std::cout << "wai:samplePDFBase" << std::endl;};

  // For adding sample dependent branches to the posteriors tree
  virtual void setMCMCBranches(TTree *outtree) {};

  protected:
  void init(double pot);
  void init(double pot, std::string mc_version);
  
  //bool gpu_rw; 
  double up_bnd; // highest energy to use (MeV)

  double getLikelihood_kernel(std::vector<double> &data);


  std::vector<double>* dataSample;
  std::vector< vector <double> >* dataSample2D;
  TH1D *dathist; // tempstore for likelihood calc
  TH2D *dathist2d;    
  
  // binned PDFs
  TH1D*_hPDF1D;
  TH2D*_hPDF2D;

  //splines
  splineBase* xsecsplines;

  TRandom3* rnd;
  bool MCthrow;
      
};

#endif
