#pragma once

#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <vector>

#include "samplePDF/Structs.h"

class samplePDFInterface
{
 public:
  virtual void reweight(double *oscpar)=0;

  virtual double getEventRate()=0;
 
  virtual std::vector< std::vector <double> > generate2D(TH2D* pdf)=0;
  virtual std::vector<double> generate()=0;

  virtual double GetLikelihood()=0;

  virtual void addData(std::vector<double> &dat)=0;
  virtual void addData(std::vector< std::vector <double> > &dat)=0;
  virtual void addData(TH1D* binneddata)=0;
  virtual void addData(TH2D* binneddata)=0;

  virtual void printPosteriors()=0; // save images of posterior histos
  
 protected:
  virtual void init(double pot)=0;
  virtual void init(double pot, std::string mc_version)=0;
  virtual double getCovLikelihood()=0;
  virtual double GetLikelihood_kernel(std::vector<double> &data)=0;
  virtual void fill1DHist()=0;
  virtual void fill2DHist()=0;
};
