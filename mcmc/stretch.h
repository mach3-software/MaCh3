#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include "TStopwatch.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom.h"
#include <map>

#include "samplePDF/samplePDFBase.h"
#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"


class stretch
{
   public:
      stretch(bool verbose = false);
      stretch(const char *name = "output.root", int nwalkers=100, bool verbose = false);
      
      virtual ~stretch();
      
      void addSamplePDF(samplePDFBase* sample);
      void addSystObj(covarianceBase* cov,bool isOsc=false);
      void runStretch();
      
      void setChainLength(int chainl);
      
      void printInitialState();
      
      void setOscOnly(bool yes){osc_only = yes;};
      bool getOscOnly(){return osc_only;};
      
      void setGpuFit(bool yes){gpu_fit = yes;};
      bool getGpuFit(){return gpu_fit;};
      
      void setA(double A){a=A;}
      double getA(){return a;}
      
   protected:
      void init(std::string name, bool verbose);
      
      // sample holder
      std::vector<samplePDFBase*> samples;
      
      //systematic holder
      std::vector<covarianceBase*> systematics;
      std::vector< std::vector < std::vector <double > > > currentpar; //first index is cov matrix, second is parameter, third is walker
      std::vector< std::vector < std::vector <double > > > proposedpar;
      
      int osc, oscbar;

      int chainLength; // number of steps in chain
      int nwalk; //number of walkers  
      double a; //tuning parameter
      int N; //dimensionality;
      
      // current state
      int step; // current step
      std::vector< double > logLCurr; // current likelihood
      std::vector< double > logLProp; // proposed likelihood
      
      float osc_llh; // oscillation covariance llh
      float *sample_llh; // store the llh breakdowns
      float *syst_llh; // systematic llh breakdowns
      
      double accProb; // current acceptance prob
      
      // benchmarking, file IO, debugging  etc
      TStopwatch clock;
      bool debug;
      int debugInterval; // how many steps to print to debug file
      ofstream debugFile;
      
      //random number
      TRandom3* random;      
      
      // output
      TFile *outputFile;
      TTree *outTree;
      int auto_save; // auto save every N steps
      bool save_nominal;
      
      // fit types
      bool osc_only; // don't fit nuisance parameters
      bool gpu_fit; // use the GPU where applicable
      

};
