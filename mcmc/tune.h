#pragma once

#include <iostream>

#include "TFile.h"
#include "TVectorF.h"
#include "TObjArray.h"

#include "covariance/covarianceBase.h"
#include "covariance/covarianceOsc.h"

class tune 
{
 public:
  tune(std::string file);
  virtual ~tune();

  void tuneParameters(covarianceBase& cov_obj, char *prefix);
  
 private:
  void init(TFile *file);
  TFile *tune_file;
  float *tune_pars;
  std::string *par_names;
  int n_pars;
};
