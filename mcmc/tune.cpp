#include "tune.h"

tune::tune(std::string name)
{
  TFile *ff = new TFile(name.c_str(), "READ");
  init(ff);
  delete ff;
}

tune::~tune(){}

void tune::init(TFile *file)
{
  TVectorF *_tune_pars = (TVectorF*)file->Get("step_size");
  TObjArray *_par_names = (TObjArray*)file->Get("par_names");
  n_pars = _tune_pars->GetNrows();

  tune_pars = new float[n_pars];
  par_names = new std::string[n_pars];
  
  for (int i = 0; i < n_pars; ++i)
    {
      tune_pars[i] = (*_tune_pars)(i);
      par_names[i] = std::string(((TObjString*)_par_names->At(i))->GetName());
    }
}

void tune::tuneParameters(covarianceBase& cov, char *prefix)
{
  for (int i = 0; i < cov.getSize(); ++i)
    {
      char name[64];
      sprintf(name, "%s%i", prefix, i);
      for (int j = 0; j < n_pars; ++j)
	{
	  if (strcmp(name, par_names[j].c_str()) == 0)
	    {
	      TF1 *p = new TF1(name,"TMath::Gaus(x,[0],[1])",-10,10);
	      p->SetParameters(0, tune_pars[j]);
	      cov.setPropFunct(i, p);
	      std::cout << "- Set " << name << " tuning" << std::endl;
	      break;
	    }
	  if(j == n_pars-1)
	    {
	      std::cout << "- couldn't find custom tuning for " << name << std::endl;
	    }
	}
    }
}
