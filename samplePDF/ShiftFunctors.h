#pragma once

#include "samplePDF/Structs.h"

class BaseFuncPar{
  public:
	BaseFuncPar(){};
	~BaseFuncPar();
	double ScaleUncertainty = 1.0;
	void SetUncertainty(double sigma){ScaleUncertainty = sigma;}
	double *Param_pos;
	void setup();
	std::vector<double> syst_pos;
	//XsecNorms4 *XsecApplicationInfo;

	fdmc_base* blah;

	virtual void Apply() =0;
};

class SKEScale : public BaseFuncPar{
  protected:
	//Have access to everything
	//skmc_base *blah;

	void Apply(){};
};

class EnergyScale : public BaseFuncPar{
  public:
   EnergyScale(){};
   ~EnergyScale();
   void Apply(){};

  protected:
   int test;

};
