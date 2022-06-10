#ifndef _splineInterface_h_
#define _splineInterface_h_

#include <TFile.h>
#include <TSpline.h>
#include <TAxis.h>
#include <vector>

class splineInterface
{
 public:
      //splineInterface(char *spline, int nutype)=0;
      //virtual ~splineInterface()=0;
  virtual void setupSplines()=0;

 protected:
  
};

#endif
