#pragma once

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
