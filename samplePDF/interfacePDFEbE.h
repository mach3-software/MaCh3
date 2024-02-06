#ifndef _interfacePDFEbE_h_
#define _interfacePDFEbE_h_
#include <vector>

class interfacePDFEbE 
{
 public:
  virtual double GetEventWeight(int sample, int event)=0;
  virtual double getDiscVar(int sample , int event, int varindx)=0;
  virtual int getNMCSamples()=0;
  virtual int getNEventsInSample(int sample)=0;
};

#endif
