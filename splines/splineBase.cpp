#include <iostream>

#include "splineBase.h"

splineBase::splineBase(const char *name, int ntype)
{
  //splinefile = new TFile(name, "READ");
  nutype = ntype;
  //  setupSplines();

}

splineBase::~splineBase()
{

}
