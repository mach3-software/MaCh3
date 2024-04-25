#pragma once

#include "yaml-cpp/yaml.h"
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <TVectorT.h>
#include <TVectorD.h>
#include <TMatrixT.h>
#include <TList.h>
#include <TObjString.h>
#include <TObjArray.h>
#include "samplePDF/Structs.h"

class covarianceParserYAML{

public:
   covarianceParserYAML();
   covarianceParserYAML(const char *filename);
   ~covarianceParserYAML();

   int GetNParameters() {return fNumPar;}
   std::vector<std::string> GetParameterNames() {return fNames;}
   
   TVectorT<double>* GetPrior(){return fPrior;}
   TVectorT<double>* GetGenerated(){return fGenerated;}
   TVectorT<double>* GetLowerBound(){return fLB;}
   TVectorT<double>* GetUpperBound(){return fUB;}
   TMatrixT<double>* GetID(){return fID;}
   TVectorT<double>* GetStepScale(){return fStepScale;}
   TVectorT<double>* GetFlatPrior(){return fFlatPrior;}
   TMatrixT<double>* GetCovMatrix(){return fCovMatrix;}


   TObjArray* GetFDSplineNames(){return fFD_spline_names;}
   TObjArray* GetFDSplineModes(){return fFD_spline_modes;}

   TObjArray* GetNDSplineNames(){return fND_spline_names;}
   
   TObjArray* GetNormModes() {return fNorm_modes;}
   TObjArray* GetNormHornCurrent() {return fNorm_horncurrents;}
   TObjArray* GetNormElements() {return fNorm_elem;}
   TObjArray* GetNormNuPDG() {return fNorm_nupdg;}
   TObjArray* GetNormProdNuPDG() {return fNorm_preoscnupdg;}
   
   TObjArray* GetKinematicPars() {return fKinematicPars;}
   TObjArray* GetKinematicBounds() {return fKinematicBounds;}

private:
   std::vector<std::string> fNames;
   int fNumPar;
   YAML::Node fYAMLDoc;

   TVectorT<double> *fPrior;
   TVectorT<double> *fGenerated;
   TVectorT<double> *fLB;
   TVectorT<double> *fUB;
   TMatrixT<double> *fID;
   TVectorT<double> *fStepScale;
   TVectorT<double> *fFlatPrior;
   TMatrixT<double> *fCovMatrix;

   TObjArray* fND_spline_names;
   

   TObjArray* fFD_spline_names;
   TObjArray* fFD_spline_modes;

   TObjArray* fNorm_modes;
   TObjArray* fNorm_horncurrents;
   TObjArray* fNorm_elem;
   TObjArray* fNorm_nupdg;
   TObjArray* fNorm_preoscnupdg;

   TObjArray* fKinematicPars;
   TObjArray* fKinematicBounds;

};

