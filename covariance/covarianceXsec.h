#pragma once

// C++ includes
#include <math.h>
#include <map>

// ROOT includes
#include "TList.h"

// MaCh3 includes
#include "covariance/covarianceBase.h"
#include "samplePDF/Structs.h"
#include "manager/YamlHelper.h"

#include "yaml-cpp/yaml.h"

class covarianceXsec : public covarianceBase {

  public:
    covarianceXsec(std::vector<std::string> FileNames, double threshold = -1, int FirstPCAdpar = -999, int LastPCAdpar = -999);
    ~covarianceXsec();

    //ETA - trying out the yaml parsing
    inline void InitXsecFromConfig();

    // Print information about the whole object once it is set
    inline void Print();

    // General Getter functions not split by detector
    // ETA - a lot of these can go... they're just duplications from the base class.
    //ETA - just return the int of the DetID, this can be removed to do a string comp at some point.
    int GetParDetID(const int i) const { return _fDetID[i];};
    //ETA - just return a string of "spline", "norm" or "functional"
    const char*  GetParamType(const int i) const {return _fParamType[i].c_str();}

    const std::vector<SplineInterpolation>& GetSplineInterpolation() const{return _fSplineInterpolationType;}
    SplineInterpolation GetParSplineInterpolation(int i) {return _fSplineInterpolationType.at(i);}

    //DB Get spline parameters depending on given DetID
    const std::vector<std::string> GetSplineParsNamesFromDetID(int DetID);
    const std::vector<std::string> GetSplineFileParsNamesFromDetID(int DetID);
    //ETA - what does this even do?
    //const std::vector<std::string> GetFDSplineFileParsNamesFromDetID(int DetID);
    //const std::vector<std::string> GetNDSplineFileParsNamesFromDetID(int DetID);
    const std::vector< std::vector<int> > GetSplineModeVecFromDetID(int DetID);
    const std::vector<int> GetSplineParsIndexFromDetID(int DetID);
    const std::vector<int> GetGlobalSystIndexFromDetID(int DetID);
    int GetNumSplineParamsFromDetID(int DetID);

    //DB Get norm/func parameters depending on given DetID
    const std::vector<XsecNorms4> GetNormParsFromDetID(int DetID);
    void SetupNormPars();
    int GetNumFuncParamsFromDetID(int DetID);
    const std::vector<std::string> GetFuncParsNamesFromDetID(int DetID);
    const std::vector<int> GetFuncParsIndexFromDetID(int DetID);

    //KS: For most covariances nominal and fparInit (prior) are the same, however for Xsec those can be different
    // For example Sigma Var are done around nominal in ND280, no idea why though...
    std::vector<double> getNominalArray() override
    {
      std::vector<double> nominal;
      for (int i = 0; i < size; i++)
      {
        nominal.push_back(_fPreFitValue.at(i));
      }
      return nominal;
    }
    inline double getNominal(const int i) override { return _fPreFitValue.at(i); };
   
    bool IsParFlux(const int i){ return isFlux[i]; }

    //KS Function to set to nominal either flux or xsec parameters
    void setXsecOnlyParameters();
    void setFluxOnlyParameters();
    
  protected:
    void initParams(const double fScale);
    void setXsecParNames();
    std::vector<bool> isFlux;

    std::vector<int> _fDetID;
    //std::vector<std::string> _fDetString;
    std::vector<std::string> _fParamType;

    //Some "usual" variables. Don't think we really need the ND/FD split
    std::vector<std::vector<int>> _fNormModes;
    std::vector<std::vector<int>> _fTargetNuclei;
    std::vector<std::vector<int>> _fNeutrinoFlavour;
    std::vector<std::vector<int>> _fNeutrinoFlavourUnosc;

    //Variables related to spline systematics
    //std::vector<std::string> _fSplineNames;
    std::vector<std::string> _fSplineNames;
    std::vector<std::vector<int>> _fSplineModes;
    std::vector<SplineInterpolation> _fSplineInterpolationType;
    std::map<int, int> _fSplineToSystIndexMap;

    //Information to be able to apply generic cuts
    std::vector<std::vector<std::string>> _fKinematicPars;
    std::vector<std::vector<std::vector<double>>> _fKinematicBounds;

    //Vector containing info for normalisation systematics
    std::vector<XsecNorms4> NormParams;
};
