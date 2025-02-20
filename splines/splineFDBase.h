#pragma once

//MaCh3 includes
#include "samplePDF/Structs.h"
#include "splines/SplineBase.h"
#include "manager/MaCh3Modes.h"

_MaCh3_Safe_Include_Start_ //{
// ROOT includes
#include "TH3F.h"
_MaCh3_Safe_Include_End_ //}

/// @brief Bin-by-bin class calculating response for spline parameters.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/05.-Splines).
/// @author Dan Barrow
/// @author Ed Atkin
/// @author Henry Wallace
class splineFDBase : public SplineBase {
  /// @todo ETA - do all of these functions and members actually need to be public?
  public:
	/// @brief Constructor
  splineFDBase(covarianceXsec *xsec_, MaCh3Modes *Modes_);
	/// @brief Destructor
	/// @todo it need some love
	virtual ~splineFDBase();

	/// @brief CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays 
	/// to store the weights; probably the best one here! Same thing but pass parameter 
	/// spline segments instead of variations
	void Evaluate();

	/// @brief add oscillation channel to spline monolith
	bool AddSample(const std::string& SampleName,
				   const std::string& DetID,
				   const std::vector<std::string>& OscChanFileNames,
				   const std::vector<std::string>& SplineVarNames);
	void TransferToMonolith();
	/// @brief Remove setup variables not needed for spline evaluations
	void cleanUpMemory();

	//Have to define this in your own class 
        virtual void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames);
	/// @brief Check if there are any repeated modes. This is used to reduce the number
	/// of modes in case many interaction modes get averaged into one spline
        std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector);
	/// @brief Return the splines which affect a given event
        std::vector< std::vector<int> > GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val);

	/// @brief
	std::vector<TAxis*> FindSplineBinning(std::string FileName, std::string SampleName);

	int CountNumberOfLoadedSplines(bool NonFlat=false, int Verbosity=0);
	int getNDim(int BinningOpt);
	std::string getDimLabel(int BinningOpt, unsigned int Axis);
	int getSampleIndex(const std::string& SampleName);
	bool isValidSplineIndex(const std::string& SampleName, int iSyst, int iOscChan, int iMode, int iVar1, int iVar2, int iVar3);

	void BuildSampleIndexingArray(const std::string& SampleName);
	void PrepForReweight();
	void getSplineCoeff_SepMany(int splineindex, M3::float_t *& xArray, M3::float_t *&manyArray);
	void PrintBinning(TAxis* Axis);
	void PrintSampleDetails(const std::string& SampleName);
	void PrintArrayDetails(const std::string& SampleName);
	
	const M3::float_t* retPointer(const int sample, const int oscchan, const int syst, const int mode,
									const int var1bin, const int var2bin, const int var3bin) const{
	  int index = indexvec[sample][oscchan][syst][mode][var1bin][var2bin][var3bin];
	  return &weightvec_Monolith[index];
	}

  protected:
	/// @brief CPU based code which eval weight for each spline
	void CalcSplineWeights() override;
	/// @brief Calc total event weight, not used by Bin-by-bin splines
	void ModifyWeights() override {return;};
	/// Pointer to covariance xsec
	covarianceXsec* xsec;

	//And now the actual member variables	
	std::vector<std::string> SampleNames;
	std::vector<int> Dimensions;
	std::vector<std::vector<std::string>> DimensionLabels;
	std::vector<std::string> DetIDs;
	std::vector<int> nSplineParams;
	std::vector<int> nOscChans;

	std::vector< std::vector< std::vector<TAxis*> > > SplineBinning;
	std::vector< std::vector<std::string> > SplineFileParPrefixNames;
	/// A vector of vectors of the spline modes that a systematic applies to
	/// This gets compared against the event mode to figure out if a syst should
	/// apply to an event or not
	std::vector< std::vector< std::vector<int> > > SplineModeVecs;
	/// @brief This holds the global spline index and is used to grab the current parameter value
	/// to evaluate splines at. Each internal vector will be of size of the number of spline
	/// systematics which affect that sample.
	std::vector< std::vector<int> > GlobalSystIndex;
    /// @brief spline interpolation types for each sample. These vectors are from
    /// a call to GetSplineInterpolationFromDetID()
    std::vector< std::vector<SplineInterpolation> > SplineInterpolationTypes;

	std::vector<std::string> UniqueSystNames;
	std::vector<int> UniqueSystIndices;
	std::vector<int> UniqueSystCurrSegment;

	/// @brief Variables related to determined which modes have splines and which piggy-back of other modes
	std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< int > > > > > > > indexvec;
	std::vector<int > coeffindexvec;
	std::vector<int> uniquecoeffindices; //Unique coefficient indices

	std::vector< TSpline3_red* > splinevec_Monolith;

	int MonolithSize;
	int MonolithIndex;
	int CoeffIndex;

	/// Need to keep track of which splines are flat and which aren't
	bool *isflatarray;
	/// x coefficients for each spline
	M3::float_t *xcoeff_arr;
	/// ybcd coefficients for each spline
	M3::float_t *manycoeff_arr;

	std::vector<M3::float_t> weightvec_Monolith;
	std::vector<int> uniquesplinevec_Monolith;

  MaCh3Modes* Modes;
  enum TokenOrdering{kSystToken,kModeToken,kVar1BinToken,kVar2BinToken,kVar3BinToken,kNTokens};
  virtual std::vector<std::string> GetTokensFromSplineName(std::string FullSplineName) = 0;
};
