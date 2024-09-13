#pragma once

//ROOT
#include "TH3F.h"

//MaCh3
#include "splines/SplineBase.h"

/// @brief Bin-by-bin class calculating response for spline parameters.
/// @see For more details, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/05.-Splines).
class splineFDBase : public SplineBase {
  /// @todo ETA - do all of these functions and members actually need to be public?
  public:
    /// @brief Constructor
    splineFDBase(covarianceXsec *xsec_ = NULL);
    /// @brief Destructor
    /// @todo it need some love
    virtual ~splineFDBase();

    /// @brief  CW: This Eval should be used when using two separate x,{y,a,b,c,d} arrays to store the weights; probably the best one here! Same thing but pass parameter spline segments instead of variations
    void Evaluate();

    //Spline Monolith things
    //Essential methods used externally
    /// @todo Move these to splineFDBase in core
    bool AddSample(std::string SampleName, int BinningOpt, int DetID, std::vector<std::string> OscChanFileNames, std::vector<std::string> SplineVarNames);
    void TransferToMonolith();
    void cleanUpMemory();

	//Have to define this in your own class 
	virtual void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames)=0;
	virtual std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector)=0;
	virtual std::vector< std::vector<int> > GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val)=0;

	std::vector<TAxis*> FindSplineBinning(std::string FileName, std::string SampleName);

	int CountNumberOfLoadedSplines(bool NonFlat=false, int Verbosity=0);
	std::string getDimLabel(int iSample, unsigned int Axis);
	std::vector<std::vector<std::string>> DimensionLabels;
    
	int getSampleIndex(std::string SampleName);
	bool isValidSplineIndex(std::string SampleName, int iSyst, int iOscChan, int iMode, int iVar1, int iVar2, int iVar3);

	void BuildSampleIndexingArray(std::string SampleName);
	void PrepForReweight();
	void getSplineCoeff_SepMany(int splineindex, _float_ *& xArray, _float_ *&manyArray);
	void PrintBinning(TAxis* Axis);
	void PrintSampleDetails(std::string SampleName);
	void PrintArrayDetails(std::string SampleName);
	void PrintArrayDimension();

	const double* retPointer(int sample, int oscchan, int syst, int mode, int var1bin, int var2bin, int var3bin){
	  int index = indexvec[sample][oscchan][syst][mode][var1bin][var2bin][var3bin];
	  return &weightvec_Monolith[index];
	}

  protected:
    /// @brief CW:Code used in step by step reweighting, Find Spline Segment for each param
    inline void FindSplineSegment() override;
    /// @brief CPU based code which eval weight for each spline
    inline void CalcSplineWeights() override;
    /// @brief Calc total event weight, not used by Bin-by-bin splines
    inline void ModifyWeights() override {return;};
    /// Pointer to covariance xsec
    covarianceXsec* xsec;

	//And now the actual member variables	
	std::vector<std::string> SampleNames;
	std::vector<int> BinningOpts;
	std::vector<int> Dimensions;
	std::vector<int> DetIDs;
	std::vector<int> nSplineParams;
	std::vector<int> nOscChans;

    /// This holds the global spline index and is used to grab the current parameter value
    /// to evaluate splines at. Each internal vector will be of size of the number of spline
    /// systematics which affect that sample.
    std::vector< std::vector<int> > GlobalSystIndex;
	std::vector< std::vector< std::vector<TAxis*> > > SplineBinning;
	/// 
	std::vector< std::vector<std::string> > SplineFileParPrefixNames;
    /// A vector of vectors of the spline modes that a systematic applies to
    /// This gets compared against the event mode to figure out if a syst should
    /// apply to an event or not
    std::vector< std::vector< std::vector<int> > > SplineModeVecs;

	int nUniqueSysts;
	std::vector<std::string> UniqueSystNames;
	std::vector<int> UniqueSystIndices;
	std::vector<int> UniqueSystNKnots;
	std::vector<int> UniqueSystCurrSegment;
	std::vector< std::vector<_float_> > UniqueSystXPts;

	/// DB Variables related to determined which modes have splines and which piggy-back of other modes
	std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< int > > > > > > > indexvec;
	std::vector<int > coeffindexvec;
    /// Unique coefficient indices
	std::vector<int>uniquecoeffindices;

    std::vector<TSpline3_red*> splinevec_Monolith;

    int MonolithSize;
    int MonolithIndex;
    int CoeffIndex;

    //Probably need to clear these arrays up at some point
    _float_ *xVarArray;
    /// Need to keep track of which splines are flat and which aren't
    bool *isflatarray;
    /// x coefficients for each spline
    _float_ *xcoeff_arr;
    /// ybcd coefficients for each spline
    _float_ *manycoeff_arr;

    std::vector<double> weightvec_Monolith;
    std::vector<int> uniquesplinevec_Monolith;
};
