#ifndef SPLINEMONOLITH_H
#define SPLINEMONOLITH_H

// C++ includes
#include <cstdlib>
#include <iomanip>

// ROOT include
#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TString.h"
#include "TIterator.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TStyle.h"

#ifdef MULTITHREAD
#include "omp.h"
#endif

// MaCh3  includes
#include "samplePDF/Structs.h"

//KS:We store ccoefficeints {y,b,c,d} in one array one by one, this is only to define it once rather then insert "4" all over the code
#define _nCoeff_ 4
//KS: For TF1 we store at most 5 coefficeints, we could make it more flexible but for now define it here to make future changes easier to track
#define _nTF1Coeff_ 5

// Make template class so we can use TF1 and TSpline3
class SMonolith {
  public:
    // Vector of TSpline3 pointers which we strip back
    SMonolith(std::vector< std::vector<TSpline3*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TSpline3_red*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TF1*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TF1_red*> > &MasterSpline);
    ~SMonolith();

    void setSplinePointers(std::vector< const double* > spline_ParsPointers) {
      splineParsPointer = spline_ParsPointers;
      for (__int__ i = 0; i < nParams; ++i) SplineInfoArray[i].splineParsPointer = spline_ParsPointers[i];
    };
    // This Eval should be used when using two separate x,{y,a,b,c,d} arrays to store the weights; probably the best one here!
    // Same thing but pass parameter spline segments instead of variations
    void Evaluate();

    // Evaluate weights on the CPU/GPU
    void Evaluate_TF1();

    // The returned gpu weights, read by the GPU
    float* cpu_weights;
    //KS: This holds the total CPU weights that gets read in samplePDFND
    float *cpu_total_weights;
    //KS: Get pointer to total weight to make fit faster wrooom!   
    const float* retPointer(const int event) {return &cpu_total_weights[event];}
  private:
    //KS: Set everything to null etc.
    inline void Initialise();
    // Function to scan through the MasterSpline of TSpline3
    void ScanMasterSpline(std::vector<std::vector<TSpline3_red*> > &MasterSpline, unsigned int &NEvents, int &MaxPoints, short int &nParams, int &nSplines, unsigned int &nKnots);
    // Function to scan through the MasterSpline of TF1
    void ScanMasterSpline(std::vector<std::vector<TF1_red*> > &MasterSpline, unsigned int &NEvents, int &MaxPoints, short int &nParams);
    // Prepare the TSpline3_red objects for the GPU
    void PrepareForGPU(std::vector<std::vector<TSpline3_red*> > &MasterSpline);
    inline void PrepareForGPU_TSpline3();

    // Array of FastSplineInfo structs: keeps information on each xsec spline for fast evaluation
    // Method identical to TSpline3::Eval(double) but faster because less operations
    FastSplineInfo *SplineInfoArray;
    inline void FindSplineSegment();
    short int *segments;
    float *vals;

    std::vector< const double* > splineParsPointer;
    
    // Prepare the TF1_red objects for the GPU
    void PrepareForGPU(std::vector<std::vector<TF1_red*> > &MasterSpline);
    // Reduced the TSpline3 to TSpline3_red
    std::vector<std::vector<TSpline3_red*> > ReduceTSpline3(std::vector<std::vector<TSpline3*> > &MasterSpline);
    // Reduced the TF1 to TF1_red
    std::vector<std::vector<TF1_red*> > ReduceTF1(std::vector<std::vector<TF1*> > &MasterSpline);
    inline void PrepareForGPU_TF1();

    //CPU based code which eval each spline
    inline void CalcSplineWeights();
    //Same but TF1
    inline void CalcSplineWeights_TF1();
    //Calc total weight
    inline void ModifyWeights();
    //Conversion from valid splines to all
    inline void ModifyWeights_GPU();

    // This loads up coefficients into two arrays: one x array and one yabcd array 
    // This should maximize our cache hits!
    inline void getSplineCoeff_SepMany(TSpline3_red* &spl, int &nPoints, float *&xArray, float *&manyArray);
    // Helper function used in the constructor, tests to see if the spline is flat
    inline bool isFlat(TSpline3_red* &spl);

    // Gets the polynomial coefficeints for TF1
    inline void getTF1Coeff(TF1_red* &spl, int &nPoints, float *&coeffs);

    //Number of events
    unsigned int NEvents;
    // Number of NIWG parameters that have splines
    short int nParams;
    // Max knots for production
    int _max_knots;
    // holds the index for good splines; don't do unsigned since starts with negative value!
    int *index_cpu;

    // Number of valid splines
    unsigned int NSplines_valid;
    // Number of total splines we can maximally have, if each event had the maximum number of splines found across all events
    unsigned int NSplines_total;

    // Number of total splines if each event had every parameter's spline
    unsigned int NSplines_total_large;

    unsigned int nKnots;
    
    // Just some pointers to memory that doesn't get allocated so we can access the GPU
    // GPU arrays to hold monolith and weights
    float *gpu_weights;

    //KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
    float *gpu_total_weights;
    unsigned int *gpu_nParamPerEvent;
    unsigned int *cpu_nParamPerEvent;

    // CPU arrays to hold monolith and weights
    float *cpu_weights_var;

    // GPU arrays to hold number of points
    short int *cpu_nPoints_arr;
    short int *gpu_nPoints_arr;
    //KS: Consider merging paramNo and nKnots into one consequitive array
    short int *cpu_paramNo_arr;
    short int *gpu_paramNo_arr;
    //KS: Number of knots per spline
    unsigned int *cpu_nKnots_arr;
    unsigned int *gpu_nKnots_arr;
    //KS: GPU arrays to hold X coefficient
    float *cpu_coeff_x;
    float *gpu_coeff_x;
    // GPU arrays to hold other coefficients
    float *cpu_coeff_many;
    float *gpu_coeff_many;
};
#endif
