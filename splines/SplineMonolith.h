#ifndef SPLINEMONOLITH_H
#define SPLINEMONOLITH_H

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>

#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <TString.h>
#include <TIterator.h>
#include <TSpline.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TStyle.h>

#ifdef MULTITHREAD
#include "omp.h"
#endif

#include "samplePDF/Structs.h"

// Make template class so we can use TF1 and TSpline3
class SMonolith {
  public:
    // Vector of TSpline3 pointers which we strip back
    SMonolith(std::vector< std::vector<TSpline3*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TSpline3_red*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TF1*> > &MasterSpline);
    SMonolith(std::vector< std::vector<TF1_red*> > &MasterSpline);
    ~SMonolith();

    // This EvalGPU should be used when using two separate x,{y,a,b,c,d} arrays to store the weights; probably the best one here!
    void EvalGPU_SepMany(float *val, bool plotWeights = false);
    // Same thing but pass parameter spline segments instead of variations
    void EvalGPU_SepMany_seg(float *val, int *segment, bool plotWeights = false);

    // Evaluate weights on the GPU
    void EvalGPU_TF1(float *val, bool plotWeight = false);

    // The returned gpu weights, read by the GPU
    float *cpu_weights;

  private:
    // Function to scan through the MasterSpline of TSpline3
    void ScanMasterSpline(std::vector<std::vector<TSpline3_red*> > &MasterSpline, unsigned int &NEvents, int &MaxPoints, int &nParams);
    // Function to scan through the MasterSpline of TF1
    void ScanMasterSpline(std::vector<std::vector<TF1_red*> > &MasterSpline, unsigned int &NEvents, int &MaxPoints, int &nParams);
    // Prepare the TSpline3_red objects for the GPU
    void PrepareForGPU(std::vector<std::vector<TSpline3_red*> > &MasterSpline);
    // Prepare the TF1_red objects for the GPU
    void PrepareForGPU(std::vector<std::vector<TF1_red*> > &MasterSpline);
    // Reduced the TSpline3 to TSpline3_red
    std::vector<std::vector<TSpline3_red*> > ReduceTSpline3(std::vector<std::vector<TSpline3*> > &MasterSpline);
    // Reduced the TF1 to TF1_red
    std::vector<std::vector<TF1_red*> > ReduceTF1(std::vector<std::vector<TF1*> > &MasterSpline);

    // This loads up coefficients into two arrays: one x array and one yabcd array
    // This should maximize our cache hits!
    inline void getSplineCoeff_SepMany(TSpline3_red* &spl, int &nPoints, float *&xArray, float *&manyArray);
    // Helper function used in the constructor, tests to see if the spline is flat
    inline bool isFlat(TSpline3_red* &spl);

    // Gets the polynomial coefficeints for TF1
    inline void getTF1Coeff(TF1_red* &spl, int &nPoints, float *&coeffs);

    // Number of NIWG parameters that have splines
    int nParams;
    // Max knots for production
    int _max_knots;
    // holds the index for good splines; don't do unsigned since starts with negative value!
    int *index;
    // Number of valid splines
    unsigned int NSplines_valid;
    // Number of total splines we can maximally have
    unsigned int NSplines_total;
    // CPU arrays to hold monolith and weights
    float *cpu_weights_var;  

    // Just some pointers to memory that doesn't get allocated so we can access the GPU
    // GPU arrays to hold monolith and weights
    float *gpu_weights;
    // GPU arrays to hold coefficients and number of points separately
    int   *gpu_nPoints_arr;
    int   *gpu_paramNo_arr;
    float *gpu_coeff_x;
    float *gpu_coeff_many;
};
#endif

