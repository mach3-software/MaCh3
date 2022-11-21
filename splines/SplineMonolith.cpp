#include "SplineMonolith.h"

// All the functions that are declared extern here exist in splines/gpuSplineUtils.cu
// They are the CUDA code which we use to do the GPU processing
// For older style (Rich/Asher era) see git commit history previous to 27 Nov 2017

extern void InitGPU_SepMany(
    float **gpu_coeff_x, 
    float **gpu_coeff_many, 
    float **gpu_weights, 

    short int** gpu_paramNo_arr,
    unsigned int** gpu_nKnots_arr,
#ifndef Weight_On_SplineBySpline_Basis
    float **cpu_total_weights, 
    float **gpu_total_weights, 
    int n_events,
    unsigned int** gpu_nParamPerEvent,
#endif    
    unsigned int coeff_array_size,
    int n_splines,
    int Eve_size);

extern void InitGPU_TF1(
    float **gpu_coeff_many, 
    short int **gpu_paramNo_arr,
    short int **gpu_nPoints_arr,
    float **gpu_weights, 

#ifndef Weight_On_SplineBySpline_Basis
    float **cpu_total_weights,
    float **gpu_total_weights, 
    int n_events,
    unsigned int** gpu_nParamPerEvent,
#endif 

    int n_splines);

extern void CopyToGPU_SepMany(
    short int *gpu_paramNo_arr,
    unsigned int *gpu_nKnots_arr,
    float *gpu_x_array,
    float *gpu_many_array,

    short int *paramNo_arr,
    unsigned int *nKnots_arr,
    float *cpu_x_array, 
    float *cpu_many_array, 
#ifndef Weight_On_SplineBySpline_Basis
    int n_events,
    unsigned int *cpu_nParamPerEvent,
    unsigned int *gpu_nParamPerEvent,
 #endif   
    int n_params, 
    int n_splines,
    short int _max_knots,
    unsigned int coeff_array_size);

extern void CopyToGPU_TF1(
    float *gpu_coeff_many, 
    short int *gpu_paramNo_arr,
    short int *gpu_nPoints_arr,

    float *cpu_coeffs,
    short int *paramNo_arr,
    short int *nPoints_arr,
#ifndef Weight_On_SplineBySpline_Basis
    int n_events,
    unsigned int *cpu_nParamPerEvent,
    unsigned int *gpu_nParamPerEvent,
#endif  
    int nParams,
    int n_splines,
    short int _max_knots);

extern void RunGPU_SepMany(
    short int* gpu_paramNo_arr,
    unsigned int* gpu_nKnots_arr,

    float *gpu_coeff_many, 

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

    float *vals,
    short int *segment);

extern void RunGPU_TF1(
    float *gpu_coeff_many, 
    short int* gpu_paramNo_arr,
    short int* gpu_nPoints_arr,

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

    float *vals);


extern void CleanupGPU_SepMany(
    short int *gpu_paramNo_arr,
    unsigned int *gpu_nKnots_arr,

    float *gpu_x_array, 
    float *gpu_many_array, 

#ifndef Weight_On_SplineBySpline_Basis
    float *gpu_total_weights,
    unsigned int *gpu_nParamPerEvent,
    float *cpu_total_weights,
#endif 
    float *gpu_weights
);

extern void CleanupGPU_TF1(
    float *gpu_coeffs,
    short int *gpu_paramNo_arr,
    short int *gpu_nPoints_arr,
     
#ifdef Weight_On_SplineBySpline_Basis
    float *gpu_total_weights,
    float *cpu_total_weights,
#endif 
float *gpu_weights
);


// *****************************************
// Uses an x array and one combined yabd array
// This should optimise cache hitting because we use the same yabd points once we've found the x point
// So make these yabd points lay right next to each other in memory
SMonolith::SMonolith(std::vector<std::vector<TSpline3*> > &MasterSpline) {
  std::cout << "Using full TSpline3, about to reduce it and send to GPU" << std::endl;
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TSpline3_red*> > ReducedSpline = ReduceTSpline3(MasterSpline);
  PrepareForGPU(ReducedSpline);
}


// *****************************************
// Constructor for the reduced TSpline3 object
SMonolith::SMonolith(std::vector<std::vector<TSpline3_red*> > &MasterSpline) {
// *****************************************
  std::cout << "-- GPUING WITH {X} and {Y,B,C,D} arrays and master spline containing TSpline3_red" << std::endl;
  PrepareForGPU(MasterSpline);
}

// *****************************************
// Uses a fifth order polynomial for most shape except 2p2h shape C/O which are two superimposed linear eq
// Reduce first
SMonolith::SMonolith(std::vector<std::vector<TF1*> > &MasterSpline) {
// *****************************************
  std::cout << "Using full TF1, about to reduce it and send to GPU" << std::endl;
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TF1_red*> > ReducedSpline = ReduceTF1(MasterSpline);
  PrepareForGPU(ReducedSpline);
}

// *****************************************
// constructor for monotone spline 
SMonolith::SMonolith(std::vector<std::vector<Monotone_Spline*> > &MasterSpline) {
// *****************************************
  std::cout << "Using Monotone spline, about to convert it to TSpline3_red and send to GPU" << std::endl;
  // Convert the Monotone pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TSpline3_red*> > ReducedSpline = ReduceMonotone(MasterSpline);
  PrepareForGPU(ReducedSpline);
}

// *****************************************
// constructor for akima spline vector
SMonolith::SMonolith(std::vector<std::vector<Akima_Spline*> > &MasterSpline) {
// *****************************************
  std::cout << "Using Akima spline, about to convert it to TSpline3_red and send to GPU" << std::endl;
  // Convert the Akima pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TSpline3_red*> > ReducedSpline = ReduceAkima(MasterSpline);
  PrepareForGPU(ReducedSpline);
}

// *****************************************
// Uses a fifth order polynomial for most shape except 2p2h shape C/O which are two superimposed linear eq
// Reduce first
SMonolith::SMonolith(std::vector<std::vector<TF1_red*> > &MasterSpline) {
// *****************************************
  std::cout << "-- GPUING WITH TF1_red" << std::endl;
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  PrepareForGPU(MasterSpline);
}

// *****************************************
// The shared initialiser from constructors of TSpline3 and TSpline3_red
void SMonolith::PrepareForGPU(std::vector<std::vector<TSpline3_red*> > &MasterSpline) {
// *****************************************

  // Scan for the max number of knots, the number of events (number of splines), and number of parameters
  unsigned int NEvents = 0;
  _max_knots = 0;
  nParams = 0;
  int nSplines = 0;
  unsigned int nKnots = 0;
  ScanMasterSpline(MasterSpline, NEvents, _max_knots, nParams, nSplines, nKnots);
  std::cout << "Found " << NEvents << " events" << std::endl;
  std::cout << "Found " << _max_knots << " knots at max" << std::endl;
  std::cout << "Found " << nParams << " parameters" << std::endl;
  std::cout << "Found " << nSplines << " maximum number of splines in an event" << std::endl;
  std::cout << "Found total " << nKnots << " knots in all splines" << std::endl;

  // Total number of events in our Spline, read from TSpline3 entries
  // Number of TSpline3 we have in total if each event had the maximal number of splines (nSplines written by ScanMasterSpline)
  NSplines_total = NEvents * nSplines;
  // Number of TSpline3 we have in total if each event has *EVERY* spline. Needed for some arrays
  NSplines_total_large = NEvents*nParams;

  unsigned int event_size_max = _max_knots * nParams;
  // Declare the {x}, {y,b,c,d} arrays for all possible splines which the event has
  // We'll filter off the flat and "disabled" (e.g. CCQE event should not have MARES spline) ones in the next for loop, but need to declare these beasts here

  // Declare the {y,b,c,d} for each knot
  // float because GPU precision (could change to double, but will incur signficant speed reduction on GPU unless you're very rich!)
  float *coeff_many = new float[nKnots*_nCoeff_]; // *4 because we store y,b,c,d parameters in this array
  //KS: For x coeff we assume that for given dial (MAQE) spacing is identical, here we are sloppy and assume each dial has the same number of knots, not a big problem
  float *coeff_x = new float[event_size_max];

  // Set all the big arrays to -999 to keep us safe...
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int j = 0; j < nKnots; j++) {
    for (int k = 0; k < _nCoeff_; k++) {
      coeff_many[j*_nCoeff_+k] = -999;
    }
  }

  for (unsigned int j = 0; j < event_size_max; j++) {
    coeff_x[j] = -999;
  }
  // Will hold what spline number a certain spline has
  short int *paramNo_big = new short int[NSplines_total];
  // Will hold map where given spline starts
  unsigned int *knotNo_big = new unsigned int[NSplines_total];

  // This holds the index of each spline
  index_cpu = new int[NSplines_total_large];
  // This holds the total CPU weights that gets read in samplePDFND
#ifdef Weight_On_SplineBySpline_Basis
   cpu_weights = new float[NSplines_total_large];
#else
  //KS: Map keeping track how many parmaeters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
  unsigned int* cpu_nParamPerEvent = new unsigned int[2*NSplines_total];
  int ParamCounter = 0;
  int ParamCounterGlobal = 0;
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int j = 0; j < 2*NSplines_total; j++) {
    cpu_nParamPerEvent[j] = -1;
  }
#endif 
  // Set index array and number of points arrays to something silly, again to be safe...
#ifdef MULTITHREAD
#pragma omp parallel
{
//KS: Gain from nowait is minimal moslty done as I am learning how to better utilise OMP
  #pragma omp for nowait
#endif
  for (unsigned int j = 0; j < NSplines_total_large; j++) {
    index_cpu[j] = -1;
  }

  #ifdef MULTITHREAD
  #pragma omp for
  #endif
  for (unsigned int j = 0; j < NSplines_total; j++) {
    paramNo_big[j] = -1;
    knotNo_big[j] = -1;
  }
#ifdef MULTITHREAD
} //end parrarel region
#endif
  // Temporary arrays to hold the coefficients for each spline
  // We get one x, one y, one b,... for each point, so only need to be _max_knots big
  //KS: Some params has less splines but this is all right main array will get proper number while this temp will be deleted
  float *x_tmp = new float[_max_knots]();
  float *many_tmp = new float[_max_knots*_nCoeff_]();
  // Number of valid splines in total for our entire ensemble of events
  NSplines_valid = 0;

  std::vector<std::vector<TSpline3_red*> >::iterator OuterIt;
  std::vector<TSpline3_red*>::iterator InnerIt;
  // Count the number of events
  unsigned int EventCounter = 0;
  unsigned int KnotCounter = 0;
  // Loop over events and extract the spline coefficients
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++EventCounter) {
    // Structure of MasterSpline is std::vector<std::vector<TSpline3*>>
    // A conventional iterator to count which parameter a given spline should be applied to
    int ParamNumber = 0;
    for (InnerIt = (*OuterIt).begin(); InnerIt != (*OuterIt).end(); ++InnerIt, ++ParamNumber) {
      //KS: how much knots each spline has
      int nPoints_tmp = 0;
      // Get a pointer to the current spline for this event
      TSpline3_red* CurrSpline = (*InnerIt);
      // If NULL we don't have this spline for the event, so move to next spline
      if (CurrSpline == NULL) continue;

      // If the number of knots are greater than 2 the spline is not a dummy and we should extract coefficients to load onto the GPU
      getSplineCoeff_SepMany(CurrSpline, nPoints_tmp, x_tmp, many_tmp);
      
      //KS: One knot means flat spline so ignore
      if (nPoints_tmp == 1) continue;
      
      (*InnerIt) = CurrSpline;
      for (int j = 0; j < _max_knots; ++j) {
        coeff_x[ParamNumber*_max_knots + j] = x_tmp[j];
      }
      //KS: Contrary to X coeff we keep for other coeef only filled knots, there is no much gain for doing so for x coeff
      for (int j = 0; j < nPoints_tmp; ++j) {
        for (int k = 0; k < _nCoeff_; k++) {
          coeff_many[KnotCounter*_nCoeff_ + j*_nCoeff_ + k] = many_tmp[j*_nCoeff_+k];
        }
      }

      // Set the parameter number for this spline
      paramNo_big[NSplines_valid] = ParamNumber;
      //KS: Fill map when each spline starts
      knotNo_big[NSplines_valid] = KnotCounter;

      KnotCounter += nPoints_tmp;
      // Set the index of the spline so we can tell apart from flat splines
      index_cpu[EventCounter*nParams + ParamNumber] = NSplines_valid;
      #ifndef Weight_On_SplineBySpline_Basis
      ParamCounter++;
      #endif
      // Incremement the counter for the number of good splines we have
      ++NSplines_valid;
    } // End the loop over the parameters in the MasterSpline
      #ifndef Weight_On_SplineBySpline_Basis
      cpu_nParamPerEvent[2*EventCounter] = ParamCounter;
      cpu_nParamPerEvent[2*EventCounter+1] = ParamCounterGlobal;
      ParamCounterGlobal += ParamCounter;
      ParamCounter = 0;
      #endif
  } // End the loop over the number of events
  delete[] many_tmp;
  delete[] x_tmp;
  
  // Now that we have looped through all events we can make the number of splines smaller
  // Going from number of events * number of points per spline * number of NIWG params to spln_counter (=valid splines)

  std::cout << "  Number of splines = " << NSplines_valid << std::endl;

  // Make array with the number of points per spline (not per spline point!)
  short int *paramNo_arr = new short int[NSplines_valid];
  //KS: And array which tells where each spline stars in a big monolit array, sort of knot map
  unsigned int *knotNo_arr = new unsigned int[NSplines_valid];

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < NSplines_valid; i++) {
    paramNo_arr[i] = paramNo_big[i];
    knotNo_arr[i] = knotNo_big[i];
    // Perform checks that all array entries have been changed from negative numbers and no number of points is greater than max knots, inputted by user
    // Don't do this for the index array since some entries there should be -1 so we know what splines to include and not include in each event for loading onto the GPU
    if (paramNo_arr[i] < 0) {
      std::cerr << "***** NEGATIVE PARAMETER NUMBER!!! ***** \n" << "On spline " << i << " " << paramNo_arr[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
    if (knotNo_arr[i] < 0) {
      std::cerr << "***** NEGATIVE Knot NUMBER!!! ***** \n" << "On spline " << i << " " << knotNo_arr[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
  }
  delete[] paramNo_big;
  delete[] knotNo_big;

  for (unsigned int j = 0; j < event_size_max; j++) {
    // Perform checks that all entries have been modified from intial values
    if (coeff_x[j] == -999) {
      std::cerr << "***** BAD X !! ***** \n" << std::endl;
      std::cerr << "Indicates some parameter doesn't have a single spline" << std::endl;
      std::cerr << "j = " << j << std::endl;
      std::cerr << __FILE__ << "::" << __LINE__ << std::endl;
      //throw;
    }
  }

  // Print some info; could probably make this to a separate function
  std::cout << "--- INITIALISED {X}, {YBCD} ARRAYS ---" << std::endl;
  std::cout << "  " << NEvents << " events with " << NSplines_valid << " splines" << std::endl;

  std::cout << "  On average " << float(NSplines_valid)/float(NEvents) << " splines per event (" << NSplines_valid << "/" << NEvents << ")" << std::endl;

  std::cout << "  Size of x array = " << double(sizeof(float)*event_size_max)/1.E6 << " MB" << std::endl;
  std::cout << "  Size of coefficient {y,b,c,d} array = " << double(sizeof(float)*nKnots*_nCoeff_)/1.E6 << " MB" << std::endl;

  std::cout << "  Size of parameter # array = " << double(sizeof(short int)*NSplines_valid)/1.E6 << " MB" << std::endl;

  std::cout << "  Total size = " << (double(sizeof(float)*nKnots*_nCoeff_)+ double(sizeof(float)*event_size_max)/1.E6+double(sizeof(short int)*NSplines_valid))/1.E6 << " MB memory on CPU to move to GPU" << std::endl;

  std::cout << "  GPU weight array (GPU->CPU every step) = " << double(sizeof(float)*NSplines_valid)/1.E6 << " MB" << std::endl;
#ifndef Weight_On_SplineBySpline_Basis  
  std::cout << "  Since you are runing Total event weight mode then GPU weight array (GPU->CPU every step) = " << double(sizeof(float)*NEvents)/1.E6 << " MB" << std::endl;
#endif  
  std::cout << "  Parameter value array (CPU->GPU every step) = " << double(sizeof(float)*nParams)/1.E6 << " MB" << std::endl;

#ifdef Weight_On_SplineBySpline_Basis
  // Make the array that holds all the returned weights from the GPU to pass to the CPU
  cpu_weights_var = new float[NSplines_valid]();
#endif
  // With the new set-up we have:   1 coefficient array of size coeff_array_size, all same size
  //                                1 coefficient array of size coeff_array_size*4, holding y,b,c,d in order (y11,b11,c11,d11; y12,b12,c12,d12;...) where ynm is n = spline number, m = spline point. Should really make array so that order is (y11,b11,c11,d11; y21,b21,c21,d21;...) because it will optimise cache hits I think; try this if you have time
  //                                return gpu_weights

  // The gpu_XY arrays don't actually need initialising, since they are only placeholders for what we'll move onto the GPU. As long as we cudaMalloc the size of the arrays correctly there shouldn't be any problems
  // Can probably make this a bit prettier but will do for now
  // Could be a lot smaller of a function...
  InitGPU_SepMany(
      &gpu_coeff_x,
      &gpu_coeff_many,
      &gpu_weights,

      &gpu_paramNo_arr,
      &gpu_nKnots_arr,
#ifndef Weight_On_SplineBySpline_Basis
      &cpu_total_weights,
      &gpu_total_weights, 
      NEvents,
                  
      &gpu_nParamPerEvent,
#endif     
      nKnots, // How many entries in coefficient array (*4 for the "many" array)
      NSplines_valid, // What's the number of splines we have (also number of entries in gpu_nPoints_arr)
      event_size_max //Knots times event number of uniqe splines           
);

  // Move number of splines and spline size to constant GPU memory; every thread does not need a copy...
  // The implementation lives in splines/gpuSplineUtils.cu
  // The GPU splines don't actually need declaring but is good for demonstration, kind of
  // fixed by passing const reference
  CopyToGPU_SepMany(
      gpu_paramNo_arr,
      gpu_nKnots_arr,
      gpu_coeff_x, 
      gpu_coeff_many,

      paramNo_arr,
      knotNo_arr,
      coeff_x,
      coeff_many,
#ifndef Weight_On_SplineBySpline_Basis
      NEvents,
      cpu_nParamPerEvent,
      gpu_nParamPerEvent,
#endif    
      nParams,
      NSplines_valid,
      _max_knots,
      nKnots);

  // Delete all the coefficient arrays from the CPU once they are on the GPU
  delete[] coeff_x;
  delete[] coeff_many;
  delete[] paramNo_arr;
  delete[] knotNo_arr;
#ifndef Weight_On_SplineBySpline_Basis
  delete[] cpu_nParamPerEvent;
#endif 
  std::cout << "Good GPU loading" << std::endl;
}

// *****************************************
// The shared initialiser from constructors of TF1 and TF1_red
void SMonolith::PrepareForGPU(std::vector<std::vector<TF1_red*> > &MasterSpline) {
// *****************************************

  // Scan for the max number of knots, the number of events (number of splines), and number of parameters
  unsigned int NEvents = 0;
  _max_knots = 0;
  nParams = 0;
  ScanMasterSpline(MasterSpline, NEvents, _max_knots, nParams);
  std::cout << "Found " << NEvents << " events" << std::endl;
  std::cout << "Found " << _max_knots << " polynomial at max" << std::endl;
  std::cout << "Found " << nParams << " parameters" << std::endl;

  // Every event maximally has nParams TF1s which we've saved
  NSplines_total = NEvents * nParams;

  // With TF1 we only save the coefficients and the order of the polynomial
  // Makes most sense to have one large monolithic array, but then it becomes impossible to tell apart a coefficient from a "number of points". So have two arrays: one of cofficients and one of number of points
  // Let's first assume all are of _max_knots size
  short int *nPoints_big = new short int[NSplines_total];
  float *coeffs_big = new float[NSplines_total*_nTF1Coeff_]; // *5 because we store a,b,c,d,e coefficients of fifth order polynomial in this array

  // Set all the big arrays to -999 to keep us safe...
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int j = 0; j < NSplines_total; j++) {
    // 1 number of points for each spline
    nPoints_big[j] = -999;
    // 5 coefficients for each spline
    for (int k = 0; k < _nTF1Coeff_; k++) {
      coeffs_big[j*_nTF1Coeff_+k] = -999;
    }
  }

  // Will hold what spline number a certain spline has
  short int *paramNo_big = new short int[NSplines_total];

  // This holds the index of each spline
  index_cpu = new int[NSplines_total];
#ifdef Weight_On_SplineBySpline_Basis
  //This holds the CPU weights for EACH SPLINE that gets read in samplePDFND
  cpu_weights = new float[NSplines_total];
#else
  //KS: Map keeping track how many parmaeters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
  unsigned int* cpu_nParamPerEvent = new unsigned int[2*NSplines_total];
  int ParamCounter = 0;
  int ParamCounterGlobal = 0;
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int j = 0; j < 2*NSplines_total; j++) {
    cpu_nParamPerEvent[j] = -1;
  }
#endif 
  // Set index array and number of points arrays to something silly, again to be safe...
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int j = 0; j < NSplines_total; j++) {
    index_cpu[j] = -1;
    paramNo_big[j] = -1;
    #ifndef Weight_On_SplineBySpline_Basis
    cpu_nParamPerEvent[j*2] = -1;
    cpu_nParamPerEvent[j*2+1] = -1;
    #endif
  }

  // Temporary arrays to hold the coefficients for each spline
  // We get one a,b,c,d,e for each spline so need 5 
  float *temp_coeffs = new float[_max_knots]();
  NSplines_valid = 0;

  std::vector<std::vector<TF1_red*> >::iterator OuterIt;
  std::vector<TF1_red*>::iterator InnerIt;
  unsigned int EventCounter = 0;
  // Loop over events and extract the spline coefficients
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++EventCounter) {
    // Structure of MasterSpline is std::vector<std::vector<TSpline3*>>
    // A conventional iterator to count which parameter a given spline should be applied to
    int ParamNumber = 0;
    for (InnerIt = (*OuterIt).begin(); InnerIt != (*OuterIt).end(); ++InnerIt, ++ParamNumber) {
      // Don't actually use this ever -- we give each spline the maximum number of points found in all splines
      int nPoints_tmp = 0;
      // Get a pointer to the current spline for this event
      TF1_red* CurrSpline = (*InnerIt);
      // If NULL we don't have this spline for the event, so move to next spline
      if (CurrSpline == NULL) continue;

      // If the number of knots are greater than 2 the spline is not a dummy and we should extract coefficients to load onto the GPU
      getTF1Coeff(CurrSpline, nPoints_tmp, temp_coeffs);
      (*InnerIt) = CurrSpline;
      for (int j = 0; j < _max_knots; ++j) {
        coeffs_big[NSplines_valid*_max_knots+j] = temp_coeffs[j];
      }
      // Save the number of points for this spline
      nPoints_big[NSplines_valid] = nPoints_tmp;

      // Set the parameter number for this spline
      paramNo_big[NSplines_valid] = ParamNumber;
      // Set the index of the spline so we can tell apart from flat splines
      index_cpu[EventCounter*nParams + ParamNumber] = NSplines_valid;
      // Incremement the counter for the number of good splines we have
      ++NSplines_valid;
      #ifndef Weight_On_SplineBySpline_Basis
      ParamCounter++;
      #endif
    } // End the loop over the parameters in the MasterSpline
      #ifndef Weight_On_SplineBySpline_Basis
      cpu_nParamPerEvent[2*EventCounter] = ParamCounter;
      ParamCounterGlobal += ParamCounter;
      cpu_nParamPerEvent[2*EventCounter+1] = ParamCounterGlobal;
      ParamCounter = 0;
      #endif
  } // End the loop over the number of events
  // Delete the temporary arrays
  delete[] temp_coeffs;

  // Now that we have looped through all events we can make the number of splines smaller
  // Going from number of events * number of points per spline * number of NIWG params to spln_counter (=valid splines)

  std::cout << "  Number of splines = " << NSplines_valid << std::endl;

  // Now declare the arrays for each point in the valid splines which the event actually has (i.e. include the splines that the event undergoes)
  // Also make array with the number of points per spline (not per spline point!)
  // float because GPU precision (could change to double, but will incur signficant speed reduction on GPU unless you're very rich!)
  short int *nPoints = new short int[NSplines_valid];
  float *coeffs = new float[NSplines_valid*_nTF1Coeff_]; // *5 because this array holds  a,b,c,d,e parameters
  short int *paramNo = new short int[NSplines_valid];

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < NSplines_valid; i++) {
    nPoints[i] = nPoints_big[i];
    paramNo[i] = paramNo_big[i];
    for (unsigned int j = 0; j < _nTF1Coeff_; ++j) {
      coeffs[i*_nTF1Coeff_+j] = coeffs_big[i*_nTF1Coeff_+j];
      if (coeffs[i*_nTF1Coeff_+j] == -999) {
        std::cerr << "***** BAD COEFFICIENT OF POLY!!! ***** \n" << "On spline " << i << " = " << coeffs[i*_nTF1Coeff_+j] << std::endl;
        std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }

    // Perform checks that all array entries have been changed from negative numbers and no number of points is greater than max knots, inputted by user
    // Don't do this for the index array since some entries there should be -1 so we know what splines to include and not include in each event for loading onto the GPU
    if (paramNo[i] < 0) {
      std::cerr << "***** NEGATIVE PARAMETER NUMBER!!! ***** \n" << "On spline " << i << " = " << paramNo[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    if (nPoints[i] < 0) {
      std::cerr << "***** NEGATIVE NUMBER OF POINTS!!! ***** \n" << "On spline " << i << " = " << nPoints[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

  }  
  // Delete allocated memory
  delete[] paramNo_big;
  delete[] nPoints_big;
  delete[] coeffs_big;

  // Print some info; could probably make this to a separate function
  std::cout << "--- INITIALISED TF1 ARRAYS ---" << std::endl;
  std::cout << "  " << NEvents << " events with " << NSplines_valid << " splines" << std::endl;

  std::cout << "  On average " << float(NSplines_valid)/float(NEvents) << " splines per event (" << NSplines_valid << "/" << NEvents << ")" << std::endl;

  std::cout << "  Size of coefficient {a,b,c,d,e} array = " << double(sizeof(float)*NSplines_valid*_nTF1Coeff_)/1.E6 << " MB" << std::endl;

  std::cout << "  Size of parameter # array = " << double(sizeof(short int)*NSplines_valid)/1.E6 << " MB" << std::endl;
  std::cout << "  Size of polynomial type array = " << double(sizeof(short int)*NSplines_valid)/1.E6 << " MB" << std::endl;

  std::cout << "  Total size = " << (double(sizeof(float)*NSplines_valid*_nTF1Coeff_)+double(2.0*sizeof(short int)*NSplines_valid))/1.E6 << " MB memory on CPU to move to GPU" << std::endl;

  std::cout << "  GPU weight array (GPU->CPU every step) = " << double(sizeof(float)*NSplines_valid)/1.E6 << " MB" << std::endl;
  std::cout << "  Parameter value array (CPU->GPU every step) = " << double(sizeof(float)*nParams)/1.E6 << " MB" << std::endl;
#ifdef Weight_On_SplineBySpline_Basis
  // Make the array that holds all the returned weights from the GPU to pass to the CPU
  cpu_weights_var = new float[NSplines_valid]();
#endif
  // With the new set-up we have:   1 coefficient array of size coeff_array_size, all same size
  //                                1 coefficient array of size coeff_array_size*4, holding y,b,c,d in order (y11,b11,c11,d11; y12,b12,c12,d12;...) where ynm is n = spline number, m = spline point. Should really make array so that order is (y11,b11,c11,d11; y21,b21,c21,d21;...) because it will optimise cache hits I think; try this if you have time
  //                                return gpu_weights

  // The gpu_XY arrays don't actually need initialising, since they are only placeholders for what we'll move onto the GPU. As long as we cudaMalloc the size of the arrays correctly there shouldn't be any problems
  // Can probably make this a bit prettier but will do for now
  // Could be a lot smaller of a function...
  InitGPU_TF1(
      &gpu_coeff_many,
      &gpu_paramNo_arr,
      &gpu_nPoints_arr,

      &gpu_weights,

#ifndef Weight_On_SplineBySpline_Basis
      &cpu_total_weights,
      &gpu_total_weights, 
      NEvents,        
      &gpu_nParamPerEvent,
#endif 
      NSplines_valid); // What's the number of splines we have (also number of entries in gpu_nPoints_arr)

  // Move number of splines and spline size to constant GPU memory; every thread does not need a copy...
  // The implementation lives in splines/gpuSplineUtils.cu
  // The GPU splines don't actually need declaring but is good for demonstration, kind of
  // fixed by passing const reference
  CopyToGPU_TF1(
      gpu_coeff_many,
      gpu_paramNo_arr,
      gpu_nPoints_arr,

      coeffs,
      paramNo,
      nPoints,
#ifndef Weight_On_SplineBySpline_Basis
      NEvents,
      cpu_nParamPerEvent,
      gpu_nParamPerEvent,
#endif  
      nParams,
      NSplines_valid,
      _max_knots);

  // Delete all the coefficient arrays from the CPU once they are on the GPU
  delete[] coeffs;
  delete[] nPoints;
  delete[] paramNo;
#ifndef Weight_On_SplineBySpline_Basis
  delete[]cpu_nParamPerEvent;
#endif  
  std::cout << "Good TF1 GPU loading" << std::endl;
}

// Need to specify template functions in header
// *****************************************
// Scan the master spline to get the maximum number of knots in any of the TSpline3*
void SMonolith::ScanMasterSpline(std::vector<std::vector<TSpline3_red*> > & MasterSpline, unsigned int &nEvents, int &MaxPoints, short int &nParams, int &nSplines, unsigned int &nKnots) {
// *****************************************

  // Need to extract: the total number of events
  //                  number of parameters
  //                  maximum number of knots
  MaxPoints = 0;
  nEvents   = 0;
  nParams   = 0;
  nSplines = 0;
  nKnots = 0;
  std::vector<std::vector<TSpline3_red*> >::iterator OuterIt;
  std::vector<TSpline3_red*>::iterator InnerIt;

  // Check the number of events
  nEvents = MasterSpline.size();

  // Maximum number of splines one event can have (scan through and find this number)
  int nMaxSplines_PerEvent = 0;

  unsigned int EventCounter = 0;
  // Loop over each parameter
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt) {
    // Check that each event has each spline saved
    if (nParams > 0) {
      int TempSize = (*OuterIt).size();
      if (TempSize != nParams) {
        std::cerr << "Found " << TempSize << " parameters for event " << EventCounter << std::endl;
        std::cerr << "but was expecting " << nParams << " since that's what I found for the previous event" << std::endl;
        std::cerr << "Somehow this event has a different number of spline parameters... Please study further!" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }
    nParams = (*OuterIt).size();

    int nSplines_SingleEvent = 0;
    // Loop over each pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt) {
      if ((*InnerIt) == NULL) continue;
      int nPoints = (*InnerIt)->GetNp();
      if (nPoints > MaxPoints) {
        MaxPoints = nPoints;
      }
      nKnots += nPoints;
      nSplines_SingleEvent++;
    }
    if (nSplines_SingleEvent > nMaxSplines_PerEvent) nMaxSplines_PerEvent = nSplines_SingleEvent;
    EventCounter++;
  }
  nSplines = nMaxSplines_PerEvent;
}

// Need to specify template functions in header
// *****************************************
// Scan the master spline to get the maximum number of knots in any of the TSpline3*
void SMonolith::ScanMasterSpline(std::vector<std::vector<TF1_red*> > & MasterSpline, unsigned int &nEvents, int &MaxPoints, short int &nParams) {
// *****************************************

  // Need to extract: the total number of events
  //                  number of parameters
  //                  maximum number of knots
  MaxPoints = 0;
  nEvents   = 0;
  nParams   = 0;

  std::vector<std::vector<TF1_red*> >::iterator OuterIt;
  std::vector<TF1_red*>::iterator InnerIt;

  // Check the number of events
  nEvents = MasterSpline.size();

  unsigned int EventCounter = 0;
  // Loop over each parameter
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt) {
    // Check that each event has each spline saved
    if (nParams > 0) {
      int TempSize = (*OuterIt).size();
      if (TempSize != nParams) {
        std::cerr << "Found " << TempSize << " parameters for event " << EventCounter << std::endl;
        std::cerr << "but was expecting " << nParams << " since that's what I found for the previous event" << std::endl;
        std::cerr << "Somehow this event has a different number of spline parameters... Please study further!" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }
    nParams = (*OuterIt).size();

    // Loop over each pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt) {
      if ((*InnerIt) == NULL) continue;
      int nPoints = (*InnerIt)->GetSize();
      if (nPoints > MaxPoints) {
        MaxPoints = nPoints;
      }
    }
    EventCounter++;
  }
}

// *****************************************
// Destructor
// Cleans up the allocated GPU memory
SMonolith::~SMonolith() {
// *****************************************

  CleanupGPU_SepMany(
      gpu_paramNo_arr,
      gpu_nKnots_arr,

      gpu_coeff_x,
      gpu_coeff_many,
      
#ifndef Weight_On_SplineBySpline_Basis
      gpu_total_weights,
      gpu_nParamPerEvent,
      cpu_total_weights,
#endif
      gpu_weights);
  
#ifdef Weight_On_SplineBySpline_Basis
  delete[] cpu_weights;
  delete[] cpu_weights_var;
#endif
  delete[] index_cpu;
}

// *********************************
// Reduce the large TSpline3 vector to TSpline3_red
std::vector<std::vector<TSpline3_red*> > SMonolith::ReduceTSpline3(std::vector<std::vector<TSpline3*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<TSpline3*> >::iterator OuterIt;
  std::vector<TSpline3*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TSpline3_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TSpline3_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      TSpline3 *spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer
      TSpline3_red *red = NULL;
      if (spline != NULL) {
        red = new TSpline3_red(spline);
        (*InnerIt) = spline;
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}

// *********************************
// convert Monotone_Spline vector to TSpline3_red
std::vector<std::vector<TSpline3_red*> > SMonolith::ReduceMonotone(std::vector<std::vector<Monotone_Spline*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<Monotone_Spline*> >::iterator OuterIt;
  std::vector<Monotone_Spline*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TSpline3_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TSpline3_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      Monotone_Spline *spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer
      TSpline3_red *red = NULL;
      if (spline != NULL) {
        red = spline->ConstructTSpline3_red();
        (*InnerIt) = spline; 
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}


// *********************************
// convert Akima_Spline vector to TSpline3_red
std::vector<std::vector<TSpline3_red*> > SMonolith::ReduceAkima(std::vector<std::vector<Akima_Spline*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<Akima_Spline*> >::iterator OuterIt;
  std::vector<Akima_Spline*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TSpline3_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TSpline3_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer 
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      Akima_Spline *spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer (which deleted TSpline3)
      TSpline3_red *red = NULL;

      if (spline != NULL) {
        red = spline->ConstructTSpline3_red(); 
        (*InnerIt) = spline;
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}


// *********************************
// Reduce the large TF1 vector to a TF1_red
std::vector<std::vector<TF1_red*> > SMonolith::ReduceTF1(std::vector<std::vector<TF1*> > &MasterSpline) {
// *********************************
  std::vector<std::vector<TF1*> >::iterator OuterIt;
  std::vector<TF1*>::iterator InnerIt;

  // The return vector
  std::vector<std::vector<TF1_red*> > ReducedVector;
  ReducedVector.reserve(MasterSpline.size());

  // Loop over each parameter
  int OuterCounter = 0;
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt, ++OuterCounter) {
    // Make the temp vector
    std::vector<TF1_red*> TempVector;
    TempVector.reserve(OuterIt->size());
    int InnerCounter = 0;
    // Loop over each TSpline3 pointer
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ++InnerCounter) {
      // Here's our delicious TSpline3 object
      TF1* spline = (*InnerIt);
      // Now make the reduced TSpline3 pointer (which deleted TSpline3)
      TF1_red* red = NULL;
      if (spline != NULL) {
        red = new TF1_red(spline);
        (*InnerIt) = spline;
      }
      // Push back onto new vector
      TempVector.push_back(red);
    } // End inner for loop
    ReducedVector.push_back(TempVector);
  } // End outer for loop
  // Now have the reduced vector
  return ReducedVector;
}

// *****************************************
// Get the spline coefficients from the TSpline3 so that we can load ONLY these onto the GPU, not the whole TSpline3 object
// This loads up coefficients into two arrays: one x array and one yabcd array
// This should maximize our cache hits!
void SMonolith::getSplineCoeff_SepMany(TSpline3_red* &spl, int &nPoints, float *& xArray, float *& manyArray) {
// *****************************************
  // Initialise all arrays to 1.0
  for (int i = 0; i < _max_knots; ++i) {
    xArray[i] = 1.0;
    for (int j = 0; j < _nCoeff_; j++) {
      manyArray[i*_nCoeff_+j] = 1.0;
    }
  }
  // Get number of points in spline
  int Np = spl->GetNp();
  // If spline is flat, set number of knots to 1.0,
  // This is used later to expedite the calcuations for flat splines
  // tmpArray[0] is number of knots
  if (isFlat(spl)) {
    nPoints = 1;
  } else {
    nPoints = Np;
    if (Np > _max_knots) {
      std::cerr << "Error, number of points is greater than saved " << _max_knots << std::endl;
      std::cerr << "This _WILL_ cause problems with GPU splines and _SHOULD_ be fixed!" << std::endl;
      std::cerr << "nPoints = " << nPoints << ", _max_knots = " << _max_knots << std::endl;
      std::cerr << __FILE__ << "::" << __LINE__ << std::endl;
      throw;
    }
  }

  // The coefficients we're writing to
  double x, y, b, c, d;
  // TSpline3 can only take doubles, not floats
  // But our GPU is slow with doubles, so need to cast to float
  for(int i = 0; i < Np; i++) {
    // Get the coefficients from the TSpline3 object
    spl->GetCoeff(i, x, y, b, c, d);
    // Write the arrays
    xArray[i] = float(x);
    manyArray[i*_nCoeff_] = float(y); // 4 because manyArray stores y,b,c,d
    manyArray[i*_nCoeff_+1] = float(b);
    manyArray[i*_nCoeff_+2] = float(c);
    manyArray[i*_nCoeff_+3] = float(d);
    
    if((xArray[i] == -999) | (manyArray[i*_nCoeff_] == -999) | (manyArray[i*4+1] == -999) | (manyArray[i*_nCoeff_+2] == -999) | (manyArray[i*_nCoeff_+3] == -999)){
      std::cerr << "*********** Bad params in getSplineCoeff_SepMany() ************"<<std::endl;
      std::cerr << "pre cast to float (x, y, b, c, d) = "<<x<<", "<<y<<", "<<b<<", "<<c<<", "<<d<<std::endl;
      std::cerr << "post cast to float (x, y, b, c, d) = "<<xArray[i]<<", "<<manyArray[i*4]<<", "<<manyArray[i*4+1]<<", "<<manyArray[i*4+2]<<", "<<manyArray[i*_nCoeff_+3]<<std::endl;
      std::cerr << "This will cause problems when preparing for GPU"<<std::endl;
      std::cerr << "***************************************************************"<<std::endl;
    }
  }
  // The structure is now xarray  ={x1,x2,x3} 
  //                      manyArr ={y1,y2,y3, b1,b2,b3, c1,c2,c3, d1,d2,d3}
#ifndef DEBUG_DUMP
  delete spl;
  spl = NULL;
#endif
}

// *****************************************
// Get the spline coefficients from the TSpline3 so that we can load ONLY these onto the GPU, not the whole TSpline3 object
// This loads up coefficients into two arrays: one x array and one yabcd array
// This should maximize our cache hits!
void SMonolith::getTF1Coeff(TF1_red* &spl, int &nPoints, float *& coeffs) {
// *****************************************

  // Initialise all arrays to 1.0
  for (int i = 0; i < _max_knots; ++i) {
    coeffs[i] = 0.0;
  }

  // Get number of points in spline
  nPoints = spl->GetSize();

  // TSpline3 can only take doubles, not floats
  // But our GPU is slow with doubles, so need to cast to float
  for (int i = 0; i < nPoints; i++) {
    coeffs[i] = spl->GetParameter(i);
  }
  // The structure is now coeffs  = {a,b,c,d,e}
#ifndef DEBUG_DUMP
  delete spl;
  spl = NULL;
#endif
}

// *****************************************
// Check if the TSpline3 object is flat; if it is, we won't load it onto the GPU
bool SMonolith::isFlat(TSpline3_red* &spl) {
// *****************************************
  int Np = spl->GetNp();
  double x, y, b, c, d;
  // Go through spline segment parameters,
  // Get y values for each spline knot,
  // Every knot must evaluate to 1.0 to create a flat spline
  for(int i = 0; i < Np; i++) {
    spl->GetCoeff(i, x, y, b, c, d);
    if (y != 1) {
      return false;
    }
  }
  return true;
}


// *****************************************
// Tell the GPU to evaluate the weights
// Load up the two x,{y,b,c,d} arrays into memory and have GPU read them with more coalescense instead of one monolithic array
// This should be used when we're using separate x,y,a,b,c,d arrays
// Also pass the segments for the parameter along with their parameter values
// This avoids doing lots of binary searches on the GPU
void SMonolith::EvalGPU_SepMany(float* vals, short int *segment, bool plotWeight) {
// *****************************************

#ifdef DEBUG
  TStopwatch clock;
  clock.Start();
#endif

  // The main call to the GPU
  RunGPU_SepMany(
      gpu_paramNo_arr,
      gpu_nKnots_arr,

      gpu_coeff_many, 

      gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
      cpu_weights_var,
#else
    gpu_total_weights,
    cpu_total_weights,
#endif
      vals,
      segment);

#ifdef DEBUG
  clock.Stop();
  std::cout << "RunGPU_SepMany " << clock.RealTime() << "s" << std::endl;
#endif

#ifdef Weight_On_SplineBySpline_Basis
  #ifdef DEBUG
  clock.Start();
  #endif

  // Multi-thread here because _numIndex is really quite large!
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int i = 0; i < NSplines_total_large; i++) {
    if (index_cpu[i] >= 0) {
      cpu_weights[i] = cpu_weights_var[index_cpu[i]];
    } else {
      cpu_weights[i] = 1.;
    }
  }

  #ifdef DEBUG
  clock.Stop();
  std::cout << "splineMonolith MP part took " << clock.RealTime() << "s" << std::endl;

  if (plotWeight == true) {
    // Weight plots; fill with gpu_weights first
    std::string name = "weightFile_GPU_new";
    int val = vals[0];
    std::stringstream ss;
    ss << val << ".root";
    TFile *weightFile = new TFile((name+ss.str()).c_str(),"RECREATE");
    TH1F *weights = new TH1F("weights","weights", 100, 0, 1.5);
    TH1F *weightsBig = new TH1F("weightsBig","weightsBig", 100, 2, 10);
    TGraph *weightsPlot = new TGraph(NSplines_valid);

    for (unsigned int i = 0; i < NSplines_valid; i++) {
      weights->Fill(cpu_weights_var[i]);
      weightsBig->Fill(cpu_weights_var[i]);
      weightsPlot->SetPoint(i, i, cpu_weights_var[i]);
    }

    weightFile->cd();
    weights->Write();
    weightsBig->Write();
    weightsPlot->Write();
    std::cout << "Wrote " << weightFile->GetName() << " to file" << std::endl;
    weightFile->Close();
  }
  #endif
#endif
}

// *****************************************
// Tell the GPU to evaluate the weights
// TF1 version
void SMonolith::EvalGPU_TF1(float* vals, bool plotWeight) {
// *****************************************

#ifdef DEBUG
  TStopwatch clock;
  clock.Start();
#endif

  // The main call to the GPU
  RunGPU_TF1(
      gpu_coeff_many, 
      gpu_paramNo_arr,
      gpu_nPoints_arr,

      gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
      cpu_weights_var,
#else
    gpu_total_weights,
    cpu_total_weights,
#endif
      vals);

#ifdef DEBUG
  clock.Stop();
  std::cout << "RunGPU_TF1 " << clock.RealTime() << "s" << std::endl;
#endif

#ifdef Weight_On_SplineBySpline_Basis
  #ifdef DEBUG
  clock.Start();
  #endif

  // Multi-thread here because _numIndex is really quite large!
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int i = 0; i < NSplines_total; i++) {
    if (index_cpu[i] >= 0) {
      cpu_weights[i] = cpu_weights_var[index_cpu[i]];
    } else {
      cpu_weights[i] = 1.;
    }
  }

  #ifdef DEBUG
  clock.Stop();
  std::cout << "splineMonolith MP part took " << clock.RealTime() << "s" << std::endl;

  if (plotWeight == true) {
    // Weight plots; fill with gpu_weights first
    std::string name = "weightFile_GPU_TF1";
    int val = vals[0];
    std::stringstream ss;
    ss << val << ".root";
    TFile *weightFile = new TFile((name+ss.str()).c_str(),"RECREATE");
    TH1F *weights = new TH1F("weights","weights", 100, 0, 1.5);
    TH1F *weightsBig = new TH1F("weightsBig","weightsBig", 100, 2, 10);
    TGraph *weightsPlot = new TGraph(NSplines_valid);

    for (unsigned int i = 0; i < NSplines_valid; i++) {
      weights->Fill(cpu_weights_var[i]);
      weightsBig->Fill(cpu_weights_var[i]);
      weightsPlot->SetPoint(i, i, cpu_weights_var[i]);
    }

    weightFile->cd();
    weights->Write();
    weightsBig->Write();
    weightsPlot->Write();
    std::cout << "Wrote " << weightFile->GetName() << " to file" << std::endl;
    weightFile->Close();
  }
  #endif
#endif
}
