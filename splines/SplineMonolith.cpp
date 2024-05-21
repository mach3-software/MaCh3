#include "SplineMonolith.h"

// All the functions that are declared extern here exist in splines/gpuSplineUtils.cu
// They are the CUDA code which we use to do the GPU processing
// For older style (Rich/Asher era) see git commit history previous to 27 Nov 2017

#ifdef CUDA
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
    unsigned int n_splines,
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

    unsigned int n_splines);

extern void CopyToGPU_SepMany(
    short int *gpu_paramNo_arr,
    unsigned int *gpu_nKnots_arr,
    float *gpu_x_array,
    float *gpu_many_array,

    std::vector<short int> cpu_paramNo_arr,
    std::vector<unsigned int> cpu_nKnots_arr,
    std::vector<float> cpu_x_array,
    std::vector<float> cpu_coeff_many,
#ifndef Weight_On_SplineBySpline_Basis
    int n_events,
    std::vector<unsigned int> cpu_nParamPerEvent,
    unsigned int *gpu_nParamPerEvent,
 #endif   
    int n_params, 
    unsigned int n_splines,
    short int _max_knots,
    unsigned int coeff_array_size);

extern void CopyToGPU_TF1(
    float *gpu_coeff_many, 
    short int *gpu_paramNo_arr,
    short int *gpu_nPoints_arr,

    std::vector<float> cpu_coeff_many,
    std::vector<short int> cpu_paramNo_arr,
    std::vector<short int> cpu_nPoints_arr,
#ifndef Weight_On_SplineBySpline_Basis
    int n_events,
    std::vector<unsigned int> cpu_nParamPerEvent,
    unsigned int *gpu_nParamPerEvent,
#endif  
    int nParams,
    unsigned int n_splines,
    short int _max_knots);

extern void RunGPU_SepMany(
    const short int* gpu_paramNo_arr,
    const unsigned int* gpu_nKnots_arr,

    const float *gpu_coeff_many,

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

    float *vals,
    short int *segment,
    const unsigned int h_n_splines);

extern void RunGPU_TF1(
    const float *gpu_coeff_many,
    const short int* gpu_paramNo_arr,
    const short int* gpu_nPoints_arr,

    float* gpu_weights, 
#ifdef Weight_On_SplineBySpline_Basis
    float* cpu_weights,
#else
    float* gpu_total_weights,
    float* cpu_total_weights,
#endif

    float *vals,
    const unsigned int h_n_splines);


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

extern void InitGPU_Segments(short int **segment);
extern void InitGPU_Vals(float **vals);
extern void CleanupGPU_Segments(short int *segment, float *vals);
#endif


// *****************************************
//Set everything to NULL or 0
void SMonolith::Initialise() {
// *****************************************

#ifdef CUDA
  MACH3LOG_INFO("Using GPU version event by event monolith");
#endif
  //KS: If true it will save spline monolith into huge ROOT file
  SaveSplineFile = false;

  nKnots = 0;
  NEvents = 0;
  _max_knots = 0;
  nParams = 0;
  NSplines_total = 0;
  NSplines_valid = 0;
  NSplines_total_large = 0;

  index_cpu = NULL;
  cpu_weights_var = NULL;
  gpu_weights = nullptr;

  cpu_total_weights = NULL;
  gpu_total_weights = nullptr;
  gpu_nParamPerEvent = nullptr;
  gpu_nPoints_arr = nullptr;
  gpu_paramNo_arr = nullptr;
  gpu_nKnots_arr = nullptr;
  gpu_coeff_x = nullptr;
  gpu_coeff_many = nullptr;
  
  SplineInfoArray = NULL;
  segments = NULL;
  vals = NULL;
  
  return;
}


// *****************************************
// Uses an x array and one combined yabd array
// This should optimise cache hitting because we use the same yabd points once we've found the x point
// So make these yabd points lay right next to each other in memory
SMonolith::SMonolith(std::vector<std::vector<TSpline3*> > &MasterSpline)
          : SplineBase() {
// *****************************************

  Initialise();
  MACH3LOG_INFO("Using full TSpline3, about to reduce it and send to GPU");
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TSpline3_red*> > ReducedSpline = ReduceTSpline3(MasterSpline);
  PrepareForGPU(ReducedSpline);
}

// *****************************************
// Constructor for the reduced TSpline3 object
SMonolith::SMonolith(std::vector<std::vector<TSpline3_red*> > &MasterSpline)
          : SplineBase() {
// *****************************************
  Initialise();
  MACH3LOG_INFO("-- GPUING WITH {X} and {Y,B,C,D} arrays and master spline containing TSpline3_red");
  PrepareForGPU(MasterSpline);
}

// *****************************************
// Uses a fifth order polynomial for most shape except 2p2h shape C/O which are two superimposed linear eq
// Reduce first
SMonolith::SMonolith(std::vector<std::vector<TF1*> > &MasterSpline)
          : SplineBase() {
// *****************************************

  Initialise();
  MACH3LOG_INFO("Using full TF1, about to reduce it and send to GPU");
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  std::vector<std::vector<TF1_red*> > ReducedSpline = ReduceTF1(MasterSpline);
  PrepareForGPU(ReducedSpline);
}

// *****************************************
// Uses a fifth order polynomial for most shape except 2p2h shape C/O which are two superimposed linear eq
// Reduce first
SMonolith::SMonolith(std::vector<std::vector<TF1_red*> > &MasterSpline)
          : SplineBase() {
// *****************************************
  Initialise();
  MACH3LOG_INFO("-- GPUING WITH TF1_red");
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  PrepareForGPU(MasterSpline);
}


// *****************************************
// Load SplineFile
SMonolith::SMonolith(std::string FileName)
          : SplineBase() {
// *****************************************
  Initialise();
  MACH3LOG_INFO("-- GPUING WITH {X} and {Y,B,C,D} arrays and master spline containing TSpline3_red");
  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  LoadSplineFile(FileName);
}

// *****************************************
// The shared initialiser from constructors of TSpline3 and TSpline3_red
void SMonolith::PrepareForGPU(std::vector<std::vector<TSpline3_red*> > &MasterSpline) {
// *****************************************

  // Scan for the max number of knots, the number of events (number of splines), and number of parameters
  int nSplines = 0;
  ScanMasterSpline(MasterSpline, NEvents, _max_knots, nParams, nSplines, nKnots);
  MACH3LOG_INFO("Found {} events", NEvents);
  MACH3LOG_INFO("Found {} knots at max", _max_knots);
  MACH3LOG_INFO("Found {} parameters", nParams);
  MACH3LOG_INFO("Found {} maximum number of splines in an event", nSplines);
  MACH3LOG_INFO("Found total {} knots in all splines", nKnots);

  // Can pass the spline segments to the GPU instead of the values
  // Make these here and only refill them for each loop, avoiding unnecessary new/delete on each reconfigure
  //KS: Since we are going to copy it each step use fancy CUDA memory allocation
#ifdef CUDA
  InitGPU_Segments(&segments);
  InitGPU_Vals(&vals);
#else
  segments = new short int[nParams]();
  vals = new float[nParams]();
#endif

  for (_int_ j = 0; j < nParams; j++)
  {
    segments[j] = 0;
    vals[j] = -999;
  }
  // Total number of events in our Spline, read from TSpline3 entries
  // Number of TSpline3 we have in total if each event had the maximal number of splines (nSplines written by ScanMasterSpline)
  NSplines_total = NEvents * nSplines;
  // Number of TSpline3 we have in total if each event has *EVERY* spline. Needed for some arrays
  NSplines_total_large = NEvents*nParams;

  unsigned int event_size_max = _max_knots * nParams;
  // Declare the {x}, {y,b,c,d} arrays for all possible splines which the event has
  // We'll filter off the flat and "disabled" (e.g. CCQE event should not have MARES spline) ones in the next for loop, but need to declare these beasts here

  // Declare the {y,b,c,d} for each knot
  // float because GPU precision (could change to double, but will incur significant speed reduction on GPU unless you're very rich!)
  cpu_coeff_many.resize(nKnots*_nCoeff_); // *4 because we store y,b,c,d parameters in this array
  //KS: For x coeff we assume that for given dial (MAQE) spacing is identical, here we are sloppy and assume each dial has the same number of knots, not a big problem
  cpu_coeff_x.resize(event_size_max);

  // Set all the big arrays to -999 to keep us safe...
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int j = 0; j < nKnots; j++) {
    for (int k = 0; k < _nCoeff_; k++) {
      cpu_coeff_many[j*_nCoeff_+k] = -999;
    }
  }

  for (unsigned int j = 0; j < event_size_max; j++) {
    cpu_coeff_x[j] = -999;
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
  //KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
  cpu_nParamPerEvent.resize(2*NSplines_total);
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
//KS: Gain from nowait is minimal mostly done as I am learning how to better utilise OMP
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
} //end parallel region
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
        cpu_coeff_x[ParamNumber*_max_knots + j] = x_tmp[j];
      }
      //KS: Contrary to X coeff we keep for other coeff only filled knots, there is no much gain for doing so for x coeff
      for (int j = 0; j < nPoints_tmp; ++j) {
        for (int k = 0; k < _nCoeff_; k++) {
          cpu_coeff_many[KnotCounter*_nCoeff_ + j*_nCoeff_ + k] = many_tmp[j*_nCoeff_+k];
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
      // Increment the counter for the number of good splines we have
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

  MACH3LOG_INFO("Number of splines = {}", NSplines_valid);

  // Make array with the number of points per spline (not per spline point!)
  cpu_paramNo_arr.resize(NSplines_valid);
  //KS: And array which tells where each spline stars in a big monolith array, sort of knot map
  cpu_nKnots_arr.resize(NSplines_valid);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < NSplines_valid; i++) {
    cpu_paramNo_arr[i] = paramNo_big[i];
    cpu_nKnots_arr[i] = knotNo_big[i];
    // Perform checks that all array entries have been changed from negative numbers and no number of points is greater than max knots, inputted by user
    // Don't do this for the index array since some entries there should be -1 so we know what splines to include and not include in each event for loading onto the GPU
    if (cpu_paramNo_arr[i] < 0) {
      std::cerr << "***** NEGATIVE PARAMETER NUMBER!!! ***** \n" << "On spline " << i << " " << cpu_paramNo_arr[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
  }
  delete[] paramNo_big;
  delete[] knotNo_big;

  int BadXCounter = 0;
  for (unsigned int j = 0; j < event_size_max; j++) {
    if (cpu_coeff_x[j] == -999) BadXCounter++;
    // Perform checks that all entries have been modified from initial values
    if (cpu_coeff_x[j] == -999 && BadXCounter < 5) {
      MACH3LOG_WARN("***** BAD X !! *****");
      MACH3LOG_WARN("Indicates some parameter doesn't have a single spline");
      MACH3LOG_WARN("j = {}");
      //throw;
    }
    if(BadXCounter == 5) MACH3LOG_WARN("There is more unutilised knots although I will stop spamming");
  }
  MACH3LOG_WARN("Found in total {} BAD X", BadXCounter);
#ifdef Weight_On_SplineBySpline_Basis
  // Make the array that holds all the returned weights from the GPU to pass to the CPU
  cpu_weights_var = new float[NSplines_valid]();
#else
  //KS: This is tricky as this variable use both by CPU and GPU, however if use CUDA we use cudaMallocHost
  #ifndef CUDA
  cpu_total_weights = new float[NEvents]();
  cpu_weights_var = new float[NSplines_valid]();
  #endif
#endif
  // Print some info; could probably make this to a separate function
  MACH3LOG_INFO("--- INITIALISED (X), (YBCD) ARRAYS ---");
  MACH3LOG_INFO("{} events with {} splines", NEvents, NSplines_valid);
  MACH3LOG_INFO("On average {:.2f} splines per event ({}/{})", float(NSplines_valid)/float(NEvents), NSplines_valid, NEvents);
  MACH3LOG_INFO("Size of x array = {:.4f} MB", double(sizeof(float)*event_size_max)/1.E6);
  MACH3LOG_INFO("Size of coefficient (y,b,c,d) array = {:.2f} MB", double(sizeof(float)*nKnots*_nCoeff_)/1.E6);
  MACH3LOG_INFO("Size of parameter # array = {:.2f} MB", double(sizeof(short int)*NSplines_valid)/1.E6);

  if(SaveSplineFile) PrepareSplineFile();

  PrepareForGPU_TSpline3();
}

// *****************************************
// Load SplineMonolith from ROOT file
void SMonolith::LoadSplineFile(std::string FileName) {
// *****************************************

  #ifdef Weight_On_SplineBySpline_Basis
  MACH3LOG_ERROR("Trying to load Monolith from file using weight by weight base, this is not supported right now, sorry");
  throw;
  #endif

  if (std::getenv("MACH3") != NULL) {
      FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
   }

  TFile *SplineFile = new TFile(FileName.c_str(), "OPEN");
  TTree *Settings = (TTree*)SplineFile->Get("Settings");
  TTree *Monolith = (TTree*)SplineFile->Get("Monolith");
  TTree *ParamInfo = (TTree*)SplineFile->Get("ParamInfo");
  TTree *XKnots = (TTree*)SplineFile->Get("XKnots");
  TTree *EventInfo = (TTree*)SplineFile->Get("EventInfo");
  TTree *FastSplineInfoTree = (TTree*)SplineFile->Get("FastSplineInfoTree");

  unsigned int NEvents_temp;
  short int nParams_temp;
  int _max_knots_temp;
  unsigned int nKnots_temp;
  unsigned int NSplines_valid_temp;

  Settings->SetBranchAddress("NEvents", &NEvents_temp);
  Settings->SetBranchAddress("nParams", &nParams_temp);
  Settings->SetBranchAddress("_max_knots", &_max_knots_temp);
  Settings->SetBranchAddress("nKnots", &nKnots_temp);
  Settings->SetBranchAddress("NSplines_valid", &NSplines_valid_temp);

  Settings->GetEntry(0);

  NEvents = NEvents_temp;
  nParams = nParams_temp;
  _max_knots = _max_knots_temp;
  nKnots = nKnots_temp;
  NSplines_valid = NSplines_valid_temp;

  unsigned int event_size_max = _max_knots * nParams;
  NSplines_total = NEvents * nParams;

  //KS: Since we are going to copy it each step use fancy CUDA memory allocation
#ifdef CUDA
  InitGPU_Segments(&segments);
  InitGPU_Vals(&vals);
#else
  segments = new short int[nParams]();
  vals = new float[nParams]();
#endif

  cpu_nParamPerEvent.resize(2*NSplines_total);
  cpu_paramNo_arr.resize(NSplines_valid);
  //KS: And array which tells where each spline stars in a big monolith array, sort of knot map
  cpu_nKnots_arr.resize(NSplines_valid);

  cpu_coeff_many.resize(nKnots*_nCoeff_); // *4 because we store y,b,c,d parameters in this array
  cpu_coeff_x.resize(event_size_max);

  //KS: This is tricky as this variable use both by CPU and GPU, however if use CUDA we use cudaMallocHost
#ifndef CUDA
  cpu_total_weights = new float[NEvents]();
  cpu_weights_var = new float[NSplines_valid]();
#endif

  float coeff = 0.;
  Monolith->SetBranchAddress("cpu_coeff_many", &coeff);
  for(unsigned int i = 0; i < nKnots*_nCoeff_; i++)
  {
    Monolith->GetEntry(i);
    cpu_coeff_many[i] = coeff;
  }

  short int paramNo_arr = 0;
  unsigned int nKnots_arr = 0;
  ParamInfo->SetBranchAddress("cpu_paramNo_arr", &paramNo_arr);
  ParamInfo->SetBranchAddress("cpu_nKnots_arr", &nKnots_arr);
  for(unsigned int i = 0; i < NSplines_valid; i++)
  {
    ParamInfo->GetEntry(i);
    cpu_paramNo_arr[i] = paramNo_arr;
    cpu_nKnots_arr[i] = nKnots_arr;
  }

  float coeff_x = 0.;
  XKnots->SetBranchAddress("cpu_coeff_x", &coeff_x);
  for(unsigned int i = 0; i < event_size_max; i++)
  {
    XKnots->GetEntry(i);
    cpu_coeff_x[i] = coeff_x;
  }

  unsigned int nParamPerEvent = 0;
  EventInfo->SetBranchAddress("cpu_nParamPerEvent", &nParamPerEvent);
  for(unsigned int i = 0; i < 2*NSplines_total; i++)
  {
    EventInfo->GetEntry(i);
    cpu_nParamPerEvent[i] = nParamPerEvent;
  }

  _int_ nPoints = 0;
  float xtemp[20];
  FastSplineInfoTree->SetBranchAddress("nPts", &nPoints);
  FastSplineInfoTree->SetBranchAddress("xPts", &xtemp);

  SplineInfoArray = new FastSplineInfo[nParams];
  for (_int_ i = 0; i < nParams; ++i) {
    FastSplineInfoTree->GetEntry(i);
    SplineInfoArray[i].nPts = -999;
    SplineInfoArray[i].xPts = NULL;
    SplineInfoArray[i].CurrSegment = 0;
    SplineInfoArray[i].splineParsPointer = NULL;

    // Fill the number of points
    SplineInfoArray[i].nPts = nPoints;
    if(nPoints == -999) continue;
    SplineInfoArray[i].xPts = new _float_[SplineInfoArray[i].nPts];
    for (_int_ k = 0; k < SplineInfoArray[i].nPts; ++k)
    {
      SplineInfoArray[i].xPts[k] = xtemp[k];
    }
  }

  SplineFile->Close();
  delete SplineFile;

  // Print some info; could probably make this to a separate function
  MACH3LOG_INFO("--- INITIALISED (X), (YBCD) ARRAYS ---");
  MACH3LOG_INFO("{} events with {} splines", NEvents, NSplines_valid);
  MACH3LOG_INFO("On average {:.2f} splines per event ({}/{})", float(NSplines_valid)/float(NEvents), NSplines_valid, NEvents);
  MACH3LOG_INFO("Size of x array = {:.4f} MB", double(sizeof(float)*event_size_max)/1.E6);
  MACH3LOG_INFO("Size of coefficient (y,b,c,d) array = {:.2f} MB", double(sizeof(float)*nKnots*_nCoeff_)/1.E6);
  MACH3LOG_INFO("Size of parameter # array = {:.2f} MB", double(sizeof(short int)*NSplines_valid)/1.E6);

  PrepareForGPU_TSpline3();
}

// *****************************************
// Save SplineMonolith into ROOT file
void SMonolith::PrepareSplineFile() {
// *****************************************

  std::string FileName = "inputs/SplineFile.root";
  if (std::getenv("MACH3") != NULL) {
      FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
   }

  TFile *SplineFile = new TFile(FileName.c_str(), "recreate");
  TTree *Settings = new TTree("Settings", "Settings");
  TTree *Monolith = new TTree("Monolith", "Monolith");
  TTree *ParamInfo = new TTree("ParamInfo", "ParamInfo");
  TTree *XKnots = new TTree("XKnots", "XKnots");
  TTree *EventInfo = new TTree("EventInfo", "EventInfo");
  TTree *FastSplineInfoTree = new TTree("FastSplineInfoTree", "FastSplineInfoTree");

  unsigned int NEvents_temp = NEvents;
  short int nParams_temp = nParams;
  int _max_knots_temp = _max_knots;
  unsigned int nKnots_temp = nKnots;
  unsigned int NSplines_valid_temp = NSplines_valid;

  Settings->Branch("NEvents", &NEvents_temp, "NEvents/i");
  Settings->Branch("nParams", &nParams_temp, "nParams/S");
  Settings->Branch("_max_knots", &_max_knots_temp, "_max_knots/I");
  Settings->Branch("nKnots", &nKnots_temp, "nKnots/i");
  Settings->Branch("NSplines_valid", &NSplines_valid_temp, "NSplines_valid/i");

  Settings->Fill();

  SplineFile->cd();
  Settings->Write();

  float coeff = 0.;
  Monolith->Branch("cpu_coeff_many", &coeff, "cpu_coeff_many/F");
  for(unsigned int i = 0; i < nKnots*_nCoeff_; i++)
  {
    coeff = cpu_coeff_many[i];
    Monolith->Fill();
  }
  SplineFile->cd();
  Monolith->Write();

  short int paramNo_arr = 0;
  unsigned int nKnots_arr = 0;
  ParamInfo->Branch("cpu_paramNo_arr", &paramNo_arr, "cpu_paramNo_arr/S");
  ParamInfo->Branch("cpu_nKnots_arr", &nKnots_arr, "cpu_nKnots_arr/i");
  for(unsigned int i = 0; i < NSplines_valid; i++)
  {
    paramNo_arr = cpu_paramNo_arr[i];
    nKnots_arr = cpu_nKnots_arr[i];

    ParamInfo->Fill();
  }
  SplineFile->cd();
  ParamInfo->Write();

  unsigned int event_size_max = _max_knots * nParams;

  float coeff_x = 0.;
  XKnots->Branch("cpu_coeff_x", &coeff_x, "cpu_coeff_x/F");
  for(unsigned int i = 0; i < event_size_max; i++)
  {
    coeff_x = cpu_coeff_x[i];
    XKnots->Fill();
  }
  SplineFile->cd();
  XKnots->Write();

  unsigned int nParamPerEvent = 0;
  EventInfo->Branch("cpu_nParamPerEvent", &nParamPerEvent, "cpu_nParamPerEvent/i");
  for(unsigned int i = 0; i < 2*NSplines_total; i++)
  {
    nParamPerEvent = cpu_nParamPerEvent[i];

    EventInfo->Fill();
  }
  SplineFile->cd();
  EventInfo->Write();

  _int_ nPoints = 0;
  float xtemp[20];
  FastSplineInfoTree->Branch("nPts", &nPoints, "nPts/I");
  FastSplineInfoTree->Branch("xPts", xtemp, "xPts[nPts]/F");

  for (_int_ i = 0; i < nParams; ++i)
  {
    nPoints = SplineInfoArray[i].nPts;

    for (_int_ k = 0; k < SplineInfoArray[i].nPts; ++k)
    {
      xtemp[k] = SplineInfoArray[i].xPts[k];
    }
    FastSplineInfoTree->Fill();
  }

  SplineFile->cd();
  FastSplineInfoTree->Write();

  delete Settings;
  delete Monolith;
  delete ParamInfo;
  delete XKnots;
  delete EventInfo;
  delete FastSplineInfoTree;
  SplineFile->Close();
  delete SplineFile;

}

// *****************************************
// The shared initialiser from constructors of TSpline3 and TSpline3_red
void SMonolith::PrepareForGPU_TSpline3() {
// *****************************************
  #ifdef CUDA
  unsigned int event_size_max = _max_knots * nParams;
  MACH3LOG_INFO("Total size = {:.2f} MB memory on CPU to move to GPU",
              (double(sizeof(float) * nKnots * _nCoeff_) + double(sizeof(float) * event_size_max) / 1.E6 +
              double(sizeof(short int) * NSplines_valid)) / 1.E6);
  MACH3LOG_INFO("GPU weight array (GPU->CPU every step) = {:.2f} MB", double(sizeof(float) * NSplines_valid) / 1.E6);
  #ifndef Weight_On_SplineBySpline_Basis
  MACH3LOG_INFO("Since you are running Total event weight mode then GPU weight array (GPU->CPU every step) = {:.2f} MB",
              double(sizeof(float) * NEvents) / 1.E6);
  #endif
  MACH3LOG_INFO("Parameter value array (CPU->GPU every step) = {:.4f} MB", double(sizeof(float) * nParams) / 1.E6);


  //CW: With the new set-up we have:   1 coefficient array of size coeff_array_size, all same size
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
      event_size_max //Knots times event number of unique splines
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

      cpu_paramNo_arr,
      cpu_nKnots_arr,
      cpu_coeff_x,
      cpu_coeff_many,
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
  cpu_coeff_x.clear();
  cpu_coeff_x.shrink_to_fit();
  cpu_coeff_many.clear();
  cpu_coeff_many.shrink_to_fit();
  cpu_paramNo_arr.clear();
  cpu_paramNo_arr.shrink_to_fit();
  cpu_nKnots_arr.clear();
  cpu_nKnots_arr.shrink_to_fit();
  #ifndef Weight_On_SplineBySpline_Basis
  cpu_nParamPerEvent.clear();
  cpu_nParamPerEvent.shrink_to_fit();
  #endif
  MACH3LOG_INFO("Good GPU loading");
#endif
  return;
}

// *****************************************
// The shared initialiser from constructors of TF1 and TF1_red
void SMonolith::PrepareForGPU(std::vector<std::vector<TF1_red*> > &MasterSpline) {
// *****************************************

  // Scan for the max number of knots, the number of events (number of splines), and number of parameters
  ScanMasterSpline(MasterSpline, NEvents, _max_knots, nParams);
  MACH3LOG_INFO("Found {} events", NEvents);
  MACH3LOG_INFO("Found {} polynomial at max", _max_knots);
  MACH3LOG_INFO("Found {} parameters", nParams);

  // Can pass the spline segments to the GPU instead of the values
  // Make these here and only refill them for each loop, avoiding unnecessary new/delete on each reconfigure
  //KS: Since we are going to copy it each step use fancy CUDA memory allocation
#ifdef CUDA
  InitGPU_Vals(&vals);
#else
  vals = new float[nParams]();
#endif
  for (_int_ j = 0; j < nParams; j++)
  {
    vals[j] = -999;
  }

  // Every event maximally has nParams TF1s which we've saved
  NSplines_total = NEvents * nParams;

  //CW: With TF1 we only save the coefficients and the order of the polynomial
  // Makes most sense to have one large monolithic array, but then it becomes impossible to tell apart a coefficient from a "number of points". So have two arrays: one of coefficients and one of number of points
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
  //KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
  cpu_nParamPerEvent.resize(2*NSplines_total);
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
      // Increment the counter for the number of good splines we have
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
  MACH3LOG_INFO("Number of splines = {}", NSplines_valid);

  // Now declare the arrays for each point in the valid splines which the event actually has (i.e. include the splines that the event undergoes)
  // Also make array with the number of points per spline (not per spline point!)
  // float because GPU precision (could change to double, but will incur significant speed reduction on GPU unless you're very rich!)
  cpu_nPoints_arr.resize(NSplines_valid);
  cpu_coeff_many.resize(NSplines_valid*_nTF1Coeff_); // *5 because this array holds  a,b,c,d,e parameters
  cpu_paramNo_arr.resize(NSplines_valid);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < NSplines_valid; i++) {
    cpu_nPoints_arr[i] = nPoints_big[i];
    cpu_paramNo_arr[i] = paramNo_big[i];
    for (unsigned int j = 0; j < _nTF1Coeff_; ++j) {
      cpu_coeff_many[i*_nTF1Coeff_+j] = coeffs_big[i*_nTF1Coeff_+j];
      if (cpu_coeff_many[i*_nTF1Coeff_+j] == -999) {
        std::cerr << "***** BAD COEFFICIENT OF POLY!!! ***** \n" << "On spline " << i << " = " << cpu_coeff_many[i*_nTF1Coeff_+j] << std::endl;
        std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }

    // Perform checks that all array entries have been changed from negative numbers and no number of points is greater than max knots, inputted by user
    // Don't do this for the index array since some entries there should be -1 so we know what splines to include and not include in each event for loading onto the GPU
    if (cpu_paramNo_arr[i] < 0) {
      std::cerr << "***** NEGATIVE PARAMETER NUMBER!!! ***** \n" << "On spline " << i << " = " << cpu_paramNo_arr[i] << std::endl;
      std::cerr << "Indicates bad reading and stripping back of splines pre-GPU" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    if (cpu_nPoints_arr[i] < 0) {
      std::cerr << "***** NEGATIVE NUMBER OF POINTS!!! ***** \n" << "On spline " << i << " = " << cpu_nPoints_arr[i] << std::endl;
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
  MACH3LOG_INFO("--- INITIALISED TF1 ARRAYS ---");
  std::cout << "  " << NEvents << " events with " << NSplines_valid << " splines" << std::endl;

  std::cout << "  On average " << float(NSplines_valid)/float(NEvents) << " splines per event (" << NSplines_valid << "/" << NEvents << ")" << std::endl;
  std::cout << "  Size of coefficient {a,b,c,d,e} array = " << double(sizeof(float)*NSplines_valid*_nTF1Coeff_)/1.E6 << " MB" << std::endl;
  std::cout << "  Size of parameter # array = " << double(sizeof(short int)*NSplines_valid)/1.E6 << " MB" << std::endl;
  std::cout << "  Size of polynomial type array = " << double(sizeof(short int)*NSplines_valid)/1.E6 << " MB" << std::endl;

  #ifdef Weight_On_SplineBySpline_Basis
  // Make the array that holds all the returned weights from the GPU to pass to the CPU
  cpu_weights_var = new float[NSplines_valid]();
#else
  //KS: This is tricky as this variable use both by CPU and GPU, however if use CUDA we use cudaMallocHost
  #ifndef CUDA
  cpu_total_weights = new float[NEvents]();
  cpu_weights_var = new float[NSplines_valid]();
  #endif
#endif

  PrepareForGPU_TF1();
}


// *****************************************
// The shared initialiser from constructors of TF1 and TF1_red
void SMonolith::PrepareForGPU_TF1() {
// *****************************************

#ifdef CUDA
  std::cout << "  Total size = " << (double(sizeof(float)*NSplines_valid*_nTF1Coeff_)+double(2.0*sizeof(short int)*NSplines_valid))/1.E6 << " MB memory on CPU to move to GPU" << std::endl;
  std::cout << "  GPU weight array (GPU->CPU every step) = " << double(sizeof(float)*NSplines_valid)/1.E6 << " MB" << std::endl;
  std::cout << "  Parameter value array (CPU->GPU every step) = " << double(sizeof(float)*nParams)/1.E6 << " MB" << std::endl;
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

      cpu_coeff_many,
      cpu_paramNo_arr,
      cpu_nPoints_arr,
#ifndef Weight_On_SplineBySpline_Basis
      NEvents,
      cpu_nParamPerEvent,
      gpu_nParamPerEvent,
#endif
      nParams,
      NSplines_valid,
      _max_knots);

  // Delete all the coefficient arrays from the CPU once they are on the GPU
  cpu_coeff_many.clear();
  cpu_coeff_many.shrink_to_fit();
  cpu_nPoints_arr.clear();
  cpu_nPoints_arr.shrink_to_fit();
  cpu_paramNo_arr.clear();
  cpu_paramNo_arr.shrink_to_fit();
  #ifndef Weight_On_SplineBySpline_Basis
  cpu_nParamPerEvent.clear();
  cpu_nParamPerEvent.shrink_to_fit();
  #endif
  MACH3LOG_INFO("Good TF1 GPU loading");
#endif

  return;
}

// Need to specify template functions in header
// *****************************************
// Scan the master spline to get the maximum number of knots in any of the TSpline3*
void SMonolith::ScanMasterSpline(std::vector<std::vector<TSpline3_red*> > & MasterSpline, unsigned int &nEvents, int &MaxPoints, short int &numParams, int &nSplines, unsigned int &numKnots) {
// *****************************************

  // Need to extract: the total number of events
  //                  number of parameters
  //                  maximum number of knots
  MaxPoints = 0;
  nEvents   = 0;
  numParams   = 0;
  nSplines = 0;
  numKnots = 0;
  std::vector<std::vector<TSpline3_red*> >::iterator OuterIt;
  std::vector<TSpline3_red*>::iterator InnerIt;

  // Check the number of events
  nEvents = MasterSpline.size();

  // Maximum number of splines one event can have (scan through and find this number)
  int nMaxSplines_PerEvent = 0;
  
  //KS: We later check that each event has the same number of splines so this is fine
  numParams = MasterSpline[0].size();
  // Initialise
  SplineInfoArray = new FastSplineInfo[numParams];
  for (_int_ i = 0; i < numParams; ++i) {
    SplineInfoArray[i].nPts = -999;
    SplineInfoArray[i].xPts = NULL;
    SplineInfoArray[i].CurrSegment = 0;
    SplineInfoArray[i].splineParsPointer = NULL;
  }

  unsigned int EventCounter = 0;
  // Loop over each parameter
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt) {
    // Check that each event has each spline saved
    if (numParams > 0) {
      int TempSize = (*OuterIt).size();
      if (TempSize != numParams) {
        MACH3LOG_ERROR("Found {} parameters for event {}", TempSize, EventCounter);
        MACH3LOG_ERROR("but was expecting {} since that's what I found for the previous event", numParams);
        MACH3LOG_ERROR("Somehow this event has a different number of spline parameters... Please study further!");
        throw;
      }
    }
    numParams = (*OuterIt).size();

    int nSplines_SingleEvent = 0;
    // Loop over each pointer
    int ij = 0;
    for (InnerIt = OuterIt->begin(); InnerIt != OuterIt->end(); ++InnerIt, ij++) {
      if ((*InnerIt) == NULL) continue;
      int nPoints = (*InnerIt)->GetNp();
      if (nPoints > MaxPoints) {
        MaxPoints = nPoints;
      }
      numKnots += nPoints;
      nSplines_SingleEvent++;
      
      // Fill the SplineInfoArray entries with information on each splinified parameter
      if (SplineInfoArray[ij].xPts == NULL)
      {
        // Fill the number of points
        SplineInfoArray[ij].nPts = (*InnerIt)->GetNp();

        // Fill the x points
        SplineInfoArray[ij].xPts = new _float_[SplineInfoArray[ij].nPts];
        for (_int_ k = 0; k < SplineInfoArray[ij].nPts; ++k)
        {
          _float_ xtemp = -999.99;
          _float_ ytemp = -999.99;
          (*InnerIt)->GetKnot(k, xtemp, ytemp);
          SplineInfoArray[ij].xPts[k] = xtemp;
        }
      }
    }

    if (nSplines_SingleEvent > nMaxSplines_PerEvent) nMaxSplines_PerEvent = nSplines_SingleEvent;
    EventCounter++;
  }
  nSplines = nMaxSplines_PerEvent;
  
  int Counter = 0;
  //KS: Sanity check that everything was set correctly
  for (_int_ i = 0; i < numParams; ++i)
  {
    const _int_ nPoints = SplineInfoArray[i].nPts;
    const _float_* xArray = SplineInfoArray[i].xPts;
    if (nPoints == -999 || xArray == NULL) {
      Counter++;
      if(Counter < 5)
      {
        MACH3LOG_WARN("SplineInfoArray[{}] isn't set yet", i);
      }
      continue;
      //throw;
    }
  }
  MACH3LOG_WARN("In total SplineInfoArray for {} hasn't been initialised", Counter);
}

// Need to specify template functions in header
// *****************************************
// Scan the master spline to get the maximum number of knots in any of the TSpline3*
void SMonolith::ScanMasterSpline(std::vector<std::vector<TF1_red*> > & MasterSpline, unsigned int &nEvents, int &MaxPoints, short int &numParams) {
// *****************************************

  // Need to extract: the total number of events
  //                  number of parameters
  //                  maximum number of knots
  MaxPoints = 0;
  nEvents   = 0;
  numParams   = 0;

  std::vector<std::vector<TF1_red*> >::iterator OuterIt;
  std::vector<TF1_red*>::iterator InnerIt;

  // Check the number of events
  nEvents = MasterSpline.size();

  unsigned int EventCounter = 0;
  // Loop over each parameter
  for (OuterIt = MasterSpline.begin(); OuterIt != MasterSpline.end(); ++OuterIt) {
    // Check that each event has each spline saved
    if (numParams > 0) {
      int TempSize = (*OuterIt).size();
      if (TempSize != numParams) {
        MACH3LOG_ERROR("Found {} parameters for event {}", TempSize, EventCounter);
        MACH3LOG_ERROR("but was expecting {} since that's what I found for the previous event", numParams);
        MACH3LOG_ERROR("Somehow this event has a different number of spline parameters... Please study further!");
        throw;
      }
    }
    numParams = (*OuterIt).size();

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

#ifdef CUDA
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
  
  //KS: Since we declared them using CUDA alloc we have to free memory using also cuda functions
  CleanupGPU_Segments(segments, vals);
#else
  if(segments != NULL) delete[] segments;
  if(vals != NULL) delete[] vals;
  if(cpu_total_weights != NULL) delete[] cpu_total_weights;
#endif

  if(SplineInfoArray != NULL) delete[] SplineInfoArray;

#ifdef Weight_On_SplineBySpline_Basis
  if(cpu_weights != NULL) delete[] cpu_weights;
  if(cpu_weights_var != NULL) delete[] cpu_weights_var;
#endif
  if(index_cpu != NULL) delete[] index_cpu;

  //KS: Those might be deleted or not depending on GPU/CPU TSpline3/TF1 DEBUG or not hence we check if not NULL
  cpu_coeff_x.clear();
  cpu_coeff_x.shrink_to_fit();
  cpu_coeff_many.clear();
  cpu_coeff_many.shrink_to_fit();
  cpu_paramNo_arr.clear();
  cpu_paramNo_arr.shrink_to_fit();
  cpu_nKnots_arr.clear();
  cpu_nKnots_arr.shrink_to_fit();
  #ifndef Weight_On_SplineBySpline_Basis
  cpu_nParamPerEvent.clear();
  cpu_nParamPerEvent.shrink_to_fit();
  #endif
  cpu_nPoints_arr.clear();
  cpu_nPoints_arr.shrink_to_fit();
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
  // This is used later to expedite the calculations for flat splines
  // tmpArray[0] is number of knots
  if (isFlat(spl)) {
    nPoints = 1;
  } else {
    nPoints = Np;
    if (Np > _max_knots) {
      MACH3LOG_ERROR("Error, number of points is greater than saved {}", _max_knots);
      MACH3LOG_ERROR("This _WILL_ cause problems with GPU splines and _SHOULD_ be fixed!");
      MACH3LOG_ERROR("nPoints = {}, _max_knots = {}", nPoints, _max_knots);
      throw;
    }
  }

  // The coefficients we're writing to
  _float_ x, y, b, c, d;
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
  #ifndef DEBUG
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
  delete spl;
  spl = NULL;
}

// *****************************************
//CW: Check if the TSpline3 object is flat; if it is, we won't load it onto the GPU
bool SMonolith::isFlat(TSpline3_red* &spl) {
// *****************************************
  int Np = spl->GetNp();
  _float_ x, y, b, c, d;
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


#ifdef CUDA
// *****************************************
// Tell the GPU to evaluate the weights
// Load up the two x,{y,b,c,d} arrays into memory and have GPU read them with more coalescence instead of one monolithic array
// This should be used when we're using separate x,y,a,b,c,d arrays
// Also pass the segments for the parameter along with their parameter values
// This avoids doing lots of binary searches on the GPU
void SMonolith::Evaluate() {
// *****************************************

  // There's a parameter mapping that goes from spline parameter to a global parameter index
  // Find the spline segments
  FindSplineSegment();

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
      segments,
      NSplines_valid);

  //KS: Normally it does nothing, in case you want to have weight for each spline it does the mapping, used mostly for debugging
  ModifyWeights_GPU();
}
#else
//If CUDA is not enabled do the same on CPU
// *****************************************
void SMonolith::Evaluate() {
// *****************************************

  // There's a parameter mapping that goes from spline parameter to a global parameter index
  // Find the spline segments
  FindSplineSegment();

  //KS: Huge MP loop over all valid splines
  CalcSplineWeights();

  //KS: Huge MP loop over all events calculating total weight
  ModifyWeights();

  return;
}
#endif

// *************************
// CW: Only need to do the binary search once per parameter, not once per event!
// Takes down the number of binary searches from 1.2M to 17, ha!
// For ROOT version see root/hist/hist/src/TSpline3.cxx TSpline3::FindX(double)
void SMonolith::FindSplineSegment() {
// *************************

  // Loop over the splines
  //KS: Tried multithreading here with 48 splines and it is faster with one thread, maybe in future multithreading will be worth revisiting
  for (_int_ i = 0; i < nParams; ++i)
  {
    const _int_ nPoints = SplineInfoArray[i].nPts;
    const _float_* xArray = SplineInfoArray[i].xPts;

    // EM: if we have a parameter that has no response for any event (i.e. all splines have just one knot), then skip it and avoid a seg fault here
    //     In principle, such parameters shouldn't really be included in the first place, but with new det syst splines this
    //     could happen if say you were just running on one FHC run, then all RHC parameters would be flat and the code below would break.
    if(xArray == NULL) continue;

    // Get the variation for this reconfigure for the ith parameter
    const _float_ xvar = *SplineInfoArray[i].splineParsPointer;
    vals[i] = xvar;

    // The segment we're interested in (klow in ROOT code)
    _int_ segment = 0;
    _int_ kHigh = nPoints-1;
    //KS: We expect new segment is very close to previous
    const _int_ PreviousSegment = SplineInfoArray[i].CurrSegment;

    //KS: It is quite probable the new segment is same as in previous step so try to avoid binary search
    if( xArray[PreviousSegment+1] > xvar && xvar >= xArray[PreviousSegment] ) segment = PreviousSegment;
    // If the variation is below the lowest saved spline point
    else if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      //CW: Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search
    } else {
      // The top point we've got
      _int_ kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point

    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;

    //CW: This way we avoid doing 1.2M+ binary searches on the GPU
    // and literally just multiply lots of numbers together on the GPU without any algorithm
    // Update the values and which segment it belongs to
    SplineInfoArray[i].CurrSegment = segment;
    segments[i] = SplineInfoArray[i].CurrSegment;

#ifdef DEBUG
    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
      std::cerr << "Found a segment which is _ABOVE_ the variation!" << std::endl;
      std::cerr << "IT SHOULD ALWAYS BE BELOW! (except when segment 0)" << std::endl;
      std::cerr << "Spline: "<< i << std::endl;

      std::cerr << "Found segment   = " << segment << std::endl;
      std::cerr << "Doing variation = " << xvar << std::endl;
      std::cerr << "x in spline     = " << SplineInfoArray[i].xPts[segment] << std::endl;
      for (_int_ j = 0; j < SplineInfoArray[j].nPts; ++j) {
        std::cerr << "    " << j << " = " << SplineInfoArray[i].xPts[j] << std::endl;
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
#endif
  } //end loop over params
}

//*********************************************************
void SMonolith::CalcSplineWeights() {
//*********************************************************
  #ifdef MULTITHREAD
  #pragma omp parallel for simd
  #endif
  for (unsigned int splineNum = 0; splineNum < NSplines_valid; ++splineNum)
  {
    //CW: Which Parameter we are accessing
    const short int Param = cpu_paramNo_arr[splineNum];

    //CW: Avoids doing costly binary search on GPU
    const short int segment = segments[Param];

    //KS: Segment for coeff_x is simply parameter*max knots + segment as each parameters has the same spacing
    const short int segment_X = Param*_max_knots+segment;

    //KS: Find knot position in out monolithical structure
    const unsigned int CurrentKnotPos = cpu_nKnots_arr[splineNum]*_nCoeff_+segment*_nCoeff_;

    // We've read the segment straight from CPU and is saved in segment_gpu
    // polynomial parameters from the monolithic splineMonolith
    const float fY = cpu_coeff_many[CurrentKnotPos];
    const float fB = cpu_coeff_many[CurrentKnotPos+1];
    const float fC = cpu_coeff_many[CurrentKnotPos+2];
    const float fD = cpu_coeff_many[CurrentKnotPos+3];
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float dx = vals[Param] - cpu_coeff_x[segment_X];

    //CW: Wooow, let's use some fancy intrinsic and pull down the processing time by <1% from normal multiplication! HURRAY
    cpu_weights_var[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
    // Or for the more "easy to read" version:
    //cpu_weights_var[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));
  }
  return;
}


#ifdef CUDA
// *****************************************
// Tell the GPU to evaluate the weights
// TF1 version
void SMonolith::Evaluate_TF1() {
// *****************************************

  // Feed the parameter variations
  for (_int_ i = 0; i < nParams; ++i) {
    // Update the values and which segment it belongs to
    vals[i] = *splineParsPointer[i];
  }
  
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
      vals,
      NSplines_valid);

  //KS: Normally it does nothing, in case you want to have weight for each spline it does the mapping, used mostly for debugging
  ModifyWeights_GPU();
}
#else
//If CUDA is not enabled do the same on CPU
// *****************************************
void SMonolith::Evaluate_TF1() {
// *****************************************

  // Feed the parameter variations
  for (_int_ i = 0; i < nParams; ++i) {
    // Update the values and which segment it belongs to
    vals[i] = *splineParsPointer[i];
  }
  
  //KS: Huge MP loop over all valid splines
  CalcSplineWeights_TF1();

  //KS: Huge MP loop over all events calculating total weight
  ModifyWeights();

  return;
}
#endif

//*********************************************************
void SMonolith::CalcSplineWeights_TF1() {
//*********************************************************

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int splineNum = 0; splineNum < NSplines_valid; ++splineNum)
  {
    // The is the variation itself (needed to evaluate variation - stored spline point = dx)
    const float x = vals[cpu_paramNo_arr[splineNum]];

    // Read the coefficients
    const float a = cpu_coeff_many[splineNum*_max_knots];
    const float b = cpu_coeff_many[splineNum*_max_knots+1];
    const float c = cpu_coeff_many[splineNum*_max_knots+2];
    const float d = cpu_coeff_many[splineNum*_max_knots+3];
    const float e = cpu_coeff_many[splineNum*_max_knots+4];

    // Match these with form in SetSplines
    // Might not be great to have this if statement: maybe split two kernels?
    if (gpu_nPoints_arr[splineNum] == 5) {
      cpu_weights_var[splineNum] = 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x + e*x*x*x*x*x;
    } else if (gpu_nPoints_arr[splineNum] == 2) {
      cpu_weights_var[splineNum] = (x<=0)*(1+a*x) + (x>0)*(1+b*x);
    } else {
      printf("Big problems, I found a nPoints array which is not 5 or 2 on GPU!\n");
    }
  }
  return;
}

//*********************************************************
//KS: Calc total event weight on CPU
void SMonolith::ModifyWeights(){
//*********************************************************
#ifndef Weight_On_SplineBySpline_Basis
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int EventNum = 0; EventNum < NEvents; ++EventNum)
  {
    float totalWeight = 1.0f; // Initialize total weight for each event

    // Extract the parameters for the current event
    const unsigned int startIndex = cpu_nParamPerEvent[2 * EventNum + 1];
    const unsigned int numParams = cpu_nParamPerEvent[2 * EventNum];

    // Compute total weight for the current event
    for (unsigned int id = 0; id < numParams; ++id)
    {
      totalWeight *= cpu_weights_var[startIndex + id];
    }

    // Store the total weight for the current event
    cpu_total_weights[EventNum] = totalWeight;
  }
#else
  //KS: Name is confusing but what it does it make a nice mapping used for debugging
  ModifyWeights_GPU();
#endif
  return;
}

//*********************************************************
//KS: Normally it does nothing, in case you want to have weight for each spline it does the mapping, used mostly for debugging
void SMonolith::ModifyWeights_GPU(){
//*********************************************************

#ifdef Weight_On_SplineBySpline_Basis
  // Multi-thread here because _numIndex is really quite large!
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int i = 0; i < NSplines_total; ++i) {
    if (index_cpu[i] >= 0) {
      cpu_weights[i] = cpu_weights_var[index_cpu[i]];
    } else {
      cpu_weights[i] = 1.;
    }
  }
#endif
  return;
}
