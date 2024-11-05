#include "SplineMonolith.h"

#ifdef CUDA
#include "splines/gpuSplineUtils.cuh"
#endif

// *****************************************
//Set everything to NULL or 0
void SMonolith::Initialise() {
// *****************************************

#ifdef CUDA
  MACH3LOG_INFO("Using GPU version event by event monolith");
  gpu_spline_handler = nullptr;
#endif

#ifndef USE_FPGA
  cpu_spline_handler = new SplineMonoStruct();
#endif

  nKnots = 0;
  nTF1coeff = 0;
  NEvents = 0;
  _max_knots = 0;
  nParams = 0;

  NSplines_valid = 0;
  NTF1_valid = 0;
  NSplines_total_large = 0;

  index_cpu = nullptr;
  index_TF1_cpu = nullptr;
  cpu_weights_var = nullptr;
  cpu_weights = nullptr;
  cpu_weights_tf1_var = nullptr;

  cpu_total_weights = nullptr;

  SplineInfoArray = nullptr;
  segments = NULL;
  vals = NULL;
  
  return;
}

// *****************************************
SMonolith::SMonolith(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline, const std::vector<RespFuncType> &SplineType, const bool SaveFlatTree)
: SplineBase() {
// *****************************************

  //KS: If true it will save spline monolith into huge ROOT file
  SaveSplineFile = SaveFlatTree;
  Initialise();
  MACH3LOG_INFO("-- GPUING WITH arrays and master spline containing TResponseFunction_red");

  // Convert the TSpline3 pointers to the reduced form and call the reduced constructor
  PrepareForGPU(MasterSpline, SplineType);
}

// *****************************************
// The shared initialiser from constructors of TSpline3 and TSpline3_red
void SMonolith::PrepareForGPU(std::vector<std::vector<TResponseFunction_red*> > &MasterSpline, const std::vector<RespFuncType> &SplineType) {
// *****************************************

  // Scan for the max number of knots, the number of events (number of splines), and number of parameters
  int maxnSplines = 0;
  ScanMasterSpline(MasterSpline,
                   NEvents,
                   _max_knots,
                   nParams,
                   maxnSplines,
                   NSplines_valid,
                   nKnots,
                   NTF1_valid,
                   nTF1coeff,
                   SplineType);

  MACH3LOG_INFO("Found {} events", NEvents);
  MACH3LOG_INFO("Found {} knots at max", _max_knots);
  MACH3LOG_INFO("Found {} parameters", nParams);
  MACH3LOG_INFO("Found {} maximum number of splines in an event", maxnSplines);
  MACH3LOG_INFO("Found total {} knots in all splines", nKnots);
  MACH3LOG_INFO("Number of splines = {}", NSplines_valid);
  MACH3LOG_INFO("Found total {} coeffs in all TF1", nTF1coeff);
  MACH3LOG_INFO("Number of TF1 = {}", NTF1_valid);

  // Can pass the spline segments to the GPU instead of the values
  // Make these here and only refill them for each loop, avoiding unnecessary new/delete on each reconfigure
  //KS: Since we are going to copy it each step use fancy CUDA memory allocation
  #ifdef CUDA
  gpu_spline_handler->InitGPU_Segments(&segments);
  gpu_spline_handler->InitGPU_Vals(&vals);
  #else
  segments = new short int[nParams]();
  vals = new float[nParams]();
  #endif

  for (_int_ j = 0; j < nParams; j++)
  {
    segments[j] = 0;
    vals[j] = -999;
  }

  // Number of objects we have in total if each event has *EVERY* spline. Needed for some arrays
  NSplines_total_large = NEvents*nParams;

  unsigned int event_size_max = _max_knots * nParams;
  // Declare the {x}, {y,b,c,d} arrays for all possible splines which the event has
  // We'll filter off the flat and "disabled" (e.g. CCQE event should not have MARES spline) ones in the next for loop, but need to declare these beasts here

  // Declare the {y,b,c,d} for each knot
  // float because GPU precision (could change to double, but will incur significant speed reduction on GPU unless you're very rich!)
  cpu_spline_handler->coeff_many.resize(nKnots*_nCoeff_); // *4 because we store y,b,c,d parameters in this array
  //KS: For x coeff we assume that for given dial (MAQE) spacing is identical, here we are sloppy and assume each dial has the same number of knots, not a big problem
  cpu_spline_handler->coeff_x.resize(event_size_max);

  // Set all the big arrays to -999 to keep us safe...
  for (unsigned int j = 0; j < event_size_max; j++) {
    cpu_spline_handler->coeff_x[j] = -999;
  }

  //CW: With TF1 we only save the coefficients and the order of the polynomial
  // Makes most sense to have one large monolithic array, but then it becomes impossible to tell apart a coefficient from a "number of points". So have two arrays: one of coefficients and one of number of points
  // Let's first assume all are of _max_knots size
  // Now declare the arrays for each point in the valid splines which the event actually has (i.e. include the splines that the event undergoes)
  // Also make array with the number of points per spline (not per spline point!)
  // float because GPU precision (could change to double, but will incur significant speed reduction on GPU unless you're very rich!)
  cpu_nPoints_arr.resize(NTF1_valid);
  cpu_coeff_TF1_many.resize(nTF1coeff); // *5 because this array holds  a,b,c,d,e parameters

  #ifdef Weight_On_SplineBySpline_Basis
  // This holds the index of each spline
  index_cpu = new int[NSplines_total_large];
  index_TF1_cpu = new int[NSplines_total_large];

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int j = 0; j < NSplines_total_large; j++) {
    index_cpu[j] = -1;
    index_TF1_cpu[j] = -1;
  }
  // This holds the total CPU weights that gets read in samplePDFND
  cpu_weights = new float[NSplines_total_large];
  #else
  //KS: Map keeping track how many parameters applies to each event, we keep two numbers here {number of splines per event, index where splines start for a given event}
  cpu_nParamPerEvent.resize(2*NEvents);
  cpu_nParamPerEvent_tf1.resize(2*NEvents);
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int j = 0; j < 2*NEvents; j++) {
    cpu_nParamPerEvent[j] = -1;
    cpu_nParamPerEvent_tf1[j] = -1;
  }
  #endif

  // Make array with the number of points per spline (not per spline point!)
  cpu_spline_handler->paramNo_arr.resize(NSplines_valid);
  //KS: And array which tells where each spline stars in a big monolith array, sort of knot map
  cpu_spline_handler->nKnots_arr.resize(NSplines_valid);
  cpu_paramNo_TF1_arr.resize(NTF1_valid);

  // Temporary arrays to hold the coefficients for each spline
  // We get one x, one y, one b,... for each point, so only need to be _max_knots big
  //KS: Some params has less splines but this is all right main array will get proper number while this temp will be deleted
  float *x_tmp = new float[_max_knots]();
  float *many_tmp = new float[_max_knots*_nCoeff_]();
  float *temp_coeffs = new float[_nTF1Coeff_]();

  // Count the number of events
  unsigned int KnotCounter = 0;
  unsigned int TF1PointsCounter = 0;
  unsigned int NSplinesCounter = 0;
  unsigned int TF1sCounter = 0;
  int ParamCounter = 0;
  int ParamCounterGlobal = 0;
  int ParamCounter_TF1 = 0;
  int ParamCounterGlobalTF1 = 0;
  // Loop over events and extract the spline coefficients
  for(unsigned int EventCounter = 0; EventCounter < MasterSpline.size(); ++EventCounter) {
    // Structure of MasterSpline is std::vector<std::vector<TSpline3*>>
    // A conventional iterator to count which parameter a given spline should be applied to
    for(unsigned int ParamNumber = 0; ParamNumber < MasterSpline[EventCounter].size(); ++ParamNumber) {

      // If NULL we don't have this spline for the event, so move to next spline
      if (MasterSpline[EventCounter][ParamNumber] == NULL) continue;

      if(SplineType[ParamNumber] == kTSpline3_red)
      {
        //KS: how much knots each spline has
        int nPoints_tmp = 0;
        // Get a pointer to the current spline for this event
        TResponseFunction_red* TespFunc = MasterSpline[EventCounter][ParamNumber];
        TSpline3_red* CurrSpline = static_cast<TSpline3_red*>(TespFunc);

        // If the number of knots are greater than 2 the spline is not a dummy and we should extract coefficients to load onto the GPU
        getSplineCoeff_SepMany(CurrSpline, nPoints_tmp, x_tmp, many_tmp);

        //KS: One knot means flat spline so ignore
        if (nPoints_tmp == 1) continue;
        for (int j = 0; j < _max_knots; ++j) {
          cpu_spline_handler->coeff_x[ParamNumber*_max_knots + j] = x_tmp[j];
        }
        //KS: Contrary to X coeff we keep for other coeff only filled knots, there is no much gain for doing so for x coeff
        for (int j = 0; j < nPoints_tmp; ++j) {
          for (int k = 0; k < _nCoeff_; k++) {
            cpu_spline_handler->coeff_many[KnotCounter*_nCoeff_ + j*_nCoeff_ + k] = many_tmp[j*_nCoeff_+k];
          }
        }
        // Set the parameter number for this spline
        cpu_spline_handler->paramNo_arr[NSplinesCounter] = ParamNumber;
        //KS: Fill map when each spline starts
        cpu_spline_handler->nKnots_arr[NSplinesCounter] = KnotCounter;
        KnotCounter += nPoints_tmp;

        #ifdef Weight_On_SplineBySpline_Basis
        // Set the index of the spline so we can tell apart from flat splines
        index_cpu[EventCounter*nParams + ParamNumber] = NSplinesCounter;
        #else
        ++ParamCounter;
        #endif
        // Increment the counter for the number of good splines we have
        ++NSplinesCounter;
      }
      else if (SplineType[ParamNumber] == kTF1_red)
      {
        // Don't actually use this ever -- we give each spline the maximum number of points found in all splines
        int nPoints_tmp = 0;
        // Get a pointer to the current spline for this event
        TF1_red* CurrSpline = dynamic_cast<TF1_red*>(MasterSpline[EventCounter][ParamNumber]);

        // If the number of knots are greater than 2 the spline is not a dummy and we should extract coefficients to load onto the GPU
        getTF1Coeff(CurrSpline, nPoints_tmp, temp_coeffs);
        for (int j = 0; j < _nTF1Coeff_; ++j) {
          cpu_coeff_TF1_many[TF1PointsCounter+j] = temp_coeffs[j];
        }
        // Save the number of points for this spline
        cpu_nPoints_arr[TF1sCounter] = nPoints_tmp;

        TF1PointsCounter += nPoints_tmp;
        // Set the parameter number for this spline
        cpu_paramNo_TF1_arr[TF1sCounter] = ParamNumber;
        #ifdef Weight_On_SplineBySpline_Basis
        // Set the index of the spline so we can tell apart from flat splines
        index_TF1_cpu[EventCounter*nParams + ParamNumber] = TF1sCounter;
        #else
        ++ParamCounter_TF1;
        #endif
        // Increment the counter for the number of good splines we have
        ++TF1sCounter;
      }
      //KS: Don't delete in debug
      #ifndef DEBUG
      delete MasterSpline[EventCounter][ParamNumber];
      MasterSpline[EventCounter][ParamNumber] = NULL;
      #endif
    } // End the loop over the parameters in the MasterSpline
    #ifndef Weight_On_SplineBySpline_Basis
    cpu_nParamPerEvent[2*EventCounter] = ParamCounter;
    cpu_nParamPerEvent[2*EventCounter+1] = ParamCounterGlobal;
    ParamCounterGlobal += ParamCounter;

    cpu_nParamPerEvent_tf1[2*EventCounter] = ParamCounter_TF1;
    cpu_nParamPerEvent_tf1[2*EventCounter+1] = ParamCounterGlobalTF1;
    ParamCounterGlobalTF1 += ParamCounter_TF1;

    ParamCounter = 0;
    ParamCounter_TF1 = 0;
    #endif
  } // End the loop over the number of events
  delete[] many_tmp;
  delete[] x_tmp;
  delete[] temp_coeffs;

  int BadXCounter = 0;
  for (unsigned int j = 0; j < event_size_max; j++) {
    if (cpu_spline_handler->coeff_x[j] == -999) BadXCounter++;
    // Perform checks that all entries have been modified from initial values
    if (cpu_spline_handler->coeff_x[j] == -999 && BadXCounter < 5) {
      MACH3LOG_WARN("***** BAD X !! *****");
      MACH3LOG_WARN("Indicates some parameter doesn't have a single spline");
      MACH3LOG_WARN("j = {}", j);
      //throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    if(BadXCounter == 5) MACH3LOG_WARN("There is more unutilised knots although I will stop spamming");
  }

  MACH3LOG_WARN("Found in total {} BAD X", BadXCounter);
  #ifdef Weight_On_SplineBySpline_Basis
  // Make the array that holds all the returned weights from the GPU to pass to the CPU
  cpu_weights_var = new float[NSplines_valid]();
  cpu_weights_tf1_var = new float[NTF1_valid]();
  #else
    //KS: This is tricky as this variable use both by CPU and GPU, however if use CUDA we use cudaMallocHost
    #ifndef CUDA
    cpu_total_weights = new float[NEvents]();
    cpu_weights_var = new float[NSplines_valid]();
    cpu_weights_tf1_var = new float[NTF1_valid]();
    #endif
  #endif

  // Print some info; could probably make this to a separate function
  PrintInitialsiation();

  if(SaveSplineFile) PrepareSplineFile();

  MoveToGPU();
}

// *****************************************
// The shared initialiser from constructors of TSpline3 and TSpline3_red
void SMonolith::MoveToGPU() {
// *****************************************
  #ifdef CUDA
  unsigned int event_size_max = _max_knots * nParams;
  MACH3LOG_INFO("Total size = {:.2f} MB memory on CPU to move to GPU",
                (double(sizeof(float) * nKnots * _nCoeff_) + double(sizeof(float) * event_size_max) / 1.E6 +
                double(sizeof(short int) * NSplines_valid)) / 1.E6);
  MACH3LOG_INFO("Total TF1 size = {:.2f} MB memory on CPU to move to GPU",
                double(sizeof(float) * NTF1_valid * _nTF1Coeff_) / 1.E6);
  MACH3LOG_INFO("GPU weight array (GPU->CPU every step) = {:.2f} MB", double(sizeof(float) * (NSplines_valid + NTF1_valid) / 1.E6));
  #ifndef Weight_On_SplineBySpline_Basis
  MACH3LOG_INFO("Since you are running Total event weight mode then GPU weight array (GPU->CPU every step) = {:.2f} MB",
                double(sizeof(float) * NEvents) / 1.E6);
  #endif
  MACH3LOG_INFO("Parameter value array (CPU->GPU every step) = {:.4f} MB", double(sizeof(float) * nParams) / 1.E6);
  //CW: With the new set-up we have:   1 coefficient array of size coeff_array_size, all same size
  //                                1 coefficient array of size coeff_array_size*4, holding y,b,c,d in order (y11,b11,c11,d11; y12,b12,c12,d12;...) where ynm is n = spline number, m = spline point. Should really make array so that order is (y11,b11,c11,d11; y21,b21,c21,d21;...) because it will optimise cache hits I think; try this if you have time
  //                                return gpu_weights

  gpu_spline_handler = new SMonolithGPU();

  // The gpu_XY arrays don't actually need initialising, since they are only placeholders for what we'll move onto the GPU. As long as we cudaMalloc the size of the arrays correctly there shouldn't be any problems
  // Can probably make this a bit prettier but will do for now
  // Could be a lot smaller of a function...
  gpu_spline_handler->InitGPU_SplineMonolith(
          #ifndef Weight_On_SplineBySpline_Basis
          &cpu_total_weights,
          NEvents,
          #endif
          nKnots, // How many entries in coefficient array (*4 for the "many" array)
          NSplines_valid, // What's the number of splines we have (also number of entries in gpu_nPoints_arr)
          NTF1_valid,
          event_size_max //Knots times event number of unique splines
  );

  // Move number of splines and spline size to constant GPU memory; every thread does not need a copy...
  // The implementation lives in splines/gpuSplineUtils.cu
  // The GPU splines don't actually need declaring but is good for demonstration, kind of
  // fixed by passing const reference
  gpu_spline_handler->CopyToGPU_SplineMonolith(
          cpu_spline_handler,

          // TFI related now
          cpu_coeff_TF1_many,
          cpu_paramNo_TF1_arr,
          #ifndef Weight_On_SplineBySpline_Basis
          NEvents,
          cpu_nParamPerEvent,
          cpu_nParamPerEvent_tf1,
          #endif
          nParams,
          NSplines_valid,
          _max_knots,
          nKnots,
          NTF1_valid);

  // Delete all the coefficient arrays from the CPU once they are on the GPU
  cpu_spline_handler->coeff_x.clear();
  cpu_spline_handler->coeff_x.shrink_to_fit();
  cpu_spline_handler->coeff_many.clear();
  cpu_spline_handler->coeff_many.shrink_to_fit();
  cpu_spline_handler->paramNo_arr.clear();
  cpu_spline_handler->paramNo_arr.shrink_to_fit();
  cpu_spline_handler->nKnots_arr.clear();
  cpu_spline_handler->nKnots_arr.shrink_to_fit();
  cpu_coeff_TF1_many.clear();
  cpu_coeff_TF1_many.shrink_to_fit();
  cpu_paramNo_TF1_arr.clear();
  cpu_paramNo_TF1_arr.shrink_to_fit();
  #ifndef Weight_On_SplineBySpline_Basis
  cpu_nParamPerEvent.clear();
  cpu_nParamPerEvent.shrink_to_fit();
  cpu_nParamPerEvent_tf1.clear();
  cpu_nParamPerEvent_tf1.shrink_to_fit();
  #endif
  delete cpu_spline_handler;
  cpu_spline_handler = nullptr;
  MACH3LOG_INFO("Good GPU loading");
  #endif
  return;
}

// Need to specify template functions in header
// *****************************************
// Scan the master spline to get the maximum number of knots in any of the TSpline3*
void SMonolith::ScanMasterSpline(std::vector<std::vector<TResponseFunction_red*> > & MasterSpline,
                                 unsigned int &nEvents,
                                 int &MaxPoints,
                                 short int &numParams,
                                 int &nSplines,
                                 unsigned int &NSplinesValid,
                                 unsigned int &numKnots,
                                 unsigned int &nTF1Valid,
                                 unsigned int &nTF1_coeff,
                                 const std::vector<RespFuncType> &SplineType) {
// *****************************************
  // Need to extract: the total number of events
  //                  number of parameters
  //                  maximum number of knots
  MaxPoints = 0;
  nEvents   = 0;
  numParams   = 0;
  nSplines = 0;
  numKnots = 0;
  NSplinesValid = 0;
  nTF1Valid = 0;
  nTF1_coeff = 0;

  // Check the number of events
  nEvents = MasterSpline.size();

  // Maximum number of splines one event can have (scan through and find this number)
  int nMaxSplines_PerEvent = 0;

  //KS: We later check that each event has the same number of splines so this is fine
  numParams = MasterSpline[0].size();
  // Initialise
  SplineInfoArray = new FastSplineInfo[numParams];

  // Loop over each parameter
  for(unsigned int EventCounter = 0; EventCounter < MasterSpline.size(); ++EventCounter) {
    // Check that each event has each spline saved
    if (numParams > 0) {
      int TempSize = MasterSpline[EventCounter].size();
      if (TempSize != numParams) {
        MACH3LOG_ERROR("Found {} parameters for event {}", TempSize, EventCounter);
        MACH3LOG_ERROR("but was expecting {} since that's what I found for the previous event", numParams);
        MACH3LOG_ERROR("Somehow this event has a different number of spline parameters... Please study further!");
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }
    }
    numParams = MasterSpline[EventCounter].size();

    int nSplines_SingleEvent = 0;
    // Loop over each pointer
    for(unsigned int ParamNumber = 0; ParamNumber < MasterSpline[EventCounter].size(); ++ParamNumber) {
      // If NULL we don't have this spline for the event, so move to next spline
      if (MasterSpline[EventCounter][ParamNumber] == NULL) continue;

      if(SplineType[ParamNumber] == kTSpline3_red)
      {
        TResponseFunction_red* TespFunc = MasterSpline[EventCounter][ParamNumber];
        TSpline3_red* CurrSpline = dynamic_cast<TSpline3_red*>(TespFunc);
        int nPoints = CurrSpline->GetNp();
        if (nPoints > MaxPoints) {
          MaxPoints = nPoints;
        }
        numKnots += nPoints;
        nSplines_SingleEvent++;

        // Fill the SplineInfoArray entries with information on each splinified parameter
        if (SplineInfoArray[ParamNumber].xPts == NULL)
        {
          // Fill the number of points
          SplineInfoArray[ParamNumber].nPts = CurrSpline->GetNp();

          // Fill the x points
          SplineInfoArray[ParamNumber].xPts = new _float_[SplineInfoArray[ParamNumber].nPts];
          for (_int_ k = 0; k < SplineInfoArray[ParamNumber].nPts; ++k)
          {
            _float_ xtemp = -999.99;
            _float_ ytemp = -999.99;
            CurrSpline->GetKnot(k, xtemp, ytemp);
            SplineInfoArray[ParamNumber].xPts[k] = xtemp;
          }
        }
        NSplinesValid++;
      }
      else if (SplineType[ParamNumber] == kTF1_red)
      {
        TResponseFunction_red* TespFunc = MasterSpline[EventCounter][ParamNumber];
        TF1_red* CurrSpline = dynamic_cast<TF1_red*>(TespFunc);
        int nPoints = CurrSpline->GetSize();
        nTF1_coeff += nPoints;
        nTF1Valid++;
      }
    }
    if (nSplines_SingleEvent > nMaxSplines_PerEvent) nMaxSplines_PerEvent = nSplines_SingleEvent;
  }
  nSplines = nMaxSplines_PerEvent;

  int Counter = 0;
  //KS: Sanity check that everything was set correctly
  for (_int_ i = 0; i < numParams; ++i)
  {
    // KS: We don't find segment for TF1, so ignore this
    if (SplineType[i] == kTF1_red) continue;

    const _int_ nPoints = SplineInfoArray[i].nPts;
    const _float_* xArray = SplineInfoArray[i].xPts;
    if (nPoints == -999 || xArray == NULL) {
      Counter++;
      if(Counter < 5) {
        MACH3LOG_WARN("SplineInfoArray[{}] isn't set yet", i);
      }
      continue;
      //throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  }
  MACH3LOG_WARN("In total SplineInfoArray for {} hasn't been initialised", Counter);
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
// Load SplineMonolith from ROOT file
void SMonolith::LoadSplineFile(std::string FileName) {
// *****************************************
  #ifdef Weight_On_SplineBySpline_Basis
  MACH3LOG_ERROR("Trying to load Monolith from file using weight by weight base, this is not supported right now, sorry");
  throw MaCh3Exception(__FILE__ , __LINE__ );
  #endif

  if (std::getenv("MACH3") != NULL) {
      FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
   }

  TFile *SplineFile = new TFile(FileName.c_str(), "OPEN");
  TTree *Settings = (TTree*)SplineFile->Get("Settings");
  TTree *Monolith = (TTree*)SplineFile->Get("Monolith");
  TTree *Monolith_TF1 = (TTree*)SplineFile->Get("Monolith_TF1");
  TTree *ParamInfo = (TTree*)SplineFile->Get("ParamInfo");
  TTree *XKnots = (TTree*)SplineFile->Get("XKnots");
  TTree *EventInfo = (TTree*)SplineFile->Get("EventInfo");
  TTree *FastSplineInfoTree = (TTree*)SplineFile->Get("FastSplineInfoTree");

  unsigned int NEvents_temp;
  short int nParams_temp;
  int _max_knots_temp;
  unsigned int nKnots_temp;
  unsigned int NSplines_valid_temp;
  unsigned int nTF1Valid_temp;
  unsigned int nTF1coeff_temp;

  Settings->SetBranchAddress("NEvents", &NEvents_temp);
  Settings->SetBranchAddress("nParams", &nParams_temp);
  Settings->SetBranchAddress("_max_knots", &_max_knots_temp);
  Settings->SetBranchAddress("nKnots", &nKnots_temp);
  Settings->SetBranchAddress("NSplines_valid", &NSplines_valid_temp);
  Settings->SetBranchAddress("NTF1_valid", &nTF1Valid_temp);
  Settings->SetBranchAddress("nTF1coeff", &nTF1coeff_temp);

  Settings->GetEntry(0);

  NEvents = NEvents_temp;
  nParams = nParams_temp;
  _max_knots = _max_knots_temp;
  nKnots = nKnots_temp;
  NSplines_valid = NSplines_valid_temp;
  NTF1_valid = nTF1Valid_temp;
  unsigned int event_size_max = _max_knots * nParams;
  nTF1coeff = nTF1coeff_temp;

  //KS: Since we are going to copy it each step use fancy CUDA memory allocation
#ifdef CUDA
  gpu_spline_handler->InitGPU_Segments(&segments);
  gpu_spline_handler->InitGPU_Vals(&vals);
#elif USE_FPGA
  queue = sycl::queue(sycl::default_selector{})
  segments = sycl::malloc_shared<short int>(nParams, queue);
  vals = new float[nParams]();
#else
  segments = new short int[nParams]();
  vals = new float[nParams]();
#endif



  cpu_nParamPerEvent.resize(2*NEvents);
  cpu_nParamPerEvent_tf1.resize(2*NEvents);
  #ifdef USE_FPGA
    cpu_spline_handler = new SplineMonoUSM(queue, event_size_max, nKnots*_nCoeff_, NSplines_valid, NSplines_valid);
  #else
    cpu_spline_handler->paramNo_arr.resize(NSplines_valid);
    //KS: And array which tells where each spline stars in a big monolith array, sort of knot map
    cpu_spline_handler->nKnots_arr.resize(NSplines_valid);

    cpu_spline_handler->coeff_many.resize(nKnots*_nCoeff_); // *4 because we store y,b,c,d parameters in this array
    cpu_spline_handler->coeff_x.resize(event_size_max);
  #endif


  cpu_coeff_TF1_many.resize(nTF1coeff);

  //KS: This is tricky as this variable use both by CPU and GPU, however if use CUDA we use cudaMallocHost
#ifndef CUDA
  cpu_total_weights = new float[NEvents]();
  #ifdef USE_FPGA
    cpu_weights_var = sycl::malloc_shared<float>(NSplines_valid, queue);
  #else
    cpu_weights_var = new float[NSplines_valid]();
  #endif

  
  cpu_weights_tf1_var = new float[NTF1_valid]();
#endif

  float coeff = 0.;
  Monolith->SetBranchAddress("cpu_coeff_many", &coeff);
  for(unsigned int i = 0; i < nKnots*_nCoeff_; i++)
  {
    Monolith->GetEntry(i);
    cpu_spline_handler->coeff_many[i] = coeff;
  }

  float coeff_tf1 = 0.;
  Monolith_TF1->SetBranchAddress("cpu_coeff_TF1_many", &coeff_tf1);
  for(unsigned int i = 0; i < nTF1coeff; i++)
  {
    Monolith_TF1->GetEntry(i);
    cpu_coeff_TF1_many[i] = coeff_tf1;
  }

  short int paramNo_arr = 0;
  unsigned int nKnots_arr = 0;
  ParamInfo->SetBranchAddress("cpu_paramNo_arr", &paramNo_arr);
  ParamInfo->SetBranchAddress("cpu_nKnots_arr", &nKnots_arr);
  for(unsigned int i = 0; i < NSplines_valid; i++)
  {
    ParamInfo->GetEntry(i);
    cpu_spline_handler->paramNo_arr[i] = paramNo_arr;
    cpu_spline_handler->nKnots_arr[i] = nKnots_arr;
  }

  float coeff_x = 0.;
  XKnots->SetBranchAddress("cpu_coeff_x", &coeff_x);
  for(unsigned int i = 0; i < event_size_max; i++)
  {
    XKnots->GetEntry(i);
    cpu_spline_handler->coeff_x[i] = coeff_x;
  }

  unsigned int nParamPerEvent = 0;
  unsigned int nParamPerEvent_tf1 = 0;

  EventInfo->SetBranchAddress("cpu_nParamPerEvent", &nParamPerEvent);
  EventInfo->SetBranchAddress("cpu_nParamPerEvent_tf1", &nParamPerEvent_tf1);
  for(unsigned int i = 0; i < 2*NEvents; i++)
  {
    EventInfo->GetEntry(i);
    cpu_nParamPerEvent[i] = nParamPerEvent;
    cpu_nParamPerEvent_tf1[i] = nParamPerEvent_tf1;
  }

  _int_ nPoints = 0;
  float xtemp[20];
  FastSplineInfoTree->SetBranchAddress("nPts", &nPoints);
  FastSplineInfoTree->SetBranchAddress("xPts", &xtemp);

  SplineInfoArray = new FastSplineInfo[nParams];
  for (_int_ i = 0; i < nParams; ++i) {
    FastSplineInfoTree->GetEntry(i);

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
  PrintInitialsiation();

  MoveToGPU();
}

// *****************************************
// Save SplineMonolith into ROOT file
void SMonolith::PrepareSplineFile() {
// *****************************************
  std::string FileName = "SplineFile.root";
  if (std::getenv("MACH3") != NULL) {
      FileName.insert(0, std::string(std::getenv("MACH3"))+"/");
   }

  TFile *SplineFile = new TFile(FileName.c_str(), "recreate");
  TTree *Settings = new TTree("Settings", "Settings");
  TTree *Monolith = new TTree("Monolith", "Monolith");
  TTree *Monolith_TF1 = new TTree("Monolith_TF1", "Monolith_TF1");
  TTree *ParamInfo = new TTree("ParamInfo", "ParamInfo");
  TTree *XKnots = new TTree("XKnots", "XKnots");
  TTree *EventInfo = new TTree("EventInfo", "EventInfo");
  TTree *FastSplineInfoTree = new TTree("FastSplineInfoTree", "FastSplineInfoTree");

  unsigned int NEvents_temp = NEvents;
  short int nParams_temp = nParams;
  int _max_knots_temp = _max_knots;
  unsigned int nKnots_temp = nKnots;
  unsigned int NSplines_valid_temp = NSplines_valid;
  unsigned int nTF1Valid_temp = NTF1_valid;
  unsigned int nTF1coeff_temp = nTF1coeff;

  Settings->Branch("NEvents", &NEvents_temp, "NEvents/i");
  Settings->Branch("nParams", &nParams_temp, "nParams/S");
  Settings->Branch("_max_knots", &_max_knots_temp, "_max_knots/I");
  Settings->Branch("nKnots", &nKnots_temp, "nKnots/i");
  Settings->Branch("NSplines_valid", &NSplines_valid_temp, "NSplines_valid/i");
  Settings->Branch("NTF1_valid", &nTF1Valid_temp, "NTF1_valid/i");
  Settings->Branch("nTF1coeff", &nTF1coeff_temp, "nTF1coeff/i");

  Settings->Fill();

  SplineFile->cd();
  Settings->Write();

  float coeff = 0.;
  Monolith->Branch("cpu_coeff_many", &coeff, "cpu_coeff_many/F");
  for(unsigned int i = 0; i < nKnots*_nCoeff_; i++)
  {
    coeff = cpu_spline_handler->coeff_many[i];
    Monolith->Fill();
  }
  SplineFile->cd();
  Monolith->Write();


  float coeff_tf1 = 0.;
  Monolith_TF1->Branch("cpu_coeff_TF1_many", &coeff_tf1, "cpu_coeff_TF1_many/F");
  for(unsigned int i = 0; i < nTF1coeff; i++)
  {
    coeff_tf1 = cpu_coeff_TF1_many[i];
    Monolith_TF1->Fill();
  }
  SplineFile->cd();
  Monolith_TF1->Write();

  short int paramNo_arr = 0;
  unsigned int nKnots_arr = 0;
  ParamInfo->Branch("cpu_paramNo_arr", &paramNo_arr, "cpu_paramNo_arr/S");
  ParamInfo->Branch("cpu_nKnots_arr", &nKnots_arr, "cpu_nKnots_arr/i");
  for(unsigned int i = 0; i < NSplines_valid; i++)
  {
    paramNo_arr = cpu_spline_handler->paramNo_arr[i];
    nKnots_arr = cpu_spline_handler->nKnots_arr[i];

    ParamInfo->Fill();
  }
  SplineFile->cd();
  ParamInfo->Write();

  unsigned int event_size_max = _max_knots * nParams;

  float coeff_x = 0.;
  XKnots->Branch("cpu_coeff_x", &coeff_x, "cpu_coeff_x/F");
  for(unsigned int i = 0; i < event_size_max; i++)
  {
    coeff_x = cpu_spline_handler->coeff_x[i];
    XKnots->Fill();
  }
  SplineFile->cd();
  XKnots->Write();

  unsigned int nParamPerEvent = 0;
  unsigned int nParamPerEvent_tf1 = 0;

  EventInfo->Branch("cpu_nParamPerEvent", &nParamPerEvent, "cpu_nParamPerEvent/i");
  EventInfo->Branch("cpu_nParamPerEvent_tf1", &nParamPerEvent_tf1, "cpu_nParamPerEvent_tf1/i");

  for(unsigned int i = 0; i < 2*NEvents; i++)
  {
    nParamPerEvent = cpu_nParamPerEvent[i];
    nParamPerEvent_tf1 = cpu_nParamPerEvent_tf1[i];
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
  delete Monolith_TF1;
  delete ParamInfo;
  delete XKnots;
  delete EventInfo;
  delete FastSplineInfoTree;
  SplineFile->Close();
  delete SplineFile;
}

// *****************************************
// Destructor
// Cleans up the allocated GPU memory
SMonolith::~SMonolith() {
  // *****************************************

  #ifdef CUDA
  gpu_spline_handler->CleanupGPU_SplineMonolith(
        #ifndef Weight_On_SplineBySpline_Basis
        cpu_total_weights
        #endif
        );

  //KS: Since we declared them using CUDA alloc we have to free memory using also cuda functions
  gpu_spline_handler->CleanupGPU_Segments(segments, vals);

  delete gpu_spline_handler;
  #else
    #ifdef USE_FPGA
      if(segments != nullptr) sycl::free(segments, queue);
    #else
      if(segments != nullptr) delete[] segments;
    #endif
  if(vals != nullptr) delete[] vals;
  if(cpu_total_weights != nullptr) delete[] cpu_total_weights;
  #endif

  if(SplineInfoArray != nullptr) delete[] SplineInfoArray;
  if(cpu_weights != nullptr) delete[] cpu_weights;
  #ifdef USE_FPGA
    if(cpu_weights_var != nullptr) sycl::free(cpu_weights_var, queue);
  #else
    if(cpu_weights_var != nullptr) delete[] cpu_weights_var;
  #endif
  if(cpu_weights_tf1_var != nullptr) delete[] cpu_weights_tf1_var;
  if(index_cpu != nullptr) delete[] index_cpu;
  if(index_TF1_cpu != nullptr) delete[] index_TF1_cpu;

  //KS: Those might be deleted or not depending on GPU/CPU TSpline3/TF1 DEBUG or not hence we check if not NULL
  #ifndef USE_FPGA
    if(cpu_spline_handler != nullptr)
    {
      cpu_spline_handler->coeff_x.clear();
      cpu_spline_handler->coeff_x.shrink_to_fit();
      cpu_spline_handler->coeff_many.clear();
      cpu_spline_handler->coeff_many.shrink_to_fit();
      cpu_spline_handler->paramNo_arr.clear();
      cpu_spline_handler->paramNo_arr.shrink_to_fit();
      cpu_spline_handler->nKnots_arr.clear();
      cpu_spline_handler->nKnots_arr.shrink_to_fit();
    }
  #endif
  cpu_coeff_TF1_many.clear();
  cpu_coeff_TF1_many.shrink_to_fit();
  cpu_paramNo_TF1_arr.clear();
  cpu_paramNo_TF1_arr.shrink_to_fit();
  #ifndef Weight_On_SplineBySpline_Basis
  cpu_nParamPerEvent.clear();
  cpu_nParamPerEvent.shrink_to_fit();
  cpu_nParamPerEvent_tf1.clear();
  cpu_nParamPerEvent_tf1.shrink_to_fit();
  #endif
  cpu_nPoints_arr.clear();
  cpu_nPoints_arr.shrink_to_fit();

  if(cpu_spline_handler != nullptr) delete cpu_spline_handler;
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
  nPoints = Np;
  if (Np > _max_knots) {
    MACH3LOG_ERROR("Error, number of points is greater than saved {}", _max_knots);
    MACH3LOG_ERROR("This _WILL_ cause problems with GPU splines and _SHOULD_ be fixed!");
    MACH3LOG_ERROR("nPoints = {}, _max_knots = {}", nPoints, _max_knots);
    throw MaCh3Exception(__FILE__ , __LINE__ );
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
    if((xArray[i] == -999) || (manyArray[i*_nCoeff_] == -999) || (manyArray[i*4+1] == -999) || (manyArray[i*_nCoeff_+2] == -999) || (manyArray[i*_nCoeff_+3] == -999)){
      MACH3LOG_ERROR("*********** Bad params in getSplineCoeff_SepMany() ************");
      MACH3LOG_ERROR("pre cast to float (x, y, b, c, d) = {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}", x, y, b, c, d);
      MACH3LOG_ERROR("pre cast to float (x, y, b, c, d) = {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}", xArray[i], manyArray[i*4], manyArray[i*4+1], manyArray[i*4+2], manyArray[i*_nCoeff_+3]);
      MACH3LOG_ERROR("This will cause problems when preparing for GPU");
      MACH3LOG_ERROR("***************************************************************");
    }
  }
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
  gpu_spline_handler->RunGPU_SplineMonolith(
    #ifdef Weight_On_SplineBySpline_Basis
          cpu_weights_var,
          cpu_weights_tf1_var,
    #else
          cpu_total_weights,
    #endif
          vals,
          segments,
          NSplines_valid,
          NTF1_valid);

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

    // Get the variation for this reconfigure for the ith parameter
    const _float_ xvar = *SplineInfoArray[i].splineParsPointer;
    vals[i] = xvar;

    // EM: if we have a parameter that has no response for any event (i.e. all splines have just one knot), then skip it and avoid a seg fault here
    //     In principle, such parameters shouldn't really be included in the first place, but with new det syst splines this
    //     could happen if say you were just running on one FHC run, then all RHC parameters would be flat and the code below would break.
    if(xArray == NULL) continue;

    // The segment we're interested in (klow in ROOT code)
    _int_ segment = 0;
    _int_ kHigh = nPoints-1;
    //KS: We expect new segment is very close to previous
    const _int_ PreviousSegment = SplineInfoArray[i].CurrSegment;

    // If the variation is below the lowest saved spline point
    if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      //CW: Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      //KS: It is quite probable the new segment is same as in previous step so try to avoid binary search, first we have to check if it is in bounds to avoid seg fault
    } else if( xArray[PreviousSegment+1] > xvar && xvar >= xArray[PreviousSegment] ) {
      segment = PreviousSegment;
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
      MACH3LOG_ERROR("Found a segment which is _ABOVE_ the variation!");
      MACH3LOG_ERROR("IT SHOULD ALWAYS BE BELOW! (except when segment 0)");
      MACH3LOG_ERROR("Spline: {}", i);
      MACH3LOG_ERROR("Found segment = {}", segment);
      MACH3LOG_ERROR("Doing variation = {}", xvar);
      MACH3LOG_ERROR("x in spline = {}", SplineInfoArray[i].xPts[segment]);
      for (_int_ j = 0; j < SplineInfoArray[j].nPts; ++j) {
        MACH3LOG_ERROR("    {} = {}", j, SplineInfoArray[i].xPts[j]);
      }
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
#endif
  } //end loop over params
}

//*********************************************************
void SMonolith::CalcSplineWeights() {
//*********************************************************
  #ifdef MULTITHREAD
  //KS: Open parallel region
  #pragma omp parallel
  {
  #endif
    //KS: First we calculate
    #ifdef MULTITHREAD
    #pragma omp for simd nowait
    #endif
    for (unsigned int splineNum = 0; splineNum < NSplines_valid; ++splineNum)
    {
      //CW: Which Parameter we are accessing
      const short int Param = cpu_spline_handler->paramNo_arr[splineNum];

      //CW: Avoids doing costly binary search on GPU
      const short int segment = segments[Param];

      //KS: Segment for coeff_x is simply parameter*max knots + segment as each parameters has the same spacing
      const short int segment_X = Param*_max_knots+segment;

      //KS: Find knot position in out monolithical structure
      const unsigned int CurrentKnotPos = cpu_spline_handler->nKnots_arr[splineNum]*_nCoeff_+segment*_nCoeff_;

      // We've read the segment straight from CPU and is saved in segment_gpu
      // polynomial parameters from the monolithic splineMonolith
      const float fY = cpu_spline_handler->coeff_many[CurrentKnotPos];
      const float fB = cpu_spline_handler->coeff_many[CurrentKnotPos+1];
      const float fC = cpu_spline_handler->coeff_many[CurrentKnotPos+2];
      const float fD = cpu_spline_handler->coeff_many[CurrentKnotPos+3];
      // The is the variation itself (needed to evaluate variation - stored spline point = dx)
      const float dx = vals[Param] - cpu_spline_handler->coeff_x[segment_X];

      //CW: Wooow, let's use some fancy intrinsic and pull down the processing time by <1% from normal multiplication! HURRAY
      cpu_weights_var[splineNum] = fmaf(dx, fmaf(dx, fmaf(dx, fD, fC), fB), fY);
      // Or for the more "easy to read" version:
      //cpu_weights_var[splineNum] = (fY+dx*(fB+dx*(fC+dx*fD)));
    }

    #ifdef MULTITHREAD
    #pragma omp for simd
    #endif
    for (unsigned int tf1Num = 0; tf1Num < NTF1_valid; ++tf1Num)
    {
      // The is the variation itself (needed to evaluate variation - stored spline point = dx)
      const float x = vals[cpu_paramNo_TF1_arr[tf1Num]];

      // Read the coefficients
      const float a = cpu_coeff_TF1_many[tf1Num*_nTF1Coeff_];
      const float b = cpu_coeff_TF1_many[tf1Num*_nTF1Coeff_+1];

      cpu_weights_tf1_var[tf1Num] = fmaf(a, x, b);
      // cpu_weights_tf1_var[tf1Num] = a*x + b;

      //cpu_weights_tf1_var[splineNum] = 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x + e*x*x*x*x*x;
    }
  #ifdef MULTITHREAD
  //KS: End parallel region
  }
  #endif
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
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (unsigned int id = 0; id < numParams; ++id) {
      totalWeight *= cpu_weights_var[startIndex + id];
    }
    //Now TF1
    // Extract the parameters for the current event
    const unsigned int startIndex_tf1 = cpu_nParamPerEvent_tf1[2 * EventNum + 1];
    const unsigned int numParams_tf1 = cpu_nParamPerEvent_tf1[2 * EventNum];

    // Compute total weight for the current event
    #ifdef MULTITHREAD
    #pragma omp simd
    #endif
    for (unsigned int id = 0; id < numParams_tf1; ++id) {
      totalWeight *= cpu_weights_tf1_var[startIndex_tf1 + id];
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
  for (unsigned int i = 0; i < NSplines_total_large; ++i) {
    if (index_cpu[i] >= 0) {
      cpu_weights[i] = cpu_weights_var[index_cpu[i]];
    } else if (index_TF1_cpu[i] >= 0) {
      cpu_weights[i] = cpu_weights_tf1_var[index_TF1_cpu[i]];
    }  else {
      cpu_weights[i] = 1.;
    }
  }
#endif
  return;
}

//*********************************************************
//KS: Print info about how much knots etc has been initialised
void SMonolith::PrintInitialsiation() {
//*********************************************************
  unsigned int event_size_max = _max_knots * nParams;

  MACH3LOG_INFO("--- INITIALISED Spline Monolith ---");
  MACH3LOG_INFO("{} events with {} splines", NEvents, NSplines_valid);
  MACH3LOG_INFO("On average {:.2f} splines per event ({}/{})", float(NSplines_valid)/float(NEvents), NSplines_valid, NEvents);
  MACH3LOG_INFO("Size of x array = {:.4f} MB", double(sizeof(float)*event_size_max)/1.E6);
  MACH3LOG_INFO("Size of coefficient (y,b,c,d) array = {:.2f} MB", double(sizeof(float)*nKnots*_nCoeff_)/1.E6);
  MACH3LOG_INFO("Size of parameter # array = {:.2f} MB", double(sizeof(short int)*NSplines_valid)/1.E6);

  MACH3LOG_INFO("On average {:.2f} TF1 per event ({}/{})", float(NTF1_valid)/float(NEvents), NTF1_valid, NEvents);
  MACH3LOG_INFO("Size of TF1 coefficient (a,b,c,d,e) array = {:.2f} MB", double(sizeof(float)*NTF1_valid*_nTF1Coeff_)/1.E6);

  return;
}

//*********************************************************
//KS: After calculations are done on GPU we copy memory to CPU. This operation is asynchronous meaning while memory is being copied some operations are being carried. Memory must be copied before actual reweight. This function make sure all has been copied.
void SMonolith::SynchroniseMemTransfer() {
//*********************************************************
  #ifdef CUDA
  SynchroniseSplines();
  #endif
  return;
}
