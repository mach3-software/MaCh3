# Making a Sample Handler {#making-a-sample-handler}

The samplePDF classes contain information about a particular sample in your analysis. Practically this involves storing the number of events predicted in each analysis bin for each sample as well as modifying this for a given set of systematic and oscillation parameter values (i.e. "reweighting" the nominal Monte-Carlo prediction). Therefore, the samplePDF class also "talks" to a lot of other classes in MaCh3 especially the covariance and spline classes. In turn the FitterBase class then gets a poisson likelihood value from comparing the MC prediction and data in the samplePDF class.

Each experiment will need to create a samplePDF experiment class which inherits from a the samplePDFFDBase class in the core repository. This is intended as a guide to help you produce your own samplePDF class for your own experiment. This is intended to be as lightweight as possible and should mainly relate to the MC format that a particular experiment is using. You need to read this information into the data structures which the samplePDFFDBase class is expecting so that it can handle the "heavy lifting" such as the reweighting of the MC prediction.

**N.B.** - Currently (April 2024), experiment samplePDF classes should inherit from the samplePDFFDBase class in core (where FD stands for Far Detector). There is a samplePDFND (Near Detector) base class but this has not been tested for multiple experiments yet. In future releases of the core code (v2) we plan to have a single samplePDF class for both ND and FD to simplify this inheritance structure.

# Constructing a samplePDF experiment class
Discussion of FD data structures; FDMC and ExpFDMC.
The main purpose of any samplePDF class is to contain all the information you need to make MC predictions of your selections for some given set of systematic parameter values. To do this, the samplePDF object needs to read in information from different input files and keep track of them in some data structures. Most of your samplePDF class will be reading in this information and putting it in the correct place.

Below is a skeleton of what the constructor of your samplePDF class will look like.
```
samplePDFSKBase::samplePDFSKBase(double pot, std::string samplecfgfile, covarianceXsec* xsec_cov)
  : samplePDFFDBase(pot, mc_version, xsec_cov)
{
  manager* SampleManager = new manager(samplecfgfile);

  ////////////////
  /// Read information about sample from sample.yaml config file
  ////////////////
  //Need to set this before call to SetXecCov as SampleDetID is used there
  ReadInformationFromSampleCfg();

  //Read input MC and spline file names for each oscillation channel
  //also need to store oscillation channel information i.e. numu->nue? or numu->numu etc. 
  ReadMCandSplineFiles(SampleManager);

  //Read any selection cuts
  GetSelectionCuts(SampleManager);

  //Set the xsec covariance
  SetXsecCov(xsec_cov);

  //Now create a spline object 
  splineFile = new splineExpBase(XsecCov);
  
  ////////////////
  //Create some structs to store the information in both for the expimernt
  //and the information for the core samplePDF FD class
  //then fill these. There is one of these for each oscillaiton channel.
  ////////////////

  for(unsigned iOscChannel=0 ; iOscChannel < nOscChannels ; iOscChannel++){
    //Fill from MC file
    SetupExpMC(&ExpStructs[iSample]);

        //Now fill the core bits in the samplePDFFDBase class
        //This is often just pointing to information in the experiment class
    SetupFDMC(&BaseClassStructs[iSample]);
  }

  //Setup the splines
  SetupSplines();

  ////////////////////////////////////
  // Setup your norm parameters as well here
  ////////////////////////////////////
  SetupNormParameters();

  //Also your func parameters if you're lucky
  //enough to have them
  SetupFuncParameters();

  //Setup pointers to all the weights you will apply
  //e.g. from xsec systsematics, detector, pot weight etc.
  SetupWeightPointers();

  //The binning here is arbitrary, now we get info from cfg so the
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D   = new TH1D("h"+histname1D+samplename,histtitle1D, 1, 0, 1);
  dathist   = new TH1D("d"+histname1D+samplename,histtitle1D, 1, 0, 1);
  _hPDF2D   = new TH2D("h"+histname2D+samplename,histtitle2D, 1, 0, 1, 1, 0, 1);
  dathist2d = new TH2D("d"+histname2D+samplename,histtitle2D, 1, 0, 1, 1, 0, 1);

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double XBinEdges[SampleXBins.size()];
  double YBinEdges[SampleYBins.size()];
  for(unsigned XBin_i = 0 ; XBin_i < SampleXBins.size() ; XBin_i++){
        XBinEdges[XBin_i] = SampleXBins[XBin_i];
  }
  //And now the YBin Edges
  for(unsigned YBin_i = 0 ; YBin_i < SampleYBins.size() ; YBin_i++){
        YBinEdges[YBin_i] = SampleYBins[YBin_i];
  }

  //Calls to functions in the core code
  set1DBinning(SampleXBins);
  set2DBinning(SampleXBins, SampleYBins);

  //This then points each event to all the splines it is associated with
  fillSplineBins();
  splineFile->cleanUpMemory();
}
```
The code above will not compile and is only illustrative only the main steps. As you can see the majority of the setup is above grabbing things from a **sample yaml file** and **passing information into an experiment struct and a core struct**. The discussion below is intended to expand on these. A discussion of the spline setup is given in another section.

# Sample Config files
As described above, a lot of information is specified in a "sample config". This is a yaml file which contains information about the MC events you would like to read in. Below is an example:
```
SampleName: "YourFavouriteSample"
Selection:
  #See below for a discussion on this
  SampleDetID: 8
#Apply cuts on particular variables if you would like to
#You will need to make sure this string maps to a variable
SelectionCuts:
  - KinematicStr: "RecoNeutrinoEnergy"
    Bounds: [ 0., 30000. ]
Binning:
  #See below for a discussion on this
  BinningOpt: 0
  XVarBins: [0.0, 100.0, 200.0, 300.0, 450.0, 600.0, 750.0, 1000.0]
  #If you want another axis as well specify this here
  YVarBins: [0.0, 180.0]
#Any bools you might want to keep track of
SampleBools:
  isrhc: no
#Paths to the MC inputs and the spline files
InputFiles:
  mtupleprefix: "inputs/PrefixFileName"
  mtuplesuffix: "_blah.root"
  splineprefix: "inputs/PrefixSplineFileName"
  splinesuffix: "_blah.root"
#This corresponds to the number of oscillation channels present
#you probably have 6 of these (if you don't have tau channels) otherwise you have the full 12
NSubSamples: 1
SubSamples:
  - name: "numu-x-numu"
    mtuplefile: "numu_x_numu"
    splinefile: "numu_x_numu"
    samplevecno: 0
    #PDG pre oscillation
    nutype: 14
    #PDG post oscillation
    oscnutype: 14
    signal: false
```
This information simply needs to be grabbed from the yaml file and stored in your samplePDFExperiment base.
**A note on SampleDetID**: SampleDetID is a little complicated. This needs to be a unique number that is a power of 2. This is used to check whether a systematic applies to a sample or not. It has to be a power of two as this is checked using bit-shifting. This will be replaced in the future.
**A note on BinningOpt**: BinningOpt is currently used to describe the number of variables you want to bin in and what they are. This is quite hard-coded and will be replaced in the future to be flexible.
**A note on KinematicStr**: we will discussion the Kinematic Variables later in the virtual functions. But you will need to supply a string which can be mapped back to get the value of a specific variable for an event so you can cut on it.

# The Far Detector Struct
One of the most important things is passing information from your MC file to the core code. All the information required for the core code to reweight MC.

The far detector struct can be seen [here](https://github.com/mach3-software/MaCh3/blob/develop/samplePDF/FDMCStruct.h). All of these variables must be filled from your MC. It is highly recommended that you also create a similar struct to hold the information of your experiment. A lot of the members of the Far Detector struct are arrays of pointers, this is so they can point to your experiment struct without storing lots of floats or double twice. Also this has the advantage that if you update a value in your samplePDFExperiment class then the far detector struct does not need updating as your are just pointing to the same part of memory. This structure means that your SetupExpMC and SetupFDMC() functions would look something like this:

```
void SetupExpMC(){
  //Grab variables from your MC file (let's say it's a ROOT TTree)
  TTree *MCTree = InputFile->Get("MyExperimentsMC");

  double TrueNeutrinoEnergy;
  double RecoNeutrinoEnergy;

  MCTree->SetBranchAddress("TrueNeutrinoEnergy", &TrueNeutrinoEnergy);
  MCTree->SetBranchAddress("RecoNeutrinoEnergy", &RecoNeutrinoEnergy);

  //Now create the "space" to store these variables in your experiment struct (ExpMC)
  struct ExpMC ExpObj = ExpMC();
  ExpObj->TrueEnu = new double[MCTree->GetEntries()];
  ExpObj->RecoEnu = new double[MCTree->GetEntries()];

  //Loop over all your MC events and fill your Exp MC struct with the info you need
  for(unsigned int iEntry = 0 ; iEntry < MCTree->GetEntries() ; ++iEntry){
    ExpObj->TrueEnu[iEntry] = TrueNeutrinoEnergy;
    ExpObj->RecoEnu[iEntry] = RecoNeutrinoEnergy;
  }
}

void SetupFDMC(){
  
  struct fdmc_base FDObj = fdmc_base();
  FDObj->x_var = new double*[nEvents];
  FDObj->rw_etrue = new double*[nEvents];

  //Loop over and setup your fdmc_base struct to point to or to copy values from you experiment struct
  for(unsigned int iEvent = 0 ; iEvent < nEvents ; ++iEvent){
    //Let's say you want your sample to be binned in reconstructed neutrino energy
    //To do this we set the x_var to point to the Reconstructed neutrino energy value in the experiment struct
    FDObj->x_var[iEvent] = &(ExpObj->RecoEnu[iEvent]);
  }
}
```
Some of the variables are setup automatically for you by the core code. For example NomXBin, NomYBin, XBin, YBin, rw_lower_xbinedge, rw_upper_xbinedge, rw_lower_lower_xbinedge and rw_upper_upper_xbinedge are set in the [FindNominalBinAndEdges1D](https://github.com/mach3-software/MaCh3/blob/develop/samplePDF/samplePDFFDBase.cpp#L1036C23-L1036C47) function (or it's 2D counterpart). So in SetupFDMC you just need to create the array and set it to some dummy value to begin with.

# Setting Up Systematics
Thankfully, a lot of this is now handled in the core code and if you have passed all the information to the fdmc_base struct you should not have to worry too much about this. The information for systematics is stored in your covarianceXsec class which parses information from a YAML file.
[**SetupNormParameters()**](https://github.com/mach3-software/MaCh3/blob/develop/samplePDF/samplePDFFDBase.cpp#L684) will check whether normalisation parameters apply to events and automatically set this up for you.
**SetupSplines()** is discussed in the spline page of the wiki [here](https://github.com/mach3-software/MaCh3/wiki/5.-Splines).
**SetupFuncParameters()** setups up functional parameters. These will be experiment specific systematics and you will need to define this function. What you really need to do for these is keep track of the parameter value which then gets passed to a custom function to apply your systematic (this could be changing a kinematic variable, it could be scaling something etc.). A basic template is below:
```
void SetupFuncParameters(){
  funcParsIndex = XsecCov->GetFuncParsIndexFromDetID(SampleDetID);
  funcParsNames = XsecCov->GetFuncParsNamesFromDetID(SampleDetID);

  // The value of the functional parameter. Initialised to -999
  func_par_pos = -999;

  int func_it = 0;
  //Now find the correct index of xsec parameters to point to
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
	std::string name = funcParsNames.at(func_it);
	if (name == "FuncParNameInYaml") {
	  func_par_pos = *it;
	}
  }
}
```
The other simple Setup function is **SetupWeightPointers()** this is meant to fill the [total_weights_pointer](https://github.com/mach3-software/MaCh3/blob/develop/samplePDF/FDMCStruct.h#L39) of the fdmc_base struct. This is a convenient way to loop over all of the weights which you might want to apply to your events from different sources (xsec, detector, oscillation). See below for a dummy example:
```
void samplePDFSKBase::SetupWeightPointers() {
for (int i = 0; i < nOscillationChannels; ++i) {
  for (int j = 0; j < MCSamples[i].nEvents; ++j) {
    MCSamples[i].ntotal_weight_pointers[j] = 7;
    MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
    MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
    MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
  }
}
```

# Pure Virtual Functions
The samplePDFFDBase class is a virtual object and cannot be initialised on its own due to it having several pure virtual functions. It is these virtual functions which must be implemented in the experiment samplePDFFDBase class. Below is a list of the pure virtual functions that exist in the [samplePDFFDBase.h](https://github.com/mach3-software/MaCh3/blob/develop/samplePDF/samplePDFFDBase.h) file. These are functions which are likely to be very experiment specific and so you need to help the core code to grab key bits of information. Some of these we have already encountered.
```
virtual void SetupWeightPointers() = 0;
// Calculate the norm weight for a given event
virtual double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) = 0;
virtual double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) = 0;
virtual std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter) = 0;
```
SetupWeightPointers has already been described above.
ReturnKinematicParameter and ReturnKinematicParameterBinning are all functions related to getting specific variables from your experiment specific samplePDF and passing it to core. It is recommended that this is done through a set of enums and strings which match these. For example:
```
enum KinematicVariables{
  kTrueNeutrinoEnergy = 0,
  kRecoNeutrinoEnergy = 1,
  kNKinematicParameters
}

inline int ReturnKinematicParameterFromString(std::string KinematicParameterStr){
  if(KinematicParameterStr.find("TrueNeutrinoEnergy") != std::npos){return kTrueNeutrinoEnergy;}
  else if(KinematicParameterStr.find("RecoNeutrinoEnergy") != std::npos){return kRecoNeutrinoEnergy;}
  
  return kNKinematicParameters;
}
```
This will then make the ReturnKinematicParaeter and ReturnKinematicParameterBinning functions simpler.
```
double ReturnKinematicParameter(double KinematicVariable, int iEvent, int iOscChannel){
  //Convert double into an enum
  KinematicVariables Var = static_cast<KinematicVariables>(KinematicVariable);
  double Val;
  switch(Var){
    case(kTrueNeutrinoEnergy):
      Val = ExpObj->TrueEnu[i];
      break;
    case(kRecoNeutrinoEnergy):
      Val = ExpObj->RecoEnu[i];
      break;
  }

  return Val;
}
```
similarly then the ReturnKinematicParameter function can be passed a string
```
double ReturnKinematicParameter(std::string KinematicString, int iEvent, int iOscChannel){
  KinematicVariables Var = static_cast<KinematicVariables>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameter(Var, iEvent, iOscChannel);
}
```
It is this second function which allows you to go from a string in your sample or xsec config to an actual value to cut on. This way the core code can "know" which variable to get by reading in from a config file. This is aimed at making the code as flexible across different experiments as possible. The core code only needs to be told how to associate a string with an actual value.

Finally the ReturnKinematicBinning function can be handy if you want to create a histogram in different variables. You can do this by defining a binning to associate with a variable.
```
std::vector<double> ReturnKinematicParameter(double KinematicVariable){
  //Convert double into an enum
  KinematicVariables Var = static_cast<KinematicVariables>(KinematicVariable);
  std::vector<double> Binning;
  switch(Var){
    case(kTrueNeutrinoEnergy):
      double BinWidth = 50;//50MeV bin widths
      for(int ibin = 0 ; ibin < 10 ; ++ibin){
        Binning.push_back(BinWidth*ibin);
      }
      break;
    case(kRecoNeutrinoEnergy):
      double BinWidth = 40;//40MeV bin widths
      for(int ibin = 0 ; ibin < 10 ; ++ibin){
        Binning.push_back(BinWidth*ibin);
      }
      break;
  }
  return Binning;
}  
```