#ifndef __MANAGER_H__
#define __MANAGER_H__

#ifndef EXIT_FAILURE
#define EXIT_FAILURE -1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

// ROOT include
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"


//******************************************************
// This class reads the configuration file config.cfg
// More documentation is available in the README
//******************************************************


#include "yaml-cpp/yaml.h"

class manager {

public:
  manager(std::string const &);
  YAML::Node const &raw();

private:
  YAML::Node config;

/*
 * Keeping this all here as we should start adding in functions here to do reading
 * of different options
 public:
  manager(char *config, bool print = true);
  virtual ~manager();

  int readConfig(char *config);

// Print and checks
  void Print();
  //bool checkSettings();

  //void SaveSettings(TFile * const OutputFile);

// Lots of yummy Get functions
  const char *GetOutputFilename() { return (const char*)output_file.c_str(); }

  int GetMCStatLLH() { return mc_stat_llh; }
  int GetUpdateW2() { return UpdateW2; }
  
  double GetPOT()       { return protons_on_target; }
  double GetNubarPOT()  { return nubar_protons_on_target; }

  int GetNSteps() { return n_steps; }
  int GetBurnInSteps() { return n_burnin_steps; }

  int GetAutoSave() { return n_AutoSave; }
  
  double GetStepScale()     { return step_scale; }
  double GetXsecStepScale() { return step_scale_xsec; }
  double GetNearDetStepScale(){ return step_scale_near_det;}
  double GetFarDetStepScale(){ return step_scale_far_det;}
  double GetOscStepScale()  { return osc_step_scale;}

  bool GetRandomStart() { return random_start; }
  bool GetProcessMCMC() { return process_mcmc; }

  bool GetDebug() { return debug_file; }
  int GetSeed()   { return random_seed; }

  double GetPCAThreshold() { return pca_threshold; }

  bool GetGPU()   { return gpu_on; }
  bool GetCPUMP() { return cpu_mp_on; }
  int GetNCPU()   { return n_cpus; }

  bool GetBinnedOsc()   { return binned_osc;}
  bool GetDoubleAngle() { return double_angle; }
  std::vector<double> GetOscBins() { return osc_bins; }

  double GetBinningOpt() { return binning_opt; }

  bool GetFakeDataFitFlag() { return fake_data; }
  bool GetFakeDataFluctuatedFitFlag() { return fake_data_fluctuated; }
  bool GetRealDataFitFlag() { return real_data; }
  bool GetAsimovFitFlag()   { return asimov_fit; }
  bool GetNDAsimovFlag()    { return nd280_asimov; }
  bool GetMakeAsimovOnly()  { return make_asimov_only; }

  bool GetToyFitFlag()  { return toy_experiment_fit; }
  int GetNtoy()         { return toy_number; }

  std::string const & GetToyFilename() const { return toy_experiment_file; }
  std::string const & GetDataFilename() const { return data_file; }
  std::string const & GetFakeDataFilename() const { return fake_data_file; }

  // ! This single function is deprecated but used in old code!
  bool GetDetSystOpt() { std::cerr << "ERROR YOU'RE USING A DEPRECATED PARAMETER!" << std::endl;
                         std::cerr << "PLEASE USE GetNearDetSystOpt or GetFarDetSystOpt INSTEAD!" << std::endl;
                         throw;
                       };

  //bool GetNearDetSystOpt() { return near_det_syst_opt; }
  //bool GetFarDetSystOpt() { return far_det_syst_opt; }
  bool GetXsecSystOpt()  { return xsec_syst_opt; }

  bool GetSaveNom() { return save_nom; }

  double GetTemp() { return temp; }

  const char* GetNDRuns() { return nd_runs.c_str(); }

  // Getters for systematic parameters
  std::vector<int> GetNearDetFlat() { return near_det_flat; }
  std::vector<int> GetNearDetFix()  { return near_det_fix; }
  std::vector<int> GetFarDetFlat() { return far_det_flat; }
  std::vector<int> GetFarDetFix()  { return far_det_fix; }
  std::vector<int> GetXsecFlat()  { return xsec_flat; }
  std::vector<int> GetXsecFix()   { return xsec_fix; }
  std::vector<int> GetOscFlat()   { return osc_flat; }
  std::vector<int> GetOscFix()    { return osc_fix; }
  //std::vector<int> GetOscFlat_antinu()  { return osc_flat_antinu; }
  //std::vector<int> GetOscFix_antinu()   { return osc_fix_antinu; }

  std::vector<double> GetXsecIndivStepSize() { return xsec_stepsize; }

  std::vector<double> GetOscParameters()        { return osc_parameters; }
  std::vector<double> GetOscParameters_antinu() { return osc_parameters_antinu; }
  std::vector<double> GetOscParameters_asimovinput() { return osc_parameters_asimovinput; }
  std::vector<int> GetOscParamEval()            { return osc_parameters_eval; }

  //Binning options for FD samples
  //DB Generic axis binniing
  std::vector<double> GetXAxisBins() { return xaxisbins;}
  std::vector<double> GetYAxisBins() { return yaxisbins;}

  //DB Adding verbosity for no print out on sample cfgs
  bool verbosity;

  std::vector<std::string> Getkinematic_cut_vars() {return kinematic_cut_vars;}
  std::vector<std::string> Getkinematic_cut_low_bnd() {return kinematic_cut_low_bnd;}
  std::vector<std::string> Getkinematic_cut_up_bnd() {return kinematic_cut_up_bnd;}

  int GetSampleNumber() {return sample_number;}

  const char* GetXsecCovMatrix() { return (const char*)xsec_cov_file.c_str(); }
  const char* GetNearDetCovMatrix() { return (const char*)near_det_cov_file.c_str(); }
  const char* GetFarDetCovMatrix() { return (const char*)far_det_cov_file.c_str(); }
  const char* GetOscCovMatrix()         { return (const char*)osc_cov_file.c_str(); }
  const char* GetOscCovMatrix_antinu()  { return (const char*)osc_cov_antinu_file.c_str(); }

  const char* GetXsecCovName() { return (const char*)xsec_cov_name.c_str(); }
  const char* GetNearDetCovName() { return (const char*)near_det_cov_name.c_str(); }
  const char* GetFarDetCovName() { return (const char*)far_det_cov_name.c_str(); }
  const char* GetOscCovName() { return (const char*)osc_cov_name.c_str(); }
  const char* GetOscCovName_antinu() { return (const char*)osc_cov_antinu_name.c_str(); }
  
  bool GetStartFromPosterior() { return start_from_pos; }
  std::string GetPosteriorFiles() { return pos_files; }
  
  std::string GetPolyFile() { return poly_file; }

  bool GetDoOsc() { return do_osc; }

  bool GetRC() { return rc_on; }

  bool GetStatOnly() { return is_stat_only; }

  bool GetGoodConfig() { return is_good_config; }

  bool GetUseBeta() { return use_beta; }
  bool GetFlipBeta() { return flip_beta; }
  bool GetApplyBetaNue() { return apply_beta_nue; }
  bool GetApplyBetaDiag(){ return apply_beta_diag; }

  bool GetFullLLH() {return full_llh;}

  std::vector<int> GetNuType() {return nu_type;}
  std::vector<int> GetOscNuType() {return oscnu_type;}
  std::vector<bool> GetSignal() {return signal;}
  
 protected:
  // fit flags
  // POT
  double protons_on_target;
  double nubar_protons_on_target;

  // Livetime
  double livetime;
  // nsteps
  int n_steps;
  int n_burnin_steps;
  
  //KS: Sets how frequent (number of steps) you want to autosave your chain, if you autosave rarely chain will be sliglhy faster, if you job wall is small you might want more freqeunt autosave to reduce number of lost steps
  int n_AutoSave;
  
  // cooling temperature(for simulated annealing)
  double temp;

  // stepsize
  //double step_size;
  // step scale for nuisance parameters
  double step_scale;
  double step_scale_xsec;
  double step_scale_near_det;
  double step_scale_far_det;
  // osc step scale
  double osc_step_scale;

  // outfile name
  std::string output_file;

  // summary and debug log file
  bool debug_file;
  bool full_llh;

  // random seed
  int random_seed;

  // GPU and CPU settings
  bool gpu_on;
  bool cpu_mp_on;
  int n_cpus;

  // toy data fit
  bool fake_data;
  // poisson fluc of asimov
  bool fake_data_fluctuated;
  // Data fit
  bool real_data;
  // Asimov fit
  bool asimov_fit;
  bool nd280_asimov;
  // Make asimov only (no fit)
  bool make_asimov_only;
  // Toy experiment fit
  bool toy_experiment_fit;
  // toy number
  int toy_number;
  
  //binned fd predictions
  bool binned_fd_pred;

  // toy exp file name
  std::string toy_experiment_file;
  std::string data_file;
  std::string fake_data_file;
  std::string pos_files;

  std::string poly_file;

  bool start_from_pos;

  //ND runs string
  std::string nd_runs;

  std::vector<double> osc_parameters ;
  std::vector<double> osc_parameters_antinu ;
  std::vector<double> osc_parameters_asimovinput ;
  std::vector<int> osc_parameters_eval ;
  std::vector<double> osc_bins ;

  // Priors
  std::vector<int> xsec_flat;
  std::vector<int> near_det_flat;
  std::vector<int> far_det_flat;
  std::vector<int> osc_flat;
  std::vector<int> osc_flat_antinu;
  
  std::vector<int> xsec_fix;
  std::vector<int> near_det_fix;
  std::vector<int> far_det_fix;
  std::vector<int> osc_fix;
  std::vector<int> osc_fix_antinu;

  std::vector<double> xsec_stepsize;

  // Binning option for splines
  double  binning_opt ;

  // Handling new Kinematic cuts
  std::vector<std::string> kinematic_cut_vars;
  std::vector<std::string> kinematic_cut_low_bnd;
  std::vector<std::string> kinematic_cut_up_bnd;

  // Binning used for samples in samplePDF object
  std::vector<double> xaxisbins;
  std::vector<double> yaxisbins;

  //General inputs
  std::string xsec_cov_file ;
  std::string near_det_cov_file ;
  std::string far_det_cov_file ;
  std::string osc_cov_file ;
  std::string osc_cov_antinu_file ;
 
  //name of the useful matrices in the above files
  std::string xsec_cov_name ;
  std::string near_det_cov_name ;
  std::string far_det_cov_name ;
  std::string osc_cov_name ;
  std::string osc_cov_antinu_name ;

  // Location of data and MC
  std::string data_location;
  std::string mc_location;

  // Use reactor constraint ?
  bool rc_on ;

  // Do oscillation (used in sigma variation)?
  bool do_osc;

  // Statistics only fit?
  bool is_stat_only;

  // Statistics only fit?
  bool use_far_syst; //option to run without far detector systematics for step size tuning

  // Near Detector systematics option (0: off, 1: on)
  bool near_det_syst_opt;
  // Far Detector systematics option (0: off, 1: on)
  bool far_det_syst_opt;
  // Cross section systematics option (0: off, 1: on)
  bool xsec_syst_opt;

  // Start chain in random position?
  bool random_start;

  // Process mcmc
  bool process_mcmc;

  bool binned_osc;
  bool double_angle;

  // Bool to hold if read config file is legitimate
  // Currently there's an EXIT_FAILURE when read fails but we could do stuff here if we're interested
  bool is_good_config;

  // Do we save nominal at ND280 (could implement at FD)
  bool save_nom;
  
  // Use beta in nuebar oscillation probability?
  bool use_beta;

  // Flip beta between 0 and 1? (If false, beta is continuous)
  bool flip_beta;

  // Apply beta to nue appearance probability (instead of nuebar appearance probability)?
  bool apply_beta_nue;

  // Apply (1/beta) to nue appearance probability and beta to nuebar appearance probability?
  bool apply_beta_diag;

  // What PCA threshold we apply for ND280
  double pca_threshold;

  // Apply Barlow Beeston
  int mc_stat_llh;

  //Whether you want to update W2 in Likelihood calcaution
  int UpdateW2;
  
  // Variables to read in sample config files
  // these are thnigs that typically need to be used in samplePDF
  std::string sample_name;
  int sample_number;
  int sample_det_id;
  double up_bnd;
  std::vector<int> sample_vecno;
  std::vector<int> nu_type;
  std::vector<int> oscnu_type;
  std::vector<bool> signal;

  */
};

#endif
