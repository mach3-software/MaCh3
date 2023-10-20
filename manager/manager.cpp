//******************************************************
// This class reads the configuration file config.cfg
// More documentation is available in the header.
//
// A brief description of each parameter is available
// in doc/executables_description/run_executable
// and the README of this folder
// look also at example_config in the main dir.
//******************************************************

#include "manager.h"
#include <iostream>

manager::manager(std::string const &filename)
    : config(YAML::LoadFile(filename)) {std::cout << "Setting config to be " << filename << std::endl; std::cout << "config is now " << config << std::endl;}

YAML::Node const & manager::raw(){ return config; }

/* Old Mananger that needs translation
 * we're moving to YAML BABY!

manager::manager(char *config, bool print) {

  verbosity = print;
  if (verbosity) std::cout << std::endl << "Reading configuration from " <<
config << std::endl;

  int success = readConfig(config);

  // Print the manager settings
  if(print){
    Print();
  }

  if (success == EXIT_FAILURE) {
    std::cerr << "Something went wrong in the manager, check above error message
and your config file" << std::endl; is_good_config = false; exit(EXIT_FAILURE);
  } else {
    is_good_config = true;
    if (verbosity) std::cout << "Succesfully read config " << config << ", now
proceeding" << std::endl;
  }

}

// Empty destructor
manager::~manager() {
}

// Read the supplied config file
// This is Getting pretty huge by now!
int manager::readConfig(char *config) {

  // Create the libconfig object
  libconfig::Config cfg;
  // Could turn on automatic type conversion in the future, but could be a bad
idea..
  // cfg.setAutoConvert(true);

  // Read the config file
  try {
    cfg.readFile(config);
  }

  // Do error reporting if readFile fails (not found or if parsing error)
  catch(const libconfig::FileIOException &fioex) {
    std::cerr << "I/O problem when reading " << config << std::endl;
    std::cerr << "Does the file exist?" << std::endl;
    return(EXIT_FAILURE);
  }

  // Catch parsing exception
  catch(const libconfig::ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() << "
- " << pex.getError() << std::endl; return(EXIT_FAILURE);
  }

  // Get the diferent entries in the config file
  // Catch any errors at the end of the try (all through libconfig)
  // Note, exception is only thrown on cfg.lookup
  // Any crucial settings should hence not be wrapped in the
if(cfg.exists("PARNAME")) try {
    // Output filename for the sampling
    if (cfg.exists("OUTPUTNAME")) {
      output_file = (const char*)cfg.lookup("OUTPUTNAME");
    } else {
      output_file = "mach3_output_file_default_name.root";
    }

    // Number of steps
    if (cfg.exists("NSTEPS")) {
      n_steps = cfg.lookup("NSTEPS");
    } else {
      n_steps = 10;
      if (verbosity) std::cout << "Number of steps not specified, using " <<
n_steps << std::endl;
    }

    //How often (number of steps) you want to autosave chain
    if (cfg.exists("AUTOSAVE")) n_AutoSave = cfg.lookup("AUTOSAVE");
    else n_AutoSave = 500;
    
    //Number of burn in steps
    if (cfg.exists("NBURNINSTEPS")){
      n_burnin_steps = cfg.lookup("NBURNINSTEPS");
    } else {
      n_burnin_steps = 1000;
      if (verbosity) std::cout << "Number of burn in steps not specified, using
" << n_burnin_steps << std::endl;
    }

    // Use Barlow Beeston likelihood in ND280?
    if (cfg.exists("MCSTAT")) {
      std::string likelihood = cfg.lookup("MCSTAT");
      if (likelihood == "Barlow-Beeston")                mc_stat_llh = TestStatistic(kBarlowBeeston);
      else if (likelihood == "IceCube")                  mc_stat_llh = TestStatistic(kIceCube);
      else if (likelihood == "Poisson")                  mc_stat_llh = TestStatistic(kPoisson);
      else if (likelihood == "Pearson")                  mc_stat_llh = TestStatistic(kPearson);
      else if (likelihood == "Dembinski-Abdelmottele")   mc_stat_llh = TestStatistic(kDembinskiAbdelmottele);
      else { 
        std::cerr << "Wrong form of test-statistic specified!" << std::endl;
        std::cerr << "You gave " << likelihood << " and I only support:" << std::endl;
        for(int i = 0; i < kNTestStatistics; i++)
        {
          std::cerr << TestStatistic_ToString(TestStatistic(i)) << std::endl;
        }
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
      if (verbosity) std::cout << "- MC stat test statistic specified, using " << TestStatistic_ToString(TestStatistic(mc_stat_llh)) << std::endl;
    } else {
      mc_stat_llh = kPoisson;
      if (verbosity) std::cout << "- MC stat test statistic not specified," << TestStatistic_ToString(TestStatistic(mc_stat_llh)) << std::endl;
    }
    
    if (cfg.exists("USE_UpdateW2")) {
      UpdateW2 = cfg.lookup("USE_UpdateW2");
    } else {
      UpdateW2 = true;
    }

    // Stepping scale for nuisance parameters
    if (cfg.exists("STEPSCALE")) {
      step_scale = cfg.lookup("STEPSCALE");
    } else {
      step_scale = 0.05;
    }

    // Stepping scale for xsec
    if (cfg.exists("STEPSCALEXSEC")) {
      step_scale_xsec = cfg.lookup("STEPSCALEXSEC");
    } else {
      step_scale_xsec = 0.05;
    }

    // Stepping scale for near detector systematics
    if (cfg.exists("STEPSCALENEARDET")) {
      step_scale_near_det = cfg.lookup("STEPSCALENEARDET");
    } else {
      step_scale_near_det = 0.04;
    }

    // Stepping scale for far detector systematics
    if (cfg.exists("STEPSCALEFARDET")) {
      step_scale_far_det = cfg.lookup("STEPSCALEFARDET");
    } else {
      step_scale_far_det = 1;
    }

    // Stepping scale for osc parameters
    if (cfg.exists("OSCSTEPSCALE")) {
      osc_step_scale = cfg.lookup("OSCSTEPSCALE");
    } else {
      osc_step_scale = 0.1;
    }

    // Neutrino POT
    if (cfg.exists("POT")) {
      protons_on_target = cfg.lookup("POT");
    } else {
      protons_on_target = 3.01E20;
      if (verbosity) std::cout << "- no POT specified, setting to 3.01E20" <<
std::endl;
    }

    // Anti-neutrino POT
    if (cfg.exists("NUBARPOT")) {
      nubar_protons_on_target = cfg.lookup("NUBARPOT");
    } else {
      nubar_protons_on_target = 3.01E20;
      if (verbosity) std::cout << "- no nubar POT specified, setting to 3.01E20"
<< std::endl;
    }

    // Eigen value threshold for PCA
    if (cfg.exists("PCA_THRESHOLD")) {
      pca_threshold = cfg.lookup("PCA_THRESHOLD");
    } else {
      pca_threshold = -1;
      if (verbosity) std::cout << "- no PCA threshold specified, setting to -1"
<< std::endl;
    }

    // Debug / Summary info
    if (cfg.exists("SUMMARY")) {
      debug_file = cfg.lookup("SUMMARY");
    } else {
      debug_file = true;
    }

    // Do we want to save nominal spectra
    if (cfg.exists("SAVENOM")) {
      save_nom = cfg.lookup("SAVENOM");
    } else {
      save_nom = false;
    }

    // posterior files (to use for posterior predictive)
    if (cfg.exists("STARTFROMPOS")) {
      start_from_pos = cfg.lookup("STARTFROMPOS");
    } else {
      start_from_pos = false;
    }

    // posterior files (to use for posterior predictive)
    if (cfg.exists("POSFILES")) {
      pos_files = (const char*)cfg.lookup("POSFILES");
    } else {
      pos_files = "EMPTY";
    }

    // DB X-axis Sample Binning
    if (cfg.exists("XAXISBINS")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["XAXISBINS"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        xaxisbins.push_back((double)setting[i]);
      }

    } else {
          std::cout << "[WARNING] - now x-axis bins specified. Using 0 and 1 as
default" << std::endl; xaxisbins.push_back(0); xaxisbins.push_back(1);
    }

    // DB Y-axis Sample Binning
    if (cfg.exists("YAXISBINS")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["YAXISBINS"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        yaxisbins.push_back((double)setting[i]);
      }

    } else {
          std::cout << "[WARNING] - now y-axis bins specified. Using 0 and 1 as
default" << std::endl; yaxisbins.push_back(0); yaxisbins.push_back(1);
    }

    // Binning option for splines
    if (cfg.exists("BINNINGOPT")) {
      binning_opt = cfg.lookup("BINNINGOPT");
    } else {
      binning_opt = 0.;
    }

    // Seed for our random number
    if (cfg.exists("SEED")) {
      random_seed = cfg.lookup("SEED");
    } else {
      random_seed = 0;
    }

    // Where is the fake data
    if (cfg.exists("FAKEDATAFILE")) {
      fake_data_file = (const char*)cfg.lookup("FAKEDATAFILE");
    } else {
      fake_data_file = "";
    }

    // Where is the data file
    if (cfg.exists("DATAFILE")) {
      data_file = (const char*)cfg.lookup("DATAFILE");
    } else {
      data_file = "";
    }

    // Fake data fit
    if (cfg.exists("FAKEDATAFIT")) {
      fake_data = cfg.lookup("FAKEDATAFIT");
    } else {
      fake_data = false;
    }

    // Asimov fit
    if (cfg.exists("ASIMOVFIT")) {
      asimov_fit = cfg.lookup("ASIMOVFIT");
    } else {
      asimov_fit = false;
    }

    // Data fit
    if (cfg.exists("REALDATAFIT")) {
      real_data = cfg.lookup("REALDATAFIT");
    } else {
      real_data = false;
    }

    // Turn on/off ND detector systematics
    if (cfg.exists("NEARDETSYSTOPT")) {
      near_det_syst_opt = (bool)cfg.lookup("NEARDETSYSTOPT");
    } else {
      near_det_syst_opt = true;
    }

    // Turn on/off Far detector systematics
    if (cfg.exists("FARDETSYSTOPT")) {
      far_det_syst_opt = (bool)cfg.lookup("FARDETSYSTOPT");
    } else {
      far_det_syst_opt = true;
    }

    // Turn off xsec systematics
    if (cfg.exists("XSECSYSTOPT")) {
      xsec_syst_opt = (bool)cfg.lookup("XSECSYSTOPT");
    } else {
      xsec_syst_opt = true;
    }

    // Set annealing temp (not used by standard!)
    if (cfg.exists("TEMP")) {
      temp = cfg.lookup("TEMP");
    } else if (cfg.exists("T")) {
      temp = cfg.lookup("T");
    } else {
      temp = -1;
    }

    // Ask if apply oscillation weight bins per bins (if not, event per event)
    if (cfg.exists("BINNEDOSC")) {
      binned_osc = (bool)cfg.lookup("BINNEDOSC");
    } else {
      binned_osc = false;
    }

    if(cfg.exists("OSCBINS")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["OSCBINS"];
      for (int i = 0; i < int(setting.getLength()); i++){
        osc_bins.push_back((double)setting[i]);
      }
    }

    // Ask if use sin2(2theta23) or sin2(2theta13)
    if (cfg.exists("DOUBLEANGLE")) {
      double_angle = (bool)cfg.lookup("DOUBLEANGLE");
    } else {
      double_angle = false;
    }

    // input oscillation parameters
    if (cfg.exists("OSCPARAM")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["OSCPARAM"];
      for (int i = 0; i < int(setting.getLength()); i++){
        osc_parameters.push_back((double)setting[i]);
      }
    }

    // input oscillation parameters if 2 sets
    if (cfg.exists("OSCPARAM_ANTINU")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["OSCPARAM_ANTINU"];
      for (int i = 0; i < int(setting.getLength()); i++){
        osc_parameters_antinu.push_back((double)setting[i]);
      }
    }

    // input oscillation parameters for Asimov fits
    // (if you don't want them to be the same as the starting values in the fit)
    if (cfg.exists("OSCPARAM_ASIMOVINPUT")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["OSCPARAM_ASIMOVINPUT"];
      for (int i = 0; i < int(setting.getLength()); i++){
        osc_parameters_asimovinput.push_back((double)setting[i]);
      }
    }

    // oscillation parameters to evaluate
    if (cfg.exists("OSCPARAMEVAL")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["OSCPARAMEVAL"];
      for (int i = 0; i < int(setting.getLength()); i++)
        osc_parameters_eval.push_back((int)setting[i]);
    }

    // Where is flux covariance
    if (cfg.exists("FLUXCOVFILE")) {
      flux_cov_file = (const char*)cfg.lookup("FLUXCOVFILE");
      if (verbosity) std::cout << "- found fluxCovMatrix = " << flux_cov_file <<
std::endl; } else { flux_cov_file =
std::string("inputs/flux_covariance_banff_13av2.root"); if (verbosity) std::cout
<< "- didn't find fluxCovMatrix! Setting default = " << flux_cov_file <<
std::endl;
    }

    // The name of the matrix to use
    if(cfg.exists("FLUXCOVNAME")) {
      flux_cov_name = (const char*)cfg.lookup("FLUXCOVNAME") ;
      if (verbosity) std::cout << "- found fluxCovMatrix name = " <<
flux_cov_name << std::endl; } else { flux_cov_name = "total_flux_cov"; if
(verbosity) std::cout << "- didn't find fluxCovMatrix! Setting default = " <<
flux_cov_name << std::endl;
    }

    // Where is xsec covariance
    if(cfg.exists("XSECCOVFILE")) {
      xsec_cov_file = (const char*)cfg.lookup("XSECCOVFILE");
      if (verbosity) std::cout << "- found xsecCovMatrix = " << xsec_cov_file <<
std::endl; } else {
          //TODO (ETA) - support not specifiying xsec cov if xsec systs are
turned off in xsec_syst_opt std::cerr << "[ERROR] - no xsec covariance matrix
specified which is required!" << std::endl; throw;
    }

    // The name of the matrix to use
    if(cfg.exists("XSECCOVNAME")) {
      xsec_cov_name = (const char*)cfg.lookup("XSECCOVNAME") ;
      if (verbosity) std::cout << "- found xsec covariance name = " <<
xsec_cov_name << std::endl; } else { xsec_cov_name = "xsec_cov"; if (verbosity)
std::cout << "- didn't find xsec covariance matrix name! Setting default = " <<
xsec_cov_name << std::endl;
    }

    // Look for specified ND280 covariance
    if (cfg.exists("NEARDETCOVFILE")) {
      near_det_cov_file = (const char*)cfg.lookup("NEARDETCOVFILE");
      if (verbosity) std::cout << "- found Near detector covariance matrix file
= " << near_det_cov_file << std::endl; } else {
          //TODO (ETA) - support not specifiying near detector covariance matrix
if near det systs are turned off in near_det_syst_opt std::cerr << "[ERROR] - no
near detector covariance matrix specified which is required!" << std::endl;
          throw;
    }

    // The name of the matrix to use
    if(cfg.exists("NEARDETCOVNAME")) {
      near_det_cov_name = (const char*)cfg.lookup("NEARDETCOVNAME") ;
      if (verbosity) std::cout << "- found Near detector covariance matrix name
= " << near_det_cov_name << std::endl; } else { near_det_cov_name =
"near_det_cov"; if (verbosity) std::cout << "- didn't find Near detector
covariance matrix name! Setting default = " << near_det_cov_name << std::endl;
    }

    // Where is Far detector covariance
    if(cfg.exists("FARDETCOVFILE")) {
      far_det_cov_file = (const char*)cfg.lookup("FARDETCOVFILE") ;
    } else {
          //TODO (ETA) - support not specifiying far detector covariance matrix
if far det systs are turned off in far_det_syst_opt std::cerr << "[ERROR] - no
far detector covariance matrix specified which is required!" << std::endl;
          throw;

    }

    // The name of the matrix to use
    if(cfg.exists("FARDETCOVNAME")) {
      far_det_cov_name = (const char*)cfg.lookup("FARDETCOVNAME") ;
      if (verbosity) std::cout << "- found Far detector covariance matrix name =
" << far_det_cov_name << std::endl; } else { far_det_cov_name = "far_det_cov";
      if (verbosity) std::cout << "- didn't find Far detector covariance matrix
name! Setting default = " << far_det_cov_name << std::endl;
    }

    // Where is osc covariance
    if(cfg.exists("OSCCOVFILE")) {
      osc_cov_file = (const char*)cfg.lookup("OSCCOVFILE") ;
      if (verbosity) std::cout << "- found oscCovMatrix = " << osc_cov_file <<
std::endl; } else { osc_cov_file =
"inputs/oscillation_covariance_6par_nondouble.root"; if (verbosity) std::cout <<
"- didn't find oscCovMatrix! Setting default = " << osc_cov_file << std::endl;
    }

    // The name of the matrix to use
    if(cfg.exists("OSCCOVNAME")) {
      osc_cov_name = (const char*)cfg.lookup("OSCCOVNAME") ;
      if (verbosity) std::cout << "- found oscCovMatrix name = " << osc_cov_name
<< std::endl; } else { osc_cov_name = "osc_cov"; if (verbosity) std::cout << "-
didn't find oscCovMatrix! Setting default = " << osc_cov_name << std::endl;
    }

    // Where is osc covariance (if we need a different one for antinu)
    if(cfg.exists("OSCCOVFILE_ANTINU")) {
      osc_cov_antinu_file = (const char*)cfg.lookup("OSCCOVFILE_ANTINU") ;
      if (verbosity) std::cout << "- found oscCovMatrix_antinu = " <<
osc_cov_antinu_file << std::endl; } else { osc_cov_antinu_file = "osc_cov";
    }

    // The name of the matrix to use
    if(cfg.exists("OSCCOVNAME_ANTINU")) {
      osc_cov_antinu_name = (const char*)cfg.lookup("OSCCOVNAME_ANTINU") ;
      if (verbosity) std::cout << "- found oscCovMatrix_antinu name = " <<
osc_cov_antinu_name << std::endl; } else { osc_cov_antinu_name = "";
    }

    // Are we using reactor constraints on 13 oscillation params
    if(cfg.exists("USERC")) {
      rc_on = cfg.lookup("USERC") ;
    } else {
      rc_on = false;
    }

    // Location of data files (requires passing manager object to samplePDF)
    if (cfg.exists("DATALOC")) {
      data_location = (const char*)(cfg.lookup("DATALOC"));
    } else if (std::getenv("MACH3_DATA") != NULL) {
      data_location = std::string(std::getenv("MACH3_DATA"));
    } else {
      data_location = "./";
    }

    // Location of MC files (requires passing manager object to samplePDF)
    if (cfg.exists("MCLOC")) {
      mc_location = (const char*)(cfg.lookup("MCLOC"));
    } else if (std::getenv("MACH3_MC") != NULL) {
      mc_location = std::string(std::getenv("MACH3_MC"));
    } else {
      mc_location = std::string("./");
    }

    // Is this statistics only fit?
    if (cfg.exists("STATONLY")) {
      is_stat_only = cfg.lookup("STATONLY");
      if (verbosity) std::cout << "- found STAT only, setting " << is_stat_only
<< std::endl; } else { is_stat_only = 0; if (verbosity) std::cout << "- didn't
find STATONLY, setting default " << is_stat_only << std::endl;
    }

    // Do oscillation?
    if (cfg.exists("DOOSC")) {
      do_osc = cfg.lookup("DOOSC");
      if (verbosity) std::cout << "- found DOOSC, setting " << do_osc <<
std::endl; } else { do_osc = 0; if (verbosity) std::cout << "- didn't find
DOOSC, setting default " << do_osc << std::endl;
    }

    // Is this statistics only fit?
    if (cfg.exists("USEFARSYST")) {
      use_far_syst = cfg.lookup("USEFARSYST");
      if (verbosity) std::cout << "- found USEFARSYST, setting " << use_far_syst
<< std::endl; } else { use_far_syst = true; if (verbosity) std::cout << "-
didn't find USEFARSYST, setting default " << use_far_syst << std::endl;
    }

    // What far detector params have flat priors ("array" in libconfig,
std::vector<int> in class) if (cfg.exists("FARDETPARAMFLAT")) { const
libconfig::Setting &root = cfg.getRoot(); const libconfig::Setting &setting =
root["FARDETPARAMFLAT"]; for (int i = 0; i < int(setting.getLength()); i++)
        far_det_flat.push_back((int)setting[i]);
    }

    // What far detector params are fixed
    if (cfg.exists("FARDETPARAMFIX")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["FARDETPARAMFIX"];
      for (int i = 0; i < int(setting.getLength()); i++)
        far_det_fix.push_back((int)setting[i]);
    }

    // What xsec params are flat
    if (cfg.exists("XSECPARAMFLAT")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["XSECPARAMFLAT"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        xsec_flat.push_back((int)setting[i]);
      }
    }

    // What xsec params are fixed
    if (cfg.exists("XSECPARAMFIX")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["XSECPARAMFIX"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        xsec_fix.push_back((int)setting[i]);
      }
    }

    // DB Xsec Param Step Sizes
    if (cfg.exists("XSECSTEPSIZE")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["XSECSTEPSIZE"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        xsec_stepsize.push_back((double)setting[i]);
      }
    } else {
      xsec_stepsize.push_back(-1);
    }

    // Use random start for the sample?
    if (cfg.exists("RANDOM_START")) {
      random_start = cfg.lookup("RANDOM_START");
    } else {
      random_start = false;
    }

    // Process MCMC in mcmc code?
    if (cfg.exists("PROCESS_MCMC")) {
      process_mcmc = cfg.lookup("PROCESS_MCMC");
    } else {
      process_mcmc = false;
    }

    // Decide whether or not to apply beta to the nuebar appearance probability
    // Default is not to use it
    if(cfg.exists("USEBETA")){
      use_beta = cfg.lookup("USEBETA");
        } else {
      use_beta = false;
        }

    // Decide whether beta should flip (true) or be continuous (false)
    // Default is continuous
    if(cfg.exists("FLIPBETA")){
      flip_beta = cfg.lookup("FLIPBETA");
        } else{
      flip_beta = false;
        }

    // Decide whether beta should be applied to nue (true) appearance
probability instead of nuebar
    // Default is nuebar (false)
    if(cfg.exists("APPLYBETANUE"))
      apply_beta_nue = cfg.lookup("APPLYBETANUE");
    else
      apply_beta_nue = false;

    // Decide wheter you should apply (1/beta) to nue appearance probability and
beta to nuebar appearance probability (true)
    // Default is to just apply beta to nuebar (false)
    if(cfg.exists("APPLYBETADIAG"))
      apply_beta_diag = cfg.lookup("APPLYBETADIAG");
    else
      apply_beta_diag = false;

    if (cfg.exists("sample_name")) {
      sample_name = (std::string) cfg.lookup("sample_name").c_str();
    } else {
      sample_name = "Sample";
    }

    if (cfg.exists("samplenumber")) {
      sample_number = (int)cfg.lookup("samplenumber");
    } else {
      sample_number = 1;
    }

    if (cfg.exists("SampleDetID")) {
      sample_det_id = cfg.lookup("SampleDetID");
    } else {
      sample_det_id = -1;
    }

    if (cfg.exists("kinematic_cut_vars")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["kinematic_cut_vars"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        kinematic_cut_vars.push_back((std::string)setting[i].c_str());
      }
    }

    if (cfg.exists("kinematic_cut_low_bnd")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["kinematic_cut_low_bnd"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        kinematic_cut_low_bnd.push_back((std::string)setting[i].c_str());
      }
    }

    if (cfg.exists("kinematic_cut_up_bnd")) {
      const libconfig::Setting &root = cfg.getRoot();
      const libconfig::Setting &setting = root["kinematic_cut_up_bnd"];
      for (int i = 0; i < int(setting.getLength()); i++) {
        kinematic_cut_up_bnd.push_back((std::string)setting[i].c_str());
      }
    }

    if (cfg.exists("up_bnd")) {
      up_bnd = cfg.lookup("up_bnd");
    } else {
      up_bnd = -999;
    }

  } // end the try

  // Catch if the config file doesn't contain some crucial parameter which we
haven't wrapped in exist catch (const libconfig::SettingNotFoundException &nfex)
{ std::cerr << "Setting " << nfex.getPath() << " not found in file " << config
<< std::endl; return(EXIT_FAILURE);
  }

  // Catch invalid names in config
  catch (const libconfig::SettingNameException &nex) {
    std::cerr << "Invalid setting name in " << nex.getPath() << std::endl;
    return(EXIT_FAILURE);
  }

  // If invalid type conversion (e.g. char* to double) - Catches if we give
FLUXCOVFILE=-1 or suttin catch (const libconfig::SettingTypeException &tyex) {
    std::cerr << "Invalid type conversion in " << tyex.getPath() << std::endl;
    return(EXIT_FAILURE);
  }

  // Get settings defined by pre-processor directives, e.g. CPU MP and GPU
#if defined(MULTITHREAD)
  cpu_mp_on = true;
  n_cpus = omp_get_max_threads();
#else
  cpu_mp_on = false;
  n_cpus = 0;
#endif

#if defined(CUDA)
  gpu_on = true;
#else
  gpu_on = false;
#endif

  // Check that all the settings make sense for us
  //if (!checkSettings()) {
    //std::cerr << "Invalid settings combination" << std::endl;
    //return(EXIT_FAILURE);
  //}

  return(EXIT_SUCCESS);
}

void manager::Print() {

  std::cout << "---------------------------------" << std::endl;
  std::cout << "General settings      " << std::endl;
  std::cout << "    Output file:      " << GetOutputFilename() << std::endl;
  std::cout << "    Seed:             " << GetSeed() << std::endl;
  std::cout << "    Debug:            " << GetDebug() << std::endl;
  std::cout << "    GPU fit:          " << GetGPU() << std::endl;
  std::cout << "    CPU MP fit:       " << GetCPUMP() << std::endl;
  std::cout << "       N cores:       " << GetNCPU() << std::endl;

  std::cout << "---------------------------------" << std::endl;
  std::cout << "MCMC settings         "    << std::endl;
  std::cout << "    N steps:          " << GetNSteps() << std::endl;
  std::cout << "    AutoSave:         " << GetAutoSave() << std::endl;
  std::cout << "    Xsec step scale:  " << GetXsecStepScale() << std::endl;
  std::cout << "    ND det step scale:" << GetNearDetStepScale() << std::endl;
  std::cout << "    Far det step scale:" << GetFarDetStepScale() << std::endl;
  std::cout << "    Osc step scale:   " << GetOscStepScale() << std::endl;
  std::cout << "    Anneal. temp:     " << GetTemp() << std::endl;
  std::cout << "    Random start:     " << GetRandomStart() << std::endl;
  std::cout << "    Process MCMC:     " << GetProcessMCMC() << std::endl;
  std::cout << "    Posterior start:  " << GetStartFromPosterior() << std::endl;
  std::cout << "    Posterior chain:  " << GetPosteriorFiles() << std::endl;

  std::cout << "---------------------------------" << std::endl;
  std::cout << "Run settings      " << std::endl;
  std::cout << "    Nnu POT:        " << GetPOT() << std::endl;
  std::cout << "    Antinu POT:    " << GetNubarPOT() << std::endl;
  std::cout << "    Save nominal:     " << GetSaveNom() << std::endl;
}

*/
