// Class for interfacing with Eb and Xsec Apps
#include "BinnedSplineMakerBase.h"

BinnedSplineMakerBase::BinnedSplineMakerBase(const std::string& Config) {
  Setting = std::make_unique<Manager>(Config.c_str());

  auto& options = Setting->raw()["settings"];
  auto ModeName = Get<std::string>(options["modes"], __FILE__ , __LINE__);
  Modes = std::make_unique<MaCh3Modes>(ModeName);
  UseKnotForWeight = GetFromManager<bool>(options["UseKnot"], true, __FILE__, __LINE__);
  int nSplineModes = 0;
  KinematicParameters = nullptr;
  ReversedKinematicParameters = nullptr;
  for (int mode = 0; mode < Modes->GetNModes(); ++mode) {
    const std::string modeName = Modes->GetMaCh3ModeName(mode);
    const std::string suffix =
    Modes->GetSplineSuffixFromMaCh3Mode(mode);

    // Reuse the spline mode if this suffix already exists.
    auto it = fSuffixToSplineMode.find(suffix);

    int splineMode;

    if (it == fSuffixToSplineMode.end()) {
      splineMode = nSplineModes++;
      fSuffixToSplineMode.emplace(suffix, splineMode);
    } else {
      splineMode = it->second;
    }

    // Map MaCh3 mode -> spline mode.
    fModeToSplineMode.emplace(mode, splineMode);

    MACH3LOG_INFO("MaCh3 mode {} ({}) -> spline suffix '{}' (mode {})", mode, modeName, suffix, splineMode);
  }
}

void BinnedSplineMakerBase::SetupBinning() {
  // Okay now we can grab the etrue bins if we provide a template file, else use a default setting.
  std::vector<double> etrue_bins = GetETrueBinning();
  auto binning_opts = Setting->raw()["binning_opts"];
  // Okay next we need to deal with reco energy bins! (var1)
  // This really really ought to be configurable!
  // HW : Should now just work, Var1 bins now set up correctly
  auto var1_binning = GetFromManager<std::vector<double>>(binning_opts["var1_bins"], {-1e15, 1e15}, __FILE__, __LINE__);
  auto var2_binning = GetFromManager<std::vector<double>>(binning_opts["var2_bins"], {-1e15, 1e15}, __FILE__, __LINE__);

  MACH3LOG_INFO("Setting binning");
  SetBinning(etrue_bins, var1_binning, var2_binning);
  MACH3LOG_INFO("Binning set");
}

void BinnedSplineMakerBase::LoadSystematics() {
  const auto& systematics = Setting->raw()["systematics"];
  std::vector<std::string> dial_names;
  std::vector<int> dial_nominal_knots;
  std::vector<std::vector<int>> dial_modes;
  std::vector<std::string> dial_output_names;

  // Next we deal with our systematics
  MACH3LOG_INFO("Dial Names :");
  for (const auto& table : systematics) {
    std::string name = table.first.as<std::string>();

    // Access the values for this systematic
    const auto& sys = systematics[name];

    // Let's user write splines with name that's not the dial name
    if(sys["output_name"]) {
      dial_output_names.push_back(sys["output_name"].as<std::string>());
    } else{
      dial_output_names.push_back(name);
    }

    int nom_knot = sys["nominal_knot"].as<int>();
    auto interaction_mode_names = sys["interaction_modes"].as<std::vector<std::string>>();

    dial_names.push_back(name);

    std::vector<int> systematic_dial_modes;
    systematic_dial_modes.reserve(interaction_mode_names.size());
    const bool include_all = std::find(interaction_mode_names.begin(), interaction_mode_names.end(), "all") != interaction_mode_names.end();

    if (include_all) {
      if (interaction_mode_names.size() == 1) {
        // Add every unique spline mode
        systematic_dial_modes.reserve(fSuffixToSplineMode.size());

        for (const auto& [suffix, spline_mode] : fSuffixToSplineMode) {
          systematic_dial_modes.push_back(spline_mode);
        }
      } else {
        throw MaCh3Exception(__FILE__, __LINE__, "'all' is currently supported only when capacity == 1");
      }
    } else {
      systematic_dial_modes.reserve(interaction_mode_names.size());

      for (const auto& suffix : interaction_mode_names) {
        const auto it = fSuffixToSplineMode.find(suffix);

        if (it == fSuffixToSplineMode.end()) {
          throw MaCh3Exception(__FILE__, __LINE__, "No spline mode found for suffix: " + suffix);
        }

        systematic_dial_modes.push_back(it->second);
      }
    }
    dial_modes.push_back(std::move(systematic_dial_modes));
    dial_nominal_knots.push_back(nom_knot);
    MACH3LOG_INFO("  {}", name);
  }

  // Okay let's grab all of our dial stuff
  systematicProperties.resize(dial_names.size());

  int numSplines = 0;
  for(size_t iSyst = 0; iSyst < dial_names.size(); iSyst++){
    systematicProperties[iSyst].systName = dial_names[iSyst];
    systematicProperties[iSyst].outputSystName = dial_output_names[iSyst];
    systematicProperties[iSyst].nominalPosition = dial_nominal_knots[iSyst];
    systematicProperties[iSyst].intModes = dial_modes[iSyst];
    numSplines += static_cast<int>(systematicProperties[iSyst].intModes.size());
  }

  // Set up the histogram vector:
  dev_tmp.resize(numSplines);
  MACH3LOG_INFO("Total splines loaded: {}", numSplines);
}

void BinnedSplineMakerBase::Initialise()
{
  const auto& binning_opts = Setting->raw()["binning_opts"];
  auto binvar_1_str = binning_opts["binvar1"].as<std::string>();
  auto binvar_2_str = binning_opts["binvar2"].as<std::string>();

  binTypes_vec.insert(binTypes_vec.end(),
                      {ReturnKinematicParameterFromString(binvar_1_str),
                       ReturnKinematicParameterFromString(binvar_2_str)});
  // Currently we can only bin in (at most) etrue, var1, var2
  if(binTypes_vec.size() != 2){
    MACH3LOG_ERROR("Wrong number of binning types, got {}, must provide 2!", binTypes_vec.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Select binning types
  if (binTypes_vec[1] != 0) {
    MACH3LOG_INFO("Using 2D splines");
  } else {
    MACH3LOG_INFO("Using 1D splines");
  }
  // Load syst info from config
  LoadSystematics();
  // Load info about knots [virtual]
  LoadKnots();
  // Actually add systematics
  AddSystematics();
  // Initialise variables and inputs [virtual]
  LoadVariables();
  // Set binning and histograms
  SetupBinning();
}

void BinnedSplineMakerBase::MakeVariations() {
  // Let's make ourselves some splines!
  MACH3LOG_INFO("Filling systematic properties");
  // Fill histograms with weight information (one bin for each etrue, var1, [var2])
  FillSystematicProperties();
  // Fit a TGraph to each bin
  MACH3LOG_INFO("Filling graphs");
  FillGraphs();
  // Fit a spline to each TGraph
  MACH3LOG_INFO("Filling splines");
  FillSplines();
  MACH3LOG_INFO("Everything's filled, writing to file!");
  WriteSplines();
}

// It's messy but it works!
void BinnedSplineMakerBase::AddSystematics()
{
  // Fill vectors before setting branch addresses, because the vectors will reallocate memory
  // while being filled. We need to know the final address to read in the weight file
  int numSysts = static_cast<int>(systematicProperties.size());
  MACH3LOG_INFO("Adding Systematics");

  int total_hist_iterator = 0;
    
  for(int iSyst = 0; iSyst < numSysts; iSyst++){
    for(size_t iMode = 0; iMode < systematicProperties[iSyst].intModes.size(); iMode++){

      if(systematicProperties[iSyst].knotArray == nullptr) {
        MACH3LOG_ERROR("Knots for {} are not defined, you should play with LoadKnots", systematicProperties[iSyst].systName);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      dev_tmp[total_hist_iterator].resize(systematicProperties[iSyst].knotArray->GetSize());
      total_hist_iterator++;

      const int spline_mode = systematicProperties[iSyst].intModes[iMode];
      std::string splineSuffix;
      for (const auto& [suffix, id] : fSuffixToSplineMode) {
        if (id == spline_mode) {
          splineSuffix = suffix;
          break;
        }
      }
      const std::string tmpname = "dev." + systematicProperties[iSyst].systName + "." + splineSuffix;
      dev_tmp_names.push_back(std::move(tmpname));
    }
  }
}


BinnedSplineMakerBase::~BinnedSplineMakerBase() {
}

void BinnedSplineMakerBase::FillGraphs() {
  int sysMode = 0;

  // loop over all systematics:
  for(size_t iSyst = 0; iSyst < systematicProperties.size(); iSyst++)
  {
    // loop over all interactions modes this systematic affects
    for(size_t iMode = 0; iMode < systematicProperties[iSyst].intModes.size(); iMode++)
    {
      int iter = 0;
      // fill the graphs with at the supplied (from the weight file) knot locations
      // Need to loop over all the bins and knots
      for(int iKnot = 0; iKnot < systematicProperties[iSyst].knotArray->GetSize(); iKnot++)
      {
        for(int xBin = 0; xBin < dev_tmp[sysMode+iMode][iKnot]->GetNbinsX(); xBin++)
        {
          for(int yBin = 0; yBin < dev_tmp[sysMode+iMode][iKnot]->GetNbinsY(); yBin++ )
          {
            for(int zBin = 0; zBin < dev_tmp[sysMode+iMode][iKnot]->GetNbinsZ(); zBin++)
            {
              int nom = systematicProperties[iSyst].nominalPosition;
              const double nomContent = dev_tmp[sysMode + iMode][nom]->GetBinContent(xBin+1,yBin+1,zBin+1);
              const double knotContent = dev_tmp[sysMode + iMode][iKnot]->GetBinContent(xBin+1,yBin+1,zBin+1);
              graphs[sysMode+iMode][xBin][yBin][zBin]->SetPoint(iter, systematicProperties[iSyst].knotArray->At(iKnot),
                                                                (nomContent > 0. ? knotContent / nomContent : 1.));
              // Don't want to fill our splines with junk!
              if(std::isnan((nomContent)) || nomContent < 0) {
                MACH3LOG_ERROR("2D: is nan at nom");
              }

              if(std::isnan(knotContent) || knotContent < 0) {
                const int spline_mode = systematicProperties[iSyst].intModes[iMode];
                std::string splineSuffix;
                for (const auto& [suffix, id] : fSuffixToSplineMode) {
                  if (id == spline_mode) {
                    splineSuffix = suffix;
                    break;
                  }
                }
                MACH3LOG_ERROR("ERROR 2D: is nan at knot a: {}", iKnot);
                MACH3LOG_ERROR("Systematic is: {} : {}", iSyst, systematicProperties[iSyst].outputSystName);
                MACH3LOG_ERROR("Are you returning a NaN? {}", std::isnan(knotContent));
                MACH3LOG_ERROR("Do you have a negative bin content? {}", knotContent);
                MACH3LOG_ERROR("Currently working on mode: {}", spline_mode);
              }
            }
          }
        }
      iter++;
      }
    }
    // To make sure we keep going over the correct modes
    sysMode += static_cast<int>(systematicProperties[iSyst].intModes.size());
  }
}

void BinnedSplineMakerBase::FillSystematicProperties()
{
  // Event loop
  for (Long64_t jentry = 0; jentry < GetEntries(); jentry++)
  {
    GetEntry(jentry);

    // Okay now we deal with our spline binning
    double Enu = ReturnKinematicParameter("TrueNeutrinoEnergy", jentry);
    double binning_par1 = ReturnKinematicParameter(binTypes_vec[0], jentry);
    double binning_par2 = ReturnKinematicParameter(binTypes_vec[1], jentry);

    auto Mode = GetMode(jentry);
    auto interaction_mode = Modes->GetModeFromGenerator(Mode);
    auto it = fModeToSplineMode.find(interaction_mode);
    if (it == fModeToSplineMode.end()) {
      MACH3LOG_ERROR("No spline mode for interaction mode: {}", std::to_string(interaction_mode));
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    auto spline_mode = it->second;
    int sysMode = 0;

    // loop over all systematics:
    for (size_t iSyst = 0; iSyst < systematicProperties.size(); iSyst++)
    {
      // loop over all interactions modes this systematic affects
      for (size_t iMode = 0; iMode < systematicProperties[iSyst].intModes.size(); iMode++)
      {
        // is this event one of these modes?
        if (systematicProperties[iSyst].intModes[iMode] == spline_mode)
        {
          for (int iKnot = 0; iKnot < systematicProperties[iSyst].knotArray->GetSize(); ++iKnot)
          {
            double wgt = 0;
            if(UseKnotForWeight){
              wgt = GetWeight(jentry, iSyst, iKnot);
            } else {
              wgt = GetWeight(jentry, iSyst, systematicProperties[iSyst].knotArray->At(iKnot));
            }

            if (std::isnan(wgt) || wgt < 0) {
              MACH3LOG_WARN("weight {}, Enu: {}, binning_par1: {}, binning_par2: {}, mode: {}",
                             wgt, Enu, binning_par1, binning_par2, Mode);
            }

            dev_tmp[sysMode + iMode][iKnot]->Fill(Enu, binning_par1, binning_par2, wgt);
          }
        }
      }
      sysMode += static_cast<int>(systematicProperties[iSyst].intModes.size());
    } // end loop over syst
  } // Event Loop
} // end function

void BinnedSplineMakerBase::FillSplines(){
  for(int iHist = 0; iHist < static_cast<int>(dev_tmp.size()); iHist++){
    for(int xBin = 0; xBin < dev_tmp[iHist][0]->GetNbinsX(); xBin++){
      std::vector< std::vector< std::unique_ptr<TSpline3> > > tmp_xbin;
      for(int yBin = 0; yBin < dev_tmp[iHist][0]->GetNbinsY(); yBin++){
        std::vector<std::unique_ptr<TSpline3> > tmp_ybin;
        for(int zBin = 0; zBin < dev_tmp[iHist][0]->GetNbinsZ(); zBin++){
          char sname[1000];
          sprintf(sname,"%s.sp.%d.%d.%d",dev_tmp_names[iHist].c_str(),xBin,yBin,zBin);
          // Literally fits a spline to each TGraph
          auto tmp_zbin = std::make_unique<TSpline3>(sname, graphs[iHist][xBin][yBin][zBin].get());
          tmp_ybin.push_back(std::move(tmp_zbin));
        }
        tmp_xbin.push_back(std::move(tmp_ybin));
      }
      if(splines.size() == 0) splines.resize(dev_tmp.size());
      splines[iHist].push_back(std::move(tmp_xbin));
    }
  }
}

// Etrue-Erec binning (used for SK splines)
// (Well any 2D binning, for 1D splines we just set rebins and nrebins = 0 and we're grand!)
void BinnedSplineMakerBase::SetBinning(const std::vector<double>& var0Bins,
                                       const std::vector<double>& var1Bins,
                                       const std::vector<double>& var2Bins)
{
  // ROOT wants number of bins - 1!
  int nVar0Bins = static_cast<int>(var0Bins.size()-1);
  int nVar1Bins = static_cast<int>(var1Bins.size()-1);
  int nVar2Bins = static_cast<int>(var2Bins.size()-1);

  MACH3LOG_INFO("VAR 0 BINS : {{ {:.2f} }}", fmt::join(var0Bins, " "));
  MACH3LOG_INFO("VAR 1 BINS : {{ {:.2f} }}", fmt::join(var1Bins, " "));
  MACH3LOG_INFO("VAR 2 BINS : {{ {:.2f} }}", fmt::join(var2Bins, " "));

  for(int iHist = 0; iHist < static_cast<int>(dev_tmp.size()); iHist++)  {
    for(int iKnot = 0; iKnot < static_cast<int>(dev_tmp[iHist].size()); iKnot++) {
      const std::string hname = "dev_tmp." + std::to_string(iKnot) + "." + std::to_string(iHist);

      // This trick is gonna save us a big loop!
      // Vectors are contiguous so double*s are able to access them like this!
      // [.data() gets a double* from an std::vector which is the format ROOT wants]
      dev_tmp[iHist][iKnot] = std::make_unique<TH3F>(hname.c_str(), hname.c_str(), nVar0Bins, var0Bins.data(), nVar1Bins,var1Bins.data(), nVar2Bins, var2Bins.data());
      dev_tmp[iHist][iKnot]->SetDirectory(nullptr);
    }
  }

  // Might as well allocate memory now !
  for(int iHist = 0; iHist < static_cast<int>(dev_tmp.size()); iHist++){
    for(int iVar0 = 0; iVar0<nVar0Bins; iVar0++) {
      std::vector< std::vector< std::unique_ptr<TGraph> > > tmp_var0_bin;
      for(int iVar1 = 0; iVar1<nVar1Bins; iVar1++) {
        std::vector< std::unique_ptr<TGraph> > tmp_var1_bin;
        for(int iVar2 = 0; iVar2 < nVar2Bins; iVar2++) {
          auto tmp_var2_bin = std::make_unique<TGraph>(); // Dummy graph for later
          tmp_var1_bin.push_back(std::move(tmp_var2_bin));
        }
        tmp_var0_bin.push_back(std::move(tmp_var1_bin));
      }
      if(graphs.size() == 0) {
        graphs.resize(dev_tmp.size());
      }
      graphs[iHist].push_back(std::move(tmp_var0_bin));
    }
  }
}

void BinnedSplineMakerBase::WriteSplines() {
  auto outputname = Setting->raw()["settings"]["outputfile"].as<std::string>();

  TFile* fout = new TFile(outputname.c_str(), "RECREATE");
  fout->cd();

  // MaCh3 doesn't really need these anymore! Keeping in for VALOR
  int rwbin = 0;

  // We  only use the hists to get binning info but not saving them breaks something somewhere!
  // YSP: I have only seen the usage of dev_tmp.0.0 to retrieve spline binning, so let's just not save any other histograms
  const std::string hname = "dev_tmp.0.0";
  auto tmp_hist = M3::Clone(dev_tmp[0][0].get());
  tmp_hist->Write(hname.c_str());

  // Write the actual splines
  for(size_t iSyst = 0; iSyst < systematicProperties.size(); iSyst++)
  {
    for(size_t iMode = 0; iMode < systematicProperties[iSyst].intModes.size(); iMode++)
    {
      for(size_t iVar0 = 0; iVar0 < splines[rwbin+iMode].size(); iVar0++)
      {
        for(size_t iVar1 = 0; iVar1 < splines[rwbin+iMode][iVar0].size(); iVar1++)
        {
          for(size_t iVar2 = 0; iVar2 < splines[rwbin+iMode][iVar0][iVar1].size(); iVar2++)
          {
            const int spline_mode = systematicProperties[iSyst].intModes[iMode];
            std::string splineSuffix;
            for (const auto& [suffix, id] : fSuffixToSplineMode) {
              if (id == spline_mode) {
                splineSuffix = suffix;
                break;
              }
            }

            TString spname = Form("dev.%s.%s.sp.%zu.%zu.%zu", systematicProperties[iSyst].outputSystName.c_str(),
                                  splineSuffix.c_str(), iVar0, iVar1, iVar2);
            auto tmp_spline = static_cast<TSpline3*>(splines[rwbin+iMode][iVar0][iVar1][iVar2]->Clone());
            tmp_spline->Write(spname);
            delete tmp_spline;
          }
        }
      }
    }
    rwbin += static_cast<int>(systematicProperties[iSyst].intModes.size());
  }

  TMacro ConfigSave = YAMLtoTMacro(Setting->raw(), "Settings");
  ConfigSave.Write();

  fout->Close();
  MACH3LOG_INFO("Splines written to {}", outputname);
}

// ************************************************
int BinnedSplineMakerBase::ReturnKinematicParameterFromString(const std::string& KinematicParameterStr) const {
// ************************************************
  auto it = KinematicParameters->find(KinematicParameterStr);
  if (it != KinematicParameters->end()) return it->second;

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameterStr);
  throw MaCh3Exception(__FILE__, __LINE__);

  return M3::_BAD_INT_;
}

// ************************************************
std::string BinnedSplineMakerBase::ReturnStringFromKinematicParameter(const int KinematicParameter) const {
// ************************************************
  auto it = ReversedKinematicParameters->find(KinematicParameter);
  if (it != ReversedKinematicParameters->end()) {
    return it->second;
  }

  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameter);
  throw MaCh3Exception(__FILE__, __LINE__);

  return "";
}
