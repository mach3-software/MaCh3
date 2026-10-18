#include "BinnedSplineMakerSK.h"

BinnedSplineMakerSK::BinnedSplineMakerSK(const std::string& Config) : BinnedSplineMakerBase(Config)
{
  is_atm = GetFromManager<bool>(Setting->raw()["settings"]["is_atm"], false, __FILE__, __LINE__);
  erec_mev = GetFromManager<bool>(Setting->raw()["binning_opts"]["erec_mev"], false, __FILE__, __LINE__);

  const auto& binning_opts = Setting->raw()["binning_opts"];
  auto binvar_1_str = binning_opts["binvar1"].as<std::string>();
  auto binvar_2_str = binning_opts["binvar2"].as<std::string>();
  binTypes_vec.insert(binTypes_vec.end(), {GetBinningType(binvar_1_str), GetBinningType(binvar_2_str)});

  Initialise();
  //Hacky way to get normalisation weight;
  auto kinematics_file = Setting->raw()["settings"]["mtuplefile"].as<std::string>();
  TFile *_sampleFile = new TFile(kinematics_file.c_str(), "READ");
  TH1D* norm = _sampleFile->Get<TH1D>("norm");
  if(!norm) {
    MACH3LOG_ERROR("Cannot find norm histogram!!, setting weight to 1");
    wgtnvv = 1.0;
  } else {
    wgtnvv = norm->GetBinContent(1)/norm->GetBinContent(2);
  }
  delete norm;
  delete _sampleFile;
}


BinnedSplineMakerSK::~BinnedSplineMakerSK() {
  if(fChain_sk) delete fChain_sk;
}

std::vector<double> BinnedSplineMakerSK::GetETrueBinning() {
  // Next job is to setup binning ...ugh!
  // Right for now we're going to stick with the horrible method used previously
  // ....should make this a config!
  int ip = GetParticleID(0);  //Let's get our sample + current information from the file name
  bool isfhc = false;
  auto kinematics_file = Setting->raw()["settings"]["mtuplefile"].as<std::string>();

  if(kinematics_file.find("fhc") != std::string::npos) isfhc = true;
  const auto& binning_opts = Setting->raw()["binning_opts"];
  auto atm_sample_num = GetFromManager<int>(binning_opts["atm_sample_num"], -1, __FILE__, __LINE__);
  // Firstly, let's grab the etrue bins!
  std::string binning_axis = "";
  if (is_atm) {
    binning_axis = "XSecBinning_Sample_" + std::to_string(atm_sample_num);
  } else {
    if(isfhc){
      switch(ip){
        case(kNueBar):
          binning_axis = "sk_numode_nueb_bins";
          break;
        case(kNue):
          binning_axis = "sk_numode_nue_bins";
          break;
        case(kNumuBar):
          binning_axis = "sk_numode_numub_bins";
          break;
        case(kNumu):
          binning_axis = "sk_numode_numu_bins";
          break;
      };
      if(kinematics_file.find("numu_x_nue") != std::string::npos) {
        binning_axis = "sk_numode_numu_bins";
      }
      if(kinematics_file.find("numubar_x_nuebar") != std::string::npos) {
        binning_axis="sk_numode_numub_bins";
      }
    } else {
      switch(ip){
        case(kNueBar):
          binning_axis = "sk_anumode_nueb_bins";
          break;
        case(kNue):
          binning_axis = "sk_anumode_nue_bins";
          break;
        case(kNumuBar):
          binning_axis = "sk_anumode_numub_bins";
          break;
        case(kNumu):
          binning_axis = "sk_anumode_numu_bins";
          break;
      };
      if(kinematics_file.find("numu_x_nue") != std::string::npos){
        binning_axis = "sk_anumode_numu_bins";
      }

      if(kinematics_file.find("numubar_x_nuebar") != std::string::npos){
        binning_axis = "sk_anumode_numub_bins";
      }
    }
  } // end is not ATM

  std::string mach3_root = getenv("MaCh3_ROOT");
  if(mach3_root == "") {
    MACH3LOG_ERROR("Sorry you've not sourced setup.sh!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  auto etrue_bin_file = GetFromManager<std::string>(binning_opts["true_e_bins_template"], "", __FILE__, __LINE__);
  std::vector<double> etrue_bins;
  if(etrue_bin_file.empty()) {
    MACH3LOG_INFO("Using a default Etrue binning since no template file was provided");
    etrue_bins = {0.0, 0.4, 0.5, 0.6, 0.8, 1.2, 2.0, 3.0, 5.0, 30.0};
  } else {
    etrue_bin_file = mach3_root +"/"+ etrue_bin_file;
    etrue_bins = ExtractBinningFromTemplate(etrue_bin_file, binning_axis, is_atm);
  }
  return etrue_bins;
}

void BinnedSplineMakerSK::LoadKnots()
{
  auto& options = Setting->raw()["settings"];
  auto weightFile = options["weightfile"].as<std::string>();
  auto kinematics_file = options["mtuplefile"].as<std::string>();
  auto kinematics_tree_name = options["mtupletree"].as<std::string>();

  // Make a TChain using the mtuple
  fChain_sk = new TChain(kinematics_tree_name.c_str());
  if(fChain_sk) {
    fChain_sk->Add(kinematics_file.c_str());
  } else {
    MACH3LOG_ERROR("Couldn't find mtuple in your kinematics file, bye!");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int numSysts = static_cast<int>(systematicProperties.size());
  MACH3LOG_INFO("Adding Systematics");
  weightArray.resize(numSysts);
  for(int iSyst = 0; iSyst < numSysts; iSyst++)
  {
    std::string sysName      = systematicProperties[iSyst].systName;
    std::string sysTreeName  = sysName + "_tree";

    TChain  *tempChain = new TChain(sysTreeName.c_str());
    tempChain->Add(weightFile.c_str());
    fChain_sk->AddFriend(sysTreeName.c_str(), weightFile.c_str());

    weightArray[iSyst]  = nullptr;
    systematicProperties[iSyst].knotArray    = nullptr;

    std::string weightBrName = sysName + "_weights";
    std::string knotBrName   = sysName + "_knots";

    fChain_sk->SetBranchAddress(weightBrName.c_str(), &(weightArray[iSyst]));
    fChain_sk->SetBranchAddress(knotBrName.c_str(), &(systematicProperties[iSyst].knotArray));
  }
  fChain_sk->GetEntry(0);
}

void BinnedSplineMakerSK::LoadVariables()
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Note that for the for the weight chain, this is done in AddSystematics()

  // Set branch addresses and branch pointers
  fCurrent_sk = -1;
  fChain_sk->SetMakeClass(1);

  fChain_sk->SetBranchAddress("iclass", &iclass);
  fChain_sk->SetBranchAddress("Ibound", &Ibound);
  fChain_sk->SetBranchAddress("erec", &E_rec);
  fChain_sk->SetBranchAddress("wgtflx", &wgtflx);
  fChain_sk->SetBranchAddress("wgtosc", &wgtosc);
  fChain_sk->SetBranchAddress("fq1rdir", fq1rdir);
  fChain_sk->SetBranchAddress("fqmome", &fqmome);
  fChain_sk->SetBranchAddress("fqmomm", &fqmomm);
  //ETA - adding in fitqun multi-ring variables for numuCC1pi samples(s)
  fChain_sk->SetBranchAddress("ipp", &ipp);
  //YSP - fiTQun multi-ring e-like fit ID
  fChain_sk->SetBranchAddress("iep", &iep);
  fChain_sk->SetBranchAddress("fqmrdir", &fqmrdir);
  fChain_sk->SetBranchAddress("fqmrmom", &fqmrmom);
  fChain_sk->SetBranchAddress("fqmreloss", &fqmreloss);
  fChain_sk->SetBranchAddress("mode", &mode);
  fChain_sk->SetBranchAddress("ipnu", ipnu);
  fChain_sk->SetBranchAddress("pnu", pnu);

  if (is_atm) {
    fChain_sk->SetBranchAddress("genmom", &genmom);
    fChain_sk->SetBranchAddress("gencz", &gencz);

    fChain_sk->SetBranchAddress("FFF", &FFF);
    fChain_sk->SetBranchAddress("oscweight3f", &oscweight3f);
    fChain_sk->SetBranchAddress("solarweight", &solarweight);
  } else {
    genmom = 0;
    gencz = 0;
  }
}

double BinnedSplineMakerSK::GetVariableValue(int binTypes) {
  double var = 0;
  constexpr double SKNuDir[3] = {0.669764, -0.742179, 0.024223};

  double pmax_lep = 0.; // momentum of most energetic lepton
  double costheta = -999.;
  double erec = 0;

  if (!is_atm) {
    switch (iclass)
    // Get kinematics for each event class [iclass in mtuples = sample type]
    {
      case k1Re:
        costheta = fq1rdir[0][1][0] * SKNuDir[0] + fq1rdir[0][1][1] * SKNuDir[1] + fq1rdir[0][1][2] * SKNuDir[2];
        pmax_lep = fqmome;
        break;

      case k1Rmu:
        costheta = fq1rdir[0][2][0] * SKNuDir[0] + fq1rdir[0][2][1] * SKNuDir[1] + fq1rdir[0][2][2] * SKNuDir[2];
        pmax_lep = fqmomm;
        break;

      case k1Re1de:
        costheta = fq1rdir[0][1][0] * SKNuDir[0] + fq1rdir[0][1][1] * SKNuDir[1] + fq1rdir[0][1][2] * SKNuDir[2];
        pmax_lep = fqmome;
        break;

      case kMRe1de:
        costheta = fqmrdir[iep][0][0] * SKNuDir[0] + fqmrdir[iep][0][1] * SKNuDir[1] + fqmrdir[iep][0][2] * SKNuDir[2];
        pmax_lep = fqmrmom[iep][0];
        break;

      case kMRMu1de:
      case kMRMu2de:
        // Okay now we do the slightly more complex task of calculating erec for mr-muon samples
        int idx_mu;
        idx_mu = 0;
        if (fqmrmom[ipp][0] < fqmrmom[ipp][1])
        {
          idx_mu = 1;
        }
        costheta = fqmrdir[ipp][idx_mu][0] * SKNuDir[0] + fqmrdir[ipp][idx_mu][1] * SKNuDir[1] + fqmrdir[ipp][idx_mu][2] * SKNuDir[2];
        pmax_lep = fqmrmom[ipp][idx_mu];
        break;

      default:
        MACH3LOG_ERROR("CAUTION : UNKNOWN ICLASS, {} THIS BREAKS ME!", iclass);
        throw MaCh3Exception(__FILE__, __LINE__);
    };
  }

  switch (binTypes)
  {
    case kErec:
      // Erec calculated in MeV
      erec = CalculateErec(iclassName(iclass), pmax_lep, costheta);
      // HW does this need a valor check?
      if ((std::abs(erec / E_rec - 1) > 1E-2)) {
        MACH3LOG_ERROR("For the sample with iClass: {}, calculated Erec of {:.2f} but in minituples it's {:.2f}", iclass, erec, E_rec);
      }

      float erec_unit; // MeV
      if (erec_mev) {
        erec_unit = 1000.0; // MeV
      } else {
        erec_unit = 1.0; // GeV
      }
      var = erec * erec_unit;
      break;
    case kMom:
      var = pmax_lep * 1000; // MeV
      break;
    case kTheta:
      var = std::acos(costheta) / TMath::Pi() * 180.;
      break;
    case kGenmom:
      var = genmom;
      break;
    case kGencz:
      var = gencz;
      break;
    case kNoBinning:
      var = 0.0; // Dummy bin
      break;
    default:
      MACH3LOG_ERROR("Can only bin parameter 1 in erec (0) or momentum (1), white you passed {}", binTypes);
      throw MaCh3Exception(__FILE__, __LINE__);
  };
  return var;
}


// loop over all knots
double BinnedSplineMakerSK::GetWeight([[maybe_unused]] Long64_t jentry, size_t iSyst, int iKnot)
{
  const double knotWeight = weightArray[iSyst]->At(iKnot);

  double wgt = wgtnvv;
  if (is_atm) {
    constexpr double livetime_s = 3244.4 / (365.25 * 500.);
    wgt *= oscweight3f * livetime_s * solarweight * FFF * knotWeight;
  } else {
    wgt *= wgtflx * wgtosc * knotWeight;
  }
  return wgt;
}


int BinnedSplineMakerSK::GetBinningType(const std::string& binName) {
  if(binName == "none")         return kNoBinning;
  else if(binName == "erec")    return kErec;
  else if(binName == "theta")   return kTheta;
  else if(binName == "mom")     return kMom;
  else if(binName == "genmom")  return kGenmom;
  else if(binName == "gencz")   return kGencz;
  else {
    MACH3LOG_ERROR("Sorry: {} is not a recognised bin type!", binName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return 0;
}

std::vector<double> BinnedSplineMakerSK::ExtractBinningFromTemplate(std::string template_file_name, std::string  binning_name, bool atm) {
  // Simple function to get bins from a template file
  TFile* template_file = new TFile(template_file_name.c_str());
  MACH3LOG_INFO("Getting bins from {}", template_file_name);
  if(template_file->IsZombie()){
    MACH3LOG_ERROR("Couldn't open template file: {}", template_file_name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  TAxis* bin_axis;
  if (atm) {
    TH3D* h = template_file->Get<TH3D>(binning_name.c_str());
    bin_axis = h->GetXaxis();
  } else {
    bin_axis = template_file->Get<TAxis>(binning_name.c_str());
    if(!bin_axis) {
      MACH3LOG_ERROR("Couldn't find: {} in {}", binning_name, template_file_name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  const double* bin_range = bin_axis->GetXbins()->GetArray();
  std::vector<double> return_bins;
  for(int i = 0; i <= bin_axis->GetNbins(); i++){
    return_bins.push_back(bin_range[i]);
  }
  return return_bins;
}
