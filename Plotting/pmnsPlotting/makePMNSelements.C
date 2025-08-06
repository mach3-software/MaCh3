#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLine.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <sstream>

void makePMNSelements(bool RCreweight = true, std::string fit_type = "AsimovA22"){

  unsigned int n_chains = 69;

  std::string chain = "T2K";
  std::string fitter = "MaCh3";

  int n_bins = 5000;

  // Does the chain include burn-in?
  bool incl_burnin;
  if (fit_type != "OA2020"){
    incl_burnin = true;
    std::cout<<"Cut burnin..."<<std::endl;
  }
  else{
    incl_burnin = false;
  }
  unsigned int n_steps_burnin = 150000;

  bool doModulus = true;
  bool doReal = false;
  bool doImag = false;
  bool doPrior = true;
  bool do2D = false;
  bool do2D_pmns = false;

  bool solarConstraint = true;
  bool sinDelCPprior = false;

  int prior_case = 0;

  bool haarPrior = false;
  bool T2Kprior = true;
  bool NOvAprior = false;

  //prior map has following structure: prior_map[prior_case][variable] = flat_in, where prior case goes from 0 to 8
  //variable goes from 0 to 2 with: 0 = delta_CP, 1 = theta13, 2 = theta23
  //flat_in goes from 0 to 2 with: 0 = flat in the variable, 1 = flat in sin(variable), 2 = flat in cos(variable)
  int prior_map[9][3] = {{0,0,0}, {1,1,1}, {1,1,2}, {1,2,1}, {2,1,1}, {2,2,1}, {2,1,2}, {1,2,2}, {2,2,2}};

  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};
  std::map<std::string, std::string> infile_str{{"AsimovA22","oa2021/AsimovA22_NonUniform_041022_reduced_smeared_reweighted.root"},{"AsimovB22","oa2021/AsimovB22_NonUniform_Full_161222_reduced_smeared_reweighted.root"},{"OA2021","oa2021/Data_NonUniform_021222_Full_reduced_reweighted_smeared.root"},{"OA2020","MaCh3-OA2020_ALL_data_RCweights_reduced_burnInCut_smearedRevision_withPionSI.root"}};
  std::map<std::string, std::string> infile_str_dune{{"FDOnly","dune/DUNE_FDOnly_reduced.root"},{"FD+ND","dune/DUNE_NuPhys_wNDDet_reduced.root"}};

  TFile* outFile;
  std::stringstream filename;

  if (chain == "NOvA" && fitter == "Aria") {
    std::cout<<"Chain: T2K+NOvA Aria"<<std::endl;
    filename<<"pmns_"<<chain<<"_"<<fitter<<"_Asimov_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors_priorCase"<<prior_case<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  if (chain == "NOvA" && fitter == "MaCh3") {
    std::cout<<"Chain: T2K+NOvA MaCh3"<<std::endl;
    filename<<"pmns_"<<chain<<"_"<<fitter<<"_DataFit_postBANFF_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors_T2KPrior_"<<bool_str[solarConstraint]<<"SolarConstraint_"<<bool_str[sinDelCPprior]<<"SinDelCPprior_fullChain_"<<n_bins<<"bins.root";//priorCase"<<prior_case<<".root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  else if (chain == "T2K" && fitter == "Aria"){
    std::cout<<"Chain: Aria T2K-only"<<std::endl;
    filename<<"pmns_"<<chain<<"_"<<fitter<<"_Asimov1_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors_priorCase"<<prior_case<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  else {
    std::cout<<"Chain: "<<chain<<";  Fit type: "<<fit_type<<std::endl;
    filename<<"pmns_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors_priorCase"<<prior_case<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  std::cout<<"Creating output file: "<<filename.str().c_str()<<std::endl;

/*
  TTree* post = (TTree*)post_file->Get("posteriors");
  TChain* post = new TChain("posteriors");
  post->Add("/gpfs/projects/McGrewGroup/kwood/chains_OA2020/dataFit/MaCh3-OA2020_dataFit_cutBurnIn_1.root");
  post->Add("/gpfs/projects/McGrewGroup/kwood/chains_OA2020/dataFit/MaCh3-OA2020_dataFit_cutBurnIn_2.root");
  post->Add("/gpfs/projects/McGrewGroup/kwood/chains_OA2020/dataFit/MaCh3-OA2020_dataFit_cutBurnIn_3.root");
  post->Add("/gpfs/projects/McGrewGroup/kwood/chains_OA2020/dataFit/MaCh3-OA2020_dataFit_cutBurnIn_4.root");
*/

  TChain* post;
  std::stringstream chainpath;

  if (chain == "NOvA" && fitter == "Aria") {
    post = new TChain("samples");
    chainpath<<"/vols/t2k/users/mp721/data/MCMCSamples_merged_novabifrost_joint_both_systs_th13_th23_dcp_dm32_masshier_noreactor_grid_nddata_overall_asimov4.root";
    post->Add(chainpath.str().c_str());
  }
  else if (chain == "NOvA" && fitter == "MaCh3") {
    //post = new TChain("posteriors");
    post = new TChain("osc_posteriors");
    for(unsigned int n = 0;n<n_chains; n++)
    {
      chainpath.str("");
      chainpath.clear();
      chainpath<<"/vols/t2k/users/ea2817/T2K_NOvA/DataFits/IndividualChains/MaCh3-NOvAT2KDataFits_09062023-BinnedOscCPU-postBANFF-woRC_cedar_ch"<<n+1<<"_red_reweighted.root";
      std::cout<<"Adding chain: "<<chainpath.str().c_str()<<std::endl;
      post->Add(chainpath.str().c_str());
    }
    std::cout<<"Done!"<<std::endl;
  }
  else if (chain == "T2K" && fitter == "Aria") {
    post = new TChain("samples");
    chainpath<<"/vols/t2k/users/mp721/data/mcmcsamples_aria_t2konly_syst_burnin_asimov1_thinned1_postubrnin20000.root";
    post->Add(chainpath.str().c_str());
  }
  else if (chain == "DUNE"){
    post = new TChain("osc_posteriors");
    chainpath<<"/vols/t2k/users/mp721/data/"<<infile_str_dune[fit_type];
    post->Add(chainpath.str().c_str());
  }
  else {
    post = new TChain("osc_posteriors");
    chainpath<<"/vols/t2k/users/mp721/data/"<<infile_str[fit_type];
    post->Add(chainpath.str().c_str());
  }
  std::cout<<"Reading input chain: "<<chainpath.str()<<std::endl;

  unsigned int nSteps = post->GetEntries();

  std::cout << "number of mcmc steps = " << nSteps << std::endl;
  std::cout << "prior case: "<< prior_case << std::endl;

  double th13, th23;
  double s2th13, s2th23, s2th12, dcp, dm2, rc;
  double s22th13, s22th23, prob, mh;
  double prior_s2th13, prior_c2th13, prior_c4th13, prior_s2th23, prior_s2th12, prior_dcp, prior_sindcp;
  double prior_th12, prior_th13, prior_wRC_th13, prior_th23;
  double rw;//, rw_prior;
  //double prior_rw_dCP[3], prior_rw_th13[3], prior_rw_th23[3], prior_rw_dCP_prior[3], prior_rw_th13_prior[3], prior_rw_th23_prior[3];
  int step;
/*
  post->SetBranchAddress("sin2th_13",&s2th13);
  post->SetBranchAddress("sin2th_23",&s2th23);
  post->SetBranchAddress("sin2th_12",&s2th12);
  post->SetBranchAddress("delta_cp",&dcp);
  post->SetBranchAddress("delm2_23",&dm2);
  post->SetBranchAddress("RCreweight",&rc);
  post->SetBranchAddress("step",&step);
*/

  if (fitter == "Aria") {
    // post->SetBranchAddress("th13",&s22th13);
    // post->SetBranchAddress("th23",&s22th23);
    // post->SetBranchAddress("theta12",&s2th12);
    post->SetBranchAddress("th13",&th13);
    post->SetBranchAddress("th23",&th23);
    post->SetBranchAddress("delta(pi)",&dcp);
    post->SetBranchAddress("dmsq32",&dm2);
    post->SetBranchAddress("logprob",&prob);
    post->SetBranchAddress("MH",&mh);
    //post->SetBranchAddress("step",&step);
  }
  else if (fitter == "MaCh3" && chain == "NOvA"){
    // post->SetBranchAddress("sin2th_13",&s2th13);
    // post->SetBranchAddress("sin2th_23",&s2th23);
    // post->SetBranchAddress("sin2th_12",&s2th12);
    // post->SetBranchAddress("delta_cp",&dcp);
    // post->SetBranchAddress("delm2_23",&dm2);
    post->SetBranchAddress("theta13",&s2th13);
    post->SetBranchAddress("theta23",&s2th23);
    post->SetBranchAddress("theta12",&s2th12);
    post->SetBranchAddress("dcp",&dcp);
    if (RCreweight){
      post->SetBranchAddress("RCreweight",&rc);
    }
    post->SetBranchAddress("step",&step);
  }
  else {
    post->SetBranchAddress("theta13",&s2th13);
    post->SetBranchAddress("theta23",&s2th23);
    post->SetBranchAddress("dcp",&dcp);
    post->SetBranchAddress("dm23",&dm2);
    post->SetBranchAddress("step",&step);
    if (chain != "DUNE"){
      post->SetBranchAddress("theta12",&s2th12);
      post->SetBranchAddress("RCreweight",&rc);
    }
  }

  

  double s13, s23, s12, sdcp, cdcp, c13, c23, c12, j;
  double ue1, ue2, ue3, umu1, umu2, umu3, utau1, utau2, utau3;
  double real_ue1, real_ue2, real_ue3, real_umu1, real_umu2, real_umu3, real_utau1, real_utau2, real_utau3;
  double imag_ue1, imag_ue2, imag_ue3, imag_umu1, imag_umu2, imag_umu3, imag_utau1, imag_utau2, imag_utau3;
  //TComplex comp_ue1, comp_ue2, comp_ue3, comp_umu1, comp_umu2, comp_umu3, comp_utau1, comp_utau2, comp_utau3;

  double s13_prior, s13_wRC_prior, s23_prior, s12_prior, sdcp_prior, cdcp_prior, c13_prior, c13_wRC_prior, c23_prior, c12_prior;
  double ue1_prior, ue2_prior, ue3_prior, umu1_prior, umu2_prior, umu3_prior, utau1_prior, utau2_prior, utau3_prior;
  double real_ue1_prior, real_ue2_prior, real_ue3_prior, real_umu1_prior, real_umu2_prior, real_umu3_prior, real_utau1_prior, real_utau2_prior, real_utau3_prior;
  double imag_ue1_prior, imag_ue2_prior, imag_ue3_prior, imag_umu1_prior, imag_umu2_prior, imag_umu3_prior, imag_utau1_prior, imag_utau2_prior, imag_utau3_prior;
  //TComplex comp_ue1_prior, comp_ue2_prior, comp_ue3_prior, comp_umu1_prior, comp_umu2_prior, comp_umu3_prior, comp_utau1_prior, comp_utau2_prior, comp_utau3_prior;

/*
  TTree* pmnsTree = new TTree("pmnsTree","pmnsTree");
  pmnsTree->Branch("ue1",&ue1);
  pmnsTree->Branch("ue2",&ue2);
  pmnsTree->Branch("ue3",&ue3);
  pmnsTree->Branch("umu1",&umu1);
  pmnsTree->Branch("umu2",&umu2);
  pmnsTree->Branch("umu3",&umu3);
  pmnsTree->Branch("utau1",&utau1);
  pmnsTree->Branch("utau2",&utau2);
  pmnsTree->Branch("utau3",&utau3);
  pmnsTree->Branch("RCweight",&rc);
*/

  TH1D* h_ue1;
  TH1D* h_ue2;
  TH1D* h_ue3;
  TH1D* h_umu1;
  TH1D* h_umu2;
  TH1D* h_umu3;
  TH1D* h_utau1;
  TH1D* h_utau2;
  TH1D* h_utau3;
  TH1D* h_ue1_prior;
  TH1D* h_ue2_prior;
  TH1D* h_ue3_prior;
  TH1D* h_umu1_prior;
  TH1D* h_umu2_prior;
  TH1D* h_umu3_prior;
  TH1D* h_utau1_prior;
  TH1D* h_utau2_prior;
  TH1D* h_utau3_prior;

  TH1D* h_real_ue1;
  TH1D* h_real_ue2;
  TH1D* h_real_ue3;
  TH1D* h_real_umu1;
  TH1D* h_real_umu2;
  TH1D* h_real_umu3;
  TH1D* h_real_utau1;
  TH1D* h_real_utau2;
  TH1D* h_real_utau3;
  TH1D* h_real_ue1_prior;
  TH1D* h_real_ue2_prior;
  TH1D* h_real_ue3_prior;
  TH1D* h_real_umu1_prior;
  TH1D* h_real_umu2_prior;
  TH1D* h_real_umu3_prior;
  TH1D* h_real_utau1_prior;
  TH1D* h_real_utau2_prior;
  TH1D* h_real_utau3_prior;

  TH1D* h_imag_ue1;
  TH1D* h_imag_ue2;
  TH1D* h_imag_ue3;
  TH1D* h_imag_umu1;
  TH1D* h_imag_umu2;
  TH1D* h_imag_umu3;
  TH1D* h_imag_utau1;
  TH1D* h_imag_utau2;
  TH1D* h_imag_utau3;
  TH1D* h_imag_ue1_prior;
  TH1D* h_imag_ue2_prior;
  TH1D* h_imag_ue3_prior;
  TH1D* h_imag_umu1_prior;
  TH1D* h_imag_umu2_prior;
  TH1D* h_imag_umu3_prior;
  TH1D* h_imag_utau1_prior;
  TH1D* h_imag_utau2_prior;
  TH1D* h_imag_utau3_prior;

  TH1D* h_dcp;
  TH1D* h_sdcp;
  TH1D* h_cdcp;
  TH1D* h_s2th12;
  TH1D* h_th13;
  TH1D* h_s2th13;
  TH1D* h_cth13;
  TH1D* h_th23;
  TH1D* h_s2th23;
  TH1D* h_cth23;

  TH1D* h_dcp_prior;
  TH1D* h_sdcp_prior;
  TH1D* h_cdcp_prior;
  TH1D* h_th13_prior;
  TH1D* h_sth13_prior;
  TH1D* h_cth13_prior;
  TH1D* h_th23_prior;
  TH1D* h_sth23_prior;
  TH1D* h_cth23_prior;

  TH2D* h_ue1_s2th12;
  TH2D* h_ue2_s2th12;
  TH2D* h_ue3_s2th12;
  TH2D* h_umu1_s2th12;
  TH2D* h_umu2_s2th12;
  TH2D* h_umu3_s2th12;
  TH2D* h_utau1_s2th12;
  TH2D* h_utau2_s2th12;
  TH2D* h_utau3_s2th12;
  TH2D* h_ue1_s2th13;
  TH2D* h_ue2_s2th13;
  TH2D* h_ue3_s2th13;
  TH2D* h_umu1_s2th13;
  TH2D* h_umu2_s2th13;
  TH2D* h_umu3_s2th13;
  TH2D* h_utau1_s2th13;
  TH2D* h_utau2_s2th13;
  TH2D* h_utau3_s2th13;
  TH2D* h_ue1_s2th23;
  TH2D* h_ue2_s2th23;
  TH2D* h_ue3_s2th23;
  TH2D* h_umu1_s2th23;
  TH2D* h_umu2_s2th23;
  TH2D* h_umu3_s2th23;
  TH2D* h_utau1_s2th23;
  TH2D* h_utau2_s2th23;
  TH2D* h_utau3_s2th23;
  TH2D* h_ue1_dcp;
  TH2D* h_ue2_dcp;
  TH2D* h_ue3_dcp;
  TH2D* h_umu1_dcp;
  TH2D* h_umu2_dcp;
  TH2D* h_umu3_dcp;
  TH2D* h_utau1_dcp;
  TH2D* h_utau2_dcp;
  TH2D* h_utau3_dcp;  

  TH2D* h_ue1_ue2;
  TH2D* h_ue1_ue3;
  TH2D* h_ue1_umu1;
  TH2D* h_ue1_umu2;
  TH2D* h_ue1_umu3;
  TH2D* h_ue1_utau1;
  TH2D* h_ue1_utau2;
  TH2D* h_ue1_utau3;
  TH2D* h_ue2_ue3;
  TH2D* h_ue2_umu1;
  TH2D* h_ue2_umu2;
  TH2D* h_ue2_umu3;
  TH2D* h_ue2_utau1;
  TH2D* h_ue2_utau2;
  TH2D* h_ue2_utau3;
  TH2D* h_ue3_umu1;
  TH2D* h_ue3_umu2;
  TH2D* h_ue3_umu3;
  TH2D* h_ue3_utau1;
  TH2D* h_ue3_utau2;
  TH2D* h_ue3_utau3;
  TH2D* h_umu1_umu2;
  TH2D* h_umu1_umu3;
  TH2D* h_umu1_utau1;
  TH2D* h_umu1_utau2;
  TH2D* h_umu1_utau3;
  TH2D* h_umu2_umu3;
  TH2D* h_umu2_utau1;
  TH2D* h_umu2_utau2;
  TH2D* h_umu2_utau3;
  TH2D* h_umu3_utau1;
  TH2D* h_umu3_utau2;
  TH2D* h_umu3_utau3;
  TH2D* h_utau1_utau2;
  TH2D* h_utau1_utau3;
  TH2D* h_utau2_utau3;

  if(doModulus) {
    h_ue1 = new TH1D("h_ue1",";|U_{e1}|",n_bins,0.,1.);
    h_ue2 = new TH1D("h_ue2",";|U_{e2}|",n_bins,0.,1.);
    h_ue3 = new TH1D("h_ue3",";|U_{e3}|",n_bins,0.,1.);
    h_umu1 = new TH1D("h_umu1",";|U_{#mu1}|",n_bins,0.,1.);
    h_umu2 = new TH1D("h_umu2",";|U_{#mu2}|",n_bins,0.,1.);
    h_umu3 = new TH1D("h_umu3",";|U_{#mu3}|",n_bins,0.,1.);
    h_utau1 = new TH1D("h_utau1",";|U_{#tau1}|",n_bins,0.,1.);
    h_utau2 = new TH1D("h_utau2",";|U_{#tau2}|",n_bins,0.,1.);
    h_utau3 = new TH1D("h_utau3",";|U_{#tau3}|",n_bins,0.,1.);
  
    h_ue1_prior = new TH1D("h_ue1_prior",";|U_{e1}|",n_bins,0.,1.);
    h_ue2_prior = new TH1D("h_ue2_prior",";|U_{e2}|",n_bins,0.,1.);
    h_ue3_prior = new TH1D("h_ue3_prior",";|U_{e3}|",n_bins,0.,1.);
    h_umu1_prior = new TH1D("h_umu1_prior",";|U_{#mu1}|",n_bins,0.,1.);
    h_umu2_prior = new TH1D("h_umu2_prior",";|U_{#mu2}|",n_bins,0.,1.);
    h_umu3_prior = new TH1D("h_umu3_prior",";|U_{#mu3}|",n_bins,0.,1.);
    h_utau1_prior = new TH1D("h_utau1_prior",";|U_{#tau1}|",n_bins,0.,1.);
    h_utau2_prior = new TH1D("h_utau2_prior",";|U_{#tau2}|",n_bins,0.,1.);
    h_utau3_prior = new TH1D("h_utau3_prior",";|U_{#tau3}|",n_bins,0.,1.);
  }

  if(doReal) {
    h_real_ue1 = new TH1D("h_real_ue1",";Re(U_{e1})",n_bins,-1.,1.);
    h_real_ue2 = new TH1D("h_real_ue2",";Re(U_{e2})",n_bins,-1.,1.);
    h_real_ue3 = new TH1D("h_real_ue3",";Re(U_{e3})",n_bins,-1.,1.);
    h_real_umu1 = new TH1D("h_real_umu1",";Re(U_{#mu1})",n_bins,-1.,1.);
    h_real_umu2 = new TH1D("h_real_umu2",";Re(U_{#mu2})",n_bins,-1.,1.);
    h_real_umu3 = new TH1D("h_real_umu3",";Re(U_{#mu3})",n_bins,-1.,1.);
    h_real_utau1 = new TH1D("h_real_utau1",";Re(U_{#tau1})",n_bins,-1.,1.);
    h_real_utau2 = new TH1D("h_real_utau2",";Re(U_{#tau2})",n_bins,-1.,1.);
    h_real_utau3 = new TH1D("h_real_utau3",";Re(U_{#tau3})",n_bins,-1.,1.);

    h_real_ue1_prior = new TH1D("h_real_ue1_prior",";Re(U_{e1})",n_bins,-1.,1.);
    h_real_ue2_prior = new TH1D("h_real_ue2_prior",";Re(U_{e2})",n_bins,-1.,1.);
    h_real_ue3_prior = new TH1D("h_real_ue3_prior",";Re(U_{e3})",n_bins,-1.,1.);
    h_real_umu1_prior = new TH1D("h_real_umu1_prior",";Re(U_{#mu1})",n_bins,-1.,1.);
    h_real_umu2_prior = new TH1D("h_real_umu2_prior",";Re(U_{#mu2})",n_bins,-1.,1.);
    h_real_umu3_prior = new TH1D("h_real_umu3_prior",";Re(U_{#mu3})",n_bins,-1.,1.);
    h_real_utau1_prior = new TH1D("h_real_utau1_prior",";Re(U_{#tau1})",n_bins,-1.,1.);
    h_real_utau2_prior = new TH1D("h_real_utau2_prior",";Re(U_{#tau2})",n_bins,-1.,1.);
    h_real_utau3_prior = new TH1D("h_real_utau3_prior",";Re(U_{#tau3})",n_bins,-1.,1.);
  }

  if(doImag) {
    h_imag_ue1 = new TH1D("h_imag_ue1",";Im(U_{e1})",n_bins,-1.,1.);
    h_imag_ue2 = new TH1D("h_imag_ue2",";Im(U_{e2})",n_bins,-1.,1.);
    h_imag_ue3 = new TH1D("h_imag_ue3",";Im(U_{e3})",n_bins,-1.,1.);
    h_imag_umu1 = new TH1D("h_imag_umu1",";Im(U_{#mu1})",n_bins,-1.,1.);
    h_imag_umu2 = new TH1D("h_imag_umu2",";Im(U_{#mu2})",n_bins,-1.,1.);
    h_imag_umu3 = new TH1D("h_imag_umu3",";Im(U_{#mu3})",n_bins,-1.,1.);
    h_imag_utau1 = new TH1D("h_imag_utau1",";Im(U_{#tau1})",n_bins,-1.,1.);
    h_imag_utau2 = new TH1D("h_imag_utau2",";Im(U_{#tau2})",n_bins,-1.,1.);
    h_imag_utau3 = new TH1D("h_imag_utau3",";Im(U_{#tau3})",n_bins,-1.,1.);

    h_imag_ue1_prior = new TH1D("h_imag_ue1_prior",";Im(U_{e1})",n_bins,-1.,1.);
    h_imag_ue2_prior = new TH1D("h_imag_ue2_prior",";Im(U_{e2})",n_bins,-1.,1.);
    h_imag_ue3_prior = new TH1D("h_imag_ue3_prior",";Im(U_{e3})",n_bins,-1.,1.);
    h_imag_umu1_prior = new TH1D("h_imag_umu1_prior",";Im(U_{#mu1})",n_bins,-1.,1.);
    h_imag_umu2_prior = new TH1D("h_imag_umu2_prior",";Im(U_{#mu2})",n_bins,-1.,1.);
    h_imag_umu3_prior = new TH1D("h_imag_umu3_prior",";Im(U_{#mu3})",n_bins,-1.,1.);
    h_imag_utau1_prior = new TH1D("h_imag_utau1_prior",";Im(U_{#tau1})",n_bins,-1.,1.);
    h_imag_utau2_prior = new TH1D("h_imag_utau2_prior",";Im(U_{#tau2})",n_bins,-1.,1.);
    h_imag_utau3_prior = new TH1D("h_imag_utau3_prior",";Im(U_{#tau3})",n_bins,-1.,1.); 
  }
  
  if(do2D_pmns){
    h_ue1_s2th12 = new TH2D("h_ue1_s2th12",";|U_{e1}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_ue2_s2th12 = new TH2D("h_ue2_s2th12",";|U_{e2}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_ue3_s2th12 = new TH2D("h_ue3_s2th12",";|U_{e3}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_umu1_s2th12 = new TH2D("h_umu1_s2th12",";|U_{#mu1}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_umu2_s2th12 = new TH2D("h_umu2_s2th12",";|U_{#mu2}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_umu3_s2th12 = new TH2D("h_umu3_s2th12",";|U_{#mu3}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_utau1_s2th12 = new TH2D("h_utau1_s2th12",";|U_{#tau1}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_utau2_s2th12 = new TH2D("h_utau2_s2th12",";|U_{#tau2}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_utau3_s2th12 = new TH2D("h_utau3_s2th12",";|U_{#tau3}|;sin^{2}(#theta_{12})",n_bins,0.,1.,n_bins,0.15,0.5);
    h_ue1_s2th13 = new TH2D("h_ue1_s2th13",";|U_{e1}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_ue2_s2th13 = new TH2D("h_ue2_s2th13",";|U_{e2}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_ue3_s2th13 = new TH2D("h_ue3_s2th13",";|U_{e3}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_umu1_s2th13 = new TH2D("h_umu1_s2th13",";|U_{#mu1}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_umu2_s2th13 = new TH2D("h_umu2_s2th13",";|U_{#mu2}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_umu3_s2th13 = new TH2D("h_umu3_s2th13",";|U_{#mu3}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_utau1_s2th13 = new TH2D("h_utau1_s2th13",";|U_{#tau1}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_utau2_s2th13 = new TH2D("h_utau2_s2th13",";|U_{#tau2}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_utau3_s2th13 = new TH2D("h_utau3_s2th13",";|U_{#tau3}|;sin^{2}(#theta_{13})",n_bins,0.,1.,n_bins,0.,0.1);
    h_ue1_s2th23 = new TH2D("h_ue1_s2th23",";|U_{e1}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_ue2_s2th23 = new TH2D("h_ue2_s2th23",";|U_{e2}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_ue3_s2th23 = new TH2D("h_ue3_s2th23",";|U_{e3}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_umu1_s2th23 = new TH2D("h_umu1_s2th23",";|U_{#mu1}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_umu2_s2th23 = new TH2D("h_umu2_s2th23",";|U_{#mu2}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_umu3_s2th23 = new TH2D("h_umu3_s2th23",";|U_{#mu3}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_utau1_s2th23 = new TH2D("h_utau1_s2th23",";|U_{#tau1}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_utau2_s2th23 = new TH2D("h_utau2_s2th23",";|U_{#tau2}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_utau3_s2th23 = new TH2D("h_utau3_s2th23",";|U_{#tau3}|;sin^{2}(#theta_{23})",n_bins,0.,1.,n_bins,0.3,0.8);
    h_ue1_dcp = new TH2D("h_ue1_dcp",";|U_{e1}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_ue2_dcp = new TH2D("h_ue2_dcp",";|U_{e2}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_ue3_dcp = new TH2D("h_ue3_dcp",";|U_{e3}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_umu1_dcp = new TH2D("h_umu1_dcp",";|U_{#mu1}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_umu2_dcp = new TH2D("h_umu2_dcp",";|U_{#mu2}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_umu3_dcp = new TH2D("h_umu3_dcp",";|U_{#mu3}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_utau1_dcp = new TH2D("h_utau1_dcp",";|U_{#tau1}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_utau2_dcp = new TH2D("h_utau2_dcp",";|U_{#tau2}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
    h_utau3_dcp = new TH2D("h_utau3_dcp",";|U_{#tau3}|;#delta_{CP}",n_bins,0.,1.,n_bins,-TMath::Pi(),TMath::Pi());
  }

  h_dcp = new TH1D("h_dcp",";#delta_{CP}",n_bins,-TMath::Pi(),TMath::Pi());
  h_s2th12 = new TH1D("h_s2th12",";sin^{2}(#theta_{12})",n_bins,0.15,0.5);
  h_s2th13 = new TH1D("h_s2th13",";sin^{2}(#theta_{13})",n_bins,0.,0.1);
  h_s2th23 = new TH1D("h_s2th23",";sin^{2}(#theta_{23})",n_bins,0.3,0.8);

  // 2D histograms for all unique pairs of |U_{alpha i}| vs |U_{beta j}| (no self-pairs)
  TH2D* h_UU[36];
  std::string U_names[9] = {"ue1", "ue2", "ue3", "umu1", "umu2", "umu3", "utau1", "utau2", "utau3"};
  std::string U_tex[9]  = {"|U_{e1}|", "|U_{e2}|", "|U_{e3}|", "|U_{#mu1}|", "|U_{#mu2}|", "|U_{#mu3}|", "|U_{#tau1}|", "|U_{#tau2}|", "|U_{#tau3}|"};
  int idx = 0;
  if (do2D){
    for(int i=0; i<9; ++i){
      for(int j=i+1; j<9; ++j){
        std::stringstream name, title;
        name << "h_" << U_names[i] << "_" << U_names[j];
        title << ";" << U_tex[i] << ";" << U_tex[j];
        std::cout << "Creating histogram: " << name.str() << std::endl;
        h_UU[idx] = new TH2D(name.str().c_str(), title.str().c_str(), n_bins, 0., 1., n_bins, 0., 1.);
        ++idx;
      }
    }
    // h_ue1_ue2 = new TH2D("h_ue1_ue2",";|U_{e1}|;|U_{e2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_ue3 = new TH2D("h_ue1_ue3",";|U_{e1}|;|U_{e3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_umu1 = new TH2D("h_ue1_umu1",";|U_{e1}|;|U_{#mu1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_umu2 = new TH2D("h_ue1_umu2",";|U_{e1}|;|U_{#mu2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_umu3 = new TH2D("h_ue1_umu3",";|U_{e1}|;|U_{#mu3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_utau1 = new TH2D("h_ue1_utau1",";|U_{e1}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_utau2 = new TH2D("h_ue1_utau2",";|U_{e1}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue1_utau3 = new TH2D("h_ue1_utau3",";|U_{e1}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_ue3 = new TH2D("h_ue2_ue3",";|U_{e2}|;|U_{e3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_umu1 = new TH2D("h_ue2_umu1",";|U_{e2}|;|U_{#mu1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_umu2 = new TH2D("h_ue2_umu2",";|U_{e2}|;|U_{#mu2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_umu3 = new TH2D("h_ue2_umu3",";|U_{e2}|;|U_{#mu3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_utau1 = new TH2D("h_ue2_utau1",";|U_{e2}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_utau2 = new TH2D("h_ue2_utau2",";|U_{e2}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue2_utau3 = new TH2D("h_ue2_utau3",";|U_{e2}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_umu1 = new TH2D("h_ue3_umu1",";|U_{e3}|;|U_{#mu1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_umu2 = new TH2D("h_ue3_umu2",";|U_{e3}|;|U_{#mu2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_umu3 = new TH2D("h_ue3_umu3",";|U_{e3}|;|U_{#mu3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_utau1 = new TH2D("h_ue3_utau1",";|U_{e3}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_utau2 = new TH2D("h_ue3_utau2",";|U_{e3}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_ue3_utau3 = new TH2D("h_ue3_utau3",";|U_{e3}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu1_umu2 = new TH2D("h_umu1_umu2",";|U_{#mu1}|;|U_{#mu2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu1_umu3 = new TH2D("h_umu1_umu3",";|U_{#mu1}|;|U_{#mu3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu1_utau1 = new TH2D("h_umu1_utau1",";|U_{#mu1}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu1_utau2 = new TH2D("h_umu1_utau2",";|U_{#mu1}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu1_utau3 = new TH2D("h_umu1_utau3",";|U_{#mu1}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu2_umu3 = new TH2D("h_umu2_umu3",";|U_{#mu2}|;|U_{#mu3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu2_utau1 = new TH2D("h_umu2_utau1",";|U_{#mu2}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu2_utau2 = new TH2D("h_umu2_utau2",";|U_{#mu2}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu2_utau3 = new TH2D("h_umu2_utau3",";|U_{#mu2}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu3_utau1 = new TH2D("h_umu3_utau1",";|U_{#mu3}|;|U_{#tau1}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu3_utau2 = new TH2D("h_umu3_utau2",";|U_{#mu3}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_umu3_utau3 = new TH2D("h_umu3_utau3",";|U_{#mu3}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_utau1_utau2 = new TH2D("h_utau1_utau2",";|U_{#tau1}|;|U_{#tau2}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_utau1_utau3 = new TH2D("h_utau1_utau3",";|U_{#tau1}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
    // h_utau2_utau3 = new TH2D("h_utau2_utau3",";|U_{#tau2}|;|U_{#tau3}|",n_bins,0.,1.,n_bins,0.,1.);
  }

  // h_dcp_prior = new TH1D("h_dcp_prior",";#delta_{CP}",100,0.,2*TMath::Pi());
  // h_sdcp_prior = new TH1D("h_sdcp_prior",";sin(#delta_{CP})",100,-1.,1.);
  // h_cdcp_prior = new TH1D("h_cdcp_prior",";cos(#delta_{CP})",100,-1.,1.);

  // h_th13_prior = new TH1D("h_th13_prior",";#theta_{13}",100,0.,2*TMath::Pi());
  // h_sth13_prior = new TH1D("h_sth13_prior",";sin(#theta_{13})",100,-1.,1.);
  // h_cth13_prior = new TH1D("h_cth13_prior",";cos(#theta_{13})",100,-1.,1.);

  // h_th23_prior = new TH1D("h_th23_prior",";#theta_{23}",100,0,2*TMath::Pi());
  // h_sth23_prior = new TH1D("h_sth23_prior",";sin(#theta_{23})",100,-1.,1.);
  // h_cth23_prior = new TH1D("h_cth23_prior",";cos(#theta_{23})",100,-1.,1.);

  TRandom3* randGen = new TRandom3();

  // nSteps = 2000000;

  for(unsigned int i = 0;i<nSteps; i++) {
  // for(unsigned int i = 0;i<50000000; i++) {

    //std::cout << "step # " << i <<std::endl;
    //if(i%100000==0) std::cout << "step # " << i <<std::endl;
    if(i%1000000==0) std::cout << "step # " << i <<std::endl;

    // std::cout<<"Try getting entry "<<i<<std::endl;
    post->GetEntry(i);

    // if (fitter == "Aria") {
    //   int step = i;
    // }

    // burn-in cut
    if (incl_burnin){
      if(step<n_steps_burnin) continue;
    }

    if (fitter == "Aria") {
      s13 = TMath::Sin(th13);
      s23 = TMath::Sin(th23);
      s12 = TMath::Sqrt(randGen->Gaus(0.307,0.013));
      sdcp = TMath::Sin(dcp*TMath::Pi());
      cdcp = TMath::Cos(dcp*TMath::Pi());
      c13 = TMath::Cos(th13);
      c12 = TMath::Sqrt(1.-TMath::Sq(s12));
      c23 = TMath::Cos(th23);
      rc = 1.;

      // double prior_rw_dCP[3] = {1.,cdcp,-sdcp};
      // double prior_rw_th13[3] = {1.,c13,-s13};
      // double prior_rw_th23[3] = {1.,c23,-s23};
      //std::cout<<s22th13<<"   "<<s13<<"     "<<c13<<std::endl;
    }
    else {
      if (!solarConstraint){
        s2th12 = randGen->Uniform(0.,1.);
      }
      s13 = TMath::Sqrt(s2th13);
      s23 = TMath::Sqrt(s2th23);
      if (chain != "DUNE"){
        s12 = TMath::Sqrt(s2th12);
      }
      else{ 
        s2th12 = randGen->Gaus(0.307,0.013);
        s12 = TMath::Sqrt(s2th12);
      }
      sdcp = TMath::Sin(dcp);
      cdcp = TMath::Cos(dcp);
      c13 = TMath::Sqrt(1.-s2th13);
      c12 = TMath::Sqrt(1.-s2th12);
      c23 = TMath::Sqrt(1.-s2th23);

      // double prior_rw_dCP[3] = {1.,cdcp,-sdcp};
      // double prior_rw_th13[3] = {1.,c13,-s13};
      // double prior_rw_th23[3] = {1.,c23,-s23};

      if (RCreweight == false){
        rc = 1.;
      }
    }
    //reweight
    //rw = prior_rw_dCP[prior_map[prior_case][0]]*prior_rw_th13[prior_map[prior_case][1]]*prior_rw_th23[prior_map[prior_case][2]]*rc;
    if (RCreweight){
      rw = rc;
    }
    else if (RCreweight && T2Kprior){
      rw = rc;
    }
    else if (haarPrior){
      rw = 2*(1-s2th13);
    }
    else if (NOvAprior){
      //rw = TMath::Abs(8*s12*c12*s13*c13*s23*c23);
      rw = TMath::Abs(4*s13*c13*s23*c23); //not transforming theta12
    }
    else{
      rw = 1.;
    }
    
    if (sinDelCPprior){
      rw = rw * TMath::Abs(cdcp);
    }

    //std::cout<<std::endl;
    //std::cout<<prior_rw_dCP[prior_map[prior_case][0]]<<std::endl;
    //std::cout<<prior_rw_dCP<<std::endl;
    //std::cout<<prior_map[prior_case]<<std::endl;
    //std::cout<<prior_map[prior_case][0]<<std::endl;
    //std::cout<<prior_rw_th13[prior_map[prior_case][1]]<<std::endl;
    //std::cout<<prior_rw_th23[prior_map[prior_case][2]]<<std::endl;
    //std::cout<<rc<<std::endl;
    //std::cout<<rw<<std::endl;
    //std::cout<<std::endl;

    real_ue1 = c12*c13;
    imag_ue1 = 0.;
    real_ue2 = s12*c13;
    imag_ue2 = 0.;
    real_ue3 = s13*cdcp;
    imag_ue3 = -s13*sdcp;

    real_umu1 = -s12*c23-c12*s23*s13*cdcp;
    imag_umu1 = -c12*s23*s13*sdcp;
    real_umu2 = c12*c23-s12*s23*s13*cdcp;
    imag_umu2 = -s12*s23*s13*sdcp;
    real_umu3 = s23*c13;
    imag_umu3 = 0.;

    real_utau1 = s12*s23-c12*c23*s13*cdcp;
    imag_utau1 = -c12*c23*s13*sdcp;
    real_utau2 = -c12*s23-s12*c23*s13*cdcp;
    imag_utau2 = -s12*c23*s13*sdcp;
    real_utau3 = c23*c13;
    imag_utau3 = 0.;

    ue1 = sqrt(real_ue1*real_ue1+imag_ue1*imag_ue1);
    ue2 = sqrt(real_ue2*real_ue2+imag_ue2*imag_ue2);
    ue3 = sqrt(real_ue3*real_ue3+imag_ue3*imag_ue3);

    umu1 = sqrt(real_umu1*real_umu1+imag_umu1*imag_umu1);
    umu2 = sqrt(real_umu2*real_umu2+imag_umu2*imag_umu2);
    umu3 = sqrt(real_umu3*real_umu3+imag_umu3*imag_umu3);

    utau1 = sqrt(real_utau1*real_utau1+imag_utau1*imag_utau1);
    utau2 = sqrt(real_utau2*real_utau2+imag_utau2*imag_utau2);
    utau3 = sqrt(real_utau3*real_utau3+imag_utau3*imag_utau3);


    //pmnsTree->Fill();

    h_dcp->Fill(dcp,rw);
    h_s2th12->Fill(s2th12,rw);
    h_s2th13->Fill(s2th13,rw);
    h_s2th23->Fill(s2th23,rw);

    if(doModulus) {
      //MP: what is rc, i.e. RCweight? --> reactor constraint used for reweighting
      h_ue1->Fill(ue1,rw);
      h_ue2->Fill(ue2,rw);
      h_ue3->Fill(ue3,rw);
      h_umu1->Fill(umu1,rw);
      h_umu2->Fill(umu2,rw);
      h_umu3->Fill(umu3,rw);
      h_utau1->Fill(utau1,rw);
      h_utau2->Fill(utau2,rw);
      h_utau3->Fill(utau3,rw);
    }

    if(doReal) {
      h_real_ue1->Fill(real_ue1,rw);
      h_real_ue2->Fill(real_ue2,rw);
      h_real_ue3->Fill(real_ue3,rw);
      h_real_umu1->Fill(real_umu1,rw);
      h_real_umu2->Fill(real_umu2,rw);
      h_real_umu3->Fill(real_umu3,rw);
      h_real_utau1->Fill(real_utau1,rw);
      h_real_utau2->Fill(real_utau2,rw);
      h_real_utau3->Fill(real_utau3,rw);
    }


    if(doImag) {
      h_imag_ue1->Fill(imag_ue1,rw);
      h_imag_ue2->Fill(imag_ue2,rw);
      h_imag_ue3->Fill(imag_ue3,rw);
      h_imag_umu1->Fill(imag_umu1,rw);
      h_imag_umu2->Fill(imag_umu2,rw);
      h_imag_umu3->Fill(imag_umu3,rw);
      h_imag_utau1->Fill(imag_utau1,rw);
      h_imag_utau2->Fill(imag_utau2,rw);
      h_imag_utau3->Fill(imag_utau3,rw);
    }

    // std::cout<<"Filling 2D histograms"<<std::endl;

    if(do2D_pmns){
      h_ue1_s2th12->Fill(ue1,s2th12,rw);
      h_ue2_s2th12->Fill(ue2,s2th12,rw);
      h_ue3_s2th12->Fill(ue3,s2th12,rw);
      h_umu1_s2th12->Fill(umu1,s2th12,rw);
      h_umu2_s2th12->Fill(umu2,s2th12,rw);
      h_umu3_s2th12->Fill(umu3,s2th12,rw);
      h_utau1_s2th12->Fill(utau1,s2th12,rw);
      h_utau2_s2th12->Fill(utau2,s2th12,rw);
      h_utau3_s2th12->Fill(utau3,s2th12,rw);
      h_ue1_s2th13->Fill(ue1,s2th13,rw);
      h_ue2_s2th13->Fill(ue2,s2th13,rw);
      h_ue3_s2th13->Fill(ue3,s2th13,rw);
      h_umu1_s2th13->Fill(umu1,s2th13,rw);
      h_umu2_s2th13->Fill(umu2,s2th13,rw);
      h_umu3_s2th13->Fill(umu3,s2th13,rw);
      h_utau1_s2th13->Fill(utau1,s2th13,rw);
      h_utau2_s2th13->Fill(utau2,s2th13,rw);
      h_utau3_s2th13->Fill(utau3,s2th13,rw);
      h_ue1_s2th23->Fill(ue1,s2th23,rw);
      h_ue2_s2th23->Fill(ue2,s2th23,rw);
      h_ue3_s2th23->Fill(ue3,s2th23,rw);
      h_umu1_s2th23->Fill(umu1,s2th23,rw);
      h_umu2_s2th23->Fill(umu2,s2th23,rw);
      h_umu3_s2th23->Fill(umu3,s2th23,rw);
      h_utau1_s2th23->Fill(utau1,s2th23,rw);
      h_utau2_s2th23->Fill(utau2,s2th23,rw);
      h_utau3_s2th23->Fill(utau3,s2th23,rw);
      h_ue1_dcp->Fill(ue1,dcp,rw);
      h_ue2_dcp->Fill(ue2,dcp,rw);
      h_ue3_dcp->Fill(ue3,dcp,rw);
      h_umu1_dcp->Fill(umu1,dcp,rw);
      h_umu2_dcp->Fill(umu2,dcp,rw);
      h_umu3_dcp->Fill(umu3,dcp,rw);
      h_utau1_dcp->Fill(utau1,dcp,rw);
      h_utau2_dcp->Fill(utau2,dcp,rw);
      h_utau3_dcp->Fill(utau3,dcp,rw);
    }

    // std::cout << "Filling 2D histograms done" << std::endl;

    if(do2D){
      // Store in an array for easy access
      double U_mod[9] = {ue1, ue2, ue3, umu1, umu2, umu3, utau1, utau2, utau3};
      idx = 0;
      for(int i=0; i<9; ++i){
        for(int j=i+1; j<9; ++j){
          h_UU[idx]->Fill(U_mod[i], U_mod[j], rw);
          ++idx;
        }
      }
    }

  }

  if(doPrior) {
    //for(unsigned int s=0;s<50000000;s++){
    for(unsigned int s=0;s<50000000;s++){
      if(s%500000==0) std::cout << "throw # " << s <<std::endl;
      //sample from prior:

      //MP: why are priors hardcoded?? --> that's how it is done for now but can be made more configurable later

      if (fitter == "Aria") {
        prior_th23 = randGen->Uniform(0.,2*TMath::Pi());
        prior_s2th12 = randGen->Gaus(0.307,0.013);
        prior_dcp = randGen->Uniform(0.,2*TMath::Pi());
        prior_sindcp = randGen->Uniform(-1.,1.);

        if (RCreweight){
          prior_th13 = randGen->Gaus(0.1476,0.0264);
          s13_prior = TMath::Sin(prior_th13);
          c13_prior = TMath::Cos(prior_th13);
        }
        else{
          prior_th13 = randGen->Uniform(0.,2*TMath::Pi());
          s13_prior = TMath::Sin(prior_th13);
          c13_prior = TMath::Cos(prior_th13);
        }
        s23_prior = TMath::Sin(prior_th23);
        s12_prior = TMath::Sqrt(prior_s2th12);
        sdcp_prior = TMath::Sin(prior_dcp);
        cdcp_prior = TMath::Cos(prior_dcp);
        c12_prior = TMath::Sqrt(1.-prior_s2th12);
        c23_prior = TMath::Cos(prior_th23);

        // double prior_rw_dCP_prior[3] = {1.,cdcp_prior,-sdcp_prior};
        // double prior_rw_th13_prior[3] = {1.,c13_prior,-s13_prior};
        // //prior_rw_th13_wRC_prior = {1.,c13_wRC_prior,-s13_wRC_prior};
        // double prior_rw_th23_prior[3] = {1.,c23_prior,-s23_prior};
      }
      else {
        if (!sinDelCPprior){
          prior_dcp = randGen->Uniform(-1.*TMath::Pi(),TMath::Pi());
          sdcp_prior = TMath::Sin(prior_dcp);
        }
        else{
          sdcp_prior = randGen->Uniform(-1.,1.);
          prior_dcp = TMath::ASin(sdcp_prior);
        }
        cdcp_prior = TMath::Cos(prior_dcp);

        if (!NOvAprior){
          prior_s2th23 = randGen->Uniform(0.,1.);

          if (RCreweight){
            prior_s2th13 = randGen->Gaus(0.0218,0.0007);
            s13_prior = TMath::Sqrt(prior_s2th13);
            c13_prior = TMath::Sqrt(1.-prior_s2th13);
          }
          else if (haarPrior){
            prior_c4th13 = randGen->Uniform(0.,1.);
            prior_c2th13 = TMath::Sqrt(prior_c4th13);
            prior_s2th13 = 1-prior_c2th13;
            c13_prior = TMath::Sqrt(prior_c2th13);
            s13_prior = TMath::Sqrt(prior_s2th13);
          }
          else{
            prior_s2th13 = randGen->Uniform(0.,1.);
            s13_prior = TMath::Sqrt(prior_s2th13);
            c13_prior = TMath::Sqrt(1.-prior_s2th13);
          }

          if (solarConstraint){
            prior_s2th12 = randGen->Gaus(0.307,0.013);
          }
          else{
            prior_s2th12 = randGen->Uniform(0.,1.);
          }

          s23_prior = TMath::Sqrt(prior_s2th23);
          s12_prior = TMath::Sqrt(prior_s2th12);

          c12_prior = TMath::Sqrt(1.-prior_s2th12);
          c23_prior = TMath::Sqrt(1.-prior_s2th23);

        }
        else{
          prior_th23 = randGen->Uniform(0.,2*TMath::Pi());

          s23_prior = TMath::Sin(prior_th23);
          c23_prior = TMath::Cos(prior_th23);

          if (RCreweight){
            prior_th13 = randGen->Gaus(0.1476,0.0264);
          }
          else{
            prior_th13 = randGen->Uniform(0.,2*TMath::Pi());
          }

          s13_prior = TMath::Sin(prior_th13);
          c13_prior = TMath::Cos(prior_th13);

          if (solarConstraint){
            prior_s2th12 = randGen->Gaus(0.307,0.013);
            s12_prior = TMath::Sqrt(prior_s2th12);
            c12_prior = TMath::Sqrt(1.-prior_s2th12);
          }
          else{
            prior_th12 = randGen->Uniform(0.,2*TMath::Pi());
            s12_prior = TMath::Sin(prior_th12);
            c12_prior = TMath::Cos(prior_th12);
          }

        }
        
      }
      //reweighting for priors doesn't include reactor constraint so far
      //rw_prior = prior_rw_dCP_prior[prior_map[prior_case][0]]*prior_rw_th13_prior[prior_map[prior_case][1]]*prior_rw_th23_prior[prior_map[prior_case][2]];
  
      real_ue1_prior = c12_prior*c13_prior;
      imag_ue1_prior = 0.;
      real_ue2_prior = s12_prior*c13_prior;
      imag_ue2_prior = 0.;
      real_ue3_prior = s13_prior*cdcp_prior;
      imag_ue3_prior = -s13_prior*sdcp_prior;
  
      real_umu1_prior = -s12_prior*c23_prior-c12_prior*s23_prior*s13_prior*cdcp_prior;
      imag_umu1_prior = -c12_prior*s23_prior*s13_prior*sdcp_prior;
      real_umu2_prior = c12_prior*c23_prior-s12_prior*s23_prior*s13_prior*cdcp_prior;
      imag_umu2_prior = -s12_prior*s23_prior*s13_prior*sdcp_prior;
      real_umu3_prior = s23_prior*c13_prior;
      imag_umu3_prior = 0.;
  
      real_utau1_prior = s12_prior*s23_prior-c12_prior*c23_prior*s13_prior*cdcp_prior;
      imag_utau1_prior = -c12_prior*c23_prior*s13_prior*sdcp_prior;
      real_utau2_prior = -c12_prior*s23_prior-s12_prior*c23_prior*s13_prior*cdcp_prior;
      imag_utau2_prior = -s12_prior*c23_prior*s13_prior*sdcp_prior;
      real_utau3_prior = c23_prior*c13_prior;
      imag_utau3_prior = 0.;
  
      ue1_prior = sqrt(real_ue1_prior*real_ue1_prior+imag_ue1_prior*imag_ue1_prior);
      ue2_prior = sqrt(real_ue2_prior*real_ue2_prior+imag_ue2_prior*imag_ue2_prior);
      ue3_prior = sqrt(real_ue3_prior*real_ue3_prior+imag_ue3_prior*imag_ue3_prior);
  
      umu1_prior = sqrt(real_umu1_prior*real_umu1_prior+imag_umu1_prior*imag_umu1_prior);
      umu2_prior = sqrt(real_umu2_prior*real_umu2_prior+imag_umu2_prior*imag_umu2_prior);
      umu3_prior = sqrt(real_umu3_prior*real_umu3_prior+imag_umu3_prior*imag_umu3_prior);
  
      utau1_prior = sqrt(real_utau1_prior*real_utau1_prior+imag_utau1_prior*imag_utau1_prior);
      utau2_prior = sqrt(real_utau2_prior*real_utau2_prior+imag_utau2_prior*imag_utau2_prior);
      utau3_prior = sqrt(real_utau3_prior*real_utau3_prior+imag_utau3_prior*imag_utau3_prior);

      //rw_prior = c13_prior;//cdcp_prior*c13_prior;

      if(doModulus) {
        h_ue1_prior->Fill(ue1_prior);//,rw_prior);
        h_ue2_prior->Fill(ue2_prior);//,rw_prior);
        h_ue3_prior->Fill(ue3_prior);//,rw_prior);
        h_umu1_prior->Fill(umu1_prior);//,rw_prior);
        h_umu2_prior->Fill(umu2_prior);//,rw_prior);
        h_umu3_prior->Fill(umu3_prior);//,rw_prior);
        h_utau1_prior->Fill(utau1_prior);//,rw_prior);
        h_utau2_prior->Fill(utau2_prior);//,rw_prior);
        h_utau3_prior->Fill(utau3_prior);//,rw_prior);
      }
      if(doReal) {
        h_real_ue1_prior->Fill(real_ue1_prior);//,rw_prior);
        h_real_ue2_prior->Fill(real_ue2_prior);//,rw_prior);
        h_real_ue3_prior->Fill(real_ue3_prior);//,rw_prior);
        h_real_umu1_prior->Fill(real_umu1_prior);//,rw_prior);
        h_real_umu2_prior->Fill(real_umu2_prior);//,rw_prior);
        h_real_umu3_prior->Fill(real_umu3_prior);//,rw_prior);
        h_real_utau1_prior->Fill(real_utau1_prior);//,rw_prior);
        h_real_utau2_prior->Fill(real_utau2_prior);//,rw_prior);
        h_real_utau3_prior->Fill(real_utau3_prior);//,rw_prior);
      }
      if(doImag) {
        h_imag_ue1_prior->Fill(imag_ue1_prior);//,rw_prior);
        h_imag_ue2_prior->Fill(imag_ue2_prior);//,rw_prior);
        h_imag_ue3_prior->Fill(imag_ue3_prior);//,rw_prior);
        h_imag_umu1_prior->Fill(imag_umu1_prior);//,rw_prior);
        h_imag_umu2_prior->Fill(imag_umu2_prior);//,rw_prior);
        h_imag_umu3_prior->Fill(imag_umu3_prior);//,rw_prior);
        h_imag_utau1_prior->Fill(imag_utau1_prior);//,rw_prior);
        h_imag_utau2_prior->Fill(imag_utau2_prior);//,rw_prior);
        h_imag_utau3_prior->Fill(imag_utau3_prior);//,rw_prior);
      }

      // h_dcp_prior->Fill(prior_dcp,rw_prior);
      // h_sdcp_prior->Fill(sdcp_prior,rw_prior);
      // h_cdcp_prior->Fill(cdcp_prior,rw_prior);
      // h_th13_prior->Fill(prior_th13,rw_prior);
      // h_sth13_prior->Fill(s13_prior,rw_prior);
      // h_cth13_prior->Fill(c13_prior,rw_prior);
      // h_th23_prior->Fill(prior_th23,rw_prior);
      // h_sth23_prior->Fill(s23_prior,rw_prior);
      // h_cth23_prior->Fill(c23_prior,rw_prior);
    }
  }

  outFile->cd();
  //pmnsTree->Write();

  if(doModulus) {
    h_ue1->Write("h_ue1");
    h_ue2->Write("h_ue2");
    h_ue3->Write("h_ue3");
    h_umu1->Write("h_umu1");
    h_umu2->Write("h_umu2");
    h_umu3->Write("h_umu3");
    h_utau1->Write("h_utau1");
    h_utau2->Write("h_utau2");
    h_utau3->Write("h_utau3");

    h_ue1_prior->Write("h_ue1_prior");
    h_ue2_prior->Write("h_ue2_prior");
    h_ue3_prior->Write("h_ue3_prior");
    h_umu1_prior->Write("h_umu1_prior");
    h_umu2_prior->Write("h_umu2_prior");
    h_umu3_prior->Write("h_umu3_prior");
    h_utau1_prior->Write("h_utau1_prior");
    h_utau2_prior->Write("h_utau2_prior");
    h_utau3_prior->Write("h_utau3_prior");
  }

  if(doReal) {
    h_real_ue1->Write("h_real_ue1");
    h_real_ue2->Write("h_real_ue2");
    h_real_ue3->Write("h_real_ue3");
    h_real_umu1->Write("h_real_umu1");
    h_real_umu2->Write("h_real_umu2");
    h_real_umu3->Write("h_real_umu3");
    h_real_utau1->Write("h_real_utau1");
    h_real_utau2->Write("h_real_utau2");
    h_real_utau3->Write("h_real_utau3");

    h_real_ue1_prior->Write("h_real_ue1_prior");
    h_real_ue2_prior->Write("h_real_ue2_prior");
    h_real_ue3_prior->Write("h_real_ue3_prior");
    h_real_umu1_prior->Write("h_real_umu1_prior");
    h_real_umu2_prior->Write("h_real_umu2_prior");
    h_real_umu3_prior->Write("h_real_umu3_prior");
    h_real_utau1_prior->Write("h_real_utau1_prior");
    h_real_utau2_prior->Write("h_real_utau2_prior");
    h_real_utau3_prior->Write("h_real_utau3_prior");
  }

  if(doImag) {
    h_imag_ue1->Write("h_imag_ue1");
    h_imag_ue2->Write("h_imag_ue2");
    h_imag_ue3->Write("h_imag_ue3");
    h_imag_umu1->Write("h_imag_umu1");
    h_imag_umu2->Write("h_imag_umu2");
    h_imag_umu3->Write("h_imag_umu3");
    h_imag_utau1->Write("h_imag_utau1");
    h_imag_utau2->Write("h_imag_utau2");
    h_imag_utau3->Write("h_imag_utau3");

    h_imag_ue1_prior->Write("h_imag_ue1_prior");
    h_imag_ue2_prior->Write("h_imag_ue2_prior");
    h_imag_ue3_prior->Write("h_imag_ue3_prior");
    h_imag_umu1_prior->Write("h_imag_umu1_prior");
    h_imag_umu2_prior->Write("h_imag_umu2_prior");
    h_imag_umu3_prior->Write("h_imag_umu3_prior");
    h_imag_utau1_prior->Write("h_imag_utau1_prior");
    h_imag_utau2_prior->Write("h_imag_utau2_prior");
    h_imag_utau3_prior->Write("h_imag_utau3_prior");
  }

  if(do2D_pmns){
    h_ue1_s2th12->Write("h_ue1_s2th12");
    h_ue2_s2th12->Write("h_ue2_s2th12");
    h_ue3_s2th12->Write("h_ue3_s2th12");
    h_umu1_s2th12->Write("h_umu1_s2th12");
    h_umu2_s2th12->Write("h_umu2_s2th12");
    h_umu3_s2th12->Write("h_umu3_s2th12");
    h_utau1_s2th12->Write("h_utau1_s2th12");
    h_utau2_s2th12->Write("h_utau2_s2th12");
    h_utau3_s2th12->Write("h_utau3_s2th12");
    h_ue1_s2th13->Write("h_ue1_s2th13");
    h_ue2_s2th13->Write("h_ue2_s2th13");
    h_ue3_s2th13->Write("h_ue3_s2th13");
    h_umu1_s2th13->Write("h_umu1_s2th13");
    h_umu2_s2th13->Write("h_umu2_s2th13");
    h_umu3_s2th13->Write("h_umu3_s2th13");
    h_utau1_s2th13->Write("h_utau1_s2th13");
    h_utau2_s2th13->Write("h_utau2_s2th13");
    h_utau3_s2th13->Write("h_utau3_s2th13");
    h_ue1_s2th23->Write("h_ue1_s2th23");
    h_ue2_s2th23->Write("h_ue2_s2th23");
    h_ue3_s2th23->Write("h_ue3_s2th23");
    h_umu1_s2th23->Write("h_umu1_s2th23");
    h_umu2_s2th23->Write("h_umu2_s2th23");
    h_umu3_s2th23->Write("h_umu3_s2th23");
    h_utau1_s2th23->Write("h_utau1_s2th23");
    h_utau2_s2th23->Write("h_utau2_s2th23");
    h_utau3_s2th23->Write("h_utau3_s2th23");
    h_ue1_dcp->Write("h_ue1_dcp");
    h_ue2_dcp->Write("h_ue2_dcp");
    h_ue3_dcp->Write("h_ue3_dcp");
    h_umu1_dcp->Write("h_umu1_dcp");
    h_umu2_dcp->Write("h_umu2_dcp");
    h_umu3_dcp->Write("h_umu3_dcp");
    h_utau1_dcp->Write("h_utau1_dcp");
    h_utau2_dcp->Write("h_utau2_dcp");
    h_utau3_dcp->Write("h_utau3_dcp");
  }

  if(do2D){
    // Write all |U_{alpha i}| vs |U_{beta j}| 2D histograms (no self-pairs)
    for(int idx=0; idx<36; ++idx){
      h_UU[idx]->Write();
    }
  }

  h_dcp->Write("h_dcp");
  h_s2th12->Write("h_s2th12");
  h_s2th13->Write("h_s2th13");
  h_s2th23->Write("h_s2th23");

  outFile->Close();

  // h_dcp_prior->Write("h_dcp_prior");
  // h_sdcp_prior->Write("h_sdcp_prior");
  // h_cdcp_prior->Write("h_cdcp_prior");
  // h_th13_prior->Write("h_th13_prior");
  // h_sth13_prior->Write("h_sth13_prior");
  // h_cth13_prior->Write("h_cth13_prior");
  // h_th23_prior->Write("h_th23_prior");
  // h_sth23_prior->Write("h_sth23_prior");
  // h_cth23_prior->Write("h_cth23_prior");
}
