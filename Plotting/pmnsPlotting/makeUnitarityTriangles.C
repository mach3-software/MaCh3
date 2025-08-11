#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TComplex.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLine.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <sstream>

void makeUnitarityTriangles(bool RCreweight = false, std::string fit_type = "FD+ND"){

  std::string chain = "DUNE";
  std::string fitter = "MaCh3";

  int n_bins = 200;

  unsigned int n_chains = 70;

  // Does the chain include burn-in?
  bool incl_burnin;
  if (fit_type != "OA2020"){
    incl_burnin = true;
    std::cout<<"Cut burnin..."<<std::endl;
  }
  else{
    incl_burnin = false;
  }
  unsigned int n_steps_burnin;  // 150000; T2K , 80000 DUNE
  if (chain != "DUNE") {
    n_steps_burnin = 150000;
  }
  else {
    n_steps_burnin = 80000;
  }

  bool doModulus = false;
  bool doReal = false;
  bool doImag = false;
  bool doPrior = false;

  bool doJarlskog = true;

  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};
  std::map<std::string, std::string> infile_str{{"AsimovA22","oa2021/AsimovA22_NonUniform_041022_reduced_smeared_reweighted.root"},{"AsimovB22","oa2021/AsimovB22_NonUniform_Full_161222_reduced_smeared_reweighted.root"},{"OA2021","oa2021/Data_NonUniform_021222_Full_reduced_reweighted_smeared.root"},{"OA2020","MaCh3-OA2020_ALL_data_RCweights_reduced_burnInCut_smearedRevision_withPionSI.root"}};
  std::map<std::string, std::string> infile_str_dune{{"FDOnly","dune/DUNE_FDOnly_reduced.root"},{"FD+ND","dune/DUNE_NuPhys_wNDDet_reduced_all_reweighted.root"}};

  TFile* outFile;
  std::stringstream filename;

  // Print out the date and time of the run
  time_t now = time(0);
  char* dt = ctime(&now);
  std::cout << "Run date and time: " << dt << std::endl;

  if (chain == "NOvA") {
    std::cout<<"Chain: T2K+NOvA"<<std::endl;
    filename<<"unitarity_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors"<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  else if (chain == "T2K" && fitter == "Aria"){
    std::cout<<"Chain: Aria T2K-only"<<std::endl;
    filename<<"unitarity_T2Konly_Aria_Asimov1_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors"<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  else {
    std::cout<<"Chain: "<<chain<<";  Fit type: "<<fit_type<<std::endl;
    filename<<"unitarity_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_"<<bool_str[doPrior]<<"Priors"<<"_"<<n_bins<<"bins.root";
    outFile = new TFile(filename.str().c_str(),"RECREATE");
  }
  std::cout<<"Outfile: "<<filename.str().c_str()<<std::endl;

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

  std::cout<<"Doing Jarlskog? "<<doJarlskog<<std::endl;

  double th13, th23;
  double s2th13, s2th23, s2th12, dcp, dm2, rc;
  double s22th13, s22th23, prob, mh;
  double prior_s2th13, prior_wRC_s2th13, prior_s2th23, prior_s2th12, prior_dcp, prior_sindcp;
  double prior_th13, prior_wRC_th13, prior_th23;
  int step;

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
  else {
    post->SetBranchAddress("theta13",&s2th13);
    post->SetBranchAddress("theta23",&s2th23);
    post->SetBranchAddress("dcp",&dcp);
    post->SetBranchAddress("dm23",&dm2);
    post->SetBranchAddress("step",&step);
    post->SetBranchAddress("theta12",&s2th12);
    if (fit_type != "FDOnly"){
      post->SetBranchAddress("RCreweight",&rc);
    }
  }

  

  double s13, s23, s12, sdcp, cdcp, c13, c23, c12, j;
  double real_ue1, real_ue2, real_ue3, real_umu1, real_umu2, real_umu3, real_utau1, real_utau2, real_utau3;
  double imag_ue1, imag_ue2, imag_ue3, imag_umu1, imag_umu2, imag_umu3, imag_utau1, imag_utau2, imag_utau3;
  TComplex ue1, ue2, ue3, umu1, umu2, umu3, utau1, utau2, utau3;

  TComplex tr_emu_num, tr_etau_num, tr_mutau_num, tr_12_num, tr_13_num, tr_23_num;
  TComplex tr_emu_denom, tr_etau_denom, tr_mutau_denom, tr_12_denom, tr_13_denom, tr_23_denom;
  TComplex tr_emu, tr_etau, tr_mutau, tr_12, tr_13, tr_23;

  double s13_prior, s13_wRC_prior, s23_prior, s12_prior, sdcp_prior, cdcp_prior, c13_prior, c13_wRC_prior, c23_prior, c12_prior;
  //double ue1_prior, ue2_prior, ue3_prior, umu1_prior, umu2_prior, umu3_prior, utau1_prior, utau2_prior, utau3_prior;
  double real_ue1_prior, real_ue2_prior, real_ue3_prior, real_umu1_prior, real_umu2_prior, real_umu3_prior, real_utau1_prior, real_utau2_prior, real_utau3_prior;
  double imag_ue1_prior, imag_ue2_prior, imag_ue3_prior, imag_umu1_prior, imag_umu2_prior, imag_umu3_prior, imag_utau1_prior, imag_utau2_prior, imag_utau3_prior;
  TComplex ue1_prior, ue2_prior, ue3_prior, umu1_prior, umu2_prior, umu3_prior, utau1_prior, utau2_prior, utau3_prior;

  TComplex tr_emu_num_prior, tr_etau_num_prior, tr_mutau_num_prior, tr_12_num_prior, tr_13_num_prior, tr_23_num_prior;
  TComplex tr_emu_denom_prior, tr_etau_denom_prior, tr_mutau_denom_prior, tr_12_denom_prior, tr_13_denom_prior, tr_23_denom_prior;
  TComplex tr_emu_prior, tr_etau_prior, tr_mutau_prior, tr_12_prior, tr_13_prior, tr_23_prior;

  double j_e1, j_e2, j_mu1, j_mu3, j_tau2, j_tau3;

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

  TH2D* h_tr_emu;
  TH2D* h_tr_etau;
  TH2D* h_tr_mutau;
  TH2D* h_tr_12;
  TH2D* h_tr_13;
  TH2D* h_tr_23;
  TH2D* h_tr_emu_prior;
  TH2D* h_tr_etau_prior;
  TH2D* h_tr_mutau_prior;
  TH2D* h_tr_12_prior;
  TH2D* h_tr_13_prior;
  TH2D* h_tr_23_prior;

  TH1D* h_rho_emu;
  TH1D* h_rho_etau;
  TH1D* h_rho_mutau;
  TH1D* h_rho_12;
  TH1D* h_rho_13;
  TH1D* h_rho_23;
  TH1D* h_rho_emu_prior;
  TH1D* h_rho_etau_prior;
  TH1D* h_rho_mutau_prior;
  TH1D* h_rho_12_prior;
  TH1D* h_rho_13_prior;
  TH1D* h_rho_23_prior;

  TH1D* h_eta_emu;
  TH1D* h_eta_etau;
  TH1D* h_eta_mutau;
  TH1D* h_eta_12;
  TH1D* h_eta_13;
  TH1D* h_eta_23;
  TH1D* h_eta_emu_prior;
  TH1D* h_eta_etau_prior;
  TH1D* h_eta_mutau_prior;
  TH1D* h_eta_12_prior;
  TH1D* h_eta_13_prior;
  TH1D* h_eta_23_prior;

  TH1D* h_j_e1;
  TH1D* h_j_e2;
  TH1D* h_j_mu1;
  TH1D* h_j_mu3;
  TH1D* h_j_tau2;
  TH1D* h_j_tau3;

  if(doModulus) {
    h_ue1 = new TH1D("h_ue1",";U_{e1}",100,0.,1.);
    h_ue2 = new TH1D("h_ue2",";U_{e2}",100,0.,1.);
    h_ue3 = new TH1D("h_ue3",";U_{e3}",100,0.,1.);
    h_umu1 = new TH1D("h_umu1",";U_{#mu1}",100,0.,1.);
    h_umu2 = new TH1D("h_umu2",";U_{#mu2}",100,0.,1.);
    h_umu3 = new TH1D("h_umu3",";U_{#mu3}",100,0.,1.);
    h_utau1 = new TH1D("h_utau1",";U_{#tau1}",100,0.,1.);
    h_utau2 = new TH1D("h_utau2",";U_{#tau2}",100,0.,1.);
    h_utau3 = new TH1D("h_utau3",";U_{#tau3}",100,0.,1.);
  
    h_ue1_prior = new TH1D("h_ue1_prior",";U_{e1}",100,0.,1.);
    h_ue2_prior = new TH1D("h_ue2_prior",";U_{e2}",100,0.,1.);
    h_ue3_prior = new TH1D("h_ue3_prior",";U_{e3}",100,0.,1.);
    h_umu1_prior = new TH1D("h_umu1_prior",";U_{#mu1}",100,0.,1.);
    h_umu2_prior = new TH1D("h_umu2_prior",";U_{#mu2}",100,0.,1.);
    h_umu3_prior = new TH1D("h_umu3_prior",";U_{#mu3}",100,0.,1.);
    h_utau1_prior = new TH1D("h_utau1_prior",";U_{#tau1}",100,0.,1.);
    h_utau2_prior = new TH1D("h_utau2_prior",";U_{#tau2}",100,0.,1.);
    h_utau3_prior = new TH1D("h_utau3_prior",";U_{#tau3}",100,0.,1.);
  }

  if(doReal) {
    h_real_ue1 = new TH1D("h_real_ue1",";Re(U_{e1})",100,-1.,1.);
    h_real_ue2 = new TH1D("h_real_ue2",";Re(U_{e2})",100,-1.,1.);
    h_real_ue3 = new TH1D("h_real_ue3",";Re(U_{e3})",100,-1.,1.);
    h_real_umu1 = new TH1D("h_real_umu1",";Re(U_{#mu1})",100,-1.,1.);
    h_real_umu2 = new TH1D("h_real_umu2",";Re(U_{#mu2})",100,-1.,1.);
    h_real_umu3 = new TH1D("h_real_umu3",";Re(U_{#mu3})",100,-1.,1.);
    h_real_utau1 = new TH1D("h_real_utau1",";Re(U_{#tau1})",100,-1.,1.);
    h_real_utau2 = new TH1D("h_real_utau2",";Re(U_{#tau2})",100,-1.,1.);
    h_real_utau3 = new TH1D("h_real_utau3",";Re(U_{#tau3})",100,-1.,1.);

    h_real_ue1_prior = new TH1D("h_real_ue1_prior",";Re(U_{e1})",100,-1.,1.);
    h_real_ue2_prior = new TH1D("h_real_ue2_prior",";Re(U_{e2})",100,-1.,1.);
    h_real_ue3_prior = new TH1D("h_real_ue3_prior",";Re(U_{e3})",100,-1.,1.);
    h_real_umu1_prior = new TH1D("h_real_umu1_prior",";Re(U_{#mu1})",100,-1.,1.);
    h_real_umu2_prior = new TH1D("h_real_umu2_prior",";Re(U_{#mu2})",100,-1.,1.);
    h_real_umu3_prior = new TH1D("h_real_umu3_prior",";Re(U_{#mu3})",100,-1.,1.);
    h_real_utau1_prior = new TH1D("h_real_utau1_prior",";Re(U_{#tau1})",100,-1.,1.);
    h_real_utau2_prior = new TH1D("h_real_utau2_prior",";Re(U_{#tau2})",100,-1.,1.);
    h_real_utau3_prior = new TH1D("h_real_utau3_prior",";Re(U_{#tau3})",100,-1.,1.);
  }

  if(doImag) {
    h_imag_ue1 = new TH1D("h_imag_ue1",";Im(U_{e1})",100,-1.,1.);
    h_imag_ue2 = new TH1D("h_imag_ue2",";Im(U_{e2})",100,-1.,1.);
    h_imag_ue3 = new TH1D("h_imag_ue3",";Im(U_{e3})",100,-1.,1.);
    h_imag_umu1 = new TH1D("h_imag_umu1",";Im(U_{#mu1})",100,-1.,1.);
    h_imag_umu2 = new TH1D("h_imag_umu2",";Im(U_{#mu2})",100,-1.,1.);
    h_imag_umu3 = new TH1D("h_imag_umu3",";Im(U_{#mu3})",100,-1.,1.);
    h_imag_utau1 = new TH1D("h_imag_utau1",";Im(U_{#tau1})",100,-1.,1.);
    h_imag_utau2 = new TH1D("h_imag_utau2",";Im(U_{#tau2})",100,-1.,1.);
    h_imag_utau3 = new TH1D("h_imag_utau3",";Im(U_{#tau3})",100,-1.,1.);

    h_imag_ue1_prior = new TH1D("h_imag_ue1_prior",";Im(U_{e1})",100,-1.,1.);
    h_imag_ue2_prior = new TH1D("h_imag_ue2_prior",";Im(U_{e2})",100,-1.,1.);
    h_imag_ue3_prior = new TH1D("h_imag_ue3_prior",";Im(U_{e3})",100,-1.,1.);
    h_imag_umu1_prior = new TH1D("h_imag_umu1_prior",";Im(U_{#mu1})",100,-1.,1.);
    h_imag_umu2_prior = new TH1D("h_imag_umu2_prior",";Im(U_{#mu2})",100,-1.,1.);
    h_imag_umu3_prior = new TH1D("h_imag_umu3_prior",";Im(U_{#mu3})",100,-1.,1.);
    h_imag_utau1_prior = new TH1D("h_imag_utau1_prior",";Im(U_{#tau1})",100,-1.,1.);
    h_imag_utau2_prior = new TH1D("h_imag_utau2_prior",";Im(U_{#tau2})",100,-1.,1.);
    h_imag_utau3_prior = new TH1D("h_imag_utau3_prior",";Im(U_{#tau3})",100,-1.,1.); 
  }

  if(doJarlskog) {
    h_j_e1 = new TH1D("h_j_e1",";J_{e1}",200,-0.05,0.05);
    h_j_e2 = new TH1D("h_j_e2",";J_{e2}",200,-0.05,0.05);
    h_j_mu1 = new TH1D("h_j_mu1",";J_{#mu1}",200,-0.05,0.05);
    h_j_mu3 = new TH1D("h_j_mu3",";J_{#mu3}",200,-0.05,0.05);
    h_j_tau2 = new TH1D("h_j_tau2",";J_{#tau2}",200,-0.05,0.05);
    h_j_tau3 = new TH1D("h_j_tau3",";J_{#tau3}",200,-0.05,0.05);
  }

  gStyle->SetOptStat(0);

  //triangles
  h_tr_emu = new TH2D("h_tr_emu",";#rho_{e#mu};#eta_{e#mu}",n_bins,-5,5,n_bins,-5,5);
  h_tr_etau = new TH2D("h_tr_etau",";#rho_{e#tau};#eta_{e#tau}",n_bins,0,2,n_bins,-1,1);
  h_tr_mutau = new TH2D("h_tr_mutau",";#rho_{#mu#tau};#eta_{#mu#tau}",n_bins,0,2,n_bins,-1,1);
  h_tr_12 = new TH2D("h_tr_12",";#rho_{12};#eta_{12}",n_bins,0,3.2,n_bins,-1.5,1.5);
  h_tr_13 = new TH2D("h_tr_13",";#rho_{13};#eta_{13}",n_bins,0,2,n_bins,-1,1);
  h_tr_23 = new TH2D("h_tr_23",";#rho_{23};#eta_{23}",n_bins,-8,8,n_bins,-8,8);

  h_tr_emu_prior = new TH2D("h_tr_emu_prior",";#rho_{e#mu};#eta_{e#mu}",n_bins,-5,5,n_bins,-5,5);
  h_tr_etau_prior = new TH2D("h_tr_etau_prior",";#rho_{e#tau};#eta_{e#tau}",n_bins,0,2,n_bins,-1,1);
  h_tr_mutau_prior = new TH2D("h_tr_mutau_prior",";#rho_{#mu#tau};#eta_{#mu#tau}",n_bins,0,2,n_bins,-1,1);
  h_tr_12_prior = new TH2D("h_tr_12_prior",";#rho_{12};#eta_{12}",n_bins,0,3.2,n_bins,-1.5,1.5);
  h_tr_13_prior = new TH2D("h_tr_13_prior",";#rho_{13};#eta_{13}",n_bins,0,2,n_bins,-1,1);
  h_tr_23_prior = new TH2D("h_tr_23_prior",";#rho_{23};#eta_{23}",n_bins,-8,8,n_bins,-8,8);

  gStyle->SetOptStat(1);

  TRandom3* randGen = new TRandom3();

  // nSteps = 50000000;

  for(unsigned int i = 0;i<nSteps; i++) {
  //for(unsigned int i = 0;i<50000000; i++) {

    //std::cout << "step # " << i <<std::endl;
    //if(i%100000==0) std::cout << "step # " << i <<std::endl;
    if(i%1000000==0) std::cout << "step # " << i <<std::endl;
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
      //std::cout<<s22th13<<"   "<<s13<<"     "<<c13<<std::endl;
    }
    else {
      s13 = TMath::Sqrt(s2th13);
      s23 = TMath::Sqrt(s2th23);
      // if (chain != "DUNE"){
      s12 = TMath::Sqrt(s2th12);
      // }
      // else{ 
      // s2th12 = randGen->Gaus(0.307,0.013);
      // s12 = TMath::Sqrt(s2th12);
      // }
      sdcp = TMath::Sin(dcp);
      cdcp = TMath::Cos(dcp);
      c13 = TMath::Sqrt(1.-s2th13);
      c12 = TMath::Sqrt(1.-s2th12);
      c23 = TMath::Sqrt(1.-s2th23);
      if (!RCreweight)
      {
        rc = 1.;
      }
    }

    //Uei
    real_ue1 = c12*c13;
    imag_ue1 = 0.;
    ue1 = TComplex(real_ue1,imag_ue1);

    real_ue2 = s12*c13;
    imag_ue2 = 0.;
    ue2 = TComplex(real_ue2,imag_ue2);

    real_ue3 = s13*cdcp;
    imag_ue3 = -s13*sdcp;
    ue3 = TComplex(real_ue3,imag_ue3);

    //Umui
    real_umu1 = -s12*c23-c12*s23*s13*cdcp;
    imag_umu1 = -c12*s23*s13*sdcp;
    umu1 = TComplex(real_umu1,imag_umu1);

    real_umu2 = c12*c23-s12*s23*s13*cdcp;
    imag_umu2 = -s12*s23*s13*sdcp;
    umu2 = TComplex(real_umu2,imag_umu2);

    real_umu3 = s23*c13;
    imag_umu3 = 0.;
    umu3 = TComplex(real_umu3,imag_umu3);

    //Utaui
    real_utau1 = s12*s23-c12*c23*s13*cdcp;
    imag_utau1 = -c12*c23*s13*sdcp;
    utau1 = TComplex(real_utau1,imag_utau1);

    real_utau2 = -c12*s23-s12*c23*s13*cdcp;
    imag_utau2 = -s12*c23*s13*sdcp;
    utau2 = TComplex(real_utau2,imag_utau2);

    real_utau3 = c23*c13;
    imag_utau3 = 0.;
    utau3 = TComplex(real_utau3,imag_utau3);

    
    //triangles
    tr_emu_num   = ue1.operator*(TComplex::Conjugate(umu1));
    tr_emu_denom = ue3.operator*(TComplex::Conjugate(umu3));
    tr_emu       = - tr_emu_num.operator/(tr_emu_denom);

    tr_etau_num   = ue2.operator*(TComplex::Conjugate(utau2));
    tr_etau_denom = ue1.operator*(TComplex::Conjugate(utau1));
    tr_etau       = - tr_etau_num.operator/(tr_etau_denom);

    tr_mutau_num   = umu3.operator*(TComplex::Conjugate(utau3));
    tr_mutau_denom = umu2.operator*(TComplex::Conjugate(utau2));
    tr_mutau       = - tr_mutau_num.operator/(tr_mutau_denom);

    tr_12_num   = ue1.operator*(TComplex::Conjugate(ue2));
    tr_12_denom = umu1.operator*(TComplex::Conjugate(umu2));
    tr_12       = - tr_12_num.operator/(tr_12_denom);

    tr_13_num   = umu1.operator*(TComplex::Conjugate(umu3));
    tr_13_denom = utau1.operator*(TComplex::Conjugate(utau3));
    tr_13       = - tr_13_num.operator/(tr_13_denom);

    tr_23_num   = utau2.operator*(TComplex::Conjugate(utau3));
    tr_23_denom = ue2.operator*(TComplex::Conjugate(ue3));
    tr_23       = - tr_23_num.operator/(tr_23_denom);

    //Jarlskog
    j_e1 = tr_mutau.Im()*TComplex::Abs(umu2)*TComplex::Abs(umu2)*TComplex::Abs(utau2)*TComplex::Abs(utau2);
    j_e2 = TComplex::Power(tr_13,-1).Im()*TComplex::Abs(umu1)*TComplex::Abs(umu1)*TComplex::Abs(umu3)*TComplex::Abs(umu3);
    j_mu1 = TComplex::Power(tr_23,-1).Im()*TComplex::Abs(utau2)*TComplex::Abs(utau2)*TComplex::Abs(utau3)*TComplex::Abs(utau3);
    j_mu3 = tr_etau.Im()*TComplex::Abs(ue1)*TComplex::Abs(ue1)*TComplex::Abs(utau1)*TComplex::Abs(utau1);
    j_tau2 = tr_emu.Im()*TComplex::Abs(ue3)*TComplex::Abs(ue3)*TComplex::Abs(umu3)*TComplex::Abs(umu3);
    j_tau3 = TComplex::Power(tr_12,-1).Im()*TComplex::Abs(ue1)*TComplex::Abs(ue1)*TComplex::Abs(ue2)*TComplex::Abs(ue2);

    if(doModulus) {
      h_ue1->Fill(TComplex::Abs(ue1),rc);
      h_ue2->Fill(TComplex::Abs(ue2),rc);
      h_ue3->Fill(TComplex::Abs(ue3),rc);
      h_umu1->Fill(TComplex::Abs(umu1),rc);
      h_umu2->Fill(TComplex::Abs(umu2),rc);
      h_umu3->Fill(TComplex::Abs(umu3),rc);
      h_utau1->Fill(TComplex::Abs(utau1),rc);
      h_utau2->Fill(TComplex::Abs(utau2),rc);
      h_utau3->Fill(TComplex::Abs(utau3),rc);
    }

    if(doReal) {
      h_real_ue1->Fill(real_ue1,rc);
      h_real_ue2->Fill(real_ue2,rc);
      h_real_ue3->Fill(real_ue3,rc);
      h_real_umu1->Fill(real_umu1,rc);
      h_real_umu2->Fill(real_umu2,rc);
      h_real_umu3->Fill(real_umu3,rc);
      h_real_utau1->Fill(real_utau1,rc);
      h_real_utau2->Fill(real_utau2,rc);
      h_real_utau3->Fill(real_utau3,rc);
    }


    if(doImag) {
      h_imag_ue1->Fill(imag_ue1,rc);
      h_imag_ue2->Fill(imag_ue2,rc);
      h_imag_ue3->Fill(imag_ue3,rc);
      h_imag_umu1->Fill(imag_umu1,rc);
      h_imag_umu2->Fill(imag_umu2,rc);
      h_imag_umu3->Fill(imag_umu3,rc);
      h_imag_utau1->Fill(imag_utau1,rc);
      h_imag_utau2->Fill(imag_utau2,rc);
      h_imag_utau3->Fill(imag_utau3,rc);
    }

    if(doJarlskog) {
      h_j_e1->Fill(j_e1,rc);
      h_j_e2->Fill(j_e2,rc);
      h_j_mu1->Fill(j_mu1,rc);
      h_j_mu3->Fill(j_mu3,rc);
      h_j_tau2->Fill(j_tau2,rc);
      h_j_tau3->Fill(j_tau3,rc);
    }

    gStyle->SetOptStat(0);

    //triangles
    h_tr_emu->Fill(tr_emu.Re(),tr_emu.Im(),rc);
    h_tr_etau->Fill(tr_etau.Re(),tr_etau.Im(),rc); 
    h_tr_mutau->Fill(tr_mutau.Re(),tr_mutau.Im(),rc);

    h_tr_12->Fill(tr_12.Re(),tr_12.Im(),rc);
    h_tr_13->Fill(tr_13.Re(),tr_13.Im(),rc);
    h_tr_23->Fill(tr_23.Re(),tr_23.Im(),rc);
  }

  if(doPrior) {
    //for(unsigned int s=0;s<50000000;s++){
    for(unsigned int s=0;s<50000000;s++){
      if(s%500000==0) std::cout << "throw # " << s <<std::endl;
      //sample from prior:

      if (fitter == "Aria") {
        prior_th13 = randGen->Uniform(0.,2*TMath::Pi());
        prior_wRC_th13 = randGen->Gaus(0.1476,0.0264);
        prior_th23 = randGen->Uniform(0.,2*TMath::Pi());
        prior_s2th12 = randGen->Gaus(0.307,0.013);
        prior_dcp = randGen->Uniform(0.,2*TMath::Pi());
        prior_sindcp = randGen->Uniform(-1.,1.);

        s13_prior = TMath::Sin(prior_th13);
        s13_wRC_prior = TMath::Sin(prior_wRC_th13);
        s23_prior = TMath::Sin(prior_th23);
        s12_prior = TMath::Sqrt(prior_s2th12);
        sdcp_prior = TMath::Sin(prior_dcp);
        cdcp_prior = TMath::Cos(prior_dcp);
        c13_prior = TMath::Cos(prior_th13);
        c13_wRC_prior = TMath::Cos(prior_wRC_th13);
        c12_prior = TMath::Sqrt(1.-prior_s2th12);
        c23_prior = TMath::Cos(prior_th23);
      }
      else {
        prior_s2th13 = randGen->Uniform(0.,1.);
        prior_wRC_s2th13 = randGen->Gaus(0.0218,0.0007);
        prior_s2th23 = randGen->Uniform(0.,1.);
        prior_s2th12 = randGen->Gaus(0.307,0.013);
        prior_dcp = randGen->Uniform(-1.*TMath::Pi(),TMath::Pi());
        prior_sindcp = randGen->Uniform(-1.,1.);

        s13_prior = TMath::Sqrt(prior_wRC_s2th13);
        s23_prior = TMath::Sqrt(prior_s2th23);
        s12_prior = TMath::Sqrt(prior_s2th12);
        sdcp_prior = TMath::Sin(prior_dcp);
        cdcp_prior = TMath::Cos(prior_dcp);
        c13_prior = TMath::Sqrt(1.-prior_s2th13);
        c12_prior = TMath::Sqrt(1.-prior_s2th12);
        c23_prior = TMath::Sqrt(1.-prior_s2th23);
      }

      //Uei
      real_ue1_prior = c12_prior*c13_prior;
      imag_ue1_prior = 0.;
      ue1_prior = TComplex(real_ue1_prior,imag_ue1_prior);
      
      real_ue2_prior = s12_prior*c13_prior;
      imag_ue2_prior = 0.;
      ue2_prior = TComplex(real_ue2_prior,imag_ue2_prior);
      
      real_ue3_prior = s13_prior*cdcp_prior;
      imag_ue3_prior = -s13_prior*sdcp_prior;
      ue3_prior = TComplex(real_ue3_prior,imag_ue3_prior);
  
      //Umui
      real_umu1_prior = -s12_prior*c23_prior-c12_prior*s23_prior*s13_prior*cdcp_prior;
      imag_umu1_prior = -c12_prior*s23_prior*s13_prior*sdcp_prior;
      umu1_prior = TComplex(real_umu1_prior,imag_umu1_prior);
      
      real_umu2_prior = c12_prior*c23_prior-s12_prior*s23_prior*s13_prior*cdcp_prior;
      imag_umu2_prior = -s12_prior*s23_prior*s13_prior*sdcp_prior;
      umu2_prior = TComplex(real_umu2_prior,imag_umu2_prior);
      
      real_umu3_prior = s23_prior*c13_prior;
      imag_umu3_prior = 0.;
      umu3_prior = TComplex(real_umu3_prior,imag_umu3_prior);
  
      //Utaui
      real_utau1_prior = s12_prior*s23_prior-c12_prior*c23_prior*s13_prior*cdcp_prior;
      imag_utau1_prior = -c12_prior*c23_prior*s13_prior*sdcp_prior;
      utau1_prior = TComplex(real_utau1_prior,imag_utau1_prior);
      
      real_utau2_prior = -c12_prior*s23_prior-s12_prior*c23_prior*s13_prior*cdcp_prior;
      imag_utau2_prior = -s12_prior*c23_prior*s13_prior*sdcp_prior;
      utau2_prior = TComplex(real_utau2_prior,imag_utau2_prior);
      
      real_utau3_prior = c23_prior*c13_prior;
      imag_utau3_prior = 0.;
      utau3_prior = TComplex(real_utau3_prior,imag_utau3_prior);

      //triangles
      tr_emu_num_prior   = ue1_prior.operator*(TComplex::Conjugate(umu1_prior));
      tr_emu_denom_prior = ue3_prior.operator*(TComplex::Conjugate(umu3_prior));
      tr_emu_prior       = - tr_emu_num_prior.operator/(tr_emu_denom_prior);

      tr_etau_num_prior   = ue2_prior.operator*(TComplex::Conjugate(utau2_prior));
      tr_etau_denom_prior = ue1_prior.operator*(TComplex::Conjugate(utau1_prior));
      tr_etau_prior       = - tr_etau_num_prior.operator/(tr_etau_denom_prior);

      tr_mutau_num_prior   = umu3_prior.operator*(TComplex::Conjugate(utau3_prior));
      tr_mutau_denom_prior = umu2_prior.operator*(TComplex::Conjugate(utau2_prior));
      tr_mutau_prior       = - tr_mutau_num_prior.operator/(tr_mutau_denom_prior);

      tr_12_num_prior   = ue1_prior.operator*(TComplex::Conjugate(ue2_prior));
      tr_12_denom_prior = umu1_prior.operator*(TComplex::Conjugate(umu2_prior));
      tr_12_prior       = - tr_12_num_prior.operator/(tr_12_denom_prior);

      tr_13_num_prior   = umu1_prior.operator*(TComplex::Conjugate(umu3_prior));
      tr_13_denom_prior = utau1_prior.operator*(TComplex::Conjugate(utau3_prior));
      tr_13_prior       = - tr_13_num_prior.operator/(tr_13_denom_prior);

      tr_23_num_prior   = utau2_prior.operator*(TComplex::Conjugate(utau3_prior));
      tr_23_denom_prior = ue2_prior.operator*(TComplex::Conjugate(ue3_prior));
      tr_23_prior       = - tr_23_num_prior.operator/(tr_23_denom_prior);

      ////////
      //PLOT//
      ////////

      gStyle->SetOptStat(1);

      if(doModulus) {
        h_ue1_prior->Fill(TComplex::Abs(ue1_prior));
        h_ue2_prior->Fill(TComplex::Abs(ue2_prior));
        h_ue3_prior->Fill(TComplex::Abs(ue3_prior));
        h_umu1_prior->Fill(TComplex::Abs(umu1_prior));
        h_umu2_prior->Fill(TComplex::Abs(umu2_prior));
        h_umu3_prior->Fill(TComplex::Abs(umu3_prior));
        h_utau1_prior->Fill(TComplex::Abs(utau1_prior));
        h_utau2_prior->Fill(TComplex::Abs(utau2_prior));
        h_utau3_prior->Fill(TComplex::Abs(utau3_prior));
      }
      if(doReal) {
        h_real_ue1_prior->Fill(real_ue1_prior);
        h_real_ue2_prior->Fill(real_ue2_prior);
        h_real_ue3_prior->Fill(real_ue3_prior);
        h_real_umu1_prior->Fill(real_umu1_prior);
        h_real_umu2_prior->Fill(real_umu2_prior);
        h_real_umu3_prior->Fill(real_umu3_prior);
        h_real_utau1_prior->Fill(real_utau1_prior);
        h_real_utau2_prior->Fill(real_utau2_prior);
        h_real_utau3_prior->Fill(real_utau3_prior);
      }
      if(doImag) {
        h_imag_ue1_prior->Fill(imag_ue1_prior);
        h_imag_ue2_prior->Fill(imag_ue2_prior);
        h_imag_ue3_prior->Fill(imag_ue3_prior);
        h_imag_umu1_prior->Fill(imag_umu1_prior);
        h_imag_umu2_prior->Fill(imag_umu2_prior);
        h_imag_umu3_prior->Fill(imag_umu3_prior);
        h_imag_utau1_prior->Fill(imag_utau1_prior);
        h_imag_utau2_prior->Fill(imag_utau2_prior);
        h_imag_utau3_prior->Fill(imag_utau3_prior);
      }

      gStyle->SetOptStat(0);

      //triangles
      h_tr_emu_prior->Fill(tr_emu_prior.Re(),tr_emu_prior.Im());
      h_tr_etau_prior->Fill(tr_etau_prior.Re(),tr_etau_prior.Im());
      h_tr_mutau_prior->Fill(tr_mutau_prior.Re(),tr_mutau_prior.Im());

      h_tr_12_prior->Fill(tr_12_prior.Re(),tr_12_prior.Im());
      h_tr_13_prior->Fill(tr_13_prior.Re(),tr_13_prior.Im());
      h_tr_23_prior->Fill(tr_23_prior.Re(),tr_23_prior.Im());
    }
  }

  outFile->cd();

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

  if(doJarlskog) {
    h_j_e1->Write("h_j_e1");
    h_j_e2->Write("h_j_e2");
    h_j_mu1->Write("h_j_mu1");
    h_j_mu3->Write("h_j_mu3");
    h_j_tau2->Write("h_j_tau2");
    h_j_tau3->Write("h_j_tau3");
  }

  //triangles
  h_tr_emu->Write("h_tr_emu");
  h_tr_etau->Write("h_tr_etau");
  h_tr_mutau->Write("h_tr_mutau");

  h_tr_12->Write("h_tr_12");
  h_tr_13->Write("h_tr_13");
  h_tr_23->Write("h_tr_23");

  //triangle priors
  h_tr_emu_prior->Write("h_tr_emu_prior");
  h_tr_etau_prior->Write("h_tr_etau_prior");
  h_tr_mutau_prior->Write("h_tr_mutau_prior");

  h_tr_12_prior->Write("h_tr_12_prior");
  h_tr_13_prior->Write("h_tr_13_prior");
  h_tr_23_prior->Write("h_tr_23_prior");

  //rho
  h_tr_emu->ProjectionX()->Write("h_rho_emu");
  h_tr_etau->ProjectionX()->Write("h_rho_etau");
  h_tr_mutau->ProjectionX()->Write("h_rho_mutau");

  h_tr_12->ProjectionX()->Write("h_rho_12");
  h_tr_13->ProjectionX()->Write("h_rho_13");
  h_tr_23->ProjectionX()->Write("h_rho_23");

  //rho priors
  h_tr_emu_prior->ProjectionX()->Write("h_rho_emu_prior");
  h_tr_etau_prior->ProjectionX()->Write("h_rho_etau_prior");
  h_tr_mutau_prior->ProjectionX()->Write("h_rho_mutau_prior");

  h_tr_12_prior->ProjectionX()->Write("h_rho_12_prior");
  h_tr_13_prior->ProjectionX()->Write("h_rho_13_prior");
  h_tr_23_prior->ProjectionX()->Write("h_rho_23_prior");

  //eta
  h_tr_emu->ProjectionY()->Write("h_eta_emu");
  h_tr_etau->ProjectionY()->Write("h_eta_etau");
  h_tr_mutau->ProjectionY()->Write("h_eta_mutau");

  h_tr_12->ProjectionY()->Write("h_eta_12");
  h_tr_13->ProjectionY()->Write("h_eta_13");
  h_tr_23->ProjectionY()->Write("h_eta_23");

  //eta priors
  h_tr_emu_prior->ProjectionY()->Write("h_eta_emu_prior");
  h_tr_etau_prior->ProjectionY()->Write("h_eta_etau_prior");
  h_tr_mutau_prior->ProjectionY()->Write("h_eta_mutau_prior");

  h_tr_12_prior->ProjectionY()->Write("h_eta_12_prior");
  h_tr_13_prior->ProjectionY()->Write("h_eta_13_prior");
  h_tr_23_prior->ProjectionY()->Write("h_eta_23_prior");
}
