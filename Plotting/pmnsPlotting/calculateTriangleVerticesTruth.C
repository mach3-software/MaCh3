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

void calculateTriangleVerticesTruth(){

  std::map<std::string, std::map<std::string, double>> pmns_pars;
  pmns_pars["AsimovA22"]["s2th12"] = 0.307;
  pmns_pars["AsimovA22"]["s2th13"] = 0.0220;
  pmns_pars["AsimovA22"]["s2th23"] = 0.561;
  pmns_pars["AsimovA22"]["dcp"] = -1.601;

  pmns_pars["AsimovB22"]["s2th12"] = 0.307;
  pmns_pars["AsimovB22"]["s2th13"] = 0.0220;
  pmns_pars["AsimovB22"]["s2th23"] = 0.450;
  pmns_pars["AsimovB22"]["dcp"] = 0.0;

  pmns_pars["FD+ND"]["s2th12"] = 0.31;       // These values are taken from NuFIT 4.0 and were used for the DUNE TDR analysis
  pmns_pars["FD+ND"]["s2th13"] = 0.0224;
  pmns_pars["FD+ND"]["s2th23"] = 0.582;
  pmns_pars["FD+ND"]["dcp"] = -2.498;

  std::string outFileName = "unitarityTriangles_truth.root";
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  // Loop over the outer map
  for (const auto& fit_type_pair : pmns_pars) {
    const std::string& fit_type = fit_type_pair.first;
    std::cout << "Fit type: " << fit_type << std::endl;

    std::cout<<"Calculating true values of PMNS elements for "<<fit_type<<std::endl;
    std::cout<<"True PMNS parameters:"<<std::endl;
    for (auto const& pmns_par : pmns_pars[fit_type]){
      std::cout << pmns_par.first << ": " << pmns_par.second << std::endl;
    }

    double s2th12 = pmns_pars[fit_type]["s2th12"];
    double s2th13 = pmns_pars[fit_type]["s2th13"];
    double s2th23 = pmns_pars[fit_type]["s2th23"];
    double dcp = pmns_pars[fit_type]["dcp"];

    double s12 = TMath::Sqrt(s2th12);
    double s13 = TMath::Sqrt(s2th13);
    double s23 = TMath::Sqrt(s2th23);
    double sdcp = TMath::Sin(dcp);
    double c12 = TMath::Sqrt(1.-s2th12);
    double c13 = TMath::Sqrt(1.-s2th13);
    double c23 = TMath::Sqrt(1.-s2th23);
    double cdcp = TMath::Cos(dcp);

    double real_ue1 = c12*c13;
    double imag_ue1 = 0.;
    double real_ue2 = s12*c13;
    double imag_ue2 = 0.;
    double real_ue3 = s13*cdcp;
    double imag_ue3 = -s13*sdcp;

    double real_umu1 = -s12*c23-c12*s23*s13*cdcp;
    double imag_umu1 = -c12*s23*s13*sdcp;
    double real_umu2 = c12*c23-s12*s23*s13*cdcp;
    double imag_umu2 = -s12*s23*s13*sdcp;
    double real_umu3 = s23*c13;
    double imag_umu3 = 0.;

    double real_utau1 = s12*s23-c12*c23*s13*cdcp;
    double imag_utau1 = -c12*c23*s13*sdcp;
    double real_utau2 = -c12*s23-s12*c23*s13*cdcp;
    double imag_utau2 = -s12*c23*s13*sdcp;
    double real_utau3 = c23*c13;
    double imag_utau3 = 0.;

    TComplex ue1, ue2, ue3, umu1, umu2, umu3, utau1, utau2, utau3;
    TComplex tr_emu_num, tr_etau_num, tr_mutau_num, tr_12_num, tr_13_num, tr_23_num;
    TComplex tr_emu_denom, tr_etau_denom, tr_mutau_denom, tr_12_denom, tr_13_denom, tr_23_denom;
    TComplex tr_emu, tr_etau, tr_mutau, tr_12, tr_13, tr_23;

    ue1 = TComplex(real_ue1,imag_ue1);
    ue2 = TComplex(real_ue2,imag_ue2);
    ue3 = TComplex(real_ue3,imag_ue3);
    umu1 = TComplex(real_umu1,imag_umu1);
    umu2 = TComplex(real_umu2,imag_umu2);
    umu3 = TComplex(real_umu3,imag_umu3);
    utau1 = TComplex(real_utau1,imag_utau1);
    utau2 = TComplex(real_utau2,imag_utau2);
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

    double eta_emu = tr_emu.Im();
    double rho_emu = tr_emu.Re();
    double eta_etau = tr_etau.Im();
    double rho_etau = tr_etau.Re();
    double eta_mutau = tr_mutau.Im();
    double rho_mutau = tr_mutau.Re();

    double eta_12 = tr_12.Im();
    double rho_12 = tr_12.Re();
    double eta_13 = tr_13.Im();
    double rho_13 = tr_13.Re();
    double eta_23 = tr_23.Im();
    double rho_23 = tr_23.Re();

    std::stringstream tree_name;
    tree_name << fit_type << "_truth";

    TTree* pmnsTree = new TTree(tree_name.str().c_str(),tree_name.str().c_str());
    pmnsTree->Branch("eta_emu",&eta_emu);
    pmnsTree->Branch("rho_emu",&rho_emu);
    pmnsTree->Branch("eta_etau",&eta_etau);
    pmnsTree->Branch("rho_etau",&rho_etau);
    pmnsTree->Branch("eta_mutau",&eta_mutau);
    pmnsTree->Branch("rho_mutau",&rho_mutau);

    pmnsTree->Branch("eta_12",&eta_12);
    pmnsTree->Branch("rho_12",&rho_12);
    pmnsTree->Branch("eta_13",&eta_13);
    pmnsTree->Branch("rho_13",&rho_13);
    pmnsTree->Branch("eta_23",&eta_23);
    pmnsTree->Branch("rho_23",&rho_23);

    pmnsTree->Branch("s2th12",&s2th12);
    pmnsTree->Branch("s2th13",&s2th13);
    pmnsTree->Branch("s2th23",&s2th23);
    pmnsTree->Branch("dcp",&dcp);

    pmnsTree->Fill();

    outFile->cd();
    pmnsTree->Write();

    // Write to CSV file
    std::ofstream csvFile;
    std::stringstream csvFileName;
    csvFileName << "unitarityTriangles_" << fit_type << "_truth.csv";
    csvFile.open(csvFileName.str());

    csvFile << "emu," << rho_emu << " + " << eta_emu << " i"<< "\n";
    csvFile << "etau," << rho_etau << " + " << eta_etau << " i"<< "\n";
    csvFile << "mutau," << rho_mutau << " + " << eta_mutau << " i"<< "\n";
    csvFile << "12," << rho_12 << " + " << eta_12 << " i"<< "\n";
    csvFile << "13," << rho_13 << " + " << eta_13 << " i"<< "\n";
    csvFile << "23," << rho_23 << " + " << eta_23 << " i"<< "\n";
    csvFile.close();

    std::cout << "CSV file " << csvFileName.str() << " has been written." << std::endl;

    std::cout << std::endl;

  }

  outFile->Close();

  std::cout<<"Output file "<<outFileName<<" has been written"<<std::endl;

 }
