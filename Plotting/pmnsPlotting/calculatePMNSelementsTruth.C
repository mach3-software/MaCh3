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

void calculatePMNSelementsTruth(){

  std::map<std::string, std::map<std::string, double>> pmns_pars;
  pmns_pars["AsimovA22"]["s2th12"] = 0.307;
  pmns_pars["AsimovA22"]["s2th13"] = 0.0220;
  pmns_pars["AsimovA22"]["s2th23"] = 0.561;
  pmns_pars["AsimovA22"]["dcp"] = -1.601;

  pmns_pars["AsimovB22"]["s2th12"] = 0.307;
  pmns_pars["AsimovB22"]["s2th13"] = 0.0220;
  pmns_pars["AsimovB22"]["s2th23"] = 0.450;
  pmns_pars["AsimovB22"]["dcp"] = 0.0;

  std::string outFileName = "pmns_truth.root";
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

    double ue1 = sqrt(real_ue1*real_ue1+imag_ue1*imag_ue1);
    double ue2 = sqrt(real_ue2*real_ue2+imag_ue2*imag_ue2);
    double ue3 = sqrt(real_ue3*real_ue3+imag_ue3*imag_ue3);

    double umu1 = sqrt(real_umu1*real_umu1+imag_umu1*imag_umu1);
    double umu2 = sqrt(real_umu2*real_umu2+imag_umu2*imag_umu2);
    double umu3 = sqrt(real_umu3*real_umu3+imag_umu3*imag_umu3);

    double utau1 = sqrt(real_utau1*real_utau1+imag_utau1*imag_utau1);
    double utau2 = sqrt(real_utau2*real_utau2+imag_utau2*imag_utau2);
    double utau3 = sqrt(real_utau3*real_utau3+imag_utau3*imag_utau3);

    std::stringstream tree_name;
    tree_name << fit_type << "_truth";

    TTree* pmnsTree = new TTree(tree_name.str().c_str(),tree_name.str().c_str());
    pmnsTree->Branch("ue1",&ue1);
    pmnsTree->Branch("ue2",&ue2);
    pmnsTree->Branch("ue3",&ue3);
    pmnsTree->Branch("umu1",&umu1);
    pmnsTree->Branch("umu2",&umu2);
    pmnsTree->Branch("umu3",&umu3);
    pmnsTree->Branch("utau1",&utau1);
    pmnsTree->Branch("utau2",&utau2);
    pmnsTree->Branch("utau3",&utau3);

    pmnsTree->Branch("ue1_real",&real_ue1);
    pmnsTree->Branch("ue2_real",&real_ue2);
    pmnsTree->Branch("ue3_real",&real_ue3);
    pmnsTree->Branch("umu1_real",&real_umu1);
    pmnsTree->Branch("umu2_real",&real_umu2);
    pmnsTree->Branch("umu3_real",&real_umu3);
    pmnsTree->Branch("utau1_real",&real_utau1);
    pmnsTree->Branch("utau2_real",&real_utau2);
    pmnsTree->Branch("utau3_real",&real_utau3);

    pmnsTree->Branch("ue1_imag",&imag_ue1);
    pmnsTree->Branch("ue2_imag",&imag_ue2);
    pmnsTree->Branch("ue3_imag",&imag_ue3);
    pmnsTree->Branch("umu1_imag",&imag_umu1);
    pmnsTree->Branch("umu2_imag",&imag_umu2);
    pmnsTree->Branch("umu3_imag",&imag_umu3);
    pmnsTree->Branch("utau1_imag",&imag_utau1);
    pmnsTree->Branch("utau2_imag",&imag_utau2);
    pmnsTree->Branch("utau3_imag",&imag_utau3);

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
    csvFileName << "pmnsElements_" << fit_type << "_truth.csv";
    csvFile.open(csvFileName.str());

    csvFile << ",1,2,3\n";
    csvFile << "e," << ue1 << "," << ue2 << "," << ue3 << "\n";
    csvFile << "mu," << umu1 << "," << umu2 << "," << umu3 << "\n";
    csvFile << "tau," << utau1 << "," << utau2 << "," << utau3 << "\n";

    csvFile.close();

    std::cout << "CSV file " << csvFileName.str() << " has been written." << std::endl;

    std::cout << std::endl;

  }

  outFile->Close();

  std::cout<<"Output file "<<outFileName<<" has been written"<<std::endl;

 }
