#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <sstream>

double* getIntervals2D(TH2D *hist, double &p68, double &p90, double &p95, double &p99 ,double &p3sig)//, TH2D &h68, TH2D &h90)
{ 
  std::cout << "getting interval" << std::endl;
  TH2D *hCopy = (TH2D*)hist->Clone("hCopy");
  TH2D *hCont68 = (TH2D*)hist->Clone("hCont68");
  TH2D *hCont90 = (TH2D*)hist->Clone("hCont90");
  TH2D *hCont95 = (TH2D*)hist->Clone("hCont95");
  TH2D *hCont99 = (TH2D*)hist->Clone("hCont99");
  TH2D *hCont3sig = (TH2D*)hist->Clone("hCont3sig");
  
  double integral = hCopy->Integral();
  double tsum = 0; 
  double cont68lev = 0;// array('d', [0.0])
  double cont90lev = 0;//array('d', [0.0])
  double cont95lev = 0;//array('d', [0.0])
  double cont99lev = 0;//array('d', [0.0])
  double cont3siglev = 0;//array('d', [0.0])
  
  std::cout << integral << " " << tsum << std::endl;
  
  while ((tsum / float(integral)) < 0.9973)
    { 
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / float(integral) < 0.68)
        { 
          cont68lev = tmax;
          hCopy->SetBinContent(bin, -1.0);
        }
      if ((tsum / float(integral) < 0.9) && (tsum / float(integral) > 0.68))
        { 
          cont90lev = tmax;
          hCopy->SetBinContent(bin, -3.0);
        }
      if ((tsum / float(integral) < 0.955) && (tsum / float(integral) > 0.9))
        { 
          cont95lev = tmax;
          hCopy->SetBinContent(bin, -5.0);
        }
      if ((tsum / float(integral) < 0.99) && (tsum / float(integral) > 0.955))
        { 
          cont99lev = tmax;
          hCopy->SetBinContent(bin, -7.0);
        }       
      if ((tsum / float(integral) < 0.9973) && (tsum / float(integral) > 0.99))
        { 
          cont3siglev = tmax;
          hCopy->SetBinContent(bin, -9.0);
        }
    
    }
  
  double quant[5];
  quant[0] = cont90lev;
  quant[1] = cont68lev;
  
  quant[2] = cont95lev;
  quant[3] = cont99lev;
  quant[4] = cont3siglev;
  
  p90 = cont90lev;
  p68 = cont68lev;
  
  p95 = cont95lev;
  p99 = cont99lev;
  p3sig = cont3siglev;
  
  std::cout << "p68 = " << p68 << ", p90 = " << p90 << std::endl;
  std::cout << "p95 = " << p95 << ", p99 = " << p99 << std::endl;
  std::cout << "p3sig = " << p3sig << std::endl;
  
  return quant;//, hCont68, hCont90;
}

TLegend* plotIntervals2D(TH2D *hist, TLegend* leg, int levOpt)//, TH2D &h68, TH2D &h90)
{
  double p68, p90, p95, p99, p3sig;
  double *cont;
  
  cont = getIntervals2D(hist, p68, p90, p95, p99, p3sig);

  double tpp68[1];
  double tpp90[1];
  double tpp95[1];
  double tpp99[1];
  double tpp3sig[1];

  tpp68[0] = p68;
  tpp90[0] = p90;
  tpp95[0] = p95;
  tpp99[0] = p99;
  tpp3sig[0] = p3sig;

  TH2D* h_tr_cont_1sig;
  TH2D* h_tr_cont_2sig;
  TH2D* h_tr_cont_3sig;

  h_tr_cont_1sig = (TH2D*)hist->Clone("h_tr_cont_1sig");
  h_tr_cont_2sig = (TH2D*)hist->Clone("h_tr_cont_2sig");
  h_tr_cont_3sig = (TH2D*)hist->Clone("h_tr_cont_3sig");

  h_tr_cont_1sig->Smooth(1);
  h_tr_cont_2sig->Smooth(1);
  h_tr_cont_3sig->Smooth(1);

  if(levOpt==1) {
    h_tr_cont_1sig->SetContour(1,tpp68);
    h_tr_cont_2sig->SetContour(1,tpp90);
    h_tr_cont_3sig->SetContour(1,tpp99);

    leg->AddEntry(h_tr_cont_1sig,"68% C.I.","l");
    leg->AddEntry(h_tr_cont_2sig,"90% C.I.","l");
    leg->AddEntry(h_tr_cont_3sig,"99% C.I.","l");
  }
  else {
    h_tr_cont_1sig->SetContour(1,tpp68);
    h_tr_cont_2sig->SetContour(1,tpp95);
    h_tr_cont_3sig->SetContour(1,tpp3sig);

    leg->AddEntry(h_tr_cont_1sig,"1#sigma C.I.","l");
    leg->AddEntry(h_tr_cont_2sig,"2#sigma C.I.","l");
    leg->AddEntry(h_tr_cont_3sig,"3#sigma C.I.","l");
  }

  // h_tr_cont_1sig->Rebin2D(2);
  // h_tr_cont_2sig->Rebin2D(2);
  // h_tr_cont_3sig->Rebin2D(2);

  // Convert HEX #00FFFF → RGB (0, 1, 1)
  int cyan = TColor::GetFreeColorIndex();
  new TColor(cyan, 0.0, 1.0, 1.0);  // Define Cyan
  
  h_tr_cont_1sig->SetLineColor(cyan);
  h_tr_cont_2sig->SetLineColor(cyan);
  h_tr_cont_3sig->SetLineColor(cyan);
  h_tr_cont_1sig->SetLineWidth(2);
  h_tr_cont_2sig->SetLineWidth(2);
  h_tr_cont_3sig->SetLineWidth(2);
  h_tr_cont_1sig->SetLineStyle(1);
  h_tr_cont_2sig->SetLineStyle(2);
  h_tr_cont_3sig->SetLineStyle(3);

  h_tr_cont_1sig->Draw("cont3 LIST same");
  h_tr_cont_2sig->Draw("cont3 LIST same");
  h_tr_cont_3sig->Draw("cont3 LIST same");

  return leg;
}


void drawUnitarityTriangles(int n_tr = 0, std::string fit_type = "OA2021", bool RCreweight = true, bool draw_hist = true, bool draw_cont = true, bool drawTruth = true, int levOpt=1){

  int n_bins = 200;

  Float_t font_size = 0.05;

  std::vector<std::string> histnames = {"h_tr_emu","h_tr_etau","h_tr_mutau","h_tr_12","h_tr_13","h_tr_23"};
  std::map<int, std::string> triangle{{0, "emu"}, {1, "etau"}, {2, "mutau"}, {3, "12"}, {4, "13"}, {5, "23"}};
  std::map<int, std::string> levOpt_str{{0,"-w1sig2sig3sigContours"},{1,"-w68per90per99perContours"}};
  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};

  std::string chain = "T2K";

  TFile* inFile;
  std::stringstream infilename;

  if (chain == "NOvA" || chain == "DUNE"){
    infilename<<"unitarity_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_woPriors_"<<n_bins<<"bins.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
  }
  else if (chain == "PDG"){
    infilename<<"unitarity_pdg_"<<n_bins<<"bins.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
    drawTruth = false;
    RCreweight = false;
    n_bins = 1000;
  }
  else {
    infilename<<"unitarity_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_woPriors_"<<n_bins<<"bins.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
  }

  std::cout<<"Input file: "<<infilename.str().c_str()<<std::endl;

  TCanvas* c = new TCanvas("c","c",1350,1200);
  gStyle->SetPadLeftMargin(.09);
  // gStyle->SetPadLeftMargin(.1);
  gStyle->SetPadRightMargin(0.17);
  // gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(.11);
  gStyle->SetPalette(kCividis);
  c->UseCurrentStyle();
  c->Draw();
  c->cd();

  TH2D* h_tr;

  TH1D* h_rho;
  TH1D* h_eta;

  double peakpos[2];
  double peakpos_proj[2];
  int peakheight_proj[2];

  gStyle->SetOptStat(0);

  h_tr = (TH2D*)inFile->Get(histnames[n_tr].c_str());
  // h_tr->Rebin2D(3);
  // Add +1 to each bin to get rid of white spaces in the plot
  for (int i = 1; i <= h_tr->GetNcells(); i++) {  // Loop over bins (ignore underflow/overflow)
    h_tr->SetBinContent(i, h_tr->GetBinContent(i) + 1E-20);
  }

  h_rho = h_tr->ProjectionX();
  h_eta = h_tr->ProjectionY();

  h_tr->Scale(1./h_tr->Integral());

  std::vector<TH1D*> projections = {h_rho,h_eta};

  Int_t MaxBin_2d = h_tr->GetMaximumBin();
  Int_t x_2d,y_2d,z_2d;
  h_tr->GetBinXYZ(MaxBin_2d,x_2d,y_2d,z_2d);

  peakpos[0] = h_rho->GetBinCenter(x_2d);
  peakpos[1] = h_eta->GetBinCenter(y_2d);

  std::cout<<"2D Peak Position: ("<<peakpos[0]<<","<<peakpos[1]<<")"<<std::endl;

  for(unsigned int l=0;l<2;l++){ 
    Int_t MaxBin = projections[l]->GetMaximumBin();

    peakpos_proj[l] = projections[l]->GetBinCenter(MaxBin);
    peakheight_proj[l] = projections[l]->GetBinContent(MaxBin);

    std::cout<<"Projection "<<l<<" - Peak Position: "<<peakpos_proj[l]<<std::endl;
    std::cout<<std::endl;
  }

  h_tr->SetContour(255);
  h_tr->Smooth(1);
  // h_tr->GetZaxis()->SetTickSize(0.);
  // h_tr->GetZaxis()->SetLabelSize(0.);
  h_tr->GetZaxis()->SetTitle("Posterior Probability   ");
  h_tr->GetZaxis()->SetTitleOffset(1.2);
  h_tr->GetZaxis()->SetMaxDigits(1);
  std::cout<<"Max Digits:"<<h_tr->GetZaxis()->GetMaxDigits()<<std::endl;
  h_tr->GetYaxis()->SetTitleOffset(0.8);
  //gPad->SetLogz(1);

  h_tr->GetYaxis()->CenterTitle();
  h_tr->GetXaxis()->CenterTitle();

  h_tr->GetXaxis()->SetLabelSize(font_size);
  h_tr->GetXaxis()->SetTitleSize(font_size);
  h_tr->GetYaxis()->SetLabelSize(font_size);
  h_tr->GetYaxis()->SetTitleSize(font_size);
  h_tr->GetZaxis()->SetLabelSize(font_size);
  h_tr->GetZaxis()->SetTitleSize(font_size);

  std::cout<<"Title size:"<<h_tr->GetXaxis()->GetTitleSize()<<std::endl;
  std::cout<<"Label size:"<<h_tr->GetXaxis()->GetLabelSize()<<std::endl;

  // Change X and Y axis tick color
  h_tr->GetXaxis()->SetAxisColor(kWhite);   // Change X-axis tick color to white
  h_tr->GetYaxis()->SetAxisColor(kWhite);   // Change Y-axis tick color to white
  h_tr->GetZaxis()->SetAxisColor(kWhite);   // Change Z-axis tick color to white
  gPad->SetTicks();

  int n_leg_entries = 0;
  if (draw_cont) {
    n_leg_entries = n_leg_entries + 3;
  }
  if (((fit_type.find("Asimov") != std::string::npos) || chain == "DUNE" ) && drawTruth) n_leg_entries++;

  TLegend* leg;
  if (n_tr != 5){
    // leg = new TLegend(0.13,0.82,0.36,0.97);
    leg = new TLegend(0.12,0.92-0.12*n_leg_entries/2,0.35,0.92);
  }
  else{
    leg = new TLegend(0.115,0.13,0.35,0.13+0.12*n_leg_entries/2);
  }

  if (draw_hist)
  {
    h_tr->Draw("COLZ");
  }

  double rho_val = 0.0;
  double eta_val = 0.0;

  std::string fileAdd_truth;

  TMarker* mrkr_truth;

  if (((fit_type.find("Asimov") != std::string::npos) || chain == "DUNE" ) && drawTruth){
    std::cout<<"Using Asimov data set-->Draw marker at true value"<<std::endl;
    TLine* lin;
    TFile* inFile = new TFile("unitarityTriangles_truth.root","READ");
    std::stringstream tree_name;
    std::stringstream rho_branch_name;
    std::stringstream eta_branch_name;
    tree_name << fit_type << "_truth";
    rho_branch_name << "rho_" << triangle[n_tr];
    eta_branch_name << "eta_" << triangle[n_tr];

    TTree* pmnsTree = (TTree*)inFile->Get(tree_name.str().c_str());
    pmnsTree->SetBranchAddress(rho_branch_name.str().c_str(),&rho_val);
    pmnsTree->SetBranchAddress(eta_branch_name.str().c_str(),&eta_val);
    pmnsTree->GetEntry(0);

    // Convert HEX #00FFFF → RGB (0, 1, 1)
    int cyan = TColor::GetFreeColorIndex();
    new TColor(cyan, 0.0, 1.0, 1.0);  // Define Cyan

    mrkr_truth = new TMarker(rho_val,eta_val,29);
    mrkr_truth->SetMarkerSize(4);
    mrkr_truth->SetMarkerColor(cyan);

    mrkr_truth->Draw("same");

    if (n_tr == 5){
      leg->AddEntry(mrkr_truth,"Truth","p");
    }

    fileAdd_truth = "_"+bool_str[drawTruth]+"Truth";
  }

  if (draw_cont)
  {
    leg = plotIntervals2D(h_tr,leg,levOpt);
    c->Update();
  }

  if (drawTruth && n_tr != 5 && ((fit_type.find("Asimov") != std::string::npos) || chain == "DUNE" )){
    leg->AddEntry(mrkr_truth,"Truth","p");
  }

  TLine* lin_base;
  TLine* lin_1;
  TLine* lin_2;
  lin_base = new TLine(0.,0.,1.,0.);
  lin_1 = new TLine(0.,0.,peakpos[0],peakpos[1]);
  lin_2 = new TLine(1.,0.,peakpos[0],peakpos[1]);
  lin_base->SetLineWidth(2);
  lin_1->SetLineWidth(2);
  lin_2->SetLineWidth(2);
  lin_base->SetLineColor(kWhite);
  lin_1->SetLineColor(kWhite);
  lin_2->SetLineColor(kWhite);
  lin_base->Draw("same");
  lin_1->Draw("same");
  lin_2->Draw("same");

  // Change legend text color
  leg->SetTextColor(kWhite);    // Change text color to white
  leg->SetTextSize(font_size*0.8);  // Change text size

  leg->SetTextSize(font_size*0.8);  // Change text size
  leg->SetFillStyle(0);  // Removes the fill color (makes it transparent)
  leg->SetBorderSize(0); // Removes the border
  leg->Draw();

  TLatex* latex = new TLatex();
  latex->SetNDC();                // Use normalized coordinates
  latex->SetTextColor(kWhite);    //light gray text
  latex->SetTextFont(leg->GetTextFont()+10); // Set font to Helvetica
  latex->SetTextSize(font_size*0.9); // Set text size
  latex->DrawLatex(0.51, 0.89, "assuming unitarity");  //#splitline{assuming}{unitarity}

  std::cout<<"Left:"<<gStyle->GetPadLeftMargin()<<std::endl;
  std::cout<<"Right:"<<gStyle->GetPadRightMargin()<<std::endl;
  std::cout<<"YTitle:"<<h_tr->GetYaxis()->GetTitleOffset()<<std::endl;

  std::string fileAdd[2];

  if (draw_hist)
  {
    fileAdd[0] = "-2dHist";
  }
  if (draw_cont)
  {
    fileAdd[1] = levOpt_str[levOpt];
  }
  std::stringstream filename;
  std::stringstream filenameX;
  std::stringstream filenameY;
  filename<<chain<<"-"<<fit_type<<"-UnitarityTriangle-"<<triangle[n_tr]<<fileAdd[0]<<fileAdd[1]<<"-"<<bool_str[RCreweight]<<"RC"<<fileAdd_truth<<".pdf";
  // filenameX<<chain<<"-"<<fit_type<<"-rho-"<<triangle[n_tr]<<"-"<<bool_str[RCreweight]<<"RC.pdf";
  // filenameY<<chain<<"-"<<fit_type<<"-eta-"<<triangle[n_tr]<<"-"<<bool_str[RCreweight]<<"RC.pdf";
  c->Print(filename.str().c_str());

  // h_rho->Draw("hist");
  // TLine* lin_rho = new TLine(peakpos_proj[0],0.,peakpos_proj[0],peakheight_proj[0]);
  // lin_rho->SetLineStyle(8);
  // lin_rho->SetLineWidth(2);
  // lin_rho->Draw("same");
  // c->Print(filenameX.str().c_str());

  // h_eta->Draw("hist");

  // TLine* lin_eta = new TLine(peakpos_proj[1],0.,peakpos_proj[1],peakheight_proj[1]);
  // lin_eta->SetLineStyle(8);
  // lin_eta->SetLineWidth(2);
  // lin_eta->Draw("same");
  // c->Print(filenameY.str().c_str());
  
}