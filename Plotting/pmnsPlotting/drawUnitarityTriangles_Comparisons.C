#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
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

void plotIntervals2D(TH2D *hist, int levOpt, int color_index = 632)//, TH2D &h68, TH2D &h90)
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
  }
  else {
    h_tr_cont_1sig->SetContour(1,tpp68);
    h_tr_cont_2sig->SetContour(1,tpp95);
    h_tr_cont_3sig->SetContour(1,tpp3sig);
  }

  // h_tr_cont_1sig->Rebin2D(2);
  // h_tr_cont_2sig->Rebin2D(2);
  // h_tr_cont_3sig->Rebin2D(2);
  
  h_tr_cont_1sig->SetLineColor(color_index);
  h_tr_cont_2sig->SetLineColor(color_index);
  h_tr_cont_3sig->SetLineColor(color_index);
  h_tr_cont_1sig->SetLineWidth(2);
  h_tr_cont_2sig->SetLineWidth(2);
  h_tr_cont_3sig->SetLineWidth(2);
  h_tr_cont_1sig->SetLineStyle(1);
  h_tr_cont_2sig->SetLineStyle(2);
  h_tr_cont_3sig->SetLineStyle(3);

  h_tr_cont_1sig->Draw("cont3 LIST same");
  h_tr_cont_2sig->Draw("cont3 LIST same");
  // h_tr_cont_3sig->Draw("cont3 LIST same");
}

std::vector<double> plotBestFitPoint(TH2D *hist, int color_index = 632){

  std::vector<double> peakpos;

  Int_t MaxBin_2d = hist->GetMaximumBin();
  Int_t x_2d,y_2d,z_2d;
  hist->GetBinXYZ(MaxBin_2d,x_2d,y_2d,z_2d);

  peakpos.push_back(hist->ProjectionX()->GetBinCenter(x_2d));
  peakpos.push_back(hist->ProjectionY()->GetBinCenter(y_2d));

  std::cout<<"2D Peak Position: ("<<peakpos[0]<<","<<peakpos[1]<<")"<<std::endl;

  TMarker* mrkr_bf = new TMarker(peakpos[0],peakpos[1],29);
  mrkr_bf->SetMarkerSize(3);
  mrkr_bf->SetMarkerColor(color_index);

  mrkr_bf->Draw();

  return peakpos;
}

TLegend* plotUnitarityTriangle(std::vector<double> peakpos, TLegend* leg_tr, std::string fit_type, int color_index = 632){
  TLine* lin_base;
  TLine* lin_1;
  TLine* lin_2;
  lin_base = new TLine(0.,0.,1.,0.);
  lin_1 = new TLine(0.,0.,peakpos[0],peakpos[1]);
  lin_2 = new TLine(1.,0.,peakpos[0],peakpos[1]);
  lin_base->SetLineWidth(2);
  lin_1->SetLineWidth(2);
  lin_2->SetLineWidth(2);
  lin_base->SetLineColor(color_index);
  lin_1->SetLineColor(color_index);
  lin_2->SetLineColor(color_index);
  lin_base->Draw("same");
  lin_1->Draw("same");
  lin_2->Draw("same");

  if (fit_type == "AsimovA22"){
    fit_type = "Asimov A22";
  }
  else if (fit_type == "AsimovB22"){
    fit_type = "Asimov B22";
  }

  leg_tr->AddEntry(lin_base,fit_type.c_str(),"l");

  return leg_tr;
}


void drawUnitarityTriangles_Comparisons(int n_tr = 0, std::string fit_type_1 = "AsimovA22", std::string fit_type_2 = "AsimovB22", bool RCreweight = false, bool draw_cont = true, int levOpt=1){

  int n_bins = 200; //100 - This is what was used for T2K TN485 , 200 - This is what was used for DUNE

  std::string chain = "DUNE";
  // int rebin_factor = 3;

  Float_t font_size = 0.05;

  bool rcComparison = true;

  std::vector<std::string> fit_types;
  if (!rcComparison){
    fit_types = {fit_type_1,fit_type_2};
  }
  else{
    fit_types = {fit_type_1,fit_type_1};
  }

  std::map<int, std::string> triangle{{0, "emu"}, {1, "etau"}, {2, "mutau"}, {3, "12"}, {4, "13"}, {5, "23"}};
  std::map<int, std::string> levOpt_str{{0,"-w1sig2sig3sigContours"},{1,"-w68per90per99perContours"}};
  std::map<int, std::string> rc_case{{0,"w/o RC"}, {1, "w/   RC"}};
  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};

  std::vector<std::string> filenames;

  for (auto iter = fit_types.begin(); iter<fit_types.end(); iter++){
    
    std::stringstream filename;
    std::cout<<"Input fit type: "<<*iter<<std::endl;
    std::cout<<"RC reweight: "<<RCreweight<<std::endl;
    if (chain == "T2K") {
      filename<<"unitarity_"<<chain<<"_"<<(*iter).c_str()<<"_"<<bool_str[RCreweight]<<"RC_woPriors_"<<n_bins<<"bins.root";
    }
    else {
      filename<<"unitarity_"<<chain<<"_"<<(*iter).c_str()<<"_"<<bool_str[RCreweight]<<"RC_woPriors_"<<n_bins<<"bins.root";
    }

    filenames.push_back(filename.str());
    if (rcComparison){
      (*iter) = rc_case[RCreweight];
      RCreweight = !RCreweight;
    }
  }

  TFile* inFile;

  TCanvas* c = new TCanvas("c","c",1200,1200);
  gStyle->SetPadLeftMargin(.09);
  // gStyle->SetPadLeftMargin(.1);
  gStyle->SetPadRightMargin(0.03);
  // gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(.11);
  c->UseCurrentStyle();
  c->Draw();
  c->cd();

  TH2D* hists;
  // TH2D* h_tr_rc;

  // TH2D* h_tr_comp;

  std::vector<double> peakpos;
  double peakpos_proj[2];
  int peakheight_proj[2];

  gStyle->SetOptStat(0);

  // std::vector<TH2D*> tr_hists = {h_tr,h_tr_rc};
  std::vector<int> hist_color_indices = {632,600};

  int max_length = std::max(fit_types[0].length(),fit_types[1].length())+1;
  // double left_border;
  // if (n_tr != 5){
  //   left_border = 0.83-max_length*0.012;
  // }
  // else{
  //   left border = 0.82-max_length*0.012;
  // }


  std::cout<<"Max length: "<<max_length<<std::endl;

  TLegend* leg;
  TLegend* leg_tr;
  TLatex* latex;
  if (n_tr != 5){
    // leg = new TLegend(0.13,0.82,0.36,0.97);
    leg = new TLegend(0.12,0.77,0.35,0.95);
    leg_tr = new TLegend(0.80-max_length*0.012,0.83,0.91,0.95);
    latex = new TLatex(0.13,0.16,"assuming unitarity");
  }
  else{
    leg = new TLegend(0.12,0.14,0.35,0.32);
    leg_tr = new TLegend(0.80-max_length*0.012,0.14,0.91,0.26);
    latex = new TLatex(0.13,0.91,"assuming unitarity");
  }

  int k = 0;

  for (auto iter = filenames.begin(); iter<filenames.end(); iter++){
    
    std::cout<<"Input file: "<<*iter<<std::endl;
    inFile = new TFile((*iter).c_str(),"READ");

    std::stringstream histname;
    histname<<"h_tr_"<<triangle[n_tr];

    hists = (TH2D*)inFile->Get(histname.str().c_str());
    // hists->Rebin2D(rebin_factor);

    // hists->Scale(1./h_tr->Integral());

    hists->SetContour(255);
    hists->Smooth(1);
    // hists->GetZaxis()->SetTickSize(0.);
    // hists->GetZaxis()->SetLabelSize(0.);
    // hists->GetZaxis()->SetTitle("Posterior Probability Density  ");
    // hists->GetZaxis()->SetTitleOffset(0.4);
    hists->GetYaxis()->SetTitleOffset(0.8);
    //gPad->SetLogz(1);

    hists->GetXaxis()->SetLabelSize(font_size);
    hists->GetXaxis()->SetTitleSize(font_size);
    hists->GetYaxis()->SetLabelSize(font_size);
    hists->GetYaxis()->SetTitleSize(font_size);

    hists->GetYaxis()->CenterTitle();
    hists->GetXaxis()->CenterTitle();

    plotIntervals2D(hists,levOpt,hist_color_indices[k]);
    peakpos = plotBestFitPoint(hists,hist_color_indices[k]);
    leg_tr = plotUnitarityTriangle(peakpos,leg_tr,fit_types[k],hist_color_indices[k]);
    c->Update();

    k++;
  }

  TMarker* hpd = new TMarker(peakpos[0],peakpos[1],29);
  hpd->SetMarkerSize(3);

  TLine* cont1;
  TLine* cont2;
  TLine* cont3;
  cont1 = new TLine(0.,0.,1.,0.);
  cont2 = new TLine(0.,0.,1.,0.);
  cont3 = new TLine(0.,0.,1.,0.);
  cont1->SetLineWidth(2);
  cont2->SetLineWidth(2);
  cont3->SetLineWidth(2);
  cont1->SetLineStyle(1);
  cont2->SetLineStyle(2);
  cont3->SetLineStyle(3);

  leg->AddEntry(hpd,"H.P.D.","p");
  leg->AddEntry(cont1,"68% C.I.","l");
  leg->AddEntry(cont2,"90% C.I.","l");
  // leg->AddEntry(cont3,"99% C.I.","l");

  // Change legend properties
  leg->SetTextSize(font_size*0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  leg_tr->SetTextSize(font_size*0.8);
  leg_tr->SetFillStyle(0);
  leg_tr->SetBorderSize(0);
  leg_tr->Draw();

  latex->SetNDC();                // Use normalized coordinates
  latex->SetTextColorAlpha(TColor::GetColor("#555555"), 0.7);    //light gray text
  latex->SetTextFont(leg->GetTextFont()+10); // Set font to Helvetica
  latex->SetTextSize(font_size*0.9); // Set text size
  latex->Draw("SAME");  //#splitline{assuming}{unitarity}

  gPad->SetTicks();
  //gPad->SetTickx();
  //gPad->SetTicky();
  gStyle->SetLineWidth(2);
  gPad->RedrawAxis();

  std::string fileAdd;

  if (draw_cont)
  {
    fileAdd = levOpt_str[levOpt];
  }
  std::stringstream filename;
  if (!rcComparison){
    filename<<"Comparison-UnitarityTriangle-"<<triangle[n_tr]<<"-"<<chain<<"-"<<fit_types[0]<<"-"<<fit_types[1]<<"-"<<bool_str[RCreweight]<<"RC"<<fileAdd<<"-"<<n_bins<<"bins.pdf";
  }
  else{
    filename<<"Comparison-UnitarityTriangle-"<<triangle[n_tr]<<"-"<<chain<<"-"<<fit_type_1<<"-rcComparison"<<fileAdd<<"-"<<n_bins<<"bins.pdf";
  }
  c->Print(filename.str().c_str());
  
}
