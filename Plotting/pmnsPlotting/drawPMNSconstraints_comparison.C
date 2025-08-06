// 0 -> draw modulus, 1-> draw real part, 2->draw imaginary part
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <TH1.h>

#include <iostream>
#include <sstream>

std::pair<double, double> findBinEdges(TH1D* hist) {
  int nbins = hist->GetNbinsX();
  double max = hist->GetMaximum();
  std::cout<<"Maximum: "<<max<<std::endl;

  int firstBin = 0, lastBin = 0;
  double firstEdge;
  double lastEdge;

  // Find the first non-empty bin
  for (int i = 1; i <= nbins; i++) { // Bins start from 1 to nbins in ROOT
      if (hist->GetBinContent(i) > max/100) {
          firstBin = i;
          break;
      }
  }

  // Find the last non-empty bin
  for (int i = nbins; i >= 1; i--) {
      if (hist->GetBinContent(i) > max/100) {
          lastBin = i;
          break;
      }
  }

  // Get bin edges
  if (firstBin > 0 && lastBin > 0) {
      firstEdge = hist->GetBinLowEdge(firstBin);
      lastEdge = hist->GetBinLowEdge(lastBin) + hist->GetBinWidth(lastBin);
      
      std::cout << "First non-empty bin lower edge: " << firstEdge << ", Content: " << hist->GetBinContent(firstBin) << std::endl;
      std::cout << "Last non-empty bin upper edge: " << lastEdge << ", Content: " << hist->GetBinContent(lastBin) << std::endl;
  } else {
      std::cout << "No non-empty bins found!" << std::endl;
  }

  return {firstEdge, lastEdge};
}

// Function to insert a space between "Asimov" and "X22"
std::string insertSpaceInAsimov(const std::string& fit_type) {
  std::string modified_fit_type = fit_type;
  if (fit_type.find("Asimov") != std::string::npos) {
      modified_fit_type.replace(6, 0, " ");
  }
  return modified_fit_type;
}

void drawPMNSconstraints_simple_comparison(int m, int n, std::string fit_type_1 = "AsimovA22", std::string fit_type_2 = "OA2021", bool RCreweight = true, int drawOpt=0, int legOpt=0){

  std::string chain = "T2K";

  bool drawPrior = false;

  bool rcComparison = false;
  bool oaComparison = false;
  bool asimovComparison = false;
  bool asimovOAComparison = false;
  bool fdNDComparison = false;
  //int drawOpt = 0; 

  if (fit_type_1 == fit_type_2){
    rcComparison = true;
  }
  else if (fit_type_1.find("OA") != std::string::npos && fit_type_2.find("OA") != std::string::npos){
    oaComparison = true;
  }
  else if (fit_type_1.find("Asimov") != std::string::npos && fit_type_2.find("Asimov") != std::string::npos){
    asimovComparison = true;
  }
  else if ((fit_type_1.find("Asimov") != std::string::npos && fit_type_2.find("OA") != std::string::npos) || (fit_type_1.find("OA") != std::string::npos && fit_type_2.find("Asimov") != std::string::npos)){
    asimovOAComparison = true;
  }
  else if (fit_type_1.find("FD") != std::string::npos && fit_type_2.find("FD") != std::string::npos){
    fdNDComparison = true;
  }

  std::vector<std::string> fit_types = {fit_type_1, fit_type_2};
  
  // std::string years = "2020-2021";
  // std::string asimovs = "A22-B22";
  // std::string asimovOA = "Asimov-A22_OA2021";

  float font_size = 0.05;

  int n_bins = 5000;
  std::vector<int> rebin_factors{25,25};

  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};
  std::map<bool, std::string> bool_str_pretty{{false, "w/o"}, {true, "w/"}};

  std::vector<std::string> filenames;
  bool RCreweight_temp = RCreweight;

  int i = 0;

  for (auto fit_type : fit_types){
    if (fit_type != "AsimovA22"){
      rebin_factors[i] = 10;
    }
    if (m == 1) {
      rebin_factors[i] = 5;
      if (n == 3 && RCreweight_temp){
        // n_bins = 5000;
        rebin_factors[i] = 2;
      }
    }
    // if (m == 1 && n == 3 && RCreweight_temp){
    //   n_bins = 5000;
    // }
    // else if (m == 1 && n == 3){
    //   n_bins = 1000;
    // }
    std::stringstream infilename;
    infilename<<"pmns_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight_temp]<<"RC_wPriors_priorCase0_5000bins.root";
    filenames.push_back(infilename.str());
    if (rcComparison){
      std::cout<<"RCreweight_temp: "<<RCreweight_temp<<std::endl;
      RCreweight_temp = !RCreweight_temp;
    }
    i++;
  }
  // if (rcComparison){
  //   filenames = {"pmns_T2K_AsimovA22_wRC_wPriors_priorCase0_200bins.root","pmns_T2K_AsimovA22_woRC_wPriors_priorCase0_200bins.root"};
  // }
  // else if (oaComparison){
  //   filenames = {"pmns_T2K_OA2020_woRC_wPriors_priorCase0_200bins.root","pmns_T2K_OA2021_woRC_wPriors_priorCase0_200bins.root"};
  // }
  // else if (asimovComparison){
  //   filenames = {"pmns_T2K_AsimovA22_woRC_wPriors_priorCase0_200bins.root","pmns_T2K_AsimovB22_woRC_wPriors_priorCase0_200bins.root"};
  // }
  // else if (asimovOAComparison){
  //   filenames = {"pmns_T2K_AsimovA22_woRC_wPriors_priorCase0_200bins.root","pmns_T2K_OA2021_woRC_wPriors_priorCase0_200bins.root"};
  // }
  // else if (fdNDComparison){
  //   filenames = {"pmns_DUNE_FDOnly_woRC_wPriors_priorCase0_200bins.root","pmns_DUNE_FD+ND_woRC_wPriors_priorCase0_200bins.root"};
  // }
  // else{
  //   filenames = {"pmns_NOvA_MaCh3_DataFit_postBANFF_woRC_wPriors_T2KPrior_woSolarConstraint.root","pmns_NOvA_MaCh3_DataFit_postBANFF_woRC_wPriors_NOvAPrior_woSolarConstraint.root","pmns_NOvA_MaCh3_DataFit_postBANFF_woRC_wPriors_haarPrior_woSolarConstraint.root"};
  // }
  //int m = 2;
  //int n = 1;

  TFile* inFile;

  TCanvas* c = new TCanvas("c","c",1200,1100);
  std::cout<<gStyle->GetPadLeftMargin()<<std::endl;
  std::cout<<gStyle->GetPadLeftMargin()<<std::endl;
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.03);
  std::cout<<gStyle->GetPadTopMargin()<<std::endl;
  std::cout<<gStyle->GetPadBottomMargin()<<std::endl;
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  c->UseCurrentStyle();
  c->Draw();
  c->cd();
  TPad* pads[3][3];
  // TH1D* hists[3][3];
  // TH1D* hists_comp[3][3];
  TH1D* hists;

  std::vector<TH1D*> hists_after;
  std::vector<TH1D*> hists_zoomed;

  std::map<int, std::string> flavour{{1, "e"}, {2, "mu"}, {3, "tau"}};
  std::map<int, std::string> draw_opt{{0, ""}, {1, "real_"}, {2, "imag_"}};
  std::map<int, std::map<int, int>> legOpt_map;

  //define legend options (corresponding to different positions--see below) for each of the elements
  if (!drawPrior){
    std::cout<<"Okay"<<std::endl;
    legOpt_map ={ {1, {{1,0},{2,0},{3,1}} },
                  {2, {{1,1},{2,0},{3,0}} },
                  {3, {{1,0},{2,0},{3,0}} }};
  }

  legOpt = legOpt_map[m][n];

  Int_t t2k_color = TColor::GetFreeColorIndex();
  TColor *c_t2k = new TColor(t2k_color, 94/255.0, 53/255.0, 151/255.0, "", 0.7);

  Int_t nova_color = TColor::GetFreeColorIndex();
  TColor *c_nova = new TColor(nova_color, 1/255.0, 150/255.0, 190/255.0, "", 0.7);

  Int_t haar_color = TColor::GetFreeColorIndex();
  TColor *c_haar = new TColor(haar_color, 248/255.0, 136/255.0, 16/255.0, "", 0.7);

  std::vector<Int_t> colors;
  std::vector<std::string> leg_entries;
  // std::vector<int> rebin_factors;

  std::string leg_fit_type_1 = insertSpaceInAsimov(fit_type_1);
  std::string leg_fit_type_2 = insertSpaceInAsimov(fit_type_2);

  colors = {t2k_color,haar_color};
  leg_entries = {leg_fit_type_1,leg_fit_type_2};
  // rebin_factors = {1,1};

  if (rcComparison){
    leg_entries = {bool_str_pretty[RCreweight]+" RC",bool_str_pretty[!RCreweight]+" RC"};
    // if (m == 1 && n == 3){
    //   rebin_factors = {1+RCreweight*1,1+!RCreweight*1};
    // }
  }
  // else{
  //   colors = {t2k_color,nova_color,haar_color};
  //   leg_entries = {"T2K Prior","NOvA Prior","Haar Prior"};
  // }

  for (int factor: rebin_factors){
    std::cout<<"Rebin factor: "<<factor<<std::endl;
  }

  gStyle->SetOptStat(0);

  TLine* lin;
  if(drawOpt==0) lin = new TLine(0.5,0.,0.5,100.);
  else lin = new TLine(0.,0.,0.,100.);
  lin->SetLineStyle(9);
  lin->SetLineColor(kGray+1);

  double y_max = 0.;

  double x_low;
  if (asimovComparison || asimovOAComparison){
    x_low = 0.64;
  }
  else {
    x_low = 0.74;
  }

  //TLegend* leg = new TLegend();
  TLegend* leg;

  if (legOpt == 0){ //Upper Left
    leg = new TLegend(0.16,0.80,0.43,0.92);
  }
  else if (legOpt == 1){ //Upper Right
    leg = new TLegend(x_low,0.80,0.97,0.92);
  }
  else if (legOpt == 2){ //Lower Left
    leg = new TLegend(0.2,0.18,0.43,0.33);
  }
  else if (legOpt == 3){ //Lower Right
    leg = new TLegend(0.74,0.18,0.97,0.33);
  }
  else if (legOpt == 4){ //Lower Middle Left
    leg = new TLegend(0.385,0.14,0.615,0.29);
  }
  else if (legOpt == 5){ //Upper Middle Left
    leg = new TLegend(0.385,0.73,0.615,0.88);
    //leg = new TLegend(0.42,0.73,0.65,0.88);
  }

  bool comp = false;
  
  double left_bin_edge = 1.;
  double right_bin_edge = 0.;

  double left_pad_edge;

  double lower_limit;
  double upper_limit;

  i = 0;

  for (auto iter = filenames.begin(); iter<filenames.end(); iter++){
    
    std::cout<<"Input file: "<<*iter<<std::endl;
    inFile = new TFile((*iter).c_str(),"READ");

    std::stringstream histname;
    if (!drawPrior){
      histname<<"h_"<<draw_opt[drawOpt]<<"u"<<flavour[m]<<n;
    }
    else{
      histname<<"h_"<<draw_opt[drawOpt]<<"u"<<flavour[m]<<n<<"_prior";
    }

    hists = (TH1D*)inFile->Get(histname.str().c_str());
    hists->Rebin(rebin_factors[i]);
    hists->Scale(1./hists->Integral());
    hists->Scale(1./hists->GetBinWidth(1));
    if(drawOpt==0) hists->GetXaxis()->SetRangeUser(0.01,0.99);
    else if(drawOpt==1) hists->GetXaxis()->SetRangeUser(-0.99,0.99);
    else if(drawOpt==2) hists->GetXaxis()->SetRangeUser(-0.19,0.19);
    //hists->GetYaxis()->SetTickSize(0.);
    // hists->SetFillColor(kBlue);
    // hists->SetLineColor(kBlue);
    hists->SetLineColor(colors[i]);
    hists->SetLineWidth(2);

    if (hists->GetMaximum()>y_max){
      y_max = hists->GetMaximum();
      std::cout<<"y_max: "<<y_max<<std::endl;
    }

    std::pair<double, double> bin_edges;
    bin_edges = findBinEdges(hists);

    if (bin_edges.first < left_bin_edge){
      left_bin_edge = bin_edges.first;
    }
    if (bin_edges.second > right_bin_edge){
      right_bin_edge = bin_edges.second;
    }

    std::cout<<"Lower edge: "<<left_bin_edge<<std::endl;
    std::cout<<"Upper edge: "<<right_bin_edge<<std::endl;

    if (1.-right_bin_edge < left_bin_edge){
      left_pad_edge = 0.165;
    }
    else{
      left_pad_edge = 0.465;
    }
    lower_limit = left_bin_edge - (right_bin_edge - left_bin_edge)*0.2;
    upper_limit = right_bin_edge + (right_bin_edge - left_bin_edge)*0.2;

    //hists->GetYaxis()->SetLimits(0.,y_max);

    leg->AddEntry(hists,leg_entries[i].c_str(),"l");

    hists_after.push_back(hists);

    // Draw the zoomed-in histogram in the inset pad
    TH1D* zoomedHist = (TH1D*)hists->Clone("zoomedHist");
    hists_zoomed.push_back(zoomedHist);

    // hists->Draw("hist SAME");

    i++;

  }

  
  for (auto iter = hists_after.begin(); iter<hists_after.end(); iter++){
    // (*iter)->GetXaxis()->SetRangeUser(0.4,0.6);
    (*iter)->GetYaxis()->SetRangeUser(0.,y_max*1.1);

    (*iter)->SetLineWidth(2);

    (*iter)->GetXaxis()->SetLabelSize(font_size);
    (*iter)->GetXaxis()->SetTitleSize(font_size);
    (*iter)->GetYaxis()->SetLabelSize(font_size);
    (*iter)->GetYaxis()->SetTitleSize(font_size);
    (*iter)->GetYaxis()->CenterTitle(true);
    (*iter)->GetYaxis()->SetTitle("Posterior Probability Density");

    (*iter)->GetXaxis()->SetTitleOffset(1.18);
    (*iter)->GetXaxis()->SetLabelOffset(0.01);
    (*iter)->GetYaxis()->SetTitleOffset(1.3);
    (*iter)->GetYaxis()->SetLabelOffset(0.01);

    (*iter)->Draw("hist SAME");
  }

  leg->SetTextSize(font_size*0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");

  gPad->SetTicks();
  //gPad->SetTickx();
  //gPad->SetTicky();
  gStyle->SetLineWidth(2);
  gPad->RedrawAxis();

  TLatex* latex = new TLatex();
  latex->SetNDC();                // Use normalized coordinates
  latex->SetTextColorAlpha(TColor::GetColor("#555555"), 0.7);    //light gray text
  latex->SetTextFont(leg->GetTextFont()+10); // Set font to Helvetica
  latex->SetTextSize(font_size*0.9); // Set text size
  latex->DrawLatex(0.65, 0.965, "assuming unitarity");  //#splitline{assuming}{unitarity}

  // Create an inset pad for the zoomed-in plot
  TPad* insetPad = new TPad("insetPad", "insetPad", left_pad_edge, 0.25, left_pad_edge+0.4, 0.65);
  insetPad->SetFillStyle(4000); // Make the inset pad transparent
  insetPad->SetLeftMargin(0.14); // was 0.17 before scaling to 1./bin_width
  insetPad->Draw();
  insetPad->cd();

  // Possible major axis divisions for the zoomed-in histogram
  std::vector<double> axis_divisions{0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
  std::vector<int> y_axis_divisions{1, 2, 5, 10, 20, 50, 100};

  if (m == 1 ){ //&& (RCreweight || rcComparison)

    int root_divisions_x = 503;
    int root_divisions_y = 503;

    for (auto div : axis_divisions){
      if ((upper_limit - lower_limit) / div < 5){
        int n_divisions_x =  ceil((upper_limit - lower_limit) / div);
        std::cout<<"Width: "<<upper_limit - lower_limit<<", Div: "<<div<<", n_divisions_x: "<<n_divisions_x<<std::endl;
        root_divisions_x = 500 + n_divisions_x;
        // x_axis_divisions_set = true;
        std::cout<<"Setting x-axis divisions to: "<<root_divisions_x<<", major axis division label: "<<div<<std::endl;
        break;
      }
    }
    for (auto div: y_axis_divisions){
      if (y_max / div < 8 ){
        int n_divisions_y =  ceil(y_max / div);
        std::cout<<"Height: "<<y_max<<", Div: "<<div<<", n_divisions_y: "<<n_divisions_y<<std::endl;
        root_divisions_y = 500 + n_divisions_y;
        // y_axis_divisions_set = true;
        std::cout<<"Setting y-axis divisions to: "<<root_divisions_y<<", major axis division label: "<<div<<std::endl;
        break;
      }
    }

    i=0;

    for (auto iter = hists_zoomed.begin(); iter<hists_zoomed.end(); iter++){
      // (*iter)->GetXaxis()->SetMinimum(0.4);
      // (*iter)->GetXaxis()->SetMaximum(0.6);
      // (*iter)->Rebin(rebin_factors[i]);
      // (*iter)->Scale(1./rebin_factors[i]);
      std::cout<<"Rebin factor: "<<rebin_factors[i]<<std::endl;

      std::cout<<"Lower limit: "<<lower_limit<<", Upper limit: "<<upper_limit<<std::endl;
      (*iter)->GetXaxis()->SetLabelOffset(0.015);
      (*iter)->GetXaxis()->SetLabelSize(font_size*1.8);
      (*iter)->GetXaxis()->SetTitleSize(0.);
      (*iter)->GetXaxis()->SetRangeUser(lower_limit,upper_limit);
      (*iter)->GetYaxis()->SetRangeUser(0.,y_max*1.1);
      (*iter)->GetYaxis()->SetLabelOffset(0.025);
      (*iter)->GetYaxis()->SetLabelSize(font_size*1.8);
      (*iter)->GetYaxis()->SetTitleSize(0.);

      (*iter)->GetXaxis()->SetNdivisions(root_divisions_x);
      (*iter)->GetYaxis()->SetNdivisions(root_divisions_y);
      
      std::cout<<"Drawing histogram."<<std::endl;
      // if (i == 0){
      //   (*iter)->Draw("hist");
      // }
      // else{
      //   (*iter)->Draw("hist SAME");
      // }
      (*iter)->Draw("hist SAME");
      i++;
    }
  }

  gPad->SetTicks();
  //gPad->SetTickx();
  //gPad->SetTicky();
  gStyle->SetLineWidth(2);
  gPad->RedrawAxis();

  std::stringstream filename;
  filename<<"Comparison-PMNSconstraints-"<<chain<<"-";
  if (rcComparison){
    filename<<"rcComparison_"<<fit_type_1<<"_";
  }
  else {
    filename<<fit_type_1<<"_"<<fit_type_2<<"_"<<bool_str[RCreweight]<<"RC_";
  }
  filename<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // std::cout<<"Check8"<<std::endl;
  // if (!drawPrior){
  //   filename<<"Comparison-PMNSconstraints-";
  // }
  // else{
  //   filename<<"Comparison-PMNSpriors-";
  // }
  // if (rcComparison){
  //   filename<<fit_type<<"-rcComparison_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // else if (oaComparison){
  //   filename<<"oaComparison_"<<years<<"_woRC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // else if (asimovComparison){
  //   filename<<"asimovComparison_"<<asimovs<<"_woRC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // else if (asimovOAComparison){
  //   filename<<"asimovOAComparison_"<<asimovOA<<"_woRC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // else if (fdNDComparison){
  //   filename<<"FDNDComparison_woRC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // else{
  //   filename<<"priorComparison-woRC_wSolarConstraint_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<".pdf";
  // }
  // std::cout<<"Check9"<<std::endl;
  c->Print(filename.str().c_str());
  // std::cout<<"Check10"<<std::endl;
  
}