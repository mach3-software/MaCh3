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
#include <fstream>
#include <sstream>
#include <cmath>

// Function to read the CSV file into a 3x3 table
std::vector<std::vector<std::string>> readCSV(const std::string& filename) {
  std::vector<std::vector<std::string>> table(3, std::vector<std::string>(3, ""));
  std::ifstream file(filename);
  if (file.is_open()) {
      std::string line;
      int row = 0;
      while (std::getline(file, line) && row < 3) {
          std::stringstream ss(line);
          std::string cell;
          int col = 0;
          while (std::getline(ss, cell, ',') && col < 3) {
              table[row][col] = cell;
              col++;
          }
          row++;
      }
      file.close();
  }
  return table;
}

// Function to write the 3x3 table to a CSV file
void writeCSV(const std::vector<std::vector<std::string>>& table, const std::string& filename) {
  std::ofstream file(filename);
  if (file.is_open()) {
      for (const auto& row : table) {
          for (size_t i = 0; i < row.size(); ++i) {
              file << row[i];
              if (i < row.size() - 1) {
                  file << ",";
              }
          }
          file << "\n";
      }
      file.close();
  } else {
      std::cerr << "Unable to open file: " << filename << std::endl;
  }
}

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

// Function to find optimal rebinning factor
int findRebinFactor(TH1D* hist){
  int n_bins = hist->GetNbinsX();
  std::cout << "Number of bins: " << n_bins << std::endl;

  std::vector<int> rebin_factors = {1, 2, 4, 5, 8, 10, 20, 25, 50};

  std::string hist_name = hist->GetName();

  int n_bins_old = 200;
  if (hist_name == "h_ue3"){
    n_bins_old = 2500;
  }
  else if (hist_name.find("h_ue") != std::string::npos){
    n_bins_old = 1000;
  }

  std::pair<double, double> bin_edges_rebin;
  bin_edges_rebin = findBinEdges(hist);

  TH1D* hist_oldBinning = (TH1D*)hist->Clone("hist_oldBinning");
  hist_oldBinning->Rebin(n_bins/n_bins_old);

  std::pair<double, double> bin_edges_oldBinning;
  bin_edges_oldBinning = findBinEdges(hist_oldBinning);

  // Calculate number of bins above threshold in the old and new histograms
  std::cout<<"Old bin edges: "<<bin_edges_oldBinning.first<<" - "<<bin_edges_oldBinning.second<<std::endl;
  std::cout<<"Old number of bins: "<<n_bins_old<<std::endl;
  std::cout<<"Number of bins between edges for old binning: "<< (bin_edges_oldBinning.second - bin_edges_oldBinning.first)*n_bins_old <<std::endl;

  int rebin_factor = 1;

  for (int i = 0; i < rebin_factors.size(); i++){
    if ( (bin_edges_rebin.second - bin_edges_rebin.first) * n_bins / rebin_factors[i] < 40 ){
      rebin_factor = rebin_factors[i];
      break;
    }
  }

  std::cout<<"Rebin factor: "<<rebin_factor<<std::endl;
  std::cout<<"Number of bins in the new histogram: "<<n_bins/rebin_factor<<std::endl;
  std::cout<<"New bin edges: "<<bin_edges_rebin.first<<" - "<<bin_edges_rebin.second<<std::endl;
  std::cout<<"Number of bins between edges for new binning: "<< (bin_edges_rebin.second - bin_edges_rebin.first)*n_bins/rebin_factor <<std::endl;

  // if (n_bins > 1000){
  //   rebin_factor = 5;
  // }
  // else if (n_bins > 500){
  //   rebin_factor = 4;
  // }
  // else if (n_bins > 200){
  //   rebin_factor = 2;
  // }
  return rebin_factor;
}

double* getInterval1D(TH1D* hist, double &p68, double &p90, double &p95, double &p99 ,double &p3sig){

  TH1D* hCopy = (TH1D*)hist->Clone("hCopy");

  double contlevel1=0.68;
  double contlevel2=0.90;
  double contlevel3=0.954;
  double contlevel4=0.99;
  double contlevel5=0.9973;

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

void drawPMNSconstraints_simple(int m, int n, int drawOpt=0, std::string fit_type = "AsimovA22", bool RCreweight = true, bool drawTruth = true, bool draw68CI = true) {

  std::string chain = "DUNE";

  bool drawPrior = true;

  if (chain == "PDG"){
    drawTruth = false;
    RCreweight = false;
    draw68CI = false;
    fit_type = "generated";
  }

  int n_bins = 5000;
  int rebin_factor = 25;
  if ((fit_type != "AsimovA22" && chain == "T2K") || (chain == "DUNE")){
    rebin_factor = 10;
  }
  if (m == 1 || (n == 3 && chain == "DUNE")) {
    rebin_factor = 5;
    if (n == 3 && RCreweight){
      // n_bins = 5000;
      rebin_factor = 2;
    }
  }

  bool drawPMNSparams = true;

  float font_size = 0.05;

  //int m = 2;
  //int n = 1;

  int legOpt = 0;
  //int drawOpt = 0;

  std::map<bool, std::string> bool_str{{false, "wo"}, {true, "w"}};

  TFile* inFile;
  std::stringstream infilename;
  std::stringstream infilenameZoomed;
  // UPDATE: Change to dynamic names wrt w/ w/o RC
  if (chain == "NOvA") {
    //inFile = new TFile("pmns_NOvA_MaCh3_DataFit_postBANFF_woRC_wPriors_woReweighting_100mSteps.root","READ");
    // inFile = new TFile("pmns_NOvA_MaCh3_DataFit_postBANFF_woRC_wPriors_T2KPrior_wSolarConstraint_woSinDelCPprior_fullChain.root","READ");
    infilename<<"pmns_NOvA_MaCh3_DataFit_postBANFF_"<<bool_str[RCreweight]<<"RC_wPriors_T2KPrior_wSolarConstraint_woSinDelCPprior_fullChain.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
  }
  else if (chain == "PDG") {
    // inFile = new TFile("pmns_pdg_200bins.root","READ");
    n_bins =  5000;
    drawPrior = false;
    infilename<<"pmns_pdg_"<<n_bins<<"bins_lowOct_dCPsin.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
  }
  else {
    //TFile* inFile = new TFile("pmns_OA2020_wRC.root","READ");
    infilename<<"pmns_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_wPriors_priorCase0_"<<n_bins<<"bins.root";
    // infilenameZoomed<<"pmns_"<<chain<<"_"<<fit_type<<"_"<<bool_str[RCreweight]<<"RC_wPriors_priorCase0_1000bins.root";
    inFile = new TFile(infilename.str().c_str(),"READ");
  }

  // TCanvas* c = new TCanvas("c","c",1350,1200);
  // gStyle->SetPadLeftMargin(.14);
  // // gStyle->SetPadLeftMargin(.1);
  // gStyle->SetPadRightMargin(0.17);
  // // gStyle->SetPadRightMargin(0.06);
  // gStyle->SetPadTopMargin(0.05);
  // gStyle->SetPadBottomMargin(.13);
  // gStyle->SetPalette(kCividis);
  // c->UseCurrentStyle();
  // c->Draw();
  // c->cd();
  TCanvas* c = new TCanvas("c","c",1200,1100);
  std::cout<<gStyle->GetPadLeftMargin()<<std::endl;
  std::cout<<gStyle->GetPadLeftMargin()<<std::endl;
  gStyle->SetPadLeftMargin(0.13); // was 0.17 before scaling to 1./bin_width
  gStyle->SetPadRightMargin(0.03);
  std::cout<<gStyle->GetPadTopMargin()<<std::endl;
  std::cout<<gStyle->GetPadBottomMargin()<<std::endl;
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  c->UseCurrentStyle();
  c->Draw();
  c->cd();
  TPad* pads[3][3];
  TH1D* hists[3][3];
  TH1D* hists_prior[3][3];

  std::map<int, std::string> flavour{{1, "e"}, {2, "mu"}, {3, "tau"}};
  std::map<int, std::string> draw_opt{{0, ""}, {1, "real_"}, {2, "imag_"}};
  std::map<int, std::map<int, int>> legOpt_map;

  //define legend options (corresponding to different positions--see below) for each of the elements
  // if (RCreweight){
  legOpt_map ={ {1, {{1,0},{2,0},{3,1}} },
                {2, {{1,1},{2,0},{3,0}} },
                {3, {{1,0},{2,0},{3,0}} }};
  // }
  // else{
  //   legOpt_map ={ {1, {{1,0},{2,0},{3,3}} },
  //                 {2, {{1,1},{2,0},{3,2}} },
  //                 {3, {{1,0},{2,0},{3,2}} }};
  // }

  int l = m;
  int i = n;

  legOpt = legOpt_map[l][i];

  gStyle->SetOptStat(0);

  Int_t post_color = TColor::GetFreeColorIndex();
  TColor *c_post = new TColor(post_color, 130/255.0, 86/255.0, 192/255.0);

  Int_t post_color_1sig = TColor::GetFreeColorIndex();
  TColor *c_1sig = new TColor(post_color_1sig, 94/255.0, 53/255.0, 151/255.0);

  Int_t prior_color = TColor::GetFreeColorIndex();
  TColor *c_prior = new TColor(prior_color, 248/255.0, 136/255.0, 16/255.0, "", 0.7);

  // for(unsigned int l=1;l<4;l++){
    // for(unsigned int i=1;i<4;i++){
      std::stringstream histname;

      if(drawOpt==0) {
        if(l==1)histname<<"h_ue"<<i;
        if(l==2)histname<<"h_umu"<<i;
        if(l==3)histname<<"h_utau"<<i;
      }
      else if(drawOpt==1) {
        if(l==1)histname<<"h_real_ue"<<i;
        if(l==2)histname<<"h_real_umu"<<i;
        if(l==3)histname<<"h_real_utau"<<i;
      }
      else if(drawOpt==2) {
        if(l==1)histname<<"h_imag_ue"<<i;
        if(l==2)histname<<"h_imag_umu"<<i;
        if(l==3)histname<<"h_imag_utau"<<i;
      }
      std::cout<<"l:"<<l<<" ; i:"<<i<<endl;
      hists[l][i] = (TH1D*)inFile->Get(histname.str().c_str());
      // if(drawOpt!=2) hists[l][i]->Rebin(5);
      // if (l != 1 && RCreweight){
      //   rebin_factor = 25; //findRebinFactor(hists[l][i]);
      // }
      hists[l][i]->Rebin(rebin_factor);
      hists[l][i]->Scale(1./hists[l][i]->Integral());
      if (chain == "NOvA") {
        // if(drawOpt==0) hists[l][i]->GetXaxis()->SetRangeUser(0.01,0.99);
        if(drawOpt==0) hists[l][i]->GetXaxis()->SetRangeUser(0.0,1.0);
        else if(drawOpt==1) hists[l][i]->GetXaxis()->SetRangeUser(-1.0,1.0);
        else if(drawOpt==2) hists[l][i]->GetXaxis()->SetRangeUser(-0.3,0.3);
      }
      else{
        if(drawOpt==0) hists[l][i]->GetXaxis()->SetRangeUser(0.01,0.99);
        else if(drawOpt==1) hists[l][i]->GetXaxis()->SetRangeUser(-0.99,0.99);
        else if(drawOpt==2) hists[l][i]->GetXaxis()->SetRangeUser(-0.19,0.19);
      }
      //hists[l][i]->GetYaxis()->SetTickSize(0.);
      hists[l][i]->SetFillColor(post_color);
      hists[l][i]->SetLineColor(post_color);
      hists[l][i]->SetLineWidth(2);

      hists[l][i]->GetXaxis()->SetLabelSize(font_size);
      hists[l][i]->GetXaxis()->SetTitleSize(font_size);
      hists[l][i]->GetYaxis()->SetLabelSize(font_size);
      hists[l][i]->GetYaxis()->SetTitleSize(font_size);
      hists[l][i]->GetYaxis()->CenterTitle(true);
      hists[l][i]->GetYaxis()->SetTitle("Posterior Probability Density");

      // std::cout<<std::endl;
      // std::cout<<hists[l][i]->GetXaxis()->GetTitleOffset()<<std::endl;
      // std::cout<<hists[l][i]->GetXaxis()->GetLabelOffset()<<std::endl;
      // std::cout<<hists[l][i]->GetYaxis()->GetTitleOffset()<<std::endl;
      // std::cout<<hists[l][i]->GetYaxis()->GetLabelOffset()<<std::endl;
      // std::cout<<std::endl;

      hists[l][i]->GetXaxis()->SetTitleOffset(1.18);
      hists[l][i]->GetXaxis()->SetLabelOffset(0.01);
      hists[l][i]->GetYaxis()->SetTitleOffset(1.3); //1.75 before scaling to 1./bin_width
      hists[l][i]->GetYaxis()->SetLabelOffset(0.01);

      std::stringstream histname_prior;
      if(drawPrior) {
        histname_prior<<histname.str()<<"_prior";
        // std::cout << histname.str() << std::endl;
        hists_prior[l][i] = (TH1D*)inFile->Get(histname_prior.str().c_str());
        // if(drawOpt!=2) hists_prior[l][i]->Rebin(5);
        hists_prior[l][i]->Rebin(rebin_factor);
        hists_prior[l][i]->Scale(1./hists_prior[l][i]->Integral());
        // hists_prior[l][i]->Scale(hists[l][i]->GetMaximum()/hists_prior[l][i]->GetMaximum());
        //hists_prior[l][i]->SetFillColor(kRed);
        hists_prior[l][i]->SetLineColorAlpha(prior_color,1.0);
        //hists_prior[l][i]->SetLineColor(kRed);
        hists_prior[l][i]->SetMarkerColorAlpha(0,0.);
        hists_prior[l][i]->SetLineWidth(2);
      }

    //   if(l==3){
    //     hists[l][i]->GetXaxis()->SetLabelSize(0.06);
    //     hists[l][i]->GetXaxis()->SetTitle("");
    //   }
    // }
  // }

  //TLegend* leg = new TLegend(0.13,0.72,0.38,0.87);
  //TLegend* leg = new TLegend(0.62,0.72,0.87,0.87);
  //TLegend* leg = new TLegend(0.13,0.15,0.38,0.3);

  int n_leg_entries = 1;
  if (drawPrior) n_leg_entries++;
  if (draw68CI) n_leg_entries++;
  if ((fit_type.find("Asimov") != std::string::npos) && drawTruth) n_leg_entries++;

  TLegend* leg;

  if (legOpt == 0){ //Upper Left
    leg = new TLegend(0.16,0.92-0.12*n_leg_entries/2,0.43,0.92); // left edge was 0.2 before scaling to 1./bin_width
  }
  else if (legOpt == 1){ //Upper Right
    leg = new TLegend(0.74,0.92-0.12*n_leg_entries/2,0.97,0.92);
  }
  else if (legOpt == 2){ //Lower Left
    leg = new TLegend(0.16,0.18,0.43,0.18+0.12*n_leg_entries/2);
  }
  else if (legOpt == 3){ //Lower Right
    leg = new TLegend(0.74,0.18,0.97,0.18+0.12*n_leg_entries/2);
  }
  else if (legOpt == 4){ //Lower Middle Left
    leg = new TLegend(0.385,0.14,0.615,0.29);
  }
  else if (legOpt == 5){ //Upper Middle Left
    leg = new TLegend(0.385,0.73,0.615,0.88);
    //leg = new TLegend(0.42,0.73,0.65,0.88);
  }

  leg->SetTextSize(font_size*0.8);  // Change text size
  leg->SetFillStyle(0);  // Removes the fill color (makes it transparent)
  leg->SetBorderSize(0); // Removes the border

  leg->AddEntry(hists[l][i],"Posterior","f");
  double bin_width = hists[l][i]->GetBinWidth(1);
  hists[l][i]->Scale(1./bin_width);
  hists[l][i]->Draw("hist");

  double p68, p90, p95, p99, p3sig;
  double *p = getInterval1D(hists[l][i], p68, p90, p95, p99, p3sig);

  int n_firstBin_p68 = hists[l][i]->FindFirstBinAbove(p68);
  int n_lastBin_p68 = hists[l][i]->FindLastBinAbove(p68);
  double x_firstBin_p68 = hists[l][i]->GetBinLowEdge(n_firstBin_p68);
  double x_lastBin_p68 = hists[l][i]->GetBinLowEdge(n_lastBin_p68) + hists[l][i]->GetBinWidth(n_lastBin_p68);

  std::cout<<"68% interval: "<<x_firstBin_p68<<" - "<<x_lastBin_p68<<std::endl;

  TH1D* hist_p68 = (TH1D*)hists[l][i]->Clone("hist_p68");
  hist_p68->SetFillColor(post_color_1sig);
  hist_p68->SetLineColor(post_color_1sig);
  for (int i = 1; i <= hist_p68->GetNbinsX(); i++){
    if (hist_p68->GetBinContent(i) < p68 ){
      hist_p68->SetBinContent(i, 0);
    }
  }
  if (draw68CI){
    leg->AddEntry(hist_p68,"68% C.I.","f");
    hist_p68->Draw("SAME HIST");
  }

  double val = 0.0;
  double height = 0.0;

  std::string fileAdd;

  if ((fit_type.find("Asimov") != std::string::npos) && drawTruth){
    std::cout<<"Using Asimov data set-->Draw line at true value"<<std::endl;
    TLine* lin;
    TFile* inFile = new TFile("pmns_truth.root","READ");
    std::stringstream tree_name;
    std::stringstream branch_name;
    tree_name << fit_type << "_truth";
    branch_name << draw_opt[drawOpt] << "u" << flavour[m] << n;

    TTree* pmnsTree = (TTree*)inFile->Get(tree_name.str().c_str());
    pmnsTree->SetBranchAddress(branch_name.str().c_str(),&val);
    pmnsTree->GetEntry(0);

    height = hists[l][i]->GetBinContent(hists[l][i]->FindBin(val));
    std::cout<<"height: "<<height<<std::endl;

    lin = new TLine(val,0.,val,height*0.99);
    lin->SetLineStyle(9);
    lin->SetLineWidth(2);
    lin->SetLineColorAlpha(kBlack,0.8);
    // lin->SetLineColor(kBlack);

    lin->Draw("SAME");

    leg->AddEntry(lin,"Truth","l");
    leg->Draw();

    fileAdd = "_"+bool_str[drawTruth]+"Truth";
  }

  if (drawPrior){
    leg->AddEntry(hists_prior[l][i],"Prior","l");
    double bin_width_prior = hists_prior[l][i]->GetBinWidth(1);
    hists_prior[l][i]->Scale(1./bin_width_prior);
    hists_prior[l][i]->Draw("SAME HIST");
    leg->Draw();
  }

  gPad->SetTicks();
  // gPad->SetTickx();
  // gPad->SetTicky();
  gStyle->SetLineWidth(2);
  gPad->RedrawAxis();

  TLatex* latex = new TLatex();
  latex->SetNDC();                // Use normalized coordinates
  latex->SetTextColorAlpha(TColor::GetColor("#555555"), 0.7);    //light gray text
  latex->SetTextFont(leg->GetTextFont()+10); // Set font to Helvetica
  latex->SetTextSize(font_size*0.9); // Set text size
  latex->DrawLatex(0.65, 0.965, "assuming unitarity");  //#splitline{assuming}{unitarity}

  // Open the zoomed-in file with finer binning
  // TFile* inFileZoomed = new TFile(infilenameZoomed.str().c_str(),"READ");

  // Draw the zoomed-in histogram in the inset pad
  TH1D* zoomedHist = (TH1D*)hists[l][i]->Clone("zoomedHist"); //inFileZoomed->Get(histname.str().c_str());
  // zoomedHist->Scale(1./zoomedHist->Integral());
  // zoomedHist->SetFillColor(post_color);
  // zoomedHist->SetLineColor(post_color);
  // zoomedHist->SetLineWidth(2);

  if (m == 1 || n == 3){
    // Possible major axis divisions for the zoomed-in histogram
    std::vector<double> axis_divisions{0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
    std::vector<int> y_axis_divisions{1, 2, 5, 10, 20, 50, 100};

    std::pair<double, double> bin_edges;
    bin_edges = findBinEdges(zoomedHist);

    std::cout<<"Lower edge: "<<bin_edges.first<<std::endl;
    std::cout<<"Upper edge: "<<bin_edges.second<<std::endl;

    double left_pad_edge;

    if (1.-bin_edges.second < bin_edges.first){
      left_pad_edge = 0.165;
    }
    else{
      left_pad_edge = 0.465;
    }
    // Create an inset pad for the zoomed-in plot
    TPad* insetPad = new TPad("insetPad", "insetPad", left_pad_edge, 0.25, left_pad_edge+0.4, 0.65);
    insetPad->SetFillStyle(4000); // Make the inset pad transparent
    insetPad->SetLeftMargin(0.14); // was 0.17 before scaling to 1./bin_width
    insetPad->Draw();
    insetPad->cd();

    double lower_limit = bin_edges.first - (bin_edges.second - bin_edges.first)*0.2;
    double upper_limit = bin_edges.second + (bin_edges.second - bin_edges.first)*0.2;

    double maxZoomed = zoomedHist->GetMaximum();

    double heightZoomed = zoomedHist->GetBinContent(zoomedHist->FindBin(val));

    zoomedHist->GetXaxis()->SetRangeUser(lower_limit,upper_limit);
    zoomedHist->GetXaxis()->SetLabelOffset(0.015);
    zoomedHist->GetXaxis()->SetLabelSize(font_size*1.8);
    zoomedHist->GetXaxis()->SetTitleSize(0.);
    zoomedHist->GetYaxis()->SetRangeUser(0.,maxZoomed*1.1);
    zoomedHist->GetYaxis()->SetLabelOffset(0.025);
    zoomedHist->GetYaxis()->SetLabelSize(font_size*1.8);
    zoomedHist->GetYaxis()->SetTitleSize(0.);

    // // Automatically adjust the axis divisions
    // bool x_axis_divisions_set = false;
    // bool y_axis_divisions_set = false;

    for (auto div : axis_divisions){
      if ((upper_limit - lower_limit) / div < 5){
        int n_divisions_x =  ceil((upper_limit - lower_limit) / div);
        std::cout<<"Width: "<<upper_limit - lower_limit<<", Div: "<<div<<", n_divisions_x: "<<n_divisions_x<<std::endl;
        int root_divisions_x = 500 + n_divisions_x;
        zoomedHist->GetXaxis()->SetNdivisions(root_divisions_x);
        // x_axis_divisions_set = true;
        std::cout<<"Setting x-axis divisions to: "<<root_divisions_x<<", major axis division label: "<<div<<std::endl;
        break;
      }
      // if ((heightZoomed / div < 8) && !y_axis_divisions_set){
      //   int n_divisions_y =  ceil(heightZoomed / div);
      //   std::cout<<"Height: "<<heightZoomed<<", Div: "<<div<<", n_divisions_y: "<<n_divisions_y<<std::endl;
      //   int root_divisions_y = 500 + n_divisions_y;
      //   zoomedHist->GetYaxis()->SetNdivisions(root_divisions_y);
      //   y_axis_divisions_set = true;
      //   std::cout<<"Setting y-axis divisions to: "<<root_divisions_y<<", major axis division label: "<<div<<std::endl;
      // }
      // if (x_axis_divisions_set && y_axis_divisions_set) {
      //   break;
      //   std::cout<<"Both axis-divisions set."<<std::endl;
      // }
    }
    for (auto div: y_axis_divisions){
      if (maxZoomed / div < 8 ){
        int n_divisions_y =  ceil(maxZoomed / div);
        std::cout<<"Height: "<<maxZoomed<<", Div: "<<div<<", n_divisions_y: "<<n_divisions_y<<std::endl;
        int root_divisions_y = 500 + n_divisions_y;
        zoomedHist->GetYaxis()->SetNdivisions(root_divisions_y);
        // y_axis_divisions_set = true;
        std::cout<<"Setting y-axis divisions to: "<<root_divisions_y<<", major axis division label: "<<div<<std::endl;
        break;
      }
    }

    // Automatically adjust the axis labels
    int ndivisions = zoomedHist->GetXaxis()->GetNdivisions();
    std::cout<<"Number of divisions: "<<ndivisions<<std::endl;

    // zoomedHist->GetXaxis()->SetNdivisions(203);

    // int nlabels = zoomedHist->GetXaxis()->GetNlabels();
    // std::cout<<"Number of labels: "<<nlabels<<std::endl;

    // auto labels = zoomedHist->GetXaxis()->GetLabels();
    // if (labels) {
    //   for (int i = 0; i < labels->GetEntries(); i++) {
    //     std::cout << labels->At(i)->GetName() << std::endl;
    //   }
    // }
    // else {
    //   std::cout << "No labels found on the X-axis." << std::endl;
    // }
    // zoomedHist->GetXaxis()->SetNdivisions(ndivisions / 2); // Show only every second label

    double bin_width_zoomed = zoomedHist->GetBinWidth(1);
    // zoomedHist->Scale(1./bin_width_zoomed);
    zoomedHist->Draw("hist");

    TH1D* hist_p68_zoomed = (TH1D*)hist_p68->Clone("hist_p68_zoomed");
    hist_p68_zoomed->SetFillColor(post_color_1sig);
    hist_p68_zoomed->SetLineColor(post_color_1sig);
    for (int i = 1; i <= hist_p68_zoomed->GetNbinsX(); i++){
      if (hist_p68_zoomed->GetBinContent(i) < p68 ){
        hist_p68_zoomed->SetBinContent(i, 0);
      }
    }
    if (draw68CI){
      hist_p68_zoomed->Draw("SAME HIST");
    }

    if ((fit_type.find("Asimov") != std::string::npos) && drawTruth){
      TLine* zoomedLin = new TLine(val, 0., val, heightZoomed*0.99);
      zoomedLin->SetLineStyle(9);
      zoomedLin->SetLineWidth(2);
      zoomedLin->SetLineColorAlpha(kBlack, 0.8);
      zoomedLin->Draw("SAME");
    }

    if (drawPrior){
      TH1D* zoomedHistPrior = (TH1D*)hists_prior[l][i]->Clone("zoomedHistPrior"); //inFileZoomed->Get(histname_prior.str().c_str());
      // zoomedHistPrior->Scale(1./zoomedHistPrior->Integral());
      zoomedHistPrior->SetLineColorAlpha(prior_color,1.0);
      zoomedHistPrior->SetLineWidth(2);

      zoomedHistPrior->GetXaxis()->SetRangeUser(lower_limit,upper_limit);
      double bin_width_prior_zoomed = zoomedHistPrior->GetBinWidth(1);
      // zoomedHistPrior->Scale(1./bin_width_prior_zoomed);
      zoomedHistPrior->Draw("SAME HIST");
    }
  }

  // Read the existing CSV file into a 3x3 table
  std::stringstream csvFilename;
  csvFilename << "pmns_intervals_" << chain << "_" << fit_type << "_" << draw_opt[drawOpt] <<  bool_str[RCreweight] << "RC" << ".csv";
  std::vector<std::vector<std::string>> intervals = readCSV(csvFilename.str().c_str());

  std::cout<<"CSV file: "<<csvFilename.str()<<std::endl;

  // Update the interval in the mth row and nth column
  intervals[m-1][n-1] = std::to_string(x_firstBin_p68) + " - " + std::to_string(x_lastBin_p68);

  // Write the updated table to the CSV file
  writeCSV(intervals, csvFilename.str().c_str());

  gPad->SetTicks();
  // gPad->SetTickx();
  // gPad->SetTicky();
  gStyle->SetLineWidth(2);
  gPad->RedrawAxis();

  //std::cout<<"Checkpoint"<<std::endl;

  int lowest_bin = hists[l][i]->FindFirstBinAbove(hists[l][i]->GetMaximum()/1000);
  int highest_bin = hists[l][i]->FindLastBinAbove(hists[l][i]->GetMaximum()/1000);
  int n_filled_bins = highest_bin - lowest_bin + 1;
  std::cout << "Lowest bin: " << lowest_bin << ", Highest bin: " << highest_bin << ", Number of filled bins: " << n_filled_bins << std::endl;

  std::stringstream filename;
  if (chain == "NOvA") {
    //filename<<"NOvA-DataFit-postBANFF-100mSteps-PMNSconstraints-"<<bool_str[RCreweight]<<"RC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<"_"<<bool_str[drawPrior]<<"Priors.pdf";
    filename<<"NOvA-DataFit-postBANFF-allSteps-PMNSconstraints-"<<bool_str[RCreweight]<<"RC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<"_legOpt"<<legOpt<<"_"<<bool_str[drawPrior]<<"Priors.pdf";
  }
  else {
    filename<<chain<<"-"<<fit_type<<"-PMNSconstraints-"<<bool_str[RCreweight]<<"RC_"<<draw_opt[drawOpt]<<flavour[m]<<n<<fileAdd<<"_legOpt"<<legOpt<<"_"<<bool_str[drawPrior]<<"Priors.pdf";
  }
  c->Print(filename.str().c_str());

  // std::vector<std::string> hist_names_pmns{"h_s2th12", "h_s2th13", "h_s2th23", "h_dcp"};

  // for (int i = 0; i < 4; i++){
  //   TH1D* hist_pmns = (TH1D*)inFile->Get(hist_names_pmns[i].c_str());
  //   hist_pmns->Rebin(5);
  //   hist_pmns->Scale(1./hist_pmns->Integral());
  //   hist_pmns->Scale(1./hist_pmns->GetBinWidth(1));
  //   hist_pmns->GetXaxis()->SetLabelSize(font_size);
  //   hist_pmns->GetXaxis()->SetTitleSize(font_size);
  //   hist_pmns->GetYaxis()->SetLabelSize(font_size);
  //   hist_pmns->GetYaxis()->SetTitleSize(font_size);
  //   hist_pmns->GetYaxis()->CenterTitle(true);
  //   hist_pmns->GetYaxis()->SetTitle("Posterior Probability Density");

  //   hist_pmns->GetXaxis()->SetTitleOffset(1.18);
  //   hist_pmns->GetXaxis()->SetLabelOffset(0.01);
  //   hist_pmns->GetYaxis()->SetTitleOffset(1.3); //1.75 before scaling to 1./bin_width
  //   hist_pmns->GetYaxis()->SetLabelOffset(0.01);

  //   if (hist_names_pmns[i] == "h_s2th13"){
  //     hist_pmns->GetXaxis()->SetRangeUser(0.01,0.03);
  //     hist_pmns->GetXaxis()->SetNdivisions(505);
  //   }

  //   hist_pmns->SetFillColor(post_color);
  //   hist_pmns->SetLineColor(post_color);

  //   TCanvas* c_pmns = new TCanvas("c_pmns","c_pmns",1200,1100);
  //   gStyle->SetPadLeftMargin(0.13); // was 0.17 before scaling to 1./bin_width
  //   gStyle->SetPadRightMargin(0.03);
  //   gStyle->SetPadTopMargin(0.03);
  //   gStyle->SetPadBottomMargin(0.14);
  //   c_pmns->UseCurrentStyle();
  //   c_pmns->Draw();
  //   c_pmns->cd();

  //   hist_pmns->Draw("hist");

  //   gPad->SetTicks();
  //   // gPad->SetTickx();
  //   // gPad->SetTicky();
  //   gStyle->SetLineWidth(2);
  //   gPad->RedrawAxis();

  //   c_pmns->Print(("pmns_"+hist_names_pmns[i]+".pdf").c_str());
  // }

}