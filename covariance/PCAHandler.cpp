#include "covariance/PCAHandler.h"


// ********************************************
PCAHandler::PCAHandler() {
// ********************************************

}

// ********************************************
PCAHandler::~PCAHandler() {
// ********************************************


}

// ********************************************
void PCAHandler::ConstructPCA(TMatrixDSym * covMatrix, const int firstPCAd, const int lastPCAd, const double eigen_thresh, int& _fNumParPCA) {
// ********************************************

  FirstPCAdpar = firstPCAd;
  LastPCAdpar = lastPCAd;
  eigen_threshold = eigen_thresh;
  MACH3LOG_INFO("PCAing parameters {} through {} inclusive", FirstPCAdpar, LastPCAdpar);
  int numunpcadpars = covMatrix->GetNrows()-(LastPCAdpar-FirstPCAdpar+1);

  TMatrixDSym submat(covMatrix->GetSub(FirstPCAdpar,LastPCAdpar,FirstPCAdpar,LastPCAdpar));

  //CW: Calculate how many eigen values this threshold corresponds to
  TMatrixDSymEigen eigen(submat);
  eigen_values.ResizeTo(eigen.GetEigenValues());
  eigen_vectors.ResizeTo(eigen.GetEigenVectors());
  eigen_values = eigen.GetEigenValues();
  eigen_vectors = eigen.GetEigenVectors();
  double sum = 0;
  // Loop over eigen values and sum them up
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    sum += eigen_values(i);
  }
  nKeptPCApars = eigen_values.GetNrows();
  //CW: Now go through again and see how many eigen values correspond to threshold
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    // Get the relative size of the eigen value
    double sig = eigen_values(i)/sum;
    // Check against the threshold
    if (sig < eigen_threshold) {
      nKeptPCApars = i;
      break;
    }
  }
  _fNumParPCA = numunpcadpars+nKeptPCApars;
  MACH3LOG_INFO("Threshold of {} on eigen values relative sum of eigen value ({}) generates {} eigen vectors, plus we have {} unpcad pars, for a total of {}", eigen_threshold, sum, nKeptPCApars, numunpcadpars, _fNumParPCA);

  //DB Create array of correct size so eigen_values can be used in CorrelateSteps
  eigen_values_master = std::vector<double>(_fNumParPCA, 1.0);
  for (int i = FirstPCAdpar; i < FirstPCAdpar+nKeptPCApars; ++i) {eigen_values_master[i] = eigen_values(i-FirstPCAdpar);}

  // Now construct the transfer matrices
  //These matrices will be as big as number of unPCAd pars plus number of eigenvalues kept
  TransferMat.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);
  TransferMatT.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);

  // Get a subset of the eigen vector matrix
  TMatrixD temp(eigen_vectors.GetSub(0, eigen_vectors.GetNrows()-1, 0, nKeptPCApars-1));

  //Make transfer matrix which is two blocks of identity with a block of the PCA transfer matrix in between
  TMatrixD temp2;
  temp2.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);

  //First set the whole thing to 0
  for(int iRow = 0; iRow < covMatrix->GetNrows(); iRow++){
    for(int iCol = 0; iCol < _fNumParPCA; iCol++){
      temp2[iRow][iCol] = 0;
    }
  }
  //Set the first identity block
  if(FirstPCAdpar != 0){
    for(int iRow = 0; iRow < FirstPCAdpar; iRow++){
      temp2[iRow][iRow] = 1;
    }
  }

  //Set the transfer matrix block for the PCAd pars
  temp2.SetSub(FirstPCAdpar,FirstPCAdpar,temp);

  //Set the second identity block
  if(LastPCAdpar != covMatrix->GetNrows()-1){
    for(int iRow = 0;iRow < (covMatrix->GetNrows()-1)-LastPCAdpar; iRow++){
      temp2[LastPCAdpar+1+iRow][FirstPCAdpar+nKeptPCApars+iRow] = 1;
    }
  }

  TransferMat = temp2;
  // Copy the contents
  TransferMatT = TransferMat;
  // And then transpose
  TransferMatT.T();

  #ifdef DEBUG_PCA
  //KS: Let's dump all useful matrices to properly validate PCA
  DebugPCA(sum, temp, submat, covMatrix->GetNrows());
  #endif
}

#ifdef DEBUG_PCA
// ********************************************
//KS: Let's dump all useful matrices to properly validate PCA
void PCAHandler::DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat, int NumPar) {
// ********************************************
  (void)submat;//This is used if DEBUG_PCA==2, this hack is to avoid compiler warnings
  TFile *PCA_Debug = new TFile("Debug_PCA.root", "RECREATE");
  PCA_Debug->cd();

  bool PlotText = true;
  //KS: If we have more than 200 plot becomes unreadable :(
  if(NumPar > 200) PlotText = false;

  TH1D* heigen_values = new TH1D("eigen_values", "Eigen Values", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_cumulative = new TH1D("heigen_cumulative", "heigen_cumulative", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_frac = new TH1D("heigen_fractional", "heigen_fractional", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  heigen_values->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative = 0;
  for(int i = 0; i < eigen_values.GetNrows(); i++)
  {
    heigen_values->SetBinContent(i+1, (eigen_values)(i));
    heigen_cumulative->SetBinContent(i+1, (eigen_values)(i)/sum + Cumulative);
    heigen_frac->SetBinContent(i+1, (eigen_values)(i)/sum);
    Cumulative += (eigen_values)(i)/sum;
  }
  heigen_values->Write("heigen_values");
  eigen_values.Write("eigen_values");
  heigen_cumulative->Write("heigen_values_cumulative");
  heigen_frac->Write("heigen_values_frac");

  TH2D* heigen_vectors = new TH2D(eigen_vectors);
  heigen_vectors->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors->Write("heigen_vectors");
  eigen_vectors.Write("eigen_vectors");

  TH2D* SubsetPCA = new TH2D(temp);
  SubsetPCA->GetXaxis()->SetTitle("Parameter in Normal Base");
  SubsetPCA->GetYaxis()->SetTitle("Parameter in Decomposed Base");

  SubsetPCA->Write("hSubsetPCA");
  temp.Write("SubsetPCA");
  TH2D* hTransferMat = new TH2D(TransferMat);
  hTransferMat->GetXaxis()->SetTitle("Parameter in Normal Base");
  hTransferMat->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  TH2D* hTransferMatT = new TH2D(TransferMatT);

  hTransferMatT->GetXaxis()->SetTitle("Parameter in Decomposed Base");
  hTransferMatT->GetYaxis()->SetTitle("Parameter in Normal Base");

  hTransferMat->Write("hTransferMat");
  TransferMat.Write("TransferMat");
  hTransferMatT->Write("hTransferMatT");
  TransferMatT.Write("TransferMatT");

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 1024, 1024);
  c1->SetBottomMargin(0.1);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetGrid();

  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // Make pretty correlation colors (red to blue)
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  double maxz = 0;
  double minz = 0;

  c1->Print("Debug_PCA.pdf[");
  TLine *EigenLine = new TLine(nKeptPCApars, 0, nKeptPCApars, heigen_cumulative->GetMaximum());
  EigenLine->SetLineColor(kPink);
  EigenLine->SetLineWidth(2);
  EigenLine->SetLineStyle(kSolid);

  TText* text = new TText(0.5, 0.5, Form("Threshold = %g", eigen_threshold));
  text->SetTextFont (43);
  text->SetTextSize (40);

  heigen_values->SetLineColor(kRed);
  heigen_values->SetLineWidth(2);
  heigen_cumulative->SetLineColor(kGreen);
  heigen_cumulative->SetLineWidth(2);
  heigen_frac->SetLineColor(kBlue);
  heigen_frac->SetLineWidth(2);

  c1->SetLogy();
  heigen_values->SetMaximum(heigen_cumulative->GetMaximum()+heigen_cumulative->GetMaximum()*0.4);
  heigen_values->Draw();
  heigen_frac->Draw("SAME");
  heigen_cumulative->Draw("SAME");
  EigenLine->Draw("Same");
  text->DrawTextNDC(0.42, 0.84,Form("Threshold = %g", eigen_threshold));

  TLegend *leg = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg->SetTextSize(0.04);
  leg->AddEntry(heigen_values, "Absolute", "l");
  leg->AddEntry(heigen_frac, "Fractional", "l");
  leg->AddEntry(heigen_cumulative, "Cumulative", "l");

  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("Same");

  c1->Print("Debug_PCA.pdf");
  c1->SetRightMargin(0.15);
  c1->SetLogy(0);
  delete EigenLine;
  delete leg;
  delete text;
  delete heigen_values;
  delete heigen_frac;
  delete heigen_cumulative;

  heigen_vectors->SetMarkerSize(0.2);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");

  TLine *Eigen_Line = new TLine(0, nKeptPCApars, LastPCAdpar-FirstPCAdpar, nKeptPCApars);
  Eigen_Line->SetLineColor(kGreen);
  Eigen_Line->SetLineWidth(2);
  Eigen_Line->SetLineStyle(kDotted);
  Eigen_Line->Draw("SAME");
  c1->Print("Debug_PCA.pdf");
  delete Eigen_Line;

  SubsetPCA->SetMarkerSize(0.2);
  minz = SubsetPCA->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) SubsetPCA->Draw("COLZ TEXT");
  else SubsetPCA->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete SubsetPCA;

  hTransferMat->SetMarkerSize(0.15);
  minz = hTransferMat->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMat->Draw("COLZ TEXT");
  else hTransferMat->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete hTransferMat;

  hTransferMatT->SetMarkerSize(0.15);
  minz = hTransferMatT->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMatT->Draw("COLZ TEXT");
  else hTransferMatT->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete hTransferMatT;


  //KS: Crosscheck against Eigen library
  #if DEBUG_PCA == 2
  Eigen::MatrixXd Submat_Eigen(submat.GetNrows(), submat.GetNcols());

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < submat.GetNrows(); i++)
  {
    for(int j = 0; j < submat.GetNcols(); j++)
    {
      Submat_Eigen(i,j) = (submat)(i,j);
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> EigenSolver;
  EigenSolver.compute(Submat_Eigen);
  Eigen::VectorXd eigen_val = EigenSolver.eigenvalues().real();
  Eigen::MatrixXd eigen_vect = EigenSolver.eigenvectors().real();
  std::vector<std::tuple<double, Eigen::VectorXd>> eigen_vectors_and_values;
  double Sum_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++)
  {
    std::tuple<double, Eigen::VectorXd> vec_and_val(eigen_val[i], eigen_vect.row(i));
    eigen_vectors_and_values.push_back(vec_and_val);
    Sum_Eigen += eigen_val[i];
  }
  std::sort(eigen_vectors_and_values.begin(), eigen_vectors_and_values.end(),
            [&](const std::tuple<double, Eigen::VectorXd>& a, const std::tuple<double, Eigen::VectorXd>& b) -> bool
            { return std::get<0>(a) > std::get<0>(b); } );
  int index = 0;
  for(auto const vect : eigen_vectors_and_values)
  {
    eigen_val(index) = std::get<0>(vect);
    eigen_vect.row(index) = std::get<1>(vect);
    index++;
  }
  TH1D* heigen_values_Eigen = new TH1D("eig_values", "Eigen Values", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_cumulative_Eigen = new TH1D("eig_cumulative", "heigen_cumulative", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_frac_Eigen = new TH1D("eig_fractional", "heigen_fractional", eigen_val.size(), 0.0, eigen_val.size());
  heigen_values_Eigen->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values_Eigen->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++)
  {
    heigen_values_Eigen->SetBinContent(i+1, eigen_val(i));
    heigen_cumulative_Eigen->SetBinContent(i+1, eigen_val(i)/sum + Cumulative_Eigen);
    heigen_frac_Eigen->SetBinContent(i+1, eigen_val(i)/sum);
    Cumulative_Eigen += eigen_val(i)/sum;
  }
  heigen_values_Eigen->SetLineColor(kRed);
  heigen_values_Eigen->SetLineWidth(2);
  heigen_cumulative_Eigen->SetLineColor(kGreen);
  heigen_cumulative_Eigen->SetLineWidth(2);
  heigen_frac_Eigen->SetLineColor(kBlue);
  heigen_frac_Eigen->SetLineWidth(2);

  c1->SetLogy();
  heigen_values_Eigen->SetMaximum(heigen_cumulative_Eigen->GetMaximum()+heigen_cumulative_Eigen->GetMaximum()*0.4);
  heigen_values_Eigen->Draw();
  heigen_cumulative_Eigen->Draw("SAME");
  heigen_frac_Eigen->Draw("SAME");

  TLegend *leg_Eigen = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg_Eigen->SetTextSize(0.04);
  leg_Eigen->AddEntry(heigen_values_Eigen, "Absolute", "l");
  leg_Eigen->AddEntry(heigen_frac_Eigen, "Fractional", "l");
  leg_Eigen->AddEntry(heigen_cumulative_Eigen, "Cumulative", "l");

  leg_Eigen->SetLineColor(0);
  leg_Eigen->SetLineStyle(0);
  leg_Eigen->SetFillColor(0);
  leg_Eigen->SetFillStyle(0);
  leg_Eigen->Draw("Same");

  c1->Print( "Debug_PCA.pdf");
  c1->SetLogy(0);
  delete heigen_values_Eigen;
  delete heigen_cumulative_Eigen;
  delete heigen_frac_Eigen;
  delete leg_Eigen;

  TH2D* heigen_vectors_Eigen = new TH2D("Eigen_Vectors", "Eigen_Vectors", eigen_val.size(), 0.0, eigen_val.size(), eigen_val.size(), 0.0, eigen_val.size());

  for(int i = 0; i < eigen_val.size(); i++)
  {
    for(int j = 0; j < eigen_val.size(); j++)
    {
      //KS: +1 because there is offset in histogram relative to TMatrix
      heigen_vectors_Eigen->SetBinContent(i+1,j+1, eigen_vect(i,j));
    }
  }
  heigen_vectors_Eigen->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors_Eigen->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors_Eigen->SetMarkerSize(0.15);
  minz = heigen_vectors_Eigen->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));

  if(PlotText) heigen_vectors_Eigen->Draw("COLZ TEXT");
  else heigen_vectors_Eigen->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");

  heigen_vectors->SetTitle("ROOT minus Eigen");
  heigen_vectors->Add(heigen_vectors_Eigen, -1.);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete heigen_vectors_Eigen;

  #endif
  delete heigen_vectors;

  c1->Print( "Debug_PCA.pdf]");
  delete c1;
  PCA_Debug->Close();
  delete PCA_Debug;
}
#endif
