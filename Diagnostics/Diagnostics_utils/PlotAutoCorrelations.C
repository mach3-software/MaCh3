#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>

void PlotAutoCorrelations(TString FileName) {

  TFile* File = new TFile(FileName.Data(),"READ");
  TDirectory* AutoCorrDir = (TDirectory*)File->Get("Auto_corr");
  TIter next(AutoCorrDir->GetListOfKeys());

  TString HistName;
  TH1D* Hist;

  TCanvas* Canvas = new TCanvas("Canv","");
  Canvas->Print("AutoCorr.pdf[");

  TObject* obj;
  while (obj = next()) {
    HistName = obj->GetName();
    std::cout << HistName << std::endl;
    Hist = (TH1D*)AutoCorrDir->Get(HistName);
    Hist->Draw();
    Canvas->Print("AutoCorr.pdf");
  }
  Canvas->Print("AutoCorr.pdf]");
  
  /*
  for (int i=0;i<1;i++) {
    std::cout << KeyList[i] << std::endl;
  }
  */
}
