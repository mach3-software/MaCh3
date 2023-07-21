#include <iostream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>


//KS: This script is used to analyse output form DiagMCMC. This script support comparing up to 4 files, there is easy way to expand it up to five or six,
//Since I doubt anyone would like to comapre more I did it like this
//.x PlotMCMC_Diagnosis.cpp("DiagMCMC_Output1.root")
//.x PlotMCMC_Diagnosis.cpp("DiagMCMC_Output1.root","DiagMCMC_Output2.root", "DiagMCMC_Output3.root", "DiagMCMC_Output4.root")

void PlotMCMC_Diagnosis(TString fname1);
void PlotMCMC_Diagnosis(TString fname1, TString fname2);
void PlotMCMC_Diagnosis(TString fname1, TString fname2, TString fname3);
void PlotMCMC_Diagnosis(TString fname1, TString fname2, TString fname3, TString fname4);
void MakePlot(TString fname1, TString fname2,TString fname3, TString fname4);
void PlotAutoCorr(TString fname1, TString fname2,TString fname3, TString fname4);
double GetMinimumInRange(TH1D *hist, double minRange, double maxRange);

TString DUMMYFILE = "KillMePlease";
bool PlotFlux = false;

void PlotMCMC_Diagnosis(TString fname1) 
{
    MakePlot(fname1, DUMMYFILE, DUMMYFILE, DUMMYFILE); 
    PlotAutoCorr(fname1, DUMMYFILE, DUMMYFILE, DUMMYFILE);
    return 0;
}

void PlotMCMC_Diagnosis(TString fname1, TString fname2) 
{
      MakePlot(fname1, fname2, DUMMYFILE ,DUMMYFILE);
      PlotAutoCorr(fname1, fname2, DUMMYFILE, DUMMYFILE);
      return 0;
}

void PlotMCMC_Diagnosis(TString fname1, TString fname2, TString fname3) 
{
      MakePlot(fname1, fname2, fname3, DUMMYFILE);
      PlotAutoCorr(fname1, fname2, fname3, DUMMYFILE);
      return 0;
}

void PlotMCMC_Diagnosis(TString fname1, TString fname2, TString fname3, TString fname4) 
{
      MakePlot(fname1, fname2, fname3, fname4);
      PlotAutoCorr(fname1, fname2, fname3, fname4);
      return 0;
}

void MakePlot(TString fname1, TString fname2,TString fname3, TString fname4)
{
    TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
    gStyle->SetOptStat(0); //Set 0 to disable statystic box
    
    TKey *key;
    TFile *infile = TFile::Open(fname1.Data());
    
    TFile *infile2 = NULL;
    if(fname2 != DUMMYFILE)infile2 = TFile::Open(fname2.Data());
    TFile *infile3 = NULL;
    if(fname3 != DUMMYFILE) infile3 = TFile::Open(fname3.Data());
    TFile *infile4 = NULL;
    if(fname4 != DUMMYFILE) infile4 = TFile::Open(fname4.Data());
            
        
    TIter next(infile->GetListOfKeys());

    while ((key = (TKey*)next())) { 
            std::string dirname = std::string(key->GetName());  
            //KS: Script will work with LogL and Batched_means, you can comment it if you are interested in it
            if( (dirname == "LogL") || (dirname == "Batched_means") ) continue;
            //KS: Trace wo longer chains is super big, the way to avoid is to plot as png but I don't like png,
            //keepp possiblilty to skip it
            //if( (dirname == "Trace") ) continue;
            infile->cd(dirname.c_str());
            TDirectory *dir;
            dir=gDirectory;
            TIter nextsub(gDirectory->GetListOfKeys());
            c1->Print(Form("%s.pdf[",dirname.c_str()), "pdf");
            TKey *subkey;
            while ((subkey = (TKey*)nextsub()))
            {
                    std::string name = std::string(subkey->GetName());
                    //KS: Skip flux params
                    if ((name.rfind("b_", 0) == 0 ) && !PlotFlux) continue;
                    name = dirname + "/" + name;
                    std::cout<<name<<std::endl;

                if (std::string(subkey->GetClassName()) != "TH1D") continue;

                    TH1D* blarb[4];
                    blarb[0] = (TH1D*)infile->Get(name.c_str())->Clone();
                    //KS: Some fixe pramas can go crazy
                    if(TMath::IsNaN(blarb[0]->GetBinContent(1)) ) continue;
                    
                    //KS: This is unfortunately hardcoded, need to find better way to write this 
                    //blarb[0]->GetListOfFunctions()->ls();
                    delete blarb[0]->GetListOfFunctions()->FindObject("Fitter");
                    blarb[0]->SetLineStyle(kSolid);
                    blarb[0]->SetLineColor(kRed);
                    blarb[0]->Draw();
                    if( dirname == "AccProb") blarb[0]->GetYaxis()->SetRangeUser(0, 1.0);
                    if( name == "AccProb/AcceptanceProbability" ) continue;
                    if(infile2 != NULL)
                    {
                        blarb[1] = (TH1D*)infile2->Get(name.c_str())->Clone();
                        delete blarb[1]->GetListOfFunctions()->FindObject("Fitter");
                        blarb[1]->SetLineStyle(kDashed);
                        blarb[1]->SetLineColor(kBlue);
                        blarb[1]->Draw("same");
                    }
                    if(infile3 != NULL)
                    {
                        blarb[2] = (TH1D*)infile3->Get(name.c_str())->Clone();
                        delete blarb[2]->GetListOfFunctions()->FindObject("Fitter");
                        blarb[2]->SetLineStyle(kDotted );
                        blarb[2]->SetLineColor(kGreen);
                        blarb[2]->Draw("same");
                    }
                    if(infile4 != NULL)
                    {
                        blarb[3] = (TH1D*)infile4->Get(name.c_str())->Clone();
                        delete blarb[3]->GetListOfFunctions()->FindObject("Fitter");
                        blarb[3]->SetLineStyle(kDashDotted );
                        blarb[3]->SetLineColor(kOrange);
                        blarb[3]->Draw("same");
                    }
                
                c1->Print(Form("%s.pdf",dirname.c_str()), "pdf");
            }
            gDirectory->cd("..");
            c1->Print(Form("%s.pdf]",dirname.c_str()), "pdf");
        }

    
    infile->Close();
    if(infile2 != NULL)infile2->Close();
    if(infile3 != NULL)infile3->Close();
    if(infile4 != NULL)infile4->Close();
}

void PlotAutoCorr(TString fname1, TString fname2, TString fname3, TString fname4) 
{
    TString fname[4];
    fname[0] = fname1; fname[1] = fname2; fname[2] = fname3; fname[3] = fname4;
    Color_t PlotColor[4]={kRed, kBlue, kGreen, kOrange};
    
    TFile *infile[4];
    infile[0] = TFile::Open(fname[0].Data());
    //KS" We need to check number of files to loop over in very lazy way
    int Nfiles = 1;
    
    if(fname[1] != DUMMYFILE){ infile[1] = TFile::Open(fname[1].Data()); Nfiles++;}
    if(fname[2] != DUMMYFILE){ infile[2] = TFile::Open(fname[2].Data()); Nfiles++;}
    if(fname[3] != DUMMYFILE){ infile[3] = TFile::Open(fname[3].Data()); Nfiles++;}
    
    TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 800,630);
    gStyle->SetOptStat(0); //Set 0 to disable statystic box
      
    c1->Print("Auto_Corr_PerFile.pdf[", "pdf");
    for(int ik = 0; ik<Nfiles; ik++)
    {
        TIter next(infile[ik]->GetListOfKeys());

        TKey* key;
        while ((key = (TKey*)next())) 
        { 
            std::string dirname = std::string(key->GetName());  

            //KS: Script We are only interested in auto corr
            if( (dirname != "Auto_corr")) continue;
            
            infile[ik]->cd(dirname.c_str());
            TDirectory *dir;
            dir=gDirectory;
            TIter nextsub(gDirectory->GetListOfKeys());

            TKey *subkey;
            bool FirstTime = true;
            while ((subkey = (TKey*)nextsub()))
            {
                    std::string name = std::string(subkey->GetName());
                    //KS: Skip flux params
                    if ((name.rfind("b_", 0) == 0 ) && !PlotFlux) continue;
                    name = dirname + "/" + name;

                if (std::string(subkey->GetClassName()) != "TH1D") continue;
                    std::cout<<name.c_str()<<std::endl;
                    TH1D* blarb = (TH1D*)infile[ik]->Get(name.c_str())->Clone();
                    //KS: Some fixe pramas can go crazy
                    if(TMath::IsNaN(blarb->GetBinContent(1))) continue;
                    //KS: This is unfortunately hardcoded, need to find better way to write this 
                    //blarb[0]->GetListOfFunctions()->ls();
                    delete blarb->GetListOfFunctions()->FindObject("Fitter");
                
                    double MinValue = GetMinimumInRange(blarb, 0, 24000);
                   
                   if(MinValue >=0.80) blarb->SetLineColor(kRed);
                   else if(MinValue >=0.40 && MinValue <0.80) blarb->SetLineColor(kOrange);
                   else if(MinValue >0.20 && MinValue <0.40) blarb->SetLineColor(kYellow);
                   else if(MinValue <=0.20) blarb->SetLineColor(kGreen);
                   blarb->GetXaxis()->UnZoom();
                   
                   if(FirstTime) blarb->SetTitle( Form("Auto_Corr_%s.pdf",fname[ik].Data()) );

                    blarb->SetLineStyle(kDashed);
                    
                    if(FirstTime) blarb->Draw();
                    if(!FirstTime) blarb->Draw("same");
                    FirstTime = false;
            }
            gDirectory->cd("..");
        }
        c1->Print("Auto_Corr_PerFile.pdf", "pdf");
    }
    c1->Print("Auto_Corr_PerFile.pdf]", "pdf");
}

//KS: function which looks for minimum in given range
double GetMinimumInRange(TH1D *hist, double minRange, double maxRange)
{
    double MinVale = 1234567890.;
    TAxis* xaxis = hist->GetXaxis();
    for(int x = 1; x <= hist->GetNbinsX(); x++) 
    {
        if(xaxis->GetBinLowEdge(x)>minRange && xaxis->GetBinUpEdge(x)<maxRange){
            if(MinVale > hist->GetBinContent(x))MinVale = hist->GetBinContent(x);
        }
    }
    return MinVale;
}
