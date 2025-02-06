#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#include "TList.h"
#include "TObjArray.h"
#pragma GCC diagnostic pop

#include "samplePDF/HistogramUtils.h"

// **************************************************
//KS: ROOT changes something with binning when moving from ROOT 5 to ROOT 6. If you open ROOT5 produced file with ROOT6 you will be missing 9 last bins
// However if you use ROOT6 and have ROOT6 file exactly the same code will work. Something have changed with how TH2Poly bins are stored in TFile
void CheckTH2PolyFileVersion(TFile *file) {
// **************************************************
  int FileROOTVersion = file->GetVersion();
  int MainFileROOTVersion = FileROOTVersion;

  // Remove last digit from number
  // till only one digit is left
  while (MainFileROOTVersion >= 10)
      MainFileROOTVersion /= 10;

  std::string SystemROOTVersion = std::string(ROOT_RELEASE);
  int MainSystemROOTVersion = SystemROOTVersion.at(0)  - '0';

  if(MainFileROOTVersion != MainSystemROOTVersion)
  {
    MACH3LOG_ERROR("File was produced with: {} ROOT version", FileROOTVersion);
    MACH3LOG_ERROR("Found: {} ROOT version in the system", SystemROOTVersion);
    MACH3LOG_ERROR("For some documentation please visit me");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// **************************************************
//WP: Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly* poly) {
// **************************************************
  double overflow = 0.;
  //TH2Polys have 9 overflow bins
  for(int iOverflow = -1; iOverflow > -10; iOverflow--)
  {
    overflow += poly->GetBinContent(iOverflow);
  }
  double IntegralUn = NoOverflowIntegral(poly) + overflow;

  return IntegralUn;
} // end function

// **************************************************
//WP: Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly) {
// **************************************************
  double integral = 0.;
  for(int i = 1; i < poly->GetNumberOfBins()+1; i++)
  {
    integral += poly->GetBinContent(i);
  }

  return integral;
} // end function

// **************************************************
//WP: Helper function for projecting TH2Poly onto the X axis
TH1D* PolyProjectionX(TObject* poly, std::string TempName, const std::vector<double>& xbins, const bool computeErrors) {
// **************************************************
  TH1D* hProjX = new TH1D((TempName+"_x").c_str(),(TempName+"_x").c_str(), int(xbins.size()-1), &xbins[0]);

  //KS: Temp Histogram to store error, use double as this is thread safe
  std::vector<double> hProjX_Error(hProjX->GetXaxis()->GetNbins() + 1, 0.0);
  double xlow, xup, frac = 0.;

  //loop over bins in the poly
  for (int i = 0; i < static_cast<TH2Poly*>(poly)->GetNumberOfBins(); i++)
  {
    //get bin and its edges
    TH2PolyBin* bin = static_cast<TH2PolyBin*>(static_cast<TH2Poly*>(poly)->GetBins()->At(i));
    xlow = bin->GetXMin();
    xup = bin->GetXMax();

    //Loop over projected bins, find fraction of poly bin in each
    for(int dx = 0; dx < int(xbins.size()); dx++)
    {
      if(xbins[dx+1] <= xlow || xbins[dx] >= xup)
      {
        frac = 0;
      }
      else if(xbins[dx] <= xlow && xbins[dx+1] >= xup)
      {
        frac = 1;
      }
      else if(xbins[dx] <= xlow && xbins[dx+1] <= xup)
      {
        frac = (xbins[dx+1]-xlow)/(xup-xlow);
      }
      else if(xbins[dx] >= xlow && xbins[dx+1] >= xup)
      {
        frac = (xup-xbins[dx])/(xup-xlow);
      }
      else if(xbins[dx] >= xlow && xbins[dx+1] <= xup)
      {
        frac = (xbins[dx+1]-xbins[dx])/(xup-xlow);
      }
      else
      {
        frac = 0;
      }
      hProjX->SetBinContent(dx+1,hProjX->GetBinContent(dx+1)+frac*bin->GetContent());
      //KS: Follow ROOT implementation and sum up the variance
      if(computeErrors)
      {
        //KS: TH2PolyBin doesn't have GetError so we have to use TH2Poly,
        //but numbering of GetBinError is different than GetBins...
        const double Temp_Err = frac*static_cast<TH2Poly*>(poly)->GetBinError(i+1)*frac*static_cast<TH2Poly*>(poly)->GetBinError(i+1);
        hProjX_Error[dx+1] += Temp_Err;
      }
    }
  }
  //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
  if(computeErrors)
  {
    for (int i = 1; i <= hProjX->GetXaxis()->GetNbins(); ++i)
    {
      const double Error = std::sqrt(hProjX_Error[i]);
      hProjX->SetBinError(i, Error);
    }
  }
  return hProjX;
} // end project poly X function

// **************************************************
//WP: Helper function for projecting TH2Poly onto the Y axis
TH1D* PolyProjectionY(TObject* poly, std::string TempName, const std::vector<double>& ybins, const bool computeErrors) {
// **************************************************
  TH1D* hProjY = new TH1D((TempName+"_y").c_str(),(TempName+"_y").c_str(),int(ybins.size()-1),&ybins[0]);
  //KS: Temp Histogram to store error, use double as this is thread safe
  std::vector<double> hProjY_Error(hProjY->GetXaxis()->GetNbins() + 1, 0.0);
  double ylow, yup, frac = 0.;

  //loop over bins in the poly
  for (int i = 0; i < static_cast<TH2Poly*>(poly)->GetNumberOfBins(); i++)
  {
    //get bin and its edges
    TH2PolyBin* bin = static_cast<TH2PolyBin*>(static_cast<TH2Poly*>(poly)->GetBins()->At(i));
    ylow = bin->GetYMin();
    yup = bin->GetYMax();

    //Loop over projected bins, find fraction of poly bin in each
    for(int dy = 0; dy < int(ybins.size()); dy++)
    {
      if(ybins[dy+1]<=ylow || ybins[dy] >= yup)
      {
        frac = 0;
      }
      else if(ybins[dy] <= ylow && ybins[dy+1] >= yup)
      {
        frac = 1;
      }
      else if(ybins[dy] <= ylow && ybins[dy+1] <= yup)
      {
        frac = (ybins[dy+1]-ylow)/(yup-ylow);
      }
      else if(ybins[dy] >= ylow && ybins[dy+1] >= yup)
      {
        frac = (yup-ybins[dy])/(yup-ylow);
      }
      else if(ybins[dy] >= ylow && ybins[dy+1] <= yup)
      {
        frac = (ybins[dy+1]-ybins[dy])/(yup-ylow);
      }
      else
      {
        frac = 0;
      }
      hProjY->SetBinContent(dy+1,hProjY->GetBinContent(dy+1)+frac*bin->GetContent());
      //KS: Follow ROOT implementation and sum up the variance
      if(computeErrors)
      {
        //KS: TH2PolyBin doesn't have GetError so we have to use TH2Poly,
        //but numbering of GetBinError is different than GetBins...
        const double Temp_Err = frac*static_cast<TH2Poly*>(poly)->GetBinError(i+1)*frac*static_cast<TH2Poly*>(poly)->GetBinError(i+1);
        hProjY_Error[dy+1] += Temp_Err;
      }
    }
  }
  //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
  if(computeErrors)
  {
    for (int i = 1; i <= hProjY->GetXaxis()->GetNbins(); ++i)
    {
      const double Error = std::sqrt(hProjY_Error[i]);
      hProjY->SetBinError(i, Error);
    }
  }
  return hProjY;
} // end project poly Y function

// ****************
//WP: Normalise a th2poly
TH2Poly* NormalisePoly(TH2Poly *Histogram) {
// ****************
  TH2Poly* HistCopy = static_cast<TH2Poly*>(Histogram->Clone());
  double IntegralWidth = PolyIntegralWidth(HistCopy);
  HistCopy = PolyScaleWidth(HistCopy, IntegralWidth);
  std::string title = std::string(HistCopy->GetName())+"_norm";
  HistCopy->SetNameTitle(title.c_str(), title.c_str());

  return HistCopy;
}

// ****************
// Normalise a TH2Poly
void NormaliseTH2Poly(TH2Poly* Histogram){
// ****************
  const double Integral = NoOverflowIntegral(Histogram);
  for(int j = 1; j < Histogram->GetNumberOfBins()+1; j++)
  {
    Histogram->SetBinContent(j, Histogram->GetBinContent(j)/Integral);
  }
}

// ****************
// Make a ratio histogram
template<class HistType>
HistType* RatioHists(HistType *NumHist, HistType *DenomHist) {
// ****************
  HistType *NumCopy = static_cast<HistType*>(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());
  NumCopy->Divide(DenomHist);

  return NumCopy;
}

// ****************
// Make a ratio th2poly
TH2Poly* RatioPolys(TH2Poly *NumHist, TH2Poly *DenomHist) {
// ****************
  TH2Poly *NumCopy = static_cast<TH2Poly*>(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());

  for(int i = 1; i < NumCopy->GetNumberOfBins()+1; ++i) {
    NumCopy->SetBinContent(i,NumHist->GetBinContent(i)/DenomHist->GetBinContent(i));
  }
  return NumCopy;
}

// ****************
TH2D* ConvertTH2PolyToTH2D(TH2Poly *poly, TH2D *h2dhist) {
// ****************
  double xlow, xup, ylow, yup;
  std::string HistTempName = poly->GetName();

  HistTempName += "_";
  //make the th2d
  TH2D *hist = static_cast<TH2D*>(h2dhist->Clone());
  hist->SetNameTitle(HistTempName.c_str(), HistTempName.c_str());

  for(int ix = 0; ix < hist->GetNbinsX() + 2; ix++) {
    for(int iy = 0; iy < hist->GetNbinsY() + 2; iy++) {
      hist->SetBinContent(ix, iy, 0);
    }
  }
  //Loop over poly bins, find the corresponding th2d and setbincontent!
  for(int i = 0; i< poly->GetNumberOfBins(); i++){
    TH2PolyBin & polybin = static_cast<TH2PolyBin &>(*poly->GetBins()->At(i));
    xlow = polybin.GetXMin();
    xup = polybin.GetXMax();
    ylow = polybin.GetYMin();
    yup = polybin.GetYMax();
    int xbin, ybin;

    xbin = hist->GetXaxis()->FindBin(xlow+(xup-xlow)/2);
    ybin = hist->GetYaxis()->FindBin(ylow+(yup-ylow)/2);

    MACH3LOG_TRACE("Poly bin {}, xlow: {}, xup: {}, ylow: {}, yup: {}. Finding bin for ({}, {}). Found Bin ({}, {}) with content {}. But Poly content: {}",
                i, xlow, xup, ylow, yup, (xlow + (xup - xlow) / 2), (ylow + (yup - ylow) / 2), xbin, ybin, polybin.GetContent(), poly->GetBinContent(i));
    hist->SetBinContent(xbin, ybin, polybin.GetContent());
  }
  return hist;
}
// ****************
TH2Poly* ConvertTH2DtoTH2Poly(TH2D* hist) {
// ****************
  // Make the x axis from the momentum of lepton
  TAxis* xaxis = hist->GetXaxis();
  // Make the y axis from the cos theta of lepton
  TAxis* yaxis = hist->GetYaxis();

  TString histname = hist->GetName();
  // Convert TH2D binning to TH2Poly
  TH2Poly* poly = new TH2Poly();
  poly->SetName(histname);
  poly->SetTitle(histname);

  // Copy axis titles
  poly->GetXaxis()->SetTitle(xaxis->GetTitle());
  poly->GetYaxis()->SetTitle(yaxis->GetTitle());

  double xmax, xmin, ymax, ymin;
  for (int iy = 1; iy <= yaxis->GetNbins(); iy++) {
    ymax = yaxis->GetBinUpEdge(iy);
    ymin = yaxis->GetBinLowEdge(iy);
    for (int ix = 1; ix <= xaxis->GetNbins(); ix++) {
      xmax = xaxis->GetBinUpEdge(ix);
      xmin = xaxis->GetBinLowEdge(ix);
      double binofx[] = {xmin, xmax, xmax, xmin};
      double binofy[] = {ymin, ymin, ymax, ymax};
      poly->AddBin(4, binofx, binofy);
    }
  }

  for (int iy = 1; iy <= yaxis->GetNbins(); iy++) {
    ymax = yaxis->GetBinUpEdge(iy);
    ymin = yaxis->GetBinLowEdge(iy);
    for (int ix = 1; ix <= xaxis->GetNbins(); ix++) {
      xmax = xaxis->GetBinUpEdge(ix);
      xmin = xaxis->GetBinLowEdge(ix);

      // Get the content of the corresponding bin in TH2D and set it in TH2Poly
      int bin = hist->GetBin(ix, iy);
      double content = hist->GetBinContent(bin);
      poly->SetBinContent(poly->FindBin((xmin + xmax) / 2, (ymin + ymax) / 2), content);
    }
  }

  return poly;
}

// ****************
//WP: Scale a TH2Poly and divide by bin width
TH2Poly* PolyScaleWidth(TH2Poly *Histogram, double scale) {
// ****************
  TH2Poly* HistCopy = static_cast<TH2Poly*>(Histogram->Clone());
  double xlow, xup, ylow, yup, area;

  for(int i = 1; i < HistCopy->GetNumberOfBins()+1; i++)
  {
    TH2PolyBin* bin = static_cast<TH2PolyBin*>(Histogram->GetBins()->At(i - 1));
    xlow = bin->GetXMin();
    xup = bin->GetXMax();
    ylow = bin->GetYMin();
    yup = bin->GetYMax();
    area = (xup-xlow)*(yup-ylow);
    HistCopy->SetBinContent(i, Histogram->GetBinContent(i)/(area*scale));
  }
  return HistCopy;
}

// ****************
//WP: Integral of TH2Poly multiplied by bin width
double PolyIntegralWidth(TH2Poly *Histogram) {
// ****************
  double integral = 0;
  double xlow, xup, ylow, yup, area;

  for(int i = 1; i < Histogram->GetNumberOfBins()+1; i++)
  {
    TH2PolyBin* bin = static_cast<TH2PolyBin*>(Histogram->GetBins()->At(i - 1));
    xlow = bin->GetXMin();
    xup = bin->GetXMax();
    ylow = bin->GetYMin();
    yup = bin->GetYMax();
    area = (xup-xlow)*(yup-ylow);
    integral += Histogram->GetBinContent(i)*area;
  }
  return integral;
}

// *********************
TH2Poly* MakePolyHist(const std::string& name, const std::vector<double>& BinArray_x, const std::vector<double>& BinArray_y) {
// *********************
  TH2Poly* poly = new TH2Poly();
  poly->SetName(name.c_str());
  poly->SetTitle(name.c_str());
  double xmax, xmin, ymax, ymin;
  for(unsigned int iy = 0; iy < BinArray_y.size()-1; iy++)
  {
    ymax = BinArray_y[iy+1];
    ymin = BinArray_y[iy];
    for(unsigned int ix = 0; ix < BinArray_x.size()-1; ix++)
    {
      xmax = BinArray_x[ix+1];
      xmin = BinArray_x[ix];
      double binofx[] = {xmin, xmax, xmax, xmin};
      double binofy[] = {ymin, ymin, ymax, ymax};
      poly->AddBin(4,binofx,binofy);
    }
  }
  return poly;
}

// *********************
//KS: Remove fitted TF1 from hist to make comparison easier
void RemoveFitter(TH1D* hist, const std::string& name) {
// *********************
  TList *listOfFunctions = hist->GetListOfFunctions();
  TF1 *fitter = dynamic_cast<TF1*>(listOfFunctions->FindObject(name.c_str()));

  listOfFunctions->Remove(fitter);
  delete fitter;
}

// ****************
// Make Poisson Fluctuation of TH1D hist
void MakeFluctuatedHistogramStandard(TH1D *FluctHist, TH1D* PolyHist, TRandom3* rand){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0);

  for (int i = 1; i <= PolyHist->GetXaxis()->GetNbins(); ++i)
  {
    // Get the posterior predictive bin content
    const double MeanContent = PolyHist->GetBinContent(i);
    // Get a Poisson fluctuation of the content
    const double Random = rand->PoissonD(MeanContent);
    // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
    FluctHist->SetBinContent(i,Random);
  }
}

// ****************
// Make Poisson Fluctuation of TH2Poly hist
void MakeFluctuatedHistogramStandard(TH2Poly *FluctHist, TH2Poly* PolyHist, TRandom3* rand) {
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0, 0.0);

  for (int i = 1; i < FluctHist->GetNumberOfBins()+1; ++i)
  {
    // Get the posterior predictive bin content
    const double MeanContent = PolyHist->GetBinContent(i);
    // Get a Poisson fluctuation of the content
    const double Random = rand->PoissonD(MeanContent);
    // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
    FluctHist->SetBinContent(i,Random);
  }
}

// ****************
// Make Poisson Fluctuation of TH1D hist
void MakeFluctuatedHistogramAlternative(TH1D* FluctHist, TH1D* PolyHist, TRandom3* rand){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0);

  const double evrate = PolyHist->Integral();
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
  const int num = rand->Poisson(evrate);
  #pragma GCC diagnostic pop
  int count = 0;
  while(count < num)
  {
    const double candidate = PolyHist->GetRandom();
    FluctHist->Fill(candidate);
    count++;
  }
}

// ****************
//KS: ROOT developers were too lazy do develop getRanom2 for TH2Poly, this implementation is based on:
// https://root.cern.ch/doc/master/classTH2.html#a883f419e1f6899f9c4255b458d2afe2e
int GetRandomPoly2(const TH2Poly* PolyHist, TRandom3* rand){
// ****************
  const int nbins = PolyHist->GetNumberOfBins();
  const double r1 = rand->Rndm();

  double* fIntegral = new double[nbins+2];
  fIntegral[0] = 0.0;

  //KS: This is custom version of ComputeIntegral, once again ROOT was lazy :(
  for (int i = 1; i < nbins+1; ++i)
  {
    fIntegral[i] = 0.0;
    const double content = PolyHist->GetBinContent(i);
    fIntegral[i] += fIntegral[i - 1] + content;
  }
  for (Int_t bin = 1; bin < nbins+1; ++bin)  fIntegral[bin] /= fIntegral[nbins];
  fIntegral[nbins+1] = PolyHist->GetEntries();

  //KS: We just return one rather then X and Y, this way we can use SetBinContent rather than Fill, which is faster
  int iBin = int(TMath::BinarySearch(nbins, fIntegral, r1));
  //KS: Have to increment because TH2Poly has stupid offset arghh
  iBin += 1;

  delete[] fIntegral;
  return iBin;
}

// ****************
// Make Poisson fluctuation of TH2Poly hist
void MakeFluctuatedHistogramAlternative(TH2Poly *FluctHist, TH2Poly* PolyHist, TRandom3* rand){
// ****************
  // Make the Poisson fluctuated hist
  FluctHist->Reset("");
  FluctHist->Fill(0.0, 0.0, 0.0);

  const double evrate = NoOverflowIntegral(PolyHist);
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
  const int num = rand->Poisson(evrate);
  #pragma GCC diagnostic pop
  int count = 0;
  while(count < num)
  {
    const int iBin = GetRandomPoly2(PolyHist, rand);
    FluctHist->SetBinContent(iBin, FluctHist->GetBinContent(iBin) + 1);
    count++;
  }
}

// *************************
TGraphAsymmErrors* MakeAsymGraph(TH1D* sigmaArrayLeft, TH1D* sigmaArrayCentr, TH1D* sigmaArrayRight, const std::string& title) {
// *************************
  TGraphAsymmErrors *var = new TGraphAsymmErrors(sigmaArrayCentr);
  var->SetNameTitle((title).c_str(), (title).c_str());

  // Need to draw TGraphs to set axes labels
  var->Draw("AP");
  var->GetXaxis()->SetTitle(sigmaArrayCentr->GetXaxis()->GetTitle());
  var->GetYaxis()->SetTitle("Number of events/bin");

  for (int m = 0; m < var->GetN(); ++m)
  {
    double xlow = sigmaArrayLeft->GetBinContent(m+1);
    double xhigh = sigmaArrayRight->GetBinContent(m+1);
    double xtemp;

    // Figure out which variation is larger so we set the error correctly
    if (xlow > xhigh)
    {
      xtemp = xlow;
      xlow = xhigh;
      xhigh = xtemp;
    }

    var->SetPointEYhigh(m, xhigh - var->GetY()[m]);
    var->SetPointEYlow(m, var->GetY()[m] - xlow);
  }
  return var;
}

// ****************
//Fast and thread safe fill of violin histogram, it assumes both histograms have the same binning
void FastViolinFill(TH2D* violin, TH1D* hist_1d){
// ****************
  for (int x = 0; x < violin->GetXaxis()->GetNbins(); ++x)
  {
    const int y = violin->GetYaxis()->FindBin(hist_1d->GetBinContent(x+1));
    violin->SetBinContent(x+1, y,  violin->GetBinContent(x+1, y)+1);
  }
}


// ****************
//DB Get the Cherenkov momentum threshold in MeV
double returnCherenkovThresholdMomentum(int PDG) {
// ****************
  constexpr double refractiveIndex = 1.334; //DB From https://github.com/fiTQun/fiTQun/blob/646cf9c8ba3d4f7400bcbbde029d5ca15513a3bf/fiTQun_shared.cc#L757
  double mass =  MaCh3Utils::GetMassFromPDG(PDG)*1e3;
  double momentumThreshold = mass/sqrt(refractiveIndex*refractiveIndex-1.);
  return momentumThreshold;
}

// **************************************************************************
// Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2){
// ***************************************************************************
  constexpr double MLep = 0.10565837;

  // Caluclate muon energy
  double ELep = sqrt((MLep*MLep)+(PLep*PLep));

  double CosTh = (2*EnuTrue*ELep - MLep*MLep - InitialQ2)/(2*EnuTrue*PLep);

  ELep = sqrt((MLep*MLep)+(PUpd*PUpd));

  // Calculate the new Q2
  double Q2Upd = -(MLep*MLep) + 2.0*EnuTrue*(ELep - PUpd*CosTh);

  return Q2Upd - InitialQ2;
}

// **************************************************************************
// Recalculate Enu after Eb shift. Takes in shifted lepton momentum, lepton angle, and binding energy change, and if nu/anu
double CalculateEnu(double PLep, double costh, double Eb, bool neutrino){
// ***************************************************************************
  double mNeff = 0.93956536 - Eb / 1000.;
  double mNoth = 0.93827203;

  if (!neutrino) {
    mNeff = 0.93827203 - Eb / 1000.;
    mNoth = 0.93956536;
  }

  double mLep = 0.10565837;
  double eLep = sqrt(PLep * PLep + mLep * mLep);

  double Enu = (2 * mNeff * eLep - mLep * mLep + mNoth * mNoth - mNeff * mNeff) /(2 * (mNeff - eLep + PLep * costh));

  return Enu;
}
