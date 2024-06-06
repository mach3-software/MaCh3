#include "Structs.h"

#include "TList.h"
#include "TObjArray.h"

namespace MaCh3Utils {

  // *****************************
  // Get the mass of a particle from the PDG
  // In GeV, not MeV!
  double GetMassFromPDG(int PDG) {
  // *****************************

    switch (abs(PDG)) {  
    case 11:
      return 0.511E-3;
      break;
    case 13:
      return 105.658E-3;
      break;
    case 15:
      return 1.77682;
      break;
    case 22:
      return 0.;
      break;
    case 211:
      return 139.57E-3;
      break;
    case 111:
      return 134.98E-3;
      break;
    case 2112:
      return 939.565E-3;
      break;
    case 2212:
      return 938.27E-3;
      break;
    //Oxygen nucleus
    case 1000080160:
      return 14.89926;
      break;
	//eta
	case 221:
	  return 547.862E-3;
	  break;
	  //K^0 (s or l)
	case 311:
	case 130:
	case 310:
	  return 497.611E-3;
	  break;
	case 321:
	  return 493.677E-3;
	  break;
	// Lamda baryon
	case 3122:
	  return 1115.683E-3;
	  break;
    case 12:
    case 14:
    case 16:
      return 0.0;
      break;
    default:
      std::cerr << "Haven't got a saved mass for PDG " << PDG << std::endl;
      std::cerr << "Please implement me! " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    } // End switch
    
    std::cerr << "Warning, didn't catch a saved mass" << std::endl;
    return 0;
  }

 
  //DB Anything added here must be of the form 2^X, where X is an integer
  //
  //DB Used to contain which DetIDs are supported
  std::unordered_map<int,int>KnownDetIDsMap({
      {0,1},    //ND
      {1,8},    //FD
      {2,16},   //SK1Rmu
      {3,32},   //Nova
      {4,64},   //Atm SubGeV e-like
      {5,128},  //Atm SubGeV mu-like 
      {6,256},  //Atm MultiGeV e-like
      {7,512},  //Atm MultiGeV mu-like
	  });
	
  int nKnownDetIDs = KnownDetIDsMap.size();

}

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
        std::cerr<<"File was produced with: "<<FileROOTVersion<<" ROOT version"<<std::endl;
        std::cerr<<"Found: "<<SystemROOTVersion<<" ROOT version in the system"<<std::endl;
        std::cerr<<"For some docuemntation please visit me"<<std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
    }
}

// **************************************************
//WP: Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly* poly) {
// **************************************************

  double overflow = 0;
  //TH2Polys have 9 overflow bins
  for(int iOverflow = -1; iOverflow > -10; iOverflow--)
    {
      overflow+=poly->GetBinContent(iOverflow);
    }
  double IntegralUn = NoOverflowIntegral(poly) + overflow;

  return IntegralUn;

} // end function

// **************************************************
//WP: Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly) {
// **************************************************

  double integral = 0;

  for(int i=1; i < poly->GetNumberOfBins()+1; i++)
  {
    integral += poly->GetBinContent(i);
  }

  return integral;

} // end function

// **************************************************
//WP: Helper function for projecting TH2Poly onto the X axis
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins, bool computeErrors) {
// **************************************************

  TH1D* hProjX = new TH1D((TempName+"_x").c_str(),(TempName+"_x").c_str(),xbins.size()-1,&xbins[0]);

  //KS: Temp Histogram to store error, use double as this is thread safe
  double *hProjX_Error = new double[hProjX->GetXaxis()->GetNbins()+1];
  for (int i = 0; i <= hProjX->GetXaxis()->GetNbins(); ++i) {hProjX_Error[i] = 0;}
  double xlow, xup, frac=0;

  //loop over bins in the poly
  for(int i = 0; i<((TH2Poly*)poly)->GetNumberOfBins(); i++)
  {
    //get bin and its edges
    TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i);
    xlow = bin->GetXMin();
    xup = bin->GetXMax();

    //Loop over projected bins, find fraction of poly bin in each
    for(int dx = 0; dx < int(xbins.size()); dx++)
    {
      if(xbins[dx+1] <= xlow || xbins[dx] >= xup)
      {
        frac = 0;
      }
      else if(xbins[dx]<=xlow && xbins[dx+1] >= xup)
      {
        frac = 1;
      }
      else if(xbins[dx]<=xlow && xbins[dx+1] <= xup)
      {
        frac = (xbins[dx+1]-xlow)/(xup-xlow);
      }
      else if(xbins[dx]>=xlow && xbins[dx+1] >= xup)
      {
        frac = (xup-xbins[dx])/(xup-xlow);
      }
      else if(xbins[dx]>=xlow && xbins[dx+1] <= xup)
      {
        frac = (xbins[dx+1]-xbins[dx])/(xup-xlow);
      }
      else
      {
        frac=0;
      }
      hProjX->SetBinContent(dx+1,hProjX->GetBinContent(dx+1)+frac*bin->GetContent());
      //KS: Follow ROOT implementation and sum up the variance
      if(computeErrors)
      {
          //KS: TH2PolyBin doesn't have GetError so we have to use TH2Poly,
          //but numbering of GetBinError is differnt than GetBins...
          double Temp_Err = frac*((TH2Poly*)poly)->GetBinError(i+1) * frac*((TH2Poly*)poly)->GetBinError(i+1);
          hProjX_Error[dx+1] += Temp_Err;
      }
    }
  }
  //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
  if(computeErrors)
  {
    for (int i = 1; i <= hProjX->GetXaxis()->GetNbins(); ++i)
    {
      double Error = TMath::Sqrt(hProjX_Error[i]);
      hProjX->SetBinError(i, Error);
    }
  }
  delete[] hProjX_Error;
  return hProjX;
} // end project poly X function

// **************************************************
//WP: Helper function for projecting TH2Poly onto the Y axis
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins, bool computeErrors) {
// **************************************************

  TH1D* hProjY = new TH1D((TempName+"_y").c_str(),(TempName+"_y").c_str(),ybins.size()-1,&ybins[0]);
  //KS: Temp Histogram to store error, use double as this is thread safe
  double *hProjY_Error = new double[hProjY->GetXaxis()->GetNbins()+1];
  for (int i = 0; i <= hProjY->GetXaxis()->GetNbins(); ++i) {hProjY_Error[i] = 0;}

  double ylow, yup, frac=0;

  //loop over bins in the poly
  for(int i = 0; i < ((TH2Poly*)poly)->GetNumberOfBins(); i++)
  {
    //get bin and its edges
    TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i);
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
        frac=0;
      }
      hProjY->SetBinContent(dy+1,hProjY->GetBinContent(dy+1)+frac*bin->GetContent());
      //KS: Follow ROOT implementation and sum up the variance
      if(computeErrors)
      {
        //KS: TH2PolyBin doesn't have GetError so we have to use TH2Poly,
        //but numbering of GetBinError is differnt than GetBins...
        double Temp_Err = frac*((TH2Poly*)poly)->GetBinError(i+1) * frac*((TH2Poly*)poly)->GetBinError(i+1);
        hProjY_Error[dy+1] += Temp_Err;
      }
    }
  }
  //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
  if(computeErrors)
  {
    for (int i = 1; i <= hProjY->GetXaxis()->GetNbins(); ++i)
    {
      double Error = TMath::Sqrt(hProjY_Error[i]);
      hProjY->SetBinError(i, Error);
    }
  }
  delete[] hProjY_Error;
  return hProjY;
} // end project poly Y function

// ****************
//WP: Normalise a th2poly
TH2Poly* NormalisePoly(TH2Poly *Histogram) {
// ****************

  TH2Poly* HistCopy = (TH2Poly*)(Histogram->Clone());
  double IntegralWidth = PolyIntegralWidth(HistCopy);
  HistCopy = PolyScaleWidth(HistCopy, IntegralWidth);
  std::string title = std::string(HistCopy->GetName())+"_norm";
  HistCopy->SetNameTitle(title.c_str(), title.c_str());

  return HistCopy;
}

// ****************
TH2D* ConvertTH2PolyToTH2D(TH2Poly *poly, TH2D *h2dhist) {
// ****************
  double xlow, xup, ylow, yup;
  std::string HistTempName = poly->GetName();

  HistTempName += "_";
  //make the th2d
  TH2D *hist = (TH2D*) h2dhist->Clone();
  hist->SetNameTitle(HistTempName.c_str(), HistTempName.c_str());

  for(int ix = 0; ix < hist->GetNbinsX() + 2; ix++) {
    for(int iy = 0; iy < hist->GetNbinsY() + 2; iy++) {
      hist->SetBinContent(ix, iy, 0);
    }
  }
  //Loop over poly bins, find the corresponding th2d and setbincontent!
  for(int i = 0; i< poly->GetNumberOfBins(); i++){
    TH2PolyBin* polybin = (TH2PolyBin*) (poly->GetBins()->At(i)->Clone());
    xlow = polybin->GetXMin();
    xup = polybin->GetXMax();
    ylow = polybin->GetYMin();
    yup = polybin->GetYMax();
    int xbin, ybin;

    xbin = hist->GetXaxis()->FindBin(xlow+(xup-xlow)/2);
    ybin = hist->GetYaxis()->FindBin(ylow+(yup-ylow)/2);

    //std::cout << "Poly bin " << i << ", xlow: " << xlow << ", xup: " << xup << ", ylow: " << ylow << ", yup: " << yup << ". Finding bin for (" << (xlow+(xup-xlow)/2) << "," << (ylow+(yup-ylow)/2) << ")" << ". Found Bin (" << xbin << "," << ybin << ") with content " << polybin->GetContent() << ". But Poly content: " << poly->GetBinContent(i) << std::endl;
    hist->SetBinContent(xbin, ybin, polybin->GetContent());
  }
  return hist;
}
// ****************
TH2Poly* ConvertTH2DToTH2Poly(TH2D* hist) {
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

  TH2Poly* HistCopy = (TH2Poly*)(Histogram->Clone());
  double xlow, xup, ylow, yup, area;

  for(int i = 1; i < HistCopy->GetNumberOfBins()+1; i++)
  {
    TH2PolyBin* bin = (TH2PolyBin*)HistCopy->GetBins()->At(i-1);
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
    TH2PolyBin* bin = (TH2PolyBin*)Histogram->GetBins()->At(i-1);
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
//KS: Remove fitted TF1 from hist to make comparison easier
void RemoveFitter(TH1D* hist, std::string name) {
  // *********************

  TList *listOfFunctions = hist->GetListOfFunctions();
  TF1 *fitter = dynamic_cast<TF1*>(listOfFunctions->FindObject(name.c_str()));

  listOfFunctions->Remove(fitter);
  delete fitter;
}

// ****************
//DB Get the Cernekov momentum threshold in MeV
double returnCherenkovThresholdMomentum(int PDG) {
// ****************
  double refractiveIndex = 1.334; //DB From https://github.com/fiTQun/fiTQun/blob/646cf9c8ba3d4f7400bcbbde029d5ca15513a3bf/fiTQun_shared.cc#L757
  double mass =  MaCh3Utils::GetMassFromPDG(PDG)*1e3;
  double momentumThreshold = mass/sqrt(refractiveIndex*refractiveIndex-1.);
  return momentumThreshold;
}

// **************************************************************************
// Recalculate Q^2 after Eb shift. Takes in shifted lepton momentum, lepton angle, and true neutrino energy
double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2){
// ***************************************************************************

  const double MLep = 0.10565837;

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
