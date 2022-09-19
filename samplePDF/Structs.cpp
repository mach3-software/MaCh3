#include "Structs.h"

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

  TString SKSampleName_toLatexString(TString String) {
    TString ReturnString;

    TObjArray *Arr = String.Tokenize("-");
    TString Str1 = ((TObjString *)(Arr->At(0)))->String();
    TString Str2 = ((TObjString *)(Arr->At(2)))->String();

    if (Str1.CompareTo("numu")==0) {ReturnString = "$\\nu_{\\mu}$";}
    else if (Str1.CompareTo("numub")==0) {ReturnString = "$\\bar{\\nu_{\\mu}}$";}
    else if (Str1.CompareTo("nue")==0) {ReturnString = "$\\nu_{e}$";}
    else if (Str1.CompareTo("nueb")==0) {ReturnString = "$\\bar{\\nu_{e}}$";}
    else {std::cout << "Something's broken in SKSampleName_toLatexString. Given:" << String << std::endl; std::cout << "Exitting.." << std::endl; throw;}

    if (Str1.CompareTo(Str2.Data())) {ReturnString += " sig";}
    return ReturnString;
  }


  // ************************
  // The mean neutrino direction, amazingly enough hard-coded for our pleasure!
  // https://www.t2k.org/nd280/physics/xsec/meetings/2015/jan142015/nudir
  extern const double ND280NuDir[3] = {-0.0128224, -0.0249785, 0.999586};
  // ************************

  // Beam direction at SK
  extern const double SKNuDir[3] = {0.669764, -0.742179, 0.024223};

  //DB Anything added here must be of the form 2^X, where X is an integer
  //
  //DB Used to contain which DetIDs are supported
  std::unordered_map<int,int>KnownDetIDsMap({
      {0,1},    //ND280
      {1,8},    //SK1Re
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
// Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
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
// Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly* poly) {
// **************************************************

  double integral=0;

  for(int i=1; i < poly->GetNumberOfBins()+1; i++)
    {
      integral += poly->GetBinContent(i);
    }

  return integral;

} // end function

// **************************************************
// Helper function for projecting TH2Poly onto the X axis
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins, bool computeErrors) {
// **************************************************

  TH1D* hProjX = new TH1D((TempName+"_x").c_str(),(TempName+"_x").c_str(),xbins.size()-1,&xbins[0]);
  //KS: Temp Histogram to store error
  TH1D* hProjX_Error = new TH1D((TempName+"_x_Err").c_str(),(TempName+"_x_Err").c_str(),xbins.size()-1,&xbins[0]);
  double xlow, xup, frac = 0;

  //loop over bins in the poly
  for(int i = 0; i < ((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      xlow = bin->GetXMin();
      xup = bin->GetXMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dx=0; dx < int(xbins.size()); dx++)
        {
          if(xbins[dx+1]<=xlow || xbins[dx]>=xup)
            {
              frac = 0;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]>=xup)
            {
              frac = 1;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]<=xup)
            {
              frac = (xbins[dx+1]-xlow)/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]>=xup)
            {
              frac = (xup-xbins[dx])/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]<=xup)
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
              //but numbering of GetBinError is differnt than GetBins...
             double Temp_Err = frac*((TH2Poly*)poly)->GetBinError(i+1) * frac*((TH2Poly*)poly)->GetBinError(i+1);
             hProjX_Error->SetBinContent(dx+1, hProjX_Error->GetBinContent(dx+1)+Temp_Err);
          }
        }
    }
    //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
    if(computeErrors)
    {
        for (int i = 1; i <= hProjX_Error->GetXaxis()->GetNbins(); ++i) 
        {
            double Error = TMath::Sqrt(hProjX_Error->GetBinContent(i));
            hProjX->SetBinError(i, Error);
        }   
    }
  delete hProjX_Error;
  return hProjX;
} // end project poly X function

// **************************************************
// Helper function for projecting TH2Poly onto the Y axis
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins, bool computeErrors) {
// **************************************************

  TH1D* hProjY = new TH1D((TempName+"_y").c_str(),(TempName+"_y").c_str(),ybins.size()-1,&ybins[0]);
  TH1D* hProjY_Error = new TH1D((TempName+"_y_Err").c_str(),(TempName+"_y_Err").c_str(),ybins.size()-1,&ybins[0]);
  double ylow, yup, frac=0;

  //loop over bins in the poly
  for(int i = 0; i < ((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      ylow = bin->GetYMin();
      yup = bin->GetYMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dy=0; dy<int(ybins.size()); dy++)
        {
          if(ybins[dy+1]<=ylow || ybins[dy]>=yup)
            {
              frac=0;
            }
          else if(ybins[dy]<=ylow && ybins[dy+1]>=yup)
            {
              frac = 1;
            }
          else if(ybins[dy]<=ylow && ybins[dy+1]<=yup)
            {
              frac=(ybins[dy+1]-ylow)/(yup-ylow);
            }
          else if(ybins[dy]>=ylow && ybins[dy+1]>=yup)
            {
              frac=(yup-ybins[dy])/(yup-ylow);
            }
          else if(ybins[dy]>=ylow && ybins[dy+1]<=yup)
            {
              frac=(ybins[dy+1]-ybins[dy])/(yup-ylow);
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
             hProjY_Error->SetBinContent(dy+1, hProjY_Error->GetBinContent(dy+1)+Temp_Err);
          }
        }
    }
    //KS: The error is sqrt(summed variance)) https://root.cern.ch/doc/master/TH2_8cxx_source.html#l02266
    if(computeErrors)
    {
        for (int i = 1; i <= hProjY_Error->GetXaxis()->GetNbins(); ++i) 
        {
            double Error = TMath::Sqrt(hProjY_Error->GetBinContent(i));
            hProjY->SetBinError(i, Error);
        }   
    }
  delete hProjY_Error;
  return hProjY;
} // end project poly Y function

// ****************
// Normalise a th2poly
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
// Scale a TH2Poly and divide by bin width
TH2Poly* PolyScaleWidth(TH2Poly *Histogram, double scale) {
// ****************

  TH2Poly* HistCopy = (TH2Poly*)(Histogram->Clone());
  double xlow, xup, ylow, yup, area;

  for(int i=1; i<HistCopy->GetNumberOfBins()+1; i++)
    {
      TH2PolyBin* bin = (TH2PolyBin*)HistCopy->GetBins()->At(i-1)->Clone();
      xlow=bin->GetXMin();
      xup=bin->GetXMax();
      ylow=bin->GetYMin();
      yup=bin->GetYMax();
      area = (xup-xlow)*(yup-ylow);
      HistCopy->SetBinContent(i, Histogram->GetBinContent(i)/(area*scale));
      delete bin;
    }

  return HistCopy;
}

// ****************
// Integral of TH2Poly multiplied by bin width
double PolyIntegralWidth(TH2Poly *Histogram) {
// ****************

  double integral=0;
  double xlow, xup, ylow, yup, area;

  for(int i=1; i<Histogram->GetNumberOfBins()+1; i++)
    {
      TH2PolyBin* bin = (TH2PolyBin*)Histogram->GetBins()->At(i-1)->Clone();
      xlow=bin->GetXMin();
      xup=bin->GetXMax();
      ylow=bin->GetYMin();
      yup=bin->GetYMax();
      area = (xup-xlow)*(yup-ylow);
      integral += Histogram->GetBinContent(i)*area;
      delete bin;
    }

  return integral;
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

/*
//DB Function used to define which mode splines are selected for which MaCh3 modes
int MaCh3Mode_to_SplineMode(int Mode) {
  int returnMode = -1;

  switch (Mode) {
  case kMaCh3_CCQE:
    returnMode = kMaCh3_CCQE;
    break;
  case kMaCh3_2p2h:
    returnMode = kMaCh3_2p2h;
    break;
  case kMaCh3_CC1pi0:
  case kMaCh3_CC1pipm:
    returnMode = kMaCh3_CC1pipm;
    break;
  case kMaCh3_CCcoh:
    returnMode = kMaCh3_CCcoh;
    break;
  case kMaCh3_CCMpi:
    returnMode = kMaCh3_CCMpi;
    break;
  case kMaCh3_CCDIS:
    returnMode = kMaCh3_CCDIS;
    break;
  case kMaCh3_NC1pi0:
    returnMode = kMaCh3_NC1pi0;
    break;
  case kMaCh3_NC1pipm:
    returnMode = kMaCh3_NC1pipm;
    break;
  case kMaCh3_NCcoh:
    returnMode = kMaCh3_NCcoh;
    break;
  case kMaCh3_NCoth:
  case kMaCh3_NCMpi:
  case kMaCh3_NCDIS:
    returnMode = kMaCh3_NCoth;
    break;
  case kMaCh3_NC1gam:
    returnMode = kMaCh3_NC1gam;
    break;
  case kMaCh3_CCMisc:
    returnMode = kMaCh3_CCMisc;
    break;
  //DB NewMaCh3Mode: Add case for which spline mode applies to new MaCh3 mode. If new mode doesn't have its own splines, indicate which mode it should piggy-back off
  default:
    std::cerr << "Mode " << Mode << " not found - Quitting!" << std::endl;
    throw;
  }

  return returnMode;
}


//DB Function used to define which MaCh3 modes are Neutral current
bool isMaCh3ModeNC(int Mode) {
  bool isNC = false;

  switch (Mode) {
  case kMaCh3_CCQE:
  case kMaCh3_2p2h:
  case kMaCh3_CC1pi0:
  case kMaCh3_CC1pipm:
  case kMaCh3_CCcoh:
  case kMaCh3_CCMpi:
  case kMaCh3_CCDIS:
  case kMaCh3_CCMisc:
    isNC = false;
    break;
  case kMaCh3_NC1pi0:
  case kMaCh3_NC1pipm:
  case kMaCh3_NCcoh:
  case kMaCh3_NCoth:
  case kMaCh3_NC1gam:
  case kMaCh3_NCMpi:
  case kMaCh3_NCDIS:
    isNC = true;
    break;
  //DB NewMaCh3Mode: Add case for whether new MaCh3 mode is neutral current or charged current
  default:
    std::cerr << "Mode " << Mode << " not found - Quitting!" << std::endl;
    throw;
  }

  return isNC;
}

//DB Function used to give by-mode colours
int MaCh3ModeColor(int Mode) {
  int Colour = -1;
  
  switch (Mode) {
  case kMaCh3_CCQE:
    Colour = kCyan;
    break;
  case kMaCh3_2p2h:
    Colour = kMagenta+2;
    break;
  case kMaCh3_CC1pi0:
    Colour = kTeal;
    break;
  case kMaCh3_CC1pipm:
    Colour = kCyan-8;
    break;
  case kMaCh3_CCcoh:
    Colour = kTeal-6;
    break;
  case kMaCh3_CCMpi:
    Colour = kGreen+1;
    break;
  case kMaCh3_CCDIS:
    Colour = kSpring+8;
    break;
  case kMaCh3_CCMisc:
    Colour = kGray+2;
    break;
  case kMaCh3_NC1pi0:
    Colour = kYellow-7;
    break;
  case kMaCh3_NC1pipm:
    Colour = kOrange+1;
    break;
  case kMaCh3_NCcoh:
    Colour = kRed+1;
    break;
  case kMaCh3_NCoth:
    Colour = kPink+7;
    break;
  case kMaCh3_NC1gam:
    Colour = kViolet+7;
    break;
  case kMaCh3_NCMpi:
    Colour = kBlue;
    break;
  case kMaCh3_NCDIS:
    Colour = kGreen-2;
    break;
  //DB NewMaCh3Mode: Add case for the new MaCh3 mode's colour in THStack hists
  default:
    std::cerr << "Mode " << Mode << " not found - Quitting!" << std::endl;
    throw;
  }

  return Colour;
}
*/
