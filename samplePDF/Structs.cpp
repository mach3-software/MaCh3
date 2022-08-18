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


// ********************
// Constructor
xsec2015::xsec2015() {
// ********************
  splMAQE = NULL;
  splCA5 = NULL;
  splMARES = NULL;
  splBGSCLLOWPPI = NULL;
  splBGSCL = NULL;
  splBYDIS = NULL;
  splBYMPI = NULL;
  splAGKYMULT = NULL;
  splFEFABS = NULL;
  splFEFCX = NULL;
  splFEFQE = NULL;
  splFEFINEL = NULL;
  splFEFCXH = NULL;
  splFEFQEH = NULL;
  splPDDC = NULL;
  splPDDO = NULL;
  splMECENULOW = NULL;
  splMECENUHIGH = NULL;
  splMECENUBARLOW = NULL;
  splMECENUBARHIGH = NULL;
  splEBC = NULL;
  splEBO = NULL;
  splSCCV = NULL;
  splSCCA = NULL;
  spl2P2HEDEP_LOWENU = NULL;
  spl2P2HEDEP_HIGHENU = NULL;
  spl2P2HEDEP_LOWENUBAR = NULL;
  spl2P2HEDEP_HIGHENUBAR = NULL;
  splISO_BKG_LOWPPI = NULL;
  splDIS_BY = NULL;
  splMPI_BY = NULL;
  splMPI_AGKY_XSEC = NULL;

  relRPA = __BAD_DOUBLE__;
}

// ************************
// Destructor
xsec2015::~xsec2015() {
// ************************
  if (splMAQE) {
    delete splMAQE;
  }
  if (splCA5) {
    delete splCA5;
  }
  if (splMARES) {
    delete splMARES;
  }
  if (splBGSCLLOWPPI) {
    delete splBGSCLLOWPPI;
  }
  if (splBGSCL) {
    delete splBGSCL;
  }
  if (splBYDIS) {
    delete splBYDIS;
  }
  if (splBYMPI) {
    delete splBYMPI;
  }
  if (splAGKYMULT) {
    delete splAGKYMULT;
  }
  if (splFEFABS) {
    delete splFEFABS;
  }
  if (splFEFCX) {
    delete splFEFCX;
  }
  if (splFEFQE) {
    delete splFEFQE;
  }
  if (splFEFINEL) {
    delete splFEFINEL;
  }
  if (splFEFCXH) {
    delete splFEFCXH;
  }
  if (splFEFQEH) {
    delete splFEFQEH;
  }
  if (splPDDC) {
    delete splPDDC;
  }
  if (splPDDO) {
    delete splPDDO;
  }
  if (splMECENULOW) {
    delete splMECENULOW;
  }
  if (splMECENUHIGH) {
    delete splMECENUHIGH;
  }
  if (splMECENUBARLOW) {
    delete splMECENUBARLOW;
  }
  if (splMECENUBARHIGH) {
    delete splMECENUBARHIGH;
  }
  if (splEBC) {
    delete splEBC;
  }
  if (splEBO) {
    delete splEBO;
  }
  if (splSCCV) {
    delete splSCCV;
  }
  if (splSCCA) {
    delete splSCCA;
  }
  if (spl2P2HEDEP_LOWENU) {
    delete spl2P2HEDEP_LOWENU;
  }
  if (spl2P2HEDEP_HIGHENU) {
    delete spl2P2HEDEP_HIGHENU;
  }
  if (spl2P2HEDEP_LOWENUBAR) {
    delete spl2P2HEDEP_LOWENUBAR;
  }
  if (spl2P2HEDEP_HIGHENUBAR) {
    delete spl2P2HEDEP_HIGHENUBAR;
  }
  if (splISO_BKG_LOWPPI) {
    delete splISO_BKG_LOWPPI;
  }
  if (splDIS_BY) {
	delete splDIS_BY;
  }
  if (splMPI_BY) {
	delete splMPI_BY;
  }
  if (splMPI_AGKY_XSEC) {
	delete splMPI_AGKY_XSEC;
  }
}

// Constructor
xsecBase::xsecBase() {
  mode    = __BAD_INT__;
  species = __BAD_INT__;
  target  = __BAD_INT__;
  Q2      = __BAD_DOUBLE__;
  Enu     = __BAD_DOUBLE__;
  weight  = __BAD_DOUBLE__;
}

void xsecBase::Print() {
  std::cout << "*** Printing cross-section info:" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "species : "<< species << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "mode    : "<< mode << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "target  : "<< target << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Enu     : "<< Enu <<  " GeV" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Q2      : "<< Q2 <<  " GeV2" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "weight  : "<< weight << std::endl;
}


// ************************
// Print the cross-section information for one given event
void xsec2015::Print() {
// ************************
  xsecBase::Print();
  std::cout << "    " << std::left << std::setw(20) << "relRPA  : "<< relRPA << std::endl;

  std::cout << "    Splines: " << std::endl;
  if (splMAQE != NULL)    std::cout << "    " << std::left << std::setw(25) << "MAQE" << std::endl;
  if (splPDDC != NULL)    std::cout << "    " << std::left << std::setw(25) << "PDDC" << std::endl;
  if (splPDDO != NULL)    std::cout << "    " << std::left << std::setw(25) << "PDDO" << std::endl;
  if (splMECENULOW != NULL)    std::cout << "    " << std::left << std::setw(25) << "MECENULOW" << std::endl;
  if (splMECENUHIGH != NULL)    std::cout << "    " << std::left << std::setw(25) << "MECENUHIGH" << std::endl;
  if (splMECENUBARLOW != NULL)    std::cout << "    " << std::left << std::setw(25) << "MECENUBARLOW" << std::endl;
  if (splMECENUBARHIGH != NULL)    std::cout << "    " << std::left << std::setw(25) << "MECENUBARHIGH" << std::endl;
  if (splCA5 != NULL)     std::cout << "    " << std::left << std::setw(25) << "CA5" << std::endl;
  if (splMARES != NULL)   std::cout << "    " << std::left << std::setw(25) << "MARES" << std::endl;
  if (splBGSCLLOWPPI != NULL)   std::cout << "    " << std::left << std::setw(25) << "I12 LOW PPI" << std::endl;
  if (splBGSCL != NULL)   std::cout << "    " << std::left << std::setw(25) << "I12" << std::endl;
  if (splBYDIS != NULL)  std::cout << "    " << std::left << std::setw(25) << "BY DIS" << std::endl;
  if (splBYMPI != NULL)  std::cout << "    " << std::left << std::setw(25) << "BY MPI" << std::endl;
  if (splAGKYMULT != NULL)  std::cout << "    " << std::left << std::setw(25) << "AGKY MULT" << std::endl;
  if (splFEFABS != NULL)  std::cout << "    " << std::left << std::setw(25) << "FSI PI ABS" << std::endl;
  if (splFEFCX != NULL)   std::cout << "    " << std::left << std::setw(25) << "FSI CEX LO" << std::endl;
  if (splFEFQE != NULL)   std::cout << "    " << std::left << std::setw(25) << "FSI INEL LO" << std::endl;
  if (splFEFINEL != NULL) std::cout << "    " << std::left << std::setw(25) << "FSI PI PROD" << std::endl;
  if (splFEFCXH != NULL)  std::cout << "    " << std::left << std::setw(25) << "FSI CEX HI" << std::endl;
  if (splFEFQEH != NULL)  std::cout << "    " << std::left << std::setw(25) << "FSI INEL HI" << std::endl;
  if (spl2P2HEDEP_LOWENU != NULL)  std::cout << "    " << std::left << std::setw(25) << "2P2H EDEP LOWENU" << std::endl;
  if (spl2P2HEDEP_HIGHENU != NULL)  std::cout << "    " << std::left << std::setw(25) << "2P2H EDEP HIGHENU" << std::endl;
  if (spl2P2HEDEP_LOWENUBAR != NULL)  std::cout << "    " << std::left << std::setw(25) << "2P2H EDEP LOWENUBAR" << std::endl;
  if (spl2P2HEDEP_HIGHENUBAR != NULL)  std::cout << "    " << std::left << std::setw(25) << "2P2H EDEP HIGHENUBAR" << std::endl;
  if (splISO_BKG_LOWPPI != NULL)  std::cout << "    " << std::left << std::setw(25) << "ISO_BKG_LOWPPI" << std::endl;
  if (splDIS_BY != NULL)  std::cout << "    " << std::left << std::setw(25) << "DIS_BY" << std::endl;
  if (splMPI_BY != NULL)  std::cout << "    " << std::left << std::setw(25) << "MPI_BY" << std::endl;
  if (splMPI_AGKY_XSEC != NULL)  std::cout << "    " << std::left << std::setw(25) << "MPI_AGKY_XSEC" << std::endl;

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

  for(int i=1; i<poly->GetNumberOfBins()+1; i++)
    {
      integral += poly->GetBinContent(i);
    }

  return integral;

} // end function

// **************************************************
// Helper function for projecting TH2Poly onto the X axis
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins) {
  // **************************************************

  TH1D* hProjX = new TH1D((TempName+"_x").c_str(),(TempName+"_x").c_str(),xbins.size()-1,&xbins[0]);
  double xlow, xup, frac=0;

  //loop over bins in the poly
  for(int i=0; i<((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      xlow=bin->GetXMin();
      xup=bin->GetXMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dx=0; dx<int(xbins.size()); dx++)
        {
          if(xbins[dx+1]<=xlow || xbins[dx]>=xup)
            {
              frac=0;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]>=xup)
            {
              frac=1;
            }
          else if(xbins[dx]<=xlow && xbins[dx+1]<=xup)
            {
              frac=(xbins[dx+1]-xlow)/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]>=xup)
            {
              frac=(xup-xbins[dx])/(xup-xlow);
            }
          else if(xbins[dx]>=xlow && xbins[dx+1]<=xup)
            {
              frac=(xbins[dx+1]-xbins[dx])/(xup-xlow);
            }
          else
            {
              frac=0;
            }
          hProjX->SetBinContent(dx+1,hProjX->GetBinContent(dx+1)+frac*bin->GetContent());
        }
    }
  return hProjX;
} // end project poly X function

// **************************************************
// Helper function for projecting TH2Poly onto the Y axis
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins) {
  // **************************************************

  TH1D* hProjY = new TH1D((TempName+"_y").c_str(),(TempName+"_y").c_str(),ybins.size()-1,&ybins[0]);
  double ylow, yup, frac=0;

  //loop over bins in the poly
  for(int i=0; i<((TH2Poly*)poly)->GetNumberOfBins(); i++)
    {
      //get bin and its edges
      TH2PolyBin* bin = (TH2PolyBin*)((TH2Poly*)poly)->GetBins()->At(i)->Clone();
      ylow=bin->GetYMin();
      yup=bin->GetYMax();

      //Loop over projected bins, find fraction of poly bin in each
      for(int dy=0; dy<int(ybins.size()); dy++)
        {
          if(ybins[dy+1]<=ylow || ybins[dy]>=yup)
            {
              frac=0;
            }
          else if(ybins[dy]<=ylow && ybins[dy+1]>=yup)
            {
              frac=1;
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
        }
    }
  return hProjY;
} // end project poly Y function

//DB Get the Cernekov momentum threshold in MeV
double returnCherenkovThresholdMomentum(int PDG) {
  double refractiveIndex = 1.334; //DB From https://github.com/fiTQun/fiTQun/blob/646cf9c8ba3d4f7400bcbbde029d5ca15513a3bf/fiTQun_shared.cc#L757
  double mass =  MaCh3Utils::GetMassFromPDG(PDG)*1e3;
  double momentumThreshold = mass/sqrt(refractiveIndex*refractiveIndex-1.);
  return momentumThreshold;
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
