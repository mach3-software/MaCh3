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
