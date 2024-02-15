#include "covarianceOsc.h"

covarianceOsc::covarianceOsc(const char* name, const char *file, TH2D *hist_dcpth13NH, TH2D *hist_dcpth13IH, TH2D *hist_23)
: covarianceBase(name, file) {

  gRandom->SetSeed(0);

  if (hist_dcpth13NH) {
    h_dcpth13NH = hist_dcpth13NH;
    std::cout << "Using delta cp and th13 correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN THESE PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_dcpth13NH = NULL;
  }

  if (hist_dcpth13IH) {
    h_dcpth13IH = hist_dcpth13IH;
    std::cout << "Using delta cp and th13 IH correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN THESE PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_dcpth13IH = NULL;
  }

  if (hist_23) {
    h_23 = hist_23;
    std::cout << "Using th23 and dm23 correlated prior from histogram" << std::endl;
    std::cout << "NOTE: THIS WILL IGNORE ANY CORRELATIONS BETWEEN th_23 PARAMETERS AND OTHER PARAMETERS" << std::endl;
  } else {
    h_23 = NULL;
  }

  //Read in osc pars from xml file
  TFile *infile = new TFile(file, "READ");
  osc_prior = (TVectorD*)infile->Get("osc_nom");
  TVectorD* osc_stepscale = (TVectorD*)infile->Get("osc_stepscale");
  TVectorD* osc_sigma = (TVectorD*)infile->Get("osc_sigma");
  TVectorD* osc_flat_prior = (TVectorD*)infile->Get("osc_flat_prior");

  TObjArray* objarr_name = (TObjArray*)(infile->Get("osc_param_names"));

  TVectorD* osc_baseline = (TVectorD*)infile->Get("osc_baseline");
  TVectorD* osc_density = (TVectorD*)infile->Get("osc_density");
  double fScale = 1.0;

  //KS: Save all neccesary information from covariance
  for(int io = 0; io <size; io++)
  {
    char name[21];
    fParNames[io] = new Char_t[21];
      
    sprintf(name, std::string(((TObjString*)objarr_name->At(io))->GetString()).c_str());
    strcpy(fParNames[io], name);
    
    fParInit[io]  = (*osc_prior)(io);
    fParCurr[io] = fParProp[io] = fParInit[io];
    fParSigma[io] = (*osc_sigma)(io);
    fIndivStepScale[io] = fScale * (*osc_stepscale)(io);
    
    //KS: Set flat prior
    if( (bool)((*osc_flat_prior)(io)) ) setEvalLikelihood(io,false);
  }
    
  L = (*osc_baseline)(0);
  density = (*osc_density)(0);

  flipdelM=false;
  reactorPrior = false;
  flipBeta=false;

  fixdm23NH = -999;
  fixdm23IH = -999;
  fixth23NH = -999;
  fixth23IH = -999;

  genPropKernels();
  randomize();
  throwNominal();

  oscpars1 = new double[10];
  
  //KS:those are constant no need to overwrite them each time
  oscpars1[6] = 2;
  oscpars1[7] = L;
  oscpars1[8] = density;
  
  Print();
  CheckOrderOfParams();
    
  infile->Close();
  delete infile;
  
  std::cout << "created oscillation parameter handler" << std::endl;
}


covarianceOsc::~covarianceOsc()
{
  //  delete[] fPropKernel;
}

// void covarianceOsc::genPosteriorHists()
// {
//   Char_t name[256];
//   posterior[0] = new TH1D("post_osc_0", "sin2th_12",100,0.5, 1);
//   posterior[1] = new TH1D("post_osc_1", "sin2th_23",100, 0.7, 1);
//   posterior[2] = new TH1D("post_osc_2", "sin2th_13",100, 0, 0.3);
//   posterior[3] = new TH1D("post_osc_3", "delm_12",100, 0, 0.00001);
//   posterior[4] = new TH1D("post_osc_4", "delm_23",100, 0, 0.004);
//   posterior[5] = new TH1D("post_osc_5", "delta_cp",100, 0, 0.3);


//   /*
//   for (int i = 0; i < size; i++)
//     {
//       posterior[i] = NULL;
//       sprintf(name, "post_%s_%03i",getName(), i);
//       posterior[i] = new TH1D(name, name, 100, 0, 1);
//     }*/
// }

// CWRET changed this from double to bool!
bool covarianceOsc::checkBounds() {

  // wrap delta cp 
  // NB: previous to ~summer 2016 steps outside of this range were instead rejected
  if(fParProp[5] > TMath::Pi()) {
    fParProp[5] = (-2.*TMath::Pi() + fParProp[5]);

  } else if (fParProp[5] < -TMath::Pi()) {

    fParProp[5] = (2.*TMath::Pi() + fParProp[5]);
  }

  // ensure osc params dont go unphysical
  if (fParProp[0] > 1.0 || fParProp[0] < 0 ||
      fParProp[1] > 1.0 || fParProp[1] < 0 ||
      fParProp[2] > 1.0 || fParProp[2] < 0 
      //|| fParProp[4] < 0.0 || fParProp[4] > 20E-3 // dont let dm32 go to IH
      //|| fabs(fParProp[4]) > 0.004 || fabs(fParProp[4]) < 0.001 
      //|| fParProp[5] < -1*TMath::Pi() || fParProp[5] > TMath::Pi()
     )
  {
    return false;
  }

  if(size==7) {
    if(fParProp[6]<0) { // Don't let beta be less than 0 (no upper limit)
      return false;
    }
  }

  return true;
}

double covarianceOsc::GetLikelihood() {

  double logL = 0.0;
  if(size==6) {
    for(int i = 0; i < covMatrix->GetNrows()-1; ++i) //delta is last, don't loop over, as it gets flat prior
    {
      for(int j = 0; j < covMatrix->GetNcols()-1; ++j) //delta is last, don't loop over, as it gets flat prior
      {
        // If parameter 23 histogram exists, use that instead of matrix
        if (h_23 && (i==1 || i==4 || j==1 || j==4))
          continue;
        // If th13-dcp histogram exists, use that instead of matrix
        if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2)) // dcp not in loop
          continue;

        if(fParEvalLikelihood[i] && fParEvalLikelihood[j])
        {
          /*std::cout << i << " " << j << " nominal i " << nominal[i] << " nominal j " << nominal[j] << std::endl;
            std::cout << "parcurr " << fParProp[i] << " invmatrix " << InvertCovMatrix[i][j]; << std::endl;
            std::cout << "diff " << fParProp[i] - nominal[i] << std::endl;
            std::cout << "likelihood " << 0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*InvertCovMatrix[i][j]; << std::endl;*/

          //check
          if (h_23 && (i==1 || i==4 || j==1 || j==4))
            std::cout << "Error: using matrix when parameter-23 histogram exists" << std::endl;
          if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2))
            std::cout << "Error: using matrix when theta13-dcp histogram exists" << std::endl;

          logL+=0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*InvertCovMatrix[i][j];;
        }
      }
    }

    // If parameter 23 histogram exists, use that instead of the matrix
    if (h_23)
    {
      if (fParEvalLikelihood[1] && fParEvalLikelihood[4])
        logL+=TMath::Log(h_23->Interpolate(fParProp[1],fParProp[4]))*-1.0;
      else if (fParEvalLikelihood[1] || fParEvalLikelihood[4])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta23 and dm23 for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }

    // If 23-parameters are fixed but you're flipping MH, evaluate the MH prior
    if (fParSigma[1]<0.0 && fParSigma[4]<0.0 && flipdelM)
    {
      if (fParProp[4]==fixdm23NH)
        logL+=TMath::Log(0.684)*-1.0;
      else if (fParProp[4]==fixdm23IH)
        logL+=TMath::Log(0.316)*-1.0;
      else
        std::cerr << "ERROR: dm23 = " << fParProp[4] << ", fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << std::endl;
    }

    // If theta13-dcp histogram exists, use that instead of the matrix
    if (h_dcpth13IH && fParProp[4]<0) // if IH
    {
      if (fParEvalLikelihood[2] && fParEvalLikelihood[5])
        logL += TMath::Log(h_dcpth13IH->Interpolate(fParProp[2],fParProp[5]))*-1.0;
      else if (fParEvalLikelihood[2] || fParEvalLikelihood[5])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }      
    else if (h_dcpth13NH) // if NH or IH histogram doesn't exist (default to NH)
    {
      if (fParEvalLikelihood[2] && fParEvalLikelihood[5])
        logL += TMath::Log(h_dcpth13NH->Interpolate(fParProp[2],fParProp[5]))*-1.0;
      else if (fParEvalLikelihood[2] || fParEvalLikelihood[5])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }
    else if (fParEvalLikelihood[5]) // Evaluate likelihood for dcp (if histogram doesn't exist, th13 prior will already have been evaluated from matrix)
      logL+=1/(2.0*TMath::Pi()); // flat prior for delta
  }

  else if(size==7)
  {
    for(int i = 0; i < covMatrix->GetNrows()-2; i++) //delta is last, don't loop over, as it gets flat prior
    {
      for(int j = 0; j < covMatrix->GetNcols()-2; j++) //delta is last, don't loop over, as it gets flat prior
      {
        // If parameter 23 histogram exists, use that instead of matrix
        if (h_23 && (i==1 || i==4 || j==1 || j==4))
          continue;
        // If th13-dcp histogram exists, use that instead of matrix
        if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2)) // dcp not in loop
          continue;

        if(fParEvalLikelihood[i] && fParEvalLikelihood[j])
        {
          /*std::cout << i << " " << j << " nominal i " << nominal[i] << " nominal j " << nominal[j] << std::endl;
            std::cout << "parcurr " << fParProp[i] << " invmatrix " << InvertCovMatrix[i][j]; << std::endl;
            std::cout << "diff " << fParProp[i] - nominal[i] << std::endl;
            std::cout << "likelihood " << 0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*InvertCovMatrix[i][j]; << std::endl;*/

          //check
          if (h_23 && (i==1 || i==4 || j==1 || j==4))
            std::cout << "Error: using matrix when parameter-23 histogram exists" << std::endl;
          if ((h_dcpth13NH || h_dcpth13IH) && (i==2 || j==2))
            std::cout << "Error: using matrix when theta13-dcp histogram exists" << std::endl;

          logL+=0.5*(fParProp[i] - nominal[i])*(fParProp[j] - nominal[j])*InvertCovMatrix[i][j];;
        }
      }
    }

    // If parameter 23 histogram exists, use that instead of the matrix
    if (h_23)
    {
      if (fParEvalLikelihood[1] && fParEvalLikelihood[4])
        logL+=TMath::Log(h_23->Interpolate(fParProp[1],fParProp[4]))*-1.0;
      else if (fParEvalLikelihood[1] || fParEvalLikelihood[4])
        std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta23 and dm23 for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
    }

    // If 23-parameters are fixed but you're flipping MH, evaluate the MH prior
    if (fParSigma[1]<0.0 && fParSigma[4]<0.0 && flipdelM)
    {
      if (fParProp[4]==fixdm23NH)
        logL+=TMath::Log(0.684)*-1.0;
      else if (fParProp[4]==fixdm23IH)
        logL+=TMath::Log(0.316)*-1.0;
      else
        std::cerr << "ERROR: dm23 = " << fParProp[4] << ", fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << std::endl;
    }
    
    // If theta13-dcp histogram exists, use that instead of the matrix
    if (h_dcpth13IH && fParProp[4]<0) // IH
      {
	if (fParEvalLikelihood[2] && fParEvalLikelihood[5])
	  logL += TMath::Log(h_dcpth13IH->Interpolate(fParProp[2],fParProp[5]))*-1.0;
	else if (fParEvalLikelihood[2] || fParEvalLikelihood[5])
	  std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
      }      
    else if (h_dcpth13NH) // if NH or only one histogram given (default to NH)
      {
	if (fParEvalLikelihood[2] && fParEvalLikelihood[5])
	  logL += TMath::Log(h_dcpth13NH->Interpolate(fParProp[2],fParProp[5]))*-1.0;
	else if (fParEvalLikelihood[2] || fParEvalLikelihood[5])
	  std::cout << "ERROR: Must have setEvalLikelihood(true) for BOTH theta13 and dcp for histogram prior to work. Not using prior for these parameters - fix it!" << std::endl;
      }
    else if (fParEvalLikelihood[5]) // Evaluate likelihood for dcp (if histogram doesn't exist, th13 prior will already have been evaluated from matrix)
      logL+=1/(2.0*TMath::Pi()); // flat prior for delta
    
    // Evaluate likelihood for beta
    // Commented out = flat prior 
    //if (fParEvalLikelihood[6])
    // logL+=1/(3.0); 
  }
  
  // Check parameter boundaries:  if bad (checkBound == false), add large logL
  //                              if good (checkBound == true), add nothing
  if (!checkBounds()) {
    logL += __LARGE_LOGL__;
  }

  //std::cout << "oscpars: " << fParProp[0] << "  " << fParProp[1] << "  " << fParProp[2] << "  " << fParProp[3] << "  " << fParProp[4] << "  " << fParProp[5] <<std::endl;
  //std::cout << "oscllh pre-RC: " << logL << std::endl;

  // reactor prior
  if (reactorPrior)
  {
    // Reactor prior from 2013 PDG: sin^2(2theta13) = 0.095 +/- 0.01
    // Reactor prior from 2014 PDG: sin^2(2theta13) = 0.093 +/- 0.008
    // Reactor prior from 2015 PDG: sin^2(2theta13) = 0.085 +/- 0.005
    // Reactor prior from 2016ca PDG: sin^2(2theta13) = 0.0857 +/- 0.0046
    // Next convert between single angle (what we have) and double angle (what the PDG gives)
    // double dblang = 4*fParProp[2]*(1-fParProp[2]);
    // double tmp = (dblang-0.0857) / 0.0046;
    
    // Reactor prior from 2018 PDG: sin^2(theta13) = 0.0212 +/- 0.0008
    // Reactor prior from 2019 PDG: sin^2(theta13) = 0.0218 +/- 0.0007
    // Reactor prior from 2021 PDG: sin^2(theta13) = 0.0220 +/- 0.0007
    // This time we don't have to convert between single<->double angle, PDG gives RC in single angle.

    // Now calculate penalty           
    double tmp = (fParProp[2]-0.0220) / 0.0007;
    //double tmp = (dblang-0.095) / 0.01;

    // this line for T2K joint fit result, NOT REACTOR CONSTRAINT!
    //double tmp = (fParProp[2] - 0.04571) / 0.01125; 

    // Finally: add to likelihood
    logL += 0.5 * tmp * tmp;
    //std::cout << "oscllh post-RC: " << logL << std::endl;
  }

  //std::cout << "dm = " << fParProp[4] << ", logl osc = " << logL << std::endl;

  return logL;
}

void covarianceOsc::throwNominal(bool nomValues)
{
  //   cout << "entering throwNominal " << nomValues << endl;

  TDecompChol chdcmp(*covMatrix);

  if(!chdcmp.Decompose())
  {
    std::cerr << "Cholesky decomposition failed for " << matrixName << std::endl;
    exit(-1);
  }

  chel = new TMatrixD(chdcmp.GetU());
  CholeskyDecomp(size, *chel);

  ThrowParms* nom_throws = new ThrowParms(/**vec*/(*osc_prior), (*covMatrix));
  nom_throws->SetSeed(0);
  nominal.clear();
  nominal.resize(size);
  if(!nomValues)
  {
    nom_throws->ThrowSet(nominal);
    for (int i = 0; i < 3; i++)//only the angles
    {
      if (nominal[i] > 1.0) nominal[i] = 1.0;
      if (nominal[i] < 0.0) nominal[i] = 0.0;
    }
    if(nominal[5]>TMath::Pi())
      nominal[5]=TMath::Pi();
    if(nominal[5]<-1*TMath::Pi())
      nominal[5]=-1*TMath::Pi();
    if(size==7)
      if(nominal[6]<0)
        nominal[6]=0.0;
  }
  else
  {
    for (int i = 0; i < int(nominal.size()); i++)
      nominal[i]=(*osc_prior)(i);
  }
  delete nom_throws;

}

double *covarianceOsc::getPropPars()
{
  for(int i = 0; i < 6; i++)
    oscpars1[i] = fParProp[i];

  //Those are constant we initalised them already
  //oscpars1[6] = 2;
  //oscpars1[7] = L;
  //oscpars1[8] = density;
  if(size==7)
    oscpars1[9]=fParProp[6];
  else
    oscpars1[9]=1;

  return oscpars1;
}

void covarianceOsc::proposeStep() {
  for (int i = 0; i < size; i++)
  {
    if (fParSigma[i] > 0.0)
    {
      /*if (i == 4) // Don't do gaussian proposal - randomly sample dm23 from a uniform distribution
        {
        fParProp[i] = random_number[0]->Uniform() * 20E-3;
        }
        else */
      fParProp[i] = fParCurr[i] + fParSigma[i] * fPropKernel[i]->GetRandom()*fStepScale*fIndivStepScale[i]; //random_number[0]->Gaus(0, fParSigma[i]);  

      if (i==4 && random_number[0]->Uniform()<0.5 && flipdelM) // flip sign for delm2_23 (after taking step)
        fParProp[i]*=-1;

      if(i==6 && flipBeta)
      {
        if(random_number[0]->Uniform()<0.5)
        {
          if(fParCurr[i]==0.0)
            fParProp[i]=1.0;
          if(fParCurr[i]==1.0)
            fParProp[i]=0.0;
        }
        else 
          fParProp[i]=fParCurr[i];
      }

    }
  }

  if(random_number[0]->Uniform()<0.5){
    // flip octant around point of maximal disappearance (0.5112)
    // this ensures we move to a parameter value which has the same oscillation probability
    fParProp[1] = 0.5112 - (fParProp[1] - 0.5112);
  }

  // if flipdelM and parameter 4 is fixed, flip between fixed parameters for two hierarchies (Note: this will flip parameters 4 and 1 - dm23 *and* theta23).
  if (fParSigma[4] < 0.0 && flipdelM)
  {
    // First: check that the NH and IH fixed values were set
    if (fixdm23NH == -999 || fixdm23IH == -999 || fixth23NH == -999 || fixth23IH == -999)
      std::cerr << "ERROR: cannot flip between fixed values of normal and inverted heirarchy if values are not assigned. You have provided: fixdm23NH = " << fixdm23NH << ", fixdm23IH = " << fixdm23IH << ", fixth23NH = " << fixth23NH << ", fixth23IH = " << fixth23IH << std::endl;

    if (fParSigma[1] > 0.0) std::cout << "ERROR: cannot use this code to flip heirarchies unless BOTH dm23 and th23 are fixed" << std::endl;

    // Check fParCurr corresponds to one of the values given
    if (!(fParCurr[4] == fixdm23NH || fParCurr[4] == fixdm23IH))
      std::cout << "ERROR: fParCurr[4] = " << fParCurr[4] << " is not equal to fixdm23NH (" << fixdm23NH << ") or fixdm23IH (" << fixdm23IH << "). Not changing this parameter." << std::endl;

    if (!(fParCurr[1] == fixth23IH || fParCurr[1] == fixth23NH))
      std::cout << "ERROR: fParCurr[1] = " << fParCurr[1] << " is not equal to fixth23NH (" << fixth23NH << ") or fixth23IH (" << fixth23IH << "). Not changing this parameter." << std::endl;

    // Flip parameters together
    double a = random_number[0]->Uniform();
    if (a<0.5) 
    {
      if (fParCurr[4] == fixdm23NH)
      {
        fParProp[4] = fixdm23IH;
        fParProp[1] = fixth23IH;
      }
      else if (fParCurr[4] == fixdm23IH) 
      {
        fParProp[4] = fixdm23NH;
        fParProp[1] = fixth23NH;
      }
    }
    else 
    {
      fParProp[4] = fParCurr[4];
      fParProp[1] = fParCurr[1];
    }
    //std::cout << a << "\t" << fParProp[4] << "\t" << fParProp[1] << std::endl;
  }

}

std::vector<double> covarianceOsc::defaultPars(bool doubled)
{
  std::vector<double> oscpars;
  if(doubled)
  {
    oscpars.push_back(0.857);
    oscpars.push_back(1.0);
    oscpars.push_back(0.098);
    oscpars.push_back(7.5E-5);
    oscpars.push_back(2.5E-3);
    oscpars.push_back(0);
  }
  else
  {
    oscpars.push_back(0.311);
    oscpars.push_back(0.5);
    oscpars.push_back(0.0251);
    oscpars.push_back(7.5e-05);
    oscpars.push_back(0.0024);
    oscpars.push_back(0);
  }
  return oscpars;
}

void covarianceOsc::setExtraBranches(TTree &tree)
{
  // set branches to save current and proposed osc pars for dm_32 and th_23
  for (int i = 0; i<size; ++i)
  {
    if (!(i==1 || i==4)) continue;

    char bit_c[1024] = "c_";
    strcat(bit_c, fParNames[i]);
    strcat(bit_c, "/D");
    char name_c[1024] = "c_";
    strcat(name_c, fParNames[i]);
    tree.Branch(name_c, (double*)&fParCurr[i], bit_c);

    char bit_p[1024] = "p_";
    strcat(bit_p, fParNames[i]);
    strcat(bit_p, "/D");
    char name_p[1024] = "p_";
    strcat(name_p, fParNames[i]);
    tree.Branch(name_p, (double*)&fParProp[i], bit_p);
  }
}

// Overload setFlipDeltaM23 to provide values for fixed NH and IH dm23 and th23
void covarianceOsc::setFlipDeltaM23(double dm23NH, double dm23IH, double th23NH, double th23IH)
{
  fixdm23NH = dm23NH;
  fixdm23IH = dm23IH;
  fixth23NH = th23NH;
  fixth23IH = th23IH;

  flipdelM = true;
}

//KS: Print all usefull informations after initialization
void covarianceOsc::Print() {
  std::cout << "Number of pars: " << size << std::endl;
  std::cout << "current " << matrixName << " parameters:" << std::endl;
  std::cout << std::left << std::setw(5) << "#" << std::setw(2) << "|" << std::setw(25) << "Name" << std::setw(2) << "|" << std::setw(10) << "Nom." << std::setw(2) << "|" << std::setw(15) << "IndivStepScale" << std::setw(2) << "|" <<std::setw(15) << "fParSigma"  << std::endl;
  for(int i = 0; i < size; i++) {
    std::cout << std::fixed << std::setprecision(5) << std::left << std::setw(5) << i << std::setw(2) << "|" << std::setw(25) << fParNames[i] << std::setw(2) << "|" << std::setw(10) << fParInit[i]<< std::setw(2) << "|" << std::setw(15) << fIndivStepScale[i] << std::setw(2) << "|" << std::setw(15)<< fParSigma[i]<< std::endl;
  }
  
  std::cout<<"Baseline: "<<L<<std::endl;
  std::cout<<"Earth Density: "<<density<<std::endl;
}


//KS: Currently prob3++/probgp requiers particular order so we need to check this is the case
void covarianceOsc::CheckOrderOfParams() 
{
    std::vector<int> wrongParam;
    bool wrongMatrix = false;
    if(strcmp( fParNames[0], "sin2th_12") != 0 ){wrongParam.push_back(0); wrongMatrix = true;};
    if(strcmp( fParNames[1], "sin2th_23") != 0 ){wrongParam.push_back(1); wrongMatrix = true;};
    if(strcmp( fParNames[2], "sin2th_13") != 0 ){wrongParam.push_back(2); wrongMatrix = true;};
    if(strcmp( fParNames[3], "delm2_12")  != 0 ){wrongParam.push_back(3); wrongMatrix = true;};
    if(strcmp( fParNames[4], "delm2_23")  != 0 ){wrongParam.push_back(4); wrongMatrix = true;};
    if(strcmp( fParNames[5], "delta_cp")  != 0 ){wrongParam.push_back(5); wrongMatrix = true;};
    if(size == 7 && strcmp( fParNames[6], "beta")  != 0 ){wrongParam.push_back(6); wrongMatrix = true;};

        
    if(wrongMatrix)
    {
        for(unsigned int i =0; i < wrongParam.size(); i++ )
        {
            std::cerr << "Osc Patameter "<<fParNames[i]<<" isn't in good order"<<std::endl;  
        }
        std::cerr << "Currently prob3++/probgp requiers particular order"<< std::endl;
        std::cerr << "Please modify XML and make new matrix with good order"<< std::endl;
        std::cerr << "Find me here "<<__FILE__ << ":" << __LINE__ << std::endl;
        throw;  
    }

}
