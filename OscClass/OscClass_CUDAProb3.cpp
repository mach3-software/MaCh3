#include "OscClass_CUDAProb3.h"

Oscillator::Oscillator(std::string ConfigName) {

  //####################################################################################
  //Set definite values

  nOscpars = 6;

  nNeutrinoSigns = 2;
  nInitialFlavours = 2;
  nFinalFlavours = 3;  //DB =2 if excluding taus, =3 if including taus

  HistogramsInitialised = false;
  nPrimaryHists = nNeutrinoSigns*nInitialFlavours*nFinalFlavours;

  //####################################################################################
  //Get FitMan parameters

  manager* osc_manager = new manager(ConfigName);

  EarthDensityFile = osc_manager->raw()["OscillatorObj"]["OscillatorEarthDensityFile"].as<std::string>();
  ArrayConfig = osc_manager->raw()["OscillatorObj"]["OscillatorArrayConfig"].as<int>();
  fFillHistograms = osc_manager->raw()["OscillatorObj"]["OscillatorFillHistograms"].as<bool>();

  IsLinear = osc_manager->raw()["OscillatorObj"]["OscillatorIsLinear"].as<bool>();
  useFineBinsPerBin = osc_manager->raw()["OscillatorObj"]["OscillatorUseFineBinsPerCoarseBin"].as<bool>();

  UseBinningTemplates = osc_manager->raw()["OscillatorObj"]["OscillatorUseBinningTemplates"].as<bool>();
  PrimaryBinningTemplateName = osc_manager->raw()["OscillatorObj"]["OscillatorPrimaryBinningTemplateName"].as<std::string>();
  SecondaryBinningTemplateName = osc_manager->raw()["OscillatorObj"]["OscillatorSecondaryBinningTemplateName"].as<std::string>();
  InputFileName = osc_manager->raw()["OscillatorObj"]["OscillatorInputFileName"].as<std::string>();

  nCoarseCosz = osc_manager->raw()["OscillatorObj"]["OscillatornCoarseCosz"].as<int>();
  lCoarseCosz = osc_manager->raw()["OscillatorObj"]["OscillatorlCoarseCosz"].as<double>();
  hCoarseCosz = osc_manager->raw()["OscillatorObj"]["OscillatorhCoarseCosz"].as<double>();

  nCoarseEnergy = osc_manager->raw()["OscillatorObj"]["OscillatornCoarseEnergy"].as<int>();
  lCoarseEnergy = osc_manager->raw()["OscillatorObj"]["OscillatorlCoarseEnergy"].as<double>();
  hCoarseEnergy = osc_manager->raw()["OscillatorObj"]["OscillatorhCoarseEnergy"].as<double>();

  nFineCosz = osc_manager->raw()["OscillatorObj"]["OscillatornFineCosz"].as<int>();
  nFineEnergy = osc_manager->raw()["OscillatorObj"]["OscillatornFineEnergy"].as<int>();
  fineCoarseRatioCosz = osc_manager->raw()["OscillatorObj"]["OscillatorfineCoarseRatioCosz"].as<double>();
  fineCoarseRatioEnergy = osc_manager->raw()["OscillatorObj"]["OscillatorfineCoarseRatioEnergy"].as<double>();

  UseProductionHeightAveraging = osc_manager->raw()["OscillatorObj"]["OscillatorUseProductionHeightAveraging"].as<bool>();
  ProductionHeightFileName = osc_manager->raw()["OscillatorObj"]["OscillatorProductionHeightFileName"].as<std::string>();
  nProductionHeightAveragingBins = osc_manager->raw()["OscillatorObj"]["OscillatornProductionHeightAveragingBins"].as<int>();
  lProductionHeightRange = osc_manager->raw()["OscillatorObj"]["OscillatorlProductionHeightRange"].as<double>();
  hProductionHeightRange = osc_manager->raw()["OscillatorObj"]["OscillatorhProductionHeightRange"].as<double>();

  UseChemicalComposition = osc_manager->raw()["OscillatorObj"]["OscillatorUseChemicalComposition"].as<bool>();

  RebinMode = osc_manager->raw()["OscillatorObj"]["OscillatorRebinMode"].as<bool>(); //I do wonder whether I need this as making this work with production height averaging is going to be a pain
  std::cout << std::endl;

  /*
  //####################################################################################

  TH2D* PrimaryHistTemplate;
  TH2D* SecondaryHistTemplate;

  if (UseBinningTemplates) {

    TFile* InFile = new TFile(InputFileName.c_str(),"READ");
    if (InFile->IsZombie()) {
      std::cout << "Oscillator Template file is not found. Given:" << InputFileName << std::endl;
      throw;
    }

    PrimaryHistTemplate = (TH2D*)InFile->Get(PrimaryBinningTemplateName.c_str());
    if (!PrimaryHistTemplate) {
      std::cout << "Primary Histogram template not found. Given:" << PrimaryBinningTemplateName << std::endl;
      throw;
    }
    PrimaryHistTemplate->SetDirectory(0);

    SecondaryHistTemplate = (TH2D*)InFile->Get(SecondaryBinningTemplateName.c_str());
    if (!SecondaryHistTemplate) {
      std::cout << "Secondary Histogram template not found. Given:" << SecondaryBinningTemplateName << std::endl;
      throw;
    }
    SecondaryHistTemplate->SetDirectory(0);

    InFile->Close();

    nCoarseCosz = PrimaryHistTemplate->GetNbinsY();
    lCoarseCosz = PrimaryHistTemplate->GetYaxis()->GetBinLowEdge(1);
    hCoarseCosz = PrimaryHistTemplate->GetYaxis()->GetBinLowEdge(nCoarseCosz+1);

    nCoarseEnergy = PrimaryHistTemplate->GetNbinsX();
    lCoarseEnergy = PrimaryHistTemplate->GetXaxis()->GetBinLowEdge(1);
    hCoarseEnergy = PrimaryHistTemplate->GetXaxis()->GetBinLowEdge(nCoarseEnergy+1);

    nFineCosz = SecondaryHistTemplate->GetNbinsY();
    nFineEnergy = SecondaryHistTemplate->GetNbinsX();

    lFineCosz = lCoarseCosz;
    hFineCosz = hCoarseCosz;
    lFineEnergy = lCoarseEnergy;
    hFineEnergy = hCoarseEnergy;

    std::cout << "Oscillator set to use normal mode" << std::endl;
    std::cout << "Initialising oscillograms in TH2D constuctor" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Using Templates -" << std::endl;
    std::cout << "Primary:" << std::setw(20) << PrimaryBinningTemplateName << std::endl;
    std::cout << "Secondary:" << std::setw(20) << SecondaryBinningTemplateName << std::endl;
    std::cout << "From file:" << std::setw(20) << InputFileName << std::endl;

  } else {
    if (!IsLinear) {
      lCoarseEnergy = log10(lCoarseEnergy);
      hCoarseEnergy = log10(hCoarseEnergy);
    }

    lFineCosz = lCoarseCosz;
    hFineCosz = hCoarseCosz;
    lFineEnergy = lCoarseEnergy;
    hFineEnergy = hCoarseEnergy;

    std::vector<double> Energy_Coarse;
    std::vector<double> Cosz_Coarse;
    std::vector<double> Energy_Fine;
    std::vector<double> Cosz_Fine;

    if (IsLinear) {
      Energy_Coarse = logspace(lCoarseEnergy,hCoarseEnergy,nCoarseEnergy);
    } else {
      Energy_Coarse = linspace(lCoarseEnergy,hCoarseEnergy,nCoarseEnergy);
    }
    Cosz_Coarse = linspace(lCoarseCosz,hCoarseCosz,nCoarseCosz);

    if (!useFineBinsPerBin) {

      if (IsLinear) {
	Energy_Fine = logspace(lFineEnergy,hFineEnergy,nFineEnergy);
      } else {
	Energy_Fine = linspace(lFineEnergy,hFineEnergy,nFineEnergy);
      }
      Cosz_Fine = linspace(lFineCosz,hFineCosz,nFineCosz);

    } else {
      if (nFineEnergy != floor(nFineEnergy)) {
	std::cerr << "nFineEnergy:" << nFineEnergy << " is not an integer" << std::endl;
        throw;
      }

      if (nFineCosz != floor(nFineCosz)) {
	std::cerr << "nFineCosz:" << nFineCosz << " is not an integer" << std::endl;
        throw;
      }

      if (nFineEnergy < 1) {
	std::cerr << "Invalid binning specified. nFineEnergy should be integer >= 1. Given:" << nFineEnergy << std::endl;
        throw;
      }

      if (nFineCosz < 1) {
	std::cerr << "Invalid binning specified. nFineCosz should be integer >= 1. Given:" << nFineCosz << std::endl;
        throw;
      }

      nFineCosz = nFineCosz*nCoarseCosz;
      nFineEnergy = nFineEnergy*nCoarseEnergy;

      Cosz_Fine = ReturnFineBinningFromCoarseBinnnig(nFineCosz,Cosz_Coarse);
      Energy_Fine = ReturnFineBinningFromCoarseBinnnig(nFineEnergy,Energy_Coarse);
    }

    if (fineCoarseRatioCosz == -1 || fineCoarseRatioEnergy == -1) {

      PrimaryHistTemplate = new TH2D("hPrimaryBinning","PrimaryBinning",Energy_Coarse.size()-1,Energy_Coarse.data(),Cosz_Coarse.size()-1,Cosz_Coarse.data());
      SecondaryHistTemplate = new TH2D("hSecondaryBinning","SecondaryBinning",Energy_Fine.size()-1,Energy_Fine.data(),Cosz_Fine.size()-1,Cosz_Fine.data());

    } else if ((fineCoarseRatioCosz > 0 && fineCoarseRatioCosz <= 1.0) && (fineCoarseRatioEnergy > 0 && fineCoarseRatioEnergy <= 1.0)) {
      nFineEnergy *= fineCoarseRatioEnergy;
      nFineCosz *= fineCoarseRatioCosz;

      PrimaryHistTemplate = new TH2D("hPrimaryBinning","PrimaryBinning",nCoarseEnergy,lCoarseEnergy,hCoarseEnergy,nCoarseCosz,lCoarseCosz,hCoarseCosz);
      SecondaryHistTemplate = new TH2D("hSecondaryBinning","SecondaryBinning",nFineEnergy,lFineEnergy,hFineEnergy,nFineCosz,lFineCosz,hFineCosz);
    } else {
      std::cerr << "fineCoarseRatioCosz is an invalid number:" << fineCoarseRatioCosz << std::endl;
      std::cerr << "fineCoarseRatioEnergy is an invalid number:" << fineCoarseRatioEnergy << std::endl;
      std::cerr << "Expected both to be either -1 or 0 > fineCoarseRatio >= 1.0" << std::endl;
      throw;
    }
  }

  if (RebinMode) {
    nMaxBin = 5000*5000;
    std::cout << "Oscillator set to use RebinMode. Max Bins:" << nMaxBin << std::endl;
    std::cout << "If for whatever reason you need to increase this, go to Oscillator::Oscillator(bool RebinMode_, TString EarthDensityFile_)" << std::endl;

    std::cout << "Production height averaging does not work for a rebinnable oscillator" << std::endl;
    UseProductionHeightAveraging = false;
    nProductionHeightAveragingBins = 0;
  } else {
    std::cout << "Oscillator set to use normal mode" << std::endl;
  }

  InitOscillogram(PrimaryHistTemplate,SecondaryHistTemplate);
  if (fFillHistograms) {
    InitialiseHistograms();
  }

  isUsingGPU();
  std::cout << std::endl;
  PrintOscillatorConfig();
  std::cout << std::endl;
  ResetSavedParams();

  DefineMiscValues();
  DefinePropagatorEnums();
  InitPropagator();
  */
}

/*
void Oscillator::RebinOscillogram(int Switcher, std::vector<double> NewBinning) {
  if (!RebinMode) {
    std::cout << "Rebin called when Oscillator not in rebin mode. Skipping.." << std::endl;
    throw;
  }

  bool CoarseChanged = false;
  bool FineChanged = false;

  //DB This would only be needed in oscillogram binning studies and in normal running this would not be needed. Hence console output to remind user to call needed function
  std::cout << "Due to the newly implemented pointer method of returning oscillogram weights, please ensure that you call samplePDFSKBase::FindEventOscBin() after Oscillogram::RebinOscillogram() to update the osc_w_pointers" << std::endl;

  if (Switcher==0) {
    if (NewBinning.size()==2) {
      if ((NewBinning[0]*NewBinning[1]) > nMaxBin) {
	std::cout << "Rebinned Binning set too large. Go to Oscillator::Oscillator(bool RebinMode_, TString EarthDensityFile_) to increase nMaxBins. Skipping.." << std::endl;
	throw;
      }

      if (NewBinning[0] > nFineCosz) {
	std::cout << "Requested coarse cosZ binning is larger than the fine cosZ binning" << std::endl;
	throw;
      }

      if (NewBinning[1] > nFineEnergy) {
	std::cout << "Requested coarse energy binning is larger than the fine energy binning" << std::endl;
        throw;
      }

      nCoarseCosz = (int)NewBinning[0];
      nCoarseEnergy = (int)NewBinning[1];
    } else {
      std::cout << "Rebin with current params not allowed. Switcher should be 0 or 3 to rebin Primary, 1 to rebin Secondary or 2 to rebin both. Switcher given:" << Switcher << std::endl;
      std::cout << "Size of binning array should be 2 for Switcher == 0, 1, 3, or 4 for Switcher == 2. Size given:" << sizeof(NewBinning)/sizeof(NewBinning[0]) << std::endl;
      throw;
    }

    CoarseChanged = true;
  }

  else if (Switcher==1) {
    if (NewBinning.size()==2) {
      if ((NewBinning[0]*NewBinning[1]) > nMaxBin) {
	std::cout << "Rebinned Binning set too large. Go to Oscillator::Oscillator(bool RebinMode_, TString EarthDensityFile_) to increase nMaxBins. Skipping.." << std::endl;
        throw;
      }

      if (NewBinning[0] < nCoarseCosz) {
	std::cout << "Requested fine cosZ binning is smaller than the coarse cosZ binning" << std::endl;
        throw;
      }

      if (NewBinning[1] < nCoarseEnergy) {
	std::cout << "Requested fine energy binning is smaller than the coarse energy binning" << std::endl;
        throw;
      }

      nFineCosz = (int)NewBinning[0];
      nFineEnergy = (int)NewBinning[1];
    } else {
      std::cout << "Rebin with current params not allowed. Switcher should be 0 or 3 to rebin Primary, 1 to rebin Secondary or 2 to rebin both. Switcher given:" << Switcher << std::endl;
      std::cout << "Size of binning array should be 2 for Switcher == 0, 1, 3, or 4 for Switcher == 2. Size given:" << sizeof(NewBinning)/sizeof(NewBinning[0]) << std::endl;
      throw;
    }

    FineChanged = true;
  }

  else if (Switcher==2) {
    if (NewBinning.size()==4) {
      if (((NewBinning[0]*NewBinning[1]) > nMaxBin)||((NewBinning[2]*NewBinning[3]) > nMaxBin)) {
	std::cout << "Rebinned Binning set too large. Go to Oscillator::Oscillator(bool RebinMode_, TString EarthDensityFile_) to increase nMaxBins. Skipping.." << std::endl;
        throw;
      }
      nCoarseCosz = (int)NewBinning[0];
      nCoarseEnergy = (int)NewBinning[1];
      nFineCosz = (int)NewBinning[2];
      nFineEnergy = (int)NewBinning[3];
    } else {
      std::cout << "Rebin with current params not allowed. Switcher should be 0 or 3 to rebin Primary, 1 to rebin Secondary or 2 to rebin both. Switcher given:" << Switcher << std::endl;
      std::cout << "Size of binning array should be 2 for Switcher == 0, 1, 3, or 4 for Switcher == 2. Size given:" << sizeof(NewBinning)/sizeof(NewBinning[0]) << std::endl;
      throw;
    }

    CoarseChanged = true;
    FineChanged = true;
  }

  else if (Switcher==3) {
    if (NewBinning.size()==2) {
      if ((NewBinning[0]<0)||(NewBinning[0]>1)) {
	std::cout << "Rebinned factorised value for CosZ outside of 0 to 1 range:" << NewBinning[0] << std::endl;
	std::cout << "Skipping.." << std::endl;
      }
      if ((NewBinning[1]<0)||(NewBinning[1]>1)) {
	std::cout << "Rebinned factorised value for Energy outside of 0 to 1 range:" << NewBinning[1] << std::endl;
	std::cout << "Skipping.." << std::endl;
      }
      int tmp_nCoarseCosz = NewBinning[0]*nFineCosz;
      int tmp_nCoarseEnergy = NewBinning[1]*nFineEnergy;
      if ((tmp_nCoarseCosz*tmp_nCoarseEnergy) > nMaxBin) {
	std::cout << "Rebinned Binning set too large. Go to Oscillator::Oscillator(bool RebinMode_, TString EarthDensityFile_) to increase nMaxBins. Skipping.." << std::endl;
        throw;
      }
      nCoarseCosz = tmp_nCoarseCosz;
      nCoarseEnergy = tmp_nCoarseEnergy;
    } else {
      std::cout << "Rebin with current params not allowed. Switcher should be 0 or 3 to rebin Primary, 1 to rebin Secondary or 2 to rebin both. Switcher given:" << Switcher << std::endl;
      std::cout << "Size of binning array should be 2 for Switcher == 0, 1, 3, or 4 for Switcher == 2. Size given:" << sizeof(NewBinning)/sizeof(NewBinning[0]) << std::endl;
      throw;
    }

    CoarseChanged = true;
  }

  else if (Switcher == 4) {
    if (NewBinning.size()==2) {
      if (NewBinning[0] != floor(NewBinning[0])) {
	std::cerr << "NewBinning[0]:" << NewBinning[0] << " is not an integer" << std::endl;
	throw;
      }

      if (NewBinning[1] != floor(NewBinning[1])) {
	std::cerr << "NewBinning[1]:" << NewBinning[1] << " is not an integer" << std::endl;
	throw;
      }

      if (NewBinning[0] < 1) {
	std::cerr << "Invalid binning specified. NewBinning[0] should be integer >= 1. Given:" << NewBinning[0] << std::endl;
	throw;
      }

      if (NewBinning[1] < 1) {
	std::cerr << "Invalid binning specified. NewBinning[0] should be integer >= 1. Given:" << NewBinning[1] << std::endl;
	throw;
      }

      nFineCosz = (int)NewBinning[0]*nCoarseCosz;
      nFineEnergy = (int)NewBinning[1]*nCoarseEnergy;

    } else {
      std::cout << "Rebin with current params not allowed. Switcher should be 0 or 3 to rebin Primary, 1 or 4 to rebin Secondary or 2 to rebin both. Switcher given:" << Switcher << std::endl;
      std::cout << "Size of binning array should be 2 for Switcher == 0, 1, 3, 4 or 4 for Switcher == 2. Size given:" << sizeof(NewBinning)/sizeof(NewBinning[0]) << std::endl;
      throw;
    }

    FineChanged = true;
  }

  else {
    std::cout << "Invalid switcher option given. Switcher given:" << Switcher << std::endl;
    throw;
  }

  ResizeOscillogram(CoarseChanged,FineChanged);

  if (fFillHistograms) {
    DeleteHistograms();
    InitialiseHistograms();
  }

  ResetSavedParams();
}

void Oscillator::ResizeOscillogram(bool CoarseChanged, bool FineChanged) {
  DeleteArrays();

  std::vector<double> Energy_Coarse;
  std::vector<double> Cosz_Coarse;
  std::vector<double> Energy_Fine;
  std::vector<double> Cosz_Fine;

  if (CoarseChanged) {
    if (IsLinear) {
      Energy_Coarse = logspace(lCoarseEnergy,hCoarseEnergy,nCoarseEnergy);
    } else {
      Energy_Coarse = linspace(lCoarseEnergy,hCoarseEnergy,nCoarseEnergy);
    }

    Cosz_Coarse = linspace(lCoarseCosz,hCoarseCosz,nCoarseCosz);
  } else {
    Energy_Coarse.resize(hPrimaryBinning->GetNbinsX()+1);
    for (int i=0;i<=hPrimaryBinning->GetNbinsX();i++) {
      Energy_Coarse[i] = hPrimaryBinning->GetXaxis()->GetBinLowEdge(i+1);
    }

    Cosz_Coarse.resize(hPrimaryBinning->GetNbinsY()+1);
    for (int i=0;i<=hPrimaryBinning->GetNbinsY();i++) {
      Cosz_Coarse[i] = hPrimaryBinning->GetYaxis()->GetBinLowEdge(i+1);
    }
  }

  if (FineChanged) {

    if (!useFineBinsPerBin) {
      if (IsLinear) {
	Energy_Fine = logspace(lFineEnergy,hFineEnergy,nFineEnergy);
      } else {
	Energy_Fine = linspace(lFineEnergy,hFineEnergy,nFineEnergy);
      }
      Cosz_Fine = linspace(lFineCosz,hFineCosz,nFineCosz);
    } else {
      Cosz_Fine = ReturnFineBinningFromCoarseBinnnig(nFineCosz,Cosz_Coarse);
      Energy_Fine = ReturnFineBinningFromCoarseBinnnig(nFineEnergy,Energy_Coarse);
    }

  } else {
    Energy_Fine.resize(hSecondaryBinning->GetNbinsX()+1);
    for (int i=0;i<=hSecondaryBinning->GetNbinsX();i++) {
      Energy_Fine[i] = hSecondaryBinning->GetXaxis()->GetBinLowEdge(i+1);
    }

    Cosz_Fine.resize(hSecondaryBinning->GetNbinsY()+1);
    for (int i=0;i<=hSecondaryBinning->GetNbinsY();i++) {
      Cosz_Fine[i] = hSecondaryBinning->GetYaxis()->GetBinLowEdge(i+1);
    }
  }

  hPrimaryBinning->SetBins(Energy_Coarse.size()-1,Energy_Coarse.data(),Cosz_Coarse.size()-1,Cosz_Coarse.data());
  hSecondaryBinning->SetBins(Energy_Fine.size()-1,Energy_Fine.data(),Cosz_Fine.size()-1,Cosz_Fine.data());

  nPrimaryBinsX = hPrimaryBinning->GetNbinsX();
  nPrimaryBinsY = hPrimaryBinning->GetNbinsY();
  nPrimaryBins = nPrimaryBinsX*nPrimaryBinsY;

  nSecondaryBinsX = hSecondaryBinning->GetNbinsX();
  nSecondaryBinsY = hSecondaryBinning->GetNbinsY();
  nSecondaryBins = nSecondaryBinsX*nSecondaryBinsY;

  std::cout << "---------------------------------" << std::endl;
  std::cout << "         " << std::setw(7) << "nBins" << " " << std::setw(10) << "LowerEdge" << " " << std::setw(10) << "HigherEdge" << std::endl;
  std::cout << "Coarse Oscillogram Binning: " << std::endl;
  if (IsLinear) {
    std::cout << "Energy = " << std::setw(7) << nCoarseEnergy << " " << std::setw(10) << lCoarseEnergy  << " " << std::setw(10) << hCoarseEnergy << std::endl;
  } else {
    std::cout << "Energy = " << std::setw(7) << nCoarseEnergy << " " << std::setw(10) << pow(10,lCoarseEnergy)  << " " << std::setw(10) << pow(10,hCoarseEnergy) << std::endl;
  }
  std::cout << "Cosz   = " << std::setw(7) << nCoarseCosz << " " << std::setw(10) << lCoarseCosz  << " " << std::setw(10) << hCoarseCosz << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "Fine Oscillogram Binning: " << std::endl;
  if (IsLinear) {
    std::cout << "Energy = " << std::setw(7) << nFineEnergy << " " << std::setw(10) << lFineEnergy  << " " << std::setw(10) << hFineEnergy << std::endl;
  } else {
    std::cout << "Energy = " << std::setw(7) << nFineEnergy << " " << std::setw(10) << pow(10,lFineEnergy)  << " " << std::setw(10) << pow(10,hFineEnergy) << std::endl;
  }
  std::cout << "Cosz   = " << std::setw(7) << nFineCosz << " " << std::setw(10) << lFineCosz  << " " << std::setw(10) << hFineCosz << std::endl;
  std::cout << "---------------------------------" << std::endl;

  if (hPrimaryOscillogram.size()!=0) {
    for (int i=0;i<nNeutrinoSigns;i++) {
      for (int j=0;j<nInitialFlavours;j++) {
	for (int k=0;k<nFinalFlavours;k++) {

	  hPrimaryOscillogram[i][j][k]->SetBins(Energy_Coarse.size()-1,Energy_Coarse.data(),Cosz_Coarse.size()-1,Cosz_Coarse.data());
	  hSecondaryOscillogram[i][j][k]->SetBins(Energy_Fine.size()-1,Energy_Fine.data(),Cosz_Fine.size()-1,Cosz_Fine.data());

	}
      }
    }
  }

  ResizeArrays();

  //DB According to LP, unique pointers should auto manage their own memory usage and once a new object is help in a unique pointer, it will 'delete' the old object
  DeletePropagator();
  InitPropagator();
  FillArrays();
}

void Oscillator::DeleteArrays() {
  for (unsigned int iPrimaryHist=0;iPrimaryHist<nPrimaryHists;iPrimaryHist++) {
    delete[] hPrimaryOscillogram_Arr[iPrimaryHist];
  }
  delete[] hPrimaryOscillogram_Arr;

  delete[] hPrimaryCounter_Arr;

  for (int iBin=0;iBin<nPrimaryBins;iBin++) {
    PrimaryBinContrib_Bin[iBin].clear();
    PrimaryBinContrib_Weight[iBin].clear();
  }

  PrimaryXBinEdges = std::vector<double>();
  PrimaryYBinEdges = std::vector<double>();

  SecondaryXBinEdges = std::vector<double>();
  SecondaryYBinEdges = std::vector<double>();

}

void Oscillator::ResizeArrays() {

  hPrimaryOscillogram_Arr = new double*[nPrimaryHists];
  for (unsigned int iPrimaryHist=0;iPrimaryHist<nPrimaryHists;iPrimaryHist++) {
    hPrimaryOscillogram_Arr[iPrimaryHist] = new double[nPrimaryBins];
  }
  hPrimaryCounter_Arr = new double[nPrimaryBins];

  PrimaryBinContrib_Bin.resize(nPrimaryBins);
  PrimaryBinContrib_Weight.resize(nPrimaryBins);

  PrimaryXBinEdges.resize(hPrimaryBinning->GetNbinsX()+1);
  PrimaryYBinEdges.resize(hPrimaryBinning->GetNbinsY()+1);
  SecondaryXBinEdges.resize(hSecondaryBinning->GetNbinsX()+1);
  SecondaryYBinEdges.resize(hSecondaryBinning->GetNbinsY()+1);

  PrimaryXBinEvalPoints.resize(hPrimaryBinning->GetNbinsX());
  PrimaryYBinEvalPoints.resize(hPrimaryBinning->GetNbinsY());
  SecondaryXBinEvalPoints.resize(hSecondaryBinning->GetNbinsX());
  SecondaryYBinEvalPoints.resize(hSecondaryBinning->GetNbinsY());

  for (int yBin=1;yBin<=hPrimaryBinning->GetNbinsY()+1;yBin++) {
    PrimaryYBinEdges[yBin-1] = hPrimaryBinning->GetYaxis()->GetBinLowEdge(yBin);
  }

  for (int yBin=1;yBin<=hSecondaryBinning->GetNbinsY()+1;yBin++) {
    SecondaryYBinEdges[yBin-1] = hSecondaryBinning->GetYaxis()->GetBinLowEdge(yBin);
  }

  for (int yBin=1;yBin<=hPrimaryBinning->GetNbinsY();yBin++) {
    PrimaryYBinEvalPoints[yBin-1] = hPrimaryBinning->GetYaxis()->GetBinCenter(yBin);
  }

  for (int yBin=1;yBin<=hSecondaryBinning->GetNbinsY();yBin++) {
    SecondaryYBinEvalPoints[yBin-1] = hSecondaryBinning->GetYaxis()->GetBinCenter(yBin);
    if (fabs(SecondaryYBinEvalPoints[yBin-1])<1e-9) SecondaryYBinEvalPoints[yBin-1] = 0.;
  }

  if (IsLinear) {
    for (int xBin=1;xBin<=hPrimaryBinning->GetNbinsX()+1;xBin++) {
      PrimaryXBinEdges[xBin-1] = hPrimaryBinning->GetXaxis()->GetBinLowEdge(xBin);
    }

    for (int xBin=1;xBin<=hSecondaryBinning->GetNbinsX()+1;xBin++) {
      SecondaryXBinEdges[xBin-1] = hSecondaryBinning->GetXaxis()->GetBinLowEdge(xBin);
    }

    for (int xBin=1;xBin<=hPrimaryBinning->GetNbinsX();xBin++) {
      PrimaryXBinEvalPoints[xBin-1] = hPrimaryBinning->GetXaxis()->GetBinCenter(xBin);
    }

    for (int xBin=1;xBin<=hSecondaryBinning->GetNbinsX();xBin++) {
      SecondaryXBinEvalPoints[xBin-1] = hSecondaryBinning->GetXaxis()->GetBinCenter(xBin);
    }

  } else {
    for (int xBin=1;xBin<=hPrimaryBinning->GetNbinsX()+1;xBin++) {
      PrimaryXBinEdges[xBin-1] = pow(10.0,hPrimaryBinning->GetXaxis()->GetBinLowEdge(xBin));
    }

    for (int xBin=1;xBin<=hSecondaryBinning->GetNbinsX()+1;xBin++) {
      SecondaryXBinEdges[xBin-1] = pow(10.0,hSecondaryBinning->GetXaxis()->GetBinLowEdge(xBin));
    }

    for (int xBin=1;xBin<=hPrimaryBinning->GetNbinsX();xBin++) {
      PrimaryXBinEvalPoints[xBin-1] = pow(10.0,hPrimaryBinning->GetXaxis()->GetBinCenter(xBin));
    }

    for (int xBin=1;xBin<=hSecondaryBinning->GetNbinsX();xBin++) {
      SecondaryXBinEvalPoints[xBin-1] = pow(10.0,hSecondaryBinning->GetXaxis()->GetBinCenter(xBin));
    }
  }

}

void Oscillator::FillArrays() {
  if (ArrayConfig == 1) {
    FillArrays_ManyContrib_Area();
  } else {
    FillArrays_Standard();
  }

  for (int iPrimaryBin=0;iPrimaryBin<nPrimaryBins;iPrimaryBin++) {
    unsigned int nContribs = PrimaryBinContrib_Bin[iPrimaryBin].size();

    double Val = 0.;
    for (unsigned int iContrib=0;iContrib<nContribs;iContrib++) {
      Val += PrimaryBinContrib_Weight[iPrimaryBin][iContrib];
    }
    if (Val == 0.) {
      std::cerr << "Coarse Bin found which has no contributions" << std::endl;
      throw;
    }

    hPrimaryCounter_Arr[iPrimaryBin] = Val;
  }
}

void Oscillator::FillArrays_ManyContrib_Area() {

  for (int x=0;x<nSecondaryBinsX;x++) {
    for (int y=0;y<nSecondaryBinsY;y++) {
      int iter = x*nSecondaryBinsY+y;

      Corner Fine_BL{SecondaryXBinEdges[x],SecondaryYBinEdges[y]};
      Corner Fine_UR{SecondaryXBinEdges[x+1],SecondaryYBinEdges[y+1]};
      Box Box_Fine{Fine_BL,Fine_UR};

      if (!IsValidBox(Box_Fine)) {
	std::cout << "Fine Binning Box is not valid!" << std::endl;
	PrintBox(Box_Fine);
	throw;
      }

      for (int yBin=0;yBin<nPrimaryBinsY;yBin++) {
	//If bottom edge of Coarse bin is above than the top edge of the Fine bin
	if (PrimaryYBinEdges[yBin] > SecondaryYBinEdges[y+1]) {
	  continue;
	}

	//If top edge of Coarse bin is lower than the bottom edge of the Fine bin
	if (PrimaryYBinEdges[yBin+1] < SecondaryYBinEdges[y]) {
	  continue;
	}

	for (int xBin=0;xBin<nPrimaryBinsX;xBin++) {

	  //If left edge of Coarse bin is above than the right edge of the Fine bin
	  if (PrimaryXBinEdges[xBin] > SecondaryXBinEdges[x+1]) {
	    continue;
	  }

	  //If right edge of Coarse bin is lower than the left edge of the Fine bin
	  if (PrimaryXBinEdges[xBin+1] < SecondaryXBinEdges[x]) {
	    continue;
	  }

	  Corner Coarse_BL{PrimaryXBinEdges[xBin],PrimaryYBinEdges[yBin]};
	  Corner Coarse_UR{PrimaryXBinEdges[xBin+1],PrimaryYBinEdges[yBin+1]};
	  Box Box_Coarse{Coarse_BL,Coarse_UR};

	  if (!IsValidBox(Box_Coarse)) {
	    std::cout << "Coarse Binning Box is not valid!" << std::endl;
	    PrintBox(Box_Coarse);
	    throw;
	  }

	  double FracOverlapped = FractionOverlapped(Box_Coarse,Box_Fine);
	  if (FracOverlapped<1e-6) FracOverlapped = 0.;

	  if (FracOverlapped>1.) {
	    PrintBox(Box_Coarse);
	    PrintBox(Box_Fine);
	    throw;
	  }

	  if (FracOverlapped != 0.) {
	    int PrimaryBin_ToAdd_X = xBin;
	    int PrimaryBin_ToAdd_Y = yBin;
	    int PrimaryBin_ToAdd = PrimaryBin_ToAdd_Y*nPrimaryBinsX+PrimaryBin_ToAdd_X;

	    PrimaryBinContrib_Bin[PrimaryBin_ToAdd].push_back(iter);
	    PrimaryBinContrib_Weight[PrimaryBin_ToAdd].push_back(FracOverlapped);
	  }

	}
      }
    }
  }

}

void Oscillator::FillArrays_Standard() {

  for (int x=0;x<nSecondaryBinsX;x++) {
    for (int y=0;y<nSecondaryBinsY;y++) {
      int iter = x*nSecondaryBinsY+y;

      double EnergyEvalPoint = SecondaryXBinEvalPoints[x];
      double CoszEvalPoint = SecondaryYBinEvalPoints[y];

      int PrimaryBin_ToAdd_X;

      if (IsLinear) {
	PrimaryBin_ToAdd_X = hPrimaryBinning->GetXaxis()->FindBin(EnergyEvalPoint)-1;
      } else {
	PrimaryBin_ToAdd_X = hPrimaryBinning->GetXaxis()->FindBin(log10(EnergyEvalPoint))-1;
      }

      int PrimaryBin_ToAdd_Y = hPrimaryBinning->GetYaxis()->FindBin(CoszEvalPoint)-1;
      int PrimaryBin_ToAdd = PrimaryBin_ToAdd_Y*nPrimaryBinsX+PrimaryBin_ToAdd_X;

      PrimaryBinContrib_Bin[PrimaryBin_ToAdd].push_back(iter);
      PrimaryBinContrib_Weight[PrimaryBin_ToAdd].push_back(1.);
    }
  }

}

void Oscillator::InitOscillogram(TH2D* PrimaryHistTemplate,TH2D* SecondaryHistTemplate) {
  std::cout << "---------------------------------" << std::endl;
  std::cout << "         " << std::setw(7) << "nBins" << " " << std::setw(10) << "LowerEdge" << " " << std::setw(10) << "HigherEdge" << std::endl;
  std::cout << "Coarse Binning: " << std::endl;
  if (IsLinear) {
    std::cout << "Energy = " << std::setw(7) << nCoarseEnergy << " " << std::setw(10) << lCoarseEnergy  << " " << std::setw(10) << hCoarseEnergy << std::endl;
  } else {
    std::cout << "Energy = " << std::setw(7) << nCoarseEnergy << " " << std::setw(10) << pow(10,lCoarseEnergy)  << " " << std::setw(10) << pow(10,hCoarseEnergy) << std::endl;
  }
  std::cout << "Cosz   = " << std::setw(7) << nCoarseCosz << " " << std::setw(10) << lCoarseCosz  << " " << std::setw(10) << hCoarseCosz << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "Fine Binning: " << std::endl;
  if (IsLinear) {
    std::cout << "Energy = " << std::setw(7) << nFineEnergy << " " << std::setw(10) << lFineEnergy  << " " << std::setw(10) << hFineEnergy << std::endl;
  } else {
    std::cout << "Energy = " << std::setw(7) << nFineEnergy << " " << std::setw(10) << pow(10,lFineEnergy)  << " " << std::setw(10) << pow(10,hFineEnergy) << std::endl;
  }
  std::cout << "Cosz   = " << std::setw(7) << nFineCosz << " " << std::setw(10) << lFineCosz  << " " << std::setw(10) << hFineCosz << std::endl;
  std::cout << "---------------------------------" << std::endl;

  hPrimaryBinning = (TH2D*)PrimaryHistTemplate->Clone("hPrimaryBinning");
  hPrimaryBinning->SetTitle("hPrimaryBinning");
  hSecondaryBinning = (TH2D*)SecondaryHistTemplate->Clone("hSecondaryBinning");
  hSecondaryBinning->SetTitle("SecondaryBinning");

  nPrimaryBinsX = hPrimaryBinning->GetNbinsX();
  nPrimaryBinsY = hPrimaryBinning->GetNbinsY();
  nPrimaryBins = nPrimaryBinsX*nPrimaryBinsY;

  nSecondaryBinsX = hSecondaryBinning->GetNbinsX();
  nSecondaryBinsY = hSecondaryBinning->GetNbinsY();
  nSecondaryBins = nSecondaryBinsX*nSecondaryBinsY;

  if (!RebinMode) {
    ProbList = new FLOAT_T[nSecondaryBins];
  }
  else {
    ProbList = new FLOAT_T[nMaxBin];
  }

  ResizeArrays();
  FillArrays();
}

void Oscillator::InitialiseHistograms() {
  hPrimaryOscillogram = std::vector< std::vector< std::vector< TH2D* > > >(nNeutrinoSigns);
  hSecondaryOscillogram = std::vector< std::vector< std::vector< TH2D* > > >(nNeutrinoSigns);

  TString NeutrinoSign;
  TString HistTitle;

  for (int i=0;i<nNeutrinoSigns;i++) {
    if (i==0) {NeutrinoSign = "#bar{#nu}";}
    else {NeutrinoSign = "#nu";}

    hPrimaryOscillogram[i].resize(nInitialFlavours);
    hSecondaryOscillogram[i].resize(nInitialFlavours);

    for (int j=0;j<nInitialFlavours;j++) {
      hPrimaryOscillogram[i][j].resize(nFinalFlavours);
      hSecondaryOscillogram[i][j].resize(nFinalFlavours);

      for (int k=0;k<nFinalFlavours;k++) {
        if      (j==0) {HistTitle = NeutrinoSign+"_{e} #rightarrow ";}
        else if (j==1) {HistTitle = NeutrinoSign+"_{#mu} #rightarrow ";}
        else           {HistTitle = NeutrinoSign+"_{#tau} #rightarrow ";}

        if      (k==0) {HistTitle += NeutrinoSign+"_{e}";}
        else if (k==1) {HistTitle += NeutrinoSign+"_{#mu}";}
        else           {HistTitle += NeutrinoSign+"_{#tau}";}

        hPrimaryOscillogram[i][j][k] = (TH2D*)hPrimaryBinning->Clone(Form("hPrimaryOscillogram_%i_%i_%i",i,j,k));
        hPrimaryOscillogram[i][j][k]->SetTitle(HistTitle+";log_{10}(E_{#nu}) (GeV);Cosine Zenith");
        hSecondaryOscillogram[i][j][k] = (TH2D*)hSecondaryBinning->Clone(Form("hSecondaryOscillogram_%i_%i_%i",i,j,k));
        hSecondaryOscillogram[i][j][k]->SetTitle(HistTitle+";log_{10}(E_{#nu}) (GeV);Cosine Zenith");

      }
    }
  }

  HistogramsInitialised = true;
}

void Oscillator::DeleteHistograms() {
  for (int i=0;i<nNeutrinoSigns;i++) {
    for (int j=0;j<nInitialFlavours;j++) {
      for (int k=0;k<nFinalFlavours;k++) {
	delete hPrimaryOscillogram[i][j][k];
      }
    }
  }
}

void Oscillator::DeleteOscillogram() {
  if (hPrimaryOscillogram.size()==nPrimaryHists) {
    for (int i=0;i<nNeutrinoSigns;i++) {
      for (int j=0;j<nInitialFlavours;j++) {
	for (int k=0;k<nFinalFlavours;k++) {
	  delete hPrimaryOscillogram[i][j][k];
	  delete hSecondaryOscillogram[i][j][k];
	}
	hPrimaryOscillogram[i][j].clear();
	hSecondaryOscillogram[i][j].clear();
      }
      hPrimaryOscillogram[i].clear();
      hSecondaryOscillogram[i].clear();
    }
    hPrimaryOscillogram.clear();
    hSecondaryOscillogram.clear();
  }

  hPrimaryBinning->Reset();
  hSecondaryBinning->Reset();

  PrimaryXBinEdges.clear();
  PrimaryYBinEdges.clear();

  SecondaryXBinEdges.clear();
  SecondaryYBinEdges.clear();

  for (int i=0;i<nPrimaryBins;i++) {
    PrimaryBinContrib_Bin[i].clear();
    PrimaryBinContrib_Weight[i].clear();
  }
  PrimaryBinContrib_Bin.clear();
  PrimaryBinContrib_Weight.clear();

  nPrimaryBinsX = -1;
  nPrimaryBinsY = -1;
  nPrimaryBins = -1;

  nSecondaryBinsX = -1;
  nSecondaryBinsY = -1;
  nSecondaryBins = -1;

  delete ProbList;

  ResetSavedParams();
}

void Oscillator::ResetSavedParams() {
  for (int i=0;i<6;i++) {
    foscpar[i] = -999;
  }
  fprodH = -999;
  fYp_Val = -1.;
  for (int i=0;i<4;i++) {
    fBinning[i] = -1;
  }
}

void Oscillator::Reset(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex) {
  int iPrimaryHist = NeutrinoSignIndex*(nInitialFlavours*nFinalFlavours)+InitialNeutrinoIndex*nFinalFlavours+FinalNeutrinoIndex;
  for (int iBin=0;iBin<nPrimaryBins;iBin++) {
    hPrimaryOscillogram_Arr[iPrimaryHist][iBin] = 0.;
  }

}

bool Oscillator::isAlreadyCalculated(double* oscpar, double prodH, double Yp_Val) {
  bool fAlreadyCalculated = true;
  for (int i=0;i<nOscpars;i++) {
    if (oscpar[i]!=foscpar[i]) {
      fAlreadyCalculated = false;
    }
  }
  if (prodH!=fprodH) {
    fAlreadyCalculated = false;
  }
  if (Yp_Val!=fYp_Val) {
    fAlreadyCalculated = false;
  }
  if (fBinning[0] != nCoarseCosz) {
    fAlreadyCalculated = false;
  }
  if (fBinning[1] != nCoarseEnergy) {
    fAlreadyCalculated = false;
  }
  if (fBinning[2] != nFineCosz) {
    fAlreadyCalculated = false;
  }
  if (fBinning[3] != nFineEnergy) {
    fAlreadyCalculated = false;
  }
  return fAlreadyCalculated;
}

void Oscillator::SaveParams(double* oscpar, double prodH, double Yp_Val) {
  for (int i=0;i<nOscpars;i++) {
    foscpar[i] = oscpar[i];
  }
  fprodH = prodH;
  fYp_Val = Yp_Val;
  fBinning[0] = nCoarseCosz;
  fBinning[1] = nCoarseEnergy;
  fBinning[2] = nFineCosz;
  fBinning[3] = nFineEnergy;
}

void  Oscillator::FillOscillogram(double* oscpar, double prodH, double Yp_Val) {
  if (isAlreadyCalculated(oscpar,prodH,Yp_Val)) {
    return;
  }

  SaveParams(oscpar,prodH,Yp_Val);

  //DB oscpar, as given from MaCh3, expresses the mixing angles in sin^2(theta). The propagator expects them in theta
  for (int iOscPar=0;iOscPar<3;iOscPar++) {
    if (oscpar[iOscPar] < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << oscpar[iOscPar] << std::endl;
      throw;
    }
  }

  FLOAT_T theta12 = asin(sqrt(oscpar[0]));
  FLOAT_T theta23 = asin(sqrt(oscpar[1]));
  FLOAT_T theta13 = asin(sqrt(oscpar[2]));
  FLOAT_T dm12sq  = oscpar[3];
  FLOAT_T dm23sq  = oscpar[4];
  FLOAT_T dcp     = oscpar[5];

  // Get how many layers are in the propagator's earth file
  nLayers = propagator->getNlayerBoundaries();
  std::vector<FLOAT_T> chemicalComposition_Dial(nLayers);
  // The electron density in the first layer is pretty well known, so don't change
  chemicalComposition_Dial[0] = 0.497;
  chemicalComposition_Dial[1] = 0.497;
  // Starting layer for different density
  int start = 2;
  if (nLayers == 6) {
    start = 3;
    chemicalComposition_Dial[2] = 0.497;
  }
  for (int i = start; i < nLayers; ++i) {
    chemicalComposition_Dial[i] = Yp_Val;
  }

  if (UseChemicalComposition) {
    propagator->setChemicalComposition(chemicalComposition_Dial);
  } else {
    propagator->setChemicalComposition(chemicalComposition_Nom);
  }
  propagator->setNeutrinoMasses(dm12sq, dm23sq);
  propagator->setProductionHeight(prodH);

  for (int i=0;i<nNeutrinoSigns;i++) {

    if (NeutrinoTypes[i]==Antineutrino) {
      // DB, Haven't really thought about it, but prob3++ sets dcp->-dcp here: https://github.com/rogerwendell/Prob3plusplus/blob/fd189e232e96e2c5ebb2f7bd3a5406b288228e41/BargerPropagator.cc#L235
      // Copying that behaviour gives same behaviour as prob3++/probGPU
      propagator->setMNSMatrix(theta12, theta13, theta23, -dcp);
    } else {
      propagator->setMNSMatrix(theta12, theta13, theta23, dcp);
    }

    propagator->calculateProbabilities(NeutrinoTypes[i]);

    for (int j=0;j<nInitialFlavours;j++) {
      for (int k=0;k<nFinalFlavours;k++) {

#ifdef DEBUG
	//DB This part could be replaced with pointers
	for (int iter=0;iter<nSecondaryBins;iter++) {
	  int yBin = int(iter%nSecondaryBinsY);
	  int xBin = int(iter/nSecondaryBinsY);
	  FLOAT_T prob = propagator->getProbability(yBin,xBin,OscChannels[j][k]);
	  if (std::isnan(prob)) {
	    std::cerr << "Prob is NAN!" << std::endl;
	    std::cerr << "In Flav:" << j << std::endl;
	    std::cerr << "Out Flav:" << k << std::endl;
	    std::cerr << "NeutrinoType:" << NeutrinoTypes[i] << std::endl;

	    std::cerr << "theta12:" << theta12 << std::endl;
	    std::cerr << "theta13:" << theta13 << std::endl;
	    std::cerr << "theta23:" << theta23 << std::endl;
	    std::cerr << "dcp:" << dcp << std::endl;

	    std::cerr << "prodH:" << prodH << std::endl;
	    std::cerr << "\n" << std::endl;

	    for (int iOscPar=0;iOscPar<3;iOscPar++) {
	      std::cerr << "iOscPar:" << iOscPar << std::endl;
	      std::cerr << "oscpar[iOscPar]:" << oscpar[iOscPar] << std::endl;
	      std::cerr << "sqrt(oscpar[iOscPar]):" << sqrt(oscpar[iOscPar]) << std::endl;
	      std::cerr << "asin(sqrt(oscpar[iOscPar])):" << asin(sqrt(oscpar[iOscPar])) << std::endl;
	      std::cerr << "\n" << std::endl;
	    }

	    throw;
	  }

	  ProbList[iter] = propagator->getProbability(yBin,xBin,OscChannels[j][k]);
	}

#else
	propagator->getProbabilityArr(ProbList, OscChannels[j][k]);

	
	//for (int iter=0;iter<nSecondaryBins;iter++) {
	  //int yBin = int(iter%nSecondaryBinsY);
	  //int xBin = int(iter/nSecondaryBinsY);
	  //ProbList[iter] = propagator->getProbability(yBin,xBin,OscChannels[j][k]);
	  //std::cout << "xBin:" << xBin << " | yBin:" << yBin << " | iter:" << iter << " | ProbList[iter]:" << ProbList[iter] << std::endl;
	  //}
#endif
        // Sometimes CUDAProb3 can return *slightly* unphysical oscillation probabilities
	for (int iter=0;iter<nSecondaryBins;iter++) {
	  ProbList[iter] = ProbList[iter] > 0 ? ProbList[iter] : 0;
          ProbList[iter] = (fabs(ProbList[iter]-1.0) < 1.0e-4) ? 1.0 : ProbList[iter];
        }

	FillPrimaryOscillogram(i,j,k);
	if (fFillHistograms) {
	  FillSecondaryHistograms(i,j,k);
	  FillPrimaryHistograms(i,j,k);
	}

      }
    }
  }

}

void Oscillator::FillPrimaryOscillogram(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex) {
  int iPrimaryHist = NeutrinoSignIndex*(nInitialFlavours*nFinalFlavours)+InitialNeutrinoIndex*nFinalFlavours+FinalNeutrinoIndex;

#ifdef MULTITHREAD
#pragma omp for
#endif
  for (int iPrimaryBin=0;iPrimaryBin<nPrimaryBins;iPrimaryBin++) {
    unsigned int nContribs = PrimaryBinContrib_Bin[iPrimaryBin].size();

    double Val = 0.;
    for (unsigned int iContrib=0;iContrib<nContribs;iContrib++) {
      Val += (PrimaryBinContrib_Weight[iPrimaryBin][iContrib]*ProbList[PrimaryBinContrib_Bin[iPrimaryBin][iContrib]])/hPrimaryCounter_Arr[iPrimaryBin];
    }
    hPrimaryOscillogram_Arr[iPrimaryHist][iPrimaryBin] = Val;
  }

}

void Oscillator::FillPrimaryHistograms(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex) {
  if (hPrimaryOscillogram.size()==0) {
    std::cout << "Primary oscillograms are not initialised so I'm not going to try and fill them" << std::endl;
    throw;
  }

  hPrimaryOscillogram[NeutrinoSignIndex][InitialNeutrinoIndex][FinalNeutrinoIndex]->Reset();
  int iPrimaryHist = NeutrinoSignIndex*(nInitialFlavours*nFinalFlavours)+InitialNeutrinoIndex*nFinalFlavours+FinalNeutrinoIndex;

  for (int iBin=0;iBin<nPrimaryBins;iBin++) {
    int xBin = int(iBin%nPrimaryBinsX)+1;
    int yBin = int(iBin/nPrimaryBinsX)+1;

    hPrimaryOscillogram[NeutrinoSignIndex][InitialNeutrinoIndex][FinalNeutrinoIndex]->SetBinContent(xBin,yBin,hPrimaryOscillogram_Arr[iPrimaryHist][iBin]);
  }
}

void Oscillator::FillSecondaryHistograms(int NeutrinoSignIndex, int InitialNeutrinoIndex, int FinalNeutrinoIndex) {
  if (hSecondaryOscillogram.size()==0) {
    std::cout << "Secondary oscillograms are not initialised so I'm not going to try and fill them" << std::endl;
    throw;
  }

  hSecondaryOscillogram[NeutrinoSignIndex][InitialNeutrinoIndex][FinalNeutrinoIndex]->Reset();

  for (int iter=0;iter<nSecondaryBins;iter++) {
    int yBin = int(iter%nSecondaryBinsY)+1;
    int xBin = int(iter/nSecondaryBinsY)+1;

    hSecondaryOscillogram[NeutrinoSignIndex][InitialNeutrinoIndex][FinalNeutrinoIndex]->SetBinContent(xBin, yBin, ProbList[iter]);
  }
}

std::vector<TH2D*> Oscillator::ReturnOscillogramArray(int fPrimary) {
  if (hPrimaryOscillogram.size()==0) {
    std::cout << "Oscillograms are not initialised so I'm not going to try and return them" << std::endl;
    throw;
  }


  std::vector<TH2D*> Array;
  TH2D* Hist;

  TString NeutrinoSign;

  for (int i=0;i<nNeutrinoSigns;i++) {
    for (int j=0;j<nInitialFlavours;j++) {
      for (int k=0;k<nFinalFlavours;k++) {

	if (fPrimary) {
	  Hist = (TH2D*)hPrimaryOscillogram[i][j][k]->Clone(Form("hPrimaryArray_%i_%i_%i",i,j,k));
	}
	else {
	  Hist = (TH2D*)hSecondaryOscillogram[i][j][k]->Clone(Form("hSecondaryArray_%i_%i_%i",i,j,k));
	}

	Array.push_back(Hist);
      }
    }
  }

  return Array;
}

double Oscillator::ReturnProb(double NeutrinoEnergy, double Cosz, int InitialFlavour, int FinalFlavour) {
  int NeutrinoSignIndex = NeutrinoSignToIndex(InitialFlavour);
  int InitialNeutrinoIndex = NeutrinoFlavourToIndex(InitialFlavour);
  int FinalNeutrinoIndex = NeutrinoFlavourToIndex(FinalFlavour);

  int iPrimaryHist = NeutrinoSignIndex*(nInitialFlavours*nFinalFlavours)+InitialNeutrinoIndex*nFinalFlavours+FinalNeutrinoIndex;

  if ((InitialNeutrinoIndex>nInitialFlavours)||(FinalNeutrinoIndex>nFinalFlavours)) {
    std::cout << "Something is wrong! Oscillation probability not calculated for neutrino flavour requested" << std::endl;
    std::cout << "NInitialFlavours:" << nInitialFlavours << " | Requested Flavour:" << InitialNeutrinoIndex << std::endl;
    std::cout << "NFinalFlavours:" << nFinalFlavours << " | Requested Flavour:" << FinalNeutrinoIndex << std::endl;
    std::exit(-1);
  }

  int xBin = -999;
  if (IsLinear) {
    xBin = hPrimaryBinning->GetXaxis()->FindBin(NeutrinoEnergy);
  } else {
    xBin = hPrimaryBinning->GetXaxis()->FindBin(log10(NeutrinoEnergy));
  }
  int yBin = hPrimaryBinning->GetYaxis()->FindBin(Cosz);

  if (xBin==-999) {
    std::cout << "xBin not set in Oscillator::ReturnProb. Quitting.." << std::endl;
    throw;
  }

  return hPrimaryOscillogram_Arr[iPrimaryHist][yBin*nPrimaryBinsX+xBin];
}

void Oscillator::isUsingGPU() {
#ifdef USE_GPU
    std::cout << "-------------------------" << std::endl;
    std::cout << "Set Oscillator to use GPU" << std::endl;
    std::cout << "-------------------------" << std::endl;
#else
    std::cout << "-------------------------" << std::endl;
    std::cout << "Set Oscillator to use CPU" << std::endl;
    std::cout << "-------------------------" << std::endl;
#endif
}

void Oscillator::SaveOscillogramsToFile(TString FileName) {
  if (hPrimaryOscillogram.size()==0) {
    std::cout << "Oscillograms are not initialised so I'm not going to try and save them" << std::endl;
    throw;
  }

  std::cout << "Saving Oscillograms to:" << FileName << std::endl;
  TFile* File = new TFile(FileName,"RECREATE");

  File->mkdir("FineOsc");
  File->cd("FineOsc");
  for (int i=0;i<nNeutrinoSigns;i++) {
    for (int j=0;j<nInitialFlavours;j++) {
      for (int k=0;k<nFinalFlavours;k++) {
	hPrimaryOscillogram[i][j][k]->Write();
      }
    }
  }

  File->mkdir("CoarseOsc");
  File->cd("CoarseOsc");
  for (int i=0;i<nNeutrinoSigns;i++) {
    for (int j=0;j<nInitialFlavours;j++) {
      for (int k=0;k<nFinalFlavours;k++) {
	hSecondaryOscillogram[i][j][k]->Write();
      }
    }
  }

  File->Close();
}

const double* Oscillator::retPointer(int GenNeutrinoFlavour, int DetNeutrinoFlavour, double NeutrinoEnergy, double TrueCZ) {
  if ((GenNeutrinoFlavour*DetNeutrinoFlavour)<0) {
    std::cout << "Generated and Detected neutrino flavour are different signs - Quitting" << std::endl;
    std::cout << "GenNeutrinoFlavour:" << GenNeutrinoFlavour << std::endl;
    std::cout << "DetNeutrinoFlavour:" << DetNeutrinoFlavour << std::endl;
    throw;
  }

  int NeutrinoSignIndex = NeutrinoSignToIndex(GenNeutrinoFlavour);
  int InitialNeutrinoIndex = NeutrinoFlavourToIndex(GenNeutrinoFlavour);
  int FinalNeutrinoIndex = NeutrinoFlavourToIndex(DetNeutrinoFlavour);

  unsigned int iPrimaryHist = NeutrinoSignIndex*(nInitialFlavours*nFinalFlavours)+InitialNeutrinoIndex*nFinalFlavours+FinalNeutrinoIndex;

  int xBin = -999;
  if (IsLinear) {
    xBin = hPrimaryBinning->GetXaxis()->FindBin(NeutrinoEnergy)-1;
  } else {
    xBin = hPrimaryBinning->GetXaxis()->FindBin(log10(NeutrinoEnergy))-1;
  }
  int yBin = hPrimaryBinning->GetYaxis()->FindBin(TrueCZ)-1;

  if (xBin==-999) {
    std::cout << "xBin not set in Oscillator::ReturnProb. Quitting.." << std::endl;
    throw;
  }

  int Bin = yBin*nPrimaryBinsX+xBin;

  if (iPrimaryHist>=nPrimaryHists) {
    std::cout << "iPrimaryHist out of range - Given:" << iPrimaryHist << std::endl;
    throw;
  }

  if ((Bin<0)||(Bin>nPrimaryBins)) {
    std::cout << "Bin out of range - Given:" << Bin << std::endl;
    throw;
  }

  return &hPrimaryOscillogram_Arr[iPrimaryHist][Bin];
}

void Oscillator::PrintBox(Box Box1) {
  std::cout << Box1.BL.x << "," << Box1.BL.y << " | " << Box1.UR.x << "," << Box1.UR.y << std::endl;
}

bool Oscillator::IsValidBox(Box Box1) {
  if (Box1.UR.y <= Box1.BL.y) return false;
  if (Box1.UR.x <= Box1.BL.x) return false;
  return true;
}

//Find Fraction of Box2 overlapping Box1
double Oscillator::FractionOverlapped(Box Box1, Box Box2) {
  double Box2_xLength = Box2.UR.x-Box2.BL.x;
  double Box2_yLength = Box2.UR.y-Box2.BL.y;
  double Box2_Area = Box2_xLength*Box2_yLength;

  double Area_Overlapped = 0;
  double xLength_Overlapped = std::min(Box1.UR.x,Box2.UR.x) - std::max(Box1.BL.x,Box2.BL.x);
  double yLength_Overlapped = std::min(Box1.UR.y,Box2.UR.y) - std::max(Box1.BL.y,Box2.BL.y);

  if (xLength_Overlapped > 0 && yLength_Overlapped > 0) {
    Area_Overlapped = xLength_Overlapped*yLength_Overlapped;
  }

  double Fraction_Area_Overlapped = Area_Overlapped/Box2_Area;

  return Fraction_Area_Overlapped;
}

void Oscillator::InitPropagator() {
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Initialising Propagator:" << std::endl;
  std::cout << " - nEnergyBins : " << nSecondaryBinsX << std::endl;
  std::cout << " - nCoszBins   : " << nSecondaryBinsY << std::endl;
  std::cout << std::endl;

  int nEnergy = nSecondaryBinsX;
  int nCosine = nSecondaryBinsY;

#ifdef USE_GPU
  propagator = std::unique_ptr<Propagator<FLOAT_T>> ( new CudaPropagatorSingle<FLOAT_T>(0,nCosine, nEnergy)); // Single-GPU propagator
#else
  int nThreads = 1;
  #pragma omp parallel
  {
#pragma omp single
    nThreads = omp_get_num_threads();
  }

  propagator = std::unique_ptr<Propagator<FLOAT_T>> ( new CpuPropagator<FLOAT_T>(nCosine, nEnergy, nThreads)); // MultiThread CPU propagator
#endif

  propagator->setEnergyList(SecondaryXBinEvalPoints);
  propagator->setCosineList(SecondaryYBinEvalPoints);

  CheckEarthDensityFile();
  std::cout << "Loading Earth density from: " << EarthDensityFile << std::endl;
  propagator->setDensityFromFile(EarthDensityFile.Data());

  if (UseProductionHeightAveraging) {
    SetProductionHeightArray();
  } else {
    std::cout << "Not using production height averaging" << std::endl;
  }


  std::cout << "--------------------------------------------" << std::endl;
}

void Oscillator::DeletePropagator() {
  delete propagator.release();
}

void Oscillator::DefineMiscValues() {
  SetProductionHeightBinEdges();

  if (!UseChemicalComposition) {
    std::cout << "Not using PREM model Chemical composition dial values" << std::endl;

    // Get how many layers are in the propagator's earth file
    nLayers = propagator->getNlayerBoundaries();
    chemicalComposition_Nom = std::vector<FLOAT_T>(nLayers);
    for (int i = 0; i < nLayers; ++i) {
      chemicalComposition_Nom[i] = 0.5;
    }
  }
}

void Oscillator::DefinePropagatorEnums() {

  NeutrinoTypes.resize(nNeutrinoSigns);
  NeutrinoTypes[0] = Antineutrino;
  NeutrinoTypes[1] = Neutrino;

  NeutrinoTypes_Names.resize(nNeutrinoSigns);
  NeutrinoTypes_Names[0] = "Antineutrino";
  NeutrinoTypes_Names[1] = "Neutrino";

  OscChannels.resize(nInitialFlavours);
  for (int i=0;i<nInitialFlavours;i++) {
    OscChannels[i].resize(nFinalFlavours);
  }
  OscChannels[0][0] = e_e;
  OscChannels[0][1] = e_m;
  OscChannels[0][2] = e_t;
  OscChannels[1][0] = m_e;
  OscChannels[1][1] = m_m;
  OscChannels[1][2] = m_t;

  OscChannels_Names.resize(nInitialFlavours);
  for (int i=0;i<nInitialFlavours;i++) {
    OscChannels_Names[i].resize(nFinalFlavours);
  }
  OscChannels_Names[0][0] = "e_e";
  OscChannels_Names[0][1] = "e_m";
  OscChannels_Names[0][2] = "e_t";
  OscChannels_Names[1][0] = "m_e";
  OscChannels_Names[1][1] = "m_m";
  OscChannels_Names[1][2] = "m_t";
}

void Oscillator::CheckEarthDensityFile() {
  // If the earth density file is null, set the default
  if (EarthDensityFile.Length() == 0) {
    std::string mach3 = std::getenv("MACH3");
    if (mach3.empty()) {
      std::cerr << "MACH3 environment needs to be set" << std::endl;
      throw;
    }
    EarthDensityFile = mach3+std::string("/CUDAProb3/models/PREM_4layer_cubic.dat");
    std::cout << "No specified Earth density file - Defaulting to: " << EarthDensityFile << "\n" << std::endl;
  }

}

void Oscillator::SetProductionHeightArray() {
  std::cout << "Loading Production Height Probabilities from: " << ProductionHeightFileName << std::endl;

  TFile* File = new TFile(ProductionHeightFileName);
  if (!File || File->IsZombie()) {
    std::cerr << "Can not find file:" << ProductionHeightFileName << std::endl;
    std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
    throw;
  }

  int nNuTypes = 2;
  int nNuFlav = 3;

  std::vector< std::vector<TString> > NuFlavNames(2);
  NuFlavNames[0].resize(3);
  NuFlavNames[1].resize(3);

  NuFlavNames[0][0] = "nue";
  NuFlavNames[0][1] = "numu";
  NuFlavNames[0][2] = "nutau";
  NuFlavNames[1][0] = "nuebar";
  NuFlavNames[1][1] = "numubar";
  NuFlavNames[1][2] = "nutaubar";

  std::vector< std::vector<TH3D*> > vecHist;
  vecHist.resize(nNuTypes);
  for (int iNuType=0;iNuType<nNuTypes;iNuType++) {
    vecHist[iNuType].resize(nNuFlav);
  }

  for (int iNuType=0;iNuType<nNuTypes;iNuType++) {
    for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
      TString HistName = "ProductionHeight_"+NuFlavNames[iNuType][iNuFlav];
      TH3D* Hist = (TH3D*)File->Get(HistName);

      if (!Hist) {
        std::cerr << HistName << " not found in File:" << ProductionHeightFileName << std::endl;
        File->ls();
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }

      vecHist[iNuType][iNuFlav] = Hist;

      if (nProductionHeightAveragingBins != vecHist[iNuType][iNuFlav]->GetNbinsZ()) {
        std::cerr << HistName << " has different number of Z bins:" << vecHist[iNuType][iNuFlav]->GetNbinsZ() << std::endl;
        std::cerr << "Expected:" << nProductionHeightAveragingBins << std::endl;
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }

      if (vecHist[iNuType][iNuFlav]->GetNbinsX() != nSecondaryBinsX) {
        std::cerr << HistName << " has different number of X bins:" << vecHist[iNuType][iNuFlav]->GetNbinsX() << std::endl;
        std::cerr << "Expected:" << nSecondaryBinsX << std::endl;
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }

      if (vecHist[iNuType][iNuFlav]->GetNbinsY() != nSecondaryBinsY) {
        std::cerr << HistName << " has different number of Y bins:" << vecHist[iNuType][iNuFlav]->GetNbinsY() << std::endl;
        std::cerr << "Expected:" << nSecondaryBinsY << std::endl;
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }
    }
  }

  int ProductionHeightProbabilitiesSize = nNuTypes * nNuFlav * nSecondaryBinsX * nSecondaryBinsY * nProductionHeightAveragingBins;
  std::vector<FLOAT_T> ProductionHeightProbabilities(ProductionHeightProbabilitiesSize);

  int index = 0;
  for (int iNuType=0;iNuType<nNuTypes;iNuType++) {
    for(int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
      for (int iSecondaryBinsX=0;iSecondaryBinsX<nSecondaryBinsX;iSecondaryBinsX++) {
        for (int iSecondaryBinsY=0;iSecondaryBinsY<nSecondaryBinsY;iSecondaryBinsY++) {
          double Total = 0.;

          for (int iProductionHeight=0;iProductionHeight<nProductionHeightAveragingBins;iProductionHeight++) {
            double dP_dh = vecHist[iNuType][iNuFlav]->GetBinContent(iSecondaryBinsX+1,iSecondaryBinsY+1,iProductionHeight+1);  //Taken from MTuple
            //double dP_dh = (1.0/(double)nProductionHeightAveragingBins)/vecHist[iNuType][iNuFlav]->GetZaxis()->GetBinWidth(iProductionHeight+1); //Flat Probability

            double dh = vecHist[iNuType][iNuFlav]->GetZaxis()->GetBinWidth(iProductionHeight+1);

            ProductionHeightProbabilities[index] = dP_dh * dh;
            Total += ProductionHeightProbabilities[index];

            index += 1;
          }

          if (fabs(Total-1.) > 1e-6) {
            std::cerr << "Probabilities integrated over production height do not sum to 1" << std::endl;
            std::cerr << "Total:" << Total << std::endl;
            for (int iProductionHeight=0;iProductionHeight<nProductionHeightAveragingBins;iProductionHeight++) {
              std::cout << "iProductionHeight:" << iProductionHeight << " | dP_dh:" << vecHist[iNuType][iNuFlav]->GetBinContent(iSecondaryBinsX+1,iSecondaryBinsY+1,iProductionHeight+1) << std::endl;
            }
            std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
            throw;
          }
        }
      }
    }
  }

  int ProductionHeightBinEdgesSize = nProductionHeightAveragingBins + 1;
  for (int iBin=1;iBin<=ProductionHeightBinEdgesSize;iBin++) {
    if (vecHist[0][0]->GetZaxis()->GetBinLowEdge(iBin) != ProductionHeightBinEdges[iBin-1]) {
      std::cerr << "Invalid Production height bin edges found in TH3Ds from File:" << ProductionHeightFileName << std::endl;
      std::cerr << "Expected:" << std::endl;
      for (int jBin=0;jBin<ProductionHeightBinEdgesSize;jBin++) {
        std::cerr << ProductionHeightBinEdges[jBin] << ", ";
      }
      std::cerr << std::endl;
      std::cerr << "Got:" << std::endl;
      for (int jBin=1;jBin<=ProductionHeightBinEdgesSize;jBin++) {
        std::cerr << vecHist[0][0]->GetZaxis()->GetBinLowEdge(jBin) << ", ";
      }
      std::cerr << std::endl;
      throw;
    }
  }

  File->Close();

  propagator->SetNumberOfProductionHeightBinsForAveraging(nProductionHeightAveragingBins);
  propagator->setProductionHeightList(ProductionHeightProbabilities,ProductionHeightBinEdges);
}

void Oscillator::SetProductionHeightBinEdges() {
  int ProductionHeightBinEdgesSize = nProductionHeightAveragingBins + 1;
  ProductionHeightBinEdges = std::vector<FLOAT_T>(ProductionHeightBinEdgesSize);

  for (int iBin=0;iBin<ProductionHeightBinEdgesSize;iBin++) {
    ProductionHeightBinEdges[iBin] = lProductionHeightRange + double(iBin)*(hProductionHeightRange-lProductionHeightRange)/double(nProductionHeightAveragingBins); //In km
  }

}

void Oscillator::PrintOscillatorConfig() {
  switch (ArrayConfig) {
    case 0:
      std::cout << "Using standard Oscillator configuration where a secondary bin only has one contribution" << std::endl;
      break;
    case 1:
      std::cout << "Using ManyContrib Oscillator configuration where a secondary bin has many contributions weighted by area" << std::endl;
      break;
    default:
      std::cout << "Unknown Oscillator Config:" << ArrayConfig << std::endl;
      throw;
  }
}

void Oscillator::SetOscillatorConfig(int ArrayConfig_) {
  bool ApplyChanges = false;
  if (ArrayConfig_ != ArrayConfig) {ApplyChanges = true;}

  ArrayConfig = ArrayConfig_;
  PrintOscillatorConfig();

  if (ApplyChanges) {
    DeleteArrays();
    ResizeArrays();
    FillArrays();
  }

}

void Oscillator::PrintBinning() {
  std::cout << "\n" << "Initialising oscillograms with ";
  if (IsLinear) {
    std::cout << "Linear Binning:" << std::endl;
  } else {
    std::cout << "Logarithmic Binning:" << std::endl;
  }

  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << std::setw(16) << "nCoarseCosz:" << " | " << std::setw(16) <<nCoarseCosz << std::endl;
  std::cout << std::setw(16) << "lCoarseCosz:" << " | " << std::setw(16) <<lCoarseCosz << std::endl;
  std::cout << std::setw(16) << "hCoarseCosz:" << " | " << std::setw(16) <<hCoarseCosz << std::endl;
  std::cout << std::endl;
  std::cout << std::setw(16) << "nCoarseEnergy:" << " | " << std::setw(16) <<nCoarseEnergy << std::endl;
  std::cout << std::setw(16) << "lCoarseEnergy:" << " | " << std::setw(16) <<lCoarseEnergy << std::endl;
  std::cout << std::setw(16) << "hCoarseEnergy:" << " | " << std::setw(16) <<hFineCosz << std::endl;
  std::cout << std::endl;
  std::cout << std::setw(16) << "nFineCosz:" << " | " << std::setw(16) <<nFineCosz << std::endl;
  std::cout << std::setw(16) << "nFineEnergy:" << " | " << std::setw(16) <<nFineEnergy << std::endl;
  std::cout << std::endl;
}

int Oscillator::GetOscillogramNBins(int Switcher) {

  switch (Switcher) {
    case 0:
      return nCoarseCosz;
    case 1:
      return nCoarseEnergy;
    case 2:
      if (useFineBinsPerBin) {
        return nFineCosz/nCoarseCosz;
      } else {
        return nFineCosz;
      }
    case 3:
      if (useFineBinsPerBin) {
        return nFineEnergy/nCoarseEnergy;
      } else {
        return nFineEnergy;
      }
    default:
      std::cerr << "Unknown option given to Ocillator::GetOscillogramNBins :" << Switcher << std::endl;
      throw;
  }

  return -1;
}

std::vector<double> Oscillator::linspace(double Emin, double Emax, int nDiv){
  if (nDiv==0) {
    std::cout << "Requested linear spacing distribution with 0 divisions" << std::endl;
    throw;
  }

  std::vector<double> linpoints(nDiv+1, 0.0);

  double step_lin = (Emax - Emin)/double(nDiv);

  double EE = Emin;

  for (int i=0; i<nDiv; i++) {
    if (fabs(EE)<1e-6) {EE = 0.;}

    linpoints[i] = EE;
    EE += step_lin;
  }

  linpoints[nDiv] = Emax;

  return linpoints;
}

std::vector<double> Oscillator::logspace(double Emin, double  Emax, int nDiv){
  if (nDiv==0) {
    std::cout << "Requested log spacing distribution with 0 divisions" << std::endl;
    throw;
  }

  std::vector<double> logpoints(nDiv+1, 0.0);
  logpoints[0]=Emin;

  if (Emin == 0.) {
    Emin = 0.01;
  }

  double Emin_log,Emax_log;
  Emin_log = log10(Emin);
  Emax_log = log10(Emax);

  double step_log = (Emax_log - Emin_log)/double(nDiv);

  double EE = Emin_log+step_log;

  for (int i=1; i<nDiv; i++) {
    logpoints[i] = pow(10.,EE);
    EE += step_log;
  }

  logpoints[nDiv]=Emax;

  return logpoints;
}

std::vector<double> Oscillator::ReturnFineBinningFromCoarseBinnnig(int nFine, std::vector<double> CoarseBinning) {

  std::vector<double> ReturnVec;
  int nCoarse = (int)CoarseBinning.size()-1;
  int nBinsPerCoarseBin = (nFine/nCoarse);

  for (int iCoarseBin=0;iCoarseBin<nCoarse;iCoarseBin++) {
    std::vector<double> tmpVec = linspace(CoarseBinning[iCoarseBin],CoarseBinning[iCoarseBin+1],nBinsPerCoarseBin);

    for (int iFineBin=0;iFineBin<nBinsPerCoarseBin;iFineBin++) {
      ReturnVec.push_back(tmpVec[iFineBin]);
    }
  }
  ReturnVec.push_back(CoarseBinning[nCoarse]);

  return ReturnVec;
}
*/
