#include "stretch.h"

// uncomment below to turn on multithreading
//#define MULTITHREAD

stretch::stretch(bool verbose)
{
  random = new TRandom3(0);

  init("output.root", verbose);
  nwalk=100;
}


stretch::stretch(const char *outfile, int nwalkers, bool verbose)
{
   random = new TRandom3(0);
   nwalk=nwalkers;
   
   init(outfile, verbose);
}

void stretch::init(std::string outfile,  bool verbose)
{
  osc_only = false; // default value
//  accCount = 0;
  auto_save = 100; // default to saving every 100 steps
  save_nominal = true;
//  init_pos = false; //starting  parameters should be thrown 
  osc=oscbar=-1;
  N=0;//zero dimensions
  a=2;//tuning parameter default value

  // setup output file
  outputFile = new TFile(outfile.c_str(), "RECREATE");
//  std::cout << gDirectory->pwd() << std::endl;
  // setup output tree
  outTree = new TTree("posteriors","Posterior Distributions");
//  outputFile->ls();

  // fit summary and debug info
  debug = verbose;

  if(debug)
    {
      outfile += ".log";
      //char *logname = ".log";
      debugInterval = 50;
      //strncat(outfile, logname, 4);
      debugFile.open(outfile.c_str());
    }

  //clear samples, systs
  samples.clear();
  systematics.clear();
}

stretch::~stretch()
{
  if(debug)
    debugFile.close();
  
  outputFile->Close();
}

void stretch::setChainLength(int L)
{
  chainLength = L;
}

void stretch::addSamplePDF(samplePDFBase *sample)
{
  samples.push_back(sample);
}

void stretch::addSystObj(covarianceBase *cov, bool isOsc)
{
  systematics.push_back(cov);
  N+=cov->getSize();

  if(save_nominal)
    {
      std::vector<double> vec = cov->getNominalArray();
      int n = vec.size();
      double *n_vec = new double[n];
      for (int i = 0; i < n; ++i)
	{
	  n_vec[i] = vec[i];
	}
      TVectorT<double> t_vec(n, n_vec);
      TString nameof = TString(cov->getName());
      nameof = nameof.Append("_nom");
      t_vec.Write(nameof);
    }
  
  std::vector< std::vector < double > > walkp;

  for(int i=0; i<cov->getSize(); i++)
  {
     std::vector<double> ww;
     for(int j=0; j<nwalk; j++)
	ww.push_back(cov->getParProp(i));
     walkp.push_back(ww);
  }

  currentpar.push_back(walkp);
  proposedpar.push_back(walkp);

  if(isOsc && osc==-1)
     osc=systematics.size()-1;
  else if(isOsc && oscbar==-1)
     oscbar=systematics.size()-1;  

}

// bool stretch::accept()
// {


//   return accept;
// }

void stretch::runStretch()
{
  // some checks
//  if(!samples.size())
//    {std::cerr << "no samples! Stopping MCMC" << std::endl; return;}


  //  if(debug)
  //     std::cout << "Passed initial checks " << osc << std::endl;

  // prepare output tree
   for(size_t i = 0; i < systematics.size(); i++)
   {
      for(size_t j=0; j<currentpar[i].size(); j++)
      {
	 outTree->Branch(systematics[i]->GetParName(j).c_str(),"std::vector<double>",&currentpar[i][j]);
      }
   }
   
//  outTree->Branch("LogL",&logLCurr,"LogL/D");
//  outTree->Branch("accProb", &accProb,"accProb/D");
  outTree->Branch("step", &step, "step/I");
  outTree->Branch("LogL","std::vector<double>",&logLCurr);

  // store individual likelihood components
//  sample_llh = new float[samples.size()];
//  syst_llh = new float[systematics.size()];
    
//   for (int i =0; i < samples.size(); ++i)
//     {
//       std::stringstream oss, oss2;
//       oss << "LogL_sample_" << i;
//       oss2 << oss.str() << "/F";
//       outTree->Branch(oss.str().c_str(), &sample_llh[i], oss2.str().c_str());

//       // For adding sample dependent branches to the posteriors tree
//       samples[i]->setMCMCBranches(outTree);
//     }

//   for (int i = 0; i < systematics.size(); ++i)
//     {
//       std::stringstream oss, oss2;
//       oss << "LogL_systematic_" << systematics[i]->getName();
//       oss2 << oss.str() << "/F";
//       outTree->Branch(oss.str().c_str(), &syst_llh[i], oss2.str().c_str());
//     }

  std::cout << "\n---------------------starting MCMC---------------------" << std::endl;
//  outputFile->cd();
//   if(debug)
//     {
//       printInitialState();
//       debugFile << "-----starting MCMC-----" << std::endl;
//     }

  // initial reweight and likelihood calc
  std::vector<double> llh_init;


  for(int k=0; k<nwalk; k++)
  {
     std::cout << k << std::endl;
     llh_init.push_back(0);

     for (size_t i = 0; i < systematics.size(); i++)
     {
	systematics[i]->throwNominal(false);
	std::vector<double> startpar = systematics[i]->getNominalArray();
	systematics[i]->setParameters(startpar);

	for(int j=0; j<systematics[i]->getSize(); j++)
	   currentpar[i][j][k]=systematics[i]->getParProp(j);
	   
	systematics[i]->throwNominal();
//	systematics[i]->printPars();

	llh_init[k] += systematics[i]->GetLikelihood();
	if(debug)
	   debugFile << "LLH after " << systematics[i]->getName() << " " << llh_init[k] << std::endl;
     }
     
     for (size_t i = 0; i < samples.size(); i++)
     {
	if (osc!=-1)
	{
	   samples[i]->reweight();//((covarianceOsc*)systematics[osc])->getPropPars());
	}
	else
	{
	   //double* fake = NULL;
	   samples[i]->reweight();
	}
	
	llh_init[k] += samples[i]->GetLikelihood();
	
	if(debug)
	   debugFile << "LLH after sample " << i << " " << llh_init[k] << std::endl;
     }
     
     logLProp.push_back(llh_init[k]);
     logLCurr.push_back(llh_init[k]);
  }
  outTree->Fill();
  
  // begin!
  clock.Start();
  for(step = 0; step < chainLength; step++)
  {
     if(step%10==0)
	std::cout << "At step: " << step << std::endl;
  
     for(int nw=0; nw<nwalk; nw++)
     {
	logLProp[nw]=0;
	bool done=false;
	int ow;
	while(!done)
	{
	   ow=random->Integer(nwalk-1);
	   if(nw!=ow)
	      done=true;
	}
	double z=pow((a-1)*random->Uniform()+1,2)/a;
	for(size_t i=0; i<systematics.size(); i++)
	{
	   std::vector<double> pars(currentpar[i].size());
	   for(int j=0; j<systematics[i]->getSize(); j++)
	   {
	      proposedpar[i][j][nw] = currentpar[i][j][ow]+z*(currentpar[i][j][nw]-currentpar[i][j][ow]);
	      if(!systematics[i]->isParameterFixed(j))
		 pars[j] = proposedpar[i][j][nw];
	      else 
		 pars[j]=currentpar[i][j][nw];
	      // May have to put something in here for special cases, e.g. mass hierarchy
	   }
	   systematics[i]->setParameters(pars);
//	   systematics[i]->printPars();
	   logLProp[nw] += systematics[i]->GetLikelihood();
	}
	for (size_t i = 0; i < samples.size(); i++)
	{
	   if (osc!=-1)
	   {
	      samples[i]->reweight();//((covarianceOsc*)systematics[osc])->getPropPars());
	   }
	   else
	   {
	      //double* fake = NULL;
	      samples[i]->reweight();
	   }
	   
	   logLProp[nw] += samples[i]->GetLikelihood();
	}
	double q = pow(z,N-1)*exp(logLCurr[nw]-logLProp[nw]);
//	std::cout << nw << "," << ow << ": " << logLCurr[nw] << "," << logLProp[nw] << "," << z << "," << q << std::endl;
	if(random->Uniform(0,1)>q)
	{
	   for(size_t i=0; i<systematics.size(); i++)
	      for(int j=0; j<systematics[i]->getSize(); j++)
		 proposedpar[i][j][nw]=currentpar[i][j][nw];
	   logLProp[nw]=logLCurr[nw];
	}

     }
     for(int nw=0; nw<nwalk; nw++)
     {
	 for(size_t i=0; i<systematics.size(); i++)
	      for(int j=0; j<systematics[i]->getSize(); j++)
		 currentpar[i][j][nw]=proposedpar[i][j][nw];
	 logLCurr[nw]=logLProp[nw];
     }
     outTree->Fill();

      // auto save 
     if (step % auto_save == auto_save - 1)
     {
	//	  if (debug)
//	    std::cout << "Current Acceptance Rate: " << step / accept << std::endl;
	outTree->AutoSave();
     }
    }
  // finished MCMC //////////////////////////////
  
  std::cout << "----------- finished MCMC --------------" << std::endl;
  clock.Stop();
  outTree->Write();
  
}
