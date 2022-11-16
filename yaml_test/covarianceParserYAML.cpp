#include "covarianceParserYAML.h"


covarianceParserYAML::covarianceParserYAML()
{

}

covarianceParserYAML::covarianceParserYAML(const char* filename)
{
   fYAMLDoc = YAML::LoadFile(filename);
   fNumPar = fYAMLDoc["parameters"].size();
   fPrior = new TVectorT<double>(fNumPar);
   fGenerated = new TVectorT<double>(fNumPar);
   fLB = new TVectorT<double>(fNumPar);
   fUB = new TVectorT<double>(fNumPar);
   fStepScale = new TVectorT<double>(fNumPar);
   fFlatPrior = new TVectorT<double>(fNumPar);
   fID = new TMatrixT<double>(fNumPar,2);
   fCovMatrix = new TMatrixT<double>(fNumPar,fNumPar);

   fFD_spline_names = new TObjArray(fNumPar);
   fFD_spline_modes = new TObjArray(fNumPar);

   fND_spline_names = new TObjArray(fNumPar);

   fNorm_modes = new TObjArray(fNumPar);
   fNorm_horncurrents = new TObjArray(fNumPar);
   fNorm_elem = new TObjArray(fNumPar);
   fNorm_nupdg = new TObjArray(fNumPar);
   fNorm_preoscnupdg = new TObjArray(fNumPar);

   fKinematicPars = new TObjArray(fNumPar);
   fKinematicBounds = new TObjArray(fNumPar);

   int i=0;
   TObjString* fd_spl_name= new TObjString("");
   TVectorT<double>* fd_modes = new TVectorT<double>(0);
   TObjString* nd_spl_name= new TObjString("");

   TVectorT<double>* modes = new TVectorT<double>(0);
   TVectorT<double>* elements = new TVectorT<double>(0);
   TVectorT<double>* horncurrents = new TVectorT<double>(0);
   TVectorT<double>* nupdg = new TVectorT<double>(0);
   TVectorT<double>* prodnupdg = new TVectorT<double>(0);
   TVectorT<double>* kinematicpar = new TVectorT<double>(0);
   TMatrixT<double>* kinematicbound = new TMatrixT<double>(0,0);

   std::vector<double> error(fNumPar);
   std::vector<std::map<std::string,double>> correlations(fNumPar);
   std::map<std::string, int> corrNames;

   for (auto const &param : fYAMLDoc["parameters"]) {
      std::cout << param["name"].as<std::string>() << std::endl;
      corrNames[param["name"].as<std::string>()]=i;
      fNames.push_back(param["name"].as<std::string>());
      (*fPrior)(i) = param["prior"].as<double>();
      (*fGenerated)(i) = param["generated"].as<double>();
      (*fLB)(i) = param["lowerbound"].as<double>();
      (*fUB)(i) = param["upperbound"].as<double>();
      (*fStepScale)(i) = param["stepscale"].as<double>();
      (*fID)(i,1) = param["detid"].as<int>();
      error[i]=param["error"].as<double>();

	  if (param["flatprior"]) {
		(*fFlatPrior)(i) = param["flatprior"].as<double>();
	  } else {
		(*fFlatPrior)(i) = 0;
	  }


	  //Now load in varaibles for spline systematics only
	  if (param["type"].as<std::string>() == "spline") {
		(*fID)(i,0) = param["splineind"].as<int>();
		*fd_spl_name = "";
		if (param["fd_spline_name"]) {
		  *fd_spl_name = param["fd_spline_name"].as<std::string>().c_str();
		}
		fFD_spline_names->AddAtAndExpand(fd_spl_name->Clone(),i);

		fd_modes->ResizeTo(0);
		if (param["fd_mode"]) {
		  std::vector<double> fd_modes_vec = param["fd_mode"].as<std::vector<double>>();
		  fd_modes->ResizeTo(fd_modes_vec.size());
		  for (std::size_t k = 0; k < fd_modes_vec.size(); k++) {
			(*fd_modes)(k)=fd_modes_vec[k];
		  }
		}
		fFD_spline_modes->AddAtAndExpand(fd_modes->Clone(),i);

		*nd_spl_name = "";
		if (param["nd_spline_name"]) {
		  *nd_spl_name = param["nd_spline_name"].as<std::string>().c_str();
		}
		fND_spline_names->AddAtAndExpand(nd_spl_name->Clone(),i);
	  } else if (param["type"].as<std::string>() == "norm") {
		(*fID)(i,0) = -1;
	  } else if (param["type"].as<std::string>() == "function") {
		(*fID)(i,0) = -2;
      } else {
		std::cout << "Wrong spline type!" << std::endl;
		exit(5);
	  }

      modes->ResizeTo(0);
	  if (param["mode"]) {
		if (!param["mode"].IsNull()) {
		  std::vector<double> modes_vec = param["mode"].as<std::vector<double>>();
		  modes->ResizeTo(modes_vec.size());
		  for (std::size_t k = 0; k < modes_vec.size(); k++) {
			(*modes)(k)=modes_vec[k];
		  }
		}
	  }
      fNorm_modes->AddAtAndExpand(modes->Clone(),i);

      elements->ResizeTo(0);
	  if (param["elements"]) {
		if(!param["elements"].IsNull()) {
		  std::vector<double> elem_vec = param["elements"].as<std::vector<double>>();
		  elements->ResizeTo(elem_vec.size());
		  for (std::size_t k = 0; k < elem_vec.size(); k++) {
			(*elements)(k)=elem_vec[k];
		  }
		}
	  }
      fNorm_elem->AddAtAndExpand(elements->Clone(),i);

      horncurrents->ResizeTo(0);
      if (param["horncurrents"]) {
		if(!param["horncurrents"].IsNull()) {
		  std::vector<double> hc_vec = param["horncurrents"].as<std::vector<double>>();
		  horncurrents->ResizeTo(hc_vec.size());
		  for(std::size_t k = 0; k < hc_vec.size(); k++) {
			(*horncurrents)(k)=hc_vec[k];
		  }
		}
      }
      fNorm_horncurrents->AddAtAndExpand(horncurrents->Clone(),i);

      nupdg->ResizeTo(0);
	  if (param["nupdg"]) {
		if(!param["nupdg"].IsNull()) {
		  std::vector<double> pdg_vec = param["nupdg"].as<std::vector<double>>();
		  nupdg->ResizeTo(pdg_vec.size());
		  for(std::size_t k = 0; k < pdg_vec.size(); k++) {
			(*nupdg)(k)=pdg_vec[k];
		  }
		}
	  }
      fNorm_nupdg->AddAtAndExpand(nupdg->Clone(),i);

      prodnupdg->ResizeTo(0);
	  if (param["prodnupdg"]) {
		if(!param["prodnupdg"].IsNull()) {
		  std::vector<double> prodnupdg_vec = param["prodnupdg"].as<std::vector<double>>();
		  prodnupdg->ResizeTo(prodnupdg_vec.size());
		  for (std::size_t k = 0; k < prodnupdg_vec.size(); k++) {
			(*prodnupdg)(k)=prodnupdg_vec[k];
		  }
		}
	  }
      fNorm_preoscnupdg->AddAtAndExpand(prodnupdg->Clone(),i);

      kinematicpar->ResizeTo(0);
      if (param["kinematicpars"]) {
		if(!param["kinematicpars"].IsNull()) {
	    std::vector<std::string> kinpar_vec = param["kinematicpars"].as<std::vector<std::string>>();
	    kinematicpar->ResizeTo(kinpar_vec.size());
	    for(std::size_t k = 0; k < kinpar_vec.size(); k++)
	       //(*kinematicpar)(k)=StringToKinematicVar(kinpar_vec[k]);
		   (*kinematicpar)(k) = -999;

		}
      }
      fKinematicPars->AddAtAndExpand(kinematicpar->Clone(),i);

	  kinematicbound->ResizeTo(0,0);
	  if(param["kinematicbounds"]) {
		if(!param["kinematicbounds"].IsNull()) {
		  std::vector<std::vector<double>> kinbound_vec = param["kinematicbounds"].as<std::vector<std::vector<double>>>();
		  kinematicbound->ResizeTo(kinbound_vec.size(),2);
		  for(std::size_t k = 0; k < kinbound_vec.size(); k++) {
			if(kinbound_vec[k].size()!=2) {
			  exit(5);
			}
			(*kinematicbound)(k,0)=kinbound_vec[k][0];
			(*kinematicbound)(k,1)=kinbound_vec[k][1];
		  }
		}
	  }
	  fKinematicBounds->AddAtAndExpand(kinematicbound->Clone(),i);
      
	  if(param["correlation"]) {
		for (auto const &corr : param["correlation"]) { 
		  correlations[i][corr.first.as<std::string>()]=corr.second.as<double>();
		}
	  }
	  i++;
   }

   for(int i=0; i<fNumPar; i++) {
	 (*fCovMatrix)(i,i)=error[i]*error[i];
	 for (auto const& [key, val] : correlations[i]) {
	   int index = -1;

	   if (corrNames.find(key) != corrNames.end()) {
		 index=corrNames[key];
	   }
	   else {
		 std::cout << "Parameter " << key << " not in list! Check your spelling?" << std::endl;
		 exit(5);
	   }

	   double corr1 = val;
	   double corr2 = 0;
	   if(correlations[index].find(fNames[i]) != correlations[index].end()) {
		 corr2 = correlations[index][fNames[i]];
		 if(corr2 != corr1) {
		   std::cout << "Correlations are not equal between " << fNames[i] << " and " << key << std::endl;
		   exit(5);
		 }
	   }
	   else {
		 std::cout << "Correlation does not appear reciprocally between " << fNames[i] << " and " << key << std::endl;
		 exit(5);
	   }

	   (*fCovMatrix)(i,index)= (*fCovMatrix)(index,i) = corr1*error[i]*error[index];
	 }
   }
}

covarianceParserYAML::~covarianceParserYAML() {

}
