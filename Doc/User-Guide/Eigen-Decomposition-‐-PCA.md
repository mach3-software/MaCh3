# Eigen-Decomposition (PCA) {#eigen-decomposition}

PCA or Eigen Decomposition is turned off by default. PCA is being enabled based on the cosntuctor. This is default no PCA consutvcosntuctorotr.
```
std::vector<std::string> xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
xsec = new covarianceXsec(xsecCovMatrixFile);
```

To enable it it should look like it
```
std::vector<std::string> xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
double PCAthreshold = 0.00001;
xsec = new covarianceXsec(xsecCovMatrixFile, PCAthreshold);
```
This threshold indicates which eigen value to remove and reduce fit dimensionality. If you set the threhsold very low you may not remove anything just to run with the decomposed matrix.

<img width="350" src="https://github.com/mach3-software/MaCh3/assets/45295406/ae078a7a-724b-4297-9b3d-e941d1dedfb1">

It is also possible to decompose only part of the matrix
```
std::vector<std::string> xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
double PCAthreshold = 0.00001;
int FirstPCAdpar = 8;
int LastPCAdpar = 108
xsec = new covarianceXsec(xsecCovMatrixFile, PCAthreshold, FirstPCAdpar, LastPCAdpar)  ;
```
This way parameters through 8-108 will be decomposed while 0-8 will be in undecomposed base. You can see transition matrix between normal and eigen base below:

<img width="350" src="https://github.com/mach3-software/MaCh3/assets/45295406/404646cc-98d1-4488-b1df-372e5ce13a26">

