// -*- c++ -*-
//
// GPU versions of mosc functions
//

#include "stdio.h"

typedef double dmArray[3];
typedef double mixArray[3][2];
__constant__ dmArray dm_device[9];
__constant__ mixArray mix_device[18];

#define mix_im_sign_t -1
#define mix_im_sign_f 1



// ERROR CHECKING ///////////////////////////////////////////

//#define CUDA_ERROR_CHECK // turn this on and off to disable error checking

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err )
    {
      fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
	       file, line, cudaGetErrorString( err ) );
      exit( -1 );
    }
#endif

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err )
    {
      fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
	       file, line, cudaGetErrorString( err ) );
      exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
      {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
      }
#endif

      return;
}

/////////////////////////////////////////////////////////////

#define elec (0)
#define muon (1)
#define tau  (2)
#define re (0)
#define im (1)

typedef enum nu_type {
  data_type,
  nue_type,
  numu_type,
  nutau_type,
  sterile_type,
  unknown_type} NuType;

typedef enum matrix_type {
  standard_type,
  barger_type} MatrixType;

//#define ZERO_CP
static int matrixtype = standard_type;

/* Flag to tell us if we're doing nu_e or nu_sterile matter effects */
// CWRET: comment out to rid compiler warning (never used)
//static NuType matterFlavor = nue_type;
//static double putMix[3][3][2]; ****************************

/* 2*sqrt(2)*Gfermi in (eV^2-cm^3)/(mole-GeV) - for e<->[mu,tau] */
//static const double tworttwoGf = 1.52588e-4;

// CWRET: comment out to rid compiler warning (never used)
/*
__host__ void setMatterFlavor(int flavor)
{
  if (flavor == nue_type) matterFlavor = nue_type;
  else if (flavor == sterile_type) matterFlavor = sterile_type;
  else {
    //fprintf(stderr, "setMatterFlavor: flavor=%d", flavor);
    //moscerr("setMatterFlavor: Illegal flavor.");
  }
}
*/

__host__ void setmix_sin(double s12,double s23,double s13,double dcp, double Mix[3][3][2])
{
  double c12,c23,c13,sd,cd = 1.0;

  if ( s12>1.0 ) s12=1.0;
  if ( s23>1.0 ) s23=1.0;
  if ( s13>1.0 ) s13=1.0;
  //if ( cd >1.0 ) cd =1.0; this is not needed

  sd = sin( dcp );
  cd = cos( dcp );

  c12 = sqrt(1.0-s12*s12);
  c23 = sqrt(1.0-s23*s23);
  c13 = sqrt(1.0-s13*s13);

  if ( matrixtype == standard_type )
    {
      Mix[0][0][re] =  c12*c13;
      Mix[0][0][im] =  0.0;
      Mix[0][1][re] =  s12*c13;
      Mix[0][1][im] =  0.0;
      Mix[0][2][re] =  s13*cd;
      Mix[0][2][im] = -s13*sd;
      Mix[1][0][re] = -s12*c23-c12*s23*s13*cd;
      Mix[1][0][im] =         -c12*s23*s13*sd;
      Mix[1][1][re] =  c12*c23-s12*s23*s13*cd;
      Mix[1][1][im] =         -s12*s23*s13*sd;
      Mix[1][2][re] =  s23*c13;
      Mix[1][2][im] =  0.0;
      Mix[2][0][re] =  s12*s23-c12*c23*s13*cd;
      Mix[2][0][im] =         -c12*c23*s13*sd;
      Mix[2][1][re] = -c12*s23-s12*c23*s13*cd;
      Mix[2][1][im] =         -s12*c23*s13*sd;
      Mix[2][2][re] =  c23*c13;
      Mix[2][2][im] =  0.0;
    }
  else
    {
      Mix[0][0][re] =  c12;
      Mix[0][0][im] =  0.0;
      Mix[0][1][re] =  s12*c23;
      Mix[0][1][im] =  0.0;
      Mix[0][2][re] =  s12*s23;
      Mix[0][2][im] =  0.0;
      Mix[1][0][re] = -s12*c13;
      Mix[1][0][im] =  0.0;
      Mix[1][1][re] =  c12*c13*c23+s13*s23*cd;
      Mix[1][1][im] =              s13*s23*sd;
      Mix[1][2][re] =  c12*c13*s23-s13*c23*cd;
      Mix[1][2][im] =             -s13*c23*sd;
      Mix[2][0][re] = -s12*s13;
      Mix[2][0][im] =  0.0;
      Mix[2][1][re] =  c12*s13*c23-c13*s23*cd;
      Mix[2][1][im] =             -c13*s23*sd;
      Mix[2][2][re] =  c12*s13*s23+c13*c23*cd;
      Mix[2][2][im] =              c13*c23*sd;
    }
}

__host__ void setmass(double dms21, double dms23, double dmVacVac[][3])
{
  double delta=5.0e-9;
  double mVac[3];

  mVac[0] = 0.0;
  mVac[1] = dms21;
  mVac[2] = dms21+dms23;

  /* Break any degeneracies */
  if (dms21==0.0) mVac[0] -= delta;
  if (dms23==0.0) mVac[2] += delta;

  dmVacVac[0][0] = dmVacVac[1][1] = dmVacVac[2][2] = 0.0;
  dmVacVac[0][1] = mVac[0]-mVac[1]; dmVacVac[1][0] = -dmVacVac[0][1];
  dmVacVac[0][2] = mVac[0]-mVac[2]; dmVacVac[2][0] = -dmVacVac[0][2];
  dmVacVac[1][2] = mVac[1]-mVac[2]; dmVacVac[2][1] = -dmVacVac[1][2];
}

/// onwards are for matter effects calcs

__device__ void get_product(double L, double E, double rho,
		 double Mix[][3][2], double dmMatVac[][3], double dmMatMat[][3],
		 int antitype, double product[][3][3][2])
{
  double fac=0.0;
  double twoEHmM[3][3][3][2];
  double tworttwoGf = 1.52588e-4;

  /* (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km) */
  /* Reverse the sign of the potential depending on neutrino type */

  //if (matterFlavor == nue_type) {

    /* If we're doing matter effects for electron neutrinos */
  if (antitype<0) fac =  tworttwoGf*E*rho; /* Anti-neutrinos */
  else        fac = -tworttwoGf*E*rho; /* Real-neutrinos */
//  }

  //printf("gpu product fac = %f\n", fac);

/*
  else if (matterFlavor == sterile_type) {

    // If we're doing matter effects for sterile neutrinos
    if (antitype<0) fac = -0.5*tworttwoGf*E*rho; // Anti-neutrinos
    else        fac =  0.5*tworttwoGf*E*rho; // Real-neutrinos
  } */

  /*printf("gpu realpart");
  printf("%f %f %f \n", Mix[0][0][re], Mix[1][0][re], Mix[2][0][re]);
  printf("%f %f %f \n", Mix[0][1][re], Mix[1][1][re], Mix[2][1][re]);
  printf("%f %f %f \n", Mix[0][2][re], Mix[1][2][re], Mix[2][2][re]);
  printf("gpu imgpart");
  printf("%f %f %f \n", Mix[0][0][im], Mix[1][0][im], Mix[2][0][im]);
  printf("%f %f %f \n", Mix[0][1][im], Mix[1][1][im], Mix[2][1][im]);
  printf("%f %f %f \n", Mix[0][2][im], Mix[1][2][im], Mix[2][2][im]);*/

  int anti_ = (antitype < 0 ? mix_im_sign_t : mix_im_sign_f);

  /* Calculate the matrix 2EH-M_j */
  for (int n=0; n<3; n++)
    {
      for (int m=0; m<3; m++)
	{
	  //#ifndef ZERO_CP
      twoEHmM[n][m][0][re] =
	-fac*(Mix[0][n][re]*Mix[0][m][re]+(anti_*Mix[0][n][im])*(anti_*Mix[0][m][im]));
      twoEHmM[n][m][0][im] =
        -fac*(Mix[0][n][re]*(anti_*Mix[0][m][im])-(anti_*Mix[0][n][im])*Mix[0][m][re]);
      twoEHmM[n][m][1][re] = twoEHmM[n][m][2][re] = twoEHmM[n][m][0][re];
      twoEHmM[n][m][1][im] = twoEHmM[n][m][2][im] = twoEHmM[n][m][0][im];
      /*#else
      twoEHmM[n][m][0][re] =
        -fac*(Mix[0][n][re]*Mix[0][m][re]);
      printf("g%i %i %f %f %f\n", n, m, fac, Mix[0][n][re], Mix[0][m][re]);
      twoEHmM[n][m][0][im] = 0 ;
      twoEHmM[n][m][1][re] = twoEHmM[n][m][2][re] = twoEHmM[n][m][0][re];
      twoEHmM[n][m][1][im] = twoEHmM[n][m][2][im] = twoEHmM[n][m][0][im];
      //#endif*/
      //printf("gpuMix[][][re] = %f , gpuMix[][][im] = %f\n", Mix[n][m][re], Mix[n][m][im]);

      if (n==m) for (int j=0; j<3; j++)
		  {
		    twoEHmM[n][m][j][re] -= dmMatVac[j][n];
		    //printf("gpu -vac %f \n", dmMatVac[j][n]);
		  }

      //if (n== 2 && m ==2) printf("twoEHmMgpu %f \n", twoEHmM[n][m][0][re]);

    }
  }

/* Calculate the product in eq.(10) of twoEHmM for j!=k */
//cudaMemset(product, 0, 3*3*3*2*sizeof(double));
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
	{
	  product[i][j][k][re] = 0;
	  product[i][j][k][im] = 0;
	}


for (int i=0; i<3; i++) {
  for (int j=0; j<3; j++) {
    for (int k=0; k<3; k++) {

      //#ifndef ZERO_CP
      product[i][j][0][re] +=
	twoEHmM[i][k][1][re]*twoEHmM[k][j][2][re] -
	twoEHmM[i][k][1][im]*twoEHmM[k][j][2][im];
      product[i][j][0][im] +=
	twoEHmM[i][k][1][re]*twoEHmM[k][j][2][im] +
	twoEHmM[i][k][1][im]*twoEHmM[k][j][2][re];
      product[i][j][1][re] +=
	twoEHmM[i][k][2][re]*twoEHmM[k][j][0][re] -
          twoEHmM[i][k][2][im]*twoEHmM[k][j][0][im];
        product[i][j][1][im] +=
          twoEHmM[i][k][2][re]*twoEHmM[k][j][0][im] +
          twoEHmM[i][k][2][im]*twoEHmM[k][j][0][re];
	product[i][j][2][re] +=
          twoEHmM[i][k][0][re]*twoEHmM[k][j][1][re] -
          twoEHmM[i][k][0][im]*twoEHmM[k][j][1][im];
        product[i][j][2][im] +=
          twoEHmM[i][k][0][re]*twoEHmM[k][j][1][im] +
          twoEHmM[i][k][0][im]*twoEHmM[k][j][1][re];
	
	/*#else
	product[i][j][0][re] +=
	  twoEHmM[i][k][1][re]*twoEHmM[k][j][2][re];
        product[i][j][1][re] +=
          twoEHmM[i][k][2][re]*twoEHmM[k][j][0][re];
        product[i][j][2][re] +=
          twoEHmM[i][k][0][re]*twoEHmM[k][j][1][re];
	  //#endif                  */
      }
    //#ifndef ZERO_CP
    product[i][j][0][re] /= (dmMatMat[0][1]*dmMatMat[0][2]);
      product[i][j][0][im] /= (dmMatMat[0][1]*dmMatMat[0][2]);
      product[i][j][1][re] /= (dmMatMat[1][2]*dmMatMat[1][0]);
      product[i][j][1][im] /= (dmMatMat[1][2]*dmMatMat[1][0]);
      product[i][j][2][re] /= (dmMatMat[2][0]*dmMatMat[2][1]);
      product[i][j][2][im] /= (dmMatMat[2][0]*dmMatMat[2][1]);

      /*#else
      product[i][j][0][re] /= (dmMatMat[0][1]*dmMatMat[0][2]);
      product[i][j][1][re] /= (dmMatMat[1][2]*dmMatMat[1][0]);
      product[i][j][2][re] /= (dmMatMat[2][0]*dmMatMat[2][1]);
      //#endif*/
    }
  }
}
/***********************************************************************
  getM
  Compute the matter-mass vector M, dM = M_i-M_j and
  and dMimj. type<0 means anti-neutrinos type>0 means "real" neutrinos
***********************************************************************/

      __device__ void getM(double Enu, double rho,
          double Mix[][3][2], double dmVacVac[][3], int antitype,
          double dmMatMat[][3], double dmMatVac[][3])
{
  int i, j, k;
  double alpha, beta, gamma, fac=0.0, arg, tmp;
  double alphaV, betaV, gammaV, argV, tmpV;
  double theta0, theta1, theta2;
  double theta0V, theta1V, theta2V;
  double mMatU[3], mMatV[3], mMat[3];
  double tworttwoGf = 1.52588e-4;

  /* Equations (22) fro Barger et.al.*/
  /* Reverse the sign of the potential depending on neutrino type */
  //if (matterFlavor == nue_type) {
  /* If we're doing matter effects for electron neutrinos */
  if (antitype<0) fac =  tworttwoGf*Enu*rho; /* Anti-neutrinos */
  else        fac = -tworttwoGf*Enu*rho; /* Real-neutrinos */
  //}
  //printf("GPU fac = %f = %f %f %f\n", fac, tworttwoGf, Enu, rho);

  //else if (matterFlavor == sterile_type) {
  /* If we're doing matter effects for sterile neutrinos */
  //if (antitype<0) fac = -0.5*tworttwoGf*Enu*rho; /* Anti-neutrinos */

  //   else        fac =  0.5*tworttwoGf*Enu*rho; /* Real-neutrinos */
  // }
  /* The strategy to sort out the three roots is to compute the vacuum
   * mass the same way as the "matter" masses are computed then to sort
   * the results according to the input vacuum masses
   */

  // if we are doing anti-nu, the imaginary part is multiplied by -1
  int anti_ = (antitype < 0 ? mix_im_sign_t : mix_im_sign_f);

  alpha  = fac + dmVacVac[0][1] + dmVacVac[0][2];
  alphaV = dmVacVac[0][1] + dmVacVac[0][2];

  //#ifndef ZERO_CP
  beta = dmVacVac[0][1]*dmVacVac[0][2] +
      fac*(dmVacVac[0][1]*(1.0 -
                           Mix[elec][1][re]*Mix[elec][1][re] -
                           (Mix[elec][1][im])*(Mix[elec][1][im])) +
           dmVacVac[0][2]*(1.0-
                           Mix[elec][2][re]*Mix[elec][2][re] -
                           (anti_*Mix[elec][2][im])*(anti_*Mix[elec][2][im])));
  betaV = dmVacVac[0][1]*dmVacVac[0][2];

  /*#else
 beta = dmVacVac[0][1]*dmVacVac[0][2] +
    fac*(dmVacVac[0][1]*(1.0 -
                         Mix[elec][1][re]*Mix[elec][1][re]) +
         dmVacVac[0][2]*(1.0-
                         Mix[elec][2][re]*Mix[elec][2][re]));
  betaV = dmVacVac[0][1]*dmVacVac[0][2];
  #endif*/

  //#ifndef ZERO_CP
  gamma = fac*dmVacVac[0][1]*dmVacVac[0][2]*
    (Mix[elec][0][re]*Mix[elec][0][re]+(Mix[elec][0][im])*(Mix[elec][0][im]));
  gammaV = 0.0;
  /*  #else
  gamma = fac*dmVacVac[0][1]*dmVacVac[0][2]*
    (Mix[elec][0][re]*Mix[elec][0][re]);
  gammaV = 0.0;
  //#endif */

  /* Compute the argument of the arc-cosine */
  tmp = alpha*alpha-3.0*beta;
  tmpV = alphaV*alphaV-3.0*betaV;
  if (tmp<0.0) {
   // fprintf(stderr, "getM: alpha^2-3*beta < 0 !\n");
    tmp = 0.0;
  }

  //printf("GPU - fac %f alpha %f beta %f gamma %f\n", fac, alpha, beta, gamma);

 /* Equation (21) */
  arg = (2.0*alpha*alpha*alpha-9.0*alpha*beta+27.0*gamma)/
    (2.0*sqrt(tmp*tmp*tmp));
  if (fabs(arg)>1.0) arg = arg/fabs(arg);
  argV = (2.0*alphaV*alphaV*alphaV-9.0*alphaV*betaV+27.0*gammaV)/
    (2.0*sqrt(tmpV*tmpV*tmpV));
  if (fabs(argV)>1.0) argV = argV/fabs(argV);

  /* These are the three roots the paper refers to */
  theta0 = acos(arg)/3.0;
  theta1 = theta0-(2.0*M_PI/3.0);
  theta2 = theta0+(2.0*M_PI/3.0);

  theta0V = acos(argV)/3.0;
  theta1V = theta0V-(2.0*M_PI/3.0);
  theta2V = theta0V+(2.0*M_PI/3.0);

  mMatU[0] = mMatU[1] = mMatU[2] = -(2.0/3.0)*sqrt(tmp);
  mMatU[0] *= cos(theta0); mMatU[1] *= cos(theta1); mMatU[2] *= cos(theta2);

  tmp = dmVacVac[0][0] - alpha/3.0;
  mMatU[0] += tmp; mMatU[1] += tmp; mMatU[2] += tmp;
  mMatV[0] = mMatV[1] = mMatV[2] = -(2.0/3.0)*sqrt(tmpV);
  mMatV[0] *= cos(theta0V); mMatV[1] *= cos(theta1V); mMatV[2] *= cos(theta2V);
  tmpV = dmVacVac[0][0] - alphaV/3.0;

  mMatV[0] += tmpV; mMatV[1] += tmpV; mMatV[2] += tmpV;

  /* Sort according to which reproduce the vaccum eigenstates */
  for (i=0; i<3; i++) {
    tmpV = fabs(dmVacVac[i][0]-mMatV[0]);
    k = 0;
    for (j=1; j<3; j++) {
      tmp = fabs(dmVacVac[i][0]-mMatV[j]);
      if (tmp<tmpV) {
        k = j;
        tmpV = tmp;
      }
    }
    mMat[i] = mMatU[k];
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      dmMatMat[i][j] = mMat[i] - mMat[j];
      dmMatVac[i][j] = mMat[i] - dmVacVac[j][0];
      //printf("gpu mmat %f\n", mMat[i]);
      //printf("%i %i %f %f\n",i,j,dmMatVac[i][j], dmMatMat[i][j]);
      //if (i == 2 && j == 2) printf("hellogpu %f %f\n", dmMatMat[i][j], dmMatVac[i][j]);
    }
  }
}

/***********************************************************************
  getA
  Calculate the transition amplitude matrix A (equation 10)
***********************************************************************/

__device__ void getA(double L, double E, double rho,
          double Mix[][3][2], double dmMatVac[][3], double dmMatMat[][3],
          int antitype, double A[3][3][2], double phase_offset)
{
  //int n, m, i, j, k;
  double /*fac=0.0,*/ arg, c, s;

  double X[3][3][2];
  double product[3][3][3][2];
  /* (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km) */
  const double LoEfac = 2.534;

  if ( phase_offset==0.0 )
    {
      //printf("%f %f\n", L, E);
      get_product(L, E, rho, Mix, dmMatVac, dmMatMat, antitype, product);
    }

  /* Make the sum with the exponential factor */
  //cudaMemset(X, 0, 3*3*2*sizeof(double));
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
	X[i][j][k] = 0;

  for (int k=0; k<3; k++)
    {
      arg = -LoEfac*dmMatVac[k][0]*L/E;
      if ( k==2 ) arg += phase_offset ;
      c = cos(arg);
      s = sin(arg);
      for (int i=0; i<3; i++)
	{
	  for (int j=0; j<3; j++)
	    {
	      ///#ifndef ZERO_CP
	      X[i][j][re] += c*product[i][j][k][re] - s*product[i][j][k][im];
	      X[i][j][im] += c*product[i][j][k][im] + s*product[i][j][k][re];
	      /*#else
	      X[i][j][re] += c*product[i][j][k][re];
	      X[i][j][im] += s*product[i][j][k][re];
	      #endif                                                                                      */
	    }
	}
    }
  //  printf("\n testy %f %f ",X[0][0][im], Mix[0][0][im]);
  /* Compute the product with the mixing matrices */
  //cudaMemset(A, 0, 3*3*2*sizeof(double));

int anti_ = (antitype < 0 ? mix_im_sign_t : mix_im_sign_f);

  for(int i=0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 2; ++k)
	{
	  //printf("gpu X = %f\n", X[i][j][k]);
	  A[i][j][k] = 0;
	}

  for (int n=0; n<3; n++) {
    for (int m=0; m<3; m++) {
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {                                                                                                                      	
	  //#ifndef ZERO_CP
          A[n][m][re] +=
	    Mix[n][i][re]*X[i][j][re]*Mix[m][j][re] +
	    (Mix[n][i][re])*X[i][j][im]*(anti_*Mix[m][j][im]) +
            (anti_*Mix[n][i][im])*X[i][j][re]*(anti_*Mix[m][j][im]) -
	    (anti_*Mix[n][i][im])*X[i][j][im]*Mix[m][j][re];
	  //printf("gpu regret %f %f %f \n",Mix[n][i][re], X[i][j][im], Mix[m][j][im]);
	  A[n][m][im] +=
	    (anti_*Mix[n][i][im])*X[i][j][im]*(anti_*Mix[m][j][im]) +
            (anti_*Mix[n][i][im])*X[i][j][re]*Mix[m][j][re] +
            Mix[n][i][re]*X[i][j][im]*Mix[m][j][re] -
            Mix[n][i][re]*X[i][j][re]*(anti_*Mix[m][j][im]);
	  /*#else
          A[n][m][re] +=
            Mix[n][i][re]*X[i][j][re]*Mix[m][j][re];
          A[n][m][im] +=
            Mix[n][i][re]*X[i][j][im]*Mix[m][j][re];
	    #endif      */
	  //printf("gpu mix %f\n", Mix[m][j][re]);
	  //printf("\n gpu %i %i %i A %f", n, m, re, A[n][m][re]);
        }
      }
    }
  }
}

////#include "mosc.cu"


static double dm[3][3];
static double mix[3][3][2];
//static double Ain[3][3][2];
static double dm21,dm32,s12,s23,s31,cd;

extern "C" __host__ double getMixVal(int x, int y, int z)
{
  return mix[x][y][z];
}

extern "C" __host__ double getT13()
{
	return dm[1][1];
}
	
__host__ void init_mixing_matrix(double dm21f,double dm32f,double s12f,double s23f,double s31f,double cdf)
{
  dm21=dm21f ;  dm32=dm32f ;
  s12=s12f   ;  s23=s23f   ; s31=s31f ;
  cd=cdf;	
// CWRET: comment out to rid compiler warning (never used)
  //setMatterFlavor(nue_type);
  setmix_sin(s12,s23,s31,cd,mix);
  setmass(dm21,dm32,dm);	
  //  cudaMalloc((void **) &device_array, size);
  //cudaMalloc((void **) &Ain,3*3*2*sizeof(double));
  //Ain[0][0][re] = Ain[1][1][re]		= Ain[2][2][re] = 1.0;

    //**********
    /*    printf("dm21,dm32   : %f %f \n",dm21,dm32);
    printf("s12,s23,s31 : %f %f %f \n",s12,s23,s31);
    printf("dm  : %f %f %f \n",dm[0][0],dm[0][1],dm[0][2]);
    printf("dm  : %f %f %f \n",dm[1][0],dm[1][1],dm[1][2]);
    printf("dm  : %f %f %f \n",dm[2][0],dm[2][1],dm[2][2]);
    ***********
    **********
    printf("mix : %f %f %f \n",mix[0][0][0],mix[0][1][0],mix[0][2][0]);
    printf("mix : %f %f %f \n",mix[1][0][0],mix[1][1][0],mix[1][2][0]);
    printf("mix : %f %f %f \n",mix[2][0][0],mix[2][1][0],mix[2][2][0]);
    printf("mix : %f %f %f \n",mix[0][0][1],mix[0][1][1],mix[0][2][1]);
    printf("mix : %f %f %f \n",mix[1][0][1],mix[1][1][1],mix[1][2][1]);
    printf("mix : %f %f %f \n",mix[2][0][1],mix[2][1][1],mix[2][2][1]);*/
    //***********
}			

// main kernel
__global__ void get_vacuum_probability(double mix_device[][3][2], int nutype, int beta, double *energy, int n, double path, double *osc_weight, double tdm21, double tdm32)
{
  double lovere ;
  double s21, s32, s31, ss21, ss32, ss31 ;
  int ista, iend ;
  double prob[3][3];
    //double prob2[3][3][2] *****************

  // index
  int idx = (blockIdx.x * blockDim.x + threadIdx.x);
  //  if (idx > n) return;

  // make more precise 20081003 rvw
  lovere= 1.26693281*(path)/(energy[idx]);
  s21 = sin(tdm21*lovere);
  s32 = sin(tdm32*lovere);
  s31 = sin((tdm21+tdm32)*lovere) ;
  ss21 = s21*s21 ;
  ss32 = s32*s32 ;
  ss31 = s31*s31 ;

  /* ista = abs(*nutype) - 1 ; */
  for ( ista=0 ; ista<3 ; ista++ )
  {
    for ( iend=0 ; iend<2 ; iend++ )
    {
      prob[ista][iend]  = mix_device[ista][0][re]*mix_device[iend][0][re]*
        mix_device[ista][1][re]*mix_device[iend][1][re]*ss21;
      prob[ista][iend] += mix_device[ista][1][re]*mix_device[iend][1][re]*
        mix_device[ista][2][re]*mix_device[iend][2][re]*ss32;
      prob[ista][iend] += mix_device[ista][2][re]*mix_device[iend][2][re]*
        mix_device[ista][0][re]*mix_device[iend][0][re]*ss31;
      if ( iend == ista )
      {
        prob[ista][iend]  = 1.0-4.0*prob[ista][iend];
      }
       else
       {
        prob[ista][iend]  = -4.0*prob[ista][iend];
      }
    }
    prob[ista][2]=1.0-prob[ista][0]-prob[ista][1];
  }

  nutype = abs(nutype);
  beta = abs(beta);

  //if ( nutype > 0 )
  double ans = prob[nutype-1][beta-1];
  osc_weight[idx] = ans;

/*  if ( nutype < 0 ) // assuming CPT!!!
    osc_weight[idx] = prob[beta-1][nutype-1];

    osc_weight[idx]= 1.2;*/
}


extern "C" __host__ double* GetVacuumProb( int Alpha, int Beta , double *energy_host, int n, double Path )
{
  // alpha -> 1:e 2:mu 3:tau
  // Energy[GeV]
  // Path[km]
  /// simple referes to the fact that in the 3 flavor analysis
  //  the solar mass term is zero

  // create a pointer to device memory
  double *energy_device;

  // specify size of array
  size_t size = n * sizeof(double);

  // CUDA function to allocate memory of size bytes to the address pointed to by device_array
  cudaMalloc((void **) &energy_device, size);

  // copy the array to be squared to the device
  cudaMemcpy(energy_device, energy_host, size, cudaMemcpyHostToDevice);

  double *osc_weights;
  cudaMalloc((void **) &osc_weights, size);

    // copy the mixing matrix to the device
  size_t mixsize = 3*3*2*sizeof(double);

  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);

  //double mix_device[3][3][2];
  mixArray *mix_device;
  	 //mix[0][0][0] = 1;
  cudaMalloc((void **) &mix_device,mixsize);
  cudaMemcpy(mix_device, mix, mixsize, cudaMemcpyHostToDevice);

  dim3 block_size;
  block_size.x = 1024;

  dim3 grid_size;
  grid_size.x = (n / block_size.x) + 1;

  //int block_size = 256;
  //int num_blocks = 1;//n/block_size;

  get_vacuum_probability<<<grid_size, block_size>>>( mix_device, Alpha, Beta, energy_device, n, Path, osc_weights,  dm21, dm32);

  //cudaThreadSynchronize();

  // copy the results back
  double *osc_weights_host = (double*)malloc(size);
  cudaMemcpy(osc_weights_host, osc_weights, size, cudaMemcpyDeviceToHost);

  cudaFree(energy_device);
  cudaFree(osc_weights);
  cudaFree(mix_device);

  return osc_weights_host;	
}

extern "C" __host__ void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta,/* double Energy_ ,*/ bool kSquared)
{

  double sin12;
  double sin13;
  double sin23;

  if (kSquared)
    {
      sin12 = sqrt(x12);
      sin13 = sqrt(x13);
      sin23 = sqrt(x23);
    }
  else
    {
      sin12 = sqrt(0.5*(1 - sqrt(1 - x12)));
      sin13 = sqrt(0.5*(1 - sqrt(1 - x13)));
      sin23 = sqrt(0.5*(1 - sqrt(1 - x23)));
    }
  // 1,2,0.5,0.1,0.1,0.1
  init_mixing_matrix(m21, m23, sin12, sin23, sin13, Delta);
	
}

//////////////////////////////////////////////////////////////////////////////////
// the following functions are DEVICE functions for the matter effects calculation
//////////////////////////////////////////////////////////////////////////////////

__device__ void clear_complex_matrix(double A[][3][2])
{
  //memset(A,0,sizeof(double)*18); // turn into a cuda fucniton
 // cudaMemset((void **) A,0,sizeof(double)*18);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
	A[i][j][k] = 0;
}

// ************************************
__device__ void copy_complex_matrix(double A[][3][2], double B[][3][2])
{
  //memcpy(B,A,sizeof(double)*18);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
	B[i][j][k] = A[i][j][k];
}

/*
  multiply complex 3x3 matrix and 3 vector
      W = A X V
 */

 __device__ void multiply_complex_matvec(double A[][3][2], double V[][2],double W[][2])
{	
  for(int i=0;i<3;i++)
{
    W[i][re] = A[i][0][re]*V[0][re]-A[i][0][im]*V[0][im]+
      A[i][1][re]*V[1][re]-A[i][1][im]*V[1][im]+
      A[i][2][re]*V[2][re]-A[i][2][im]*V[2][im] ;
    W[i][im] = A[i][0][re]*V[0][im]+A[i][0][im]*V[0][re]+
      A[i][1][re]*V[1][im]+A[i][1][im]*V[1][re]+
      A[i][2][re]*V[2][im]+A[i][2][im]*V[2][re] ;
  }
}

__device__ void multiply_complex_matrix(double A[][3][2], double B[][3][2], double C[][3][2])
{
  int i,j,k;

  for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
	{
	  for (k=0; k<3; k++)
	    {
	      C[i][j][re] += A[i][k][re] * B[k][j][re] - A[i][k][im] * B[k][j][im];
	      C[i][j][im] += A[i][k][im] * B[k][j][re] + A[i][k][re] * B[k][j][im];
	    }
	}
    }
}


// want to output flavor composition of
// pure mass eigenstate, state
__device__ void convert_from_mass_eigenstate( int state, int flavor, double pure[][2], double mix[3][3][2] )
{
  int    i,j; //,k; *********************************
  double mass    [3][2];
  double conj    [3][3][2];
  int    lstate  = state - 1;
  int    factor  = ( flavor > 0 ? -1. : 1. );

  // need the conjugate for neutrinos but not for
  // anti-neutrinos
  for (i=0; i<3; i++) {
    mass[i][0] = ( lstate == i ? 1.0 : 0. );
    mass[i][1] = (                     0. );
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      conj[i][j][re] =        mix[i][j][re];
      conj[i][j][im] = factor*mix[i][j][im];
    }
  }
  multiply_complex_matvec(conj, mass, pure);
}



__device__ void get_transition_matrix(int nutypei,double Enuf,double rhof,double Lenf,double Aout[][3][2],double phase_offsetf, double mix[3][3][2], double dm[3][3])
{
  int nutype; //, make_average ;**************
  double Enu, rho, Len ;
  double dmMatVac[3][3], dmMatMat[3][3];
  double phase_offset;
  nutype=nutypei;
  Enu=Enuf ;
  rho=rhof ;
  Len=Lenf ;
  phase_offset = phase_offsetf ;
  /*   propagate_mat(Ain,rho,Len,Enu,mix,dm,nutype,Aout);    */
  //printf("GPU - %f %f %i \n", Enu, rho, nutype);
  getM(Enu, rho, mix, dm, nutype, dmMatMat, dmMatVac);
  getA(Len, Enu, rho, mix, dmMatVac, dmMatMat, nutype, Aout,phase_offset);

  //for (int i = 0; i < 3; ++i)
    //printf("mix %f %f %f\n", mix[0][i][0], mix[1][i][0], mix[2][i][0]);
  //printf("gpu dm %f %f %f\n", dm[0][0], dm[1][0], dm[2][0]);
  //printf("Aout gpu %f %f %f\n", Aout[0][0][0], Aout[1][0][0], Aout[2][0][0]);

  //  Aout[0][0][0] =  dm[0][0];
}


// the colonel! (kernel...)
__global__ void propagateLinear(int Alpha, int Beta, double Path, double Density, /*double Mix[3][3][2], double dm[3][3],*/ double *Energy, double *osc_w, int n)
{
  // here we go
  bool kUseMassEigenstates = false; // quick hack for now
  int idx = (blockIdx.x * blockDim.x + threadIdx.x);
  if (idx < n)
    {
  double Probability[3][3];

  int i,j;

  double TransitionMatrix[3][3][2];
  //  double TransitionProduct[3][3][2];
  //  double TransitionTemp[3][3][2];
  double RawInputPsi[3][2];
  double OutputPsi[3][2];

  get_transition_matrix( Alpha,
			 Energy[idx],               // in GeV
			 Density * 0.5,
			 Path,          // in km
			 TransitionMatrix,     // Output transition matrix
			 0.0,
			 mix_device,
			 dm_device);

    //copy_complex_matrix( TransitionMatrix , TransitionProduct );

  for ( i = 0 ; i < 3 ; i++ ){
    for ( j = 0 ; j < 3 ; j++ ){
      Probability[i][j]=0;
    }
  }

  for ( i = 0 ; i < 3 ; i++ )
    {
    for ( j = 0 ; j < 3 ; j++ )
    	{       RawInputPsi[j][0] = 0.0; RawInputPsi[j][1] = 0.0;   }

      if( kUseMassEigenstates )
	convert_from_mass_eigenstate( i+1, Alpha,  RawInputPsi, mix_device );
      else
	RawInputPsi[i][0] = 1.0;

      multiply_complex_matvec( TransitionMatrix /*Product*/, RawInputPsi, OutputPsi );



      Probability[i][0] += OutputPsi[0][0] * OutputPsi[0][0] + OutputPsi[0][1]*OutputPsi[0][1];
      Probability[i][1] += OutputPsi[1][0] * OutputPsi[1][0] + OutputPsi[1][1]*OutputPsi[1][1];
      Probability[i][2] += OutputPsi[2][0] * OutputPsi[2][0] + OutputPsi[2][1]*OutputPsi[2][1];

    }

  // now do the part that getprob usually does
  int In = abs( Alpha );
  int Out = abs( Beta );
  osc_w[idx] = Probability[In-1][Out-1];
  
    }
}

// this kernel calculates the weights for different combinations of P(alpha -> beta)
__global__ void propagateLinearAll(int *Alpha, int *Beta, double Path, double Density, /*double Mix[3][3][2], double dm[3][3],*/ double *Energy, double *osc_w, int n)
{
  bool kUseMassEigenstates = false; // quick hack for now
  int idx = (blockIdx.x * blockDim.x + threadIdx.x);
  if (idx < n)
    {
      double Probability[3][3];
      int i,j;
      double TransitionMatrix[3][3][2];
      double RawInputPsi[3][2];
      double OutputPsi[3][2];
      get_transition_matrix( Alpha[idx],
			     Energy[idx],               // in GeV
			     Density * 0.5,
			     Path,          // in km
			     TransitionMatrix,     // Output transition matrix
			     0.0,
			     mix_device,
			     dm_device);

      for ( i = 0 ; i < 3 ; i++ )
	{
	  for ( j = 0 ; j < 3 ; j++ )
	    {       RawInputPsi[j][0] = 0.0; RawInputPsi[j][1] = 0.0;   }
	
	  if( kUseMassEigenstates )
	    convert_from_mass_eigenstate( i+1, Alpha[idx],  RawInputPsi, mix_device );
	  else
	    RawInputPsi[i][0] = 1.0;
	  multiply_complex_matvec( TransitionMatrix /*Product*/, RawInputPsi, OutputPsi );
	
	  Probability[i][0] += OutputPsi[0][0] * OutputPsi[0][0] + OutputPsi[0][1]*OutputPsi[0][1];
	  Probability[i][1] += OutputPsi[1][0] * OutputPsi[1][0] + OutputPsi[1][1]*OutputPsi[1][1];
	  Probability[i][2] += OutputPsi[2][0] * OutputPsi[2][0] + OutputPsi[2][1]*OutputPsi[2][1];
	}
      int In = abs( Alpha[idx] );
      int Out = abs( Beta[idx] );
      osc_w[idx] = Probability[In-1][Out-1];
    }
}

extern "C" __host__ void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw)
{
  // copy DM matrix
  size_t dmsize = 3*3*sizeof(double);
  typedef double dmArray[3];
  dmArray *d = (dmArray*)malloc(dmsize);
  memcpy(d, &dm, dmsize);

  cudaMemcpyToSymbol(dm_device, dm, dmsize);
  //dmArray *dm_device;
  //cudaMalloc((void **) &dm_device, dmsize);
  //  cudaMemcpy(dm_device, dm, dmsize, cudaMemcpyHostToDevice);

  // copy mns matrix to device
  size_t mixsize = 3*3*2*sizeof(double);
  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);

  cudaMemcpyToSymbol(mix_device, mix, mixsize);
  //mixArray *mix_device;
  //  cudaMalloc((void **) &mix_device,mixsize);
  //cudaMemcpy(mix_device, mix, mixsize, cudaMemcpyHostToDevice);

  // copy energy array to device
  size_t size = n * sizeof(double);
  double *energy_device = NULL;

  cudaMalloc((void **) &energy_device, size);
  cudaMemcpy(energy_device, Energy, size, cudaMemcpyHostToDevice);

  // allocate output memory space on the device
  double *osc_weights;
  cudaMalloc((void **) &osc_weights, size);

  dim3 block_size;
  block_size.x = 128;

  dim3 grid_size;
  grid_size.x = (n / block_size.x) + 1;

  propagateLinear<<<grid_size, block_size>>>(Alpha, Beta, Path, Density, /*mix_device, dm_device,*/ energy_device, osc_weights, n);
  //CudaCheckError();

  // copy the results back
  cudaMemcpy(oscw, osc_weights, size, cudaMemcpyDeviceToHost);

  cudaFree(energy_device);
  cudaFree(osc_weights);
  //  cudaFree(mix_device);
  //cudaFree(dm_device);
  free(m);
  free(d);
}

/*extern "C"*/ __host__ void GetProbAll(int *d_alpha, int *d_beta, double Path, double Density, double *energy_device, int n, double *osc_weights, double *oscw)
{
  // copy the dm matrix
  /*   size_t dmsize = 3*3*sizeof(double);
  typedef double dmArray[3];
  dmArray *d = (dmArray*)malloc(dmsize);
  memcpy(d, &dm, dmsize);
  dmArray *dm_device;
  cudaMalloc((void **) &dm_device, dmsize);
  cudaMemcpy(dm_device, dm, dmsize, cudaMemcpyHostToDevice);

  // copy mns matrix to device
  size_t mixsize = 3*3*2*sizeof(double);
  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);
  mixArray *mix_device;
  cudaMalloc((void **) &mix_device,mixsize);
  cudaMemcpy(mix_device, mix, mixsize, cudaMemcpyHostToDevice);*/
  //CudaCheckError();

  // copy DM matrix
  size_t dmsize = 3*3*sizeof(double);
  typedef double dmArray[3];
  dmArray *d = (dmArray*)malloc(dmsize);
  memcpy(d, &dm, dmsize);
  cudaMemcpyToSymbol(dm_device, dm, dmsize);

  // copy mns matrix to device
  size_t mixsize = 3*3*2*sizeof(double);
  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);

  cudaMemcpyToSymbol(mix_device, mix, mixsize);


  // from here //
  // copy energy array to device
  //size_t size = n * sizeof(double);
  //  double *energy_device = NULL;

  //  cudaMalloc((void **) &energy_device, size);
  //  cudaMemcpy(energy_device, Energy, size, cudaMemcpyHostToDevice);

  // allocate output memory space on the device
  /*double *osc_weights = NULL;
    cudaMalloc((void **) &osc_weights, size);*/

  // copy alpha and beta
  /*int *d_alpha = NULL;
  int *d_beta = NULL;
  cudaMalloc((void **) &d_alpha, size);
  cudaMalloc((void **) &d_beta, size);
  cudaMemcpy(d_alpha, Alpha, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_beta, Beta, size, cudaMemcpyHostToDevice);*/
  // to here //

  //dim3 block_size;
  //block_size.x
    int blocksize = 128;//512;

  //dim3 grid_size;
  int gridsize = (n / blocksize) + 1;
  //grid_size.x = (n / block_size.x) + 1;

  //printf("gridsize %i total things %i", gridsize, n);
  //printf("A WEIGHT : %f : mix %f %i %i", oscw[100], mix[0][0][0], gridsize, blocksize);

  propagateLinearAll<<<gridsize, blocksize>>>(d_alpha, d_beta, Path, Density, /*mix_device, dm_device,*/ energy_device, osc_weights,n);
  //CudaCheckError();

  // copy the results back
  cudaMemcpy(oscw, osc_weights, n * sizeof(double), cudaMemcpyDeviceToHost);
  //CudaCheckError();


  //cudaFree(mix_device);
  //cudaFree(dm_device);

  // from here //

  // to here //

  free(m);
  free(d);
}

/*extern "C"*/ __host__ void probInitGPU(double **energy_device, double **osc_weights, int **d_alpha, int **d_beta, int n)
{
  cudaMalloc((void **) energy_device, n * sizeof(double));
  CudaCheckError();

  cudaMalloc((void **) osc_weights, n * sizeof(double));
  CudaCheckError();

  cudaMalloc((void **) d_alpha, n * sizeof(int));
  CudaCheckError();

  cudaMalloc((void **) d_beta, n * sizeof(int));
  CudaCheckError();
}

/*extern "C"*/ __host__ void probCopyToGPU(double *gpu_energy, double *energy, int*gpu_alpha, int *alpha, int *gpu_beta, int *beta, int n)
{
  //cudaThreadSynchronize();
  printf("- copying %i events", n);
  cudaMemcpy(gpu_energy, energy, n*sizeof(double), cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_alpha, alpha, n*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaMemcpy(gpu_beta, beta, n*sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();

  cudaDeviceSynchronize();
}

/*extern "C"*/ __host__ void probCleanupGPU(double *gpu_energy, double *gpu_osc_weights, int *gpu_alpha, int *gpu_beta)
{
  cudaFree(gpu_energy);
  cudaFree(gpu_osc_weights);
  cudaFree(gpu_alpha);
  cudaFree(gpu_beta);
}

////////////////////////////////
// the functions below are for earth matter fx
////////////////////////////////

// CONSTANT MEMORY FOR EARTH PROFILE ///////////////////////
#define maxLayers 5 // this is a temporary measure
static bool earth_profile_set = false; // do we need to set the profile?
__constant__ double device_Rhos[maxLayers];
__constant__ double device_Radii[maxLayers];
__constant__ double device_CosLimit[maxLayers];
__constant__ double device_YpMap[maxLayers];
//__constant__ double device_REarth;
////////////////////////////////////////////////////////////

// this function must be called once at the beginning to load the earth profile into constant memory
__host__ void LoadEarthProfile()
{
  // hard coded density
  // this is from official SK atm result
  int MaxDepth = maxLayers;

  double _Rhos[MaxDepth];
  double _YpMap[MaxDepth];
  double _Radii[MaxDepth];
  double _CosLimit[MaxDepth];

  // they are reversed
  _Radii[0] = 6371.0; _Rhos[0] = 3.3;  _YpMap[0] = 0.497;
  _Radii[1] = 5701.0; _Rhos[1] = 5.0;  _YpMap[1] = 0.497;
  _Radii[2] = 3480.0; _Rhos[2] = 11.3; _YpMap[2] = 0.497;
  _Radii[3] = 1220.0; _Rhos[3] = 13.0; _YpMap[3] = 0.468;
  _Radii[4] = 0.0;    _Rhos[4] = 13.0; _YpMap[4] = 0.468;

  double REarth = _Radii[0];

  /*_density[ 0 ]       =  13.0 ;
  _density[ 1220.0 ]  =  13.0 ;
  _density[ 3480.0 ]  =  11.3 ;
  _density[ 5701.0 ]  =  5.0 ;
  _density[ 6371.0 ]  =  3.3 ;*/

  ///
  double x;
  //_CosLimit.clear();

  // first element of _Radii is largest radius!
  for (int i = 0; i < /*(int) _Radii.size()*/ MaxDepth; i++)
    {
      // Using a cosine threshold instead! //
      x = -1 * sqrt( 1 - (_Radii[i] * _Radii[i] / ( REarth*REarth)) );
      if ( i  == 0 ) x = 0;
      _CosLimit[i /*_Radii[i]*/] = x;
    }
  ///

  // copy to constant memory
  CudaSafeCall( cudaMemcpyToSymbol(device_Rhos, _Rhos, maxLayers * sizeof(double)) );
  CudaSafeCall( cudaMemcpyToSymbol(device_YpMap, _YpMap, maxLayers * sizeof(double)) );
  CudaSafeCall( cudaMemcpyToSymbol(device_Radii, _Radii, maxLayers * sizeof(double)) );
  CudaSafeCall( cudaMemcpyToSymbol(device_CosLimit, _CosLimit, maxLayers * sizeof(double)) );
  //CudaSafeCall( cudaMemcpyToSymbol(device_REarth, REarth, sizeof(double)) );
}

__host__ void SetChemicalComposition(int nChemCompVals, double* ChemCompVals) {
  cudaFree(device_YpMap);

  if (nChemCompVals != maxLayers) {
    printf("Incompatible number of layers!\n");
    printf("Given:%i \n",nChemCompVals);
    printf("Expecting:%i \n",maxLayers);
    exit(-1);
  }

  int MaxDepth = maxLayers;
  double _YpMap[MaxDepth];
  _YpMap[0] = ChemCompVals[0];
  _YpMap[1] = ChemCompVals[1];
  _YpMap[2] = ChemCompVals[2];
  _YpMap[3] = ChemCompVals[3];
  _YpMap[4] = ChemCompVals[4];

  CudaSafeCall( cudaMemcpyToSymbol(device_YpMap, _YpMap, maxLayers * sizeof(double)));
}

__host__ void ResetChemicalComposition() {
  int nChemCompVals = 5;
  double ChemCompVals[nChemCompVals];
  
  ChemCompVals[0] = 0.497;
  ChemCompVals[1] = 0.497;
  ChemCompVals[2] = 0.497;
  ChemCompVals[3] = 0.468;
  ChemCompVals[4] = 0.468;
  
  SetChemicalComposition(nChemCompVals,ChemCompVals);
}

__host__ void SetChemicalComposition_PremModel(double Yp_Val) {
  int nChemCompVals = 5;
  double ChemCompVals[nChemCompVals];

  ChemCompVals[0] = 0.497;
  ChemCompVals[1] = 0.497;
  ChemCompVals[2] = 0.497;
  ChemCompVals[3] = Yp_Val;
  ChemCompVals[4] = Yp_Val;

  SetChemicalComposition(nChemCompVals,ChemCompVals);
}

__device__ void SetDensityProfile(double CosineZ, double PathLength , double ProductionHeight, int &Layers, double *_TraverseRhos, double *_TraverseDistance, double *_TraverseYpMap)
{
  double REarth = device_Radii[0];// * 1.0e5;
  int i;
  int MaxLayer;
  double km2cm = 1.0E5;
  double TotalEarthLength =  -2.0 * CosineZ * REarth * km2cm; // in [cm]
  double CrossThis, CrossNext;

  //printf("gpu total earth length in cm %f \n", TotalEarthLength);


  // path through air
  _TraverseRhos[0] = 0.0;
  _TraverseYpMap[0] = 0.0;
  _TraverseDistance[0] =  (PathLength) - TotalEarthLength ;
  //printf("setting first element %f - %f = %f \n", PathLength, TotalEarthLength, _TraverseDistance[0] );

  if ( CosineZ >= 0 )
    {
      _TraverseDistance[0] =  PathLength;
      Layers = 1;
      return;
    }

  Layers = 0;

  for (int i = 0; i < maxLayers; ++i)
    {

      if (CosineZ < device_CosLimit[i])
	{
	  Layers++;
	}
    }

  MaxLayer = Layers;

  // the zeroth layer is the air!
  //#pragma unroll
  for ( i = 0 ; i< MaxLayer ; i++ )
    {
      _TraverseRhos[i+1] = device_Rhos[i];
      _TraverseYpMap[i+1] = device_YpMap[i];

      CrossThis = 2.0 * sqrt( device_Radii[i] * device_Radii[i]         - REarth*REarth*( 1 -CosineZ*CosineZ ) );
      CrossNext = 2.0 * sqrt( device_Radii[i+1] * device_Radii[i+1]      - REarth*REarth*( 1 -CosineZ*CosineZ ) );

      if( i < MaxLayer-1 )
	_TraverseDistance[i+1]  =  0.5*( CrossThis-CrossNext ) * km2cm;
      else
	_TraverseDistance[i+1]  =  CrossThis * km2cm;

      //  if (i == 0)
      //	printf("gpu crossthis %f crossnext %f \n", CrossThis, CrossNext);

      // assumes azimuthal symmetry
      if( i < MaxLayer )
	{
	  _TraverseRhos    [ 2*MaxLayer - i ] = _TraverseRhos[i];
	  _TraverseYpMap    [ 2*MaxLayer - i ] = _TraverseYpMap[i];
	  _TraverseDistance[ 2*MaxLayer - i ] = _TraverseDistance[i];
	  //  printf("honkbeep gpu %i %f\n", 2*MaxLayer - i, _TraverseDistance[i]);
	}

      //printf("GPU:: %i rhos %f dists %f\n", i, _TraverseRhos[i], _TraverseDistance[i]);
    }
  Layers = 2*MaxLayer;
}


__global__ void propagate(int in_flav, int out_flav, double *enu_device, double *cosz_device, double prod_h, const int n, double *out_device)
{
  bool kUseMassEigenstates = false; // quick hack for now
  int idx = (blockIdx.x * blockDim.x + threadIdx.x);

  if (idx < n)
    {
      int alpha = in_flav;//[idx];

      double _TraverseDistance[2*maxLayers+1];
      double _TraverseRhos[2*maxLayers+1];
      double _TraverseYpMap[2*maxLayers+1];

      for (int i = 0; i < 2*maxLayers+1; ++i)
	{
	  _TraverseDistance[i] = 0;
	  _TraverseRhos[i] = 0;
	  _TraverseYpMap[i] = 0.;
	}

      int Layers;

      // first we need to fill the above 3 variables
      // using the earth profile

      prod_h *= 1e5;
      double rearth = device_Radii[0] * 1.0e5;
      double costh = cosz_device[idx];
      double PathLength = sqrt( (rearth + prod_h ) * (rearth + prod_h)
				- (rearth*rearth)*( 1 - costh * costh)) - rearth * costh;
      //printf("GPU pathlength: %f : ph %f rearth %f cz %f \n", PathLength, prod_h, rearth, costh);

      SetDensityProfile( costh, PathLength, prod_h, Layers, _TraverseRhos, _TraverseDistance, _TraverseYpMap);

      // now we can proceed as usual
      // with the addition of multiplying the transition matrices for each layer
      double Probability[3][3];
      int i,j;

      double TransitionMatrix[3][3][2];
      double TransitionProduct[3][3][2];
      double TransitionTemp[3][3][2];
      double RawInputPsi[3][2];
      double OutputPsi[3][2];

      // loop over the layers traversed
      //#pragma unroll
      for ( i = 0; i < Layers ; i++ )
	{
	  //printf("Layer %i, dist: %f, density: %f\n", i, _TraverseDistance[i], _TraverseRhos[i]);
	  get_transition_matrix( alpha,
				 enu_device[idx],    // in GeV
				 _TraverseRhos[i] * _TraverseYpMap[i], //density_convert,
				 _TraverseDistance[i] / 1.0e5,   // in km
				 //Density * 0.5,
				 //Path,          // in km
				 TransitionMatrix,     // Output transition matrix
				 0.0,
				 mix_device,
				 dm_device);

	  // MAKE SURE THESE FUNCTIONS EXIST AND WORK
	  if ( i == 0 )
	    copy_complex_matrix( TransitionMatrix , TransitionProduct );

	  if ( i >0 )
	    {
	      clear_complex_matrix( TransitionTemp );
	      multiply_complex_matrix( TransitionMatrix, TransitionProduct, TransitionTemp );
	      copy_complex_matrix( TransitionTemp, TransitionProduct );
	    }
	  //printf("%f %f\n", TransitionMatrix[0][0][0], TransitionTemp[0][0][0]);
	}

      //copy_complex_matrix( TransitionMatrix , TransitionProduct );

      for ( i = 0 ; i < 3 ; i++ )
	{
	  for ( j = 0 ; j < 3 ; j++ )
	    {       RawInputPsi[j][0] = 0.0; RawInputPsi[j][1] = 0.0;   }
	
	  if( kUseMassEigenstates )
	    convert_from_mass_eigenstate( i+1, in_flav,  RawInputPsi, mix_device );
	  else
	    RawInputPsi[i][0] = 1.0;
	
	  multiply_complex_matvec( TransitionProduct, RawInputPsi, OutputPsi );
	
	  Probability[i][0] += OutputPsi[0][0] * OutputPsi[0][0] + OutputPsi[0][1]*OutputPsi[0][1];
	  Probability[i][1] += OutputPsi[1][0] * OutputPsi[1][0] + OutputPsi[1][1]*OutputPsi[1][1];
	  Probability[i][2] += OutputPsi[2][0] * OutputPsi[2][0] + OutputPsi[2][1]*OutputPsi[2][1];
	}

      // now do the part that getprob usually does
      int In = abs( alpha );
      int Out = abs( out_flav );

      //printf("%i %i %f\n", In-1, Out-1, Probability[In-1][Out-1]);
      out_device[idx] = Probability[In-1][Out-1];
      //printf("%i %i %f: layers %i\n", In, Out, Probability[In-1][Out-1], Layers);
    }
}

// DB Function changed to take single out_flav instead of array
extern "C" __host__ void GetProbAtm(int in_flav, int out_flav, double *enu, double *cosz, double prod_h, double Yp_Val, const int n, double *out)
{
  if (!earth_profile_set)
    {
      LoadEarthProfile();
      earth_profile_set = true;
    }
  SetChemicalComposition_PremModel(Yp_Val);

  size_t dmsize = 3*3*sizeof(double);
  typedef double dmArray[3];
  dmArray *d = (dmArray*)malloc(dmsize);
  memcpy(d, &dm, dmsize);

  cudaMemcpyToSymbol(dm_device, dm, dmsize);

  // copy mns matrix to device
  size_t mixsize = 3*3*2*sizeof(double);
  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);
  cudaMemcpyToSymbol(mix_device, mix, mixsize);
  // copy energy array to device
  size_t size = n * sizeof(double);
  double *energy_device = NULL;

  cudaMalloc((void **) &energy_device, size);
  cudaMemcpy(energy_device, enu, size, cudaMemcpyHostToDevice);

  // copy the cosz to device
  double *cosz_device = NULL;

  cudaMalloc((void **) &cosz_device, size);
  cudaMemcpy(cosz_device, cosz, size, cudaMemcpyHostToDevice);

  // DB this function was previously expecting an array of final flavours where I only give it a single flavour. Hence don't need to copy array
  /*
  // copy the output flavours
  int *out_flav_device;
  cudaMalloc((void **) &out_flav_device, size);
  printf("%i %i %i %i",out_flav_device, out_flav, size, cudaMemcpyHostToDevice);
  cudaMemcpy(out_flav_device, out_flav, size, cudaMemcpyHostToDevice);
  */

  // make an output array
  double *out_device;
  cudaMalloc((void **) &out_device, size);
  CudaCheckError();

  // now we have copied the mixing and mass matrices, we need to
  // launch the kernel
  int blocksize = 128;//512;

  int gridsize = (n / blocksize) + 1;
  //propagate<<<gridsize, blocksize>>>(in_flav, out_flav_device, energy_device, cosz_device, prod_h, n, out_device);
  propagate<<<gridsize, blocksize>>>(in_flav, out_flav, energy_device, cosz_device, prod_h, n, out_device);
  CudaCheckError();

  // copy the results back
  cudaMemcpy(out, out_device, n * sizeof(double), cudaMemcpyDeviceToHost);
  CudaCheckError();

  cudaFree(energy_device);
  cudaFree(cosz_device);
  cudaFree(out_device);

  free(m);
  free(d);
}

__host__ void probInitGPU_atm(double *energy_host, double *cosz_host, int *out_flav_host, const int N, double **energy_device, double **cosz_device, int **out_flav_device, double **out_device)
{
  // allocate memory on gpu
  cudaMalloc((void **) energy_device, N * sizeof(double));
  cudaMalloc((void **) cosz_device, N * sizeof(double));
  cudaMalloc((void **) out_flav_device, N * sizeof(int));
  CudaCheckError();

  cudaMemcpy(*energy_device, energy_host, N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(*cosz_device, cosz_host, N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(*out_flav_device, out_flav_host, N * sizeof(int), cudaMemcpyHostToDevice);
  CudaCheckError();

  // make an output array

  cudaMalloc((void **) out_device, N * sizeof(double));
}

__host__ void probCleanupGPU_atm( double *energy_device, double *cosz_device, int *out_flav_device, double *out_device)
{
  cudaFree(energy_device);
  cudaFree(cosz_device);
  cudaFree(out_flav_device);
  cudaFree(out_device);
}

__host__ void GetProb(int in_flav, int out_flav_device, double *energy_device, double *cosz_device, double prod_h, double Yp_Val, const int n, double *out, double *out_device)
{
  if (!earth_profile_set)
    {
      LoadEarthProfile();
      earth_profile_set = true;
    }
  SetChemicalComposition_PremModel(Yp_Val);

  size_t dmsize = 3*3*sizeof(double);
  typedef double dmArray[3];
  dmArray *d = (dmArray*)malloc(dmsize);
  memcpy(d, &dm, dmsize);

  cudaMemcpyToSymbolAsync(dm_device, dm, dmsize);

  // copy mns matrix to device

  size_t mixsize = 3*3*2*sizeof(double);
  typedef double mixArray[3][2];
  mixArray *m = (mixArray*)malloc(mixsize);
  memcpy(m, &mix, mixsize);
  cudaMemcpyToSymbolAsync(mix_device, mix, mixsize);

  // launch the kernel

  int blocksize = 64;//64;//256;//256;//512;

  int gridsize = (n / blocksize) + 1;

  propagate<<<gridsize, blocksize>>>(in_flav, out_flav_device, energy_device, cosz_device, prod_h, n, out_device);


  cudaMemcpyAsync(out, out_device, n * sizeof(double), cudaMemcpyDeviceToHost);
  CudaCheckError();

  free(m);
  free(d);
}

