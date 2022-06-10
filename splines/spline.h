#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

class spline
{
 public:
  spline(int N, int tfDelta, int tfXmin, int tfXmax, double tfX[], double tfY[], double tfB[], double tfC[], double tfD[] )
    {
      
      fNp = N; 
      fKstep = 0; // sort this out
      fDelta = tfDelta;
      fXmin = tfXmin;
      fXmax = tfXmax;
      fX = tfX;
      fY = tfY;
      fB = tfB;
      fC = tfC;
      fD = tfD;
    };
  
  ~spline()
    {
    };
  
  CUDA_CALLABLE_MEMBER double Eval(double x)
  {
    int klow=0;
    if(x<=fXmin) klow=0;
    else if(x>=fXmax) klow=fNp-1;
    else
      {
	if(fKstep) 
	  {
	    // Equidistant knots, use histogramming
	    klow = int((x-fXmin)/fDelta);
	    if (klow < fNp-1) klow = fNp-1;
	  } 
	else
	  {
	    int khig=fNp/*-1*/, khalf;
	    // Non equidistant knots, binary search
	    while(khig-klow>1)
	      if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
	      else khig=khalf;
	  }
      }
    if(klow >= fNp-1) klow = fNp -2;
    
    // Evaluate now
    double dx=x-fX[klow];
    return ( fY[klow] + dx * ( fB[klow] + dx * ( fC[klow] + dx * fD[klow] ) ) );
  };
  
 protected:
  int fNp;
  int fKstep;
  
  double fDelta;
  double fXmin;
  double fXmax;
  
  double *fX;
  double *fY;
  double *fB;
  double *fC;
  double *fD;
};

