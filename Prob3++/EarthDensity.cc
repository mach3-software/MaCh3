#include "EarthDensity.h"
#include <iostream>
#include <sstream>

EarthDensity::EarthDensity( ) 
{
   //cout << "EarthDensity::EarthDensity Using Default density profile  " << endl;

   // radius: [ km ]  density  [ g/cm^3 ]
   _density[ 0 ]       =  13.0 ;
   _density[ 1220.0 ]  =  13.0 ;
   _density[ 3480.0 ]  =  11.3 ;
   _density[ 5701.0 ]  =  5.0 ;
   _density[ 6371.0 ]  =  3.3 ;

   // chemical composition
   _YpMap[ 0 ]       =  0.468 ;
   _YpMap[ 1220.0 ]  =  0.468 ;
   _YpMap[ 3480.0 ]  =  0.497 ;
   _YpMap[ 5701.0 ]  =  0.497 ;
   _YpMap[ 6371.0 ]  =  0.497 ;

   _TraverseDistance  = NULL;
   _TraverseRhos      = NULL;
   _Yp                = NULL;
   init();

       
}


EarthDensity::EarthDensity( const char * file )
{
    _TraverseDistance  = NULL;
    _TraverseRhos      = NULL;
    _Yp                = NULL;
    LoadDensityProfile( file );
}

void EarthDensity::LoadDensityProfile( const char * file )
{
        ifstream PREM_dat;
        double r_dist;          // radial distance -- map key //
        double rho;             // density at that distance -- map value //
	double lyp;

	DensityFileName = file;

        PREM_dat.open(DensityFileName.c_str());
        if(! PREM_dat)
        {
                cerr<<"EarthDensity::Load ERROR OPENING " << DensityFileName << endl;
                exit(1);
        }
        else
           cout << "Loading Density profile from: " << DensityFileName << endl;

	// count the number of columns to separate IO between reading 
        // files like 
        //   radius rho     and  
        //   radius rho  Yp and  
	std::string line;
        getline( PREM_dat , line);         
	std::stringstream s;
        s << line;                   
        int count = 0;
        double dummy; 
        while ( s >> dummy ) count++;

        // rewind the file 
        PREM_dat.clear();
        PREM_dat.seekg(0);

        while( !PREM_dat.eof( ) )
        {
                if ( r_dist > REarth ) REarth = r_dist;

                if( count == 2 ) 
		  { 
		    PREM_dat >> r_dist >> rho ;
		    _density[r_dist] = rho ;
		    _YpMap  [r_dist] = 0.5 ;
		  }
                else
		  { 
		    PREM_dat >> r_dist >> rho >> lyp ;
		    _density[r_dist] = rho ;
		    _YpMap  [r_dist] = lyp ;
		  }
        }
        PREM_dat.close();

        // must be re-initialized after a profile is loaded
        init();        
}


void EarthDensity::init()
{

	Load();
	ComputeMinLengthToLayers();
}


///// Really need to clean this bit up, slow and bulky!
void EarthDensity::SetDensityProfile( double CosineZ, double PathLength , double ProductionHeight )
{
   int i;
   int MaxLayer;
   double km2cm = 1.0e5;
   double TotalEarthLength =  -2.0*CosineZ*REarth*km2cm; // in [cm]
   double CrossThis, CrossNext;

   map<double, double>::iterator _i;


   // path through air
   _TraverseRhos[0] = 0.0 ;
   _Yp[0]           = 0.0 ;
   _TraverseDistance[0] =  PathLength - TotalEarthLength ;

   if( CosineZ >= 0 )
   {  
       _TraverseDistance[0] =  PathLength;
       Layers = 1; 
       return; 
   }

   	
   Layers = 0;
   for ( _i = _CosLimit.begin(); _i != _CosLimit.end() ; _i++ )
      if( CosineZ < _i->second )
          Layers++;	


   MaxLayer = Layers;

   // the zeroth layer is the air!
   for ( i = 0 ; i< MaxLayer ; i++ ) 
   {
 
      _TraverseRhos[i+1]      = _Rhos[i];
      _Yp          [i+1]      = _Yps [i];
      
      CrossThis = 2.0*sqrt( _Radii[i]*_Radii[i]  - REarth*REarth*( 1 -CosineZ*CosineZ ) );
      CrossNext = 2.0*sqrt( _Radii[i+1]*_Radii[i+1]      - REarth*REarth*( 1 -CosineZ*CosineZ ) );

     if( i < MaxLayer-1 )
      _TraverseDistance[i+1]  =  0.5*( CrossThis-CrossNext )*km2cm;
     else
      _TraverseDistance[i+1]  =  CrossThis*km2cm;
    
     // assumes azimuthal symmetry    
     if( i < MaxLayer ){
        _TraverseRhos    [ 2*MaxLayer - i ] = _TraverseRhos[i];
        _TraverseDistance[ 2*MaxLayer - i ] = _TraverseDistance[i];
	_Yp              [ 2*MaxLayer - i ] = _Yp[i]; 
     }
   }

   Layers = 2*MaxLayer; 

}


// now using Zenith angle to compute minimum conditions...20050620 rvw
void EarthDensity::ComputeMinLengthToLayers()
{
	double x;

	_CosLimit.clear();

	// first element of _Radii is largest radius!
	for(int i=0; i < (int) _Radii.size() ; i++ )
	{
		// Using a cosine threshold instead! //
		x = -1* sqrt( 1 - (_Radii[i] * _Radii[i] / ( REarth*REarth)) );
		if ( i  == 0 ) x = 0;
		_CosLimit[ _Radii[i] ] = x;
	}

}




void EarthDensity::Load( )
{

        int MaxDepth = 0;

	//map<double, double>::iterator _i;
	map<double, double>::reverse_iterator _i;


      if( _TraverseRhos     != NULL ) delete [] _TraverseRhos;
      if( _TraverseDistance != NULL ) delete [] _TraverseDistance;
      if( _Yp               != NULL ) delete [] _Yp;

      //_i = _density.end();
      //_i--;
      _i = _density.rbegin();
      REarth = _i->first;  // Earth is 6371.0 [km]


       
	// to get the densities in order of decreasing radii
       for( _i = _density.rbegin() ; _i != _density.rend() ; ++_i )
       {
	  _Rhos.push_back( _i->second );
   	  _Radii.push_back( _i->first  );
	  MaxDepth++;
       } 

       for( _i = _YpMap.rbegin() ; _i != _YpMap.rend() ; ++_i )
	 _Yps  .push_back( _i->second  );


      _TraverseRhos      = new double [ 2*MaxDepth + 1 ];
      _TraverseDistance  = new double [ 2*MaxDepth + 1 ];
      _Yp                = new double [ 2*MaxDepth + 1 ];

      return;
                                                                                                                                                             
}

void EarthDensity::SetChemicalComposition(int nChemCompVals, double* ChemCompVals) {

  if ((unsigned)nChemCompVals != _Yps.size()) {
    printf("Incompatible number of layers!\n");
    printf("Given:%i \n",nChemCompVals);
    printf("Expecting:%i \n",(int)_Yps.size());
    exit(-1);
  }

  _YpMap[ 0 ]       =  ChemCompVals[0];
  _YpMap[ 1220.0 ]  =  ChemCompVals[1];
  _YpMap[ 3480.0 ]  =  ChemCompVals[2];
  _YpMap[ 5701.0 ]  =  ChemCompVals[3];
  _YpMap[ 6371.0 ]  =  ChemCompVals[4];

  _Yps[ 0 ]       =  ChemCompVals[0];
  _Yps[ 1 ]       =  ChemCompVals[1];
  _Yps[ 2 ]       =  ChemCompVals[2];
  _Yps[ 3 ]       =  ChemCompVals[3];
  _Yps[ 4 ]       =  ChemCompVals[4];
}

void EarthDensity::SetChemicalComposition_PremModel(double Yp_Val) {
  int nChemCompVals = 5;
  double ChemCompVals[nChemCompVals];

  ChemCompVals[0] = 0.497;
  ChemCompVals[1] = 0.497;
  ChemCompVals[2] = 0.497;
  ChemCompVals[3] = Yp_Val;
  ChemCompVals[4] = Yp_Val;
  
  SetChemicalComposition(nChemCompVals,ChemCompVals);
}

EarthDensity::~EarthDensity( )
{
   if( _TraverseRhos     != NULL ) delete [] _TraverseRhos;
   if( _TraverseDistance != NULL ) delete [] _TraverseDistance;
   if( _Yp               != NULL ) delete [] _Yp; 
}


