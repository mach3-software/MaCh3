#ifndef _EarthDensity_
#define _EarthDensity_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iostream>

//
//  EarthDensity is an object designed to represent the density profile 
//  of the earth, or other planet. It can read in a radial density profile
//  via a user-specified text file. Once the density
//  profile is loaded for given neutrino trajectory, the distance
//  across each layer it will traverse is computed. These are later
//  available for code wishing to propagate the neutrino. In essense
//  a density profile specific to a neutrino's path is created
//
// 20050816 rvw 
// 20081001 rvw (update)

//  User-specified density profiles must contain two columns of floating point
//  numbers, the first is the radial distance [km] from the sphere center, the second
//  is the density [g/cm^3] for 
//  0.      x_0
//  r_1     x_1
//  ..      ..
//  r_n     x_n
//  the last entry should contain the radius of the sphere.
//  each x_i represents the density up to and including r_i
//  the entry for zero radial density must be included. 


using namespace std;

class EarthDensity
{
	public:
                 // default contstructor for the Earth, and a radial density profile 
                 // as specified by the SK 3f paper: PRD.74.032002 (2006) 
                 EarthDensity( );
                 
                 // constructor for a user-specified density profile, see PREM.dat 
		 EarthDensity( const char * );
		~EarthDensity( );

                void init();
	
		// Load the Density profile for a given neutrino path:
		//
		// 	Cosine Zenith,
	        //	Path Length - 
		//      Production Height - amount of vacuum to travese before matter 
		//		corresponds to height in atmosphere for atmospheric neutrinos
		//		or distance from neutrion source for extraplanetary nu's
		void SetDensityProfile( double,double,double );


		// Read in radii and densities
		void Load();

		// if one wants to use a non earth sized sphere...
		void   SetEarthRadiuskm( double x ) { REarth = x; };
		double GetEarthRadiuskm( ) { return  REarth; };

		// The next three functions are only available after a call to 
		// SetDensityProfile...
		
		// number of layers the current neutrino sees
		int get_LayersTraversed( ) {return Layers;};

		// self-explanatory
		double get_DistanceAcrossLayer( int i) { return _TraverseDistance[i];};
		double get_DensityInLayer( int i) { return _TraverseRhos[i];};
		double get_YpInLayer( int i)      { return _Yp[i];          };

		//DB Set Chemical composition in Prem model
		void SetChemicalComposition(int nChemCompVals, double* ChemCompVals);
		void SetChemicalComposition_PremModel(double Yp_Val);

                // return total path length through the sphere, including vacuum layers
                double get_Pathlength(){ 
			double Sum = 0.0;
			for( int i=0 ; i < Layers ; i++ )
				Sum += _TraverseDistance[i];
			return Sum;

		}	

                void LoadDensityProfile( const char * );
		
	private:
		// Computes the minimum pathlenth ( zenith angle ) a track needs to cross
		// each of the radial layers of the earth
		void ComputeMinLengthToLayers( );

                string DensityFileName;
			
		map<double, double>	_CosLimit;
                map<double, double>	_density;
		map<double, double>     _YpMap; 

		vector< double >	_Radii;
		vector< double >	_Rhos;
		vector< double >        _Yps ;

		double * _TraverseDistance;
		double * _TraverseRhos;
		double * _Yp;

		double REarth;
		int Layers;			
};


#endif



