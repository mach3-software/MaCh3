#ifndef _BargerPropagator_
#define _BargerPropagator_

#include "EarthDensity.h"
#include "NeutrinoPropagator.h"

#include "mosc3.h"
#include "mosc.h"


#include <iostream>
#include <cstdlib>

// A neutrino propagator class based on the the 1980 Barger paper PRD.22.11, Dec. 1 1980
// The underlying oscillation code is written in mosc*
//
// Capable of computing oscillation probabilities through constant denisity matter
// or through sphere's of varying radial density.


class BargerPropagator : public NeutrinoPropagator
{
  public:
      
      BargerPropagator( );
      // const char specifies an input radial density profile, c.f. PREM.dat
      BargerPropagator( bool );
      BargerPropagator( const char *);
     ~BargerPropagator(      );
      
      // main driving routine for computing oscillations through a sphere
      // called after SetMNS(...) !!
      // specify neutrino type:   +int : neutrino   -int: anti-neutrino 
      virtual void propagate( int );
      
      // driving routine for oscillations through linear media of contstant density 
      // called after SetMNS(...)
      // specify neutrino type:   +int : neutrino   -int: anti-neutrino 
      // specify Path length in the matter
      // specify density of the matter
      virtual void propagateLinear( int , double, double );
      
      // driving routine for oscillations in vaccuum!
      // called after SetMNS(...)
      // specify           nu_in, nu_out, Energy [GeV] , pathlength [km]  
      virtual double GetVacuumProb( int , int   , double , double );
      
      // determines the pathlength and density profile for a neutrino propagating through a sphere
      // specify, cosine of zenith angle   -1 : upward going 0: horizontal +1: downward going
      // specify production height in the the atmosphere [km]
      // specify if the profile withing EarthDensity object should be recomputed, default is true
      virtual void DefinePath( double, double, bool kSetProfile = true  );
      
      // determine the neutrino oscillation parameters
      // This routine must be called _before_ propagate* routines!
      // Specify the neutrino oscillation parameters, and energy
      // the final boolean specifies which form  of mixing angle is input
      //            x12   ,  x13   ,  x23   ,  dm21  ,  dm32  ,  d_cp  , Energy [GeV],  T: sin^2(x) F: sin^2(2x)
      // 
      //  The last argument is the neutrino type  nu > 0 : neutrinos  nu < 0 : antineutrinos
      //  The default is to define the MNS matrix for neutrino propagation
      //  This type must agree with the type used in the call to propagate() and propagateLinear()
      virtual void SetMNS( double , double , double , double , double , double , double , bool, int kNuType = 1 );
      
      // for changing the conversion factor from matter density to electron density
      void SetDensityConversion( double x ) { density_convert = x; }
      
      //DB Passthrough information for setting chemical composition of Earth's core
      void SetChemicalComposition_PremModel(double Yp_Val=0.468) {
	Earth->SetChemicalComposition_PremModel(Yp_Val);
      }
      
      // return oscillation probabilities nu_in -> nu_out
      // nu_ - 1:e 2:mu 3:tau  -1:e_bar -2:mu_bar -3:tau_bar
      double GetProb( int nuIn, int nuOut ){ 
     	int In = abs( nuIn );
	int Out = abs( nuOut );
	return Probability[In-1][Out-1];
      };
      
      
      // miscellaneuos
      double GetPathLength(){return Earth->get_Pathlength();}
      void SetPathLength( double x ){ PathLength = x;}
      void SetEnergy    ( double x ){ Energy     = x;}   

      virtual void SetMatterPathLength();
      virtual void SetAirPathLength(double);
      
      
      // Specify weather oscillition probabilities are computed from neutrino mass eigenstates
      // of from neutrino flavor eigen states   T: mass  F: flavor
      void UseMassEigenstates( bool x ) { kUseMassEigenstates = x ;}
      
      void SetWarningSuppression( bool x = true ) { kSuppressWarnings = x ; }

      // Specify how the user inputs the atmospheric neutrino mixing mass 
      // true (default mode) means the mixing input for SetMNS corresponds to
      //    NH: m32 
      //    IH: m31   
      //  That is, in this mode the code will correct the input value of the 
      //  atmospheric mass splitting parameter by the solar mass splitting if 
      //  the input is negative (corresponding to IH input).
      //  If SetOneMassScaleMode(false) is called this correction is not 
      //  performed, and the user is responsible for supplying the correct 
      //  value of Dm23 for oscillations in both hierarchies.  
      void SetOneMassScaleMode  ( bool x = true)  { kOneDominantMass  = x ; }

  protected:
      void init();
      
      void ClearProbabilities();
      
      double Probability[3][3];
      
      EarthDensity * Earth;		
      
      double REarth;
      double ProductionHeight;
      double PathLength;
      double AirPathLength;
      double MatterPathLength;
      double CosineZenith;
      double Energy;	
      double density_convert;
      
      
      bool kAntiMNSMatrix     ;
      bool kSuppressWarnings  ; 
      bool kUseMassEigenstates;
      bool kOneDominantMass   ;  
      
};
      
#endif
      
