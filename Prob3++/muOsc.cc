#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

double oscfunction(double* x, double* p) 
{
   //parameters
   //0=sin^2(2theta12)
   //1=sin^2(2theta23)
   //2=sin^2(theta13)
   //3=delta m_12
   //4=delta m_32
   //5=delta_cp
   //6=input species
   //7=output species
   //8=L
   //9=density
   
    cout << "E: " << x[0] << endl;
    cout << "Parameters: " << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "," << p[4] << "," << p[5] 
	 << "," << int(p[6]) << "," << int(p[7]) << "," << p[8] << "," << p[9] << endl;


   BargerPropagator   * bNu; 
   
   bNu = new BargerPropagator( );
   bNu->UseMassEigenstates( false );
   bNu->SetMNS( p[0],  p[2], p[1], p[3], p[4], p[5] , x[0], false);
   bNu->propagateLinear( int(p[6]), p[8], p[9] );
	
   double prob = bNu->GetProb(int(p[6]), int(p[7])); 

   cout << "Probability: " << prob << endl;

   delete bNu;
   return prob;
}

using namespace std;
int main(int argc, char * argv[] )
{

	double DM2 = 2.5e-3;
	double Theta23 = 1.0;
	double Theta13 = 0.1;
	double dm2 = 7.9e-5;
	double Theta12 = 0.825;
	double L = 295.335;
	double delta_cp = 0.0;
	double density = 2.7;

	std::cout << "Using          " << std::endl
		  << "      DM2      " <<  DM2      << std::endl
		  << "      Theta23  " <<  Theta23  << std::endl
		  << "      Theta13  " <<  Theta13  << std::endl
		  << "      dm2      " <<  dm2      << std::endl
		  << "      Theta12  " <<  Theta12  << std::endl;

		
	TF1* mutomu = new TF1("mutomu",oscfunction,0.01,10,10);
	mutomu->SetParameters(Theta12,Theta23,Theta13,dm2,DM2,delta_cp,2,2,L,density);
	mutomu->SetNpx(10000);

	TF1* mutoe = new TF1("mutoe",oscfunction,0.01,10,10);
	mutoe->SetParameters(Theta12,Theta23,Theta13,dm2,DM2,delta_cp,2,1,L,density);
	mutoe->SetNpx(10000);

	TF1* mutotau = new TF1("mutotau",oscfunction,0.01,10,10);
	mutotau->SetParameters(Theta12,Theta23,Theta13,dm2,DM2,delta_cp,2,3,L,density);
	mutotau->SetNpx(10000);

	TF1* mutoe_inverted = new TF1("mutoe_inverted",oscfunction,0.01,10,10);
	mutoe_inverted->SetParameters(Theta12,Theta23,Theta13,dm2,-1*DM2,delta_cp,2,1,L,density);
	mutoe_inverted->SetNpx(10000);

	cout << "mutomu: " << mutomu->Eval(0.6)
	     << "\tmutoe: " << mutoe->Eval(0.6)
	     << "\tmutotau: " << mutotau->Eval(0.6)
	     << "\tmutoe_inverted: " << mutoe_inverted->Eval(0.6) << endl;
	


	TFile *tmp = new TFile("muOsc.root", "recreate");
	tmp->cd();
	mutomu->Write();
	mutoe->Write();
	mutotau->Write();
	mutoe_inverted->Write();

	tmp->Write();
	tmp->Close();
		
	cout << endl<<"Done Cowboy!" << endl;
	return 0;
}

