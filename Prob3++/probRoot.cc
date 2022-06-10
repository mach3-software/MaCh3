#include <iostream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH2D.h"

int main(int argc, char * argv[] )
{

   double cosineZ, energy;
   double e_start, e_end, e_step, cz_start, cz_end, cz_step;
   int i, j ;

//// Types of Averaging	
   int NBinsEnergy = 200;
   int Dm12ZenithNbin = 200; 


// Zenith Angle Range
   double Dm12ZenithEdge[Dm12ZenithNbin+1];
   cz_start = -1.001;
   cz_end = 1.101;
   cz_step = ( cz_end - cz_start)/double(Dm12ZenithNbin);


// Energy Range
   double EnergyBins[NBinsEnergy+1];
   e_start = 0.110000001;
   e_end = 300.;
   e_step = log10(e_end/e_start)/double(NBinsEnergy);

   
/// Oscillation Parameters
   bool kSquared = false;   // are we using sin^2(x) variables?

   int    kNuBar  = 1;
   double DM2     = 2.5e-3;
   double Theta23 = 1.0;
   double Theta13 = 0.10;
   double dm2     = 7.9e-5;
   double Theta12 = 0.825;

   std::cout << "Using          " << std::endl
             << "      DM2      " <<  DM2      << std::endl
             << "      Theta23  " <<  Theta23  << std::endl
             << "      Theta13  " <<  Theta13  << std::endl
             << "      dm2      " <<  dm2      << std::endl
             << "      Theta12  " <<  Theta12  << std::endl;

   std::cout << "From "
	     << " [ " << e_start << " - " << e_end << " ] GeV " << endl;
   
// Methods to Compute Probability
   NeutrinoPropagator * myNu;
   BargerPropagator   * bNu; 

   bNu = new BargerPropagator( );
   bNu->UseMassEigenstates( false ); // this is default (use flavor eigenstates)

   // use the standard barger
   myNu = bNu;

   double Entry = e_start;
   for(i=0; i<NBinsEnergy; i++ )
        {
                Entry = e_start*pow( 10.0 , double(i)*e_step );
                EnergyBins[i] = Entry;
   }
   EnergyBins[NBinsEnergy] = EnergyBins[NBinsEnergy-1]*1.001;

   
   Dm12ZenithEdge[0]= cz_start*0.9999;
   for ( i=1; i<Dm12ZenithNbin ; i++ )
   	Dm12ZenithEdge[i] = Dm12ZenithEdge[0] + double(i)*cz_step;

   Dm12ZenithEdge[Dm12ZenithNbin] = Dm12ZenithEdge[Dm12ZenithNbin-1]*1.001;

   
   TH2D *NuEToNuE3f 	=  new TH2D("NuEToNuE3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuEToNuMu3f 	=  new TH2D("NuEToNuMu3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuEToNuTau3f 	=  new TH2D("NuEToNuTau3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuEToNuX3f 	=  new TH2D("NuEToNuX3f","3 Flavor P_{#nu_{e}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   

   
   TH2D *NuMuToNuE3f 	=  new TH2D("NuMuToNuE3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuMuToNuMu3f 	=  new TH2D("NuMuToNuMu3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuMuToNuTau3f 	=  new TH2D("NuMuToNuTau3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuMuToNuX3f 	=  new TH2D("NuMuToNuX3f","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);

   
   TH2D *NuTauToNuE3f 	=  new TH2D("NuTauToNuE3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{e}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuTauToNuMu3f 	=  new TH2D("NuTauToNuMu3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#mu}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuTauToNuTau3f 	=  new TH2D("NuTauToNuTau3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);
   TH2D *NuTauToNuX3f 	=  new TH2D("NuTauToNuX3f","3 Flavor P_{#nu_{#tau}#rightarrow#nu_{x}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);


   
   TH2D *NuMuToNuTau2f 	=  new TH2D("NuMuToNuTau2f","2 Flavor P_{#nu_{#mu}#rightarrow#nu_{#tau}}",
   					NBinsEnergy  -1 , EnergyBins, Dm12ZenithNbin -1, Dm12ZenithEdge);

//------------------- End of Raw Probability Plots	


   double total = 0.0;
   // fill the probabilities
   for ( i = 0 ; i <= NBinsEnergy ; i ++ ) 
   {

     energy = e_start*pow(10.0, double(i)*e_step);
     for ( j = 0 ; j <= Dm12ZenithNbin ; j++ )
     { 

        cosineZ = cz_start + double(j)*cz_step;
  
        myNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, 0.0 , energy, kSquared); // Compare with 3Flavor05.ps//
        myNu->DefinePath( cosineZ, 42.00  );
        myNu->propagate( 1*kNuBar );
  
        total = 0.0;
	total += myNu->GetProb(1,1);
	total += myNu->GetProb(1,2);
	total += myNu->GetProb(1,3);

	if( fabs( 1.00 - total ) > 1.0e-7 ) abort();
	
	
        NuEToNuE3f 	->Fill( energy, cosineZ, myNu->GetProb(1,1) ) ;
        NuEToNuMu3f 	->Fill( energy, cosineZ, myNu->GetProb(1,2) ) ;
        NuEToNuTau3f 	->Fill( energy, cosineZ, myNu->GetProb(1,3) ) ;
        NuEToNuX3f 	->Fill( energy, cosineZ, 1.0 - myNu->GetProb(1,1) ) ;

        NuMuToNuE3f 	->Fill( energy, cosineZ, myNu->GetProb(2,1) ) ;
        NuMuToNuMu3f 	->Fill( energy, cosineZ, myNu->GetProb(2,2) ) ;
        NuMuToNuTau3f 	->Fill( energy, cosineZ, myNu->GetProb(2,3) ) ;
        NuMuToNuX3f 	->Fill( energy, cosineZ, 1.0 - myNu->GetProb(2,2) ) ;

        NuTauToNuE3f 	->Fill( energy, cosineZ, myNu->GetProb(3,1) ) ;
        NuTauToNuMu3f 	->Fill( energy, cosineZ, myNu->GetProb(3,2) ) ;
        NuTauToNuTau3f 	->Fill( energy, cosineZ, myNu->GetProb(3,3) ) ;
        NuTauToNuX3f 	->Fill( energy, cosineZ, 1.0 - myNu->GetProb(3,3) ) ;
     }// End of Cosine Z Looping //
   
   } // End Energy Loop //


   TFile *tmp = new TFile("RawProb.root", "recreate");
   tmp->cd();

   NuEToNuE3f    ->Write();
   NuEToNuMu3f   ->Write();
   NuEToNuTau3f  ->Write();
   NuEToNuX3f    ->Write();
   
   NuMuToNuE3f   ->Write();
   NuMuToNuMu3f  ->Write();
   NuMuToNuTau3f ->Write();
   NuMuToNuX3f   ->Write();

   NuTauToNuE3f  ->Write();
   NuTauToNuMu3f ->Write();
   NuTauToNuTau3f->Write();
   NuTauToNuX3f  ->Write();

   NuMuToNuTau2f ->Write();

   tmp->Write();
   tmp->Close();
   	
   std::cout << std::endl<<"Done Cowboy!" << std::endl;
   return 0;
}

