//#!/bin/csh -f 
/// \file
/// \ingroup tutorial_pythia
/// pythia8 basic example
///
/// to run, do:
///
/// ~~~{.cpp}
///  root > .x pythia8.C
/// ~~~
///
/// \macro_code
///

#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TVector3.h"

using namespace fastjet;  //fastjet is a library that builds jets 
using namespace std;      //std means standard. set of data types and commands prepackaged with the language
class MyInfo: public PseudoJet::UserInfoBase {
  public:
    MyInfo(double charge) : _charge(charge){}
    double pdg_charge() const {return _charge;}
    double _charge;
};

bool verbosity = 1; //printing out 
void pythia8Jets(Int_t nev  = 10000, Int_t ndeb = 1) // nev= number of events . We dont know what deb means lol
{


// Load libraries
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");

// Histograms
   TH1F* etaH = new TH1F("etaH", "Pseudorapidity", 120, -12., 12.);
   TH1F* ptH  = new TH1F("ptH",  "pt",              50,   0., 10.);
   TH1D* h= new TH1D("h", "1st Histogram",30, -3 ,3); //this is mj's histogram for jet charge.
	 h->SetXTitle("Jet charge"); //x axis for h histogram
	 h->SetYTitle("jet number");//y axis for h histogram 

  // choose a jet definition
  double R = 0.5; //in LHCb radius is fixed
  JetDefinition jet_def(antikt_algorithm, R);
  vector<PseudoJet> parts;

// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 1000);//Tparticle a dynamic particle class created by event generators and used during
// Create pythia8 object
   TPythia8* pythia8 = new TPythia8();

#if PYTHIA_VERSION_INTEGER == 8235
   // Pythia 8.235 is known to cause crashes:
   printf("ABORTING PYTHIA8 TUTORIAL!\n");
   printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
   printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
   return;
#endif

// Configure
   pythia8->ReadString("HardQCD:qq2qq = on") ; // turing on all head on collisions from QCD
   pythia8->ReadString("Random:setSeed = on");//for simulation, we use random number genarator 
   // use a reproducible seed: always the same results for the tutorial.
   pythia8->ReadString("Random:seed = 42");
   pythia8->ReadString("PhaseSpace:pTHatMin = 35.0");


// Initialize

   pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);

// Event loop
   for (Int_t iev = 0; iev < nev; iev++) {
      pythia8->GenerateEvent();
      
      //if (iev < ndeb) pythia8->EventListing();
    //  if (verbosity) { pythia8->EventListing(); }
      pythia8->ImportParticles(particles,"All");
      Int_t np = particles->GetEntriesFast();
      parts.clear();
     bool is_u_quark = false;
TParticle* part4 = (TParticle*) particles->At(4);
TParticle* part5 = (TParticle*) particles->At(5);

 Int_t pdg4 = part4->GetPdgCode();
 Int_t pdg5 = part5->GetPdgCode();
 if (pdg4!=pdg5)continue;
if (pdg4!=2)continue;
// Particle loop
      for (Int_t ip = 0; ip < np; ip++) {
	
         TParticle* part = (TParticle*) particles->At(ip); // creating the T particle
         Int_t ist = part->GetStatusCode(); //very inportant probably
	//Int_t status= particles->At(ip)->status();      
// 	bool is_u_quark = false;
     // cout<<ist<<",";
     //

/*
 Int_t pdg = part->GetPdgCode();
	//cout<<"about to eneter if";
	if(ist==0)

{
//cout<< " hellloooo this is running!!!!!";
//	Int_t pdg = part->GetPdgCode();
	if (pdg ==2){
		cout<<"found u quark";
  		is_u_quark=true ;
}

}*/
          // Positive codes are final particles.
         if (ist <= 0) continue; //  we get it from previous line and  greater than 0 means it's not an intermidiate particle
         Int_t pdg = part->GetPdgCode();
         Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
         //if (charge == 0.) continue;
         Float_t eta = part->Eta();
         Float_t pt  = part->Pt();
	 if (pt < 0.2) continue;
	 parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
	 parts.back().set_user_info(new MyInfo(charge));
         etaH->Fill(eta);
         if (pt > 0.) ptH->Fill(pt, 1./(2. * pt));
      }
//if(!is_u_quark) continue;

      //cout << " parts size : " << parts.size() << endl;

      ClusterSequence cs(parts, jet_def); // building the jets in the events
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // 
      //cout << " num jets : " << jets.size() << endl;
	int NumJets = 0;
      // Jet loop

      for (unsigned i = 0; i < jets.size(); i++) //
      {
          if (jets[i].pt() < 5) continue;
          cout << "jet num : " <<  i << endl;
          PseudoJet jet = jets[i];
	  TVector3 jet3(jet.px(), jet.py(), jet.pz());

          vector<PseudoJet> constituents = jet.constituents(); //give me jet daughter
	      double jetcharge = 0;
            double kappa = 1.0;

          for (unsigned j = 0; j < constituents.size(); j++) // knowing what each jet daugther property for us to maybe build jet charge
	  {
	    double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();
	    PseudoJet con = constituents[j];
            TVector3 con3(con.px(), con.py(), con.pz());
  // jetdtrs created
//      	    double jetcharge = 0;
  //
//    	    double kappa = 1.0;
 
      //if(jetdtrs.size()!=jetdtrscharge.size()) continue;
      //for(int i = 0; i < jetdtrs.size(); i++){
       	   jetcharge+=pow(con.pt(), kappa)*dtrcharge/3.; 
	   //jetcharge/=pow(jets[i].pt(), kappa);
		// we are trying to make Histograms
		//h= new TH1D("h", "1st Histogram",15, 0 ,15);

	//	for (int i=0; i<1000; i++){
  //      	h->Fill(jetcharge));

	//	}


//	h->Draw();

//	}

	//   h1_hfcharge->Fill(hfcharge/3.);
     	 //  h1_jetcharge->Fill(jetcharge);
	//   if(truthLevel) h1_hfchargetrue->Fill(hfchargetrue.at(i));
	  
//            double z = (con3.Dot(jet3))/(jet3.Mag2()); // i dont want this i want jet charge
	  /*  if (verbosity==1)
              {
	        cout << "    constituent " << j << "'s mass: " << con.m() << " (px,py,pz,E) : ( " << con.px() << ", " << con.py() << ", " << con.pz() << ", " << con.e() << " )" << endl;
                cout << " z = : " << z << endl;
              }*/
	  }
  
//if(!is_u_quark)continue;
	jetcharge/=pow(jets[i].pt(), kappa);
        NumJets++;	
	 h->Fill(jetcharge);
	 
      }
    
   }
    //h->Draw();
   pythia8->PrintStatistics();
   TCanvas c1("c1", "c1");
   h->Draw();
   c1.SaveAs("plot.pdf");

}

