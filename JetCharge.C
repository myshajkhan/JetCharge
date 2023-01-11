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
/// \author Andreas Morsch

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
#include "JetConf.h"
using namespace fastjet;  //fastjet is a library that builds jets
using namespace std;

//This class is to store the charge of the particles in FastJet
class MyInfo: public PseudoJet::UserInfoBase {
public:
  MyInfo(double charge, int pdgid) : _charge(charge), _pdg_id(pdgid){}
  double pdg_charge() const {return _charge;}
  double _charge;
  //mj inserting
  int pdg_id() const {return _pdg_id;}
  int _pdg_id;
  //mj done
};
class MyPDG: public PseudoJet::UserInfoBase {
public:
  MyPDG(int pdg) : _pdg(pdg){}
  int pdg_id() const {return _pdg;}
  int _pdg;
};


bool verbosity = 0; //printing out


void JetCharge(Int_t nev  = 100, int mode = 0,
              bool UE_switch = true, bool maydecay = true
              )
{

  TString str_UE = "";
  TString str_maydecay = "";
  TString str_mode = "";

  if(UE_switch) str_UE = "_UE";
  if(!maydecay) str_maydecay = "_noDecay";
  if(mode == 0) str_mode = "_all";
  else if(mode == 4) str_mode = "_c";
  else if(mode == 5) str_mode = "_b";
  else if(mode == 45) str_mode = "_bc";

  // choose a jet definition
  double R = 0.5; //in LHCb radius is fixed
  double kappa = 0.7;


  TString extension = TString("jetcharge")+Form("_ev_%d", nev)+
                      Form("_R0%d", int(R*10))+Form("_k0%d", int(kappa*10))+
                      str_UE+
                      str_maydecay+str_mode; //Name of output file
  // Load libraries
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  // Histograms

  int binlow = -3; //lower edge of jet charge histogram
  int binhigh = 3; // upper edge of jet charge histogram

  TFile f("hists/"+ extension+".root", "RECREATE");

  TH1F* etaH = new TH1F("etaH", "Pseudorapidity", 120, -12., 12.); // pseudorapidity of jets
  TH1F* ptH  = new TH1F("ptH",  "pt",              50,   0., 10.); // pT of jets


  // Histgrams for jet charge
  TH1D* h1_u_jetcharge_pT = new TH1D("pT u_jetcharge", "", 60, binlow, binhigh); // up quark histo
  TH1D* h1_ubar_jetcharge_pT= new TH1D("pT ubar_jetcharge", "", 60, binlow, binhigh); // up bar quark histo
  TH1D* h1_d_jetcharge_pT = new TH1D("pT d_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_dbar_jetcharge_pT = new TH1D("pT dbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_c_jetcharge_pT= new TH1D("pT c_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_cbar_jetcharge_pT = new TH1D("pT cbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_b_jetcharge_pT= new TH1D("pT b_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_bbar_jetcharge_pT = new TH1D("pT bbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_u_jetcharge_z = new TH1D("z u_jetcharge", "", 60, binlow, binhigh); // up quark histo
  TH1D* h1_ubar_jetcharge_z= new TH1D("z ubar_jetcharge", "", 60, binlow, binhigh); // up bar quark histo
  TH1D* h1_d_jetcharge_z = new TH1D("z d_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_dbar_jetcharge_z = new TH1D("z dbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_c_jetcharge_z= new TH1D("z c_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_cbar_jetcharge_z = new TH1D("z cbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_b_jetcharge_z= new TH1D(" z b_jetcharge", "", 60, binlow, binhigh); // down quark histo
  TH1D* h1_bbar_jetcharge_z = new TH1D(" zbbar_jetcharge", "", 60, binlow, binhigh); // down bar histo
  TH1D* h1_frag_charge = new TH1D("even_frag_charge", "", 6, -2.5, 3.5);
  TH2* pTvsz = new TH2F("ptvsz", "h2 title", 30, 0, 2.0, 30, 0.0, 2.0);




  JetDefinition jet_def(antikt_algorithm, R); // Anti-kT algorithm with radius R
  vector<PseudoJet> parts; // vector of PseudoJets to store the particles. Each PseudoJet is a four-vector (px, py, pz, e)

  // Array of particles
  TClonesArray* particles = new TClonesArray("TParticle", 1000);
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
  if(mode == 0) pythia8->ReadString("HardQCD:all = on") ;
  else if (mode == 4) pythia8->ReadString("HardQCD:hardccbar = on");
  else if (mode == 5) pythia8->ReadString("HardQCD:hardbbbar = on");
  else if (mode == 45){
    pythia8->ReadString("HardQCD:hardbbbar = on");
    pythia8->ReadString("HardQCD:hardccbar = on");
  }
   // turing on all head on collisions from QCD


  pythia8->ReadString("Random:setSeed = on");
  // use a reproducible seed: always the same results for the tutorial.
  pythia8->ReadString("Random:seed = 42");
  pythia8->ReadString("PhaseSpace:pTHatMin = 20.0");

  // Turn off underlying event
  if(!UE_switch){
    pythia8->ReadString("PartonLevel:MPI = off"); // MultipartonInteractions
    pythia8->ReadString("PartonLevel:ISR = off"); // Initial State Radiation
    pythia8->ReadString("PartonLevel:FSR = off"); // Final State Radiation
  }
  if(!maydecay){
    pythia8->ReadString("511:mayDecay = off");
    pythia8->ReadString("521:mayDecay = off");
    pythia8->ReadString("10511:mayDecay = off");
    pythia8->ReadString("10521:mayDecay = off");
    pythia8->ReadString("513:mayDecay = off");
    pythia8->ReadString("523:mayDecay = off");
    pythia8->ReadString("-511:mayDecay = off");
    pythia8->ReadString("-521:mayDecay = off");
    //pythia8->ReadString("-10511:mayDecay = off");
    //pythia8->ReadString("-10521:mayDecay = off");
    pythia8->ReadString("-513:mayDecay = off");
    pythia8->ReadString("-523:mayDecay = off");

    pythia8->ReadString("411:mayDecay = off");
    pythia8->ReadString("421:mayDecay = off");
    //pythia8->ReadString("10411:mayDecay = off");
    //pythia8->ReadString("10421:mayDecay = off");
    pythia8->ReadString("413:mayDecay = off");
    pythia8->ReadString("423:mayDecay = off");
    pythia8->ReadString("-411:mayDecay = off");
    pythia8->ReadString("-421:mayDecay = off");
    //pythia8->ReadString("-10411:mayDecay = off");
    //pythia8->ReadString("-10421:mayDecay = off");
    pythia8->ReadString("-413:mayDecay = off");
    pythia8->ReadString("-423:mayDecay = off");

  }



  // Initialize

  pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */); //proton proton collision at 13 TeV

  // Event loop
  for (Int_t iev = 0; iev < nev; iev++) {
    pythia8->GenerateEvent();

    //if (iev < ndeb) pythia8->EventListing();
    //if (verbosity) { pythia8->EventListing(); }
    pythia8->ImportParticles(particles,"All");
    Int_t np = particles->GetEntriesFast();
    parts.clear();
    // Particle loop
    TLorentzVector parton1, parton2; // Four-vectors for outgoing quarks/gluons
    int pdg_parton1, pdg_parton2; // pdg code for outgoing quarks/gluons
    for (Int_t ip = 0; ip < np; ip++) {
      TParticle* part = (TParticle*) particles->At(ip); // creating the T particle
      Int_t ist = part->GetStatusCode(); //not inportant probably

      //if the particle is the first outgoing quark/gluon, store its info
      if(ip == 4){
        parton1.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy() );
        pdg_parton1 = part->GetPdgCode();
      }

      //if the particle is the second outgoing quark/gluon, store its info
      if(ip == 5){
        parton2.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy() );
        pdg_parton2 = part->GetPdgCode();
      } 
      // Positive codes are final particles.
      if (ist <= 0) continue; //  we get it from previous line and  greater than 0 means it's not an intermidiate particle
      Int_t pdg = part->GetPdgCode();
      if(abs(pdg/100 %10) == 5){
      h1_frag_charge->Fill(pdg/10 %10);
      }
      Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      //if (charge == 0.) continue;
      Float_t eta = part->Eta();
      Float_t pt  = part->Pt();
      if (pt < 0.2) continue;
      parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
      //cout<<"("<<pdg << "," << charge<<"),";
      parts.back().set_user_info(new MyInfo(charge, pdg));
  
      etaH->Fill(eta);
      if (pt > 0.) ptH->Fill(pt, 1./(2. * pt));
    }

    
    ClusterSequence cs(parts, jet_def); // building the jets in the events
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // Store all the jets in the event
    if(jets.size()<2) continue; //Throw away events with less than two jets
    //cout << " num jets : " << jets.size() << endl;
    int NumJets = 0; // Not important for now
    // Jet loop
    int pdg_jet = -99; // variable to store which flavor the jet came from (up jet? down jet? etc)
    bool has_parton1 = false; // boolean to tell whether the jet contained the first outgoing quark/gluon
    bool has_parton2 = false;// boolean to tell whether the jet contained the second outgoing quark/gluon

    for (unsigned i = 0; i < jets.size(); i++) //Loop over the jets!
    {
      if (jets[i].pt() < 20) continue; // get rid of low pT jets
      //if (jets[i].eta() < 2 || jets[i].eta() > 4.5) continue;
      if(fabs(jets[i].eta()) > 5 || fabs(jets[i].eta()) < 2) continue;
      //cout << "jet num : " <<  i << endl;
      PseudoJet jet = jets[i];
      TLorentzVector jetvec(jet.px(), jet.py(), jet.pz(), jet.e()); // Four vector of the jet
      if(jetvec.DeltaR(parton1) < 0.1) has_parton1 = true;  //Check whether the direction of the jet matches the first outgoing parton
      else if(jetvec.DeltaR(parton2) < 0.1) has_parton2 = true; //Check whether the direction of the jet matches the second outgoing parton
      TVector3 pJet(jet.px(),jet.py(),jet.pz());// defining pjet for the denominator of
      double pJetmag=pJet.Mag();//magnitude of pJet
      


      if(has_parton1&&(!has_parton2)){// If this is true, then the jet came from the first outgoing parton
        pdg_jet = pdg_parton1;
      }
      else if(has_parton2 && (!has_parton1)){// If this is true, then the jet came from the second outgoing parton
        pdg_jet = pdg_parton2;
      }
    
      vector<PseudoJet> constituents = jet.constituents(); //give me jet daughters
      if( constituents.size() >2){
      double jetcharge_pT = 0;
      double jetcharge_z = 0;
      double PT = 0;
      double Z = 0;
      double jetcharge_z_norm=0;
      for (unsigned j = 0; j < constituents.size(); j++) // knowing what each jet daugther property for us to maybe build jet charge
      {
      
       
        double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();
        PseudoJet con = constituents[j];
        TVector3 con3(con.px(), con.py(), con.pz());
	      double pjetCon= pJet.Dot(con3); //the numerator of z
	      double z= pjetCon/pow(pJet.Mag(),2);
        
        jetcharge_z+=pow(z, kappa)*dtrcharge/3.; //compute jet charge
        jetcharge_z_norm+= pow(z,kappa);
        //double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();
        //PseudoJet con = constituents[j];
        //TVector3 con3(con.px(), con.py(), con.pz());
        jetcharge_pT+=pow(con.pt(), kappa)*dtrcharge/3.; //compute jet charge

       // double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();
       // PseudoJet con = constituents[j];
        //TVector3 con3(con.px(), con.py(), con.pz());
        //double pjetCon= pJet.Dot(con3); //the numerator of z
        //double z= pjetCon/pow(pJet.Mag(),2);

        PT=con.pt() /jets[i].pt();
        pTvsz->Fill(z,PT);
      }
      //jetcharge_z/=pow(pJetmag, kappa); //normalize by jet pT
  
      jetcharge_pT/=pow(jets[i].pt(), kappa); //normalize by jet pT
       

      //Fill histograms
      if(pdg_jet == 1) h1_d_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == -1) h1_dbar_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == 2) h1_u_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == -2) h1_ubar_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == 4) h1_c_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == -4) h1_cbar_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == 5) h1_b_jetcharge_z->Fill(jetcharge_z);
      else if(pdg_jet == -5) h1_bbar_jetcharge_z->Fill(jetcharge_z);
      if(pdg_jet == 1) h1_d_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == -1) h1_dbar_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == 2) h1_u_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == -2) h1_ubar_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == 4) h1_c_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == -4) h1_cbar_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == 5) h1_b_jetcharge_pT->Fill(jetcharge_pT);
      else if(pdg_jet == -5) h1_bbar_jetcharge_pT->Fill(jetcharge_pT);
      
      }
        
     
   if( constituents.size() >2){
     

      for (unsigned j = 0; j < constituents.size(); j++) // knowing what each jet daugther property for us to maybe build jet charge
      {


       
      }
   }

    }

  }
 

    
   
    //The following doesn't matter much, skip ahead!!
    int ican=-1,iframe=-1,itext=-1;
    TCanvas *ccan[1000];
    TH1F	*frame[1000];
    TLatex	*text[1000];
    for (int i=0;i<1000;i++){
      text[i]	= new TLatex();
      text[i]->SetNDC(kTRUE);
      text[i]->SetTextSize(0.06);
    }
    TLatex Tl;
    Tl.SetNDC(kTRUE);
    Tl.SetTextSize(0.04);
    //
    //gStyle->SetOptStat(0);
    //gStyle->SetPaperSize(TStyle::kUSLetter);
    //gStyle->SetPadBottomMargin(0.08);
    //gStyle->SetPadTopMargin(0.005);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadRightMargin(0.13);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetTitleXSize(0.055);
    gStyle->SetTitleYSize(0.055);
    gStyle->SetTitleOffset(0.85,"X");
    gStyle->SetTitleOffset(1.2,"Y");
    gStyle->SetStatW(0.2);
    //gStyle->SetPalette(1);
    //gStyle->SetErrorX(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetStatStyle(0);
    //gStyle->SetLineWidth(3);

    //---- paint...
    char buf[100];
    char bufb[100];
    TString rootfile;
    TString plotfile;
    TString plotfilePDF;
    TString plotfileO;
    TString plotfileC;
    //TString OutputFileBase	= outbase+outinfo;
    rootfile	= TString("./hists/") + extension + TString(".root");
    plotfile	= TString("./plots/")   + extension + TString(".ps");
    plotfilePDF	= TString("./plots/")   + extension + TString(".pdf");
    plotfileO	= plotfilePDF + TString("(");
    plotfileC	= plotfilePDF + TString("]");
    //c->SaveAs("plots/"+extension+".pdf");
    
  
    
    // Start new page!!!!
   
    ++ican;
    sprintf(buf,"ccan%d",ican);
    ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
    ccan[ican]->SetFillColor(10);
    //gPad->SetLeftMargin(0.16);
    //gPad->SetBottomMargin(0.06);
    ccan[ican]->cd(); ccan[ican]->Divide(2,2,0.0001,0.0001);
    ccan[ican]->cd(1);
    //plot stuff here!!!!
 
    h1_u_jetcharge_pT->SetXTitle("Jet Charge");
    h1_u_jetcharge_pT->SetLineColor(kBlue);
    h1_u_jetcharge_pT->SetLineStyle(kSolid);

    h1_ubar_jetcharge_pT->SetLineColor(kGreen);
    h1_ubar_jetcharge_pT->SetLineStyle(kDashed);

    h1_u_jetcharge_pT->Draw("HIST same");
    h1_ubar_jetcharge_pT->Draw("HIST same");
    h1_u_jetcharge_z->SetXTitle("Jet Charge");
    h1_u_jetcharge_z->SetLineColor(kBlack);
    h1_u_jetcharge_z->SetLineStyle(kSolid);

    h1_ubar_jetcharge_z->SetLineColor(kRed);
    h1_ubar_jetcharge_z->SetLineStyle(kDashed);

    h1_u_jetcharge_z->Draw("HIST same");
    h1_ubar_jetcharge_z->Draw("HIST same");
    

    ccan[ican]->cd(2);
   

    h1_d_jetcharge_pT->SetLineColor(kBlue);
    h1_d_jetcharge_pT->SetLineStyle(kSolid);

    h1_dbar_jetcharge_pT->SetLineColor(kGreen);
    h1_dbar_jetcharge_pT->SetLineStyle(kDashed);

    h1_d_jetcharge_pT->Draw("HIST same");
    h1_dbar_jetcharge_pT->Draw("HIST same");
    h1_d_jetcharge_z->SetLineColor(kBlack);
    h1_d_jetcharge_z->SetLineStyle(kSolid);

    h1_dbar_jetcharge_z->SetLineColor(kRed);
    h1_dbar_jetcharge_z->SetLineStyle(kDashed);

    h1_d_jetcharge_z->Draw("HIST same");
    h1_dbar_jetcharge_z->Draw("HIST same");


    ccan[ican]->cd(3);
    //plot stuff here!!!!

    
    h1_c_jetcharge_pT->SetXTitle("Jet Charge");
    h1_c_jetcharge_pT->SetLineColor(kBlue);
    h1_c_jetcharge_pT->SetLineStyle(kSolid);

    h1_cbar_jetcharge_pT->SetLineColor(kGreen);
    h1_cbar_jetcharge_pT->SetLineStyle(kDashed);

    h1_c_jetcharge_pT->Draw("HIST same");
    h1_cbar_jetcharge_pT->Draw("HIST same");

    h1_c_jetcharge_z->SetXTitle("Jet Charge");
    h1_c_jetcharge_z->SetLineColor(kBlack);
    h1_c_jetcharge_z->SetLineStyle(kSolid);

    h1_cbar_jetcharge_z->SetLineColor(kRed);
    h1_cbar_jetcharge_z->SetLineStyle(kDashed);

    h1_c_jetcharge_z->Draw("HIST same");
    h1_cbar_jetcharge_z->Draw("HIST same");


    ccan[ican]->cd(4);

    
    h1_b_jetcharge_pT->SetLineColor(kBlue);
    h1_b_jetcharge_pT->SetLineStyle(kSolid);

    h1_bbar_jetcharge_pT->SetLineColor(kGreen);
    h1_bbar_jetcharge_pT->SetLineStyle(kDashed);

    h1_b_jetcharge_pT->Draw("HIST same");
    h1_bbar_jetcharge_pT->Draw("HIST same");

    h1_b_jetcharge_z->SetLineColor(kBlack);
    h1_b_jetcharge_z->SetLineStyle(kSolid);

    h1_bbar_jetcharge_z->SetLineColor(kRed);
    h1_bbar_jetcharge_z->SetLineStyle(kDashed);

    h1_b_jetcharge_z->Draw("HIST same");
    h1_bbar_jetcharge_z->Draw("HIST same");
    ccan[ican]->cd();ccan[ican]->Update();
    if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
       else { ccan[ican]->Print(plotfilePDF.Data()); }
    
       

       ++ican;
      sprintf(buf,"ccan%d",ican);
      ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
      ccan[ican]->SetFillColor(10);
      ccan[ican]->cd(); ccan[ican]->Divide(2,2,0.0001,0.0001);
      ccan[ican]->cd(1);
    //plot stuff here!!!!
    h1_u_jetcharge_z->SetXTitle("Jet Charge");
    h1_u_jetcharge_z->SetLineColor(kBlack);
    h1_u_jetcharge_z->SetLineStyle(kSolid);

    h1_ubar_jetcharge_z->SetLineColor(kRed);
    h1_ubar_jetcharge_z->SetLineStyle(kDashed);

    h1_u_jetcharge_z->Draw("HIST same");
    h1_ubar_jetcharge_z->Draw("HIST same");
    
    

    ccan[ican]->cd(2);
    h1_d_jetcharge_z->SetLineColor(kBlack);
    h1_d_jetcharge_z->SetLineStyle(kSolid);

    h1_dbar_jetcharge_z->SetLineColor(kRed);
    h1_dbar_jetcharge_z->SetLineStyle(kDashed);

    h1_d_jetcharge_z->Draw("HIST same");
    h1_dbar_jetcharge_z->Draw("HIST same");

   


    ccan[ican]->cd(3);
    //plot stuff here!!!!
    h1_c_jetcharge_z->SetXTitle("Jet Charge");
    h1_c_jetcharge_z->SetLineColor(kBlack);
    h1_c_jetcharge_z->SetLineStyle(kSolid);

    h1_cbar_jetcharge_z->SetLineColor(kRed);
    h1_cbar_jetcharge_z->SetLineStyle(kDashed);

    h1_c_jetcharge_z->Draw("HIST same");
    h1_cbar_jetcharge_z->Draw("HIST same");
    
   

    ccan[ican]->cd(4);
    h1_b_jetcharge_z->SetLineColor(kBlack);
    h1_b_jetcharge_z->SetLineStyle(kSolid);

    h1_bbar_jetcharge_z->SetLineColor(kRed);
    h1_bbar_jetcharge_z->SetLineStyle(kDashed);

    h1_b_jetcharge_z->Draw("HIST same");
    h1_bbar_jetcharge_z->Draw("HIST same");
    
   

    ccan[ican]->cd();ccan[ican]->Update();
    if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
       else { ccan[ican]->Print(plotfilePDF.Data()); }
    
    



      ++ican;
    sprintf(buf,"ccan%d",ican);
    ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
    ccan[ican]->SetFillColor(10);
    //gPad->SetLeftMargin(0.16);
    //gPad->SetBottomMargin(0.06);
    ccan[ican]->cd(); ccan[ican]->Divide(2,2,0.0001,0.0001);
   


    ccan[ican]->cd(1);
    h1_frag_charge->Draw();



    ccan[ican]->cd(2);
    pTvsz->Draw();

    ccan[ican]->cd();ccan[ican]->Update();
    if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
       else { ccan[ican]->Print(plotfilePDF.Data()); }





    f.Write();
    f.Close();


    if (ican>-1){
    cout<<" You plotted "<<ican+1<<" canvasses......."<<endl;
    ccan[ican]->Print(plotfileC.Data());
    }


}
