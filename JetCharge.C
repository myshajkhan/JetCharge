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
#include "JetConf.h"
using namespace fastjet;  //fastjet is a library that builds jets
using namespace std;

//This class is to store the charge of the particles in FastJet
class MyInfo: public PseudoJet::UserInfoBase {  // what is a pseudoJet
public:
  MyInfo(double charge) : _charge(charge){}
  double pdg_charge() const {return _charge;}
  double _charge;
 // mj copied this from the internet because they want to know the pdg code of jet dtrs that's making the jet charge 0
  MyInfo(int id) : _pdg_id(id){}
  int pdg_id() const {return _pdg_id;}
  int _pdg_id;
};



 


double kappa = 1;
bool verbosity = 0; //printing out


void JetCharge(Int_t nev  = 1000000, Int_t ndeb = 1) // nev= number of events . We dont know what deb means lol
{


  TString extension = "jetcharge"; //Name of output file
  // Load libraries
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");
  // Histograms

  int binlow = -2; //lower edge of jet charge histogram
  int binhigh = 2; // upper edge of jet charge histogram

  TFile f("hists/"+ extension+".root", "RECREATE");

  TH1F* etaH = new TH1F("etaH", "Pseudorapidity", 120, -12., 12.); // pseudorapidity of jets
  TH1F* ptH  = new TH1F("ptH",  "pt",              50,   0., 10.); // pT of jets


  // Histgrams for jet charge
  TH1D* h1_u_jetcharge = new TH1D("h1_u_jetcharge", "", 300, binlow, binhigh); // up quark histo
  TH1D* h1_ubar_jetcharge = new TH1D("h1_ubar_jetcharge", "", 300, binlow, binhigh); // up bar quark histo
  TH1D* h1_d_jetcharge = new TH1D("h1_d_jetcharge", "", 300, binlow, binhigh); // down quark histo
  TH1D* h1_dbar_jetcharge = new TH1D("h1_dbar_jetcharge", "", 300, binlow, binhigh); // down bar histo
  TH1D* h1_g_jetcharge = new TH1D("h1_g_jetcharge", "", 300, binlow, binhigh); // gluon histo
  TH1D* h1_s_jetcharge = new TH1D("h1_s_jetcharge", "", 300, binlow, binhigh); // strange bar histo
  TH1D* h1_sbar_jetcharge = new TH1D("h1_sbar_jetcharge", "", 300, binlow, binhigh); // strage bar histo
  TH1D* h1_b_jetcharge = new TH1D("h1_b_jetcharge", "", 300, binlow, binhigh); // down bar histo
  TH1D* h1_bbar_jetcharge = new TH1D("h1_bbar_jetcharge", "", 300, binlow, binhigh); // gluon histo
  TH1D* h1_c_jetcharge = new TH1D("h1_c_jetcharge", "", 300, binlow, binhigh); // strange bar histo
  TH1D* h1_cbar_jetcharge = new TH1D("h1_cbar_jetcharge", "", 300, binlow, binhigh); // strage bar histo
  // choose a jet definition
  double R = 0.5; //in LHCb radius is fixed
  JetDefinition jet_def(antikt_algorithm, R); // Anti-kT algorithm with radius R
  vector<PseudoJet> parts; // vector of PseudoJets to store the particles. Each PseudoJet is a four-vector (px, py, pz, e)

  // Array of particles
  TClonesArray* particles = new TClonesArray("TParticle", 10000);////Tparticle a dynamic particle class created by event generators and used during

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
 



  pythia8->ReadString("HardQCD:all = on") ; // turing on all head on collisions from QCD
//  pythia8->ReadString("HardQCD:hardbbbar = on") ; // turing on bbbar head on collisions from QCD
//  pythia8->ReadString("HardQCD:hardccbar = on") ; // turing on ccbar head on collisions from QCD
  pythia8->ReadString("Random:setSeed = on");
  // use a reproducible seed: always the same results for the tutorial.
  pythia8->ReadString("Random:seed = 42");
  pythia8->ReadString("PhaseSpace:pTHatMin = 20.0");

  // Turn off underlying event
 pythia8->ReadString("PartonLevel:MPI = off"); // MultipartonInteractions  what are these 3 lines??
 pythia8->ReadString("PartonLevel:ISR = off"); // Initial State Radiation
 pythia8->ReadString("PartonLevel:FSR = off"); // Final State Radiation
 pythia8->ReadString("511:mayDecay = off");
 pythia8->ReadString("521:mayDecay = off");
 pythia8->ReadString("10511:mayDecay = off");
 pythia8->ReadString("10521:mayDecay = off");
 pythia8->ReadString("513:mayDecay = off");
 pythia8->ReadString("523:mayDecay = off");
 pythia8->ReadString("-511:mayDecay = off");
 pythia8->ReadString("-521:mayDecay = off");
 pythia8->ReadString("-10511:mayDecay = off");
 pythia8->ReadString("-10521:mayDecay = off");
 pythia8->ReadString("-513:mayDecay = off");
 pythia8->ReadString("-523:mayDecay = off");
 pythia8->ReadString("411:mayDecay = off");
 pythia8->ReadString("421:mayDecay = off");
 pythia8->ReadString("10411:mayDecay = off");
 pythia8->ReadString("10421:mayDecay = off");
 pythia8->ReadString("413:mayDecay = off");
 pythia8->ReadString("423:mayDecay = off");
 pythia8->ReadString("-411:mayDecay = off");
 pythia8->ReadString("-421:mayDecay = off");
 pythia8->ReadString("-10411:mayDecay = off");
 pythia8->ReadString("-10421:mayDecay = off");
 pythia8->ReadString("-413:mayDecay = off");
 pythia8->ReadString("-423:mayDecay = off");

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
/*       if (ip < 10) {
 	cout << "ip: " << ip << " PDG Code: " << part->GetPdgCode() << endl;
       }
*/
      //if the particle is the first outgoing quark/gluon, store its info
      if(ip == 4){ // 4 and 5 generated by the collision.  
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
 
      Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      //if (charge == 0.) continue;
      Float_t eta = part->Eta();
      Float_t pt  = part->Pt();
   //   if (pt < 0.2) continue;
      parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
      parts.back().set_user_info(new MyInfo(charge));
      etaH->Fill(eta);
      if (pt > 0.) ptH->Fill(pt, 1./(2. * pt));
    }

    ClusterSequence cs(parts, jet_def); // building the jets in the events
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // Store all the jets in the event
    if(jets.size()<2 || jets.size()>5) continue; //Throw away events with less than two jets
    cout << " num jets : " << jets.size() << endl;
    int NumJets = 0; // Not important for now
    // Jet loop
    int pdg_jet = -99; // variable to store which flavor the jet came from (up jet? down jet? etc)
    bool has_parton1 = false; // boolean to tell whether the jet contained the first outgoing quark/gluon
    bool has_parton2 = false;// boolean to tell whether the jet contained the second outgoing quark/gluon
//cout<< "jets size" << jets.size() <<endl;
    for (unsigned i = 0; i < jets.size(); i++) //Loop over the jets!
    {
      if (jets[i].pt() < 20) continue; // get rid of low pT jets
	
      //cout << "jet num : " <<  i << endl;
      PseudoJet jet = jets[i];
      TLorentzVector jetvec(jet.px(), jet.py(), jet.pz(), jet.e()); // Four vector of the jet
      if(checkdR(jetvec, parton1) < 0.1) has_parton1 = true;  //Check whether the direction of the jet matches the first outgoing parton
      else if(checkdR(jetvec, parton2) < 0.1) has_parton2 = true; //Check whether the direction of the jet matches the second outgoing parton


      if(has_parton1&&(!has_parton2)){// If this is true, then the jet came from the first outgoing parton
        pdg_jet = pdg_parton1;

      }
      else if(has_parton2 && (!has_parton1)){// If this is true, then the jet came from the second outgoing parton
        pdg_jet = pdg_parton2;
      }

      vector<PseudoJet> constituents = jet.constituents(); //give me jet daughters
     if (constituents.size()<2) continue;
     
      double jetcharge = 0;
    //  double kappa = 0.3;

      for (unsigned j = 0; j < constituents.size(); j++) // knowing what each jet daugther property for us to maybe build jet charge
      {
 PseudoJet con = constituents[j];
	     
    
	      double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();

  
        TVector3 con3(con.px(), con.py(), con.pz());
     //   double  dtrid = constituents.at(j).user_info<MyInfo>().pdg_id();
  	jetcharge+=pow(con.pt(), kappa)*dtrcharge/3.; //compute jet charge

      }

      jetcharge/=pow(jets[i].pt(), kappa); //normalize by jet pT
      //cout<<jetcharge<<",";
      //Fill histograms
 
if ( jetcharge == 0){
 
 	
        cout<< "particles of jets "<<  constituents.size() << endl;
 for (unsigned j = 0; j < constituents.size(); j++){
 PseudoJet con = constituents[j];
	 //	 int  dtrid = constituents.at(j).user_info<MyInfo>().pdg_id();
//	 cout << "pdg of  constituents " << dtrid << endl; 
//         int  dtrid = constituents.at(j).user_info<MyInfo>().pdg_id();
//  double  dtrcharge = constituents.at(j).user_info<MyInfo>().pdg_charge();
//  	 cout << " charge of  constituents " <<  dtrcharge << endl;
      for (Int_t ip = 0; ip < np; ip++) {
         TParticle* part = (TParticle*) particles->At(ip); // creating the T particle
         Int_t ist = part->GetStatusCode();
         if (ist <= 0) continue;

         // check the 4vector of cons against 4vec of pythia particles to find a possible match
         if (con.px()==part->Px()&& con.py()==part->Py() && con.pz()==part->Pz()&& con.e()==part->Energy()){
                 //this condition will tell us the constituent we are looking at will match the T particle. Bc we can use Tparticle to get the pdg code for dtrs
        //       cout<< "this is working T~T"<< endl; }
         cout << "pdg for dtrs " << part->GetPdgCode() << endl;
         int cPid;
         cPid=part->GetPdgCode();
        cout << "charge of dtrs " << TDatabasePDG::Instance()->GetParticle(cPid)->Charge()<< endl;
      }

        }
cout<< endl; 

}
}

  
 

      if(pdg_jet == 1) h1_d_jetcharge->Fill(jetcharge);//pdg code for down quark
      else if(pdg_jet == -1) h1_dbar_jetcharge->Fill(jetcharge);// pdg code for downbar quark
      else if(pdg_jet == 2) h1_u_jetcharge->Fill(jetcharge);// pdg code for u
      else if(pdg_jet == -2) h1_ubar_jetcharge->Fill(jetcharge);//pdg code for ubar
      else if(pdg_jet == 21) h1_g_jetcharge->Fill(jetcharge);// pdg code for gluon
      else if(pdg_jet == 3) h1_s_jetcharge->Fill(jetcharge);//pdg code for s
      else if(pdg_jet == -3) h1_sbar_jetcharge->Fill(jetcharge);//pdg code for sbar
      else if(pdg_jet == 4) h1_c_jetcharge->Fill(jetcharge);//pdg code for c
      else if(pdg_jet == -4) h1_cbar_jetcharge->Fill(jetcharge);// pdg code for cbar
      else if(pdg_jet == 5) h1_b_jetcharge->Fill(jetcharge);//pdg code for b
      else if(pdg_jet == -5) h1_bbar_jetcharge->Fill(jetcharge);//pdg code for bbar

//cout<< "pdg code" << pdg_jet << endl;
   
      //else if(pdg_jet == -2) h1_ubar_jetcharge->Fill(jetcharge);
      //h->Fill(jetcharge);

    }

  
  }


  //if normalazation is on, you will see error bars thus you will not get the boxy histogram boys
  //say no to boxy boys
  //Normalize histograms such that their integral = 1
  h1_u_jetcharge->Scale(1./h1_u_jetcharge->Integral());
  h1_ubar_jetcharge->Scale(1./h1_ubar_jetcharge->Integral());
  h1_d_jetcharge->Scale(1./h1_d_jetcharge->Integral());
  h1_dbar_jetcharge->Scale(1./h1_dbar_jetcharge->Integral());
  h1_g_jetcharge->Scale(1./h1_g_jetcharge->Integral());
  h1_s_jetcharge->Scale(1./h1_s_jetcharge->Integral());
  h1_sbar_jetcharge->Scale(1./h1_sbar_jetcharge->Integral());
  h1_c_jetcharge->Scale(1./h1_c_jetcharge->Integral());
  h1_cbar_jetcharge->Scale(1./h1_cbar_jetcharge->Integral());
  h1_b_jetcharge->Scale(1./h1_b_jetcharge->Integral());
  h1_bbar_jetcharge->Scale(1./h1_bbar_jetcharge->Integral());




  pythia8->PrintStatistics();

// The following makes pdf plots

      //---- paint setup...
    //
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
    gStyle->SetOptStat(1);
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




    //c1->SaveAs(Form("plots_misc/Misc_%s"+extension+".pdf", file_name.c_str()));
    //cout<<"...outbase   = "<<outbase.Data()<<endl;
    //cout<<"...rootfile  = "<<rootfile.Data()<<endl;
    //cout<<"...plotfile  = "<<plotfile.Data()<<endl;


    //
    // Start new page!!!!
  /*  //
    ++ican;
    sprintf(buf,"ccan%d",ican);
    ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
    ccan[ican]->SetFillColor(10);
    //gPad->SetLeftMargin(0.16);
    //gPad->SetBottomMargin(0.06);
    ccan[ican]->cd(); ccan[ican]->Divide(1,1,0.0001,0.0001);
    ccan[ican]->cd(1);
    */
    // mj radomly writing stuff
    //
  // new page
++ican;
sprintf(buf,"ccan%d",ican);
ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
ccan[ican]->SetFillColor(10);
//gPad->SetLeftMargin(0.16);
//gPad->SetBottomMargin(0.06);
ccan[ican]->cd(); ccan[ican]->Divide(2,2,0.0001,0.0001);
  ccan[ican]->cd(1);
    h1_c_jetcharge->SetXTitle("Jet charge for c");
    h1_c_jetcharge->SetYTitle("Counts");
    h1_c_jetcharge->SetLineColor(kBlack);
    h1_c_jetcharge->SetLineStyle(kSolid);
    h1_c_jetcharge->Draw();
  ccan[ican]->cd(2);
    h1_cbar_jetcharge->Draw();
    h1_cbar_jetcharge->SetXTitle("Jet charge for cbar");
    h1_cbar_jetcharge->SetYTitle("Counts");
    h1_cbar_jetcharge->SetLineColor(kBlack);
    h1_cbar_jetcharge->SetLineStyle(kSolid);
  ccan[ican]->cd(3);
    h1_b_jetcharge->Draw();
    h1_b_jetcharge->SetXTitle("Jet charge for b");
    h1_b_jetcharge->SetYTitle("Counts");
    h1_b_jetcharge->SetLineColor(kBlack);
    h1_b_jetcharge->SetLineStyle(kSolid);
   ccan[ican]->cd(4);
    h1_bbar_jetcharge->Draw();
    h1_bbar_jetcharge->SetXTitle("Jet charge for bbar");
    h1_bbar_jetcharge->SetYTitle("Counts");
    h1_bbar_jetcharge->SetLineColor(kBlack);
    h1_bbar_jetcharge->SetLineStyle(kSolid);
                                          
      ccan[ican]->cd();ccan[ican]->Update();
if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
     else { ccan[ican]->Print(plotfilePDF.Data()); }

 
    // new page
++ican;
sprintf(buf,"ccan%d",ican);
ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
ccan[ican]->SetFillColor(10);
//gPad->SetLeftMargin(0.16);
//gPad->SetBottomMargin(0.06);
ccan[ican]->cd(); ccan[ican]->Divide(2,2/*,0.0001,0.0001*/);
  ccan[ican]->cd(1);
  cout << " where did this go " << endl;
    h1_u_jetcharge->SetXTitle("Jet charge for u");
    h1_u_jetcharge->SetYTitle("Counts");
    h1_u_jetcharge->SetLineColor(kBlack);
    h1_u_jetcharge->SetLineStyle(kSolid);
    h1_u_jetcharge->Draw();
  ccan[ican]->cd(2);
    h1_ubar_jetcharge->Draw();
    h1_ubar_jetcharge->SetXTitle("Jet charge for ubar");
    h1_ubar_jetcharge->SetYTitle("Counts");
    h1_ubar_jetcharge->SetLineColor(kBlack);
    h1_ubar_jetcharge->SetLineStyle(kSolid);
  ccan[ican]->cd(3);
    h1_d_jetcharge->Draw();
    h1_d_jetcharge->SetXTitle("Jet charge for d");
    h1_d_jetcharge->SetYTitle("Counts");
    h1_d_jetcharge->SetLineColor(kBlack);
    h1_d_jetcharge->SetLineStyle(kSolid);
  ccan[ican]->cd(4);
    h1_dbar_jetcharge->Draw();
    h1_dbar_jetcharge->SetXTitle("Jet charge for dbar");
    h1_dbar_jetcharge->SetYTitle("Counts");
    h1_dbar_jetcharge->SetLineColor(kBlack); 
  h1_dbar_jetcharge->SetLineStyle(kSolid);


  /*auto legend1 = new TLegend(0.6,0.7,0.8,0.9);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetFillColor(3);
    legend1->AddEntry(kappa, "kappa");
    legend1->AddEntry(h1_ubar_jetcharge, "ubar");
    legend1->AddEntry(h1_d_jetcharge, "d");
    legend1->AddEntry(h1_dbar_jetcharge, "dbar");
     legend1->AddEntry(h1_g_jetcharge, "g");
    legend1->AddEntry(h1_s_jetcharge, "s");
   legend1->AddEntry(h1_sbar_jetcharge, "sbar");
  // legend1->AddEntry(nev, "events");
    legend1->Draw("same");
*/
/*
TLatex tt;
  
char cc[256] = "kappa";
kappa = 0.3;
printf(cc, kappa);
tt.DrawLatex(0.5,0.1,cc,kappa);
*/

     ccan[ican]->cd();ccan[ican]->Update();
if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
else { ccan[ican]->Print(plotfilePDF.Data()); }

    /*
  f.Write();
  f.Close();
    if (ican>-1){
    cout<<" You plotted "<<ican+1<<" canvasses......."<<endl;
    ccan[ican]->Print(plotfileC.Data());
    }
*/
/*
//newpage
++ican;
sprintf(buf,"ccan%d",ican);
ccan[ican] = new TCanvas(buf,buf,30*ican,30*ican,800,(8.5/11.)*800);
ccan[ican]->SetFillColor(10);
//gPad->SetLeftMargin(0.16);
//gPad->SetBottomMargin(0.06);
ccan[ican]->cd(); ccan[ican]->Divide(2,2,0.0001,0.0001);
  ccan[ican]->cd(1);
    h1_g_jetcharge->SetStats(0);
    h1_g_jetcharge->Draw();
    h1_g_jetcharge->SetXTitle("Jet charge for g");
    h1_g_jetcharge->SetYTitle("Counts");
    h1_g_jetcharge->SetLineColor(kBlack);
    h1_g_jetcharge->SetLineStyle(kSolid);
  ccan[ican]->cd(2);
    h1_s_jetcharge->SetStats(0);
    h1_s_jetcharge->Draw();
    h1_s_jetcharge->SetXTitle("Jet charge for s");
    h1_s_jetcharge->SetYTitle("Counts");
    h1_s_jetcharge->SetLineColor(kBlack);
    h1_s_jetcharge->SetLineStyle(kSolid);
  ccan[ican]->cd(3);
    h1_sbar_jetcharge->SetStats(0);
    h1_sbar_jetcharge->Draw();
    h1_sbar_jetcharge->SetXTitle("Jet charge for sbar");
    h1_sbar_jetcharge->SetYTitle("Counts");
    h1_sbar_jetcharge->SetLineColor(kBlack);
    h1_sbar_jetcharge->SetLineStyle(kSolid);
*/
  

     
 /*
    //plot stuff here!!!!
    h1_u_jetcharge->SetXTitle("Jet Charge");
    h1_u_jetcharge->SetLineColor(kBlack);
    h1_u_jetcharge->SetLineStyle(kSolid);
    h1_ubar_jetcharge->SetLineColor(kBlack);
    h1_ubar_jetcharge->SetLineStyle(kDashed);
    h1_d_jetcharge->SetLineColor(kRed);
    h1_d_jetcharge->SetLineStyle(kSolid);
    h1_dbar_jetcharge->SetLineColor(kRed);
    h1_dbar_jetcharge->SetLineStyle(kDashed);
    h1_g_jetcharge->SetLineColor(kGreen);
    h1_g_jetcharge->SetLineStyle(kSolid);
    h1_s_jetcharge->SetLineColor(kMagenta);
    h1_s_jetcharge->SetLineStyle(kSolid);
    h1_sbar_jetcharge->SetLineColor(kMagenta);
    h1_sbar_jetcharge->SetLineStyle(kDashed);
    h1_u_jetcharge->Draw("HIST same");
    h1_ubar_jetcharge->Draw("HIST same");
    h1_d_jetcharge->Draw("HIST same");
    h1_dbar_jetcharge->Draw("HIST same");
    h1_g_jetcharge->Draw("HIST same");
    h1_s_jetcharge->Draw("HIST same");
    h1_sbar_jetcharge->Draw("HIST same");
*/
 /*   // add a legend!!
    auto legend1 = new TLegend(0.6,0.7,0.8,0.9);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetFillColor(3);
    legend1->AddEntry(h1_u_jetcharge, "u");
    legend1->AddEntry(h1_ubar_jetcharge, "ubar");
    legend1->AddEntry(h1_d_jetcharge, "d");
    legend1->AddEntry(h1_dbar_jetcharge, "dbar");
     legend1->AddEntry(h1_g_jetcharge, "g");
    legend1->AddEntry(h1_s_jetcharge, "s");
   legend1->AddEntry(h1_sbar_jetcharge, "sbar");
  // legend1->AddEntry(nev, "events");
    legend1->Draw("same");
*/    
  


/*
    ccan[ican]->cd();ccan[ican]->Update();
    if (ican==0){ ccan[ican]->Print(plotfileO.Data()); }
       else { ccan[ican]->Print(plotfilePDF.Data()); }
*/       
    //
//    f.Write();
 //   f.Close();

/*   
TLatex tt;
  
char cc[256] = "";
sprintf(cc,"kappa",kappa);
tt.DrawLatex(0.5,0.1,cc);
*/

    if (ican>-1){
    cout<<" You plotted "<<ican+1<<" canvasses......."<<endl;
    ccan[ican]->Print(plotfileC.Data());
    }




}




