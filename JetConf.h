#include <vector>
#include <string>
#include <iostream>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TRandom3.h>
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/SoftDrop.hh"

//using namespace fastjet;

//double dot_pseudo(PseudoJet pj1, PseudoJet pj2);
//void decluster_jet(PseudoJet mom, PseudoJet hfhad, vector<double> &kTs, vector<double> &Erads,
 //                   vector<double> &thetas, vector<double> &deltars);

double checkphi(double phi);
double folddphi(double dphi);
double checkdR(TLorentzVector v1, TLorentzVector v2);

bool pt_cut = true;
bool Erad_cut = true;
bool kt_cut = true;
bool pt_track_cut = false;
bool removeHFDTRs = true;
bool followHardest = true;
bool fullReconstruction = false;
bool truthLevel = false;
bool ghostCut = true;
bool bdtCut = true;
bool softDrop = false;
bool chargedJetCut = false;
bool PrimaryVertexCut = true;

const float zjetdphicut = 7. * TMath::Pi() / 8.;
double B_mass = 5.2;
double mass_num;
double BDTbc_cut = 0.1;
double BDTbcudsg_cut = 0.0;
double ghostProb = 0.5;
double jetradius = 0.5;
//double Erad_low = 110.; //GeV
//double Erad_high = 120.;
double etaMin = 2.5;
double etaMax = 4.0;
double LambdaQCD = 0.2;
double ptMin = 20;
double ptMax = 150;
double sd_z_cut = 0.10;
double sd_beta  = 0.0;

TF1 *f1_theta_Erad = new TF1("f1_theta_Erad", "log(x/[0])", 1, 800);
TF1 *f1_kt_Erad = new TF1("f1_kt_Erad", "log(x/[0])", 1, 800);


const int NumHists = 5;
const int NumProj = 16;
int MinErad = 100;
int MaxErad = 200;
int EradStep = 20;
int MinEradProj = 20;
int MaxEradProj = MinEradProj+EradStep*NumProj;
int StepErad = (MaxErad - MinErad)/NumHists;
int h2_erad_bins = 100;
int h2_erad_low = 0;
int h2_erad_high = 400;
int R_bins = 4;
const int NumRbins = 5;
int RbinsArr[NumRbins] = {5, 6, 7, 8, 9};
double R_min = 2.2;
double R_max = 6.5;

const int pTbinsize = 4;
float pTres_binedges[pTbinsize+1] = {10, 20, 30, 70, 150};

/*
double dot_pseudo(PseudoJet pj1, PseudoJet pj2)
{
    return pj1.px()*pj2.px()+pj1.py()*pj2.py()+pj1.pz()*pj2.pz();
}

void decluster_jet(PseudoJet mom, PseudoJet hfhad, vector<double> &kTs, vector<double> &Erads,
                    vector<double> &thetas, vector<double> &deltars){
  PseudoJet dtr1, dtr2;
  double Erad, Esoft, pt_leading,kT, theta, deltar;
  bool check_splitting;

  while(mom.has_parents(dtr1, dtr2))
  {
      // if(followHardest && removeHFDTRs && (flavor == 4 || flavor == 5)){
      //   if(mom.contains(hfhad)) hardest_is_hf++;
      //   else hardest_not_hf++;
      // }
      //if(!mom.contains(hfhad)) break;

      //cout<<mom.contains(gdtr1)<<", ";
      //cout<<"dtr loop";
      double mag1 = dtr1.modp();
      double mag2 = dtr2.modp();
      Erad = mom.e();
      //cout<<"mag1 = "<< mag1 <<", mag2 = "<< mag2<<endl;
      //cout<<dot_pseudo(dtr1, dtr2)<<endl;
      //TLorentzVector p1(dtr1.px(), dtr1.py(), dtr1.pz(), dtr1.e());
      //TLorentzVector p2(dtr2.px(), dtr2.py(), dtr2.pz(), dtr2.e());
      //double theta_tlv = p1.Angle(p2.Vect());
      //cout<<mag1<<","<<mag2<<",";
      double theta = acos(dot_pseudo(dtr1, dtr2)/(mag1*mag2));
      double deltar = sqrt(pow(dtr1.rap()-dtr2.rap(), 2)+pow(dtr1.phi()-dtr2.phi(), 2));
      //cout<<theta - theta_tlv<<endl;
      //printf("theta = %f \n", theta);
      if(followHardest) check_splitting = (dtr1.e() > dtr2.e());
      else check_splitting = (dtr1.contains(hfhad)|| abs(dtr1.e() - hfhad.e()) < 1e-4);
      if(check_splitting)
      {
        Esoft = dtr2.e();
        //Erad = dtr1.e();
        pt_leading = dtr1.pt();
        kT = Esoft*theta;
        mom.reset(dtr1);

      }
      else
      {
        Esoft = dtr1.e();
        // Erad = dtr2.e();
        pt_leading = dtr2.pt();
        kT = Esoft*theta;
        mom.reset(dtr2);
      }

      if (kt_cut)
      {
        //if (kT < LambdaQCD) cout<<Erad<<", ";
        if (kT < LambdaQCD) continue; //Non-perturbative cut
      }
      kTs.push_back(kT);
      Erads.push_back(Erad);
      thetas.push_back(theta);
      deltars.push_back(deltar);

      //h1_erad->Fill(Erad);
      //h2_thetaErad->Fill(Erad, log(1./theta));
      //h1_theta->Fill(theta);

      // Begin prong cuts


      // if (pt_track_cut)
      // {
      //   //cout<<pt_leading<<",";
      //   if (pt_leading < pt_track) continue;
      // }
      // End prong cuts
      //theta_counter++;



    }
}
*/
double checkphi(double phi)
{
  float returnphi = phi;

  if(phi > 1. * TMath::Pi())
    returnphi -= 2 * TMath::Pi();
  else if(phi < -1 * TMath::Pi() )
    returnphi += 2 * TMath::Pi();

  return returnphi;
}


double folddphi(double dphi)
{
  double returndphi = dphi;
  if(dphi < 0)
    returndphi += 2. * TMath::Pi();
  if(dphi > 2. * TMath::Pi())
    returndphi -= 2. * TMath::Pi();
  if( dphi > TMath::Pi() )
    returndphi = TMath::Pi() - (dphi - TMath::Pi());

  return returndphi;

}
double checkdR(TLorentzVector v1, TLorentzVector v2)
{

  double dr = -99;

  double dphi = checkphi(v1.Phi()) - checkphi(v2.Phi());
  double deta = v1.Eta() - v2.Eta();
  dphi = checkphi(dphi);

  dr = sqrt(dphi * dphi + deta * deta);

  return dr;



}
void BinLogX(TH1*h)
{

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = TMath::Log10(axis->GetXmin());
   Axis_t to = TMath::Log10(axis->GetXmax());
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
}
