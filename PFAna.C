#define PFAna_cxx
#include "PFAna.h"
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

#include <map>

const double deltaRmax = 0.4;//0.2;//0.4;

void PFAna::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PFAna.C
//      root> PFAna t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  // disable all branches

   //PFAnaTree->Draw("PFCand_pt/genPart_pt","PFCand_dRToGenPart<0.16 && nPFCand==1 && PFCand_mass>0 && nSimHitEE<60 && !PFCand_has_trk","same")
   
   fChain->SetBranchStatus("genPart_pt",1);
   fChain->SetBranchStatus("genPart_eta",1);
   fChain->SetBranchStatus("genPart_phi",1);
   fChain->SetBranchStatus("genPart_energy",1);
   fChain->SetBranchStatus("genPart_extrapolated_HCAL_ieta",1);
   fChain->SetBranchStatus("genPart_extrapolated_HCAL_iphi",1);

   fChain->SetBranchStatus("nSimHitEE",1);
   fChain->SetBranchStatus("SimHitEE_eta",1);
   fChain->SetBranchStatus("SimHitEE_phi",1);
   fChain->SetBranchStatus("SimHitEE_energy",1);

   fChain->SetBranchStatus("nSimHitHBHE",1);
   fChain->SetBranchStatus("SimHitHBHE_eta",1);
   fChain->SetBranchStatus("SimHitHBHE_phi",1);
   fChain->SetBranchStatus("SimHitHBHE_depth",1);
   fChain->SetBranchStatus("SimHitHBHE_energy",1);
   fChain->SetBranchStatus("SimHitHBHE_samplingFactor",1);
   fChain->SetBranchStatus("SimHitHBHE_ieta",1);
   fChain->SetBranchStatus("SimHitHBHE_iphi",1);

   fChain->SetBranchStatus("nPFRecHitEE",1);
   fChain->SetBranchStatus("PFRecHitEE_eta",1);
   fChain->SetBranchStatus("PFRecHitEE_phi",1);
   fChain->SetBranchStatus("PFRecHitEE_energy",1);

   fChain->SetBranchStatus("nPFRecHitHBHE",1);
   fChain->SetBranchStatus("PFRecHitHBHE_eta",1);
   fChain->SetBranchStatus("PFRecHitHBHE_phi",1);
   fChain->SetBranchStatus("PFRecHitHBHE_depth",1);
   fChain->SetBranchStatus("PFRecHitHBHE_energy",1);
   fChain->SetBranchStatus("PFRecHitHBHE_ieta",1);
   fChain->SetBranchStatus("PFRecHitHBHE_iphi",1);
   
   fChain->SetBranchStatus("nPFClusterECAL",1);
   fChain->SetBranchStatus("PFClusterECAL_eta",1);
   fChain->SetBranchStatus("PFClusterECAL_phi",1);
   fChain->SetBranchStatus("PFClusterECAL_energy",1);

   fChain->SetBranchStatus("nPFClusterHCAL",1);
   fChain->SetBranchStatus("PFClusterHCAL_eta",1);
   fChain->SetBranchStatus("PFClusterHCAL_phi",1);
   fChain->SetBranchStatus("PFClusterHCAL_layer",1);
   fChain->SetBranchStatus("PFClusterHCAL_energy",1);

   TLorentzVector p4g, p4s, p4rc;

   // Create a ROOT file to store the TTree
   TFile *file = new TFile("PFAnaTree.root", "RECREATE");

   // Create a TTree named "s"
   TTree *tree = new TTree("s", "PFAna TTree");

   // Variables for branches
   //Float_t genPart_pt, genPart_eta, genPart_phi, genPart_energy;
   Int_t genPart_ieta;
   Int_t genPart_iphi;

   const int nDepth = 8;
   Float_t genPart_fe;
   Float_t genPart_sum;
   Float_t genPart_depth[nDepth];

   Float_t recHit_fe;
   Float_t recHit_sum;
   Float_t recHit_depth[nDepth];

   Float_t randomCone_fe;
   Float_t randomCone_sum;
   Float_t randomCone_depth[nDepth];
   
   Float_t pfCluster_fe;
   Float_t pfCluster_sum;
   Float_t pfCluster_depth[nDepth];

   Float_t randomCone2_fe;
   Float_t randomCone2_sum;
   Float_t randomCone2_depth[nDepth];

   // Create branches
   tree->Branch("genPart_pt", &genPart_pt, "genPart_pt/F");
   tree->Branch("genPart_eta", &genPart_eta, "genPart_eta/F");
   tree->Branch("genPart_phi", &genPart_phi, "genPart_phi/F");
   tree->Branch("genPart_energy", &genPart_energy, "genPart_energy/F");
   tree->Branch("genPart_ieta", &genPart_ieta, "genPart_ieta/I");
   tree->Branch("genPart_iphi", &genPart_iphi, "genPart_iphi/I");

   tree->Branch("genPart_fe", &genPart_fe, "genPart_fe/F");
   tree->Branch("genPart_sum", &genPart_sum, "genPart_sum/F");
   tree->Branch("genPart_depth", genPart_depth, "genPart_depth[8]/F");

   tree->Branch("recHit_fe", &recHit_fe, "recHit_fe/F");
   tree->Branch("recHit_sum", &recHit_sum, "recHit_sum/F");
   tree->Branch("recHit_depth", recHit_depth, "recHit_depth[8]/F");

   tree->Branch("randomCone_fe", &randomCone_fe, "randomCone_fe/F");
   tree->Branch("randomCone_sum", &randomCone_sum, "randomCone_sum/F");
   tree->Branch("randomCone_depth", randomCone_depth, "randomCone_depth[8]/F");

   tree->Branch("pfCluster_fe", &pfCluster_fe, "pfCluster_fe/F");
   tree->Branch("pfCluster_sum", &pfCluster_sum, "pfCluster_sum/F");
   tree->Branch("pfCluster_depth", pfCluster_depth, "pfCluster_depth[8]/F");

   tree->Branch("randomCone2_fe", &randomCone2_fe, "randomCone2_fe/F");
   tree->Branch("randomCone2_sum", &randomCone2_sum, "randomCone2_sum/F");
   tree->Branch("randomCone2_depth", randomCone2_depth, "randomCone2_depth[8]/F");

   TFile *fout = new TFile("PFAna.root","RECREATE");

   TH1D *hfe = new TH1D("hfe","",100,0,1);

   TH1D *hnrechit_hh = new TH1D("hnrechit_hh","",140,0,140);
   TH1D *hnrechite_hh = new TH1D("hnrechite_hh","",140,0,140);
   TH1D *hnrechith_hh = new TH1D("hnrechith_hh","",140,0,140);
   TH1D *hnrechitrc_hh = new TH1D("hnrechitrc_hh","",140,0,140);
   TH1D *hnrechitrce_hh = new TH1D("hnrechitrce_hh","",140,0,140);
   TH1D *hnrechitrch_hh = new TH1D("hnrechitrch_hh","",140,0,140);

   TH1D *hntwr_hh = new TH1D("hntwr_hh","",20,0,20);
   TH1D *hntwrrc_hh = new TH1D("hntwrrc_hh","",20,0,20);
      
   TH1D *hnpf_hh = new TH1D("hnpf_hh","",20,0,20);
   TH1D *hnpfe_hh = new TH1D("hnpfe_hh","",20,0,20);
   TH1D *hnpfh_hh = new TH1D("hnpfh_hh","",20,0,20);
   TH1D *hnpfrc_hh = new TH1D("hnpfrc_hh","",20,0,20);
   TH1D *hnpfrce_hh = new TH1D("hnpfrce_hh","",20,0,20);
   TH1D *hnpfrch_hh = new TH1D("hnpfrch_hh","",20,0,20);

   TH1D *hgensum = new TH1D("hgensum","",400,-100,300);
   TH1D *hgensum_hh = new TH1D("hgensum_hh","",400,-100,300);
   TH1D *hgensum_eh = new TH1D("hgensum_eh","",400,-100,300);
   TH1D *hgensum_ee = new TH1D("hgensum_ee","",400,-100,300);

   TH1D *hrecsum = new TH1D("hrecsum","",400,-100,300);
   TH1D *hrecsum_hh = new TH1D("hrecsum_hh","",400,-100,300);
   TH1D *hrecsum_eh = new TH1D("hrecsum_eh","",400,-100,300);
   TH1D *hrecsum_ee = new TH1D("hrecsum_ee","",400,-100,300);

   TH1D *hpfsum = new TH1D("hpfsum","",400,-100,300);
   TH1D *hpfsum_hh = new TH1D("hpfsum_hh","",400,-100,300);
   TH1D *hpfsum_eh = new TH1D("hpfsum_eh","",400,-100,300);
   TH1D *hpfsum_ee = new TH1D("hpfsum_ee","",400,-100,300);

   TH1D *hrcsum = new TH1D("hrcsum","",400,-100,300);
   TH1D *hrcsum_hh = new TH1D("hrcsum_hh","",400,-100,300);
   TH1D *hrcsum_eh = new TH1D("hrcsum_eh","",400,-100,300);
   TH1D *hrcsum_ee = new TH1D("hrcsum_ee","",400,-100,300);

   // RC-subtracted sums
   TH1D *hrecrcsum = new TH1D("hrecrcsum","",400,-100,300);
   TH1D *hrecrcsum_hh = new TH1D("hrecrcsum_hh","",400,-100,300);
   TH1D *hrecrcsum_eh = new TH1D("hrecrcsum_eh","",400,-100,300);
   TH1D *hrecrcsum_ee = new TH1D("hrecrcsum_ee","",400,-100,300);

   TH1D *hpfrcsum = new TH1D("hpfrcsum","",400,-100,300);
   TH1D *hpfrcsum_hh = new TH1D("hpfrcsum_hh","",400,-100,300);
   TH1D *hpfrcsum_eh = new TH1D("hpfrcsum_eh","",400,-100,300);
   TH1D *hpfrcsum_ee = new TH1D("hpfrcsum_ee","",400,-100,300);

   // Tower energies
   TH1D *htwr = new TH1D("htwr","",400,-100,300);
   TH1D *htwr_hh = new TH1D("htwr_hh","",400,-100,300);
   TH1D *htwr_eh = new TH1D("htwr_eh","",400,-100,300);
   TH1D *htwr_ee = new TH1D("htwr_ee","",400,-100,300);

   TH1D *htwrrc = new TH1D("htwrrc","",400,-100,300);
   TH1D *htwrrc_hh = new TH1D("htwrrc_hh","",400,-100,300);
   TH1D *htwrrc_eh = new TH1D("htwrrc_eh","",400,-100,300);
   TH1D *htwrrc_ee = new TH1D("htwrrc_ee","",400,-100,300);

   // Leading Towers
   TH1D *htwr1_hh = new TH1D("htwr1_hh","",400,-100,300);
   TH1D *htwr2_hh = new TH1D("htwr2_hh","",400,-100,300);
   TH1D *htwr3_hh = new TH1D("htwr3_hh","",400,-100,300);
   TH1D *htwr4_hh = new TH1D("htwr4_hh","",400,-100,300);
   TH1D *htwr12_hh = new TH1D("htwr12_hh","",400,-100,300);
   TH1D *htwr123_hh = new TH1D("htwr123_hh","",400,-100,300);

   TH1D *htwr1rc_hh = new TH1D("htwr1rc_hh","",400,-100,300);
   TH1D *htwr2rc_hh = new TH1D("htwr2rc_hh","",400,-100,300);
   TH1D *htwr3rc_hh = new TH1D("htwr3rc_hh","",400,-100,300);
   TH1D *htwr4rc_hh = new TH1D("htwr4rc_hh","",400,-100,300);
   TH1D *htwr12rc_hh = new TH1D("htwr12rc_hh","",400,-100,300);
   TH1D *htwr123rc_hh = new TH1D("htwr123rc_hh","",400,-100,300);
   
   // Leading PFClusters
   TH1D *hpf1_hh = new TH1D("hpf1_hh","",400,-100,300);
   TH1D *hpf2_hh = new TH1D("hpf2_hh","",400,-100,300);
   TH1D *hpf3_hh = new TH1D("hpf3_hh","",400,-100,300);
   TH1D *hpf4_hh = new TH1D("hpf4_hh","",400,-100,300);
   TH1D *hpf12_hh = new TH1D("hpf12_hh","",400,-100,300);
   TH1D *hpf123_hh = new TH1D("hpf123_hh","",400,-100,300);

   TH1D *hpf1rc_hh = new TH1D("hpf1rc_hh","",400,-100,300);
   TH1D *hpf2rc_hh = new TH1D("hpf2rc_hh","",400,-100,300);
   TH1D *hpf3rc_hh = new TH1D("hpf3rc_hh","",400,-100,300);
   TH1D *hpf4rc_hh = new TH1D("hpf4rc_hh","",400,-100,300);
   TH1D *hpf12rc_hh = new TH1D("hpf12rc_hh","",400,-100,300);
   TH1D *hpf123rc_hh = new TH1D("hpf123rc_hh","",400,-100,300);
   
   /*
   TH2D *h2gendepth = new TH2D("h2gendepth","",8,0,8,100,0,100);
   TH2D *h2recdepth = new TH2D("h2recdepth","",8,0,8,100,0,100);
   TH2D *h2rcdepth = new TH2D("h2rcdept","",8,0,8,100,0,100);
   
   TH2D *h2gendepth_hh = new TH2D("h2gendepth_hh","",8,0,8,100,0,100);
   TH2D *h2recdepth_hh = new TH2D("h2recdepth_hh","",8,0,8,100,0,100);
   TH2D *h2rcdepth_hh = new TH2D("h2rcdepth_hh","",8,0,8,100,0,100);
   
   TH2D *h2gendepth_eh = new TH2D("h2gendepth_eh","",8,0,8,100,0,100);
   TH2D *h2recdepth_eh = new TH2D("h2recdepth_eh","",8,0,8,100,0,100);
   TH2D *h2rcdepth_eh = new TH2D("h2rcdepth_eh","",8,0,8,100,0,100);
   
   TH2D *h2gendepth_ee = new TH2D("h2gendepth_ee","",8,0,8,100,0,100);
   TH2D *h2recdepth_ee = new TH2D("h2recdepth_ee","",8,0,8,100,0,100);
   TH2D *h2rcdepth_ee = new TH2D("h2rcdepth_ee","",8,0,8,100,0,100);
   */

   // Studies of pT cut per depth
   double gensum(0), gensum_hh(0), gensum_eh(0), gensum_ee(0);
   double recsum(0), recsum_hh(0), recsum_eh(0), recsum_ee(0);
   double rcsum(0), rcsum_hh(0), rcsum_eh(0), rcsum_ee(0);
   double recrcsum(0), recrcsum_hh(0), recrcsum_eh(0), recrcsum_ee(0);
   TH2D *h2gendepth = new TH2D("h2gendepth","",8,0,8,50,0,5);
   TH2D *h2recdepth = new TH2D("h2recdepth","",8,0,8,50,0,5);
   TH2D *h2rcdepth = new TH2D("h2rcdepth","",8,0,8,50,0,5);
   TH2D *h2recrcdepth = new TH2D("h2recrcdepth","",8,0,8,50,0,5);
   
   TH2D *h2gendepth_hh = new TH2D("h2gendepth_hh","",8,0,8,50,0,5);
   TH2D *h2recdepth_hh = new TH2D("h2recdepth_hh","",8,0,8,50,0,5);
   TH2D *h2rcdepth_hh = new TH2D("h2rcdepth_hh","",8,0,8,50,0,5);
   TH2D *h2recrcdepth_hh = new TH2D("h2recrcdepth_hh","",8,0,8,50,0,5);
   
   TH2D *h2gendepth_eh = new TH2D("h2gendepth_eh","",8,0,8,50,0,5);
   TH2D *h2recdepth_eh = new TH2D("h2recdepth_eh","",8,0,8,50,0,5);
   TH2D *h2rcdepth_eh = new TH2D("h2rcdepth_eh","",8,0,8,50,0,5);
   TH2D *h2recrcdepth_eh = new TH2D("h2recrcdepth_eh","",8,0,8,50,0,5);
   
   TH2D *h2gendepth_ee = new TH2D("h2gendepth_ee","",8,0,8,50,0,5);
   TH2D *h2recdepth_ee = new TH2D("h2recdepth_ee","",8,0,8,50,0,5);
   TH2D *h2rcdepth_ee = new TH2D("h2rcdepth_ee","",8,0,8,50,0,5);
   TH2D *h2recrcdepth_ee = new TH2D("h2recrcdepth_ee","",8,0,8,50,0,5);

   // Energy depositions per depth
   TProfile *pgendepth = new TProfile("pgendepth","",8,0,8);
   TProfile *precdepth = new TProfile("precdepth","",8,0,8);
   TProfile *prcdepth = new TProfile("prcdepth","",8,0,8);
   TProfile *precrcdepth = new TProfile("precrcdepth","",8,0,8);
   
   TProfile *pgendepth_hh = new TProfile("pgendepth_hh","",8,0,8);
   TProfile *precdepth_hh = new TProfile("precdepth_hh","",8,0,8);
   TProfile *prcdepth_hh = new TProfile("prcdepth_hh","",8,0,8);
   TProfile *precrcdepth_hh = new TProfile("precrcdepth_hh","",8,0,8);

   TProfile *pgendepth_eh = new TProfile("pgendepth_eh","",8,0,8);
   TProfile *precdepth_eh = new TProfile("precdepth_eh","",8,0,8);
   TProfile *prcdepth_eh = new TProfile("prcdepth_eh","",8,0,8);
   TProfile *precrcdepth_eh = new TProfile("precrcdepth_eh","",8,0,8);

   TProfile *pgendepth_ee = new TProfile("pgendepth_ee","",8,0,8);
   TProfile *precdepth_ee = new TProfile("precdepth_ee","",8,0,8);
   TProfile *prcdepth_ee = new TProfile("prcdepth_ee","",8,0,8);
   TProfile *precrcdepth_ee = new TProfile("precrcdepth_ee","",8,0,8);


   // Studies of how energy is distributed around gen particle in ieta-iphi
   double sumetaphigen(0), sumetaphirec(0), sumetaphirc(0);
   TH2D *h2etaphigen = new TH2D("h2etaphigen","",14,-7.5,6.5,14,-7.5,6.5);
   TH2D *h2etaphirec = new TH2D("h2etaphirec","",14,-7.5,6.5,14,-7.5,6.5);
   TH2D *h2etaphirc = new TH2D("h2etaphirc","",14,-7.5,6.5,14,-7.5,6.5);
   TH2D *h2etaphigen_1 = new TH2D("h2etaphigen_1","",14,-7.5,6.5,14,-7.5,6.5);
   TH2D *h2etaphirec_1 = new TH2D("h2etaphirec_1","",14,-7.5,6.5,14,-7.5,6.5);
   TH2D *h2etaphirc_1 = new TH2D("h2etaphirc_1","",14,-7.5,6.5,14,-7.5,6.5);
   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%10000==0) cout << "." << flush;

      // Task 1: find gen particle in ieta=28
      genPart_ieta = genPart_extrapolated_HCAL_ieta;
      genPart_iphi = genPart_extrapolated_HCAL_iphi;
      if (genPart_ieta!=28) continue;
      p4g.SetPtEtaPhiM(genPart_pt, genPart_eta, genPart_phi, 0);
      p4rc.SetPtEtaPhiM(genPart_pt, genPart_eta, genPart_phi+TMath::Pi(), 0);
            
      // Task 2: find energy deposits to ECAL (depth=0) and each of HCAL depths
      genPart_sum = 0;
      for (int i = 0; i != nDepth; ++i) {
	genPart_depth[i] = 0;
      } // for Depth

      // ECAL endcap
      for (int i = 0; i != nSimHitEE; ++i) {
	double erg = (*SimHitEE_energy)[i];
	double pts = erg / cosh((*SimHitEE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*SimHitEE_eta)[i], (*SimHitEE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  genPart_sum += erg;
	  genPart_depth[0] += erg;
	}
      } // for SimhitEE

      // HCAL endcap
      for (int i = 0; i != nSimHitHBHE; ++i) {
	double k = (*SimHitHBHE_samplingFactor)[i];
	double erg = k * (*SimHitHBHE_energy)[i];
	double pts = erg / cosh((*SimHitHBHE_eta)[i]);
	int depth = (*SimHitHBHE_depth)[i];
	p4s.SetPtEtaPhiM(pts, (*SimHitHBHE_eta)[i], (*SimHitHBHE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  genPart_sum += erg;
	  genPart_depth[depth] += erg;
	}
      } // for SimHitHBHE

      // Task 2b: repeat histogramming for H, EH split
      genPart_fe = (genPart_sum!=0 ? genPart_depth[0] / genPart_sum : 0.);
      bool ishh = (genPart_fe<0.05);
      bool iseh = (genPart_fe>=0.05 && genPart_sum-genPart_depth[0]>2.5);
      bool isee = (genPart_fe>=0.05 && genPart_sum-genPart_depth[0]<=2.5);
      
      // ECAL endcap
      for (int i = 0; i != nSimHitEE; ++i) {
	double erg = (*SimHitEE_energy)[i];
	double pts = erg / cosh((*SimHitEE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*SimHitEE_eta)[i], (*SimHitEE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  h2gendepth->Fill(0., erg, erg); gensum += erg;
	  if (ishh) { h2gendepth_hh->Fill(0., erg, erg); gensum_hh += erg; }
	  if (iseh) { h2gendepth_eh->Fill(0., erg, erg); gensum_eh += erg; }
	  if (isee) { h2gendepth_ee->Fill(0., erg, erg); gensum_ee += erg; }
	}
      } // for SimhitEE

      // HCAL endcap
      for (int i = 0; i != nSimHitHBHE; ++i) {
	double k = (*SimHitHBHE_samplingFactor)[i];
	double erg = k * (*SimHitHBHE_energy)[i];
	double pts = erg / cosh((*SimHitHBHE_eta)[i]);
	int depth = (*SimHitHBHE_depth)[i];
	p4s.SetPtEtaPhiM(pts, (*SimHitHBHE_eta)[i], (*SimHitHBHE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  //genPart_sum += erg;
	  h2gendepth->Fill(depth, erg, erg); gensum += erg;
	  if (ishh) { h2gendepth_hh->Fill(depth, erg, erg); gensum_hh += erg; }
	  if (iseh) { h2gendepth_eh->Fill(depth, erg, erg); gensum_eh += erg; }
	  if (isee) { h2gendepth_ee->Fill(depth, erg, erg); gensum_ee += erg; }
	}
	int dieta = (*SimHitHBHE_ieta)[i]-genPart_ieta;
	int diphi = (*SimHitHBHE_iphi)[i]-genPart_iphi;
	diphi = (diphi>36 ? diphi-72 : (diphi<-36 ? diphi+72 : diphi));
	if (ishh && p4g.DeltaR(p4s)<deltaRmax) {
	  h2etaphigen->Fill(dieta, diphi, erg);
	  if (depth==1) h2etaphigen_1->Fill(dieta, diphi, erg);
	  //if (fabs(dieta)<6 && fabs(diphi)<6)
	  sumetaphigen += erg;
	}
      } // for SimHitHBHE
      
      
      // Task 3: find reconstructed ECAL (depth=0) and HCAL rechits per depth
      // Bonus: find "random cone" energy for noise (and PU) subtraction
      // Bonus2: build towers out of HCAL RecHits
      int nRH(0), nRHh(0), nRHe(0), nRHrc(0), nRHrce(0), nRHrch(0);
      recHit_sum = 0;
      randomCone_sum = 0;
      map<int, map<int, double> > mt;
      map<int, map<int, double> > mtrc;
      for (int i = 0; i != nDepth; ++i) {
	recHit_depth[i] = 0;
	randomCone_depth[i] = 0;
      } // for nDepth

      // ECAL endcap
      for (int i = 0; i != nPFRecHitEE; ++i) {
	double erg = (*PFRecHitEE_energy)[i];
	double pts = erg / cosh((*PFRecHitEE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFRecHitEE_eta)[i], (*PFRecHitEE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  ++nRH; ++nRHe;
	  recHit_sum += erg;
	  recHit_depth[0] += erg;
	  h2recdepth->Fill(0., erg, erg); recsum += erg;
	  if (ishh) { h2recdepth_hh->Fill(0., erg, erg); recsum_hh += erg; }
	  if (iseh) { h2recdepth_eh->Fill(0., erg, erg); recsum_eh += erg; }
	  if (isee) { h2recdepth_ee->Fill(0., erg, erg); recsum_ee += erg; }
	  //
	  h2recrcdepth->Fill(0., erg, erg); recrcsum += erg;
	  if (ishh) { h2recrcdepth_hh->Fill(0., erg, erg); recrcsum_hh += erg; }
	  if (iseh) { h2recrcdepth_eh->Fill(0., erg, erg); recrcsum_eh += erg; }
	  if (isee) { h2recrcdepth_ee->Fill(0., erg, erg); recrcsum_ee += erg; }
	}
	if (p4rc.DeltaR(p4s)<deltaRmax) {
	  ++nRHrc; ++nRHrce;
	  randomCone_sum += erg;
	  randomCone_depth[0] += erg;
	  h2rcdepth->Fill(0., erg, erg); rcsum += erg;
	  if (ishh) { h2rcdepth_hh->Fill(0., erg, erg); rcsum_hh += erg; }
	  if (iseh) { h2rcdepth_eh->Fill(0., erg, erg); rcsum_eh += erg; }
	  if (isee) { h2rcdepth_ee->Fill(0., erg, erg); rcsum_ee += erg; }
	  //
	  h2recrcdepth->Fill(0., erg, -erg); recrcsum -= erg;
	  if (ishh) { h2recrcdepth_hh->Fill(0., erg,-erg); recrcsum_hh -= erg; }
	  if (iseh) { h2recrcdepth_eh->Fill(0., erg,-erg); recrcsum_eh -= erg; }
	  if (isee) { h2recrcdepth_ee->Fill(0., erg,-erg); recrcsum_ee -= erg; }
	}
      } // for PFRecHitEE

      // HCAL endcap
      for (int i = 0; i != nPFRecHitHBHE; ++i) {
	double erg = (*PFRecHitHBHE_energy)[i];
	double pts = erg / cosh((*PFRecHitHBHE_eta)[i]);
	int depth = (*PFRecHitHBHE_depth)[i];
	int ieta = (*PFRecHitHBHE_ieta)[i];
	int iphi = (*PFRecHitHBHE_iphi)[i];
	p4s.SetPtEtaPhiM(pts, (*PFRecHitHBHE_eta)[i], (*PFRecHitHBHE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  ++nRH; ++nRHh;
	  recHit_sum += erg;
	  recHit_depth[depth] += erg;
	  h2recdepth->Fill(depth, erg, erg); recsum += erg;
	  if (ishh) { h2recdepth_hh->Fill(depth, erg, erg); recsum_hh += erg; }
	  if (iseh) { h2recdepth_eh->Fill(depth, erg, erg); recsum_eh += erg; }
	  if (isee) { h2recdepth_ee->Fill(depth, erg, erg); recsum_ee += erg; }
	  //
	  h2recrcdepth->Fill(depth, erg, erg); recrcsum += erg;
	  if (ishh) { h2recrcdepth_hh->Fill(depth,erg,erg); recrcsum_hh+=erg; }
	  if (iseh) { h2recrcdepth_eh->Fill(depth,erg,erg); recrcsum_eh+=erg; }
	  if (isee) { h2recrcdepth_ee->Fill(depth,erg,erg); recrcsum_ee+=erg; }
	  // Build tower
	  mt[ieta][iphi] += erg;
	}
	int dieta = ieta - -genPart_ieta;
	int diphi = iphi - genPart_iphi;
	diphi = (diphi>36 ? diphi-72 : (diphi<-36 ? diphi+72 : diphi));
	if (ishh && p4g.DeltaR(p4s)<deltaRmax) {
	  h2etaphirec->Fill(dieta, diphi, erg);
	  if (depth==1) h2etaphirec_1->Fill(dieta, diphi, erg);
	  //if (fabs(dieta)<6 && fabs(diphi)<6)
	  sumetaphirec += erg;
	}
	
	if (p4rc.DeltaR(p4s)<deltaRmax) {
	  ++nRHrc; ++nRHrch;
	  randomCone_sum += erg;
	  randomCone_depth[depth] += erg;
	  h2rcdepth->Fill(depth, erg, erg); rcsum += erg;
	  if (ishh) { h2rcdepth_hh->Fill(depth, erg, erg); rcsum_hh += erg; }
	  if (iseh) { h2rcdepth_eh->Fill(depth, erg, erg); rcsum_eh += erg; }
	  if (isee) { h2rcdepth_ee->Fill(depth, erg, erg); rcsum_ee += erg; }
	  //
	  h2recrcdepth->Fill(depth, erg, -erg); recrcsum -= erg;
	  if (ishh) { h2recrcdepth_hh->Fill(depth,erg,-erg); recrcsum_hh-=erg; }
	  if (iseh) { h2recrcdepth_eh->Fill(depth,erg,-erg); recrcsum_eh-=erg; }
	  if (isee) { h2recrcdepth_ee->Fill(depth,erg,-erg); recrcsum_ee-=erg; }
	  // Build tower
	  mtrc[ieta][iphi] += erg;
	}
	int iphirc = genPart_iphi+36;
	if (iphirc>72) iphirc -= 72;
	int diphirc = (*PFRecHitHBHE_iphi)[i]-iphirc;
	diphirc = (diphirc>36 ? diphirc-72 : (diphirc<-36 ? diphirc+72 : diphirc));
	if (ishh && p4rc.DeltaR(p4s)<deltaRmax) {
	  h2etaphirc->Fill(dieta, diphirc, erg);
	  if (depth==1) h2etaphirc_1->Fill(dieta, diphirc, erg);
	  //if (fabs(dieta)<6 && fabs(diphirc)<6)
	  sumetaphirc += erg;
	}
      } // for PFRecHitHBHE

      // Task 3b: find the leading and subleading towers
      int nTWR(0), nTWRrc(0);
      double eTWRh(0), eTWRh2(0), eTWRh3(0);
      double eTWRrc(0), eTWRrc2(0), eTWRrc3(0);
      typedef map<int, map<int, double> >::const_iterator IT;
      typedef map<int, double>::const_iterator JT;
      for (IT it = mt.begin(); it != mt.end(); ++it) {
	map<int, double> const& mt2 = it->second;
	for (JT jt = mt2.begin(); jt != mt2.end(); ++jt) {

	  // Find leading and subleading clusters
	  double erg = jt->second;
	  ++nTWR;
	  htwr->Fill(erg);
	  if (ishh) htwr_hh->Fill(erg);
	  if (iseh) htwr_eh->Fill(erg);
	  if (isee) htwr_ee->Fill(erg);
	  if (erg>eTWRh) {
	    eTWRh3 = eTWRh2;
	    eTWRh2 = eTWRh;
	    eTWRh = erg;
	  }
	  else if (erg>eTWRh2) {
	    eTWRh3 = eTWRh2;
	    eTWRh2 = erg;
	  }
	  else if (erg>eTWRh3) {
	    eTWRh3 = erg;
	  }
	} // for JT
      } // for IT

      // Same for random cone
      for (IT it = mtrc.begin(); it != mtrc.end(); ++it) {
	map<int, double> const& mt2 = it->second;
	for (JT jt = mt2.begin(); jt != mt2.end(); ++jt) {
	  
	  // Find leading and subleading clusters
	  double erg = jt->second;
	  ++nTWRrc;
	  htwrrc->Fill(erg);
	  if (ishh) htwrrc_hh->Fill(erg);
	  if (iseh) htwrrc_eh->Fill(erg);
	  if (isee) htwrrc_ee->Fill(erg);
	  if (erg>eTWRrc) {
	    eTWRrc3 = eTWRrc2;
	    eTWRrc2 = eTWRrc;
	    eTWRrc = erg;
	  }
	  else if (erg>eTWRrc2) {
	    eTWRrc3 = eTWRrc2;
	    eTWRrc2 = erg;
	  }
	  else if (erg>eTWRrc3) {
	    eTWRrc3 = erg;
	  }
	} // for JT
      } // for IT
      
      // Task 4: find ECAL (depth=0) and HCAL PF clusters per depth
      // Bonus: find "random cone" energy for noise (and PU) subtraction
      int nPF(0), nPFh(0), nPFe(0), nPFrc(0), nPFrce(0), nPFrch(0);
      int iPFh(-1), iPFh2(-1), iPFh3(-1);
      int iPFrc(-1), iPFrc2(-1), iPFrc3(-1);
      double ePFh(0), ePFh2(0), ePFh3(0);
      double ePFrc(0), ePFrc2(0), ePFrc3(0);
      pfCluster_sum = 0;
      randomCone2_sum = 0;
      for (int i = 0; i != nDepth; ++i) {
	pfCluster_depth[i] = 0;
	randomCone2_depth[i] = 0;
      } // for nDepth

      // ECAL endcap
      for (int i = 0; i != nPFClusterECAL; ++i) {
	double pts = (*PFClusterECAL_energy)[i] / cosh((*PFClusterECAL_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFClusterECAL_eta)[i], (*PFClusterECAL_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  ++nPF; ++nPFe;
	  pfCluster_sum += (*PFClusterECAL_energy)[i];
	  pfCluster_depth[0] += (*PFClusterECAL_energy)[i];
	}
	if (p4rc.DeltaR(p4s)<deltaRmax) {
	  ++nPFrc; ++nPFrce;
	  randomCone2_sum += (*PFClusterECAL_energy)[i];
	  randomCone2_depth[0] += (*PFClusterECAL_energy)[i];
	}
      } // for PFClusterECAL

      // HCAL endcap
      for (int i = 0; i != nPFClusterHCAL; ++i) {
	double erg = (*PFClusterHCAL_energy)[i];
	double pts = erg / cosh((*PFClusterHCAL_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFClusterHCAL_eta)[i], (*PFClusterHCAL_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<deltaRmax) {
	  ++nPF; ++nPFh;
	  pfCluster_sum += erg;
	  pfCluster_depth[(*PFClusterHCAL_layer)[i]] += erg;

	  // Find leading and subleading clusters
	  if (erg>ePFh) {
	    iPFh3 = iPFh2; ePFh3 = ePFh2;
	    iPFh2 = iPFh;  ePFh2 = ePFh;
	    iPFh = i;      ePFh = erg;
	  }
	  else if (erg>ePFh2) {
	    iPFh3 = iPFh2; ePFh3 = ePFh2;
	    iPFh2 = i;     ePFh2 = erg;
	  }
	  else if (erg>ePFh3) {
	    iPFh3 = i;     ePFh3 = erg;
	  }
	}
	if (p4rc.DeltaR(p4s)<deltaRmax) {
	  ++nPFrc; ++nPFrch;
	  randomCone2_sum += erg;
	  randomCone2_depth[(*PFClusterHCAL_layer)[i]] += erg;

	  // Find leading and subleading clusters
	  if (erg>ePFrc) {
	    iPFrc3 = iPFrc2; ePFrc3 = ePFrc2;
	    iPFrc2 = iPFrc;  ePFrc2 = ePFrc;
	    iPFrc = i;      ePFrc = erg;
	  }
	  else if (erg>ePFrc2) {
	    iPFrc3 = iPFrc2; ePFrc3 = ePFrc2;
	    iPFrc2 = i;     ePFrc2 = erg;
	  }
	  else if (erg>ePFrc3) {
	    iPFrc3 = i;     ePFrc3 = erg;
	  }
	}
      } // for PFClusterHCAL

      //genPart_fe = (genPart_sum!=0 ? genPart_depth[0] / genPart_sum : 0.);
      recHit_fe = (recHit_sum!=0 ? recHit_depth[0] / recHit_sum : 0.);
      randomCone_fe = (randomCone_sum!=0 ? randomCone_depth[0] / randomCone_sum : 0.);
      pfCluster_fe = (pfCluster_sum!=0 ? pfCluster_depth[0] / pfCluster_sum : 0.);
      randomCone2_fe = (randomCone2_sum!=0 ? randomCone2_depth[0] / randomCone2_sum : 0.);
      
      // Fill the tree
      tree->Fill();

      hfe->Fill(genPart_fe);

      if (ishh) {
	hnrechit_hh->Fill(nRH);
	hnrechite_hh->Fill(nRHe);
	hnrechith_hh->Fill(nRHh);
	hnrechitrc_hh->Fill(nRHrc);
	hnrechitrce_hh->Fill(nRHrce);
	hnrechitrch_hh->Fill(nRHrch);
	
	hntwr_hh->Fill(nTWR);
	hntwrrc_hh->Fill(nTWRrc);
	
	hnpf_hh->Fill(nPF);
	hnpfe_hh->Fill(nPFe);
	hnpfh_hh->Fill(nPFh);
	hnpfrc_hh->Fill(nPFrc);
	hnpfrce_hh->Fill(nPFrce);
	hnpfrch_hh->Fill(nPFrch);
      }
	
      hgensum->Fill(genPart_sum);
      hrecsum->Fill(recHit_sum);
      hpfsum->Fill(pfCluster_sum);
      //hrcsum->Fill(randomCone2_sum);
      hrcsum->Fill(randomCone_sum);
      hrecrcsum->Fill(recHit_sum-randomCone_sum);
      hpfrcsum->Fill(pfCluster_sum-randomCone2_sum);
      
      for (int i = 0; i != nDepth; ++i) {
	//h2gendepth->Fill(i, genPart_depth[i]);
	//h2recdepth->Fill(i, recHit_depth[i]);
	//h2rcdepth->Fill(i, randomCone_depth[i]);
	
	pgendepth->Fill(i, genPart_depth[i]);
	precdepth->Fill(i, recHit_depth[i]);
	prcdepth->Fill(i, randomCone_depth[i]);
	precrcdepth->Fill(i, recHit_depth[i]-randomCone_depth[i]);
      } // for nDepth
      
      //if (recHit_fe<0.05) {
      //if (genPart_fe<0.05) {
      if (ishh) {
	hgensum_hh->Fill(genPart_sum);
	hrecsum_hh->Fill(recHit_sum);
	hpfsum_hh->Fill(pfCluster_sum);
	hrcsum_hh->Fill(randomCone_sum);
	hrecrcsum_hh->Fill(recHit_sum-randomCone_sum);
	hpfrcsum_hh->Fill(pfCluster_sum-randomCone2_sum);

	htwr1_hh->Fill(eTWRh);
	htwr2_hh->Fill(eTWRh2);
	htwr3_hh->Fill(eTWRh3);
	htwr4_hh->Fill(recHit_sum-(eTWRh+eTWRh2+eTWRh3));
	htwr12_hh->Fill(eTWRh+eTWRh2);
	htwr123_hh->Fill(eTWRh+eTWRh2+eTWRh3);

	htwr1rc_hh->Fill(eTWRrc);
	htwr2rc_hh->Fill(eTWRrc2);
	htwr3rc_hh->Fill(eTWRrc3);
	htwr4rc_hh->Fill(randomCone_sum-(eTWRrc+eTWRrc2+eTWRrc3));
	htwr12rc_hh->Fill(eTWRrc+eTWRrc2);
	htwr123rc_hh->Fill(eTWRrc+eTWRrc2+eTWRrc3);
	
	hpf1_hh->Fill(ePFh);
	hpf2_hh->Fill(ePFh2);
	hpf3_hh->Fill(ePFh3);
	hpf4_hh->Fill(pfCluster_sum-(ePFh+ePFh2+ePFh3));
	hpf12_hh->Fill(ePFh+ePFh2);
	hpf123_hh->Fill(ePFh+ePFh2+ePFh3);

	hpf1rc_hh->Fill(ePFrc);
	hpf2rc_hh->Fill(ePFrc2);
	hpf3rc_hh->Fill(ePFrc3);
	hpf4rc_hh->Fill(randomCone2_sum-(ePFrc+ePFrc2+ePFrc3));
	hpf12rc_hh->Fill(ePFrc+ePFrc2);
	hpf123rc_hh->Fill(ePFrc+ePFrc2+ePFrc3);

	for (int i = 0; i != nDepth; ++i) {
	  //h2gendepth_hh->Fill(i, genPart_depth[i]);
	  //h2recdepth_hh->Fill(i, recHit_depth[i]);
	  //h2rcdepth_hh->Fill(i, randomCone_depth[i]);

	  pgendepth_hh->Fill(i, genPart_depth[i]);
	  precdepth_hh->Fill(i, recHit_depth[i]);
	  prcdepth_hh->Fill(i, randomCone_depth[i]);
	  precrcdepth_hh->Fill(i, recHit_depth[i]-randomCone_depth[i]);
	} // for nDepth
      }
      //else if (genPart_fe<0.95) {
      // Less than 5 GeV (<10%) deposited in HCAL
      //else if (genPart_sum-genPart_depth[0]<2.5) {
      if (isee) {
      	hgensum_ee->Fill(genPart_sum);
	hrecsum_ee->Fill(recHit_sum);
	hpfsum_ee->Fill(pfCluster_sum);
	hrcsum_ee->Fill(randomCone_sum);
	hrecrcsum_ee->Fill(recHit_sum-randomCone_sum);
	hpfrcsum_ee->Fill(pfCluster_sum-randomCone2_sum);
      
	for (int i = 0; i != nDepth; ++i) {
	  //h2gendepth_ee->Fill(i, genPart_depth[i]);
	  //h2recdepth_ee->Fill(i, recHit_depth[i]);
	  //h2rcdepth_ee->Fill(i, randomCone_depth[i]);

	  pgendepth_ee->Fill(i, genPart_depth[i]);
	  precdepth_ee->Fill(i, recHit_depth[i]);
	  prcdepth_ee->Fill(i, randomCone_depth[i]);
	  precrcdepth_ee->Fill(i, recHit_depth[i]-randomCone_depth[i]);
	} // for nDepth
      }
      //else {
      if (iseh) {
      	hgensum_eh->Fill(genPart_sum);
	hrecsum_eh->Fill(recHit_sum);
	hpfsum_eh->Fill(pfCluster_sum);
	hrcsum_eh->Fill(randomCone_sum);
	hrecrcsum_eh->Fill(recHit_sum-randomCone_sum);
	hpfrcsum_eh->Fill(pfCluster_sum-randomCone2_sum);
      
	for (int i = 0; i != nDepth; ++i) {
	  //h2gendepth_eh->Fill(i, genPart_depth[i]);
	  //h2recdepth_eh->Fill(i, recHit_depth[i]);
	  //h2rcdepth_eh->Fill(i, randomCone_depth[i]);

	  pgendepth_eh->Fill(i, genPart_depth[i]);
	  precdepth_eh->Fill(i, recHit_depth[i]);
	  prcdepth_eh->Fill(i, randomCone_depth[i]);
	  precrcdepth_eh->Fill(i, recHit_depth[i]-randomCone_depth[i]);
	} // for nDepth
      }
	
   } // for jentry

   // Write the TTree to the file
   file->cd();
   tree->Write();
   
   // Close the file
   file->Close();

   h2gendepth->Scale(1./gensum);//h2gendepth->Integral());
   h2gendepth_hh->Scale(1./gensum_hh);//h2gendepth_hh->Integral());
   h2gendepth_eh->Scale(1./gensum_eh);//h2gendepth_eh->Integral());
   h2gendepth_ee->Scale(1./gensum_ee);//h2gendepth_ee->Integral());

   h2recdepth->Scale(1./recsum);//h2recdepth->Integral());
   h2recdepth_hh->Scale(1./recsum_hh);//h2recdepth_hh->Integral());
   h2recdepth_eh->Scale(1./recsum_eh);//h2recdepth_eh->Integral());
   h2recdepth_ee->Scale(1./recsum_ee);//h2recdepth_ee->Integral());

   h2rcdepth->Scale(1./rcsum);//h2rcdepth->Integral());
   h2rcdepth_hh->Scale(1./rcsum_hh);//h2rcdepth_hh->Integral());
   h2rcdepth_eh->Scale(1./rcsum_eh);//h2rcdepth_eh->Integral());
   h2rcdepth_ee->Scale(1./rcsum_ee);//h2rcdepth_ee->Integral());

   h2recrcdepth->Scale(1./recrcsum);
   h2recrcdepth_hh->Scale(1./recrcsum_hh);
   h2recrcdepth_eh->Scale(1./recrcsum_eh);
   h2recrcdepth_ee->Scale(1./recrcsum_ee);

   h2etaphigen->Scale(1./sumetaphigen);
   h2etaphirec->Scale(1./sumetaphirec);
   // Normalize depth1 with sum of all depths
   h2etaphigen_1->Scale(1./sumetaphigen);
   h2etaphirec_1->Scale(1./sumetaphirec);
   // Normalize RC with signal cone energy to have same scale
   h2etaphirc->Scale(1./sumetaphirec);//sumetaphirc);
   h2etaphirc_1->Scale(1./sumetaphirec);//sumetaphirc);
   
   fout->cd();
   fout->Write();
   fout->Close();
   
   cout << endl;
   cout << "Histograms have been written to 'PFAna.root'." << endl;
   cout << "TTree 's' has been created and written to 'PFAnaTree.root'." << endl;
    
} // PFAna::Loop
