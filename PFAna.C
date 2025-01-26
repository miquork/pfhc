#define PFAna_cxx
#include "PFAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>


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

   fChain->SetBranchStatus("nPFRecHitEE",1);
   fChain->SetBranchStatus("PFRecHitEE_eta",1);
   fChain->SetBranchStatus("PFRecHitEE_phi",1);
   fChain->SetBranchStatus("PFRecHitEE_energy",1);

   fChain->SetBranchStatus("nPFRecHitHBHE",1);
   fChain->SetBranchStatus("PFRecHitHBHE_eta",1);
   fChain->SetBranchStatus("PFRecHitHBHE_phi",1);
   fChain->SetBranchStatus("PFRecHitHBHE_depth",1);
   fChain->SetBranchStatus("PFRecHitHBHE_energy",1);

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
   TH1D *hgensum = new TH1D("hgensum","",100,0,100);
   TH1D *hgensum_h = new TH1D("hgensum_h","",100,0,100);
   TH1D *hgensum_e = new TH1D("hgensum_e","",100,0,100);
   TH1D *hrecsum = new TH1D("hrecsum","",100,0,100);
   TH1D *hrecsum_h = new TH1D("hrecsum_h","",100,0,100);
   TH1D *hrecsum_e = new TH1D("hrecsum_e","",100,0,100);
   TH1D *hpfsum = new TH1D("hpfsum","",100,0,100);
   TH1D *hpfsum_h = new TH1D("hpfsum_h","",100,0,100);
   TH1D *hpfsum_e = new TH1D("hpfsum_e","",100,0,100);
   
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
	double pts = (*SimHitEE_energy)[i] / cosh((*SimHitEE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*SimHitEE_eta)[i], (*SimHitEE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<0.4) {
	  genPart_sum += (*SimHitEE_energy)[i];
	  genPart_depth[0] += (*SimHitEE_energy)[i];
	}
      } // for SimhitEE

      // HCAL endcap
      for (int i = 0; i != nSimHitHBHE; ++i) {
	double k = (*SimHitHBHE_samplingFactor)[i];
	double pts = k * (*SimHitHBHE_energy)[i] / cosh((*SimHitHBHE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*SimHitHBHE_eta)[i], (*SimHitHBHE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<0.4) {
	  genPart_sum += k * (*SimHitHBHE_energy)[i];
	  genPart_depth[(*SimHitHBHE_depth)[i]] += k * (*SimHitHBHE_energy)[i];
	}
      } // for SimHitHBHE

      // Task 3: find reconstructed ECAL (depth=0) and HCAL rechits per depth
      // Bonus: find "random cone" energy for noise (and PU) subtraction
      recHit_sum = 0;
      randomCone_sum = 0;
      for (int i = 0; i != nDepth; ++i) {
	recHit_depth[i] = 0;
	randomCone_depth[i] = 0;
      } // for nDepth

      // ECAL endcap
      for (int i = 0; i != nPFRecHitEE; ++i) {
	double pts = (*PFRecHitEE_energy)[i] / cosh((*PFRecHitEE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFRecHitEE_eta)[i], (*PFRecHitEE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<0.4) {
	  recHit_sum += (*PFRecHitEE_energy)[i];
	  recHit_depth[0] += (*PFRecHitEE_energy)[i];
	}
	if (p4rc.DeltaR(p4s)<0.4) {
	  randomCone_sum += (*PFRecHitEE_energy)[i];
	  randomCone_depth[0] += (*PFRecHitEE_energy)[i];
	}
      } // for PFRecHitEE

      // HCAL endcap
      for (int i = 0; i != nPFRecHitHBHE; ++i) {
	double pts = (*PFRecHitHBHE_energy)[i] / cosh((*PFRecHitHBHE_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFRecHitHBHE_eta)[i], (*PFRecHitHBHE_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<0.4) {
	  recHit_sum += (*PFRecHitHBHE_energy)[i];
	  recHit_depth[(*PFRecHitHBHE_depth)[i]] += (*PFRecHitHBHE_energy)[i];
	}
	if (p4rc.DeltaR(p4s)<0.4) {
	  randomCone_sum += (*PFRecHitHBHE_energy)[i];
	  randomCone_depth[(*PFRecHitHBHE_depth)[i]] += (*PFRecHitHBHE_energy)[i];
	}
      } // for PFRecHitHBHE


      // Task 4: find ECAL (depth=0) and HCAL PF clusters per depth
      // Bonus: find "random cone" energy for noise (and PU) subtraction
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
	if (p4g.DeltaR(p4s)<0.4) {
	  pfCluster_sum += (*PFClusterECAL_energy)[i];
	  pfCluster_depth[0] += (*PFClusterECAL_energy)[i];
	}
	if (p4rc.DeltaR(p4s)<0.4) {
	  randomCone2_sum += (*PFClusterECAL_energy)[i];
	  randomCone2_depth[0] += (*PFClusterECAL_energy)[i];
	}
      } // for PFClusterECAL

      // HCAL endcap
      for (int i = 0; i != nPFClusterHCAL; ++i) {
	double pts = (*PFClusterHCAL_energy)[i] / cosh((*PFClusterHCAL_eta)[i]);
	p4s.SetPtEtaPhiM(pts, (*PFClusterHCAL_eta)[i], (*PFClusterHCAL_phi)[i], 0.);
	if (p4g.DeltaR(p4s)<0.4) {
	  pfCluster_sum += (*PFClusterHCAL_energy)[i];
	  pfCluster_depth[(*PFClusterHCAL_layer)[i]] += (*PFClusterHCAL_energy)[i];
	}
	if (p4rc.DeltaR(p4s)<0.4) {
	  randomCone2_sum += (*PFClusterHCAL_energy)[i];
	  randomCone2_depth[(*PFClusterHCAL_layer)[i]] += (*PFClusterHCAL_energy)[i];
	}
      } // for PFClusterHCAL

      genPart_fe = (genPart_sum!=0 ? genPart_depth[0] / genPart_sum : 0.);
      recHit_fe = (recHit_sum!=0 ? recHit_depth[0] / recHit_sum : 0.);
      randomCone_fe = (randomCone_sum!=0 ? randomCone_depth[0] / randomCone_sum : 0.);
      pfCluster_fe = (pfCluster_sum!=0 ? pfCluster_depth[0] / pfCluster_sum : 0.);
      randomCone2_fe = (randomCone2_sum!=0 ? randomCone2_depth[0] / randomCone2_sum : 0.);
      
      // Fill the tree
      tree->Fill();

      hfe->Fill(genPart_fe);
      
      hgensum->Fill(genPart_sum);
      hrecsum->Fill(recHit_sum);
      hpfsum->Fill(pfCluster_sum);
      if (genPart_fe<0.05) {
	hgensum_h->Fill(genPart_sum);
	hrecsum_h->Fill(recHit_sum);
	hpfsum_h->Fill(pfCluster_sum);
      }
      else {
      	hgensum_e->Fill(genPart_sum);
	hrecsum_e->Fill(recHit_sum);
	hpfsum_e->Fill(pfCluster_sum);
      }
	
   } // for jentry

   // Write the TTree to the file
   file->cd();
   tree->Write();
   
   // Close the file
   file->Close();

   fout->cd();
   fout->Write();
   fout->Close();
   
   cout << endl;
   cout << "Histograms have been written to 'PFAna.root'." << endl;
   cout << "TTree 's' has been created and written to 'PFAnaTree.root'." << endl;
    
} // PFAna::Loop
