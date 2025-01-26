//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 25 17:08:53 2025 by ROOT version 6.28/06
// from TTree PFAnaTree/PFAnaTree
// found on file: MergedNtuple_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root
//////////////////////////////////////////////////////////

#ifndef PFAna_h
#define PFAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class PFAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       run;
   ULong64_t       evt;
   ULong64_t       lumiBlock;
   ULong64_t       time;
   Float_t         genVertex_x;
   Float_t         genVertex_y;
   Float_t         genVertex_z;
   Float_t         genPart_pt;
   Float_t         genPart_eta;
   Float_t         genPart_phi;
   Float_t         genPart_mass;
   Float_t         genPart_energy;
   Float_t         genPart_p;
   Int_t           genPart_pdgId;
   Int_t           genPart_charge;
   Bool_t          genPart_extrapolated_okECAL;
   Bool_t          genPart_extrapolated_okHCAL;
   Int_t           genPart_extrapolated_EB_ieta;
   Int_t           genPart_extrapolated_EB_iphi;
   Float_t         genPart_extrapolated_EB_eta;
   Float_t         genPart_extrapolated_EB_phi;
   Int_t           genPart_extrapolated_PFRecHitEB_Idx;
   UInt_t          genPart_extrapolated_PFRecHitEB_detId;
   Int_t           genPart_extrapolated_EE_ix;
   Int_t           genPart_extrapolated_EE_iy;
   Float_t         genPart_extrapolated_EE_eta;
   Float_t         genPart_extrapolated_EE_phi;
   Int_t           genPart_extrapolated_PFRecHitEE_Idx;
   UInt_t          genPart_extrapolated_PFRecHitEE_detId;
   Int_t           genPart_extrapolated_EHCAL_ieta;
   Int_t           genPart_extrapolated_EHCAL_iphi;
   Int_t           genPart_extrapolated_EHCAL_depth;
   Float_t         genPart_extrapolated_EHCAL_eta;
   Float_t         genPart_extrapolated_EHCAL_phi;
   Int_t           genPart_extrapolated_EHCAL_subdetId;
   Int_t           genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx;
   UInt_t          genPart_extrapolated_EHCAL_PFRecHitHBHE_detId;
   Int_t           genPart_extrapolated_EHCAL_PFRecHitHF_Idx;
   UInt_t          genPart_extrapolated_EHCAL_PFRecHitHF_detId;
   Int_t           genPart_extrapolated_HCAL_ieta;
   Int_t           genPart_extrapolated_HCAL_iphi;
   Int_t           genPart_extrapolated_HCAL_depth;
   Float_t         genPart_extrapolated_HCAL_eta;
   Float_t         genPart_extrapolated_HCAL_phi;
   Int_t           genPart_extrapolated_HCAL_subdetId;
   Int_t           genPart_extrapolated_HCAL_PFRecHitHBHE_Idx;
   UInt_t          genPart_extrapolated_HCAL_PFRecHitHBHE_detId;
   Int_t           genPart_extrapolated_HCAL_PFRecHitHF_Idx;
   UInt_t          genPart_extrapolated_HCAL_PFRecHitHF_detId;
   Float_t         genPart_extrapolated_pointECAL_x;
   Float_t         genPart_extrapolated_pointECAL_y;
   Float_t         genPart_extrapolated_pointECAL_z;
   Float_t         genPart_extrapolated_pointECAL_eta;
   Float_t         genPart_extrapolated_pointECAL_phi;
   Float_t         genPart_extrapolated_pointHCAL_x;
   Float_t         genPart_extrapolated_pointHCAL_y;
   Float_t         genPart_extrapolated_pointHCAL_z;
   Float_t         genPart_extrapolated_pointHCAL_eta;
   Float_t         genPart_extrapolated_pointHCAL_phi;
   Int_t           nPFCand;
   vector<float>   *PFCand_pt;
   vector<float>   *PFCand_eta;
   vector<float>   *PFCand_phi;
   vector<float>   *PFCand_mass;
   vector<float>   *PFCand_energy;
   vector<int>     *PFCand_pdgId;
   vector<int>     *PFCand_isChgHadIso;
   vector<unsigned int> *PFCand_birthId;
   vector<float>   *PFCand_dRToGenPart;
   vector<int>     *PFCand_has_trk;
   vector<float>   *PFCand_trk_pt;
   vector<float>   *PFCand_trk_p;
   vector<float>   *PFCand_trk_eta;
   vector<float>   *PFCand_trk_phi;
   vector<int>     *PFCand_trk_algo;
   vector<int>     *PFCand_trk_charge;
   vector<float>   *PFCand_trk_ptError;
   vector<float>   *PFCand_trk_dxy;
   vector<float>   *PFCand_trk_dz;
   vector<float>   *PFCand_trk_dxyError;
   vector<float>   *PFCand_trk_dzError;
   vector<float>   *PFCand_trk_chi2;
   vector<float>   *PFCand_trk_ndof;
   vector<float>   *PFCand_trk_vx;
   vector<float>   *PFCand_trk_vy;
   vector<float>   *PFCand_trk_vz;
   vector<int>     *PFCand_trk_numberOfValidPixelHits;
   vector<int>     *PFCand_trk_numberOfValidTrackerHits;
   vector<int>     *PFCand_trk_pixelHitOK;
   vector<int>     *PFCand_trk_trackerHitOK;
   vector<float>   *PFCand_trk_positionAtECALEntrance_x;
   vector<float>   *PFCand_trk_positionAtECALEntrance_y;
   vector<float>   *PFCand_trk_positionAtECALEntrance_z;
   vector<float>   *PFCand_trk_positionAtECALEntrance_eta;
   vector<float>   *PFCand_trk_positionAtECALEntrance_phi;
   vector<float>   *PFCand_trk_closestECAL_eta;
   vector<float>   *PFCand_trk_closestECAL_phi;
   vector<int>     *PFCand_trk_closestECAL_detId;
   vector<float>   *PFCand_trk_closestECAL_EBieta;
   vector<float>   *PFCand_trk_closestECAL_EBiphi;
   vector<float>   *PFCand_trk_closestECAL_EEix;
   vector<float>   *PFCand_trk_closestECAL_EEiy;
   vector<float>   *PFCand_trk_closestHCAL_eta;
   vector<float>   *PFCand_trk_closestHCAL_phi;
   vector<int>     *PFCand_trk_closestHCAL_detId;
   vector<float>   *PFCand_trk_closestHCAL_ieta;
   vector<float>   *PFCand_trk_closestHCAL_iphi;
   vector<float>   *PFCand_ecalEnergy;
   vector<float>   *PFCand_hcalEnergy;
   vector<float>   *PFCand_caloEnergy;
   vector<float>   *PFCand_rawEcalEnergy;
   vector<float>   *PFCand_rawHcalEnergy;
   vector<float>   *PFCand_rawCaloEnergy;
   vector<float>   *PFCand_hoEnergy;
   vector<float>   *PFCand_rawHoEnergy;
   vector<float>   *PFCand_pS1Energy;
   vector<float>   *PFCand_p21Energy;
   vector<float>   *PFCand_hcalDepth1EFrac;
   vector<float>   *PFCand_hcalDepth2EFrac;
   vector<float>   *PFCand_hcalDepth3EFrac;
   vector<float>   *PFCand_hcalDepth4EFrac;
   vector<float>   *PFCand_hcalDepth5EFrac;
   vector<float>   *PFCand_hcalDepth6EFrac;
   vector<float>   *PFCand_hcalDepth7EFrac;
   vector<int>     *PFCand_PFBlock_Idx;
   vector<unsigned int> *PFCand_nPFTrack;
   vector<unsigned int> *PFCand_nPFClusterHCAL;
   vector<unsigned int> *PFCand_nPFClusterECAL;
   vector<unsigned int> *PFCand_nPFClusterPS;
   vector<unsigned int> *PFCand_nPFClusterHF;
   vector<unsigned int> *PFCand_nPFClusterSC;
   vector<unsigned int> *PFCand_nPFClusterGSF;
   vector<unsigned int> *PFCand_nPFClusterBREM;
   vector<vector<int> > *PFCand_PFClusterHCAL_Idx;
   vector<vector<int> > *PFCand_PFClusterECAL_Idx;
   vector<vector<int> > *PFCand_PFClusterPS_Idx;
   vector<vector<int> > *PFCand_PFClusterHF_Idx;
   vector<unsigned int> *PFCand_nPFTrackInBlock;
   vector<unsigned int> *PFCand_nPFClusterHCALInBlock;
   vector<unsigned int> *PFCand_nPFClusterECALInBlock;
   vector<unsigned int> *PFCand_nPFClusterPSInBlock;
   vector<unsigned int> *PFCand_nPFClusterHFInBlock;
   vector<unsigned int> *PFCand_nPFClusterSCInBlock;
   vector<unsigned int> *PFCand_nPFClusterGSFInBlock;
   vector<unsigned int> *PFCand_nPFClusterBREMInBlock;
   vector<vector<int> > *PFCand_PFClusterHCALInBlock_Idx;
   vector<vector<int> > *PFCand_PFClusterECALInBlock_Idx;
   vector<vector<int> > *PFCand_PFClusterPSInBlock_Idx;
   vector<vector<int> > *PFCand_PFClusterHFInBlock_Idx;
   Int_t           Idx_ClosestPFCandChgHad;
   Int_t           Idx_ClosestPFCandChgHadIso;
   Int_t           Idx_ClosestPFCandNeuHad;
   Int_t           Idx_ClosestPFCandPhoton;
   Int_t           nPFRecHitHBHE;
   vector<float>   *PFRecHitHBHE_energy;
   vector<int>     *PFRecHitHBHE_ieta;
   vector<int>     *PFRecHitHBHE_iphi;
   vector<float>   *PFRecHitHBHE_eta;
   vector<float>   *PFRecHitHBHE_phi;
   vector<int>     *PFRecHitHBHE_depth;
   vector<float>   *PFRecHitHBHE_cutThreshold;
   vector<float>   *PFRecHitHBHE_hcalRespCorr;
   vector<unsigned int> *PFRecHitHBHE_detId;
   Int_t           nPFRecHitEB;
   vector<float>   *PFRecHitEB_energy;
   vector<int>     *PFRecHitEB_ieta;
   vector<int>     *PFRecHitEB_iphi;
   vector<int>     *PFRecHitEB_tower_ieta;
   vector<int>     *PFRecHitEB_tower_iphi;
   vector<float>   *PFRecHitEB_approxEta;
   vector<float>   *PFRecHitEB_eta;
   vector<float>   *PFRecHitEB_phi;
   vector<float>   *PFRecHitEB_cutThreshold;
   vector<unsigned int> *PFRecHitEB_detId;
   Int_t           nPFRecHitEE;
   vector<float>   *PFRecHitEE_energy;
   vector<int>     *PFRecHitEE_ix;
   vector<int>     *PFRecHitEE_iy;
   vector<float>   *PFRecHitEE_eta;
   vector<float>   *PFRecHitEE_phi;
   vector<float>   *PFRecHitEE_cutThreshold;
   vector<unsigned int> *PFRecHitEE_detId;
   Int_t           nSimHitHBHE;
   vector<float>   *SimHitHBHE_energy;
   vector<float>   *SimHitHBHE_energyEM;
   vector<float>   *SimHitHBHE_energyHAD;
   vector<float>   *SimHitHBHE_samplingFactor;
   vector<int>     *SimHitHBHE_nCaloHits;
   vector<int>     *SimHitHBHE_ieta;
   vector<int>     *SimHitHBHE_iphi;
   vector<int>     *SimHitHBHE_depth;
   vector<float>   *SimHitHBHE_eta;
   vector<float>   *SimHitHBHE_phi;
   vector<float>   *SimHitHBHE_detId;
   vector<int>     *SimHitHBHE_subdetId;
   vector<float>   *SimHitHBHE_energy_25ns;
   vector<float>   *SimHitHBHE_energyEM_25ns;
   vector<float>   *SimHitHBHE_energyHAD_25ns;
   Int_t           nSimHitEB;
   vector<float>   *SimHitEB_energy;
   vector<float>   *SimHitEB_energyEM;
   vector<float>   *SimHitEB_energyHAD;
   vector<int>     *SimHitEB_nCaloHits;
   vector<int>     *SimHitEB_ieta;
   vector<int>     *SimHitEB_iphi;
   vector<float>   *SimHitEB_eta;
   vector<float>   *SimHitEB_phi;
   vector<float>   *SimHitEB_detId;
   vector<int>     *SimHitEB_subdetId;
   Int_t           nSimHitEE;
   vector<float>   *SimHitEE_energy;
   vector<float>   *SimHitEE_energyEM;
   vector<float>   *SimHitEE_energyHAD;
   vector<int>     *SimHitEE_nCaloHits;
   vector<int>     *SimHitEE_ix;
   vector<int>     *SimHitEE_iy;
   vector<float>   *SimHitEE_eta;
   vector<float>   *SimHitEE_phi;
   vector<float>   *SimHitEE_detId;
   vector<int>     *SimHitEE_subdetId;
   Int_t           nPFClusterECAL;
   vector<float>   *PFClusterECAL_pt;
   vector<float>   *PFClusterECAL_energy;
   vector<float>   *PFClusterECAL_correctedEnergy;
   vector<float>   *PFClusterECAL_eta;
   vector<float>   *PFClusterECAL_phi;
   vector<int>     *PFClusterECAL_layer;
   vector<unsigned int> *PFClusterECAL_seedhit_detId;
   vector<int>     *PFClusterECAL_nhits;
   vector<vector<unsigned int> > *PFClusterECAL_hits_detId;
   vector<vector<float> > *PFClusterECAL_hits_fraction;
   vector<vector<int> > *PFClusterECAL_hits_PFRecHitEB_Idx;
   vector<vector<int> > *PFClusterECAL_hits_PFRecHitEE_Idx;
   vector<int>     *PFClusterECAL_key;
   Int_t           nPFClusterPS;
   vector<float>   *PFClusterPS_pt;
   vector<float>   *PFClusterPS_energy;
   vector<float>   *PFClusterPS_correctedEnergy;
   vector<float>   *PFClusterPS_eta;
   vector<float>   *PFClusterPS_phi;
   vector<int>     *PFClusterPS_layer;
   vector<unsigned int> *PFClusterPS_seedhit_detId;
   vector<int>     *PFClusterPS_nhits;
   vector<vector<unsigned int> > *PFClusterPS_hits_detId;
   vector<vector<float> > *PFClusterPS_hits_fraction;
   vector<int>     *PFClusterPS_key;
   Int_t           nPFClusterHCAL;
   vector<float>   *PFClusterHCAL_pt;
   vector<float>   *PFClusterHCAL_energy;
   vector<float>   *PFClusterHCAL_correctedEnergy;
   vector<float>   *PFClusterHCAL_eta;
   vector<float>   *PFClusterHCAL_phi;
   vector<int>     *PFClusterHCAL_layer;
   vector<unsigned int> *PFClusterHCAL_seedhit_detId;
   vector<int>     *PFClusterHCAL_nhits;
   vector<vector<unsigned int> > *PFClusterHCAL_hits_detId;
   vector<vector<float> > *PFClusterHCAL_hits_fraction;
   vector<vector<int> > *PFClusterHCAL_hits_PFRecHitHBHE_Idx;
   vector<int>     *PFClusterHCAL_key;

   // List of branches
   //TBranch        *b_orun;   //!
   TBranch        *b_genVertex_x;   //!
   TBranch        *b_genVertex_y;   //!
   TBranch        *b_genVertex_z;   //!
   TBranch        *b_genPart_pt;   //!
   TBranch        *b_genPart_eta;   //!
   TBranch        *b_genPart_phi;   //!
   TBranch        *b_genPart_mass;   //!
   TBranch        *b_genPart_energy;   //!
   TBranch        *b_genPart_p;   //!
   TBranch        *b_genPart_pdgId;   //!
   TBranch        *b_genPart_charge;   //!
   TBranch        *b_genPart_extrapolated_okECAL;   //!
   TBranch        *b_genPart_extrapolated_okHCAL;   //!
   TBranch        *b_genPart_extrapolated_EB_ieta;   //!
   TBranch        *b_genPart_extrapolated_EB_iphi;   //!
   TBranch        *b_genPart_extrapolated_EB_eta;   //!
   TBranch        *b_genPart_extrapolated_EB_phi;   //!
   TBranch        *b_genPart_extrapolated_PFRecHitEB_Idx;   //!
   TBranch        *b_genPart_extrapolated_PFRecHitEB_detId;   //!
   TBranch        *b_genPart_extrapolated_EE_ix;   //!
   TBranch        *b_genPart_extrapolated_EE_iy;   //!
   TBranch        *b_genPart_extrapolated_EE_eta;   //!
   TBranch        *b_genPart_extrapolated_EE_phi;   //!
   TBranch        *b_genPart_extrapolated_PFRecHitEE_Idx;   //!
   TBranch        *b_genPart_extrapolated_PFRecHitEE_detId;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_ieta;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_iphi;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_depth;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_eta;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_phi;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_subdetId;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_PFRecHitHBHE_detId;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_PFRecHitHF_Idx;   //!
   TBranch        *b_genPart_extrapolated_EHCAL_PFRecHitHF_detId;   //!
   TBranch        *b_genPart_extrapolated_HCAL_ieta;   //!
   TBranch        *b_genPart_extrapolated_HCAL_iphi;   //!
   TBranch        *b_genPart_extrapolated_HCAL_depth;   //!
   TBranch        *b_genPart_extrapolated_HCAL_eta;   //!
   TBranch        *b_genPart_extrapolated_HCAL_phi;   //!
   TBranch        *b_genPart_extrapolated_HCAL_subdetId;   //!
   TBranch        *b_genPart_extrapolated_HCAL_PFRecHitHBHE_Idx;   //!
   TBranch        *b_genPart_extrapolated_HCAL_PFRecHitHBHE_detId;   //!
   TBranch        *b_genPart_extrapolated_HCAL_PFRecHitHF_Idx;   //!
   TBranch        *b_genPart_extrapolated_HCAL_PFRecHitHF_detId;   //!
   TBranch        *b_genPart_extrapolated_pointECAL_x;   //!
   TBranch        *b_genPart_extrapolated_pointECAL_y;   //!
   TBranch        *b_genPart_extrapolated_pointECAL_z;   //!
   TBranch        *b_genPart_extrapolated_pointECAL_eta;   //!
   TBranch        *b_genPart_extrapolated_pointECAL_phi;   //!
   TBranch        *b_genPart_extrapolated_pointHCAL_x;   //!
   TBranch        *b_genPart_extrapolated_pointHCAL_y;   //!
   TBranch        *b_genPart_extrapolated_pointHCAL_z;   //!
   TBranch        *b_genPart_extrapolated_pointHCAL_eta;   //!
   TBranch        *b_genPart_extrapolated_pointHCAL_phi;   //!
   TBranch        *b_nPFCand;   //!
   TBranch        *b_PFCand_pt;   //!
   TBranch        *b_PFCand_eta;   //!
   TBranch        *b_PFCand_phi;   //!
   TBranch        *b_PFCand_mass;   //!
   TBranch        *b_PFCand_energy;   //!
   TBranch        *b_PFCand_pdgId;   //!
   TBranch        *b_PFCand_isChgHadIso;   //!
   TBranch        *b_PFCand_birthId;   //!
   TBranch        *b_PFCand_dRToGenPart;   //!
   TBranch        *b_PFCand_has_trk;   //!
   TBranch        *b_PFCand_trk_pt;   //!
   TBranch        *b_PFCand_trk_p;   //!
   TBranch        *b_PFCand_trk_eta;   //!
   TBranch        *b_PFCand_trk_phi;   //!
   TBranch        *b_PFCand_trk_algo;   //!
   TBranch        *b_PFCand_trk_charge;   //!
   TBranch        *b_PFCand_trk_ptError;   //!
   TBranch        *b_PFCand_trk_dxy;   //!
   TBranch        *b_PFCand_trk_dz;   //!
   TBranch        *b_PFCand_trk_dxyError;   //!
   TBranch        *b_PFCand_trk_dzError;   //!
   TBranch        *b_PFCand_trk_chi2;   //!
   TBranch        *b_PFCand_trk_ndof;   //!
   TBranch        *b_PFCand_trk_vx;   //!
   TBranch        *b_PFCand_trk_vy;   //!
   TBranch        *b_PFCand_trk_vz;   //!
   TBranch        *b_PFCand_trk_numberOfValidPixelHits;   //!
   TBranch        *b_PFCand_trk_numberOfValidTrackerHits;   //!
   TBranch        *b_PFCand_trk_pixelHitOK;   //!
   TBranch        *b_PFCand_trk_trackerHitOK;   //!
   TBranch        *b_PFCand_trk_positionAtECALEntrance_x;   //!
   TBranch        *b_PFCand_trk_positionAtECALEntrance_y;   //!
   TBranch        *b_PFCand_trk_positionAtECALEntrance_z;   //!
   TBranch        *b_PFCand_trk_positionAtECALEntrance_eta;   //!
   TBranch        *b_PFCand_trk_positionAtECALEntrance_phi;   //!
   TBranch        *b_PFCand_trk_closestECAL_eta;   //!
   TBranch        *b_PFCand_trk_closestECAL_phi;   //!
   TBranch        *b_PFCand_trk_closestECAL_detId;   //!
   TBranch        *b_PFCand_trk_closestECAL_EBieta;   //!
   TBranch        *b_PFCand_trk_closestECAL_EBiphi;   //!
   TBranch        *b_PFCand_trk_closestECAL_EEix;   //!
   TBranch        *b_PFCand_trk_closestECAL_EEiy;   //!
   TBranch        *b_PFCand_trk_closestHCAL_eta;   //!
   TBranch        *b_PFCand_trk_closestHCAL_phi;   //!
   TBranch        *b_PFCand_trk_closestHCAL_detId;   //!
   TBranch        *b_PFCand_trk_closestHCAL_ieta;   //!
   TBranch        *b_PFCand_trk_closestHCAL_iphi;   //!
   TBranch        *b_PFCand_ecalEnergy;   //!
   TBranch        *b_PFCand_hcalEnergy;   //!
   TBranch        *b_PFCand_caloEnergy;   //!
   TBranch        *b_PFCand_rawEcalEnergy;   //!
   TBranch        *b_PFCand_rawHcalEnergy;   //!
   TBranch        *b_PFCand_rawCaloEnergy;   //!
   TBranch        *b_PFCand_hoEnergy;   //!
   TBranch        *b_PFCand_rawHoEnergy;   //!
   TBranch        *b_PFCand_pS1Energy;   //!
   TBranch        *b_PFCand_p21Energy;   //!
   TBranch        *b_PFCand_hcalDepth1EFrac;   //!
   TBranch        *b_PFCand_hcalDepth2EFrac;   //!
   TBranch        *b_PFCand_hcalDepth3EFrac;   //!
   TBranch        *b_PFCand_hcalDepth4EFrac;   //!
   TBranch        *b_PFCand_hcalDepth5EFrac;   //!
   TBranch        *b_PFCand_hcalDepth6EFrac;   //!
   TBranch        *b_PFCand_hcalDepth7EFrac;   //!
   TBranch        *b_PFCand_PFBlock_Idx;   //!
   TBranch        *b_PFCand_nPFTrack;   //!
   TBranch        *b_PFCand_nPFClusterHCAL;   //!
   TBranch        *b_PFCand_nPFClusterECAL;   //!
   TBranch        *b_PFCand_nPFClusterPS;   //!
   TBranch        *b_PFCand_nPFClusterHF;   //!
   TBranch        *b_PFCand_nPFClusterSC;   //!
   TBranch        *b_PFCand_nPFClusterGSF;   //!
   TBranch        *b_PFCand_nPFClusterBREM;   //!
   TBranch        *b_PFCand_PFClusterHCAL_Idx;   //!
   TBranch        *b_PFCand_PFClusterECAL_Idx;   //!
   TBranch        *b_PFCand_PFClusterPS_Idx;   //!
   TBranch        *b_PFCand_PFClusterHF_Idx;   //!
   TBranch        *b_PFCand_nPFTrackInBlock;   //!
   TBranch        *b_PFCand_nPFClusterHCALInBlock;   //!
   TBranch        *b_PFCand_nPFClusterECALInBlock;   //!
   TBranch        *b_PFCand_nPFClusterPSInBlock;   //!
   TBranch        *b_PFCand_nPFClusterHFInBlock;   //!
   TBranch        *b_PFCand_nPFClusterSCInBlock;   //!
   TBranch        *b_PFCand_nPFClusterGSFInBlock;   //!
   TBranch        *b_PFCand_nPFClusterBREMInBlock;   //!
   TBranch        *b_PFCand_PFClusterHCALInBlock_Idx;   //!
   TBranch        *b_PFCand_PFClusterECALInBlock_Idx;   //!
   TBranch        *b_PFCand_PFClusterPSInBlock_Idx;   //!
   TBranch        *b_PFCand_PFClusterHFInBlock_Idx;   //!
   TBranch        *b_Idx_ClosestPFCandChgHad;   //!
   TBranch        *b_Idx_ClosestPFCandChgHadIso;   //!
   TBranch        *b_Idx_ClosestPFCandNeuHad;   //!
   TBranch        *b_Idx_ClosestPFCandPhoton;   //!
   TBranch        *b_nPFRecHitHBHE;   //!
   TBranch        *b_PFRecHitHBHE_energy;   //!
   TBranch        *b_PFRecHitHBHE_ieta;   //!
   TBranch        *b_PFRecHitHBHE_iphi;   //!
   TBranch        *b_PFRecHitHBHE_eta;   //!
   TBranch        *b_PFRecHitHBHE_phi;   //!
   TBranch        *b_PFRecHitHBHE_depth;   //!
   TBranch        *b_PFRecHitHBHE_cutThreshold;   //!
   TBranch        *b_PFRecHitHBHE_hcalRespCorr;   //!
   TBranch        *b_PFRecHitHBHE_detId;   //!
   TBranch        *b_nPFRecHitEB;   //!
   TBranch        *b_PFRecHitEB_energy;   //!
   TBranch        *b_PFRecHitEB_ieta;   //!
   TBranch        *b_PFRecHitEB_iphi;   //!
   TBranch        *b_PFRecHitEB_tower_ieta;   //!
   TBranch        *b_PFRecHitEB_tower_iphi;   //!
   TBranch        *b_PFRecHitEB_approxEta;   //!
   TBranch        *b_PFRecHitEB_eta;   //!
   TBranch        *b_PFRecHitEB_phi;   //!
   TBranch        *b_PFRecHitEB_cutThreshold;   //!
   TBranch        *b_PFRecHitEB_detId;   //!
   TBranch        *b_nPFRecHitEE;   //!
   TBranch        *b_PFRecHitEE_energy;   //!
   TBranch        *b_PFRecHitEE_ix;   //!
   TBranch        *b_PFRecHitEE_iy;   //!
   TBranch        *b_PFRecHitEE_eta;   //!
   TBranch        *b_PFRecHitEE_phi;   //!
   TBranch        *b_PFRecHitEE_cutThreshold;   //!
   TBranch        *b_PFRecHitEE_detId;   //!
   TBranch        *b_nSimHitHBHE;   //!
   TBranch        *b_SimHitHBHE_energy;   //!
   TBranch        *b_SimHitHBHE_energyEM;   //!
   TBranch        *b_SimHitHBHE_energyHAD;   //!
   TBranch        *b_SimHitHBHE_samplingFactor;   //!
   TBranch        *b_SimHitHBHE_nCaloHits;   //!
   TBranch        *b_SimHitHBHE_ieta;   //!
   TBranch        *b_SimHitHBHE_iphi;   //!
   TBranch        *b_SimHitHBHE_depth;   //!
   TBranch        *b_SimHitHBHE_eta;   //!
   TBranch        *b_SimHitHBHE_phi;   //!
   TBranch        *b_SimHitHBHE_detId;   //!
   TBranch        *b_SimHitHBHE_subdetId;   //!
   TBranch        *b_SimHitHBHE_energy_25ns;   //!
   TBranch        *b_SimHitHBHE_energyEM_25ns;   //!
   TBranch        *b_SimHitHBHE_energyHAD_25ns;   //!
   TBranch        *b_nSimHitEB;   //!
   TBranch        *b_SimHitEB_energy;   //!
   TBranch        *b_SimHitEB_energyEM;   //!
   TBranch        *b_SimHitEB_energyHAD;   //!
   TBranch        *b_SimHitEB_nCaloHits;   //!
   TBranch        *b_SimHitEB_ieta;   //!
   TBranch        *b_SimHitEB_iphi;   //!
   TBranch        *b_SimHitEB_eta;   //!
   TBranch        *b_SimHitEB_phi;   //!
   TBranch        *b_SimHitEB_detId;   //!
   TBranch        *b_SimHitEB_subdetId;   //!
   TBranch        *b_nSimHitEE;   //!
   TBranch        *b_SimHitEE_energy;   //!
   TBranch        *b_SimHitEE_energyEM;   //!
   TBranch        *b_SimHitEE_energyHAD;   //!
   TBranch        *b_SimHitEE_nCaloHits;   //!
   TBranch        *b_SimHitEE_ix;   //!
   TBranch        *b_SimHitEE_iy;   //!
   TBranch        *b_SimHitEE_eta;   //!
   TBranch        *b_SimHitEE_phi;   //!
   TBranch        *b_SimHitEE_detId;   //!
   TBranch        *b_SimHitEE_subdetId;   //!
   TBranch        *b_nPFClusterECAL;   //!
   TBranch        *b_PFClusterECAL_pt;   //!
   TBranch        *b_PFClusterECAL_energy;   //!
   TBranch        *b_PFClusterECAL_correctedEnergy;   //!
   TBranch        *b_PFClusterECAL_eta;   //!
   TBranch        *b_PFClusterECAL_phi;   //!
   TBranch        *b_PFClusterECAL_layer;   //!
   TBranch        *b_PFClusterECAL_seedhit_detId;   //!
   TBranch        *b_PFClusterECAL_nhits;   //!
   TBranch        *b_PFClusterECAL_hits_detId;   //!
   TBranch        *b_PFClusterECAL_hits_fraction;   //!
   TBranch        *b_PFClusterECAL_hits_PFRecHitEB_Idx;   //!
   TBranch        *b_PFClusterECAL_hits_PFRecHitEE_Idx;   //!
   TBranch        *b_PFClusterECAL_key;   //!
   TBranch        *b_nPFClusterPS;   //!
   TBranch        *b_PFClusterPS_pt;   //!
   TBranch        *b_PFClusterPS_energy;   //!
   TBranch        *b_PFClusterPS_correctedEnergy;   //!
   TBranch        *b_PFClusterPS_eta;   //!
   TBranch        *b_PFClusterPS_phi;   //!
   TBranch        *b_PFClusterPS_layer;   //!
   TBranch        *b_PFClusterPS_seedhit_detId;   //!
   TBranch        *b_PFClusterPS_nhits;   //!
   TBranch        *b_PFClusterPS_hits_detId;   //!
   TBranch        *b_PFClusterPS_hits_fraction;   //!
   TBranch        *b_PFClusterPS_key;   //!
   TBranch        *b_nPFClusterHCAL;   //!
   TBranch        *b_PFClusterHCAL_pt;   //!
   TBranch        *b_PFClusterHCAL_energy;   //!
   TBranch        *b_PFClusterHCAL_correctedEnergy;   //!
   TBranch        *b_PFClusterHCAL_eta;   //!
   TBranch        *b_PFClusterHCAL_phi;   //!
   TBranch        *b_PFClusterHCAL_layer;   //!
   TBranch        *b_PFClusterHCAL_seedhit_detId;   //!
   TBranch        *b_PFClusterHCAL_nhits;   //!
   TBranch        *b_PFClusterHCAL_hits_detId;   //!
   TBranch        *b_PFClusterHCAL_hits_fraction;   //!
   TBranch        *b_PFClusterHCAL_hits_PFRecHitHBHE_Idx;   //!
   TBranch        *b_PFClusterHCAL_key;   //!

   PFAna(TTree *tree=0);
   virtual ~PFAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PFAna_cxx
PFAna::PFAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/piongun/MergedNtuple_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/piongun/MergedNtuple_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../data/piongun/MergedNtuple_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root:/pfAnaNtuplizer");
      dir->GetObject("PFAnaTree",tree);

   }
   Init(tree);
}

PFAna::~PFAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PFAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PFAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PFAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PFCand_pt = 0;
   PFCand_eta = 0;
   PFCand_phi = 0;
   PFCand_mass = 0;
   PFCand_energy = 0;
   PFCand_pdgId = 0;
   PFCand_isChgHadIso = 0;
   PFCand_birthId = 0;
   PFCand_dRToGenPart = 0;
   PFCand_has_trk = 0;
   PFCand_trk_pt = 0;
   PFCand_trk_p = 0;
   PFCand_trk_eta = 0;
   PFCand_trk_phi = 0;
   PFCand_trk_algo = 0;
   PFCand_trk_charge = 0;
   PFCand_trk_ptError = 0;
   PFCand_trk_dxy = 0;
   PFCand_trk_dz = 0;
   PFCand_trk_dxyError = 0;
   PFCand_trk_dzError = 0;
   PFCand_trk_chi2 = 0;
   PFCand_trk_ndof = 0;
   PFCand_trk_vx = 0;
   PFCand_trk_vy = 0;
   PFCand_trk_vz = 0;
   PFCand_trk_numberOfValidPixelHits = 0;
   PFCand_trk_numberOfValidTrackerHits = 0;
   PFCand_trk_pixelHitOK = 0;
   PFCand_trk_trackerHitOK = 0;
   PFCand_trk_positionAtECALEntrance_x = 0;
   PFCand_trk_positionAtECALEntrance_y = 0;
   PFCand_trk_positionAtECALEntrance_z = 0;
   PFCand_trk_positionAtECALEntrance_eta = 0;
   PFCand_trk_positionAtECALEntrance_phi = 0;
   PFCand_trk_closestECAL_eta = 0;
   PFCand_trk_closestECAL_phi = 0;
   PFCand_trk_closestECAL_detId = 0;
   PFCand_trk_closestECAL_EBieta = 0;
   PFCand_trk_closestECAL_EBiphi = 0;
   PFCand_trk_closestECAL_EEix = 0;
   PFCand_trk_closestECAL_EEiy = 0;
   PFCand_trk_closestHCAL_eta = 0;
   PFCand_trk_closestHCAL_phi = 0;
   PFCand_trk_closestHCAL_detId = 0;
   PFCand_trk_closestHCAL_ieta = 0;
   PFCand_trk_closestHCAL_iphi = 0;
   PFCand_ecalEnergy = 0;
   PFCand_hcalEnergy = 0;
   PFCand_caloEnergy = 0;
   PFCand_rawEcalEnergy = 0;
   PFCand_rawHcalEnergy = 0;
   PFCand_rawCaloEnergy = 0;
   PFCand_hoEnergy = 0;
   PFCand_rawHoEnergy = 0;
   PFCand_pS1Energy = 0;
   PFCand_p21Energy = 0;
   PFCand_hcalDepth1EFrac = 0;
   PFCand_hcalDepth2EFrac = 0;
   PFCand_hcalDepth3EFrac = 0;
   PFCand_hcalDepth4EFrac = 0;
   PFCand_hcalDepth5EFrac = 0;
   PFCand_hcalDepth6EFrac = 0;
   PFCand_hcalDepth7EFrac = 0;
   PFCand_PFBlock_Idx = 0;
   PFCand_nPFTrack = 0;
   PFCand_nPFClusterHCAL = 0;
   PFCand_nPFClusterECAL = 0;
   PFCand_nPFClusterPS = 0;
   PFCand_nPFClusterHF = 0;
   PFCand_nPFClusterSC = 0;
   PFCand_nPFClusterGSF = 0;
   PFCand_nPFClusterBREM = 0;
   PFCand_PFClusterHCAL_Idx = 0;
   PFCand_PFClusterECAL_Idx = 0;
   PFCand_PFClusterPS_Idx = 0;
   PFCand_PFClusterHF_Idx = 0;
   PFCand_nPFTrackInBlock = 0;
   PFCand_nPFClusterHCALInBlock = 0;
   PFCand_nPFClusterECALInBlock = 0;
   PFCand_nPFClusterPSInBlock = 0;
   PFCand_nPFClusterHFInBlock = 0;
   PFCand_nPFClusterSCInBlock = 0;
   PFCand_nPFClusterGSFInBlock = 0;
   PFCand_nPFClusterBREMInBlock = 0;
   PFCand_PFClusterHCALInBlock_Idx = 0;
   PFCand_PFClusterECALInBlock_Idx = 0;
   PFCand_PFClusterPSInBlock_Idx = 0;
   PFCand_PFClusterHFInBlock_Idx = 0;
   PFRecHitHBHE_energy = 0;
   PFRecHitHBHE_ieta = 0;
   PFRecHitHBHE_iphi = 0;
   PFRecHitHBHE_eta = 0;
   PFRecHitHBHE_phi = 0;
   PFRecHitHBHE_depth = 0;
   PFRecHitHBHE_cutThreshold = 0;
   PFRecHitHBHE_hcalRespCorr = 0;
   PFRecHitHBHE_detId = 0;
   PFRecHitEB_energy = 0;
   PFRecHitEB_ieta = 0;
   PFRecHitEB_iphi = 0;
   PFRecHitEB_tower_ieta = 0;
   PFRecHitEB_tower_iphi = 0;
   PFRecHitEB_approxEta = 0;
   PFRecHitEB_eta = 0;
   PFRecHitEB_phi = 0;
   PFRecHitEB_cutThreshold = 0;
   PFRecHitEB_detId = 0;
   PFRecHitEE_energy = 0;
   PFRecHitEE_ix = 0;
   PFRecHitEE_iy = 0;
   PFRecHitEE_eta = 0;
   PFRecHitEE_phi = 0;
   PFRecHitEE_cutThreshold = 0;
   PFRecHitEE_detId = 0;
   SimHitHBHE_energy = 0;
   SimHitHBHE_energyEM = 0;
   SimHitHBHE_energyHAD = 0;
   SimHitHBHE_samplingFactor = 0;
   SimHitHBHE_nCaloHits = 0;
   SimHitHBHE_ieta = 0;
   SimHitHBHE_iphi = 0;
   SimHitHBHE_depth = 0;
   SimHitHBHE_eta = 0;
   SimHitHBHE_phi = 0;
   SimHitHBHE_detId = 0;
   SimHitHBHE_subdetId = 0;
   SimHitHBHE_energy_25ns = 0;
   SimHitHBHE_energyEM_25ns = 0;
   SimHitHBHE_energyHAD_25ns = 0;
   SimHitEB_energy = 0;
   SimHitEB_energyEM = 0;
   SimHitEB_energyHAD = 0;
   SimHitEB_nCaloHits = 0;
   SimHitEB_ieta = 0;
   SimHitEB_iphi = 0;
   SimHitEB_eta = 0;
   SimHitEB_phi = 0;
   SimHitEB_detId = 0;
   SimHitEB_subdetId = 0;
   SimHitEE_energy = 0;
   SimHitEE_energyEM = 0;
   SimHitEE_energyHAD = 0;
   SimHitEE_nCaloHits = 0;
   SimHitEE_ix = 0;
   SimHitEE_iy = 0;
   SimHitEE_eta = 0;
   SimHitEE_phi = 0;
   SimHitEE_detId = 0;
   SimHitEE_subdetId = 0;
   PFClusterECAL_pt = 0;
   PFClusterECAL_energy = 0;
   PFClusterECAL_correctedEnergy = 0;
   PFClusterECAL_eta = 0;
   PFClusterECAL_phi = 0;
   PFClusterECAL_layer = 0;
   PFClusterECAL_seedhit_detId = 0;
   PFClusterECAL_nhits = 0;
   PFClusterECAL_hits_detId = 0;
   PFClusterECAL_hits_fraction = 0;
   PFClusterECAL_hits_PFRecHitEB_Idx = 0;
   PFClusterECAL_hits_PFRecHitEE_Idx = 0;
   PFClusterECAL_key = 0;
   PFClusterPS_pt = 0;
   PFClusterPS_energy = 0;
   PFClusterPS_correctedEnergy = 0;
   PFClusterPS_eta = 0;
   PFClusterPS_phi = 0;
   PFClusterPS_layer = 0;
   PFClusterPS_seedhit_detId = 0;
   PFClusterPS_nhits = 0;
   PFClusterPS_hits_detId = 0;
   PFClusterPS_hits_fraction = 0;
   PFClusterPS_key = 0;
   PFClusterHCAL_pt = 0;
   PFClusterHCAL_energy = 0;
   PFClusterHCAL_correctedEnergy = 0;
   PFClusterHCAL_eta = 0;
   PFClusterHCAL_phi = 0;
   PFClusterHCAL_layer = 0;
   PFClusterHCAL_seedhit_detId = 0;
   PFClusterHCAL_nhits = 0;
   PFClusterHCAL_hits_detId = 0;
   PFClusterHCAL_hits_fraction = 0;
   PFClusterHCAL_hits_PFRecHitHBHE_Idx = 0;
   PFClusterHCAL_key = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("run", &run, &b_orun);
   //fChain->SetBranchAddress("evt", &evt, &b_orun);
   //fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_orun);
   //fChain->SetBranchAddress("time", &time, &b_orun);
   fChain->SetBranchAddress("genVertex_x", &genVertex_x, &b_genVertex_x);
   fChain->SetBranchAddress("genVertex_y", &genVertex_y, &b_genVertex_y);
   fChain->SetBranchAddress("genVertex_z", &genVertex_z, &b_genVertex_z);
   fChain->SetBranchAddress("genPart_pt", &genPart_pt, &b_genPart_pt);
   fChain->SetBranchAddress("genPart_eta", &genPart_eta, &b_genPart_eta);
   fChain->SetBranchAddress("genPart_phi", &genPart_phi, &b_genPart_phi);
   fChain->SetBranchAddress("genPart_mass", &genPart_mass, &b_genPart_mass);
   fChain->SetBranchAddress("genPart_energy", &genPart_energy, &b_genPart_energy);
   fChain->SetBranchAddress("genPart_p", &genPart_p, &b_genPart_p);
   fChain->SetBranchAddress("genPart_pdgId", &genPart_pdgId, &b_genPart_pdgId);
   fChain->SetBranchAddress("genPart_charge", &genPart_charge, &b_genPart_charge);
   fChain->SetBranchAddress("genPart_extrapolated_okECAL", &genPart_extrapolated_okECAL, &b_genPart_extrapolated_okECAL);
   fChain->SetBranchAddress("genPart_extrapolated_okHCAL", &genPart_extrapolated_okHCAL, &b_genPart_extrapolated_okHCAL);
   fChain->SetBranchAddress("genPart_extrapolated_EB_ieta", &genPart_extrapolated_EB_ieta, &b_genPart_extrapolated_EB_ieta);
   fChain->SetBranchAddress("genPart_extrapolated_EB_iphi", &genPart_extrapolated_EB_iphi, &b_genPart_extrapolated_EB_iphi);
   fChain->SetBranchAddress("genPart_extrapolated_EB_eta", &genPart_extrapolated_EB_eta, &b_genPart_extrapolated_EB_eta);
   fChain->SetBranchAddress("genPart_extrapolated_EB_phi", &genPart_extrapolated_EB_phi, &b_genPart_extrapolated_EB_phi);
   fChain->SetBranchAddress("genPart_extrapolated_PFRecHitEB_Idx", &genPart_extrapolated_PFRecHitEB_Idx, &b_genPart_extrapolated_PFRecHitEB_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_PFRecHitEB_detId", &genPart_extrapolated_PFRecHitEB_detId, &b_genPart_extrapolated_PFRecHitEB_detId);
   fChain->SetBranchAddress("genPart_extrapolated_EE_ix", &genPart_extrapolated_EE_ix, &b_genPart_extrapolated_EE_ix);
   fChain->SetBranchAddress("genPart_extrapolated_EE_iy", &genPart_extrapolated_EE_iy, &b_genPart_extrapolated_EE_iy);
   fChain->SetBranchAddress("genPart_extrapolated_EE_eta", &genPart_extrapolated_EE_eta, &b_genPart_extrapolated_EE_eta);
   fChain->SetBranchAddress("genPart_extrapolated_EE_phi", &genPart_extrapolated_EE_phi, &b_genPart_extrapolated_EE_phi);
   fChain->SetBranchAddress("genPart_extrapolated_PFRecHitEE_Idx", &genPart_extrapolated_PFRecHitEE_Idx, &b_genPart_extrapolated_PFRecHitEE_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_PFRecHitEE_detId", &genPart_extrapolated_PFRecHitEE_detId, &b_genPart_extrapolated_PFRecHitEE_detId);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_ieta", &genPart_extrapolated_EHCAL_ieta, &b_genPart_extrapolated_EHCAL_ieta);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_iphi", &genPart_extrapolated_EHCAL_iphi, &b_genPart_extrapolated_EHCAL_iphi);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_depth", &genPart_extrapolated_EHCAL_depth, &b_genPart_extrapolated_EHCAL_depth);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_eta", &genPart_extrapolated_EHCAL_eta, &b_genPart_extrapolated_EHCAL_eta);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_phi", &genPart_extrapolated_EHCAL_phi, &b_genPart_extrapolated_EHCAL_phi);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_subdetId", &genPart_extrapolated_EHCAL_subdetId, &b_genPart_extrapolated_EHCAL_subdetId);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx", &genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx, &b_genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_PFRecHitHBHE_detId", &genPart_extrapolated_EHCAL_PFRecHitHBHE_detId, &b_genPart_extrapolated_EHCAL_PFRecHitHBHE_detId);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_PFRecHitHF_Idx", &genPart_extrapolated_EHCAL_PFRecHitHF_Idx, &b_genPart_extrapolated_EHCAL_PFRecHitHF_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_EHCAL_PFRecHitHF_detId", &genPart_extrapolated_EHCAL_PFRecHitHF_detId, &b_genPart_extrapolated_EHCAL_PFRecHitHF_detId);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_ieta", &genPart_extrapolated_HCAL_ieta, &b_genPart_extrapolated_HCAL_ieta);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_iphi", &genPart_extrapolated_HCAL_iphi, &b_genPart_extrapolated_HCAL_iphi);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_depth", &genPart_extrapolated_HCAL_depth, &b_genPart_extrapolated_HCAL_depth);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_eta", &genPart_extrapolated_HCAL_eta, &b_genPart_extrapolated_HCAL_eta);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_phi", &genPart_extrapolated_HCAL_phi, &b_genPart_extrapolated_HCAL_phi);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_subdetId", &genPart_extrapolated_HCAL_subdetId, &b_genPart_extrapolated_HCAL_subdetId);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_PFRecHitHBHE_Idx", &genPart_extrapolated_HCAL_PFRecHitHBHE_Idx, &b_genPart_extrapolated_HCAL_PFRecHitHBHE_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_PFRecHitHBHE_detId", &genPart_extrapolated_HCAL_PFRecHitHBHE_detId, &b_genPart_extrapolated_HCAL_PFRecHitHBHE_detId);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_PFRecHitHF_Idx", &genPart_extrapolated_HCAL_PFRecHitHF_Idx, &b_genPart_extrapolated_HCAL_PFRecHitHF_Idx);
   fChain->SetBranchAddress("genPart_extrapolated_HCAL_PFRecHitHF_detId", &genPart_extrapolated_HCAL_PFRecHitHF_detId, &b_genPart_extrapolated_HCAL_PFRecHitHF_detId);
   fChain->SetBranchAddress("genPart_extrapolated_pointECAL_x", &genPart_extrapolated_pointECAL_x, &b_genPart_extrapolated_pointECAL_x);
   fChain->SetBranchAddress("genPart_extrapolated_pointECAL_y", &genPart_extrapolated_pointECAL_y, &b_genPart_extrapolated_pointECAL_y);
   fChain->SetBranchAddress("genPart_extrapolated_pointECAL_z", &genPart_extrapolated_pointECAL_z, &b_genPart_extrapolated_pointECAL_z);
   fChain->SetBranchAddress("genPart_extrapolated_pointECAL_eta", &genPart_extrapolated_pointECAL_eta, &b_genPart_extrapolated_pointECAL_eta);
   fChain->SetBranchAddress("genPart_extrapolated_pointECAL_phi", &genPart_extrapolated_pointECAL_phi, &b_genPart_extrapolated_pointECAL_phi);
   fChain->SetBranchAddress("genPart_extrapolated_pointHCAL_x", &genPart_extrapolated_pointHCAL_x, &b_genPart_extrapolated_pointHCAL_x);
   fChain->SetBranchAddress("genPart_extrapolated_pointHCAL_y", &genPart_extrapolated_pointHCAL_y, &b_genPart_extrapolated_pointHCAL_y);
   fChain->SetBranchAddress("genPart_extrapolated_pointHCAL_z", &genPart_extrapolated_pointHCAL_z, &b_genPart_extrapolated_pointHCAL_z);
   fChain->SetBranchAddress("genPart_extrapolated_pointHCAL_eta", &genPart_extrapolated_pointHCAL_eta, &b_genPart_extrapolated_pointHCAL_eta);
   fChain->SetBranchAddress("genPart_extrapolated_pointHCAL_phi", &genPart_extrapolated_pointHCAL_phi, &b_genPart_extrapolated_pointHCAL_phi);
   fChain->SetBranchAddress("nPFCand", &nPFCand, &b_nPFCand);
   fChain->SetBranchAddress("PFCand_pt", &PFCand_pt, &b_PFCand_pt);
   fChain->SetBranchAddress("PFCand_eta", &PFCand_eta, &b_PFCand_eta);
   fChain->SetBranchAddress("PFCand_phi", &PFCand_phi, &b_PFCand_phi);
   fChain->SetBranchAddress("PFCand_mass", &PFCand_mass, &b_PFCand_mass);
   fChain->SetBranchAddress("PFCand_energy", &PFCand_energy, &b_PFCand_energy);
   fChain->SetBranchAddress("PFCand_pdgId", &PFCand_pdgId, &b_PFCand_pdgId);
   fChain->SetBranchAddress("PFCand_isChgHadIso", &PFCand_isChgHadIso, &b_PFCand_isChgHadIso);
   fChain->SetBranchAddress("PFCand_birthId", &PFCand_birthId, &b_PFCand_birthId);
   fChain->SetBranchAddress("PFCand_dRToGenPart", &PFCand_dRToGenPart, &b_PFCand_dRToGenPart);
   fChain->SetBranchAddress("PFCand_has_trk", &PFCand_has_trk, &b_PFCand_has_trk);
   fChain->SetBranchAddress("PFCand_trk_pt", &PFCand_trk_pt, &b_PFCand_trk_pt);
   fChain->SetBranchAddress("PFCand_trk_p", &PFCand_trk_p, &b_PFCand_trk_p);
   fChain->SetBranchAddress("PFCand_trk_eta", &PFCand_trk_eta, &b_PFCand_trk_eta);
   fChain->SetBranchAddress("PFCand_trk_phi", &PFCand_trk_phi, &b_PFCand_trk_phi);
   fChain->SetBranchAddress("PFCand_trk_algo", &PFCand_trk_algo, &b_PFCand_trk_algo);
   fChain->SetBranchAddress("PFCand_trk_charge", &PFCand_trk_charge, &b_PFCand_trk_charge);
   fChain->SetBranchAddress("PFCand_trk_ptError", &PFCand_trk_ptError, &b_PFCand_trk_ptError);
   fChain->SetBranchAddress("PFCand_trk_dxy", &PFCand_trk_dxy, &b_PFCand_trk_dxy);
   fChain->SetBranchAddress("PFCand_trk_dz", &PFCand_trk_dz, &b_PFCand_trk_dz);
   fChain->SetBranchAddress("PFCand_trk_dxyError", &PFCand_trk_dxyError, &b_PFCand_trk_dxyError);
   fChain->SetBranchAddress("PFCand_trk_dzError", &PFCand_trk_dzError, &b_PFCand_trk_dzError);
   fChain->SetBranchAddress("PFCand_trk_chi2", &PFCand_trk_chi2, &b_PFCand_trk_chi2);
   fChain->SetBranchAddress("PFCand_trk_ndof", &PFCand_trk_ndof, &b_PFCand_trk_ndof);
   fChain->SetBranchAddress("PFCand_trk_vx", &PFCand_trk_vx, &b_PFCand_trk_vx);
   fChain->SetBranchAddress("PFCand_trk_vy", &PFCand_trk_vy, &b_PFCand_trk_vy);
   fChain->SetBranchAddress("PFCand_trk_vz", &PFCand_trk_vz, &b_PFCand_trk_vz);
   fChain->SetBranchAddress("PFCand_trk_numberOfValidPixelHits", &PFCand_trk_numberOfValidPixelHits, &b_PFCand_trk_numberOfValidPixelHits);
   fChain->SetBranchAddress("PFCand_trk_numberOfValidTrackerHits", &PFCand_trk_numberOfValidTrackerHits, &b_PFCand_trk_numberOfValidTrackerHits);
   fChain->SetBranchAddress("PFCand_trk_pixelHitOK", &PFCand_trk_pixelHitOK, &b_PFCand_trk_pixelHitOK);
   fChain->SetBranchAddress("PFCand_trk_trackerHitOK", &PFCand_trk_trackerHitOK, &b_PFCand_trk_trackerHitOK);
   fChain->SetBranchAddress("PFCand_trk_positionAtECALEntrance_x", &PFCand_trk_positionAtECALEntrance_x, &b_PFCand_trk_positionAtECALEntrance_x);
   fChain->SetBranchAddress("PFCand_trk_positionAtECALEntrance_y", &PFCand_trk_positionAtECALEntrance_y, &b_PFCand_trk_positionAtECALEntrance_y);
   fChain->SetBranchAddress("PFCand_trk_positionAtECALEntrance_z", &PFCand_trk_positionAtECALEntrance_z, &b_PFCand_trk_positionAtECALEntrance_z);
   fChain->SetBranchAddress("PFCand_trk_positionAtECALEntrance_eta", &PFCand_trk_positionAtECALEntrance_eta, &b_PFCand_trk_positionAtECALEntrance_eta);
   fChain->SetBranchAddress("PFCand_trk_positionAtECALEntrance_phi", &PFCand_trk_positionAtECALEntrance_phi, &b_PFCand_trk_positionAtECALEntrance_phi);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_eta", &PFCand_trk_closestECAL_eta, &b_PFCand_trk_closestECAL_eta);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_phi", &PFCand_trk_closestECAL_phi, &b_PFCand_trk_closestECAL_phi);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_detId", &PFCand_trk_closestECAL_detId, &b_PFCand_trk_closestECAL_detId);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_EBieta", &PFCand_trk_closestECAL_EBieta, &b_PFCand_trk_closestECAL_EBieta);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_EBiphi", &PFCand_trk_closestECAL_EBiphi, &b_PFCand_trk_closestECAL_EBiphi);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_EEix", &PFCand_trk_closestECAL_EEix, &b_PFCand_trk_closestECAL_EEix);
   fChain->SetBranchAddress("PFCand_trk_closestECAL_EEiy", &PFCand_trk_closestECAL_EEiy, &b_PFCand_trk_closestECAL_EEiy);
   fChain->SetBranchAddress("PFCand_trk_closestHCAL_eta", &PFCand_trk_closestHCAL_eta, &b_PFCand_trk_closestHCAL_eta);
   fChain->SetBranchAddress("PFCand_trk_closestHCAL_phi", &PFCand_trk_closestHCAL_phi, &b_PFCand_trk_closestHCAL_phi);
   fChain->SetBranchAddress("PFCand_trk_closestHCAL_detId", &PFCand_trk_closestHCAL_detId, &b_PFCand_trk_closestHCAL_detId);
   fChain->SetBranchAddress("PFCand_trk_closestHCAL_ieta", &PFCand_trk_closestHCAL_ieta, &b_PFCand_trk_closestHCAL_ieta);
   fChain->SetBranchAddress("PFCand_trk_closestHCAL_iphi", &PFCand_trk_closestHCAL_iphi, &b_PFCand_trk_closestHCAL_iphi);
   fChain->SetBranchAddress("PFCand_ecalEnergy", &PFCand_ecalEnergy, &b_PFCand_ecalEnergy);
   fChain->SetBranchAddress("PFCand_hcalEnergy", &PFCand_hcalEnergy, &b_PFCand_hcalEnergy);
   fChain->SetBranchAddress("PFCand_caloEnergy", &PFCand_caloEnergy, &b_PFCand_caloEnergy);
   fChain->SetBranchAddress("PFCand_rawEcalEnergy", &PFCand_rawEcalEnergy, &b_PFCand_rawEcalEnergy);
   fChain->SetBranchAddress("PFCand_rawHcalEnergy", &PFCand_rawHcalEnergy, &b_PFCand_rawHcalEnergy);
   fChain->SetBranchAddress("PFCand_rawCaloEnergy", &PFCand_rawCaloEnergy, &b_PFCand_rawCaloEnergy);
   fChain->SetBranchAddress("PFCand_hoEnergy", &PFCand_hoEnergy, &b_PFCand_hoEnergy);
   fChain->SetBranchAddress("PFCand_rawHoEnergy", &PFCand_rawHoEnergy, &b_PFCand_rawHoEnergy);
   fChain->SetBranchAddress("PFCand_pS1Energy", &PFCand_pS1Energy, &b_PFCand_pS1Energy);
   fChain->SetBranchAddress("PFCand_p21Energy", &PFCand_p21Energy, &b_PFCand_p21Energy);
   fChain->SetBranchAddress("PFCand_hcalDepth1EFrac", &PFCand_hcalDepth1EFrac, &b_PFCand_hcalDepth1EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth2EFrac", &PFCand_hcalDepth2EFrac, &b_PFCand_hcalDepth2EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth3EFrac", &PFCand_hcalDepth3EFrac, &b_PFCand_hcalDepth3EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth4EFrac", &PFCand_hcalDepth4EFrac, &b_PFCand_hcalDepth4EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth5EFrac", &PFCand_hcalDepth5EFrac, &b_PFCand_hcalDepth5EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth6EFrac", &PFCand_hcalDepth6EFrac, &b_PFCand_hcalDepth6EFrac);
   fChain->SetBranchAddress("PFCand_hcalDepth7EFrac", &PFCand_hcalDepth7EFrac, &b_PFCand_hcalDepth7EFrac);
   fChain->SetBranchAddress("PFCand_PFBlock_Idx", &PFCand_PFBlock_Idx, &b_PFCand_PFBlock_Idx);
   fChain->SetBranchAddress("PFCand_nPFTrack", &PFCand_nPFTrack, &b_PFCand_nPFTrack);
   fChain->SetBranchAddress("PFCand_nPFClusterHCAL", &PFCand_nPFClusterHCAL, &b_PFCand_nPFClusterHCAL);
   fChain->SetBranchAddress("PFCand_nPFClusterECAL", &PFCand_nPFClusterECAL, &b_PFCand_nPFClusterECAL);
   fChain->SetBranchAddress("PFCand_nPFClusterPS", &PFCand_nPFClusterPS, &b_PFCand_nPFClusterPS);
   fChain->SetBranchAddress("PFCand_nPFClusterHF", &PFCand_nPFClusterHF, &b_PFCand_nPFClusterHF);
   fChain->SetBranchAddress("PFCand_nPFClusterSC", &PFCand_nPFClusterSC, &b_PFCand_nPFClusterSC);
   fChain->SetBranchAddress("PFCand_nPFClusterGSF", &PFCand_nPFClusterGSF, &b_PFCand_nPFClusterGSF);
   fChain->SetBranchAddress("PFCand_nPFClusterBREM", &PFCand_nPFClusterBREM, &b_PFCand_nPFClusterBREM);
   fChain->SetBranchAddress("PFCand_PFClusterHCAL_Idx", &PFCand_PFClusterHCAL_Idx, &b_PFCand_PFClusterHCAL_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterECAL_Idx", &PFCand_PFClusterECAL_Idx, &b_PFCand_PFClusterECAL_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterPS_Idx", &PFCand_PFClusterPS_Idx, &b_PFCand_PFClusterPS_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterHF_Idx", &PFCand_PFClusterHF_Idx, &b_PFCand_PFClusterHF_Idx);
   fChain->SetBranchAddress("PFCand_nPFTrackInBlock", &PFCand_nPFTrackInBlock, &b_PFCand_nPFTrackInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterHCALInBlock", &PFCand_nPFClusterHCALInBlock, &b_PFCand_nPFClusterHCALInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterECALInBlock", &PFCand_nPFClusterECALInBlock, &b_PFCand_nPFClusterECALInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterPSInBlock", &PFCand_nPFClusterPSInBlock, &b_PFCand_nPFClusterPSInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterHFInBlock", &PFCand_nPFClusterHFInBlock, &b_PFCand_nPFClusterHFInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterSCInBlock", &PFCand_nPFClusterSCInBlock, &b_PFCand_nPFClusterSCInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterGSFInBlock", &PFCand_nPFClusterGSFInBlock, &b_PFCand_nPFClusterGSFInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterBREMInBlock", &PFCand_nPFClusterBREMInBlock, &b_PFCand_nPFClusterBREMInBlock);
   fChain->SetBranchAddress("PFCand_PFClusterHCALInBlock_Idx", &PFCand_PFClusterHCALInBlock_Idx, &b_PFCand_PFClusterHCALInBlock_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterECALInBlock_Idx", &PFCand_PFClusterECALInBlock_Idx, &b_PFCand_PFClusterECALInBlock_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterPSInBlock_Idx", &PFCand_PFClusterPSInBlock_Idx, &b_PFCand_PFClusterPSInBlock_Idx);
   fChain->SetBranchAddress("PFCand_PFClusterHFInBlock_Idx", &PFCand_PFClusterHFInBlock_Idx, &b_PFCand_PFClusterHFInBlock_Idx);
   fChain->SetBranchAddress("Idx_ClosestPFCandChgHad", &Idx_ClosestPFCandChgHad, &b_Idx_ClosestPFCandChgHad);
   fChain->SetBranchAddress("Idx_ClosestPFCandChgHadIso", &Idx_ClosestPFCandChgHadIso, &b_Idx_ClosestPFCandChgHadIso);
   fChain->SetBranchAddress("Idx_ClosestPFCandNeuHad", &Idx_ClosestPFCandNeuHad, &b_Idx_ClosestPFCandNeuHad);
   fChain->SetBranchAddress("Idx_ClosestPFCandPhoton", &Idx_ClosestPFCandPhoton, &b_Idx_ClosestPFCandPhoton);
   fChain->SetBranchAddress("nPFRecHitHBHE", &nPFRecHitHBHE, &b_nPFRecHitHBHE);
   fChain->SetBranchAddress("PFRecHitHBHE_energy", &PFRecHitHBHE_energy, &b_PFRecHitHBHE_energy);
   fChain->SetBranchAddress("PFRecHitHBHE_ieta", &PFRecHitHBHE_ieta, &b_PFRecHitHBHE_ieta);
   fChain->SetBranchAddress("PFRecHitHBHE_iphi", &PFRecHitHBHE_iphi, &b_PFRecHitHBHE_iphi);
   fChain->SetBranchAddress("PFRecHitHBHE_eta", &PFRecHitHBHE_eta, &b_PFRecHitHBHE_eta);
   fChain->SetBranchAddress("PFRecHitHBHE_phi", &PFRecHitHBHE_phi, &b_PFRecHitHBHE_phi);
   fChain->SetBranchAddress("PFRecHitHBHE_depth", &PFRecHitHBHE_depth, &b_PFRecHitHBHE_depth);
   fChain->SetBranchAddress("PFRecHitHBHE_cutThreshold", &PFRecHitHBHE_cutThreshold, &b_PFRecHitHBHE_cutThreshold);
   fChain->SetBranchAddress("PFRecHitHBHE_hcalRespCorr", &PFRecHitHBHE_hcalRespCorr, &b_PFRecHitHBHE_hcalRespCorr);
   fChain->SetBranchAddress("PFRecHitHBHE_detId", &PFRecHitHBHE_detId, &b_PFRecHitHBHE_detId);
   fChain->SetBranchAddress("nPFRecHitEB", &nPFRecHitEB, &b_nPFRecHitEB);
   fChain->SetBranchAddress("PFRecHitEB_energy", &PFRecHitEB_energy, &b_PFRecHitEB_energy);
   fChain->SetBranchAddress("PFRecHitEB_ieta", &PFRecHitEB_ieta, &b_PFRecHitEB_ieta);
   fChain->SetBranchAddress("PFRecHitEB_iphi", &PFRecHitEB_iphi, &b_PFRecHitEB_iphi);
   fChain->SetBranchAddress("PFRecHitEB_tower_ieta", &PFRecHitEB_tower_ieta, &b_PFRecHitEB_tower_ieta);
   fChain->SetBranchAddress("PFRecHitEB_tower_iphi", &PFRecHitEB_tower_iphi, &b_PFRecHitEB_tower_iphi);
   fChain->SetBranchAddress("PFRecHitEB_approxEta", &PFRecHitEB_approxEta, &b_PFRecHitEB_approxEta);
   fChain->SetBranchAddress("PFRecHitEB_eta", &PFRecHitEB_eta, &b_PFRecHitEB_eta);
   fChain->SetBranchAddress("PFRecHitEB_phi", &PFRecHitEB_phi, &b_PFRecHitEB_phi);
   fChain->SetBranchAddress("PFRecHitEB_cutThreshold", &PFRecHitEB_cutThreshold, &b_PFRecHitEB_cutThreshold);
   fChain->SetBranchAddress("PFRecHitEB_detId", &PFRecHitEB_detId, &b_PFRecHitEB_detId);
   fChain->SetBranchAddress("nPFRecHitEE", &nPFRecHitEE, &b_nPFRecHitEE);
   fChain->SetBranchAddress("PFRecHitEE_energy", &PFRecHitEE_energy, &b_PFRecHitEE_energy);
   fChain->SetBranchAddress("PFRecHitEE_ix", &PFRecHitEE_ix, &b_PFRecHitEE_ix);
   fChain->SetBranchAddress("PFRecHitEE_iy", &PFRecHitEE_iy, &b_PFRecHitEE_iy);
   fChain->SetBranchAddress("PFRecHitEE_eta", &PFRecHitEE_eta, &b_PFRecHitEE_eta);
   fChain->SetBranchAddress("PFRecHitEE_phi", &PFRecHitEE_phi, &b_PFRecHitEE_phi);
   fChain->SetBranchAddress("PFRecHitEE_cutThreshold", &PFRecHitEE_cutThreshold, &b_PFRecHitEE_cutThreshold);
   fChain->SetBranchAddress("PFRecHitEE_detId", &PFRecHitEE_detId, &b_PFRecHitEE_detId);
   fChain->SetBranchAddress("nSimHitHBHE", &nSimHitHBHE, &b_nSimHitHBHE);
   fChain->SetBranchAddress("SimHitHBHE_energy", &SimHitHBHE_energy, &b_SimHitHBHE_energy);
   fChain->SetBranchAddress("SimHitHBHE_energyEM", &SimHitHBHE_energyEM, &b_SimHitHBHE_energyEM);
   fChain->SetBranchAddress("SimHitHBHE_energyHAD", &SimHitHBHE_energyHAD, &b_SimHitHBHE_energyHAD);
   fChain->SetBranchAddress("SimHitHBHE_samplingFactor", &SimHitHBHE_samplingFactor, &b_SimHitHBHE_samplingFactor);
   fChain->SetBranchAddress("SimHitHBHE_nCaloHits", &SimHitHBHE_nCaloHits, &b_SimHitHBHE_nCaloHits);
   fChain->SetBranchAddress("SimHitHBHE_ieta", &SimHitHBHE_ieta, &b_SimHitHBHE_ieta);
   fChain->SetBranchAddress("SimHitHBHE_iphi", &SimHitHBHE_iphi, &b_SimHitHBHE_iphi);
   fChain->SetBranchAddress("SimHitHBHE_depth", &SimHitHBHE_depth, &b_SimHitHBHE_depth);
   fChain->SetBranchAddress("SimHitHBHE_eta", &SimHitHBHE_eta, &b_SimHitHBHE_eta);
   fChain->SetBranchAddress("SimHitHBHE_phi", &SimHitHBHE_phi, &b_SimHitHBHE_phi);
   fChain->SetBranchAddress("SimHitHBHE_detId", &SimHitHBHE_detId, &b_SimHitHBHE_detId);
   fChain->SetBranchAddress("SimHitHBHE_subdetId", &SimHitHBHE_subdetId, &b_SimHitHBHE_subdetId);
   fChain->SetBranchAddress("SimHitHBHE_energy_25ns", &SimHitHBHE_energy_25ns, &b_SimHitHBHE_energy_25ns);
   fChain->SetBranchAddress("SimHitHBHE_energyEM_25ns", &SimHitHBHE_energyEM_25ns, &b_SimHitHBHE_energyEM_25ns);
   fChain->SetBranchAddress("SimHitHBHE_energyHAD_25ns", &SimHitHBHE_energyHAD_25ns, &b_SimHitHBHE_energyHAD_25ns);
   fChain->SetBranchAddress("nSimHitEB", &nSimHitEB, &b_nSimHitEB);
   fChain->SetBranchAddress("SimHitEB_energy", &SimHitEB_energy, &b_SimHitEB_energy);
   fChain->SetBranchAddress("SimHitEB_energyEM", &SimHitEB_energyEM, &b_SimHitEB_energyEM);
   fChain->SetBranchAddress("SimHitEB_energyHAD", &SimHitEB_energyHAD, &b_SimHitEB_energyHAD);
   fChain->SetBranchAddress("SimHitEB_nCaloHits", &SimHitEB_nCaloHits, &b_SimHitEB_nCaloHits);
   fChain->SetBranchAddress("SimHitEB_ieta", &SimHitEB_ieta, &b_SimHitEB_ieta);
   fChain->SetBranchAddress("SimHitEB_iphi", &SimHitEB_iphi, &b_SimHitEB_iphi);
   fChain->SetBranchAddress("SimHitEB_eta", &SimHitEB_eta, &b_SimHitEB_eta);
   fChain->SetBranchAddress("SimHitEB_phi", &SimHitEB_phi, &b_SimHitEB_phi);
   fChain->SetBranchAddress("SimHitEB_detId", &SimHitEB_detId, &b_SimHitEB_detId);
   fChain->SetBranchAddress("SimHitEB_subdetId", &SimHitEB_subdetId, &b_SimHitEB_subdetId);
   fChain->SetBranchAddress("nSimHitEE", &nSimHitEE, &b_nSimHitEE);
   fChain->SetBranchAddress("SimHitEE_energy", &SimHitEE_energy, &b_SimHitEE_energy);
   fChain->SetBranchAddress("SimHitEE_energyEM", &SimHitEE_energyEM, &b_SimHitEE_energyEM);
   fChain->SetBranchAddress("SimHitEE_energyHAD", &SimHitEE_energyHAD, &b_SimHitEE_energyHAD);
   fChain->SetBranchAddress("SimHitEE_nCaloHits", &SimHitEE_nCaloHits, &b_SimHitEE_nCaloHits);
   fChain->SetBranchAddress("SimHitEE_ix", &SimHitEE_ix, &b_SimHitEE_ix);
   fChain->SetBranchAddress("SimHitEE_iy", &SimHitEE_iy, &b_SimHitEE_iy);
   fChain->SetBranchAddress("SimHitEE_eta", &SimHitEE_eta, &b_SimHitEE_eta);
   fChain->SetBranchAddress("SimHitEE_phi", &SimHitEE_phi, &b_SimHitEE_phi);
   fChain->SetBranchAddress("SimHitEE_detId", &SimHitEE_detId, &b_SimHitEE_detId);
   fChain->SetBranchAddress("SimHitEE_subdetId", &SimHitEE_subdetId, &b_SimHitEE_subdetId);
   fChain->SetBranchAddress("nPFClusterECAL", &nPFClusterECAL, &b_nPFClusterECAL);
   fChain->SetBranchAddress("PFClusterECAL_pt", &PFClusterECAL_pt, &b_PFClusterECAL_pt);
   fChain->SetBranchAddress("PFClusterECAL_energy", &PFClusterECAL_energy, &b_PFClusterECAL_energy);
   fChain->SetBranchAddress("PFClusterECAL_correctedEnergy", &PFClusterECAL_correctedEnergy, &b_PFClusterECAL_correctedEnergy);
   fChain->SetBranchAddress("PFClusterECAL_eta", &PFClusterECAL_eta, &b_PFClusterECAL_eta);
   fChain->SetBranchAddress("PFClusterECAL_phi", &PFClusterECAL_phi, &b_PFClusterECAL_phi);
   fChain->SetBranchAddress("PFClusterECAL_layer", &PFClusterECAL_layer, &b_PFClusterECAL_layer);
   fChain->SetBranchAddress("PFClusterECAL_seedhit_detId", &PFClusterECAL_seedhit_detId, &b_PFClusterECAL_seedhit_detId);
   fChain->SetBranchAddress("PFClusterECAL_nhits", &PFClusterECAL_nhits, &b_PFClusterECAL_nhits);
   fChain->SetBranchAddress("PFClusterECAL_hits_detId", &PFClusterECAL_hits_detId, &b_PFClusterECAL_hits_detId);
   fChain->SetBranchAddress("PFClusterECAL_hits_fraction", &PFClusterECAL_hits_fraction, &b_PFClusterECAL_hits_fraction);
   fChain->SetBranchAddress("PFClusterECAL_hits_PFRecHitEB_Idx", &PFClusterECAL_hits_PFRecHitEB_Idx, &b_PFClusterECAL_hits_PFRecHitEB_Idx);
   fChain->SetBranchAddress("PFClusterECAL_hits_PFRecHitEE_Idx", &PFClusterECAL_hits_PFRecHitEE_Idx, &b_PFClusterECAL_hits_PFRecHitEE_Idx);
   fChain->SetBranchAddress("PFClusterECAL_key", &PFClusterECAL_key, &b_PFClusterECAL_key);
   fChain->SetBranchAddress("nPFClusterPS", &nPFClusterPS, &b_nPFClusterPS);
   fChain->SetBranchAddress("PFClusterPS_pt", &PFClusterPS_pt, &b_PFClusterPS_pt);
   fChain->SetBranchAddress("PFClusterPS_energy", &PFClusterPS_energy, &b_PFClusterPS_energy);
   fChain->SetBranchAddress("PFClusterPS_correctedEnergy", &PFClusterPS_correctedEnergy, &b_PFClusterPS_correctedEnergy);
   fChain->SetBranchAddress("PFClusterPS_eta", &PFClusterPS_eta, &b_PFClusterPS_eta);
   fChain->SetBranchAddress("PFClusterPS_phi", &PFClusterPS_phi, &b_PFClusterPS_phi);
   fChain->SetBranchAddress("PFClusterPS_layer", &PFClusterPS_layer, &b_PFClusterPS_layer);
   fChain->SetBranchAddress("PFClusterPS_seedhit_detId", &PFClusterPS_seedhit_detId, &b_PFClusterPS_seedhit_detId);
   fChain->SetBranchAddress("PFClusterPS_nhits", &PFClusterPS_nhits, &b_PFClusterPS_nhits);
   fChain->SetBranchAddress("PFClusterPS_hits_detId", &PFClusterPS_hits_detId, &b_PFClusterPS_hits_detId);
   fChain->SetBranchAddress("PFClusterPS_hits_fraction", &PFClusterPS_hits_fraction, &b_PFClusterPS_hits_fraction);
   fChain->SetBranchAddress("PFClusterPS_key", &PFClusterPS_key, &b_PFClusterPS_key);
   fChain->SetBranchAddress("nPFClusterHCAL", &nPFClusterHCAL, &b_nPFClusterHCAL);
   fChain->SetBranchAddress("PFClusterHCAL_pt", &PFClusterHCAL_pt, &b_PFClusterHCAL_pt);
   fChain->SetBranchAddress("PFClusterHCAL_energy", &PFClusterHCAL_energy, &b_PFClusterHCAL_energy);
   fChain->SetBranchAddress("PFClusterHCAL_correctedEnergy", &PFClusterHCAL_correctedEnergy, &b_PFClusterHCAL_correctedEnergy);
   fChain->SetBranchAddress("PFClusterHCAL_eta", &PFClusterHCAL_eta, &b_PFClusterHCAL_eta);
   fChain->SetBranchAddress("PFClusterHCAL_phi", &PFClusterHCAL_phi, &b_PFClusterHCAL_phi);
   fChain->SetBranchAddress("PFClusterHCAL_layer", &PFClusterHCAL_layer, &b_PFClusterHCAL_layer);
   fChain->SetBranchAddress("PFClusterHCAL_seedhit_detId", &PFClusterHCAL_seedhit_detId, &b_PFClusterHCAL_seedhit_detId);
   fChain->SetBranchAddress("PFClusterHCAL_nhits", &PFClusterHCAL_nhits, &b_PFClusterHCAL_nhits);
   fChain->SetBranchAddress("PFClusterHCAL_hits_detId", &PFClusterHCAL_hits_detId, &b_PFClusterHCAL_hits_detId);
   fChain->SetBranchAddress("PFClusterHCAL_hits_fraction", &PFClusterHCAL_hits_fraction, &b_PFClusterHCAL_hits_fraction);
   fChain->SetBranchAddress("PFClusterHCAL_hits_PFRecHitHBHE_Idx", &PFClusterHCAL_hits_PFRecHitHBHE_Idx, &b_PFClusterHCAL_hits_PFRecHitHBHE_Idx);
   fChain->SetBranchAddress("PFClusterHCAL_key", &PFClusterHCAL_key, &b_PFClusterHCAL_key);
   Notify();
}

Bool_t PFAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PFAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PFAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PFAna_cxx
