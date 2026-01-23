// Purpose: derive PF hadron calibration (PFHC) using powerlaw approach
//          parameterize vs <rawEcal+rawHcal> and rawEcal/(rawEcal+rawHcal)
//          |eta| binning to follow that of JEC L2Relative

#define piongun_cxx
#include "piongun.h"
#include "PFEnergyCalibrationFromMikko.h"
#include "PFEnergyCalibrationFromMikko.cc"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH3D.h"

#include <iostream>

// Filter maximum E to compare merged sample to old
bool filterMaxP = false;
double maxP = 500;//5000;

// Filter with p>0 and/or charge!=0
// Without this |eta|<2.322 high E end seems bad
// (Nuclear interactions? With cut efficiency at high E is low)
// Trye turning both off for 2026v2
bool filterP = false;//true;//false;//v1
bool filterC = true;//false;//true;

// Subtract Random Cone (necessary for withPU, optional for noPU)
bool subRC = false;//true;

// Testing corrections from Conrado stored in tuple
bool usePFHC = false;//true;//false;
bool usePFEC = false;

// Testing corrections in PFEnergyCalibrationFromMikko.cc + piongun*.txt
bool applyPFEC_Charged = false;        // use charged/true energy
bool applyPFEC_Neutral = false;//true;//false; // use calorimeter energy

// Make plost for depth1 studies
bool doDepth1 = true;

void piongun::Loop()
{
//   In a ROOT session, you can do:
//      root> .L piongun.C
//      root> piongun t
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

  PFEnergyCalibration *pfec = new PFEnergyCalibration();

  // Sanity check test
  double etestraw(1.), htestraw(1.);
  double ttest(1000.), etestcorr(etestraw), htestcorr(htestraw), etatest(1.7);
  pfec->energyEmHad(ttest, etestcorr, htestcorr, etatest, 0.);
  double corrtest = (etestcorr+htestcorr) / (etestraw + htestraw);
  cout << "piongun.C:" << endl;
  cout << " ttest = " << ttest
       << " etestraw = " << etestraw
       << " htestraw = " << htestraw
       << " etatest = " << etatest << endl
       << " etestcorr = " << etestcorr
       << " htestcorr = " << htestcorr
       << " corrtest = " << corrtest
       << endl << flush;
  
  fChain->SetBranchStatus("*",0);
  //fChain->SetBranchStatus("true",1); // v1, 20241031 // 2025
  fChain->SetBranchStatus("genP",1); // v2 // 2025_v2
  //fChain->SetBranchStatus("eta",1); // v1, 20241031 // 2025
  fChain->SetBranchStatus("genEta",1); // v2 // 2025_v2, closure
  //fChain->SetBranchStatus("rawEcal",1);
  fChain->SetBranchStatus("ecal",1);
  //fChain->SetBranchStatus("rawHcal",1);
  fChain->SetBranchStatus("hcal",1);
  if (subRC) fChain->SetBranchStatus("rcEcal",1);
  if (subRC) fChain->SetBranchStatus("rcHcal",1);
  //fChain->SetBranchStatus("ho",1);
  //if (filterP) fChain->SetBranchStatus("p",1);
  if (filterP) fChain->SetBranchStatus("trkP",1);

  if (filterC || subRC) fChain->SetBranchStatus("charge",1); // v4
  
  if (usePFHC) fChain->SetBranchStatus("PFHC_energy",1);
  if (usePFEC) fChain->SetBranchStatus("PFEC_energy",1);

  if (doDepth1) fChain->SetBranchStatus("hcalDepthFractions",1);
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  // Regular L2Relative eta binning
  double veta[] =
    {0, 0.087, 0.174, 0.261, 0.348, 0.435,
     0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
     1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
     2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
     4.363, 4.538, 4.716, 4.889, 5.191};
  const int neta = sizeof(veta) / sizeof(veta[0]) - 1;

  // Inclusive jets pT binning adapted to single particle gun (x2 to x4)
  // Add extra threshold at 200 GeV, move 507 to 500 to match sample edges
  double vpt[] =
    {0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0,
     5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     //97, 133, 174, 220, 272, 330, 395, 507, 592, 686, 790, 905, 1032,
     97, 133, 174, 200., 220, 272, 330, 395, 500., 592, 686, 790, 905, 1032,
     1327, 1684, 2116, 2640, 3273, 4037, 5000};
  double vptf[] =
    {0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
     5, 6, 7, 8, 9, 10., 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 200., 220, 245, 272, 300, 330, 362, 395, 430, 468,
     500., 592, 686, 790, 905, 1032,
     1327, 1684, 2116, 2640, 3273, 4037, 5000.};
  //97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  //507, 592, 686, 790, 905, 1032, 1327, 1684, 2116, 2640, 3273, 4037, 5000};
  //507, 548, 592, 638, 686, 737, 790, 846, 905, 967, //1000,
  // 1032, 1172, 1327, 1497, 1684, 1890, 2116, 2366,
  // 2640, 2941, 3273, 3637, 4037, 4477, 5000};
  //1032, 1101, 1172, 1248,
  //1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
  //2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
  //5000};//4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  double npt = sizeof(vpt) / sizeof(vpt[0]) - 1;
  double nptf = sizeof(vptf) / sizeof(vptf[0]) - 1;

  /*
  const int nfe = 20;
  double vfe[nfe+1];
  for (int i = 0; i != nfe+1; ++i) {
    vfe[i] = 1./nfe*i;
  }
  vfe[nfe] = 1+1e-4; // ensure fE==1 does not go to overflow
  */
  const int nfe = 40;
  double vfe[nfe+1];
  for (int i = 0; i != 10; ++i) {
    vfe[i] = 0.1/10*i;
  }
  for (int i = 10; i != 30; ++i) {
    vfe[i] = 0.1+0.8/20*(i-10);
  }
  for (int i = 30; i != 40+1; ++i) {
    vfe[i] = 0.9+0.1/10*(i-30);
  }
  vfe[nfe] = 1+1e-4; // ensure fE==1 does not go to overflow

  // Response bins
  const int nr = 400;
  double vr[nr+1];
  for (int i = 0; i != nr+1; ++i) vr[i] = i*(4.-0.)/nr;
  
  TDirectory *curdir = gDirectory;
  TFile *fout = new TFile("piongun.root","RECREATE");

  // MIP energy threshold in ECAL for considering hadron an H-hadron
  double e_mip = 1.;

  // General particle distribution for uniformity checks
  TH2D *h2np = new TH2D("h2np",";p_{gen};#eta_{gen};",npt,vpt,10*32,0,3.2);
  TH2D *h2npt = new TH2D("h2npt",";p_{T,gen};#eta_{gen};",npt,vpt,10*32,0,3.2);

  // Efficiencies and H-hadron fractions
  TProfile2D *p2e_trk, *p2e_cal, *p2e_sel, *p2e_had, *p2e_has;
  p2e_trk = new TProfile2D("p2e_trk",";p_{T,gen} (GeV);#eta_{gen};"
			   "Tracking efficiency only",npt,vpt, 10*32, 0,3.2);
  p2e_cal = new TProfile2D("p2e_cal",";p_{T,gen} (GeV);#eta_{gen};"
			   "Calorimeter efficiency only",npt,vpt, 10*32, 0,3.2);
  p2e_sel = new TProfile2D("p2e_sel",";p_{T,gen} (GeV);#eta_{gen};"
			   "Selection efficiency",npt,vpt, 10*32, 0,3.2);
  p2e_had = new TProfile2D("p2e_had",";p_{T,gen} (GeV);#eta_{gen};"
			   "H-hadron fraction",npt,vpt, 10*32, 0,3.2);
  p2e_has = new TProfile2D("p2e_has",";p_{T,gen} (GeV);#eta_{gen};"
			   "H-hadron fraction after efficiency",
			   npt,vpt, 10*32, 0,3.2);

  TProfile2D *p2r_all, *p2r_hhh, *p2r_ehh;
  p2r_all = new TProfile2D("p2r_all",";p_{T,gen} (GeV);#eta_{gen};"
			   "Hadron response",
			   npt,vpt, 10*32, 0,3.2);
  p2r_hhh = new TProfile2D("p2r_hhh",";p_{T,gen} (GeV);#eta_{gen};"
			   "H-hadron response",
			   npt,vpt, 10*32, 0,3.2);
  p2r_ehh = new TProfile2D("p2r_ehh",";p_{T,gen} (GeV);#eta_{gen};"
			   "EH-hadron response",
			   npt,vpt, 10*32, 0,3.2);
  
  TProfile2D *p2e, *p2h, *p2r_a, *p2r_h, *p2r_e;
  p2e = new TProfile2D("p2e",";p_{T,gen} (GeV);#eta_{gen};Efficiency",
		       npt,vpt, neta,veta);
  p2h = new TProfile2D("p2h",";p_{T,gen} (GeV);#eta_{gen};H fraction",
		       npt,vpt, neta,veta);
  p2r_a = new TProfile2D("p2r_a",";p_{T,gen} (GeV);|#eta_{gen}|;Response (A)",
  		       npt,vpt, neta,veta);
  p2r_h = new TProfile2D("p2r_h",";p_{T,gen} (GeV);|#eta_{gen}|;Response (H)",
			 npt,vpt, neta,veta);
  p2r_e = new TProfile2D("p2r_e",";p_{T,gen} (GeV);|#eta_{gen}|;Response (EH)",
			 npt,vpt, neta,veta);

  // For closure tests
  TProfile2D *p2c_a, *p2c_h, *p2c_e;
  p2c_a = new TProfile2D("p2c_a",";p_{T,gen} (GeV);|#eta_{gen}|;"
			 "Corrrection^{-1} (A)",npt,vpt, neta,veta);
  p2c_h = new TProfile2D("p2c_h",";p_{T,gen} (GeV);|#eta_{gen}|;"
			 "Correction^{-1} (H)",npt,vpt, neta,veta);
  p2c_e = new TProfile2D("p2c_e",";p_{T,gen} (GeV);|#eta_{gen}|;"
			 "Correction^{-1} (EH)",npt,vpt, neta,veta);

  TProfile *pe_bb0, *pe_bb, *pe_ec1, *pe_ec2;
  pe_bb0 = new TProfile("pe_bb0",";p_{T,gen} (GeV);Efficiency (BB0)",npt,vpt);
  pe_bb = new TProfile("pe_bb",";p_{T,gen} (GeV);Efficiency (BB)",npt,vpt);
  pe_ec1 = new TProfile("pe_ec1",";p_{T,gen} (GeV);Efficiency (EC1)",npt,vpt);
  pe_ec2 = new TProfile("pe_ec2",";p_{T,gen} (GeV);Efficiency (EC2)",npt,vpt);
  TProfile2D *p2rf_bb0, *p2rf_bb, *p2rf_ec1, *p2rf_ec2;
  //TProfile2D *p2rfmip;
  p2rf_bb0 = new TProfile2D("p2rf_bb0",";p_{T,gen} (GeV);f_{ECAL,raw};"
			    "Response (BB0)",npt,vpt, nfe,vfe);
  p2rf_bb = new TProfile2D("p2rf_bb",";p_{T,gen} (GeV);f_{ECAL,raw};"
			   "Response (BB)",npt,vpt, nfe,vfe);
  p2rf_ec1 = new TProfile2D("p2rf_ec1",";p_{T,gen} (GeV);f_{ECAL,raw};"
			   "Response (EC1)",npt,vpt, nfe,vfe);
  p2rf_ec2 = new TProfile2D("p2rf_ec2",";p_{T,gen} (GeV);f_{ECAL,raw};"
			   "Response (EC2)",npt,vpt, nfe,vfe);
  //p2rfmip = new TProfile2D("p2rfmip",";p_{T,gen} (GeV);f_{ECAL,raw};Response",
  //			   npt,vpt, nfe,vfe);
  TH2D *h2rf_bb0, *h2rf_bb, *h2rf_ec1, *h2rf_ec2;
  h2rf_bb0 = new TH2D("h2rf_bb0",";p_{T,gen} (GeV);f_{ECAL,raw};N_{ptcl} (BB0)",
		      npt,vpt, nfe,vfe);
  h2rf_bb = new TH2D("h2rf_bb",";p_{T,gen} (GeV);f_{ECAL,raw};N_{ptcl} (BB)",
		     npt,vpt, nfe,vfe);
  h2rf_ec1 = new TH2D("h2rf_ec1",";p_{T,gen} (GeV);f_{ECAL,raw};N_{ptcl} (EC1)",
		      npt,vpt, nfe,vfe);
  h2rf_ec2 = new TH2D("h2rf_ec2",";p_{T,gen} (GeV);f_{ECAL,raw};N_{ptcl} (EC2)",
		      npt,vpt, nfe,vfe);

  TProfile3D *p3rf;
  p3rf = new TProfile3D("p3rf",";|#eta_{gen}|;p_{T,gen} (GeV);f_{ECAL,raw};"
			"Response",neta,veta, npt,vpt, nfe,vfe);

  // For closure tests
  TProfile3D *p3cf;
  p3cf = new TProfile3D("p3cf",";|#eta_{gen}|;p_{T,gen} (GeV);f_{ECAL,raw}"
			"Correction^{-1}",neta,veta, npt,vpt, nfe,vfe);

  // For resolution studies
  TH3D *h3r, *h3c;
  h3r = new TH3D("h3r",";|#eta_{gen}|;p_{T,gen} (GeV);Response;",
		 neta,veta, npt,vpt, nr,vr);
  h3c = new TH3D("h3c",";|#eta_{gen}|;p_{T,gen} (GeV);Correction^{-1};",
		 neta,veta, npt,vpt, nr,vr);


  // For depth1 studies
  const int np1 = 8;
  TProfile2D *p2r(0), *p2re(0), *p2rh(0);
  TProfile *vpf1[np1];
  TProfile *vph1[np1];
  TProfile *vpe1[np1];
  TProfile *vpr1[np1];
  if (doDepth1) {
    fout->mkdir("depth1");
    fout->cd("depth1");

    p2r = new TProfile2D("p2r",";p_{T,gen};f_{depth1};Response",nptf,vptf, nfe,vfe);
    p2re = new TProfile2D("p2re",";p_{T,gen};f_{depth1};Response (e>0)",nptf,vptf, nfe,vfe);
    p2rh = new TProfile2D("p2rh",";p_{T,gen};f_{depth1};Response (e=0)",nptf,vptf, nfe,vfe);
      
    for (int i = 0; i != np1; ++i) {
      vpf1[i] = new TProfile(Form("pf1_%d",i),";p_{T,gen};Depth1 category fraction;",nptf,vptf);
    }
    for (int i = 0; i != np1; ++i) {
      vph1[i] = new TProfile(Form("ph1_%d",i),";p_{T,gen};Depth1 HCAL fraction;",nptf,vptf);
    }
    for (int i = 0; i != np1; ++i) {
      vpe1[i] = new TProfile(Form("pe1_%d",i),";p_{T,gen};Depth1 energy fraction;",nptf,vptf);
    }
    for (int i = 0; i != np1; ++i) {
      vpr1[i] = new TProfile(Form("pr1_%d",i),";p_{T,gen};Depth1 category response;",nptf,vptf);
    } // for i
    fout->cd();
  } // doDepth1
    
  // Store new PFEC as a friend treee
  TFile *outFriend(0);
  TTree *friendTree(0);
  Float_t PFECC_energy(0.0);
  Float_t PFECN_energy(0.0);
  if (applyPFEC_Charged || applyPFEC_Neutral) {
    outFriend = new TFile("piongunFriend.root","RECREATE");
    friendTree = new TTree("piongunFriend","piongunFriend");
  }
  if (applyPFEC_Charged) {
    friendTree->Branch("PFECC_energy",&PFECC_energy,"PFECC_energy/F");
  }
  if (applyPFEC_Neutral) {
    friendTree->Branch("PFECN_energy",&PFECN_energy,"PFECN_energy/F");
  }

  curdir->cd();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    rcEcal = rcHcal = 0; // safety for noPU without these branches
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if (jentry%100000==0) cout << "." << flush;

    double abseta = fabs(genEta);
    double genPt = genP/cosh(genEta);
    //double fe = ((rawHcal+rawEcal)>0 ? rawEcal / (rawEcal+rawHcal) :
    //		 (jentry%3==0 ? 0 : 1)); // v1
    double fe = ((rawHcal+rawEcal)>0 ? rawEcal / (rawEcal+rawHcal) : 1); // v2
    double resp = (genP>0 ? (rawEcal+rawHcal) / genP : 0);
    if (charge==0) // 2025v2
      resp = (genP>0 ? ((rawEcal-rcEcal)+(rawHcal-rcHcal)) / genP : 0);
    //double eff = ((rawHcal+rawEcal)>0 ? 1 : 0);
    // Patch Conrado Munoz Diaz's tuples for p==0 in tracker coverage
    // Patch V2 tup[les with !filterP for missing p altogether
    double eff = ((rawHcal+rawEcal)>0 && genP>0 &&
		  (!usePFHC || (pfhcE>0 && pfhcE<13000.)) &&
		  (!filterMaxP || (genP<maxP)) &&
		  (!filterP || (p>0 || abseta>2.6)) &&
		  //(!filterP || (p>0 || abseta>2.5)) &&
		  //&& genP>0 ? 1 : 0);
		  //&& genP>0 && charge!=0 ? 1 : 0);
		  //&& genP>0 && (charge!=0 || abseta>2.65) ? 1 : 0);
		  //&& genP>0 && (charge!=0 || abseta>2.5) ? 1 : 0);
		  //(!filterC || (charge!=0 || abseta>2.322)) ? 1 : 0); // W25
    		  //(!filterC  || ((charge!=0 && abseta<2.6) ||
		  //		 (charge==0 && abseta>=2.6)))  // 2025v2
		  (!filterC  || ((charge!=0 && abseta<2.6) ||
				 (charge==0 && abseta>=2.6)))  // 2026v2
		  ? 1 : 0);
    //(!filterC || (charge!=0 || abseta>2.65)) ? 1 : 0); // W24
    bool ish = (fe<0.01);// as in drawPionGun.C
    bool ise = (fe>0.2 && fe<0.8);// as in drawPionGun.C
    assert(!(ish && ise));

    // Apply PFEC
    double corr(1.);
    if (usePFHC) {
      if (eff>0 && charge!=0)
	corr = pfhcE / (rawEcal+rawHcal);
      if (eff>0 && ((rawEcal-rcEcal)+(rawHcal-rcHcal)>0) && charge==0)
	corr = pfhcE / (rawEcal+rawHcal-rcEcal-rcHcal);
    }
    if (usePFEC) {
      if (eff>0 && charge!=0)
	corr = pfecE / (rawEcal+rawHcal);
      if (eff>0 && ((rawEcal-rcEcal)+(rawHcal-rcHcal)>0) && charge==0)
	corr = pfecE / (rawEcal+rawHcal-rcEcal-rcHcal);
    }
    // Predict calorimeter response based on track/true E (genP here)
    if (applyPFEC_Charged) {
      double corrEcal(charge!=0 ? rawEcal : rawEcal-rcEcal);
      double corrHcal(charge!=0 ? rawHcal : rawHcal-rcHcal);
      if (eff>0) pfec->energyEmHad(genP, corrEcal, corrHcal, genEta, 0.);
      if (eff>0 && charge!= 0)
	corr = (corrEcal+corrHcal)/(rawEcal+rawHcal);
      if (eff>0 && ((rawEcal-rcEcal)+(rawHcal-rcHcal)>0) && charge==0)
	corr = (corrEcal+corrHcal)/(rawEcal+rawHcal-rcEcal-rcHcal);
      PFECC_energy = (corrEcal+corrHcal);
    }
    // Predict calo response based on measured E (-1 => use corrEcal+corrHcal)
    if (applyPFEC_Neutral) {
      double corrEcal(charge!=0 ? rawEcal : rawEcal-rcEcal);
      double corrHcal(charge!=0 ? rawHcal : rawHcal-rcHcal);
      if (eff>0) pfec->energyEmHad(-1, corrEcal, corrHcal, genEta, 0.);
      if (eff>0 && charge!=0)
	corr = (corrEcal+corrHcal)/(rawEcal+rawHcal);
      if (eff>0 && ((rawEcal-rcEcal)+(rawHcal-rcHcal)>0) && charge==0)
	corr = (corrEcal+corrHcal)/(rawEcal+rawHcal-rcEcal-rcHcal);
      PFECN_energy = (corrEcal+corrHcal);
    }
    if (applyPFEC_Charged || applyPFEC_Neutral) friendTree->Fill();
    resp *= corr;

    h2np->Fill(genP, abseta);
    h2npt->Fill(genPt, abseta);

    p2e_trk->Fill(genPt, abseta, charge!=0);
    p2e_cal->Fill(genPt, abseta, (rawEcal+rawHcal)>0);
    p2e_sel->Fill(genPt, abseta, eff);
    p2e_had->Fill(genPt, abseta, ish ? 1 : 0);
    if (eff>0) p2e_has->Fill(genPt, abseta, ish ? 1 : 0);
    if (eff>0) p2r_all->Fill(genPt, abseta, resp);
    if (eff>0 && ish) p2r_hhh->Fill(genPt, abseta, resp);
    if (eff>0 && ise) p2r_ehh->Fill(genPt, abseta, resp);
    
    p2e->Fill(genPt, abseta, eff);
    if (abseta<0.522) pe_bb0->Fill(genPt, eff);
    if (abseta<1.479) pe_bb->Fill(genPt, eff);
    if (abseta>1.479 && abseta<2.5) pe_ec1->Fill(genPt, eff);
    if (abseta>2.5)   pe_ec2->Fill(genPt, eff);
    
    if (eff>0) {
      p2h->Fill(genPt, abseta, ish ? 1 : 0);
      p2r_a->Fill(genPt, abseta, resp);
      if (ish) p2r_h->Fill(genPt, abseta, resp);
      if (ise) p2r_e->Fill(genPt, abseta, resp);

      p2c_a->Fill(genPt, abseta, 1./corr);
      if (ish) p2c_h->Fill(genPt, abseta, 1./corr);
      if (ise) p2c_e->Fill(genPt, abseta, 1./corr);
      
      if (abseta<0.522) {
	p2rf_bb0->Fill(genPt, fe, resp);
	h2rf_bb0->Fill(genPt, fe);
      }
      //if (abseta<1.131) {
      if (abseta<1.479) {
	p2rf_bb->Fill(genPt, fe, resp);
	h2rf_bb->Fill(genPt, fe);
      }
      //if (abseta>1.566 && abseta<2.5) {
      if (abseta>1.479 && abseta<2.5) {
	p2rf_ec1->Fill(genPt, fe, resp);
	h2rf_ec1->Fill(genPt, fe);
      }
      if (abseta>2.5) {
	p2rf_ec2->Fill(genPt, fe, resp);
	h2rf_ec2->Fill(genPt, fe);
      }

      //if (resp>0 && resp<3.)
      p3rf->Fill(abseta, genPt, fe, resp);
      p3cf->Fill(abseta, genPt, fe, 1./corr);

      h3r->Fill(abseta, genPt, resp);
      h3c->Fill(abseta, genPt, 1./corr);
    } // eff>0


    // Depth1 studies
    if (doDepth1 && charge!=0 && genP>0 && fabs(genEta)<1.305) {

      double f1 = hcalDepthFractions[0];
      double e1 = hcalDepthFractions[0] * rawHcal / genP;
      double rp = (rawEcal+rawHcal) / genP;
      int cat(-1);
      if (rawEcal==0  && rawHcal>0  && f1==1)        cat = 0;
      if (rawEcal==0  && rawHcal>0  && f1>0 && f1<1) cat = 1;
      if (rawEcal==0  && rawHcal>0  && f1==0)        cat = 2;
      if (rawEcal>0   && rawHcal>0  && f1==1)        cat = 3;
      if (rawEcal>0   && rawHcal>0  && f1>0 && f1<1) cat = 4;
      if (rawEcal>0   && rawHcal>0  && f1==0)        cat = 5;
      if (rawEcal>0   && rawHcal==0)                 cat = 6;
      if (rawEcal==0  && rawHcal==0)                 cat = 7;
      assert(cat>=0);

      double f = (rawEcal+rawHcal>0 ? f1 * rawHcal / (rawHcal+rawEcal) : 0.);
      p2r->Fill(genP, f, rp);
      if (rawEcal>0)   p2re->Fill(genP, f, rp);
      if (rawEcal==0)  p2rh->Fill(genP, f, rp);
      
      for (int i = 0; i != np1; ++i) {
	vpf1[i]->Fill(genP, i==cat ? 1 : 0);
	if (i==cat) {
	  vph1[i]->Fill(genP, f1);
	  vpe1[i]->Fill(genP, e1);
	  vpr1[i]->Fill(genP, rp);
	}
      } // for i
    } // doDepth1
  } // for jentry

  fout->Write();
  fout->Close();

  if (applyPFEC_Charged || applyPFEC_Neutral) {
    friendTree->Write();
    outFriend->Close();
  }
}
