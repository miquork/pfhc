#include "piongun.h"

#include "drawPiongun.C"

//R__LOAD_LIBRARY(piongun.C+g);
R__LOAD_LIBRARY(piongun_C.so);
R__LOAD_LIBRARY(drawPiongun_C.so);
R__LOAD_LIBRARY(drawPiongunResolution_C.so);
R__LOAD_LIBRARY(drawPiongunEfficiency_C.so);

void mk_piongun() {

  TChain *c = new TChain("s");
  //c->AddFile("../data/piongun/SinglePion_0p2to500_Run3Winter24_fromConrado.root"); // original derivation
  //c->AddFile("../data/piongun/NTuple_0_500_PFHC_PFEC_24MC24Corrections_Closures.root"); // closure tests

  //c->AddFile("/Users/manvouti/Downloads/2024_smallSample/2024_smallSample_merged.root");
  //c->AddFile("/Users/voutila/Downloads/2024_smallSampleII/2024_smallSampleII_merged.root");

  //c->AddFile("../data/piongun/Conrado-20241031_smallSample.root");
  //c->AddFile("../data/piongun/2025_E0p2to200GeV.root");
  //c->AddFile("../data/piongun/2025_E200to500GeV.root");

  // New PFHC derivation for 2026, first GTv9
  /*
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to200_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E200to500_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E500to5000_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  */
  // Original GTv7 reference for PFHC25
  /*
  string tag = "PFHC25_GTv7";
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to10_Winter25_NoPU_GTv7.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to200_Winter25_NoPU_GTv7.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E200to500_Winter25_NoPU_GTv7.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E500to5000_Winter25_NoPU_GTv7.root");
  */
  /*
  string tag = "PFHC26_GTv7_2025MOYv4";
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to10_Winter25_NoPU_GTv7_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to200_Winter25_NoPU_GTv7_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E200to500_Winter25_NoPU_GTv7_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E500to5000_Winter25_NoPU_GTv7_HcalRespCorrs2025MOYv4.root");
  */

  // PFHC26 version v1
  /*
  //string tag = "PFHC26_GTv9_2025MOYv4";
  string tag = "PFHC26v3";
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to10_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to200_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E200to500_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E500to5000_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4.root");
  */

  // With power law corrections
  string tag = "PFHC26v3";
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to10_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4_PFHCPowerLaw.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E0p2to200_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4_PFHCPowerLaw.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E200to500_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4_PFHCPowerLaw.root");
  c->AddFile("../data/piongun/FikriPFHC26/MERGED/Ntuples_Pi_E500to5000_Winter25_NoPU_GTv9_HcalRespCorrs2025MOYv4_PFHCPowerLaw.root");

  // New withRC Winter25v2 (beamspot+HE scale)
  //c->AddFile("../data/piongun/2025_0p2to5000GeV_withCorrections_and_PowerLaw_NoPU_bfix_v2.root"); // baseline v2, "HB1X", has hcal>0 cut
  //c->AddFile("../data/piongun/2025_0p2to5000GeV_withCorrections_and_PowerLaw_NoPU_bfix.root"); // baseline v1, no hcal>0 cut? // Last baseline from 2025
  //c->AddFile("../data/piongun/2025_0p2to500GeV_withCorrections_randomCone.root"); // withPU?
  //
  // Winter25v2 (beamspot+HE scale)
  //c->AddFile("../data/piongun/2025_0p2to5000GeV_withCorrections_and_PowerLaw_NoPU.root"); // old
  //c->AddFile("../data/piongun/2025_0p2to5000GeV_withCorrections_and_PowerLaw_withPU.root"); // old

  // HB cuts test
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_5xPFThres_HcalCuts_compressed.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_4xPFThres_HcalCuts_compressed.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_3xPFThres_HcalCuts_compressed.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_2xPFThres_HcalCuts_compressed.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_1p5xPFThres_HcalCuts_compressed.root");
  
  // EE cuts test
  //c->AddFile("../data/piongun/PFHC_0to200_EEcuts_4sigma.root"); // old
  //c->AddFile("../data/piongun/PFHC_0to200_EEcuts_0sigma.root"); // old
  // Proper updated files:
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_zeroEcalThres_reduced.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_1sigmaEcalThres_reduced.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_2sigmaEcalThres_reduced.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_3sigmaEcalThres_reduced.root");
  //c->AddFile("../data/piongun/2025_Merged_NoPU_pionGun_4sigmaEcalThres_reduced.root");
  
  // Winter25v1
  //c->AddFile("../data/piongun/2025_E0p2to200GeV_v2.root");
  //c->AddFile("../data/piongun/2025_E200to500GeV_v2.root");
  // Winter25 PFHC closure
  //c->AddFile("../data/piongun/2025_0p2to200GeV_withCorrections.root");
  //c->AddFile("../data/piongun/2025_200to500GeV_withCorrections.root");

  // Winter24 refence for resolution
  //c->AddFile("../data/piongun/2024_E0p2to200GeV_v10.root");
  //c->AddFile("../data/piongun/2024_E200to500GeV_v9.root");
  
  piongun pg(c);
  pg.Loop();

  const char *ct = tag.c_str();
  gROOT->ProcessLine(Form(".! cp -pi piongun.root piongun_%s.root",ct));

  gROOT->ProcessLine(".L drawPiongun.C+g");
  drawPiongun("piongun.root",tag);

  gROOT->ProcessLine(Form(".! cp -pi drawPiongun.root drawPiongun_%s.root",ct));

  gROOT->ProcessLine(".L drawPiongunResolution.C+g");
  drawPiongunResolution();//"piongun.root","","","",tag);

  //gROOT->ProcessLine(Form(".! cp -pi drawPiongunResolution.root drawPiongunResolution_%s.root",ct));

  gROOT->ProcessLine(".L drawPiongunEfficiency.C+g");
  drawPiongunEfficiency();

  
} // mk_piongun
