#include "piongun.h"

#include "drawPiongun.C"

//R__LOAD_LIBRARY(piongun.C+g);
R__LOAD_LIBRARY(piongun_C.so);

void mk_piongun() {

  TChain *c = new TChain("s");
  //c->AddFile("../data/piongun/SinglePion_0p2to500_Run3Winter24_fromConrado.root"); // original derivation
  //c->AddFile("../data/piongun/NTuple_0_500_PFHC_PFEC_24MC24Corrections_Closures.root"); // closure tests

  //c->AddFile("/Users/manvouti/Downloads/2024_smallSample/2024_smallSample_merged.root");
  //c->AddFile("/Users/voutila/Downloads/2024_smallSampleII/2024_smallSampleII_merged.root");

  //c->AddFile("../data/piongun/Conrado-20241031_smallSample.root");
  //c->AddFile("../data/piongun/2025_E0p2to200GeV.root");
  //c->AddFile("../data/piongun/2025_E200to500GeV.root");

  // Winter25v2 (beamspot+HE scale)
  //c->AddFile("2025_0p2to5000GeV_withCorrections_and_PowerLaw_NoPU.root");
  //c->AddFile("2025_0p2to5000GeV_withCorrections_and_PowerLaw_withPU.root");
  
  // Winter25v1
  c->AddFile("../data/piongun/2025_E0p2to200GeV_v2.root");
  c->AddFile("../data/piongun/2025_E200to500GeV_v2.root");
  // Winter25 PFHC closure
  //c->AddFile("../data/piongun/2025_0p2to200GeV_withCorrections.root");
  //c->AddFile("../data/piongun/2025_200to500GeV_withCorrections.root");

  // Winter24 refence for resolution
  //c->AddFile("../data/piongun/2024_E0p2to200GeV_v10.root");
  //c->AddFile("../data/piongun/2024_E200to500GeV_v9.root");
  
  piongun pg(c);
  pg.Loop();

  gROOT->ProcessLine(".L drawPiongun.C+g");
  drawPiongun();
  
} // mk_piongun
