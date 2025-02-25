#include "PFAna.h"

//#include "drawPFAna.C"

//R__LOAD_LIBRARY(PFAna.C+g);
R__LOAD_LIBRARY(PFAna_C.so);

void mk_PFAna() {

  TChain *c = new TChain("pfAnaNtuplizer/PFAnaTree");
  //c->AddFile("../data/piongun/MergedNtuple_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root");
  c->AddFile("../data/piongun/MergedNtuple_Premix_SinglePionMinus_E-50_Eta-2p5to3p0_Run1to100.root");
  
  PFAna pf(c);
  pf.Loop();

  //gROOT->ProcessLine(".L drawPFAna.C+g");
  //drawPFAna();
  
} // mk_PFAna
