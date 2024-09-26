#include "piongun.h"

//R__LOAD_LIBRARY(piongun.C+g);
R__LOAD_LIBRARY(piongun_C.so);

void mk_piongun() {

  TChain *c = new TChain("s");
  //c->AddFile("../data/piongun/SinglePion_0p2to500_Run3Winter24_fromConrado.root"); // original derivation
  c->AddFile("../data/piongun/NTuple_0_500_PFHC_PFEC_24MC24Corrections_Closures.root"); // closure tests
  
  piongun pg(c);
  pg.Loop();
  
} // mk_piongun
