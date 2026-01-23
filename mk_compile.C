// Run this script to compile CondFormats libraries. After this can easily run 
// root -l -b -q mk_GamHistosFill.C
// using R__LOAD_LIBRARY to load *.so
{

  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // For Gamjet code (v6.30/04)
  //gROOT->ProcessLine(".L pfhadronsRun3.C+g");
  //gROOT->ProcessLine(".L PFEnergyCalibrationFromMikko.cc+g");
  gROOT->ProcessLine(".L piongun.C++g"); // ++g to ensure updating #include .txt

  gROOT->ProcessLine(".L drawPiongun.C+g");
  gROOT->ProcessLine(".L drawPiongunResolution.C+g");
  gROOT->ProcessLine(".L drawPiongunEfficiency.C+g");

  gROOT->ProcessLine(".L PFAna.C+g");

}
