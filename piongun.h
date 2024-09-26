//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 18 11:15:03 2024 by ROOT version 6.30/04
// from TTree s/ PFCalibration
// found on file: SinglePion_0p2to500_Run3Winter24_fromCondado.root
//////////////////////////////////////////////////////////

#ifndef piongun_h
#define piongun_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class piongun {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
  /*
   Float_t         genP; //true;  // Q: what is true? genE, genP?
   Float_t         p;  // Q: is p trkP?
   Float_t         rawEcal; //ecal;
   Float_t         rawHcal; //hcal;
  //Float_t         rawHcal; //hcal+ho;
   Float_t         ho; // new 
   Float_t         genEta; //eta; // Q: is eta genEta?
   Float_t         phi;
  */
   Int_t           charge;
   vector<float>   *dr;
   vector<float>   *Eecal;
   vector<float>   *Ehcal;
   vector<float>   *pfcID;
   vector<int>     *pfcs;
   vector<float>   *correcal;
   vector<float>   *corrhcal;
   Float_t         Ccorrecal;
   Float_t         Ccorrhcal;
   ULong64_t       run;
   ULong64_t       evt;
   ULong64_t       lumiBlock;
   ULong64_t       time;

   // New tuples from Conrado change variable types
   Float_t         p; // not there, just placeholder
   Double_t        genP;//true;
   Double_t        rawEcal;//ecal;
   Double_t        rawHcal;//hcal;
   Double_t        genEta;//eta;
   Double_t        phi;
  //Double_t        PFHC_energy;
  //Double_t        PFEC_energy;
  //Double_t        PFHC_closure;
  //Double_t        PFEC_closure;
   Double_t          pfecE;
   Double_t          pfhcE;

   // List of branches
   TBranch        *b_true;   //!
   TBranch        *b_p;   //!
   TBranch        *b_ecal;   //!
   TBranch        *b_hcal;   //!
   TBranch        *b_ho;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_dr;   //!
   TBranch        *b_Eecal;   //!
   TBranch        *b_Ehcal;   //!
   TBranch        *b_pfcID;   //!
   TBranch        *b_pfcs;   //!
   TBranch        *b_correcal;   //!
   TBranch        *b_corrhcal;   //!
   TBranch        *b_Ccorrecal;   //!
   TBranch        *b_Ccorrhcal;   //!
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_time;   //!

   TBranch        *b_pfecE;
   TBranch        *b_pfhcE;

   piongun(TTree *tree=0);
   virtual ~piongun();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef piongun_cxx
piongun::piongun(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SinglePion_0p2to500_Run3Winter24_fromCondado.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SinglePion_0p2to500_Run3Winter24_fromCondado.root");
      }
      f->GetObject("s",tree);

   }
   Init(tree);
}

piongun::~piongun()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t piongun::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t piongun::LoadTree(Long64_t entry)
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

void piongun::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dr = 0;
   Eecal = 0;
   Ehcal = 0;
   pfcID = 0;
   pfcs = 0;
   correcal = 0;
   corrhcal = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // //fChain->SetBranchAddress("true", &true, &b_true);
   fChain->SetBranchAddress("true", &genP, &b_true);
   //fChain->SetBranchAddress("p", &p, &b_p);
   //fChain->SetBranchAddress("ecal", &ecal, &b_ecal);
   fChain->SetBranchAddress("ecal", &rawEcal, &b_ecal);
   //fChain->SetBranchAddress("hcal", &hcal, &b_hcal);
   fChain->SetBranchAddress("hcal", &rawHcal, &b_hcal);
   //fChain->SetBranchAddress("ho", &ho, &b_ho);
   // //fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("eta", &genEta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   /*
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("dr", &dr, &b_dr);
   fChain->SetBranchAddress("Eecal", &Eecal, &b_Eecal);
   fChain->SetBranchAddress("Ehcal", &Ehcal, &b_Ehcal);
   fChain->SetBranchAddress("pfcID", &pfcID, &b_pfcID);
   fChain->SetBranchAddress("pfcs", &pfcs, &b_pfcs);
   fChain->SetBranchAddress("correcal", &correcal, &b_correcal);
   fChain->SetBranchAddress("corrhcal", &corrhcal, &b_corrhcal);
   fChain->SetBranchAddress("Ccorrecal", &Ccorrecal, &b_Ccorrecal);
   fChain->SetBranchAddress("Ccorrhcal", &Ccorrhcal, &b_Ccorrhcal);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("time", &time, &b_time);
   */

   fChain->SetBranchAddress("PFEC_energy", &pfecE, &b_pfecE);
   fChain->SetBranchAddress("PFHC_energy", &pfhcE, &b_pfhcE);
   
   Notify();
}

Bool_t piongun::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void piongun::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t piongun::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef piongun_cxx
