// Purpose: Exploratory study of HCAL depth fractions to compare with IsoTrack
//          1) select E=40-60 GeV range (wrt E=50 GeV for IsoTrack)
//          2) compare EH vs H hadrons (vs H for IsoTrack)
//          3) compare tracker vs gen truth
//          4) directly compare to HCAL on the same footing
#include "TChain.h"
#include "TH2D.h"

#include <iostream>
#include <vector>

#include "tdrstyle_mod22.C"

void fillHistos();

void hcalDepthFractions() {
  
  fillHistos();
}


void fillHistos() {
  
  cout << "Calling fillHistos" << endl << flush;
  
  TChain *c = new TChain("s");
  //c->AddFile("../data/piongun/2025_E0p2to200GeV.root");
  //c->AddFile("../data/piongun/2025_E200to500GeV.root");
  c->AddFile("../data/piongun/2025_E0p2to200GeV_v2.root");
  //c->AddFile("../data/piongun/2025_E200to500GeV_v2.root");
  
  Float_t        genP;     // true;
  Float_t        trkP;     // p;
  Float_t        rawEcal;  // ecal;
  Float_t        rawHcal;  // hcal;
  Float_t        rawHFEm;  // hfem;
  Float_t        rawHFHad; // hfhaf;
  Float_t        genEta;   // eta;
  Float_t        genPhi;   // phi
  Int_t          charge;
  const int ndepth = 7;
  Float_t        hcalDepthFractions[ndepth];

  TBranch        *b_genP;
  cout << "Set branch addresses." << endl << flush;
  //c->SetBranchAddress("true", &genP, &b_genP);
  c->SetBranchAddress("genP", &genP, &b_genP);
  //c->SetBranchAddress("p", &trkP);
  c->SetBranchAddress("trkP", &trkP);
  c->SetBranchAddress("ecal", &rawEcal);
  c->SetBranchAddress("hcal", &rawHcal);
  c->SetBranchAddress("hfem", &rawHFEm);
  c->SetBranchAddress("hfhad", &rawHFHad);
  c->SetBranchAddress("eta", &genEta);
  c->SetBranchAddress("phi", &genPhi);
  c->SetBranchAddress("charge", &charge);
  c->SetBranchAddress("hcalDepthFractions", hcalDepthFractions);

  cout << "Setup histograms." << endl << flush;
  const int nb = 120;//240;//60;
  TH2D *h2 = new TH2D("h2",";#eta;p (GeV)",nb,-3,3,100,0,500);
  TProfile* vd[ndepth];
  TProfile* vd2[ndepth];
  TProfile *pcal, *pecal, *phcal, *phcal2, *psum, *psum2;
  pcal = new TProfile("pcal",";#eta;(ECAL+HCAL)/E_{true}", nb,-3,3);
  pecal = new TProfile("pecal",";#eta;ECAL/E_{true}", nb,-3,3);
  phcal = new TProfile("phcal",";#eta;HCAL/E_{true}", nb,-3,3);
  phcal2 = new TProfile("phcal2",";#eta;HCAL/E_{true}", nb,-3,3);
  psum = new TProfile("psum",";#eta;#sum{HCAL}/E_{true}", nb,-3,3);
  psum2 = new TProfile("psum2",";#eta;#sum{HCAL}/E_{true}", nb,-3,3);

  TH2D *h2sum = new TH2D("h2sum",";#eta;#sum{HCAL}/E_{true}", nb,-3,3, 100,0,1);
  
  for (int depth = 0; depth != ndepth; ++depth) {
    vd[depth] = new TProfile(Form("pdepth%d",depth+1),
			     ";#eta;E_{depth}/E_{true}", nb,-3,3);
    vd2[depth] = new TProfile(Form("pdepth%d_h",depth+1),
			      ";#eta;E_{depth}/E_{true}", nb,-3,3);
  }

  cout << "Count entries." << endl << flush;
  int nentries = c->GetEntries();

  cout << "Loop over " << nentries << " entries." << endl << flush;
  for (Long64_t jentry=0; jentry<nentries; ++jentry) {

    if (jentry%100000==0) cout << "." << flush;

    // Speedup for energy range selection
    Long64_t ientry = c->LoadTree(jentry);
    b_genP->GetEntry(ientry);
    if (genP<=40 || genP>=60) continue;

    // Load rest of the tree
    c->GetEntry(jentry);

    bool pass = (charge!=0 || fabs(genEta)>2.322);
    if (!pass) continue;

    h2->Fill(genEta, genP);
    //pcal->Fill(genEta, (rawEcal + rawHcal) / genP);
    //pecal->Fill(genEta, rawHcal / genP);
    //phcal->Fill(genEta, rawHcal / genP);
    // Check that fractions properly filled
    double sum(0);
    for (int depth = 0; depth != ndepth; ++depth) {
      sum += hcalDepthFractions[depth];
    }
    //if (sum==0) continue;

    // First enegy sums
    pcal->Fill(genEta, (rawEcal + rawHcal + rawHFEm + rawHFHad) / genP);
    pecal->Fill(genEta, rawEcal / genP);
    phcal->Fill(genEta, (rawHcal + rawHFEm + rawHFHad) / genP);
    if (rawEcal<2.) {
      phcal2->Fill(genEta, (rawHcal + rawHFEm + rawHFHad) / genP);
    }

    // Don't do fractions for hcal=0 cases
    if (sum==0) continue;

    // Then fractions
    for (int depth = 0; depth != ndepth; ++depth) {
      //vd[depth]->Fill(genEta, hcalDepthFractions[depth] * rawHcal / genP);
      vd[depth]->Fill(genEta, hcalDepthFractions[depth]);
      //sum += hcalDepthFractions[depth];
      if (rawEcal<2.) {
	//vd2[depth]->Fill(genEta, hcalDepthFractions[depth] * rawHcal / genP);
	vd2[depth]->Fill(genEta, hcalDepthFractions[depth]);
      }
    }
    psum->Fill(genEta, sum);
    h2sum->Fill(genEta, sum);
    if (rawEcal<2.) psum2->Fill(genEta, sum);
  } // for jentry

  //h2->Draw("COLZ");
  setTDRStyle();
  TH1D *h = tdrHist("h","E_{meas}/E_{true}  or  E_{depth}/E_{meas}",0,1.45,
		    "#eta",-3,3);
  extraText = "Private";
  lumi_136TeV = "NoPU PionGun for PFHC, Winter25";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  h->GetYaxis()->SetTitleOffset(1.20);
  
  TLegend *leg1 = tdrLeg(0.75,0.90-0.03*7,1.00,0.90);
  leg1->SetTextSize(0.035);
  TLegend *leg2 = tdrLeg(0.42,0.90-0.035*6,0.67,0.90);
  leg2->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  //tex->DrawLatex(0.18,0.74,"40#leqE_{true}#leq60 GeV");
  tex->DrawLatex(0.18,0.74,"40<E_{true}<60 GeV");
  
  int marker[ndepth] =
    {kOpenTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp,
     kOpenStar, kOpenDiamond, kOpenCross};
  int color[ndepth] =
    {kBlue, kOrange+1, kGreen+1, kRed,
     kYellow+1, kOrange+3, kGray+2};
  double size[ndepth] =
    {0.7, 0.5, 0.5, 0.7,
     1.0, 1.0, 1.0};

  tdrDraw(phcal2,"HIST][",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(pcal,"HIST][",kNone,kBlack,kDashed,-1,kNone,0);
  tdrDraw(phcal,"HP][",kFullCircle,kBlack,kSolid,-1,kNone,0,0.5);
  tdrDraw(pecal,"HP][",kOpenCircle,kBlack,kDotted,-1,kNone,0,0.5);

  tdrDraw(psum,"HP][",kFullStar,kBlue-9,kDotted,-1,kNone,0,0.5);
  tdrDraw(psum2,"HP][",kFullStar,kRed-9,kDotted,-1,kNone,0,0.5);

  for (int i = 0; i != ndepth; ++i) {
    tdrDraw(vd2[i],"HIST][",kNone,color[i],kSolid,-1,kNone,0);
    tdrDraw(vd[i],"HP][",marker[i],color[i],kSolid,-1,kNone,0,size[i]*0.5);
    leg1->AddEntry(vd[i],Form("Depth %d",i+1),"PLE");
  }

  leg2->AddEntry(phcal2,"H/T (H-had.)","L");
  leg2->AddEntry(pcal,"(E+H)/T (all had.)","L");
  leg2->AddEntry(phcal,"H/T (all had.)","PL");
  leg2->AddEntry(pecal,"E/T (all had.)","PL");  
  leg2->AddEntry(vd2[1],"H_{2}/H (H-had.)","L");
  leg2->AddEntry(vd[1],"H_{2}/H (all had.)","PLE");

  gPad->RedrawAxis();

  c1->SaveAs("pdf/hcalDepthFractions.pdf");

  //h->GetXaxis()->SetRangeUser(1.5,3.0);
  h->GetXaxis()->SetRangeUser(2.0,3.0);
  c1->Update();
  c1->SaveAs("pdf/hcalDepthFractions_HE.pdf");
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  h2sum->Draw("COLZ");
}
