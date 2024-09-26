// Purpose: draw piongun resuts
//          - maps of reconstruction efficiency vs pT, |eta|
//          - confusion rate between E and EH hadrons vs pT, |eta|
//          - power law fits to response vs pT in bins of |eta|
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"

#include "tdrstyle_mod22.C"

#include <fstream>

// Conrado Munoz Diaz's tuple format
bool useCMD = true;

void drawPiongun(string file = "", bool isClosure = false) {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/vsEta");

  TFile *f = new TFile(file=="" ? "piongun.root" : file.c_str(), "READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  // Load results produced by [mk_]piongun.C
  //////////////////////////////////////////
  
  TProfile2D *p2e, *p2h, *p2r_e, *p2r_h, *p2r_a;
  p2e = (TProfile2D*)f->Get("p2e"); assert(p2e);
  p2h = (TProfile2D*)f->Get("p2h"); assert(p2h);
  p2r_e = (TProfile2D*)f->Get("p2r_e"); assert(p2r_e);
  p2r_h = (TProfile2D*)f->Get("p2r_h"); assert(p2r_h);
  p2r_a = (TProfile2D*)f->Get("p2r_a"); assert(p2r_a);

  TProfile2D *p2c_e, *p2c_h, *p2c_a;
  p2c_e = (TProfile2D*)f->Get("p2c_e"); assert(p2c_e);
  p2c_h = (TProfile2D*)f->Get("p2c_h"); assert(p2c_h);
  p2c_a = (TProfile2D*)f->Get("p2c_a"); assert(p2c_a);

  TProfile2D *p2rf_bb, *p2rf_ec1, *p2rf_ec2;
  p2rf_bb = (TProfile2D*)f->Get("p2rf_bb0"); assert(p2rf_bb);
  p2rf_ec1 = (TProfile2D*)f->Get("p2rf_ec1"); assert(p2rf_ec1);
  p2rf_ec2 = (TProfile2D*)f->Get("p2rf_ec2"); assert(p2rf_ec2);
  
  TProfile *pe_bb, *pe_ec1, *pe_ec2;
  pe_bb = (TProfile*)f->Get("pe_bb0"); assert(pe_bb);
  pe_ec1 = (TProfile*)f->Get("pe_ec1"); assert(pe_ec1);
  pe_ec2 = (TProfile*)f->Get("pe_ec2"); assert(pe_ec2);
  
  TH2D *h2rf_bb, *h2rf_ec1, *h2rf_ec2;
  h2rf_bb = (TH2D*)f->Get("h2rf_bb0"); assert(h2rf_bb);
  h2rf_ec1 = (TH2D*)f->Get("h2rf_ec1"); assert(h2rf_ec1);
  h2rf_ec2 = (TH2D*)f->Get("h2rf_ec2"); assert(h2rf_ec2);

  
  // Draw hadron detection efficiency in 2D
  /////////////////////////////////////////
  
  TH1D *h1e = tdrHist("h1e","|#eta_{gen}|",0,3.139,
		      "p_{T,gen} (GeV)",0.2,300);
  extraText = "Private";
  lumi_136TeV = "Winter24 piongun";
  TCanvas *c1e = tdrCanvas("c1e",h1e,8,11,kRectangular);
  gPad->SetLogx();
  gPad->SetRightMargin(0.15);

  p2e->UseCurrentStyle();
  p2e->GetZaxis()->SetTitleOffset(0.85);
  tdrDraw(p2e,"COLZ");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.2,1.479,300,1.479);
  l->DrawLine(0.2,2.500,300,2.500);
  
  gPad->RedrawAxis();
  gPad->Update();
  c1e->SaveAs("pdf/drawPionGun_eff2D.pdf");

  
  // Draw hadron detection efficiency in 1D
  /////////////////////////////////////////
  
  TH1D *h1e1 = tdrHist("h1e1","Efficiency",0,1.2,
		      "p_{T,gen} (GeV)",0.2,300);
  TCanvas *c1e1 = tdrCanvas("c1e1",h1e1,8,11,kSquare);
  gPad->SetLogx();

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1);
  l->DrawLine(2.5,0,2.5,1);
  l->DrawLine(0.7,0,0.7,0.5);
  l->DrawLine(0.2,0.84,300,0.84);
  l->DrawLine(0.2,1.00,300,1.00);
  
  tdrDraw(pe_bb,"Pz",kNone,kBlue);
  tdrDraw(pe_ec1,"Pz",kNone,kGreen+2);
  tdrDraw(pe_ec2,"Pz",kNone,kRed);

  TLegend *leg2e = tdrLeg(0.60,0.65-0.05*3,0.85,0.65);
  leg2e->AddEntry(pe_bb,"BB","PLE");
  leg2e->AddEntry(pe_ec1,"EC1","PLE");
  leg2e->AddEntry(pe_ec2,"EC2","PLE");
  
  gPad->RedrawAxis();
  c1e1->SaveAs("pdf/drawPionGun_eff1D.pdf");


  // Draw E/H identification rate (H-hadron fraction) in 2D
  /////////////////////////////////////////////////////////x
  
  TH1D *h2h = tdrHist("h1h","|#eta_{gen}|",0,3.139,
		      "p_{T,gen} (GeV)",0.2,300);
  TCanvas *c2h = tdrCanvas("c2h",h2h,8,11,kRectangular);
  gPad->SetLogx();
  gPad->SetRightMargin(0.15);

  p2h->UseCurrentStyle();
  p2h->GetZaxis()->SetTitleOffset(0.85);
  p2h->GetZaxis()->SetTitle("H-hadron fraction");
  tdrDraw(p2h,"COLZ");

  TF1 *f1 = new TF1("f2","[0]+[1]*pow(log(x),[2])",1.0,30);
  f1->SetParameters(1.305,0.8,0.5);
  f1->SetLineColor(kGray+2);
  f1->Draw("SAME");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.2,1.479,300,1.479);
  l->DrawLine(0.2,2.500,300,2.500);
  
  gPad->RedrawAxis();
  gPad->Update();
  c2h->SaveAs("pdf/drawPionGun_Hfrac2D.pdf");

  
  // Draw response in 2D for H, EH and A(ll) hadrons
  //////////////////////////////////////////////////
  
  TCanvas *c3r = new TCanvas("c3r","c3r",1800,600);
  c3r->Divide(3,1,0,0);

  c3r->cd(1);
  gPad->SetLogx();

  p2r_e->UseCurrentStyle();
  p2r_e->GetXaxis()->SetMoreLogLabels();
  p2r_e->GetXaxis()->SetNoExponent();
  p2r_e->GetXaxis()->SetRangeUser(0.7,300.);
  p2r_e->GetYaxis()->SetRangeUser(0,3.139);
  p2r_e->GetZaxis()->SetRangeUser(0,1.1);
  if (isClosure) p2r_e->GetZaxis()->SetRangeUser(0.8,1.2);
  tdrDraw(p2r_e,"COL");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.75,0.92,"E-hadrons");
  
  gPad->RedrawAxis();

  c3r->cd(2);
  gPad->SetLogx();

  p2r_h->UseCurrentStyle();
  p2r_h->GetXaxis()->SetMoreLogLabels();
  p2r_h->GetXaxis()->SetNoExponent();
  p2r_h->GetXaxis()->SetRangeUser(0.7,300.);
  p2r_h->GetYaxis()->SetRangeUser(0,3.139);
  p2r_h->GetZaxis()->SetRangeUser(0,1.1);
  if (isClosure) p2r_h->GetZaxis()->SetRangeUser(0.8,1.2);
  tdrDraw(p2r_h,"COL");

  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.70,0.92,"H-hadrons");
  
  gPad->RedrawAxis();

  c3r->cd(3);
  gPad->SetLogx();
  gPad->SetRightMargin(0.20);

  p2r_a->UseCurrentStyle();
  p2r_a->GetXaxis()->SetMoreLogLabels();
  p2r_a->GetXaxis()->SetNoExponent();
  p2r_a->GetXaxis()->SetRangeUser(0.7,300.);
  p2r_a->GetYaxis()->SetRangeUser(0,3.139);
  p2r_a->GetZaxis()->SetRangeUser(0,1.1);
  p2r_a->GetZaxis()->SetTitle("Response");
  if (isClosure) {
    p2r_a->GetZaxis()->SetRangeUser(0.8,1.2);
    p2r_a->GetZaxis()->SetTitle("Response closure");
  }
  tdrDraw(p2r_a,"COLZ");
  
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");

  tex->DrawLatex(0.52,0.92,"All hadrons");

  gPad->RedrawAxis();
  
  c3r->SaveAs("pdf/drawPionGun_resp2D_x3.pdf");


  // Draw (correction)^{-1} in 2D for H, EH and A(ll) hadrons
  // (Sanity check that it agrees with original response)
  ///////////////////////////////////////////////////////////
  
  TCanvas *c3c = new TCanvas("c3c","c3c",1800,600);
  c3c->Divide(3,1,0,0);

  c3c->cd(1);
  gPad->SetLogx();

  p2c_e->UseCurrentStyle();
  p2c_e->GetXaxis()->SetMoreLogLabels();
  p2c_e->GetXaxis()->SetNoExponent();
  p2c_e->GetXaxis()->SetRangeUser(0.7,300.);
  p2c_e->GetYaxis()->SetRangeUser(0,3.139);
  p2c_e->GetZaxis()->SetRangeUser(0,1.1);
  tdrDraw(p2c_e,"COL");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.75,0.92,"E-hadrons");
  
  gPad->RedrawAxis();

  c3c->cd(2);
  gPad->SetLogx();

  p2c_h->UseCurrentStyle();
  p2c_h->GetXaxis()->SetMoreLogLabels();
  p2c_h->GetXaxis()->SetNoExponent();
  p2c_h->GetXaxis()->SetRangeUser(0.7,300.);
  p2c_h->GetYaxis()->SetRangeUser(0,3.139);
  p2c_h->GetZaxis()->SetRangeUser(0,1.1);
  tdrDraw(p2c_h,"COL");

  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.70,0.92,"H-hadrons");
  
  gPad->RedrawAxis();

  c3c->cd(3);
  gPad->SetLogx();
  gPad->SetRightMargin(0.20);

  p2c_a->UseCurrentStyle();
  p2c_a->GetXaxis()->SetMoreLogLabels();
  p2c_a->GetXaxis()->SetNoExponent();
  p2c_a->GetXaxis()->SetRangeUser(0.7,300.);
  p2c_a->GetYaxis()->SetRangeUser(0,3.139);
  p2c_a->GetZaxis()->SetRangeUser(0,1.1);
  p2c_a->GetZaxis()->SetTitle("Correction^{-1}");
  tdrDraw(p2c_a,"COLZ");
  
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  f1->Draw("SAME");

  tex->DrawLatex(0.52,0.92,"All hadrons");

  gPad->RedrawAxis();
  
  c3c->SaveAs("pdf/drawPionGun_invcorr2D_x3.pdf");

  
  // Draw response in 2D fECAL vs pT for BB, EC1 and EC2
  //////////////////////////////////////////////////////
  
  TCanvas *c4f = new TCanvas("c4f","c4f",1800,600);
  c4f->Divide(3,1,0,0);

  c4f->cd(1);
  gPad->SetLogx();

  p2rf_bb->UseCurrentStyle();
  p2rf_bb->GetXaxis()->SetMoreLogLabels();
  p2rf_bb->GetXaxis()->SetNoExponent();
  p2rf_bb->GetXaxis()->SetRangeUser(0.7,300.);
  p2rf_bb->GetYaxis()->SetRangeUser(0,3.139);
  p2rf_bb->GetZaxis()->SetRangeUser(0,1.1);
  if (isClosure) p2rf_bb->GetZaxis()->SetRangeUser(0.8,1.2);
  tdrDraw(p2rf_bb,"COL");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.);
  
  TF1 *f1mip = new TF1("f1mip","max([0]/x,[1])",0.7,300.);
  f1mip->SetParameters(1.,0.01);
  f1mip->SetLineColor(kGray+2);
  f1mip->Draw("SAME");

  tex->DrawLatex(0.75,0.92,"BB hadrons");  

  gPad->RedrawAxis();

  c4f->cd(2);
  gPad->SetLogx();

  p2rf_ec1->UseCurrentStyle();
  p2rf_ec1->GetXaxis()->SetMoreLogLabels();
  p2rf_ec1->GetXaxis()->SetNoExponent();
  p2rf_ec1->GetXaxis()->SetRangeUser(0.7,300.);
  p2rf_ec1->GetYaxis()->SetRangeUser(0,3.139);
  p2rf_ec1->GetZaxis()->SetRangeUser(0,1.1);
  if (isClosure) p2rf_ec1->GetZaxis()->SetRangeUser(0.8,1.2);
  tdrDraw(p2rf_ec1,"COL");

  l->DrawLine(2.5,0,2.5,1.);
  f1mip->Draw("SAME");
  
  tex->DrawLatex(0.70,0.92,"EC1 hadrons");
  
  gPad->RedrawAxis();

  c4f->cd(3);
  gPad->SetLogx();
  gPad->SetRightMargin(0.20);

  p2rf_ec2->UseCurrentStyle();
  p2rf_ec2->GetXaxis()->SetMoreLogLabels();
  p2rf_ec2->GetXaxis()->SetNoExponent();
  p2rf_ec2->GetXaxis()->SetRangeUser(0.7,300.);
  p2rf_ec2->GetYaxis()->SetRangeUser(0,3.139);
  p2rf_ec2->GetZaxis()->SetRangeUser(0,1.1);
  p2rf_ec2->GetZaxis()->SetTitle("Response");
  if (isClosure) {
    p2rf_ec2->GetZaxis()->SetRangeUser(0.8,1.2);
    p2rf_ec2->GetZaxis()->SetTitle("Response closure");
  }
  tdrDraw(p2rf_ec2,"COLZ");

  l->DrawLine(2.5,0,2.5,1.);
  f1mip->Draw("SAME");
  
  tex->DrawLatex(0.52,0.92,"EC2 hadrons");

  gPad->RedrawAxis();
  
  c4f->SaveAs("pdf/drawPionGun_fecal2D_x3.pdf");


  // Draw response vs fECAL in 1D (for barrel)
  ////////////////////////////////////////////
  
  TH1D *h4f1 = tdrHist("h4f1","(Fraction or) Response",
		       0.0, isClosure ? 1.7 : 1.3,
		       "f_{ECAL,raw}",0.,1.0);
  if (isClosure) h4f1->SetYTitle("(Fraction or) Response closure");
  TCanvas *c4f1 = tdrCanvas("c4f1",h4f1,8,11,kSquare);
  gPad->RedrawAxis();

  TF1 *f1rf = new TF1("f1rf","[0]*(1-x)+[1]*x",0,1);
  TF1 *f1rfm = new TF1("f1rfm","[0]+[1]*(x-0.5)",0,1);
  
  double vx[] = {5, 10, 20, 40, 80, 160, 320};
  const int nx = sizeof(vx)/sizeof(vx[0]);
  int color[] = {kBlue, kGreen+2, kYellow+2, kOrange+1, kRed, kBlack, kGray+1};
  const int nc = sizeof(color)/sizeof(color[0]);
  
  TLegend *leg4f1 = tdrLeg(0.35,0.90-0.035*nx,0.60,0.90);
  leg4f1->SetTextSize(0.035);

  tex->DrawLatex(0.68,0.85,"Barrel:");
  tex->DrawLatex(0.68,0.80,"Point Resp.,");
  tex->DrawLatex(0.68,0.75,"Hist Frac.");
  tex->DrawLatex(0.68,0.70,"[0.1,0.9]#times10");
  
  for (int i = 0; i != nx; ++i) {
    double pt = vx[i];
    int j = p2rf_bb->GetXaxis()->FindBin(pt);
    double ptmin = p2rf_bb->GetXaxis()->GetBinLowEdge(j);
    double ptmax = p2rf_bb->GetXaxis()->GetBinLowEdge(j+1);
    TH1D *h = p2rf_bb->ProjectionY(Form("h%1.0f",pt),j,j,"o");

    tdrDraw(h,"Pz",kNone,color[i]);

    double fmin = max(f1mip->Eval(pt)*1.5,0.10);
    double fmax = min(0.80,1-2.0/pt);
    f1rf->SetRange(fmin,fmax);
    f1rf->SetLineColor(color[i]);

    f1rf->SetParameters(0.75,0.60);
    h->Fit(f1rf,"QRN");

    f1rf->DrawClone("SAME");
    f1rf->SetRange(0,1);
    f1rf->SetLineStyle(kDotted);
    f1rf->DrawClone("SAME");
    f1rf->SetLineStyle(kSolid);
    
    f1rfm->SetRange(fmin,fmax);
    f1rfm->SetLineColor(color[i]);
    f1rf->SetParameters(0.675,-0.075);
    h->Fit(f1rfm,"QRN");
    f1rfm->DrawClone("SAME");

    TH1D *hf = h2rf_bb->ProjectionY(Form("hf%1.0f",pt),j,j,"o");
    hf->Scale(1./hf->Integral());
    tdrDraw(hf,"HIST",kNone,color[i],kSolid,-1,1001,color[i]);
    hf->SetFillColorAlpha(color[i],0.1);

    for (int k = 1; k != hf->GetNbinsX()+1; ++k) {
      if (hf->GetBinWidth(k)>0.02) {
	hf->SetBinContent(k, hf->GetBinContent(k)*10);
      }
    }
    leg4f1->AddEntry(hf,Form("[%1.3g,%1.3g] GeV",ptmin,ptmax),"PLEF");
  } // for i

  gPad->RedrawAxis();
  c4f1->SaveAs("pdf/drawPionGun_fecal1D.pdf");

  /*
  // Repeat systematically for all the pT bins, plot [0] and [1]
  TH1D *hre = p2rf_bb->ProjectionX("hre",0,-1,"o"); hre->Reset();
  TH1D *hrh = p2rf_bb->ProjectionX("hrh",0,-1,"o"); hrh->Reset();
  TH1D *hr = p2rf_bb->ProjectionX("hr",0,-1,"o"); hr->Reset();
  for (int i = 1; i != p2rf_bb->GetNbinsX()+1; ++i) {
    double pt = p2rf_bb->GetXaxis()->GetBinLowEdge(i);
    TH1D *hrf = p2rf_bb->ProjectionY(Form("hrf%1.0f",pt),i,i,"o");
    
    double fmin = max(f1mip->Eval(pt)*1.5,0.10);
    double fmax = min(0.80,1-2.0/pt);
    f1rf->SetRange(fmin,fmax);
    f1rf->SetParameters(0.75,-0.15);
    f1rfm->SetRange(fmin,fmax);
    f1rfm->SetParameters(0.675,-0.075);
    if (hrf->Integral()>0) {
      hrf->Fit(f1rf,"QRN");
      hrf->Fit(f1rfm,"QRN");

      double k = sqrt(max(1.,f1rf->GetChisquare()/max(1,f1rf->GetNDF())));
      hrh->SetBinContent(i, f1rf->GetParameter(0));
      hrh->SetBinError(i, k*f1rf->GetParError(0));
      hre->SetBinContent(i, f1rf->GetParameter(1));
      hre->SetBinError(i, k*f1rf->GetParError(1));
      hr->SetBinContent(i, f1rfm->GetParameter(0));
      hr->SetBinError(i, k*f1rfm->GetParError(0));
    }
  }

  TProfile *prhh = p2rf_bb->ProfileX("prhh",1,1,"o"); // ref.
  int ieta05 = p2r_h->GetYaxis()->FindBin(0.522-0.05); // def.
  TProfile *prhh2 = p2r_h->ProfileX("prhh2",1,ieta05,"o"); // def.
  TProfile *pree2 = p2r_e->ProfileX("pree2",1,ieta05,"o"); // def.
  int j1 = p2rf_bb->GetYaxis()->FindBin(0.20);
  int j2 = p2rf_bb->GetYaxis()->FindBin(0.80);
  TProfile *pree = p2rf_bb->ProfileX("pree",j1,j2,"o"); // ref.
  
  TH1D *h4rf = tdrHist("h4rf","Response",0,1.3,"p_{T,gen} (GeV)",0.2,1000.);
  if (isClosure) h4rf->SetYTitle("Corrected response");
  if (isClosure) h4rf->GetYaxis()->SetRangeUser(0.8+1e-5,1.5-1e-5);
  TCanvas *c4rf = tdrCanvas("c4rf",h4rf,8,11,kSquare);
  gPad->SetLogx();

  tdrDraw(prhh,"Pz",kFullCircle,kRed, kSolid,-1,kNone,0,0.7);
  tdrDraw(prhh2,"Pz",kOpenCircle,kRed, kSolid,-1,kNone,0,0.7);
  tdrDraw(pree,"Pz",kFullCircle,kBlue, kSolid,-1,kNone,0,0.7);
  tdrDraw(pree2,"Pz",kOpenCircle,kBlue, kSolid,-1,kNone,0,0.7);

  tdrDraw(hr,"Pz",kNone,kBlack);
  tdrDraw(hre,"Pz",kNone,kBlue);
  tdrDraw(hrh,"Pz",kNone,kOrange+1);


  double emin = 3.5;
  double eref = 500;
  double etaref = 0.5;
  cout << "H-component of EH" << endl;
  TF1 *f1rh = new TF1("f1rh","[0]+[1]*pow(x,[2])",emin,eref/cosh(etaref));
  f1rh->SetParameters(1,-1,-0.3);
  f1rh->SetParLimits(0,0.9,1.2);
  f1rh->SetParLimits(2,-0.5,-0.1);
  hrh->Fit(f1rh,"RN");
  f1rh->SetLineColor(kOrange+1);
  f1rh->DrawClone("SAME");
  f1rh->SetLineStyle(kDotted);
  f1rh->SetRange(0.2,1000.);
  f1rh->Draw("SAME");

  cout << "E-component of EH" << endl;
  TF1 *f1re = new TF1("f1re","[0]+[1]*pow(x,[2])",emin,eref/cosh(etaref));
  f1re->SetParameters(1,-1,-0.3);
  f1re->SetParLimits(0,0.9,1.3);//1.2);
  f1re->SetParLimits(2,-0.5,-0.1);
  hre->Fit(f1re,"RN");
  f1re->SetLineColor(kBlue);
  f1re->DrawClone("SAME");
  f1re->SetLineStyle(kDotted);
  f1re->SetRange(0.2,1000.);
  f1re->Draw("SAME");

  cout << "EH-hadrons" << endl;
  TF1 *f1r = new TF1("f1r","[0]+[1]*pow(x,[2])",emin,eref/cosh(etaref));
  f1r->SetParameters(1,-1,-0.3);
  f1r->SetParLimits(0,0.9,1.2);
  f1r->SetParLimits(2,-0.5,-0.1);
  hr->Fit(f1r,"RN");
  f1r->SetLineColor(kBlack);
  f1r->DrawClone("SAME");
  f1r->SetLineStyle(kDotted);
  f1r->SetRange(0.2,1000.);
  f1r->Draw("SAME");

  cout << "H-hadrons" << endl;
  TF1 *f1rhh = new TF1("f1rhh","[0]+[1]*pow(x,[2])",emin,eref/cosh(etaref));
  f1rhh->SetParameters(1,-1,-0.3);
  f1rhh->SetParLimits(0,0.9,1.2);
  f1rhh->SetParLimits(2,-0.5,-0.1);
  prhh->Fit(f1rhh,"RN");
  f1rhh->SetLineColor(kRed);
  f1rhh->DrawClone("SAME");
  f1rhh->SetLineStyle(kDotted);
  f1rhh->SetRange(0.2,1000.);
  f1rhh->Draw("SAME");
  
  gPad->RedrawAxis();
  c4rf->SaveAs("pdf/drawPionGun_respHE.pdf");
*/

  ///////////////////////////////////////////////////////////
  // Complete 3D analysis in eta(x), pT(y), fe(z); 30 bins //
  ///////////////////////////////////////////////////////////
  
  TCanvas *c5 = new TCanvas("c5","c5",6*300,5*300);
  c5->Divide(6,5,0,0);
  TCanvas *c5f = new TCanvas("c5f","c5f",6*300,5*300);
  c5f->Divide(6,5,0,0);
  const int neta = 6*5;

  // Set limits to fit parameters c,a,m (also used for automatic plot ranges)
  // Functional form is c*(1 - a*pT^{m-1})
  // Expectations are c~1, a~0.5 (but seems a~1) and m~0.7 (or ~0.85?)
  double minc = (isClosure ? 0.80 : 0.90);
  double maxc = (isClosure ? 1.20 : 1.40);
  double mina = (isClosure ?   -1 : 0.45);
  double maxa = (isClosure ?   +1 : 1.20);
  double minm = 0.0;
  double maxm = 0.90;

  double refa = (isClosure ? 0 : 1) ;
  
  c5f->cd(neta);
  TLegend *leg5f = tdrLeg(0.05,0.90,0.55,0.90);

  // Load full eta,pT,f_ECAL 3D map of single-pion response
  TProfile3D *p3 = (TProfile3D*)f->Get("p3rf"); assert(p3);
  
  // Map to store graphs of the fit results
  map<string, map<int, TGraphErrors*> > mg;

  // Vector to store fit results
  vector< map<string, double[4]> > vm(p3->GetNbinsX());
  
  // Fit shape vs f_ECAL: now linear, but maybe quadratic in the future
  // (HCAL+ECAL response is slightly banana shaped, esp. at high pT)
  TF1 *f1f = new TF1("f1f","[0]*(1-x)+[1]*x",0,1);

  for (int ieta = 1; ieta != p3->GetNbinsX()+1; ++ieta) {
    double eta = p3->GetXaxis()->GetBinCenter(ieta);
    double deta = 0.5*p3->GetXaxis()->GetBinWidth(ieta);

    // There are hadrons up to (eta>3.139), but these are partially in HF
    if (eta>2.964) continue;

    // Select eta bin
    p3->GetXaxis()->SetRange(ieta,ieta);
    TProfile2D *p2 = p3->Project3DProfile("zy");
    p2->SetName(Form("p2yz_%d",ieta));
					  
    // All, H-hadrons and "average" EH-hadrons directly vs pT
    // Need "o" (=original) option to not mess up variable binned axis

    // Range 0,-1 projects all fe bins
    TProfile *pa = p2->ProfileX(Form("pa_%d",ieta),0,-1,"o");

    // Range [1,1] effectively picks fe<0.01, so below MIP cut up to 100 GeV
    // at higher pT still useful to cut out delta rays(?) above MIP levels
    // there is strong fe dependence at low fe due to H-EH mixture so cut tight
    TProfile *ph = p2->ProfileX(Form("ph_%d",ieta),1,1,"o");

    // Use less biased EH hadrons in the central range fe [0,2.0.8]:
    // 0<fe<0.2 is mix of H and EH where E was lost, and very variable vs fe
    // fe>0.8 is sometimes ok, sometimes with H that was lost, so not linear
    int i1 = p2->GetYaxis()->FindBin(0.2);
    int i2 = p2->GetYaxis()->FindBin(0.8);
    TProfile *pe = p2->ProfileX(Form("pe_%d",ieta),i1,i2,"o");

    // EH-hadrons differentially vs fe. Hand-pick linear region E>MIP+[0.2,0.8]
    ///////////////////////////////////////////////////////////////////////////
    
    c5f->cd(ieta);
    TH1D *hf = tdrHist(Form("h5f_%d",ieta),"(rawECAL+rawHCAL)/genP",
		       0.2+1e-4,1.2-1e-5,
		       "rawEcal/(rawEcal+rawHcal)",0,1);
    if (useCMD) hf->SetYTitle("(ecal+hcal)/true");
    if (useCMD) hf->SetXTitle("ecal/(ecal+hcal)");
    if (useCMD && isClosure) hf->SetYTitle("corrected (ecal+hcal)/true");
    if (useCMD && isClosure) hf->SetXTitle("raw ecal/(ecal+hcal)");
    if (isClosure) hf->GetYaxis()->SetRangeUser(0.8+1e-5,1.5-1e-5);
    hf->Draw();

    tex->SetTextSize(0.045*1.5);
    tex->DrawLatex(0.50,0.87,Form("%1.3f#leq|#eta|<%1.3f",eta-deta,eta+deta));

    TH1D *hee = p2->ProjectionX(Form("hee_%d",ieta),i1,i2,"o");
    hee->Reset();
    TH1D *heh = p2->ProjectionX(Form("heh_%d",ieta),i1,i2,"o");
    heh->Reset();
    TH1D *hea = p2->ProjectionX(Form("hea_%d",ieta),i1,i2,"o");
    hea->Reset();
    f1f->SetParameters(0.90,0.90); // high pT starting guess

    // Reverse loop so lower pT uses better previous fit values as starter
    for (int ipt = p2->GetNbinsX(); ipt !=0; --ipt) {

      double pt = p2->GetXaxis()->GetBinCenter(ipt);
      double ptmin = p2->GetXaxis()->GetBinLowEdge(ipt);
      double ptmax = p2->GetXaxis()->GetBinLowEdge(ipt+1);
	    
      TProfile *pfe =p2->ProfileY(Form("pfe_%d_%d",ieta,ipt),ipt,ipt,"o");
      double e_mip = 1. * 1.5; // GeV
      double h_thr = 2.0; // GeV
      f1f->SetRange(min(max(0.1,e_mip/pt),0.3), max(0.7,min(0.9,1-h_thr/pt)));
      int i50 = pfe->GetXaxis()->FindBin(0.50);
      if (pfe->GetBinError(i50)!=0) {

	TFitResultPtr fp1f = pfe->Fit(f1f,"QRNS"); // S to return fit result

	// Get the covariance matrix
	TMatrixDSym covMatrix = fp1f->GetCovarianceMatrix();

	// Get the uncertainties on the parameters
	double sigma_p0 = sqrt(covMatrix(0, 0));
	double sigma_p1 = sqrt(covMatrix(1, 1));
	double cov_p0p1 = covMatrix(0, 1);

	// Set the value of x
	double x = 0.5;

	// Derivatives of the function at x = 0.5
	double df_dp0 = (1 - x);
	double df_dp1 = x;

	// Error propagation (to x=0.5)
	double error = sqrt( df_dp0 * df_dp0 * sigma_p0 * sigma_p0 +
			     df_dp1 * df_dp1 * sigma_p1 * sigma_p1 +
			     2 * df_dp0 * df_dp1 * cov_p0p1 );
	double errmin = 0;//0.001;
	
	double rh = f1f->GetParameter(0);
	double re = f1f->GetParameter(1);

	if (f1f->GetNDF()>0) {
	  double k  = sqrt(f1f->GetChisquare()/max(1,f1f->GetNDF()));
	  heh->SetBinContent(ipt, f1f->GetParameter(0));
	  heh->SetBinError(ipt, max(errmin,k*f1f->GetParError(0)));
	  hee->SetBinContent(ipt, f1f->GetParameter(1));
	  hee->SetBinError(ipt, max(errmin,k*f1f->GetParError(1)));

	  // 50-50 point for better 'm' later
	  hea->SetBinContent(ipt, 0.5*(rh+re));
	  hea->SetBinError(ipt, max(errmin,k*error));
	} // good fit
	
	// Draw fe for a few reference pTs starting from 5 to 500 GeV
	int ipt5 = p2->GetXaxis()->FindBin(5.);
	if (pt>=5 && pt <=500 && (ipt-ipt5)%4==0) {
	  c5f->cd(ieta);
	    
	  tdrDraw(pfe,"Pz",kNone,color[((ipt-ipt5)/4)%nc],kSolid,-1);
	  if (ieta==1) {
	    leg5f->AddEntry(pfe,Form("[%1.0f,%1.0f] GeV",ptmin,ptmax),"PLE");
	    leg5f->SetY1NDC(leg5f->GetY1NDC()-0.05*1.5);
	  }
	  
	  if (f1f->GetNDF()>0) {
	    f1f->SetLineStyle(kSolid);
	    f1f->SetLineColor(color[(ipt-ipt5)/4%nc]);
	    f1f->DrawClone("SAME");
	    f1f->SetLineStyle(kDotted);
	    f1f->SetRange(0,1);
	    f1f->DrawClone("SAME");
	  } // good fit
	} // find in setipt
	
      } // pfe50!=0

    } // for ipt

    // A(ll), H, EH and E/H components of HE vs pT
    //////////////////////////////////////////////
    
    c5->cd(ieta);
    gPad->SetLogx();

    TH1D *h = tdrHist(Form("h5_%d",ieta),"(rawECAL+rawHCAL)/genP",0+1e-5,1.3,
		      "p_{T,gen} (GeV)",0.2,1000.-1e-3);
    if (useCMD) h->SetYTitle("(ecal+hcal)/true");
    if (useCMD) h->SetXTitle("true/cosh(eta) (GeV)");
    if (useCMD && isClosure) h->SetYTitle("corrected (ecal+hcal)/true");
    if (isClosure) h->GetYaxis()->SetRangeUser(0.8+1e-5,1.5-1e-5);
    h->Draw();

    //TLatex *tex = new TLatex();
    tex->SetTextSize(0.045*1.5);
    tex->SetNDC();

    tex->DrawLatex(0.50,0.87,Form("%1.3f#leq|#eta|<%1.3f",eta-deta,eta+deta));
    
    tdrDraw(hea,"Pz",kOpenCircle,kBlue-9, kSolid,-1,kNone,0, 0.7);
    tdrDraw(pa,"Pz",kOpenCircle,kGray+1, kSolid,-1,kNone,0, 0.7);
    tdrDraw(pe,"Pz",kFullCircle,kBlue-9, kSolid,-1,kNone,0, 0.7);

    tdrDraw(ph,"Pz",kFullCircle,kRed, kSolid,-1,kNone,0, 0.7);
    tdrDraw(hee,"Pz",kOpenDiamond,kMagenta+2, kSolid,-1,kNone,0, 1.0);
    tdrDraw(heh,"Pz",kOpenDiamond,kOrange+2, kSolid,-1,kNone,0, 1.0);

    double ptmin = 5.; // broadly safe, expect parts of EC
    double ptmax = 500./cosh(eta);
    double fixm_a = 0;//0.60; // 0 for free
    double fixm_h = 0;//0.70; // 0 for free
    double fixm_e = 0;//0.70; // 0 for free
    bool fixm_toe = false;//true;
    bool fixa = false;//true;
    
    // Adjust minimum pT range for E hadrons
    double ptmin_e = ptmin;
    //if (eta<1.479) ptmin_e = 3.5; // go bit lower
    if (eta>2.043) ptmin_e = 8; // go bit higher
    
    // NB: parameter 'a' should be vs E, but works better with pT??

    // All hadrons
    TFitResultPtr fp1a;
    TF1 *f1a = new TF1(Form("f1a_%d",ieta),
		       "max([3],[0]*(1-[1]*pow(x,[2]-1)))",
		       ptmin,ptmax);
    f1a->SetParameters(1,refa,0.60,0.50);
    f1a->SetParLimits(0,minc,maxc);
    f1a->SetParLimits(1,mina,maxa);
    if (fixa) f1a->FixParameter(1,1.0);
    f1a->SetParLimits(2,minm,maxm);
    if (fixm_a) f1a->FixParameter(2,fixm_a);
    f1a->FixParameter(3,0.50);

    f1a->SetLineColor(kGray+1);
    fp1a = pa->Fit(f1a,"QRNS");
    TF1 *f1a_0 = (TF1*)f1a->DrawClone("SAME");
    f1a->SetRange(0.2,1000.);
    f1a->SetLineStyle(kDotted);
    //f1a->DrawClone("SAME");

    // H-hadrons
    TF1 *f1h = new TF1(Form("f1h_%d",ieta),
		       "max([3],[0]*(1-[1]*pow(x,[2]-1)))",
		       ptmin,ptmax);
    f1h->SetParameters(1.1,refa,0.70,0.65);
    f1h->SetParLimits(0,minc,maxc);
    f1h->SetParLimits(1,mina,maxa);
    if (fixa) f1h->FixParameter(1,0.55);
    f1h->SetParLimits(2,minm,maxm);
    if (fixm_h && eta<1.479) f1h->FixParameter(2,fixm_h);
    f1h->FixParameter(3,0.65);
    if (eta>1.392) f1h->FixParameter(3,0.75);
    if (eta>1.740) f1h->FixParameter(3,0.65);
    if (eta>2.172) f1h->FixParameter(3,0.60);
    if (eta>2.500) f1h->FixParameter(3,0.50);
    if (eta>2.853) f1h->FixParameter(3,0.40);

    f1h->SetLineColor(kRed);
    ph->Fit(f1h,"QRNS");
    TF1 *f1h_0 = (TF1*)f1h->DrawClone("SAME");
    f1h->SetRange(0.2,1000.);
    f1h->SetLineStyle(kDotted);
    f1h->DrawClone("SAME");

    TF1 *f1e = new TF1(Form("f1e_%d",ieta),
		       "max([3],[0]*(1-[1]*pow(x,[2]-1)))",
		       ptmin_e,ptmax);
    f1e->SetParameters(1.1,refa,0.70,0.50);
    f1e->SetParLimits(0,minc,maxc);
    f1e->SetParLimits(1,mina,maxa);
    if (fixa) f1e->FixParameter(1,1.0);
    f1e->SetParLimits(2,minm,maxm);
    if (fixm_e) f1e->FixParameter(2,fixm_e);
    f1e->FixParameter(3,0.50);

    f1e->SetLineColor(kBlue-9);
    pe->Fit(f1e,"QRNS");
    TF1 *f1e_0 = (TF1*)f1e->DrawClone("SAME");
    f1e->SetRange(0.2,1000.);
    f1e->SetLineStyle(kDotted);
    //f1e->DrawClone("SAME");

    TF1 *f1ea = new TF1(Form("f1ea_%d",ieta),
			"max([3],[0]*(1-[1]*pow(x,[2]-1)))",
			ptmin_e,ptmax);
    f1ea->SetParameters(1.1,refa,0.70,0.50);
    f1ea->SetParLimits(0,minc,maxc);
    f1ea->SetParLimits(1,mina,maxa);
    if (fixa) f1ea->FixParameter(1,1.0);
    f1ea->SetParLimits(2,minm,maxm);
    if (fixm_e) f1ea->FixParameter(2,fixm_e);
    f1ea->FixParameter(3,0.40);
    if (eta>2.322) f1ea->FixParameter(3,0.50);

    f1ea->SetLineColor(kBlue-9);
    hea->Fit(f1ea,"QRNS");
    TF1 *f1ea_0 = (TF1*)f1ea->DrawClone("SAME");
    f1ea->SetRange(0.2,1000.);
    f1ea->SetLineStyle(kDotted);
    //f1ea->DrawClone("SAME");

    TF1 *f1ee = new TF1(Form("f1ee_%d",ieta),
			"max([3],[0]*(1-[1]*pow(x,[2]-1)))",
			ptmin,ptmax);
    f1ee->SetParameters(1.15,refa,0.70,0.35);
    f1ee->SetParLimits(0,minc,maxc);
    f1ee->SetParLimits(1,mina,maxa);
    if (fixa) f1ee->FixParameter(1,1.3);
    f1ee->SetParLimits(2,minm,maxm);
    if (fixm_toe) f1ee->FixParameter(2,f1ea->GetParameter(2));
    f1ee->FixParameter(3,0.35);
    if (eta>1.479) f1ee->FixParameter(3,0.25);

    f1ee->SetLineColor(kMagenta+2);
    hee->Fit(f1ee,"QRNS");
    TF1 *f1ee_0 = (TF1*)f1ee->DrawClone("SAME");
    f1ee->SetRange(0.2,1000.);
    f1ee->SetLineStyle(kDotted);
    f1ee->DrawClone("SAME");

    TF1 *f1eh = new TF1(Form("f1eh_%d",ieta),
			"max([3],[0]*(1-[1]*pow(x,[2]-1)))",
			ptmin,ptmax);
    f1eh->SetParameters(1.0,refa,0.70,0.50);
    f1eh->SetParLimits(0,minc,maxc);
    f1eh->SetParLimits(1,mina,maxa);
    if (fixa) f1eh->FixParameter(1,0.85);
    f1eh->SetParLimits(2,minm,maxm);
    if (fixm_toe) f1eh->FixParameter(2,f1ea->GetParameter(2));
    f1eh->FixParameter(3,0.50);
    if (eta>1.479) f1eh->FixParameter(3,0.70);
    if (eta>1.566) f1eh->FixParameter(3,0.60);
    if (eta>2.322) f1eh->FixParameter(3,0.65);
    if (eta>2.500) f1eh->FixParameter(3,0.70);
    if (eta>2.650) f1eh->FixParameter(3,0.65);
    if (eta>2.853) f1eh->FixParameter(3,0.50);

    f1eh->SetLineColor(kOrange+2);
    heh->Fit(f1eh,"QRNS");
    TF1 *f1eh_0 = (TF1*)f1eh->DrawClone("SAME");
    f1eh->SetRange(0.2,1000.);
    f1eh->SetLineStyle(kDotted);
    f1eh->DrawClone("SAME");

    gPad->RedrawAxis();

    // Draw legend in the last empty pad
    TLegend *legi(0);
    if (ieta==1) {
      c5->cd(neta);

      TLegend *leg = tdrLeg(0.05,0.95-2*10*0.045,0.80,0.95);
      leg->SetTextSize(1.5*0.045);
      leg->AddEntry(pa,"All hadrons","PLE");
      leg->AddEntry(pa,"Split by f_{ECAL}:","");
      leg->AddEntry(ph,"H (f_{ECAL}<0.01)","PLE");
      leg->AddEntry(pe,"E (0.20<f_{ECAL}<0.80)","PLE");
      leg->AddEntry(pe,"Fit H#times(1-f_{ECAL})+E#timesf_{ECAL}:","");
      leg->AddEntry(hea,"E (f_{ECAL}#rightarrow0.5)","PLE");
      leg->AddEntry(heh,"H of E (f_{ECAL}#rightarrow0)","PLE");
      leg->AddEntry(hee,"E of E (f_{ECAL}#rightarrow1)","PLE");
      leg->AddEntry(f1a,"Fit c #times (1 - a #times p_{T}^{m-1})","L");
      leg->AddEntry(f1a,"(p_{T}>5 GeV to E<500 GeV)","");

      legi = (TLegend*)leg->Clone();
      legi->SetY1NDC(0.90);
      legi->SetY2NDC(0.90-0.05*9);
      legi->SetX1NDC(0.20);
      legi->SetX2NDC(0.45);
    }

    // Redraw fe-dependence after pT-fits as a sanity check
    TF1 *f1f2 = new TF1("f1f2","[0]*(1-x)+[1]*x",0,1);
    for (int ipt = p2->GetNbinsX(); ipt !=0; --ipt) {

      double pt = p2->GetXaxis()->GetBinCenter(ipt);

	
      // Draw fe for a few reference pTs
      int ipt5 = p2->GetXaxis()->FindBin(5.);
      if (pt>=5 && pt <=500 && (ipt-ipt5)%4==0) {
	c5f->cd(ieta);

	if (hee->GetBinContent(hee->FindBin(pt))>0.) {
	  f1f2->SetParameters(f1eh->Eval(pt), f1ee->Eval(pt));
	  f1f2->SetLineStyle(kSolid);
	  f1f2->SetLineColor(color[(ipt-ipt5)/4%nc] + 1);
	  f1f2->SetLineStyle(kDashed);
	  f1f2->DrawClone("SAME");
	}
      } // find in setipt
	
    } // for ipt

    
    // Store fit results for later re-analysis
    // Shift |eta| a bit in the graphs to avoid error bars overlapping
    double dh = -0.01;
    double de = +0.01;
    double deh = -0.005;
    double dee = +0.005;

    for (int i = 0; i != 4; ++i) {
      
      if (mg["h"][i]==0) mg["h"][i] = new TGraphErrors(0);
      if (mg["e"][i]==0) mg["e"][i] = new TGraphErrors(0);
      if (mg["a"][i]==0) mg["a"][i] = new TGraphErrors(0);
      if (mg["ea"][i]==0) mg["ea"][i] = new TGraphErrors(0);
      if (mg["ee"][i]==0) mg["ee"][i] = new TGraphErrors(0);
      if (mg["eh"][i]==0) mg["eh"][i] = new TGraphErrors(0);
      if (i<3) {
	mg["h"][i]->SetPoint(ieta-1, eta+dh, f1h->GetParameter(i));
	mg["h"][i]->SetPointError(ieta-1, deta, f1h->GetParError(i));
	mg["e"][i]->SetPoint(ieta-1, eta+de, f1e->GetParameter(i));
	mg["e"][i]->SetPointError(ieta-1, deta, f1e->GetParError(i));
	mg["a"][i]->SetPoint(ieta-1, eta, f1a->GetParameter(i));
	mg["a"][i]->SetPointError(ieta-1, deta, f1a->GetParError(i));
	
	mg["ea"][i]->SetPoint(ieta-1, eta+dee, f1ea->GetParameter(i));
	mg["ea"][i]->SetPointError(ieta-1, deta, f1ea->GetParError(i));
	mg["ee"][i]->SetPoint(ieta-1, eta+dee, f1ee->GetParameter(i));
	mg["ee"][i]->SetPointError(ieta-1, deta, f1ee->GetParError(i));
	mg["eh"][i]->SetPoint(ieta-1, eta+deh, f1eh->GetParameter(i));
	mg["eh"][i]->SetPointError(ieta-1, deta, f1eh->GetParError(i));

	vm[ieta-1]["eta"][i] = (eta-deta+i*deta);
	vm[ieta-1]["hh"][i] = f1h->GetParameter(i);
	vm[ieta-1]["eh"][i] = f1eh->GetParameter(i);
	vm[ieta-1]["ee"][i] = f1ee->GetParameter(i);
      }
      else {
	double h_chi2 = f1h->GetChisquare();
	int h_ndf = max(1,f1h->GetNDF());
	double k = 1.;//1./3.;
	mg["h"][i]->SetPoint(ieta-1, eta, k*h_chi2/h_ndf);
	mg["h"][i]->SetPointError(ieta-1, deta, k*1./sqrt(h_ndf));
	double e_chi2 = f1e->GetChisquare();
	int e_ndf = max(1,f1e->GetNDF());
	mg["e"][i]->SetPoint(ieta-1, eta, k*e_chi2/e_ndf);
	mg["e"][i]->SetPointError(ieta-1, deta, k*1./sqrt(e_ndf));
	double a_chi2 = f1a->GetChisquare();
	int a_ndf = max(1,f1a->GetNDF());
	mg["a"][i]->SetPoint(ieta-1, eta, k*a_chi2/a_ndf);
	mg["a"][i]->SetPointError(ieta-1, deta, k*1./sqrt(a_ndf));

	double ee_chi2 = f1ee->GetChisquare();
	int ee_ndf = max(1,f1ee->GetNDF());
	mg["ee"][i]->SetPoint(ieta-1, eta, k*ee_chi2/ee_ndf);
	mg["ee"][i]->SetPointError(ieta-1, deta, k*1./sqrt(ee_ndf));
	double ea_chi2 = f1ea->GetChisquare();
	int ea_ndf = max(1,f1ea->GetNDF());
	mg["ea"][i]->SetPoint(ieta-1, eta, k*ea_chi2/ea_ndf);
	mg["ea"][i]->SetPointError(ieta-1, deta, k*1./sqrt(ea_ndf));
	double eh_chi2 = f1ee->GetChisquare();
	int eh_ndf = max(1,f1eh->GetNDF());
	mg["eh"][i]->SetPoint(ieta-1, eta, k*eh_chi2/ee_ndf);
	mg["eh"][i]->SetPointError(ieta-1, deta, k*1./sqrt(eh_ndf));

	vm[ieta-1]["eta"][i] = ieta-1;
	vm[ieta-1]["hh"][i] = f1h->GetParameter(i);
	vm[ieta-1]["eh"][i] = f1eh->GetParameter(i);
	vm[ieta-1]["ee"][i] = f1ee->GetParameter(i);
      }

      // Store also response vs |eta| from fits at a few given pTs
      for (int ix = 0; ix != nx; ++ix) {
	int pt = vx[ix];
	if (mg["h_eta"][ix]==0) mg["h_eta"][ix] = new TGraphErrors(0);
	if (mg["e_eta"][ix]==0) mg["e_eta"][ix] = new TGraphErrors(0);
	if (mg["a_eta"][ix]==0) mg["a_eta"][ix] = new TGraphErrors(0);
	if (mg["ea_eta"][ix]==0) mg["ea_eta"][ix] = new TGraphErrors(0);
	if (mg["ee_eta"][ix]==0) mg["ee_eta"][ix] = new TGraphErrors(0);
	if (mg["eh_eta"][ix]==0) mg["eh_eta"][ix] = new TGraphErrors(0);
	double eps = 0.01*pt;
	mg["h_eta"][ix]->SetPoint(ieta-1, eta+dh, f1h->Eval(pt));
	mg["e_eta"][ix]->SetPoint(ieta-1, eta+de, f1e->Eval(pt));
	mg["a_eta"][ix]->SetPoint(ieta-1, eta, f1a->Eval(pt));
	mg["ea_eta"][ix]->SetPoint(ieta-1, eta+dee, f1ea->Eval(pt));
	mg["ee_eta"][ix]->SetPoint(ieta-1, eta+dee, f1ee->Eval(pt));
	mg["eh_eta"][ix]->SetPoint(ieta-1, eta+deh, f1eh->Eval(pt));

	// Calculate error (to be done)
	/*
	mg["h_eta"][ix]->SetPointError(ieta-1, deta, f1h->IntegralError(pt-eps,pt+eps) / (2.*eps));
      	mg["e_eta"][ix]->SetPointError(ieta-1, deta, f1e->IntegralError(pt-eps,pt+eps) / (2.*eps));
    	mg["a_eta"][ix]->SetPointError(ieta-1, deta, f1a->IntegralError(pt-eps,pt+eps) / (2.*eps));
  	mg["ee_eta"][ix]->SetPointError(ieta-1, deta, f1ee->IntegralError(pt-eps,pt+eps) / (2.*eps));
	mg["eh_eta"][ix]->SetPointError(ieta-1, deta, f1eh->IntegralError(pt-eps,pt+eps) / (2.*eps));
	*/
      } // for ix
      
    } // for i

    // Store results also in an individual canvas for better eta bin control
    if (true) { // single canvas plotting

      double eps = 1e-4;
      TH1D *h5i = tdrHist(Form("h5i_%d",ieta),"(rawEcal+rawHcal)/genP",
			  0.+eps,1.3-eps,
			  "p_{T,gen} (GeV)",0.2,1000-eps);
      if (useCMD) h5i->SetYTitle("(ecal+hcal)/true");
      if (useCMD) h5i->SetXTitle("true/cosh(eta) (GeV)");
      TCanvas *c5i = tdrCanvas(Form("c5_%d",ieta),h5i,8,11,kSquare);
      gPad->SetLogx();

      tex->SetTextSize(0.045);
      tex->DrawLatex(0.50,0.87,Form("%1.3f#leq|#eta|<%1.3f",eta-deta,eta+deta));
      
      f1a->Draw("SAME");  f1a_0->Draw("SAME");
      f1e->Draw("SAME");  f1e_0->Draw("SAME");
      f1h->Draw("SAME");  f1h_0->Draw("SAME");
      f1ea->Draw("SAME"); f1ea_0->Draw("SAME");
      f1ee->Draw("SAME"); f1ee_0->Draw("SAME");
      f1eh->Draw("SAME"); f1eh_0->Draw("SAME");

      tdrDraw(hea,"Pz",kOpenCircle,kBlue-9, kSolid,-1,kNone,0, 0.7);
      tdrDraw(pa,"Pz",kOpenCircle,kGray+1, kSolid,-1,kNone,0, 0.7);
      tdrDraw(pe,"Pz",kFullCircle,kBlue-9, kSolid,-1,kNone,0, 0.7);

      tdrDraw(ph,"Pz",kFullCircle,kRed, kSolid,-1,kNone,0, 0.7);
      tdrDraw(hee,"Pz",kOpenDiamond,kMagenta+2, kSolid,-1,kNone,0, 1.0);
      tdrDraw(heh,"Pz",kOpenDiamond,kOrange+2, kSolid,-1,kNone,0, 1.0);

      //legi->Draw();

      c5i->SaveAs(Form("pdf/vsEta/drawPionGun_respHE_eta_%04.0f_%04.0f.pdf",
		       1000*(eta-deta),1000.*(eta+deta)));
      
    } // if single canvas

  } // for ieta

  c5->SaveAs("pdf/drawPionGun_respHE_3D.pdf");
  c5f->SaveAs("pdf/drawPionGun_respFE_3D.pdf");

  
  // Draw fit results vs |eta|
  TH1D *h6 = tdrHist("h6","Parameter",0.,2.0,"|#eta_{gen}|",0,3.139);
  TCanvas *c6 = tdrCanvas("c6",h6,8,11,kSquare);

  tdrDraw(mg["h"][0],"Pz",kFullSquare,kRed,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["h"][1],"Pz",kFullCircle,kRed+1,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["h"][2],"Pz",kFullDiamond,kRed+2,kSolid,-1,kNone,0, 1.0);

  tdrDraw(mg["e"][0],"Pz",kFullSquare,kBlue,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["e"][1],"Pz",kFullCircle,kBlue+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["e"][2],"Pz",kFullDiamond,kBlue+2,kSolid,-1,kNone,0, 0.8);

  tdrDraw(mg["a"][0],"Pz",kFullSquare,kGray+1,kSolid,-1,kNone,0, 0.4);
  tdrDraw(mg["a"][1],"Pz",kFullCircle,kGray+2,kSolid,-1,kNone,0, 0.4);
  tdrDraw(mg["a"][2],"Pz",kFullDiamond,kGray+3,kSolid,-1,kNone,0, 0.7);

  tdrDraw(mg["ee"][0],"Pz",kOpenSquare,kMagenta,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["ee"][1],"Pz",kOpenCircle,kMagenta+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["ee"][2],"Pz",kOpenDiamond,kMagenta+2,kSolid,-1,kNone,0, 0.8);

  tdrDraw(mg["eh"][0],"Pz",kOpenSquare,kOrange,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["eh"][1],"Pz",kOpenCircle,kOrange+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["eh"][2],"Pz",kOpenDiamond,kOrange+2,kSolid,-1,kNone,0, 0.8);

  gPad->RedrawAxis();
  c6->SaveAs("pdf/drawPionGun_parsVsEta.pdf");


  // c range is 0.9,1.3, so 0.4 -> 0.6
  TH1D *h6c = tdrHist("h6c","Parameter 'c'",
		      minc,minc+(maxc-minc)*1.75,"|#eta_{gen}|",0,3.139);
  TCanvas *c6c = tdrCanvas("c6c",h6c,8,11,kSquare);

  tex->SetTextSize(0.045);
  tex->DrawLatex(0.19,0.75,"Ratio to EM scale");
  tex->DrawLatex(0.19,0.70,"c #approx R_{#pi^{+}}(50 GeV)/EM");
  
  tdrDraw(mg["h"][0],"Pz",kFullSquare,kRed,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["e"][0],"Pz",kFullSquare,kBlue,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["a"][0],"Pz",kFullSquare,kGray+2,kSolid,-1,kNone,0, 0.4);
  tdrDraw(mg["ea"][0],"Pz",kOpenSquare,kBlue+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["ee"][0],"Pz",kOpenSquare,kMagenta+2,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["eh"][0],"Pz",kOpenSquare,kOrange+1,kSolid,-1,kNone,0, 0.5);

  TLegend *leg6c = tdrLeg(0.55,0.90-0.045*6,0.80,0.90);
  leg6c->AddEntry(mg["a"][0],"All hadrons","PLE");
  leg6c->AddEntry(mg["h"][0],"H (fE<0.01)","PLE");
  leg6c->AddEntry(mg["e"][0],"E (0.20<fE<0.80)","PLE");
  leg6c->AddEntry(mg["ea"][0],"E (fE#rightarrow0.5)","PLE");
  leg6c->AddEntry(mg["eh"][0],"H of E (fE#rightarrow0)","PLE");
  leg6c->AddEntry(mg["ee"][0],"E of E (fE#rightarrow1)","PLE");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(1.479,0.9,1.479,1.5);
  l->DrawLine(2.500,0.9,2.500,1.3);
  l->SetLineStyle(kDashed);
  l->DrawLine(0.,1.0,3.139,1.0);
  l->SetLineColor(kRed);
  l->DrawLine(0.,1.1,3.139,1.1);
  
  gPad->RedrawAxis();
  c6c->SaveAs("pdf/drawPionGun_par1cVsEta.pdf");

  // m range is 0,0.9, so 0.9 -> 1.35
  TH1D *h6m = tdrHist("h6m","Parameter 'm'",
		      minm,minm+(maxm-minm)*1.75,"|#eta_{gen}|",0,3.139);
  TCanvas *c6m = tdrCanvas("c6m",h6m,8,11,kSquare);

  // https://www.osti.gov/servlets/purl/923647
  tex->DrawLatex(0.19,0.30,"Nuclear interactions");
  tex->DrawLatex(0.19,0.20,"m = 1- #frac{ln(1/(1-F_{#pi^{0}}))}{ln(n)}");
  
  tdrDraw(mg["h"][2],"Pz",kFullDiamond,kRed+2,kSolid,-1,kNone,0, 1.0);
  tdrDraw(mg["e"][2],"Pz",kFullDiamond,kBlue+2,kSolid,-1,kNone,0, 0.8);
  tdrDraw(mg["a"][2],"Pz",kFullDiamond,kGray+3,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["ea"][2],"Pz",kOpenDiamond,kBlue+1,kSolid,-1,kNone,0, 0.8);
  tdrDraw(mg["ee"][2],"Pz",kOpenDiamond,kMagenta+2,kSolid,-1,kNone,0, 0.8);
  tdrDraw(mg["eh"][2],"Pz",kOpenDiamond,kOrange+2,kSolid,-1,kNone,0, 0.8);

  TLegend *leg6m = tdrLeg(0.55,0.90-0.045*6,0.80,0.90);
  leg6m->AddEntry(mg["a"][2],"All hadrons","PLE");
  leg6m->AddEntry(mg["h"][2],"H (fE<0.01)","PLE");
  leg6m->AddEntry(mg["e"][2],"E (0.20<fE<0.80)","PLE");
  leg6m->AddEntry(mg["ea"][2],"E (fE#rightarrow0.5)","PLE");
  leg6m->AddEntry(mg["eh"][2],"H of E (fE#rightarrow0)","PLE");
  leg6m->AddEntry(mg["ee"][2],"E of E (fE#rightarrow1)","PLE");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(1.479,0.0,1.479,1.35);
  l->DrawLine(2.500,0.0,2.500,0.9);
  // Core parameters for nuclear interactions: Fpi0 and n
  double Fpi0 = 1./3.+0.05; // Fraction of pi0 in secondaries ~1./3 + rho0
  double n = 5; // Typical number of nuclear secondaries 5-6
  double m = 1-log(1/(1-Fpi0))/log(n);
  l->SetLineStyle(kDashed);
  l->SetLineColor(kRed);
  l->DrawLine(0,m,3.139,m);
  
  gPad->RedrawAxis();
  c6m->SaveAs("pdf/drawPionGun_par3mVsEta.pdf");
  
  
  // a range is 0.3,1.5, so 1.2 -> 1.8
  TH1D *h6a = tdrHist("h6a","Parameter 'a'",
		      mina,mina+(maxa-mina)*1.75,"|#eta_{gen}|",0,3.139);
  TCanvas *c6a = tdrCanvas("c6a",h6a,8,11,kSquare);

  // https://www.osti.gov/servlets/purl/923647
  tex->DrawLatex(0.19,0.75,"Calorimeter h/e");
  tex->DrawLatex(0.19,0.70,"a = (1-h/e)#timesE_{0}^{1-m}");
  
  tdrDraw(mg["h"][1],"Pz",kFullCircle,kRed+1,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["e"][1],"Pz",kFullCircle,kBlue+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["a"][1],"Pz",kFullCircle,kGray+2,kSolid,-1,kNone,0, 0.4);
  tdrDraw(mg["ea"][1],"Pz",kOpenCircle,kBlue+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["ee"][1],"Pz",kOpenCircle,kMagenta+1,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["eh"][1],"Pz",kOpenCircle,kOrange+1,kSolid,-1,kNone,0, 0.5);

  TLegend *leg6a = tdrLeg(0.55,0.90-0.045*6,0.80,0.90);
  leg6a->AddEntry(mg["a"][1],"All hadrons","PLE");
  leg6a->AddEntry(mg["h"][1],"H (fE<0.01)","PLE");
  leg6a->AddEntry(mg["e"][1],"E (0.20<fE<0.80)","PLE");
  leg6a->AddEntry(mg["ea"][1],"E (fE#rightarrow0.5)","PLE");
  leg6a->AddEntry(mg["eh"][1],"H of E (fE#rightarrow0)","PLE");
  leg6a->AddEntry(mg["ee"][1],"E of E (fE#rightarrow1)","PLE");

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(1.479,0.3,1.479,2.1);
  l->DrawLine(2.500,0.3,2.500,1.5);

  // Calorimeter parameters for h/e and E0
  // https://indico.cern.ch/event/31463/contributions/726204/attachments/603326/830281/08-AndrisSkuja-HCALstatus_LeHCReport.pdf?#page=5
  // (see also pi+/pi- ratio, pi+/p ratio)
  double he_h = 1./2.;//should be: 1./1.4;
  double e0_h = 1.0;
  double a_h = (1-he_h)/pow(e0_h,m-1);
  double he_e = 1/4.;//should be: 1./3.
  double e0_e = 1.0;//2.5;
  double a_e = (1-he_e)/pow(e0_e,m-1);
  l->SetLineStyle(kDashed);
  l->SetLineColor(kRed);
    l->DrawLine(0,a_h,3.139,a_h);
  l->SetLineColor(kBlue);
  l->DrawLine(0,a_e,3.139,a_e);
  
  gPad->RedrawAxis();
  c6a->SaveAs("pdf/drawPionGun_par2aVsEta.pdf");


  // Fit chisquare for quality control
  TH1D *h7 = tdrHist("h7","#chi^{2} / NDF",0,15,"|#eta_{gen}|",0,3.139);
  TCanvas *c7 = tdrCanvas("c7",h7,8,11,kSquare);

  tdrDraw(mg["h"][3],"Pz",kFullStar,kRed,kSolid,-1,kNone,0, 0.7);
  tdrDraw(mg["e"][3],"Pz",kFullStar,kBlue,kSolid,-1,kNone,0, 0.5);
  tdrDraw(mg["a"][3],"Pz",kFullStar,kGray+2,kSolid,-1,kNone,0, 0.4);
  tdrDraw(mg["ea"][3],"Pz",kOpenStar,kBlue+2,kSolid,-1,kNone,0, 1.0);
  tdrDraw(mg["ee"][3],"Pz",kOpenStar,kMagenta+2,kSolid,-1,kNone,0, 1.0);
  tdrDraw(mg["eh"][3],"Pz",kOpenStar,kOrange+1,kSolid,-1,kNone,0, 0.5);

  TLegend *leg7 = tdrLeg(0.35,0.90-0.045*6,0.60,0.90);
  leg7->AddEntry(mg["a"][3],"All hadrons","PLE");
  leg7->AddEntry(mg["h"][3],"H (fE<0.01)","PLE");
  leg7->AddEntry(mg["e"][3],"E (0.20<fE<0.80)","PLE");
  leg7->AddEntry(mg["ea"][3],"E (fE#rightarrow0.5)","PLE");
  leg7->AddEntry(mg["eh"][3],"H of E (fE#rightarrow0)","PLE");
  leg7->AddEntry(mg["ee"][3],"E of E (fE#rightarrow1)","PLE");
  
  gPad->RedrawAxis();
  c7->SaveAs("pdf/drawPionGun_chi2VsEta.pdf");


  // Responvs vs |eta| at a few given pT for stability monitoring
  for (int ix = 0; ix != nx; ++ix) {
    TH1D *h8 = tdrHist(Form("h8_%d",ix),"Response",0,2,"|#eta_{gen}|",0,3.139);
    if (isClosure) h8->SetYTitle("Corrected response");
    if (isClosure) h8->GetYaxis()->SetRangeUser(0.8+1e-8,1.5-1e-5);
    TCanvas *c8 = tdrCanvas(Form("c8_%d",ix),h8,8,11,kSquare);

    tdrDraw(mg["h_eta"][ix],"Pz",kFullCircle,kRed,kSolid,-1,kNone,0, 0.7);
    tdrDraw(mg["e_eta"][ix],"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0, 0.5);
    tdrDraw(mg["a_eta"][ix],"Pz",kFullCircle,kGray+2,kSolid,-1,kNone,0, 0.4);
    tdrDraw(mg["ea_eta"][ix],"Pz",kOpenCircle,kBlue+1,kSolid,-1,kNone,0,0.7);
    tdrDraw(mg["ee_eta"][ix],"Pz",kOpenCircle,kMagenta+2,kSolid,-1,kNone,0,0.7);
    tdrDraw(mg["eh_eta"][ix],"Pz",kOpenCircle,kOrange+1,kSolid,-1,kNone,0,0.6);
    
    TLegend *leg8 = tdrLeg(0.35,0.90-0.045*6,0.60,0.90);
    leg8->AddEntry(mg["a_eta"][ix],"All hadrons","PLE");
    leg8->AddEntry(mg["h_eta"][ix],"H (fE<0.01)","PLE");
    leg8->AddEntry(mg["e_eta"][ix],"E (0.20<fE<0.80)","PLE");
    leg8->AddEntry(mg["ea_eta"][ix],"E (fE#rightarrow0.5)","PLE");
    leg8->AddEntry(mg["eh_eta"][ix],"H of E (fE#rightarrow0)","PLE");
    leg8->AddEntry(mg["ee_eta"][ix],"E of E (fE#rightarrow1)","PLE");
    
    gPad->RedrawAxis();
    c8->SaveAs(Form("pdf/drawPionGun_respVsEta_pt%1.0f.pdf",vx[ix]));
  } // for ix 


  /////////////////////////////////////////////////////////////////////////////
  // Print fit results in text file for implementing in PFEnergyCalibration //
  /////////////////////////////////////////////////////////////////////////////

  ofstream fout("piongun.txt");

  // Count eta bins
  int nb(0);
  for (unsigned int i = 0; i != vm.size() && vm[i]["eta"][2]!=0; ++i) {
    ++nb;
  }

  fout << "// Code auto-generated by pfhc/drawPionGun.C" << endl;
  fout << "// Hadron response |eta| bins" << endl;
  fout << Form("  veta[%d] =",nb+1) << endl;
  for (unsigned int i = 0; i != vm.size() && vm[i]["eta"][2]!=0; ++i) {
    fout << Form("%s%1.3f", i==0 ? "    {" : ", ", vm[i]["eta"][0]);
  }
  fout << Form(", %1.3f};", vm[nb]["eta"][0]) << endl << endl;
  
  fout << "// H-hadron response (c, a, m, R_h,min)" << endl;
  //fout << Form("  phh[%d][4] = {\n",nb);
  fout << Form("array< array<double, 4>, %d> vhh = {{\n",nb);
  for (unsigned int i = 0; i != vm.size() && vm[i]["eta"][2]!=0; ++i) {
    for (int j = 0; j != 4; ++j) {
      fout << Form("%s%6.4f", j==0 ? "    {" : ", ", vm[i]["hh"][j]);
    }
    fout << Form("}, // [%1.3f,%1.3f]\n",vm[i]["eta"][0],vm[i]["eta"][2]);
  }
  fout << "  }};"<<endl<<endl;

  fout << "// EH-hadron H-component response (c, a, m, R_eh,min)" << endl;
  //fout << Form("  peh[%d][4] = {\n",nb);
  fout << Form("array< array<double, 4>, %d> veh = {{\n",nb);
  for (unsigned int i = 0; i != vm.size() && vm[i]["eta"][2]!=0; ++i) {
    for (int j = 0; j != 4; ++j) {
      fout << Form("%s%6.4f", j==0 ? "    {" : ", ", vm[i]["eh"][j]);
    }
    fout << Form("}, // [%1.3f,%1.3f]\n",vm[i]["eta"][0],vm[i]["eta"][2]);
  }
  fout << "  }};"<<endl<<endl;

  fout << "// EH-hadron E-component response (c, a, m, R_ee,min)" << endl;
  //fout << Form("  pee[%d][4] = {\n",nb);
  fout << Form("array< array<double, 4>, %d> vee = {{\n",nb);
  for (unsigned int i = 0; i != vm.size() && vm[i]["eta"][2]!=0; ++i) {
    for (int j = 0; j != 4; ++j) {
      fout << Form("%s%6.4f", j==0 ? "    {" : ", ", vm[i]["ee"][j]);
    }
    fout << Form("}, // [%1.3f,%1.3f]\n",vm[i]["eta"][0],vm[i]["eta"][2]);
  }
  fout << "  }};"<<endl;

  
} // drawPiongun
