// Purpose: Draw tracking, calorimeter, and selection efficiency in 2D
//          Also draw H-hadron fraction before and after selection
//          Expectation is that PFHC makes sense when efficiencies are high
//          and H-hadron fraction is close to expected 30-40%
//          Otherwise we're mixing up H and HE hadrons, keeping only
//          ECAL or HCAL component, or cutting out calo tails
#include "TFile.h"
#include "TH2D.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"

void normalizeByP(TH2D* h2);
TH2D* getResolution(TProfile2D *p2, string name);

void drawPiongunEfficiency() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/drawPiongunEfficiency");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/drawPiongunEfficiency");

  TFile *f = new TFile("piongun.root","READ");
  assert(f && !f->IsZombie());

  TH2D *h2p = (TH2D*)f->Get("h2np"); assert(h2p);
  TH2D *h2pn = (TH2D*)h2p->Clone("h2pn");
  normalizeByP(h2pn);
  h2p->Scale(1./10000.);

  TProfile2D *p2t = (TProfile2D*)f->Get("p2e_trk"); assert(p2t);
  TProfile2D *p2c = (TProfile2D*)f->Get("p2e_cal"); assert(p2c);
  TProfile2D *p2s = (TProfile2D*)f->Get("p2e_sel"); assert(p2s);
  TProfile2D *p2h = (TProfile2D*)f->Get("p2e_had"); assert(p2h);
  TProfile2D *p2hh = (TProfile2D*)f->Get("p2e_has"); assert(p2hh);

  TProfile2D *p2ra = (TProfile2D*)f->Get("p2r_all"); assert(p2ra);
  TProfile2D *p2rh = (TProfile2D*)f->Get("p2r_hhh"); assert(p2rh);
  TProfile2D *p2re = (TProfile2D*)f->Get("p2r_ehh"); assert(p2re);

  TH2D *h2sa = getResolution(p2ra,"h2sa");
  TH2D *h2sh = getResolution(p2rh,"h2sh");
  TH2D *h2se = getResolution(p2re,"h2se");
  
  
  TH1D *h1p = tdrHist("h1p","#eta_{gen}",0,3.0,"p_{gen}",0.1,5000.);
  lumi_136TeV = "Simulation";
  extraText = "Private";
  TCanvas *c1p = tdrCanvas("c1p",h1p,8,11,kSquare);
  h1p->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  h2p->Draw("COLZ SAME");
  h2p->SetZTitle("Events per bin (per 10k)");
  h2p->GetZaxis()->SetRangeUser(0,2);
  h2p->GetZaxis()->SetTitleOffset(1.5);
  
  gPad->RedrawAxis();
  gPad->Update();
  c1p->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c1p.pdf");

  
  TH1D *h1pn = tdrHist("h1pn","#eta_{gen}",0,3.0,"p_{gen}",0.1,5000.);
  TCanvas *c1pn = tdrCanvas("c1pn",h1pn,8,11,kSquare);
  h1pn->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  h2pn->Draw("COLZ SAME");
  h2pn->SetZTitle("Events per GeV (relative to sample)");
  h2pn->GetZaxis()->SetRangeUser(0,2);
  h2pn->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c1pn->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c1pn.pdf");

  
  TH1D *h2t = tdrHist("h2t","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c2t = tdrCanvas("c2t",h2t,8,11,kSquare);
  h2t->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2t->Draw("COLZ SAME");
  p2t->SetZTitle("Tracking efficiency");
  p2t->GetZaxis()->SetRangeUser(0,1);
  p2t->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c2t->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c2t.pdf");

  
  TH1D *h2c = tdrHist("h2c","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c2c = tdrCanvas("c2c",h2c,8,11,kSquare);
  h2c->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2c->Draw("COLZ SAME");
  p2c->SetZTitle("Calorimeter efficiency");
  p2c->GetZaxis()->SetRangeUser(0,1);
  p2c->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c2c->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c2c.pdf");

  
  TH1D *h2s = tdrHist("h2s","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c2s = tdrCanvas("c2s",h2s,8,11,kSquare);
  h2s->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2s->Draw("COLZ SAME");
  p2s->SetZTitle("Selection efficiency");
  p2s->GetZaxis()->SetRangeUser(0,1);
  p2s->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c2s->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c2s.pdf");


  TH1D *h3h = tdrHist("h3h","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c3h = tdrCanvas("c3h",h3h,8,11,kSquare);
  h3h->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2h->Draw("COLZ SAME");
  p2h->SetZTitle("H-hadron fraction (before selection)");
  p2h->GetZaxis()->SetRangeUser(0,1);
  p2h->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c3h->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c3h.pdf");

  
  TH1D *h3hh = tdrHist("h3hh","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c3hh = tdrCanvas("c3hh",h3hh,8,11,kSquare);
  h3hh->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2hh->Draw("COLZ SAME");
  p2hh->SetZTitle("H-hadron fraction (after selection)");
  p2hh->GetZaxis()->SetRangeUser(0,1);
  p2hh->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c3hh->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c3hh.pdf");


  TH1D *h4a = tdrHist("h4a","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c4a = tdrCanvas("c4a",h4a,8,11,kSquare);
  h4a->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2ra->Draw("COLZ SAME");
  p2ra->SetZTitle("Hadron response");
  p2ra->GetZaxis()->SetRangeUser(0,1.3);
  p2ra->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c4a->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c4a.pdf");

  
  TH1D *h4h = tdrHist("h4h","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c4h = tdrCanvas("c4h",h4h,8,11,kSquare);
  h4h->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2rh->Draw("COLZ SAME");
  p2rh->SetZTitle("H-hadron response");
  p2rh->GetZaxis()->SetRangeUser(0,1.3);
  p2rh->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c4h->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c4h.pdf");

  
  TH1D *h4e = tdrHist("h4e","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c4e = tdrCanvas("c4e",h4e,8,11,kSquare);
  h4e->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  p2re->Draw("COLZ SAME");
  p2re->SetZTitle("EH-hadron response");
  p2re->GetZaxis()->SetRangeUser(0,1.3);
  p2re->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c4e->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c4e.pdf");


  TH1D *h5a = tdrHist("h5a","|#eta_{gen}|",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c5a = tdrCanvas("c5a",h5a,8,11,kSquare);
  h5a->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  h2sa->Draw("COLZ SAME");
  h2sa->SetZTitle("Hadron resolution");
  h2sa->GetZaxis()->SetRangeUser(0,1.0);
  h2sa->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c5a->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c5a.pdf");

  
  TH1D *h5h = tdrHist("h5h","|#eta_{gen}|",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c5h = tdrCanvas("c5h",h5h,8,11,kSquare);
  h5h->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  h2sh->Draw("COLZ SAME");
  h2sh->SetZTitle("H-hadron resolution");
  h2sh->GetZaxis()->SetRangeUser(0,1.0);
  h2sh->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c5h->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c5h.pdf");

  
  TH1D *h5e = tdrHist("h5e","#eta_{gen}",0,3.0,"p_{T,gen}",0.1,5000.);
  TCanvas *c5e = tdrCanvas("c5e",h5e,8,11,kSquare);
  h5e->GetXaxis()->SetMoreLogLabels(kFALSE);
  gPad->SetRightMargin(0.15);
  gPad->SetLogx();
  h2se->Draw("COLZ SAME");
  h2se->SetZTitle("EH-hadron resolution");
  h2se->GetZaxis()->SetRangeUser(0,1.0);
  h2se->GetZaxis()->SetTitleOffset(1.5);

  gPad->RedrawAxis();
  gPad->Update();
  c5e->SaveAs("pdf/drawPiongunEfficiency/drawPiongunEfficiency_c5e.pdf");
  
  
} // drawPionGunEfficiency


void normalizeByP(TH2D *h2) {
  // Events [0.2,10], [0.2,200], [200,500], [500, 5000]
  //double vn[] = {11.8e6, 5.5e6, 5.4e6, 12.3e6};
  //double vn[] = {11.8e6*2.5, 5.5e6, 5.4e6*1.05, 12.3e6*1.05*1.10};
  double vn[] = {19997000., 7955000., 7919000., 20000000.};
  vn[0] = vn[0]+vn[1]*9.8/199.8;
  vn[1] = vn[1]*(1-9.8/199.8);
  double vw[] = {9.8, 180., 300., 4500.};
  double k = 2./1e4;
  //double nb = (h2->GetNbinsX() * h2->GetNbinsY() * 3.0/3.2);
  for (int i = 1 ; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1 ; j != h2->GetNbinsY()+1; ++j) {
      double p = h2->GetXaxis()->GetBinLowEdge(i);
      double dp = h2->GetXaxis()->GetBinWidth(i);
      double w = (p<10 ? vw[0] : p<200 ? vw[1] : p<500 ? vw[2] : vw[3]);
      double n = (p<10 ? vn[0] : p<200 ? vn[1] : p<500 ? vn[2] : vn[3]);
      h2->SetBinContent(i, j, h2->GetBinContent(i, j) / dp * (w/9.8) * (vn[0]/n) * k);
    } // for j
  } // for i
} // for normalizeByP


TH2D *getResolution(TProfile2D *p2, string name) {

  Option_t *opt = p2->GetErrorOption();
  p2->SetErrorOption("S");
  TH2D *h2 = p2->ProjectionXY(name.c_str(),"e");
  p2->SetErrorOption(opt);
  
  for (int i = 1 ; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1 ; j != h2->GetNbinsY()+1; ++j) {
      double s = h2->GetBinError(i,j);
      double r = p2->GetBinContent(i,j);
      double e = p2->GetBinError(i,j);
      if (s>0 && r>0) {
	h2->SetBinContent(i, j, s / r);
	h2->SetBinError(i, j, e / s);
      }
    } // for j
  } // for i

  return h2;
} // getResolution
