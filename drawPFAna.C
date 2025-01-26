#include "tdrstyle_mod22.C"

void drawPFAna() {

  setTDRStyle();
  
  TH1D *h = tdrHist("h","Fraction (%)",0,3.5,"Energy (GeV)",0,100);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.33,0.88,"i#eta=28, E=50 GeV, #pi^{-}, H-hadron: f_{E}<0.05");
  
  TLegend *leg = tdrLeg(0.62,0.85-0.035*9,0.87,0.85);
  leg->SetTextSize(0.035);
  
  TFile *f = new TFile("PFAna.root","READ");

  TH1D *hgensum = (TH1D*)f->Get("hgensum"); assert(hgensum);
  TH1D *hgensum_h = (TH1D*)f->Get("hgensum_h"); assert(hgensum_h);
  TH1D *hgensum_e = (TH1D*)f->Get("hgensum_e"); assert(hgensum_e);
  
  TH1D *hrecsum = (TH1D*)f->Get("hrecsum"); assert(hrecsum);
  TH1D *hrecsum_h = (TH1D*)f->Get("hrecsum_h"); assert(hrecsum_h);
  TH1D *hrecsum_e = (TH1D*)f->Get("hrecsum_e"); assert(hrecsum_e);

  TH1D *hpfsum = (TH1D*)f->Get("hpfsum"); assert(hpfsum);
  TH1D *hpfsum_h = (TH1D*)f->Get("hpfsum_h"); assert(hpfsum_h);
  TH1D *hpfsum_e = (TH1D*)f->Get("hpfsum_e"); assert(hpfsum_e);

  hgensum_h->Scale(100./hgensum->Integral());
  hgensum_e->Scale(100./hgensum->Integral());
  hgensum->Scale(100./hgensum->Integral());

  hrecsum_h->Scale(100./hrecsum->Integral());
  hrecsum_e->Scale(100./hrecsum->Integral());
  hrecsum->Scale(100./hrecsum->Integral());

  hpfsum_h->Scale(100./hpfsum->Integral());
  hpfsum_e->Scale(100./hpfsum->Integral());
  hpfsum->Scale(100./hpfsum->Integral());
  
  tdrDraw(hgensum,"HIST",kNone,kBlack,kSolid,-1,1001,kGray);
  tdrDraw(hgensum_e,"HIST",kNone,kBlue,kSolid,-1,1001,kRed-1);
  tdrDraw(hgensum_h,"HIST",kNone,kRed,kSolid,-1,1001,kBlue-9);

  hgensum->SetFillColorAlpha(kGray,0.20);
  hgensum_e->SetFillColorAlpha(kBlue-9,0.50);
  hgensum_h->SetFillColorAlpha(kRed-9,0.50);

  tdrDraw(hrecsum,"HIST",kNone,kBlack,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_e,"HIST",kNone,kBlue+1,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_h,"HIST",kNone,kRed+1,kSolid,-1,kNone,0,0.7,2);

  tdrDraw(hpfsum,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_e,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_h,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7,1);

  leg->AddEntry(hgensum,"SimHits","FL");
  leg->AddEntry(hgensum_h,"SimHits, H","FL");
  leg->AddEntry(hgensum_e,"SimHits, EH","FL");
  leg->AddEntry(hrecsum,"RecHits","FL");
  leg->AddEntry(hrecsum_h,"RecHits, H","FL");
  leg->AddEntry(hrecsum_e,"RecHits, EH","FL");
  leg->AddEntry(hpfsum,"PFClusters","PLE");
  leg->AddEntry(hpfsum_h,"PFClusters, H","PLE");
  leg->AddEntry(hpfsum_e,"PFClusters, EH","PLE");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawPFAna.pdf");
}
