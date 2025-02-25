#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"

#include "tdrstyle_mod22.C"

const double dRmax = 0.2;//0.4;

void drawPFAna() {

  setTDRStyle();
  
  //TH1D *h = tdrHist("h","Fraction (%/GeV)",0,3.5,"Energy (GeV)",0,100);
  TH1D *h = tdrHist("h","Fraction (%/GeV)",0,3.5,"Energy (GeV)",-30,160);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.33,0.88,"i#eta=28, E=50 GeV, #pi^{-}, H-hadron: f_{E}^{rec}<0.05");
  //tex->DrawLatex(0.33,0.88,"i#eta=28, E=50 GeV, #pi^{-}, H-hadron: f_{E}^{sim}<0.05");
  //tex->DrawLatex(0.53,0.84,"E-hadron: E_{HCAL}^{sim} < 2.5 GeV");
  tex->DrawLatex(0.50,0.84,"EH-hadron: E_{HCAL}^{sim} > 2.5 GeV");
  //tex->DrawLatex(0.64,0.80,"#DeltaR(Rec,Sim)<0.4");
  tex->DrawLatex(0.50,0.80,Form("#DeltaR(Rec,Sim)<%1.2g minus RC",dRmax));
  
  TLegend *leg = tdrLeg(0.63,0.78-0.035*12,0.88,0.78);
  leg->SetTextSize(0.035);
  
  TFile *f = new TFile("PFAna.root","READ");

  TH1D *hgensum = (TH1D*)f->Get("hgensum"); assert(hgensum);
  TH1D *hgensum_hh = (TH1D*)f->Get("hgensum_hh"); assert(hgensum_hh);
  TH1D *hgensum_eh = (TH1D*)f->Get("hgensum_eh"); assert(hgensum_eh);
  TH1D *hgensum_ee = (TH1D*)f->Get("hgensum_ee"); assert(hgensum_ee);

  /*
  TH1D *hrecsum = (TH1D*)f->Get("hrecsum"); assert(hrecsum);
  TH1D *hrecsum_hh = (TH1D*)f->Get("hrecsum_hh"); assert(hrecsum_hh);
  TH1D *hrecsum_eh = (TH1D*)f->Get("hrecsum_eh"); assert(hrecsum_eh);
  TH1D *hrecsum_ee = (TH1D*)f->Get("hrecsum_ee"); assert(hrecsum_ee);
  
  TH1D *hpfsum = (TH1D*)f->Get("hpfsum"); assert(hpfsum);
  TH1D *hpfsum_hh = (TH1D*)f->Get("hpfsum_hh"); assert(hpfsum_hh);
  TH1D *hpfsum_eh = (TH1D*)f->Get("hpfsum_eh"); assert(hpfsum_eh);
  TH1D *hpfsum_ee = (TH1D*)f->Get("hpfsum_ee"); assert(hpfsum_ee);
  */

  // Variants for withPU-RC
  TH1D *hrecsum = (TH1D*)f->Get("hrecrcsum"); assert(hrecsum);
  TH1D *hrecsum_hh = (TH1D*)f->Get("hrecrcsum_hh"); assert(hrecsum_hh);
  TH1D *hrecsum_eh = (TH1D*)f->Get("hrecrcsum_eh"); assert(hrecsum_eh);
  TH1D *hrecsum_ee = (TH1D*)f->Get("hrecrcsum_ee"); assert(hrecsum_ee);

  TH1D *hpfsum = (TH1D*)f->Get("hpfrcsum"); assert(hpfsum);
  TH1D *hpfsum_hh = (TH1D*)f->Get("hpfrcsum_hh"); assert(hpfsum_hh);
  TH1D *hpfsum_eh = (TH1D*)f->Get("hpfrcsum_eh"); assert(hpfsum_eh);
  TH1D *hpfsum_ee = (TH1D*)f->Get("hpfrcsum_ee"); assert(hpfsum_ee);

  TH1D *hrcsum = (TH1D*)f->Get("hrcsum"); assert(hrcsum);
  TH1D *hrcsum_hh = (TH1D*)f->Get("hrcsum_hh"); assert(hrcsum_hh);
  TH1D *hrcsum_eh = (TH1D*)f->Get("hrcsum_eh"); assert(hrcsum_eh);
  TH1D *hrcsum_ee = (TH1D*)f->Get("hrcsum_ee"); assert(hrcsum_ee);

  // Make backup copies before modifying
  // so we don't later accidentally get modified copy from open file
  hgensum = (TH1D*)hgensum->Clone("hgensum_cp");
  hgensum_hh = (TH1D*)hgensum_hh->Clone("hgensum_hh_cp");
  hgensum_eh = (TH1D*)hgensum_eh->Clone("hgensum_eh_cp");
  hgensum_ee = (TH1D*)hgensum_ee->Clone("hgensum_ee_cp");

  hrecsum = (TH1D*)hrecsum->Clone("hrecsum_cp");
  hrecsum_hh = (TH1D*)hrecsum_hh->Clone("hrecsum_hh_cp");
  hrecsum_eh = (TH1D*)hrecsum_eh->Clone("hrecsum_eh_cp");
  hrecsum_ee = (TH1D*)hrecsum_ee->Clone("hrecsum_ee_cp");

  hpfsum = (TH1D*)hpfsum->Clone("hpfsum_cp");
  hpfsum_hh = (TH1D*)hpfsum_hh->Clone("hpfsum_hh_cp");
  hpfsum_eh = (TH1D*)hpfsum_eh->Clone("hpfsum_eh_cp");
  hpfsum_ee = (TH1D*)hpfsum_ee->Clone("hpfsum_ee_cp");

  hrcsum = (TH1D*)hrcsum->Clone("hrcsum_cp");
  hrcsum_hh = (TH1D*)hrcsum_hh->Clone("hrcsum_hh_cp");
  hrcsum_eh = (TH1D*)hrcsum_eh->Clone("hrcsum_eh_cp");
  hrcsum_ee = (TH1D*)hrcsum_ee->Clone("hrcsum_ee_cp");
  
  
  hgensum_hh->Scale(100./hgensum->Integral());
  hgensum_eh->Scale(100./hgensum->Integral());
  hgensum_ee->Scale(100./hgensum->Integral());
  hgensum->Scale(100./hgensum->Integral());

  hrecsum_hh->Scale(100./hrecsum->Integral());
  hrecsum_eh->Scale(100./hrecsum->Integral());
  hrecsum_ee->Scale(100./hrecsum->Integral());
  hrecsum->Scale(100./hrecsum->Integral());

  hpfsum_hh->Scale(100./hpfsum->Integral());
  hpfsum_eh->Scale(100./hpfsum->Integral());
  hpfsum_ee->Scale(100./hpfsum->Integral());
  hpfsum->Scale(100./hpfsum->Integral());

  hrcsum_hh->Scale(100./hrcsum->Integral());
  hrcsum_eh->Scale(100./hrcsum->Integral());
  hrcsum_ee->Scale(100./hrcsum->Integral());
  hrcsum->Scale(100./hrcsum->Integral());
  
  tdrDraw(hgensum,"HIST",kNone,kBlack,kSolid,-1,1001,kGray);
  tdrDraw(hgensum_ee,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray+1);
  tdrDraw(hgensum_eh,"HIST",kNone,kBlue,kSolid,-1,1001,kRed-1);
  tdrDraw(hgensum_hh,"HIST",kNone,kRed,kSolid,-1,1001,kBlue-9);

  hgensum->SetFillColorAlpha(kGray,0.20);
  hgensum_ee->SetFillColorAlpha(kGray+1,0.50);
  hgensum_eh->SetFillColorAlpha(kBlue-9,0.50);
  hgensum_hh->SetFillColorAlpha(kRed-9,0.50);

  tdrDraw(hrecsum,"HIST",kNone,kBlack,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_ee,"HIST",kNone,kGray+3,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_eh,"HIST",kNone,kBlue+1,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_hh,"HIST",kNone,kRed+1,kSolid,-1,kNone,0,0.7,2);

  tdrDraw(hpfsum,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_ee,"Pz",kFullCircle,kGray+2,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_eh,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_hh,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7,1);

  leg->AddEntry(hgensum,"SimHits","FL");
  leg->AddEntry(hgensum_hh,"SimHits, H","FL");
  leg->AddEntry(hgensum_eh,"SimHits, EH","FL");
  leg->AddEntry(hgensum_ee,"SimHits, E","FL");
  leg->AddEntry(hrecsum,"RecHits","FL");
  leg->AddEntry(hrecsum_hh,"RecHits, H","FL");
  leg->AddEntry(hrecsum_eh,"RecHits, EH","FL");
  leg->AddEntry(hrecsum_ee,"RecHits, E","FL");
  leg->AddEntry(hpfsum,"PFClusters","PLE");
  leg->AddEntry(hpfsum_hh,"PFClusters, H","PLE");
  leg->AddEntry(hpfsum_eh,"PFClusters, EH","PLE");
  leg->AddEntry(hpfsum_ee,"PFClusters, E","PLE");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawPFAna.pdf");
  

  TH1D *h_2 = tdrHist("h_2","Energy per depth (GeV)",0,22,"Depth",0,8);
  //TH1D *h_2 = tdrHist("h_2","Energy per depth (GeV)",0,40,"Depth",0,8);
  TH1D *hd_2 = tdrHist("hd_2","Rec / Sim",0.65,1.2,"Depth",0,8);
  //TH1D *hd_2 = tdrHist("hd_2","Rec / Sim",0.,12.0,"Depth",0,8);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c2 = tdrDiCanvas("c2",h_2,hd_2,8,11);

  TProfile *pg = (TProfile*)f->Get("pgendepth"); assert(pg);
  TProfile *pg_hh = (TProfile*)f->Get("pgendepth_hh"); assert(pg_hh);
  TProfile *pg_eh = (TProfile*)f->Get("pgendepth_eh"); assert(pg_eh);
  /*
  TProfile *pr = (TProfile*)f->Get("precdepth"); assert(pr);
  TProfile *pr_hh = (TProfile*)f->Get("precdepth_hh"); assert(pr_hh);
  TProfile *pr_eh = (TProfile*)f->Get("precdepth_eh"); assert(pr_eh);
  */
  // Variant designed for withPU-RC
  TProfile *pr = (TProfile*)f->Get("precrcdepth"); assert(pr);
  TProfile *pr_hh = (TProfile*)f->Get("precrcdepth_hh"); assert(pr_hh);
  TProfile *pr_eh = (TProfile*)f->Get("precrcdepth_eh"); assert(pr_eh);
  
  c2->cd(1);

  TLegend *leg2 = tdrLeg(0.63,0.70-0.05*6,0.88,0.70);
  leg2->SetTextSize(0.045);

  tex->SetTextSize(0.044);
  tex->DrawLatex(0.33,0.85,"i#eta=28, E=50 GeV, #pi^{-}, H-had.: f_{E}^{sim}<0.05");
  tex->DrawLatex(0.53,0.80,"EH-had.: E_{HCAL}^{sim} > 2.5 GeV");
  //tex->DrawLatex(0.64,0.75,"#DeltaR(Rec,Sim)<0.4");
  tex->DrawLatex(0.50,0.75,Form("#DeltaR(Rec,Sim)<%1.2g minus RC",dRmax));
  
  tdrDraw(pg,"HIST",kNone,kBlack,kSolid,-1,1001,kGray);
  //tdrDraw(pg_ee,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray+1);
  tdrDraw(pg_eh,"HIST",kNone,kBlue,kSolid,-1,1001,kRed-1);
  tdrDraw(pg_hh,"HIST",kNone,kRed,kSolid,-1,1001,kBlue-9);

  pg->SetFillColorAlpha(kGray,0.50);
  //pg_ee->SetFillColorAlpha(kGray+1,0.50);
  pg_eh->SetFillColorAlpha(kBlue-9,0.50);
  pg_hh->SetFillColorAlpha(kRed-9,0.50);

  tdrDraw(pr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7,1);
  //tdrDraw(pr_ee,"Pz",kFullCircle,kGray+2,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(pr_eh,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(pr_hh,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7,1);

  leg2->AddEntry(pg,"SimHits","FL");
  leg2->AddEntry(pg_hh,"SimHits, H","FL");
  leg2->AddEntry(pg_eh,"SimHits, EH","FL");

  leg2->AddEntry(pr,"RecHits","PLE");
  leg2->AddEntry(pr_hh,"RecHits, H","PLE");
  leg2->AddEntry(pr_eh,"RecHits, EH","PLE");
  
  gPad->RedrawAxis();

  c2->cd(2);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0,1,8,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,1.1,8,1.1);
  l->DrawLine(0,0.8,8,0.8);

  TH1D *hr = pr->ProjectionX("hr");
  TH1D *hr_hh = pr_hh->ProjectionX("hr_hh");
  TH1D *hr_eh = pr_eh->ProjectionX("hr_eh");

  hr->Divide(pg);
  hr_hh->Divide(pg_hh);
  hr_eh->Divide(pg_eh);

  tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hr_eh,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hr_hh,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7,1);

  gPad->RedrawAxis();

  c2->SaveAs("pdf/drawPFAna_depths.pdf");


//TH1D *h_3 = tdrHist("hd_3","Fractional energy loss (%)",0.,10.,
//TH1D *h_3 = tdrHist("hd_3","Fractional energy loss (%)",-5,5.,
  TH1D *h_3 = tdrHist("hd_3","Fractional energy loss (%)",-5,10.,
		      //"HCalPFcut (GeV)",0,2);
		      "HCalPFcut (GeV)",0,5);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);

  /*
  TH2D *h2r = (TH2D*)f->Get("h2recdepth"); assert(h2r);
  TH2D *h2r_hh = (TH2D*)f->Get("h2recdepth_hh"); assert(h2r_hh);
  TH2D *h2r_eh = (TH2D*)f->Get("h2recdepth_eh"); assert(h2r_eh);
  */
  TH2D *h2r = (TH2D*)f->Get("h2recrcdepth"); assert(h2r);
  TH2D *h2r_hh = (TH2D*)f->Get("h2recrcdepth_hh"); assert(h2r_hh);
  TH2D *h2r_eh = (TH2D*)f->Get("h2recrcdepth_eh"); assert(h2r_eh);

  TH1D *ht1 = h2r->ProjectionY("ht1",2,2);
  TH1D *ht1_hh = h2r_hh->ProjectionY("ht1_hh",2,2);
  TH1D *ht1_eh = h2r_hh->ProjectionY("ht1_eh",2,2);
  for (int i = 1; i != ht1->GetNbinsX()+1; ++i) {
    ht1->SetBinContent(i, ht1->GetBinContent(i-1) + ht1->GetBinContent(i));
    ht1_hh->SetBinContent(i, ht1_hh->GetBinContent(i-1) +
			  ht1_hh->GetBinContent(i));
    ht1_eh->SetBinContent(i, ht1->GetBinContent(i-1) +
			  ht1->GetBinContent(i));
  } // for i
  ht1->Scale(100.);
  ht1_hh->Scale(100.);
  ht1_eh->Scale(100.);

  TH1D *ht2 = h2r->ProjectionY("ht2",3,3);
  TH1D *ht2_hh = h2r_hh->ProjectionY("ht2_hh",3,3);
  TH1D *ht2_eh = h2r_hh->ProjectionY("ht2_eh",3,3);
  for (int i = 1; i != ht2->GetNbinsX()+1; ++i) {
    ht2->SetBinContent(i, ht2->GetBinContent(i-1) + ht2->GetBinContent(i));
    ht2_hh->SetBinContent(i, ht2_hh->GetBinContent(i-1) +
			  ht2_hh->GetBinContent(i));
    ht2_eh->SetBinContent(i, ht2->GetBinContent(i-1) +
			  ht2->GetBinContent(i));
  } // for i
  ht2->Scale(100.);
  ht2_hh->Scale(100.);
  ht2_eh->Scale(100.);

  tdrDraw(ht1,"Pz",kFullSquare,kBlack,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(ht1_eh,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(ht1_hh,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7,1);
  
  tdrDraw(ht2,"Pz",kOpenSquare,kBlack,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(ht2_eh,"Pz",kOpenCircle,kBlue,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(ht2_hh,"Pz",kOpenCircle,kRed,kSolid,-1,kNone,0,0.7,1);

  double emax = 4.0;//3.0;//1.2;
  double emax_hh = 5.0;//2.0;

  ht1->GetXaxis()->SetRangeUser(0,emax);
  ht1_eh->GetXaxis()->SetRangeUser(0,emax);
  ht1_hh->GetXaxis()->SetRangeUser(0,emax_hh);
  
  ht2->GetXaxis()->SetRangeUser(0,emax);
  ht2_eh->GetXaxis()->SetRangeUser(0,emax);
  ht2_hh->GetXaxis()->SetRangeUser(0,emax_hh);

  
  TLegend *leg3 = tdrLeg(0.65,0.85-6*0.05,0.90,0.85);
  leg3->AddEntry(ht1_eh,"Depth1, EH","PLE");
  leg3->AddEntry(ht2_eh,"Depth2, EH","PLE");
  leg3->AddEntry(ht1,"Depth1, all","PLE");
  leg3->AddEntry(ht2,"Depth2, all","PLE");
  leg3->AddEntry(ht1_hh,"Depth1, H","PLE");
  leg3->AddEntry(ht2_hh,"Depth2, H","PLE");

  l->DrawLine(0.8,0,0.8,10);
  l->DrawLine(0,0,h_3->GetXaxis()->GetXmax(),0);
  
  
  c3->SaveAs("pdf/drawPFAna_HcalPFcut.pdf");

  
  TH1D *h_4 = tdrHist("hd_4","Fraction (% / GeV)",0,5.5,
		      "Energy (GeV)",-30,160);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c4 = tdrCanvas("c4",h_4,8,11,kSquare);

  TH1D *hgensum_hh_raw = (TH1D*)f->Get("hgensum_hh"); assert(hgensum_hh_raw);
  TH1D *hrecsum_hh_raw = (TH1D*)f->Get("hrecsum_hh"); assert(hrecsum_hh_raw);
  TH1D *hpfsum_hh_raw = (TH1D*)f->Get("hpfsum_hh"); assert(hpfsum_hh_raw);
  TH1D *hrcsum_hh_raw = (TH1D*)f->Get("hrcsum_hh"); assert(hrcsum_hh_raw);

  TH1D *hrecsum_hh_sub = (TH1D*)f->Get("hrecrcsum_hh"); assert(hrecsum_hh_sub);
  TH1D *hpfsum_hh_sub = (TH1D*)f->Get("hpfrcsum_hh"); assert(hpfsum_hh_sub);
  
  hrecsum_hh_raw->Scale(100./hgensum_hh_raw->Integral());
  hpfsum_hh_raw->Scale(100./hgensum_hh_raw->Integral());
  hrcsum_hh_raw->Scale(100./hgensum_hh_raw->Integral());

  hrecsum_hh_sub->Scale(100./hgensum_hh_raw->Integral());
  hpfsum_hh_sub->Scale(100./hgensum_hh_raw->Integral());
  
  hgensum_hh_raw->Scale(100./hgensum_hh_raw->Integral());
  
  tdrDraw(hgensum_hh_raw,"HIST",kNone,kRed,kSolid,-1,1001,kBlue-9);
  hgensum_hh_raw->SetFillColorAlpha(kRed-9,0.50);
  tdrDraw(hrecsum_hh_raw,"HIST",kNone,kGreen+1,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hrecsum_hh_sub,"HIST",kNone,kGreen+2,kSolid,-1,kNone,0,0.7,2);
  tdrDraw(hpfsum_hh_raw,"HPz",kOpenCircle,kRed,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hpfsum_hh_sub,"HPz",kFullCircle,kRed+1,kSolid,-1,kNone,0,0.7,1);
  tdrDraw(hrcsum_hh_raw,"HIST",kNone,kBlue,kSolid,-1,kNone,0,0.7,1);


  TLegend *leg4 = tdrLeg(0.55,0.85-6*0.05,0.80,0.85);
  leg4->AddEntry(hgensum_hh_raw,"SimHits","FL");
  leg4->AddEntry(hrecsum_hh_raw,"RecHits","FL");
  leg4->AddEntry(hpfsum_hh_raw,"PFClusters","PLE");
  leg4->AddEntry(hrcsum_hh_raw,"RandomCone","FL");
  
  gPad->RedrawAxis();

  c4->SaveAs("pdf/drawPFAna_RecoLevelsForH.pdf");


  double eps = 1e-4;
  TH1D *h_5 = tdrHist("hd_5","#Deltai#phi",-7.5,6.5,//-6.5,6.5,
		      "#Deltai#eta",-3.5-eps,1.5+eps);//-6.5,4.5);
  extraText = "Private";
  lumi_136TeV = "Fikri's private production";
  TCanvas *c5 = tdrCanvas("c5",h_5,8,11,kSquare);
  gPad->SetRightMargin(0.15);

  // Is already H hadrons here
  TH2D *h2etaphi = (TH2D*)f->Get("h2etaphigen"); assert(h2etaphi);
  //TH2D *h2etaphi = (TH2D*)f->Get("h2etaphirec"); assert(h2etaphi);

  h2etaphi->Scale(100.);
  h2etaphi->RebinY(2);
  h2etaphi->GetZaxis()->SetTitle("Fraction (%)");
  //h2etaphi->GetZaxis()->SetRangeUser(0,68);//0.68);
  h2etaphi->GetXaxis()->SetRangeUser(-3.5,1.5);//-6.5,1.5);
  h2etaphi->GetZaxis()->SetRangeUser(0,11);
  //h2etaphi->Draw("COLZ SAME");
  gStyle->SetPaintTextFormat("1.2g");
  h2etaphi->Draw("COLZ TEXT SAME");

  gPad->RedrawAxis();
  gPad->Update();
  c5->SaveAs("pdf/drawPFAna_h2etaphi.pdf");
}
