// Draw resolution from corrected piongun results

// Include necessary headers
#include <iostream>
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"

// Include your TDR style macros (adjust the path if necessary)
#include "tdrstyle_mod22.C"

void drawPiongunResolution(string file = "", string file2 = "") {
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/vsEta");

  
  // Open the ROOT file containing p3rf
  TFile *f = new TFile(file=="" ? "piongun.root" : file.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file " << file << "!" << std::endl;
    return;
  }

  TFile *f2 = new TFile(file2=="" ? "piongun_pfhc_v2.root" : file2.c_str(),
			"READ");
  if (!f2 || f2->IsZombie()) {
    std::cerr << "Error opening file " << file2 << "!" << std::endl;
    return;
  }
  
  // Get the TProfile3D object (abseta, genpt, fe, resp)
  TProfile3D *p3rf = (TProfile3D*)f->Get("p3rf");
  if (!p3rf) {
    std::cerr << "Error: p3rf not found in file " << file << "!" << std::endl;
    return;
  }

  TProfile3D *p3rf2 = (TProfile3D*)f2->Get("p3rf");
  if (!p3rf2) {
    std::cerr << "Error: p3rf not found in file " << file2 << "!" << std::endl;
    return;
  }

  
  // Step 1: 2D scan of resolution vs eta and pT (EH, H and all)
  //////////////////////////////////////////////////////////////
  TProfile2D *p2a = p3rf->Project3DProfile("xy");
  p2a->SetName("p2a"); // any fe
  
  p3rf->GetZaxis()->SetRangeUser(0,0.01); // fe<0.01, bin 1
  TProfile2D *p2h = p3rf->Project3DProfile("xy");
  p2h->SetName("p2h");
  
  p3rf->GetZaxis()->SetRangeUser(0.2,0.8); // 0.1<fe<0.9, bins 11-30
  TProfile2D *p2e = p3rf->Project3DProfile("xy");
  p2e->SetName("p2e");
  

  TH2D *h2ae = p2a->ProjectionXY("h2ae","e");
  TH2D *h2he = p2h->ProjectionXY("h2he","e");
  TH2D *h2ee = p2e->ProjectionXY("h2ee","e");

  p2a->SetErrorOption("s");
  p2h->SetErrorOption("s");
  p2e->SetErrorOption("s");
  TH2D *h2a = p2a->ProjectionXY("h2a","e");
  TH2D *h2h = p2h->ProjectionXY("h2h","e");
  TH2D *h2e = p2e->ProjectionXY("h2e","e");

  for (int i = 1; i != h2a->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2a->GetNbinsY()+1; ++j) {
      h2a->SetBinContent(i, j, h2a->GetBinError(i, j));
      h2h->SetBinContent(i, j, h2h->GetBinError(i, j));
      h2e->SetBinContent(i, j, h2e->GetBinError(i, j));

      h2a->SetBinError(i, j, 0.);
      h2h->SetBinError(i, j, 0.);
      h2e->SetBinError(i, j, 0.);
    /*      
      if (h2a->GetBinError(i)!=0)
	h2a->SetBinContent(i, j, h2ae->GetBinError(i)/h2a->GetBinError(i));
      if (h2h->GetBinError(i)!=0)
	h2h->SetBinContent(i, j, h2he->GetBinError(i)/h2h->GetBinError(i));
      if (h2e->GetBinError(i)!=0)
	h2e->SetBinContent(i, j, h2ee->GetBinError(i)/h2e->GetBinError(i));
    */
    } // for j
  } // for i

  TLine *l = new TLine();
  l->SetLineColor(kGray+2);
  l->SetLineStyle(kDashed);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  TF1 *f1 = new TF1("f1","[0]+[1]*pow(log(x),[2])",1.0,30);
  f1->SetParameters(1.305,0.8,0.5);
  f1->SetLineColor(kGray+2);
  
  TCanvas *c3s = new TCanvas("c3s","c3s",1800,600);
  c3s->Divide(3,1,0,0);
  
  c3s->cd(1);
  gPad->SetLogx();
  
  h2e->UseCurrentStyle();
  h2e->GetXaxis()->SetMoreLogLabels();
  h2e->GetXaxis()->SetNoExponent();
  h2e->GetXaxis()->SetRangeUser(0.7,300.);
  h2e->GetYaxis()->SetRangeUser(0,3.139);
  h2e->GetZaxis()->SetRangeUser(0,1.0);

  tdrDraw(h2e,"COL");
  
  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  l->SetLineStyle(kDotted);
  l->DrawLine(0.7,1.305,300,1.305);
  l->DrawLine(0.7,1.566,300,1.566);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.75,0.92,"E-hadrons");
  
  gPad->RedrawAxis();
  
  c3s->cd(2);
  gPad->SetLogx();

  h2h->UseCurrentStyle();
  h2h->GetXaxis()->SetMoreLogLabels();
  h2h->GetXaxis()->SetNoExponent();
  h2h->GetXaxis()->SetRangeUser(0.7,300.);
  h2h->GetYaxis()->SetRangeUser(0,3.139);
  h2h->GetZaxis()->SetRangeUser(0,1.0);

  tdrDraw(h2h,"COL");

  l->SetLineStyle(kSolid);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  l->SetLineStyle(kDotted);
  l->DrawLine(0.7,1.305,300,1.305);
  l->DrawLine(0.7,1.566,300,1.566);
  f1->Draw("SAME");
  
  tex->DrawLatex(0.70,0.92,"H-hadrons");
  
  gPad->RedrawAxis();

  c3s->cd(3);
  gPad->SetLogx();
  gPad->SetRightMargin(0.20);

  h2a->UseCurrentStyle();
  h2a->GetXaxis()->SetMoreLogLabels();
  h2a->GetXaxis()->SetNoExponent();
  h2a->GetXaxis()->SetRangeUser(0.7,300.);
  h2a->GetYaxis()->SetRangeUser(0,3.139);
  h2a->GetZaxis()->SetRangeUser(0,1.0);
  h2a->GetZaxis()->SetTitle("Resolution (RMS)");

  tdrDraw(h2a,"COLZ");

  l->SetLineStyle(kSolid);
  l->DrawLine(3.5,0,3.5,1.479);
  l->DrawLine(2.5,1.479,2.5,3.139);
  l->DrawLine(0.7,1.479,300,1.479);
  l->DrawLine(0.7,2.500,300,2.500);
  l->SetLineStyle(kDotted);
  l->DrawLine(0.7,1.305,300,1.305);
  l->DrawLine(0.7,1.566,300,1.566);
  f1->Draw("SAME");

  tex->DrawLatex(0.52,0.92,"All hadrons");

  gPad->RedrawAxis();
  
  c3s->SaveAs("pdf/drawPiongunResolution_res2D_x3.pdf");
    

  
  // Return fe cut to nominal
  p3rf->GetZaxis()->SetRangeUser(0,-1);
  
  // Set the abseta range to |eta| < 1.305
  p3rf->GetXaxis()->SetRangeUser(0.0, 1.305);
  //p3rf->GetXaxis()->SetRangeUser(2.5, 2.964);
  
  // Project p3rf onto genPt (y) and fe (z) axes
  TProfile2D *p2_genPt_fe_resp = p3rf->Project3DProfile("zy");//"yz");
  
  // Now project onto genPt axis, integrating over fe
  // ProfileX() operates along the x-axis (which corresponds to genPt in p2_genPt_fe_resp)
  //p2_genPt_fe_resp->GetYaxis()->SetRangeUser(0.0, 0.01);
  TProfile *p_genPt_resp = p2_genPt_fe_resp->ProfileX("p_all", 0, -1, "s");
  TProfile *p_genPt_resp_hh = p2_genPt_fe_resp->ProfileX("p_hh", 1, 1, "s");
  TProfile *p_genPt_resp_eh = p2_genPt_fe_resp->ProfileX("p_eh", 10, 30, "s");
  
  // Ensure that the error option is set to "s" to store RMS as bin errors
  p_genPt_resp->SetErrorOption("s");
  
  // Create a histogram to store RMS vs genPt
  int nGenPtBins = p_genPt_resp->GetNbinsX();
  TH1D *hRMSvsGenPt = new TH1D("hRMSvsGenPt", ";genPt [GeV];Resolution (RMS of resp)", nGenPtBins, p_genPt_resp->GetXaxis()->GetXbins()->GetArray());
  TH1D *hRMSvsGenPt_hh = new TH1D("hRMSvsGenPt_hh", ";genPt [GeV];Resolution (RMS of resp)", nGenPtBins, p_genPt_resp->GetXaxis()->GetXbins()->GetArray());
  TH1D *hRMSvsGenPt_eh = new TH1D("hRMSvsGenPt_eh", ";genPt [GeV];Resolution (RMS of resp)", nGenPtBins, p_genPt_resp->GetXaxis()->GetXbins()->GetArray());
  
  // Loop over genPt bins to extract the RMS from bin errors
  for (int iGenPt = 1; iGenPt <= nGenPtBins; ++iGenPt) {
    double RMS = p_genPt_resp->GetBinError(iGenPt);
    double RMS_hh = p_genPt_resp_hh->GetBinError(iGenPt);
    double RMS_eh = p_genPt_resp_eh->GetBinError(iGenPt);
    hRMSvsGenPt->SetBinContent(iGenPt, RMS);
    hRMSvsGenPt_hh->SetBinContent(iGenPt, RMS_hh);
    hRMSvsGenPt_eh->SetBinContent(iGenPt, RMS_eh);
    //hRMSvsGenPt->SetBinError(iGenPt, 0); // Optional: set error bars if needed
  }
  
  // Draw the histogram using TDR style
  TH1D *h = tdrHist("h","RMS",0,1.0,"p_{T,gen} (GeV)",0.2,1000.);
  lumi_136TeV = "Winter24 PionGun";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1", h, 8, 11, kSquare);
  gPad->SetLogx();
  
  tdrDraw(hRMSvsGenPt, "Pz", kFullCircle, kBlack);
  tdrDraw(hRMSvsGenPt_hh, "Pz", kOpenSquare, kRed);
  tdrDraw(hRMSvsGenPt_eh, "Pz", kOpenCircle, kBlue);
  
  TF1 *f1a = new TF1("f1a","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
		     6.0,500.);
  f1a->SetParameters(1,1,0.1);
  hRMSvsGenPt->Fit(f1a,"QRNW");
  f1a->SetLineColor(kBlack);
  f1a->Draw("SAME");
  
  TF1 *f1hh = new TF1("f1hh","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
		      4.0,500.);
  f1hh->SetParameters(1,1,0.1);
  hRMSvsGenPt_hh->Fit(f1hh,"QRNW");
  f1hh->SetLineColor(kRed);
  f1hh->Draw("SAME");
  
  
  TF1 *f1eh = new TF1("f1eh","sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
		      8.0,500.);
  f1eh->SetParameters(1,1,0.1);
  hRMSvsGenPt_eh->Fit(f1eh,"QRNW");
  f1eh->SetLineColor(kBlue);
  f1eh->Draw("SAME");
  
  
  // Save the canvas as a PDF
  c1->SaveAs("pdf/drawPiongunResolution.pdf");


  /////////////////////////////////////////////////////////////////////
  // Complete 2D analysis in eta(x), pT(y) for A(ll), H, EH; 30 bins //
  /////////////////////////////////////////////////////////////////////
  
  TCanvas *c3 = new TCanvas("c5","c5",6*300,5*300); // 1800x1500
  c3->Divide(6,5,0,0);
  TCanvas *c3r = new TCanvas("c5r","c5r",6*300,5*300); // 1800x1500
  c3r->Divide(6,5,0,0);
  const int neta = 6*5;

  TProfile3D *p3 = p3rf; // reuse code from drawPiongun.C
  TProfile3D *p3b = p3rf2;
  TLegend *legi(0);
  for (int ieta = 1; ieta != p3->GetNbinsX()+1; ++ieta) {
    double eta = p3->GetXaxis()->GetBinCenter(ieta);
    double deta = 0.5*p3->GetXaxis()->GetBinWidth(ieta);

    // There are hadrons up to (eta>3.139), but these are partially in HF
    if (eta>2.964) continue;

    // Return fe cut to nominal
    p3->GetZaxis()->SetRangeUser(0,-1);
    p3b->GetZaxis()->SetRangeUser(0,-1);

    // Select eta bin
    p3->GetXaxis()->SetRange(ieta,ieta);
    p3b->GetXaxis()->SetRange(ieta,ieta);

    // Project p3rf onto genPt (y) and fe (z) axes
    TProfile2D *p2 = p3->Project3DProfile("zy");
    p2->SetName(Form("p2yz_%d",ieta));
    TProfile2D *p2b = p3b->Project3DProfile("zy");
    p2b->SetName(Form("p2byz_%d",ieta));
					  
    // All, H-hadrons and "average" EH-hadrons directly vs pT
    // Need "o" (=original) option to not mess up variable binned axis
    // Need "s" option to use RMS as error so we can extract JER

    // Range 0,-1 projects all fe bins
    TProfile *pa = p2->ProfileX(Form("pa_%d",ieta),0,-1,"os");
    TProfile *pab = p2b->ProfileX(Form("pab_%d",ieta),0,-1,"os");

    // Range [1,1] effectively picks fe<0.01, so below MIP cut up to 100 GeV
    // at higher pT still useful to cut out delta rays(?) above MIP levels
    // there is strong fe dependence at low fe due to H-EH mixture so cut tight
    TProfile *ph = p2->ProfileX(Form("ph_%d",ieta),1,1,"os");
    TProfile *phb = p2b->ProfileX(Form("phb_%d",ieta),1,1,"os");

    // Use less biased EH hadrons in the central range fe [0,2.0.8]:
    // 0<fe<0.2 is mix of H and EH where E was lost, and very variable vs fe
    // fe>0.8 is sometimes ok, sometimes with H that was lost, so not linear
    int i1 = p2->GetYaxis()->FindBin(0.2);
    int i2 = p2->GetYaxis()->FindBin(0.8);
    TProfile *pe = p2->ProfileX(Form("pe_%d",ieta),i1,i2,"os");
    TProfile *peb = p2b->ProfileX(Form("peb_%d",ieta),i1,i2,"os");

    // Ensure that the error option is set to "s" to store RMS as bin errors
    //pa->SetErrorOption("s");
    //pab->SetErrorOption("s");
    //ph->SetErrorOption("s");
    //phb->SetErrorOption("s");
    //pe->SetErrorOption("s");
    //peb->SetErrorOption("s");
  
    // Create a histogram to store RMS vs genPt
    int nx = pa->GetNbinsX();
    const double *vx = pa->GetXaxis()->GetXbins()->GetArray();
    TH1D *ha, *hh, *he;
    ha = new TH1D(Form("ha_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
		  nx, vx);
    hh = new TH1D(Form("hh_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
		  nx, vx);
    he = new TH1D(Form("he_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
		  nx, vx);

    TH1D *hab, *hhb, *heb;
    hab = new TH1D(Form("hab_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
			nx, vx);
    hhb = new TH1D(Form("hhb_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
			nx, vx);
    heb = new TH1D(Form("heb_%d",ieta),";p_{T,gen} (GeV);Resolution (RMS)",
			nx, vx);
    
    // Loop over genPt bins to extract the RMS from bin errors
    for (int i = 1; i != nx; ++i) {
      double RMS_a = pa->GetBinError(i);
      double RMS_h = ph->GetBinError(i);
      double RMS_e = pe->GetBinError(i);
      ha->SetBinContent(i, RMS_a);
      hh->SetBinContent(i, RMS_h);
      he->SetBinContent(i, RMS_e);

      double RMS_ab = pab->GetBinError(i);
      double RMS_hb = phb->GetBinError(i);
      double RMS_eb = peb->GetBinError(i);
      hab->SetBinContent(i, RMS_ab);
      hhb->SetBinContent(i, RMS_hb);
      heb->SetBinContent(i, RMS_eb);
      // TBD: uncertainties
    }
  
    // Draw the histogram using TDR style
    c3->cd(ieta);
    gPad->SetLogx();

    double eps = 1e-3;
    TH1D *h = tdrHist(Form("h_%d",ieta),"Resolution (RMS)",0+eps,1.0-eps,
		      "p_{T,gen} (GeV)",0.2*(1+eps),1000.*(1-eps));
    h->Draw();

    tex->SetTextSize(0.045*1.5);
    tex->DrawLatex(0.50,0.87,Form("%1.3f#leq|#eta|<%1.3f",eta-deta,eta+deta));

    tdrDraw(hhb, "HIST", kNone, kRed, kSolid, -1, kNone);
    tdrDraw(heb, "HIST", kNone, kBlue, kSolid, -1, kNone);
    tdrDraw(hab, "HIST", kNone, kBlack, kSolid, -1, kNone);
    
    tdrDraw(hh, "Pz", kOpenSquare, kRed-9, kSolid, -1, kNone, 0, 0.8);
    tdrDraw(he, "Pz", kOpenCircle, kBlue-9, kSolid, -1, kNone, 0, 0.8);
    tdrDraw(ha, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.8);
  
    TF1 *f1hh = new TF1(Form("f1hh_%d",ieta),
			"sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
			4.0,500.);
    f1hh->SetParameters(1,1,0.1);
    hh->Fit(f1hh,"QRNW");
    f1hh->SetLineColor(kRed-9);
    TF1 *f1hh_0 = (TF1*)f1hh->DrawClone("SAME");
    f1hh->SetLineStyle(kDotted);
    f1hh->SetRange(0.2,1000.);
    f1hh->Draw("SAME");
    
    TF1 *f1eh = new TF1("f1eh",
			"sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
			8.0,500.);
    f1eh->SetParameters(1,1,0.1);
    he->Fit(f1eh,"QRNW");
    f1eh->SetLineColor(kBlue-9);
    TF1 *f1eh_0 = (TF1*)f1eh->DrawClone("SAME");
    f1eh->SetLineStyle(kDotted);
    f1eh->SetRange(0.2,1000.);
    f1eh->Draw("SAME");

    TF1 *f1a = new TF1(Form("f1a_%d",ieta),
		       "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
		       6.0, 500.);
    f1a->SetParameters(1,1,0.1);
    ha->Fit(f1a,"QRNW");
    f1a->SetLineColor(kBlack);
    TF1 *f1a_0 = (TF1*)f1a->DrawClone("SAME");
    f1a->SetLineStyle(kDotted);
    f1a->SetRange(0.2,1000.);
    f1a->Draw("SAME");

    gPad->RedrawAxis();
    
    // Draw legend in the last empty pad
    if (ieta==1) {
      c3->cd(neta);

      TLegend *leg = tdrLeg(0.05,0.95-2*6*0.045,0.80,0.95);
      leg->SetTextSize(1.5*0.045);
      leg->AddEntry(ha,"All hadrons","PLE");
      leg->AddEntry(hh,"H hadrons","PLE");
      leg->AddEntry(he,"EH hadrons","PLE");
      leg->AddEntry(hab,"All (PFHC)","F");
      leg->AddEntry(hhb,"H (PFHC)","F");
      leg->AddEntry(heb,"EH (PFHC)","F");
    }

    // Store results also in an individual canvas for better eta bin control
    if (true) { // single canvas plotting

      double eps = 1e-3;
      TH1D *h3i = tdrHist(Form("h3i_%d",ieta),"Resolution (RMS)",0+eps,1-eps,
			  "p_{T,gen} (GeV)",0.2*(1+eps),1000.*(1-eps));
      TCanvas *c3i = tdrCanvas(Form("c3_%d",ieta),h3i,8,11,kSquare);
      gPad->SetLogx();

      tex->SetTextSize(0.045);
      tex->DrawLatex(0.50,0.87,Form("%1.3f#leq|#eta|<%1.3f",eta-deta,eta+deta));

      tdrDraw(hhb, "HIST", kNone, kRed, kSolid, -1, kNone);
      tdrDraw(heb, "HIST", kNone, kBlue, kSolid, -1, kNone);
      tdrDraw(hab, "HIST", kNone, kBlack, kSolid, -1, kNone);

      f1hh_0->Draw("SAME");
      f1eh_0->Draw("SAME");
      f1a_0->Draw("SAME");
      f1hh->Draw("SAME");
      f1eh->Draw("SAME");
      f1a->Draw("SAME");
      
      tdrDraw(hh, "Pz", kOpenSquare, kRed-9, kSolid, -1, kNone, 0, 0.8);
      tdrDraw(he, "Pz", kOpenCircle, kBlue-9, kSolid, -1, kNone, 0, 0.8);
      tdrDraw(ha, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.8);


      // For single pads
      if (ieta==1) {
	legi = tdrLeg(0.60,0.80-6*0.045,0.85,0.80);
	legi->AddEntry(ha,"All hadrons","PLE");
	legi->AddEntry(hh,"H hadrons","PLE");
	legi->AddEntry(he,"EH hadrons","PLE");
	legi->AddEntry(hab,"All (PFHC)","F");
	legi->AddEntry(hhb,"H (PFHC)","F");
	legi->AddEntry(heb,"EH (PFHC)","F");
      }
      legi->Draw();

      gPad->RedrawAxis();
      c3i->SaveAs(Form("pdf/vsEta/drawPiongunResolution_eta_%04.0f_%04.0f.pdf",
		       1000*(eta-deta),1000.*(eta+deta)));
      
    } // if single canvas
    
  } // for ieta
  
  c3->SaveAs("pdf/drawPiongunResolution_2D.pdf");
  
  // Clean up
  //f->Close();
  //delete f;
}
