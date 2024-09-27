// Draw response shape for all hadrons in each eta, pT bin
#include "TFile.h"
#include "TH3D.h"

#include "tdrstyle_mod22.C"

void drawPiongunShapes() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();
  
  TFile *f = new TFile("piongun_nhcorrfast_v2.root","READ");
  assert(f && !f->IsZombie());

  TFile *f2 = new TFile("piongun_pfhc_v2.root","READ");
  assert(f2 && !f2->IsZombie());

  curdir->cd();

  TH3D *h3, *h3b;
  h3 = (TH3D*)f->Get("h3r");
  h3b = (TH3D*)f2->Get("h3r");

  // Step 1. Single plot of one |eta|, pT bin
  
  double eps = 1e-4;
  //double minf = 1e-6; // 1e-6 for hight pT wide bins, 1e-4 otherwise
  double minf = 1e-4; // 1e-6 for hight pT wide bins, 1e-4 otherwise
  TH1D *h = tdrHist("h","Fraction",minf,0.07-eps,"Response",
		    //-0.05+eps,2-eps); // for high pT
		    -0.05+eps,4-eps); // for low pT
  extraText = "Private";
  lumi_136TeV = "Winter24 SinglePionGun";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  int ix1 = h3->GetXaxis()->FindBin(2.5);//2.65);//0.);
  double x1 = h3->GetXaxis()->GetBinLowEdge(ix1);
  int ix2 = h3->GetXaxis()->FindBin(2.65)-1;//2.853)-1;//1.305)-1;
  double x2 = h3->GetXaxis()->GetBinLowEdge(ix2+1);
  int iy1 = h3->GetYaxis()->FindBin(1.0);//64);//2.0);//400.);//1.4);//2.2);//499.);//2.75);
  double y1 = h3->GetYaxis()->GetBinLowEdge(iy1);
  int iy2 = h3->GetYaxis()->FindBin(1.5)-1;//74)-1;//2.5)-1;//500.);
  double y2 = h3->GetYaxis()->GetBinLowEdge(iy2+1);

  h3->GetXaxis()->SetRange(ix1,ix2);
  h3->GetYaxis()->SetRange(iy1,iy2);
  TH1D *h1 = (TH1D*)h3->Project3D("z");
  h1->SetName("h1");

  h3b->GetXaxis()->SetRange(ix1,ix2);
  h3b->GetYaxis()->SetRange(iy1,iy2);
  TH1D *h1b = (TH1D*)h3b->Project3D("z");
  h1b->SetName("h1b");
  
  if (h1->Integral()!=0) h1->Scale(1./h1->Integral());
  if (h1b->Integral()!=0) h1b->Scale(1./h1b->Integral());
  
  tdrDraw(h1b,"HIST",kNone,kBlue+1,kSolid,-1,1000,kBlue-9);
  tdrDraw(h1,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.5);

  if (x1!=0)
    tex->DrawLatex(0.65,0.85,Form("%1.3f#leq|#eta|<%1.3f",x1,x2));
  else
    tex->DrawLatex(0.65,0.85,Form("|#eta| < %1.3f",x2));

  //tex->DrawLatex(0.65,0.80,Form("%1.0f#leqp_{T,gen}<%1.0f GeV",y1,y2));
  if (min(y1,y2)<10.)
    tex->DrawLatex(0.65,0.80,Form("[%1.1f,%1.1f] GeV",y1,y2));
  else
    tex->DrawLatex(0.65,0.80,Form("[%1.0f,%1.0f] GeV",y1,y2));

  TLegend *leg = tdrLeg(0.66,0.78-0.05*2,0.91,0.78);
  leg->AddEntry(h1,"New PFEC","PLE");
  leg->AddEntry(h1b,"Old PFHC","F");
  
  gPad->RedrawAxis();

  //c1->SaveAs("pdf/drawPiongunShapes_liny.pdf");
  c1->SaveAs(Form("pdf/drawPiongunShapes_liny_eta_%04.0f_%04.0f_pt_%04.0f_%04.0f.pdf",1000.*x1,1000.*x2,y1,y2));
  c1->SetLogy();
  //c1->SaveAs("pdf/drawPiongunShapes_logy.pdf");
  c1->SaveAs(Form("pdf/drawPiongunShapes_logy_eta_%04.0f_%04.0f_pt_%04.0f_%04.0f.pdf",1000.*x1,1000.*x2,y1,y2));
  
  int ieta1 = h3->GetXaxis()->FindBin(2.5);//2.65);//0.);
  double eta1 = h3->GetXaxis()->GetBinLowEdge(ieta1);
  int ieta2 = h3->GetXaxis()->FindBin(2.65)-1;//2.853)-1;//1.305)-1;
  double eta2 = h3->GetXaxis()->GetBinLowEdge(ieta2+1);
  
  const int npt = h3->GetNbinsY();
  h3->GetXaxis()->SetRange(ieta1,ieta2);
  h3->GetYaxis()->SetRange(0,-1);
  h3->GetZaxis()->SetRange(0,-1);
  TH1D *hpt = (TH1D*)h3->Project3D("y");
  hpt->SetName(Form("hpt_%d",ieta1));
  hpt->Reset();
  TH1D *hptr = (TH1D*)h3->Project3D("y");
  hptr->SetName(Form("hptr_%d",ieta1));
  hptr->Reset();

  h3b->GetXaxis()->SetRange(ieta1,ieta2);
  h3b->GetYaxis()->SetRange(0,-1);
  h3b->GetZaxis()->SetRange(0,-1);
  TH1D *hptb = (TH1D*)h3b->Project3D("y");
  hptb->SetName(Form("hptb_%d",ieta1));
  hptb->Reset();
  TH1D *hptrb = (TH1D*)h3b->Project3D("y");
  hptrb->SetName(Form("hptrb_%d",ieta1));
  hptrb->Reset();

  for (int ipt = 0; ipt != npt; ++ipt) {

    h3->GetXaxis()->SetRange(ieta1,ieta2);
    h3->GetYaxis()->SetRange(ipt,ipt);
    h3->GetZaxis()->SetRange(0,-1);
    TH1D *h1 = (TH1D*)h3->Project3D("z");
    h1->SetName(Form("h1_%d_%d",ieta1,ipt));

    h1->GetXaxis()->SetRangeUser(0.1,10);//1.9);
    hpt->SetBinContent(ipt, h1->GetRMS());
    hpt->SetBinError(ipt, 1e-5);
    hptr->SetBinContent(ipt, h1->GetMean());
    hptr->SetBinError(ipt, 1e-5);
    
    h3b->GetXaxis()->SetRange(ieta1,ieta2);
    h3b->GetYaxis()->SetRange(ipt,ipt);
    h3b->GetZaxis()->SetRange(0,-1);
    TH1D *h1b = (TH1D*)h3b->Project3D("z");
    h1b->SetName(Form("h1b_%d_%d",ieta1,ipt));

    h1b->GetXaxis()->SetRangeUser(0.1,10);//1.9);
    hptb->SetBinContent(ipt, h1b->GetRMS());
    hptb->SetBinError(ipt, 1e-5);
    hptrb->SetBinContent(ipt, h1b->GetMean());
    hptrb->SetBinError(ipt, 1e-5);
  } // for ipt

  

  TH1D *h2 = tdrHist("h2","Resolution (RMS)",0+eps,1-eps,
		     "p_{T,gen} (GeV)",0.2*(1+eps),1000*(1-eps));
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);
  gPad->SetLogx();

  tdrDraw(hptb,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hpt,"PEz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);

  if (eta1!=0)
    tex->DrawLatex(0.50,0.85,Form("%1.3f#leq|#eta|<%1.3f",eta1,eta2));
  else
    tex->DrawLatex(0.50,0.85,Form("|#eta| < %1.3f",eta2));

  TLegend *leg2 = tdrLeg(0.50,0.80-0.05*2,0.75,0.80);
  leg2->AddEntry(hptr,"New PFEC","PLE");
  leg2->AddEntry(hptrb,"Old PFHC","F");
  
  gPad->RedrawAxis();
  c2->SaveAs(Form("pdf/vsEta/drawPionGunShapes_eta_%04.0f_%04.0f.pdf",
		  1000.*eta1,1000.*eta2));

  
  TH1D *h_3 = tdrHist("h_3","Response (mean)",0.+eps,3.0-eps,
		     "p_{T,gen} (GeV)",0.2*(1+eps),1000*(1-eps));
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);
  gPad->SetLogx();

  l->DrawLine(0.2,1,1000.,1);
  
  tdrDraw(hptrb,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hptr,"PEz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);

  if (eta1!=0)
    tex->DrawLatex(0.50,0.85,Form("%1.3f#leq|#eta|<%1.3f",eta1,eta2));
  else
    tex->DrawLatex(0.50,0.85,Form("|#eta| < %1.3f",eta2));

  TLegend *leg3 = tdrLeg(0.50,0.80-0.05*2,0.75,0.80);
  leg3->AddEntry(hptr,"New PFEC","PLE");
  leg3->AddEntry(hptrb,"Old PFHC","F");
  
  gPad->RedrawAxis();
  c3->SaveAs(Form("pdf/vsEta/drawPionGunShapes_eta_%04.0f_%04.0f_resp.pdf",
		  1000.*eta1,1000.*eta2));
  
} // void drawPiongunShapes
