#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize, char *text) {
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetTextFont(52); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void drawCMS(Double_t x,Double_t y, bool prel) {

  float cmsTextSize = 0.07;
  float extraOverCmsTextSize = 0.68;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  TLatex l;
  l.SetTextSize(cmsTextSize); 
  l.SetTextFont(61); 
  l.SetTextAngle(0);
  l.SetNDC();
  l.SetTextColor(1);
  l.DrawLatex(x,y,"CMS");

  TLatex lp;
  lp.SetTextSize(extraTextSize); 
  lp.SetTextFont(52); 
  lp.SetNDC();
  lp.SetTextColor(1);
  float offset = 0.11;
  lp.DrawLatex(x+offset,y,"Simulation");
  
  TLatex ll;
  ll.SetTextSize(extraTextSize); 
  ll.SetTextFont(42); 
  ll.SetTextAngle(0);
  ll.SetNDC();
  ll.SetTextColor(1);
  ll.DrawLatex(0.68,y,"35.9 fb^{-1} (13 TeV)");

}


void drawCombineMatrix() {
  
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetTextFont(42);
  gStyle->SetPalette(kLightTemperature);

  
  TString var[2] = {"pt","y"};
  TString level[2] = {"","_PL"};
  
  TCanvas c;

  /*
  double length[2] {0.00, 1.00};
  //double red[2] = {0.99, 0.32};
  //double green[2] = {0.99, 0.42};
  //double blue[2] = {0.99, 0.9};
  double red[2] = {0.98, 0.0};
  double green[2] = {0.98, 0.0};
  double blue[2] = {1.0, 1.0};
  TColor::CreateGradientColorTable(2,length,red,green,blue,256);
  gStyle->SetNumberContours(256);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPaintTextFormat(".1f");
  */

  for (int iv=0; iv<2; iv++) {
    for (int il=0; il<2; il++) {
      
      TFile* file0 = TFile::Open("response_"+var[iv]+"_el"+level[il]+".root");
      TFile* file1 = TFile::Open("response_"+var[iv]+"_mu"+level[il]+".root");

      TH2D* h0 = (TH2D*) file0->Get("response_final");
      TH2D* h1 = (TH2D*) file1->Get("response_final");

      h0->Add(h1);


      // normalize so that for each bin of true top quark pt(eta), the bins in measured top pt(eta) add up to 100%
      int nbins = h0->GetNbinsX();
      for (int ir = 1; ir < nbins+1; ir++){
	double rowsum = h0->Integral(1,nbins,ir,ir);
	for (int ic = 1; ic < nbins+1; ic++){
	  //double normval = rowsum > 0.0 ? h0->GetBinContent(ic,ir) / rowsum * 100.0 : 0.0;
	  double normval = rowsum > 0.0 ? h0->GetBinContent(ic,ir) / rowsum : 0.0;
	  h0->SetBinContent(ic,ir,normval);
	}
      }

      if (var[iv] == "pt") {
	if (level[il] != "") {
	  h0->GetXaxis()->SetTitle("Detector-level top jet p_{T} [GeV]");
	  h0->GetYaxis()->SetTitle("Particle-level top jet p_{T} [GeV]");
	}
	else {
	  h0->GetXaxis()->SetTitle("Detector-level top jet p_{T} [GeV]");
	  h0->GetYaxis()->SetTitle("Parton-level top quark p_{T} [GeV]");
	}
	
	h0->SetAxisRange(401.0,1199.0,"X");
	h0->SetAxisRange(401.0,1199.0,"Y");
      }
      else {
	if (level[il] != "") {
	  h0->GetXaxis()->SetTitle("Detector-level top jet rapidity");
	  h0->GetYaxis()->SetTitle("Particle-level top jet rapidity");
	}
	else {
	  h0->GetXaxis()->SetTitle("Detector-level top jet rapidity");
	  h0->GetYaxis()->SetTitle("Parton-level top quark rapidity");
	}
      }

      h0->GetXaxis()->SetTitleOffset(1.2);
      h0->GetYaxis()->SetTitleOffset(1.3);

      h0->SetAxisRange(0.01,1,"Z");

      h0->Draw("text colz");
      //h0->Draw("text,same");

      drawCMS(0.10,0.92, false);
      mySmallText(0.15,0.82,1,0.035,"l+jets events");

      c.SaveAs("UnfoldingPlots/combined_responseMatrix_"+var[iv]+level[il]+".pdf");

      file0->TFile::Close();
      file1->TFile::Close();

    }    
  }

}
