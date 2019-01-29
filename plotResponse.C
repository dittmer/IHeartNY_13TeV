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

void drawCMS(Double_t x,Double_t y, bool prel) {

  float cmsTextSize = 0.065;
  float extraOverCmsTextSize = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  TLatex l;
  l.SetTextSize(cmsTextSize); 
  l.SetTextFont(61); 
  l.SetTextAngle(0);
  l.SetNDC();
  l.SetTextColor(1);
  l.DrawLatex(x,y,"CMS");

  if (prel) {
    TLatex lp;
    lp.SetTextSize(extraTextSize); 
    lp.SetTextFont(52); 
    lp.SetNDC();
    lp.SetTextColor(1);
    float offset = 0.09;
    lp.DrawLatex(x+offset,y,"Preliminary");
  }
}


void plotResponse(TString which="pt") {

  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.15);
  double length[2] {0.00, 1.00};
  double red[2] = {0.9, 0.0};
  double green[2] = {0.9, 0.0};
  double blue[2] = {1.0, 1.0};
  TColor::CreateGradientColorTable(2,length,red,green,blue,256);
  gStyle->SetNumberContours(256);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPaintTextFormat(".1f");
  gStyle->SetOptStat(0);

  TCanvas c;

  TString channels[2] = {"mu","el"};
  //TString matrices[2] = {"response_pt_split_TH2","response2_pt_split_TH2"};
  if (which!="pt" && which!="y") return;
  TString matrices[2] = {"response_"+which+"_TH2","response_"+which+"_TH2_PL"};

  // Louise version
  TString name_TTbarNom = "PLnew";
  TString name_TTbarNom_p2 = "v2_PLnew";
  TString name_TTbar_m700to1000 = "m700to1000_PLnew";
  TString name_TTbar_m1000toInf = "m1000toInf_PLnew";
  // Susan version
  //TString name_TTbarNom = "fullTruth_PLnew";
  //TString name_TTbarNom_p2 = "fullTruth_PLnew_p2";
  //TString name_TTbar_m700to1000 = "fullTruth_m700to1000_PLnew";
  //TString name_TTbar_m1000toInf = "fullTruth_m1000toInf_PLnew";


  for (int ich = 0; ich < 2; ich++){
    TFile* f_m0to700_p1 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_"+name_TTbarNom+"_"+channels[ich]+"_nom_post.root");
    TFile* f_m0to700_p2 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_"+name_TTbarNom_p2+"_"+channels[ich]+"_nom_post.root");
    TFile* f_m700to1000 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_"+name_TTbar_m700to1000+"_"+channels[ich]+"_nom_post.root");
    TFile* f_m1000toInf = TFile::Open("histfiles_full2016/hists_PowhegPythia8_"+name_TTbar_m1000toInf+"_"+channels[ich]+"_nom_post.root");

    for (int im = 0; im < 2; im++){
      TH2F* h_m0to700_p1 = (TH2F*) f_m0to700_p1->Get(matrices[im]);
      TH2F* h_m0to700_p2 = (TH2F*) f_m0to700_p2->Get(matrices[im]);
      TH2F* h_m700to1000 = (TH2F*) f_m700to1000->Get(matrices[im]);
      TH2F* h_m1000toInf = (TH2F*) f_m1000toInf->Get(matrices[im]);
      
      TH2F* h_m0to700 = (TH2F*) h_m0to700_p1->Clone();
      h_m0to700->Add(h_m0to700_p2);
      h_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
      h_m700to1000->Scale(831.76 * 35867.0 * 0.0967 / 38578334.0);
      h_m1000toInf->Scale(831.76 * 35867.0 * 0.0256 / 24495211.0);
    
      TH2F* h_sum = (TH2F*) h_m0to700->Clone();
      h_sum->Add(h_m700to1000);
      h_sum->Add(h_m1000toInf);

      // Plot response matrix
      h_sum->Draw("colz");

      //TString append = (matrices[im].Contains("response2")) ? "_ptFull" : "";
      TString append = "";
      if (matrices[im].Contains("response2")) append="_ptFull";
      if (matrices[im].Contains("PL")) append="_PL";
      //c.SaveAs("UnfoldingPlots/unfold_responseMatrix_"+which+"_"+channels[ich]+append+".pdf");  

      // Plot normalized response matrix
      for (int ir = 1; ir < h_sum->GetNbinsY()+1; ir++){
	double rowsum = h_sum->Integral(1,h_sum->GetNbinsX(),ir,ir);
	for (int ic = 1; ic < h_sum->GetNbinsX()+1; ic++){
	  double norm = (rowsum > 0.0) ? h_sum->GetBinContent(ic,ir) / rowsum * 100.0 : 0.0;
	  h_sum->SetBinContent(ic,ir,norm);
	}
      }

      //float maxZ = (matrices[im].Contains("response2")) ? 85.0 : 60.0;
      //h_sum->SetMaximum(maxZ);
      h_sum->Draw("colz");
      c.SaveAs("UnfoldingPlots/unfold_responseMatrix_"+which+"_"+channels[ich]+append+"_norm.pdf");  
    }
  }
}
