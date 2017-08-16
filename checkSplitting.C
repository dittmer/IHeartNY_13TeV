#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void checkSplitting(){

  gStyle->SetOptStat(0);
  
  TString channels[2] = {"mu","el"};
  const int nHIST = 5;
  TString hists[nHIST] = {"genTTbarMass","genTTbarMass_full","hardTTbarMass","ptGenTop","ptRecoTop_split"};

  // Values from mInc sample
  float e_m700to1000 = 0.0967;
  float e_m1000toInf = 0.0256;

  //Values from McM
  //float e_m700to1000 = 0.0921;
  //float e_m1000toInf = 0.02474;
  
  for (int ich = 0; ich < 2; ich++){
    
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.08);

    TFile* f_mInc_p1    = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_mInc_"+channels[ich]+"_nom_post.root");
    TFile* f_mInc_p2    = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_mInc_"+channels[ich]+"_nom_post.root");
    TFile* f_m0to700_p1 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+channels[ich]+"_nom_post.root");
    TFile* f_m0to700_p2 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+channels[ich]+"_nom_post.root");
    TFile* f_m700to1000 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_m700to1000_"+channels[ich]+"_nom_post.root");
    TFile* f_m1000toInf = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_m1000toInf_"+channels[ich]+"_nom_post.root");

    for (int ih = 0; ih < nHIST; ih++){
      TH1F* h_mInc_p1    = (TH1F*) f_mInc_p1->Get(hists[ih]);
      TH1F* h_mInc_p2    = (TH1F*) f_mInc_p2->Get(hists[ih]);
      TH1F* h_m0to700_p1 = (TH1F*) f_m0to700_p1->Get(hists[ih]);
      TH1F* h_m0to700_p2 = (TH1F*) f_m0to700_p2->Get(hists[ih]);
      TH1F* h_m700to1000 = (TH1F*) f_m700to1000->Get(hists[ih]);
      TH1F* h_m1000toInf = (TH1F*) f_m1000toInf->Get(hists[ih]);

      TH1F* h_mInc = (TH1F*) h_mInc_p1->Clone();
      h_mInc->Add(h_mInc_p2);
      h_mInc->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));

      TH1F* h_m0to700 = (TH1F*) h_m0to700_p1->Clone();
      h_m0to700->Add(h_m0to700_p2);
      h_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));

      h_m700to1000->Scale(831.76 * 35867.0 * e_m700to1000 / 38578334.0);
      h_m1000toInf->Scale(831.76 * 35867.0 * e_m1000toInf / 24495211.0);

      TH1F* h_sum = (TH1F*) h_m0to700->Clone();
      h_sum->Add(h_m700to1000);
      h_sum->Add(h_m1000toInf);
      TH1F* h_ratio = (TH1F*) h_mInc->Clone();
      h_ratio->Divide(h_sum);

      THStack* h_split = new THStack();
      h_m0to700->SetFillColor(2);
      h_m0to700->SetLineWidth(1);
      h_m700to1000->SetFillColor(3);
      h_m700to1000->SetLineWidth(1);
      h_m1000toInf->SetFillColor(4);
      h_m1000toInf->SetLineWidth(1);
      h_split->Add(h_m0to700);
      h_split->Add(h_m700to1000);
      h_split->Add(h_m1000toInf);

      TCanvas* c = new TCanvas("c","c",900,800);
      TPad* p1 = new TPad("p1","p1",0,0.3,1,1);
      p1->SetTopMargin(0.08);
      p1->SetBottomMargin(0.05);
      p1->SetNumber(1);
      TPad* p2 = new TPad("p2","p2",0,0,1,0.3);
      p2->SetNumber(2);
      p2->SetTopMargin(0.05);
      p2->SetBottomMargin(0.35);
    
      p1->Draw();
      p2->Draw();
      p1->cd();
      
      h_ratio->GetXaxis()->SetTitle(h_mInc->GetXaxis()->GetTitle());
    
      h_mInc->GetYaxis()->SetLabelSize(0.05);
      h_mInc->GetYaxis()->SetTitleSize(0.05);
      h_mInc->GetYaxis()->SetTitleOffset(1.2);
      h_mInc->GetXaxis()->SetTitle("");
    
      h_mInc->Draw("");
      h_split->Draw("hist,same");
      h_mInc->Draw("same");

      // legend
      TLegend* leg;
      leg = new TLegend(0.65,0.65,0.91,0.92);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->AddEntry(h_mInc, "Inclusive", "p");
      leg->AddEntry(h_m0to700, "m0to700", "f");
      leg->AddEntry(h_m700to1000, "m700to1000", "f");
      leg->AddEntry(h_m1000toInf, "m1000toInf", "f");
      leg->Draw();

      p2->cd();
      p2->SetGridy();
      h_ratio->Draw();
      h_ratio->SetMaximum(1.1);
      h_ratio->SetMinimum(0.9);
      h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
      h_ratio->GetYaxis()->SetTitle("Inclusive / Stitched");
      h_ratio->GetXaxis()->SetLabelSize(0.1);
      h_ratio->GetYaxis()->SetLabelSize(0.1);
      h_ratio->GetXaxis()->SetTitleOffset(1.0);
      h_ratio->GetYaxis()->SetTitleOffset(0.4);
      h_ratio->GetXaxis()->SetTitleSize(0.1);
      h_ratio->GetYaxis()->SetTitleSize(0.09);

      c->SaveAs("Plots/compareStitched_"+hists[ih]+"_"+channels[ich]+".pdf");
    }

    // Now compare response matrices
    TH2F* h_mInc_p1    = (TH2F*) f_mInc_p1->Get("response_pt_split_TH2");
    TH2F* h_mInc_p2    = (TH2F*) f_mInc_p2->Get("response_pt_split_TH2");
    TH2F* h_m0to700_p1 = (TH2F*) f_m0to700_p1->Get("response_pt_split_TH2");
    TH2F* h_m0to700_p2 = (TH2F*) f_m0to700_p2->Get("response_pt_split_TH2");
    TH2F* h_m700to1000 = (TH2F*) f_m700to1000->Get("response_pt_split_TH2");
    TH2F* h_m1000toInf = (TH2F*) f_m1000toInf->Get("response_pt_split_TH2");
    
    TH2F* h_mInc = (TH2F*) h_mInc_p1->Clone();
    h_mInc->Add(h_mInc_p2);
    h_mInc->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
    
    TH2F* h_m0to700 = (TH2F*) h_m0to700_p1->Clone();
    h_m0to700->Add(h_m0to700_p2);
    h_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
    
    h_m700to1000->Scale(831.76 * 35867.0 * e_m700to1000 / 38578334.0);
    h_m1000toInf->Scale(831.76 * 35867.0 * e_m1000toInf / 24495211.0);
    
    TH2F* h_sum = (TH2F*) h_m0to700->Clone();
    h_sum->Add(h_m700to1000);
    h_sum->Add(h_m1000toInf);

    TH2F* h_mInc_err = (TH2F*) h_mInc->Clone();
    TH2F* h_sum_err = (TH2F*) h_sum->Clone();
    for (int ibiny = 1; ibiny < h_mInc->GetNbinsY()+1; ibiny++){
      for (int ibinx = 1; ibinx < h_mInc->GetNbinsX()+1; ibinx++){
	h_mInc_err->SetBinContent(ibinx,ibiny,h_mInc->GetBinError(ibinx,ibiny)/h_mInc->GetBinContent(ibinx,ibiny));
	h_sum_err->SetBinContent(ibinx,ibiny,h_sum->GetBinError(ibinx,ibiny)/h_sum->GetBinContent(ibinx,ibiny));
	
      }
    }
    
    for (int ibiny = 1; ibiny < h_mInc->GetNbinsY()+1; ibiny++){
      double rowsum_mInc = h_mInc->Integral(1,h_mInc->GetNbinsX(),ibiny,ibiny);
      double rowsum = h_sum->Integral(1,h_mInc->GetNbinsX(),ibiny,ibiny);
      for (int ibinx = 1; ibinx < h_mInc->GetNbinsX()+1; ibinx++){
	double norm_mInc = rowsum_mInc > 0.0 ? h_mInc->GetBinContent(ibinx,ibiny) / rowsum_mInc * 100.0 : 0.0;
	double norm = rowsum > 0.0 ? h_sum->GetBinContent(ibinx,ibiny) / rowsum * 100.0 : 0.0;
	h_mInc->SetBinContent(ibinx,ibiny,norm_mInc);
	h_sum->SetBinContent(ibinx,ibiny,norm);
      }
    }

    TH2F* dummy = (TH2F*) h_mInc->Clone();
    dummy->Reset();
    
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

    TCanvas c2;
    h_mInc->SetMaximum(60.);
    h_mInc->Draw("colz");
    dummy->Draw("same");
    c2.SaveAs("UnfoldingPlots/responseMatrix_mInc_"+channels[ich]+".pdf");

    h_sum->SetMaximum(60.);
    h_sum->Draw("colz");
    dummy->Draw("same");
    c2.SaveAs("UnfoldingPlots/responseMatrix_stitched_"+channels[ich]+".pdf");
    
    h_mInc_err->Draw("colz");
    dummy->Draw("same");
    c2.SaveAs("UnfoldingPlots/responseMatrixErr_mInc_"+channels[ich]+".pdf");
    
    h_sum_err->Draw("colz");
    dummy->Draw("same");
    c2.SaveAs("UnfoldingPlots/responseMatrixErr_stitched_"+channels[ich]+".pdf");
  }
}
