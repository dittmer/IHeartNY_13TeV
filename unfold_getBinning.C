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

int get400bin(TH1* pt_hist){
  for (int ib = 1; ib < pt_hist->GetNbinsX(); ib++){
    if (pt_hist->GetXaxis()->GetBinLowEdge(ib) == 400.0) return ib;
  }
  return -1;
}
void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize,char *text); 

void getPurityStabilityEfficiency(TH2D* h_response, TString which, int ibinlow, int ibinhigh, double &purity, double &stability, double &efficiency){

  int ibin_low = (which == "pt") ? get400bin(h_response) : 1;
  int ibin_high = h_response->GetNbinsX();

  double ngen = h_response->Integral(ibin_low,ibin_high,ibinlow,ibinhigh);
  double ngenfull = h_response->Integral(0,ibin_high+1,ibinlow,ibinhigh);
  double nrec = h_response->Integral(ibinlow,ibinhigh,ibin_low,ibin_high);
  double nrecfull = h_response->Integral(ibinlow,ibinhigh,0,ibin_high+1);
  double nrecgen = h_response->Integral(ibinlow,ibinhigh,ibinlow,ibinhigh);
  
  if (ngen == 0 || nrec == 0) {
    purity = stability = efficiency = -1.0;
    return;
  }
  
  stability = nrecgen/ngen;
  purity = nrecgen/nrec;
  efficiency = nrecfull/ngenfull;
  
  return;
}

//Fill h_output with contents of h_input
void rebinTH2(TH2D* h_input, TH2D* h_output){
  int outbinx = 0;
  int inbinx_first = 0;
  for (int ibx = 0; ibx < h_input->GetXaxis()->GetNbins()+2; ibx++){
    if (((float)h_input->GetXaxis()->GetBinUpEdge(ibx) == (float)h_output->GetXaxis()->GetBinUpEdge(outbinx)) ||
	(ibx == h_input->GetXaxis()->GetNbins()+1 && outbinx == h_output->GetXaxis()->GetNbins()+1)){ //Overflow bins will not have same upper edge
      int outbiny = 0;
      int inbiny_first = 0;
      for (int iby = 0; iby < h_input->GetYaxis()->GetNbins()+2; iby++){
	if (((float)h_input->GetYaxis()->GetBinUpEdge(iby) == (float)h_output->GetYaxis()->GetBinUpEdge(outbiny)) ||
	(iby == h_input->GetYaxis()->GetNbins()+1 && outbiny == h_output->GetYaxis()->GetNbins()+1)){
	  double err;
	  double sum = h_input->IntegralAndError(inbinx_first,ibx,inbiny_first,iby,err);
	  h_output->SetBinContent(outbinx,outbiny,sum);
	  h_output->SetBinError(outbinx,outbiny,err);
	  outbiny += 1;
	  inbiny_first = iby+1;
	}
      }
      outbinx += 1;
      inbinx_first = ibx+1;
    }
  }
  h_output->GetXaxis()->SetTitle(h_input->GetXaxis()->GetTitle());
  h_output->GetYaxis()->SetTitle(h_input->GetYaxis()->GetTitle());
  return;
}

void unfold_getBinning(TString channel, TString which, bool doPL = false) {
  
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  
  bool drawAll = false; // if set to true, this script will draw also response matrices with fine binning, full pt range, etc. else just final one


  TString append = (doPL) ? "_PL" : ""; // particle level response matrix 

  // Louise version
  //TString name_TTbarNom = "PLnew";
  //TString name_TTbarNom_p2 = "v2_PLnew";
  //TString name_TTbar_m700to1000 = "m700to1000_PLnew";
  //TString name_TTbar_m1000toInf = "m1000toInf_PLnew";
  // Susan version
  TString name_TTbarNom = "fullTruth";
  TString name_TTbarNom_p2 = "fullTruth_p2";
  TString name_TTbar_m700to1000 = "fullTruth_m700to1000";
  TString name_TTbar_m1000toInf = "fullTruth_m1000toInf";
  TString DIR = "histfiles_full2016_latest";

  TFile* f_ttbar_m0to700_p1 = TFile::Open(DIR+"/hists_PowhegPythia8_"+name_TTbarNom+"_"+channel+"_nom_post.root");
  TFile* f_ttbar_m0to700_p2 = TFile::Open(DIR+"/hists_PowhegPythia8_"+name_TTbarNom_p2+"_"+channel+"_nom_post.root");
  TFile* f_ttbar_m700to1000 = TFile::Open(DIR+"/hists_PowhegPythia8_"+name_TTbar_m700to1000+"_"+channel+"_nom_post.root");
  TFile* f_ttbar_m1000toInf = TFile::Open(DIR+"/hists_PowhegPythia8_"+name_TTbar_m1000toInf+"_"+channel+"_nom_post.root");   	

  TH2D* h_response_m0to700_p1 = (TH2D*) f_ttbar_m0to700_p1->Get("response_fine_"+which+"_TH2"+append);
  TH2D* h_response_m0to700_p2 = (TH2D*) f_ttbar_m0to700_p2->Get("response_fine_"+which+"_TH2"+append);
  TH2D* h_response_m700to1000 = (TH2D*) f_ttbar_m700to1000->Get("response_fine_"+which+"_TH2"+append);
  TH2D* h_response_m1000toInf = (TH2D*) f_ttbar_m1000toInf->Get("response_fine_"+which+"_TH2"+append);
  h_response_m0to700_p1->Sumw2();
  h_response_m0to700_p2->Sumw2();
  h_response_m700to1000->Sumw2();
  h_response_m1000toInf->Sumw2();
  TH2D* h_response_m0to700 = (TH2D*) h_response_m0to700_p1->Clone();
  h_response_m0to700->Add(h_response_m0to700_p2);
  h_response_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
  h_response_m700to1000->Scale(831.76 * 35867.0 * 0.0967 / 38578334.0);
  h_response_m1000toInf->Scale(831.76 * 35867.0 * 0.0256 / 24495211.0);
  TH2D* h_response = (TH2D*) h_response_m0to700->Clone();
  h_response->Add(h_response_m700to1000);
  h_response->Add(h_response_m1000toInf);


  //Find optimal binning
  std::vector<float> binedges;
  std::vector<float> v_stability;
  std::vector<float> v_purity;
  std::vector<float> v_efficiency;
  std::vector<float> v_res_gaus;
  std::vector<float> v_res_rms;
  std::vector<float> v_stat_unc;
  
  if (which == "pt") binedges.push_back(400.0);
  else binedges.push_back(h_response->GetXaxis()->GetBinLowEdge(1)); //Should be -2.4

  int lastbin = (which == "pt") ? get400bin(h_response) - 1 : 0;
  for (int ib = lastbin+1; ib < h_response->GetNbinsX()+1; ib++){

    // Get purity / stability / efficiency
    double purity, stability, efficiency;
    getPurityStabilityEfficiency(h_response, which, lastbin+1, ib, purity, stability, efficiency);

    // Get resolution
    float binwidth = (h_response->GetXaxis()->GetBinUpEdge(ib) - h_response->GetXaxis()->GetBinLowEdge(lastbin+1));
    TH1D* h_slice = h_response->ProjectionX(Form("slice_%i",(int)binedges.size()),lastbin+1,ib,"e"); //TODO make sure this is correct projection
    float res_gaus = 2.0 * binwidth; //Just to make sure it's too high if following fails
    if (h_slice->Integral() > 0.0){
      TF1* mygaus = new TF1("mygaus","gaus",h_response->GetXaxis()->GetBinLowEdge(1),h_response->GetXaxis()->GetBinUpEdge(h_response->GetNbinsX()));
      mygaus->SetParameters(h_slice->GetMaximum(),h_slice->GetMean(),h_slice->GetRMS());
      h_slice->Fit("mygaus","Q");
      res_gaus = mygaus->GetParameter(2);
    }
    float res_RMS = (h_slice->Integral() > 0.0) ? h_slice->GetRMS() : 2.0 * binwidth;

    double err_stat;
    double stat = h_slice->IntegralAndError(1,h_response->GetNbinsX(),err_stat);
    float stat_unc = (float) err_stat / stat;


    // Quality condition: we want bins with purity / stability > 0.5, resolution < binwidth
    bool condition = false;
    if ( (which=="y") && ( (ib%2 == 0 && purity > 0.5 && stability > 0.5 && res_gaus < binwidth && res_RMS < binwidth) || ib == h_response->GetNbinsX()) ) condition=true;
    if ( (which=="pt") && (ib == 20 || ib == 35 || ib == 50 || ib == 70 || ib == 90 || ib == 115 || ib == 170 || ib == 330) ) condition=true;
    if (condition) { 
      binedges.push_back(h_response->GetXaxis()->GetBinUpEdge(ib));
      v_stability.push_back(stability);
      v_purity.push_back(purity);
      v_efficiency.push_back(efficiency);
      v_res_gaus.push_back(res_gaus / binwidth);
      v_res_rms.push_back(res_RMS / binwidth);
      v_stat_unc.push_back(stat_unc);

      lastbin = ib;
      cout << "Added bin edge at " << h_response->GetXaxis()->GetBinUpEdge(ib) << ", lastbin is now " << lastbin << endl;
      cout << "Purity " << purity << " stability " << stability << " res_gaus " << res_gaus/binwidth << " res_RMS " << res_RMS/binwidth << endl;
    }
  }
  
  //Create new response matrix with this binning
  //Convert vector into array
  const int nbins_final = binedges.size() - 1;
  float binedges_final[nbins_final+1];
  for (int ii = 0; ii < nbins_final+1; ii++){
    binedges_final[ii] = binedges.at(ii);
  }
  
  TH2D* h_response_final = new TH2D("response_final","",nbins_final,binedges_final,nbins_final,binedges_final);
  rebinTH2(h_response,h_response_final);
  
  // Get purity / stability / efficiency
  TH1F* h_stability = new TH1F("stability","",nbins_final,binedges_final);
  TH1F* h_purity = (TH1F*) h_stability->Clone("purity");
  TH1F* h_efficiency = (TH1F*) h_stability->Clone("efficiency");
  TH1F* h_res_gaus = (TH1F*) h_stability->Clone("res_gaus");
  TH1F* h_res_rms = (TH1F*) h_stability->Clone("res_rms");
  TH1F* h_stat_unc = (TH1F*) h_stability->Clone("stat_unc");
    
  // fill hists
  for (int ib=1; ib < nbins_final+1; ib++) {
    h_stability->SetBinContent(ib, v_stability.at(ib-1));
    h_purity->SetBinContent(ib, v_purity.at(ib-1));
    h_efficiency->SetBinContent(ib, v_efficiency.at(ib-1));
    h_res_gaus->SetBinContent(ib, v_res_gaus.at(ib-1));
    h_res_rms->SetBinContent(ib, v_res_rms.at(ib-1));
    h_stat_unc->SetBinContent(ib, v_stat_unc.at(ib-1));
  }


  // ----------------------------------------------------------------------------------------------------------------
  // Plot purity / stability / efficiency
  // ----------------------------------------------------------------------------------------------------------------

  h_purity->SetAxisRange(0,1.2,"Y");
  h_purity->GetYaxis()->SetTitleOffset(1.1);
  h_purity->GetYaxis()->SetTitle("");  
  h_purity->GetYaxis()->SetTitle("Fractional");
  if (which == "pt") h_purity->GetXaxis()->SetTitle("Top p_{T} [GeV]");
  else h_purity->GetXaxis()->SetTitle("Top rapidity");
  h_purity->GetXaxis()->SetRangeUser(400.0,1199.0);

  TCanvas c;
  h_purity->SetLineColor(4); 
  h_purity->SetLineWidth(2);
  h_purity->Draw("hist");
  h_stability->SetLineColor(2); 
  h_stability->SetLineWidth(2);
  h_stability->Draw("same,hist");
  h_efficiency->SetLineColor(8); 
  h_efficiency->SetLineWidth(2);
  h_efficiency->Draw("same,hist");

  TLegend* leg;
  if (which=="y") leg = new TLegend(0.62,0.4,0.87,0.6);
  else leg = new TLegend(0.62,0.7,0.87,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_purity, " Purity", "l");
  leg->AddEntry(h_stability, " Stability", "l");
  leg->AddEntry(h_efficiency, " Efficiency", "l");
  leg->Draw();

  if (doPL) mySmallText(0.2,0.8,1,0.04,"Particle level");
  else mySmallText(0.2,0.8,1,0.04,"Parton level");

  c.SaveAs("UnfoldingPlots/purity-stability_"+channel+"_"+which+append+".pdf");


  // ------------------------------------------------------------------------------------
  //Plot resolution / statistical uncertainty
  // ------------------------------------------------------------------------------------

  h_res_gaus->SetAxisRange(0,2.0,"Y");
  h_res_gaus->GetYaxis()->SetTitleOffset(1.1);
  h_res_gaus->GetYaxis()->SetTitle("");  
  h_res_gaus->GetYaxis()->SetTitle("Fractional");
  if (which == "pt") h_res_gaus->GetXaxis()->SetTitle("Top p_{T} [GeV]");
  else h_res_gaus->GetXaxis()->SetTitle("Top rapidity");
  h_res_gaus->GetXaxis()->SetRangeUser(400.0,1199.0);
  
  TCanvas c2;
  h_res_gaus->SetLineColor(4); 
  h_res_gaus->SetLineWidth(2);
  h_res_gaus->Draw("hist");
  h_res_rms->SetLineColor(2); 
  h_res_rms->SetLineWidth(2);
  h_res_rms->Draw("same,hist");
  h_stat_unc->SetLineColor(8); 
  h_stat_unc->SetLineWidth(2);
  h_stat_unc->Draw("same,hist");

  TLegend* leg2 = new TLegend(0.62,0.7,0.87,0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h_res_gaus, "Resolution (Gaus)", "l");
  leg2->AddEntry(h_res_rms, "Resolution (RMS)", "l");
  leg2->AddEntry(h_stat_unc, "Stat. Unc.", "l");
  leg2->Draw();

  if (doPL) mySmallText(0.2,0.8,1,0.04,"Particle level");
  else mySmallText(0.2,0.8,1,0.04,"Parton level");

  c2.SaveAs("UnfoldingPlots/resolution_"+channel+"_"+which+append+".pdf");


  // -------------------------------------------------------------------------------------
  // plot response matrices 
  // -------------------------------------------------------------------------------------

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

  TCanvas cr;

  //Draw fine-binned response matrix
  h_response->GetXaxis()->SetTitleOffset(0.8);
  h_response->GetYaxis()->SetTitleOffset(0.8);
  h_response->Draw("colz");
  if (drawAll) cr.SaveAs("UnfoldingPlots/unfold_fine_responseMatrix_full_"+which+"_"+channel+"_nom"+append+".pdf");

  // normalize so that for each bin of true top quark pt(eta), the bins in measured top pt(eta) add up to 100%
  for (int ir = 1; ir < h_response->GetNbinsY()+1; ir++){
    double rowsum = h_response->Integral(1,h_response->GetNbinsX()+1,ir,ir);
    for (int ic = 1; ic < h_response->GetNbinsX()+1; ic++){
      double normval = rowsum > 0.0 ? h_response->GetBinContent(ic,ir) / rowsum * 100.0 : 0.0;
      h_response->SetBinContent(ic,ir,normval);
    }
  }
  

  if (drawAll) {
    if (which == "pt"){
      h_response->SetAxisRange(401.0,1199.0,"X");
      h_response->SetAxisRange(401.0,1199.0,"Y");
      h_response->Draw("colz");
      cr.SaveAs("UnfoldingPlots/unfold_fine_responseMatrix_zoom_"+which+"_"+channel+"_nom"+append+".pdf");
    }
    else {
      h_response->Draw("colz");
      cr.SaveAs("UnfoldingPlots/unfold_fine_responseMatrix_"+which+"_"+channel+"_nom"+append+".pdf");
    }
  }

    
  // draw final binning response matrix
  if (drawAll) {
    //h_response_final->GetXaxis()->SetTitleOffset(0.8);
    //h_response_final->GetYaxis()->SetTitleOffset(1.2);
    h_response_final->Draw("colz");
    cr.SaveAs("UnfoldingPlots/unfold_responseMatrix_full_"+which+"_"+channel+"_nom"+append+".pdf");
  }

  // normalize so that for each bin of true top quark pt(eta), the bins in measured top pt(eta) add up to 100%
  for (int ir = 1; ir < nbins_final+1; ir++){
    double rowsum = h_response_final->Integral(1,nbins_final,ir,ir);
    for (int ic = 1; ic < nbins_final+1; ic++){
      double normval = rowsum > 0.0 ? h_response_final->GetBinContent(ic,ir) / rowsum * 100.0 : 0.0;
      h_response_final->SetBinContent(ic,ir,normval);
    }
  }


  if (which == "pt"){
    h_response_final->SetAxisRange(401.0,1199.0,"X");
    h_response_final->SetAxisRange(401.0,1199.0,"Y");
    h_response_final->Draw("colz");
    cr.SaveAs("UnfoldingPlots/unfold_responseMatrix_zoom_"+which+"_"+channel+"_nom"+append+".pdf");
  }
  else {
    h_response_final->Draw("colz");
    cr.SaveAs("UnfoldingPlots/unfold_responseMatrix_"+which+"_"+channel+"_nom"+append+".pdf");
  }

}
void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize, char *text) {
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
