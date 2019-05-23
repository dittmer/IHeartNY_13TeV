#include "makePlots.h"

#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "TExec.h"
#include "TPaletteAxis.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

#include <iostream>
#include <iomanip>


void setStyle() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(000000);
  
  gStyle->SetTitleFont(43);
  gStyle->SetTitleFont(43, "XYZ");
  gStyle->SetTitleSize(28, "XYZ");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.0, "Y");
  gStyle->SetLabelFont(43, "XYZ");
  gStyle->SetLabelSize(24, "XYZ");

  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(0.);

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetPalette(1);
  
  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  
}

void myLargeText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.06); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void myText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.043); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void myItalicText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.05); 
  l.SetTextFont(52); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void mySmallText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.042); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void drawCMS(Double_t x,Double_t y, bool prel) {

  float cmsTextSize = 0.08;
  float extraOverCmsTextSize = 0.68;
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

  TLatex ll;
  ll.SetTextSize(extraTextSize); 
  ll.SetTextFont(42); 
  ll.SetTextAngle(0);
  ll.SetNDC();
  ll.SetTextColor(1);
  ll.DrawLatex(0.71,y,"35.9 fb^{-1} (13 TeV)");

}


// -------------------------------------------------------------------------------------
// make pretty plots
// -------------------------------------------------------------------------------------

void makePlots(TString DIR, TString DIRqcd, TString channel, TString var, TString region, bool useQCDMC = false, bool unBlind = false, bool usePost = false, TString split = "") {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+region;
  
  // get histograms
  SummedHist* diboson; 
  SummedHist* zjets;
  SummedHist* wjets;
  SummedHist* singletop;
  SummedHist* ttbar;
  SummedHist* ttbar_nonSemiLep;
  SummedHist* data;

  SummedHist* diboson2; 
  SummedHist* zjets2;
  SummedHist* wjets2;
  SummedHist* singletop2;
  SummedHist* ttbar2;
  SummedHist* ttbar_nonSemiLep2;
  SummedHist* data2;

  if (channel=="comb") {
    diboson = getDiboson( DIR, var, region, "mu", false, "nom", usePost, split );
    zjets = getZJets( DIR, var, region, "mu", false, "nom", usePost, split );
    wjets  = getWJets( DIR, var, region, "mu", false, "nom", usePost, split );
    singletop = getSingleTop( DIR, var, region, "mu", false, "nom", usePost, split );
    ttbar = getTTbar( DIR, var, region, "mu", false, "nom", usePost, split );
    ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, region, "mu", false, "nom", usePost, split );
    data = getData( DIR, var, region, "mu", false, split);

    diboson2 = getDiboson( DIR, var, region, "el", false, "nom", usePost, split );
    zjets2 = getZJets( DIR, var, region, "el", false, "nom", usePost, split );
    wjets2  = getWJets( DIR, var, region, "el", false, "nom", usePost, split );
    singletop2 = getSingleTop( DIR, var, region, "el", false, "nom", usePost, split );
    ttbar2 = getTTbar( DIR, var, region, "el", false, "nom", usePost, split );
    ttbar_nonSemiLep2 = getTTbarNonSemiLep( DIR, var, region, "el", false, "nom", usePost, split );
    data2 = getData( DIR, var, region, "el", false, split);
  }
  else {
    diboson = getDiboson( DIR, var, region, channel, false, "nom", usePost, split );
    zjets = getZJets( DIR, var, region, channel, false, "nom", usePost, split );
    wjets  = getWJets( DIR, var, region, channel, false, "nom", usePost, split );
    singletop = getSingleTop( DIR, var, region, channel, false, "nom", usePost, split );
    ttbar = getTTbar( DIR, var, region, channel, false, "nom", usePost, split );
    ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, region, channel, false, "nom", usePost, split );
    data = getData( DIR, var, region, channel, false, split);
  }

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  
  TH1F* h_qcd;
  TH1F* h_qcd2;

  if (useQCDMC) {
    if (channel=="comb") {
      SummedHist* qcd = getQCDMC(DIR, var, region, "mu", false, "nom", usePost, split);
      h_qcd = (TH1F*) qcd->hist();
      SummedHist* qcd2 = getQCDMC(DIR, var, region, "el", false, "nom", usePost, split);
      h_qcd2 = (TH1F*) qcd2->hist();
    }
    else {
      SummedHist* qcd = getQCDMC(DIR, var, region, channel, false, "nom", usePost, split);
      h_qcd = (TH1F*) qcd->hist();
    }
  }
  else {
    if (channel=="comb") {
      h_qcd = (TH1F*) getQCDData( DIR, DIRqcd, var, region, "mu", "nom", usePost, split);
      h_qcd2 = (TH1F*) getQCDData( DIR, DIRqcd, var, region, "el", "nom", usePost, split);
    }
    else {
      h_qcd = (TH1F*) getQCDData( DIR, DIRqcd, var, region, channel, "nom", usePost, split); // Currently taking QCD from data sideband with normalization from MC in signal region
    }
  }
  TH1F* h_diboson = (TH1F*) diboson->hist();
  TH1F* h_zjets = (TH1F*) zjets->hist();
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  TH1F* h_singletop = (TH1F*) singletop->hist();
  TH1F* h_data = (TH1F*) data->hist();

  TH1F* h_diboson2;
  TH1F* h_zjets2;
  TH1F* h_wjets2;
  TH1F* h_ttbar2;
  TH1F* h_ttbar_nonSemiLep2;
  TH1F* h_singletop2;
  TH1F* h_data2;

  if (channel=="comb") {
    h_diboson2 = (TH1F*) diboson2->hist();
    h_zjets2 = (TH1F*) zjets2->hist();
    h_wjets2 = (TH1F*) wjets2->hist();
    h_ttbar2 = (TH1F*) ttbar2->hist();
    h_ttbar_nonSemiLep2 = (TH1F*) ttbar_nonSemiLep2->hist();
    h_singletop2 = (TH1F*) singletop2->hist();
    h_data2 = (TH1F*) data2->hist();
  }


  // -----------------------------------------------------------------------------------------------------
  // Apply post-fit normalizations if desired
  // The scale factors are determined from the post-fit nuisance parameters, not the post-fit event yields

  if (usePost && !(var.Contains("Raw"))){
    if (channel == "comb") {
      h_qcd->Scale(1.02);
      h_qcd2->Scale(1.32);
    }
    else if (channel == "mu") h_qcd->Scale(1.02);
    else h_qcd->Scale(1.32);                 
    h_diboson->Scale(1.01);
    h_zjets->Scale(0.94);
    h_singletop->Scale(1.16);
    h_ttbar->Scale(0.81);
    h_ttbar_nonSemiLep->Scale(0.81);
    if (channel=="comb") {
      h_diboson2->Scale(1.01);
      h_zjets2->Scale(0.94);
      h_singletop2->Scale(1.16);
      h_ttbar2->Scale(0.81);
      h_ttbar_nonSemiLep2->Scale(0.81);
    }

    if (var == "ak8jetPt" && region == "1t1b") {
      if (channel=="comb") {
	SummedHist* wjetsL = getWJets( DIR, var, region, "mu", false, "nom", usePost, "l");
	TH1F* h_wjetsL = (TH1F*) wjetsL->hist();
	h_wjets->Add(h_wjetsL,-1.0);
	h_wjets->Scale(0.98);
	h_wjetsL->Scale(0.77);
	h_wjets->Add(h_wjetsL);

	SummedHist* wjetsL2 = getWJets( DIR, var, region, "el", false, "nom", usePost, "l");
	TH1F* h_wjetsL2 = (TH1F*) wjetsL2->hist();
	h_wjets2->Add(h_wjetsL2,-1.0);
	h_wjets2->Scale(0.98);
	h_wjetsL2->Scale(0.77);
	h_wjets2->Add(h_wjetsL2);
      }
      else {
	SummedHist* wjetsL = getWJets( DIR, var, region, channel, false, "nom", usePost, "l");
	TH1F* h_wjetsL = (TH1F*) wjetsL->hist();
	h_wjets->Add(h_wjetsL,-1.0);
	h_wjets->Scale(0.98);
	h_wjetsL->Scale(0.77);
	h_wjets->Add(h_wjetsL);
      }
    }
    else {
      h_wjets->Scale(0.84); //WJets is 1/3 HF, 2/3 LF --> combine above scale factors
      if (channel=="comb") h_wjets2->Scale(0.84);
    }
  }


  // -------------------------------------------------------------------------------------
  // for combined electron+muon channel, merge the two histograms here

  if (channel=="comb") {
    h_data->Add(h_data2);
    h_qcd->Add(h_qcd2);
    h_diboson->Add(h_diboson2);
    h_zjets->Add(h_zjets2);
    h_wjets->Add(h_wjets2);
    h_singletop->Add(h_singletop2);
    h_ttbar_nonSemiLep->Add(h_ttbar_nonSemiLep2);
    h_ttbar->Add(h_ttbar2);

    h_data2->Delete();
    h_qcd2->Delete();
    h_diboson2->Delete();
    h_zjets2->Delete();
    h_wjets2->Delete();
    h_singletop2->Delete();
    h_ttbar_nonSemiLep2->Delete();
    h_ttbar2->Delete();
  }


  // -------------------------------------------------------------------------------------
  // various hist plotting edits

  int rebinby = 1;
  if ((hist.Contains("ht") && !hist.Contains("htLep")) || hist.Contains("ak4jetPt")) rebinby = 3;
  else if (hist.Contains("ak8jetPt")) rebinby = 2; //4;
  else if (hist.Contains("ak8jetTau") || hist.Contains("ak4jetEta") || hist.Contains("Phi") || hist.Contains("dR")) rebinby = 5;
  else if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))) rebinby = 2;
  
  if (h_qcd) h_qcd->Rebin(rebinby);
  if (h_diboson) h_diboson->Rebin(rebinby);
  if (h_zjets) h_zjets->Rebin(rebinby);
  if (h_wjets) h_wjets->Rebin(rebinby);
  if (h_singletop) h_singletop->Rebin(rebinby);
  if (h_ttbar_nonSemiLep) h_ttbar_nonSemiLep->Rebin(rebinby);
  if (h_ttbar) h_ttbar->Rebin(rebinby);
  h_data->Rebin(rebinby);

  h_data->UseCurrentStyle();
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);

  // create stack & summed histogram for ratio plot
  THStack* h_stack = new THStack();    
  if (h_qcd) h_stack->Add(h_qcd);
  if (h_diboson) h_stack->Add(h_diboson);
  if (h_zjets) h_stack->Add(h_zjets);
  if (h_wjets) h_stack->Add(h_wjets);
  if (h_singletop) h_stack->Add(h_singletop);
  if (h_ttbar_nonSemiLep) h_stack->Add(h_ttbar_nonSemiLep);
  if (h_ttbar) h_stack->Add(h_ttbar);

  TH1F* h_totalbkg = (TH1F*) h_ttbar->Clone("totalbkg_"+hist);
  if (h_ttbar_nonSemiLep) h_totalbkg->Add(h_ttbar_nonSemiLep);
  if (h_diboson) h_totalbkg->Add(h_diboson);
  if (h_zjets) h_totalbkg->Add(h_zjets);
  if (h_wjets) h_totalbkg->Add(h_wjets);
  if (h_singletop) h_totalbkg->Add(h_singletop);
  if (h_qcd) h_totalbkg->Add(h_qcd);

  // -------------------------------------------------------------------------------------
  // poisson errors
  h_data->SetBinErrorOption(TH1::kPoisson);
  for (int ibindata=0; ibindata<(int)h_data->GetNbinsX(); ibindata++) {
    double err_low = h_data->GetBinErrorLow(ibindata);
    double err_up = h_data->GetBinErrorUp(ibindata);
  }
  // -------------------------------------------------------------------------------------
  
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {

    float stat_tt, stat_tt_non, stat_st, stat_wj, stat_db, stat_zj, stat_qcd = 0.0;
    stat_tt = h_ttbar->GetBinError(ib+1);
    if (h_ttbar_nonSemiLep) stat_tt_non = h_ttbar_nonSemiLep->GetBinError(ib+1);
    if (h_singletop) stat_st = h_singletop->GetBinError(ib+1);
    if (h_diboson) stat_db = h_diboson->GetBinError(ib+1);
    if (h_zjets) stat_zj = h_zjets->GetBinError(ib+1);
    if (h_wjets) stat_wj = h_wjets->GetBinError(ib+1);
    if (h_qcd) stat_qcd = h_qcd->GetBinError(ib+1);

    float binbkg = sqrt( stat_tt*stat_tt + // Stat error only for now
			 stat_tt_non*stat_tt_non + 
			 stat_st*stat_st + 
			 stat_wj*stat_wj +
			 stat_zj*stat_zj +
			 stat_db*stat_db +
			 stat_qcd*stat_qcd 
			 );
    
    h_totalbkg->SetBinError(ib+1, binbkg);
  }

  if (hist.Contains("ak8jetPt")) h_data->GetXaxis()->SetRangeUser(400.,1200.);
  if (hist.Contains("met") && (channel == "mu" || channel == "comb")) h_data->GetXaxis()->SetRangeUser(30.,250.);
  if (hist.Contains("met") && channel == "el") h_data->GetXaxis()->SetRangeUser(50.,250.);
  
  TH1F* h_ratio;
  TH1F* h_ratio2;
  h_ratio = (TH1F*) h_data->Clone("ratio_"+hist);  // Data / MC
  h_ratio->Divide(h_totalbkg);
  h_ratio2 = (TH1F*) h_totalbkg->Clone("ratio2_"+hist); // Uncertainty on Data / MC
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {
    h_ratio2->SetBinContent(ib+1, 1.0);
    float tmperr = h_totalbkg->GetBinError(ib+1);
    float tmpcount = h_totalbkg->GetBinContent(ib+1);
    
    if (tmpcount==0) continue;
    h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
  }

  float mymax = max(h_data->GetMaximum(),h_totalbkg->GetMaximum());
  if (hist.Contains("Pt")) h_data->SetAxisRange(0,mymax*1.1,"Y");
  else h_data->SetAxisRange(0,mymax*1.4,"Y");

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1;
  TPad* p2;
  if (region == "Pre" || unBlind){
    p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.28,1,1);
    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.02);
    p1->SetNumber(1);
    p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.28);
    p2->SetNumber(2);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.35);
    
    p1->Draw();
    p2->Draw();
    p1->cd();
    //if (hist.Contains("Pt")) p1->SetLogy();
    
    h_ratio2->SetMarkerSize(0);
    h_ratio2->SetLineColor(0);
    h_ratio2->SetFillColor(15);
    h_ratio2->SetFillStyle(1001);
    
    h_ratio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    
    h_data->GetXaxis()->SetLabelSize(0);
    h_data->GetYaxis()->SetLabelSize(26);
    h_data->GetYaxis()->SetTitleSize(32);
    h_data->GetYaxis()->SetTitleOffset(1.4);
    h_data->GetXaxis()->SetTitle("");
    
    if (hist.Contains("ak8jetY")) h_data->GetYaxis()->SetTitle("Events / 0.2");
    else if (hist.Contains("ak8jetPt")) h_data->GetYaxis()->SetTitle("Events / 20 GeV");

    h_data->Draw("LE0P");
  }

  h_totalbkg->UseCurrentStyle();
  h_totalbkg->SetFillColor(0);
  h_totalbkg->SetLineWidth(1);
  h_totalbkg->SetLineColor(1);
  
  if (region == "Pre" || unBlind) h_totalbkg->Draw("hist,same");
  else h_totalbkg->Draw("hist");
  
  h_stack->Draw("hist,same");
  if (region == "Pre" || unBlind) h_data->Draw("LE0P,same");
  h_data->Draw("axis,same");

  float xmin = 0.73;
  float ymin = 0.48;

  float xwidth = 0.20;
  float ywidth = 0.40;

  if (region != "Pre" && !unBlind){
    ymin = 0.57;
    xwidth = 0.18;
    ywidth = 0.34;
  }

  //Legend top left
  if (var.Contains("ak8jetTau") || (var.Contains("CSV") && (region == "0t1b" || region == "1t1b")) || ((var == "ak8jetMass" || var == "ak8jetSDm01" || var == "ak8jetSDmass" || var == "ak8jetSubjetMaxCSV") && (region == "1t0b" || region == "1t1b"))) xmin = 0.16;

  //Legend bottom center
  if (var.Contains("Phi")) {
    xmin = 0.40;
    ymin = 0.12;
    if (region != "Pre" && !unBlind) ymin = 0.25;
    if ((region == "0t1b" || region == "1t0b" || region == "1t1b" ) && var == "ak4jetPhi") ymin = 0.45;
  }

  //Legend center center
  if (var.Contains("Phi") && channel == "el" && (region == "1t0b" || region == "1t1b")){
    ymin = 0.35;
    if (var == "ak4jetPhi") ymin = 0.31;
  }

  //Legend top center
  if (var.Contains("CSV") && (region == "0t" || region == "Pre")) xmin = 0.40;
  
  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  if (region == "Pre" || unBlind) leg->SetTextSize(0.045);
  else leg->SetTextSize(0.037);
  if (region == "Pre" || unBlind) leg->AddEntry(h_data, "Data", "pe");
  leg->AddEntry(h_ttbar, "t#bar{t} signal", "f");
  leg->AddEntry(h_ttbar_nonSemiLep, "t#bar{t} other", "f");
  leg->AddEntry(h_singletop, "Single t", "f");
  leg->AddEntry(h_wjets, "W+jets", "f");
  leg->AddEntry(h_zjets, "Z+jets", "f");
  leg->AddEntry(h_diboson, "Diboson", "f");
  leg->AddEntry(h_qcd, "Multijet" , "f");
  leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
  leg->Draw();

  //myText(0.10,0.94,1,"#intLdt = 35.9 fb^{-1}");
  //myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");
  drawCMS(0.15,0.92, false);
  if (channel=="comb") myLargeText(0.6,0.82,1,"l+jets");
  

  // plot ratio part
  if (region == "Pre" || unBlind){
    p2->cd();
    p2->SetGridy();
    h_ratio->UseCurrentStyle();
    h_ratio->SetMarkerStyle(8);
    h_ratio->SetMarkerSize(1);
    h_ratio->Draw("le0p");
    h_ratio2->Draw("same,e2");
    h_ratio->Draw("le0p,same");
    h_ratio->SetMaximum(1.75);
    h_ratio->SetMinimum(0.25);
    h_ratio->GetYaxis()->SetNdivisions(4,4,0,true);
    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetXaxis()->SetLabelSize(28);
    h_ratio->GetYaxis()->SetLabelSize(26);
    h_ratio->GetXaxis()->SetTitleOffset(3.5);
    h_ratio->GetYaxis()->SetTitleOffset(1.2);
    h_ratio->GetXaxis()->SetTitleSize(32);
    h_ratio->GetYaxis()->SetTitleSize(32);
  }

  // save output
  TString append = "";
  if (usePost) append = "_post";
  TString outname = "Plots/"+channel+"_"+hist+append+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_stack->Delete();
  leg->Delete();
}

void troubleshootQCD(TString DIR, TString var, TString channel) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+"Pre";
  
  // get histograms
  SummedHist* wjets = getWJets( DIR, var, "Pre", channel, true, "nom", false );
  SummedHist* singletop = getSingleTop( DIR, var, "Pre", channel, true, "nom", false );
  SummedHist* ttbar = getTTbar( DIR, var, "Pre", channel, true, "nom", false );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, "Pre", channel, true, "nom", false );
  SummedHist* qcd = getQCDMC(DIR, var, "Pre", channel, true, "nom", false); 
  SummedHist* data = getData( DIR, var, "Pre", channel, true);

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  TH1F* h_singletop = (TH1F*) singletop->hist();
  TH1F* h_qcd = (TH1F*) qcd->hist();
  TH1F* h_data = (TH1F*) data->hist();

  // -------------------------------------------------------------------------------------
  // various hist plotting edits
  
  if (hist.Contains("ak8jetTau") || hist.Contains("lepAbsEta")){
    if (h_qcd) h_qcd->Rebin(5);
    if (h_wjets) h_wjets->Rebin(5);
    if (h_singletop) h_singletop->Rebin(5);
    if (h_ttbar_nonSemiLep) h_ttbar_nonSemiLep->Rebin(5);
    if (h_ttbar) h_ttbar->Rebin(5);
    h_data->Rebin(5);
  }
  
  else if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))){
    if (h_qcd) h_qcd->Rebin(2);
    if (h_wjets) h_wjets->Rebin(2);
    if (h_singletop) h_singletop->Rebin(2);
    if (h_ttbar_nonSemiLep) h_ttbar_nonSemiLep->Rebin(2);
    if (h_ttbar) h_ttbar->Rebin(2);
    h_data->Rebin(2);
  }

  h_data->UseCurrentStyle();
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);

  // create stack & summed histogram for ratio plot
  THStack* h_stack = new THStack();    
  if (h_qcd) h_stack->Add(h_qcd);
  if (h_wjets) h_stack->Add(h_wjets);
  if (h_singletop) h_stack->Add(h_singletop);
  if (h_ttbar_nonSemiLep) h_stack->Add(h_ttbar_nonSemiLep);
  if (h_ttbar) h_stack->Add(h_ttbar);

  TH1F* h_totalbkg = (TH1F*) h_ttbar->Clone("totalbkg_"+hist);
  if (h_ttbar_nonSemiLep) h_totalbkg->Add(h_ttbar_nonSemiLep);
  if (h_wjets) h_totalbkg->Add(h_wjets);
  if (h_singletop) h_totalbkg->Add(h_singletop);
  if (h_qcd) h_totalbkg->Add(h_qcd);

  // -------------------------------------------------------------------------------------
  // poisson errors
  h_data->SetBinErrorOption(TH1::kPoisson);
  for (int ibindata=0; ibindata<(int)h_data->GetNbinsX(); ibindata++) {
    double err_low = h_data->GetBinErrorLow(ibindata);
    double err_up = h_data->GetBinErrorUp(ibindata);
  }
  // -------------------------------------------------------------------------------------
  
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {

    float stat_tt, stat_tt_non, stat_st, stat_wj, stat_qcd = 0.0;
    stat_tt = h_ttbar->GetBinError(ib+1);
    if (h_ttbar_nonSemiLep) stat_tt_non = h_ttbar_nonSemiLep->GetBinError(ib+1);
    if (h_singletop) stat_st = h_singletop->GetBinError(ib+1);
    if (h_wjets) stat_wj = h_wjets->GetBinError(ib+1);
    if (h_qcd) stat_qcd = h_qcd->GetBinError(ib+1);

    float binbkg = sqrt( stat_tt*stat_tt + // Stat error only for now
			 stat_tt_non*stat_tt_non + 
			 stat_st*stat_st + 
			 stat_wj*stat_wj + 
			 stat_qcd*stat_qcd 
			 );
    
    h_totalbkg->SetBinError(ib+1, binbkg);
  }

  if (hist.Contains("ak8jetPt")) h_data->GetXaxis()->SetRangeUser(400.,1200.);
  
  TH1F* h_ratio;
  TH1F* h_ratio2;
  h_ratio = (TH1F*) h_data->Clone("ratio_"+hist);  // Data / MC
  h_ratio->Divide(h_totalbkg);
  h_ratio2 = (TH1F*) h_totalbkg->Clone("ratio2_"+hist); // Uncertainty on Data / MC
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {
    h_ratio2->SetBinContent(ib+1, 1.0);
    float tmperr = h_totalbkg->GetBinError(ib+1);
    float tmpcount = h_totalbkg->GetBinContent(ib+1);
    
    if (tmpcount==0) continue;
    h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
  }

  float mymax = max(h_data->GetMaximum(),h_totalbkg->GetMaximum());
  h_data->SetAxisRange(0,mymax*1.05,"Y");

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
  p1->SetTopMargin(0.08);
  p1->SetBottomMargin(0.05);
  p1->SetNumber(1);
  TPad* p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
  p2->SetNumber(2);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.35);
    
  p1->Draw();
  p2->Draw();
  p1->cd();
    
  h_ratio2->SetMarkerSize(0);
  h_ratio2->SetLineColor(0);
  h_ratio2->SetFillColor(15);
  h_ratio2->SetFillStyle(1001);
  
  h_ratio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());

  h_data->GetYaxis()->SetLabelSize(26);
  h_data->GetYaxis()->SetTitleSize(32);
  h_data->GetYaxis()->SetTitleOffset(1.4);
  h_data->GetXaxis()->SetTitle("");
    
  h_data->Draw("LE0P");

  h_totalbkg->UseCurrentStyle();
  h_totalbkg->SetFillColor(0);
  h_totalbkg->SetLineWidth(1);
  h_totalbkg->SetLineColor(1);
  
  h_totalbkg->Draw("hist,same");
  h_stack->Draw("hist,same");
  h_data->Draw("LE0P,same");

  float xmin = 0.71;
  float ymin = 0.52;
  float xwidth = 0.20;
  float ywidth = 0.38;

  //Legend top left
  if (var == "ak8jetTau32" || var == "ak8jetTau21") xmin = 0.16;

  //Legend bottom center
  if (var.Contains("Phi")) {
    xmin = 0.40;
    ymin = 0.12;
  }

  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_data, "Data", "pe");
  leg->AddEntry(h_ttbar, "t#bar{t} signal", "f");
  leg->AddEntry(h_ttbar_nonSemiLep, "t#bar{t} other", "f");
  leg->AddEntry(h_singletop, "Single t", "f");
  leg->AddEntry(h_wjets, "W+jets", "f");
  leg->AddEntry(h_qcd, "Multijet" , "f");
  leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
  leg->Draw();

  myText(0.10,0.94,1,"#intLdt = 35.9 fb^{-1}");
  myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");

  // plot ratio part
  p2->cd();
  p2->SetGridy();
  h_ratio->UseCurrentStyle();
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(1);
  h_ratio->Draw("le0p");
  h_ratio2->Draw("same,e2");
  h_ratio->Draw("le0p,same");
  h_ratio->SetMaximum(1.8);
  h_ratio->SetMinimum(0.2);
  h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
  h_ratio->GetYaxis()->SetTitle("Data / MC");
  h_ratio->GetXaxis()->SetLabelSize(26);
  h_ratio->GetYaxis()->SetLabelSize(26);
  h_ratio->GetXaxis()->SetTitleOffset(2.8);
  h_ratio->GetYaxis()->SetTitleOffset(1.4);
  h_ratio->GetXaxis()->SetTitleSize(32);
  h_ratio->GetYaxis()->SetTitleSize(32);

  // save output
  TString outname = "Plots/Sideband_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_stack->Delete();
  leg->Delete();
}

void compareShapes(TString DIR, TString DIRqcd, TString channel, TString var, TString region) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+region;
  
  // get histograms
  SummedHist* wjets = getWJets( DIR, var, region, channel, false, "nom", false );
  SummedHist* ttbar = getTTbar( DIR, var, region, channel, false, "nom", false );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, region, channel, false, "nom", false );
  //SummedHist* singletop = getSingleTop( DIR, var, region, channel, false, "nom", false );

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  
  TH1F* h_qcd = (TH1F*) getQCDData( DIR, DIRqcd, var, region, channel, "nom", false); // Currently taking QCD from data sideband with no triangular cut,
                                                                                      // with normalization from MC in signal region
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  h_ttbar->Add(h_ttbar_nonSemiLep);
  //TH1F* h_singletop = (TH1F*) singletop->hist();

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  
  if (h_qcd){
    h_qcd->SetLineColor(h_qcd->GetFillColor());
    h_qcd->SetFillColor(0);
    h_qcd->Scale(1./h_qcd->Integral());
  }
  h_wjets->SetLineColor(h_wjets->GetFillColor());
  h_wjets->SetFillColor(0);
  h_wjets->Scale(1./h_wjets->Integral());
  h_ttbar->SetLineColor(h_ttbar->GetFillColor());
  h_ttbar->SetFillColor(0);
  h_ttbar->Scale(1./h_ttbar->Integral());
  //h_singletop->SetLineColor(h_singletop->GetFillColor());
  //h_singletop->SetFillColor(0);
  //h_singletop->Scale(1./h_singletop->Integral());

  if (hist.Contains("ak4jetEta") || hist.Contains("ak8jetSDmass") || hist.Contains("ak8jetTau")){
    if (h_qcd) h_qcd->Rebin(5);
    h_wjets->Rebin(5);
    h_ttbar->Rebin(5);
    //h_singletop->Rebin(5);
  }

  else if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))){
    if (h_qcd) h_qcd->Rebin(2);
    h_wjets->Rebin(2);
    h_ttbar->Rebin(2);
    //h_singletop->Rebin(2);
  }


  h_ttbar->SetMaximum(h_ttbar->GetMaximum()*1.5);
  
  h_ttbar->Draw("hist");
  h_wjets->Draw("hist,same");
  if (h_qcd) h_qcd->Draw("hist,same");
  //h_singletop->Draw("hist,same");

  float xmin = 0.70;
  float ymin = 0.70;

  if (hist.Contains("Tau32") || hist.Contains("SDmass")) xmin = 0.2;

  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+0.2,ymin+0.2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.042);
  leg->AddEntry(h_ttbar, "t#bar{t}", "l");
  leg->AddEntry(h_wjets, "W+jets", "l");
  if (h_qcd) leg->AddEntry(h_qcd, "Multijet" , "l");
  //leg->AddEntry(h_singletop, "Single t", "l");
  leg->Draw();

  // save output
  TString outname = "Plots/compShapes_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  leg->Delete();
  if (h_qcd) h_qcd->Delete();
  h_ttbar->Delete();
  h_ttbar_nonSemiLep->Delete();
  h_wjets->Delete();
  //h_singletop->Delete();
  delete wjets;
  delete ttbar;
  delete ttbar_nonSemiLep;
  //delete singletop;
}


void makeCombineInputs(TString DIR, TString DIRqcd, TString whichQCD) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  const int nchannels = 2;
  TString channels[nchannels] = {"mu","el"};
  const int nhist = 13;
  TString histnames[nhist] = {"lepEta","lepEta","ak4jetEta","ak4jetEta","ak8jetTau21","ak8jetTau32","ak8jetSDmass","lepAbsEta","lepAbsEta","ak4jetAbsEta","ak4jetAbsEta","ak4jetCSV","counts"};
  TString regions[nhist] = {"0t","1t0b","0t","1t0b","0t","1t0b","1t1b","0t","1t0b","0t","1t0b","0t",""};
  int rebinby[nhist] = {5,5,5,5,5,5,2,5,5,5,5,5,1};
  float lowbounds[nhist] = {-2.5,-2.5,-2.5,-2.5,0.1,0.2,100.0,0.0,0.0,0.0,0.0,0.0,-0.5};
  float highbounds[nhist] = {2.5,2.5,2.5,2.5,1.0,0.84,220.0,2.5,2.5,2.5,2.5,1.0,2.5};
  const int nbins = 3;
  TString binnames[nbins] = {"","barrel","endcap"};
  const int nsys = 11;
  TString sysnames[nsys] = {"nom","lepUp","lepDown","JECUp","JECDown","JERUp","JERDown","BTagUp","BTagDown","TopTagUp","TopTagDown"};
  
  TH1F* h_qcd[nchannels][nhist][nbins][nsys];
  TH1F* h_diboson[nchannels][nhist][nbins][nsys];
  TH1F* h_zjets[nchannels][nhist][nbins][nsys];
  TH1F* h_wjets[nchannels][nhist][nbins][nsys];
  TH1F* h_wjetsL[nchannels][nhist][nbins][nsys];
  TH1F* h_wjetsHF[nchannels][nhist][nbins][nsys];
  TH1F* h_ttbar[nchannels][nhist][nbins][nsys];
  TH1F* h_singletop[nchannels][nhist][nbins][nsys];
  TH1F* h_singletop_tW[nchannels][nhist][nbins][nsys];
  TH1F* h_singletop_other[nchannels][nhist][nbins][nsys];
  TH1F* h_data[nchannels][nhist][nbins];
  
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){
      for (int ib = 0; ib < nbins; ib++){
	for (int kk = 0; kk < nsys; kk++){
	  TString append = (sysnames[kk] == "nom") ? "" : "_"+(sysnames[kk].Contains("TopTag") ? binnames[ib]+sysnames[kk] : sysnames[kk]);
	  append.ReplaceAll("lep",channels[ii]+"SF");
	  
	  // get histograms
	  SummedHist* diboson = getDiboson( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib]);
	  SummedHist* zjets = getZJets( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib]);
	  SummedHist* wjets = getWJets( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib]);
	  SummedHist* wjetsL  = getWJets( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, "l"+(binnames[ib] == "" ? "" : "_"+binnames[ib]));
	  SummedHist* singletop = getSingleTop( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib] );
	  SummedHist* singletop_tW = getSingleTop_tW( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib] );
	  SummedHist* ttbar = getTTbar( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib] );
	  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib] );

	  // -------------------------------------------------------------------------------------
	  // get the TH1F versions
	  
	  TH1F* tmp_diboson = (TH1F*) diboson->hist()->Clone("Diboson"+append);
	  TH1F* tmp_zjets = (TH1F*) zjets->hist()->Clone("ZJets"+append);
	  TH1F* tmp_wjets = (TH1F*) wjets->hist()->Clone("WJets"+append);
	  TH1F* tmp_wjetsL = (TH1F*) wjetsL->hist()->Clone("WJetsL"+append);
	  TH1F* tmp_wjetsHF = (TH1F*) tmp_wjets->Clone("WJetsHF"+append);
	  tmp_wjetsHF->Add(tmp_wjetsL,-1.0);
	  TH1F* tmp_ttbar = (TH1F*) ttbar->hist()->Clone("TTbar"+append);
	  tmp_ttbar->Add((TH1F*) ttbar_nonSemiLep->hist());
	  TH1F* tmp_singletop = (TH1F*) singletop->hist()->Clone("SingleTop"+append);
	  TH1F* tmp_singletop_tW = (TH1F*) singletop_tW->hist()->Clone("ST_tW"+append);
	  TH1F* tmp_singletop_other = (TH1F*) tmp_singletop->Clone("ST_other"+append);
	  tmp_singletop_other->Add(tmp_singletop_tW,-1.0);
	  
	  // Get QCD
	  TH1F* tmp_qcd;
	  if (whichQCD == "MC"){
	    SummedHist* qcd = getQCDMC( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk], false, binnames[ib] );
	    tmp_qcd = (TH1F*) qcd->hist()->Clone("QCD"+append);
	  }
	  if (whichQCD == "data"){
	    tmp_qcd = (TH1F*) getQCDData( DIR, DIRqcd, histnames[jj], regions[jj], channels[ii], sysnames[kk], false, binnames[ib])->Clone("QCD"+append);
	  }	    
	  
	  tmp_qcd->Rebin(rebinby[jj]);
	  tmp_diboson->Rebin(rebinby[jj]);
	  tmp_zjets->Rebin(rebinby[jj]);
	  tmp_wjets->Rebin(rebinby[jj]);
	  tmp_wjetsL->Rebin(rebinby[jj]);
	  tmp_wjetsHF->Rebin(rebinby[jj]);
	  tmp_ttbar->Rebin(rebinby[jj]);
	  tmp_singletop->Rebin(rebinby[jj]);
	  tmp_singletop_tW->Rebin(rebinby[jj]);
	  tmp_singletop_other->Rebin(rebinby[jj]);
	  
	  h_qcd[ii][jj][ib][kk] = adjustRange(tmp_qcd,lowbounds[jj],highbounds[jj]);
	  h_diboson[ii][jj][ib][kk] = adjustRange(tmp_diboson,lowbounds[jj],highbounds[jj]);
	  h_zjets[ii][jj][ib][kk] = adjustRange(tmp_zjets,lowbounds[jj],highbounds[jj]);
	  h_wjets[ii][jj][ib][kk] = adjustRange(tmp_wjets,lowbounds[jj],highbounds[jj]);
	  h_wjetsL[ii][jj][ib][kk] = adjustRange(tmp_wjetsL,lowbounds[jj],highbounds[jj]);
	  h_wjetsHF[ii][jj][ib][kk] = adjustRange(tmp_wjetsHF,lowbounds[jj],highbounds[jj]);
	  h_singletop[ii][jj][ib][kk] = adjustRange(tmp_singletop,lowbounds[jj],highbounds[jj]);
	  h_singletop_tW[ii][jj][ib][kk] = adjustRange(tmp_singletop_tW,lowbounds[jj],highbounds[jj]);
	  h_singletop_other[ii][jj][ib][kk] = adjustRange(tmp_singletop_other,lowbounds[jj],highbounds[jj]);
	  h_ttbar[ii][jj][ib][kk] = adjustRange(tmp_ttbar,lowbounds[jj],highbounds[jj]);
	  
	  tmp_qcd->Delete();
	  tmp_diboson->Delete();
	  tmp_zjets->Delete();
	  tmp_wjets->Delete();
	  tmp_wjetsL->Delete();
	  tmp_wjetsHF->Delete();
	  tmp_singletop->Delete();
	  tmp_singletop_tW->Delete();
	  tmp_singletop_other->Delete();
	  tmp_ttbar->Delete();
	  
	  h_qcd[ii][jj][ib][kk]->SetName(((TString)h_qcd[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_diboson[ii][jj][ib][kk]->SetName(((TString)h_diboson[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_zjets[ii][jj][ib][kk]->SetName(((TString)h_zjets[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_wjets[ii][jj][ib][kk]->SetName(((TString)h_wjets[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_wjetsL[ii][jj][ib][kk]->SetName(((TString)h_wjetsL[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_wjetsHF[ii][jj][ib][kk]->SetName(((TString)h_wjetsHF[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	  h_singletop_other[ii][jj][ib][kk]->SetName(((TString)h_singletop_other[ii][jj][ib][kk]->GetName()).ReplaceAll("TopTag","TopMisTag"));
	}

	// Get data
	SummedHist* data = getData( DIR, histnames[jj], regions[jj], channels[ii], false, binnames[ib]);
	TH1F* tmp_data = (TH1F*) data->hist()->Clone("data_obs");
	tmp_data->Rebin(rebinby[jj]);
	h_data[ii][jj][ib] = adjustRange(tmp_data,lowbounds[jj],highbounds[jj]);
	tmp_data->Delete();
      }
    }
  }

  // Write input file
  TFile* fout = new TFile("combineInputs_"+whichQCD+".root","recreate");
  TDirectory* topdir = fout->GetDirectory("");
  TDirectory* dirs[nchannels][nhist][nbins];
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){
      for (int ib = 0; ib < nbins; ib++){
	dirs[ii][jj][ib] = fout->mkdir(histnames[jj]+regions[jj]+(binnames[ib] != "" ? "_"+binnames[ib] : "")+"_"+channels[ii]);
	dirs[ii][jj][ib]->cd();
	for (int kk = 0; kk < nsys; kk++){
	  h_qcd[ii][jj][ib][kk]->Write();
	  h_diboson[ii][jj][ib][kk]->Write();
	  h_zjets[ii][jj][ib][kk]->Write();
	  h_wjets[ii][jj][ib][kk]->Write();
	  h_wjetsL[ii][jj][ib][kk]->Write();
	  h_wjetsHF[ii][jj][ib][kk]->Write();
	  h_singletop[ii][jj][ib][kk]->Write();
	  h_singletop_tW[ii][jj][ib][kk]->Write();
	  h_singletop_other[ii][jj][ib][kk]->Write();
	  h_ttbar[ii][jj][ib][kk]->Write();
	}
	h_data[ii][jj][ib]->Write();
	topdir->cd();
      }
    }
  }
  fout->Close();

  // Make comparison plots
  int colors[nsys-1] = {2,2,3,3,4,4,5,5,6,6};
  int styles[nsys-1] = {1,2,1,2,1,2,1,2,1,2};
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){

      TCanvas* c = new TCanvas("c_"+histnames[jj]+regions[jj]+"_"+channels[ii],"c_"+histnames[jj]+regions[jj]+"_"+channels[ii],900,800);

      h_ttbar[ii][jj][0][0]->SetFillColor(0);
      h_ttbar[ii][jj][0][0]->SetMaximum(1.5*h_ttbar[ii][jj][0][0]->GetMaximum());
      h_ttbar[ii][jj][0][0]->Draw("hist");
    
      TLegend* leg;
      float xmin = 0.66;
      if ((regions[jj] == "1t1b" || regions[jj] == "1t0b") && !(histnames[jj].Contains("lepAbsEta"))) xmin = 0.2;
      leg = new TLegend(xmin,0.65,xmin+0.2,0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      leg->AddEntry(h_ttbar[ii][jj][0][0], "TTbar nominal", "l");

      for (int kk = 1; kk < nsys; kk++){
	h_ttbar[ii][jj][0][kk]->SetFillColor(0);
	h_ttbar[ii][jj][0][kk]->SetLineColor(colors[kk-1]);
	h_ttbar[ii][jj][0][kk]->SetLineStyle(styles[kk-1]);
	h_ttbar[ii][jj][0][kk]->Draw("hist,same");
	if (kk%2 == 1) leg->AddEntry(h_ttbar[ii][jj][0][kk],sysnames[kk],"l");
      }

      leg->Draw();

      TString outname = "Plots/sysVar_"+channels[ii]+"_"+histnames[jj]+regions[jj]+".pdf";
      c->SaveAs(outname);
    }
  }
  
  // Cleanup
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj ++){
      for (int ib = 0; ib < nbins; ib++){
	for (int kk = 0; kk < nsys; kk++){
	  h_qcd[ii][jj][ib][kk]->Delete();
	  h_diboson[ii][jj][ib][kk]->Delete();
	  h_zjets[ii][jj][ib][kk]->Delete();
	  h_wjets[ii][jj][ib][kk]->Delete();
	  h_wjetsL[ii][jj][ib][kk]->Delete();
	  h_wjetsHF[ii][jj][ib][kk]->Delete();
	  h_singletop[ii][jj][ib][kk]->Delete();
	  h_singletop_tW[ii][jj][ib][kk]->Delete();
	  h_singletop_other[ii][jj][ib][kk]->Delete();
	  h_ttbar[ii][jj][ib][kk]->Delete();
	}
	h_data[ii][jj][ib]->Delete();
      }
    }
  }
}

void makeQCDComp(TString DIR, TString DIRqcd, TString channel, TString var) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+"Pre";
  
  // get histograms
  SummedHist* sigQCDMC = getQCDMC(DIR, var, "Pre", channel, false, "nom", false);
  SummedHist* sideQCDMC = getQCDMC(DIRqcd, var, "Pre", channel, true, "nom", false);
  TH1F* h_sigQCDMC = (TH1F*) sigQCDMC->hist();
  TH1F* h_sideQCDMC = (TH1F*) sideQCDMC->hist();

  TH1F* h_sideQCDData = (TH1F*) getQCDData(DIR,DIRqcd,var,"Pre",channel,"nom", false);

  h_sigQCDMC->SetLineColor(2);
  h_sigQCDMC->SetMarkerColor(2);
  h_sideQCDMC->SetLineColor(4);
  h_sideQCDMC->SetMarkerColor(4);
  h_sideQCDData->SetLineColor(1);
  h_sideQCDData->SetMarkerColor(1);
  h_sigQCDMC->Sumw2();
  h_sideQCDMC->Sumw2();
  h_sideQCDData->Sumw2();

  int rebinby = 1;
  if (hist.Contains("dR")) rebinby = 5;
  else if ((hist.Contains("ht") && !hist.Contains("htLep")) || hist.Contains("jetPt")) rebinby = 6;
  else if (hist.Contains("Phi")) rebinby = 7;
  else if (hist.Contains("Tau")) rebinby = 10;
  else if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))) rebinby = 5;

  h_sigQCDMC->Rebin(rebinby);
  h_sideQCDMC->Rebin(rebinby);
  h_sideQCDData->Rebin(rebinby);
  
  if (var == "ht"){
    TCanvas* c1 = new TCanvas("c1_"+hist,"c1_"+hist,900,600);
    h_sigQCDMC->SetFillColor(0);
    h_sigQCDMC->Draw("hist");
    c1->SaveAs("Plots/QCDtroubleshoot_sig_"+channel+"_"+hist+".pdf");
    h_sideQCDMC->SetFillColor(0);
    h_sideQCDMC->Draw("hist");
    c1->SaveAs("Plots/QCDtroubleshoot_side_"+channel+"_"+hist+".pdf");
  }    
  
  h_sigQCDMC->Scale(1./h_sigQCDMC->Integral());
  h_sideQCDMC->Scale(1./h_sideQCDMC->Integral());
  h_sideQCDData->Scale(1./h_sideQCDData->Integral());

  if (hist.Contains("ak8jetPt")) h_sideQCDData->GetXaxis()->SetRangeUser(400.,1200.);

  TH1F* h_ratio = (TH1F*) h_sideQCDData->Clone();
  TH1F* h_ratio2 = (TH1F*) h_sigQCDMC->Clone();
  TH1F* h_ratio3 = (TH1F*) h_sideQCDMC->Clone();
  h_ratio->Divide(h_sideQCDData);
  h_ratio2->Divide(h_sideQCDData);
  h_ratio3->Divide(h_sideQCDData);
  
  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
  p1->SetTopMargin(0.08);
  p1->SetBottomMargin(0.05);
  p1->SetNumber(1);
  TPad* p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
  p2->SetNumber(2);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.35);

  p1->Draw();
  p2->Draw();
  p1->cd();

  h_ratio->GetXaxis()->SetTitle(h_sigQCDMC->GetXaxis()->GetTitle());

  h_sideQCDData->GetXaxis()->SetTitle("");
  h_sideQCDData->SetMaximum(h_sideQCDData->GetMaximum()*1.4);
  
  h_sideQCDData->Draw();
  h_sideQCDMC->Draw("same");
  h_sigQCDMC->Draw("same");

  float xmin = 0.64;
  float xwidth = 0.2;
  float ymin = 0.70;
  float ywidth = 0.2;

  if (var.Contains("AbsEta") || var.Contains("Tau")) xmin = 0.2;
  if (var.Contains("ak4jetCSV")) xmin = 0.4;
  
  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_sideQCDData, "Data-driven QCD", "l");
  leg->AddEntry(h_sideQCDMC, "QCD MC, sideband", "l");
  leg->AddEntry(h_sigQCDMC, "QCD MC, sig. region", "l");
  leg->Draw();

  // plot ratio part
  p2->cd();
  //h_ratio->UseCurrentStyle();
  h_ratio->Draw();
  h_ratio2->Draw("same");
  h_ratio3->Draw("same");
  h_ratio->SetMaximum(2.0);
  h_ratio->SetMinimum(0.0);
  h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
  h_ratio->GetYaxis()->SetTitle("MC / Data-driven");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleOffset(1.0);
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleSize(0.1);

  // save output
  TString outname = "Plots/compQCD_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_sideQCDData->Delete();
  h_sideQCDMC->Delete();
  h_sigQCDMC->Delete();
  h_ratio->Delete();
  h_ratio2->Delete();
  h_ratio3->Delete();
}

// -------------------------------------------------------------------------------------
// print cutflow latex table
// -------------------------------------------------------------------------------------
void makeTable(TString DIR, TString DIRqcd, TString channel, bool inSideband, bool useQCDMC, bool unBlind, bool usePost = false) {

  const int nhist = 4;
  TString what = "lepPhi";
  TString regions[nhist] = {"Pre","0t","1t0b","1t1b"};

  // counts for cutflow
  float count_tt_semiLep[nhist] = {0};
  float count_tt_nonsemilep[nhist] = {0};
  float count_singletop[nhist] = {0};
  float count_diboson[nhist] = {0};
  float count_zjets[nhist] = {0};
  float count_wjets[nhist] = {0};
  float count_qcd[nhist] = {0};
  float count_tot[nhist] = {0};
  float count_data[nhist] = {0};
  
  // errors for cutflow
  float err_tt_semiLep[nhist] = {0};
  float err_tt_nonsemilep[nhist] = {0};
  float err_singletop[nhist] = {0};
  float err_diboson[nhist] = {0};
  float err_zjets[nhist] = {0};
  float err_wjets[nhist] = {0};
  float err_qcd[nhist] = {0};
  float err_tot[nhist] = {0};

  // Get count and error for each sample
  for (int ii = 0; ii < nhist; ii ++){
    // get histograms
    SummedHist* diboson = getDiboson(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* zjets = getZJets(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* wjets = getWJets(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* singletop = getSingleTop(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* ttbar = getTTbar(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep(DIR, what, regions[ii], channel, inSideband, "nom", usePost );
    SummedHist* data = getData(DIR, what, regions[ii], channel, inSideband);

    TH1F* h_diboson = (TH1F*) diboson->hist();
    TH1F* h_zjets = (TH1F*) zjets->hist();
    TH1F* h_wjets = (TH1F*) wjets->hist();
    TH1F* h_ttbar_semiLep = (TH1F*) ttbar->hist();
    TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
    TH1F* h_singletop = (TH1F*) singletop->hist();
    TH1F* h_data = (TH1F*) data->hist();

    TH1F* h_qcd;
    if(inSideband || useQCDMC) {
      SummedHist* qcd = getQCDMC(DIRqcd, what, regions[ii], channel, inSideband, "nom", usePost);
      h_qcd = (TH1F*) qcd->hist();
    }
    else h_qcd = (TH1F*) getQCDData(DIR,DIRqcd,what, regions[ii],channel,"nom", usePost);

    if (usePost){
      if (channel == "mu") h_qcd->Scale(0.71);
      else h_qcd->Scale(0.91);                 
      h_diboson->Scale(0.99);
      h_zjets->Scale(0.67);
      h_singletop->Scale(1.32);
      h_ttbar_semiLep->Scale(0.80);
      h_ttbar_nonSemiLep->Scale(0.80);
      h_wjets->Scale(0.97);
    }

    // error on pre-fit counts
    for (int ib=0; ib<h_ttbar_semiLep->GetNbinsX(); ib++) {
      if (h_ttbar_semiLep) err_tt_semiLep[ii]    += h_ttbar_semiLep->GetBinError(ib+1)*h_ttbar_semiLep->GetBinError(ib+1);
      if (h_ttbar_nonSemiLep) err_tt_nonsemilep[ii] += h_ttbar_nonSemiLep->GetBinError(ib+1)*h_ttbar_nonSemiLep->GetBinError(ib+1);
      if (h_singletop) err_singletop[ii]     += h_singletop->GetBinError(ib+1)*h_singletop->GetBinError(ib+1);
      if (h_diboson) err_diboson[ii]         += h_diboson->GetBinError(ib+1)*h_diboson->GetBinError(ib+1);
      if (h_zjets) err_zjets[ii]         += h_zjets->GetBinError(ib+1)*h_zjets->GetBinError(ib+1);
      if (h_wjets) err_wjets[ii]         += h_wjets->GetBinError(ib+1)*h_wjets->GetBinError(ib+1);
      if (h_qcd) err_qcd[ii]           += h_qcd->GetBinError(ib+1)*h_qcd->GetBinError(ib+1);
    }

    err_tt_semiLep[ii]    = sqrt(err_tt_semiLep[ii]);
    err_tt_nonsemilep[ii] = sqrt(err_tt_nonsemilep[ii]);
    err_singletop[ii]     = sqrt(err_singletop[ii]);
    err_diboson[ii]       = sqrt(err_diboson[ii]);
    err_zjets[ii]         = sqrt(err_zjets[ii]);
    err_wjets[ii]         = sqrt(err_wjets[ii]);
    err_qcd[ii]           = sqrt(err_qcd[ii]);
    
    err_tot[ii] = err_tt_semiLep[ii]*err_tt_semiLep[ii] + err_tt_nonsemilep[ii]*err_tt_nonsemilep[ii] + err_singletop[ii]*err_singletop[ii] + err_diboson[ii]*err_diboson[ii] + err_zjets[ii]*err_zjets[ii] + err_wjets[ii]*err_wjets[ii] + err_qcd[ii]*err_qcd[ii];
    err_tot[ii] = sqrt(err_tot[ii]);

    if (h_ttbar_semiLep) count_tt_semiLep[ii] = h_ttbar_semiLep->GetSum();
    if (h_ttbar_nonSemiLep) count_tt_nonsemilep[ii] = h_ttbar_nonSemiLep->GetSum();
    if (h_singletop) count_singletop[ii] = h_singletop->GetSum();
    if (h_diboson) count_diboson[ii] = h_diboson->GetSum();
    if (h_zjets) count_zjets[ii] = h_zjets->GetSum();
    if (h_wjets) count_wjets[ii] = h_wjets->GetSum();
    if (h_qcd) count_qcd[ii] = h_qcd->GetSum();
    count_tot[ii] = count_tt_semiLep[ii] + count_tt_nonsemilep[ii] + count_singletop[ii] + count_diboson[ii] + count_zjets[ii] + count_wjets[ii] + count_qcd[ii];
    if (h_data) count_data[ii] = h_data->GetSum();
    
  }

  // Print table

  cout << endl << "--------------------------------------------------" << endl;
  if (channel == "el") cout << "*** " << DIR << ", electron+jets channel ***" ;
  else cout << "*** " << DIR << ", muon+jets channel ***" ;
  if (inSideband) cout << " (sideband counts) ***";
  cout << endl;
  cout         << "--------------------------------------------------" << endl;
  //cout << setprecision(5);
  cout << "Sample                &     Preselection      &          0t          &         1t0b         &         1t1b         \\\\ " << endl;
  cout << "\\hline \\hline" << endl;
  cout << "\\ttbar (signal)       & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_tt_semiLep[ii]) << " $\\pm$ " << setw(4) << (int)round(err_tt_semiLep[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\ttbar (non-semilep)  & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_tt_nonsemilep[ii]) << " $\\pm$ " << setw(4) << (int)round(err_tt_nonsemilep[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "Single top            & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_singletop[ii]) << " $\\pm$ " << setw(4) << (int)round(err_singletop[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "W+jets                & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_wjets[ii]) << " $\\pm$ " << setw(4) << (int)round(err_wjets[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "Z+jets                & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_zjets[ii]) << " $\\pm$ " << setw(4) << (int)round(err_zjets[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "Diboson               & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_diboson[ii]) << " $\\pm$ " << setw(4) << (int)round(err_diboson[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "QCD                   & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_qcd[ii]) << " $\\pm$ " << setw(4) << (int)round(err_qcd[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline" << endl;
  cout << "Total                 & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(5) << (int)round(count_tot[ii]) << " $\\pm$ " << setw(4) << (int)round(err_tot[ii]);
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline \\hline" << endl;
  cout << "Data                  &       ";
  for (int ii = 0; ii < nhist; ii++){
    if (ii == 0 || unBlind || inSideband) cout << setw(5) << count_data[ii];
    else cout << "  N/A  ";
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << "      &      ";
  }
}

void plotWJetsSplit(TString var, TString region, TString channel){
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+region;
    
  // get histograms
  SummedHist* wjets_bb = getWJets( "histfiles_full2016/", var, region, channel, false, "nom", false, "bb");
  SummedHist* wjets_b  = getWJets( "histfiles_full2016/", var, region, channel, false, "nom", false, "b");
  SummedHist* wjets_cc = getWJets( "histfiles_full2016/", var, region, channel, false, "nom", false, "cc");
  SummedHist* wjets_c  = getWJets( "histfiles_full2016/", var, region, channel, false, "nom", false, "c");
  SummedHist* wjets_l  = getWJets( "histfiles_full2016/", var, region, channel, false, "nom", false, "l");

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  TH1F* h_bb = (TH1F*) wjets_bb->hist();
  TH1F* h_b  = (TH1F*) wjets_b->hist();
  TH1F* h_cc = (TH1F*) wjets_cc->hist();
  TH1F* h_c  = (TH1F*) wjets_c->hist();
  TH1F* h_l  = (TH1F*) wjets_l->hist();

  int rebinfac;
  if (hist.Contains("lepEta") || hist.Contains("ak8jetTau") || hist.Contains("ak4jetEta")) rebinfac = 5;
  else rebinfac = 2;
  
  h_bb->Rebin(rebinfac);
  h_b->Rebin(rebinfac);
  h_cc->Rebin(rebinfac);
  h_c->Rebin(rebinfac);
  h_l->Rebin(rebinfac);

  // create stack & summed histogram for ratio plot
  THStack* h_stack = new THStack();    
  h_stack->Add(h_bb);
  h_stack->Add(h_b);
  h_stack->Add(h_cc);
  h_stack->Add(h_c);
  h_stack->Add(h_l);

  TH1F* h_totalHF = (TH1F*) h_bb->Clone("total_"+hist);
  h_totalHF->Add(h_b);
  h_totalHF->Add(h_cc);
  h_totalHF->Add(h_c);
  h_totalHF->Scale(1.0/h_totalHF->Integral());

  // -------------------------------------------------------------------------------------
  
  TH1F* h_ratio;
  TH1F* h_ratio2;
  h_ratio = (TH1F*) h_l->Clone("ratio_"+hist);  // Data / MC
  h_ratio->Scale(1.0/h_ratio->Integral());
  h_ratio->Divide(h_totalHF);
  
  h_ratio2 = (TH1F*) h_totalHF->Clone("ratio2_"+hist); // Uncertainty on Data / MC
  for (int ib=0; ib<h_totalHF->GetNbinsX(); ib++) {
    h_ratio2->SetBinContent(ib+1, 1.0);
    float tmperr = h_totalHF->GetBinError(ib+1);
    float tmpcount = h_totalHF->GetBinContent(ib+1);
    
    if (tmpcount==0) continue;
    h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
  }

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
  p1->SetTopMargin(0.08);
  p1->SetBottomMargin(0.05);
  p1->SetNumber(1);
  TPad* p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
  p2->SetNumber(2);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.35);
  
  p1->Draw();
  p2->Draw();
  p1->cd();
  
  h_ratio2->SetMarkerSize(0);
  h_ratio2->SetLineColor(0);
  h_ratio2->SetFillColor(15);
  h_ratio2->SetFillStyle(1001);
  
  h_ratio->GetXaxis()->SetTitle(h_totalHF->GetXaxis()->GetTitle());
  h_stack->Draw("hist");

  float xmin = 0.71;
  float ymin = 0.52;

  float xwidth = 0.20;
  float ywidth = 0.38;

  //Legend top left
  if (var == "ak8jetTau32" || var == "ak8jetTau21" || var == "ak8jetSDmass") xmin = 0.16;

  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_bb,     "W+bb", "f");
  leg->AddEntry(h_b,      "W+b", "f");
  leg->AddEntry(h_cc,     "W+cc", "f");
  leg->AddEntry(h_c,      "W+c", "f");
  leg->AddEntry(h_l,      "W+l" , "f");
  leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
  leg->Draw();

  myText(0.10,0.94,1,"#intLdt = 35.9 fb^{-1}");
  myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");

  // plot ratio part
  p2->cd();
  p2->SetGridy();
  h_ratio->UseCurrentStyle();
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(1);
  h_ratio->Draw("le0p");
  h_ratio2->Draw("same,e2");
  h_ratio->Draw("le0p,same");
  h_ratio->SetMaximum(1.4);
  h_ratio->SetMinimum(0.6);
  h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
  h_ratio->GetYaxis()->SetTitle("l / HF");
  h_ratio->GetXaxis()->SetLabelSize(26);
  h_ratio->GetYaxis()->SetLabelSize(26);
  h_ratio->GetXaxis()->SetTitleOffset(2.8);
  h_ratio->GetYaxis()->SetTitleOffset(1.4);
  h_ratio->GetXaxis()->SetTitleSize(32);
  h_ratio->GetYaxis()->SetTitleSize(32);

  // save output
  TString outname = "Plots/splitWJets_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_stack->Delete();
  leg->Delete();

}

void calcBtagSF(TString channel, TString sys){
  SummedHist* ttbar            = getTTbar("histfiles_full2016/", "bTagSFvsPt", "", channel, false, sys, false);
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep("histfiles_full2016/", "bTagSFvsPt", "", channel, false, sys, false);
  SummedHist* singletop        = getSingleTop("histfiles_full2016/", "bTagSFvsPt", "", channel, false, sys, false);
  SummedHist* wjets            = getWJets("histfiles_full2016/", "bTagSFvsPt", "", channel, false, sys, false);
  SummedHist* qcd              = getQCDMC("histfiles_full2016/", "bTagSFvsPt", "", channel, false, sys, false);

  TH2F* h_ttbar = (TH2F*) ttbar->hist();
  TH2F* h_ttbar_nonSemiLep = (TH2F*) ttbar_nonSemiLep->hist();
  TH2F* h_singletop = (TH2F*) singletop->hist();
  TH2F* h_wjets = (TH2F*) wjets->hist();
  TH2F* h_qcd = (TH2F*) qcd->hist();

  TH2F* h_total = (TH2F*) h_ttbar->Clone();
  h_total->Add(h_ttbar_nonSemiLep);
  h_total->Add(h_singletop);
  h_total->Add(h_wjets);
  h_total->Add(h_qcd);

  double avgSF = h_total->GetMean(1);
  cout << "Average SF in " << channel << " channel for " << sys << " variation is " << avgSF << endl;
  
}

void plot2D(TString sample, TString hist){
  
  SummedHist* hsum;
  if (sample == "ttbar") hsum = getTTbar("histfiles_full2016/", hist, "", "el", false, "nom", false);
  if (sample == "qcd") hsum = getQCDMC("histfiles_full2016/", hist, "", "el", false, "nom", false);
  if (sample == "data") hsum = getData("histfiles_full2016/", hist, "", "el", false);

  TH2F* h_tmp = (TH2F*) hsum->hist();
  cout << "Number of events in " << hist << " is " << h_tmp->GetSum() << endl;

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c","c",900,700);
  c->SetLeftMargin(0.18);
  c->SetRightMargin(0.12);
  c->SetBottomMargin(0.12);
  h_tmp->GetZaxis()->SetLabelSize(0.04);
  h_tmp->Draw("COLZ");
  gPad->Update();
  TPaletteAxis *ttbar_palette = (TPaletteAxis*)h_tmp->GetListOfFunctions()->FindObject("palette");
  ttbar_palette->SetX1NDC(0.88);
  ttbar_palette->SetX2NDC(0.93);
  ttbar_palette->SetY1NDC(0.12);
  ttbar_palette->SetY2NDC(0.9);
  gPad->Modified();
  gPad->Update();
  c->SaveAs("Plots/"+hist+"_"+sample+".pdf");
}

void combineResults(TString channel, TString fit) {

  TH1::AddDirectory(kFALSE); 
  setStyle();

  const int nhist = 3;
  TString what[nhist] = {"ak4jetEta0t","ak4jetEta1t0b","ak8jetSDmass1t1b"};
  //TString what[nhist] = {"ak4jetAbsEta0t","ak4jetAbsEta1t0b","ak8jetSDmass1t1b"};
  //TString what[nhist] = {"lepEta0t","lepEta1t0b","ak8jetSDmass1t1b"};
  //TString what[nhist] = {"lepEta0t_barrel","lepEta0t_endcap","lepEta1t0b_barrel","lepEta1t0b_endcap","ak8jetSDmass1t1b_barrel","ak8jetSDmass1t1b_endcap"};
  //TString what[nhist] = {"lepAbsEta0t","lepAbsEta1t0b","ak8jetSDmass1t1b"};
  //TString what[nhist] = {"lepAbsEta0t_barrel","lepAbsEta0t_endcap","lepAbsEta1t0b_barrel","lepAbsEta1t0b_endcap","ak8jetSDmass1t1b_barrel","ak8jetSDmass1t1b_endcap"};
  //TString what[nhist] = {"ak8jetTau210t","ak8jetTau321t0b","ak8jetSDmass1t1b"};
  //TString what[nhist] = {"counts"};
  //TString what[nhist] = {"counts_barrel","counts_endcap"};
  const int ncats = 7;
  TString cats[ncats] = {"TTbar","SingleTop","WJets","ZJets","Diboson","QCD","total"};
  bool mergewjets = true;
  bool mergeST = true;
  TString whichQCD = "data";

  // counts and errors for cutflow
  float counts[2][nhist][ncats] = {0};
  float errs[2][nhist][ncats] = {0};
  float counts_data[nhist] = {0};

  // Get prefit hists and counts
  TFile* f_prefit = TFile::Open("combineInputs_"+whichQCD+".root");
  TFile* f_postfit = TFile::Open( "fitDiagnostics"+fit+".root" );
  RooArgSet* norm_post = (RooArgSet*) f_postfit->Get("norm_fit_s");

  for (int ih = 0; ih < nhist; ih++){
    TH1F* h_data    = (TH1F*) f_prefit->Get(what[ih]+"_"+channel+"/data_obs");
    counts_data[ih] = h_data->GetSum();

    TH1F* hists[2][ncats];
    for (int ic = 0; ic < ncats-1; ic++){
      if (cats[ic] == "QCD" && channel == "el" && what[ih] == "ak8jetSDmass1t1b_endcap") continue;
      hists[0][ic] = (TH1F*) f_prefit->Get(what[ih]+"_"+channel+"/"+cats[ic]);
      if (ic == 0) hists[0][ncats-1] = (TH1F*) hists[0][ic]->Clone();
      else hists[0][ncats-1]->Add(hists[0][ic]);
    }

    for (int ic = 0; ic < ncats; ic++){
      if (cats[ic] == "QCD" && channel == "el" && what[ih] == "ak8jetSDmass1t1b_endcap") continue;
      hists[1][ic]     = (TH1F*) hists[0][ic]->Clone();

      // These have correct bin contents, sys errors; wrong binning and no stat. errors
      TH1F* h_tmp;
      if (mergewjets && cats[ic] == "WJets") {
	h_tmp = (TH1F*) f_postfit->Get("shapes_fit_s/"+what[ih]+"_"+channel+"/"+cats[ic]+"L");
	h_tmp->Add((TH1F*) f_postfit->Get("shapes_fit_s/"+what[ih]+"_"+channel+"/"+cats[ic]+"HF"));
      }
      else if (mergeST && cats[ic] == "SingleTop") {
	h_tmp = (TH1F*) f_postfit->Get("shapes_fit_s/"+what[ih]+"_"+channel+"/ST_tW");
	h_tmp->Add((TH1F*) f_postfit->Get("shapes_fit_s/"+what[ih]+"_"+channel+"/ST_other"));
      }
      else h_tmp = (TH1F*) f_postfit->Get("shapes_fit_s/"+what[ih]+"_"+channel+"/"+cats[ic]);

      counts[0][ih][ic] = hists[0][ic]->GetSum();
      counts[1][ih][ic] = h_tmp->GetSum();

      // Get systematic postfit errors from 'norm_fit_s' object -- will add stat. later
      if (ic != ncats-1){
	if (mergewjets && cats[ic] == "WJets"){
	  RooRealVar* thisnorm1 = (RooRealVar*) norm_post->find(what[ih]+"_"+channel+"/"+cats[ic]+"L");
	  RooRealVar* thisnorm2 = (RooRealVar*) norm_post->find(what[ih]+"_"+channel+"/"+cats[ic]+"HF");
	  errs[1][ih][ic] = std::pow(thisnorm1->getError(),2)+std::pow(thisnorm2->getError(),2);
	}
	else if (mergeST && cats[ic] == "SingleTop"){
	  RooRealVar* thisnorm1 = (RooRealVar*) norm_post->find(what[ih]+"_"+channel+"/ST_tW");
	  RooRealVar* thisnorm2 = (RooRealVar*) norm_post->find(what[ih]+"_"+channel+"/ST_other");
	  errs[1][ih][ic] = std::pow(thisnorm1->getError(),2)+std::pow(thisnorm2->getError(),2);
	}
	else {
	  RooRealVar* thisnorm = (RooRealVar*) norm_post->find(what[ih]+"_"+channel+"/"+cats[ic]);
	  errs[1][ih][ic] = std::pow(thisnorm->getError(),2);
	}
      }

      // Construct postfit hists with correct binning, uncertainties
      for (int ib=0; ib<h_data->GetNbinsX(); ib++) {
	hists[1][ic]->SetBinContent(ib+1,h_tmp->GetBinContent(ib+1));
	
	// Set posterior error as sum in quadrature of statistical uncertainty (from prefit)
	// scaled by post/pre norm ratio and systematic uncertainty (from postfit)
	// This is not quite right if there are large shape shifts in the nominal, but approximately correct
	hists[1][ic]->SetBinError(ib+1,sqrt(std::pow(h_tmp->GetBinError(ib+1),2)
					    +std::pow(hists[0][ic]->GetBinError(ib+1)*counts[1][ih][ic]/counts[0][ih][ic],2)));
	
	// Prefit uncertainties must be taken from histograms, not norm_prefit -- norm_prefit has 'prefit' systematic uncertainties
	errs[0][ih][ic] += std::pow(hists[0][ic]->GetBinError(ib+1),2);
      }

      // Add postfit statistical uncertainties (prefit statistical scaled by post/pre ratio)
      if (ic != ncats-1){
	errs[1][ih][ic] += errs[0][ih][ic] * std::pow(counts[1][ih][ic]/counts[0][ih][ic],2);
	errs[1][ih][ncats-1] += errs[1][ih][ic];
	errs[0][ih][ncats-1] += errs[0][ih][ic];
	errs[1][ih][ic] = sqrt(errs[1][ih][ic]);
	errs[0][ih][ic] = sqrt(errs[0][ih][ic]);
      }
    }

    if (channel == "el" && what[ih] == "ak8jetSDmass1t1b_endcap"){
      errs[0][ih][3] = 0.0;
      errs[1][ih][3] = 0.0;
      counts[0][ih][3] = 0.0;
      counts[1][ih][3] = 0.0;
    }

    errs[0][ih][ncats-1] = sqrt(errs[0][ih][ncats-1]);
    errs[1][ih][ncats-1] = sqrt(errs[1][ih][ncats-1]);

    TString histtitle = h_data->GetXaxis()->GetTitle();

    // Now we've gotten all the counts and errors, move on to plotting
    int colors[ncats-1] = {kRed+1,6,kGreen-3,kYellow,kAzure-9,kViolet-6};
    for (int ff = 0; ff < 2; ff++){
      THStack* h_stack = new THStack();
      for (int ic = 0; ic < ncats-1; ic++){
	if (hists[ff][5-ic]) {
	  hists[ff][5-ic]->SetFillColor(colors[5-ic]);
	  h_stack->Add(hists[ff][5-ic]);
	}
      }

      h_data->SetBinErrorOption(TH1::kPoisson);
      
      TH1F* h_ratio;
      TH1F* h_ratio2;
      h_ratio = (TH1F*) h_data->Clone("ratio_"+what[ih]+"_"+channel);  // Data / MC
      h_ratio->Sumw2();
      h_ratio->Divide(hists[ff][ncats-1]);
      
      h_ratio2 = (TH1F*) hists[ff][ncats-1]->Clone("ratio2_"+what[ih]+"_"+channel); // Uncertainty on Data / MC
      h_ratio2->Sumw2();
      for (int ib=0; ib<h_ratio2->GetNbinsX(); ib++) {
	float tmperr = h_ratio2->GetBinError(ib+1);
	float tmpcount = h_ratio2->GetBinContent(ib+1);
	h_ratio2->SetBinContent(ib+1, 1.0);
	if (tmpcount <= 0.00001) h_ratio2->SetBinError(ib+1,0.0);
	else h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
      }

      float mymax = 0;
      if (ff==0) mymax = max(h_data->GetMaximum(),hists[ff][ncats-1]->GetMaximum());
      else mymax = hists[ff][ncats-1]->GetMaximum();

      h_data->SetAxisRange(0,mymax*1.4,"Y");


      // -------------------------------------------------------------------------------------
      // fix axis labels

      if (what[ih].Contains("ak4jetEta")) 
	h_data->GetYaxis()->SetTitle("Events / 0.5");
      else if (what[ih].Contains("ak8jetSDmass1t1b"))
	h_data->GetYaxis()->SetTitle("Events / 10 GeV");


      // -------------------------------------------------------------------------------------
      // plotting!

      TCanvas* c = new TCanvas("c","c",900,800);
      TPad* p1 = new TPad("datamcp1_"+what[ih]+"_"+channel,"datamcp1_"+what[ih]+"_"+channel,0,0.28,1,1);
      p1->SetTopMargin(0.1);
      p1->SetBottomMargin(0.02);
      p1->SetNumber(1);
      TPad* p2 = new TPad("datamcp2_"+what[ih]+"_"+channel,"datamcp2_"+what[ih]+"_"+channel,0,0,1,0.28);
      p2->SetNumber(2);
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(0.35);
      
      p1->Draw();
      p2->Draw();
      p1->cd();
      
      hists[ff][ncats-1]->UseCurrentStyle();
      hists[ff][ncats-1]->SetFillColor(0);
      hists[ff][ncats-1]->SetLineWidth(1);
      hists[ff][ncats-1]->SetLineColor(1);

      h_ratio2->SetMarkerSize(0);
      h_ratio2->SetLineColor(0);
      h_ratio2->SetFillColor(15);
      h_ratio2->SetFillStyle(1001);
      
      h_ratio->GetXaxis()->SetTitle(histtitle);
      
      h_data->UseCurrentStyle();
      h_data->SetMarkerColor(1);
      h_data->SetMarkerStyle(8);
      h_data->SetMarkerSize(1);
      h_data->GetXaxis()->SetLabelSize(0);
      h_data->GetYaxis()->SetLabelSize(26);
      h_data->GetYaxis()->SetTitleSize(32);
      h_data->GetYaxis()->SetTitleOffset(1.4);
      h_data->GetXaxis()->SetTitle("");
      
      h_data->Draw("LE0P");
      hists[ff][ncats-1]->Draw("hist,same");
      h_stack->Draw("hist,same");
      h_data->Draw("LE0P,same");
      h_data->Draw("axis,same");

      float xmin = 0.73;
      float ymin = 0.48;
      float xwidth = 0.20;
      float ywidth = 0.40;
      
      if (what[ih] == "ak8jetSDmass1t1b" || what[ih].Contains("ak8jetTau")) xmin = 0.18;

      // legend
      TLegend* leg;
      leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->AddEntry(h_data, "Data", "pe");
      leg->AddEntry(hists[ff][0], "t#bar{t}", "f");
      leg->AddEntry(hists[ff][1], "Single t", "f");
      leg->AddEntry(hists[ff][2], "W+jets", "f");
      leg->AddEntry(hists[ff][3], "Z+jets", "f");
      leg->AddEntry(hists[ff][4], "Diboson", "f");
      if (hists[ff][5]) leg->AddEntry(hists[ff][5], "Multijet" , "f");
      leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
      leg->Draw();

      //myText(0.10,0.94,1,"#intLdt = 35.9 fb^{-1}");
      //myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");
      drawCMS(0.15,0.92, false);

      float whereput=0.2;
      if (what[ih] == "ak8jetSDmass1t1b") whereput=0.82;
      if (channel=="el") myLargeText(whereput,0.82,1,"e+jets");
      else myLargeText(whereput,0.82,1,"#mu+jets");


      // plot ratio part
      p2->cd();
      p2->SetGridy();
      h_ratio->UseCurrentStyle();
      h_ratio->SetMarkerStyle(8);
      h_ratio->SetMarkerSize(1);
      h_ratio->Draw("le0p");
      h_ratio2->Draw("same,e2");
      h_ratio->Draw("le0p,same");
      h_ratio->SetMaximum(1.75);
      h_ratio->SetMinimum(0.25);
      h_ratio->GetYaxis()->SetNdivisions(4,4,0,true);
      h_ratio->GetYaxis()->SetTitle("Data / MC");
      h_ratio->GetXaxis()->SetLabelSize(28);
      h_ratio->GetYaxis()->SetLabelSize(26);
      h_ratio->GetXaxis()->SetTitleOffset(3.5);
      h_ratio->GetYaxis()->SetTitleOffset(1.2);
      h_ratio->GetXaxis()->SetTitleSize(32);
      h_ratio->GetYaxis()->SetTitleSize(32);

      // save output
      TString preorpost = "pre";
      if (ff == 1) preorpost = "post";
      TString outname = "Plots/"+channel+"_"+what[ih]+"_"+fit+"_"+preorpost+".pdf";
      c->SaveAs(outname);
      c->Close();
      
    } // End pre / post plotting loop
  } // End histogram loop

  // Now plot correlation matrix
  TH2F* h_mcorr_raw = (TH2F*) f_postfit->Get("covariance_fit_s"); //Note that this has the per-bin stat uncertainties if autoMCStats is on
  int countnuis = 0;
  for (int ibin = 1; ibin < h_mcorr_raw->GetNbinsX()+1; ibin++){
    TString binlabel = h_mcorr_raw->GetXaxis()->GetBinLabel(ibin);
    if (!(binlabel.Contains("bin"))) countnuis++;
  }
  TH2F* h_mcorr = new TH2F("covariance_final","",countnuis,0,countnuis,countnuis,0,countnuis);
  int mybinx = 0;
  for (int ibinx = 1; ibinx < h_mcorr_raw->GetNbinsX()+1; ibinx++){
    TString binlabelx = h_mcorr_raw->GetXaxis()->GetBinLabel(ibinx);
    if (binlabelx.Contains("bin")) continue;
    h_mcorr->GetXaxis()->SetBinLabel(mybinx+1,binlabelx);
    h_mcorr->GetYaxis()->SetBinLabel(h_mcorr->GetNbinsX()-mybinx,binlabelx);
    int mybiny = 0;
    for (int ibiny = 1; ibiny < h_mcorr_raw->GetNbinsY()+1; ibiny++){
      TString binlabely = h_mcorr_raw->GetYaxis()->GetBinLabel(ibiny);
      if (binlabely.Contains("bin")) continue;
      h_mcorr->SetBinContent(mybinx+1,mybiny+1,h_mcorr_raw->GetBinContent(ibinx,ibiny));
      mybiny++;
    }
    mybinx++;
  }
    
  TCanvas* c2 = new TCanvas("c2","c2",900,700);
  c2->SetLeftMargin(0.18);
  c2->SetRightMargin(0.12);
  c2->SetBottomMargin(0.1);
  h_mcorr->GetZaxis()->SetLabelSize(0.04);
  h_mcorr->Draw("COLZ");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)h_mcorr->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.88);
  palette->SetX2NDC(0.93);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  gPad->Modified();
  gPad->Update();
  c2->SaveAs("Plots/correlation_"+fit+".pdf");
  c2->Close();
  
  // Finally make tables
  for (int ff = 0; ff < 2; ff++){
    cout << endl << "--------------------------------------------------" << endl;
    if (ff == 1){
      if (channel == "el") cout << "*** Postfit event counts, electron+jets channel ***" ;
      else cout << "*** Postfit event counts, muon+jets channel ***" ;
    }
    else {
      if (channel == "el") cout << "*** Prefit event counts, electron+jets channel ***" ;
      else cout << "*** Prefit event counts, muon+jets channel ***" ;
    }
    cout << endl;
    cout << fit << " fit in " << channel << " channel" << endl;
    cout         << "--------------------------------------------------" << endl;
    cout << setprecision(5);
    cout << "Sample      &          0t          &         1t0b         &         1t1b         \\\\ " << endl;
    cout << "\\hline \\hline" << endl;
    for (int ic = 0; ic < ncats; ic++){
      if (ic == ncats-1) cout << "\\hline" << endl;
      cout << setw(11) << cats[ic] << " & ";
      for (int ih = 0; ih < nhist; ih++){
	cout << setw(5) << (int)round(counts[ff][ih][ic]) << " $\\pm$ " << setw(4) << (int)round(errs[ff][ih][ic]);
	if (ih == nhist - 1) cout << " \\\\ " << endl;
	else cout << " & ";
      }
    }
    cout << "\\hline \\hline" << endl;
    cout << "Data        &       ";
    for (int ih = 0; ih < nhist; ih++){
      cout << setw(5) << counts_data[ih];
      if (ih == nhist - 1) cout << " \\\\ " << endl;
      else cout << "      &      ";
    } //end data line
  }
}

void make2DScans(TString DIR, TString channel) {

  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hists[7] = {"2DisoPt15","2DisoPt20","2DisoPt25","2DisoPt30","2DisoPt35","2DisoPt40","2DisoPt45"};
  
  // get histograms
  for (int ii = 0; ii < 7; ii ++){
    SummedHist* signal = getTTbar(DIR, hists[ii], "", channel, false, "nom", false);
    TH2F* h_signal = (TH2F*) signal->hist();
    TH2F* h_bkg = (TH2F*) getBackground( DIR, hists[ii], "", channel, false);

    //Construct test statistic per bin
    TH2F* h_sigma = (TH2F*) h_signal->Clone();
    h_sigma->Reset();
    int nbinsx = h_sigma->GetNbinsX();
    int nbinsy = h_sigma->GetNbinsY();
    for (int jj = 1; jj < nbinsx; jj++){
      for (int kk = 1; kk < nbinsy; kk++){
	float n_sig = h_signal->GetSum() - h_signal->Integral(1,jj,1,kk); //Number of signal events
	float n_bkg = h_bkg->GetSum() - h_bkg->Integral(1,jj,1,kk); //Number of background events
	float sigma = n_sig / sqrt(n_sig + n_bkg);
	h_sigma->SetBinContent(jj,kk,sigma);
      }
    }
    
    // Plot
    TCanvas* c = new TCanvas("c_"+hists[ii],"c_"+hists[ii],900,700);
    c->SetRightMargin(1.2);
    h_signal->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_signal.pdf");
    h_bkg->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_bkg.pdf");
    h_sigma->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_sigma.pdf");

  }
}

void make1DScans(TString channel) {

  TH1::AddDirectory(kFALSE); 
  setStyle();

  const int nHIST = 6;
  //TString hists[5] = {"metPt","ht","htLep","lepBJetdR","lepTJetdR"};
  TString hists[nHIST] = {"metPtPre","lepMETdPhiRaw","ak4METdPhiRaw","ak8METdPhiRaw","eMetTriRaw","jetMetTriRaw"};
  
  // get histograms
  for (int ii = 0; ii < nHIST; ii ++){
    if (channel == "mu" && hists[ii].Contains("Tri")) continue;
    
    SummedHist* signal = getTTbar("histfiles_full2016/", hists[ii], "", channel, false, "nom", false);
    TH1F* h_signal = (TH1F*) signal->hist(true);
    cout << "Got signal" << endl;
    //TH1F* h_bkg = (TH1F*) getBackground( "histfiles_full2016/", hists[ii], "", channel, false);
    SummedHist* bkg = getQCDMC("histfiles_full2016/", hists[ii], "", channel, false, "nom", false);
    TH1F* h_bkg = (TH1F*) bkg->hist(true);
    cout << "Got background" << endl;

    //Construct test statistic per bin
    TH1F* h_sigma = (TH1F*) h_signal->Clone();
    h_sigma->Reset();
    int Nbins = h_sigma->GetNbinsX();
    for (int jj = 1; jj < Nbins+1; jj++){
      float n_sig = h_signal->Integral(jj,Nbins+1); //Number of signal events
      float n_bkg = h_bkg->Integral(jj,Nbins+1); //Number of background events
      if (hists[ii] == "lepBJetdRPre"){
	n_sig = h_signal->Integral(0,jj-1); //Number of signal events
	n_bkg = h_bkg->Integral(0,jj-1); //Number of background events
      }

      float sigma = 0.0;
      if (n_sig > 0.0 && n_bkg > 0.0) sigma = n_sig / sqrt(n_sig + n_bkg);
      h_sigma->SetBinContent(jj,sigma);
    }
    
    // Plot
    TCanvas* c = new TCanvas("c_"+hists[ii],"c_"+hists[ii],900,700);
    h_sigma->GetYaxis()->SetTitle("S / sqrt(S+B)");
    h_sigma->GetYaxis()->SetTitleOffset(1.0);
    h_sigma->Draw("hist");
    c->SaveAs("Plots/"+channel+"_"+hists[ii]+"_sigma.pdf");

  }
}

void makeEffPlots(TString channel, TString which) {
  setStyle();
  TH1::AddDirectory(kFALSE); 

  const int niso = 6;
  TString isos[niso] = {"2DisoPt25","2DisoPt45","2DisoB2G","2DisoIHNY","MiniIso10","MiniIso20"};
  int colors[niso] = {1,2,3,4,6,7};
  TH1F* h_signal_eff[niso];
  TH1F* h_bkg_eff[niso];
  
  for (int ii = 0; ii < niso; ii ++){
    // Get signal
    SummedHist* signal_pre = getTTbar("histfiles_mMu_mEl_Loose/", which, "Pre", channel, false, "nom", false ); //These should be for no iso
    SummedHist* signal_post = getTTbar("histfiles_mMu_mEl_"+isos[ii]+"/", which, "Pre", channel, false, "nom", false ); //These should be with iso
    TH1F* h_signal_pre = (TH1F*) signal_pre->hist(true);
    TH1F* h_signal_post = (TH1F*) signal_post->hist(true);

    // Get QCD
    SummedHist* bkg_pre = getQCDMC("histfiles_mMu_mEl_Loose/",which, "Pre",channel,false,"nom", false);
    SummedHist* bkg_post = getQCDMC("histfiles_mMu_mEl_"+isos[ii]+"/",which, "Pre",channel,false,"nom",false);
    TH1F* h_bkg_pre = (TH1F*) bkg_pre->hist(true);
    TH1F* h_bkg_post = (TH1F*) bkg_post->hist(true);

    h_signal_pre->Rebin(2);
    h_bkg_pre->Rebin(2);
    h_signal_post->Rebin(2);
    h_bkg_post->Rebin(2);
     
    //Construct signal efficiency
    h_signal_eff[ii] = (TH1F*) h_signal_post->Clone();
    h_signal_eff[ii]->Divide(h_signal_pre);
    if (which == "ak8jetPt") h_signal_eff[ii]->GetXaxis()->SetRangeUser(400.,1200.);
    if (which == "leadLepPt") h_signal_eff[ii]->GetXaxis()->SetRangeUser(50.,500.);
    h_bkg_eff[ii] = (TH1F*) h_bkg_post->Clone();
    h_bkg_eff[ii]->Divide(h_bkg_pre);
    if (which == "ak8jetPt") h_bkg_eff[ii]->GetXaxis()->SetRangeUser(400.,1200.);
    if (which == "leadLepPt") h_bkg_eff[ii]->GetXaxis()->SetRangeUser(50.,500.);
  }

  // Plot
  TCanvas* c = new TCanvas("c","c",900,700);

  float xcor = 0.2;
  if (which == "leadLepPt") xcor = 0.7;  
  TLegend* leg1 = new TLegend(xcor,0.15,xcor+0.2,0.45);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.045);
  leg1->AddEntry(h_signal_eff[0], isos[0], "l");

  h_signal_eff[0]->SetLineColor(colors[0]);
  h_signal_eff[0]->GetYaxis()->SetRangeUser(0.5,1.0);
  h_signal_eff[0]->GetYaxis()->SetTitle("Signal (semilep ttbar) efficiency");
  h_signal_eff[0]->GetXaxis()->SetTitleOffset(1.1);
  h_signal_eff[0]->Draw("hist");
  
  for (int jj = 1; jj < niso; jj ++){
    h_signal_eff[jj]->SetLineColor(colors[jj]);
    h_signal_eff[jj]->Draw("hist,same");
    leg1->AddEntry(h_signal_eff[jj], isos[jj], "l");
  }
  leg1->Draw();
  c->SaveAs("Plots/eff_sig_"+which+"_"+channel+".pdf");

  TLegend* leg2 = new TLegend(xcor,0.58,xcor+0.2,0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.045);
  leg2->AddEntry(h_bkg_eff[0], isos[0], "l");

  h_bkg_eff[0]->SetLineColor(colors[0]);
  h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,1.0);
  h_bkg_eff[0]->GetYaxis()->SetTitle("Background (QCD) efficiency");
  h_bkg_eff[0]->GetXaxis()->SetTitleOffset(1.1);
  if (channel == "mu" && which == "ak8jetPt")  h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,0.1);
  h_bkg_eff[0]->Draw("hist");
  
  for (int jj = 1; jj < niso; jj ++){
    h_bkg_eff[jj]->SetLineColor(colors[jj]);
    h_bkg_eff[jj]->Draw("hist,same");
    leg2->AddEntry(h_bkg_eff[jj], isos[jj], "l");
  }
  leg2->Draw();
  c->SaveAs("Plots/eff_bkg_"+which+"_"+channel+".pdf");
}

void makeEffPlotsFinal() {
  setStyle();
  TH1::AddDirectory(kFALSE); 

  TString channels[2] = {"mu","el"};
  TString hists[2] = {"ak8jetPt","leadLepPt"};
  
  for (int ii = 0; ii < 2; ii++){
    TH1F* h_signal_eff[2];
    TH1F* h_bkg_eff[2];
    for (int jj = 0; jj < 2; jj++){
      // get histograms
      SummedHist* signal_pre = getTTbar("histfiles_mMu_mEl_Loose/", hists[ii], "Pre", channels[jj], false, "nom", false);
      SummedHist* signal_post = getTTbar("histfiles_mMu_tEl_MiniIso10/", hists[ii], "Pre", channels[jj], false, "nom", false);
      TH1F* h_signal_pre = (TH1F*) signal_pre->hist(true);
      TH1F* h_signal_post = (TH1F*) signal_post->hist(true);
      
      //This version uses QCD only
      SummedHist* bkg_pre = getQCDMC("histfiles_mMu_mEl_Loose/",hists[ii],"Pre",channels[jj],false,"nom",false);
      SummedHist* bkg_post = getQCDMC("histfiles_mMu_tEl_MiniIso10/",hists[ii],"Pre",channels[jj],false,"nom",false);
      TH1F* h_bkg_pre = (TH1F*) bkg_pre->hist(true);
      TH1F* h_bkg_post = (TH1F*) bkg_post->hist(true);

      h_signal_pre->Rebin(2);
      h_bkg_pre->Rebin(2);
      h_signal_post->Rebin(2);
      h_bkg_post->Rebin(2);
     
      //Construct signal efficiency
      h_signal_eff[jj] = (TH1F*) h_signal_post->Clone();
      h_signal_eff[jj]->SetName("signal_eff_"+channels[jj]);
      h_signal_eff[jj]->GetYaxis()->SetTitle("Signal Efficiency");
      h_signal_eff[jj]->Divide(h_signal_post,h_signal_pre,1.0,1.0,"B");
      if (hists[ii] == "ak8jetPt") h_signal_eff[jj]->GetXaxis()->SetRangeUser(400.,1200.);
      if (hists[ii] == "leadLepPt") h_signal_eff[jj]->GetXaxis()->SetRangeUser(50.,500.);
      h_bkg_eff[jj] = (TH1F*) h_bkg_post->Clone();
      h_bkg_eff[jj]->SetName("bkg_eff_"+channels[jj]);
      h_bkg_eff[jj]->GetYaxis()->SetTitle("Background Efficiency");
      h_bkg_eff[jj]->Divide(h_bkg_post,h_bkg_pre,1.0,1.0,"B");
      if (hists[ii] == "ak8jetPt") h_bkg_eff[jj]->GetXaxis()->SetRangeUser(400.,1200.);
      if (hists[ii] == "leadLepPt") h_bkg_eff[jj]->GetXaxis()->SetRangeUser(50.,500.);
    }//End channel loop
    
    // Plot
    TCanvas* c = new TCanvas("c","c",900,700);

    float xcor = 0.2;
    if (hists[ii] == "leadLepPt") xcor = 0.7;  
    TLegend* leg1 = new TLegend(xcor,0.15,xcor+0.2,0.45);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.045);
    leg1->AddEntry(h_signal_eff[0], "Eff_sig, mu", "l");
    leg1->AddEntry(h_signal_eff[1], "Eff_sig, el", "l");
    
    h_signal_eff[0]->SetLineColor(1);
    h_signal_eff[0]->GetYaxis()->SetRangeUser(0.5,1.0);
    h_signal_eff[0]->GetXaxis()->SetTitleOffset(1.1);
    h_signal_eff[0]->Draw("e");
  
    h_signal_eff[1]->SetLineColor(2);
    h_signal_eff[1]->Draw("e,same");
    leg1->Draw();
    c->SaveAs("Plots/eff_sig_final_"+hists[ii]+".pdf");

    TLegend* leg2 = new TLegend(xcor,0.58,xcor+0.2,0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.045);
    leg2->AddEntry(h_bkg_eff[0], "Eff_bkg, mu", "l");
    leg2->AddEntry(h_bkg_eff[1], "Eff_bkg, el", "l");

    h_bkg_eff[0]->SetLineColor(1);
    h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,0.5);
    h_bkg_eff[0]->GetXaxis()->SetTitleOffset(1.1);
    h_bkg_eff[0]->Draw("e");
  
    h_bkg_eff[1]->SetLineColor(2);
    h_bkg_eff[1]->Draw("e,same");
    leg2->Draw();
    c->SaveAs("Plots/eff_bkg_final_"+hists[ii]+".pdf");
  }
}

void drawROCCurve(TString channel){

  setStyle();
  
  SummedHist* Iso2DScanPoints_sig = getTTbar("histfiles_mMu_mEl_Loose/","2DisoScanPoints","",channel,false,"nom",false);
  SummedHist* Iso2DScanPoints_bkg = getQCDMC("histfiles_mMu_mEl_Loose/","2DisoScanPoints","",channel,false,"nom",false);
  SummedHist* MiniIsoScanPoints_sig = getTTbar("histfiles_mMu_mEl_Loose/","MiniIsoScanPoints","",channel,false,"nom",false);
  SummedHist* MiniIsoScanPoints_bkg = getQCDMC("histfiles_mMu_mEl_Loose/","MiniIsoScanPoints","",channel,false,"nom",false);
  SummedHist* Iso2DScanPoints_sig_tight = getTTbar("histfiles_tMu_tEl/","2DisoScanPoints","",channel,false,"nom",false);
  SummedHist* Iso2DScanPoints_bkg_tight = getQCDMC("histfiles_tMu_tEl/","2DisoScanPoints","",channel,false,"nom",false);
  SummedHist* MiniIsoScanPoints_sig_tight = getTTbar("histfiles_tMu_tEl/","MiniIsoScanPoints","",channel,false,"nom",false);
  SummedHist* MiniIsoScanPoints_bkg_tight = getQCDMC("histfiles_tMu_tEl/","MiniIsoScanPoints","",channel,false,"nom",false);
    
  TH1F* h_Iso2DScanPoints_sig = (TH1F*) Iso2DScanPoints_sig->hist();
  TH1F* h_Iso2DScanPoints_bkg = (TH1F*) Iso2DScanPoints_bkg->hist();
  TH1F* h_MiniIsoScanPoints_sig = (TH1F*) MiniIsoScanPoints_sig->hist();
  TH1F* h_MiniIsoScanPoints_bkg = (TH1F*) MiniIsoScanPoints_bkg->hist();
  TH1F* h_Iso2DScanPoints_sig_tight = (TH1F*) Iso2DScanPoints_sig_tight->hist();
  TH1F* h_Iso2DScanPoints_bkg_tight = (TH1F*) Iso2DScanPoints_bkg_tight->hist();
  TH1F* h_MiniIsoScanPoints_sig_tight = (TH1F*) MiniIsoScanPoints_sig_tight->hist();
  TH1F* h_MiniIsoScanPoints_bkg_tight = (TH1F*) MiniIsoScanPoints_bkg_tight->hist();
  
  const int n2D = 168;
  const int nMini = 30;
  
  float Sig_eff_2D[n2D];
  float Bkg_eff_2D[n2D];
  float Sig_eff_2D_tight[n2D];
  float Bkg_eff_2D_tight[n2D];
  float Mini_sig_eff[nMini];
  float Mini_bkg_eff[nMini];
  float Mini_sig_eff_tight[nMini];
  float Mini_bkg_eff_tight[nMini];
  
  for (int i=0; i<n2D; i++) {
    Sig_eff_2D[i] = h_Iso2DScanPoints_sig->GetBinContent(i+2) / h_Iso2DScanPoints_sig->GetBinContent(1);
    Bkg_eff_2D[i] = h_Iso2DScanPoints_bkg->GetBinContent(i+2) / h_Iso2DScanPoints_bkg->GetBinContent(1);
    Sig_eff_2D_tight[i] = h_Iso2DScanPoints_sig_tight->GetBinContent(i+2) / h_Iso2DScanPoints_sig_tight->GetBinContent(1);
    Bkg_eff_2D_tight[i] = h_Iso2DScanPoints_bkg_tight->GetBinContent(i+2) / h_Iso2DScanPoints_bkg_tight->GetBinContent(1);
  }
  
  for (int i=0; i<nMini; i++) {
    Mini_sig_eff[i] = h_MiniIsoScanPoints_sig->GetBinContent(i+2) / h_MiniIsoScanPoints_sig->GetBinContent(1);
    Mini_bkg_eff[i] = h_MiniIsoScanPoints_bkg->GetBinContent(i+2) / h_MiniIsoScanPoints_bkg->GetBinContent(1);
    Mini_sig_eff_tight[i] = h_MiniIsoScanPoints_sig_tight->GetBinContent(i+2) / h_MiniIsoScanPoints_sig_tight->GetBinContent(1);
    Mini_bkg_eff_tight[i] = h_MiniIsoScanPoints_bkg_tight->GetBinContent(i+2) / h_MiniIsoScanPoints_bkg_tight->GetBinContent(1);
  }
    
  TGraph* g_eff2D = new TGraph(n2D,Sig_eff_2D,Bkg_eff_2D);
  TGraph* g_effMini = new TGraph(nMini,Mini_sig_eff,Mini_bkg_eff);
  TGraph* g_eff2D_tight = new TGraph(n2D,Sig_eff_2D_tight,Bkg_eff_2D_tight);
  TGraph* g_effMini_tight = new TGraph(nMini,Mini_sig_eff_tight,Mini_bkg_eff_tight);

  TH2F* h_dummy_mu = new TH2F("dummy_mu", "; efficiency (ttbar); efficiency (QCD)",50,0.5,1.0,100,0.0,0.1);
  TH2F* h_dummy_el = new TH2F("dummy_el", "; efficiency (ttbar); efficiency (QCD)",50,0.5,1.0,100,0.0,1.0);
  
  g_eff2D->SetMarkerStyle(8);
  g_effMini->SetMarkerColor(2);
  g_effMini->SetLineColor(2);
  g_effMini->SetMarkerStyle(8);
  g_eff2D_tight->SetMarkerColor(3);
  g_eff2D_tight->SetLineColor(3);
  g_eff2D_tight->SetMarkerStyle(8);
  g_effMini_tight->SetMarkerColor(4);
  g_effMini_tight->SetLineColor(4);
  g_effMini_tight->SetMarkerStyle(8);
  
  TCanvas* c = new TCanvas("c","c",900,700);
  if (channel == "mu") h_dummy_mu->Draw();
  if (channel == "el") h_dummy_el->Draw();
  g_eff2D->Draw("same,ep");
  g_effMini->Draw("same,ep");
  g_eff2D_tight->Draw("same,ep");
  g_effMini_tight->Draw("same,ep");
  
  TLegend* lh = new TLegend(0.20,0.60,0.45,0.8);
  lh->SetFillStyle(0);
  lh->SetBorderSize(0);
  lh->SetTextSize(0.045);
  lh->AddEntry(g_eff2D,"2D iso, medium ID","lp");
  lh->AddEntry(g_effMini,"MiniIso, medium ID","lp");
  lh->AddEntry(g_eff2D_tight,"2D iso, tight ID","lp");
  lh->AddEntry(g_effMini_tight,"MiniIso, tight ID","lp");
  lh->SetTextFont(42);
  lh->Draw();	
  
  c->SaveAs("Plots/ROCCurve_"+channel+".pdf");
  
}

void checkElID(){
  //Compare signal and sideband shapes for electron pt, eta; AK4 eta
  TString histnames[3] = {"elPt","elEta","ak4jetEta"};
  TString sidebands[16] = {"1a","1b","1c","1d","2a","2b","2c","2d","3a","3b","3c","3d","4a","4b","4c","4d"};
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 4; jj++){
      SummedHist* signal = getQCDMC("histfiles_full2016/", histnames[ii]+"PassSel", "", "el", false, "nom", false, "", true);
      TH1F* h_sig = (TH1F*) signal->hist();
      SummedHist* sidebandA = getQCDMC("histfiles_full2016/", histnames[ii]+"PassQCD"+sidebands[4*jj+0], "", "el", false, "nom", false, "", true);
      TH1F* h_sideA = (TH1F*) sidebandA->hist();
      SummedHist* sidebandB = getQCDMC("histfiles_full2016/", histnames[ii]+"PassQCD"+sidebands[4*jj+1], "", "el", false, "nom", false, "", true);
      TH1F* h_sideB = (TH1F*) sidebandB->hist();
      SummedHist* sidebandC = getQCDMC("histfiles_full2016/", histnames[ii]+"PassQCD"+sidebands[4*jj+2], "", "el", false, "nom", false, "", true);
      TH1F* h_sideC = (TH1F*) sidebandC->hist();
      SummedHist* sidebandD = getQCDMC("histfiles_full2016/", histnames[ii]+"PassQCD"+sidebands[4*jj+3], "", "el", false, "nom", false, "", true);
      TH1F* h_sideD = (TH1F*) sidebandD->hist();
      
      cout << "*********************" << endl;
      cout << "Signal: " << h_sig->GetSum() << endl;
      cout << "Sideband " << sidebands[4*jj+0] << ": " << h_sideA->GetSum() << endl;
      cout << "Sideband " << sidebands[4*jj+1] << ": " << h_sideB->GetSum() << endl;
      cout << "Sideband " << sidebands[4*jj+2] << ": " << h_sideC->GetSum() << endl;
      cout << "Sideband " << sidebands[4*jj+3] << ": " << h_sideD->GetSum() << endl;
      cout << "*********************" << endl;
      
      h_sig->Rebin(5);
      h_sideA->Rebin(5);
      h_sideB->Rebin(5);
      h_sideC->Rebin(5);
      h_sideD->Rebin(5);
      
      h_sig->Scale(1.0/h_sig->Integral());
      h_sideA->Scale(1.0/h_sideA->Integral());
      h_sideB->Scale(1.0/h_sideB->Integral());
      h_sideC->Scale(1.0/h_sideC->Integral());
      h_sideD->Scale(1.0/h_sideD->Integral());
      
      TH1F* h_ratioA = (TH1F*) h_sig->Clone();
      h_ratioA->Divide(h_sideA);
      TH1F* h_ratioB = (TH1F*) h_sig->Clone();
      h_ratioB->Divide(h_sideB);
      TH1F* h_ratioC = (TH1F*) h_sig->Clone();
      h_ratioC->Divide(h_sideC);
      TH1F* h_ratioD = (TH1F*) h_sig->Clone();
      h_ratioD->Divide(h_sideD);
      
      // Plot
      TCanvas* c = new TCanvas("c_"+histnames[ii],"c_"+histnames[ii],900,800);
      TPad* p1 = new TPad("p1_"+histnames[ii],"p1_"+histnames[ii],0,0.3,1,1);
      p1->SetTopMargin(0.08);
      p1->SetBottomMargin(0.05);
      p1->SetNumber(1);
      TPad* p2 = new TPad("p2_"+histnames[ii],"p2_"+histnames[ii],0,0,1,0.3);
      p2->SetNumber(2);
      p2->SetTopMargin(0.05);
      p2->SetBottomMargin(0.35);
      
      p1->Draw();
      p2->Draw();
      p1->cd();
      
      h_ratioA->GetXaxis()->SetTitle(h_sig->GetXaxis()->GetTitle());
      h_sig->GetXaxis()->SetTitle("");
      
      h_sig->SetFillColor(0);
      h_sig->SetLineColor(1);
      h_sideA->SetFillColor(0);
      h_sideA->SetLineColor(kRed);
      h_sideB->SetFillColor(0);
      h_sideB->SetLineColor(kBlue);
      h_sideC->SetFillColor(0);
      h_sideC->SetLineColor(kGreen);
      h_sideD->SetFillColor(0);
      h_sideD->SetLineColor(kOrange);
      h_ratioA->SetFillColor(0);
      h_ratioA->SetLineColor(kRed);
      h_ratioB->SetFillColor(0);
      h_ratioB->SetLineColor(kBlue);
      h_ratioC->SetFillColor(0);
      h_ratioC->SetLineColor(kGreen);
      h_ratioD->SetFillColor(0);
      h_ratioD->SetLineColor(kOrange);
      
      h_sig->Draw("hist");
      h_sideA->Draw("hist,same");
      h_sideB->Draw("hist,same");
      h_sideC->Draw("hist,same");
      h_sideD->Draw("hist,same");
      
      TLegend* leg = new TLegend(0.45,0.6,0.75,0.8);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->AddEntry(h_sig, "Signal region", "l");
      leg->AddEntry(h_sideA, "Sideband "+sidebands[4*jj+0], "l");
      leg->AddEntry(h_sideB, "Sideband "+sidebands[4*jj+1], "l");
      leg->AddEntry(h_sideC, "Sideband "+sidebands[4*jj+2], "l");
      leg->AddEntry(h_sideD, "Sideband "+sidebands[4*jj+3], "l");
      leg->Draw();
      
      p2->cd();
      h_ratioA->Draw("hist");
      h_ratioB->Draw("hist,same");
      h_ratioC->Draw("hist,same");
      h_ratioD->Draw("hist,same");
      h_ratioA->SetMaximum(2.0);
      h_ratioA->SetMinimum(0.0);
      h_ratioA->GetYaxis()->SetNdivisions(2,4,0,false);
      h_ratioA->GetYaxis()->SetTitle("Signal / Sideband");
      h_ratioA->GetXaxis()->SetLabelSize(0.1);
      h_ratioA->GetYaxis()->SetLabelSize(0.1);
      h_ratioA->GetXaxis()->SetTitleOffset(1.0);
      h_ratioA->GetYaxis()->SetTitleOffset(0.5);
      h_ratioA->GetXaxis()->SetTitleSize(0.12);
      h_ratioA->GetYaxis()->SetTitleSize(0.1);

      TString itos = Form("%d",jj+1);
      c->SaveAs("Plots/compQCD_elID_"+histnames[ii]+"_side"+itos+".pdf");
    }
  }
}

void compLepQCD(TString DIR, TString var, bool isSide) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();
  
  TString hist = var+"Pre";
  
  // get histograms
  SummedHist* muQCD = getQCDMC(DIR, var, "Pre", "mu", isSide, "nom", false);
  SummedHist* elQCD = getQCDMC(DIR, var, "Pre", "el", isSide, "nom", false);
  TH1F* h_muQCD = (TH1F*) muQCD->hist();
  TH1F* h_elQCD = (TH1F*) elQCD->hist();
  
  h_muQCD->SetLineColor(2);
  h_muQCD->SetMarkerColor(2);
  h_elQCD->SetLineColor(4);
  h_elQCD->SetMarkerColor(4);
  
  if (hist.Contains("ak8jetSDsubjetMaxCSV") || hist.Contains("ak4jetEta") || hist.Contains("lepAbsEta")){
    h_muQCD->Rebin(5);
    h_elQCD->Rebin(5);
  }
  
  else if (var != "nAK4jet" && var != "nAK8jet" && var != "nBjet" && var != "nTjet"){
    h_muQCD->Rebin(2);
    h_elQCD->Rebin(2);
  }
  
  h_muQCD->Scale(1./h_muQCD->Integral());
  h_elQCD->Scale(1./h_elQCD->Integral());

  TH1F* h_ratio = (TH1F*) h_muQCD->Clone();
  TH1F* h_ratio2 = (TH1F*) h_elQCD->Clone();
  h_ratio->Divide(h_muQCD);
  h_ratio2->Divide(h_muQCD);
  
  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
  p1->SetTopMargin(0.08);
  p1->SetBottomMargin(0.05);
  p1->SetNumber(1);
  TPad* p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
  p2->SetNumber(2);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.35);

  p1->Draw();
  p2->Draw();
  p1->cd();

  h_ratio->GetXaxis()->SetTitle(h_muQCD->GetXaxis()->GetTitle());

  h_muQCD->Draw();
  h_elQCD->Draw("same");

  float xmin = 0.66;
  float xmax = 0.85;
  float ymin = 0.65;
  float ymax = 0.88;

  if (var.Contains("AbsEta")) {
    xmin = 0.2;
    xmax = 0.4;
  }
  
  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmax,ymax);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_elQCD, "QCD, e+jets", "l");
  leg->AddEntry(h_muQCD, "QCD, #mu+jets", "l");
  leg->Draw();

  // plot ratio part
  p2->cd();
  h_ratio->Draw();
  h_ratio2->Draw("same");
  h_ratio->SetMaximum(2.0);
  h_ratio->SetMinimum(0.0);
  h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
  h_ratio->GetYaxis()->SetTitle("MC / Data-driven");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleOffset(1.0);
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleSize(0.1);

  // save output
  TString outname = "Plots/compLepQCD_"+hist;
  if (isSide) outname += "_side.pdf";
  else outname += ".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_elQCD->Delete();
  h_muQCD->Delete();
  h_ratio->Delete();
  h_ratio2->Delete();
}



