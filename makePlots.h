#ifndef makePlots_h
#define makePlots_h

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <THStack.h>
#include <TColor.h>
#include <TFile.h>
#include <TROOT.h>
#include <Rtypes.h>
#include <vector>
#include <iostream>


// -------------------------------------------------------------------------------------
// various declarations
// -------------------------------------------------------------------------------------

const double LUM[2] = {35867.0,35867.0}; //pb-1
// TODO: use correct lumis for B2G datasets (rather than lumi for total golden 2016 dataset)

// -------------------------------------------------------------------------------------
// helper class for summed, weighted histograms (e.g. single top)
// -------------------------------------------------------------------------------------

class SummedHist {
 public : 
 SummedHist( TString const & name, int color ) : name_(name), color_(color) {
    summedHist_ = 0;
  };
  
  // return the summed histogram if created, else create it (summing up histograms in vector hists_) 
  TH1* hist(bool noFill = false) { 
    if (summedHist_ != 0) {
      if (noFill) summedHist_->SetFillColor(0);
      return summedHist_; 
    }
    else if (hists_.size() == 0) {
      return 0; 
    } 
    else {
      summedHist_ = (TH1*)hists_[0]->Clone();
      summedHist_->SetName( name_ );
      if (!noFill) summedHist_->SetFillColor( color_ );
      for (unsigned int j = 1; j<hists_.size(); ++j) {
	if (hists_[j]->Integral() > 0.0) summedHist_->Add( hists_[j], 1.0 );
      }
      return summedHist_; 
    };
  }
  
  // return the vector of input histograms  
  std::vector<TH1*> const & hists() const {
    return hists_;
  }
  
  // add histogram to the vector of inputs
  void push_back( TH1 const * ihist, double norm ) {
    TH1* clone = (TH1*) ihist->Clone();
    TString iname( name_ );
    iname += hists_.size();
    clone->SetName( iname );
    if (clone->Integral() > 0.0) clone->Scale( norm );
    hists_.push_back( clone );
    norms_.push_back( norm );
  };
  

 protected : 
  
  std::vector<TH1*> hists_;
  std::vector<double> norms_;
  
  TString name_; 
  int color_;
  TH1* summedHist_; 
  
};

TH1* getHist(TString filename, TString histname, TString region, TString split = ""){
  TH1::AddDirectory(kFALSE); 

  TString append = "";
  if (split != "") append = "_"+split;

  TFile* infile = TFile::Open( filename );
  TH1* hist;
  if (histname == "counts"){
    Double_t err_pre, err_1b, err_1t, err_1t1b;
    TH1F* h_pre = (TH1F*) infile->Get("ak4jetEtaPre"+append);
    TH1F* h_1t = (TH1F*) infile->Get("ak4jetEta1t"+append);
    TH1F* h_1t1b = (TH1F*) infile->Get("ak4jetEta1t1b"+append);
    hist = new TH1F("count","",3,-0.5,2.5);
    hist->SetBinContent(1,h_pre->IntegralAndError(0,h_pre->GetNbinsX()+1,err_pre) - h_1t->IntegralAndError(0,h_1t->GetNbinsX()+1,err_1t));
    hist->SetBinError(1,sqrt(pow(err_pre,2)+pow(err_1t,2)));
    hist->SetBinContent(2,h_1t->IntegralAndError(0,h_1t->GetNbinsX()+1,err_1t) - h_1t1b->IntegralAndError(0,h_1t1b->GetNbinsX()+1,err_1t1b));
    hist->SetBinError(2,sqrt(pow(err_1t,2)+pow(err_1t1b,2)));
    hist->SetBinContent(3,h_1t1b->IntegralAndError(0,h_1t1b->GetNbinsX()+1,err_1t1b));
    hist->SetBinError(3,err_1t1b);
    delete h_pre;
    delete h_1t;
    delete h_1t1b;
  }
  else {
    if (region == "1t1b" || region == "1t" || region == "1b" || region == "Pre" || region == "") {
      hist = (TH1*) infile->Get(histname+region+append);
    }
    if (region == "1t0b"){
      hist = (TH1*) infile->Get(histname+"1t"+append);
      TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b"+append);
      hist->Add(h_tmp,-1.0);
      delete h_tmp;
    }
    if (region == "0t1b"){
      hist = (TH1*) infile->Get(histname+"1b"+append);
      TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b"+append);
      hist->Add(h_tmp,-1.0);
      delete h_tmp;
    }
    if (region == "0t"){
      hist = (TH1*) infile->Get(histname+"Pre"+append);
      TH1* h_tmp = (TH1*) infile->Get(histname+"1t"+append);
      hist->Add(h_tmp,-1.0);
      delete h_tmp;
    }
    if (region == "0b"){
      hist = (TH1*) infile->Get(histname+"Pre"+append);
      TH1* h_tmp = (TH1*) infile->Get(histname+"1b"+append);
      hist->Add(h_tmp,-1.0);
      delete h_tmp;
    }
    if (region == "0t0b"){
      hist = (TH1*) infile->Get(histname+"Pre"+append);
      TH1* h_tmp1 = (TH1*) infile->Get(histname+"1t"+append);
      TH1* h_tmp2 = (TH1*) infile->Get(histname+"1b"+append);
      TH1* h_tmp3 = (TH1*) infile->Get(histname+"1t1b"+append);
      hist->Add(h_tmp3);
      hist->Add(h_tmp2,-1.0);
      hist->Add(h_tmp1,-1.0);
      delete h_tmp1;
      delete h_tmp2;
      delete h_tmp3;
    }
  }
  delete infile;
  return hist;  
}

// -------------------------------------------------------------------------------------
// W+jets
// -------------------------------------------------------------------------------------

SummedHist * getWJets( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  const int nwjets = 7;
  
  TString wjets_names[nwjets] = {
    "WJets_HT100to200",
    "WJets_HT200to400",
    "WJets_HT400to600",
    "WJets_HT600to800",
    "WJets_HT800to1200",
    "WJets_HT1200to2500",
    "WJets_HT2500toInf",
  };
  
  double wjets_norms[nwjets] = {
    1345.0 * 1.21 * LUM[ich] / 39617787.,  // Note: cross sections are from AN-15-107, may need updating
    359.7 * 1.21 * LUM[ich] / 19914590.,  
    48.91 * 1.21 * LUM[ich] / 5796237.,  
    12.05 * 1.21 * LUM[ich] / 14822888.,
    5.501 * 1.21 * LUM[ich] / 6200954.,
    1.329 * 1.21 * LUM[ich] / 6324934.,
    0.03216 * 1.21 * LUM[ich] / 2384260.,
  };

  int plotcolor = kGreen-3;
  if (split == "bb") plotcolor = kRed+1;
  if (split == "b") plotcolor = kRed-7;
  if (split == "cc") plotcolor = 6;
  if (split == "c") plotcolor = kGreen-3;
  if (split == "l") plotcolor = kYellow;
  
  SummedHist* wjets = new SummedHist( histname, plotcolor );
  
  for (int i=0 ; i<nwjets; i++) {
    TString iname = DIR + "hists_" + wjets_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
    wjets->push_back( hist, wjets_norms[i] );
  }
  
  return wjets;
  
}

// -------------------------------------------------------------------------------------
// Diboson (WW, WZ, ZZ)
// -------------------------------------------------------------------------------------

SummedHist * getDiboson( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  const int ndiboson = 3;
  
  TString diboson_names[ndiboson] = {
    "WW",
    "WZ",
    "ZZ",
  };
  
  double diboson_norms[ndiboson] = {
    118.7 * LUM[ich] / 6987124., //Xsec / Nevents from Louise -- source?
    44.9 * LUM[ich] / 2995828.,
    15.4 * LUM[ich] / 990064.,
  };

  int plotcolor = kViolet-6;
  
  SummedHist* diboson = new SummedHist( histname, plotcolor );
  
  for (int i=0 ; i<ndiboson; i++) {
    TString iname = DIR + "hists_" + diboson_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
    diboson->push_back( hist, diboson_norms[i] );
  }
  
  return diboson;
  
}

// -------------------------------------------------------------------------------------
// ZJets
// -------------------------------------------------------------------------------------

SummedHist * getZJets( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  double zjets_norm = 5765 * LUM[ich] / 42923575.; //Nevents from Louise, xsec from AN-17-003

  int plotcolor = kAzure-9;
  
  SummedHist* zjets = new SummedHist( histname, plotcolor );
  TString iname = DIR + "hists_ZJets_" + channel + "_" + syst + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region,split);
  zjets->push_back( hist, zjets_norm);
  
  return zjets;
  
}

// -------------------------------------------------------------------------------------
// single top
// -------------------------------------------------------------------------------------

SummedHist * getSingleTop( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nsingletop = 5;
  
  TString singletop_names[nsingletop] = { 
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
    "SingleTop_s",
  };
  
  double singletop_norms[nsingletop] = {
    136.02 * LUM[ich] / 67240808., //Event counts from Louise, xsec from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    80.95 * LUM[ich] / 38811017.,  
    19.3 * LUM[ich] / 8629641.,    // Xsec from AN-16-020 -- NoFullyHadronicDecays
    19.3 * LUM[ich] / 8681541.,
    10.32 * 0.322 * LUM[ich] / 9651642. // BR is needed because s-channel sample is leptonic final state only
  };
  
  SummedHist* singletop = new SummedHist( histname, 6 );
  
  for (int i=0; i<nsingletop; i++) {
    TString iname = DIR + "hists_" + singletop_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
    singletop->push_back( hist, singletop_norms[i] );
  }
  
  return singletop;
  
}


// -------------------------------------------------------------------------------------
// single top, tW channel only
// -------------------------------------------------------------------------------------

SummedHist * getSingleTop_tW( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nsingletop_tW = 2;
  
  TString singletop_tW_names[nsingletop_tW] = { 
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
  };
  
  double singletop_tW_norms[nsingletop_tW] = {
    19.3 * LUM[ich] / 8629641., //Event counts from Louise, xsec from AN-16-020 -- NoFullyHadronicDecays
    19.3 * LUM[ich] / 8681541.,
  };
  
  SummedHist* singletop_tW = new SummedHist( histname, 6 );
  
  for (int i=0; i<nsingletop_tW; i++) {
    TString iname = DIR + "hists_" + singletop_tW_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
    singletop_tW->push_back( hist, singletop_tW_norms[i] );
  }
  
  return singletop_tW;
  
}


// -------------------------------------------------------------------------------------
// non-semileptonic ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbarNonSemiLep( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "" ) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  TString ttbar_name = "PowhegPythia8_nonsemilep";
  double ttbar_norm = 831.76 * LUM[ich] / 77229341.; //TODO: technically this is incorrect, as there is a job missing -- but can't tell what files the job corresponds to
  
  SummedHist* ttbar = new SummedHist( histname, kRed-7);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region,split);
  ttbar->push_back( hist, ttbar_norm );
  
  return ttbar;
  
}


// -------------------------------------------------------------------------------------
// signal ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbar( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {
  
  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  TString ttbar_name = (usePost && !isQCD) ? "PowhegPythia8_fullTruth_mInc" : "PowhegPythia8_fullTruth";
  double ttbar_norm = 831.76 * LUM[ich] / 77229341.;
  
  SummedHist* ttbar = new SummedHist( histname, kRed+1);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region,split);
  ttbar->push_back( hist, ttbar_norm );
  
  return ttbar;
  
}

// -------------------------------------------------------------------------------------
// QCD
// -------------------------------------------------------------------------------------

SummedHist * getQCDMC( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "", bool elID = false) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nqcd = 5;
  
  TString qcd_names[nqcd] = {
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
  };
  
  double qcd_norms[nqcd] = {
    32100. * LUM[ich] / 18929951.,  // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    6831. * LUM[ich] / 15629253.,   // Technically part of this is missing, but can't tell which part -- ~1% though
    1207. * LUM[ich] / 4767100., 
    119.9 * LUM[ich] / 3970819.,
    25.24 * LUM[ich] / 1991645.,
  };
  
  SummedHist* qcd = new SummedHist( histname, kYellow );
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    if (elID) iname = DIR + "hists_" + qcd_names[i] + "_elID.root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
    qcd->push_back( hist, qcd_norms[i] );
  }
  
  return qcd;
  
}

SummedHist * getData( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString split = "") {
  
  TString append = "";
  if (isQCD) append = "_qcd";
    
  SummedHist* data = new SummedHist( histname, 0 );
  
  TString iname = DIR + "hists_Data_" + channel + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region,split);
  data->push_back( hist, 1.0 );
  
  return data;
  
}

TH1 * getQCDData(TString sigDIR, TString sideDIR, TString histname, TString region, TString channel, TString syst, bool usePost, TString split = "") {

  if (histname == "counts") {
    SummedHist* qcd = getQCDMC(sigDIR,histname,region,channel,false,syst,usePost);
    return qcd->hist();
  }

  else {
    SummedHist* diboson = getDiboson( sideDIR, histname, region, channel, true, syst, usePost,split);
    SummedHist* zjets = getZJets( sideDIR, histname, region, channel, true, syst, usePost,split);
    SummedHist* wjets = getWJets( sideDIR, histname, region, channel, true, syst, usePost,split);
    SummedHist* singletop = getSingleTop( sideDIR, histname, region, channel, true, syst, usePost,split );
    SummedHist* ttbar = getTTbar( sideDIR, histname, region, channel, true, syst, usePost,split );
    SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( sideDIR, histname, region, channel, true, syst, usePost,split );
    SummedHist* data = getData( sideDIR, histname, region, channel, true,split);
    
    TH1* h_diboson = (TH1*) diboson->hist();
    TH1* h_zjets = (TH1*) zjets->hist();
    TH1* h_wjets = (TH1*) wjets->hist();
    TH1* h_singletop = (TH1*) singletop->hist();
    TH1* h_ttbar = (TH1*) ttbar->hist();
    TH1* h_ttbar_nonSemiLep = (TH1*) ttbar_nonSemiLep->hist();
    TH1* h_data = (TH1*) data->hist();
    
    TH1* h_qcd = (TH1*) h_data->Clone("QCD");
    h_qcd->Add(h_diboson,-1.0);
    h_qcd->Add(h_zjets,-1.0);
    h_qcd->Add(h_wjets,-1.0);
    h_qcd->Add(h_singletop,-1.0);
    h_qcd->Add(h_ttbar,-1.0);
    h_qcd->Add(h_ttbar_nonSemiLep,-1.0);
      
    for (int ii = 1; ii < h_qcd->GetNbinsX()+1; ii++){
      if (h_qcd->GetBinContent(ii) < 0.0) h_qcd->SetBinContent(ii,0.0);
    }
    SummedHist* qcd_mc = getQCDMC(sigDIR,histname,region,channel,false,syst,usePost,split);
    TH1F* h_qcd_mc = (TH1F*) qcd_mc->hist();
    float n_qcd = 0.0;
    n_qcd = h_qcd_mc->GetSum();
    if (h_qcd->Integral() > 0.0) h_qcd->Scale(n_qcd / h_qcd->Integral());
    h_qcd->SetFillColor(kYellow);
    
    //wjets->Delete();
    //singletop->Delete();
    //ttbar->Delete();
    //ttbar_nonSemiLep->Delete();
    //data->Delete();
    
    return h_qcd;
  }
}

TObject * getBackground( TString DIR, TString histname, TString region, TString channel, bool isQCD ) {

  TString append = "";
  if (isQCD) append = "_qcd";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  const int nbkg = 22;
  
  TString bkg_names[nbkg] = {
    "PowhegPythia8_nonsemilep",
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
    "SingleTop_s",
    "WJets_HT100to200",
    "WJets_HT200to400",
    "WJets_HT400to600",
    "WJets_HT600to800",
    "WJets_HT800to1200",
    "WJets_HT1200to2500",
    "WJets_HT2500toInf",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
    "WW",
    "WZ",
    "ZZ",
    "ZJets"
  };
  
  double bkg_norms[nbkg] = {
    831.76 * LUM[ich] / 77229341.,       //TTbar nonsignal
    136.02 * LUM[ich] / 67240808., //SingleTop
    80.95 * LUM[ich] / 38811017.,  
    19.3 * LUM[ich] / 8629641.,    
    19.3 * LUM[ich] / 8681541.,
    10.32 * 0.322 * LUM[ich] / 9651642.,
    1345.0 * 1.21 * LUM[ich] / 39617787., // WJets
    359.7 * 1.21 * LUM[ich] / 19914590.,  
    48.91 * 1.21 * LUM[ich] / 5796237.,  
    12.05 * 1.21 * LUM[ich] / 14822888.,
    5.501 * 1.21 * LUM[ich] / 6200954.,
    1.329 * 1.21 * LUM[ich] / 6324934.,
    0.03216 * 1.21 * LUM[ich] / 2384260.,
    32100. * LUM[ich] / 18929951., //QCD 
    6831. * LUM[ich] / 15629253.,  
    1207. * LUM[ich] / 4767100., 
    119.9 * LUM[ich] / 3970819.,
    25.24 * LUM[ich] / 1991645.,
    118.7 * LUM[ich] / 6987124., //WW
    44.9 * LUM[ich] / 2995828.,  //WZ
    15.4 * LUM[ich] / 990064.,   //ZZ
    5765 * LUM[ich] / 42923575., //ZJets
  };
  
  TH2F* bkg_hists[nbkg];
  
  for (int i=0; i<nbkg; i++) {
    TString iname = DIR + "hists_" + bkg_names[i] + "_" + channel + "_nom" + append + ".root";
    bkg_hists[i] = (TH2F*) getHist(iname,histname,region);
    if (bkg_hists[i]->Integral() > 0.0) bkg_hists[i]->Scale(bkg_norms[i]);
  }
  
  TH2F* bkg = (TH2F*) bkg_hists[0]->Clone();
  for (int j=1; j<nbkg; j++){
    bkg->Add(bkg_hists[j]);
  }

  return bkg;
  
}

TH1F* adjustRange(TH1F* h_input, float xlow, float xhigh){
  if (!h_input) return 0;
  const int xbins = h_input->GetXaxis()->GetNbins();
  double edges[xbins];
  h_input->GetXaxis()->GetLowEdge(edges);
  int lowbin = 1;
  int highbin = xbins;
  for (int ii = 0; ii < xbins; ii++){
    if (edges[ii] <= xlow) lowbin = ii+1;
    if (edges[xbins-1-ii] >= xhigh) highbin = xbins-ii-1;
  }
  const int xbins_new = highbin - lowbin + 1;
  TString xlabel = h_input->GetXaxis()->GetTitle();
  TString ylabel = h_input->GetYaxis()->GetTitle();
  TH1F* h_output = new TH1F(h_input->GetName(),";"+xlabel+";"+ylabel,xbins_new,edges[lowbin-1],edges[highbin-1]+h_input->GetBinWidth(1));
  for (int ii = 0; ii < xbins_new; ii++){
    h_output->SetBinContent(ii+1,max(0.0001,h_input->GetBinContent(lowbin+ii)));
    h_output->SetBinError(ii+1,max(0.0001,h_input->GetBinError(lowbin+ii)));
  }
  return h_output;
}


#endif
