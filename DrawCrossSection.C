#include "CMS_lumi.C"
#include "tdrStyle.C"

void DrawCrossSection(TString LEVEL,TString VAR,TString XTITLE,TString YTITLE,TString YTITLE_NORM,float YMIN,float YMAX,float RMIN,float RMAX,float YMIN_NORM,float YMAX_NORM,float RMIN_NORM,float RMAX_NORM,bool LOGY)
{
  setTDRStyle();

  gROOT->ForceStyle();
  
  const int NTHEORY = 2;
  
  /*
  TString THEORY_ALIAS[NTHEORY] = {"PowhegPythia8","amcatnloPythia8","PowhegHerwigpp"};
  int THEORY_COLOR[NTHEORY]  = {kRed,kBlue,kGreen+2};
  int THEORY_MSTYLE[NTHEORY] = {20,25,26};
  int THEORY_LSTYLE[NTHEORY] = {1,1,1};
  int THEORY_FSTYLE[NTHEORY] = {3004,3003,3005};
  */
  TString THEORY_ALIAS[NTHEORY] = {"PowhegPythia8","PowhegHerwigpp"};
  int THEORY_COLOR[NTHEORY]  = {kRed,kGreen+2};
  int THEORY_MSTYLE[NTHEORY] = {20,26};
  int THEORY_LSTYLE[NTHEORY] = {1,1};
  int THEORY_FSTYLE[NTHEORY] = {3004,3005};
  
  TFile *inf;
  TH1F *hTheory[NTHEORY],*hTheoryNorm[NTHEORY],*hTheoryHIST[NTHEORY],*hTheoryNormHIST[NTHEORY];
  TH1F *hTheoryRatio[NTHEORY],*hTheoryRatioNorm[NTHEORY],*hTheoryRatioHIST[NTHEORY],*hTheoryRatioNormHIST[NTHEORY];
  TH1F *hSignal,*hSignalNorm,*hSignalRatio,*hSignalRatioNorm,*hRelUncMinus,*hRelUncPlus,*hRelUncMinusNorm,*hRelUncPlusNorm;

  if ((LEVEL == "Parton") || (LEVEL == "Gen") || (LEVEL == "Particle")) {
    inf = TFile::Open("CrossSection_"+LEVEL+"_"+VAR+".root");
    hSignal          = (TH1F*)inf->Get("CrossSection_"+LEVEL+"_Nominal");
    hSignalNorm      = (TH1F*)inf->Get("NormCrossSection_"+LEVEL+"_Nominal");
    hSignalRatio     = (TH1F*)inf->Get("CrossSectionRatio");
    hSignalRatioNorm = (TH1F*)inf->Get("NormCrossSectionRatio");
    hRelUncMinus     = (TH1F*)inf->Get("RelUncMinus");
    hRelUncPlus      = (TH1F*)inf->Get("RelUncPlus");
    hRelUncMinusNorm = (TH1F*)inf->Get("RelUncNormMinus");
    hRelUncPlusNorm  = (TH1F*)inf->Get("RelUncNormPlus");
  } 
  else {
    inf = TFile::Open("CrossSection_Parton_"+VAR+".root");
    hSignal          = (TH1F*)inf->Get("CrossSection_Fiducial_Nominal");
    hSignalNorm      = (TH1F*)inf->Get("NormCrossSection_Fiducial_Nominal");
    hSignalRatio     = (TH1F*)inf->Get("FiducialCrossSectionRatio");
    hSignalRatioNorm = (TH1F*)inf->Get("FiducialNormCrossSectionRatio");
    hRelUncMinus     = (TH1F*)inf->Get("FiducialRelUncMinus");
    hRelUncPlus      = (TH1F*)inf->Get("FiducialRelUncPlus");
    hRelUncMinusNorm = (TH1F*)inf->Get("FiducialRelUncNormMinus");
    hRelUncPlusNorm  = (TH1F*)inf->Get("FiducialRelUncNormPlus");
  }  

  TH1F *hTotalRelUncMinus     = (TH1F*)hRelUncMinus->Clone("TotalRelUncMinus");
  TH1F *hTotalRelUncMinusNorm = (TH1F*)hRelUncMinusNorm->Clone("TotalRelUncMinusNorm");
  TH1F *hTotalRelUncPlus      = (TH1F*)hRelUncPlus->Clone("TotalRelUncPlus");
  TH1F *hTotalRelUncPlusNorm  = (TH1F*)hRelUncPlusNorm->Clone("TotalRelUncPlusNorm");

  hTotalRelUncMinus->SetFillColor(kGray);
  hTotalRelUncMinusNorm->SetFillColor(kGray);
  hTotalRelUncPlus->SetFillColor(kGray);
  hTotalRelUncPlusNorm->SetFillColor(kGray);
  float vx[20],vexl[20],vexh[20],vy[20],veyl[20],veyh[20];

  for(int i=0;i<hTotalRelUncMinus->GetNbinsX();i++) {
    float eStat = hSignal->GetBinError(i+1);
    float y = hSignal->GetBinContent(i+1);
    float erSystMinus = hTotalRelUncMinus->GetBinContent(i+1);
    float erSystPlus = hTotalRelUncPlus->GetBinContent(i+1);
    float erTotMinus = sqrt(TMath::Power(eStat/y,2)+TMath::Power(erSystMinus,2));
    float erTotPlus = sqrt(TMath::Power(eStat/y,2)+TMath::Power(erSystPlus,2));
    hTotalRelUncMinus->SetBinContent(i+1,-erTotMinus);
    hTotalRelUncPlus->SetBinContent(i+1,erTotPlus);

    vx[i]   = hSignal->GetBinCenter(i+1);
    vexl[i] = 0.5*hSignal->GetBinWidth(i+1);
    vexh[i] = 0.5*hSignal->GetBinWidth(i+1);
    vy[i]   = y;
    veyl[i] = y*erTotMinus;
    veyh[i] = y*erTotPlus;
  }

  TGraphAsymmErrors *gSignalWithUnc = new TGraphAsymmErrors(hSignal->GetNbinsX(),vx,vy,vexl,vexh,veyl,veyh);
  gSignalWithUnc->SetFillColor(kGray);
  gSignalWithUnc->SetMarkerSize(0);

  for(int i=0;i<hTotalRelUncMinusNorm->GetNbinsX();i++) {
    float eStat = hSignalNorm->GetBinError(i+1);
    float y = hSignalNorm->GetBinContent(i+1);
    float erSystMinus = hTotalRelUncMinusNorm->GetBinContent(i+1);
    float erSystPlus = hTotalRelUncPlusNorm->GetBinContent(i+1);
    float erTotMinus = sqrt(TMath::Power(eStat/y,2)+TMath::Power(erSystMinus,2));
    float erTotPlus = sqrt(TMath::Power(eStat/y,2)+TMath::Power(erSystPlus,2));
    hTotalRelUncMinusNorm->SetBinContent(i+1,-erTotMinus);
    hTotalRelUncPlusNorm->SetBinContent(i+1,erTotPlus);

    vx[i]   = hSignalNorm->GetBinCenter(i+1);
    vexl[i] = 0.5*hSignalNorm->GetBinWidth(i+1);
    vexh[i] = 0.5*hSignalNorm->GetBinWidth(i+1);
    vy[i]   = y;
    veyl[i] = y*erTotMinus;
    veyh[i] = y*erTotPlus;
  }

  TGraphAsymmErrors *gSignalWithUncNorm = new TGraphAsymmErrors(hSignalNorm->GetNbinsX(),vx,vy,vexl,vexh,veyl,veyh);
  gSignalWithUncNorm->SetFillColor(kGray);
  gSignalWithUncNorm->SetMarkerSize(0);

  for(int ith=0;ith<NTHEORY;ith++) {
    if ((LEVEL == "Parton") || (LEVEL == "Gen") || (LEVEL == "Particle")) {
      hTheory[ith] = (TH1F*)inf->Get("CrossSection_"+THEORY_ALIAS[ith]);
      hTheoryNorm[ith] = (TH1F*)inf->Get("NormCrossSection_"+THEORY_ALIAS[ith]);
      hTheoryRatio[ith] = (TH1F*)inf->Get("RatioOverData_"+THEORY_ALIAS[ith]);
      hTheoryRatioNorm[ith] = (TH1F*)inf->Get("NormRatioOverData_"+THEORY_ALIAS[ith]);
    }
    else {
      hTheory[ith] = (TH1F*)inf->Get("Fiducial_CrossSection_"+THEORY_ALIAS[ith]);
      hTheoryNorm[ith] = (TH1F*)inf->Get("Fiducial_NormCrossSection_"+THEORY_ALIAS[ith]);
      hTheoryRatio[ith] = (TH1F*)inf->Get("Fiducial_RatioOverData_"+THEORY_ALIAS[ith]);
      hTheoryRatioNorm[ith] = (TH1F*)inf->Get("Fiducial_NormRatioOverData_"+THEORY_ALIAS[ith]);
    }
    hTheory[ith]->SetLineColor(THEORY_COLOR[ith]);
    hTheoryNorm[ith]->SetLineColor(THEORY_COLOR[ith]);
    hTheoryRatio[ith]->SetLineColor(THEORY_COLOR[ith]);
    hTheoryRatioNorm[ith]->SetLineColor(THEORY_COLOR[ith]);
    hTheory[ith]->SetLineWidth(2);
    hTheoryNorm[ith]->SetLineWidth(2);
    hTheoryRatio[ith]->SetLineWidth(2);
    hTheoryRatioNorm[ith]->SetLineWidth(2);
    hTheory[ith]->SetLineStyle(THEORY_LSTYLE[ith]);
    hTheoryNorm[ith]->SetLineStyle(THEORY_LSTYLE[ith]);
    hTheoryRatio[ith]->SetLineStyle(THEORY_LSTYLE[ith]);
    hTheoryRatioNorm[ith]->SetLineStyle(THEORY_LSTYLE[ith]);
    hTheory[ith]->SetMarkerColor(THEORY_COLOR[ith]);
    hTheoryNorm[ith]->SetMarkerColor(THEORY_COLOR[ith]);
    hTheoryRatio[ith]->SetMarkerColor(THEORY_COLOR[ith]);
    hTheoryRatioNorm[ith]->SetMarkerColor(THEORY_COLOR[ith]);
    hTheory[ith]->SetFillColor(THEORY_COLOR[ith]);
    hTheoryNorm[ith]->SetFillColor(THEORY_COLOR[ith]);
    hTheoryRatio[ith]->SetFillColor(THEORY_COLOR[ith]);
    hTheoryRatioNorm[ith]->SetFillColor(THEORY_COLOR[ith]);
    hTheory[ith]->SetFillStyle(THEORY_FSTYLE[ith]);
    hTheoryNorm[ith]->SetFillStyle(THEORY_FSTYLE[ith]);
    hTheoryRatio[ith]->SetFillStyle(THEORY_FSTYLE[ith]);
    hTheoryRatioNorm[ith]->SetFillStyle(THEORY_FSTYLE[ith]);
    hTheory[ith]->SetMarkerSize(0);
    hTheoryNorm[ith]->SetMarkerSize(0);
    hTheoryRatio[ith]->SetMarkerSize(0);
    hTheoryRatioNorm[ith]->SetMarkerSize(0);
    hTheory[ith]->SetMarkerStyle(THEORY_MSTYLE[ith]);
    hTheoryNorm[ith]->SetMarkerStyle(THEORY_MSTYLE[ith]);
    hTheoryRatio[ith]->SetMarkerStyle(THEORY_MSTYLE[ith]);
    hTheoryRatioNorm[ith]->SetMarkerStyle(THEORY_MSTYLE[ith]);

    hTheoryHIST[ith] = (TH1F*)hTheory[ith]->Clone(TString(hTheory[ith]->GetName())+"HIST");
    hTheoryNormHIST[ith] = (TH1F*)hTheoryNorm[ith]->Clone(TString(hTheoryNorm[ith]->GetName())+"HIST");
    hTheoryRatioHIST[ith] = (TH1F*)hTheoryRatio[ith]->Clone(TString(hTheoryRatio[ith]->GetName())+"HIST");
    hTheoryRatioNormHIST[ith] = (TH1F*)hTheoryRatioNorm[ith]->Clone(TString(hTheoryRatioNorm[ith]->GetName())+"HIST");
     
    hTheoryHIST[ith]->SetFillStyle(0);
    hTheoryNormHIST[ith]->SetFillStyle(0);
    hTheoryRatioHIST[ith]->SetFillStyle(0);
    hTheoryRatioNormHIST[ith]->SetFillStyle(0);
  }

  TLegend *leg = new TLegend(0.55,0.68,0.9,0.92);
  if (LEVEL == "Parton") {
    leg->SetHeader("Parton phase space");
  }
  else if (LEVEL == "Gen" || LEVEL == "Particle"){
    leg->SetHeader("Particle phase space");
  }
  else {
    leg->SetHeader("Fiducial phase space");
  }
  leg->AddEntry(hSignal,"Data","PEL");
  leg->AddEntry(hTotalRelUncMinus,"Total unc.","F");
  leg->AddEntry(hTheory[0],"Powheg+Pythia8","F");
  //leg->AddEntry(hTheory[1],"aMC@NLO+Pythia8","F");
  leg->AddEntry(hTheory[1],"Powheg+Herwigpp","F");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  //TPaveText *pave = new TPaveText(0.2,0.85,0.4,0.87,"NDC");
  TPaveText *pave = new TPaveText(0.2,0.88,0.4,0.90,"NDC");
  //pave->AddText("Hadronic channel");
  pave->AddText("l+jets channel");
  pave->SetBorderSize(0);
  pave->SetFillColor(0);
  pave->SetTextFont(42);
  pave->SetTextSize(0.04);

  TString CAN_TITLE = "CrossSection_"+LEVEL+"_"+VAR;
  
  TCanvas *can = new TCanvas(CAN_TITLE,CAN_TITLE,800,704);
  if (LOGY) {
    gPad->SetLogy();
    //gPad->SetLogx();
  }
  
  can->cd(1)->SetBottomMargin(0.35);
  
  hSignal->GetXaxis()->SetTitleSize(0);
  hSignal->GetXaxis()->SetLabelSize(0);
  
  hSignal->GetYaxis()->SetTitle(YTITLE);
  hSignal->GetYaxis()->SetRangeUser(YMIN,YMAX);

  hSignal->SetMarkerColor(kBlack);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerSize(1.2);
  hSignal->SetLineWidth(2.0);
  hSignal->Draw();
  gSignalWithUnc->Draw("E2 same");
  hTheory[0]->Draw("E2 same");
  hTheory[1]->Draw("E2 same");
  //hTheory[2]->Draw("E2 same");
  hTheoryHIST[0]->Draw("hist same");
  hTheoryHIST[1]->Draw("hist same");
  //hTheoryHIST[2]->Draw("hist same");
  hSignal->Draw("same");
  leg->Draw("same");
  pave->Draw("same");

  gPad->RedrawAxis();

  TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
  pad->SetTopMargin(0.67);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  pad->SetGridy();
  //gPad->SetLogx();

  hRelUncMinus->SetFillColor(kYellow);
  hRelUncPlus->SetFillColor(kYellow);
  hSignalRatio->SetLineWidth(2);
  hSignalRatio->SetMarkerStyle(20);
  hSignalRatio->SetMarkerSize(1.2);

  hTotalRelUncMinus->GetYaxis()->SetRangeUser(RMIN,RMAX);
  hTotalRelUncMinus->GetXaxis()->SetTitle(XTITLE);
  hTotalRelUncMinus->GetYaxis()->SetTitleSize(0.04);
  hTotalRelUncMinus->GetYaxis()->SetLabelSize(0.03);
  hTotalRelUncMinus->GetYaxis()->SetNdivisions(505);
  hTotalRelUncMinus->GetYaxis()->SetTitle("The./data-1");
  hTotalRelUncMinus->Draw("hist");
  hTotalRelUncPlus->Draw("hist same");
  //hRelUncMinus->Draw("hist same");
  //hRelUncPlus->Draw("hist same");
  
  hTheoryRatio[0]->Draw("E2 same");
  hTheoryRatio[1]->Draw("E2 same");
  //hTheoryRatio[2]->Draw("E2 same");
  hTheoryRatioHIST[0]->Draw("hist same");
  hTheoryRatioHIST[1]->Draw("hist same");
  //hTheoryRatioHIST[2]->Draw("hist same");
  hSignalRatio->Draw("same E");
  
  gPad->RedrawAxis();

  CMS_lumi(can,4,0);

  can->Print("PaperPlots/"+TString(can->GetName())+".pdf");

  TCanvas *canNorm = new TCanvas(CAN_TITLE+"_Norm",CAN_TITLE+"_Norm",800,704);
  canNorm->cd(1)->SetBottomMargin(0.35);
  if (LOGY) {
    gPad->SetLogy();
    //gPad->SetLogx();
  }
  
  hSignalNorm->GetXaxis()->SetTitleSize(0);
  hSignalNorm->GetXaxis()->SetLabelSize(0);
  hSignalNorm->GetYaxis()->SetRangeUser(YMIN_NORM,YMAX_NORM);
  hSignalNorm->GetXaxis()->SetTitle(XTITLE);
  hSignalNorm->GetYaxis()->SetTitle(YTITLE_NORM);
  hSignalNorm->SetMarkerColor(kBlack);
  hSignalNorm->SetMarkerStyle(20);
  hSignalNorm->SetMarkerSize(1.2);
  hSignalNorm->SetLineWidth(2.0);
  hSignalNorm->Draw();
  gSignalWithUncNorm->Draw("E2 same");
  hTheoryNorm[0]->Draw("E2 same");
  hTheoryNorm[1]->Draw("E2 same");
  //hTheoryNorm[2]->Draw("E2 same");
  hTheoryNormHIST[0]->Draw("hist same");
  hTheoryNormHIST[1]->Draw("hist same");
  //hTheoryNormHIST[2]->Draw("hist same");
  hSignalNorm->Draw("same");
  leg->Draw("same");
  pave->Draw("same");

  gPad->RedrawAxis();

  TPad *padNorm = new TPad("padNorm","padNorm",0.,0.,1.,1.);
  padNorm->SetTopMargin(0.67);
  padNorm->SetFillColor(0);
  padNorm->SetFillStyle(0);
  padNorm->Draw();
  padNorm->cd(0);
  padNorm->SetGridy();
  //gPad->SetLogx();

  hRelUncMinusNorm->SetFillColor(kYellow);
  hRelUncPlusNorm->SetFillColor(kYellow);
  hSignalRatioNorm->SetLineWidth(2);
  hSignalRatioNorm->SetMarkerColor(kBlack);
  hSignalRatioNorm->SetMarkerStyle(20);
  hSignalRatioNorm->SetMarkerSize(1.2);
  
  hTotalRelUncMinusNorm->GetYaxis()->SetRangeUser(RMIN_NORM,RMAX_NORM);
  hTotalRelUncMinusNorm->GetXaxis()->SetTitle(XTITLE);
  hTotalRelUncMinusNorm->GetYaxis()->SetTitleSize(0.04);
  hTotalRelUncMinusNorm->GetYaxis()->SetLabelSize(0.03);
  hTotalRelUncMinusNorm->GetYaxis()->SetNdivisions(505);
  hTotalRelUncMinusNorm->GetYaxis()->SetTitle("The./data-1");
  hTotalRelUncMinusNorm->Draw("hist");
  hTotalRelUncPlusNorm->Draw("hist same");
  //hRelUncMinusNorm->Draw("hist same");
  //hRelUncPlusNorm->Draw("hist same");
  
  hTheoryRatioNorm[0]->Draw("E2 same");
  hTheoryRatioNorm[1]->Draw("E2 same");
  //hTheoryRatioNorm[2]->Draw("E2 same");
  hTheoryRatioNormHIST[0]->Draw("hist same");
  hTheoryRatioNormHIST[1]->Draw("hist same");
  //hTheoryRatioNormHIST[2]->Draw("hist same");
  hSignalRatioNorm->Draw("same E");
  gPad->RedrawAxis();

  CMS_lumi(canNorm,4,0);

  canNorm->Print("PaperPlots/"+TString(canNorm->GetName())+".pdf");
}
