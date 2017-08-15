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

void getCorr(bool check = false){
  TString channels[2] = {"mu","el"};
  const int nSYS = 21;
  TString sysnames[nSYS] = {"nom","puUp","puDown","JECUp","JECDown","JERUp","JERDown","BTagUp","BTagDown","TopTagUp","TopTagDown","lepUp","lepDown","PDFUp","PDFDown","Q2Up","Q2Down","ISRUp","ISRDown","FSRUp","FSRDown"};
  
  TH1D* h_weight[2][nSYS][2];
  
  for (int ich = 0; ich < 2; ich++){
    for (int isys = 0; isys < nSYS; isys++){
      
      TH2D* h_response;
      TH2D* h_response_noweight;
      TH1F* h_ptGenTop;
      
      if (sysnames[isys].Contains("ISR") || sysnames[isys].Contains("FSR")){
	TFile* f = TFile::Open("histfiles_full2016/hists_PowhegPythia8_"+sysnames[isys]+"_fullTruth_"+channels[ich]+"_"+sysnames[isys]+"_post.root");
	h_response = (TH2D*) f->Get("response2_pt_TH2");
	h_response_noweight = (TH2D*) f->Get("response2_pt_TH2_noweight");
      }
      
      else{ //Stitch final response matrix
	TFile* f_m0to700_p1 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+channels[ich]+"_"+sysnames[isys]+"_post.root");
	TFile* f_m0to700_p2 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+channels[ich]+"_"+sysnames[isys]+"_post.root");
	TFile* f_m700to1000 = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_m700to1000_"+channels[ich]+"_"+sysnames[isys]+"_post.root");
	TFile* f_m1000toInf = TFile::Open("histfiles_full2016/hists_PowhegPythia8_fullTruth_m1000toInf_"+channels[ich]+"_"+sysnames[isys]+"_post.root");
	TH2D* h_response_m0to700_p1 = (TH2D*) f_m0to700_p1->Get("response2_pt_TH2");
	TH2D* h_response_m0to700_p2 = (TH2D*) f_m0to700_p2->Get("response2_pt_TH2");
	TH2D* h_response_m0to700 = (TH2D*) h_response_m0to700_p1->Clone();
	h_response_m0to700->Add(h_response_m0to700_p2);
	h_response_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
	TH2D* h_response_m700to1000 = (TH2D*) f_m700to1000->Get("response2_pt_TH2");
	h_response_m700to1000->Scale(831.76 * 35867.0 * 0.0967 / 38578334.0);
	TH2D* h_response_m1000toInf = (TH2D*) f_m1000toInf->Get("response2_pt_TH2");
	h_response_m1000toInf->Scale(831.76 * 35867.0 * 0.0256 / 24495211.0);
	h_response = (TH2D*) h_response_m0to700->Clone();
	h_response->Add(h_response_m700to1000);
	h_response->Add(h_response_m1000toInf);

	TH2D* h_response_noweight_m0to700_p1 = (TH2D*) f_m0to700_p1->Get("response2_pt_TH2_noweight");
	TH2D* h_response_noweight_m0to700_p2 = (TH2D*) f_m0to700_p2->Get("response2_pt_TH2_noweight");
	TH2D* h_response_noweight_m0to700 = (TH2D*) h_response_noweight_m0to700_p1->Clone();
	h_response_noweight_m0to700->Add(h_response_noweight_m0to700_p2);
	h_response_noweight_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
	TH2D* h_response_noweight_m700to1000 = (TH2D*) f_m700to1000->Get("response2_pt_TH2_noweight");
	h_response_noweight_m700to1000->Scale(831.76 * 35867.0 * 0.0967 / 38578334.0);
	TH2D* h_response_noweight_m1000toInf = (TH2D*) f_m1000toInf->Get("response2_pt_TH2_noweight");
	h_response_noweight_m1000toInf->Scale(831.76 * 35867.0 * 0.0256 / 24495211.0);
	h_response_noweight = (TH2D*) h_response_noweight_m0to700->Clone();
	h_response_noweight->Add(h_response_noweight_m700to1000);
	h_response_noweight->Add(h_response_noweight_m1000toInf);

	TH1F* h_ptGenTop_m0to700_p1 = (TH1F*) f_m0to700_p1->Get("ptGenTop");
	TH1F* h_ptGenTop_m0to700_p2 = (TH1F*) f_m0to700_p2->Get("ptGenTop");
	TH1F* h_ptGenTop_m0to700 = (TH1F*) h_ptGenTop_m0to700_p1->Clone();
	h_ptGenTop_m0to700->Add(h_ptGenTop_m0to700_p2);
	h_ptGenTop_m0to700->Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.));
	TH1F* h_ptGenTop_m700to1000 = (TH1F*) f_m700to1000->Get("ptGenTop");
	h_ptGenTop_m700to1000->Scale(831.76 * 35867.0 * 0.0967 / 38578334.0);
	TH1F* h_ptGenTop_m1000toInf = (TH1F*) f_m1000toInf->Get("ptGenTop");
	h_ptGenTop_m1000toInf->Scale(831.76 * 35867.0 * 0.0256 / 24495211.0);
	h_ptGenTop = (TH1F*) h_ptGenTop_m0to700->Clone();
	h_ptGenTop->Add(h_ptGenTop_m700to1000);
	h_ptGenTop->Add(h_ptGenTop_m1000toInf);
      }
      
      //Projection of h_response onto X is ptRecoTop
      //Projection of h_response onto Y WOULD BE ptGenTop if misses were weighted correctly
      for (int ih = 0; ih < 2; ih++){
	h_weight[ich][isys][ih] = (TH1D*) h_response->ProjectionY();
	h_weight[ich][isys][ih]->Reset();
	h_weight[ich][isys][ih]->SetName(Form("corrFac%i_miss_",ih+1)+channels[ich]+"_"+sysnames[isys]);
	for (int iy = 2-ih; iy < h_response->GetYaxis()->GetNbins()+ih; iy++){
	  float sum_weight = h_response->Integral(2-ih,h_response->GetXaxis()->GetNbins()-1+ih,iy,iy);
	  float sum_noweight = h_response_noweight->Integral(2-ih,h_response_noweight->GetXaxis()->GetNbins()-1+ih,iy,iy);
	  float underflow = h_response_noweight->Integral(0,1-ih,iy,iy) + h_response_noweight->Integral(h_response_noweight->GetXaxis()->GetNbins()+ih,h_response_noweight->GetXaxis()->GetNbins()+1,iy,iy);
	  float weight = underflow > 0.0 ? 1.0 + (sum_noweight - sum_weight) / underflow : 1.0;
	  h_weight[ich][isys][ih]->SetBinContent(iy,weight);
	}
      }

      if (check && isys == 0){
	cout << "Bin h_response ptGenTop" << endl;
	for (int ibiny = 1; ibiny < h_response->GetNbinsY()-1; ibiny++){
	  cout << "[" << h_response->GetYaxis()->GetBinLowEdge(ibiny+1) << "," << h_response->GetYaxis()->GetBinUpEdge(ibiny+1) << "] " << h_response->Integral(0,h_response->GetNbinsX()+1,ibiny+1,ibiny+1) << " " << h_ptGenTop->GetBinContent(ibiny) << endl;
	}
      }
    }
  }
  TFile* fout = new TFile("corrFac_miss.root","recreate");
  for (int ich = 0; ich < 2; ich++){
    for (int isys = 0; isys < nSYS; isys++){
      for (int ih = 0; ih < 2; ih++){
	h_weight[ich][isys][ih]->Write();
      }
    }
  }
  return;
}
