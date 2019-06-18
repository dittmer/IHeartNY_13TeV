// Script to run makeHists.C
// Syntax is makeHists(TString INDIR, TString OUTDIR, TString sample, TString channel, bool isData = false, bool isSignal = false, TString lepID = "Medium", TString iso = "None", bool doHemiCuts = false, float metCut = 0.0, bool doTriangular = false, bool isQCD = false, TString systematic = "nom", int oddOrEven = 0, bool usePost = false)

#include "TROOT.h"
#include "TSystem.h"

using namespace std;

#include "makeHists.C"
#include "BTagCalibrationStandalone.cpp"

//R__LOAD_LIBRARY(RooUnfold/libRooUnfold.so)

void runMakeHists(TString toMake = "prefit"){


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


  //gROOT->ProcessLine(".include RooUnfold/src");

  const int nBKG = 22;
  TString bkgMCnames[nBKG] = {
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

  const int nSYS = 13;
  TString sysnames[nSYS] = {"nom","puUp","puDown","JECUp","JECDown","JERUp","JERDown","lepUp","lepDown","BTagUp","BTagDown","TopTagUp","TopTagDown"};
  const int nTHSYS = 4;
  TString thsysnames[nTHSYS] = {"PDFUp","PDFDown","Q2Up","Q2Down"};
  const int nSAMPLES = 10;
  TString thsamples[nSAMPLES] = {"ISRUp","ISRDown","FSRUp","FSRDown","TuneUp","TuneDown","HdampUp","HdampDown","ErdOn","Herwig"};
  const int nISO = 7;
  TString isoWPs[nISO] = {"MiniIso10","MiniIso20","2DisoPt25","2DisoPt45","2DisoB2G","2DisoIHNY","Loose"};

  // ----------------------------------------------
  // Make histfiles for lepton optimization studies

  if (toMake == "all" || toMake == "lepOpt"){

    for (int ii = 0; ii < nISO; ii++){
      // Run signal
      makeHists("skimTrees_full2016/","histfiles_full2016_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
      makeHists("skimTrees_full2016/","histfiles_full2016_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_"+name_TTbarNom,"el",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
      
      //Run QCD
      for (int jj = nBKG-6; jj < nBKG; jj ++){
	makeHists("skimTrees_full2016/","histfiles_full2016_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"mu",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
	makeHists("skimTrees_full2016/","histfiles_full2016_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"el",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);    
      }
    }

    // Run signal
    makeHists("skimTrees_full2016/","histfiles_full2016_tMu_tEl","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
    makeHists("skimTrees_full2016/","histfiles_full2016_tMu_tEl","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
    
    //Run QCD
    for (int jj = nBKG-6; jj < nBKG; jj ++){
      makeHists("skimTrees_full2016/","histfiles_full2016_tMu_tEl",bkgMCnames[jj],"mu",false,false,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
      makeHists("skimTrees_full2016/","histfiles_full2016_tMu_tEl",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);    
    }
  }

  // ----------------------------------------------------------------------
  // Make histfiles with final lepton selection, pre-selection-optimization

  if (toMake == "all" || toMake == "selOpt"){
    //run data
    makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10","Data_mu","mu",true,false,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
    makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10","Data_el","el",true,false,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    
    // Run signal
    makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
    makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    
    //run other MCs
    for (int jj = 0; jj < nBKG; jj ++){
      makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
      makeHists("skimTrees_full2016/","histfiles_full2016_mMu_tEl_MiniIso10",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    }
  }
    
  // ----------------------------------------------
  // Make histfiles for plotting / fit

  if (toMake == "all" || toMake == "prefit" || toMake == "postfit"){

    bool postTopTagSF = false;
    int nSysToRun = nSYS;
    if (toMake == "postfit") {
      postTopTagSF = true;
      nSysToRun = 1;
    }

    //run data
    makeHists("skimTrees_full2016/","histfiles_full2016_latest","Data_mu","mu",true,false,"Medium","MiniIso10",true,35.0,false,false,"nom",0,false); //Data
    makeHists("skimTrees_full2016/","histfiles_full2016_latest","Data_el","el",true,false,"Tight","MiniIso10",true,50.0,true,false,"nom",0,false);
    
    makeHists("skimTrees_full2016/","histfiles_full2016_latest","Data_mu","mu",true,false,"Medium","MiniIso10",true,35.0,false,true,"nom",0,false); //QCD
    makeHists("skimTrees_full2016/","histfiles_full2016_latest","Data_el","el",true,false,"Medium","MiniIso10",true,50.0,true,true,"nom",0,false);

    for (int ii = 0; ii < nSysToRun; ii++){
 
      // Run signal
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); //Signal
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,true,sysnames[ii],0,postTopTagSF); //QCD
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Medium","MiniIso10",true,50.0,true,true,sysnames[ii],0,postTopTagSF);

      //run other MCs
      for (int jj = 0; jj < nBKG; jj++){
    	makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); //Signal
    	makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);    
    	makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,true,sysnames[ii],0,postTopTagSF); //QCD
    	makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"el",false,false,"Medium","MiniIso10",true,50.0,true,true,sysnames[ii],0,postTopTagSF);    
      }
    }
  }
  
  if (toMake == "unfold") {

    cout << "unfold" << endl;

    bool postTopTagSF = true;
    /*
    // non-stiched sample for plots    
    makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,false);
    makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF,false);

    for (int ii = 0; ii < nSYS; ii++){

      cout << "systematic " << sysnames[ii] << endl;
      
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF,true); //Signal
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF,true); //Odd
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF,true); //Even
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF,true);
      
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF);

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF);
    }

    for (int ii = 0; ii < nTHSYS; ii++){

      cout << "theory systematic " << thsysnames[ii] << endl;      

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_p2","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m700to1000","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF);

      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth_m1000toInf","el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF);
    }

    // Theory variant samples
    // Note these are the original samples! Not from newJER area
    for (int ii = 0; ii < nSAMPLES; ii++){
      makeHists("skimTrees_full2016/","histfiles_full2016_latest","PowhegPythia8_"+thsamples[ii]+"_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsamples[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016_latest","PowhegPythia8_"+thsamples[ii]+"_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsamples[ii],0,postTopTagSF);
    }

    // TTbar for QCD
    makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,false,true,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest","PowhegPythia8_fullTruth","el",false,true,"Medium","MiniIso10",true,50.0,true,true,"nom",0,postTopTagSF);
*/
    // Backgrounds
    for (int jj = 0; jj < nBKG; jj++){
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); //Signal
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);    
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,true,"nom",0,postTopTagSF); //QCD
      makeHists("skimTrees_full2016_newJER/","histfiles_full2016_latest",bkgMCnames[jj],"el",false,false,"Medium","MiniIso10",true,50.0,true,true,"nom",0,postTopTagSF);    
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------------------------

  if (toMake == "unfoldParticle") {

    cout << "unfoldParticle" << endl;

    bool postTopTagSF = true;

    // non-stiched sample for plots    
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,false);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF,false);

    for (int ii = 0; ii < nSYS; ii++){

      cout << "systematic " << sysnames[ii] << endl;
      
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF,true); //Signal
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF,true); //Odd
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF,true); //Even
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF,true);
      
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF);

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,sysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,sysnames[ii],2,postTopTagSF);
    }

    for (int ii = 0; ii < nTHSYS; ii++){

      cout << "theory systematic " << thsysnames[ii] << endl;

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF,true); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF,true);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF,true);

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF);

      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],0,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],1,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsysnames[ii],2,postTopTagSF); 
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],1,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsysnames[ii],2,postTopTagSF);
    }

    // Theory variant samples
    for (int ii = 0; ii < nSAMPLES; ii++){
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+thsamples[ii]+"_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,thsamples[ii],0,postTopTagSF);
      makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+thsamples[ii]+"_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,thsamples[ii],0,postTopTagSF);
    }

    makeHists("skimTrees_full2016/","histfiles_full2016","TTJets_amcatnloFXFX_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","TTJets_amcatnloFXFX_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);


    //Data, nominal
    makeHists("skimTrees_full2016/","histfiles_full2016","Data_mu","mu",true,false,"Medium","MiniIso10",true,35.0,false,false,"nom",0,false);
    makeHists("skimTrees_full2016/","histfiles_full2016","Data_el","el",true,false,"Tight","MiniIso10",true,50.0,true,false,"nom",0,false);

    //Data for QCD
    makeHists("skimTrees_full2016/","histfiles_full2016","Data_mu","mu",true,false,"Medium","MiniIso10",true,35.0,false,true,"nom",0,false);
    makeHists("skimTrees_full2016/","histfiles_full2016","Data_el","el",true,false,"Medium","MiniIso10",true,50.0,true,true,"nom",0,false);

    // TTbar for QCD
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,false,"Medium","MiniIso10",true,35.0,false,true,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,false,"Medium","MiniIso10",true,50.0,true,true,"nom",0,postTopTagSF);

    // Backgrounds
    for (int jj = 0; jj < nBKG; jj++){
      makeHists("skimTrees_full2016/","histfiles_full2016",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); //Signal
      makeHists("skimTrees_full2016/","histfiles_full2016",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);    
      makeHists("skimTrees_full2016/","histfiles_full2016",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,false,true,"nom",0,postTopTagSF); //QCD
      makeHists("skimTrees_full2016/","histfiles_full2016",bkgMCnames[jj],"el",false,false,"Medium","MiniIso10",true,50.0,true,true,"nom",0,postTopTagSF);    
    }
  }

  // ------------------------------------------------------------------------------------------------------------------------------------------

  if (toMake == "unfoldParticleNom") {
    
    bool postTopTagSF = true;

    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,true); //Signal
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",1,postTopTagSF,true); //Odd
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",2,postTopTagSF,true); //Even
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,true); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",1,postTopTagSF,true); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",2,postTopTagSF,true); 

    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",1,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m700to1000,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",2,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",1,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbar_m1000toInf,"mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",2,postTopTagSF); 

    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",1,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",2,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",1,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_"+name_TTbarNom_p2,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",2,postTopTagSF,true);

    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",1,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m700to1000,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",2,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",1,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_fullTruth_"+name_TTbar_m1000toInf,"el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",2,postTopTagSF);
  }

  // ------------------------------------------------------------------------------------------------------------------------------------------

  if (toMake == "test") {
    bool postTopTagSF = true;
    //bool postTopTagSF = false;

    /*
    makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_fullTruth_PL","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_fullTruth_PL_p2","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,true);
    makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m700to1000_PL","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 
    makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PL","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 
    */

    //makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PL","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",1,postTopTagSF); 
    //makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PL","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",2,postTopTagSF); 
    //makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 

    //makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF); 
    //makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_m1000toInf_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF); 
    // makeHists("skimTrees_full2016/","histfiles_full2016_debug","PowhegPythia8_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF,true); //Signal

    makeHists("skimTrees_full2016/","histfiles_full2016","TTJets_amcatnloFXFX_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","TTJets_amcatnloFXFX_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);

    /*
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_mtop1715_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_mtop1715_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);

    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_mtop1735_PLnew","mu",false,true,"Medium","MiniIso10",true,35.0,false,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_full2016/","histfiles_full2016","PowhegPythia8_mtop1735_PLnew","el",false,true,"Tight","MiniIso10",true,50.0,true,false,"nom",0,postTopTagSF);
    */
  }

}
