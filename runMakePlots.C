// Script to run makePlots.cc
// Syntax is makePlots(TString channel, TString var, TString region)

#include "TSystem.h"

#include "makePlots.cc"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void runMakePlots(TString toPlot = "final"){

  //gSystem->CompileMacro("makePlots.cc");
  //cout << "Compilation successful!" << endl;

  gROOT->ProcessLine("gErrorIgnoreLevel = 1;");

  TString channels[2] = {"mu","el"};
  const int nregion = 4;
  TString regions[nregion] = {"Pre","0t","1t0b","1t1b"};
  const int nhist = 24;
  TString hists[nhist] = {"metPt",
			  "ht",
			  "htLep",
			  "ak4jetPt",
			  "ak4jetEta",
			  "ak4jetPhi",
			  "ak4jetMass",
			  "ak4jetCSV",
			  "ak8jetPt",
			  "ak8jetEta",
			  "ak8jetPhi",
			  "ak8jetY",
			  "ak8jetMass",
			  "ak8jetTau32",
			  "ak8jetTau21",
			  "ak8jetSDmass",
			  "ak8jetSubjetMaxCSV",
			  "lepPt",
			  "lepEta",
			  "lepAbsEta",
			  "lepSignEta",
			  "lepPhi",
			  "lepBJetdR",
			  //"lepTJetdR",
			  //"lepBJetPtRel",
			  "lepMiniIso"
  };

  if (toPlot == "all" || toPlot == "lepSelOpt"){
    makeEffPlots("mu","ak8jetPt");
    makeEffPlots("el","ak8jetPt");
    makeEffPlots("mu","leadLepPt");
    makeEffPlots("el","leadLepPt");
    drawROCCurve("mu");
    drawROCCurve("el");
    //makeEffPlotsFinal();
    
    cout << endl << "Finished with lepton selection optimization plots!" << endl << endl;
  }

  if (toPlot == "all" || toPlot == "selOpt"){

    make1DScans("mu");
    make1DScans("el");

    /*
    TString opthists[5] = {"metPt","ht","htLep","lepBJetdR","lepTJetdR"};
    
    for (int ii = 0; ii < 5; ii++){
      makePlots("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","mu",opthists[ii],"Pre",true,false,false);
      makePlots("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","el",opthists[ii],"Pre",true,false,false);
    }

    makeTable("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","el",false,true,false);
    
    cout << endl << "Finished with selection optimization plots!" << endl << endl;
    */
  }

  if (toPlot=="louise") {
    /*
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetPt","0t",false,true,true);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetPt","1t0b",false,true,true);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetPt","1t1b",false,true,true);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetY","0t",false,true,true);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetY","1t0b",false,true,true);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","comb","ak8jetY","1t1b",false,true,true);
    */
    combineResults("mu","comb1");
    combineResults("el","comb1");
  }

  if (toPlot == "all" || toPlot == "final"){

    bool unBlind = true;

    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/",channels[ii],hists[kk],regions[jj],false,unBlind,false);
	}
      }
    }
    
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","nAK4jet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","nAK4jet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","nBjet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","nBjet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","nAK8jet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","nAK8jet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","nTjet","Pre",false,unBlind,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","nTjet","Pre",false,unBlind,false);

    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","elPtRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","elEtaRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","elMiniIsoRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","elMiniIsoRaw_b","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","elMiniIsoRaw_e","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","muPtRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","muEtaRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","muMiniIsoRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","muMiniIsoRaw_b","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","mu","muMiniIsoRaw_e","",true,true,false);

    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","lepMETdPhiRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","ak4METdPhiRaw","",true,true,false);
    makePlots("histfiles_full2016_latest/","histfiles_full2016_latest/","el","ak8METdPhiRaw","",true,true,false);

    //plot2D("ttbar","bTagSFvsPt");
    //plot2D("qcd","bTagSFvsPt");
    //plot2D("ttbar","bTagSFvsCSV");
    //plot2D("qcd","bTagSFvsCSV");
    
    cout << endl << "Finished with regular plots!" << endl << endl;

    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nhist; jj++){
	makeQCDComp("histfiles_full2016/","histfiles_full2016/",channels[ii],hists[jj]);
      }
      makeQCDComp("histfiles_full2016/","histfiles_full2016/",channels[ii],"nAK4jet");
      makeQCDComp("histfiles_full2016/","histfiles_full2016/",channels[ii],"nBjet");
      makeQCDComp("histfiles_full2016/","histfiles_full2016/",channels[ii],"nAK8jet");
      makeQCDComp("histfiles_full2016/","histfiles_full2016/",channels[ii],"nTjet");
    }

    //for (int ii = 0; ii < nhist; ii++){
    //  compLepQCD("histfiles_full2016/",hists[ii],true);
    //  compLepQCD("histfiles_full2016/",hists[ii],false);
    //}
    
    //cout << endl << "Finished with QCD comparison plots!" << endl << endl;
    /*
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  compareShapes("histfiles_full2016/","histfiles_full2016/",channels[ii],hists[kk],regions[jj]);
	}
      }
    }
    */
    //cout << endl << "Finished making shape comparisons!" << endl << endl;

    makeCombineInputs("histfiles_full2016_latest/","histfiles_full2016_latest/","MC");
    makeCombineInputs("histfiles_full2016_latest/","histfiles_full2016_latest/","data");
    
    cout << endl << "Finished making combine inputs!" << endl << endl;

    makeTable("histfiles_full2016_latest/","histfiles_full2016_latest/","mu",true,true,unBlind);
    cout << endl;
    makeTable("histfiles_full2016_latest/","histfiles_full2016_latest/","el",true,true,unBlind);
    cout << endl;
    makeTable("histfiles_full2016_latest/","histfiles_full2016_latest/","mu",false,false,unBlind);
    cout << endl;
    makeTable("histfiles_full2016_latest/","histfiles_full2016_latest/","el",false,false,unBlind);
  }

  if (toPlot == "post"){
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  makePlots("histfiles_full2016/","histfiles_full2016/",channels[ii],hists[kk],regions[jj],false,true,true);
	}
      }
    }
    
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","nAK4jet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","nAK4jet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","nBjet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","nBjet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","nAK8jet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","nAK8jet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","nTjet","Pre",false,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","nTjet","Pre",false,true,true);

    makePlots("histfiles_full2016/","histfiles_full2016/","el","elPtRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","elEtaRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","elMiniIsoRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","elMiniIsoRaw_b","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","el","elMiniIsoRaw_e","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","muPtRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","muEtaRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","muMiniIsoRaw","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","muMiniIsoRaw_b","",true,true,true);
    makePlots("histfiles_full2016/","histfiles_full2016/","mu","muMiniIsoRaw_e","",true,true,true);
    
    cout << endl << "Finished with postfit plots!" << endl << endl;

    makeTable("histfiles_full2016/","histfiles_full2016/","mu",false,false,true,true);
    cout << endl;
    makeTable("histfiles_full2016/","histfiles_full2016/","el",false,false,true,true);

  }

  if (toPlot == "tablepost"){
    makeTable("histfiles_full2016/","histfiles_full2016/","mu",false,false,true,true);
    cout << endl;
    makeTable("histfiles_full2016/","histfiles_full2016/","el",false,false,true,true);
  }

  if (toPlot == "combine"){
    combineResults("mu","mu1");
    combineResults("el","el1");
    combineResults("mu","comb1");
    combineResults("el","comb1");
  }

  if (toPlot == "troubleshoot"){
    for (int ii = 1; ii < 2; ii++){
      for (int jj = 0; jj < nhist; jj++){
      	troubleshootQCD("histfiles_full2016/",hists[jj],channels[ii]);
      }
    }
  }

  if (toPlot == "BtagSF"){
    for (int ii = 0; ii < 2; ii++){
      calcBtagSF(channels[ii],"nom");
      calcBtagSF(channels[ii],"BTagUp");
      calcBtagSF(channels[ii],"BTagDown");
    }
  }

  if (toPlot == "elID"){
    TString vars[9] = {"dEtaIn","dPhiIn","Full5x5siee","HoE","ooEmooP","MissHits","ConVeto","Dxy","Dz"};
    for (int ii = 0; ii < 9; ii++){
      plot2D("qcd","elPtVs"+vars[ii]+"Raw");
      plot2D("qcd","elEtaVs"+vars[ii]+"Raw");
      plot2D("qcd","elPtVs"+vars[ii]+"Nm1");
      plot2D("qcd","elEtaVs"+vars[ii]+"Nm1");
    }
    checkElID();
  }
}
