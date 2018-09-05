#!/usr/bin/env python

import math

import sys

from ROOT import gRandom, TH1, TH1D, TH1F, TH2F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD, TGraph, TUnfold, Double, TSpline, TSpline3, TUnfoldDensity, TUnfoldSys, TUnfold, TAttLine, TStyle, THStack, TMatrixT, TMatrixD

from array import array
import string

gROOT.Macro("rootlogon.C")
gROOT.SetBatch(True)

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0)

gStyle.SetTitleFont(43)
gStyle.SetTitleFont(43, "XYZ")
gStyle.SetTitleSize(30, "XYZ")
gStyle.SetLabelFont(43, "XYZ")
gStyle.SetLabelSize(24, "XYZ")

gStyle.SetPadTopMargin(0.07)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.18)

gROOT.ProcessLine('gErrorIgnoreLevel = kError;')

# ------------------
# Helper function
# ------------------

def myText(x, y, color, text) :
  l = TLatex()
  l.SetTextSize(0.043) 
  l.SetTextFont(42)
  l.SetNDC()
  l.SetTextColor(color)
  l.DrawLatex(x,y,text)

# -------------------------------------------------------------------------------------
# Script for doing RooUnfold on the ttbar differential cross secion
# -------------------------------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()

# -------------------------------------------------------------------------------------
# input options
# -------------------------------------------------------------------------------------

parser.add_option('--toUnfold', metavar='F', type='string', action='store',
                  default='pt',
                  dest='toUnfold',
                  help='Distribution to unfold (pt or y)')

parser.add_option('--level', metavar='F', type='string', action='store',
                  default='gen',
                  dest='level',
                  help='Level to unfold (gen or part)')

# -------------------------------------------------------------------------------------
# load options 
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []
  
print "TUnfold version is " + str(TUnfold.GetTUnfoldVersion())
TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# Define helper functions
# -------------------------------------------------------------------------------------
def noNegBins( hist ):
    for i in xrange(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i) < 0.0 :
            hist.SetBinContent(i,0.0)

# Reweight response underflow (effective anti-tag SF)
def antiTagWeight(h_true, h_response):
    for ii in xrange(1,h_response.GetNbinsY()+1):
        rowsum = h_response.Integral(1,h_response.GetNbinsX()+1,ii,ii)
        binweight = (h_true.GetBinContent(ii) - rowsum)/h_response.GetBinContent(0,ii)
        h_response.SetBinContent(0,ii,h_response.GetBinContent(0,ii)*binweight)
        
    
# Combine underflow / overflow
def convertForTUnfold(h_response):
    h_response.SetBinContent(0,0,0)
    h_response.SetBinContent(h_response.GetNbinsX()+1,h_response.GetNbinsY()+1,0)
    for ii in xrange(1,h_response.GetNbinsY()+1):
        underflow = h_response.GetBinContent(0,                       ii)
        overflow  = h_response.GetBinContent(h_response.GetNbinsX()+1,ii)
        h_response.SetBinContent(0,ii,underflow+overflow)
    for jj in xrange(1,h_response.GetNbinsX()+1):
        underflow = h_response.GetBinContent(jj,0)
        overflow  = h_response.GetBinContent(jj,h_response.GetNbinsY()+1)
        h_response.SetBinContent(jj,0,underflow+overflow)

# Subtract fakes from input
def removeFakes(h_input, h_response):
    for ii in xrange(1,h_response.GetNbinsX()+1):
        if h_response.Integral(ii,ii,0,h_response.GetNbinsY()+1) > 0.0 :
            fakefraction = h_response.GetBinContent(ii,0) / h_response.Integral(ii,ii,0,h_response.GetNbinsY()+1)
            h_input.SetBinContent(ii,h_input.GetBinContent(ii)*(1.0-fakefraction))
            h_input.SetBinError(ii,h_input.GetBinError(ii)*(1.0-fakefraction))
            h_response.SetBinContent(ii,0,0)
            h_response.SetBinError(ii,0,0)


def drawCMS( x1, y1, size = 0.056, status = "Work In Progress") :
    cmsTextSize = size
    extraTextSize = 0.76 * size
    
    l = TLatex()
    l.SetTextSize(cmsTextSize)
    l.SetTextFont(61)
    l.SetTextAngle(0)
    l.SetNDC()
    l.SetTextColor(1)

    lp = TLatex()
    lp.SetTextSize(extraTextSize) 
    lp.SetTextFont(52)
    lp.SetNDC()
    lp.SetTextColor(1)

    l.DrawLatex(x1,y1,"CMS")
    if status is not None:
        lp.DrawLatex(x1+0.1,y1,status)

class Background :
    def __init__(self, name, norm, err, color=1):
        self.name = name
        self.hist = None
        self.norm = norm
        self.err  = err
        self.color = color
    def addHist(self, hist):
        if self.hist is None:
            self.hist = hist.Clone()
        else:
            self.hist.Add(hist)
    
    
# -------------------------------------------------------------------------------------
# Define normalization constants
# -------------------------------------------------------------------------------------

lum = 35867.0

PowhegPythia8_norm         = 831.76         * lum / 77229341.
SingleTop_t_t_norm         = 136.02         * lum / 67240808.
SingleTop_tbar_t_norm      = 80.95          * lum / 38811017.
SingleTop_t_tW_norm        = 19.3           * lum / 8629641.
SingleTop_tbar_tW_norm     = 19.3           * lum / 8681541.
SingleTop_s_norm           = 10.32 * 0.322  * lum / 9651642.
WJets_HT100to200_norm      = 1345.0 * 1.21  * lum / 39617787.
WJets_HT200to400_norm      = 359.7 * 1.21   * lum / 19914590.
WJets_HT400to600_norm      = 48.91 * 1.21   * lum / 5796237.
WJets_HT600to800_norm      = 12.05 * 1.21   * lum / 14822888.
WJets_HT800to1200_norm     = 5.501 * 1.21   * lum / 6200954.
WJets_HT1200to2500_norm    = 1.329 * 1.21   * lum / 6324934.
WJets_HT2500toInf_norm     = 0.03216 * 1.21 * lum / 2384260.
ZJets_norm                 = 5765           * lum / 42923575.
WW_norm                    = 118.7          * lum / 6987124.
WZ_norm                    = 44.9           * lum / 2995828.
ZZ_norm                    = 15.4           * lum / 990064.
QCD_HT500to700_norm        = 32100.         * lum / 18929951.
QCD_HT700to1000_norm       = 6831.          * lum / 15629253.  
QCD_HT1000to1500_norm      = 1207.          * lum / 4767100. 
QCD_HT1500to2000_norm      = 119.9          * lum / 3970819.
QCD_HT2000toInf_norm       = 25.24          * lum / 1991645.

response_name = "response_"+options.toUnfold+"_split_TH2"
hMeas_name    = options.toUnfold+"RecoTop_split"
hTrue_name    = options.toUnfold+"GenTop"

if options.level == "part":
    response_name += "_PL"
    hTrue_name = hTrue_name.replace("Gen","Part")

labelstring1 = "quark" if (options.level == "gen") else "jet"
labelstring2 = "p_{T} [GeV]" if (options.toUnfold == "pt") else "rapidity"

# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

response = {}
thisMeas = {}
thisTrue = {}
thisExpect = {}
Hres_sys = {}
backgrounds = {}

channels = ["mu","el"]
sysnames = ["JEC","JER","BTag","TopTag","lep","pu","PDF","Q2"]
thsysnames = ["ISR","FSR","Tune","Hdamp","ErdOn","Herwig"]
allsysnames = sysnames+thsysnames

longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","Pileup","PDF Uncertainty","#mu_{R}, #mu_{F} scales","ISR","FSR","Tune","ME-PS matching","Color reconnection","Parton shower"]
variants = ["Up","Down"]

for channel in channels:
    f_data = TFile("histfiles_full2016/hists_Data_"+channel+".root")
    f_QCD  = TFile("histfiles_full2016/hists_Data_"+channel+"_qcd.root")

    f_ttbar_m0to700_p1 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_PL_"+channel+"_nom_post.root")
    f_ttbar_m0to700_p2 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_PL_p2_"+channel+"_nom_post.root")
    f_ttbar_m700to1000 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_m700to1000_PL_"+channel+"_nom_post.root")
    f_ttbar_m1000toInf = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_m1000toInf_PL_"+channel+"_nom_post.root")
    
    f_ttbar_nonsemilep   = TFile("histfiles_full2016/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_post.root")
    f_T_t                = TFile("histfiles_full2016/hists_SingleTop_t_t_"+channel+"_nom_post.root")
    f_Tbar_t             = TFile("histfiles_full2016/hists_SingleTop_tbar_t_"+channel+"_nom_post.root")
    f_T_tW               = TFile("histfiles_full2016/hists_SingleTop_t_tW_"+channel+"_nom_post.root")
    f_Tbar_tW            = TFile("histfiles_full2016/hists_SingleTop_tbar_tW_"+channel+"_nom_post.root")
    f_T_s                = TFile("histfiles_full2016/hists_SingleTop_s_"+channel+"_nom_post.root")    
    f_WJets_HT100to200   = TFile("histfiles_full2016/hists_WJets_HT100to200_"+channel+"_nom_post.root")
    f_WJets_HT200to400   = TFile("histfiles_full2016/hists_WJets_HT200to400_"+channel+"_nom_post.root")
    f_WJets_HT400to600   = TFile("histfiles_full2016/hists_WJets_HT400to600_"+channel+"_nom_post.root")
    f_WJets_HT600to800   = TFile("histfiles_full2016/hists_WJets_HT600to800_"+channel+"_nom_post.root")
    f_WJets_HT800to1200  = TFile("histfiles_full2016/hists_WJets_HT800to1200_"+channel+"_nom_post.root")
    f_WJets_HT1200to2500 = TFile("histfiles_full2016/hists_WJets_HT1200to2500_"+channel+"_nom_post.root")
    f_WJets_HT2500toInf  = TFile("histfiles_full2016/hists_WJets_HT2500toInf_"+channel+"_nom_post.root")
    f_ZJets              = TFile("histfiles_full2016/hists_ZJets_"+channel+"_nom_post.root")
    f_WW                 = TFile("histfiles_full2016/hists_WW_"+channel+"_nom_post.root")
    f_WZ                 = TFile("histfiles_full2016/hists_WZ_"+channel+"_nom_post.root")
    f_ZZ                 = TFile("histfiles_full2016/hists_ZZ_"+channel+"_nom_post.root")
    
    f_qcd_ttbar              = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+channel+"_nom_qcd_post.root")
    f_qcd_ttbar_nonsemilep   = TFile("histfiles_full2016/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_qcd_post.root")
    f_qcd_T_t                = TFile("histfiles_full2016/hists_SingleTop_t_t_"+channel+"_nom_qcd_post.root")
    f_qcd_Tbar_t             = TFile("histfiles_full2016/hists_SingleTop_tbar_t_"+channel+"_nom_qcd_post.root")
    f_qcd_T_tW               = TFile("histfiles_full2016/hists_SingleTop_t_tW_"+channel+"_nom_qcd_post.root")
    f_qcd_Tbar_tW            = TFile("histfiles_full2016/hists_SingleTop_tbar_tW_"+channel+"_nom_qcd_post.root")
    f_qcd_T_s                = TFile("histfiles_full2016/hists_SingleTop_s_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT100to200   = TFile("histfiles_full2016/hists_WJets_HT100to200_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT200to400   = TFile("histfiles_full2016/hists_WJets_HT200to400_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT400to600   = TFile("histfiles_full2016/hists_WJets_HT400to600_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT600to800   = TFile("histfiles_full2016/hists_WJets_HT600to800_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT800to1200  = TFile("histfiles_full2016/hists_WJets_HT800to1200_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT1200to2500 = TFile("histfiles_full2016/hists_WJets_HT1200to2500_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT2500toInf  = TFile("histfiles_full2016/hists_WJets_HT2500toInf_"+channel+"_nom_qcd_post.root")
    f_qcd_ZJets              = TFile("histfiles_full2016/hists_ZJets_"+channel+"_nom_qcd_post.root")
    f_qcd_WW                 = TFile("histfiles_full2016/hists_WW_"+channel+"_nom_qcd_post.root")
    f_qcd_WZ                 = TFile("histfiles_full2016/hists_WZ_"+channel+"_nom_qcd_post.root")
    f_qcd_ZZ                 = TFile("histfiles_full2016/hists_ZZ_"+channel+"_nom_qcd_post.root")
    
    # -------------------------------------------------------------------------------------
    # Get response matrix
    # -------------------------------------------------------------------------------------
    response_m0to700_p1 = f_ttbar_m0to700_p1.Get(response_name)
    response_m0to700_p2 = f_ttbar_m0to700_p2.Get(response_name)
    response_m700to1000 = f_ttbar_m700to1000.Get(response_name)
    response_m1000toInf = f_ttbar_m1000toInf.Get(response_name)
    response_m0to700_p1.Sumw2()
    response_m0to700_p2.Sumw2()
    response_m700to1000.Sumw2()
    response_m1000toInf.Sumw2()
    response_m0to700 = response_m0to700_p1.Clone()
    response_m0to700.Add(response_m0to700_p2)
    response_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
    response_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
    response_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
    response[channel] = response_m0to700.Clone()
    response[channel].Add(response_m700to1000)
    response[channel].Add(response_m1000toInf)

    # -------------------------------------------------------------------------------------
    # Get systematic variations
    # -------------------------------------------------------------------------------------
    for sysname in sysnames:
        for var in variants:
            f_ttbar_sys_m0to700_p1 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_PL_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m0to700_p2 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_PL_p2_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m700to1000 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_m700to1000_PL_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m1000toInf = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_m1000toInf_PL_"+channel+"_"+sysname+var+"_post.root")    

            response_sys_m0to700_p1 = f_ttbar_sys_m0to700_p1.Get(response_name)
            response_sys_m0to700_p2 = f_ttbar_sys_m0to700_p2.Get(response_name)
            response_sys_m700to1000 = f_ttbar_sys_m700to1000.Get(response_name)
            response_sys_m1000toInf = f_ttbar_sys_m1000toInf.Get(response_name)
            response_sys_m0to700_p1.Sumw2()
            response_sys_m0to700_p2.Sumw2()
            response_sys_m700to1000.Sumw2()
            response_sys_m1000toInf.Sumw2()
            response_sys_m0to700 = response_sys_m0to700_p1.Clone()
            response_sys_m0to700.Add(response_sys_m0to700_p2)
            response_sys_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
            response_sys_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
            response_sys_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
            response_sys = response_sys_m0to700.Clone()
            response_sys.Add(response_sys_m700to1000)
            response_sys.Add(response_sys_m1000toInf)

            true_sys_m0to700_p1 = f_ttbar_sys_m0to700_p1.Get(hTrue_name)
            true_sys_m0to700_p2 = f_ttbar_sys_m0to700_p2.Get(hTrue_name)
            true_sys_m700to1000 = f_ttbar_sys_m700to1000.Get(hTrue_name)
            true_sys_m1000toInf = f_ttbar_sys_m1000toInf.Get(hTrue_name)
            true_sys_m0to700_p1.Sumw2()
            true_sys_m0to700_p2.Sumw2()
            true_sys_m700to1000.Sumw2()
            true_sys_m1000toInf.Sumw2()
            true_sys_m0to700 = true_sys_m0to700_p1.Clone()
            true_sys_m0to700.Add(true_sys_m0to700_p2)
            true_sys_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
            true_sys_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
            true_sys_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
            true_sys = true_sys_m0to700.Clone()
            true_sys.Add(true_sys_m700to1000)
            true_sys.Add(true_sys_m1000toInf)
            
            antiTagWeight(true_sys,response_sys)
            convertForTUnfold(response_sys)
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinError(ibin,0,0)
                
            Hres_sys[sysname+var+"_"+channel] = response_sys
            
    for thsysname in thsysnames:
        for var in variants:
            f_ttbar_sys = TFile("histfiles_full2016/hists_PowhegPythia8_"+thsysname+var+"_fullTruth_"+channel+"_"+thsysname+var+"_post.root")
            response_sys = f_ttbar_sys.Get(response_name)
            response_sys.Sumw2()
            
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinContent(ibin,0,0)

            Hres_sys[thsysname+var+"_"+channel] = response_sys

    for thsysname in thsysnames:
        if thsysname is "ErdOn" or thsysname is "Herwig":
            f_ttbar_sys = TFile("histfiles_full2016/hists_PowhegPythia8_"+thsysname+"_fullTruth_"+channel+"_"+thsysname+"_post.root")
            response_sys = f_ttbar_sys.Get(response_name)
            response_sys.Sumw2()
            true_sys = f_ttbar_sys.Get(hTrue_name)
            antiTagWeight(true_sys,response_sys)
            convertForTUnfold(response_sys)
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinContent(ibin,0,0)

            Hres_sys[thsysname+"_"+channel] = response_sys

        else :
            for var in variants:
                f_ttbar_sys = TFile("histfiles_full2016/hists_PowhegPythia8_"+thsysname+var+"_fullTruth_"+channel+"_"+thsysname+var+"_post.root")
                response_sys = f_ttbar_sys.Get(response_name)
                response_sys.Sumw2()
                true_sys = f_ttbar_sys.Get(hTrue_name)
                antiTagWeight(true_sys,response_sys)
                convertForTUnfold(response_sys)
                for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                    response_sys.SetBinContent(ibin,0,0)
                    response_sys.SetBinContent(ibin,0,0)

                Hres_sys[thsysname+var+"_"+channel] = response_sys


    # -------------------------------------------------------------------------------------
    # read & normalize histograms
    # -------------------------------------------------------------------------------------

    thisMeas[channel] = f_data.Get(hMeas_name).Clone()
    thisMeas[channel].Sumw2()

    thisTrue_m0to700_p1 = f_ttbar_m0to700_p1.Get(hTrue_name)
    thisTrue_m0to700_p2 = f_ttbar_m0to700_p2.Get(hTrue_name)
    thisTrue_m700to1000 = f_ttbar_m700to1000.Get(hTrue_name)
    thisTrue_m1000toInf = f_ttbar_m1000toInf.Get(hTrue_name)
    thisTrue_m0to700_p1.Sumw2()
    thisTrue_m0to700_p2.Sumw2()
    thisTrue_m700to1000.Sumw2()
    thisTrue_m1000toInf.Sumw2()
    thisTrue_m0to700 = thisTrue_m0to700_p1.Clone()
    thisTrue_m0to700.Add(thisTrue_m0to700_p2)
    thisTrue_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
    thisTrue_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
    thisTrue_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
    thisTrue[channel] = thisTrue_m0to700.Clone()
    thisTrue[channel].Add(thisTrue_m700to1000)
    thisTrue[channel].Add(thisTrue_m1000toInf)

    thisExpect_m0to700_p1 = f_ttbar_m0to700_p1.Get(hMeas_name)
    thisExpect_m0to700_p2 = f_ttbar_m0to700_p2.Get(hMeas_name)
    thisExpect_m700to1000 = f_ttbar_m700to1000.Get(hMeas_name)
    thisExpect_m1000toInf = f_ttbar_m1000toInf.Get(hMeas_name)
    thisExpect_m0to700_p1.Sumw2()
    thisExpect_m0to700_p2.Sumw2()
    thisExpect_m700to1000.Sumw2()
    thisExpect_m1000toInf.Sumw2()
    thisExpect_m0to700 = thisExpect_m0to700_p1.Clone()
    thisExpect_m0to700.Add(thisExpect_m0to700_p2)
    thisExpect_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
    thisExpect_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
    thisExpect_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
    thisExpect[channel] = thisExpect_m0to700.Clone()
    thisExpect[channel].Add(thisExpect_m700to1000)
    thisExpect[channel].Add(thisExpect_m1000toInf)    

    noNegBins(thisMeas[channel])
    noNegBins(thisTrue[channel])
    noNegBins(thisExpect[channel])
    
    thisMeas[channel].SetName("recolevel") 
    thisTrue[channel].SetName("truthlevel")
    
    nbinsTrue = thisTrue[channel].GetNbinsX()
    nbinsMeas = thisMeas[channel].GetNbinsX()

    hMeas_tt_nonsemi          = f_ttbar_nonsemilep.Get(hMeas_name)
    hMeas_T_t                 = f_T_t.Get(hMeas_name)
    hMeas_Tbar_t              = f_Tbar_t.Get(hMeas_name)
    hMeas_T_tW                = f_T_tW.Get(hMeas_name)
    hMeas_Tbar_tW             = f_Tbar_tW.Get(hMeas_name)
    hMeas_T_s                 = f_T_s.Get(hMeas_name)
    hMeas_WJets_HT100to200    = f_WJets_HT100to200.Get(hMeas_name)
    hMeas_WJets_HT200to400    = f_WJets_HT200to400.Get(hMeas_name)
    hMeas_WJets_HT400to600    = f_WJets_HT400to600.Get(hMeas_name)
    hMeas_WJets_HT600to800    = f_WJets_HT600to800.Get(hMeas_name)
    hMeas_WJets_HT800to1200   = f_WJets_HT800to1200.Get(hMeas_name)
    hMeas_WJets_HT1200to2500  = f_WJets_HT1200to2500.Get(hMeas_name)
    hMeas_WJets_HT2500toInf   = f_WJets_HT2500toInf.Get(hMeas_name)
    hMeas_WJetsL_HT100to200   = f_WJets_HT100to200.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT200to400   = f_WJets_HT200to400.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT400to600   = f_WJets_HT400to600.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT600to800   = f_WJets_HT600to800.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT800to1200  = f_WJets_HT800to1200.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT1200to2500 = f_WJets_HT1200to2500.Get(hMeas_name+'_l')
    hMeas_WJetsL_HT2500toInf  = f_WJets_HT2500toInf.Get(hMeas_name+'_l')
    hMeas_ZJets               = f_ZJets.Get(hMeas_name)
    hMeas_WW                  = f_WW.Get(hMeas_name)
    hMeas_WZ                  = f_WZ.Get(hMeas_name)
    hMeas_ZZ                  = f_ZZ.Get(hMeas_name)
    
    hMeas_tt_nonsemi         .Sumw2()
    hMeas_T_t                .Sumw2()
    hMeas_Tbar_t             .Sumw2()
    hMeas_T_tW               .Sumw2()
    hMeas_Tbar_tW            .Sumw2()
    hMeas_T_s                .Sumw2()
    hMeas_WJets_HT100to200   .Sumw2()
    hMeas_WJets_HT200to400   .Sumw2()
    hMeas_WJets_HT400to600   .Sumw2()
    hMeas_WJets_HT600to800   .Sumw2()
    hMeas_WJets_HT800to1200  .Sumw2()
    hMeas_WJets_HT1200to2500 .Sumw2()
    hMeas_WJets_HT2500toInf  .Sumw2()
    hMeas_WJetsL_HT100to200  .Sumw2()
    hMeas_WJetsL_HT200to400  .Sumw2()
    hMeas_WJetsL_HT400to600  .Sumw2()
    hMeas_WJetsL_HT600to800  .Sumw2()
    hMeas_WJetsL_HT800to1200 .Sumw2()
    hMeas_WJetsL_HT1200to2500.Sumw2()
    hMeas_WJetsL_HT2500toInf .Sumw2()
    hMeas_ZJets              .Sumw2()
    hMeas_WW                 .Sumw2()
    hMeas_WZ                 .Sumw2()
    hMeas_ZZ                 .Sumw2()
    
    hMeas_tt_nonsemi.Scale(PowhegPythia8_norm)
    hMeas_T_t.Scale( SingleTop_t_t_norm)
    hMeas_Tbar_t.Scale( SingleTop_tbar_t_norm)
    hMeas_T_tW.Scale(SingleTop_t_tW_norm)
    hMeas_Tbar_tW.Scale(SingleTop_tbar_tW_norm)
    hMeas_T_s.Scale(SingleTop_s_norm)
    hMeas_WJets_HT100to200.Scale(WJets_HT100to200_norm)
    hMeas_WJets_HT200to400.Scale(WJets_HT200to400_norm)
    hMeas_WJets_HT400to600.Scale(WJets_HT400to600_norm)
    hMeas_WJets_HT600to800.Scale(WJets_HT600to800_norm)
    hMeas_WJets_HT800to1200.Scale(WJets_HT800to1200_norm)
    hMeas_WJets_HT1200to2500.Scale(WJets_HT1200to2500_norm)
    hMeas_WJets_HT2500toInf.Scale(WJets_HT2500toInf_norm)
    hMeas_WJetsL_HT100to200.Scale(WJets_HT100to200_norm)
    hMeas_WJetsL_HT200to400.Scale(WJets_HT200to400_norm)
    hMeas_WJetsL_HT400to600.Scale(WJets_HT400to600_norm)
    hMeas_WJetsL_HT600to800.Scale(WJets_HT600to800_norm)
    hMeas_WJetsL_HT800to1200.Scale(WJets_HT800to1200_norm)
    hMeas_WJetsL_HT1200to2500.Scale(WJets_HT1200to2500_norm)
    hMeas_WJetsL_HT2500toInf.Scale(WJets_HT2500toInf_norm)
    hMeas_ZJets.Scale(ZJets_norm)
    hMeas_WW.Scale(WW_norm)
    hMeas_WZ.Scale(WZ_norm)
    hMeas_ZZ.Scale(ZZ_norm)
    
    hMeas_qcd                    = f_QCD.Get(hMeas_name)
    hMeas_qcd_tt                 = f_qcd_ttbar.Get(hMeas_name)
    hMeas_qcd_tt_nonsemi         = f_qcd_ttbar_nonsemilep.Get(hMeas_name)
    hMeas_qcd_T_t                = f_qcd_T_t.Get(hMeas_name)
    hMeas_qcd_Tbar_t             = f_qcd_Tbar_t.Get(hMeas_name)
    hMeas_qcd_T_tW               = f_qcd_T_tW.Get(hMeas_name)
    hMeas_qcd_Tbar_tW            = f_qcd_Tbar_tW.Get(hMeas_name)
    hMeas_qcd_T_s                = f_qcd_T_s.Get(hMeas_name)
    hMeas_qcd_WJets_HT100to200   = f_qcd_WJets_HT100to200.Get(hMeas_name)
    hMeas_qcd_WJets_HT200to400   = f_qcd_WJets_HT200to400.Get(hMeas_name)
    hMeas_qcd_WJets_HT400to600   = f_qcd_WJets_HT400to600.Get(hMeas_name)
    hMeas_qcd_WJets_HT600to800   = f_qcd_WJets_HT600to800.Get(hMeas_name)
    hMeas_qcd_WJets_HT800to1200  = f_qcd_WJets_HT800to1200.Get(hMeas_name)
    hMeas_qcd_WJets_HT1200to2500 = f_qcd_WJets_HT1200to2500.Get(hMeas_name)
    hMeas_qcd_WJets_HT2500toInf  = f_qcd_WJets_HT2500toInf.Get(hMeas_name)
    hMeas_qcd_ZJets              = f_qcd_ZJets.Get(hMeas_name)
    hMeas_qcd_WW                 = f_qcd_WW.Get(hMeas_name)
    hMeas_qcd_WZ                 = f_qcd_WZ.Get(hMeas_name)
    hMeas_qcd_ZZ                 = f_qcd_ZZ.Get(hMeas_name)
    
    hMeas_qcd_tt_nonsemi         .Sumw2()
    hMeas_qcd_T_t                .Sumw2()
    hMeas_qcd_Tbar_t             .Sumw2()
    hMeas_qcd_T_tW               .Sumw2()
    hMeas_qcd_Tbar_tW            .Sumw2()
    hMeas_qcd_T_s                .Sumw2()
    hMeas_qcd_WJets_HT100to200   .Sumw2()
    hMeas_qcd_WJets_HT200to400   .Sumw2()
    hMeas_qcd_WJets_HT400to600   .Sumw2()
    hMeas_qcd_WJets_HT600to800   .Sumw2()
    hMeas_qcd_WJets_HT800to1200  .Sumw2()
    hMeas_qcd_WJets_HT1200to2500 .Sumw2()
    hMeas_qcd_WJets_HT2500toInf  .Sumw2()
    hMeas_qcd_ZJets              .Sumw2()
    hMeas_qcd_WW                 .Sumw2()
    hMeas_qcd_WZ                 .Sumw2()
    hMeas_qcd_ZZ                 .Sumw2()
    
    hMeas_qcd_tt.Scale(PowhegPythia8_norm)
    hMeas_qcd_tt_nonsemi.Scale(PowhegPythia8_norm)
    hMeas_qcd_T_t.Scale( SingleTop_t_t_norm)
    hMeas_qcd_Tbar_t.Scale( SingleTop_tbar_t_norm)
    hMeas_qcd_T_tW.Scale(SingleTop_t_tW_norm)
    hMeas_qcd_Tbar_tW.Scale(SingleTop_tbar_tW_norm)
    hMeas_qcd_T_s.Scale(SingleTop_s_norm)
    hMeas_qcd_WJets_HT100to200.Scale(WJets_HT100to200_norm)
    hMeas_qcd_WJets_HT200to400.Scale(WJets_HT200to400_norm)
    hMeas_qcd_WJets_HT400to600.Scale(WJets_HT400to600_norm)
    hMeas_qcd_WJets_HT600to800.Scale(WJets_HT600to800_norm)
    hMeas_qcd_WJets_HT800to1200.Scale(WJets_HT800to1200_norm)
    hMeas_qcd_WJets_HT1200to2500.Scale(WJets_HT1200to2500_norm)
    hMeas_qcd_WJets_HT2500toInf.Scale(WJets_HT2500toInf_norm)
    hMeas_qcd_ZJets.Scale(ZJets_norm)
    hMeas_qcd_WW.Scale(WW_norm)
    hMeas_qcd_WZ.Scale(WZ_norm)
    hMeas_qcd_ZZ.Scale(ZZ_norm)

    backgrounds[channel] = []
    bkg_TT_nonsemi = Background("t#bar{t} other",0.80,0.08,632-7)
    bkg_TT_nonsemi.addHist(hMeas_tt_nonsemi)
    backgrounds[channel].append(bkg_TT_nonsemi)
    
    bkg_SingleTop = Background("Single t",1.32,0.57,6)
    for hist in [hMeas_T_t,hMeas_Tbar_t, hMeas_T_tW, hMeas_Tbar_tW, hMeas_T_s] :
        bkg_SingleTop.addHist(hist)
    backgrounds[channel].append(bkg_SingleTop)

    bkg_WJetsL = Background("W+jets LF",0.78,0.18,416)
    bkg_WJetsHF = Background("W+jets HF",1.05,0.30,416-3)
    for hist in [hMeas_WJets_HT100to200,hMeas_WJets_HT200to400,hMeas_WJets_HT400to600,hMeas_WJets_HT600to800,hMeas_WJets_HT800to1200,hMeas_WJets_HT1200to2500,hMeas_WJets_HT2500toInf] :
        bkg_WJetsHF.addHist( hist )
    for hist in [hMeas_WJetsL_HT100to200,hMeas_WJetsL_HT200to400,hMeas_WJetsL_HT400to600,hMeas_WJetsL_HT600to800,hMeas_WJetsL_HT800to1200,hMeas_WJetsL_HT1200to2500,hMeas_WJetsL_HT2500toInf] :
        bkg_WJetsL.addHist( hist )
        hist.Scale(-1.0)
        bkg_WJetsHF.addHist( hist )
    backgrounds[channel].append(bkg_WJetsL)
    backgrounds[channel].append(bkg_WJetsHF)

    bkg_ZJets = Background("Z+jets",0.67,0.24,860-9)
    bkg_ZJets.addHist(hMeas_ZJets)
    backgrounds[channel].append(bkg_ZJets)

    bkg_Diboson = Background("Diboson",0.99,0.31,880-6)
    for hist in [hMeas_WW, hMeas_WZ, hMeas_ZZ]:
        bkg_Diboson.addHist(hist)
    backgrounds[channel].append(bkg_Diboson)

    hMeas_QCD = hMeas_qcd.Clone()
    for hist in [hMeas_qcd_tt, hMeas_qcd_tt_nonsemi, hMeas_qcd_T_t,hMeas_qcd_Tbar_t,hMeas_qcd_T_tW,hMeas_qcd_Tbar_tW, hMeas_qcd_T_s, hMeas_qcd_WJets_HT100to200,hMeas_qcd_WJets_HT200to400,hMeas_qcd_WJets_HT400to600, hMeas_qcd_WJets_HT600to800,hMeas_qcd_WJets_HT800to1200,hMeas_qcd_WJets_HT1200to2500,hMeas_qcd_WJets_HT2500toInf, hMeas_qcd_ZJets, hMeas_qcd_WW, hMeas_qcd_WZ, hMeas_qcd_ZZ] :
        hMeas_QCD.Add(hist,-1.0)
    noNegBins(hMeas_QCD)
    
    # -------------------------------
    # Normalize QCD to MC prediction
    # -------------------------------

    f_QCD_HT500to700   = TFile("histfiles_full2016/hists_QCD_HT500to700_"+channel+"_nom_post.root")
    f_QCD_HT700to1000  = TFile("histfiles_full2016/hists_QCD_HT700to1000_"+channel+"_nom_post.root")
    f_QCD_HT1000to1500 = TFile("histfiles_full2016/hists_QCD_HT1000to1500_"+channel+"_nom_post.root")
    f_QCD_HT1500to2000 = TFile("histfiles_full2016/hists_QCD_HT1500to2000_"+channel+"_nom_post.root")
    f_QCD_HT2000toInf  = TFile("histfiles_full2016/hists_QCD_HT2000toInf_"+channel+"_nom_post.root")

    hNorm_QCD_HT500to700   = f_QCD_HT500to700.Get(hMeas_name)
    hNorm_QCD_HT700to1000  = f_QCD_HT700to1000.Get(hMeas_name)
    hNorm_QCD_HT1000to1500 = f_QCD_HT1000to1500.Get(hMeas_name)
    hNorm_QCD_HT1500to2000 = f_QCD_HT1500to2000.Get(hMeas_name)
    hNorm_QCD_HT2000toInf  = f_QCD_HT2000toInf.Get(hMeas_name)
    
    hNorm_QCD_HT500to700   .Sumw2()
    hNorm_QCD_HT700to1000  .Sumw2()
    hNorm_QCD_HT1000to1500 .Sumw2()
    hNorm_QCD_HT1500to2000 .Sumw2()
    hNorm_QCD_HT2000toInf  .Sumw2()
    
    hNorm_QCD_HT500to700.Scale(QCD_HT500to700_norm)
    hNorm_QCD_HT700to1000.Scale(QCD_HT700to1000_norm)
    hNorm_QCD_HT1000to1500.Scale(QCD_HT1000to1500_norm)
    hNorm_QCD_HT1500to2000.Scale(QCD_HT1500to2000_norm)
    hNorm_QCD_HT2000toInf.Scale(QCD_HT2000toInf_norm)
    
    QCD_norm = 0.0
    for hist in [hNorm_QCD_HT500to700,hNorm_QCD_HT700to1000,hNorm_QCD_HT1000to1500,hNorm_QCD_HT1500to2000,hNorm_QCD_HT2000toInf]:
        QCD_norm += hist.Integral()
    
    hMeas_QCD.Scale(QCD_norm / hMeas_QCD.Integral())

    if channel is "mu":
        bkg_QCD = Background("muQCD",0.71,0.75,400)
    else:
        bkg_QCD = Background("elQCD",0.90,0.65,400)
    bkg_QCD.addHist(hMeas_QCD)
    backgrounds[channel].append(bkg_QCD)

    # -------------------------------
    # Construct total background hist, with posterior normalizations
    # -------------------------------
    
    hBkg = thisMeas[channel].Clone()
    hBkg.Reset()

    for background in backgrounds[channel]:
        hBkg.Add(background.hist,background.norm)

    # -------------------------------
    # Troubleshoot -- print and plot comparison of data vs. prediction
    # -------------------------------

    print "Bin             & Data             & $tt nonsemi     & SingleTop      & WJetsL         & WJetsHF        & ZJets          & Diboson        & QCD           & Fake fraction "
    for ibin in xrange(1,nbinsMeas+1):
        print string.ljust("$[{:.1f},{:.1f}]$".format(thisMeas[channel].GetXaxis().GetBinLowEdge(ibin), thisMeas[channel].GetXaxis().GetBinUpEdge(ibin)),16) + "&" + string.rjust("{:.1f}".format(thisMeas[channel].GetBinContent(ibin)),6) + " $\pm$ " + string.rjust("{:.1f}".format(thisMeas[channel].GetBinError(ibin)),4) + " & " + " & ".join(" $\pm$ ".join([string.rjust("{:.1f}".format(background.hist.GetBinContent(ibin)*background.norm),4),string.rjust("{:.1f}".format(background.hist.GetBinError(ibin)*background.norm),3)]) for background in backgrounds[channel]) + " & " + "{:.3f}".format(response[channel].GetBinContent(ibin,0) / response[channel].Integral(ibin,ibin,0,response[channel].GetYaxis().GetNbins()+1))

    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.14)
    
    c7 = TCanvas("c7","",900,800)
    hStack = THStack("hStack","")
    for background in backgrounds[channel]:
        background.hist.Scale(background.norm) #Will need to undo this later
        background.hist.SetFillColor(background.color)
        hStack.Add(background.hist)
    thisExpect[channel].Scale(0.80)
    thisExpect[channel].SetFillColor(632+1)
    hStack.Add(thisExpect[channel])
    hMC = hBkg.Clone()
    hMC.Add(thisExpect[channel])
    thisMeas[channel].SetBinErrorOption(TH1.kPoisson)

    h_ratio = thisMeas[channel].Clone()
    h_ratio.Divide(hMC)
    h_ratio2 = hMC.Clone()
    for ib in xrange(0,hMC.GetNbinsX()):
        tmpval = 0.0 if h_ratio2.GetBinContent(ib+1) == 0.0 else h_ratio2.GetBinError(ib+1)/h_ratio2.GetBinContent(ib+1)
        h_ratio2.SetBinError(ib+1,tmpval)
        h_ratio2.SetBinContent(ib+1,1.0)

    p1 = TPad("p1","",0,0.3,1,1)
    p1.SetTopMargin(0.08)
    p1.SetBottomMargin(0.05)
    p1.SetNumber(1)
    p2 = TPad("p2","",0,0,1,0.3)
    p2.SetTopMargin(0.05)
    p2.SetBottomMargin(0.35)
    p2.SetNumber(2)
    
    p1.Draw()
    p2.Draw()
    p1.cd()

    thisMeas[channel].SetMaximum(max(thisMeas[channel].GetMaximum(),hMC.GetMaximum())*1.3)
    thisMeas[channel].SetMarkerStyle(8)
    thisMeas[channel].SetMarkerSize(1)
    #thisMeas[channel].GetYaxis().SetLabelSize(26)
    thisMeas[channel].SetTitle("")
    
    thisMeas[channel].Draw("LE0P")  
    hStack.Draw("hist,same")
    thisMeas[channel].Draw("LE0P,same")

    leg7 = TLegend(0.71,0.52,0.91,0.92)
    leg7.SetBorderSize(0)
    leg7.SetFillStyle(0)
    leg7.SetTextFont(42)
    leg7.SetTextSize(0.045)
    leg7.AddEntry(thisMeas[channel],"Data","pe")
    leg7.AddEntry(thisExpect[channel],"t#bar{t} signal","f")
    for background in backgrounds[channel] :
        if "QCD" in background.name :
            leg7.AddEntry(background.hist,"Multijet","f")
        else :
            leg7.AddEntry(background.hist,background.name,"f")
    leg7.AddEntry(h_ratio2, "MC Stat. Unc.", "f")
    leg7.Draw()

    myText(0.10,0.94,1,"#intLdt = 35.9 fb^{-1}")
    myText(0.80,0.94,1,"#sqrt{s} = 13 TeV")

    p2.cd()
    p2.SetGridy()
    
    h_ratio.SetMarkerStyle(8)
    h_ratio.SetMarkerSize(1)
    h_ratio.SetMaximum(1.8)
    h_ratio.SetMinimum(0.2)
    h_ratio.GetYaxis().SetNdivisions(2,4,0,False)
    h_ratio.GetYaxis().SetTitle("Data / MC")
    h_ratio.GetXaxis().SetTitle(thisMeas[channel].GetXaxis().GetTitle())
    #h_ratio.GetXaxis().SetLabelSize(26)
    #h_ratio.GetYaxis().SetLabelSize(26)
    #h_ratio.GetXaxis().SetTitleOffset(2.8)
    #h_ratio.GetYaxis().SetTitleOffset(1.4)
    #h_ratio.GetXaxis().SetTitleSize(32)
    #h_ratio.GetYaxis().SetTitleSize(32)

    h_ratio2.SetMarkerSize(0)
    h_ratio2.SetLineColor(0)
    h_ratio2.SetFillColor(15)
    h_ratio2.SetFillStyle(1001)

    h_ratio.Draw("le0p")
    h_ratio2.Draw("same,e2")
    h_ratio.Draw("le0p,same")
    
    c7.SaveAs("UnfoldingPlots/compare_"+options.toUnfold+"RecoTop_"+channel+".pdf")

    for background in backgrounds[channel]:
        background.hist.Scale(1.0/background.norm)

    # -------------------------------
    # Convert for TUnfold
    # -------------------------------
    
    # Get 'ttbar' part of thisMeas 
    thisMeas[channel].Add(hBkg,-1.0)

    antiTagWeight(thisTrue[channel],response[channel])
    convertForTUnfold(response[channel])
    removeFakes(thisMeas[channel],response[channel])
            
    # Now add back background contribution to thisMeas
    thisMeas[channel].Add(hBkg)

# -------------------------------------------------------------------------------------
# Construct combined response matrices / inputs
# -------------------------------------------------------------------------------------
thisMeas["comb"] = thisMeas["mu"].Clone()
thisMeas["comb"].Add(thisMeas["el"])

thisTrue["comb"] = thisTrue["mu"].Clone()
thisTrue["comb"].Add(thisTrue["el"])

thisExpect["comb"] = thisExpect["mu"].Clone()
thisExpect["comb"].Add(thisExpect["el"])

backgrounds["comb"] = []
for ibkg in xrange(0,len(backgrounds["mu"])-1):
    bkg_comb = Background(backgrounds["mu"][ibkg].name, backgrounds["mu"][ibkg].norm, backgrounds["mu"][ibkg].err, backgrounds["mu"][ibkg].color)
    bkg_comb.addHist(backgrounds["mu"][ibkg].hist)
    bkg_comb.addHist(backgrounds["el"][ibkg].hist)
    backgrounds["comb"].append(bkg_comb)
backgrounds["comb"].append(backgrounds["mu"][len(backgrounds["mu"])-1])
backgrounds["comb"].append(backgrounds["el"][len(backgrounds["el"])-1])

response["comb"] = response["mu"].Clone()
response["comb"].Add(response["el"])

for sysname in allsysnames:
    for var in variants:
        Hres_sys[sysname+var+"_comb"] = Hres_sys[sysname+var+"_mu"].Clone()
        Hres_sys[sysname+var+"_comb"].Add(Hres_sys[sysname+var+"_el"])

# -------------------------------------------------------------------------------------
# Do unfolding
# -------------------------------------------------------------------------------------

channels = ["mu","el","comb"]

for channel in channels:
    unfold = {}
    for var in variants:
        unfold[var] = TUnfoldDensity(response[channel],TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintArea, TUnfoldDensity.kDensityModeBinWidth)
        unfold[var].SetInput(thisMeas[channel])

        # Add systematic uncertainties
        for sysname in allsysnames:
            if sysname == "ErdOn" or sysname == "Herwig":
                unfold[var].AddSysError(Hres_sys[sysname+"_"+channel],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            else :
                unfold[var].AddSysError(Hres_sys[sysname+var+"_"+channel],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            
        # Subtract backgrounds with appropriate scale / uncertainty
        for background in backgrounds[channel]:
            unfold[var].SubtractBackground(background.hist,background.name,background.norm,background.err)

        # Do unfolding
        unfold[var].DoUnfold(0)

        print "chi**2=" + str(unfold[var].GetChi2A()) + "+" + str(unfold[var].GetChi2L()) + " / " + str(unfold[var].GetNdf())
    
    # unfolded distribution (histogram)
    thisReco = unfold["Up"].GetOutput("reco") #Same for both, so we just pick one
    thisReco.Sumw2()

    '''
    xvec = TMatrixD(nbinsTrue,1)
    for ibin in xrange(0,nbinsTrue):
        xvec[ibin][0] = thisReco.GetBinContent(ibin+1)

    yvec = TMatrixD(nbinsMeas,1)
    for ibin in xrange(0,nbinsMeas):
        tmp = thisMeas[channel].GetBinContent(ibin+1)
        for background in backgrounds[channel]:
            tmp -= background.hist.GetBinContent(ibin+1)
        yvec[ibin][0] = tmp
    '''
    # -------------------------------------------------------------------------------------
    #Plot error breakdown
    # -------------------------------------------------------------------------------------
    
    #Statistical -- input and unfolding matrix (GetEmatrixSysUncorr() and GetEmatrixInput())
    h_STAT = thisTrue[channel].Clone("stat")
    h_STAT.Reset()

    h_EXP = thisTrue[channel].Clone("stat")
    h_EXP.Reset()
    
    h_TH = thisTrue[channel].Clone("stat")
    h_TH.Reset()

    h_TOT = thisTrue[channel].Clone("tot")
    h_TOT.Reset()
    
    # Individual error sources
    h_SYS = {}
    for sysname in allsysnames:
        h_SYS[sysname] = thisTrue[channel].Clone(sysname)
        h_SYS[sysname].Reset()
        
    h_BKG = thisTrue[channel].Clone()
    h_BKG.Reset()

    h_BSCALE = thisTrue[channel].Clone()
    h_BSCALE.Reset()
    
    h_BSTAT = thisTrue[channel].Clone()
    h_BSTAT.Reset()
    
    h_INPUT = thisTrue[channel].Clone()
    h_INPUT.Reset()
    
    h_MATRIX = thisTrue[channel].Clone()
    h_MATRIX.Reset()
    
    h_LUMI = thisTrue[channel].Clone()
    h_LUMI.Reset()

    #Get actual error matrices / hists
    hErrInput = unfold["Up"].GetEmatrixInput("mErrInput")
    hErrStat = unfold["Up"].GetEmatrixSysUncorr("mErrStat")
    
    # Construct average systematic errors; get FSR covariance matrix for later
    hErrSys = {}
    hCovSysFSR = hErrInput.Clone()
    hCovSysFSR.Reset()
    for sysname in allsysnames:
        hErrSysUp = unfold["Up"].GetDeltaSysSource(sysname,"hErrSysUp_"+sysname)
        hErrSysDn =  unfold["Down"].GetDeltaSysSource(sysname,"hErrSysDown_"+sysname)
        hErrSys[sysname] = hErrSysUp.Clone()
        hErrSys[sysname].Reset()
        for ibinx in xrange(1,nbinsTrue+1):
            hErrSys[sysname].SetBinContent(ibinx,(pow(hErrSysUp.GetBinContent(ibinx),2)+pow(hErrSysDn.GetBinContent(ibinx),2))/2.0)
            if sysname is "FSR":
                for ibiny in xrange(1,nbinsTrue+1):
                    hCovSysFSR.SetBinContent(ibinx,ibiny,(hErrSysUp.GetBinContent(ibinx)*hErrSysUp.GetBinContent(ibiny)+hErrSysDn.GetBinContent(ibinx)*hErrSysDn.GetBinContent(ibiny))/2.0)

    hErrBkgStat = {}
    hErrBkgScale = {}
    for background in backgrounds[channel]:
        hErrBkgStat[background.name] = unfold["Up"].GetEmatrixSysBackgroundUncorr(background.name,"mErrBkgStat_"+background.name)
        hErrBkgScale[background.name] = unfold["Up"].GetDeltaSysBackgroundScale(background.name,"mErrBkgScale_"+background.name)

    # Average total covariance matrices
    hErrTotUp = unfold["Up"].GetEmatrixTotal("mErrTotalUp")
    hErrTotDn = unfold["Down"].GetEmatrixTotal("mErrTotalDown")
    hErrTot = hErrTotUp.Clone()
    hErrTot.Add(hErrTotDn)
    hErrTot.Scale(0.5)
    
    # Correct FSR
    if "FSR" in allsysnames:
        hErrSys["FSR"].Scale(0.5); #Scale uncertainty from FSR down by sqrt(2) --> scale uncertainty squared by 0.5
        hErrTot.Add(hCovSysFSR,-0.5) #Contribution from FSR to total covariance decreases by 0.5       
    
    # Fill uncertainty histograms
    for ibin in xrange(1,nbinsTrue+1):
        h_STAT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)+hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_INPUT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_MATRIX.SetBinContent(ibin,math.sqrt(hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_TOT.SetBinContent(ibin,math.sqrt(hErrTot.GetBinContent(ibin,ibin)/pow(thisReco.GetBinContent(ibin),2)+0.025*0.025))
        h_LUMI.SetBinContent(ibin,0.025)
        tot_exp = 0.0
        tot_th = 0.0
        for sysname in allsysnames:
            h_SYS[sysname].SetBinContent(ibin,math.sqrt(hErrSys[sysname].GetBinContent(ibin))/thisReco.GetBinContent(ibin))
        for sysname in ["PDF","Q2","ISR","FSR","Tune","Hdamp","ErdOn","Herwig"]:
            tot_th += hErrSys[sysname].GetBinContent(ibin)
        for sysname in ["JEC","JER","BTag","TopTag","lep","pu"]:
            tot_exp += hErrSys[sysname].GetBinContent(ibin)
        tot_bkg_stat = 0.0
        tot_bkg_scale = 0.0
        for background in backgrounds[channel]:
            tot_bkg_stat += hErrBkgStat[background.name].GetBinContent(ibin,ibin)
            tot_bkg_scale += pow(hErrBkgScale[background.name].GetBinContent(ibin),2)
        h_BSTAT.SetBinContent(ibin,math.sqrt(tot_bkg_stat)/thisReco.GetBinContent(ibin))
        h_BSCALE.SetBinContent(ibin,math.sqrt(tot_bkg_scale)/thisReco.GetBinContent(ibin))
        h_BKG.SetBinContent(ibin,math.sqrt(tot_bkg_stat+tot_bkg_scale)/thisReco.GetBinContent(ibin))
        tot_exp += tot_bkg_stat
        tot_exp += tot_bkg_scale
        h_EXP.SetBinContent(ibin,math.sqrt(tot_exp)/thisReco.GetBinContent(ibin))
        h_TH.SetBinContent(ibin,math.sqrt(tot_th)/thisReco.GetBinContent(ibin))

    # Plot error breakdown
    c6 = TCanvas("c6", "", 800, 600)
    c6.SetTopMargin(0.08)
    c6.SetRightMargin(0.05)
    c6.SetBottomMargin(0.14)
    c6.SetLeftMargin(0.16)

    h_TOT.GetXaxis().SetTitle("Top "+labelstring1+" "+labelstring2)
    h_TOT.GetYaxis().SetTitle("Uncertainty [%]")
    if options.toUnfold == "pt":
        h_TOT.SetAxisRange(400,1199,"X")
    
    h_TOT.GetYaxis().SetTitleSize(0.055)    
    h_TOT.GetYaxis().SetTitleOffset(1.1)
    h_TOT.GetYaxis().SetLabelSize(0.045)
    
    h_TOT.GetXaxis().SetTitleSize(0.05)
    h_TOT.GetXaxis().SetTitleOffset(1.2)
    h_TOT.GetXaxis().SetLabelSize(0.0455)
    
    c6.cd()
    
    h_TOT.SetFillColor(17)
    h_TOT.SetFillStyle(3344)
    h_TOT.SetLineColor(16)
    h_TOT.SetLineWidth(2)
    
    h_STAT.SetLineColor(1)
    h_STAT.SetLineWidth(2)
    h_STAT.SetMarkerColor(1)
    h_STAT.SetMarkerStyle(20)
    
    h_LUMI.SetLineColor(40)
    h_LUMI.SetLineWidth(2)
    h_LUMI.SetMarkerColor(40)
    h_LUMI.SetMarkerStyle(34)

    colors = [632,600,617,417,432,4,1,419,600,882,632,600,617,2]
    markers = [20,21,22,23,33,26,24,25,27,32,23,33,26,24]
    for isys in xrange(0,len(allsysnames)):
        h_SYS[allsysnames[isys]].SetLineColor(colors[isys])
        h_SYS[allsysnames[isys]].SetLineWidth(2)
        h_SYS[allsysnames[isys]].SetMarkerColor(colors[isys])
        h_SYS[allsysnames[isys]].SetMarkerStyle(markers[isys])
        
    h_BKG.SetLineColor(29)
    h_BKG.SetLineWidth(2)
    h_BKG.SetMarkerColor(29)
    h_BKG.SetMarkerStyle(22)
  
    leg6 = TLegend(0.2,0.39,0.45,0.88)
    leg6.AddEntry(h_TOT,"Total syst. uncertainty","f")
    leg6.AddEntry(h_STAT,"Input stat. unc.","lp")
    leg6.AddEntry(h_LUMI,"Int. luminosity","lp")
    for isys in xrange(0,len(allsysnames)):
        leg6.AddEntry(h_SYS[allsysnames[isys]],longnames[isys],"lp")
    leg6.AddEntry(h_BKG,"Backgrounds","lp")
        
    leg6.SetFillStyle(0)
    leg6.SetBorderSize(0)
    leg6.SetTextSize(0.04)
    leg6.SetTextFont(42)

    h_TOT.GetYaxis().SetRangeUser(0.0,1.0)
    h_TOT.Draw("hist")
    h_STAT.Draw("ep,same")
    h_LUMI.Draw("ep,same")
    for sysname in allsysnames:
        h_SYS[sysname].Draw("ep,same")
    h_BKG.Draw("ep,same")
        
    leg6.Draw() 

    drawCMS(0.17,0.935)
    
    c6.SaveAs("UnfoldingPlots/unfold_relative_uncertainties_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")

    # -------------------------------------------------------------------------------------
    # Troubleshoot stat. unc.
    # -------------------------------------------------------------------------------------

    # Get Vyy (manually)
    '''
    Vyy_input = TMatrixD(nbinsMeas,nbinsMeas)
    for ibin in xrange(0,nbinsMeas): #rows
        Vyy_input[ibin][ibin] = pow(thisMeas[channel].GetBinError(ibin+1),2)

    Vyy_bstat = TMatrixD(nbinsMeas,nbinsMeas)
    Vyy_bscale = TMatrixD(nbinsMeas,nbinsMeas)

    for ibinx in xrange(0,nbinsMeas): #rows
        sumbkg = 0.0
        row = '[{:.1f},{:.1f}] '.format(thisMeas[channel].GetXaxis().GetBinLowEdge(ibinx+1),thisMeas[channel].GetXaxis().GetBinUpEdge(ibinx+1))
        for background in backgrounds[channel]:
            sumbkg += pow(background.norm*background.hist.GetBinError(ibinx+1),2)
            row += ' & {:.1f}'.format(100.0*background.norm*background.hist.GetBinError(ibinx+1)/(thisMeas[channel].GetBinContent(ibinx+1)-hBkg.GetBinContent(ibinx+1)))
        Vyy_bstat[ibinx][ibinx] = sumbkg
        row += ' & {:.1f}'.format(100.0*math.sqrt(sumbkg)/(thisMeas[channel].GetBinContent(ibinx+1)-hBkg.GetBinContent(ibinx+1)))
            
        for ibiny in xrange(0,nbinsMeas): #cols
            sumbkg2 = 0.0
            for background in backgrounds[channel]:
                sumbkg2 += pow(background.err,2)*background.hist.GetBinContent(ibinx+1)*background.hist.GetBinContent(ibiny+1)
                if ibiny == ibinx:
                    row += ' & {:.1f}'.format(100.0*background.err*math.sqrt(background.hist.GetBinContent(ibinx+1)*background.hist.GetBinContent(ibiny+1))/(thisMeas[channel].GetBinContent(ibinx+1)-hBkg.GetBinContent(ibinx+1)))
            Vyy_bscale[ibinx][ibiny] = sumbkg2
            if ibiny == ibinx:
                row += ' & {:.1f}'.format(100.0*math.sqrt(sumbkg2)/(thisMeas[channel].GetBinContent(ibinx+1)-hBkg.GetBinContent(ibinx+1)))
                row += ' & {:.1f}'.format(100.0*math.sqrt(sumbkg+sumbkg2)/(thisMeas[channel].GetBinContent(ibinx+1)-hBkg.GetBinContent(ibinx+1)))
        print row

    Vyy = TMatrixD(Vyy_input)
    Vyy += Vyy_bstat
    Vyy += Vyy_bscale

    VyyInv = TMatrixD(Vyy)
    VyyInv.Invert()

    #Get A (normalized response matrix)
    A = TMatrixD(nbinsMeas,nbinsTrue)
    for irow in xrange(0,nbinsMeas):
        for icol in xrange(0,nbinsTrue):
            A[irow][icol] = response[channel].GetBinContent(irow+1,icol+1) / response[channel].Integral(0,nbinsMeas+1,icol+1,icol+1)
            
    #Manually calculate Vxx
    AT = TMatrixD(nbinsTrue,nbinsMeas)
    AT.Transpose(A)
    
    VxxInv1 = TMatrixD(nbinsTrue,nbinsMeas)
    VxxInv1.Mult(AT,VyyInv)
    VxxInv = TMatrixD(nbinsTrue,nbinsTrue)
    VxxInv.Mult(VxxInv1,A)
    Vxx = TMatrixD(VxxInv)
    Vxx.Invert()

    #Plot matrices
    hVyy_input  = TH2F("hVyy_input","",nbinsMeas,0,nbinsMeas,nbinsMeas,0,nbinsMeas)
    hVyy_bstat  = TH2F("hVyy_bstat","",nbinsMeas,0,nbinsMeas,nbinsMeas,0,nbinsMeas)
    hVyy_bscale = TH2F("hVyy_bscale","",nbinsMeas,0,nbinsMeas,nbinsMeas,0,nbinsMeas)
    hVyy        = TH2F("hVyy","",nbinsMeas,0,nbinsMeas,nbinsMeas,0,nbinsMeas)
    hVyyInv     = TH2F("hVyyInv","",nbinsMeas,0,nbinsMeas,nbinsMeas,0,nbinsMeas)
    hVxxInv     = TH2F("hVxxInv","",nbinsTrue,0,nbinsTrue,nbinsTrue,0,nbinsTrue)
    hVxx        = TH2F("hVxx","",nbinsTrue,0,nbinsTrue,nbinsTrue,0,nbinsTrue)
    
    for ibinx in xrange(0,nbinsMeas):
        hVyy_input.SetBinContent(ibinx+1,ibinx+1,Vyy_input[ibinx][ibinx])
        hVyy_bstat.SetBinContent(ibinx+1,ibinx+1,Vyy_bstat[ibinx][ibinx])
        for ibiny in xrange(0,nbinsMeas):
            hVyy_bscale.SetBinContent(ibinx+1,ibiny+1,Vyy_bscale[ibinx][ibiny])
            hVyy.SetBinContent(ibinx+1,ibiny+1,Vyy[ibinx][ibiny])
            hVyyInv.SetBinContent(ibinx+1,ibiny+1,VyyInv[ibinx][ibiny])
    for ibinx in xrange(0,nbinsTrue):
        for ibiny in xrange(0,nbinsTrue):
            hVxxInv.SetBinContent(ibinx+1,ibiny+1,VxxInv[ibinx][ibiny])
            hVxx.SetBinContent(ibinx+1,ibiny+1,Vxx[ibinx][ibiny])
                    
    gStyle.SetPadTopMargin(0.06)
    gStyle.SetPadRightMargin(0.15)
    gStyle.SetPadBottomMargin(0.06)
    gStyle.SetPadLeftMargin(0.06)

    gStyle.SetPalette(55)
    
    c9 = TCanvas()
    hVyy_input.SetMaximum(1000.)
    hVyy_input.Draw("colz")
    c9.SaveAs("UnfoldingPlots/covariance_input_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")
    hVyy_bstat.SetMaximum(1000.)
    hVyy_bstat.Draw("colz")
    c9.SaveAs("UnfoldingPlots/covariance_bstat_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")
    hVyy_bscale.SetMaximum(1000.)
    hVyy_bscale.Draw("colz")
    c9.SaveAs("UnfoldingPlots/covariance_bscale_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")    
    #hVyy.SetMaximum(1000.)
    hVyy.Draw("colz")
    c9.SaveAs("UnfoldingPlots/covariance_totalVyy_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")    
    hVyyInv.Draw("colz")
    c9.SaveAs("UnfoldingPlots/VyyInv_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")    
    hVxxInv.Draw("colz")
    c9.SaveAs("UnfoldingPlots/VxxInv_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")    
    hVxx.Draw("colz")
    c9.SaveAs("UnfoldingPlots/Vxx_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")    
    
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadBottomMargin(0.16)
    gStyle.SetPadLeftMargin(0.18)

    #Plot comparisons of pre / post unfolding relative uncertainties
    hVyyInput = thisMeas[channel].Clone()
    hVyyInput.Reset()
    hVyyBstat = thisMeas[channel].Clone()
    hVyyBstat.Reset()
    hVyyBscale = thisMeas[channel].Clone()
    hVyyBscale.Reset()
    hVyyTotal = thisMeas[channel].Clone()
    hVyyTotal.Reset()
    hVxxDiag = thisTrue[channel].Clone()
    hVxxDiag.Reset()

    for ibin in xrange(0,nbinsMeas):
        hVyyInput.SetBinContent(ibin+1,math.sqrt(Vyy_input[ibin][ibin])/(thisMeas[channel].GetBinContent(ibin+1)-hBkg.GetBinContent(ibin+1)))
        hVyyBstat.SetBinContent(ibin+1,math.sqrt(Vyy_bstat[ibin][ibin])/(thisMeas[channel].GetBinContent(ibin+1)-hBkg.GetBinContent(ibin+1)))
        hVyyBscale.SetBinContent(ibin+1,math.sqrt(Vyy_bscale[ibin][ibin])/(thisMeas[channel].GetBinContent(ibin+1)-hBkg.GetBinContent(ibin+1)))
        hVyyTotal.SetBinContent(ibin+1,math.sqrt(Vyy[ibin][ibin])/(thisMeas[channel].GetBinContent(ibin+1)-hBkg.GetBinContent(ibin+1)))
    
    for ibin in xrange(0,nbinsTrue):
        hVxxDiag.SetBinContent(ibin+1,math.sqrt(Vxx[ibin][ibin])/thisReco.GetBinContent(ibin+1))

    hVyyInput.SetLineColor(2)
    hVyyBstat.SetLineColor(2)
    hVyyBscale.SetLineColor(2)
    hVyyTotal.SetLineColor(2)

    c8 = TCanvas()
    h_INPUT.SetMinimum(0.0)
    h_INPUT.SetTitle(";top "+labelstring2+";Relative unc.")
    h_INPUT.Draw("hist")
    hVyyInput.Draw("hist,same")
    leg8 = TLegend(0.3,0.8,0.5,0.9)
    leg8.SetBorderSize(0)
    leg8.SetFillStyle(0)
    leg8.SetTextSize(0.04)
    leg8.AddEntry(hVyyInput,"Measured","l")
    leg8.AddEntry(h_INPUT,"Unfolded","l")
    leg8.Draw()
    c8.SaveAs("UnfoldingPlots/compareUnc_input_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")

    h_BSTAT.SetMinimum(0.0)
    h_BSTAT.Draw("hist")
    h_BSTAT.SetTitle(";top "+labelstring2+";Relative unc.")
    hVyyBstat.Draw("hist,same")
    leg8.Draw()
    c8.SaveAs("UnfoldingPlots/compareUnc_bstat_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")

    h_BSCALE.SetMinimum(0.0)
    h_BSCALE.Draw("hist")
    h_BSCALE.SetTitle(";top "+labelstring2+";Relative unc.")
    hVyyBscale.Draw("hist,same")
    leg8.Draw()
    c8.SaveAs("UnfoldingPlots/compareUnc_bscale_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")

    hVxxDiag.SetMinimum(0.0)
    hVxxDiag.SetLineColor(1)
    hVxxDiag.SetLineWidth(2)
    hVxxDiag.SetTitle(";top "+labelstring2+";Relative unc.")
    hVxxDiag.Draw("hist")
    hVyyTotal.Draw("hist,same")
    leg8.Draw()
    c8.SaveAs("UnfoldingPlots/compareUnc_total_"+options.toUnfold+"_"+options.level+"_"+channel+".pdf")
    '''
    # -------------------------------------------------------------------------------------
    # Plot covariance matrix
    # -------------------------------------------------------------------------------------
    gStyle.SetPadRightMargin(0.12)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadBottomMargin(0.15)

    diags = []
    for ibin in xrange(0,nbinsTrue):
        tmp = math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) if (math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) is not 0) else 1.0
        diags.append(tmp)
    for ibinx in xrange(1,nbinsTrue+1):
        for ibiny in xrange(1,nbinsTrue+1):
            tmp = hErrTot.GetBinContent(ibinx,ibiny) / diags[ibinx-1] / diags[ibiny-1]
            hErrTot.SetBinContent(ibinx,ibiny,tmp)

    c4 = TCanvas()
    if options.toUnfold == "pt":
        hErrTot.GetXaxis().SetRangeUser(400.,1199.)
        hErrTot.GetYaxis().SetRangeUser(400.,1199.)
    hErrTot.GetZaxis().SetRangeUser(-1.0,1.0)
    hErrTot.SetTitle(";unfolded top "+labelstring1+" "+labelstring2+";unfolded top "+labelstring1+" "+labelstring2)

    hErrTot.Draw("colz")

    drawCMS(0.16,0.945)

    c4.SaveAs("UnfoldingPlots/covariance_"+options.toUnfold+"_"+options.level+"_"+channel+"_data.pdf")

    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadBottomMargin(0.16)
    gStyle.SetPadLeftMargin(0.18)
    
    # -------------------------------------------------------------------------------------
    # Translate to cross section (not events) in bins of pt N/L/BR)
    # -------------------------------------------------------------------------------------

    BR = 0.172
    if channel is "el":
        BR = 0.173
    if channel is "comb":
        BR = 0.172 + 0.173

    thisTrue[channel].Scale(1.0/(lum*BR)) # true @ parton level
    thisMeas[channel].Scale(1.0/(lum*BR)) # measured @ reco level
    thisReco.Scale(1.0/(lum*BR)) # unfolded to parton level
    
    # -------------------------------------------------------------------------------------
    # Adjust for bin width
    # -------------------------------------------------------------------------------------
    
    for ibin in range(1, nbinsTrue+1 ) :
        
        width = thisTrue[channel].GetBinWidth(ibin)
        
        thisTrue[channel].SetBinContent(ibin, thisTrue[channel].GetBinContent(ibin) / width )
        thisTrue[channel].SetBinError(ibin, thisTrue[channel].GetBinError(ibin) / width )
        
        thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
        thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )
        
    for ibin in range(1, nbinsMeas+1) :
        
        width = thisMeas[channel].GetBinWidth(ibin)
        
        thisMeas[channel].SetBinContent(ibin,  thisMeas[channel].GetBinContent(ibin) / width )
        thisMeas[channel].SetBinError(ibin,  thisMeas[channel].GetBinError(ibin) / width )
        
    # -------------------------------
    # Make table of unfolding (w/ unc.)
    # -------------------------------
    
    print labelstring2+"  & data  & stat. [%] & exp. [%] & th. [%] & lumi [%] & total [%] & PowhegPythia8"
    for ibin in xrange(1,nbinsTrue+1):
        print string.ljust("$[{:.1f},{:.1f}]$".format(thisTrue[channel].GetXaxis().GetBinLowEdge(ibin), thisTrue[channel].GetXaxis().GetBinUpEdge(ibin)),16) + " & " + string.rjust("{:.2f}".format(thisReco.GetBinContent(ibin)*1000.0),6) + " & " + string.rjust("{:.2f}".format(h_STAT.GetBinContent(ibin)*100.),4) + " & " + string.rjust("{:.2f}".format(h_EXP.GetBinContent(ibin)*100.),4) + " & " + string.rjust("{:.2f}".format(h_TH.GetBinContent(ibin)*100.),4) + " & 2.50 & " + string.rjust("{:.2f}".format(h_TOT.GetBinContent(ibin)*100.),4) + " & " + string.rjust("{:.2f}".format(thisTrue[channel].GetBinContent(ibin)*1000.0),4)        

    # -------------------------------------------------------------------------------------
    # draw parton-level unfolding
    # -------------------------------------------------------------------------------------
    
    ## ratio of unfolded data to generator-level
    hFrac = thisTrue[channel].Clone()
    hFrac.SetName("hFrac")
    hFrac.SetTitle(";Top "+labelstring1+" "+labelstring2+";Theory/Data")
    hFrac.Divide(thisReco)
    
    # Convert h_TOT, h_STAT into error bands
    for ibin in xrange(1,nbinsTrue+1):
        h_TOT.SetBinError(ibin,h_TOT.GetBinContent(ibin))
        h_TOT.SetBinContent(ibin,1.0)
        h_STAT.SetBinError(ibin,h_STAT.GetBinContent(ibin))
        h_STAT.SetBinContent(ibin,1.0)
        
    h_TOT.SetMarkerSize(0)
    h_TOT.SetLineColor(0)
    h_TOT.SetFillColor(860-9)
    h_TOT.SetFillStyle(1001)
    
    h_STAT.SetMarkerSize(0)
    h_STAT.SetLineColor(0)
    h_STAT.SetFillColor(18)
    h_STAT.SetFillStyle(1001)
    
    # Plot
    c1 = TCanvas("c", "c", 700, 700)
    pad1 =  TPad("pad1","pad1",0,0.3,1,1)
    pad1.SetBottomMargin(0.05)
    pad1.Draw()
    pad1.cd()
    
    thisReco.SetMarkerStyle(21)
    thisMeas[channel].SetMarkerStyle(25)

    if options.toUnfold == "pt":
        thisReco.GetXaxis().SetRangeUser(400.,1199.)
        thisTrue[channel].GetXaxis().SetRangeUser(400.,1199.)
        thisMeas[channel].GetXaxis().SetRangeUser(400.,1199.)
    
    if options.toUnfold == "pt":
        xsec_title = ";;d#sigma/dp_{T} [fb/GeV]"
    else:
        xsec_title = ";;d#sigma/dy [fb]"        
    
    thisReco.SetTitle(xsec_title)
    thisReco.GetYaxis().SetTitleOffset(1.2)
    thisReco.SetMinimum(0.0)
    thisReco.SetMaximum(max(thisTrue[channel].GetMaximum(),thisReco.GetMaximum())*1.15)
    thisReco.Draw()
    thisTrue[channel].Draw('hist same')
    thisMeas[channel].Draw('same')
    thisTrue[channel].UseCurrentStyle()
    thisTrue[channel].SetLineColor(4)
    thisTrue[channel].GetYaxis().SetTitleSize(25)
    thisTrue[channel].GetXaxis().SetLabelSize(0)
    
    leg = TLegend(0.5, 0.5, 0.9, 0.75)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.SetBorderSize(0)
    
    tt = TLatex()
    tt.SetNDC()
    tt.SetTextFont(42)
    leg.AddEntry( thisReco, 'Unfolded data', 'p')
    leg.AddEntry( thisTrue[channel], 'Generated (Powheg)', 'l')
    leg.AddEntry( thisMeas[channel], 'Measured data', 'p')
    leg.AddEntry( h_STAT, 'Stat. uncertainty','f')
    leg.AddEntry( h_TOT, 'Stat. #oplus syst. uncertainties','f')
    
    leg.Draw()
    drawCMS(0.19,0.945)
    
    text1 = TLatex()
    text1.SetNDC()
    text1.SetTextFont(42)
    text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 35.9 fb^{-1}, #sqrt{s} = 13 TeV}")
    
    c1.cd()
    pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.4)
    pad2.Draw()
    pad2.cd()
    pad2.SetGridy()
    hFrac.SetMaximum(1.8)
    hFrac.SetMinimum(0.2)
    hFrac.UseCurrentStyle()
    hFrac.GetYaxis().SetTitleSize(25)
    hFrac.GetYaxis().SetTitleOffset(2.0)
    hFrac.GetXaxis().SetTitleOffset(4.0)
    hFrac.GetXaxis().SetLabelSize(25)
    hFrac.GetYaxis().SetNdivisions(4,4,0,False)
    
    hFrac.Draw("hist")
    h_TOT.Draw("same,e2")
    h_STAT.Draw("same,e2")
    hFrac.Draw("same,hist")
    if options.toUnfold == "pt":
        hFrac.GetXaxis().SetRangeUser(400., 1199.)
    
    c1.Update()
    
    c1.SaveAs("UnfoldingPlots/closure_"+options.toUnfold+"_"+options.level+"_"+channel+"_data.pdf")

    # Finally, plot response matrices
    gStyle.SetPadRightMargin(0.15)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadBottomMargin(0.15)

    for ibiny in xrange(1,nbinsTrue+1):
        tmp = response[channel].Integral(0,nbinsMeas+1,ibiny,ibiny)
        for ibinx in xrange(1,nbinsMeas+1):
            response[channel].SetBinContent(ibinx,ibiny,response[channel].GetBinContent(ibinx,ibiny)/tmp)    

    c5 = TCanvas()
    response[channel].SetTitle(";reconstructed top jet "+labelstring2+";true top "+labelstring1+" "+labelstring2)
    if options.toUnfold == "pt" :
        response[channel].GetZaxis().SetRangeUser(0.0,0.035)
    response[channel].Draw("colz")

    c5.SaveAs("UnfoldingPlots/responseMatrix_"+options.toUnfold+"_"+options.level+"_"+channel+"_nom.pdf")

    '''
    xvec.Print()
    yvec.Print()

    for systematic in ["ISRUp","ISRDown","FSRUp","FSRDown","TuneUp","TuneDown","HdampUp","HdampDown","ErdOn"]:
        for ibiny in xrange(1,nbinsTrue+1):
            tmp = Hres_sys[systematic+"_"+channel].Integral(0,nbinsMeas+1,ibiny,ibiny)
            for ibinx in xrange(1,nbinsMeas+1):
                Hres_sys[systematic+"_"+channel].SetBinContent(ibinx,ibiny,Hres_sys[systematic+"_"+channel].GetBinContent(ibinx,ibiny)/tmp)
        Hres_sys[systematic+"_"+channel].Add(response[channel],-1.0)
        Hres_sys[systematic+"_"+channel].GetZaxis().SetRangeUser(-0.015,0.015)
        Hres_sys[systematic+"_"+channel].Draw("colz")
        c5.SaveAs("UnfoldingPlots/responseMatrix_"+options.toUnfold+"_"+options.level+"_"+channel+"_"+systematic+".pdf")

        dA = TMatrixD(nbinsMeas,nbinsTrue)
        for irow in xrange(0,nbinsMeas):
            for icol in xrange(0,nbinsTrue):
                dA[irow][icol] = Hres_sys[systematic+"_"+channel].GetBinContent(irow+1,icol+1)
            
        #Manually calculate Vxx
        dAT = TMatrixD(nbinsTrue,nbinsMeas)
        dAT.Transpose(dA)
        
        dx1 = TMatrixD(nbinsMeas,1)
        dx1.Mult(dA,xvec)
        dx2 = TMatrixD(nbinsMeas,1)
        dx2.Mult(VyyInv,dx1)
        dx3 = TMatrixD(nbinsTrue,1)
        dx3.Mult(AT,dx2)

        dx4 = TMatrixD(nbinsMeas,1)
        dx4.Mult(A,xvec)
        dx5 = TMatrixD(yvec)
        dx5 -= dx4
        dx6 = TMatrixD(nbinsMeas,1)
        dx6.Mult(VyyInv,dx5)
        dx7 = TMatrixD(nbinsTrue,1)
        dx7.Mult(dAT,dx6)

        dx7 -= dx3

        dx = TMatrixD(nbinsTrue,1)
        dx.Mult(Vxx,dx7)
        dx.Print()

    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadBottomMargin(0.16)
    gStyle.SetPadLeftMargin(0.18)
    '''
