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

def mySmallText(x, y, color, text) :
  l = TLatex()
  l.SetTextSize(0.035) 
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

parser.add_option('--norm', metavar='M', action='store_true',
                  default=False,
                  dest='norm',
                  help='Do normalized unfolding')

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
        taggedsum = h_response.Integral(1,h_response.GetNbinsX()+1,ii,ii)
        untagged = h_true.GetBinContent(ii) - taggedsum
        h_response.SetBinContent(0,ii,untagged)

# Combine underflow / overflow
def convertForTUnfold(h_response):
    for ii in xrange(0,h_response.GetNbinsY()+2):
        underflow = h_response.GetBinContent(0,                       ii)
        overflow  = h_response.GetBinContent(h_response.GetNbinsX()+1,ii)
        h_response.SetBinContent(0,ii,underflow+overflow)
        h_response.SetBinContent(h_response.GetNbinsX()+1,ii,0)
    for jj in xrange(0,h_response.GetNbinsX()+2):
        underflow = h_response.GetBinContent(jj,0)
        overflow  = h_response.GetBinContent(jj,h_response.GetNbinsY()+1)
        h_response.SetBinContent(jj,0,underflow+overflow)
        h_response.SetBinContent(jj,h_response.GetNbinsY()+1,0)
    h_response.SetBinContent(0,0,0)
    h_response.SetBinError(0,0,0)

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
    extraTextSize = 0.68 * size

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

    l2 = TLatex()
    l2.SetTextSize(extraTextSize)
    l2.SetTextFont(42)
    l2.SetTextAngle(0)
    l2.SetNDC()
    l2.SetTextColor(1)
    l2.DrawLatex(0.63,y1,"(2016) 35.9 fb^{-1} (13 TeV)")


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

normname = ""
if options.norm:
    normname = "_norm"


# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

response = {}
thisMeas = {}
thisTrue = {}
thisExpect = {}
Hres_sys = {}
backgrounds = {}

trueHerwig = {}
plotmcnlo=False # need to rerun to account for negative weights... i (=louise) forgot about that...
if plotmcnlo:
    trueMCNLO = {}

channels = ["mu","el"]
sysnames = ["JEC","JER","BTag","TopTag","lep","pu","PDF","Q2"]
thsysnames = ["ISR","FSR","Tune","Hdamp","ErdOn","Herwig"]
allsysnames = sysnames+thsysnames

longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","Pileup","PDF Uncertainty","#mu_{R}, #mu_{F} scales","ISR","FSR","Tune","ME-PS matching","Color reconnection","Parton shower"]
variants = ["Up","Down"]

DIR = "histfiles_full2016_latest"

#Louise version
#name_TTbarNom = "PLnew"
#name_TTbarNom_p2 = "v2_PLnew"
#name_TTbar_m700to1000 = "m700to1000_PLnew"
#name_TTbar_m1000toInf = "m1000toInf_PLnew"
#name_TTbarVar = "PLnew"
#Susan version
name_TTbarNom = "fullTruth"
name_TTbarNom_p2 = "fullTruth_p2"
name_TTbar_m700to1000 = "fullTruth_m700to1000"
name_TTbar_m1000toInf = "fullTruth_m1000toInf"
name_TTbarVar = "fullTruth_PLnew"

for channel in channels:
    f_data = TFile(DIR+"/hists_Data_"+channel+".root")
    f_QCD  = TFile(DIR+"/hists_Data_"+channel+"_qcd.root")

    f_ttbar_m0to700_p1 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbarNom+"_"+channel+"_nom_post.root")
    f_ttbar_m0to700_p2 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbarNom_p2+"_"+channel+"_nom_post.root")
    f_ttbar_m700to1000 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbar_m700to1000+"_"+channel+"_nom_post.root")
    f_ttbar_m1000toInf = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbar_m1000toInf+"_"+channel+"_nom_post.root")

    f_ttbar_nonsemilep   = TFile(DIR+"/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_post.root")
    f_T_t                = TFile(DIR+"/hists_SingleTop_t_t_"+channel+"_nom_post.root")
    f_Tbar_t             = TFile(DIR+"/hists_SingleTop_tbar_t_"+channel+"_nom_post.root")
    f_T_tW               = TFile(DIR+"/hists_SingleTop_t_tW_"+channel+"_nom_post.root")
    f_Tbar_tW            = TFile(DIR+"/hists_SingleTop_tbar_tW_"+channel+"_nom_post.root")
    f_T_s                = TFile(DIR+"/hists_SingleTop_s_"+channel+"_nom_post.root")    
    f_WJets_HT100to200   = TFile(DIR+"/hists_WJets_HT100to200_"+channel+"_nom_post.root")
    f_WJets_HT200to400   = TFile(DIR+"/hists_WJets_HT200to400_"+channel+"_nom_post.root")
    f_WJets_HT400to600   = TFile(DIR+"/hists_WJets_HT400to600_"+channel+"_nom_post.root")
    f_WJets_HT600to800   = TFile(DIR+"/hists_WJets_HT600to800_"+channel+"_nom_post.root")
    f_WJets_HT800to1200  = TFile(DIR+"/hists_WJets_HT800to1200_"+channel+"_nom_post.root")
    f_WJets_HT1200to2500 = TFile(DIR+"/hists_WJets_HT1200to2500_"+channel+"_nom_post.root")
    f_WJets_HT2500toInf  = TFile(DIR+"/hists_WJets_HT2500toInf_"+channel+"_nom_post.root")
    f_ZJets              = TFile(DIR+"/hists_ZJets_"+channel+"_nom_post.root")
    f_WW                 = TFile(DIR+"/hists_WW_"+channel+"_nom_post.root")
    f_WZ                 = TFile(DIR+"/hists_WZ_"+channel+"_nom_post.root")
    f_ZZ                 = TFile(DIR+"/hists_ZZ_"+channel+"_nom_post.root")

    f_qcd_ttbar              = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbarNom+"_"+channel+"_nom_qcd_post.root")
    f_qcd_ttbar_nonsemilep   = TFile(DIR+"/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_qcd_post.root")
    f_qcd_T_t                = TFile(DIR+"/hists_SingleTop_t_t_"+channel+"_nom_qcd_post.root")
    f_qcd_Tbar_t             = TFile(DIR+"/hists_SingleTop_tbar_t_"+channel+"_nom_qcd_post.root")
    f_qcd_T_tW               = TFile(DIR+"/hists_SingleTop_t_tW_"+channel+"_nom_qcd_post.root")
    f_qcd_Tbar_tW            = TFile(DIR+"/hists_SingleTop_tbar_tW_"+channel+"_nom_qcd_post.root")
    f_qcd_T_s                = TFile(DIR+"/hists_SingleTop_s_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT100to200   = TFile(DIR+"/hists_WJets_HT100to200_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT200to400   = TFile(DIR+"/hists_WJets_HT200to400_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT400to600   = TFile(DIR+"/hists_WJets_HT400to600_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT600to800   = TFile(DIR+"/hists_WJets_HT600to800_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT800to1200  = TFile(DIR+"/hists_WJets_HT800to1200_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT1200to2500 = TFile(DIR+"/hists_WJets_HT1200to2500_"+channel+"_nom_qcd_post.root")
    f_qcd_WJets_HT2500toInf  = TFile(DIR+"/hists_WJets_HT2500toInf_"+channel+"_nom_qcd_post.root")
    f_qcd_ZJets              = TFile(DIR+"/hists_ZJets_"+channel+"_nom_qcd_post.root")
    f_qcd_WW                 = TFile(DIR+"/hists_WW_"+channel+"_nom_qcd_post.root")
    f_qcd_WZ                 = TFile(DIR+"/hists_WZ_"+channel+"_nom_qcd_post.root")
    f_qcd_ZZ                 = TFile(DIR+"/hists_ZZ_"+channel+"_nom_qcd_post.root")

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
            f_ttbar_sys_m0to700_p1 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbarNom+"_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m0to700_p2 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbarNom_p2+"_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m700to1000 = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbar_m700to1000+"_"+channel+"_"+sysname+var+"_post.root")
            f_ttbar_sys_m1000toInf = TFile(DIR+"/hists_PowhegPythia8_"+name_TTbar_m1000toInf+"_"+channel+"_"+sysname+var+"_post.root")    

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
        if thsysname is "ErdOn" or thsysname is "Herwig":
            f_ttbar_sys = TFile(DIR+"/hists_PowhegPythia8_"+thsysname+"_"+name_TTbarVar+"_"+channel+"_"+thsysname+"_post.root")
            response_sys = f_ttbar_sys.Get(response_name)
            response_sys.Sumw2()
            true_sys = f_ttbar_sys.Get(hTrue_name)

            antiTagWeight(true_sys,response_sys)
            convertForTUnfold(response_sys)
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinContent(ibin,0,0)

            Hres_sys[thsysname+"_"+channel] = response_sys

            ## additional theory comparison 
            if thsysname is "Herwig":
                trueHerwig[channel] = f_ttbar_sys.Get(hTrue_name)
                trueHerwig[channel].Sumw2()
                trueHerwig[channel].Scale(831.76 * 35867.0 / 59174465) #Number events in miniAOD datasets

        else :
            for var in variants:
                f_ttbar_sys = TFile(DIR+"/hists_PowhegPythia8_"+thsysname+var+"_"+name_TTbarVar+"_"+channel+"_"+thsysname+var+"_post.root")
                response_sys = f_ttbar_sys.Get(response_name)
                response_sys.Sumw2()
                true_sys = f_ttbar_sys.Get(hTrue_name)

                antiTagWeight(true_sys,response_sys)
                convertForTUnfold(response_sys)
                for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                    response_sys.SetBinContent(ibin,0,0)
                    response_sys.SetBinContent(ibin,0,0)

                Hres_sys[thsysname+var+"_"+channel] = response_sys

    ## additional theory comparison 
    if plotmcnlo:
        f_mcnlo = TFile(DIR+"/hists_TTJets_amcatnloFXFX_"+name_TTbarVar+"_"+channel+"_nom_post.root")
        trueMCNLO[channel] = f_mcnlo.Get(hTrue_name)
        trueMCNLO[channel].Sumw2()
        trueMCNLO[channel].Scale(831.76 * 35867.0 / 44350533) #Number events in miniAOD datasets


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
    bkg_TT_nonsemi = Background("t#bar{t} other",0.81,0.09,632-7)
    bkg_TT_nonsemi.addHist(hMeas_tt_nonsemi)
    backgrounds[channel].append(bkg_TT_nonsemi)

    bkg_SingleTop = Background("Single t",1.16,0.37,6)
    for hist in [hMeas_T_t,hMeas_Tbar_t, hMeas_T_tW, hMeas_Tbar_tW, hMeas_T_s] :
        bkg_SingleTop.addHist(hist)
    backgrounds[channel].append(bkg_SingleTop)

    bkg_WJetsL = Background("W+jets LF",0.77,0.21,416)
    bkg_WJetsHF = Background("W+jets HF",0.98,0.32,416-3)
    for hist in [hMeas_WJets_HT100to200,hMeas_WJets_HT200to400,hMeas_WJets_HT400to600,hMeas_WJets_HT600to800,hMeas_WJets_HT800to1200,hMeas_WJets_HT1200to2500,hMeas_WJets_HT2500toInf] :
        bkg_WJetsHF.addHist( hist )
    for hist in [hMeas_WJetsL_HT100to200,hMeas_WJetsL_HT200to400,hMeas_WJetsL_HT400to600,hMeas_WJetsL_HT600to800,hMeas_WJetsL_HT800to1200,hMeas_WJetsL_HT1200to2500,hMeas_WJetsL_HT2500toInf] :
        bkg_WJetsL.addHist( hist )
        hist.Scale(-1.0)
        bkg_WJetsHF.addHist( hist )
    backgrounds[channel].append(bkg_WJetsL)
    backgrounds[channel].append(bkg_WJetsHF)

    bkg_ZJets = Background("Z+jets",0.94,0.30,860-9)
    bkg_ZJets.addHist(hMeas_ZJets)
    backgrounds[channel].append(bkg_ZJets)

    bkg_Diboson = Background("Diboson",1.01,0.30,880-6)
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

    f_QCD_HT500to700   = TFile(DIR+"/hists_QCD_HT500to700_"+channel+"_nom_post.root")
    f_QCD_HT700to1000  = TFile(DIR+"/hists_QCD_HT700to1000_"+channel+"_nom_post.root")
    f_QCD_HT1000to1500 = TFile(DIR+"/hists_QCD_HT1000to1500_"+channel+"_nom_post.root")
    f_QCD_HT1500to2000 = TFile(DIR+"/hists_QCD_HT1500to2000_"+channel+"_nom_post.root")
    f_QCD_HT2000toInf  = TFile(DIR+"/hists_QCD_HT2000toInf_"+channel+"_nom_post.root")

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
        bkg_QCD = Background("muQCD",1.02,1.23,400)
    else:
        bkg_QCD = Background("elQCD",1.32,0.85,400)
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

trueHerwig["comb"] = trueHerwig["mu"].Clone()
trueHerwig["comb"].Add(trueHerwig["el"])

if plotmcnlo:
    trueMCNLO["comb"] = trueMCNLO["mu"].Clone()
    trueMCNLO["comb"].Add(trueMCNLO["el"])

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
    if sysname == "ErdOn" or sysname == "Herwig":
        Hres_sys[sysname+"_comb"] = Hres_sys[sysname+"_mu"].Clone()
        Hres_sys[sysname+"_comb"].Add(Hres_sys[sysname+"_el"])
    else :
        for var in variants:
            Hres_sys[sysname+var+"_comb"] = Hres_sys[sysname+var+"_mu"].Clone()
            Hres_sys[sysname+var+"_comb"].Add(Hres_sys[sysname+var+"_el"])

# -------------------------------------------------------------------------------------
# Do unfolding
# -------------------------------------------------------------------------------------

channels = ["mu","el","comb"]

# for writing output root file for using with Kostas' plotting script
LEVEL = "Parton"
if options.level == "part":
    LEVEL = "Particle"

fout = TFile("CrossSection_"+LEVEL+"_"+options.toUnfold+normname+"_TMP.root", "RECREATE")
fout.cd()


for channel in channels:
    unfold = {}
    for var in variants:
        unfold[var] = TUnfoldDensity(response[channel],TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintNone, TUnfoldDensity.kDensityModeBinWidth)
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

    # -------------------------------------------------------------------------------------
    # Plot error breakdown
    # -------------------------------------------------------------------------------------
    h_DUMMY = thisTrue[channel].Clone()
    h_DUMMY.Reset()
    
    #Statistical uncertainties
    h_INPUT = h_DUMMY.Clone()
    h_MATRIX = h_DUMMY.Clone()
    hErrInput = unfold["Up"].GetEmatrixInput("mErrInput")
    hErrStat  = unfold["Up"].GetEmatrixSysUncorr("mErrStat")
    for ibin in xrange(1,nbinsTrue+1):
      h_INPUT. SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)))
      h_MATRIX.SetBinContent(ibin,math.sqrt(hErrStat.GetBinContent(ibin,ibin)))

    #Background uncertainties
    hErrBkgStat = {}
    hErrBkgScale = {}
    for background in backgrounds[channel]:
      hErrBkgStat[background.name] = h_DUMMY.Clone()
      hErrBkgScale[background.name] = h_DUMMY.Clone()
      tmpErrBkgStat  = unfold["Up"].GetEmatrixSysBackgroundUncorr(background.name,"mErrBkgStat_"+background.name)
      tmpErrBkgScale = unfold["Up"].GetDeltaSysBackgroundScale(   background.name,"mErrBkgScale_"+background.name)
      for ibin in xrange(1,nbinsTrue+1):
        hErrBkgStat[background.name] .SetBinContent(ibin,math.sqrt(tmpErrBkgStat .GetBinContent(ibin,ibin)))
        hErrBkgScale[background.name].SetBinContent(ibin,          tmpErrBkgScale.GetBinContent(ibin,ibin))
        
    #Systematic uncertainties
    hErrSys = {}
    for variation in ["Up","Down"]:
      hErrSys[variation] = {}
      for sysname in allsysnames:
        hErrSys[variation][sysname] = unfold[variation].GetDeltaSysSource(sysname,"hErrSys"+variation+sysname)
        if sysname == "FSR" : hErrSys[variation][sysname].Scale(1.0/math.sqrt(2.0))
        
    #Combined uncertainties (+lumi)
    h_LUMI   = h_DUMMY.Clone()
    h_STAT   = h_DUMMY.Clone()
    h_BSCALE = h_DUMMY.Clone()
    h_BSTAT  = h_DUMMY.Clone()
    h_BKG    = h_DUMMY.Clone()
    h_EXP    = h_DUMMY.Clone()
    h_TH     = h_DUMMY.Clone()
    h_TOT    = h_DUMMY.Clone()
    h_SYS    = {}
    for key in hErrSys["Up"].keys():
      h_SYS[key] = h_DUMMY.Clone()

    #Construct relative uncertainties to plot, as appropriate
    hSingleSource = [h_INPUT, h_MATRIX] + hErrBkgStat.values() + hErrBkgScale.values() + hErrSys["Up"].values() + hErrSys["Down"].values()

    if options.norm :
      for hist in hSingleSource:
        hist.Add(thisReco)
        hist.Scale(1.0/hist.Integral(1,nbinsTrue+1))
        
      thisTrue[channel]  .Scale(1.0/thisTrue[channel]  .Integral(1,nbinsTrue+1))
      thisMeas[channel]  .Scale(1.0/thisMeas[channel]  .Integral(1,nbinsTrue+1))
      thisReco           .Scale(1.0/thisReco           .Integral(1,nbinsTrue+1))
      trueHerwig[channel].Scale(1.0/trueHerwig[channel].Integral(1,nbinsTrue+1))
      if plotmcnlo:
          trueMCNLO[channel] .Scale(1.0/trueMCNLO[channel].Integral(1,nbinsTrue+1))

      for hist in hSingleSource:
        hist.Add(thisReco,-1.0)
        hist.Divide(thisReco)
        for ibin in xrange(hist.GetNbinsX()):
          hist.SetBinContent(ibin+1,abs(hist.GetBinContent(ibin+1)))
          hist.SetBinError(ibin+1,0)
          
    else :
      for hist in hSingleSource:
        #hist.Divide(thisReco)
        for ibin in xrange(hist.GetNbinsX()):
          hist.SetBinContent(ibin+1,hist.GetBinContent(ibin+1)/thisReco.GetBinContent(ibin+1))
          #hist.SetBinError(ibin+1,0)
          
    for ibin in xrange(1,nbinsTrue+1):
      h_STAT.SetBinContent(ibin,math.sqrt(pow(h_INPUT.GetBinContent(ibin),2)+pow(h_MATRIX.GetBinContent(ibin),2)))

      tot_bkg_stat  = 0.0
      tot_bkg_scale = 0.0
      for background in backgrounds[channel]:
        tot_bkg_stat  += pow(hErrBkgStat[background.name] .GetBinContent(ibin),2)
        tot_bkg_scale += pow(hErrBkgScale[background.name].GetBinContent(ibin),2)
      h_BSTAT .SetBinContent(ibin,math.sqrt(tot_bkg_stat))
      h_BSCALE.SetBinContent(ibin,math.sqrt(tot_bkg_scale))
      h_BKG   .SetBinContent(ibin,math.sqrt(tot_bkg_stat+tot_bkg_scale))

      tot_exp = 0.0
      tot_th  = 0.0

      for sysname in ["PDF","Q2","ISR","FSR","Tune","Hdamp","ErdOn","Herwig"]:
        h_SYS[sysname].SetBinContent(ibin,math.sqrt((pow(hErrSys["Up"][sysname].GetBinContent(ibin),2)+pow(hErrSys["Down"][sysname].GetBinContent(ibin),2))/2.0))
        tot_th += (pow(hErrSys["Up"][sysname].GetBinContent(ibin),2)+pow(hErrSys["Down"][sysname].GetBinContent(ibin),2))/2.0
      for sysname in ["JEC","JER","BTag","TopTag","lep","pu"]:
        h_SYS[sysname].SetBinContent(ibin,math.sqrt((pow(hErrSys["Up"][sysname].GetBinContent(ibin),2)+pow(hErrSys["Down"][sysname].GetBinContent(ibin),2))/2.0))
        tot_exp += (pow(hErrSys["Up"][sysname].GetBinContent(ibin),2)+pow(hErrSys["Down"][sysname].GetBinContent(ibin),2))/2.0
      tot_exp += tot_bkg_stat
      tot_exp += tot_bkg_scale
      h_EXP.SetBinContent(ibin,math.sqrt(tot_exp))
      h_TH.SetBinContent(ibin,math.sqrt(tot_th))

      if options.norm:
        h_TOT.SetBinContent(ibin,math.sqrt(tot_exp+tot_th+pow(h_STAT.GetBinContent(ibin),2)))
      else:
        h_TOT.SetBinContent(ibin,math.sqrt(tot_exp+tot_th+pow(h_STAT.GetBinContent(ibin),2)+0.025*0.025))
        h_LUMI.SetBinContent(ibin,0.025)                        

    # WHAT TO DO ABOUT COVARIANCE MATRICES FOR NORM?
    #
    # Probably best to construct new
    # This is straightforward for systematic uncertainties, but how to deal with statistical ones?
    if not options.norm:
      hErrTotUp = unfold["Up"].GetEmatrixTotal("mErrTotalUp")
      hErrTotDn = unfold["Down"].GetEmatrixTotal("mErrTotalDown")
      hErrTot = hErrTotUp.Clone()
      hErrTot.Add(hErrTotDn)
      hErrTot.Scale(0.5)

      hCovSysFSR = hErrInput.Clone()
      hCovSysFSR.Reset()
      for ibinx in xrange(1,nbinsTrue+1):
        for ibiny in xrange(1,nbinsTrue+1):
          hCovSysFSR.SetBinContent(ibinx,ibiny,(hErrSys["Up"]["FSR"].GetBinContent(ibinx)*hErrSys["Up"]["FSR"].GetBinContent(ibiny)+hErrSys["Down"]["FSR"].GetBinContent(ibinx)*hErrSys["Down"]["FSR"].GetBinContent(ibiny))/2.0)
      
      # Correct FSR
      if "FSR" in allsysnames:
        hErrTot.Add(hCovSysFSR,-0.5) #Contribution from FSR to total covariance decreases by 0.5       
            
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
    if not options.norm:
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
    if not options.norm:
        h_LUMI.Draw("ep,same")
    for sysname in allsysnames:
        h_SYS[sysname].Draw("ep,same")
    h_BKG.Draw("ep,same")
        
    leg6.Draw() 

    drawCMS(0.17,0.935)
    
    c6.SaveAs("UnfoldingPlots/unfold_relative_uncertainties_"+options.toUnfold+"_"+options.level+"_"+channel+normname+".pdf")


    ##################################################################################################
    # prettified paper version of uncertainty plot 
    ##################################################################################################
    
    cerr = TCanvas("cerr", "cerr", 800, 600)
    cerr.SetTopMargin(1.0)
    cerr.SetRightMargin(0.05)
    cerr.SetBottomMargin(0.12)
    cerr.SetLeftMargin(0.16)
    cerr.cd()

    # ----------------------------------------------------------
    # all uncertainties     
    # ----------------------------------------------------------

    sysleg = TLegend(0.2,0.39,0.45,0.88)
    sysleg.AddEntry(h_TOT,"Total syst. uncertainty","f")
    sysleg.AddEntry(h_STAT,"Input stat. unc.","lp")
    if not options.norm:
        sysleg.AddEntry(h_LUMI,"Int. luminosity","lp")
    for isys in xrange(0,len(allsysnames)):
        sysleg.AddEntry(h_SYS[allsysnames[isys]],longnames[isys],"lp")
    sysleg.AddEntry(h_BKG,"Backgrounds","lp")
        
    sysleg.SetFillStyle(0)
    sysleg.SetBorderSize(0)
    sysleg.SetTextSize(0.03)
    sysleg.SetTextFont(42)

    h_TOT.GetXaxis().SetTitleOffset(1.0)

    h_TOT.GetYaxis().SetRangeUser(0.0,1.0)
    h_TOT.Draw("hist")
    h_STAT.Draw("ep,same")
    if not options.norm:
        h_LUMI.Draw("ep,same")
    for sysname in allsysnames:
        h_SYS[sysname].Draw("ep,same")
    h_BKG.Draw("ep,same")
        
    sysleg.Draw() 

    drawCMS(0.17,0.92,0.065,"")
    
    cerr.SaveAs("UnfoldingPlots/xsec_uncertainties_"+options.toUnfold+"_"+options.level+"_"+channel+normname+".pdf")


    # ----------------------------------------------------------
    # merged categories
    # ----------------------------------------------------------

    # h_TOT    => total uncertainty shaded band     
    # h_STAT   => statistical uncertainty
    # h_LUMI   => lum 
    # h_BKG    => backgrounds 

    #sysnames = ["JEC","JER","BTag","TopTag","lep","pu","PDF","Q2"]
    #thsysnames = ["ISR","FSR","Tune","Hdamp","ErdOn","Herwig"]
    #longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","Pileup","PDF Uncertainty","#mu_{R}, #mu_{F} scales","ISR","FSR","Tune","ME-PS matching","Color reconnection","Parton shower"]

    # h_sys_jet     : JEC, JER, BTag
    # h_sys_toptag  : TopTag 
    # h_sys_other   : bkg, lep, pu, h_LUMI
    # h_sys_hardscatter  : PDF, Q2
    # h_sys_partonshower : ISR, FSR, Tune, Hdamp, EdrOn, Herwig

    h_sys_stat = h_STAT.Clone("sys_stat")
    h_sys_stat.Scale(100.0)
    h_sys_toptag = h_SYS["TopTag"].Clone("sys_toptag")
    h_sys_toptag.Scale(100.0) 

    h_sys_jet = h_SYS["JEC"].Clone("sys_jet")
    h_sys_jet.Reset()
    h_sys_other = h_SYS["JEC"].Clone("sys_other")
    h_sys_other.Reset()
    h_sys_hardscatter = h_SYS["JEC"].Clone("sys_hardscatter")
    h_sys_hardscatter.Reset()
    h_sys_partonshower = h_SYS["JEC"].Clone("sys_partonshower")
    h_sys_partonshower.Reset()
    
    for ibin in xrange(1,nbinsTrue+1):
        
        bin_jet = 0;
        bin_other = 0;
        bin_hardscatter = 0;
        bin_partonshower = 0;
        
        for sysname in allsysnames:
        
            if sysname == "JEC" or sysname == "JER" or sysname == "BTag" :
                bin_jet += pow(h_SYS[sysname].GetBinContent(ibin),2)
            if sysname == "lep" or sysname == "pu":
                bin_other += pow(h_SYS[sysname].GetBinContent(ibin),2)
            if sysname == "PDF" or sysname == "Q2":
                bin_hardscatter += pow(h_SYS[sysname].GetBinContent(ibin),2)
            if sysname == "ISR" or sysname == "FSR" or sysname == "Tune" or sysname == "Hdamp" or sysname == "ErdOn" or sysname == "Herwig":
                bin_partonshower += pow(h_SYS[sysname].GetBinContent(ibin),2)

        bin_other += pow(h_BKG.GetBinContent(ibin),2) 
        if not options.norm:
            bin_other += pow(h_LUMI.GetBinContent(ibin),2) 
        
        h_sys_jet.SetBinContent(ibin,math.sqrt(bin_jet)*100.0)
        h_sys_other.SetBinContent(ibin,math.sqrt(bin_other)*100.0)
        h_sys_hardscatter.SetBinContent(ibin,math.sqrt(bin_hardscatter)*100.0)
        h_sys_partonshower.SetBinContent(ibin,math.sqrt(bin_partonshower)*100.0)


    msysleg = TLegend(0.2,0.52,0.45,0.8)
    msysleg.AddEntry(h_sys_stat,"Stat. uncertainty","f")
    msysleg.AddEntry(h_sys_jet,"JES+JER+b tagging","l")
    msysleg.AddEntry(h_sys_toptag,"t tagging","l")
    msysleg.AddEntry(h_sys_other,"Other experimental","l")
    msysleg.AddEntry(h_sys_partonshower,"Parton shower","l")
    msysleg.AddEntry(h_sys_hardscatter,"Hard scattering","l")
        
    msysleg.SetFillStyle(0)
    msysleg.SetBorderSize(0)
    msysleg.SetTextSize(0.03)
    msysleg.SetTextFont(42)
    
    h_sys_stat.SetFillColor(920)
    h_sys_stat.SetFillStyle(1001)
    h_sys_stat.SetLineWidth(0)

    h_sys_stat.GetXaxis().SetTitle("Top "+labelstring1+" "+labelstring2)
    h_sys_stat.GetYaxis().SetTitle("Relative uncertainty [%]")

    h_sys_jet.SetLineColor(2)
    h_sys_toptag.SetLineColor(4)
    h_sys_other.SetLineColor(6)
    h_sys_partonshower.SetLineColor(416)
    h_sys_hardscatter.SetLineColor(800)

    h_sys_jet.SetLineWidth(3)
    h_sys_toptag.SetLineWidth(3)
    h_sys_other.SetLineWidth(3)
    h_sys_partonshower.SetLineWidth(3)
    h_sys_hardscatter.SetLineWidth(3)

    
    h_sys_stat.GetXaxis().SetTitleOffset(1.0)
    h_sys_stat.GetYaxis().SetTitleOffset(1.1)
    h_sys_stat.GetYaxis().SetRangeUser(0.0,100.0)
    h_sys_stat.Draw("hist")

    h_sys_jet.Draw("hist same")
    h_sys_toptag.Draw("hist same")
    h_sys_other.Draw("hist same")
    h_sys_partonshower.Draw("hist same")
    h_sys_hardscatter.Draw("hist same")
        
    msysleg.Draw() 

    drawCMS(0.17,0.92,0.065,"")
    if options.level=="part":
        mySmallText(0.21,0.84,1,"Particle level (l+jets channel)")
    else :
        mySmallText(0.21,0.84,1,"Parton level (l+jets channel)")
        
    if options.norm:
        mySmallText(0.63,0.84,1,"Normalized cross section")
    else:
        mySmallText(0.63,0.84,1,"Absolute cross section")

    
    cerr.SaveAs("UnfoldingPlots/xsec_mergedUncertainties_"+options.toUnfold+"_"+options.level+"_"+channel+normname+".pdf")

    
    ##################################################################################################

    
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
    if not options.norm:
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

    if not options.norm:
        thisTrue[channel].Scale(1.0/(lum*BR)) # true @ parton level
        thisMeas[channel].Scale(1.0/(lum*BR)) # measured @ reco level
        thisReco.Scale(1.0/(lum*BR)) # unfolded to parton level

        trueHerwig[channel].Scale(1.0/(lum*BR)) # true @ parton level
        if plotmcnlo:
            trueMCNLO[channel].Scale(1.0/(lum*BR)) # true @ parton level

        
    # -------------------------------------------------------------------------------------
    # Adjust for bin width
    # -------------------------------------------------------------------------------------
    
    for ibin in range(1, nbinsTrue+1 ) :
        
        width = thisTrue[channel].GetBinWidth(ibin)
        
        thisTrue[channel].SetBinContent(ibin, thisTrue[channel].GetBinContent(ibin) / width )
        thisTrue[channel].SetBinError(ibin, thisTrue[channel].GetBinError(ibin) / width )
        
        thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
        thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )

        trueHerwig[channel].SetBinContent(ibin, trueHerwig[channel].GetBinContent(ibin) / width )
        trueHerwig[channel].SetBinError(ibin, trueHerwig[channel].GetBinError(ibin) / width )

        if plotmcnlo:
            trueMCNLO[channel].SetBinContent(ibin, trueMCNLO[channel].GetBinContent(ibin) / width )
            trueMCNLO[channel].SetBinError(ibin, trueMCNLO[channel].GetBinError(ibin) / width )

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
    dataNoUnc = thisReco.Clone("dataNoUnc")
    for ibin in xrange(1,nbinsTrue+1):
        dataNoUnc.SetBinError(ibin,0)
    hFrac.Divide(dataNoUnc)

    ## same but for Herwig
    hFracHerwig = trueHerwig[channel].Clone()
    hFracHerwig.SetName("hFracHerwig")
    hFracHerwig.SetTitle(";Top "+labelstring1+" "+labelstring2+";Theory/Data")
    hFracHerwig.Divide(dataNoUnc)

    ## and MC@NLO
    if plotmcnlo:
        hFracMCNLO = trueMCNLO[channel].Clone()
        hFracMCNLO.SetName("hFracMCNLO")
        hFracMCNLO.SetTitle(";Top "+labelstring1+" "+labelstring2+";Theory/Data")
        hFracMCNLO.Divide(dataNoUnc)

    if channel == "comb":
        hRelUncMinus = h_TOT.Clone("RelUncMinus")
        hRelUncPlus = h_TOT.Clone("RelUncPlus")
    
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

        trueHerwig[channel].GetXaxis().SetRangeUser(400.,1199.)
        if plotmcnlo:
            trueMCNLO[channel].GetXaxis().SetRangeUser(400.,1199.)
        
    if options.norm: 
        if options.toUnfold == "pt":
            xsec_title = ";;1/#sigma d#sigma/dp_{T} [1/GeV]"
        else:
            xsec_title = ";;1/#sigma d#sigma/dy []"        
    else:
        if options.toUnfold == "pt":
            xsec_title = ";;d#sigma/dp_{T} [pb/GeV]"
        else:
            xsec_title = ";;d#sigma/dy [pb]"        
    
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
    
    c1.SaveAs("UnfoldingPlots/closure_"+options.toUnfold+"_"+options.level+"_"+channel+normname+"_data.pdf")


    ##################################################################################################
    # pretty unfolded plot for paper
    ##################################################################################################
    
    # histos 
    # hFrac              => ratio 
    # h_TOT, h_STAT      => error bands for ratio 
    # thisReco           => measurement
    # thisTrue[channel]  => powheg+pythia

    # canvas
    paper_c = TCanvas("paper_c", "paper_c", 700, 700)
    paper_pad1 =  TPad("pad1","pad1",0,0.28,1,1)
    paper_pad1.SetTopMargin(0.1)
    paper_pad1.SetBottomMargin(0.02)
    paper_pad1.Draw()
    paper_pad1.cd()

    # title etc
    thisReco.GetXaxis().SetLabelSize(0)
    thisReco.GetYaxis().SetTitleSize(32)
    thisReco.GetYaxis().SetTitleOffset(1.4)
    thisReco.SetMarkerColor(1)
    thisReco.SetMarkerStyle(8)
    thisReco.SetLineColor(1)
    
    thisTrue[channel].SetLineColor(2) 
    thisTrue[channel].SetLineWidth(2) 

    # draw
    thisReco.Draw()
    thisTrue[channel].Draw('hist, same')
    thisReco.Draw("axis,same")
    
    # legend
    paper_leg = TLegend(0.5, 0.6, 0.9, 0.85)
    paper_leg.SetFillStyle(0)
    paper_leg.SetTextFont(42)
    paper_leg.SetTextSize(0.045)
    paper_leg.SetBorderSize(0)
    
    paper_tt = TLatex()
    paper_tt.SetNDC()
    paper_tt.SetTextFont(42)
    paper_leg.AddEntry( thisReco, 'Data', 'lep')
    paper_leg.AddEntry( thisTrue[channel], 'Powheg+Pythia8', 'l')
    paper_leg.AddEntry( h_STAT, 'Stat. uncertainty','f')
    paper_leg.AddEntry( h_TOT, 'Stat. #oplus syst. uncertainties','f')
    paper_leg.Draw()

    drawCMS(0.18,0.92,0.08,"")

    # ratio part of plot
    paper_c.cd()
    paper_pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
    paper_pad2.SetTopMargin(0.02)
    paper_pad2.SetBottomMargin(0.35)
    paper_pad2.Draw()
    paper_pad2.cd()
    paper_pad2.SetGridy()

    hFrac.SetLineColor(2)
    hFrac.SetLineWidth(2)
    hFrac.GetYaxis().SetTitleSize(26)
    hFrac.GetYaxis().SetTitleOffset(1.4)
    hFrac.GetXaxis().SetTitleOffset(3.3)
    hFrac.SetMaximum(1.9)
    hFrac.SetMinimum(0.1)
    hFrac.GetYaxis().SetNdivisions(4,4,0,True)

    hFrac.Draw("hist")
    h_TOT.Draw("same,e2")
    h_STAT.Draw("same,e2")
    hFrac.Draw("same,hist")
    hFrac.Draw("axis,same")
    
    paper_c.SaveAs("UnfoldingPlots/xsec_"+options.toUnfold+"_"+options.level+"_"+channel+normname+".pdf")

    if channel == "comb":

        for ibin in xrange(1,nbinsTrue+1):
            hFrac.SetBinContent(ibin, hFrac.GetBinContent(ibin)-1)
            h_STAT.SetBinContent(ibin, h_STAT.GetBinContent(ibin)-1)
            hFracHerwig.SetBinContent(ibin, hFracHerwig.GetBinContent(ibin)-1)
            if plotmcnlo:
                hFracMCNLO.SetBinContent(ibin, hFracMCNLO.GetBinContent(ibin)-1)
                    
        if options.norm: 
            thisReco.SetName("NormCrossSection_"+LEVEL+"_Nominal")
            h_STAT.SetName("NormCrossSectionRatio")
            hRelUncMinus.SetName("RelUncNormMinus")
            hRelUncPlus.SetName("RelUncNormPlus")
            thisTrue[channel].SetName("NormCrossSection_PowhegPythia8")
            hFrac.SetName("NormRatioOverData_PowhegPythia8")
            trueHerwig[channel].SetName("NormCrossSection_PowhegHerwigpp")
            hFracHerwig.SetName("NormRatioOverData_PowhegHerwigpp")
            if plotmcnlo:
                trueMCNLO[channel].SetName("NormCrossSection_amcatnloPythia8")
                hFracMCNLO.SetName("NormRatioOverData_amcatnloPythia8")
        else:
            thisReco.SetName("CrossSection_"+LEVEL+"_Nominal")
            h_STAT.SetName("CrossSectionRatio")
            thisTrue[channel].SetName("CrossSection_PowhegPythia8")
            hFrac.SetName("RatioOverData_PowhegPythia8")
            trueHerwig[channel].SetName("CrossSection_PowhegHerwigpp")
            hFracHerwig.SetName("RatioOverData_PowhegHerwigpp")
            if plotmcnlo:
                trueMCNLO[channel].SetName("CrossSection_amcatnloPythia8")
                hFracMCNLO.SetName("RatioOverData_amcatnloPythia8")
        thisReco.Write()
        h_STAT.Write()
        hRelUncMinus.Write()
        hRelUncPlus.Write()
        thisTrue[channel].Write()
        hFrac.Write()
        trueHerwig[channel].Write()
        hFracHerwig.Write()
        if plotmcnlo:
            trueMCNLO[channel].Write()
            hFracMCNLO.Write()

        fout.Close()
    
    
    ##################################################################################################

    
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
