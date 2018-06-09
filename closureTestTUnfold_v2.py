#!/usr/bin/env python

import math

# -------------------------------------------------------------------------------------
# Script for doing RooUnfold on the ttbar differential cross secion
# -------------------------------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()

# -------------------------------------------------------------------------------------
# input options
# -------------------------------------------------------------------------------------

parser.add_option('--lepType', metavar='F', type='string', action='store',
                  default='muon',
                  dest='lepType',
                  help='Lepton type (ele or muon)')

parser.add_option('--toUnfold', metavar='F', type='string', action='store',
                  default='pt',
                  dest='toUnfold',
                  help='Distribution to unfold (pt or y)')

parser.add_option('--level', metavar='F', type='string', action='store',
                  default='gen',
                  dest='level',
                  help='Level to unfold (gen or part)')

parser.add_option('--type', metavar='F', type='string', action='store',
                  default='full',
                  dest='type',
                  help='')

parser.add_option('--toy', metavar='F', type='string', action='store',
                  default='',
                  dest='toy',
                  help='')

parser.add_option('--tauMode', metavar='F', type='string', action='store',
                  default='scanTau',
                  dest='tauMode',
                  help='ScanTau, LCurve, or Tau0')

parser.add_option('--regMode', metavar='F', type='string', action='store',
                  default='curvature',
                  dest='regMode',
                  help='curvature or derivative')

parser.add_option('--fullRange', metavar='F', action='store_true',
                  default=False,
                  dest='fullRange',
                  help='Pt range [350,2000]')

parser.add_option('--areaConstraint', metavar='F', action='store_true',
                  default=False,
                  dest='areaConstraint',
                  help='Add area constraint')

parser.add_option('--doSys', metavar='F', action='store_true',
                  default=False,
                  dest='doSys',
                  help='Add sys uncertainties')

# -------------------------------------------------------------------------------------
# load options & check they make sense 
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

if options.toUnfold != "pt" and options.fullRange :
    print 'Full range option only works with pt unfolding! Exiting...'
    exit(0)

# -------------------------------------------------------------------------------------
# set plot style
# -------------------------------------------------------------------------------------

import sys

from ROOT import gRandom, TH1, TH1D, TH1F, TH2F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD, TGraph, TUnfold, Double, TSpline, TSpline3, TUnfoldDensity, TUnfoldSys, TUnfold, TAttLine, TStyle, TMatrixT, TMatrixD
from array import array

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

print "TUnfold version is " + str(TUnfold.GetTUnfoldVersion())

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
        print 'Bin (' + str(h_response.GetNbinsX()+1) + ',' + str(ii) + ') ' + str(overflow) + ' added to bin (0,'+str(ii)+') ' + str(underflow)
    for jj in xrange(1,h_response.GetNbinsX()+1):
        underflow = h_response.GetBinContent(jj,0)
        overflow  = h_response.GetBinContent(jj,h_response.GetNbinsY()+1)
        h_response.SetBinContent(jj,0,underflow+overflow)
        print 'Bin (' + str(jj) + ',0) ' + str(overflow) + ' added to bin ('+str(jj)+','+str(h_response.GetNbinsY()+1)+') ' + str(underflow)

# Subtract fakes from input
def removeFakes(h_input, h_response):
    for ii in xrange(1,h_response.GetNbinsX()+1):
        if h_response.Integral(ii,ii,0,response.GetNbinsY()+1) > 0.0 :
            fakefraction = h_response.GetBinContent(ii,0) / h_response.Integral(ii,ii,0,h_response.GetNbinsY()+1)
            print 'Bin [' + str(h_response.GetXaxis().GetBinLowEdge(ii)) + ',' + str(h_response.GetXaxis().GetBinUpEdge(ii)) + '] fakefraction ' + str(fakefraction)
            h_input.SetBinContent(ii,h_input.GetBinContent(ii)*(1.0-fakefraction))
            h_input.SetBinError(ii,h_input.GetBinError(ii)*(1.0-fakefraction))
            h_response.SetBinContent(ii,0,0)
            h_response.SetBinError(ii,0,0)


# -------------------------------------------------------------------------------------
# Set unfolding options
# -------------------------------------------------------------------------------------
    
lum = 35867.0

response_name = "response_"+options.toUnfold+"_split_TH2"
hMeas_name    = options.toUnfold+"RecoTop_split"
hTrue_name    = options.toUnfold+"GenTop"

if options.level == "part":
    response_name += "_PL"
    hTrue_name.replace("Gen","Part")

if options.fullRange:
    response_name.replace("response","response2")
    hMeas_name.replace("Top","Top2")
    hTrue_name.replace("Top","Top2")

labelstring = "quark" if (options.level == "gen") else "jet"

theRegMode = TUnfold.kRegModeCurvature
if options.regMode == "derivative" :
    print 'Doing derivative regularization'
    theRegMode = TUnfold.kRegModeDerivative

theAreaConstraint = TUnfold.kEConstraintNone
if options.areaConstraint :
    theAreaConstraint = TUnfold.kEConstraintArea

suffix = ""
if options.areaConstraint:
    suffix += "_area"
elif options.fullRange:
    suffix += "_ptFull"
elif options.tauMode == "Tau0" :
    suffix = ""
else :
    suffix += ("_" + options.regMode)

if options.doSys:
    suffix += "_sys"
# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

DIR="histfiles_full2016"
PL="_PL"

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"

if options.type == "full":
    f_ttbar_m0to700_p1_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_nom_post.root")
    f_ttbar_m0to700_p2_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_nom_post.root")
    f_ttbar_m700to1000_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_nom_post.root")
    f_ttbar_m1000toInf_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_nom_post.root")
    f_ttbar_m0to700_p1_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_nom_post.root")
    f_ttbar_m0to700_p2_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_nom_post.root")
    f_ttbar_m700to1000_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_nom_post.root")
    f_ttbar_m1000toInf_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_nom_post.root")    
else:
    f_ttbar_m0to700_p1_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_nom_even_post.root")
    f_ttbar_m0to700_p2_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_nom_even_post.root")
    f_ttbar_m700to1000_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_nom_even_post.root")
    f_ttbar_m1000toInf_input    = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_nom_even_post.root")
    f_ttbar_m0to700_p1_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_nom_odd_post.root")
    f_ttbar_m0to700_p2_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_nom_odd_post.root")
    f_ttbar_m700to1000_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_nom_odd_post.root")
    f_ttbar_m1000toInf_response = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_nom_odd_post.root")
    

# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------
response_m0to700_p1 = f_ttbar_m0to700_p1_response.Get(response_name)
response_m0to700_p2 = f_ttbar_m0to700_p2_response.Get(response_name)
response_m700to1000 = f_ttbar_m700to1000_response.Get(response_name)
response_m1000toInf = f_ttbar_m1000toInf_response.Get(response_name)
response_m0to700_p1.Sumw2()
response_m0to700_p2.Sumw2()
response_m700to1000.Sumw2()
response_m1000toInf.Sumw2()
response_m0to700 = response_m0to700_p1.Clone()
response_m0to700.Add(response_m0to700_p2)
response_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
response_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
response_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
response = response_m0to700.Clone()
response.Add(response_m700to1000);
response.Add(response_m1000toInf);
    
TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# Get systematic variations, if using
# -------------------------------------------------------------------------------------
Hres_sys = {}
variants = ['Up']

if options.doSys:
    sysnames = ['JEC','JER','BTag','TopTag','lep','pu','PDF','Q2']
    thsysnames = ['ISR','FSR','Tune','Hdamp','ErdOn']
    allsysnames = sysnames+thsysnames
    longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","Pileup","PDF Uncertainty","#mu_{R}, #mu_{F} scales","ISR","FSR","Tune","ME-PS matching","Color reconnection"]
    variants = ['Up','Down']

    for sysname in sysnames:
        for var in variants:
            if options.type == "full":
                f_ttbar_m0to700_p1_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_"+sysname+var+"_post.root")
                f_ttbar_m0to700_p2_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_"+sysname+var+"_post.root")
                f_ttbar_m700to1000_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_"+sysname+var+"_post.root")
                f_ttbar_m1000toInf_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_"+sysname+var+"_post.root")
            else:
                f_ttbar_m0to700_p1_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_"+muOrEl+"_"+sysname+var+"_odd_post.root")
                f_ttbar_m0to700_p2_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth"+PL+"_p2_"+muOrEl+"_"+sysname+var+"_odd_post.root")
                f_ttbar_m700to1000_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m700to1000"+PL+"_"+muOrEl+"_"+sysname+var+"_odd_post.root")
                f_ttbar_m1000toInf_sys = TFile(DIR+"/hists_PowhegPythia8_fullTruth_m1000toInf"+PL+"_"+muOrEl+"_"+sysname+var+"_odd_post.root")
                
            response_sys_m0to700_p1 = f_ttbar_m0to700_p1_sys.Get(response_name)
            response_sys_m0to700_p2 = f_ttbar_m0to700_p2_sys.Get(response_name)
            response_sys_m700to1000 = f_ttbar_m700to1000_sys.Get(response_name)
            response_sys_m1000toInf = f_ttbar_m1000toInf_sys.Get(response_name)
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
            response_sys.Add(response_sys_m700to1000);
            response_sys.Add(response_sys_m1000toInf);

            true_sys_m0to700_p1 = f_ttbar_m0to700_p1_sys.Get(hTrue_name)
            true_sys_m0to700_p2 = f_ttbar_m0to700_p2_sys.Get(hTrue_name)
            true_sys_m700to1000 = f_ttbar_m700to1000_sys.Get(hTrue_name)
            true_sys_m1000toInf = f_ttbar_m1000toInf_sys.Get(hTrue_name)
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
            true_sys.Add(true_sys_m700to1000);
            true_sys.Add(true_sys_m1000toInf);

            antiTagWeight(true_sys,response_sys)
            convertForTUnfold(response_sys)
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinError(ibin,0,0)            
            
            if options.type == "half":
                response_sys.Scale(2.0)
                
            Hres_sys[sysname+var] = response_sys

    for thsysname in thsysnames:
        if thsysname is "ErdOn":
            f_ttbar_sys = TFile(DIR+"/hists_PowhegPythia8_"+thsysname+"_fullTruth_"+muOrEl+"_"+thsysname+"_post.root")
            response_sys = f_ttbar_sys.Get(response_name)
            response_sys.Sumw2()
            true_sys = f_ttbar_sys.Get(hTrue_name)
            antiTagWeight(true_sys,response_sys)
            convertForTUnfold(response_sys)
            for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                response_sys.SetBinContent(ibin,0,0)
                response_sys.SetBinContent(ibin,0,0)

            Hres_sys[thsysname] = response_sys

        else :
            for var in variants:
                f_ttbar_sys = TFile(DIR+"/hists_PowhegPythia8_"+thsysname+var+"_fullTruth_"+muOrEl+"_"+thsysname+var+"_post.root")
                response_sys = f_ttbar_sys.Get(response_name)
                response_sys.Sumw2()
                true_sys = f_ttbar_sys.Get(hTrue_name)
                antiTagWeight(true_sys,response_sys)
                convertForTUnfold(response_sys)
                for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
                    response_sys.SetBinContent(ibin,0,0)
                    response_sys.SetBinContent(ibin,0,0)

                Hres_sys[thsysname+var] = response_sys

# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

'Getting input/true'
hTrue_name_orig = hTrue_name
if options.toy == "up":
    hMeas_name = hMeas_name.replace("Top","TopMod")
    hTrue_name = hTrue_name.replace("Top","TopMod")
elif options.toy == "dn":
    hMeas_name = hMeas_name.replace("Top","TopModDown")
    hTrue_name = hTrue_name.replace("Top","TopModDown")

thisMeas_m0to700_p1 = f_ttbar_m0to700_p1_input.Get(hMeas_name)
thisMeas_m0to700_p2 = f_ttbar_m0to700_p2_input.Get(hMeas_name)
thisMeas_m700to1000 = f_ttbar_m700to1000_input.Get(hMeas_name)
thisMeas_m1000toInf = f_ttbar_m1000toInf_input.Get(hMeas_name)
thisMeas_m0to700_p1.Sumw2()
thisMeas_m0to700_p2.Sumw2()
thisMeas_m700to1000.Sumw2()
thisMeas_m1000toInf.Sumw2()
thisMeas_m0to700 = thisMeas_m0to700_p1.Clone()
thisMeas_m0to700.Add(thisMeas_m0to700_p2)
thisMeas_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
thisMeas_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
thisMeas_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
thisMeas = thisMeas_m0to700.Clone()
thisMeas.Add(thisMeas_m700to1000);
thisMeas.Add(thisMeas_m1000toInf);

thisTrue_m0to700_p1 = f_ttbar_m0to700_p1_input.Get(hTrue_name)
thisTrue_m0to700_p2 = f_ttbar_m0to700_p2_input.Get(hTrue_name)
thisTrue_m700to1000 = f_ttbar_m700to1000_input.Get(hTrue_name)
thisTrue_m1000toInf = f_ttbar_m1000toInf_input.Get(hTrue_name)
thisTrue_m0to700_p1.Sumw2()
thisTrue_m0to700_p2.Sumw2()
thisTrue_m700to1000.Sumw2()
thisTrue_m1000toInf.Sumw2()
thisTrue_m0to700 = thisTrue_m0to700_p1.Clone()
thisTrue_m0to700.Add(thisTrue_m0to700_p2)
thisTrue_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
thisTrue_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
thisTrue_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
thisTrue = thisTrue_m0to700.Clone()
thisTrue.Add(thisTrue_m700to1000);
thisTrue.Add(thisTrue_m1000toInf);

if options.type == "half":
    thisMeas.Scale(2.0)
    thisTrue.Scale(2.0)
    
noNegBins(thisMeas)
noNegBins(thisTrue)
    
thisMeas.SetName("recolevel") 
thisTrue.SetName("truthlevel")

nbinsMeas = thisMeas.GetNbinsX()
nbinsTrue = thisTrue.GetNbinsX()

# -------------------------------------------------------------------------------------
# Convert for TUnfold
# -------------------------------------------------------------------------------------

refTrue_m0to700_p1 = f_ttbar_m0to700_p1_response.Get(hTrue_name_orig)
refTrue_m0to700_p2 = f_ttbar_m0to700_p2_response.Get(hTrue_name_orig)
refTrue_m700to1000 = f_ttbar_m700to1000_response.Get(hTrue_name_orig)
refTrue_m1000toInf = f_ttbar_m1000toInf_response.Get(hTrue_name_orig)
refTrue_m0to700_p1.Sumw2()
refTrue_m0to700_p2.Sumw2()
refTrue_m700to1000.Sumw2()
refTrue_m1000toInf.Sumw2()
refTrue_m0to700 = refTrue_m0to700_p1.Clone()
refTrue_m0to700.Add(refTrue_m0to700_p2)
refTrue_m0to700.Scale(831.76 * 35867.0 / (77229341. + 78006311. * 1191. / 1192.))
refTrue_m700to1000.Scale(831.76 * 35867.0 * 0.0967 / 38578334.0)
refTrue_m1000toInf.Scale(831.76 * 35867.0 * 0.0256 / 24495211.0)
refTrue = refTrue_m0to700.Clone()
refTrue.Add(refTrue_m700to1000);
refTrue.Add(refTrue_m1000toInf);

antiTagWeight(refTrue,response)
convertForTUnfold(response)
removeFakes(thisMeas,response)

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
#
# RUN PSEUDO EXPERIMENTS !!!
#
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# generate toys
# -------------------------------------------------------------------------------------

hToy = thisMeas.Clone() 
hToy.SetName("toys")
hToy.Reset() 

hToy_i = [] 

ntoys = 1000

ntot = thisMeas.Integral()
nentries = int(thisMeas.GetEntries())

for itoy in xrange(0,ntoys) :
    hToy_tmp = hToy.Clone() 
    hToy_tmp.SetName("toys"+str(itoy))
     
    hToy_tmp.FillRandom(thisMeas,nentries)
    hToy_tmp.Scale(ntot/hToy_tmp.Integral())
    for ibin in xrange(1,thisMeas.GetNbinsX()+1):
        hToy_tmp.SetBinError(ibin,thisMeas.GetBinError(ibin))
    noNegBins(hToy_tmp)
    hToy_i.append(hToy_tmp)


# -------------------------------------------------------------------------------------
# UNFOLDING FOR TOYS
# -------------------------------------------------------------------------------------

hBias = thisTrue.Clone()
hBias.Reset()
hBias2 = hBias.Clone()
hErr_stat = hBias.Clone()
hErr_stat2 = hBias.Clone()
hErr_tot = hBias.Clone()
hErr_tot2 = hBias.Clone()

for itoy in xrange(0,ntoys) :
    unfold_tmp = TUnfoldDensity(response,TUnfold.kHistMapOutputVert, theRegMode, theAreaConstraint, TUnfoldDensity.kDensityModeBinWidth)
    unfold_tmp.SetInput(hToy_i[itoy])
    if options.doSys:
        for sysname in allsysnames:
            if sysname == "ErdOn" :
                unfold_tmp.AddSysError(Hres_sys[sysname],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            else :
                unfold_tmp.AddSysError(Hres_sys[sysname+"Up"],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            
    if options.tauMode == "LCurve" :
        logTauX_tmp = TSpline3()
        logTauY_tmp = TSpline3()
        lCurve_tmp = TGraph()
        iBest_tmp = unfold_tmp.ScanLcurve(30,0.,0.,lCurve_tmp,logTauX_tmp,logTauY_tmp)
    elif options.tauMode == "ScanTau" :
        scanResult_tmp = TSpline3()
        iBest_tmp = unfold_tmp.ScanTau(100,0.0001,0.1,scanResult_tmp,TUnfoldDensity.kEScanTauRhoAvgSys) #CAN CHANGE
    else :
        unfold_tmp.DoUnfold(0)

    hReco_tmp = unfold_tmp.GetOutput("tmp_output")
    hReco_tmp.Add(thisTrue,-1.0)
    hBias.Add(hReco_tmp)
    hReco_tmp.Multiply(hReco_tmp)
    hBias2.Add(hReco_tmp)

    hErrInput_tmp = unfold_tmp.GetEmatrixInput("mErrInput")
    hErrStat_tmp = unfold_tmp.GetEmatrixSysUncorr("mErrStat")

    if options.doSys :
        hErrTot_tmp = unfold_tmp.GetEmatrixTotal("mErrTotalTmp")    
        for ibin in xrange(0,nbinsTrue):
            hErr_stat.SetBinContent(ibin+1,hErr_stat.GetBinContent(ibin+1)+math.sqrt(hErrInput_tmp.GetBinContent(ibin+1,ibin+1)+hErrStat_tmp.GetBinContent(ibin+1,ibin+1)))
            hErr_stat2.SetBinContent(ibin+1,hErr_stat2.GetBinContent(ibin+1)+hErrInput_tmp.GetBinContent(ibin+1,ibin+1)+hErrStat_tmp.GetBinContent(ibin+1,ibin+1))
            hErr_tot.SetBinContent(ibin+1,hErr_tot.GetBinContent(ibin+1)+math.sqrt(hErrTot_tmp.GetBinContent(ibin+1,ibin+1)))
            hErr_tot2.SetBinContent(ibin+1,hErr_tot2.GetBinContent(ibin+1)+hErrTot_tmp.GetBinContent(ibin+1,ibin+1))
    else : #If not including systematic uncertainties, look at input stat. unc. vs matrix stat. unc. instead
        for ibin in xrange(0,nbinsTrue):
            hErr_stat.SetBinContent(ibin+1,hErr_stat.GetBinContent(ibin+1)+math.sqrt(hErrInput_tmp.GetBinContent(ibin+1,ibin+1)))
            hErr_stat2.SetBinContent(ibin+1,hErr_stat2.GetBinContent(ibin+1)+hErrInput_tmp.GetBinContent(ibin+1,ibin+1))
            hErr_tot.SetBinContent(ibin+1,hErr_tot.GetBinContent(ibin+1)+math.sqrt(hErrStat_tmp.GetBinContent(ibin+1,ibin+1)))
            hErr_tot2.SetBinContent(ibin+1,hErr_tot2.GetBinContent(ibin+1)+hErrStat_tmp.GetBinContent(ibin+1,ibin+1))
            
hBias.Scale(1.0/ntoys)
hErr_stat.Scale(1.0/ntoys)
hErr_tot.Scale(1.0/ntoys)
    
for ibin in xrange(0,nbinsTrue):
    hBias.SetBinError(ibin+1,math.sqrt(max(hBias2.GetBinContent(ibin+1)/ntoys - pow(hBias.GetBinContent(ibin+1),2),0.0000001))/thisTrue.GetBinContent(ibin+1))
    hErr_stat.SetBinError(ibin+1,math.sqrt(max(hErr_stat2.GetBinContent(ibin+1)/ntoys - pow(hErr_stat.GetBinContent(ibin+1),2),0.0000001))/thisTrue.GetBinContent(ibin+1))
    hErr_tot.SetBinError(ibin+1,math.sqrt(max(hErr_tot2.GetBinContent(ibin+1)/ntoys - pow(hErr_tot.GetBinContent(ibin+1),2),0.0000001))/thisTrue.GetBinContent(ibin+1))
    hBias.SetBinContent(ibin+1,hBias.GetBinContent(ibin+1)/thisTrue.GetBinContent(ibin+1))
    hErr_stat.SetBinContent(ibin+1,hErr_stat.GetBinContent(ibin+1)/thisTrue.GetBinContent(ibin+1))
    hErr_tot.SetBinContent(ibin+1,hErr_tot.GetBinContent(ibin+1)/thisTrue.GetBinContent(ibin+1))

ccc = TCanvas()
gPad.SetGridy()
hBias.GetYaxis().SetTitle("Relative bias / unc.")
if options.toUnfold == "pt":
    hBias.GetXaxis().SetTitle("Top "+labelstring+" p_{T} (GeV)")
    hBias.SetAxisRange(400,1199,"X")
else :
    hBias.GetXaxis().SetTitle("Top "+labelstring+" rapidity")
hBias.SetAxisRange(-0.5,0.5,"Y")
hBias.SetLineColor(1)
hErr_stat.SetLineColor(2)
hErr_tot.SetLineColor(4)
hBias.Draw("e")
hErr_stat.Draw("e,same")
hErr_tot.Draw("e,same")
hBias.Draw("e,same")
leg4 = TLegend(0.65, 0.25, 0.9, 0.45)
leg4.SetFillStyle(0)
leg4.SetTextFont(42)
leg4.SetTextSize(0.045)
leg4.SetBorderSize(0)
leg4.AddEntry(hBias, 'Bias', 'l')
if options.doSys :
    leg4.AddEntry(hErr_stat, 'Stat. Unc.', 'l')
    leg4.AddEntry(hErr_tot, 'Total Unc.', 'l')
else :
    leg4.AddEntry(hErr_stat, 'Input Stat. Unc.', 'l')
    leg4.AddEntry(hErr_tot, 'Matrix Stat. Unc.', 'l')
leg4.Draw()
ccc.SaveAs("UnfoldingPlots/bias_vs"+options.toUnfold+options.toy+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.type+suffix+".pdf")

# -------------------------------------------------------------------------------------
# Done with toys, doing actual unfolding
# -------------------------------------------------------------------------------------
unfold = {}
for var in variants:
    unfold[var] = TUnfoldDensity(response,TUnfold.kHistMapOutputVert, theRegMode, theAreaConstraint, TUnfoldDensity.kDensityModeBinWidth)
    unfold[var].SetInput(thisMeas)
    if options.doSys:
        for sysname in allsysnames:
            if sysname == "ErdOn" :
                unfold[var].AddSysError(Hres_sys[sysname],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            else :
                unfold[var].AddSysError(Hres_sys[sysname+var],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
            
    varFlag = ''
    if len(variants) > 1:
        varFlag = '_'+var+'Unc'

    # -------------------------------------------------------------------------------------
    # Do manual scan to get bias / stat. unc. / total unc.
    # -------------------------------------------------------------------------------------
    h_BiasVsTau    = TH1F("biasVsTau"   ,";log(#tau);Bias"      ,21,-4.0,-1.0)
    h_StatUncVsTau = TH1F("statUncVsTau",";log(#tau);Stat. Unc.",21,-4.0,-1.0)
    h_TotUncVsTau  = TH1F("totUncVsTau" ,";log(#tau);Total Unc.",21,-4.0,-1.0)
    for ii in xrange(0,21):
        thistau = pow(10,-4.0+3.0/21.0*(ii+0.5))
        unfold[var].DoUnfold(thistau)
        thisReco_tmp = unfold[var].GetOutput("reco")
        bias = 0.0
        statUnc = 0.0
        totalUnc = 0.0
        for ibin in xrange(1,nbinsTrue+1):
            bias += abs((thisReco_tmp.GetBinContent(ibin) - thisTrue.GetBinContent(ibin))/thisReco_tmp.GetBinContent(ibin))
            if options.doSys :
                statUnc += math.sqrt(unfold[var].GetEmatrixInput("inUnc").GetBinContent(ibin,ibin) + unfold[var].GetEmatrixSysUncorr("matUnc").GetBinContent(ibin,ibin))/thisReco_tmp.GetBinContent(ibin)
                totalUnc += math.sqrt(unfold[var].GetEmatrixTotal("totUnc").GetBinContent(ibin,ibin))/thisReco_tmp.GetBinContent(ibin)
            else :
                statUnc += math.sqrt(unfold[var].GetEmatrixInput("inUnc").GetBinContent(ibin,ibin))/thisReco_tmp.GetBinContent(ibin)
                totalUnc += math.sqrt(unfold[var].GetEmatrixSysUncorr("totUnc").GetBinContent(ibin,ibin))/thisReco_tmp.GetBinContent(ibin)                
        h_BiasVsTau.SetBinContent(ii+1,bias/float(nbinsTrue))
        h_StatUncVsTau.SetBinContent(ii+1,statUnc/float(nbinsTrue))
        h_TotUncVsTau.SetBinContent(ii+1,totalUnc/float(nbinsTrue))
    c5 = TCanvas()
    h_BiasVsTau.SetLineColor(1)
    h_StatUncVsTau.SetLineColor(2)
    h_TotUncVsTau.SetLineColor(4)
    h_BiasVsTau.SetLineWidth(2)
    h_StatUncVsTau.SetLineWidth(2)
    h_TotUncVsTau.SetLineWidth(2)
    h_BiasVsTau.GetYaxis().SetRangeUser(0.0,0.4)
    leg5 = TLegend(0.4, 0.7, 0.6, 0.9)
    leg5.SetFillStyle(0)
    leg5.SetTextFont(42)
    leg5.SetTextSize(0.045)
    leg5.SetBorderSize(0)
    leg5.AddEntry(h_BiasVsTau, 'Bias', 'l')
    if options.doSys:
        leg5.AddEntry(h_StatUncVsTau, 'Stat. Unc.', 'l')
        leg5.AddEntry(h_TotUncVsTau, 'Total Unc.', 'l')
    else :
        leg5.AddEntry(h_StatUncVsTau, 'Input Stat. Unc.', 'l')
        leg5.AddEntry(h_TotUncVsTau, 'Matrix Stat. Unc.', 'l')        
    h_BiasVsTau.Draw("hist")
    h_StatUncVsTau.Draw("hist,same")
    h_TotUncVsTau.Draw("hist,same")
    leg5.Draw()
    c5.SaveAs("UnfoldingPlots/CoarseTauScan_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+varFlag+suffix+".pdf")

    # -------------------------------------------------------------------------------------
    # Plot tau scan
    # -------------------------------------------------------------------------------------
    if options.tauMode == "LCurve" :
        c2 = TCanvas("c2", "c2", 700, 700)
        c2.cd()
        c2.UseCurrentStyle()
        logTauX = TSpline3()
        logTauY = TSpline3()
        lCurve = TGraph()
        bestLCurve = TGraph(1)
        iBest = unfold[var].ScanLcurve(30,0.,0.,lCurve,logTauX,logTauY)
        Tau = Double(0)
        x = Double(0)
        y = Double(0)
        logTauX.GetKnot(iBest,Tau,x)
        logTauY.GetKnot(iBest,Tau,y)
        bestLCurve.SetPoint(1,x,y)
        bestLCurve.SetMarkerColor(2)
        lCurve.GetXaxis().SetTitle("log L_{1}")
        lCurve.GetYaxis().SetTitle("log L_{2} / #tau^{2}")
        lCurve.Draw()
        bestLCurve.Draw("*")
        tl2 = TLatex()
        tl2.SetNDC()
        tl2.SetTextFont(42)
        legend = "log(#tau) = %.3e" % Tau
        tl2.DrawLatex(0.55,0.8,legend)
        c2.SaveAs("UnfoldingPlots/TauScan_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+varFlag+suffix+".pdf")

    elif options.tauMode == "ScanTau" :
        c2 = TCanvas("c2", "c2", 700, 700)
        c2.cd()
        c2.UseCurrentStyle()
        bestTau = TGraph(1)
        scanResult = TSpline3()
        iBest = unfold[var].ScanTau(100,0.0001,0.1,scanResult,TUnfoldDensity.kEScanTauRhoAvgSys)
        Tau = Double(0)
        rho = Double(0)
        scanResult.GetKnot(iBest,Tau,rho)
        bestTau.SetPoint(1,Tau,rho)
        bestTau.SetMarkerColor(2)
        scanResult.Draw("P")
        bestTau.Draw("*")
        tl2 = TLatex()
        tl2.SetNDC()
        tl2.SetTextFont(42)
        legend = "log(#tau) = %.3e" % Tau
        tl2.DrawLatex(0.55,0.8,legend)
        c2.SaveAs("UnfoldingPlots/TauScan_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+varFlag+suffix+".pdf")

    else:
        unfold[var].DoUnfold(0)

    print "chi**2=" + str(unfold[var].GetChi2A()) + "+" + str(unfold[var].GetChi2L()) + " / " + str(unfold[var].GetNdf())
    
# unfolded distribution (histogram)
# Use "Up" because uncertainty variants do not affect nominal
thisReco = unfold["Up"].GetOutput("reco")
thisReco.Sumw2()

# -------------------------------------------------------------------------------------
# Plot error breakdown
# -------------------------------------------------------------------------------------

#Statistical -- input and unfolding matrix (GetEmatrixSysUncorr() and GetEmatrixInput())
h_STAT = thisTrue.Clone("stat")
h_STAT.Reset()

h_TOT = thisTrue.Clone("tot")
h_TOT.Reset()

#Get actual error matrices / hists
hErrInput = unfold["Up"].GetEmatrixInput("mErrInput")
hErrStat = unfold["Up"].GetEmatrixSysUncorr("mErrStat")

if not options.doSys :
    hErrTot = unfold["Up"].GetEmatrixTotal("mErrTotalUp")

    for ibin in xrange(1,nbinsTrue+1):
        h_STAT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)+hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_TOT.SetBinContent (ibin,math.sqrt(hErrTot.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))

if options.doSys : #Doesn't seem necessary otherwise

    # Individual error sources
    h_SYS = {}
    for sysname in allsysnames:
        h_SYS[sysname] = thisTrue.Clone(sysname)
        h_SYS[sysname].Reset()

    h_INPUT = thisTrue.Clone()
    h_INPUT.Reset()
    
    h_MATRIX = thisTrue.Clone()
    h_MATRIX.Reset()
    
    h_LUMI = thisTrue.Clone()
    h_LUMI.Reset()
        
    # Construct average systematic error matrices
    hErrSys = {}
    hCovSysFSR = hErrInput.Clone()
    hCovSysFSR.Reset()
    for sysname in allsysnames:
        hErrSysUp = unfold["Up"].GetDeltaSysSource(sysname,"hErrSysUp_"+sysname)
        hErrSysDn =  unfold["Down"].GetDeltaSysSource(sysname,"hErrSysDown_"+sysname)
        hErrSys[sysname] = hErrSysUp.Clone()
        hErrSys[sysname].Reset()
        for ibin in xrange(1,nbinsTrue+1):
            hErrSys[sysname].SetBinContent(ibin,(pow(hErrSysUp.GetBinContent(ibin),2)+pow(hErrSysDn.GetBinContent(ibin),2))/2.0)
            if sysname is "FSR":
                for ibiny in xrange(1,nbinsTrue+1):
                    hCovSysFSR.SetBinContent(ibin,ibiny,(hErrSysUp.GetBinContent(ibin)*hErrSysUp.GetBinContent(ibiny)+hErrSysDn.GetBinContent(ibin)*hErrSysDn.GetBinContent(ibiny))/2.0);

    # Average total covariance matrices
    hErrTotUp = unfold["Up"].GetEmatrixTotal("mErrTotalUp")
    hErrTotDn = unfold["Down"].GetEmatrixTotal("mErrTotalDown")
    hErrTot = hErrTotUp.Clone()
    hErrTot.Add(hErrTotDn)
    hErrTot.Scale(0.5)
    
    # Correct FSR
    if "FSR" in allsysnames:
        hErrSys["FSR"].Scale(1.0/math.sqrt(2.0)); #Scale uncertainty by sqrt(2)
        hErrTot.Add(hCovSysFSR,1.0/math.sqrt(2.0)-1.0)        
    
    # Fill uncertainty histograms
    for ibin in xrange(1,nbinsTrue+1):
        h_STAT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)+hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_TOT.SetBinContent(ibin,math.sqrt(hErrTot.GetBinContent(ibin,ibin)/pow(thisReco.GetBinContent(ibin),2)+0.025*0.025))
        h_INPUT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_MATRIX.SetBinContent(ibin,math.sqrt(hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
        h_LUMI.SetBinContent(ibin,0.025)
        for sysname in allsysnames:
            h_SYS[sysname].SetBinContent(ibin,math.sqrt(hErrSys[sysname].GetBinContent(ibin))/thisReco.GetBinContent(ibin))

    # Plot error breakdown
    c6 = TCanvas("c6", "", 800, 600)
    c6.SetTopMargin(0.08)
    c6.SetRightMargin(0.05)
    c6.SetBottomMargin(0.14)
    c6.SetLeftMargin(0.16)

    if options.toUnfold == "pt" :
        h_TOT.GetXaxis().SetTitle("Top "+labelstring+" p_{T} (GeV)")
        h_TOT.SetAxisRange(400,1199,"X")
    else :
        h_TOT.GetXaxis().SetTitle("Top "+labelstring+" rapidity")
        
    h_TOT.GetYaxis().SetTitle("Uncertainty [%]")
    
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
    
    h_INPUT.SetLineColor(1)
    h_INPUT.SetLineWidth(2)
    h_INPUT.SetMarkerColor(1)
    h_INPUT.SetMarkerStyle(20)
    
    h_MATRIX.SetLineColor(419)
    h_MATRIX.SetLineWidth(2)
    h_MATRIX.SetMarkerColor(419)
    h_MATRIX.SetMarkerStyle(21)
    
    h_LUMI.SetLineColor(40)
    h_LUMI.SetLineWidth(2)
    h_LUMI.SetMarkerColor(40)
    h_LUMI.SetMarkerStyle(22)
    
    colors = [632,600,617,417,432,4,1,419,600,882,632,600,617]
    markers = [20,21,22,23,33,26,24,25,27,32,23,33,26]

    for isys in xrange(0,len(allsysnames)):
        h_SYS[allsysnames[isys]].SetLineColor(colors[isys])
        h_SYS[allsysnames[isys]].SetLineWidth(2)
        h_SYS[allsysnames[isys]].SetMarkerColor(colors[isys])
        h_SYS[allsysnames[isys]].SetMarkerStyle(markers[isys])

    leg6 = TLegend(0.2,0.39,0.45,0.88)
    leg6.AddEntry(h_TOT,"Total syst. uncertainty","f")
    leg6.AddEntry(h_INPUT,"Input stat. unc.","lp")
    leg6.AddEntry(h_MATRIX,"Response stat. unc.","lp")
    leg6.AddEntry(h_LUMI,"Int. luminosity","lp")
    for isys in xrange(0,9):
        leg6.AddEntry(h_SYS[allsysnames[isys]],longnames[isys],"lp")

    leg6.SetFillStyle(0);
    leg6.SetBorderSize(0);
    leg6.SetTextSize(0.04);
    leg6.SetTextFont(42);

    h_TOT.GetYaxis().SetRangeUser(0.0,1.0)
    h_TOT.Draw("hist")
    h_INPUT.Draw("ep,same")
    h_MATRIX.Draw("ep,same")
    h_LUMI.Draw("ep,same")
    for sysname in allsysnames:
        h_SYS[sysname].Draw("ep,same")

    leg6.Draw(); 

    c6.SaveAs("UnfoldingPlots/closure_relative_uncertainties_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+suffix+".pdf")

    # -------------------------------------------------------------------------------------
    # Troubleshoot stat. unc.
    # -------------------------------------------------------------------------------------
    # Get inverse of Vyy
    # TUnfold method has bug for version 17.1 -- doing this manually
    #Vyy = TMatrixD(nbinsMeas,nbinsMeas)
    #for ibin in xrange(0,nbinsMeas): #rows
    #    Vyy[ibin][ibin] = pow(thisMeas.GetBinError(ibin+1),2)

    #VyyInv = TMatrixD(Vyy)
    #VyyInv.Invert()

    #Get A (normalized response matrix)
    #A = TMatrixD(nbinsMeas,nbinsTrue)
    #for irow in xrange(0,nbinsMeas):
    #    for icol in xrange(0,nbinsTrue):
    #        A[irow][icol] = response.GetBinContent(irow+1,icol+1) / response.Integral(0,nbinsMeas+1,icol+1,icol+1)

    #Manually calculate Vxx
    #AT = TMatrixD(nbinsTrue,nbinsMeas)
    #AT.Transpose(A)
    
    #EInv1 = TMatrixD(nbinsTrue,nbinsMeas)
    #EInv1.Mult(AT,VyyInv)
    #EInv = TMatrixD(nbinsTrue,nbinsTrue)
    #EInv.Mult(EInv1,A)
    #E = TMatrixD(EInv)
    #E.Invert()

    # Plot relative stat. unc., and compare w/ data
    #relUnc1 = thisMeas.Clone()
    #for ibin in xrange(1,nbinsMeas+1): #rows
    #    relUnc1.SetBinContent(ibin,thisMeas.GetBinError(ibin) / thisMeas.GetBinContent(ibin))

    #c8 = TCanvas()
    #h_INPUT.SetTitle(";Top p_{T};Relative stat. unc.")
    #h_INPUT.SetMaximum(max(h_INPUT.GetMaximum(),relUnc1.GetMaximum())*1.2)
    #h_INPUT.Draw("hist")
    #relUnc1.SetLineColor(2)
    #relUnc1.Draw("hist,same")
    #leg8 = TLegend(0.2,0.6,0.5,0.8)
    #leg8.SetBorderSize(0)
    #leg8.SetTextSize(0.04)
    #leg8.AddEntry(relUnc1,"Stat. unc. on measurement","l")
    #leg8.AddEntry(h_INPUT,"Propagated stat. unc.","l")
    #leg8.Draw()
    #c8.SaveAs("UnfoldingPlots/compare_rel_unc_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+suffix+".pdf")

    # -------------------------------------------------------------------------------------
    # Plot covariance matrix
    # -------------------------------------------------------------------------------------
    
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
    hErrTot.Draw("colz")
    c4.SaveAs("UnfoldingPlots/closure_covariance_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+suffix+".pdf")

# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------

BR = 0.172
if options.lepType=="ele":
    BR = 0.173

thisTrue.Scale(1.0/(lum*BR)) # true @ parton level
thisMeas.Scale(1.0/(lum*BR)) # measured @ reco level
thisReco.Scale(1.0/(lum*BR)) # unfolded to parton level

# -------------------------------------------------------------------------------------
# Adjust for bin width
# -------------------------------------------------------------------------------------

for ibin in range(1, nbinsTrue+1 ) :

    width = thisTrue.GetBinWidth(ibin)
    
    thisTrue.SetBinContent(ibin, thisTrue.GetBinContent(ibin) / width )
    thisTrue.SetBinError(ibin, thisTrue.GetBinError(ibin) / width )
    
    thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
    thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )
    
for ibin in range(1, nbinsMeas+1) :
        
    width = thisMeas.GetBinWidth(ibin)
    
    thisMeas.SetBinContent(ibin,  thisMeas.GetBinContent(ibin) / width )
    thisMeas.SetBinError(ibin,  thisMeas.GetBinError(ibin) / width )

# -------------------------------------------------------------------------------------
# draw parton-level unfolding
# -------------------------------------------------------------------------------------

## ratio of unfolded data to generator-level
hFrac = thisReco.Clone()
hFrac.SetName("hFrac")
if options.toUnfold == "pt":
    hFrac.SetTitle(";Top "+labelstring+" p_{T} (GeV);Unfolded/True")
else :
    hFrac.SetTitle(";Top "+labelstring+" rapidity;Unfolded/True")    
hFrac.Divide(thisTrue)

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
thisMeas.SetMarkerStyle(25)

if options.toUnfold == "pt":
    thisReco.GetXaxis().SetRangeUser(400.,1199.)
    thisTrue.GetXaxis().SetRangeUser(400.,1199.)
    thisMeas.GetXaxis().SetRangeUser(400.,1199.)
    thisReco.SetTitle(";;d#sigma/dp_{T} [fb/GeV]")
else:
    thisReco.SetTitle(";;d#sigma/dy [fb]")

thisReco.GetYaxis().SetTitleOffset(1.2)
thisReco.SetMinimum(0.0)
max = thisTrue.GetMaximum()
max2 = thisReco.GetMaximum()
if max2 > max:
	max = max2
thisReco.SetAxisRange(0,max*1.15,"Y")
thisReco.Draw()
thisTrue.Draw('hist same')
thisMeas.Draw('same')
thisTrue.UseCurrentStyle()
thisTrue.SetLineColor(4)
thisTrue.GetYaxis().SetTitleSize(25)
thisTrue.GetXaxis().SetLabelSize(0)

leg = TLegend(0.5, 0.5, 0.9, 0.75)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)

tt = TLatex()
tt.SetNDC()
tt.SetTextFont(42)
leg.AddEntry( thisReco, 'Unfolded MC (Powheg)', 'p')
leg.AddEntry( thisTrue, 'Generated (Powheg)', 'l')
leg.AddEntry( thisMeas, 'Reco-level (Powheg)', 'p')
leg.AddEntry( h_STAT, 'Stat. uncertainty','f');
leg.AddEntry( h_TOT, 'Stat. #oplus syst. uncertainties','f');

tt.DrawLatex(0.55,0.45, "MC closure test")
leg.Draw()

# write histograms to file
thisReco.SetName("UnfoldedMC")

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

c1.SaveAs("UnfoldingPlots/closure_"+options.toUnfold+"_"+options.level+"_"+options.tauMode+"_"+options.lepType+"_"+options.toy+"_"+options.type+"_result"+suffix+".pdf")

