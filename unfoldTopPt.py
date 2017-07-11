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

# -------------------------------------------------------------------------------------
# load options & set plot style
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

import sys

from ROOT import gRandom, TH1, TH1D, TH1F, TH2F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD, TGraph, TUnfold, Double, TSpline, TSpline3, TUnfoldDensity, TUnfoldSys, TUnfold, TAttLine, TStyle
from array import array
import string

gROOT.Macro("rootlogon.C")
gROOT.SetBatch(True)

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0);

gStyle.SetTitleFont(43)
gStyle.SetTitleFont(43, "XYZ")
gStyle.SetTitleSize(30, "XYZ")
gStyle.SetLabelFont(43, "XYZ")
gStyle.SetLabelSize(24, "XYZ")

gStyle.SetPadTopMargin(0.07);
gStyle.SetPadRightMargin(0.05);
gStyle.SetPadBottomMargin(0.16);
gStyle.SetPadLeftMargin(0.18);
  
gSystem.Load("RooUnfold/libRooUnfold.so")

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold

print "TUnfold version is " + str(TUnfold.GetTUnfoldVersion())

# -------------------------------------------------------------------------------------
# Define helper functions
# -------------------------------------------------------------------------------------

#def makeResponse( response ) :
#    Hres = response.HresponseNoOverflow()
#    vtru = response.Vtruth()
#    etru = response.Etruth()
#        
#    for j in xrange(1,response.GetNbinsTruth()+1) :
#        ntru = 0.0
#        errtru = 0.0
#        for i in xrange(1,response.GetNbinsMeasured()+1) :
#            ntru += Hres.GetBinContent(i,j)
#            errtru += pow(Hres.GetBinError(i,j),2)
#        Hres.SetBinContent(0, j, vtru[j-1]-ntru)
#        Hres.SetBinError(0, j, math.sqrt(etru[j-1]-errtru)) #Normally would add, but we want underflow + normal bin errors to sum to total error
#        
#    return Hres

#def removeFakes( hMeas, response ):
#    if response.FakeEntries() :
#        fakes = response.Vfakes()
#        fac = response.Vmeasured().Sum() # Measured, from response matrix
#        if fac != 0.0 :
#            measVec = RooUnfoldResponse.H2V(hMeas,response.GetNbinsMeasured(),0) #Actual measured input
#            fac = measVec.Sum() / fac 
#            for i in xrange(1,response.GetNbinsMeasured()+1) :
#                hMeas.SetBinContent(i,hMeas.GetBinContent(i) - fac * fakes[i-1])

def noNegBins( hist ):
    for i in xrange(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i) < 0.0 :
            hist.SetBinContent(i,0.0)

# -------------------------------------------------------------------------------------
# Define normalization constants
# -------------------------------------------------------------------------------------

lum = 35867.0

PowhegPythia8_norm         = 831.76         * lum / 77229341.
SingleTop_t_t_norm         = 136.02         * lum / 5993676.
SingleTop_tbar_t_norm      = 80.95          * lum / 3928063.
SingleTop_t_tW_norm        = 35.9           * lum / 992024.
SingleTop_tbar_tW_norm     = 35.9           * lum / 998276.
SingleTop_s_norm           = 10.32          * lum / 1000000.
WJets_HT100to200_norm      = 1345.0 * 1.21  * lum / 39617787.
WJets_HT200to400_norm      = 359.7 * 1.21   * lum / 19914590.
WJets_HT400to600_norm      = 48.91 * 1.21   * lum / 5796237.
WJets_HT600to800_norm      = 12.05 * 1.21   * lum / 14822888.
WJets_HT800to1200_norm     = 5.501 * 1.21   * lum / 6200954.
WJets_HT1200to2500_norm    = 1.329 * 1.21   * lum / 6324934.
WJets_HT2500toInf_norm     = 0.03216 * 1.21 * lum / 2384260.
QCD_HT500to700_norm        = 32100.         * lum / 18929951.
QCD_HT700to1000_norm       = 6831.          * lum / 15629253.  
QCD_HT1000to1500_norm      = 1207.          * lum / 4767100. 
QCD_HT1500to2000_norm      = 119.9          * lum / 3970819.
QCD_HT2000toInf_norm       = 25.24          * lum / 1991645.

response_name = "response_pt_split_TH2"
hMeas_name = "ptRecoTop_split"
hTrue_name = "ptGenTop"

# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"

f_data = TFile("histfiles_full2016/hists_Data_"+muOrEl+".root")
f_QCD  = TFile("histfiles_full2016/hists_Data_"+muOrEl+"_qcd.root")

f_ttbar        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_post.root")
f_ttbar_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_post.root")

f_ttbar_nonsemilep   = TFile("histfiles_full2016/hists_PowhegPythia8_nonsemilep_"+muOrEl+"_nom_post.root")
f_T_t                = TFile("histfiles_full2016/hists_SingleTop_t_t_"+muOrEl+"_nom_post.root")
f_Tbar_t             = TFile("histfiles_full2016/hists_SingleTop_tbar_t_"+muOrEl+"_nom_post.root")
f_T_tW               = TFile("histfiles_full2016/hists_SingleTop_t_tW_"+muOrEl+"_nom_post.root")
f_Tbar_tW            = TFile("histfiles_full2016/hists_SingleTop_tbar_tW_"+muOrEl+"_nom_post.root")
f_T_s                = TFile("histfiles_full2016/hists_SingleTop_s_"+muOrEl+"_nom_post.root")    
f_WJets_HT100to200   = TFile("histfiles_full2016/hists_WJets_HT100to200_"+muOrEl+"_nom_post.root")
f_WJets_HT200to400   = TFile("histfiles_full2016/hists_WJets_HT200to400_"+muOrEl+"_nom_post.root")
f_WJets_HT400to600   = TFile("histfiles_full2016/hists_WJets_HT400to600_"+muOrEl+"_nom_post.root")
f_WJets_HT600to800   = TFile("histfiles_full2016/hists_WJets_HT600to800_"+muOrEl+"_nom_post.root")
f_WJets_HT800to1200  = TFile("histfiles_full2016/hists_WJets_HT800to1200_"+muOrEl+"_nom_post.root")
f_WJets_HT1200to2500 = TFile("histfiles_full2016/hists_WJets_HT1200to2500_"+muOrEl+"_nom_post.root")
f_WJets_HT2500toInf  = TFile("histfiles_full2016/hists_WJets_HT2500toInf_"+muOrEl+"_nom_post.root")
    
f_qcd_ttbar              = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_qcd_post.root")
f_qcd_ttbar_nonsemilep   = TFile("histfiles_full2016/hists_PowhegPythia8_nonsemilep_"+muOrEl+"_nom_qcd_post.root")
f_qcd_T_t                = TFile("histfiles_full2016/hists_SingleTop_t_t_"+muOrEl+"_nom_qcd_post.root")
f_qcd_Tbar_t             = TFile("histfiles_full2016/hists_SingleTop_tbar_t_"+muOrEl+"_nom_qcd_post.root")
f_qcd_T_tW               = TFile("histfiles_full2016/hists_SingleTop_t_tW_"+muOrEl+"_nom_qcd_post.root")
f_qcd_Tbar_tW            = TFile("histfiles_full2016/hists_SingleTop_tbar_tW_"+muOrEl+"_nom_qcd_post.root")
f_qcd_T_s                = TFile("histfiles_full2016/hists_SingleTop_s_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT100to200   = TFile("histfiles_full2016/hists_WJets_HT100to200_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT200to400   = TFile("histfiles_full2016/hists_WJets_HT200to400_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT400to600   = TFile("histfiles_full2016/hists_WJets_HT400to600_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT600to800   = TFile("histfiles_full2016/hists_WJets_HT600to800_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT800to1200  = TFile("histfiles_full2016/hists_WJets_HT800to1200_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT1200to2500 = TFile("histfiles_full2016/hists_WJets_HT1200to2500_"+muOrEl+"_nom_qcd_post.root")
f_qcd_WJets_HT2500toInf  = TFile("histfiles_full2016/hists_WJets_HT2500toInf_"+muOrEl+"_nom_qcd_post.root")

# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

response = f_ttbar.Get(response_name)
response.Sumw2()
response2 = f_ttbar_2.Get(response_name)
response2.Sumw2()
Hres.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
Hres_2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
Hres.Add(Hres_2)

TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# Get systematic variations
# -------------------------------------------------------------------------------------

Hres_sys = {}
sysnames = ['JECUp','JECDown','JERUp','JERDown','BTagUp','BTagDown','TopTagUp','TopTagDown','lepUp','lepDown','PDFUp','PDFDown','Q2Up','Q2Down','ASUp','ASDown']

for sysname in sysnames:
    f_ttbar_sys        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+sysname+"_post.root")
    f_ttbar_sys_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_"+sysname+"_post.root")

    response_sys = f_ttbar_sys.Get(response_name)
    response_sys.Sumw2()
    response2_sys = f_ttbar_sys_2.Get(response_name)
    response_sys.Sumw2()
    Hres_tmp.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
    Hres_tmp_2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
    Hres_tmp.Add(Hres_tmp_2)
        
    Hres_sys[sysname] = Hres_tmp

# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

thisMeas = f_data.Get(hMeas_name).Clone()
thisMeas.Sumw2()

thisTrue = f_ttbar.Get(hTrue_name).Clone()
thisTrue.Sumw2()
thisTrue2 = f_ttbar_2.Get(hTrue_name).Clone()
thisTrue2.Sumw2()
thisTrue.Add(thisTrue2)
thisTrue.Scale(lum * 831.76 / (77229341. + 78006311. * 1191. / 1192.));

noNegBins(thisMeas)
noNegBins(thisTrue)
    
thisMeas.SetName("recolevel") 
thisTrue.SetName("truthlevel")

nbinsTrue = thisTrue.GetNbinsX()
nbinsMeas = thisMeas.GetNbinsX()

hMeas_tt_nonsemi         = f_ttbar_nonsemilep.Get(hMeas_name)
hMeas_T_t                = f_T_t.Get(hMeas_name)
hMeas_Tbar_t             = f_Tbar_t.Get(hMeas_name)
hMeas_T_tW               = f_T_tW.Get(hMeas_name)
hMeas_Tbar_tW            = f_Tbar_tW.Get(hMeas_name)
hMeas_T_s                = f_T_s.Get(hMeas_name)
hMeas_WJets_HT100to200   = f_WJets_HT100to200.Get(hMeas_name)
hMeas_WJets_HT200to400   = f_WJets_HT200to400.Get(hMeas_name)
hMeas_WJets_HT400to600   = f_WJets_HT400to600.Get(hMeas_name)
hMeas_WJets_HT600to800   = f_WJets_HT600to800.Get(hMeas_name)
hMeas_WJets_HT800to1200  = f_WJets_HT800to1200.Get(hMeas_name)
hMeas_WJets_HT1200to2500 = f_WJets_HT1200to2500.Get(hMeas_name)
hMeas_WJets_HT2500toInf  = f_WJets_HT2500toInf.Get(hMeas_name)

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

hMeas_SingleTop = hMeas_T_t.Clone()
for hist in [hMeas_Tbar_t, hMeas_T_tW, hMeas_Tbar_tW, hMeas_T_s] :
    hMeas_SingleTop.Add( hist )
    
hMeas_WJets = hMeas_WJets_HT100to200.Clone()
for hist in [hMeas_WJets_HT200to400,hMeas_WJets_HT400to600,hMeas_WJets_HT600to800,hMeas_WJets_HT800to1200,hMeas_WJets_HT1200to2500,hMeas_WJets_HT2500toInf] :
    hMeas_WJets.Add( hist )
    
hMeas_QCD = hMeas_qcd.Clone()
for hist in [hMeas_qcd_tt, hMeas_qcd_tt_nonsemi,
             hMeas_qcd_T_t,hMeas_qcd_Tbar_t,hMeas_qcd_T_tW,hMeas_qcd_Tbar_tW, hMeas_qcd_T_s,
             hMeas_qcd_WJets_HT100to200,hMeas_qcd_WJets_HT200to400,hMeas_qcd_WJets_HT400to600,
             hMeas_qcd_WJets_HT600to800,hMeas_qcd_WJets_HT800to1200,hMeas_qcd_WJets_HT1200to2500,hMeas_qcd_WJets_HT2500toInf] :
    hMeas_QCD.Add(hist,-1.0)
noNegBins(hMeas_QCD)
    
# -------------------------------
# Normalize QCD to MC prediction
# -------------------------------

f_QCD_HT500to700   = TFile("histfiles_full2016/hists_QCD_HT500to700_"+muOrEl+"_nom_post.root")
f_QCD_HT700to1000  = TFile("histfiles_full2016/hists_QCD_HT700to1000_"+muOrEl+"_nom_post.root")
f_QCD_HT1000to1500 = TFile("histfiles_full2016/hists_QCD_HT1000to1500_"+muOrEl+"_nom_post.root")
f_QCD_HT1500to2000 = TFile("histfiles_full2016/hists_QCD_HT1500to2000_"+muOrEl+"_nom_post.root")
f_QCD_HT2000toInf  = TFile("histfiles_full2016/hists_QCD_HT2000toInf_"+muOrEl+"_nom_post.root")

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

# -------------------------------
# Troubleshoot
# -------------------------------

print 'Bin             & Data             & $tt nonsemi     & SingleTop      & WJets          & QCD           & Fake fraction '
for ibin in xrange(1,nbinsMeas+1):
    print string.ljust('$[{:.0f},{:.0f}]$'.format(thisMeas.GetXaxis().GetBinLowEdge(ibin), thisMeas.GetXaxis().GetBinUpEdge(ibin)),16) + '&' + string.rjust('{:.1f}'.format(thisMeas.GetBinContent(ibin)),6)        + ' $\pm$ ' + string.rjust('{:.1f}'.format(thisMeas.GetBinError(ibin)),4) + ' & ' + string.rjust('{:.1f}'.format(hMeas_tt_nonsemi.GetBinContent(ibin)*0.78),4) + ' $\pm$ ' + string.rjust('{:.1f}'.format(hMeas_tt_nonsemi.GetBinError(ibin)*0.78),3) + ' & ' + string.rjust('{:.1f}'.format(hMeas_SingleTop.GetBinContent(ibin)*1.34),4)  + ' $\pm$ ' + string.rjust('{:.1f}'.format(hMeas_SingleTop.GetBinError(ibin)*1.34),3) + ' & ' + string.rjust('{:.1f}'.format(hMeas_WJets.GetBinContent(ibin)*0.96),4) + ' $\pm$ ' + string.rjust('{:.1f}'.format(hMeas_WJets.GetBinError(ibin)*0.96),3) + ' & ' + string.rjust('{:.1f}'.format(hMeas_QCD.GetBinContent(ibin)*0.77),4) + ' $\pm$ ' + string.rjust('{:.1f}'.format(hMeas_QCD.GetBinError(ibin)*0.77),3) + ' & ' + '{:.3f}'.format(response.Vfakes()[ibin-1] / response.Vmeasured()[ibin-1])

print 'Troubleshoot 1: bin 2 has content ' + str(thisMeas.GetBinContent(2)) + '+-' + str(thisMeas.GetBinError(2))

# -------------------------------
# Construct total background hist, with posterior normalizations
# -------------------------------

hBkg = hMeas_tt_nonsemi.Clone()
hBkg.Scale(0.78)
hBkg.Add(hMeas_SingleTop,1.34)
hBkg.Add(hMeas_WJets,0.96)
if options.lepType is "mu":
    hBkg.Add(hMeas_QCD,0.77)
else:
    hBkg.Add(hMeas_QCD,0.74)

# -------------------------------
# Convert for TUnfold
# -------------------------------

# Get 'ttbar' part of thisMeas 
thisMeas.Add(hBkg,-1.0)

for ibin in xrange(1,nbinsMeas+1):
    # Remove fake part in each bin, using fake fractions from response matrix
    fakefraction = response.GetBinContent(ibin,0) / response.Integral(ibin,ibin,0,response.GetYaxis().GetNbins()+1)
    thisMeas.SetBinContent(ibin,thisMeas.GetBinContent(ibin)*(1.0-fakefraction))
    thisMeas.SetBinError(ibin,thisMeas.GetBinError(ibin)*(1.0-fakefraction))
    # And then clear underflow bins, since TUnfold does not expect them
    response.SetBinContent(ibin,0,0)
    response.SetBinError(ibin,0,0)

# Now add back background contribution to thisMeas
thisMeas.Add(hBkg)

# -------------------------------------------------------------------------------------
# Do unfolding
# -------------------------------------------------------------------------------------

unfold = TUnfoldDensity(Hres,TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintArea, TUnfoldDensity.kDensityModeBinWidth)
unfold.SetInput(thisMeas)

# Add systematic uncertainties
for sysname in sysnames:
    unfold.AddSysError(Hres_sys[sysname],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)

# Subtract backgrounds with appropriate scale / uncertainty
unfold.SubtractBackground(hMeas_tt_nonsemi,"TT_nonsemi",0.78,0.09)
unfold.SubtractBackground(hMeas_SingleTop,"SingleTop",1.34,0.58)
unfold.SubtractBackground(hMeas_WJets,"WJets",0.96,0.42)
if options.lepType == "mu":
    unfold.SubtractBackground(hMeas_QCD,"QCD",0.77,0.70)
else :
    unfold.SubtractBackground(hMeas_QCD,"QCD",0.74,0.67)

# Do unfolding
unfold.DoUnfold(0)

print "chi**2=" + str(unfold.GetChi2A()) + "+" + str(unfold.GetChi2L()) + " / " + str(unfold.GetNdf())

# unfolded distribution (histogram)
thisReco = unfold.GetOutput("reco")
thisReco.Sumw2()

# -------------------------------------------------------------------------------------
#Plot error breakdown
# -------------------------------------------------------------------------------------

systs = ['JEC','JER','BTag','TopTag','lep','PDF','Q2','AS']
backgrounds = ['TT_nonsemi','SingleTop','WJets','QCD']


#Statistical -- input and unfolding matrix (GetEmatrixSysUncorr() and GetEmatrixInput())
h_STAT = thisTrue.Clone("stat")
h_STAT.Reset()
#Experimental -- experimental systematics, plus background normalizations (GetEmatrixSysSource() for relevant, then GetEmatrixSysBackgroundUncorr() + GetEmatrixSysBackgroundScale()
h_EXP = thisTrue.Clone("exp")
h_EXP.Reset()
#Theory -- theory systematics (GetEmatrixSysSource() for relevant)
h_TH = thisTrue.Clone("th")
h_TH.Reset()
#Total is GetEmatrixTotal()
h_TOT = thisTrue.Clone("tot")
h_TOT.Reset()

h_SYS = {}
for syst in systs:
    h_SYS[syst] = thisTrue.Clone(sysname)
    h_SYS[syst].Reset()

h_BKG = thisTrue.Clone(sysname)
h_BKG.Reset()

h_LUMI = thisTrue.Clone()
h_LUMI.Reset()

#Get actual error matrices / hists
hErrInput = unfold.GetEmatrixInput("mErrInput")
hErrStat = unfold.GetEmatrixSysUncorr("mErrStat")
hErrTot = unfold.GetEmatrixTotal("mErrTot")

hErrBkgStat = {}
hErrBkgScale = {}
for background in backgrounds:
    hErrBkgStat[background] = unfold.GetEmatrixSysBackgroundUncorr(background,"mErrBkgStat_"+background)
    hErrBkgScale[background] = unfold.GetDeltaSysBackgroundScale(background,"hErrBkgScale_"+background)

hErrSys = {}
for sysname in sysnames:
    hErrSys[sysname] = unfold.GetDeltaSysSource(sysname,"hErrSys_"+sysname)
  
# Fill uncertainty histograms
for ibin in xrange(1,nbinsTrue+1):
    h_STAT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)+hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
    h_TOT.SetBinContent(ibin,math.sqrt(hErrTot.GetBinContent(ibin,ibin)/pow(thisReco.GetBinContent(ibin),2)+0.026*0.026))
    h_LUMI.SetBinContent(ibin,0.026)
    for syst in systs:
        h_SYS[syst].SetBinContent(ibin,(abs(hErrSys[syst+'Up'].GetBinContent(ibin))+abs(hErrSys[syst+'Down'].GetBinContent(ibin)))/ (2*thisReco.GetBinContent(ibin)))
    tot_bkg = 0.0
    for background in backgrounds:
        tot_bkg += (hErrBkgStat[background].GetBinContent(ibin,ibin)+pow(hErrBkgScale[background].GetBinContent(ibin),2))/pow(thisReco.GetBinContent(ibin),2)
    h_BKG.SetBinContent(ibin,math.sqrt(tot_bkg))

# Plot error breakdown
c6 = TCanvas("c6", "", 800, 600)
c6.SetTopMargin(0.08)
c6.SetRightMargin(0.05)
c6.SetBottomMargin(0.14)
c6.SetLeftMargin(0.16)

h_TOT.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
h_TOT.GetYaxis().SetTitle("Uncertainty [%]")
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

colors = [632,600,617,417,432,801,864,906]
markers = [20,21,22,23,33,24,25,26]
for isys in xrange(0,len(systs)):
    h_SYS[systs[isys]].SetLineColor(colors[isys])
    h_SYS[systs[isys]].SetLineWidth(2)
    h_SYS[systs[isys]].SetMarkerColor(colors[isys])
    h_SYS[systs[isys]].SetMarkerStyle(markers[isys])

h_BKG.SetLineColor(419)
h_BKG.SetLineWidth(2)
h_BKG.SetMarkerColor(419)
h_BKG.SetMarkerStyle(32)

longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","PDF Uncertainty","#mu_{R}, #mu_{F} scales","#alpha_{s} scale"]
  
leg6 = TLegend(0.2,0.65,0.45,0.88)
leg6.AddEntry(h_TOT,"Total syst. uncertainty","f")
leg6.AddEntry(h_STAT,"Statistical uncertainty","lp")
leg6.AddEntry(h_LUMI,"Int. luminosity","lp")
for isys in xrange(0,3):
    leg6.AddEntry(h_SYS[systs[isys]],longnames[isys],"lp")

leg6.SetFillStyle(0);
leg6.SetBorderSize(0);
leg6.SetTextSize(0.04);
leg6.SetTextFont(42);

leg66 = TLegend(0.55,0.65,0.8,0.88);
for isys in xrange(3,8):
    leg66.AddEntry(h_SYS[systs[isys]],longnames[isys],"lp")
leg66.AddEntry(h_BKG,"Background normalization","lp")

leg66.SetFillStyle(0)
leg66.SetBorderSize(0)
leg66.SetTextSize(0.04)
leg66.SetTextFont(42)

h_TOT.GetYaxis().SetRangeUser(0.0,1.0)
h_TOT.Draw("hist")
h_STAT.Draw("ep,same")
h_LUMI.Draw("ep,same")
for syst in systs:
    h_SYS[syst].Draw("ep,same")
h_BKG.Draw("ep,same")

leg6.Draw(); 
leg66.Draw(); 

c6.SaveAs("UnfoldingPlots/unfold_relative_uncertainties_"+options.lepType+".pdf")


# Fetch / fill error matrices
diags = []
for ibin in xrange(0,nbinsTrue):
    tmp = math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) if (math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) is not 0) else 1.0
    diags.append(tmp)
sum_weight = 0.0
sum_err = 0.0
for ibinx in xrange(1,nbinsTrue+1):
    for ibiny in xrange(1,nbinsTrue+1):
        tmp = hErrTot.GetBinContent(ibinx,ibiny) / diags[ibinx-1] / diags[ibiny-1]
        hErrTot.SetBinContent(ibinx,ibiny,tmp)

c4 = TCanvas()
hErrTot.GetXaxis().SetRangeUser(400.,1199.)
hErrTot.GetYaxis().SetRangeUser(400.,1199.)
hErrTot.GetZaxis().SetRangeUser(-1.0,1.0)
hErrTot.Draw("colz")
c4.SaveAs("UnfoldingPlots/covariance_"+options.lepType+"_data.pdf")

# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------
# TODO: should fix BR

thisTrue.Scale(1.0/(lum*0.438/3.)) # true @ parton level
thisMeas.Scale(1.0/(lum*0.438/3.)) # measured @ reco level
thisReco.Scale(1.0/(lum*0.438/3.)) # unfolded to parton level

print 'Troubleshoot 4: bin 2 has content ' + str(thisMeas.GetBinContent(2)) + '+-' + str(thisMeas.GetBinError(2))

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

print 'Troubleshoot 5: bin 2 has content ' + str(thisMeas.GetBinContent(2)) + '+-' + str(thisMeas.GetBinError(2))

# -------------------------------------------------------------------------------------
# draw parton-level unfolding
# -------------------------------------------------------------------------------------

## ratio of unfolded data to generator-level

hFrac = thisReco.Clone()
hFrac.SetName("hFrac")
hFrac.SetTitle(";Top quark p_{T} (GeV);Data/MC")
hFrac.Divide(thisTrue)

c1 = TCanvas("c", "c", 700, 700)
pad1 =  TPad("pad1","pad1",0,0.3,1,1)
pad1.SetBottomMargin(0.05);
pad1.Draw();
pad1.cd();

thisReco.SetMarkerStyle(21)
thisMeas.SetMarkerStyle(25);

thisReco.GetXaxis().SetRangeUser(400.,1199.)
thisTrue.GetXaxis().SetRangeUser(400.,1199.)
thisMeas.GetXaxis().SetRangeUser(400.,1199.)

xsec_title = ";;d#sigma/dp_{T} [fb/GeV]"

thisReco.SetTitle(xsec_title)
thisReco.GetYaxis().SetTitleOffset(1.2)
thisReco.SetMinimum(0.0)
max = thisTrue.GetMaximum()
max2 = thisReco.GetMaximum()
if max2 > max:
	max = max2
thisReco.SetAxisRange(0,max*1.15,"Y")
thisReco.Draw()
thisTrue.Draw('hist same')
print 'Troubleshoot 6: bin 2 has content ' + str(thisMeas.GetBinContent(2)) + '+-' + str(thisMeas.GetBinError(2))
thisMeas.Draw('same')
thisTrue.UseCurrentStyle()
thisTrue.SetLineColor(4);
thisTrue.GetYaxis().SetTitleSize(25)
thisTrue.GetXaxis().SetLabelSize(0)

leg = TLegend(0.5, 0.55, 0.9, 0.75)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)

tt = TLatex()
tt.SetNDC()
tt.SetTextFont(42)
leg.AddEntry( thisReco, 'Unfolded data', 'p')
leg.AddEntry( thisTrue, 'Generated (Powheg)', 'l')
leg.AddEntry( thisMeas, 'Measured data', 'p')
leg.Draw()

text1 = TLatex()
text1.SetNDC()
text1.SetTextFont(42)
text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 35.9 fb^{-1}, #sqrt{s} = 13 TeV}")

c1.cd();
pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
pad2.SetTopMargin(0.05);
pad2.SetBottomMargin(0.4);
pad2.Draw();
pad2.cd();
pad2.SetGridy()
hFrac.SetMaximum(1.8)
hFrac.SetMinimum(0.2)
hFrac.UseCurrentStyle()
hFrac.GetYaxis().SetTitleSize(25)
hFrac.GetYaxis().SetTitleOffset(2.0)
hFrac.GetXaxis().SetTitleOffset(4.0)
hFrac.GetXaxis().SetLabelSize(25)
hFrac.GetYaxis().SetNdivisions(4,4,0,False)

hFrac.Draw("e")
hFrac.GetXaxis().SetRangeUser(400., 1199.)

c1.Update()

c1.SaveAs("UnfoldingPlots/closure_"+options.lepType+"_data.pdf")

