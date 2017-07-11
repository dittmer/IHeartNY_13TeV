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

parser.add_option('--type', metavar='F', type='string', action='store',
                  default='full',
                  dest='type',
                  help='')

parser.add_option('--toy', metavar='F', type='string', action='store',
                  default='',
                  dest='toy',
                  help='')

parser.add_option('--regType', metavar='F', type='string', action='store',
                  default='scanTau',
                  dest='regType',
                  help='scanTau or LCurve')

# -------------------------------------------------------------------------------------
# load options & set plot style
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

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
  
#gSystem.Load("RooUnfold/libRooUnfold.so")

#from ROOT import RooUnfoldResponse
#from ROOT import RooUnfold
#from ROOT import RooUnfoldBayes
#from ROOT import RooUnfoldSvd
#from ROOT import RooUnfoldTUnfold

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
#        for i in xrange(1,response.GetNbinsMeasured()+1) :
#            ntru += Hres.GetBinContent(i,j)
#        Hres.SetBinContent(0, j, vtru[j-1]-ntru)
#        Hres.SetBinError(0, j, etru[j-1])
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

lum = 35867.0

response_name = "response_pt_split_TH2"
hMeas_name    = "ptRecoTop_split"
hMeasUp_name  = "ptRecoTopMod_split"
hMeasDn_name  = "ptRecoTopModDown_split"
hTrue_name    = "ptGenTop"

# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"

f_ttbar        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_post.root")
f_ttbar_even   = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_even_post.root")
f_ttbar_odd    = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_odd_post.root")
f_ttbar_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_post.root")
f_ttbar_even_2 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_even_post.root")
f_ttbar_odd_2  = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_odd_post.root")

# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

if options.type == "full":
    response = f_ttbar.Get(response_name)
    response.Sumw2()
    response.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
    response2 = f_ttbar_2.Get(response_name)
    response2.Sumw2()
    response2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
    response.Add(response2)
elif options.type == "half" :
    response = f_ttbar_odd.Get(response_name)
    response.Sumw2()
    response.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
    response2 = f_ttbar_odd_2.Get(response_name)
    response2.Sumw2()
    response2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
    response.Add(response2)
    response.Scale(2.0)
elif options.type == "each" : #unfold p1 with p2
    response = f_ttbar_2.Get(response_name)
    response.Sumw2()
else :
    response = f_ttbar_odd.Get(response_name)
    response.Sumw2()
    response.Scale(2.0)
    
TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# Get systematic variations
# -------------------------------------------------------------------------------------
Hres_sys = {}
sysnames = ['JEC','JER','BTag','TopTag','lep','PDF','Q2','AS']
vars = ['Up','Down']

for sysname in sysnames:
    for var in vars:
        f_ttbar_sys        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+sysname+var+"_post.root")
        f_ttbar_sys_odd    = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+sysname+var+"_odd_post.root")
        f_ttbar_sys_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_"+sysname+var+"_post.root")
        f_ttbar_sys_odd_2  = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_"+sysname+var+"_odd_post.root")

        if options.type == "full":
            response_sys = f_ttbar_sys.Get(response_name)
            response_sys.Sumw2()
            response_sys.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
            response_sys2 = f_ttbar_sys_2.Get(response_name)
            response_sys2.Sumw2()
            response_sys2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
            response_sys.Add(response_sys2)
        elif options.type == "half" :
            response_sys = f_ttbar_sys_odd.Get(response_name)
            response_sys.Sumw2()
            response_sys.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
            response_sys2 = f_ttbar_sys_odd_2.Get(response_name)
            response_sys2.Sumw2()
            response_sys2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
            response_sys.Add(response_sys2)
            response_sys.Scale(2.0)
        elif options.type == "each" : #unfold p1 with p2
            response_sys = f_ttbar_sys_2.Get(response_name)
            response_sys.Sumw2()
        else :
            response_sys = f_ttbar_sys_odd.Get(response_name)
            response_sys.Sumw2()
            response_sys.Scale(2.0)
            
        for ibin in xrange(1,response_sys.GetXaxis().GetNbins()+1):
            response_sys.SetBinContent(ibin,0,0)
            response_sys.SetBinError(ibin,0,0)
        
        Hres_sys[sysname+var] = response_sys

    # -------------------------------------------------------------------------------------
    # Troubleshoot systematic variations
    # -------------------------------------------------------------------------------------
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadRightMargin(0.12)
    
    hDiffSysUp = Hres_sys[sysname+"Up"].Clone()
    hDiffSysUp.Add(response,-1.0)
    hDiffSysDn = Hres_sys[sysname+"Down"].Clone()
    hDiffSysDn.Add(response,-1.0)

    c7 = TCanvas()
    hDiffSysUp.GetXaxis().SetRangeUser(400.,1199.)
    hDiffSysUp.GetYaxis().SetRangeUser(400.,1199.)
    hDiffSysUp.Draw("colz")
    c7.SaveAs("UnfoldingPlots/compareSysResponse_"+sysname+"Up_"+options.lepType+"_"+options.type+".pdf")

    hDiffSysDn.GetXaxis().SetRangeUser(400.,1199.)
    hDiffSysDn.GetYaxis().SetRangeUser(400.,1199.)
    hDiffSysDn.Draw("colz")
    c7.SaveAs("UnfoldingPlots/compareSysResponse_"+sysname+"Dn_"+options.lepType+"_"+options.type+".pdf")

    gStyle.SetPadLeftMargin(0.18)
    gStyle.SetPadRightMargin(0.05)


# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

if options.type == "full":
    if options.toy == "up" :
        thisMeas = f_ttbar.Get(hMeasUp_name).Clone()
        thisTrue = f_ttbar.Get(hTrue_name+"Mod").Clone()
        thisMeas2 = f_ttbar_2.Get(hMeasUp_name).Clone()
        thisTrue2 = f_ttbar_2.Get(hTrue_name+"Mod").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar.Get(hMeasDn_name).Clone() 
        thisTrue = f_ttbar.Get(hTrue_name+"ModDown").Clone() 
        thisMeas2 = f_ttbar_2.Get(hMeasDn_name).Clone() 
        thisTrue2 = f_ttbar_2.Get(hTrue_name+"ModDown").Clone() 
    else : 
        thisMeas = f_ttbar.Get(hMeas_name).Clone() 
        thisTrue = f_ttbar.Get(hTrue_name).Clone()
        thisMeas2 = f_ttbar_2.Get(hMeas_name).Clone() 
        thisTrue2 = f_ttbar_2.Get(hTrue_name).Clone()
    thisMeas.Sumw2()
    thisTrue.Sumw2()
    thisMeas2.Sumw2()
    thisTrue2.Sumw2()
    thisMeas.Add(thisMeas2)
    thisTrue.Add(thisTrue2)
    thisMeas.Scale(lum * 831.76 / (77229341. + 78006311. * 1191. / 1192.))
    thisTrue.Scale(lum * 831.76 / (77229341. + 78006311. * 1191. / 1192.))

elif options.type == "half":
    if options.toy == "up" :
        thisMeas = f_ttbar_even.Get(hMeasUp_name).Clone()
        thisTrue = f_ttbar_even.Get(hTrue_name+"Mod").Clone()
        thisMeas2 = f_ttbar_even_2.Get(hMeasUp_name).Clone()
        thisTrue2 = f_ttbar_even_2.Get(hTrue_name+"Mod").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar_even.Get(hMeasDn_name).Clone() 
        thisTrue = f_ttbar_even.Get(hTrue_name+"ModDown").Clone() 
        thisMeas2 = f_ttbar_even_2.Get(hMeasDn_name).Clone() 
        thisTrue2 = f_ttbar_even_2.Get(hTrue_name+"ModDown").Clone() 
    else : 
        thisMeas = f_ttbar_even.Get(hMeas_name).Clone() 
        thisTrue = f_ttbar_even.Get(hTrue_name).Clone() 
        thisMeas2 = f_ttbar_even_2.Get(hMeas_name).Clone() 
        thisTrue2 = f_ttbar_even_2.Get(hTrue_name).Clone()
    thisMeas.Sumw2()
    thisTrue.Sumw2()
    thisMeas2.Sumw2()
    thisTrue2.Sumw2()
    thisMeas.Add(thisMeas2)
    thisTrue.Add(thisTrue2)
    thisMeas.Scale(lum * 831.76 * 2.0 / (77229341. + 78006311. * 1191. / 1192.))
    thisTrue.Scale(lum * 831.76 * 2.0 / (77229341. + 78006311. * 1191. / 1192.))
    
elif options.type == "each" : #unfold p1 with p2
    if options.toy == "up" :
        thisMeas = f_ttbar.Get(hMeasUp_name).Clone()
        thisTrue = f_ttbar.Get(hTrue_name+"Mod").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar.Get(hMeasDn_name).Clone() 
        thisTrue = f_ttbar.Get(hTrue_name+"ModDown").Clone() 
    else : 
        thisMeas = f_ttbar.Get(hMeas_name).Clone() 
        thisTrue = f_ttbar.Get(hTrue_name).Clone()
    thisMeas.Sumw2()
    thisTrue.Sumw2()
    thisMeas.Scale( lum * 831.76 / 77229341.)
    thisTrue.Scale( lum * 831.76 / 77229341.)

else :
    if options.toy == "up" :
        thisMeas = f_ttbar_even.Get(hMeasUp_name).Clone()
        thisTrue = f_ttbar_even.Get(hTrue_name+"Mod").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar_even.Get(hMeasDn_name).Clone() 
        thisTrue = f_ttbar_even.Get(hTrue_name+"ModDown").Clone() 
    else : 
        thisMeas = f_ttbar_even.Get(hMeas_name).Clone() 
        thisTrue = f_ttbar_even.Get(hTrue_name).Clone()
    thisMeas.Sumw2()
    thisTrue.Sumw2()
    thisMeas.Scale( lum * 831.76 * 2.0 / 77229341.)
    thisTrue.Scale( lum * 831.76 * 2.0 / 77229341.)

noNegBins(thisMeas)
noNegBins(thisTrue)
    
thisMeas.SetName("recolevel") 
thisTrue.SetName("truthlevel")

nbinsMeas = thisMeas.GetNbinsX()
nbinsTrue = thisTrue.GetNbinsX()

# -------------------------------------------------------------------------------------
# Convert for TUnfold
# -------------------------------------------------------------------------------------

for ibin in xrange(1,response.GetXaxis().GetNbins()+1):
    fakefraction = response.GetBinContent(ibin,0) / response.Integral(ibin,ibin,0,response.GetYaxis().GetNbins()+1)
    thisMeas.SetBinContent(ibin,thisMeas.GetBinContent(ibin)*(1.0-fakefraction))
    thisMeas.SetBinError(ibin,thisMeas.GetBinError(ibin)*(1.0-fakefraction))
    response.SetBinContent(ibin,0,0)
    response.SetBinError(ibin,0,0)

    print 'Bin ' + str(ibin) + ' has content ' + str(thisMeas.GetBinContent(ibin)) + ' +- ' + str(thisMeas.GetBinError(ibin))

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

'''
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

hDiff_bin = []
hDiffIni_bin = []
h_tau = TH1F("tau_toys",";log(Tau);Toys",100,-4.0,-1.0)

for ibin in xrange(0,nbinsTrue) :
    hDiff_tmp = TH1F("diff_bin"+str(ibin),";(truth-unfolded)/truth; number of events",100,-1,1)
    hDiff_bin.append(hDiff_tmp)
for ibin in xrange(0,nbinsMeas) :
    hDiffIni_tmp = TH1F("diff_bin"+str(ibin),";(meas-toy)/meas; number of events",100,-1,1)
    hDiffIni_bin.append(hDiffIni_tmp)
        
for itoy in xrange(0,ntoys) :
    unfold_tmp = TUnfoldDensity(response,TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintArea, TUnfoldDensity.kDensityModeBinWidth)
    unfold_tmp.SetInput(hToy_i[itoy])
    for sysname in sysnames:
        unfold_tmp.AddSysError(Hres_sys[sysname+"Up"],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
    
    if options.regType == "LCurve" :
        logTauX_tmp = TSpline3()
        logTauY_tmp = TSpline3()
        lCurve_tmp = TGraph()
        iBest_tmp = unfold_tmp.ScanLcurve(30,0.,0.,lCurve_tmp,logTauX_tmp,logTauY_tmp)
        logTau_tmp = Double(0)
        x_tmp = Double(0)
        logTauX_tmp.GetKnot(iBest_tmp,logTau_tmp,x_tmp)
        h_tau.Fill(logTau_tmp)
    elif options.regType == "ScanTau" :
        scanResult_tmp = TSpline3()
        iBest_tmp = unfold_tmp.ScanTau(100,0.0001,0.1,scanResult_tmp,TUnfoldDensity.kEScanTauRhoAvgSys) #CAN CHANGE
        Tau_tmp = Double(0)
        rho_tmp = Double(0)
        scanResult_tmp.GetKnot(iBest_tmp,Tau_tmp,rho_tmp)
        h_tau.Fill(Tau_tmp)
    #else :
    unfold_tmp.DoUnfold(0)
    
    hReco_tmp = unfold_tmp.GetOutput("tmp_output")

    for ibin in range(0, nbinsTrue):
        if thisTrue.GetBinContent(ibin+1) > 0.0:
            hDiff_bin[ibin].Fill( (thisTrue.GetBinContent(ibin+1) - hReco_tmp.GetBinContent(ibin+1)) / thisTrue.GetBinContent(ibin+1))
    for ibin in xrange(0,nbinsMeas):
        if thisMeas.GetBinContent(ibin+1) > 0.0:
            hDiffIni_bin[ibin].Fill((thisMeas.GetBinContent(ibin+1) - hToy_i[itoy].GetBinContent(ibin+1)) / thisMeas.GetBinContent(ibin+1))

for ibin in xrange(0,nbinsTrue) :
    c = TCanvas()
    hDiff_bin[ibin].Draw("same")
    c.SaveAs("UnfoldingPlots/pull"+options.toy+"_"+options.regType+"_"+options.lepType+"_bin"+str(ibin)+"_"+options.type+".pdf")

if options.regType is not "None" :
    c2 = TCanvas()
    h_tau.Draw()
    c2.SaveAs("UnfoldingPlots/TauFromToys_"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type+".pdf")

hBias_pt = thisTrue.Clone() 
hBias_pt.SetName("biasVsPt_"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type)
hBias_pt.Reset()

for ibin in xrange(0,nbinsTrue) :
    mean = hDiff_bin[ibin].GetMean()
    err = hDiff_bin[ibin].GetRMS()

    hBias_pt.SetBinContent(ibin+1, mean)
    hBias_pt.SetBinError(ibin+1, err)
            
ccc = TCanvas()
gPad.SetGridy()
hBias_pt.GetYaxis().SetTitle("Bias")
hBias_pt.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
hBias_pt.SetAxisRange(400,1199,"X")
hBias_pt.SetAxisRange(-0.5,0.5,"Y")
hBias_pt.Draw()
ccc.SaveAs("UnfoldingPlots/bias_vspt"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type+".pdf")

hBiasIni_pt = thisMeas.Clone() 
hBiasIni_pt.SetName("initialBiasVsPt_"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type)
hBiasIni_pt.Reset()

for ibin in xrange(0,nbinsMeas) :
    mean = hDiffIni_bin[ibin].GetMean()
    err = hDiffIni_bin[ibin].GetRMS()

    hBiasIni_pt.SetBinContent(ibin+1, mean)
    hBiasIni_pt.SetBinError(ibin+1, err)
            
c3 = TCanvas()
gPad.SetGridy()
hBiasIni_pt.GetYaxis().SetTitle("Input Bias")
hBiasIni_pt.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
#hBiasIni_pt.SetAxisRange(400,1199,"X")
hBiasIni_pt.SetAxisRange(-0.5,0.5,"Y")
hBiasIni_pt.Draw()
c3.SaveAs("UnfoldingPlots/bias_vspt"+options.toy+"_"+options.lepType+"_"+options.type+".pdf")
'''

# -------------------------------------------------------------------------------------
# Done with toys, doing actual unfolding
# -------------------------------------------------------------------------------------

unfold = TUnfoldDensity(response,TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintArea, TUnfoldDensity.kDensityModeBinWidth)
unfold.SetInput(thisMeas)
for sysname in sysnames:
    unfold.AddSysError(Hres_sys[sysname+"Down"],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)

if options.regType == "LCurve" :
    logTauX = TSpline3()
    logTauY = TSpline3()
    lCurve = TGraph()
    bestLCurve = TGraph(1)
    iBest = unfold.ScanLcurve(30,0.,0.,lCurve,logTauX,logTauY)
    Tau = Double(0)
    x = Double(0)
    y = Double(0)
    logTauX.GetKnot(iBest,Tau,x)
    logTauY.GetKnot(iBest,Tau,y)
    bestLCurve.SetPoint(1,x,y)
    bestLCurve.SetMarkerColor(2)

elif options.regType == "ScanTau" :
    bestTau = TGraph(1)
    scanResult = TSpline3()
    iBest = unfold.ScanTau(100,0.0001,0.1,scanResult,TUnfoldDensity.kEScanTauRhoAvgSys)
    Tau = Double(0)
    rho = Double(0)
    scanResult.GetKnot(iBest,Tau,rho)
    bestTau.SetPoint(1,Tau,rho)
    bestTau.SetMarkerColor(2)

#else :
unfold.DoUnfold(0)

print "chi**2=" + str(unfold.GetChi2A()) + "+" + str(unfold.GetChi2L()) + " / " + str(unfold.GetNdf())

# unfolded distribution (histogram)
thisReco = unfold.GetOutput("reco")
thisReco.Sumw2()

# -------------------------------------------------------------------------------------
#Plot error breakdown
# -------------------------------------------------------------------------------------

#Statistical -- input and unfolding matrix (GetEmatrixSysUncorr() and GetEmatrixInput())
h_STAT = thisTrue.Clone("stat")
h_STAT.Reset()
#Total is GetEmatrixTotal() (+lumi, since missing)
h_TOT = thisTrue.Clone("tot")
h_TOT.Reset()

# Individual error sources
h_SYS = {}
for sysname in sysnames:
    h_SYS[sysname] = thisTrue.Clone(sysname)
    h_SYS[sysname].Reset()

h_INPUT = thisTrue.Clone()
h_INPUT.Reset()

h_MATRIX = thisTrue.Clone()
h_MATRIX.Reset()

h_LUMI = thisTrue.Clone()
h_LUMI.Reset()

#Get actual error matrices / hists
hErrInput = unfold.GetEmatrixInput("mErrInput")
hErrStat = unfold.GetEmatrixSysUncorr("mErrStat") #This might be buggy!!!
hErrTot = unfold.GetEmatrixTotal("mErrTot")

hErrSys = {}
for sysname in sysnames:
    hErrSys[sysname] = unfold.GetDeltaSysSource(sysname,"hErrSys_"+sysname)
  
# Fill uncertainty histograms
for ibin in xrange(1,nbinsTrue+1):
    h_STAT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin)+hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
    h_INPUT.SetBinContent(ibin,math.sqrt(hErrInput.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
    h_MATRIX.SetBinContent(ibin,math.sqrt(hErrStat.GetBinContent(ibin,ibin))/thisReco.GetBinContent(ibin))
    h_TOT.SetBinContent(ibin,math.sqrt(hErrTot.GetBinContent(ibin,ibin)/pow(thisReco.GetBinContent(ibin),2)+0.026*0.026))
    h_LUMI.SetBinContent(ibin,0.026)
    for sysname in sysnames:
        h_SYS[sysname].SetBinContent(ibin,hErrSys[sysname].GetBinContent(ibin)/thisReco.GetBinContent(ibin))

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

h_INPUT.SetLineColor(1)
h_INPUT.SetLineWidth(2)
h_INPUT.SetMarkerColor(1)
h_INPUT.SetMarkerStyle(20)

h_MATRIX.SetLineColor(419)
h_MATRIX.SetLineWidth(2)
h_MATRIX.SetMarkerColor(419)
h_MATRIX.SetMarkerStyle(32)

h_LUMI.SetLineColor(40)
h_LUMI.SetLineWidth(2)
h_LUMI.SetMarkerColor(40)
h_LUMI.SetMarkerStyle(34)

colors = [632,600,617,417,432,801,864,906]
markers = [20,21,22,23,33,24,25,26]
for isys in xrange(0,len(sysnames)):
    h_SYS[sysnames[isys]].SetLineColor(colors[isys])
    h_SYS[sysnames[isys]].SetLineWidth(2)
    h_SYS[sysnames[isys]].SetMarkerColor(colors[isys])
    h_SYS[sysnames[isys]].SetMarkerStyle(markers[isys])

longnames = ["Jet energy scale","Jet energy resolution","b tagging efficiency","t tagging efficiency","Lepton ID","PDF Uncertainty","#mu_{R}, #mu_{F} scales","#alpha_{s} scale"]
  
leg6 = TLegend(0.2,0.65,0.45,0.88)
leg6.AddEntry(h_TOT,"Total syst. uncertainty","f")
leg6.AddEntry(h_INPUT,"Input stat. unc.","lp")
leg6.AddEntry(h_MATRIX,"Response stat. unc.","lp")
leg6.AddEntry(h_LUMI,"Int. luminosity","lp")
for isys in xrange(0,2):
    leg6.AddEntry(h_SYS[sysnames[isys]],longnames[isys],"lp")

leg6.SetFillStyle(0);
leg6.SetBorderSize(0);
leg6.SetTextSize(0.04);
leg6.SetTextFont(42);

leg66 = TLegend(0.55,0.65,0.8,0.88);
for isys in xrange(2,8):
    leg66.AddEntry(h_SYS[sysnames[isys]],longnames[isys],"lp")

leg66.SetFillStyle(0)
leg66.SetBorderSize(0)
leg66.SetTextSize(0.04)
leg66.SetTextFont(42)

h_TOT.GetYaxis().SetRangeUser(0.0,1.0)
h_TOT.Draw("hist")
h_INPUT.Draw("ep,same")
h_MATRIX.Draw("ep,same")
h_LUMI.Draw("ep,same")
for sysname in sysnames:
    h_SYS[sysname].Draw("ep,same")

leg6.Draw(); 
leg66.Draw(); 

c6.SaveAs("UnfoldingPlots/closure_relative_uncertainties_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

# -------------------------------------------------------------------------------------
# Troubleshoot stat. unc.
# -------------------------------------------------------------------------------------
# Get inverse of Vyy
# TUnfold method has bug for version 17.1 -- doing this manually
Vyy = TMatrixD(nbinsMeas,nbinsMeas)
for ibin in xrange(0,nbinsMeas): #rows
    Vyy[ibin][ibin] = pow(thisMeas.GetBinError(ibin+1),2)

VyyInv = TMatrixD(Vyy)
VyyInv.Invert()

#Get A (normalized response matrix)
A = TMatrixD(nbinsMeas,nbinsTrue)
for irow in xrange(0,nbinsMeas):
    for icol in xrange(0,nbinsTrue):
        A[irow][icol] = response.GetBinContent(irow+1,icol+1) / response.Integral(0,nbinsMeas+1,icol+1,icol+1)

#Manually calculate Vxx
AT = TMatrixD(nbinsTrue,nbinsMeas)
AT.Transpose(A)

EInv1 = TMatrixD(nbinsTrue,nbinsMeas)
EInv1.Mult(AT,VyyInv)
EInv = TMatrixD(nbinsTrue,nbinsTrue)
EInv.Mult(EInv1,A)
E = TMatrixD(EInv)
E.Invert()

# Plot relative stat. unc., and compare w/ data
relUnc1 = thisMeas.Clone()
for ibin in xrange(1,nbinsMeas+1): #rows
    relUnc1.SetBinContent(ibin,thisMeas.GetBinError(ibin) / thisMeas.GetBinContent(ibin))

c8 = TCanvas()
h_INPUT.SetTitle(";Top p_{T};Relative stat. unc.")
h_INPUT.SetMaximum(max(h_INPUT.GetMaximum(),relUnc1.GetMaximum())*1.2)
h_INPUT.Draw("hist")
relUnc1.SetLineColor(2)
relUnc1.Draw("hist,same")
leg8 = TLegend(0.2,0.6,0.5,0.8)
leg8.SetBorderSize(0)
leg8.SetTextSize(0.04)
leg8.AddEntry(relUnc1,"Stat. unc. on measurement","l")
leg8.AddEntry(h_INPUT,"Propagated stat. unc.","l")
leg8.Draw()
c8.SaveAs("UnfoldingPlots/compare_rel_unc_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

# -------------------------------------------------------------------------------------
#Plot covariance matrix
# -------------------------------------------------------------------------------------

diags = []
for ibin in xrange(0,nbinsTrue):
    tmp = math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) if (math.sqrt(hErrTot.GetBinContent(ibin+1,ibin+1)) is not 0) else 1.0
    diags.append(tmp)
sum_weight = 0.0
sum_err = 0.0
for ibinx in xrange(1,nbinsTrue+1):
    for ibiny in xrange(1,nbinsTrue+1):
        tmp = hErrTot.GetBinContent(ibinx,ibiny) / diags[ibinx-1] / diags[ibiny-1]
        weight = abs(hErrTot.GetBinContent(ibinx,ibiny))
        hErrTot.SetBinContent(ibinx,ibiny,tmp)
        if ibinx >= ibiny:
            sum_weight += tmp * weight
            sum_err += weight

c4 = TCanvas()
hErrTot.GetXaxis().SetRangeUser(400.,1199.)
hErrTot.GetYaxis().SetRangeUser(400.,1199.)
hErrTot.GetZaxis().SetRangeUser(-1.0,1.0)
hErrTot.Draw("colz")
c4.SaveAs("UnfoldingPlots/closure_covariance_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

'''
# Do manual scan to get bias / stat. unc. / total unc.
h_BiasVsTau    = TH1F("biasVsTau"   ,";log(#tau);Bias"      ,21,-4.0,-1.0)
h_StatUncVsTau = TH1F("statUncVsTau",";log(#tau);Stat. Unc.",21,-4.0,-1.0)
h_TotUncVsTau  = TH1F("totUncVsTau" ,";log(#tau);Total Unc.",21,-4.0,-1.0)
for ii in xrange(0,21):
    thistau = pow(10,-4.0+3.0/21.0*(ii+0.5))
    unfold.DoUnfold(thistau)
    thisReco_tmp = unfold.GetOutput("reco")
    bias = 0.0
    statUnc = 0.0
    totalUnc = 0.0
    for ibin in xrange(1,nbinsTrue+1):
        bias += pow((thisReco_tmp.GetBinContent(ibin) - thisTrue.GetBinContent(ibin))/thisReco_tmp.GetBinContent(ibin),2)
        statUnc += pow(thisReco_tmp.GetBinError(ibin)/thisReco_tmp.GetBinContent(ibin),2)
        totalUnc += unfold.GetEmatrixTotal("totUnc").GetBinContent(ibin,ibin)/pow(thisReco_tmp.GetBinContent(ibin),2)
    h_BiasVsTau.SetBinContent(ii+1,math.sqrt(bias))
    h_StatUncVsTau.SetBinContent(ii+1,math.sqrt(statUnc))
    h_TotUncVsTau.SetBinContent(ii+1,math.sqrt(totalUnc))
c5 = TCanvas()
h_BiasVsTau.SetLineColor(1)
h_StatUncVsTau.SetLineColor(2)
h_TotUncVsTau.SetLineColor(3)
h_BiasVsTau.GetYaxis().SetRangeUser(0.0,1.5)
leg5 = TLegend(0.4, 0.7, 0.6, 0.9)
leg5.SetFillStyle(0)
leg5.SetTextFont(42)
leg5.SetTextSize(0.045)
leg5.SetBorderSize(0)
leg5.AddEntry(h_BiasVsTau, 'Bias', 'l')
leg5.AddEntry(h_StatUncVsTau, 'Stat. Unc.', 'l')
leg5.AddEntry(h_TotUncVsTau, 'Total Unc.', 'l')
h_BiasVsTau.Draw("hist")
h_StatUncVsTau.Draw("hist,same")
h_TotUncVsTau.Draw("hist,same")
leg5.Draw()
c5.SaveAs("UnfoldingPlots/CoarseTauScan_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")
'''

# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------
# TODO: should fix BR

thisTrue.Scale(1.0/(lum*0.438/3.)) # true @ parton level
thisMeas.Scale(1.0/(lum*0.438/3.)) # measured @ reco level
thisReco.Scale(1.0/(lum*0.438/3.)) # unfolded to parton level

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
hFrac.SetTitle(";Top quark p_{T} (GeV);Unfolded/True")
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
hFrac.GetXaxis().SetRangeUser(400., 1199.)

c1.Update()

c1.SaveAs("UnfoldingPlots/closure_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+"_result.pdf")

# -------------------------------------------------------------------------------------
# Plot L-curve scan
# -------------------------------------------------------------------------------------

if not options.regType == "None" :
    c2 = TCanvas("c2", "c2", 700, 700)
    c2.cd()
    c2.UseCurrentStyle()
    if options.regType == "LCurve" :
        lCurve.GetXaxis().SetTitle("log L_{1}")
        lCurve.GetYaxis().SetTitle("log L_{2} / #tau^{2}")
        lCurve.Draw()
        bestLCurve.Draw("*")
    else :
        #dummy = TGraph()
        #dummy.GetXaxis.SetRangeUser(-4.0,-1.0)
        #dummy.GetXaxis.SetTitle("log(#tau)")
        #dummy.GetYaxis.SetRangeUser(0.5,1.0)
        #dummy.GetYaxis.SetTitle("Correlation")
        #dummy.Draw()
        scanResult.Draw("P")
        bestTau.Draw("*")
        
    tl2 = TLatex()
    tl2.SetNDC()
    tl2.SetTextFont(42)
    legend = "log(#tau) = %.3e" % Tau
    tl2.DrawLatex(0.55,0.8,legend)

    c2.SaveAs("UnfoldingPlots/TauScan_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

