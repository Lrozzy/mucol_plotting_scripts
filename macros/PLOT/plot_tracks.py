import os
import logging
from ROOT import TH1D, TFile, TTree, TColor, TCanvas, TLegend, TLine, TLatex
from ROOT import kBlack, kBlue, kRed, kGray
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import fabs

gROOT.SetBatch(True)

# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="ntup_tracks.root", help="Name of the ROOT file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="/data", help="Name of the output folder")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFileSig = TFile(options.inFile, "READ")

# Define features here
names = ['track_d0', 'track_z0', 'track_pT', 'track_phi', 'track_theta',
         'track_nholes', 'track_Rprod', 'track_nhits', 'track_chi2ndf']
vars = ['track d_{0} [mm]', 'track z_{0} [mm]', 'track p_{T} [GeV]', 'track #phi [rad]', 'track #theta [rad]',
        'Number of holes', 'track r_{prod} [mm]', 'Number of hits', 'Track #chi^{2}/ndf']

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

for (name, var) in zip(names, vars):

    c1 = TCanvas()

    h = fFileSig.Get(name)
    h.SetLineColor(kRed)
    h.SetLineWidth(2)
    h.SetTitle("")
    h.GetYaxis().SetTitle("Number of tracks")
    h.GetYaxis().SetTitleOffset(1.13)
    h.GetYaxis().SetTitleSize(0.045)
    h.GetYaxis().SetLabelSize(0.045)
    h.GetXaxis().SetTitleOffset(1.08)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.045)
    if 'holes' in name:
        h.GetXaxis().SetNdivisions(5)
    if 'chi2' in name:
        h.GetXaxis().SetRangeUser(0., 40.)
    if 'z0' in name:
        h.GetXaxis().SetRangeUser(-5., 5.)
    if 'd0' in name:
        h.GetXaxis().SetRangeUser(-2., 2.)

    h.GetXaxis().SetTitle(var)
    h.Scale(1., "width")

    h_LLP = fFileSig.Get("LLP"+name)
    h_LLP.SetLineColor(kBlue)
    h_LLP.SetLineWidth(2)
    h_LLP.Scale(1., "width")

    h.SetMinimum(0.1)
    if 'z0' in name:
        h.SetMinimum(100)
    if 'd0' in name:
        h.SetMinimum(100)
    if 'theta' in name:
        h.SetMinimum(20)
        h.SetMaximum(50000)
    if 'nholes' in name:
        h.SetMinimum(20)
        h.SetMaximum(50000)
    if 'nhits' in name:
        h.SetMinimum(20)
        h.SetMaximum(50000)
    if 'chi2' in name:
        h.SetMaximum(80000)

    h.Draw("HIST")
    h_LLP.Draw("HISTSAME")

    if 'holes' in name:
        splitLineA = TLine(1, 20, 1, 40000)
        splitLineA.SetLineWidth(2)
        splitLineA.SetLineStyle(2)
        splitLineA.SetLineColor(kGray+2)
        splitLineA.Draw("same")
    if 'chi2' in name:
        splitLineA = TLine(5, 0, 5, 10000)
        splitLineA.SetLineWidth(2)
        splitLineA.SetLineStyle(2)
        splitLineA.SetLineColor(kGray+2)
        splitLineA.Draw("same")
    if 'nhits' in name:
        splitLineA = TLine(6, 20, 6, 20000)
        splitLineA.SetLineWidth(2)
        splitLineA.SetLineStyle(2)
        splitLineA.SetLineColor(kGray+2)
        splitLineA.Draw("same")

    '''    if 'z0' in name:
        splitLineA = TLine(-2, 0, -2, 0.1)
        splitLineA.SetLineWidth(2)
        splitLineA.SetLineStyle(2)
        splitLineA.SetLineColor(kGray+2)
        splitLineA.Draw("same")
        splitLineB = TLine(2, 0, 2, 0.1)
        splitLineB.SetLineWidth(2)
        splitLineB.SetLineStyle(2)
        splitLineB.SetLineColor(kGray+2)
        splitLineB.Draw("same")
    if 'd0' in name:
        splitLineA = TLine(-0.5, 0, -0.5, 0.1)
        splitLineA.SetLineWidth(2)
        splitLineA.SetLineStyle(2)
        splitLineA.SetLineColor(kGray+2)
        splitLineA.Draw("same")
        splitLineB = TLine(0.5, 0, 0.5, 0.1)
        splitLineB.SetLineWidth(2)
        splitLineB.SetLineStyle(2)
        splitLineB.SetLineColor(kGray+2)
        splitLineB.Draw("same")
    '''

    gPad.RedrawAxis()
    c1.SetLogy()
    if 'theta' in name:
        leg = TLegend(.15, .74, .9, .89)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.045)
        leg.SetNColumns(2)
        leg.AddEntry(h, "Standard tracks", "L")
        leg.AddEntry(h_LLP, "LL tracks", "L")
        leg.Draw()
    else:
        leg = TLegend(.6, .7, .9, .85)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.045)
        leg.AddEntry(h, "Standard tracks", "L")
        leg.AddEntry(h_LLP, "LL tracks", "L")
        leg.Draw()

    t4 = TLatex()
    t4.SetTextFont(42)
    t4.SetTextColor(1)
    t4.SetTextSize(0.04)
    t4.SetTextAlign(12)
    t4.SetNDC()
    t4.DrawLatex(
        0.15, 0.95, '#sqrt{s} = 3 TeV #mu^{+}#mu^{-} collisions')
    # t4.DrawLatex(
    #    0.15, 0.95, '#sqrt{s} = 3 TeV #mu^{+}#mu^{-} collisions, #sqrt{s} = 1.5 TeV BIB overlay')

    c1.SaveAs(options.outFolder+"/"+name+".png")

fFileSig.Close()
