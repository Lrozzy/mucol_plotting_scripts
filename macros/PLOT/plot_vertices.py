import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLine, TLatex, TMath
from ROOT import kBlack, kBlue, kRed, kGray, kGreen, kWhite
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
                  default="histos_timing.root", help="Name of the ROOT histo file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default="/data", help="Name of the output folder")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")

# Define features here
names = ['vertex_Lxy',
         'vertex_Lz',
         'vertex_score',
         'vertex_minD',
         'vertex_mass',
         'vertex_ntracks']

vars = ['vertex R [mm]',
        'vertex z [mm]',
        'vertex matching score s',
        'vertex distance to sbottom decay position [mm]',
        'vertex mass [GeV]',
        'Number of vertex tracks']

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

for (name, var) in zip(names, vars):

    c1 = TCanvas()

    h = fFile.Get(name)
    h.SetLineColor(kRed)
    h.SetLineWidth(2)
    h.SetTitle("")
    h.GetYaxis().SetTitle("Number of vertices")
    h.GetYaxis().SetTitleOffset(1.13)
    h.GetYaxis().SetTitleSize(0.045)
    h.GetYaxis().SetLabelSize(0.045)
    h.GetXaxis().SetTitleOffset(1.08)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetXaxis().SetLabelSize(0.045)

    h.GetXaxis().SetTitle(var)
    h.Scale(1., "width")

    if "score" in name:
        h.GetXaxis().SetRange(1, h.GetNbinsX() + 1)

    h.SetMinimum(0.1)

    h.Draw("HIST")

    gPad.RedrawAxis()

    if "ntracks" in name or "mass" in name:
        c1.SetLogy()

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

fFile.Close()
