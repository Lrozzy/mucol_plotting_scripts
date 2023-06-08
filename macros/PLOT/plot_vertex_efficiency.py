import os
import logging
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLine, TLatex, TMath, TEfficiency
from ROOT import kBlack, kBlue, kRed, kGray, kGreen, kWhite
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import fabs
from array import array

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
h_reco_R = fFile.Get('vertex_Lxy')
h_reco_R.Sumw2()
h_tot_R = fFile.Get('truth_sbottom_Lxy')
h_tot_R.Sumw2()

eff_vs_R = TEfficiency(h_reco_R, h_tot_R)
eff_vs_R.SetLineColor(kRed)
eff_vs_R.SetLineWidth(2)

h_reco_Lz = fFile.Get('vertex_Lz')
h_tot_Lz = fFile.Get('truth_sbottom_Lz')

eff_vs_Lz = TH1D(h_reco_Lz)
eff_vs_Lz.Divide(h_tot_Lz)

gStyle.SetOptStat(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

c1 = TCanvas()

arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
h_frame = TH1D('framec1', 'framec1', len(arrBins_R)-1, arrBins_R)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Vertex reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
#h_frame.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame.GetXaxis().SetTitle("Vertex r [mm]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.2)

h_frame.Draw()
eff_vs_R.Draw("SAME")

postick = [31, 51, 74, 102, 127, 340, 554, 819, 1153, 1486]
labtick = ["VXD0-1", "VXD2-3", "VXD4-5", "VXD6-7",
           "IT0", "IT1", "IT2", "OT0", "OT1", "OT2"]
lines = []

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(kGray+3)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetTextAngle(90)

for ipos, pos in enumerate(postick):
    lines.append(TLine(pos, 0., pos, 1.))
    lines[ipos].SetLineStyle(2)
    lines[ipos].SetLineWidth(2)
    lines[ipos].SetLineColor(kGray+3)
    lines[ipos].Draw()
    t3.DrawLatex(pos+8, 0.8, labtick[ipos])


gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.15)

c1.SaveAs(options.outFolder+"/Efficiency_vs_R.png")

c2 = TCanvas()

eff_vs_Lz.SetLineColor(kBlack)
eff_vs_Lz.SetLineWidth(2)
eff_vs_Lz.SetTitle("")
eff_vs_Lz.GetYaxis().SetTitle("Vertex reconstruction efficiency")
eff_vs_Lz.GetYaxis().SetTitleOffset(1.2)
eff_vs_Lz.GetXaxis().SetTitleOffset(1.2)
#eff_vs_Lz.GetXaxis().SetRangeUser(-0.24, 0.5)
eff_vs_Lz.GetXaxis().SetTitle("Vertex Lz [mm]")
eff_vs_Lz.SetMinimum(0.)
eff_vs_Lz.SetMaximum(1.5)

eff_vs_Lz.Draw("E0")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.15)

c2.SaveAs(options.outFolder+"/Efficiency_vs_Lz.png")

fFile.Close()
