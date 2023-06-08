import os
import logging
from ROOT import TH1D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine
from ROOT import kBlack, kBlue, kRed, kYellow, kGray
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import fabs


def check_output_directory(output_path):
    # checks if output directory exists; if not, mkdir
    if not os.path.exists(str(output_path)):
        os.makedirs(output_path)


# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="histos_occupancy.root", help="Name of the ROOT histo file")
parser.add_option('-o', '--outFile', help='--outFile TrackerOccupancy.pdf',
                  type=str, default='TrackerOccupancy.pdf')
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")

# Define features here
h_all = fFile.Get('h_nhits')
h_time = fFile.Get('h_ntimehits')

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

c1 = TCanvas("", "", 800, 600)

h_all.SetLineColor(kBlue+1)
h_all.SetFillColor(kBlue+1)
h_all.SetTitle("")
h_all.GetYaxis().SetTitle("Average number of hits / cm^{  2}")
h_all.GetYaxis().SetTitleOffset(1.4)
h_all.SetMinimum(0.0001)
h_all.GetXaxis().SetNdivisions(10)
h_all.GetXaxis().SetLabelSize(0.04)
h_all.GetXaxis().SetTitleOffset(1.3)
h_all.GetXaxis().SetTitle("Tracking Detector Layer")
h_all.Draw("HIST")

h_time.SetLineColor(kYellow+1)
h_time.SetFillColor(kYellow+1)
h_time.Draw("HISTSAME")

gPad.RedrawAxis()

detector_boundaries = [8, 16, 24, 27, 34, 41, 44, 48]
for bound in detector_boundaries:
    splitLine = TLine(bound, 0, bound, 600)
    splitLine.SetLineWidth(2)
    splitLine.SetLineStyle(2)
    splitLine.SetLineColor(kGray+1)
    splitLine.DrawClone("same")

leg = TLegend(.6, .64, .9, .78)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(h_all, "No time window", "F")
leg.AddEntry(h_time, "Time window [-3 #sigma_{t}, 5 #sigma_{t}]", "F")
leg.Draw()

# t1 = TLatex()
# t1.SetTextFont(42)
# t1.SetTextColor(1)
# t1.SetTextSize(0.04)
# t1.SetTextAlign(12)
# t1.SetNDC()
# t1.DrawLatex(0.6, 0.85, '#bf{Muon Collider}')
#
# t1_2 = TLatex()
# t1_2.SetTextFont(42)
# t1_2.SetTextColor(1)
# t1_2.SetTextSize(0.04)
# t1_2.SetTextAlign(12)
# t1_2.SetNDC()
# t1_2.DrawLatex(0.6, 0.8, '#it{Simulation}')

t2 = TLatex()
t2.SetTextFont(42)
t2.SetTextColor(1)
t2.SetTextSize(0.035)
t2.SetTextAlign(12)
t2.SetNDC()
t2.DrawLatex(0.25, 0.85, '#sigma_{t}^{VXD}   = 30 ps')

t2_2 = TLatex()
t2_2.SetTextFont(42)
t2_2.SetTextColor(1)
t2_2.SetTextSize(0.035)
t2_2.SetTextAlign(12)
t2_2.SetNDC()
t2_2.DrawLatex(0.25, 0.79, '#sigma_{t}^{IT, OT} = 60 ps')

t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(0.25, 0.73, 'B_{solenoid} = 5 T')

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.15, 0.94, 'Background hits overlay in [-360, 480] ps range')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c1.SaveAs(options.outFile)

fFile.Close()
