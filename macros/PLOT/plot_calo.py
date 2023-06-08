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
                  default="histos_calo_v1.root", help="Name of the ROOT histo file")
parser.add_option("-a", "--Alt",   dest='Alt',
                  default="histos_calo_v0B.root", help="Name of the ROOT histo file")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")
fAlt = TFile(options.Alt, "READ")

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# ECAL barrel
h_barrel = fFile.Get('ECAL_barrel_depth')
h_barrel_alt = fAlt.Get('ECAL_barrel_depth')

c1 = TCanvas("", "", 800, 600)

h_barrel.SetLineColor(kBlue+1)
h_barrel.SetLineWidth(2)
h_barrel.SetTitle("")
h_barrel.GetYaxis().SetTitle("Energy density [MeV / cm^{  2}]")
h_barrel.GetYaxis().SetTitleOffset(1.7)
h_barrel.SetMinimum(0.000001)
h_barrel.SetMaximum(0.002)
h_barrel.GetXaxis().SetNdivisions(10)
h_barrel.GetXaxis().SetLabelSize(0.04)
h_barrel.GetXaxis().SetTitleOffset(1.3)
h_barrel.GetXaxis().SetTitle("Hit Radial Depth in ECAL")
h_barrel.Draw("HIST")

h_barrel_alt.SetLineColor(kRed+1)
h_barrel_alt.SetLineWidth(2)
h_barrel_alt.Draw("HISTSAME")

gPad.RedrawAxis()

leg = TLegend(.6, .64, .9, .78)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(h_barrel, "MuColl v1", "L")
leg.AddEntry(h_barrel_alt, "MuColl10 v0B", "L")
leg.Draw()

t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.24, 0.94, 'Background hits overlay in [0, 10] ns range')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c1.SaveAs("ECALBarrelOccupancy.pdf")

# HCAL barrel
h_barrel = fFile.Get('HCAL_barrel_depth')
h_barrel_alt = fAlt.Get('HCAL_barrel_depth')

c1 = TCanvas("", "", 800, 600)

h_barrel.SetLineColor(kBlue+1)
h_barrel.SetLineWidth(2)
h_barrel.SetTitle("")
h_barrel.GetYaxis().SetTitle("Energy density [MeV / cm^{  2}]")
h_barrel.GetYaxis().SetTitleOffset(1.4)
h_barrel.SetMinimum(0.000001)
h_barrel.SetMaximum(0.01)
h_barrel.GetXaxis().SetNdivisions(10)
h_barrel.GetXaxis().SetLabelSize(0.04)
h_barrel.GetXaxis().SetTitleOffset(1.3)
h_barrel.GetXaxis().SetTitle("Hit Radial Depth in ECAL")
h_barrel.Draw("HIST")

h_barrel_alt.SetLineColor(kRed+1)
h_barrel_alt.SetLineWidth(2)
h_barrel_alt.Draw("HISTSAME")

gPad.RedrawAxis()

leg = TLegend(.6, .64, .9, .78)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(h_barrel, "MuColl v1", "L")
leg.AddEntry(h_barrel_alt, "MuColl10 v0B", "L")
leg.Draw()

t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(.62, 0.82, 'B_{solenoid} = 5 T')

t3 = TLatex()
t3.SetTextFont(42)
t3.SetTextColor(1)
t3.SetTextSize(0.035)
t3.SetTextAlign(12)
t3.SetNDC()
t3.DrawLatex(0.24, 0.94, 'Background hits overlay in [0, 10] ns range')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.035)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(0.82, 0.94, '#sqrt{s} = 10 TeV')

c1.SaveAs("HCALBarrelOccupancy.pdf")

fFile.Close()
