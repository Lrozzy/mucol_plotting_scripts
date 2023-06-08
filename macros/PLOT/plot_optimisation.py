import os
import logging
from stringprep import c22_specials
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLine, TLatex, TMath, TEfficiency
from ROOT import kBlack, kBlue, kRed, kGray, kGreen, kWhite, kAzure
from ROOT import gStyle, gPad, gROOT
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
                  default=".", help="Name of the output folder")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")

gStyle.SetOptStat(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

# Define features here
h_seed_R = fFile.Get('seed_Rprod')
h_seed_R.Sumw2()
h_seedLLP_R = fFile.Get('seedLLP_Rprod')
h_seedLLP_R.Sumw2()
h_track_R = fFile.Get('track_Rprod')
h_track_R.Sumw2()
h_trackLLP_R = fFile.Get('trackLLP_Rprod')
h_trackLLP_R.Sumw2()
h_merge_R = fFile.Get('merged_Rprod')
h_merge_R.Sumw2()
h_reco_R = fFile.Get('refitted_Rprod')
h_reco_R.Sumw2()
h_tot_R = fFile.Get('truth_Rprod')
h_tot_R.Sumw2()

eff_seed = TEfficiency(h_seed_R, h_tot_R)
eff_seedLLP = h_seedLLP_R.Clone()
eff_seedLLP.Divide(h_tot_R)
#eff_seedLLP = TEfficiency(h_seedLLP_R, h_tot_R)
eff_trk = TEfficiency(h_track_R, h_tot_R)
eff_trkLLP = TEfficiency(h_trackLLP_R, h_tot_R)
eff_merge = TEfficiency(h_merge_R, h_tot_R)
eff_std = TEfficiency(h_reco_R, h_tot_R)

eff_seed.SetLineColor(kGreen+1)
eff_seed.SetLineWidth(2)
eff_seedLLP.SetLineColor(kGreen+2)
eff_seedLLP.SetLineWidth(2)
eff_trk.SetLineColor(kAzure)
eff_trk.SetLineWidth(2)
eff_trkLLP.SetLineColor(kAzure+1)
eff_trkLLP.SetLineWidth(2)
eff_merge.SetLineColor(kGray+2)
eff_merge.SetLineWidth(2)
eff_std.SetLineColor(kBlack)
eff_std.SetLineWidth(2)

c1 = TCanvas()

arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
h_frame = TH1D('framec1', 'framec1', len(arrBins_R)-1, arrBins_R)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetRangeUser(0., 170.)
h_frame.GetXaxis().SetTitle("r_{prod} [mm]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_seed.Draw("E0SAME")
eff_seedLLP.Draw("E0SAME")
eff_trk.Draw("E0SAME")
eff_trkLLP.Draw("E0SAME")
eff_merge.Draw("E0SAME")
eff_std.Draw("E0SAME")

postick = [31, 51, 74, 102, 127, 340, 554, 819, 1153, 1486]
labtick = ["VXD0/1", "VXD2/3", "VXD4/5", "VXD6/7",
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
    t3.DrawLatex(pos+4, 0.02, labtick[ipos])

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

leg = TLegend(.61, .35, .92, .62)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
leg.AddEntry(eff_seed, "Seed tracks", "L")
leg.AddEntry(eff_seedLLP, "Seed LLP tracks", "L")
leg.AddEntry(eff_trk, "Tracks", "L")
leg.AddEntry(eff_trkLLP, "LLP tracks", "L")
leg.AddEntry(eff_merge, "Merged tracks", "L")
leg.AddEntry(eff_std, "Refitted tracks", "L")
leg.Draw()

c1.SaveAs(options.outFolder+"/Track_Efficiency_vs_R.pdf")

'''
c2 = TCanvas()

# Define features here
h_reco = fFile.Get('track_pT')
h_reco.Sumw2()
h_tot = fFile.Get('truth_pT')
h_tot.Sumw2()

eff_std = TEfficiency(h_reco, h_tot)
eff_std.SetLineColor(kBlue)
eff_std.SetLineWidth(2)

arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
h_frame = TH1D('framec2', 'framec2', len(arrBins_pT)-1, arrBins_pT)

h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitle("Truth particle p_{T} [GeV]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_std.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

c2.SaveAs(options.outFolder+"/Track_Efficiency_vs_pT.pdf")

c3 = TCanvas()
# Define features here
h_reco = fFile.Get('track_theta')
h_reco.Sumw2()
h_tot = fFile.Get('truth_theta')
h_tot.Sumw2()

eff_std = TEfficiency(h_reco, h_tot)
eff_std.SetLineColor(kBlue)
eff_std.SetLineWidth(2)

arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))
h_frame = TH1D('framec3', 'framec3', len(arrBins_theta)-1, arrBins_theta)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
#h_frame.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame.GetXaxis().SetTitle("Truth particle #theta [rad]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_std.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

c3.SaveAs(options.outFolder+"/Track_Efficiency_vs_theta.pdf")

fFile.Close()


gROOT.SetBatch(True)

# Options
parser = OptionParser()
parser.add_option("-i", "--inFile",   dest='inFile',
                  default="histos_timing.root", help="Name of the ROOT histo file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default=".", help="Name of the output folder")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")

gStyle.SetOptStat(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

# Define features here
h_reco_R = fFile.Get('track_Rprod')
h_reco_R.Sumw2()
h_tot_R = fFile.Get('truth_Rprod')
h_tot_R.Sumw2()

eff_std = TEfficiency(h_reco_R, h_tot_R)
eff_std.SetLineColor(kBlue)
eff_std.SetLineWidth(2)

c1 = TCanvas()

arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
h_frame = TH1D('framec1', 'framec1', len(arrBins_R)-1, arrBins_R)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetRangeUser(0., 170.)
h_frame.GetXaxis().SetTitle("r_{prod} [mm]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_std.Draw("E0SAME")

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
gPad.SetRightMargin(0.1)

c1.SaveAs(options.outFolder+"/Track_Efficiency_vs_R.pdf")

c2 = TCanvas()

# Define features here
h_reco = fFile.Get('track_pT')
h_reco.Sumw2()
h_tot = fFile.Get('truth_pT')
h_tot.Sumw2()

eff_std = TEfficiency(h_reco, h_tot)
eff_std.SetLineColor(kBlue)
eff_std.SetLineWidth(2)

arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
h_frame = TH1D('framec2', 'framec2', len(arrBins_pT)-1, arrBins_pT)

h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitle("Truth particle p_{T} [GeV]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_std.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

c2.SaveAs(options.outFolder+"/Track_Efficiency_vs_pT.pdf")

c3 = TCanvas()
# Define features here
h_reco = fFile.Get('track_theta')
h_reco.Sumw2()
h_tot = fFile.Get('truth_theta')
h_tot.Sumw2()

eff_std = TEfficiency(h_reco, h_tot)
eff_std.SetLineColor(kBlue)
eff_std.SetLineWidth(2)

arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))
h_frame = TH1D('framec3', 'framec3', len(arrBins_theta)-1, arrBins_theta)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Track reconstruction efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
#h_frame.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame.GetXaxis().SetTitle("Truth particle #theta [rad]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_std.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

c3.SaveAs(options.outFolder+"/Track_Efficiency_vs_theta.pdf")
'''

fFile.Close()
