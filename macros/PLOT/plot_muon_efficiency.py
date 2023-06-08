import os
import logging
from stringprep import c22_specials
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLine, TLatex, TMath, TEfficiency
from ROOT import kBlack, kBlue, kRed, kGray, kGreen, kWhite
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
                  default="ntup_muonPFO.root", help="Name of the ROOT histo file")
parser.add_option("-o", "--outFolder",   dest='outFolder',
                  default=".", help="Name of the output folder")
(options, args) = parser.parse_args()

# Init the logger for the script and various modules
fFile = TFile(options.inFile, "READ")

gStyle.SetOptStat(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

# Define features here
h_reco = fFile.Get('muon_pT')
h_reco.Sumw2()
h_trk = fFile.Get('trk_pT')
h_trk.Sumw2()
h_tot = fFile.Get('truthMu_pT')
h_tot.Sumw2()

eff_rec = TEfficiency(h_reco, h_tot)
eff_rec.SetLineColor(kRed)
eff_rec.SetLineWidth(2)

eff_trk = TEfficiency(h_trk, h_tot)
eff_trk.SetLineColor(kBlue)
eff_trk.SetLineWidth(2)

eff_max = TEfficiency(h_reco, h_trk)
eff_max.SetLineColor(kGreen+2)
eff_max.SetLineWidth(2)

c1 = TCanvas()

arrBins = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                      3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50.))
h_frame = TH1D('framec1', 'framec1', len(arrBins)-1, arrBins)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
# h_frame.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame.GetXaxis().SetTitle("p_{T} [GeV]")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_rec.Draw("E0SAME")
eff_trk.Draw("E0SAME")
eff_max.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

leg = TLegend(.5, .2, .9, .4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.AddEntry(eff_rec, "Total efficiency", "L")
leg.AddEntry(eff_trk, "Tracking efficiency", "L")
leg.AddEntry(eff_max, "Tracked only", "L")
leg.Draw()

c1.SaveAs(options.outFolder+"/Muon_Efficiency_pT.pdf")

# Define features here
h_reco2 = fFile.Get('muon_theta')
h_reco2.Sumw2()
h_trk2 = fFile.Get('trk_theta')
h_trk2.Sumw2()
h_tot2 = fFile.Get('truthMu_theta')
h_tot2.Sumw2()

eff_rec2 = TEfficiency(h_reco2, h_tot2)
eff_rec2.SetLineColor(kRed)
eff_rec2.SetLineWidth(2)

eff_trk2 = TEfficiency(h_trk2, h_tot2)
eff_trk2.SetLineColor(kBlue)
eff_trk2.SetLineWidth(2)

eff_max2 = h_reco2
eff_max2.Divide(h_trk2)
eff_max2.SetLineColor(kGreen+2)
eff_max2.SetLineWidth(2)

c2 = TCanvas()

arrBins2 = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                       90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))
h_frame = TH1D('framec1', 'framec1', len(arrBins2)-1, arrBins2)

h_frame.SetLineColor(kBlack)
h_frame.SetLineWidth(2)
h_frame.SetTitle("")
h_frame.GetYaxis().SetTitle("Efficiency")
h_frame.GetYaxis().SetTitleOffset(1.2)
h_frame.GetXaxis().SetTitleOffset(1.2)
# h_frame.GetXaxis().SetRangeUser(-0.24, 0.5)
h_frame.GetXaxis().SetTitle("Muon #theta")
h_frame.SetMinimum(0.)
h_frame.SetMaximum(1.05)

h_frame.Draw()
eff_rec2.Draw("E0SAME")
eff_trk2.Draw("E0SAME")
eff_max2.Draw("E0SAME")

gPad.RedrawAxis()
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.15)
gPad.SetRightMargin(0.1)

leg = TLegend(.5, .2, .9, .4)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.AddEntry(eff_rec2, "Total efficiency", "L")
leg.AddEntry(eff_trk2, "Tracking efficiency", "L")
leg.AddEntry(eff_max2, "Tracked only", "L")
leg.Draw()

c2.SaveAs(options.outFolder+"/Muon_Efficiency_theta.pdf")

fFile.Close()
