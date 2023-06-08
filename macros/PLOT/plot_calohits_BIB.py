import os
import logging
from ROOT import TH1D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from optparse import OptionParser
import itertools
from math import *


def compose_plot(histname, titleX, titleY, fFile_10, fFile_50, fFile_100, fFile_500, fFile_1000, fFile_5000, rangeX=None):

    c = TCanvas("", "", 800, 600)

    hist_10 = fFile_10.Get(histname)
    hist_10.SetDirectory(0)
    hist_50 = fFile_50.Get(histname)
    hist_50.SetDirectory(0)
    hist_100 = fFile_100.Get(histname)
    hist_100.SetDirectory(0)
    hist_500 = fFile_500.Get(histname)
    hist_500.SetDirectory(0)
    hist_1000 = fFile_1000.Get(histname)
    hist_1000.SetDirectory(0)
    hist_5000 = fFile_5000.Get(histname)
    hist_5000.SetDirectory(0)

    hist_10.SetLineColor(kBlue+1)
    hist_10.SetLineWidth(2)
    hist_10.SetTitle("")
    hist_10.GetYaxis().SetTitle(titleY)
    hist_10.GetYaxis().SetTitleOffset(1.7)
    if rangeX:
        hist_10.GetXaxis().SetRangeUser(rangeX[0], rangeX[1])
    hist_10.GetXaxis().SetLabelSize(0.04)
    hist_10.GetXaxis().SetTitleOffset(1.3)
    hist_10.GetXaxis().SetTitle(titleX)
    hist_10.DrawNormalized("HIST")

    hist_50.SetLineColor(kRed+1)
    hist_50.SetLineWidth(2)
    hist_50.DrawNormalized("HISTSAME")

    hist_100.SetLineColor(kGreen+1)
    hist_100.SetLineWidth(2)
    hist_100.DrawNormalized("HISTSAME")

    hist_500.SetLineColor(kOrange+1)
    hist_500.SetLineWidth(2)
    hist_500.DrawNormalized("HISTSAME")

    hist_1000.SetLineColor(kGray+1)
    hist_1000.SetLineWidth(2)
    hist_1000.DrawNormalized("HISTSAME")

    hist_5000.SetLineColor(kYellow+1)
    hist_5000.SetLineWidth(2)
    hist_5000.DrawNormalized("HISTSAME")

    gPad.RedrawAxis()

    leg = TLegend(.6, .54, .9, .82)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(hist_10, "10 GeV", "L")
    leg.AddEntry(hist_50, "50 GeV", "L")
    leg.AddEntry(hist_100, "100 GeV", "L")
    leg.AddEntry(hist_500, "500 GeV", "L")
    leg.AddEntry(hist_1000, "1000 GeV", "L")
    leg.AddEntry(hist_5000, "5000 GeV", "L")

    return c, leg


def check_output_directory(output_path):
    # checks if output directory exists; if not, mkdir
    if not os.path.exists(str(output_path)):
        os.makedirs(output_path)


# Load files, hardcoding FTW...
fFile_10 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_10_BIB.root", "READ")
fFile_50 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_50_BIB.root", "READ")
fFile_100 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_100_BIB.root", "READ")
fFile_500 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_500_BIB.root", "READ")
fFile_1000 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_1000_BIB.root", "READ")
fFile_5000 = TFile(
    "/Users/fmeloni/MuonCollider/MuCData/ntup_photons/ntup_photonGun_5000_BIB.root", "READ")

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# ECAL barrel
c1, leg = compose_plot("ECAL_hit_layer", "Calorimeter layer [0.6 X_{0}]", "Energy deposit per layer [%]", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 50])

leg.Draw()

t2_3 = TLatex()
t2_3.SetTextFont(42)
t2_3.SetTextColor(1)
t2_3.SetTextSize(0.035)
t2_3.SetTextAlign(12)
t2_3.SetNDC()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
t2_3.DrawLatex(.62, 0.5, '3.32 X_{0} before calorimeter')
t2_3.DrawLatex(.67, 0.45, 'Digitised hits')
t2_3.DrawLatex(.67, 0.4, 'Cell threshold 50 keV')

c1.SaveAs("ECAL_eloss_BIB.pdf")

# Gun energy
c2, leg = compose_plot("truth_E", "True energy [GeV]", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 5100])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c2.SaveAs("photonGun_energy_BIB.pdf")

# Gun theta
c3, leg = compose_plot("truth_theta", "True #theta", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, pi])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c3.SaveAs("photonGun_theta_BIB.pdf")

# PFO energy
c4, leg = compose_plot("matchedEMpfo_E", "Photon PFO energy [GeV]", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 5100])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c4.SaveAs("photonGun_PFOenergy_BIB.pdf")

# PFO theta
c5, leg = compose_plot("matchedEMpfo_theta", "Photon PFO #theta", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, pi])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c5.SaveAs("photonGun_PFOtheta_BIB.pdf")

# PFO type
c6, leg = compose_plot("pfo_type", "PFO PID", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 2200])
c6.SetLogy()
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c6.SaveAs("photonGun_PFOtype_BIB.pdf")

# PFO multiplicity
c7, leg = compose_plot("Npfo", "Number of PFOs", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 10])
c7.SetLogy()
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c7.SaveAs("photonGun_nPFO_BIB.pdf")

# ECAL sum energy
c8, leg = compose_plot("ECAL_sumE", "Total ECAL energy [GeV]", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 5100])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c8.SaveAs("photonGun_ECAL_sumE_BIB.pdf")

# HCAL sum energy
c8, leg = compose_plot("HCAL_sumE", "Total HCAL energy [GeV]", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [0, 500])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c8.SaveAs("photonGun_HCAL_sumE_BIB.pdf")

# delta E sum E
c9, leg = compose_plot("delta_E_sumE", "sumE - true energy [GeV]", "Events", fFile_10, fFile_50,
                       fFile_100, fFile_500, fFile_1000, fFile_5000, [-1000, 10])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c9.SaveAs("photonGun_deltaE_sumE_BIB.pdf")

# delta E sum E
c10, leg = compose_plot("deltaEM_E", "PFO Energy - true energy [GeV]", "Events", fFile_10, fFile_50,
                        fFile_100, fFile_500, fFile_1000, fFile_5000, [-1000, 10])
leg.Draw()
t2_3.DrawLatex(.62, 0.84, 'Photon particle gun')
c10.SaveAs("photonGun_deltaE_PFO_BIB.pdf")

fFile_10.Close()
fFile_50.Close()
fFile_100.Close()
fFile_500.Close()
fFile_1000.Close()
fFile_5000.Close()
