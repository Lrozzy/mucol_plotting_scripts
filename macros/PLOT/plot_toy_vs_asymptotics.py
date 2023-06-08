#!/usr/bin/env python

from ROOT import TGraph, TH1D, TCanvas, gPad, gROOT, gStyle, TLatex, TLegend, TLine, TMath
from ROOT import kBlack, kWhite, kGray, kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan, kTeal, kGreen, kSpring, kYellow, kOrange
import json
import math
import pyhf
import os
import logging
from array import array

from optparse import OptionParser
import itertools
from math import fabs

import pandas as pd
import numpy as np
import ROOT


def GetSigma(p):

    # double pres limit:
    if p > (1.0-1e-16):
        return -7.4
    # double pres limit:
    if (p < 1e-16):
        return 7.4
    # convert p-value in standard deviations("nsigma")
    nsigma = 0
    if (p > 1.0e-16):
        nsigma = ROOT.RooStats.PValueToSignificance(p)
    else:
        nsigma = -1

    return nsigma


def graphFromFile(dirname, process, selection):
    x_arr = array('d')
    y_arr = array('d')
    good = 0
    # iterate over files in directory
    for filename in os.listdir(dirname):
        f = os.path.join(dirname, filename)
        # checking if it is a file
        if (os.path.isfile(f)) and (process in filename) and (selection in filename):
            file = open(f)
            lines = file.readlines()
            # Get the signal yields first
            for line in lines:
                x_arr.append(float(line.split(" ")[0]))
                y_arr.append(GetSigma(float(line.split(" ")[1])))
                good = good+1

    if len(x_arr)>0:
        df = pd.DataFrame({'x':x_arr, 'y':y_arr})
        df.sort_values('x', axis=0, inplace=True)
        x_arr =  df[['x']].to_numpy()
        y_arr =  df[['y']].to_numpy()

        graph = TGraph(good, x_arr, y_arr)
    else:
        graph = TGraph()
    
    return graph


# Options
parser = OptionParser()
parser.add_option("-f", "--inputFolder",   dest='inputFolder',
                  default="toys", help="Path to folder containing the inputs")
parser.add_option("-o", "--outFile",   dest='outFile',
                  default="p0_vs_lumi", help="Output file name")
(options, args) = parser.parse_args()

gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

c1 = TCanvas()
c1.SetLogx()

frame = TH1D(" ", " ", 100, 1., 100000.)
frame.GetXaxis().SetTitle("Integrated luminosity [fb^{-1}]")
frame.GetXaxis().SetTitleSize(0.05)
frame.GetYaxis().SetTitle("Discovery significance")
frame.GetYaxis().SetTitleSize(0.05)
frame.GetXaxis().SetTitleOffset(1.2)
frame.GetYaxis().SetTitleOffset(0.9)
frame.GetYaxis().SetLabelSize(0.045)
frame.GetXaxis().SetLabelSize(0.045)
frame.SetMaximum(6.)
frame.Draw()

graph_1LW = graphFromFile(options.inputFolder+"_toys", "hino", "1T")
graph_1LW.SetLineColor(kGreen+3)
graph_1LW.SetLineWidth(2)
graph_1LW.SetLineStyle(6)
graph_1LW.Draw("LSAME")

graph_1LH = graphFromFile(options.inputFolder, "hino", "1T")
graph_1LH.SetLineColor(kOrange+7)
graph_1LH.SetLineWidth(2)
graph_1LH.SetLineStyle(6)
graph_1LH.Draw("LSAME")

graph_2LW = graphFromFile(options.inputFolder+"_toys", "hino", "2T")
graph_2LW.SetLineColor(kGreen+3)
graph_2LW.SetLineWidth(2)
graph_2LW.Draw("LSAME")

graph_2LH = graphFromFile(options.inputFolder, "hino", "2T")
graph_2LH.SetLineColor(kOrange+7)
graph_2LH.SetLineWidth(2)
graph_2LH.Draw("LSAME")

leg = TLegend(.18, .5, .4, .75)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)

leg.AddEntry(
    graph_2LW, "SR^{#gamma}_{2t} #tilde{W}, 2.7 TeV, #tau = 0.2 ns", "FL")
leg.AddEntry(
    graph_2LH, "SR^{#gamma}_{2t} #tilde{H}, 1.1 TeV, #tau = 0.02 ns", "FL")
leg.AddEntry(
    graph_1LW, "SR^{#gamma}_{1t} #tilde{W}, 2.7 TeV, #tau = 0.2 ns", "FL")
leg.AddEntry(
    graph_1LH, "SR^{#gamma}_{1t} #tilde{H}, 1.1 TeV, #tau = 0.02 ns", "FL")
leg.Draw()

splitLineA = TLine(1., 5., 100000., 5.)
splitLineA.SetLineWidth(2)
splitLineA.SetLineStyle(2)
splitLineA.SetLineColor(kGray+2)
splitLineA.Draw("same")

t1 = TLatex()
t1.SetTextFont(42)
t1.SetTextColor(kGray+2)
t1.SetTextSize(0.04)
t1.SetTextAlign(12)
t1.SetNDC()
t1.DrawLatex(0.18, 0.81, '5#sigma')

t4 = TLatex()
t4.SetTextFont(42)
t4.SetTextColor(1)
t4.SetTextSize(0.04)
t4.SetTextAlign(12)
t4.SetNDC()
t4.DrawLatex(
    0.15, 0.95, '#sqrt{s} = 10 TeV #mu^{+}#mu^{-} collisions, #sqrt{s} = 1.5 TeV BIB overlay')

c1.SaveAs(options.outFile + ".pdf")
