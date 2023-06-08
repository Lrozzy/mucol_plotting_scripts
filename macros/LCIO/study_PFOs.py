from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()

# declare histograms
arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))

h_truthMu_Rprod = TH1D('truthMu_Rprod', 'truthMu_Rprod',
                       len(arrBins_R)-1, arrBins_R)  # mm
h_truthMu_pT = TH1D('truthMu_pT', 'truthMu_pT', len(arrBins_pT)-1, arrBins_pT)
h_truthMu_theta = TH1D('truthMu_theta', 'truthMu_theta',
                       len(arrBins_theta)-1, arrBins_theta)

h_truthEl_Rprod = TH1D('truthEl_Rprod', 'truthEl_Rprod',
                       len(arrBins_R)-1, arrBins_R)  # mm
h_truthEl_pT = TH1D('truthEl_pT', 'truthEl_pT', len(arrBins_pT)-1, arrBins_pT)
h_truthEl_theta = TH1D('truthEl_theta', 'truthEl_theta',
                       len(arrBins_theta)-1, arrBins_theta)

h_muon_pT = TH1D('muon_pT', 'muon_pT', len(arrBins_pT)-1, arrBins_pT)
h_muon_theta = TH1D('muon_theta', 'muon_theta',
                    len(arrBins_theta)-1, arrBins_theta)
h_muon_Rprod = TH1D('muon_Rprod',
                    'muon_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_ele_pT = TH1D('ele_pT', 'ele_pT', len(arrBins_pT)-1, arrBins_pT)
h_ele_theta = TH1D('ele_theta', 'ele_theta',
                   len(arrBins_theta)-1, arrBins_theta)
h_ele_Rprod = TH1D('ele_Rprod',
                   'ele_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truthMu_Rprod, h_truthMu_pT, h_truthMu_theta,
               h_truthEl_Rprod, h_truthEl_pT, h_truthEl_theta,
               h_muon_pT, h_muon_theta, h_muon_Rprod,
               h_ele_pT, h_ele_theta, h_ele_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)


# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

Bfield = 3.56  # T

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Processing event " + str(ievt))

    pfoCollection = event.getCollection('PandoraPFOs')
    mcpCollection = event.getCollection('MCParticle')

    for mcp in mcpCollection:

        charge = mcp.getCharge()
        status = mcp.getGeneratorStatus()

        if fabs(charge) > 0:
            if fabs(mcp.getPDG()) == 13:
                vx = mcp.getVertex()
                rprod = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

                if tlv.Perp() > 1 and not mcp.isDecayedInTracker():
                    h_truthMu_Rprod.Fill(rprod)
                    h_truthMu_pT.Fill(tlv.Perp())
                    h_truthMu_theta.Fill(tlv.Theta())

                    for pfo in pfoCollection:
                        if fabs(pfo.getType()) == 13:
                            dp3 = pfo.getMomentum()
                            tlv_pfo = TLorentzVector()
                            tlv_pfo.SetPxPyPzE(
                                dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                            if tlv_pfo.DeltaR(tlv) < 1.2:
                                h_muon_Rprod.Fill(rprod)
                                h_muon_pT.Fill(tlv_pfo.Perp())
                                h_muon_theta.Fill(tlv_pfo.Theta())

            if fabs(mcp.getPDG()) == 11:
                vx = mcp.getVertex()
                rprod = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

                if tlv.Perp() > 1 and not mcp.isDecayedInTracker():
                    h_truthEl_Rprod.Fill(rprod)
                    h_truthEl_pT.Fill(tlv.Perp())
                    h_truthEl_theta.Fill(tlv.Theta())

                    for pfo in pfoCollection:
                        if fabs(pfo.getType()) == 11:
                            dp3 = pfo.getMomentum()
                            tlv_pfo = TLorentzVector()
                            tlv_pfo.SetPxPyPzE(
                                dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                            if tlv_pfo.DeltaR(tlv) < 1.2:
                                h_ele_Rprod.Fill(rprod)
                                h_ele_pT.Fill(tlv.Perp())
                                h_ele_theta.Fill(tlv.Theta())

reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_PFO.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
