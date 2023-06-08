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
h_truth_E = TH1D('truth_E', 'truth_E', 120, 0., 6000.)
h_truth_theta = TH1D('truth_theta', 'truth_theta', 100, 0., pi)

# Low-level sim hit distributions
h_ECAL_simhit_E = TH1D('ECAL_simhit_E', 'ECAL_simhit_E', 100, 0, 20)  # GeV
h_ECAL_simhit_layer = TH1D(
    'ECAL_simhit_layer', 'ECAL_simhit_layer', 100, 0, 100)
h_ECAL_simhit_layer_ele = TH1D(
    'ECAL_simhit_layer_ele', 'ECAL_simhit_layer_ele', 100, 0, 100)
h_ECAL_simhit_layer_gamma = TH1D(
    'ECAL_simhit_layer_gamma', 'ECAL_simhit_layer_gamma', 100, 0, 100)
h_ECAL_simhit_layer_other = TH1D(
    'ECAL_simhit_layer_other', 'ECAL_simhit_layer_other', 100, 0, 100)

h_HCAL_simhit_E = TH1D('HCAL_simhit_E', 'HCAL_simhit_E', 100, 0, 20)  # GeV
h_HCAL_simhit_layer = TH1D(
    'HCAL_simhit_layer', 'HCAL_simhit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta,
               h_ECAL_simhit_E, h_HCAL_simhit_E,
               h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac,
               h_ECAL_simhit_layer, h_ECAL_simhit_layer_ele, h_ECAL_simhit_layer_gamma, h_ECAL_simhit_layer_other,
               h_HCAL_simhit_layer
               ]

for histo in histos_list:
    histo.SetDirectory(0)

to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # loop over all events in the file
    for ievt, event in enumerate(reader):

        if ievt % 1 == 0:
            print(" ")
            print("Processing event " + str(ievt))

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        h_truth_E.Fill(mcpCollection[0].getEnergy())
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        print("True Photon", mcpCollection[0].getEnergy(), tlv.Theta())

        # Fill the simhit-level histos and aggregated energy
        ECAL_sumE = 0.
        ECALsimhitCollection = event.getCollection('ECalBarrelCollection')
        encoding = ECALsimhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        print("N ECAL simhits:", len(ECALsimhitCollection))

        for simhit in ECALsimhitCollection:
            cellID = int(simhit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_ECAL_simhit_E.Fill(simhit.getEnergy())
            h_ECAL_simhit_layer.Fill(layer, simhit.getEnergy())
            ECAL_sumE = ECAL_sumE + simhit.getEnergy()

            E_ele = 0
            E_gamma = 0
            E_other = 0
            for c in range(0, simhit.getNMCParticles()):
                if abs(simhit.getPDGCont(c)) == 11:
                    E_ele = E_ele + simhit.getEnergyCont(c)
                elif abs(simhit.getPDGCont(c)) == 22:
                    E_gamma = E_gamma + simhit.getEnergyCont(c)
                else:
                    E_other = E_other + simhit.getEnergyCont(c)
            h_ECAL_simhit_layer_ele.Fill(layer, E_ele)
            h_ECAL_simhit_layer_gamma.Fill(layer, E_gamma)
            h_ECAL_simhit_layer_other.Fill(layer, E_other)

        h_ECAL_sumE.Fill(ECAL_sumE)

        HCAL_sumE = 0.
        HCALsimhitCollection = event.getCollection('HCalBarrelCollection')
        encoding = ECALsimhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        print("N HCAL simhits:", len(HCALsimhitCollection))

        for simhit in HCALsimhitCollection:
            cellID = int(simhit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_HCAL_simhit_E.Fill(simhit.getEnergy())
            h_HCAL_simhit_layer.Fill(layer, simhit.getEnergy())
            HCAL_sumE = HCAL_sumE + simhit.getEnergy()

        h_HCAL_sumE.Fill(HCAL_sumE)

        print(ECAL_sumE, HCAL_sumE)
        h_sumE.Fill(ECAL_sumE+HCAL_sumE)

        if ECAL_sumE+HCAL_sumE > 0:
            h_EMfrac.Fill(ECAL_sumE/(ECAL_sumE+HCAL_sumE))
        else:
            h_EMfrac.Fill(0)

    reader.close()

# write histograms
output_file = TFile(options.outDir + "simhit.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
