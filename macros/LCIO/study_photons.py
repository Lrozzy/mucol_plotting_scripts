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

h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_EMpfo_E = TH1D('EMpfo_E', 'EMpfo_E', 120, 0., 6000.)
h_HADpfo_E = TH1D('HADpfo_E', 'HADpfo_E', 120, 0., 6000.)
h_matchedEMpfo_E = TH1D('matchedEMpfo_E', 'matchedEMpfo_E', 120, 0., 6000.)
h_matchedEMpfo_theta = TH1D('matchedEMpfo_theta', 'matchedEMpfo_theta',
                            100, 0., pi)
h_matchedHADpfo_E = TH1D('matchedHADpfo_E', 'matchedHADpfo_E', 120, 0., 6000.)
h_matchedHADpfo_theta = TH1D('matchedHADpfo_theta', 'matchedHADpfo_theta',
                             100, 0., pi)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)

h_deltaEM_E = TH1D('deltaEM_E', 'deltaEM_E', 250, -1000, 1000)
h_deltaHAD_E = TH1D('deltaHAD_E', 'deltaHAD_E', 250, -1000, 1000)
h_delta_E_sumE = TH1D('delta_E_sumE', 'delta_E_sumE', 250, -1000, 1000)

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

# Low-level digitised hit distributions
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV
h_EMfrac_PFO = TH1D('EMfrac_PFO', 'EMfrac_PFO', 100, 0, 1)  # GeV

# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta,
               h_EMpfo_E,
               h_HADpfo_E,
               h_matchedEMpfo_E, h_matchedEMpfo_theta,
               h_matchedHADpfo_E, h_matchedHADpfo_theta,
               h_deltaEM_E, h_deltaHAD_E, h_delta_E_sumE,
               h_Npfo, h_pfo_type,
               h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
               h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
               h_ECAL_simhit_E, h_HCAL_simhit_E,
               h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac, h_EMfrac_PFO,
               h_ECAL_hit_layer, h_HCAL_hit_layer,
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
        '''
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus() == 1 and len(mcp.getParents()) == 0:
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                h_truth_E.Fill(mcp.getEnergy())
                h_truth_theta.Fill(tlv.Theta())
        '''

        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')
        h_Npfo.Fill(len(pfoCollection))

        # Match true pfo with closest reco PFO in deltaR
        matchedEM_E = -1.
        matchedEM_theta = -1.
        matchedHAD_E = -1.
        matchedHAD_theta = -1.
        allEM_E = 0.
        allHAD_E = 0.

        minDREM = 999999.
        minDRHAD = 999999.

        for pfo in pfoCollection:
            h_pfo_type.Fill(abs(pfo.getType()))
            dp3 = pfo.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            print(pfo.getType())

            if abs(pfo.getType()) == 22:
                allEM_E = allEM_E + pfo.getEnergy()
            elif abs(pfo.getType()) == 2112:
                allHAD_E = allHAD_E + pfo.getEnergy()

            dR = tlv_pfo.DeltaR(tlv)

            if dR < minDREM and abs(pfo.getType()) == 22:
                minDREM = dR
                matchedEM_E = pfo.getEnergy()
                matchedEM_theta = tlv_pfo.Theta()
            if dR < minDRHAD and abs(pfo.getType()) == 2112:
                minDRHAD = dR
                matchedHAD_E = pfo.getEnergy()
                matchedHAD_theta = tlv_pfo.Theta()

        print(allEM_E, allHAD_E, matchedEM_E, matchedHAD_E)

        h_EMpfo_E.Fill(allEM_E)
        h_HADpfo_E.Fill(allHAD_E)

        if matchedEM_E > 0:
            h_matchedEMpfo_E.Fill(matchedEM_E)
            h_matchedEMpfo_theta.Fill(matchedEM_theta)
            h_deltaEM_E.Fill(matchedEM_E-mcpCollection[0].getEnergy())
        if matchedHAD_E > 0:
            h_matchedHADpfo_E.Fill(matchedHAD_E)
            h_matchedHADpfo_theta.Fill(matchedHAD_theta)
            h_deltaHAD_E.Fill(matchedHAD_E-mcpCollection[0].getEnergy())

        if allHAD_E+allEM_E > 0:
            h_EMfrac_PFO.Fill(allEM_E/(allHAD_E+allEM_E))

        # Fill the simhit-level histos and aggregated energy
        '''
        ECALsimhitCollection = event.getCollection('ECalBarrelCollection')
        encoding = ECALsimhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        for simhit in ECALsimhitCollection:
            cellID = int(simhit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_ECAL_simhit_E.Fill(simhit.getEnergy())
            h_ECAL_simhit_layer.Fill(layer, simhit.getEnergy())

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

        HCALsimhitCollection = event.getCollection('HCalBarrelCollection')
        encoding = ECALsimhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        for simhit in HCALsimhitCollection:
            cellID = int(simhit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_HCAL_simhit_E.Fill(simhit.getEnergy())
            h_HCAL_simhit_layer.Fill(layer, simhit.getEnergy())
        '''

        # Fill the hit-level histos and aggregated energy
        ECAL_sumE = 0.
        ECALhitCollection = event.getCollection('ECALBarrel')

        encoding = ECALhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        for hit in ECALhitCollection:
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_ECAL_hit_time.Fill(hit.getTime())
            h_ECAL_hit_E.Fill(hit.getEnergy())
            h_ECAL_hit_layer.Fill(layer, hit.getEnergy())

            ECAL_sumE = ECAL_sumE + hit.getEnergy()
            pos = hit.getPosition()
            h_ECAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))
        h_ECAL_sumE.Fill(ECAL_sumE)

        HCAL_sumE = 0.
        HCALhitCollection = event.getCollection('HCALBarrel')

        encoding = ECALhitCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        for hit in HCALhitCollection:
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder["layer"].value()

            h_HCAL_hit_time.Fill(hit.getTime())
            h_HCAL_hit_E.Fill(hit.getEnergy())
            h_HCAL_hit_layer.Fill(layer, hit.getEnergy())
            HCAL_sumE = HCAL_sumE + hit.getEnergy()
            pos = hit.getPosition()
            h_HCAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))
        h_HCAL_sumE.Fill(HCAL_sumE)

        print(ECAL_sumE, HCAL_sumE)

        h_sumE.Fill(ECAL_sumE+HCAL_sumE)

        if ECAL_sumE+HCAL_sumE > 0:
            h_EMfrac.Fill(ECAL_sumE/(ECAL_sumE+HCAL_sumE))
        else:
            h_EMfrac.Fill(0)
        h_delta_E_sumE.Fill(ECAL_sumE+HCAL_sumE-mcpCollection[0].getEnergy())

    reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_pfoPFO.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
