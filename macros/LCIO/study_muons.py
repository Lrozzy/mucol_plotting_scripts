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
                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))

h_truthMu_Rprod = TH1D('truthMu_Rprod', 'truthMu_Rprod',
                       len(arrBins_R)-1, arrBins_R)  # mm
h_truthMu_pT = TH1D('truthMu_pT', 'truthMu_pT', len(arrBins_pT)-1, arrBins_pT)
h_truthMu_theta = TH1D('truthMu_theta', 'truthMu_theta',
                       len(arrBins_theta)-1, arrBins_theta)

h_trk_pT = TH1D('trk_pT', 'trk_pT', len(arrBins_pT)-1, arrBins_pT)
h_trk_theta = TH1D('trk_theta', 'trk_theta',
                   len(arrBins_theta)-1, arrBins_theta)
h_trk_Rprod = TH1D('trk_Rprod',
                   'trk_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_muon_pT = TH1D('muon_pT', 'muon_pT', len(arrBins_pT)-1, arrBins_pT)
h_muon_theta = TH1D('muon_theta', 'muon_theta',
                    len(arrBins_theta)-1, arrBins_theta)
h_muon_Rprod = TH1D('muon_Rprod',
                    'muon_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truthMu_Rprod, h_truthMu_pT, h_truthMu_theta,
               h_trk_pT, h_trk_theta, h_trk_Rprod,
               h_muon_pT, h_muon_theta, h_muon_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)


# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# Bfield = 3.56  # T
Bfield = 5  # T

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Processing event " + str(ievt))

    pfoCollection = event.getCollection('PandoraPFOs')
    trkCollection = event.getCollection('SiTracks_Refitted')
    #relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
    #relation = UTIL.LCRelationNavigator(relationCollection)

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

                goodtheta = False
                if tlv.Theta() > 30.*TMath.Pi()/180. and tlv.Theta() < 150.*TMath.Pi()/180.:
                    goodtheta = True

                if tlv.Perp() > 1 and not mcp.isDecayedInTracker() and goodtheta:
                    h_truthMu_Rprod.Fill(rprod)
                    h_truthMu_pT.Fill(tlv.Perp())
                    h_truthMu_theta.Fill(tlv.Theta())

                    #tracks = relation.getRelatedToObjects(mcp)
                    # if len(tracks) > 0:
                    #    track = tracks[0]
                    #    h_trk_Rprod.Fill(rprod)
                    #    h_trk_pT.Fill(tlv.Perp())
                    #    h_trk_theta.Fill(tlv.Theta())

                    for pfo in pfoCollection:
                        if fabs(pfo.getType()) == 13:
                            dp3 = pfo.getMomentum()
                            tlv_pfo = TLorentzVector()
                            tlv_pfo.SetPxPyPzE(
                                dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                            pfotracks = pfo.getTracks()
                            if len(pfotracks) > 0:
                                trk = pfotracks[0]

                            if track == trk:
                                h_muon_Rprod.Fill(rprod)
                                h_muon_pT.Fill(tlv.Perp())
                                h_muon_theta.Fill(tlv.Theta())

reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_muonPFO.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
