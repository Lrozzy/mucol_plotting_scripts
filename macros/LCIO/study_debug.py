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

h_truth_Rprod = TH1D('truth_Rprod', 'truth_Rprod',
                     len(arrBins_R)-1, arrBins_R)  # mm
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)

h_muon_pT = TH1D('muon_pT', 'muon_pT', len(arrBins_pT)-1, arrBins_pT)
h_muon_theta = TH1D('muon_theta', 'muon_theta',
                    len(arrBins_theta)-1, arrBins_theta)
h_muon_Rprod = TH1D('muon_Rprod',
                    'muon_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truth_Rprod, h_truth_pT, h_truth_theta,
               h_muon_pT, h_muon_theta, h_muon_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)


# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

Bfield = 3.56  # T
pt_print_cut = 20.

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt > 9:
        break
        # print("Processing event " + str(ievt))
    '''
    tracks = event.getCollection("Tracks")
    print("Tracks (" + str(len(tracks)) + ")")
    for itrack, track in enumerate(tracks):
        pt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())
        hits = track.getTrackerHits()
        if pt > pt_print_cut:
            print(" " + str(pt) + " " + str(theta) +
                  " " + str(phi) + " NA " + str(len(hits)))

    tracks = event.getCollection("Tracks_LLP")
    print("TracksLLP (" + str(len(tracks)) + ")")
    for itrack, track in enumerate(tracks):
        pt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())
        hits = track.getTrackerHits()
        if pt > pt_print_cut:
            print(" " + str(pt) + " " + str(theta) +
                  " " + str(phi) + " NA " + str(len(hits)))

    tracks = event.getCollection("SelectedTracks")
    print("Selected Tracks (" + str(len(tracks)) + ")")
    for itrack, track in enumerate(tracks):
        pt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())
        hits = track.getTrackerHits()
        if pt > pt_print_cut:
            print(" " + str(pt) + " " + str(theta) +
                  " " + str(phi) + " NA " + str(len(hits)))

    tracks = event.getCollection("SelectedTracks_LLP")
    print("Selected Tracks_LLP (" + str(len(tracks)) + ")")
    for itrack, track in enumerate(tracks):
        pt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())
        hits = track.getTrackerHits()
        if pt > pt_print_cut:
            print(" " + str(pt) + " " + str(theta) +
                  " " + str(phi) + " NA " + str(len(hits)))

    rtracks = event.getCollection("SiTracks_Refitted")
    relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
    relation = UTIL.LCRelationNavigator(relationCollection)
    print("Refitted Tracks (" + str(len(tracks)) + ")")
    for itrack, track in enumerate(rtracks):
        pt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())
        mcpvec = relation.getRelatedFromObjects(track)
        hits = track.getTrackerHits()
        if pt > pt_print_cut:
            if len(mcpvec) > 0:
                print(" " + str(pt) + " " + str(theta) +
                      " " + str(phi) + " " + str(mcpvec[0].getPDG()) + " " + str(len(hits)))
            else:
                print(" " + str(pt) + " " + str(theta) +
                      " " + str(phi) + " NA " + str(len(hits)))
    '''

    print("PFOs")
    # find the last stops
    pfoCollection = event.getCollection('PandoraPFOs')
    for pfo in pfoCollection:
        dp3 = pfo.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

        if tlv.Perp() > pt_print_cut:
            print(" " + str(tlv.Perp()) + " " + str(pfo.getType()))
        '''
        if fabs(pfo.getType()) == 13:
            tracks = pfo.getTracks()
            for trk in tracks:
                pt = 0.3 * Bfield / fabs(trk.getOmega() * 1000.)
                print("  Trk " + str(pt))
        '''
    print(" ")

reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_muonPFO.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
