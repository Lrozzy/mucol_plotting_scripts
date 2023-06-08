from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TVector3, TTree, TMath
from math import *
from optparse import OptionParser

#########################
# parameters

Bfield = 3.56  # T
doHits = False

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos',
                  type=str, default='histos')
(options, args) = parser.parse_args()


def getOriginPID(mcp):
    # Look for sbottom mothers
    origin_PDGid = 0
    momVec = mcp.getParents()
    while (len(momVec) > 0 and fabs(origin_PDGid) != 1000005):
        mc_mother = momVec[0]
        origin_PDGid = mc_mother.getPDG()
        momVec = mc_mother.getParents()

    return origin_PDGid

#########################
# declare histograms


arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))

h_truth_Rprod = TH1D('truth_Rprod', 'truth_Rprod',
                     len(arrBins_R)-1, arrBins_R)  # mm
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)
h_truth_phi = TH1D('truth_phi', 'truth_phi', 20, -TMath.Pi(), TMath.Pi())

h_track_pT = TH1D('track_pT', 'track_pT', len(arrBins_pT)-1, arrBins_pT)
h_track_phi = TH1D('track_phi', 'track_phi', 20, -TMath.Pi(), TMath.Pi())
h_track_theta = TH1D('track_theta', 'track_theta',
                     len(arrBins_theta)-1, arrBins_theta)
h_track_Rprod = TH1D('track_Rprod',
                     'track_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_trackLLP_pT = TH1D('trackLLP_pT', 'trackLLP_pT',
                     len(arrBins_pT)-1, arrBins_pT)
h_trackLLP_phi = TH1D('trackLLP_phi', 'trackLLP_phi',
                      20, -TMath.Pi(), TMath.Pi())
h_trackLLP_theta = TH1D('trackLLP_theta', 'trackLLP_theta',
                        len(arrBins_theta)-1, arrBins_theta)
h_trackLLP_Rprod = TH1D('trackLLP_Rprod',
                        'trackLLP_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_merged_pT = TH1D('merged_pT', 'merged_pT',
                   len(arrBins_pT)-1, arrBins_pT)
h_merged_phi = TH1D('merged_phi', 'merged_phi',
                    20, -TMath.Pi(), TMath.Pi())
h_merged_theta = TH1D('merged_theta', 'merged_theta',
                      len(arrBins_theta)-1, arrBins_theta)
h_merged_Rprod = TH1D('merged_Rprod',
                      'merged_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_seed_pT = TH1D('seed_pT', 'seed_pT', len(arrBins_pT)-1, arrBins_pT)
h_seed_phi = TH1D('seed_phi', 'seed_phi', 20, -TMath.Pi(), TMath.Pi())
h_seed_theta = TH1D('seed_theta', 'seed_theta',
                    len(arrBins_theta)-1, arrBins_theta)
h_seed_Rprod = TH1D('seed_Rprod',
                    'seed_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_seedLLP_pT = TH1D('seedLLP_pT', 'seedLLP_pT', len(arrBins_pT)-1, arrBins_pT)
h_seedLLP_phi = TH1D('seedLLP_phi', 'seedLLP_phi', 20, -TMath.Pi(), TMath.Pi())
h_seedLLP_theta = TH1D('seedLLP_theta', 'seedLLP_theta',
                       len(arrBins_theta)-1, arrBins_theta)
h_seedLLP_Rprod = TH1D('seedLLP_Rprod',
                       'seedLLP_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_refitted_pT = TH1D('refitted_pT', 'refitted_pT',
                     len(arrBins_pT)-1, arrBins_pT)
h_refitted_phi = TH1D('refitted_phi', 'refitted_phi',
                      20, -TMath.Pi(), TMath.Pi())
h_refitted_theta = TH1D('refitted_theta', 'refitted_theta',
                        len(arrBins_theta)-1, arrBins_theta)
h_refitted_Rprod = TH1D('refitted_Rprod',
                        'refitted_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truth_Rprod, h_truth_pT, h_truth_theta, h_truth_phi,
               h_track_pT, h_track_phi, h_track_theta, h_track_Rprod,
               h_trackLLP_pT, h_trackLLP_phi, h_trackLLP_theta, h_trackLLP_Rprod,
               h_seed_pT, h_seed_phi, h_seed_theta, h_seed_Rprod,
               h_seedLLP_pT, h_seedLLP_phi, h_seedLLP_theta, h_seedLLP_Rprod,
               h_merged_pT, h_merged_phi, h_merged_theta, h_merged_Rprod,
               h_refitted_pT, h_refitted_phi, h_refitted_theta, h_refitted_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)


tree = TTree("tracks_tree", "tracks_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt = array('d', [0])
pt_truth = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
d0 = array('d', [0])
z0 = array('d', [0])
chi2 = array('d', [0])
ndf = array('i', [0])
nhits = array('i', [0])
hits_layer = array('i', [-1, -1, -1, -1, -1, -1, -1, -
                   1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
hits_detector = array('i', [-1, -1, -1, -1, -1, -1, -
                      1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
nholes = array('i', [0])
r_truth = array('d', [0])
isLLP = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt,  'var/D')
tree.Branch("pTtruth",  pt_truth,  'var/D')
tree.Branch("phi", phi, 'var/D')
tree.Branch("theta", theta, 'var/D')
tree.Branch("d0", d0, 'var/D')
tree.Branch("z0", z0, 'var/D')
tree.Branch("chi2", chi2, 'var/D')
tree.Branch("ndf", ndf, 'var/I')
tree.Branch("nhits", nhits, 'var/I')
tree.Branch("hits_layer", hits_layer, 'myArrint[20]/I')
tree.Branch("hits_detector", hits_detector, 'myArrint[20]/I')
tree.Branch("nholes", nholes, 'var/I')
tree.Branch("r_truth", r_truth, 'var/D')
tree.Branch("isLLP", isLLP, 'var/I')

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    # Look at particles from MC
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

                    origin_PDGid = 0
                    momVec = mcp.getParents()
                    if len(momVec) > 0:
                        mc_mother = momVec[0]
                        origin_PDGid = mc_mother.getPDG()

                    if fabs(origin_PDGid) == 1000013 or fabs(origin_PDGid) == 2000013:
                        h_truth_Rprod.Fill(rprod)
                        h_truth_pT.Fill(tlv.Perp())
                        h_truth_phi.Fill(tlv.Phi())
                        h_truth_theta.Fill(tlv.Theta())

                    relationCollection = event.getCollection(
                        'MCParticle_SeedTracks')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    seed = relation.getRelatedToObjects(mcp)
                    for track in seed:
                        h_seed_Rprod.Fill(rprod)
                        h_seed_pT.Fill(tlv.Perp())
                        h_seed_phi.Fill(tlv.Phi())
                        h_seed_theta.Fill(tlv.Theta())

                    relationCollection = event.getCollection(
                        'MCParticle_SeedTracks_LLP')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    seedLLP = relation.getRelatedToObjects(mcp)
                    for track in seedLLP:
                        h_seedLLP_Rprod.Fill(rprod)
                        h_seedLLP_pT.Fill(tlv.Perp())
                        h_seedLLP_phi.Fill(tlv.Phi())
                        h_seedLLP_theta.Fill(tlv.Theta())

                    relationCollection = event.getCollection(
                        'MCParticle_Tracks')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    tracks = relation.getRelatedToObjects(mcp)
                    for track in tracks:
                        h_track_Rprod.Fill(rprod)
                        h_track_pT.Fill(tlv.Perp())
                        h_track_phi.Fill(tlv.Phi())
                        h_track_theta.Fill(tlv.Theta())

                        pt[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                        phi[0] = track.getPhi()
                        theta[0] = TMath.Pi()/2-atan(track.getTanLambda())
                        d0[0] = track.getD0()
                        z0[0] = track.getZ0()
                        chi2[0] = track.getChi2()
                        ndf[0] = track.getNdf()
                        nhits[0] = len(track.getTrackerHits())
                        nholes[0] = int(track.getdEdxError())  # BADHACK
                        pt_truth[0] = tlv.Perp()
                        r_truth[0] = rprod
                        isLLP[0] = 0

                        # reset array
                        for i in range(len(hits_layer)):
                            hits_layer[i] = -1
                            hits_detector[i] = -1
                        
                        if doHits:
                            # setting decoder
                            vertexHitsCollection = event.getCollection('VBTrackerHits')
                            encoding = vertexHitsCollection.getParameters(
                            ).getStringVal(EVENT.LCIO.CellIDEncoding)
                            decoder = UTIL.BitField64(encoding)

                            trkhitsCollection = track.getTrackerHits()
                            for ihit, hit in enumerate(trkhitsCollection):
                                cellID = int(hit.getCellID0())
                                decoder.setValue(cellID)
                                hits_layer[ihit] = decoder['layer'].value()
                                hits_detector[ihit] = decoder["system"].value()

                        tree.Fill()

                    relationCollection = event.getCollection(
                        'MCParticle_Tracks_LLP')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    trackLLP = relation.getRelatedToObjects(mcp)
                    for track in trackLLP:
                        h_trackLLP_Rprod.Fill(rprod)
                        h_trackLLP_pT.Fill(tlv.Perp())
                        h_trackLLP_phi.Fill(tlv.Phi())
                        h_trackLLP_theta.Fill(tlv.Theta())

                        pt[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                        phi[0] = track.getPhi()
                        theta[0] = TMath.Pi()/2-atan(track.getTanLambda())
                        d0[0] = track.getD0()
                        z0[0] = track.getZ0()
                        chi2[0] = track.getChi2()
                        ndf[0] = track.getNdf()
                        nhits[0] = len(track.getTrackerHits())
                        nholes[0] = int(track.getdEdxError())  # BADHACK
                        pt_truth[0] = tlv.Perp()
                        r_truth[0] = rprod
                        isLLP[0] = 1
                        
                        # reset array
                        for i in range(len(hits_layer)):
                            hits_layer[i] = -1
                            hits_detector[i] = -1

                        if doHits:
                            # setting decoder
                            vertexHitsCollection = event.getCollection('VBTrackerHits')
                            encoding = vertexHitsCollection.getParameters(
                            ).getStringVal(EVENT.LCIO.CellIDEncoding)
                            decoder = UTIL.BitField64(encoding)

                            trkhitsCollection = track.getTrackerHits()
                            for ihit, hit in enumerate(trkhitsCollection):
                                cellID = int(hit.getCellID0())
                                decoder.setValue(cellID)
                                hits_layer[ihit] = decoder['layer'].value()
                                hits_detector[ihit] = decoder["system"].value()

                        tree.Fill()

                    relationCollection = event.getCollection(
                        'MCParticle_MergedTracks')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    tracks = relation.getRelatedToObjects(mcp)
                    for track in tracks:
                        h_merged_Rprod.Fill(rprod)
                        h_merged_pT.Fill(tlv.Perp())
                        h_merged_phi.Fill(tlv.Phi())
                        h_merged_theta.Fill(tlv.Theta())

                    relationCollection = event.getCollection(
                        'MCParticle_SiTracks_Refitted')
                    relation = UTIL.LCRelationNavigator(relationCollection)
                    tracks = relation.getRelatedToObjects(mcp)
                    for track in tracks:
                        h_refitted_Rprod.Fill(rprod)
                        h_refitted_pT.Fill(tlv.Perp())
                        h_refitted_phi.Fill(tlv.Phi())
                        h_refitted_theta.Fill(tlv.Theta())

reader.close()

# write histograms
output_file = TFile(options.outFile+"_hist.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()

# write tree
output_file2 = TFile(options.outFile+"_ntup.root", 'RECREATE')
tree.Write()
output_file2.Close()
