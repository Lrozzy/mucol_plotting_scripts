from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TFile, TTree, TLorentzVector, TVector3, TMath
import ROOT as r
from math import *
from optparse import OptionParser

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-t', '--trackCollection', help='--trackCollection SelectedTracks',
                  type=str, default='SelectedTracks')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()


def getOriginPID(mcp, stopID):
    # Look for sbottom mothers
    origin_PDGid = 0
    momVec = mcp.getParents()
    while (len(momVec) > 0 and fabs(origin_PDGid) != stopID):
        mc_mother = momVec[0]
        origin_PDGid = mc_mother.getPDG()
        momVec = mc_mother.getParents()

    return origin_PDGid

#########################


tree = TTree("dummy_tree", "tracks_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt = array('d', [0])
d0 = array('d', [0])
z0 = array('d', [0])
mass = array('d', [0])
vxpt = array('d', [0])
vxr = array('d', [0])
ntracks = array('i', [0])
hasmuon = array('i', [0])
target = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt,  'var/D')
tree.Branch("d0", d0, 'var/D')
tree.Branch("z0", z0, 'var/D')
tree.Branch("vxpt",  vxpt,  'var/D')
tree.Branch("vxmass",  mass,  'var/D')
tree.Branch("vxr",  vxr,  'var/D')
tree.Branch("ntracks", ntracks, 'var/I')
tree.Branch("hasmuon", hasmuon, 'var/I')
tree.Branch("target", target, 'var/I')

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    mcpCollection = event.getCollection('MCParticle')
    tracks = event.getCollection('SelectedTracks')
    svCollection = event.getCollection('MySVCollection')
    T2VRelCollection = event.getCollection('Tracks_To_Vertex')
    t2vrelation = UTIL.LCRelationNavigator(T2VRelCollection)

    good_PID = [2, 3, 4, 5]

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        '''
        if fabs(pdg) == 13:
            daughters = mcp.getDaughters()
            if len(daughters) > 0:
                for d in daughters:
                    if d.getPDG() in good_PID:'''
        if fabs(pdg) == 1000005:
            daughters = mcp.getDaughters()
            if len(daughters) > 0:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        # clear all vectors
                        pt[0] = -1.
                        d0[0] = 0
                        z0[0] = 0
                        mass[0] = -1.
                        ntracks[0] = 0
                        vxpt[0] = -1.
                        vxr[0] = -1.
                        hasmuon[0] = 0
                        # set target
                        #target[0] = 0
                        target[0] = 1
                        # get direction of sbottom
                        #momentum = d.getMomentum()
                        momentum = d.getMomentum()
                        tlv_mc = TLorentzVector()
                        tlv_mc.SetPxPyPzE(
                            momentum[0], momentum[1], momentum[2], d.getEnergy())
                        # get all tracks in cone around sbottom
                        for itrack, track in enumerate(tracks):
                            ptt = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                            theta = TMath.Pi()/2-atan(track.getTanLambda())
                            eta = -1.*TMath.Log(tan(theta/2))
                            phi = track.getPhi()

                            tlv_trk = TLorentzVector()
                            tlv_trk.SetPtEtaPhiM(ptt, eta, phi, 0.139570)

                            if fabs(tlv_mc.DeltaR(tlv_trk)) < 0.4:
                                if ptt > pt[0]:
                                    pt[0] = ptt
                                    d0[0] = track.getD0()
                                    z0[0] = track.getZ0()

                        for v in svCollection:
                            vpos = v.getPosition()

                            v3_pos = TVector3(vpos[0], vpos[1], vpos[2])
                            v4_pos = TLorentzVector(v3_pos, 4.8)
                            if fabs(tlv_mc.DeltaR(v4_pos)) < 0.4:
                                tlv_vertex = TLorentzVector()
                                tlv_vertex.SetPtEtaPhiM(0., 0., 0., 0.)

                                tracksAtVertex = t2vrelation.getRelatedToObjects(
                                    v)

                                for trk in tracksAtVertex:
                                    tlv_track = TLorentzVector()
                                    ptv = 0.3 * Bfield / \
                                        fabs(trk.getOmega() * 1000.)
                                    phi = trk.getPhi()
                                    theta = TMath.Pi()/2-atan(trk.getTanLambda())
                                    eta = -1.*TMath.Log(tan(theta/2))

                                    tlv_track.SetPtEtaPhiM(
                                        ptv, eta, phi, 0.139570)
                                    tlv_vertex = tlv_vertex+tlv_track

                                mass[0] = tlv_vertex.M()
                                vxr[0] = sqrt(vpos[0]*vpos[0]+vpos[1]*vpos[1])
                                vxpt[0] = tlv_vertex.Perp()
                                ntracks[0] = len(tracksAtVertex)

                        foundmu = False
                        for mcp2 in mcpCollection:
                            pdg2 = mcp2.getPDG()
                            # if fabs(pdg2) == 13 and getOriginPID(mcp2, d.getPDG()) == d.getPDG():
                            if fabs(pdg2) == 13 and getOriginPID(mcp2, 1000005) == 1000005:
                                momentum2 = mcp2.getMomentum()
                                tlv_mc2 = TLorentzVector()
                                tlv_mc2.SetPxPyPzE(
                                    momentum[0], momentum[1], momentum[2], mcp2.getEnergy())
                                if fabs(tlv_mc.DeltaR(tlv_mc2)) < 0.4 and tlv_mc2.Perp() > 5.:
                                    foundmu = True

                        hasmuon[0] = foundmu

                        tree.Fill()

reader.close()

# write histograms
output_file = TFile(options.outDir + "ml_" +
                    options.trackCollection + ".root", 'RECREATE')
tree.Write()
output_file.Close()
