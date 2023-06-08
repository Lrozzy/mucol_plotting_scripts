from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
import os
from array import array

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile test.slcio',
                  type=str, default='test.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_vertex.root',
                  type=str, default='ntup_vertex.root')
(options, args) = parser.parse_args()

Bfield = 3.56  # T


def isFromSbottom(mcp):
    # Look for sbottom mothers
    origin_PDGid = 0
    momVec = mcp.getParents()
    while (len(momVec) > 0 and fabs(origin_PDGid) != 1000005):
        mc_mother = momVec[0]
        origin_PDGid = mc_mother.getPDG()
        momVec = mc_mother.getParents()

    return origin_PDGid

#########################


def find_sbottom_decay(list_particles):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if pdg == 1000005:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found sbottom with no daughters')
                continue
            else:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        result.append(mcp)
    return result


def find_asbottom_decay(list_particles):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if pdg == -1000005:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found sbottom with no daughters')
                continue
            else:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        result.append(mcp)
    return result


#########################
# declare histograms
truth_sbottom_Lxy = TH1D(
    'truth_sbottom_Lxy', 'truth_sbottom_Lxy', 50, 0., 500.)  # mm
truth_sbottom_Lz = TH1D(
    'truth_sbottom_Lz', 'truth_sbottom_Lz', 50, 0., 500.)  # mm
vertex_Lxy = TH1D('vertex_Lxy', 'vertex_Lxy', 50, 0., 500.)  # mm
vertex_Lz = TH1D('vertex_Lz', 'vertex_Lz', 50, 0., 500.)  # mm
vertex_S = TH1D('vertex_score', 'vertex_score', 20, 0., 1.)  # unitless
vertex_minD = TH1D('vertex_minD', 'vertex_minD', 30, 0., 150.)  # mm
vertex_mass = TH1D('vertex_mass', 'vertex_mass', 20, 0., 200.)  # GeV
vertex_ntracks = TH1D('vertex_ntracks', 'vertex_ntracks', 25, 0., 25.)  # GeV

histos_list = [truth_sbottom_Lxy, truth_sbottom_Lz,
               vertex_Lxy, vertex_Lz, vertex_S, vertex_minD, vertex_mass, vertex_ntracks]

for histo in histos_list:
    histo.SetDirectory(0)

#########################
tree = TTree("vertex_tree", "vertex_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
x = array('d', [0])
y = array('d', [0])
z = array('d', [0])
mass = array('d', [0])
matchscore = array('d', [0])
mindistance = array('d', [0])
ntracks = array('i', [0])
pdgID = array('i', [0])

tree.Branch("x",  x,  'var/D')
tree.Branch("y",  y,  'var/D')
tree.Branch("z",  z,  'var/D')
tree.Branch("mass",  mass,  'var/D')
tree.Branch("matchscore",  matchscore,  'var/D')
tree.Branch("mindistance",  mindistance,  'var/D')
tree.Branch("ntracks", ntracks, 'var/I')
tree.Branch("pdgID", pdgID, 'var/I')

#########################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Processing event " + str(ievt))

    # get mc particle collection and loop over it
    mcpCollection = event.getCollection('MCParticle')

    MC2TRelCollection = event.getCollection('MCParticle_Tracks_LLP')
    mc2trelation = UTIL.LCRelationNavigator(MC2TRelCollection)

    T2VRelCollection = event.getCollection('Tracks_To_Vertex')
    t2vrelation = UTIL.LCRelationNavigator(T2VRelCollection)

    sbottom_decays = find_sbottom_decay(mcpCollection)
    asbottom_decays = find_asbottom_decay(mcpCollection)

    if len(sbottom_decays) > 1:
        print("Bad stuff!")
        continue

    if len(asbottom_decays) > 1:
        print("Bad stuff!")
        continue

    for sbot in sbottom_decays:
        end = sbot.getEndpoint()
        truth_sbottom_Lxy.Fill(sqrt(end[0]*end[0]+end[1]*end[1]))
        truth_sbottom_Lz.Fill(end[2])

    for asbot in asbottom_decays:
        end = asbot.getEndpoint()
        truth_sbottom_Lxy.Fill(sqrt(end[0]*end[0]+end[1]*end[1]))
        truth_sbottom_Lz.Fill(end[2])

    # get mc particle collection and loop over it
    svCollection = event.getCollection('MySVCollection')
    for v in svCollection:
        vpos = v.getPosition()

        tlv_vertex = TLorentzVector()
        tlv_vertex.SetPtEtaPhiM(0., 0., 0., 0.)

        tracksAtVertex = t2vrelation.getRelatedToObjects(v)
        #print("Found " + str(len(tracksAtVertex)) + " tracks")

        s_sbot = 0
        s_asbot = 0
        den = 0
        s_score = 0

        for trk in tracksAtVertex:

            tlv_track = TLorentzVector()
            pt = 0.3 * Bfield / fabs(trk.getOmega() * 1000.)

            phi = trk.getPhi()
            theta = TMath.Pi()/2-atan(trk.getTanLambda())
            eta = -1.*TMath.Log(tan(theta/2))

            tlv_track.SetPtEtaPhiM(pt, eta, phi, 0.139570)
            tlv_vertex = tlv_vertex+tlv_track

            den = den + pt

            relparts = mc2trelation.getRelatedFromObjects(trk)
            #print(" Found " + str(len(relparts)) + " related particles")

            if relparts:
                mcp = mc2trelation.getRelatedFromObjects(trk)[0]
                origin = isFromSbottom(mcp)

                if origin == 1000005:
                    s_sbot = s_sbot+pt
                if origin == -1000005:
                    s_asbot = s_asbot+pt

        min_distance = 999999999.
        s_score = s_sbot/den
        pdg = 1000005

        if s_asbot > s_sbot:
            s_score = s_asbot/den
            pdg = -1000005

            # check match with sbottom decays
            for sbot in asbottom_decays:
                end = sbot.getEndpoint()
                distX = vpos[0] - end[0]
                distY = vpos[1] - end[1]
                distZ = vpos[2] - end[2]
                dist = sqrt(distX*distX+distY*distY+distZ*distZ)
                if dist < min_distance:
                    min_distance = dist

        else:
            # check match with sbottom decays
            for sbot in sbottom_decays:
                end = sbot.getEndpoint()
                distX = vpos[0] - end[0]
                distY = vpos[1] - end[1]
                distZ = vpos[2] - end[2]
                dist = sqrt(distX*distX+distY*distY+distZ*distZ)
                if dist < min_distance:
                    min_distance = dist

        x[0] = vpos[0]
        y[0] = vpos[1]
        z[0] = vpos[2]
        mass[0] = tlv_vertex.M()
        matchscore[0] = s_score
        mindistance[0] = min_distance
        ntracks[0] = len(tracksAtVertex)
        pdgID[0] = pdg

        tree.Fill()

        vertex_S.Fill(s_score)
        vertex_minD.Fill(min_distance)
        vertex_mass.Fill(tlv_vertex.M())
        vertex_ntracks.Fill(len(tracksAtVertex))
        if min_distance < 5 and s_score > 0.5:
            vertex_Lxy.Fill(sqrt(vpos[0]*vpos[0]+vpos[1]*vpos[1]))
            vertex_Lz.Fill(vpos[2])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
for histo in histos_list:
    histo.Write()
output_file.Close()
