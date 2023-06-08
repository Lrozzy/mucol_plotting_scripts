from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath
from math import *
from optparse import OptionParser
import os
import fnmatch

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile test.slcio',
                  type=str, default='test.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos.root',
                  type=str, default='histos.root')
(options, args) = parser.parse_args()

Bfield = 3.57  # T
output_file_name = 'histos_vx_eff'

#########################


def find_sbottom_decay(list_particles):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if fabs(pdg) == 1000005:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found sbottom with no daughters')
                continue
            else:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        result.append(mcp)
    return result


def find_Rhad_decay(list_particles, good_PIDs):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if fabs(pdg) in good_PIDs:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found sbottom with no daughters')
                continue
            else:
                for d in daughters:
                    if fabs(d.getPDG()) == 1000005:
                        result.append(mcp)
    return result


#########################
# declare histograms
truth_Lxy = TH1D('truth_Lxy', 'truth_Lxy', 10, 0., 100.)  # mm
truth_Lxy.SetDirectory(0)

truth_Lz = TH1D('truth_Lz', 'truth_Lz', 10, 0., 100.)  # mm
truth_Lz.SetDirectory(0)

vertex_Lxy = TH1D('vertex_Lxy', 'vertex_Lxy', 10, 0., 100.)  # mm
vertex_Lxy.SetDirectory(0)

vertex_Lz = TH1D('vertex_Lz', 'vertex_Lz', 10, 0., 100.)  # mm
vertex_Lz.SetDirectory(0)
#########################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    print(" ")
    print("New event, #" + str(ievt))

    # get mc particle collection and loop over it
    mcpCollection = event.getCollection('MCParticle')

    sbottom_decays = find_sbottom_decay(mcpCollection)
    nsbottom = len(sbottom_decays)

    good_pid = [1005321, 1000522, 1000512]
    Rhad_decays = find_Rhad_decay(mcpCollection, good_pid)
    nRhad = len(Rhad_decays)

    print('Nsbottom '+str(nsbottom) + ' NRhad ' + str(nRhad))

    for sbot in sbottom_decays:
        end = sbot.getEndpoint()
        print(str(sbot.getPDG()) + " decay at " + str(end[0]) + "  " +
              str(end[1]) + "  " + str(end[2]))
        truth_Lxy.Fill(sqrt(end[0]*end[0]+end[1]*end[1]))
        truth_Lz.Fill(end[2])

    for rhad in Rhad_decays:
        end = rhad.getEndpoint()
        print(str(rhad.getPDG()) + " decay at " + str(end[0]) + "  " +
              str(end[1]) + "  " + str(end[2]))

    # get mc particle collection and loop over it
    svCollection = event.getCollection('MySVCollection')
    for v in svCollection:
        vpos = v.getPosition()
        print("DV at " + str(vpos[0]) + "  " +
              str(vpos[1]) + "  " + str(vpos[2]))
        vertex_Lxy.Fill(sqrt(vpos[0]*vpos[0]+vpos[1]*vpos[1]))
        vertex_Lz.Fill(vpos[2])

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
truth_Lxy.Write()
truth_Lz.Write()
vertex_Lxy.Write()
vertex_Lz.Write()
output_file.Close()
