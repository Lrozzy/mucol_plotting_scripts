from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TColor, TH1D, TH2D, TFile, TCanvas, gROOT, gStyle, TLegend, TVector3, TMath, TRandom3, TTree
from math import *
from optparse import OptionParser
import os
import fnmatch
from array import array

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_timing.root',
                  type=str, default='ntup_timing.root')
(options, args) = parser.parse_args()

speedoflight = 299792458/1000000  # mm/ns
'''
allCollections = [
    "InnerTrackerBarrelCollection",
    "VertexBarrelCollection"]
'''
allCollections = [
    "VertexBarrelCollection",
    "VertexEndcapCollection",
    "InnerTrackerBarrelCollection",
    "InnerTrackerEndcapCollection",
    "OuterTrackerBarrelCollection",
    "OuterTrackerEndcapCollection"
]

tree = TTree("hits_tree", "hits_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
time_v = array('d', [0])
detector_v = array('i', [0])
layer_v = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("time",  time_v,  'var/D')
tree.Branch("detector", detector_v, 'var/I')
tree.Branch("layer", layer_v, 'var/I')

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()

reader.open(options.inFile)
for ievt, event in enumerate(reader):

    rndm = TRandom3(ievt)

    if ievt % 10 == 0:
        print("Processing "+str(ievt))

    # do all first
    for coll in allCollections:

        # setting decoder
        hitsCollection = event.getCollection(coll)
        encoding = hitsCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        # looping over hit collections
        for ihit, hit in enumerate(hitsCollection):
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            detector = decoder["system"].value()
            layer = decoder["layer"].value()

            position = hit.getPosition()  # mm
            d = sqrt(position[0]*position[0] + position[1]
                     * position[1] + position[2]*position[2])
            tof = d/speedoflight

            resolution = 0.03
            if detector > 2:
                resolution = 0.06

            corrected_time = hit.getTime()*(1.+rndm.Gaus(0., resolution)) - tof

            try:
                mcp = hit.getMCParticle()
                pdg = mcp.getPDG()
            except:
                pdg = 0

            if fabs(pdg) == 13:  # if muon
                time_v[0] = corrected_time
                detector_v[0] = detector
                layer_v[0] = layer

                tree.Fill()

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
output_file.Close()
