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
parser.add_option('-t', '--trackCollection', help='--trackCollection SelectedTracks',
                  type=str, default='SelectedTracks')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 1 == 0:
        print("Processing event " + str(ievt))

    # find the last stops
    jetCollection = event.getCollection('Jets')
    for jet in jetCollection:
        dp3 = jet.getMomentum()

        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

        print(str(tlv.Perp()) + " " + str(len(jet.getParticles())))
        for constituent in jet.getParticles():
            ids = constituent.getParticleIDs()
            print(str(constituent.getCharge()) + " " +
                  str(len(constituent.getTracks())))

            if len(ids) > 0:
                print(str(ids[0]))

reader.close()
