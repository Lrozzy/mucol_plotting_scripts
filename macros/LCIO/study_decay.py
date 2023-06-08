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
parser.add_option('-i', '--inFile', help='--inFile /data/stop4body/stop4b.slcio',
                  type=str, default='/data/stop4body/stop4b.slcio')
(options, args) = parser.parse_args()
#########################

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    print(" ")
    print("New event, #" + str(ievt))

    # get mc particle collection and loop over it
    try:
        # find the last stops
        mcpCollection = event.getCollection('MCParticle')
    except:
        print("Exception for event " + str(ievt) +
              ": no MC particle collection found!")
        mcpCollection = 0

    good_pid = [1000005, 1005321, 1000522, 1000512]
    for c in mcpCollection:

        if fabs(c.getPDG()) == 1000005:
            parents = c.getParents()
            # print(fabs(parents[0].getPDG()))

            if fabs(parents[0].getPDG()) == 21:

                daughters = c.getDaughters()
                thepid = c.getPDG()
                thenextpid = c.getPDG()

                pos = c.getVertex()
                nextpos = c.getVertex()
                end = c.getEndpoint()
                nextend = c.getEndpoint()

                while len(daughters) > 0:

                    dau_list = []
                    next = []
                    thepid = thenextpid
                    end = nextend
                    pos = nextpos

                    for d in daughters:
                        dau_list.append(d.getPDG())
                        if fabs(d.getPDG()) in good_pid:
                            next = d.getDaughters()
                            thenextpid = d.getPDG()
                            nextend = d.getEndpoint()
                            nextpos = d.getVertex()

                    daughters = next
                    print(str(thepid) + " --> " + str(dau_list) + " (" +
                          str(round(pos[0], 2)) + ", " + str(round(pos[1], 2)) +
                          ", " + str(round(pos[2], 2)) + ") --> (" +
                          str(round(end[0], 2)) + ", " + str(round(end[1], 2)) +
                          ", " + str(round(end[2], 2)) + ")")

reader.close()
