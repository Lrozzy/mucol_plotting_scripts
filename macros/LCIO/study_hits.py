from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TColor, TH1D, TH2D, TFile, TCanvas, gROOT, gStyle, TLegend, TVector3, TMath
from math import *
from optparse import OptionParser
import os
import fnmatch

#########################
# parameters

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_hits.root',
                  type=str, default='histos_hits.root')
(options, args) = parser.parse_args()

hit_dtheta_closest = TH1D('hit_dtheta_closest',
                          'hit_dtheta_closest', 1500, 0., 30.)
C1_hit_dtheta_closest = TH1D(
    'C1_hit_dtheta_closest', 'C1_hit_dtheta_closest', 1500, 0., 30.)

hit_dtheta_closest_L0 = TH1D('hit_dtheta_closest_L0',
                          'hit_dtheta_closest_L0', 1500, 0., 30.)
C1_hit_dtheta_closest_L0 = TH1D(
    'C1_hit_dtheta_closest_L0', 'C1_hit_dtheta_closest_L0', 1500, 0., 30.)

hit_dtheta_closest_L2 = TH1D('hit_dtheta_closest_L2',
                          'hit_dtheta_closest_L2', 1500, 0., 30.)
C1_hit_dtheta_closest_L2 = TH1D(
    'C1_hit_dtheta_closest_L2', 'C1_hit_dtheta_closest_L2', 1500, 0., 30.)

hit_dtheta_closest_L4 = TH1D('hit_dtheta_closest_L4',
                             'hit_dtheta_closest_L4', 1500, 0., 30.)
C1_hit_dtheta_closest_L4 = TH1D(
    'C1_hit_dtheta_closest_L4', 'C1_hit_dtheta_closest_L4', 1500, 0., 30.)

hit_dtheta_closest_L6 = TH1D('hit_dtheta_closest_L6',
                             'hit_dtheta_closest_L6', 1500, 0., 30.)
C1_hit_dtheta_closest_L6 = TH1D(
    'C1_hit_dtheta_closest_L6', 'C1_hit_dtheta_closest_L6', 1500, 0., 30.)

hit_dphi_closest = TH1D('hit_dphi_closest', 'hit_dphi_closest', 1000, 0., 50.)
C1_hit_dphi_closest = TH1D('C1_hit_dphi_closest',
                           'C1_hit_dphi_closest', 1000, 0., 50.)

hit_dphi_closest_L0 = TH1D('hit_dphi_closest_L0', 'hit_dphi_closest_L0', 1000, 0., 50.)
C1_hit_dphi_closest_L0 = TH1D('C1_hit_dphi_closest_L0',
                           'C1_hit_dphi_closest_L0', 1000, 0., 50.)

hit_dphi_closest_L2 = TH1D('hit_dphi_closest_L2', 'hit_dphi_closest_L2', 1000, 0., 50.)
C1_hit_dphi_closest_L2 = TH1D('C1_hit_dphi_closest_L2',
                           'C1_hit_dphi_closest_L2', 1000, 0., 50.)

hit_dphi_closest_L4 = TH1D('hit_dphi_closest_L4', 'hit_dphi_closest_L4', 1000, 0., 50.)
C1_hit_dphi_closest_L4 = TH1D('C1_hit_dphi_closest_L4',
                           'C1_hit_dphi_closest_L4', 1000, 0., 50.)

hit_dphi_closest_L6 = TH1D('hit_dphi_closest_L6', 'hit_dphi_closest_L6', 1000, 0., 50.)
C1_hit_dphi_closest_L6 = TH1D('C1_hit_dphi_closest_L6',
                           'C1_hit_dphi_closest_L6', 1000, 0., 50.)

hit_mindR = TH1D(
    'hit_mindR', 'hit_mindR', 100, 0., 0.1)
C1_hit_mindR = TH1D(
    'C1_hit_mindR', 'C1_hit_mindR', 100, 0., 0.1)

hit_dtheta_vs_dphi = TH2D(
    'hit_dtheta_vs_dphi', 'hit_dtheta_vs_dphi', 1500, 0., 30., 1000, 0., 50.)
C1_hit_dtheta_vs_dphi = TH2D(
    'C1_hit_dtheta_vs_dphi', 'C1_hit_dtheta_vs_dphi', 1500, 0., 30., 1000, 0., 50.)

h_C1_z = TH1D('h_C1_z', 'h_C1_z', 100, -10., 10.)

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

for ievt, event in enumerate(reader):

    if ievt%100==0:
        print("Event " + str(ievt))

    # setting decoder
    vertexHitsCollection = event.getCollection('VXDTrackerHits')
    encoding = vertexHitsCollection.getParameters(
    ).getStringVal(EVENT.LCIO.CellIDEncoding)
    decoder = UTIL.BitField64(encoding)

    relationCollection=event.getCollection('VXDTrackerHitRelations')
    relation = UTIL.LCRelationNavigator(relationCollection)

    # loop once over all the hits and split them in layers for the delta theta check
    hits_by_layer = []
    for ilay in range(0, 8):
        Li_hits = []
        for imod in range(0, 32):
            mi_hits = []
            Li_hits.append(mi_hits)
        hits_by_layer.append(Li_hits)
        
    for hit in vertexHitsCollection:
        cellID = int(hit.getCellID0())
        decoder.setValue(cellID)
        layer = decoder['layer'].value()
        module = decoder['module'].value()
        position = hit.getPosition()  # mm
        hits_by_layer[layer][module].append(hit)

    # loop by layer
    for ilay in reversed(range(0, 8)):

        for imod in range(0, 32):
            the_layer = hits_by_layer[ilay][imod]

            # looping over hit collections
            for ihit, hit in enumerate(the_layer):

                if ihit%1000 == 0 and ihit>0:
                    print("Looking at hit " + str(ihit))

                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder['layer'].value()

                #particles propagate outward
                if (layer == 1) or (layer == 3) or (layer == 5) or (layer == 7):
                    continue

                position = hit.getPosition()  # mm
                r = sqrt(position[0]*position[0] + position[1]*position[1])
                z = position[2]
                tvec = TVector3(position[0], position[1], position[2])
                hit_theta = tvec.Theta()

                min_dR = 998.
                dtheta_close = -1000
                dphi_close = -1000

                the_other_layer = hits_by_layer[1][imod]
                if layer == 2:
                    the_other_layer = hits_by_layer[3][imod]
                if layer == 4:
                    the_other_layer = hits_by_layer[5][imod]
                if layer == 6:
                    the_other_layer = hits_by_layer[7][imod]

                for hit2 in the_other_layer:
                    position2 = hit2.getPosition()  # mm
                    tvec2 = TVector3(position2[0], position2[1], position2[2])
                    hit_theta2 = tvec2.Theta()

                    dtheta = fabs(hit_theta2-hit_theta)
                    dphi = tvec.DeltaPhi(tvec2)

                    dR = sqrt(dphi*dphi+dtheta*dtheta)

                    if dR < min_dR:
                        min_dR = dR
                        dtheta_close = dtheta
                        dphi_close = dphi

                has_rel = len(relation.getRelatedToObjects(hit))
                if has_rel > 0:

                    simhit = relation.getRelatedToObjects(hit)[0]

                    try:
                        mcp = simhit.getMCParticle()
                        pdg = mcp.getPDG()
                    except:
                        mcp = 0
                        pdg = 0
                else:
                    mcp = 0
                    pdg = 0

                if len(the_other_layer)>0:

                    dtheta_close = fabs(dtheta_close*1000.)
                    dphi_close = fabs(dphi_close*1000.)
    
                    if fabs(pdg) == 1000024:  # if chargino
                    
                        z_of_C1 = mcp.getVertex()[2]
                        h_C1_z.Fill(z_of_C1)
    
                        #look only within 3 sigmas from z==0 (sim has sigma_z = 1.5*mm)
                        #if fabs(z_of_C1)<4.5:
                        C1_hit_dtheta_closest.Fill(dtheta_close)
                        C1_hit_dphi_closest.Fill(dphi_close)
                        C1_hit_mindR.Fill(min_dR)
                        C1_hit_dtheta_vs_dphi.Fill(dtheta_close, dphi_close)
                        if layer == 0:
                            C1_hit_dtheta_closest_L0.Fill(dtheta_close)
                            C1_hit_dphi_closest_L0.Fill(dphi_close)
                        if layer == 2:
                            C1_hit_dtheta_closest_L2.Fill(dtheta_close)
                            C1_hit_dphi_closest_L2.Fill(dphi_close)
                        if layer == 4:
                            C1_hit_dtheta_closest_L4.Fill(dtheta_close)
                            C1_hit_dphi_closest_L4.Fill(dphi_close)
                        if layer == 6:
                            C1_hit_dtheta_closest_L6.Fill(dtheta_close)
                            C1_hit_dphi_closest_L6.Fill(dphi_close)
    
                    hit_dtheta_closest.Fill(dtheta_close)
                    hit_dphi_closest.Fill(dphi_close)
                    hit_mindR.Fill(min_dR)
                    hit_dtheta_vs_dphi.Fill(dtheta_close, dphi_close)
                    if layer == 0:
                        hit_dtheta_closest_L0.Fill(dtheta_close)
                        hit_dphi_closest_L0.Fill(dphi_close)
                    if layer == 2:
                        hit_dtheta_closest_L2.Fill(dtheta_close)
                        hit_dphi_closest_L2.Fill(dphi_close)
                    if layer == 4:
                        hit_dtheta_closest_L4.Fill(dtheta_close)
                        hit_dphi_closest_L4.Fill(dphi_close)
                    if layer == 6:
                        hit_dtheta_closest_L6.Fill(dtheta_close)
                        hit_dphi_closest_L6.Fill(dphi_close)

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')

hit_mindR.Write()
C1_hit_mindR.Write()
hit_dtheta_vs_dphi.Write()
C1_hit_dtheta_vs_dphi.Write()

hit_dtheta_closest.Write()
C1_hit_dtheta_closest.Write()
hit_dphi_closest.Write()
C1_hit_dphi_closest.Write()
hit_dtheta_closest_L0.Write()
C1_hit_dtheta_closest_L0.Write()
hit_dphi_closest_L0.Write()
C1_hit_dphi_closest_L0.Write()
hit_dtheta_closest_L2.Write()
C1_hit_dtheta_closest_L2.Write()
hit_dphi_closest_L2.Write()
C1_hit_dphi_closest_L2.Write()
hit_dtheta_closest_L4.Write()
C1_hit_dtheta_closest_L4.Write()
hit_dphi_closest_L4.Write()
C1_hit_dphi_closest_L4.Write()
hit_dtheta_closest_L6.Write()
C1_hit_dtheta_closest_L6.Write()
hit_dphi_closest_L6.Write()
C1_hit_dphi_closest_L6.Write()

h_C1_z.Write()

output_file.Close()
