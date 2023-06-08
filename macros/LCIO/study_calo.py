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
parser.add_option('-o', '--outFile', help='--outFile histos_calo.root',
                  type=str, default='histos_calo.root')
(options, args) = parser.parse_args()

ecal_hit_E = TH1D('ECAL_hitEnergy', 'ECAL_hitEnergy', 100, 0., 0.1)
hcal_hit_E = TH1D('HCAL_hitEnergy', 'HCAL_hitEnergy', 100, 0., 0.5)
ecal_hit_r_barrel = TH1D('ECAL_barrel', 'ECAL_barrel', 100, 1590, 1850)
ecal_hit_rz = TH2D('ECAL_rz', 'ECAL_rz', 1000, -2800, 2800, 1000, 0, 2000)
ecal_hit_depth_barrel = TH1D(
    'ECAL_barrel_depth', 'ECAL_barrel_depth', 100, 0, 300)
ecal_hit_r_endcap = TH1D('ECAL_endcap', 'ECAL_endcap', 100, 310, 1700)
hcal_hit_r_barrel = TH1D('HCAL_barrel', 'HCAL_barrel', 100, 1852, 3852)
hcal_hit_r_endcap = TH1D('HCAL_endcap', 'HCAL_endcap', 100, 307, 3246)
hcal_hit_depth_barrel = TH1D(
    'HCAL_barrel_depth', 'HCAL_barrel_depth', 100, 0, 2000)

# Low-level sim hit distributions
h_ECAL_simhit_E = TH1D('ECAL_simhit_E', 'ECAL_simhit_E', 100, 0, 20)  # GeV
h_ECAL_simhit_layer = TH1D(
    'ECAL_simhit_layer', 'ECAL_simhit_layer', 100, 0, 100)
h_ECAL_simhit_layer_ele = TH1D(
    'ECAL_simhit_layer_ele', 'ECAL_simhit_layer_ele', 100, 0, 100)
h_ECAL_simhit_layer_gamma = TH1D(
    'ECAL_simhit_layer_gamma', 'ECAL_simhit_layer_gamma', 100, 0, 100)
h_ECAL_simhit_layer_other = TH1D(
    'ECAL_simhit_layer_other', 'ECAL_simhit_layer_other', 100, 0, 100)

h_HCAL_simhit_E = TH1D('HCAL_simhit_E', 'HCAL_simhit_E', 100, 0, 20)  # GeV
h_HCAL_simhit_layer = TH1D(
    'HCAL_simhit_layer', 'HCAL_simhit_layer', 100, 0, 100)

# Low-level digitised hit distributions
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV

# Histo list for writing to outputs
histos_list = [h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
               h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
               h_ECAL_simhit_E, h_HCAL_simhit_E,
               h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac, h_EMfrac_PFO,
               h_ECAL_hit_layer, h_HCAL_hit_layer,
               h_ECAL_simhit_layer, h_ECAL_simhit_layer_ele, h_ECAL_simhit_layer_gamma, h_ECAL_simhit_layer_other,
               h_HCAL_simhit_layer,
               ecal_hit_E, hcal_hit_E, 
               ecal_hit_r_barrel, ecal_hit_depth_barrel, ecal_hit_rz, ecal_hit_r_endcap,
               hcal_hit_r_barrel, hcal_hit_r_endcap, hcal_hit_depth_barrel,
               ]

for histo in histos_list:
    histo.SetDirectory(0)


#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)
# loop over all events in the file

doTruth = True
r_min_ecal = 1590
r_min_hcal = 1852

for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Event " + str(ievt))

    # ECAL barrel
    hitsCollection = event.getCollection("ECALBarrel")

    if doTruth:
        relationCollection = event.getCollection('RelationCaloHit')
        relation = UTIL.LCRelationNavigator(relationCollection)
        for ihit, hit in enumerate(hitsCollection):
            sim_vec = relation.getRelatedToObjects(hit)
            print(len(sim_vec))
            '''
            for simhit in sim_vec:
                mystring = ""
                for ipart in range(0, simhit.getNMCParticles()):
                    mcpart = simhit.getPDGCont(ipart)
                    if mcpart != 0:
                        mystring = mystring + " " + str(mcpart)
                print(mystring)
            '''

    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        ecal_hit_E.Fill(hit.getEnergy())
        # [1/cm^2], I hope
        density_weight = hit.getEnergy()*100./(2.*pi*r*2210.*2.)
        ecal_hit_depth_barrel.Fill(r-r_min_ecal, density_weight)
        ecal_hit_r_barrel.Fill(r, density_weight)
        ecal_hit_rz.Fill(hit.getPosition()[2], r)

    # ECAL endcap
    hitsCollection = event.getCollection("ECALEndcap")
    for ihit, hit in enumerate(hitsCollection):
        ecal_hit_E.Fill(hit.getEnergy())
        ecal_hit_r_endcap.Fill(sqrt(hit.getPosition()[
            0]*hit.getPosition()[0] + hit.getPosition()[1]*hit.getPosition()[1]))

    # HCAL barrel
    hitsCollection = event.getCollection("HCALBarrel")
    for ihit, hit in enumerate(hitsCollection):
        r = sqrt(hit.getPosition()[0]*hit.getPosition()
                 [0] + hit.getPosition()[1]*hit.getPosition()[1])
        hcal_hit_E.Fill(hit.getEnergy())
        # [1/cm^2], I hope
        density_weight = hit.getEnergy()*100./(2.*pi*r*2569.*2.)
        hcal_hit_r_barrel.Fill(r, density_weight)
        hcal_hit_depth_barrel.Fill(r-r_min_hcal, density_weight)

    # HCAL endcap
    hitsCollection = event.getCollection("HCALEndcap")
    for ihit, hit in enumerate(hitsCollection):
        hcal_hit_E.Fill(hit.getEnergy())
        hcal_hit_r_endcap.Fill(sqrt(hit.getPosition()[
            0]*hit.getPosition()[0] + hit.getPosition()[1]*hit.getPosition()[1]))

reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
