from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
from math import *
from optparse import OptionParser
import json

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_tracks',
                  type=str, default='ntup_tracks')
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


def isGoodTrack(track):

    pt = round(0.3 * Bfield / fabs(track.getOmega() * 1000.), 2)
    if pt < 1 : return False

    chi2 = track.getChi2()
    ndf = track.getNdf()
    if chi2/ndf > 10: return False

    hits = track.getTrackerHits()
    numhits = len(hits)
    #numholes = int(track.getdEdxError())  # BADHACK

    if numhits < 6: return False

    return True  

#########################
# declare histograms

h_truth_rxy   = TH1D('truth_rxy'  , 'truth_rxy' , 100, -10, 10)  # prod rxy, mm
h_truth_pt    = TH1D('truth_pt'   , 'truth_pt'    , 100,   0, 1000) # GeV
h_truth_theta = TH1D('truth_theta', 'truth_theta' ,  50,  0.4, 2.7) # approx 30-150 deg
h_truth_phi   = TH1D('truth_phi'  , 'truth_phi'   ,  50, -3.5, 3.5)

h_track_d0      = TH1D('track_d0'     , 'track_d0'     , 100,   -5,  5)
h_track_z0      = TH1D('track_z0'     , 'track_z0'     , 100,  -20, 20)
h_track_pt      = TH1D('track_pt'     , 'track_pt'     , 100,    0, 1000) 
h_track_phi     = TH1D('track_phi'    , 'track_phi'    ,  50, -3.5, 3.5)
h_track_theta   = TH1D('track_theta'  , 'track_theta'  ,  50,  0.4, 2.7)
h_track_nholes  = TH1D('track_nholes' , 'track_nholes' ,  20,    0, 20)
h_track_nhits   = TH1D('track_nhits'  , 'track_nhits'  ,  20,    0, 20)
h_track_chi2ndf = TH1D('track_chi2ndf', 'track_chi2ndf', 100,    0, 100)

h_ntrack_all    = TH1D('ntrack_all'   , 'ntrack_all'  , 100,  0, 5000)  # mm
h_ntrack_good   = TH1D('ntrack_good'  , 'ntrack_good' , 100,  0, 5000)  # mm

h_track_hit_time  = TH1D('track_hit_time', 'track_hit_time', 100, -1, 1)

histos_list = [h_truth_rxy, h_truth_pt, h_truth_theta, h_truth_phi, 
               h_track_d0, h_track_z0, h_track_pt, h_track_phi, h_track_theta, 
               h_track_nholes, h_track_nhits, h_track_chi2ndf, 
               h_ntrack_all, h_ntrack_good, h_track_hit_time]

for histo in histos_list:
    histo.SetDirectory(0)

#########################

# if you want to save data for later analysis 
# create lists of lists where the first index is the event, the 2nd is the object 
# we'll output these to a json
# then you should be able read them in with 
# awkward_array = ak.Array(python_lists)
truth_pt    = [] 
truth_phi   = [] 
truth_theta = [] 

track_pt    = [] 
track_phi   = [] 
track_theta = [] 
track_d0    = []  
track_z0    = []  
track_d0err = []  
track_z0err = []  
track_chi2  = []  
track_ndof  = []  
track_nhits = [] 

track_hit_t = []
track_hit_x = []
track_hit_y = []
track_hit_z = []

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))


    # more info here: 
    # https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1MCParticle.html
    truthParticles = event.getCollection('MCParticle')

    # initialize truth muon info for this event
    tru_pt    = [] 
    tru_phi   = [] 
    tru_theta = [] 
    for itruth,truthParticle in enumerate(truthParticles): 

        if abs(truthParticle.getPDG()) != 13 : continue # select muons 
        if truthParticle.getGeneratorStatus()!=1 : continue # select status=1

        # get any information you need
        p3 = truthParticle.getMomentum() 
        p4 = TLorentzVector()
        p4.SetPxPyPzE(p3[0], p3[1], p3[2], truthParticle.getEnergy())

        vertex = truthParticle.getVertex()
        rxy = (vertex[0]**2 + vertex[1]**2)**0.5

        charge   = truthParticle.getCharge() 
        
        if ievent==0: 
            print("TruthMuon,",itruth)
            print("  (pt,theta,phi),",p4.Pt(), p4.Theta(), p4.Phi() ) 

        # fill histograms
        h_truth_rxy  .Fill( rxy )
        h_truth_pt   .Fill( p4.Pt() )  
        h_truth_phi  .Fill( p4.Phi() ) 
        h_truth_theta.Fill( p4.Theta() ) 

        # append to this event's list
        tru_pt   .append( p4.Pt()    ) 
        tru_theta.append( p4.Theta() ) 
        tru_phi  .append( p4.Phi()   ) 


    hit_relations = []
    IBTrackerHitsRelations = event.getCollection('IBTrackerHitsRelations')
    hit_relations.append(IBTrackerHitsRelations)
    IETrackerHitsRelations = event.getCollection('IETrackerHitsRelations')
    hit_relations.append(IETrackerHitsRelations)
    OBTrackerHitsRelations = event.getCollection('OBTrackerHitsRelations')
    hit_relations.append(OBTrackerHitsRelations)
    OETrackerHitsRelations = event.getCollection('OETrackerHitsRelations')
    hit_relations.append(OETrackerHitsRelations)
    VBTrackerHitsRelations = event.getCollection('VBTrackerHitsRelations')
    hit_relations.append(VBTrackerHitsRelations)
    VETrackerHitsRelations = event.getCollection('VETrackerHitsRelations')
    hit_relations.append(VETrackerHitsRelations)

    hit_collections = []
    IBTrackerHits = event.getCollection('IBTrackerHits')
    hit_collections.append(IBTrackerHits)
    IETrackerHits = event.getCollection('IETrackerHits')
    hit_collections.append(IETrackerHits)
    OBTrackerHits = event.getCollection('OBTrackerHits')
    hit_collections.append(OBTrackerHits)
    OETrackerHits = event.getCollection('OETrackerHits')
    hit_collections.append(OETrackerHits)
    VBTrackerHits = event.getCollection('VBTrackerHits')
    hit_collections.append(VBTrackerHits)
    VETrackerHits = event.getCollection('VETrackerHits')
    hit_collections.append(VETrackerHits)

    # filling standard tracks
    tracks = event.getCollection('SiTracks')
    nTracksGood = 0

    # initialize this event's list 
    trk_pt   = [] 
    trk_phi  = []
    trk_theta = []
    trk_d0    = []
    trk_z0    = []
    trk_d0err = []
    trk_z0err = []
    trk_chi2  = []
    trk_ndof  = []
    trk_nhits = []
    trk_hit_t = []
    trk_hit_x = []
    trk_hit_y = []
    trk_hit_z = []
    for itrack, track in enumerate(tracks):
    
        if not isGoodTrack(track) : continue # apply some very basic requirements
        nTracksGood+=1
    
        # get any info you need
        pt    = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi   = track.getPhi()
        theta = TMath.Pi()/2-atan(track.getTanLambda())

        d0 = track.getD0()
        z0 = track.getZ0()
        d0err = track.getCovMatrix()[0]
        z0err = track.getCovMatrix()[9]
        chi2 = track.getChi2()
        ndof = track.getNdf()
    
        hits = track.getTrackerHits()
        nhits = len(hits)

        if ievent == 0: 
            print("Track, ",itrack)
            print("  (pt,theta,phi),",pt,phi,theta)

        # now decode hits
        encoding = hit_collections[0].getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        hit_x = []
        hit_y = []
        hit_z = []
        hit_t = []
        for ihit, hit in enumerate(hits): 

            position = hit.getPosition()
            x = position[0]
            y = position[1]
            z = position[2]
            rxy = (x*x+y*y)**0.5
            dedx = hit.getdEdx()
            time = hit.getTime()
            if ievent == 0 : 
                print("  (ihit,rxy,z,time), {} {:.2f} {:.2f} {:.3f}".format(ihit, rxy, z, time) ) 

            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            detector = decoder["system"].value()
            layer = decoder["layer"].value() 

            if ievent == 0 : 
                print("     (detector,layer),",detector,layer)


            h_track_hit_time.Fill(time)

            hit_x.append(x)
            hit_y.append(y)
            hit_z.append(z)
            hit_t.append(time)


        # fill track histograms
        h_track_pt   .Fill(pt)
        h_track_phi  .Fill(phi)
        h_track_theta.Fill(theta)
        h_track_d0   .Fill(d0)
        h_track_z0   .Fill(z0)
        h_track_nhits  .Fill(nhits)
        h_track_chi2ndf.Fill(chi2/float(ndof))

        # save track info to this event's lists 
        trk_pt   .append(pt) 
        trk_phi  .append(phi) 
        trk_theta.append(theta) 
        trk_d0.append(d0) 
        trk_z0.append(z0) 
        trk_d0err.append(d0err)
        trk_z0err.append(z0err)
        trk_chi2 .append(chi2) 
        trk_ndof .append(ndof) 
        trk_nhits.append(nhits)

        trk_hit_x.append(hit_x)
        trk_hit_y.append(hit_y)
        trk_hit_z.append(hit_z)
        trk_hit_t.append(hit_t)

    # fill any event level histograms
    h_ntrack_all .Fill( len(tracks) ) 
    h_ntrack_good.Fill( nTracksGood ) 

    # append the event's lists to final output 
    truth_pt   .append(tru_pt)
    truth_phi  .append(tru_phi)
    truth_theta.append(tru_theta)

    track_pt   .append(trk_pt)
    track_phi  .append(trk_phi)
    track_theta.append(trk_theta) 

    track_d0   .append(trk_d0)
    track_z0   .append(trk_z0)
    track_d0err.append(trk_d0err)
    track_z0err.append(trk_z0err)
    track_chi2 .append(trk_chi2)
    track_ndof .append(trk_ndof)
    track_nhits.append(trk_nhits)

    track_hit_x.append(trk_hit_x)
    track_hit_y.append(trk_hit_y)
    track_hit_z.append(trk_hit_z)
    track_hit_t.append(trk_hit_t)


reader.close()

# write histograms

output_root = TFile(options.outFile+".root","RECREATE")

for histo in histos_list:
    histo.Write()

output_root.Close()

# write output
output = { "truth_pt"   : truth_pt,
           "truth_phi"  : truth_phi,
           "truth_theta": truth_theta,

           "track_pt"   : track_pt,
           "track_phi"  : track_phi,
           "track_theta": track_theta,
           "track_d0"   : track_d0,
           "track_z0"   : track_z0,
           "track_d0err" : track_d0err,
           "track_z0err" : track_z0err,
           "track_chi2" : track_chi2,
           "track_ndof" : track_ndof,
           "track_nhits" : track_nhits,

           "track_hit_t"    : track_hit_t,
           "track_hit_x"    : track_hit_x,
           "track_hit_y"    : track_hit_y,
           "track_hit_z"    : track_hit_z,
        }

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
            json.dump(output, fp)
