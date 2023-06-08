from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch
import itertools

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_truth.root',
                  type=str, default='ntup_truth.root')
(options, args) = parser.parse_args()

#########################

min_dphi = TH1D('min_dphi', 'min_dphi', 100, 0., TMath.Pi())  # GeV

tree = TTree("truth_tree", "truth_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt_vec = array('d', [0])
pz_vec = array('d', [0])
E_vec = array('d', [0])
m_vec = array('d', [0])
charge_vec = array('i', [0])
beta_vec = array('d', [0])
gamma_vec = array('d', [0])
phi_vec = array('d', [0])
theta_vec = array('d', [0])
z_truth_vec = array('d', [0])
r_truth_vec = array('d', [0])
zend_truth_vec = array('d', [0])
rend_truth_vec = array('d', [0])
pdgID_vec = array('i', [0])
status_vec = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt_vec,  'var/D')
tree.Branch("pz",  pz_vec,  'var/D')
tree.Branch("E",  E_vec,  'var/D')
tree.Branch("m",  m_vec,  'var/D')
tree.Branch("charge",  charge_vec,  'var/I')
tree.Branch("beta",  beta_vec,  'var/D')
tree.Branch("gamma",  gamma_vec,  'var/D')
tree.Branch("phi", phi_vec, 'var/D')
tree.Branch("theta", theta_vec, 'var/D')
tree.Branch("z_truth", z_truth_vec, 'var/D')
tree.Branch("r_truth", r_truth_vec, 'var/D')
tree.Branch("zend_truth", zend_truth_vec, 'var/D')
tree.Branch("rend_truth", rend_truth_vec, 'var/D')
tree.Branch("pdgID", pdgID_vec, 'var/I')
tree.Branch("status", status_vec, 'var/I')

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

count = 0

# loop over all events in the file
for ievt, event in enumerate(reader):

    # if ievt > 0:
    #    print("Interrupting at " + str(ievt))
    #    break

    # grab mc collection
    mcpCollection = event.getCollection('MCParticle')

    if ievt % 100 == 0:
        print("Processing event " + str(ievt))
        print(str(len(mcpCollection)) + " truth particles")

    interesting_particles_list = []

    for c in mcpCollection:
        if c.getGeneratorStatus() == 1:
            pdg = c.getPDG()
            if fabs(pdg) == 22 and c.getEnergy() > 20.:
                interesting_particles_list.append(c)

    n_ph100 = 0
    for c in interesting_particles_list:
        vx = c.getVertex()
        end = c.getEndpoint()
        dp3 = c.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], c.getEnergy())

        if c.getEnergy() > 100 and fabs(tlv.Eta()) < 2.5:
            n_ph100 = n_ph100 + 1

        pt_vec[0] = tlv.Perp()
        pz_vec[0] = dp3[2]
        E_vec[0] = c.getEnergy()
        m_vec[0] = c.getMass()
        charge_vec[0] = int(c.getCharge())
        beta_vec[0] = tlv.Beta()
        gamma_vec[0] = tlv.Gamma()
        phi_vec[0] = tlv.Phi()
        theta_vec[0] = tlv.Theta()
        z_truth_vec[0] = vx[2]
        r_truth_vec[0] = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
        zend_truth_vec[0] = end[2]
        rend_truth_vec[0] = sqrt(end[0]*end[0]+end[1]*end[1])
        pdgID_vec[0] = c.getPDG()
        status_vec[0] = c.getGeneratorStatus()

        tree.Fill()

    mindphi = 1000.
    for pair in itertools.combinations(interesting_particles_list, 2):
        dp3A = pair[0].getMomentum()
        tlvA = TLorentzVector()
        tlvA.SetPxPyPzE(dp3A[0], dp3A[1], dp3A[2], pair[0].getEnergy())

        dp3B = pair[1].getMomentum()
        tlvB = TLorentzVector()
        tlvB.SetPxPyPzE(dp3B[0], dp3B[1], dp3B[2], pair[1].getEnergy())

        newdphi = tlvA.DeltaPhi(tlvB)
        if newdphi < mindphi:
            mindphi = newdphi

    #print(str(len(interesting_particles_list)) + " " + str(mindphi))
    min_dphi.Fill(mindphi)

    if n_ph100 == 2 and fabs(mindphi) > 2:
        count = count+1

reader.close()

print("Good events: " + str(count))

# write outputs
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
min_dphi.Write()
output_file.Close()
