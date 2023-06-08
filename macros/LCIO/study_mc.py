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

Bfield = 3.57  # T
output_file_name = 'histos_mc'

#########################


def find_C1_in_decay_chain(list_particles):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if fabs(pdg) == 1000005:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found chargino with no daughters')
                continue
            else:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        result.append(mcp)
    return result


#########################
# declare histograms
photon_pt = TH1D('photon_pt', 'photon_pt', 100, 0, 5000)  # GeV
photon_pt.SetDirectory(0)
stop_beta = TH1D('stop_beta', 'stop_beta', 100, 0, 1)  # GeV
stop_beta.SetDirectory(0)
stop_theta = TH1D('stop_theta', 'stop_theta', 100, 0., 180.)  # GeV
stop_theta.SetDirectory(0)
stop_gamma = TH1D('stop_gamma', 'stop_gamma', 60, 0, 6)  # GeV
stop_gamma.SetDirectory(0)
stop_gammabeta = TH1D('stop_gammabeta', 'stop_gammabeta', 120, 0, 6)  # GeV
stop_gammabeta.SetDirectory(0)
bquark_pt = TH1D('bquark_pt', 'bquark_pt', 100, 0, 5)  # GeV
bquark_pt.SetDirectory(0)
mc_pt = TH1D('mc_pt', 'mc_pt', 100, 0, 5000)  # GeV
mc_pt.SetDirectory(0)
mc_pz = TH1D('mc_pz', 'mc_pz', 200, -10000, 10000)  # GeV
mc_pz.SetDirectory(0)
mc_E = TH1D('mc_E', 'mc_E', 200, -10000, 10000)  # GeV
mc_E.SetDirectory(0)
mc_m = TH1D('mc_m', 'mc_m', 200, -10000, 10000)  # GeV
mc_m.SetDirectory(0)
mc_charge = TH1D('mc_charge', 'mc_charge', 4, -2, 2)  # GeV
mc_charge.SetDirectory(0)
mc_genStatus = TH1D('mc_genStatus', 'mc_genStatus', 40, 0,
                    40)  # it should be in the 21-23 range
mc_genStatus.SetDirectory(0)
mc_endPoint = TH1D('mc_endPoint', 'mc_endPoint', 100, 0., 1500.)  # mm
mc_endPoint.SetDirectory(0)

B_endPoint = TH1D('B_endPoint', 'B_endPoint', 40, 0., 20.)  # mm
B_endPoint.SetDirectory(0)

mc_radialEndPoint = TH1D(
    'mc_radialEndPoint', 'mc_radialEndPoint', 100, 0., 1500.)  # mm
mc_radialEndPoint.SetDirectory(0)
mc_pdg = TH1D('mc_pdg', 'mc_pdg', 2, -2., 2.)
mc_pdg.SetDirectory(0)
nC1perEvent = TH1D('nC1perEvent', 'nC1perEvent', 4, 0, 4)
nC1perEvent.SetDirectory(0)
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
        for part in mcpCollection:
            if fabs(part.getPDG()) == 22:
                p = part.getMomentum()  # GeV
                pt = sqrt(p[0]*p[0] +
                          p[1]*p[1])
                photon_pt.Fill(pt)
            if fabs(part.getPDG()) == 5:
                p = part.getMomentum()  # GeV
                pt = sqrt(p[0]*p[0] +
                          p[1]*p[1])
                bquark_pt.Fill(pt)

        result = find_C1_in_decay_chain(mcpCollection)
        nC1 = len(result)

        good_pid = [1000005, 1005321, 1000522, 1000512]

        for c in mcpCollection:
            #        for c in result:
            if fabs(c.getPDG()) in good_pid:
                end = c.getEndpoint()
                endpoint = sqrt(
                    end[0]*end[0] + end[1]*end[1] + end[2]*end[2])
                pos = c.getVertex()
                print(str(c.getPDG()) + " at " + str(pos[0]) + "  " +
                      str(pos[1]) + "  " + str(pos[2]))
                print("decay at " + str(end[0]) + "  " +
                      str(end[1]) + "  " + str(end[2]))

                daughters = c.getDaughters()
                dau_list = []

                for d in daughters:
                    dau_list.append(d.getPDG())
                    if fabs(d.getPDG()) == 5:
                        # found the b, now we have to figure out where it's going
                        b_daughter_pgdid = 5
                        the_part = d
                        n_step = 0
                        while b_daughter_pgdid == 5:
                            #print("Step "+str(n_step))
                            b_daughters = the_part.getDaughters()
                            for k in b_daughters:
                                b_daughter_pgdid = fabs(k.getPDG())
                                end = k.getEndpoint()
                                endpoint = sqrt(
                                    end[0]*end[0] + end[1]*end[1] + end[2]*end[2])
                                momentum = k.getMomentum()  # GeV
                                pt = sqrt(
                                    momentum[0]*momentum[0] + momentum[1]*momentum[1])
                                mass = k.getMass()
                                # print(" " + str(k.getPDG()) +
                                #      " " + str(endpoint) + " " + str(pt) + " " + str(mass))
                                if fabs(b_daughter_pgdid) > 500 and fabs(b_daughter_pgdid) < 600:
                                    B_endPoint.Fill(endpoint)
                            the_part = k
                            n_step = n_step + 1

                print("daughters " + str(dau_list))

                endpoint = c.getEndpoint()
                endpoint_l = sqrt(
                    endpoint[0]*endpoint[0] + endpoint[1]*endpoint[1] + endpoint[2]*endpoint[2])
                if endpoint_l > 0:
                    endpoint_r = sqrt(
                        endpoint[0]*endpoint[0] + endpoint[1]*endpoint[1])
                    momentum = c.getMomentum()  # GeV
                    pt = sqrt(momentum[0]*momentum[0] +
                              momentum[1]*momentum[1])
                    pz = momentum[2]
                    momentum = c.getMomentum()
                    tlv = TLorentzVector()
                    tlv.SetPxPyPzE(momentum[0], momentum[1],
                                   momentum[2], c.getEnergy())
                    stop_beta.Fill(tlv.Beta())
                    stop_gamma.Fill(tlv.Gamma())
                    stop_gammabeta.Fill(tlv.Gamma()*tlv.Beta())
                    theta = tlv.Theta()*180./TMath.Pi()
                    stop_theta.Fill(theta)
                    mc_pt.Fill(pt)
                    mc_pz.Fill(pz)
                    mc_E.Fill(c.getEnergy())
                    mc_m.Fill(c.getMass())
                    mc_charge.Fill(c.getCharge())
                    mc_genStatus.Fill(c.getGeneratorStatus())
                    mc_radialEndPoint.Fill(endpoint_r)  # mm
                    mc_endPoint.Fill(endpoint_l)  # mm
                    c_pdg = c.getPDG()
                    if c_pdg > 0:
                        mc_pdg.Fill(1)
                    else:
                        mc_pdg.Fill(-1)
        nC1perEvent.Fill(nC1)
    except:
        print("Exception for event " + str(ievt) +
              ": no MC particle collection found!")
        mcpCollection = 0
reader.close()

# write histograms
output_file_name = output_file_name + '.root'
output_file = TFile(output_file_name, 'RECREATE')
mc_pt.Write()
mc_pz.Write()
mc_E.Write()
mc_m.Write()
mc_charge.Write()
mc_genStatus.Write()
mc_endPoint.Write()
mc_radialEndPoint.Write()
mc_pdg.Write()
nC1perEvent.Write()
photon_pt.Write()
bquark_pt.Write()
stop_beta.Write()
stop_gamma.Write()
stop_gammabeta.Write()
stop_theta.Write()
B_endPoint.Write()

output_file.Close()
