import ROOT
from optparse import OptionParser


def normalise_hist(hist):
    if hist.Integral() > 0:
        hist.Scale(1./hist.Integral())


def make_plots(title, sig_std, sig_llp, bg_std, bg_llp):
    c = ROOT.TCanvas(" ", " ", 800, 600)
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.15)

    max = 0
    if sig_std.GetMaximum() > max:
        max = sig_std.GetMaximum()
    sig_std.SetLineWidth(2)

    if sig_llp.GetMaximum() > max:
        max = sig_llp.GetMaximum()
    sig_llp.SetLineWidth(2)
    sig_llp.SetLineColor(ROOT.kRed)

    if bg_std.GetMaximum() > max:
        max = bg_std.GetMaximum()
    bg_std.SetLineWidth(2)
    bg_std.SetLineStyle(2)

    if bg_llp.GetMaximum() > max:
        max = bg_llp.GetMaximum()
    bg_llp.SetLineWidth(2)
    bg_llp.SetLineStyle(2)
    bg_llp.SetLineColor(ROOT.kRed)

    sig_std.SetMaximum(1.05*max)
    sig_std.DrawCopy("HIST")
    sig_llp.DrawCopy("HISTsame")
    bg_std.DrawCopy("HISTsame")
    bg_llp.DrawCopy("HISTsame")

    # Save the canvas
    c.SaveAs(title+".pdf")


parser = OptionParser()
parser.add_option('-s', '--inSig', help='--inSig ../../MuCData/',
                  type=str, default='../../MuCData/ntup_smuons.root')
parser.add_option('-b', '--inBg', help='--inBg ../../MuCData/',
                  type=str, default='../../MuCData/ntup_bg.root')
parser.add_option('-t', '--inTree', help='--inTree tracks_tree',
                  type=str, default='tracks_tree')
(options, args) = parser.parse_args()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTickX(1)

# DataFrames for signal tracks
d_sig = ROOT.RDataFrame(options.inTree, options.inSig).Define(
    'chi2ndf', 'chi2 / ndf')
d_sig_std = d_sig.Filter('isLLP==0', '!isLLP').Filter('r_truth', 'r_truth<127')
d_sig_llp = d_sig.Filter('isLLP==1', 'isLLP').Filter('r_truth', 'r_truth<127')

filt_sig_std = d_sig_std.Filter('nhits > 2', 'nhits > 2')
filt_sig_llp = d_sig_llp.Filter('nhits > 5', 'nhits > 5')

# DataFrames for background tracks
d_bg = ROOT.RDataFrame(options.inTree, options.inBg).Define(
    'chi2ndf', 'chi2 / ndf')
d_bg_std = d_bg.Filter('isLLP==0', '!isLLP')
d_bg_llp = d_bg.Filter('isLLP==1', 'isLLP')

filt_bg_std = d_bg_std.Filter('nhits > 2', 'nhits > 2')
filt_bg_llp = d_bg_llp.Filter('nhits > 5', 'nhits > 5')

# Book histograms
hist_list = []

h_d0_sig_std = filt_sig_std.Histo1D(("d0", "", 60, -30, 30), "d0")
hist_list.append(h_d0_sig_std)
h_d0_sig_llp = filt_sig_llp.Histo1D(("d0", "", 60, -30, 30), "d0")
hist_list.append(h_d0_sig_llp)
h_d0_bg_std = filt_bg_std.Histo1D(("d0", "", 60, -30, 30), "d0")
hist_list.append(h_d0_bg_std)
h_d0_bg_llp = filt_bg_llp.Histo1D(("d0", "", 60, -30, 30), "d0")
hist_list.append(h_d0_bg_llp)

h_z0_sig_std = filt_sig_std.Histo1D(("z0", "", 60, -30, 30), "z0")
hist_list.append(h_z0_sig_std)
h_z0_sig_llp = filt_sig_llp.Histo1D(("z0", "", 60, -30, 30), "z0")
hist_list.append(h_z0_sig_llp)
h_z0_bg_std = filt_bg_std.Histo1D(("z0", "", 60, -30, 30), "z0")
hist_list.append(h_z0_bg_std)
h_z0_bg_llp = filt_bg_llp.Histo1D(("z0", "", 60, -30, 30), "z0")
hist_list.append(h_z0_bg_llp)

h_nhits_sig_std = filt_sig_std.Histo1D(("nhits", "", 20, 0, 20), "nhits")
hist_list.append(h_nhits_sig_std)
h_nhits_sig_llp = filt_sig_llp.Histo1D(("nhits", "", 20, 0, 20), "nhits")
hist_list.append(h_nhits_sig_llp)
h_nhits_bg_std = filt_bg_std.Histo1D(("nhits", "", 20, 0, 20), "nhits")
hist_list.append(h_nhits_bg_std)
h_nhits_bg_llp = filt_bg_llp.Histo1D(("nhits", "", 20, 0, 20), "nhits")
hist_list.append(h_nhits_bg_llp)

h_chi2ndf_sig_std = filt_sig_std.Histo1D(("chi2ndf", "", 20, 0, 20), "chi2ndf")
hist_list.append(h_chi2ndf_sig_std)
h_chi2ndf_sig_llp = filt_sig_llp.Histo1D(("chi2ndf", "", 20, 0, 20), "chi2ndf")
hist_list.append(h_chi2ndf_sig_llp)
h_chi2ndf_bg_std = filt_bg_std.Histo1D(("chi2ndf", "", 20, 0, 20), "chi2ndf")
hist_list.append(h_chi2ndf_bg_std)
h_chi2ndf_bg_llp = filt_bg_llp.Histo1D(("chi2ndf", "", 20, 0, 20), "chi2ndf")
hist_list.append(h_chi2ndf_bg_llp)

h_rtruth_sig_std = filt_sig_std.Histo1D(("r_truth", "", 20, 0, 20), "r_truth")
hist_list.append(h_rtruth_sig_std)
h_rtruth_sig_llp = filt_sig_llp.Histo1D(("r_truth", "", 20, 0, 20), "r_truth")
hist_list.append(h_rtruth_sig_llp)
h_rtruth_bg_std = filt_bg_std.Histo1D(("r_truth", "", 20, 0, 20), "r_truth")
hist_list.append(h_rtruth_bg_std)
h_rtruth_bg_llp = filt_bg_llp.Histo1D(("r_truth", "", 20, 0, 20), "r_truth")
hist_list.append(h_rtruth_bg_llp)

# Normalise everything
for hist in hist_list:
    normalise_hist(hist)

# Plot the histograms side by side on a canvas
make_plots("d0", h_d0_sig_std, h_d0_sig_llp, h_d0_bg_std, h_d0_bg_llp)
make_plots("z0", h_z0_sig_std, h_z0_sig_llp, h_z0_bg_std, h_z0_bg_llp)
make_plots("nhits", h_nhits_sig_std, h_nhits_sig_llp,
           h_nhits_bg_std, h_nhits_bg_llp)
make_plots("chi2ndf", h_chi2ndf_sig_std, h_chi2ndf_sig_llp,
           h_chi2ndf_bg_std, h_chi2ndf_bg_llp)
make_plots("rtruth", h_rtruth_sig_std, h_rtruth_sig_llp,
           h_rtruth_bg_std, h_rtruth_bg_llp)

#augmented1 = filtered2.Define('b3', 'b1 / b2')
#filtered3 = augmented1.Filter('b3 < .5','Cut3')

print('Signal !LLP stats')
report_sig_std = filt_sig_std.Report()
report_sig_std.Print()

print('Signal LLP stats')
report_sig_llp = filt_sig_llp.Report()
report_sig_llp.Print()

print('Background !LLP stats')
report_bg_std = filt_bg_std.Report()
report_bg_std.Print()

print('Background LLP stats')
report_bg_llp = filt_bg_llp.Report()
report_bg_llp.Print()
