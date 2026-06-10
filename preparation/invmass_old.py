
from ROOT import TFile, TH1F, TDirectory, TCanvas, TDirectory, TF1, TPad, gStyle, TLine
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

NORM_LOW_MASS = 3.8 # 
NORM_HIGH_MASS = 3.9 #

def normalise_histograms_and_compute_subtraction(outdir:TDirectory, mode:str, rebin:int=1, suffix:str=''):

    h_same = file_same.Get(f'invmass/hInvariantMass{mode}{suffix}')
    h_same.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_same.SetName(f'hSameEvent{suffix}')
    if rebin > 1:
        h_same.Rebin(rebin)

    h_mixed = file_mixed.Get(f'invmass/hInvariantMass{mode}{suffix}')
    h_mixed.SetDirectory(0)  # Detach from file to avoid issues with deletion
    h_mixed.SetName(f'hMixedEvent{suffix}')
    if rebin > 1:
        h_mixed.Rebin(rebin)


    low_bin = h_same.FindBin(NORM_LOW_MASS)
    high_bin = h_same.FindBin(NORM_HIGH_MASS)
    normalization_factor = h_same.Integral(low_bin, high_bin) / h_mixed.Integral(low_bin, high_bin)
    h_mixed.Scale(normalization_factor)

    h_subtracted = h_same.Clone(f'hInvariantMassSubtracted{suffix}')
    h_subtracted.Add(h_mixed, -1.)

    if suffix != '':
        return h_same, h_mixed, h_subtracted
        
    outdir.cd()
    for hist in [h_same, h_mixed, h_subtracted]:
        hist.SetTitle(';#it{m}_{p - ^{3}He} (GeV/#it{c}^{2}); Counts (a.u.)')
        hist.Write()

    return h_same, h_mixed, h_subtracted

def invariant_mass_subtraction_corrected(h_same, h_mixed, h_correction, outdir:TDirectory, mode:str, suffix:str=''):
    """
    Compute the correlation function for the centrality integrated histograms.
    """
    
    h_mixed_corrected = h_mixed.Clone(f'hMixedEventCorrected{suffix}')
    for ibin in range(1, h_mixed_corrected.GetNbinsX()+1):
        
        invmass_value = h_mixed.GetBinCenter(ibin)
        mixed_value = h_mixed.GetBinContent(ibin)
        correction_value = h_correction.GetBinContent(h_correction.FindBin(invmass_value))

        h_mixed_corrected.SetBinContent(ibin, mixed_value*correction_value)

    h_subtracted_corrected = h_same.Clone(f'hInvariantMassSubtractedCorrected{suffix}')
    h_subtracted_corrected.Add(h_mixed_corrected, -1.)

    outdir.cd()
    for hist in [h_mixed_corrected, h_subtracted_corrected]:
        hist.SetTitle(';#it{m}_{p - ^{3}He} (GeV/#it{c}^{2}); Counts (a.u.)')
        hist.Write()
    h_correction.Write(f'hCorrection{suffix}')

    return h_mixed_corrected, h_subtracted_corrected

def nsigma_plot(h_subtracted, outdir:TDirectory, mode:str, suffix:str=''):

    h_nsigma = h_subtracted.Clone(f'hNsigma{suffix}')
    for ibin in range(1, h_nsigma.GetNbinsX()+1):
        value = h_subtracted.GetBinContent(ibin)
        error = h_subtracted.GetBinError(ibin)
        h_nsigma.SetBinContent(ibin, value/error if error > 0. else 0.)
        h_nsigma.SetBinError(ibin, 0)
    
    outdir.cd()
    h_nsigma.SetTitle(';#it{m}_{p - ^{3}He} (GeV/#it{c}^{2}); n#sigma')
    h_nsigma.Write()

    return h_nsigma

def plot_invmass_over_nsigma(h_invmass, h_nsigma, pdf_path:str, x_limits:list):

    gStyle.SetOptStat(0)
    y_portion = 0.3
    canvas = TCanvas('cInvMassOverNsigma', '', 800, 800)
    
    canvas.cd()
    upper_pad = TPad('upper_pad', '', 0, y_portion, 1, 1.)
    upper_pad.Draw()

    h_invmass.GetXaxis().SetRangeUser(*x_limits)
    set_root_object(h_invmass, marker_color=418, marker_style=20, line_color=418)
    upper_pad.cd()
    h_invmass.Draw('e0')

    upper_pad.Update()
    line_upper = TLine(3.751, upper_pad.GetUymin(), 3.751, upper_pad.GetUymax())
    set_root_object(line_upper, line_color=2, line_style=2)
    line_upper.Draw('same')

    canvas.cd()
    lower_pad = TPad('lower_pad', '', 0, 0.05, 1, y_portion)
    lower_pad.SetBottomMargin(0.22)
    lower_pad.Draw()

    h_nsigma.GetXaxis().SetRangeUser(*x_limits)
    set_root_object(h_nsigma, marker_color=797, marker_style=20, x_title_size=0.1, y_title_size=0.1,
                    x_title_offset=1, y_title_offset=0.3, x_label_size=0.1, y_label_size=0.1)
    lower_pad.cd()
    h_nsigma.Draw('p0')

    lower_pad.Update()
    line_lower = TLine(3.751, lower_pad.GetUymin(), 3.751, lower_pad.GetUymax())
    set_root_object(line_lower, line_color=2, line_style=2)
    line_lower.Draw('same')

    canvas.SaveAs(pdf_path)
    del canvas
        

if __name__ == '__main__':
    
    #infile_same = 'output/same_event.root'
    #infile_mixed = 'output/mixed_event.root'
    #_outfile = 'output/correlation.root'

    infile_same = 'checks/same_event_hadronpid_pass1_pass4_nohe3pcut.root'
    infile_mixed = 'checks/rotation_mixing_hadronpid_pass1_pass4_nohe3pcut.root'
    _outfile = 'checks/invmass_hadronpid_pass1_pass4_nohe3pcut.root'

    file_same = TFile.Open(infile_same)
    file_mixed = TFile.Open(infile_mixed)
    outfile = TFile.Open(_outfile, 'RECREATE')   

    for mode in ['Matter', 'Antimatter']:
    
        outdir = outfile.mkdir(f'InvariantMass{mode}')

        h_sames, h_mixeds = [], []

        for suffix in ['', 'PtHadronUnder1p0', 'PtHadronUnder1p1', 'PtHadronUnder1p2',
                       'PLiUnder2', 'PLiUnder3', 'PLiUnder4', 'PLiUnder5',
                       'PLi2', 'PLi3', 'PLi4', 'PLi5']:
            
            pdf_path = f'checks/invmass_hadronpid_pass1_pass4_nclstpc{mode}{suffix}.pdf'
            h_same, h_mixed, h_subtracted = normalise_histograms_and_compute_subtraction(outdir, mode, rebin=4, suffix=suffix)
            h_correction = load_hist('/home/galucia/Lithium4/femto/models/CATS_converted.root', 'hHe3_p_Coul_InvMass')
            __, h_subtracted_corrected = invariant_mass_subtraction_corrected(h_same, h_mixed, h_correction, outdir, mode, suffix)
            h_nsigma = nsigma_plot(h_subtracted_corrected, outdir, mode, suffix)
            if suffix == '':
                plot_invmass_over_nsigma(h_subtracted_corrected, h_nsigma, pdf_path, [3.747, 3.777])

    outfile.Close()
