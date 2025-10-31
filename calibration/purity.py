import yaml
import numpy as np
import pandas as pd
from ROOT import TFile, TGraphErrors, TCanvas, TPaveText, TF1, TH2F, RDataFrame, TChain
from ROOT import RooRealVar, RooCrystalBall, RooGenericPdf, RooAddPdf, RooDataHist, RooArgList, RooChebychev, RooGaussian

from torchic.utils.root import set_root_object
from torchic.core.graph import create_graph

import sys
sys.path.append('..')
from preparation.prepare import prepare_selections, prepare_input_tchain
from utils.plot_utils import write_params_to_text, PdfHandler
from utils.roopdf_utils import init_roopdf


def prepare_rdataframe(chain_data: TChain, base_selection: str, selection: str):
      
      # Recalibration
      #.Redefine('fNSigmaTPCHe3', 'ComputeNsigmaTPCHe(std::abs(fInnerParamTPCHe3), fSignalTPCHe3)') \
    
    rdf = RDataFrame(chain_data) 
    if 'fNSigmaTPCHadPr' in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTPCHad', 'fNSigmaTPCHadPr')
    if 'fNSigmaTOFHadPr' not in rdf.GetColumnNames():
        rdf = rdf.Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(std;;abs(fPtHad), fMassTOFHad)')
    
    rdf = rdf.Define('fSignedPtHad', 'fPtHad') \
      .Define('fSignHe3', 'fPtHe3/std::abs(fPtHe3)') \
      .Redefine('fPtHe3', 'std::abs(fPtHe3)') \
      .Redefine('fPtHad', 'std::abs(fPtHad)') \
      .Redefine('fPtHe3', '(fPIDtrkHe3 == 7) || (fPIDtrkHe3 == 8) || (fPtHe3 > 2.5) ? fPtHe3 : CorrectPidTrkHe(fPtHe3)') \
      .Define('fSignedPtHe3', 'fPtHe3 * fSignHe3') \
      .Redefine('fInnerParamTPCHe3', 'fInnerParamTPCHe3 * 2') \
      .Define('fNSigmaTOFHad', 'ComputeNsigmaTOFPr(fPtHad, fMassTOFHad)') \
      .Define('fClusterSizeCosLamHe3', 'ComputeAverageClusterSize(fItsClusterSizeHe3) / cosh(fEtaHe3)') \
      .Define('fClusterSizeCosLamHad', 'ComputeAverageClusterSize(fItsClusterSizeHad) / cosh(fEtaHad)') \
      .Define('fNSigmaITSHe3', 'ComputeNsigmaITSHe(fPtHe3 * std::cosh(fEtaHe3), fClusterSizeCosLamHe3)') \
      .Define('fNSigmaITSHad', 'ComputeNsigmaITSPr(fPtHad * std::cosh(fEtaHad), fClusterSizeCosLamHad)') \
      .Define('fNSigmaDCAxyHe3', 'ComputeNsigmaDCAxyHe(fPtHe3, fDCAxyHe3)') \
      .Define('fNSigmaDCAzHe3', 'ComputeNsigmaDCAzHe(fPtHe3, fDCAzHe3)') \
      .Define('fNSigmaDCAxyHad', 'ComputeNsigmaDCAxyPr(fPtHad, fDCAxyHad)') \
      .Define('fNSigmaDCAzHad', 'ComputeNsigmaDCAzPr(fPtHad, fDCAzHad)') \
      .Filter(base_selection).Filter(selection) 
    
    return rdf

def build_signal_model(x):
    # signal function
    signal_pars = {
        'mean': RooRealVar('mean', '#mu', 0., -2, 2),
        'sigma': RooRealVar('sigma', '#sigma', 1, 0.3, 2),
        'aL': RooRealVar('aL', '#alpha', 1., 0.7, 10.),
        'nL': RooRealVar('nL', 'n', 1., 0.1, 30.)
    }
    signal = RooCrystalBall('signal', 'signal', x, 
                            signal_pars['mean'], signal_pars['sigma'], 
                            signal_pars['aL'], signal_pars['nL'], doubleSided=True)
    return signal, signal_pars

def build_background_model(x, function:str='pol1', suffix:str=''):

    if function == 'exp':
        bkg_exp_pars = {
            'alpha': RooRealVar('alpha', 'alpha', 5., 0., 1.e5),
            'offset': RooRealVar('offset', 'offset', -5., -1.e5, 0.)
        }
        bkg_exp = RooGenericPdf(f'bkg{suffix}', 'bkg', f'exp(- alpha * ({x.GetName()} - offset))', 
                                [x, *bkg_exp_pars.values()])
        return bkg_exp, bkg_exp_pars

    elif function == 'gaus':

        bkg_gaus_pars = {
            'mean': RooRealVar('mean_gaus', 'mean_gaus', 0., -6, -2),
            'sigma': RooRealVar('sigma_gaus', 'sigma_gaus', 0.5, 1e-3, 1e3)
        }
        bkg_gaus = RooGaussian(f'bkg{suffix}', 'bkg', x, *bkg_gaus_pars.values())
        return bkg_gaus, bkg_gaus_pars

    elif function == 'pol1':
        bkg_pol1_pars = {
            'p0': RooRealVar('p0', 'c_{0}', 100, 0., 1e5),
            'p1': RooRealVar('p1', 'c_{1}', 0.1, -1.e5, 1.e5),
        }
        bkg_pol1 = RooGenericPdf(f'bkg{suffix}', 'bkg', f'p0 + p1 * {x.GetName()}', [x, *bkg_pol1_pars.values()])
        return bkg_pol1, bkg_pol1_pars
    
    elif function == 'pol0':
        bkg_pol0_pars = {
        }
        bkg_pol1 = RooGenericPdf(f'bkg{suffix}', 'bkg', f'1', [x, *bkg_pol0_pars.values()])
        return bkg_pol1, bkg_pol0_pars

    else:
        raise ValueError('Accepted values for the function parameter are "exp", "gaus", "pol1", "pol0')

def fit_sidebands_root(h_data, bkg, bkg_pars, fit_range, sig_window, x, outdir):
    """
    Perform a two-step fit:
      1. Fit sidebands with TF1 to get normalization
      2. Transfer parameters to RooFit background model
    
    Parameters
    ----------
    h_data : TH1
        Data histogram
    invmass : RooRealVar
        Mass variable
    bkg : RooAbsPdf
        Background model
    bkg_pars : dict
        Dictionary of background parameters (RooRealVar)
    sig_window : tuple
        Signal window (min, max) to exclude from sideband fit
    """
    
    bkg_name = bkg.GetName()
    if 'exp' in bkg_name:
        tf1_bkg = TF1(f'tf1_{bkg_name}', '[0]*exp(-[2] * (x - [1]))', fit_range[0], fit_range[1])
        tf1_bkg.SetParameters(100, -5, 1)  # Initial guesses
        tf1_bkg.SetParLimits(0, 0, 1e6)
    elif 'pol1' in bkg_name:
        tf1_bkg = TF1(f'tf1_{bkg_name}', '[0] + [1]*x', fit_range[0], fit_range[1])
        tf1_bkg.SetParameters(100, 0)
        tf1_bkg.SetParLimits(0, 0, 1e6)
    elif 'pol0' in bkg_name:
        tf1_bkg = TF1(f'tf1_{bkg_name}', '[0]', fit_range[0], fit_range[1])
        tf1_bkg.SetParameters(0.1)
        tf1_bkg.SetParLimits(0, 0, 1e6)
    elif 'gaus' in bkg_name:
        tf1_bkg = TF1(f'tf1_{bkg_name}', 'gaus', fit_range[0], fit_range[1])
        tf1_bkg.SetParameters(100, 0, 0)
        tf1_bkg.SetParLimits(0, 0, 1e6)
    else:
        raise ValueError(f'Unsupported background function for TF1 fit: {bkg_name}')
    
    h_sidebands = h_data.Clone(f'{h_data.GetName()}_sidebands')
    x_low_max = h_sidebands.GetXaxis().FindBin(sig_window[0])
    x_high_min = h_sidebands.GetXaxis().FindBin(sig_window[1])
    for ibin in range(x_low_max, x_high_min + 1):
        h_sidebands.SetBinContent(ibin, 0.)
        h_sidebands.SetBinError(ibin, 0.)
    
    h_sidebands.Fit(tf1_bkg, 'RQN')  # R=range, Q=quiet, N=no draw
    
    if 'exp' in bkg_name:
        bkg_pars[f'offset'].setVal(tf1_bkg.GetParameter(1))
        bkg_pars[f'alpha'].setVal(tf1_bkg.GetParameter(2))
    elif 'pol1' in bkg_name:
        bkg_pars['p0'].setVal(tf1_bkg.GetParameter(0))
        bkg_pars['p1'].setVal(tf1_bkg.GetParameter(1))
    elif 'gaus' in bkg_name:
        bkg_pars['mean'].setVal(tf1_bkg.GetParameter(1))
        bkg_pars['sigma'].setVal(tf1_bkg.GetParameter(2))
    
    for par in bkg_pars.values():
        par.setConstant(True)
    
    N_bkg_full = tf1_bkg.Integral(fit_range[0], fit_range[1]) / h_data.GetBinWidth(1)
    
    bkg_normalisation = RooRealVar('bkg_normalisation', '#it{N}_{bkg}', N_bkg_full)
    bkg_normalisation.setConstant(True)
    
    canvas = TCanvas(f'sidebands_{bkg_name}_fit_{x:.2f}', '')
    canvas.SetLogy(True)
    h_sidebands.SetMinimum(0.1)
    h_sidebands.Draw('E')
    tf1_bkg.Draw('same')
    
    text = TPaveText(0.7, 0.5, .85, 0.55, 'ndc')
    text.SetFillColor(0)
    text.SetBorderSize(0)
    chi2_ndf = tf1_bkg.GetChisquare()/tf1_bkg.GetNDF() if tf1_bkg.GetNDF() > 0 else 999.
    text.AddText(f'#chi^{{2}} / ndf = {chi2_ndf:.2f}')
    text.AddText(f'#it{{N}}_{{bkg}} = {N_bkg_full:.2f}')
    text.Draw('same')
    
    outdir.cd()
    canvas.Write()
    
    return bkg_normalisation, chi2_ndf

def _fit_purity_slice(signal_pdf, signal_pars, hist, x, bkg_funcs,
                      pt_range, nsigma_integral_range, outdir_bkg_fits):

    pt = 0.5 * (pt_range[1] + pt_range[0])

    signal_window = (-4, 4)
    fit_range = (-6, 6)
    bkgs, bkgs_pars, bkgs_normalisation, bkgs_chi2 = [], [], [], []
    
    for bkg_func in bkg_funcs:

        bkg, bkg_pars = build_background_model(x, function=bkg_func, suffix=f'_{bkg_func}')
        bkg_normalisation, chi2 = fit_sidebands_root(hist, bkg, bkg_pars, fit_range, signal_window, pt, outdir_bkg_fits)
        
        if chi2 < 0:
            continue
        bkgs.append(bkg)
        bkgs_pars.append(bkg_pars)
        bkgs_chi2.append(chi2)
        bkgs_normalisation.append(bkg_normalisation)

    print('check after bkg studies')

    index_min = min(range(len(bkgs_chi2)), key=bkgs_chi2.__getitem__)
    bkg_pdf = bkgs[index_min]
    bkg_pars = bkgs_pars[index_min]
    bkg_normalisation = bkgs_normalisation[index_min] if len(bkgs_normalisation) > 0 else RooRealVar('bkg_normalisation', 'bkg_normalisation', 1., 0., 1e6)

    signal_normalisation = RooRealVar('signal_normalisation', 'signal_normalisation', 1., 0., 1e6)
    model = RooAddPdf('model', 'signal + bkg', [signal_pdf, bkg_pdf], [signal_normalisation, bkg_normalisation])
    #signal_normalisation = RooRealVar('signal_normalisation', 'signal_normalisation', 0.5, 0., 1.)
    #model = RooAddPdf('model', 'signal + bkg', [signal_pdf, bkg_pdf], [signal_normalisation])

    dh = RooDataHist(f'dh_{pt:.2f}', 'dh', [x], Import=hist)
    model.fitTo(dh, PrintLevel=1, Save=True)

    frame = x.frame(Title=f'{pt_range[0]:.2f} < #it{{p}}_{{T}} < {pt_range[1]:.2f} GeV/#it{{c}}')
    dh.plotOn(frame)
    model.plotOn(frame, LineColor=2)
    model.plotOn(frame, Components={signal_pdf}, LineColor=3, LineStyle='--')
    model.plotOn(frame, Components={bkg_pdf}, LineColor=4, LineStyle='--')
    #for ipdf, bkg_pdf in enumerate(bkg_pdfs, start=4):
    #    model.plotOn(frame, Components={bkg_pdf}, LineColor=ipdf, LineStyle='--')

    x.setRange('integral_range', nsigma_integral_range[0], nsigma_integral_range[1])
    integral_signal = signal_pdf.createIntegral(x, x, 'integral_range').getVal() * signal_normalisation.getVal()
    integral_bkg = bkg_pdf.createIntegral(x, x, 'integral_range').getVal() * bkg_normalisation.getVal()

    ifit_results = {'pt': np.abs(pt),
                    'int_sig': integral_signal, 'sig_mean': signal_pars['mean'].getVal(), 
                    'sig_sigma': signal_pars['sigma'].getVal(), 
                    'sig_mean_err': signal_pars['sigma'].getVal()/np.sqrt(integral_signal),
                    'sig_sigma_err': signal_pars['sigma'].getError(),
                    'purity': integral_signal / (integral_signal + integral_bkg)}

    return frame, ifit_results

def draw_results(h2_nsigma: TH2F, fit_results: pd.DataFrame, pt_min, pt_max, pdf_canvas: PdfHandler, out_dir: TFile, sign: str):
    '''
        Draw the results of the fit
    '''
    fit_results['pt_err'] = h2_nsigma.GetXaxis().GetBinWidth(1)/2.
    graph_signal = create_graph(fit_results, 'pt', 'int_sig', 'pt_err', 0., name=f'g_signal_{sign}', title=';#it{p}_{T} (GeV/c); int (n#sigma_{TPC})')

    graph_signal_mean = create_graph(fit_results, 'pt', 'sig_mean', 'pt_err', 0., name=f'g_signal_mean_{sign}', title=';#it{p}_{T} (GeV/c); #mu(n#sigma_{TPC})')
    fit_mean = TF1('fit_mean', '[0]', pt_min, pt_max)
    graph_signal_mean.Fit(fit_mean, 'RSM+')
    graph_signal_mean.SetMarkerStyle(20)
    pdf_canvas.draw_save_and_clear(graph_signal_mean, draw_option='ap', logy=False)

    graph_signal_sigma = create_graph(fit_results, 'pt', 'sig_sigma', 'pt_err', 0., name=f'g_signal_sigma_{sign}', title=';#it{p}_{T} (GeV/c); #sigma(n#sigma_{TPC})')
    fit_sigma = TF1('fit_sigma', 'pol0', pt_min, pt_max)
    graph_signal_sigma.Fit(fit_sigma, 'RSM+')
    graph_signal_sigma.SetMarkerStyle(20)
    pdf_canvas.draw_save_and_clear(graph_signal_sigma, draw_option='ap', logy=False)

    graph_purity = create_graph(fit_results, 'pt', 'purity', 'pt_err', 0., name=f'g_purity_{sign}', title=';#it{p}_{T} (GeV/c); purity')

    out_dir.cd()
    graph_signal.Write(f'g_signal_{sign}')
    graph_signal_mean.Write(f'g_signal_mean_{sign}')
    graph_signal_sigma.Write(f'g_signal_sigma_{sign}')
    graph_purity.Write(f'g_purity_{sign}')

def fit_purity(h2_nsigma: TH2F, output_file: TFile, particle:str, detector:str='TPC'):
    
    tmp_pt_min, tmp_pt_max = (1.6, 3.5) if particle == 'He3' else (0.4, 4)
    nsigma_min, nsigma_max = -6, 6
    nsigma_low_int, nsigma_high_int = -2, 2
    
    nsigma = RooRealVar('nsigma', f'#it{{n}}#sigma_{{{detector}}}', nsigma_min, nsigma_max)

    signal, signal_pars = build_signal_model(nsigma)
    bkg_funcs = [#'pol1',
        'pol0', 'exp', 'gaus'] if (particle == 'Had' and detector == 'TPC') else ['pol0']

    detector_dir = output_file.mkdir(detector)
    out_dir = detector_dir.mkdir('fits')
    outdir_bkg_fits = detector_dir.mkdir('bkg_fits')

    for sign in ['matter', 'antimatter']:

        output_pdf_path = output_file.GetFile().GetName().split('.')[0] + f'_{particle}_{detector}_{sign}.pdf'
        with PdfHandler(output_pdf_path) as pdf_canvas:
        
            fit_results = None

            if sign == 'antimatter':
                pt_min, pt_max = -tmp_pt_max, -tmp_pt_min
            else:
                pt_min, pt_max = tmp_pt_min, tmp_pt_max

            for pt_bin in range(h2_nsigma.GetXaxis().FindBin(pt_min), h2_nsigma.GetXaxis().FindBin(pt_max)):

                pt = h2_nsigma.GetXaxis().GetBinCenter(pt_bin)
                h_nsigma = h2_nsigma.ProjectionY(f'h_nsigma{detector}_{pt:.2f}', pt_bin, pt_bin, 'e')
                #h_nsigma.Rebin(4)
                if h_nsigma.GetEntries() < 100:
                    continue

                pt_low_edge = h2_nsigma.GetXaxis().GetBinLowEdge(pt_bin)
                pt_high_edge = h2_nsigma.GetXaxis().GetBinLowEdge(pt_bin+1)
                nsigma_frame, ifit_results = _fit_purity_slice(signal, signal_pars, h_nsigma, nsigma, bkg_funcs,
                                                               (pt_low_edge, pt_high_edge), (nsigma_low_int, nsigma_high_int),
                                                               outdir_bkg_fits)
                
                if fit_results is None:
                    fit_results = pd.DataFrame.from_dict([ifit_results])
                else:
                    fit_results = pd.concat([fit_results, pd.DataFrame.from_dict([ifit_results])], ignore_index=True)

                text = write_params_to_text({**signal_pars#, **bkg_pol1_pars
                                             }.values(), coordinates=[0.4, 0.15, 0.63, 0.4])
                #text.AddText(f'#chi^{{2}} / NDF = {nsigma_frame.chiSquare():.2f}')
                nsigma_frame.addObject(text)
                nsigma_frame.SetMinimum(1)
                pdf_canvas.draw_and_save(nsigma_frame, logy=True)

                out_dir.cd()
                pdf_canvas.canvas.Write(f'nsigma_frame_{sign}_{pt:.2f}')
                pdf_canvas.clear()

            draw_results(h2_nsigma, fit_results, pt_min, pt_max, pdf_canvas, detector_dir, sign)
                
    h2_nsigma.Write()

if __name__ == '__main__':
    
    output_file = TFile('output/purity.root', 'recreate')

    config_file = '/home/galucia/Lithium4/preparation/config/config_prepare_purity.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection = prepare_selections(config)
    chain_data = prepare_input_tchain(config)
    rdf = prepare_rdataframe(chain_data, base_selection, selection)

    for particle in ['Had', 'He3']:
        
        outdir = output_file.mkdir(particle)
        
        x_variable, y_variable_tpc, y_variable_tof = ('fSignedPtHe3', 'fNSigmaTPCHe3', '') \
                    if particle == 'He3' else ('fSignedPtHad', 'fNSigmaTPCHad', 'fNSigmaTOFHad')
        
        h2_nsigma_tpc = rdf.Histo2D(('h2NsigmaTpcVsPt', ';#it{p}_{T} (GeV/#it{c}); #it{n}#sigma_{TPC};',
                                     60, -6, 6, 60, -6, 6), 
                                     x_variable, y_variable_tpc)
        fit_purity(h2_nsigma_tpc, outdir, particle=particle, detector='TPC')
        
        if particle == 'Had':
            h2_nsigma_tof = rdf.Filter('std::abs(fNSigmaTPCHad)< 2').Histo2D(('h2NsigmaTofVsPt', 
                                         ';#it{p}_{T} (GeV/#it{c}); #it{n}#sigma_{TOF};',
                                         60, -6, 6, 60, -6, 6), x_variable, y_variable_tof)
            fit_purity(h2_nsigma_tof, outdir, particle=particle, detector='TOF')
    
    output_file.Close()
