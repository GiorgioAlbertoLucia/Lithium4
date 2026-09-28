import numpy as np

from ROOT import (
    TFile, TDirectory, TMarker, TGraphErrors,
    RooWorkspace, TGraph, TH2D, TCanvas, RooMsgService, TGraph
)
import ROOT
from uncertainties import ufloat

from core.signal_fitter import SignalFitter
from core.bkg_fitter import BkgFitter
from core.model_fitter import ModelFitter

from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.root import set_alice_global_style, set_root_object, init_legend

import faulthandler
faulthandler.enable()

MEASUREMENTS = {
    'Matter': {
        '010': {'R': ufloat(6.12, 0.11), 'A': ufloat(0.157, 0.077)},
        '1050': {'R': ufloat(4.95, 0.05), 'A': ufloat(0.490, 0.073)}
    },
    'Antimatter': {
        '010': {'R': ufloat(6.12, 0.11), 'A': ufloat(0.278, 0.101)},
        '1050': {'R': ufloat(4.95, 0.05), 'A': ufloat(0.434, 0.089)}
    },
    'Both': {
        '010': {'R': ufloat(6.12, 0.11), 'A': ufloat(0.205, 0.061)},
        '1050': {'R': ufloat(4.95, 0.05), 'A': ufloat(0.465, 0.056)}
    }
}
RADIUS_VALUES = np.arange(3.5, 8.0 + 0.05, 0.1)
SIGNAL_NORMALISATIONS = np.arange(0.0, 1.0 + 0.005, 0.005)

CENTRALITIES_TO_SCAN = ['010', '1050']

MODES_TO_SCAN = ['', 'Matter', 'Antimatter']

INPUT_MODELS_PATH = 'models/lambda_models_scan_LL.root'
DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/correlation_PbPb_hadronpid.root'
SYSTEMATICS_FILE_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/systematic_uncertainties.root'
SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution_proper_sill.root', 'hCkHist')

SELECTION = 'Default'

def make_chi2_lines(delta_chi2, chi2_map: TH2D, chi2_min: float) -> TGraph:
    target = chi2_min + delta_chi2

    g_lower = TGraph()
    g_upper = TGraph()

    ilower = 0
    iupper = 0

    for ix in range(1, chi2_map.GetNbinsX() + 1):

        radius = chi2_map.GetXaxis().GetBinCenter(ix)

        crossings = []

        for iy in range(1, chi2_map.GetNbinsY()):

            chi2_1 = chi2_map.GetBinContent(ix, iy)
            chi2_2 = chi2_map.GetBinContent(ix, iy + 1)

            if chi2_1 <= 0 or chi2_2 <= 0:
                continue

            if (chi2_1 - target) * (chi2_2 - target) <= 0:
                y1 = chi2_map.GetYaxis().GetBinCenter(iy)
                y2 = chi2_map.GetYaxis().GetBinCenter(iy + 1)
                if chi2_2 != chi2_1:
                    y_cross = y1 + (target - chi2_1) * (y2 - y1) / (chi2_2 - chi2_1)
                else:
                    y_cross = 0.5 * (y1 + y2)
                crossings.append(y_cross)

        crossings.sort()

        if len(crossings) >= 2:
            g_lower.SetPoint(ilower, radius, crossings[0]); ilower += 1
            g_upper.SetPoint(iupper, radius, crossings[-1]); iupper += 1
        elif len(crossings) == 1:
            # tangent point: both edges meet here, so the tip closes
            g_lower.SetPoint(ilower, radius, crossings[0]); ilower += 1
            g_upper.SetPoint(iupper, radius, crossings[0]); iupper += 1

    return g_lower, g_upper

def get_data_histograms(mode: str, centrality: str):

    sign = 'Both' if mode == '' else mode

    if centrality in ['1050', '080', '1080']:
        h_data_name = f'Correlation{mode}/{SELECTION}/hCorrelationDirectComputation{centrality}'
        h_mixed_event_name = f'Correlation{mode}/{SELECTION}/hNormalisedMixedEventDirectComputation{centrality}'
        h_same_event_name = f'Correlation{mode}/{SELECTION}/hSameEventDirectComputation{centrality}'

    else:
        h_data_name = f'Correlation{mode}/{SELECTION}/hCorrelation{centrality}'
        h_mixed_event_name = f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent{centrality}'
        h_same_event_name = f'Correlation{mode}/{SELECTION}/hSameEvent{centrality}'

    h_systematics_name = f'Correlation{mode}/{centrality}/hCorrelationSyst{centrality}'

    h_data = load_hist(DATA_INPUT_PATH, h_data_name)
    h_mixed_event = load_hist(DATA_INPUT_PATH, h_mixed_event_name)
    h_same_event = load_hist(DATA_INPUT_PATH, h_same_event_name)
    h_systematics = load_hist(SYSTEMATICS_FILE_PATH, h_systematics_name)

    return h_data, h_same_event, h_mixed_event, h_systematics

def radius_scan(outfile: TDirectory, mode: str, centrality: str):

    sign = 'Both' if mode == '' else mode
    print(tc.GREEN + f'\nScanning radius for mode={mode}, centrality={centrality}' + tc.RESET)

    h_correlation_function, h_same_event, h_mixed_event, h_systematics = get_data_histograms(mode, centrality)

    is_first_bin_empty = h_correlation_function.GetBinContent(1) < 1e-12
    KSTAR_MIN, KSTAR_MAX = (0.02, 0.4) if is_first_bin_empty else (0.01, 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)

    scan_dir = outfile.mkdir(f'{sign}_{centrality}')

    chi2_graph = TGraph()
    chi2_graph.SetName('gChi2VsRadius')
    chi2_graph.SetTitle(f'{sign} {centrality};#it{{R}}_{{s}} (fm);#chi^{{2}}')

    for ir, radius in enumerate(RADIUS_VALUES):

        radius_str = f'{radius:.2f}'
        scan_dir_radius = scan_dir.mkdir(f'R_{radius_str}')
        print(tc.CYAN + f'  R = {radius_str} fm' + tc.RESET)

        h_bkg_name = f'{sign}/{centrality}/R_{radius_str}/hLambdaSigmaCorrectedCk_Smeared_R_{radius_str}'
        h_bkg = load_hist(INPUT_MODELS_PATH, h_bkg_name)

        workspace = RooWorkspace(f'roows_R_{radius_str}')

        signal_fitter = SignalFitter('signal', kstar_spec, scan_dir_radius, workspace)
        signal_fitter.init_signal('from_kde', h_signal, rho=3)
        signal_fitter.title = '^{4}Li'
        signal_fitter.save_to_workspace()
        
        bkg_fitter = BkgFitter('bkg', kstar_spec, scan_dir_radius, workspace)
        bkg_fitter.init_bkg('from_kde', h_bkg, rho=0.1)
        bkg_fitter.title = 'Coulomb + strong interaction'
        bkg_fitter.save_to_workspace()

        model_fitter = ModelFitter('model', kstar_spec, scan_dir_radius,
                                   ['signal_pdf'], ['bkg_pdf'], workspace,
                                   extended=True, title='^{4}Li + interaction')

        model_fitter.fractions['signal_pdf'].setRange(0., 1.)
        model_fitter.fractions['signal_pdf'].setVal(0.3)
        model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')

        model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')

        model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
        model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4),
            range_name='bkg_fit_range', save_normalisation_value=True)

        sign_label = ('p#minus^{3}He' if mode == 'Matter' 
                    else '#bar{p}#minus^{3}#bar{He}' if mode == 'Antimatter'
                    else 'p#minus^{3}He #oplus #bar{p}#minus^{3}#bar{He}')

        model_fitter.fit_model(h_correlation_function, 
                               signal_name='signal_pdf', norm_range='bkg_fit_range', data_label=sign_label)
        model_fitter.save_to_workspace()
        
        h_available_bkgs_lambda = [load_hist(INPUT_MODELS_PATH, f'{sign}/{centrality}/hLambdaSigmaCorrectedCk_Smeared{variation_str}') 
                                   for variation_str in ['', '_higher', '_lower']]
        
        chi2 = model_fitter.compute_chi2_new(h_correlation_function, h_systematics, h_bkg, 
                                            [h_available_bkgs_lambda], suffix=f'_{sign}_{centrality}_R_{radius_str}')
        chi2_graph.SetPoint(ir, radius, chi2)

        # -----------------------------------------------------------
        # Clean up
        # -----------------------------------------------------------
        del workspace, signal_fitter, bkg_fitter, model_fitter, h_bkg

    # =================================================================
    # Write scan results
    # =================================================================

    scan_dir.cd()

    chi2_graph.Write()
    tree.Write()

    print(tc.GREEN + f'Finished scan for {sign}, centrality {centrality}' + tc.RESET)

def radius_normalisation_scan(outfile: TDirectory, mode: str, centrality: str):

    sign = 'Both' if mode == '' else mode
    print(tc.GREEN + f'\nScanning radius and signal normalisation for mode={mode}, centrality={centrality}' + tc.RESET)

    h_correlation_function, h_same_event, h_mixed_event, h_systematics = \
        get_data_histograms(mode, centrality)

    is_first_bin_empty = h_correlation_function.GetBinContent(1) < 1e-12
    KSTAR_MIN, KSTAR_MAX = (0.02, 0.4) if is_first_bin_empty else (0.01, 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')

    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    scan_dir = outfile.mkdir(f'{sign}_{centrality}')

    chi2_map = TH2D('hChi2VsRadiusAndNormalisation', f'{sign} {centrality};#it{{R}}_{{s}} (fm);#it{{A}}_{{^{{4}}Li}};#chi^{{2}}',
                    len(RADIUS_VALUES), RADIUS_VALUES[0] - 0.05, RADIUS_VALUES[-1] + 0.05,
                    len(SIGNAL_NORMALISATIONS), SIGNAL_NORMALISATIONS[0] - 0.0025, SIGNAL_NORMALISATIONS[-1] + 0.0025)
    chi2_min = np.inf
    best_radius = None
    best_signal_norm = None

    for ir, radius in enumerate(RADIUS_VALUES):

        radius_str = f'{radius:.2f}'
        print(tc.CYAN + f'  R = {radius_str} fm' + tc.RESET)

        h_bkg_name = (f'{sign}/{centrality}/R_{radius_str}/hLambdaSigmaCorrectedCk_Smeared_R_{radius_str}')
        h_bkg = load_hist(INPUT_MODELS_PATH, h_bkg_name)

        workspace = RooWorkspace(f'roows_R_{radius_str}')

        signal_fitter = SignalFitter('signal', kstar_spec, scan_dir, workspace)
        signal_fitter.init_signal('from_kde', h_signal, rho=3)
        signal_fitter.title = '^{4}Li'
        signal_fitter.save_to_workspace()

        bkg_fitter = BkgFitter('bkg', kstar_spec, scan_dir, workspace)
        bkg_fitter.init_bkg('from_kde', h_bkg, rho=0.1)
        bkg_fitter.title = 'Coulomb + strong interaction'
        bkg_fitter.save_to_workspace()

        model_fitter = ModelFitter('model', kstar_spec, scan_dir,
                                   ['signal_pdf'], ['bkg_pdf'], workspace,
                                   extended=True, title='^{4}Li + interaction')

        model_fitter.fractions['signal_pdf'].setRange(0., 1.)
        model_fitter.fractions['signal_pdf'].setVal(0.3)
        model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
        model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')

        model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
        model_fitter.prefit_background(h_correlation_function, 
                                       range_limits=(0.2, 0.4), range_name='bkg_fit_range', 
                                       save_normalisation_value=True)
        sign_label = ('p#minus^{3}He' if mode == 'Matter'
                     else '#bar{p}#minus#bar{^{3}He}' if mode == 'Antimatter'
                     else 'p#minus^{3}He #oplus #bar{p}#minus#bar{^{3}He}')

        for signal_norm in SIGNAL_NORMALISATIONS:

            model_fitter.fractions['signal_pdf'].setVal(signal_norm)
            suffix = f'_R_{radius_str}_A_{signal_norm:.3f}'
            chi2 = model_fitter.compute_chi2_stat_only(h_correlation_function, h_systematics, suffix=suffix,
                                             kstar_max_chi2=0.23)
            #h_chi2 = load_hist('output/radius_scan.root', f'{sign}_{centrality}/model/chi2_model{suffix}')
            #chi2 = h_chi2.GetBinContent(h_chi2.FindBin(0.4))

            bin_x = chi2_map.GetXaxis().FindBin(radius)
            bin_y = chi2_map.GetYaxis().FindBin(signal_norm)

            chi2_map.SetBinContent(bin_x, bin_y, chi2)
            
            if chi2 < chi2_min:
                chi2_min = chi2
                best_radius = radius
                best_signal_norm = signal_norm
        
        del workspace, signal_fitter, bkg_fitter, model_fitter, h_bkg

    min_bin = chi2_map.GetMinimumBin()
    chi2_min = chi2_map.GetBinContent(min_bin)

    g_1sigma_lower, g_1sigma_upper = make_chi2_lines(2.30, chi2_map, chi2_min)
    g_2sigma_lower, g_2sigma_upper = make_chi2_lines(6.18, chi2_map, chi2_min)

    for graph in [g_1sigma_lower, g_1sigma_upper, g_2sigma_lower, g_2sigma_upper]:
        set_root_object(graph, line_color=ROOT.kBlack, line_width=2, 
                        line_style=1)
        
    canvas = TCanvas(f'cChi2VsRadiusAndNormalisation_{sign}_{centrality}', '', 900, 750)
    chi2_map.Draw('COLZ')
    for graph in [g_1sigma_lower, g_1sigma_upper, g_2sigma_lower, g_2sigma_upper]:
        graph.Draw('SAME')

    marker = TMarker(best_radius, best_signal_norm, 20)
    set_root_object(marker, marker_color=ROOT.kBlack, marker_width=1.5, marker_style=20)
    marker.Draw('SAME')
    
    marker_measurement = TGraphErrors()
    marker_measurement.SetPoint(0, MEASUREMENTS[sign][centrality]['R'].n, MEASUREMENTS[sign][centrality]['A'].n)
    marker_measurement.SetPointError(0, MEASUREMENTS[sign][centrality]['R'].s, MEASUREMENTS[sign][centrality]['A'].s)
    set_root_object(marker_measurement, marker_color=ROOT.kRed, marker_width=1.5, marker_style=20,
                    line_color=ROOT.kRed, line_width=2, line_style=1, fill_style=0)
    marker_measurement.Draw('SAME P5')
    
    legend = init_legend(0.68, 0.74, 0.8, 0.88, fill_style=1001, fill_color=ROOT.kWhite)
    legend.AddEntry(marker, f'#chi^{{2}}_{{min}}', 'p')
    legend.AddEntry(marker_measurement, '#it{R}_{s} fixed', 'p')
    legend.Draw('SAME')

    scan_dir.cd()
    chi2_map.Write()
    canvas.Write()

    print(tc.GREEN + f'Finished 2D scan for {sign}, centrality {centrality} (chi2_min = {chi2_min:.2f})' + tc.RESET)

if __name__ == '__main__':

    set_alice_global_style()
    RooMsgService.instance().setGlobalKillBelow(5) # 3 = WARNING, 4 = ERROR, 5 = FATAL
    ROOT.RooFit.PrintLevel(-1)
        
    outfile = TFile('output/radius_scan.root', 'recreate')

    for mode in MODES_TO_SCAN:
        for centrality in CENTRALITIES_TO_SCAN:
            #radius_scan(outfile, mode, centrality)
            radius_normalisation_scan(outfile, mode, centrality)

    outfile.Close()
    print(tc.CYAN + 'Output written to output/radius_scan.root' + tc.RESET)