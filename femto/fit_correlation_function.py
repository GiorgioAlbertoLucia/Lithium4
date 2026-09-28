import numpy as np

from ROOT import TFile, TDirectory, gStyle, \
                 RooWorkspace, TH1F, TCanvas, TLatex

from core.signal_fitter import SignalFitter
from core.bkg_fitter import BkgFitter
from core.model_fitter import ModelFitter
from core.plot_correlation_over_nsigma import plot_correlation_over_nsigma
from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc
from torchic.utils.root import silence_roofit

from dataclasses import dataclass, field
from typing import Dict, List

import argparse
from core.config_loader import load_yaml

HIST_BKG_NAME = None
DATA_INPUT_PATH = None
INPUT_SUFFIX = None
SELECTION = None
SELECTION_SUFFIX = None

CENTRALITY_BINS = None

AVAILABLE_BKGS = None
AVAILABLE_SIGNALS = None
AVAILABLE_MODEL_PARAM_VARIATIONS = None
SYSTEMATICS_FILE_PATH = None

@dataclass
class FitOptions:
    ground_state_only: bool = False
    use_smoothening: bool = True
    finer_binning: bool = True
    lambda_to_one: bool = False
    lambda_to_zero: bool = False
    lambda_variation: float = 10
    match_ratio: bool = False
    LL: bool = False
    coulomb: bool = False

    _SUFFIX_MAP: Dict[str, str] = field(default_factory=lambda : {
        'ground_state_only': '_ground_state_only',
        'use_smoothening': '_smoothened',
        'finer_binning': '_finer_binning',
        'lambda_to_one': '_dummy_fraction_to_one',
        'lambda_to_zero': '_dummy_fraction_to_zero',
        'LL': '_LL',
        'coulomb': '_Coulomb',
    }, repr=False)

    def _lambda_variation_suffix(self) -> str:
        suffix = ''
        if self.lambda_variation: suffix += f'_{self.lambda_variation}'
        if self.match_ratio: suffix += '_match_ratio'
        return suffix

    def suffix(self) -> str:
        '''Builds the combined filename suffix from whichever flags are True.'''
        base = ''.join(s for name, s in self._SUFFIX_MAP.items() if getattr(self, name))
        return base + self._lambda_variation_suffix()

    def lambda_suffix(self) -> str:
        '''Builds the combined filename suffix for lambda variations.'''
        return self._lambda_variation_suffix()

def _lambda_hist_names(mode_dir: str, template: str) -> list:
    '''template uses {mode_dir} and {cent} placeholders, e.g.
       "{mode_dir}/{cent}/hLambdaSigmaCorrectedCk_Smeared"'''
    return [template.format(mode_dir=mode_dir, cent=c) for c in CENTRALITY_BINS]

def prepare_centrality_dict(mode: str, opts: FitOptions):
    if mode != 'Matter' and mode != 'Antimatter' and mode != '':
        raise ValueError('Supported modes are "Matter", "Antimatter" ans "" (for inclusive).')

    mode_dir = 'Both' if mode == '' else mode
    prefix = INPUT_SUFFIX if INPUT_SUFFIX != 'PbPb' else 'LHC25_PbPb_pass1'

    hist_names = _lambda_hist_names(mode_dir, '{mode_dir}/{cent}/hLambdaSigmaCorrectedCk_Smeared')
    bkg_paths = [f'models/{prefix}_lambda_models.root'] * 8
        
    data_input_path = [DATA_INPUT_PATH] * 8
    h_data_names = [(f'Correlation{mode}/{SELECTION}/hCorrelation{cent}' if (cent != '1050' and cent != '080' and cent != '1080') else 
                     f'Correlation{mode}/{SELECTION}/hCorrelationDirectComputation{cent}') for cent in CENTRALITY_BINS]
        
    mixed_event_input_path = [DATA_INPUT_PATH] * 8
    h_mixed_event_names = [f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent{cent}' if (cent != '1050' and cent != '080' and cent != '1080') else
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEventDirectComputation{cent}' for cent in CENTRALITY_BINS]
    h_same_event_names = [f'Correlation{mode}/{SELECTION}/hSameEvent{cent}' if (cent != '1050' and cent != '080' and cent != '1080') else
                          f'Correlation{mode}/{SELECTION}/hSameEventDirectComputation{cent}' for cent in CENTRALITY_BINS]

    return {
        'name': [f'{mode}{cent}' for cent in CENTRALITY_BINS],
        'bkg_input_path': bkg_paths,
        'h_bkg_name': hist_names,
        'data_input_path': data_input_path,
        'h_data_name': h_data_names,
        'mixed_event_input_path': mixed_event_input_path,
        'h_same_event_name': h_same_event_names,
        'h_mixed_event_name': h_mixed_event_names,
    }

def build_bkg_envelope(h_bkg_variations: List[List[TH1F]], h_ref: TH1F, name_prefix: str = 'h_bkg_envelope') -> tuple:
    '''
    Build the lower and upper envelope histograms across all background variations
    (possibly with different binning), sampled on h_ref's binning via interpolation.
    Returns (h_bkg_low, h_bkg_high), suitable for compute_chi2_stat_only.
    '''
    all_variations = [h for group in h_bkg_variations for h in group]
    if not all_variations:
        return None, None

    h_bkg_low = h_ref.Clone(f'{name_prefix}_low')
    h_bkg_low.Reset()
    h_bkg_high = h_ref.Clone(f'{name_prefix}_high')
    h_bkg_high.Reset()

    for ibin in range(1, h_ref.GetNbinsX() + 1):
        kstar_value = h_ref.GetBinCenter(ibin)
        values = [h_var.Interpolate(kstar_value) for h_var in all_variations]
        h_bkg_low.SetBinContent(ibin, min(values))
        h_bkg_high.SetBinContent(ibin, max(values))

    return h_bkg_low, h_bkg_high

def fitting_routine(outfile:TDirectory, bkg_input_path:str, data_input_path:str, mixed_event_input_path:str, 
                    h_bkg_name:str, h_data_name:str, h_same_event_name:str, h_mixed_event_name:str,
                    output_pdf:str, mode:str='', centrality:str='', use_smoothening:bool=False, 
                    run_variations:bool=False, run_sigma_variations:bool=False, run_model_param_variations:bool=False,
                   variations_bkg_path:str=''):

    workspace = RooWorkspace('roows')
    
    print(tc.CYAN + f'Loading histograms for mode {mode}, centrality {centrality}:' + tc.RESET)
    print(tc.CYAN + f'Background histogram: {h_bkg_name} from {bkg_input_path}' + tc.RESET)
    print(tc.CYAN + f'Signal histogram: {SIGNAL_HIST_LOAD_INFO.hist_name} from {SIGNAL_HIST_LOAD_INFO.hist_file_path}' + tc.RESET)
    print(tc.CYAN + f'Correlation function histogram: {h_data_name} from {data_input_path}' + tc.RESET)
    print(tc.CYAN + f'Mixed event histogram: {h_mixed_event_name} from {mixed_event_input_path}' + tc.RESET)
    
    if not SIGNAL_HIST_LOAD_INFO.hist_file_path or not SIGNAL_HIST_LOAD_INFO.hist_name:
        raise ValueError('Signal histogram file path and name must be provided in SIGNAL_HIST_LOAD_INFO.')
    if not SYSTEMATICS_FILE_PATH:
        raise ValueError('Systematics file path must be provided in SYSTEMATICS_FILE_PATH.')
    if not INPUT_SUFFIX:
        raise ValueError('INPUT_SUFFIX must be provided in the configuration.')

    
    sign = 'Both' if mode == '' else mode
    h_bkg = load_hist(bkg_input_path, h_bkg_name)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = load_hist(data_input_path, h_data_name)
    h_mixed_event = load_hist(mixed_event_input_path, h_mixed_event_name)
    h_same_event = load_hist(data_input_path, h_same_event_name)
    h_systematics_name = f'Correlation{mode}/{centrality}/hCorrelationSyst{centrality}'
    h_systematics = load_hist(SYSTEMATICS_FILE_PATH, h_systematics_name)

    is_first_bin_empty = h_correlation_function.GetBinContent(1) < 1e-12
    KSTAR_MIN, KSTAR_MAX = (0.02, 0.4) if is_first_bin_empty else (0., 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    print(f'{h_signal=}')
    signal_fitter.init_signal(signal_init_mode, h_signal, rho=3)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg, rho=0.1) #(0.05 if '010' not in h_data_name else 0.1)) #, extended=True)
    bkg_fitter.title = 'Coulomb + strong interaction' 
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, 
                               extended=True, title='^{4}Li + interaction')
    
    model_fitter.REFERENCE_KSTAR_VALUE_FOR_BKG_NORMALIZATION = 0.31 # GeV/c - arbitrary value in the region where the background is expected to be dominant to perform the normalisation
    model_fitter.REFERENCE_KSTAR_VALUE_FOR_SIGNAL_NORMALIZATION = 0.07 # GeV/c - arbitrary value in the region where the signal is expected to be dominant to perform the normalisation
    model_fitter.KSTAR_MAX_SIGNIFICANCE = 0.23 # GeV/c - arbitrary value in the region where the signal is expected to be dominant to perform the significance calculation

    model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
    model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')
    
    sign_label = 'p#minus^{3}He' if mode == 'Matter' else ('#bar{p}#minus^{3}#bar{He}' if mode == 'Antimatter' 
                                                           else 'p#minus^{3}He #oplus #bar{p}#minus^{3}#bar{He}')
    
    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4), range_name='bkg_fit_range',
                                   save_normalisation_value=True) #, use_chi2_method=False)
    model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf', norm_range='bkg_fit_range',
                           data_label=sign_label)
    signal_normalisation, signal_normalisation_error = model_fitter.fractions['signal_pdf'].getVal(), model_fitter.fractions['signal_pdf'].getError()
    print(f'Signal normalisation (fraction of signal in the correlation function): {signal_normalisation:.3f} ± {signal_normalisation_error:.3f}')
    
    bkg_normalisation_value = model_fitter.get_bkg_value_at_reference_kstar()
    model_fitter.save_to_workspace()
    
    h_available_bkgs_lambdaR = [load_hist(variations_bkg_path, f'{sign}/{centrality}/{bkg_rel_name}') for bkg_rel_name in AVAILABLE_BKGS] if AVAILABLE_BKGS is not None and len(AVAILABLE_BKGS) > 0 else []
    h_available_bkgs_pars = [load_hist(variations_bkg_path, f'{sign}/{centrality}/{model_param_rel_name}') for model_param_rel_name in AVAILABLE_MODEL_PARAM_VARIATIONS] if AVAILABLE_BKGS is not None and len(AVAILABLE_BKGS) > 0 else []
    h_bkg_variations = [h_available_bkgs_lambdaR, h_available_bkgs_pars]
    h_bkg_low, h_bkg_high = build_bkg_envelope(h_bkg_variations, h_bkg) if len(h_bkg_variations) > 0 else (None, None)
    
    model_fitter.compute_chi2_stat_only(h_correlation_function, h_systematics, h_bkg_low, h_bkg_high)
    model_fitter.compute_chi2_new(h_correlation_function, h_systematics, h_bkg, h_bkg_variations, suffix=f'_{sign}_{centrality}')
    
    model_fitter.compute_raw_yield(h_same_event, h_mixed_event, 'signal_pdf', 'bkg_pdf')
    #plot_correlation_over_nsigma(outfile, output_pdf, [KSTAR_MIN, KSTAR_MAX], mode, centrality)
    plot_correlation_over_nsigma(outfile, output_pdf, [0.001, KSTAR_MAX], sign, centrality, 
                                 #use_systematics=False
                                 available_bkgs=AVAILABLE_BKGS, bkg_file_path=variations_bkg_path, 
                                 normalisation_value=bkg_normalisation_value,
                                 use_systematics=True if h_systematics is not None else False
                                 )
    del workspace, signal_fitter, bkg_fitter, model_fitter

    if run_variations and AVAILABLE_BKGS is not None and len(AVAILABLE_BKGS) > 0:
        h_raw_yields = TH1F('hRawYieldVariations', 'Raw yield variations;Raw yield;Counts', 1600, -200, 1400)
        raw_yields = []
        h_raw_yields_radii = TH1F('hRawYieldVariationsRadii', 'Raw yield variations;Raw yield;Counts', 1600, -200, 1400)
        h_raw_yields_lambda = TH1F('hRawYieldVariationsLambda', 'Raw yield variations;Raw yield;Counts', 1600, -200, 1400)
        prefix = INPUT_SUFFIX if INPUT_SUFFIX != 'PbPb' else 'LHC25_PbPb_pass1'

        for i, bkg_rel_name in enumerate(AVAILABLE_BKGS):
            var_bkg_name = f'{sign}/{centrality}/{bkg_rel_name}'
            var_dir = outfile.mkdir(f'variations/var_{i}')
            var_workspace = RooWorkspace('roows_var')

            var_signal_fitter = SignalFitter('signal', kstar_spec, var_dir, var_workspace)
            var_signal_fitter.init_signal(signal_init_mode, h_signal) #, rho=0.1)
            var_signal_fitter.title = '^{4}Li'
            var_signal_fitter.save_to_workspace()

            h_var_bkg = load_hist(variations_bkg_path, var_bkg_name)
            var_bkg_fitter = BkgFitter('bkg', kstar_spec, var_dir, var_workspace)
            var_bkg_fitter.init_bkg(bkg_init_mode, h_var_bkg, rho=0.1) #rho=0.1)
            var_bkg_fitter.title = 'Coulomb + strong interaction'
            var_bkg_fitter.save_to_workspace()

            var_model_fitter = ModelFitter('model', kstar_spec, var_dir, ['signal_pdf'], ['bkg_pdf'], var_workspace,
                                           extended=True, title='^{4}Li + interaction')
            var_model_fitter.fractions['signal_pdf'].setRange(0., 1.)
            var_model_fitter.fractions['signal_pdf'].setVal(0.3)
            var_model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
            var_model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')

            var_model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
            var_model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4),
                                               range_name='bkg_fit_range', save_normalisation_value=True)
            var_model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf',
                                       norm_range='bkg_fit_range', data_label=sign_label)
            var_model_fitter.save_to_workspace()
            var_model_fitter.compute_chi2_stat_only(h_correlation_function, h_systematics, h_bkg_low, h_bkg_high)
            var_raw_yield = var_model_fitter.compute_raw_yield(h_same_event, h_mixed_event, 'signal_pdf', 'bkg_pdf')

            h_raw_yields.Fill(var_raw_yield)
            raw_yields.append(var_raw_yield)
            if 'nominal' in bkg_rel_name: # nominal radius
                h_raw_yields_lambda.Fill(var_raw_yield)

            if 'higher' not in bkg_rel_name and 'lower_' not in bkg_rel_name and 'upper_lower' not in bkg_rel_name and 'nominal_lower' not in bkg_rel_name: # skip lambda variations for the radii variations
                h_raw_yields_radii.Fill(var_raw_yield)
            
            del var_workspace, var_signal_fitter, var_bkg_fitter, var_model_fitter
        
        if run_sigma_variations and AVAILABLE_SIGNALS is not None and len(AVAILABLE_SIGNALS) > 0:
            h_sigma_raw_yields = TH1F('hRawYieldSigmaVariations', 'Raw yield variations;Raw yield;Counts', 1600, -200, 1400)
            sign = 'Both' if mode == '' else mode
            
            for i, signal_name in enumerate(AVAILABLE_SIGNALS):
                var_dir = outfile.mkdir(f'sigma_variations/var_{i}')
                var_workspace = RooWorkspace('roows_var')
            
                var_signal_fitter = SignalFitter('signal', kstar_spec, var_dir, var_workspace)
                h_var_signal = load_hist(SIGNAL_HIST_LOAD_INFO.hist_file_path, signal_name)
                var_signal_fitter.init_signal(signal_init_mode, h_var_signal) #, rho=0.1)
                var_signal_fitter.title = '^{4}Li'
                var_signal_fitter.save_to_workspace()
            
                var_bkg_fitter = BkgFitter('bkg', kstar_spec, var_dir, var_workspace)
                var_bkg_fitter.init_bkg(bkg_init_mode, h_bkg, rho=0.1)
                var_bkg_fitter.title = 'Coulomb + strong interaction'
                var_bkg_fitter.save_to_workspace()
            
                var_model_fitter = ModelFitter('model', kstar_spec, var_dir, ['signal_pdf'], ['bkg_pdf'], var_workspace,
                                               extended=True, title='^{4}Li + interaction')
                var_model_fitter.fractions['signal_pdf'].setRange(0., 1.)
                var_model_fitter.fractions['signal_pdf'].setVal(0.3)
                signal_suffix = "#sigma + 10%" if 'SigmaUp' in signal_name else ("#sigma - 10%" if 'SigmaDown' in signal_name else "")
                var_model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}'+f' ({signal_suffix})')
                var_model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')
            
                var_model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
                var_model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4),
                                                   range_name='bkg_fit_range', save_normalisation_value=True)
                var_model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf',
                                           norm_range='bkg_fit_range', data_label=sign_label)
                var_model_fitter.save_to_workspace()
                var_model_fitter.compute_chi2_stat_only(h_correlation_function, h_systematics, h_bkg_low, h_bkg_high)
                var_raw_yield = var_model_fitter.compute_raw_yield(h_same_event, h_mixed_event, 'signal_pdf', 'bkg_pdf')
            
                h_sigma_raw_yields.Fill(var_raw_yield)
                
                del var_workspace, var_signal_fitter, var_bkg_fitter, var_model_fitter

        if run_model_param_variations and AVAILABLE_MODEL_PARAM_VARIATIONS is not None and len(AVAILABLE_MODEL_PARAM_VARIATIONS) > 0:
            h_raw_yields_model_params = TH1F('hRawYieldModelParamVariations', 'Raw yield variations;Raw yield;Counts', 1600, -200, 1400)

            for i, model_param_rel_name in enumerate(AVAILABLE_MODEL_PARAM_VARIATIONS):
                model_param_name = f'{sign}/{centrality}/{model_param_rel_name}'
                var_dir = outfile.mkdir(f'model_param_variations/var_{i}')
                var_workspace = RooWorkspace('roows_var')

                var_signal_fitter = SignalFitter('signal', kstar_spec, var_dir, var_workspace)
                var_signal_fitter.init_signal(signal_init_mode, h_signal)
                var_signal_fitter.title = '^{4}Li'
                var_signal_fitter.save_to_workspace()

                h_var_bkg = load_hist(variations_bkg_path, model_param_name)
                var_bkg_fitter = BkgFitter('bkg', kstar_spec, var_dir, var_workspace)
                var_bkg_fitter.init_bkg(bkg_init_mode, h_var_bkg, rho=0.1)
                var_bkg_fitter.title = 'Coulomb + strong interaction'
                var_bkg_fitter.save_to_workspace()

                var_model_fitter = ModelFitter('model', kstar_spec, var_dir, ['signal_pdf'], ['bkg_pdf'], var_workspace,
                                               extended=True, title='^{4}Li + interaction')
                var_model_fitter.fractions['signal_pdf'].setRange(0., 1.)
                var_model_fitter.fractions['signal_pdf'].setVal(0.3)
                var_model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
                var_model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')

                var_model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
                var_model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4),
                                                   range_name='bkg_fit_range', save_normalisation_value=True)
                var_model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf',
                                           norm_range='bkg_fit_range', data_label=sign_label)
                var_model_fitter.save_to_workspace()
                var_model_fitter.compute_chi2_stat_only(h_correlation_function, h_systematics, h_bkg_low, h_bkg_high)
                var_raw_yield = var_model_fitter.compute_raw_yield(h_same_event, h_mixed_event, 'signal_pdf', 'bkg_pdf')

                h_raw_yields_model_params.Fill(var_raw_yield)

                del var_workspace, var_signal_fitter, var_bkg_fitter, var_model_fitter

            canvas = TCanvas('cRawYieldVariations', 'Raw yield variations', 800, 600)
            h_raw_yields.Draw()
            std_dev = np.std(raw_yields, ddof=1)
            text = TLatex(0.15, 0.85, f'#sigma = {std_dev:.2f}')
            text.SetNDC()
            text.SetTextSize(0.04)
            text.Draw()
            
            outfile.cd()
            h_raw_yields.Write()
            h_sigma_raw_yields.Write()
            h_raw_yields_radii.Write()
            h_raw_yields_lambda.Write()
            h_raw_yields_model_params.Write()   
            canvas.Write()

if __name__ == '__main__':

    gStyle.SetOptStat(0)
    silence_roofit()  # Suppress RooFit messages for cleaner output
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='configs/fit_correlation_PbPb.yaml',
                        help='Path to YAML config file')
    args, _ = parser.parse_known_args()

    cfg = load_yaml(args.config)

    HIST_BKG_NAME = cfg['data'].get('bkg_hist_name', 'hHe3_p_Coul_CF')
    DATA_INPUT_PATH = cfg['data']['input_path']
    INPUT_SUFFIX = cfg['data']['suffix']
    SELECTION = cfg['data']['selection']
    SELECTION_SUFFIX = f'_{SELECTION}' if SELECTION != 'Default' else ''
    INPUT_MODELS_PATH = cfg['input_models_path']

    AVAILABLE_BKGS = cfg.get('available_bkgs', None)
    AVAILABLE_SIGNALS = cfg.get('available_signals', None)
    AVAILABLE_MODEL_PARAM_VARIATIONS = cfg.get('available_model_param_variations', None)
    SYSTEMATICS_FILE_PATH = cfg['systematics_file']
    CENTRALITY_BINS = cfg['centrality_bins']
    
    opts = FitOptions(**cfg['fit_options'])
    
    suffix = opts.suffix()
    lambda_suffix = opts.lambda_suffix()
    
    SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution_proper_sill.root', 'hCkHist' if not opts.ground_state_only else  'hCkHist_GroundStateOnly')

    outfile = TFile(f'output/{INPUT_SUFFIX}{SELECTION_SUFFIX}_fit_correlation_function_hadronpid_{suffix}.root', 'recreate')

    for mode in ['', 'Matter', 'Antimatter']:
        
        CENTRALITIES = prepare_centrality_dict(mode, opts)
        
        for (name, bkg_input_path, h_bkg_name, data_input, h_data_name, 
             mixed_event_input_path, h_same_event_name, h_mixed_event_name) in zip(
                 CENTRALITIES['name'], CENTRALITIES['bkg_input_path'], CENTRALITIES['h_bkg_name'], 
                 CENTRALITIES['data_input_path'], CENTRALITIES['h_data_name'], 
                 CENTRALITIES['mixed_event_input_path'], CENTRALITIES['h_same_event_name'], 
                 CENTRALITIES['h_mixed_event_name']):

            #if '050' in name and '3' not in name:
            #    continue
            
            outdir = outfile.mkdir(name)
            centrality = name.replace(mode, '') if mode != '' else name

            print('\n\n', tc.GREEN + f'Fitting correlation function for mode {mode}, centrality {centrality}' + tc.RESET)
            prefix = INPUT_SUFFIX if INPUT_SUFFIX != 'PbPb' else 'LHC25_PbPb_pass1'
            
            #INPUT_MODELS_PATH = f'models/{prefix}_lambda_models{lambda_to_one_suffix}{lambda_to_zero_suffix}.root'
            
            output_pdf = f'figures/{INPUT_SUFFIX}{SELECTION_SUFFIX}/fit_correlation_function_hadronpid_{name}{suffix}.pdf' 
            
            fitting_routine(outdir, bkg_input_path=INPUT_MODELS_PATH, data_input_path=data_input, mixed_event_input_path=mixed_event_input_path, 
                            h_bkg_name=h_bkg_name, h_data_name=h_data_name, h_same_event_name=h_same_event_name, h_mixed_event_name=h_mixed_event_name,
                            output_pdf=output_pdf, mode=mode, centrality=centrality,
                            use_smoothening=opts.use_smoothening, run_variations=True, run_sigma_variations=True, run_model_param_variations=True, 
                            variations_bkg_path=INPUT_MODELS_PATH)
    
    print('Output written to', tc.CYAN+outfile.GetName()+tc.RESET)
    outfile.Close()
