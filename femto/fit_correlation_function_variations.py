import os
from ROOT import TFile, TDirectory, TH1F, gStyle, RooWorkspace

from core.signal_fitter import SignalFitter
from core.bkg_fitter import BkgFitter
from core.model_fitter import ModelFitter
from core.plot_correlation_over_nsigma import plot_correlation_over_nsigma
from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

# ---------------------------------------------------------------------------
# Global configuration
# ---------------------------------------------------------------------------

SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution_proper_sill.root', 'hCkHist')
DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_refined_dca.root'

# Single .root file that contains the background variation histograms
# organised as: <folder>/<hist_name>  (e.g. nominal/hLambdaSigmaCorrectedCk_matched_smeared)
BKG_VARIATIONS_FILE = 'models/lambda_models.root'

# The three sub-folders inside BKG_VARIATIONS_FILE
BKG_VARIATION_FOLDERS = ['nominal', 'upper', 'lower']

# The three histogram names inside each folder
BKG_VARIATION_HISTS = [
    'hLambdaSigmaCorrectedCk_matched_for_smearing_matched_smeared',
    'hLambdaSigmaCorrectedCk_Higher_matched_for_smearing_matched_smeared',
    'hLambdaSigmaCorrectedCk_Lower_matched_for_smearing_matched_smeared',
]

# ---------------------------------------------------------------------------
# Helper: build the 9-entry variations list
# ---------------------------------------------------------------------------

def get_bkg_variations():
    """Return list of (label, h_bkg_name) for all 9 folder×hist combinations."""
    variations = []
    for folder in BKG_VARIATION_FOLDERS:
        for hist in BKG_VARIATION_HISTS:
            # e.g. "nominal_Higher_matched_for_smearing_matched_smeared"
            suffix = hist.split('hLambdaSigmaCorrectedCk_')[1]
            label = f'{folder}_{suffix}'
            h_bkg_name = f'{folder}/{hist}'
            variations.append((label, h_bkg_name))
    return variations  # 9 entries


# ---------------------------------------------------------------------------
# Centrality configuration (data / mixed-event side only — bkg comes from variations)
# ---------------------------------------------------------------------------

def prepare_centrality_dict(mode: str):
    if mode not in ('Matter', 'Antimatter', ''):
        raise ValueError('Supported modes are "Matter", "Antimatter" and "" (inclusive).')

    data_input_path = [DATA_INPUT_PATH] #* 5
    h_data_names = [
        #f'Correlation{mode}/Default/hCorrelation050',
        f'Correlation{mode}/Default/hCorrelation010',
        #f'Correlation{mode}/Default/hCorrelation1030',
        #f'Correlation{mode}/Default/hCorrelation3050',
        #f'Correlation{mode}/Default/hCorrelationDirectComputation1050',
    ]
    mixed_event_input_path = [DATA_INPUT_PATH] #* 5
    h_mixed_event_names = [
        #f'Correlation{mode}/Default/hNormalisedMixedEvent050',
        f'Correlation{mode}/Default/hNormalisedMixedEvent010',
        #f'Correlation{mode}/Default/hNormalisedMixedEvent1030',
        #f'Correlation{mode}/Default/hNormalisedMixedEvent3050',
        #f'Correlation{mode}/Default/hNormalisedMixedEventDirectComputation1050',
    ]

    print(tc.BLUE + f'Preparing centrality dict for mode: {mode}' + tc.RESET)

    return {
        'name':                   [#f'{mode}050', 
                                   f'{mode}010', #f'{mode}1030', f'{mode}3050', f'{mode}1050'
                                   ],
        'data_input_path':        data_input_path,
        'h_data_name':            h_data_names,
        'mixed_event_input_path': mixed_event_input_path,
        'h_mixed_event_name':     h_mixed_event_names,
    }


# ---------------------------------------------------------------------------
# Core fitting routine — returns (raw_yield, raw_yield_err)
# ---------------------------------------------------------------------------

def fitting_routine(outdir: TDirectory,
                    bkg_input_path: str, h_bkg_name: str,
                    data_input_path: str, mixed_event_input_path: str,
                    h_data_name: str, h_mixed_event_name: str,
                    output_pdf: str,
                    mode: str = '', centrality: str = '',
                    use_smoothening: bool = True):

    workspace = RooWorkspace('roows')

    h_bkg                  = load_hist(bkg_input_path, h_bkg_name, debug=True)
    h_signal               = load_hist(SIGNAL_HIST_LOAD_INFO, debug=True)
    h_correlation_function = load_hist(data_input_path, h_data_name, debug=True)
    h_mixed_event          = load_hist(mixed_event_input_path, h_mixed_event_name, debug=True)

    KSTAR_MIN, KSTAR_MAX = (
        (0.02, 0.4) if ('010' in h_data_name and mode == 'Antimatter') else (0.01, 0.4)
    )
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')

    init_mode = 'from_kde' if use_smoothening else 'from_mc'

    signal_fitter = SignalFitter('signal', kstar_spec, outdir, workspace)
    signal_fitter.init_signal(init_mode, h_signal)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outdir, workspace)
    bkg_fitter.init_bkg(init_mode, h_bkg, rho=(0.1 if '010' in h_data_name else 0.05))
    bkg_fitter.title = 'Full model'
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outdir, ['signal_pdf'], ['bkg_pdf'],
                               workspace, extended=True)
    model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
    model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Full model}')

    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function,
                                   range_limits=(0.2, 0.4),
                                   range_name='bkg_fit_range',
                                   save_normalisation_value=True)
    model_fitter.fit_model(h_correlation_function,
                           signal_name='signal_pdf',
                           norm_range='bkg_fit_range')
    model_fitter.save_to_workspace()
    model_fitter.compute_chi2(h_correlation_function)

    raw_yield = model_fitter.compute_raw_yield(h_mixed_event, 'signal_pdf', 'bkg_pdf')

    plot_correlation_over_nsigma(outdir, output_pdf, [KSTAR_MIN, KSTAR_MAX], mode, centrality)

    del workspace, signal_fitter, bkg_fitter, model_fitter

    return raw_yield


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    gStyle.SetOptStat(0)
    use_smoothening = True

    os.makedirs('figures/variations', exist_ok=True)

    outfile = TFile('output/fit_correlation_variations_smeared_lambda.root', 'recreate')

    bkg_variations = get_bkg_variations()   # list of (label, h_bkg_name), length 9
    n_variations   = len(bkg_variations)    # == 9

    for mode in ['Matter', 'Antimatter']:

        CENTRALITIES = prepare_centrality_dict(mode)
        sign_bkg_variations = bkg_variations.copy()  # will be modified with mode/centrality-specific paths

        for (name,
             data_input,        h_data_name,
             mixed_event_input, h_mixed_event_name) in zip(
                CENTRALITIES['name'],
                CENTRALITIES['data_input_path'],  CENTRALITIES['h_data_name'],
                CENTRALITIES['mixed_event_input_path'], CENTRALITIES['h_mixed_event_name']):

            centrality = name.replace(mode, '') if mode != '' else name
            for idx, (label, bkg_name) in enumerate(sign_bkg_variations):
                tmp_bkg_name = f'{mode}/{centrality}/{bkg_name}'
                sign_bkg_variations[idx] = (label, tmp_bkg_name)
                print(tc.BLUE + f'Checking background variation: {sign_bkg_variations[idx]}' + tc.RESET)

            print(tc.BLUE + f'\n=== [{mode}] centrality: {centrality} ===' + tc.RESET)

            # Summary histogram: one bin per background variation
            h_raw_yield = TH1F(
                f'hRawYield_{name}',
                f'Raw yield variations {name};#it{{N}}_{{raw}} (^{{4}}Li);Counts',
                500, 0, 500
            )
            
            # Create folder structure:  <name>/variations/<label>/
            cent_dir      = outfile.mkdir(name)
            variations_dir = cent_dir.mkdir('variations')

            for i, (label, h_bkg_name) in enumerate(sign_bkg_variations):

                print(tc.BLUE +
                      f'  variation {i + 1}/{n_variations}: {label}' +
                      tc.RESET)

                var_subdir = variations_dir.mkdir(label)
                output_pdf = f'figures/variations/fit_{name}_{label}.pdf'

                raw_yield = fitting_routine(
                    var_subdir,
                    bkg_input_path        = BKG_VARIATIONS_FILE,
                    h_bkg_name            = h_bkg_name,
                    data_input_path       = data_input,
                    mixed_event_input_path= mixed_event_input,
                    h_data_name           = h_data_name,
                    h_mixed_event_name    = h_mixed_event_name,
                    output_pdf            = output_pdf,
                    mode                  = mode,
                    centrality            = centrality,
                    use_smoothening       = use_smoothening,
                )

                h_raw_yield.Fill(raw_yield)

            # Write the summary histogram to the centrality's top-level folder
            cent_dir.cd()
            h_raw_yield.Write()

    outfile.Close()
    print(tc.BLUE + '\nDone. Output: output/fit_correlation_variations_smeared_lambda.root' + tc.RESET)