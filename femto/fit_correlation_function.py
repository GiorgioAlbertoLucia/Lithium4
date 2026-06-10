from ROOT import TFile, TDirectory, gStyle, \
                 RooWorkspace

from core.signal_fitter import SignalFitter
from core.bkg_fitter import BkgFitter
from core.model_fitter import ModelFitter
from core.plot_correlation_over_nsigma import plot_correlation_over_nsigma
from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

#SIGNAL_HIST_LOAD_INFO = HistLoadInfo('output/sampling_check.root', 'hCkHist')
#SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution.root', 'hCkHist')
#SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution_finer_binning.root', 'hCkHist')
SIGNAL_HIST_LOAD_INFO = HistLoadInfo('models/li4_contribution_proper_sill.root', 'hCkHist')
HIST_BKG_NAME = 'hHe3_p_Coul_CF'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/output/correlation.root'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_nclstpc.root'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut_offlinetpc.root'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_23_24_25_pass1_pass4_pass1_reject_multiples_pidintrk.root'
DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/correlation_PbPb_hadronpid.root'
INPUT_SUFFIX = 'PbPb'
SELECTION = 'Default' # 'Default', 'PLiGreaterThan3'
SELECTION_SUFFIX = f'_{SELECTION}' if SELECTION != 'Default' else ''

def prepare_centrality_dict(mode:str, use_systematics:bool=False, upper_radius:bool=False, lower_radius:bool=False,
                            plus_10_percent:bool=False, minus_10_percent:bool=False, smeared_lambda:bool=False):
    if mode != 'Matter' and mode != 'Antimatter' and mode != '':
        raise ValueError('Supported modes are "Matter", "Antimatter" ans "" (for inclusive).')
    
    mode_dir = 'Both' if mode == '' else mode
    prefix = INPUT_SUFFIX if INPUT_SUFFIX != 'PbPb' else 'LHC25_PbPb_pass1'

    if smeared_lambda:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_Smeared', 
                      f'{mode_dir}/010/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1030/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/3050/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/5080/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/080/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1080/hLambdaSigmaCorrectedCk_Smeared'
                      ]
        bkg_paths = [f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root', 
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root']
    elif plus_10_percent:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_Smeared', 
                      f'{mode_dir}/010/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1030/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/3050/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/5080/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/080/hLambdaSigmaCorrectedCk_Smeared',
                      f'{mode_dir}/1080/hLambdaSigmaCorrectedCk_Smeared'
                      ]
        bkg_paths = [f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root', 
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root']
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_LowerLambda', 
                      f'{mode_dir}/010/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/1030/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/3050/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/5080/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/080/hLambdaSigmaCorrectedCk_LowerLambda',
                      f'{mode_dir}/1080/hLambdaSigmaCorrectedCk_LowerLambda',
                      ]
        bkg_paths = [f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root', 
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root'] # Dummy last one
    elif upper_radius:
        hist_names = [f'{mode_dir}/050/upper/hLambdaSigmaCorrectedCk_Smeared_upper', 
                      f'{mode_dir}/010/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/1030/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/3050/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/1050/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/5080/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/080/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      f'{mode_dir}/1080/upper/hLambdaSigmaCorrectedCk_Smeared_upper',
                      ]
        bkg_paths = [f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root', 
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root']
    elif lower_radius:
        hist_names = [f'{mode_dir}/050/lower/hLambdaSigmaCorrectedCk_Smeared_lower', 
                      f'{mode_dir}/010/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/1030/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/3050/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/1050/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/5080/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/080/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      f'{mode_dir}/1080/lower/hLambdaSigmaCorrectedCk_Smeared_lower',
                      ]
        bkg_paths = [f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root', 
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root',
                     f'models/{prefix}_lambda_models.root', f'models/{prefix}_lambda_models.root']
    else:
        hist_names = [HIST_BKG_NAME] * 5
        bkg_paths = [f'models/{prefix}_CATS_converted.root', f'models/{prefix}_CATS_cent0_10_converted.root',  
                     f'models/{prefix}_CATS_cent10_30_converted.root', f'models/{prefix}_CATS_cent30_50_converted.root',
                     f'models/{prefix}_pHe3_square_well_1050_GeV.root',]
    
    data_input_path = [DATA_INPUT_PATH] * 8 if not use_systematics else \
                      ['/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root'] * 5
    h_data_names = [f'Correlation{mode}/{SELECTION}/hCorrelation050', 
                    f'Correlation{mode}/{SELECTION}/hCorrelation010', 
                    f'Correlation{mode}/{SELECTION}/hCorrelation1030', 
                    f'Correlation{mode}/{SELECTION}/hCorrelation3050',
                    f'Correlation{mode}/{SELECTION}/hCorrelationDirectComputation1050', 
                    f'Correlation{mode}/{SELECTION}/hCorrelation5080', 
                    f'Correlation{mode}/{SELECTION}/hCorrelationDirectComputation080', 
                    f'Correlation{mode}/{SELECTION}/hCorrelationDirectComputation1080'
                    ] if not use_systematics else \
                    [f'{mode}/hCorrelation050StatAndSyst', 
                     f'{mode}/hCorrelation010StatAndSyst', 
                     f'{mode}/hCorrelation1030StatAndSyst', 
                     f'{mode}/hCorrelation3050StatAndSyst',
                     f'{mode}/hCorrelationDirectComputation1050StatAndSyst',
                     f'{mode}/hCorrelation5080StatAndSyst',
                     f'{mode}/hCorrelationDirectComputation080StatAndSyst',
                     f'{mode}/hCorrelationDirectComputation1080StatAndSyst',
                    ]
    
    mixed_event_input_path = [DATA_INPUT_PATH] * 8
    h_mixed_event_names = [f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent050', 
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent010', 
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent1030',
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent3050',
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEventDirectComputation1050',
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEvent5080',
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEventDirectComputation080',
                           f'Correlation{mode}/{SELECTION}/hNormalisedMixedEventDirectComputation1080',
                          ]
        
    print(tc.BLUE + f'Preparing centrality dict for mode {mode}' + tc.RESET)
    print(tc.BLUE + f'Bkg paths: {bkg_paths}' + tc.RESET)
    print(tc.BLUE + f'Hist names: {hist_names}' + tc.RESET)
    
    return {
        'name': [f'{mode}050', 
                 f'{mode}010', f'{mode}1030', f'{mode}3050',
                 f'{mode}1050', f'{mode}5080', f'{mode}080', f'{mode}1080'
                 ],
        'bkg_input_path': bkg_paths,
        'h_bkg_name': hist_names,
        'data_input_path': data_input_path,
        'h_data_name': h_data_names,
        'mixed_event_input_path': mixed_event_input_path,
        'h_mixed_event_name': h_mixed_event_names,
    }

def fitting_routine(outfile:TDirectory, bkg_input_path:str, data_input_path:str, mixed_event_input_path:str, 
                    h_bkg_name:str, h_data_name:str, h_mixed_event_name:str,
                    output_pdf:str, mode:str='', centrality:str='',
                    use_smoothening:bool=False):

    workspace = RooWorkspace('roows')
    
    print(tc.CYAN + f'Loading histograms for mode {mode}, centrality {centrality}:' + tc.RESET)
    print(tc.CYAN + f'Background histogram: {h_bkg_name} from {bkg_input_path}' + tc.RESET)
    print(tc.CYAN + f'Signal histogram: {SIGNAL_HIST_LOAD_INFO.hist_name} from {SIGNAL_HIST_LOAD_INFO.hist_file_path}' + tc.RESET)
    print(tc.CYAN + f'Correlation function histogram: {h_data_name} from {data_input_path}' + tc.RESET)
    print(tc.CYAN + f'Mixed event histogram: {h_mixed_event_name} from {mixed_event_input_path}' + tc.RESET)
    
    h_bkg = load_hist(bkg_input_path, h_bkg_name)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = load_hist(data_input_path, h_data_name)
    h_mixed_event = load_hist(mixed_event_input_path, h_mixed_event_name)

    is_first_bin_empty = h_correlation_function.GetBinContent(1) < 1e-12
    KSTAR_MIN, KSTAR_MAX = (0.02, 0.4) if is_first_bin_empty else (0.01, 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    print(f'{h_signal=}')
    signal_fitter.init_signal(signal_init_mode, h_signal)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg, rho=0.1) #(0.05 if '010' not in h_data_name else 0.1)) #, extended=True)
    bkg_fitter.title = 'Coulomb + strong interaction' 
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, 
                               extended=True, title='^{4}Li + interaction')

    model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
    model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Coulomb + strong}')
    
    sign_label = 'p#minus^{3}He' if mode == 'Matter' else '#bar{p}#minus^{3}#bar{He}'
    
    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4), range_name='bkg_fit_range',
                                   save_normalisation_value=True) #, use_chi2_method=False)
    model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf', norm_range='bkg_fit_range',
                           data_label=sign_label)
    model_fitter.save_to_workspace()
    model_fitter.compute_chi2(h_correlation_function)
    model_fitter.compute_raw_yield(h_mixed_event, 'signal_pdf', 'bkg_pdf')
    plot_correlation_over_nsigma(outfile, output_pdf, [KSTAR_MIN, KSTAR_MAX], mode, centrality)

    del workspace, signal_fitter, bkg_fitter, model_fitter

if __name__ == '__main__':

    gStyle.SetOptStat(0)
    use_systematics = False
    use_smoothening = True
    plus_10_percent = False
    minus_10_percent = False
    finer_binning = True
    
    smeared_lambda = True
    upper_radius = False
    lower_radius = False

    systematics_suffix = '_with_systematics' if use_systematics else ''
    smoothening_suffix = '_smoothened' if use_smoothening else ''
    plus_10_percent_suffix = '_plus_10_percent' if plus_10_percent else ''
    minus_10_percent_suffix = '_minus_10_percent' if minus_10_percent else ''
    finer_binning_suffix = '_finer_binning' if finer_binning else ''
    smeared_lambda_suffix = '_smeared_lambda' if smeared_lambda else ''
    upper_radius_suffix = '_upper_radius' if upper_radius else ''
    lower_radius_suffix = '_lower_radius' if lower_radius else ''

    outfile = TFile(f'output/{INPUT_SUFFIX}{SELECTION_SUFFIX}_fit_correlation_function_hadronpid_{systematics_suffix}{upper_radius_suffix}{smoothening_suffix}{lower_radius_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.root', 'recreate')

    for mode in ['', 'Matter', 'Antimatter']:
        
        CENTRALITIES = prepare_centrality_dict(mode, use_systematics=use_systematics, upper_radius=upper_radius, lower_radius=lower_radius,
                                               plus_10_percent=plus_10_percent, minus_10_percent=minus_10_percent, smeared_lambda=smeared_lambda)
        
        for name, bkg_input_path, h_bkg_name, data_input, h_data_name, mixed_event_input_path, h_mixed_event_name in zip(CENTRALITIES['name'], CENTRALITIES['bkg_input_path'], 
                                                                 CENTRALITIES['h_bkg_name'], CENTRALITIES['data_input_path'], CENTRALITIES['h_data_name'], CENTRALITIES['mixed_event_input_path'], CENTRALITIES['h_mixed_event_name']):

            #if '050' in name and '3' not in name:
            #    continue
            
            outdir = outfile.mkdir(name)
            centrality = name.replace(mode, '') if mode != '' else name
            output_pdf = f'figures/{INPUT_SUFFIX}{SELECTION_SUFFIX}/fit_correlation_function_hadronpid_{name}{systematics_suffix}{upper_radius_suffix}{smoothening_suffix}{lower_radius_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.pdf' 

            print('\n\n', tc.GREEN + f'Fitting correlation function for mode {mode}, centrality {centrality}' + tc.RESET)
            fitting_routine(outdir, bkg_input_path=bkg_input_path, data_input_path=data_input, mixed_event_input_path=mixed_event_input_path, 
                            h_bkg_name=h_bkg_name, h_data_name=h_data_name, h_mixed_event_name=h_mixed_event_name,
                            output_pdf=output_pdf, mode=mode, centrality=centrality,
                            use_smoothening=use_smoothening)
    
    outfile.Close()
