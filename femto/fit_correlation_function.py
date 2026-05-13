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
DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/output/PbPb/correlation_PbPb_hadronpid.root'
INPUT_SUFFIX = 'PbPb'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_23_24_25_pass1_pass4_pass1_reject_multiples_pidintrk.root'

def prepare_centrality_dict(mode:str, use_lambda_parameters:bool=False, use_strong_interaction:bool=False, 
                            use_systematics:bool=False, use_rebinning:bool=False, consider_sigma:bool=False,
                            plus_10_percent:bool=False, minus_10_percent:bool=False, smeared_lambda:bool=False):
    if mode != 'Matter' and mode != 'Antimatter' and mode != '':
        raise ValueError('Supported modes are "Matter", "Antimatter" ans "" (for inclusive).')
    
    mode_dir = 'Both' if mode == '' else mode

    if smeared_lambda:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_Smeared', 
                      #f'{mode_dir}/010/hLambdaSigmaCorrectedCk_Smeared'
                      f'{mode_dir}/010/nominal/hLambdaSigmaCorrectedCk_matched_for_smearing_matched_smeared'
                      ] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_Smeared']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one

    elif plus_10_percent:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_HigherLambda', f'{mode_dir}/010/hLambdaSigmaCorrectedCk_HigherLambda'] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_HigherLambda']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one
    elif minus_10_percent:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk_LowerLambda', f'{mode_dir}/010/hLambdaSigmaCorrectedCk_LowerLambda'] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaSigmaCorrectedCk_LowerLambda']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one
    elif consider_sigma:
        hist_names = [f'{mode_dir}/050/hLambdaSigmaCorrectedCk', f'{mode_dir}/010/hLambdaSigmaCorrectedCk'] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaSigmaCorrectedCk']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one
    elif use_rebinning:
        hist_names = [f'{mode_dir}/050/hLambdaCorrectedCk_Matched', f'{mode_dir}/010/hLambdaCorrectedCk_Matched'] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaCorrectedCk_Matched']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one
    elif use_lambda_parameters:
        hist_names = [f'{mode_dir}/050/hLambdaCorrectedCk', f'{mode_dir}/010/hLambdaCorrectedCk'] + [HIST_BKG_NAME] * 2 + [f'{mode_dir}/1050/hLambdaCorrectedCk']
        bkg_paths = ['models/lambda_models.root', 'models/lambda_models.root', 
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/lambda_models.root'] # Dummy last one
    elif use_strong_interaction:    
        hist_names = [HIST_BKG_NAME] * 5
        bkg_paths = ['models/pHe3_square_well_050_GeV.root', 'models/pHe3_square_well_010_GeV.root',
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/pHe3_square_well_1050_GeV.root',]
    else:
        hist_names = [HIST_BKG_NAME] * 5
        bkg_paths = ['models/CATS_converted.root', 'models/CATS_cent0_10_converted.root',  
                     'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root',
                     'models/pHe3_square_well_1050_GeV.root',]
    
    data_input_path = [DATA_INPUT_PATH] * 5 if not use_systematics else \
                      ['/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root'] * 5
    h_data_names = [f'Correlation{mode}/Default/hCorrelation050', f'Correlation{mode}/Default/hCorrelation010', 
                    f'Correlation{mode}/Default/hCorrelation1030', f'Correlation{mode}/Default/hCorrelation3050',
                    f'Correlation{mode}/Default/hCorrelationDirectComputation1050',
                    ] if not use_systematics else \
                    [f'{mode}/hCorrelation050StatAndSyst', f'{mode}/hCorrelation010StatAndSyst', 
                    f'{mode}/hCorrelation1030StatAndSyst', f'{mode}/hCorrelation3050StatAndSyst',
                    f'{mode}/hCorrelationDirectComputation1050StatAndSyst',
                    ]
    
    mixed_event_input_path = [DATA_INPUT_PATH] * 5
    h_mixed_event_names = [f'Correlation{mode}/Default/hNormalisedMixedEvent050', f'Correlation{mode}/Default/hNormalisedMixedEvent010', 
                    f'Correlation{mode}/Default/hNormalisedMixedEvent1030', f'Correlation{mode}/Default/hNormalisedMixedEvent3050',
                    f'Correlation{mode}/Default/hNormalisedMixedEventDirectComputation1050',
                    ]
        
    print(tc.BLUE + f'Preparing centrality dict for mode {mode} with lambda parameters: {use_lambda_parameters}, strong interaction: {use_strong_interaction}' + tc.RESET)
    print(tc.BLUE + f'Bkg paths: {bkg_paths}' + tc.RESET)
    print(tc.BLUE + f'Hist names: {hist_names}' + tc.RESET)
    
    return {
        'name': [f'{mode}050', 
                 f'{mode}010', f'{mode}1030', f'{mode}3050',
                 f'{mode}1050'
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

    h_bkg = load_hist(bkg_input_path, h_bkg_name)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = load_hist(data_input_path, h_data_name)
    h_mixed_event = load_hist(mixed_event_input_path, h_mixed_event_name)

    #KSTAR_MIN, KSTAR_MAX = ((0.01, 0.4) if ('010' not in h_data_name or mode != 'Antimatter' ) 
    #                        else (0.02, 0.4))
    KSTAR_MIN, KSTAR_MAX = (0.02, 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    signal_fitter.init_signal(signal_init_mode, h_signal)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg, rho=(0.05 if '010' not in h_data_name else 0.1)) #, extended=True)
    #bkg_fitter.title = 'Full model' 
    bkg_fitter.title = 'Coulomb + strong interaction' 
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, 
                               extended=True, title='^{4}Li + interaction')

    model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].SetTitle('#it{A}_{^{4}Li}')
    #model_fitter.fractions['bkg_pdf'].SetTitle('#it{A}_{Full model}')
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
    use_lambda_parameters = False
    use_strong_interaction = False
    use_systematics = False
    use_smoothening = True
    use_rebinning = False
    consider_sigma = False
    plus_10_percent = False
    minus_10_percent = False
    finer_binning = True
    smeared_lambda = True

    suffix_lambda = '_lambda_parameters' if use_lambda_parameters else ''
    suffix_strong = '_strong_interaction' if use_strong_interaction else ''
    systematics_suffix = '_with_systematics' if use_systematics else ''
    smoothening_suffix = '_smoothened' if use_smoothening else ''
    rebinning_suffix = '_rebinning' if use_rebinning else ''
    sigma_suffix = '_consider_sigma' if consider_sigma else ''
    plus_10_percent_suffix = '_plus_10_percent' if plus_10_percent else ''
    minus_10_percent_suffix = '_minus_10_percent' if minus_10_percent else ''
    finer_binning_suffix = '_finer_binning' if finer_binning else ''
    smeared_lambda_suffix = '_smeared_lambda' if smeared_lambda else ''

    #outfile = TFile('output/fit_correlation_function_all_pass1_pass4_nclstpc.root', 'recreate')
    #outfile = TFile(f'output/fit_correlation_function_hadronpid_pass1_pass4_nohe3pcut_offlinetpc{suffix_lambda}{suffix_strong}.root', 'recreate')
    outfile = TFile(f'output/{INPUT_SUFFIX}_fit_correlation_function_hadronpid_{suffix_lambda}{suffix_strong}{systematics_suffix}{sigma_suffix}{smoothening_suffix}{rebinning_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.root', 'recreate')
    #outfile = TFile(f'output/fit_correlation_function_23_24_25{suffix_lambda}{suffix_strong}{systematics_suffix}{sigma_suffix}{smoothening_suffix}{rebinning_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.root', 'recreate')

    #for mode in ['', 'Matter', 'Antimatter']:
    for mode in ['Matter', 'Antimatter']:
        
        CENTRALITIES = prepare_centrality_dict(mode, use_lambda_parameters=use_lambda_parameters, use_strong_interaction=use_strong_interaction,
                                               use_systematics=use_systematics, use_rebinning=use_rebinning, consider_sigma=consider_sigma,
                                               plus_10_percent=plus_10_percent, minus_10_percent=minus_10_percent)
        
        for name, bkg_input_path, h_bkg_name, data_input, h_data_name, mixed_event_input_path, h_mixed_event_name in zip(CENTRALITIES['name'], CENTRALITIES['bkg_input_path'], 
                                                                 CENTRALITIES['h_bkg_name'], CENTRALITIES['data_input_path'], CENTRALITIES['h_data_name'], CENTRALITIES['mixed_event_input_path'], CENTRALITIES['h_mixed_event_name']):
        
            outdir = outfile.mkdir(name)
            centrality = name.replace(mode, '') if mode != '' else name
            #output_pdf = f'figures/fit_correlation_function_all_pass1_pass4_nclstpc_{name}.pdf'
            #output_pdf = f'figures/fit_correlation_function_hadronpid_pass1_pass4_nohe3pcut_offlinetpc_{name}{suffix_lambda}{suffix_strong}.pdf'
            output_pdf = f'figures/{INPUT_SUFFIX}/fit_correlation_function_hadronpid_{name}{suffix_lambda}{suffix_strong}{systematics_suffix}{sigma_suffix}{smoothening_suffix}{rebinning_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.pdf' 
            #output_pdf = f'figures/fit_correlation_function_23_24_25{name}{suffix_lambda}{suffix_strong}{systematics_suffix}{sigma_suffix}{smoothening_suffix}{rebinning_suffix}{plus_10_percent_suffix}{minus_10_percent_suffix}{finer_binning_suffix}{smeared_lambda_suffix}.pdf' 

            fitting_routine(outdir, bkg_input_path=bkg_input_path, data_input_path=data_input, mixed_event_input_path=mixed_event_input_path, 
                            h_bkg_name=h_bkg_name, h_data_name=h_data_name, h_mixed_event_name=h_mixed_event_name,
                            output_pdf=output_pdf, mode=mode, centrality=centrality,
                            use_smoothening=use_smoothening)
    
    outfile.Close()
