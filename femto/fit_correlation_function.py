from ROOT import TFile, TDirectory, gStyle, \
                 RooWorkspace

from core.signal_fitter import SignalFitter
from core.bkg_fitter import BkgFitter
from core.model_fitter import ModelFitter
from core.plot_correlation_over_nsigma import plot_correlation_over_nsigma
from torchic.core.histogram import load_hist, AxisSpec, HistLoadInfo
from torchic.utils.terminal_colors import TerminalColors as tc

#SIGNAL_HIST_LOAD_INFO = HistLoadInfo('output/sampling_check.root', 'hKstar')
SIGNAL_HIST_LOAD_INFO = HistLoadInfo('output/sampling_check.root', 'hCkHist')
HIST_BKG_NAME = 'hHe3_p_Coul_CF'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/output/correlation.root'
#DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_all_pass1_pass4_nclstpc.root'
DATA_INPUT_PATH = '/home/galucia/Lithium4/preparation/checks/correlation_hadronpid_pass1_pass4_nohe3pcut.root'

def prepare_centrality_dict(mode:str):
    if mode != 'Matter' and mode != 'Antimatter':
        raise ValueError('Supported modes are "Matter", "Antimatter"')
    return {
        'name': [f'{mode}050', f'{mode}010', f'{mode}1030', f'{mode}3050'],
        'bkg_input_path': ['models/CATS_converted.root', 'models/CATS_cent0_10_converted.root', 'models/CATS_cent10_30_converted.root', 'models/CATS_cent30_50_converted.root'],
        'h_data_name': [f'Correlation{mode}/hCorrelation050', f'Correlation{mode}/hCorrelation010', f'Correlation{mode}/hCorrelation1030', f'Correlation{mode}/hCorrelation3050']
    }

def fitting_routine(outfile:TDirectory, bkg_input_path:str, data_input_path:str, h_data_name:str, output_pdf:str):

    workspace = RooWorkspace('roows')

    h_bkg = load_hist(bkg_input_path, HIST_BKG_NAME)
    h_signal = load_hist(SIGNAL_HIST_LOAD_INFO)
    h_correlation_function = load_hist(data_input_path, h_data_name)

    KSTAR_MIN, KSTAR_MAX = 0.0, 0.4
    kstar_spec = AxisSpec(100, 0.0, 0.4, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_fitter.init_signal('from_mc', h_signal)
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_fitter.init_bkg('from_mc', h_bkg) #, extended=True)
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, extended=True)
    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4), range_name='bkg_fit_range',
                                   save_normalisation_value=True) #, use_chi2_method=False)
    model_fitter.fit_model(h_correlation_function, norm_range='bkg_fit_range')
    model_fitter.save_to_workspace()
    model_fitter.compute_chi2(h_correlation_function)
    plot_correlation_over_nsigma(outfile, output_pdf, [KSTAR_MIN, KSTAR_MAX])

    del workspace, signal_fitter, bkg_fitter, model_fitter

if __name__ == '__main__':

    gStyle.SetOptStat(0)

    #outfile = TFile('output/fit_correlation_function_all_pass1_pass4_nclstpc.root', 'recreate')
    outfile = TFile('output/fit_correlation_function_hadronpid_pass1_pass4_nohe3pcut.root', 'recreate')

    for mode in ['Matter', 'Antimatter']:
        
        CENTRALITIES = prepare_centrality_dict(mode)
        
        for name, bkg_input_path, h_data_name in zip(CENTRALITIES['name'], CENTRALITIES['bkg_input_path'], CENTRALITIES['h_data_name']):
        
            outdir = outfile.mkdir(name)
            #output_pdf = f'figures/fit_correlation_function_all_pass1_pass4_nclstpc_{name}.pdf'
            output_pdf = f'figures/fit_correlation_function_hadronpid_pass1_pass4_nohe3pcut_{name}.pdf'
            fitting_routine(outdir, bkg_input_path=bkg_input_path, data_input_path=DATA_INPUT_PATH, h_data_name=h_data_name, output_pdf=output_pdf)
    
    outfile.Close()
