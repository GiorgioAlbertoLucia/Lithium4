import sys
import ROOT
from ROOT import TFile, RooWorkspace, TH1F
from tqdm import tqdm

import gc
from torchic.core.histogram import load_hist, AxisSpec
from torchic.utils.root import silence_roofit

sys.path.append('..')
from femto.core.signal_fitter import SignalFitter
from femto.core.bkg_fitter import BkgFitter
from femto.core.model_fitter import ModelFitter


ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)

def fitting_routine(outfile, workspace: RooWorkspace, 
                    h_signal:TH1F,
                    h_bkg:TH1F,
                    h_data:TH1F, h_same_event:TH1F, h_mixed_event:TH1F, 
                    sign:str, centrality:str, use_smoothening:bool=False, 
                    tree_name:str='tree'):
    h_correlation_function = h_data
    h_signal_local = h_signal.Clone()
    h_bkg_local = h_bkg.Clone()

    KSTAR_MIN, KSTAR_MAX = (0.01, 0.4) if h_correlation_function.GetBinWidth(1) > 1e-5 else (0.02, 0.4)
    kstar_spec = AxisSpec(100, KSTAR_MIN, KSTAR_MAX, 'kstar', '#it{k}* (GeV/#it{c})')
    
    signal_fitter = SignalFitter('signal', kstar_spec, outfile, workspace)
    signal_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    signal_fitter.init_signal(signal_init_mode, h_signal_local, tree_name=f'signal_{tree_name}') #, extended=True)
    signal_fitter.title = '^{4}Li'
    signal_fitter.save_to_workspace()

    bkg_fitter = BkgFitter('bkg', kstar_spec, outfile, workspace)
    bkg_init_mode = 'from_mc' if not use_smoothening else 'from_kde'
    bkg_fitter.init_bkg(bkg_init_mode, h_bkg_local, rho=(0.05 if '010' not in centrality else 0.1), tree_name=f'bkg_{tree_name}') #, extended=True)
    bkg_fitter.title = 'Full model' 
    bkg_fitter.save_to_workspace()

    model_fitter = ModelFitter('model', kstar_spec, outfile, ['signal_pdf'], ['bkg_pdf'], workspace, extended=True)
    #model_fitter.fractions['signal_pdf'].setRange(0., 1.)
    model_fitter.load_data(h_correlation_function, h_correlation_function.GetName())
    model_fitter.prefit_background(h_correlation_function, range_limits=(0.2, 0.4), range_name='bkg_fit_range',
                                   save_normalisation_value=True) #, use_chi2_method=False)
    model_fitter.fractions['signal_pdf'].setVal(0.3)
    model_fitter.fractions['signal_pdf'].setRange(-1e4, 1e4)
    model_fitter.fit_model(h_correlation_function, signal_name='signal_pdf', norm_range='bkg_fit_range')
    model_fitter.save_to_workspace()
    #model_fitter.compute_chi2(h_correlation_function)
    xvar = workspace.obj(kstar_spec.name)
    xvar.setRange(KSTAR_MIN, 0.39)
    raw_yield_value = model_fitter.compute_raw_yield(h_same_event, h_mixed_event, 'signal_pdf', 'bkg_pdf')
    
    import faulthandler
    faulthandler.enable()

    signal_fitter.cleanup(keep_histograms=True)
    bkg_fitter.cleanup(keep_histograms=True)
    model_fitter.cleanup(keep_histograms=True)
    del signal_fitter, bkg_fitter, model_fitter

    for obj in (h_signal_local, h_bkg_local):
        obj.ResetBit(ROOT.kMustCleanup)
        del obj

    gc.collect()
    

    return raw_yield_value

def ratio_systematic_routine():

    N_ITERATIONS = 100
    CENTRALITY = '1050'  # combined 10-50% class, built from 1030 + 3050

    infile = TFile.Open("/home/galucia/Lithium4/preparation/output/hist_systematics_with_upper_limit.root")
    outFile = TFile.Open("output/ratio_with_1050.root", "RECREATE")
    workspace = RooWorkspace('roows')

    SIGNAL_PATH = '/home/galucia/Lithium4/femto/models/li4_contribution_proper_sill.root'
    NOMINAL_SIGNAL_HIST = load_hist(SIGNAL_PATH, 'hCkHist')

    BKG_PATH = '/home/galucia/Lithium4/femto/models/lambda_models_LL_10.root'
    # NOTE: check that this key exists for the combined 10-50% class in your bkg model file.
    # If only per-decade keys exist (e.g. '1030', '3050'), pick/derive the appropriate one here.
    BKG_CENTRALITY_KEY = CENTRALITY

    h_ratios = TH1F('hRatios_1050', ';N_{^{4}#bar{Li}}^{raw} / N_{^{4}Li}^{raw};', 200, 0, 2)
    outDirCentrality = outFile.mkdir(CENTRALITY)

    NOMINAL_BKG_HISTS = {
        sign: load_hist(BKG_PATH, f'{sign}/{BKG_CENTRALITY_KEY}/hLambdaSigmaCorrectedCk')
        for sign in ['Matter', 'Antimatter']
    }

    for iter in tqdm(range(N_ITERATIONS), desc=f'Processing ratio - Centrality {CENTRALITY}'):

        outDirIter = outDirCentrality.mkdir(f'iter_{iter}')
        raw_yields = {}

        for sign in ['Matter', 'Antimatter']:

            inDir = infile.Get(sign)

            h_same_1030 = inDir.Get(f'1030/iter_{iter}/hSame')
            h_mixed_normalised_1030 = inDir.Get(f'1030/iter_{iter}/hMixedNormalised')
            h_same_3050 = inDir.Get(f'3050/iter_{iter}/hSame')
            h_mixed_normalised_3050 = inDir.Get(f'3050/iter_{iter}/hMixedNormalised')

            h_same = h_same_1030.Clone(f'hSame_{sign}_{CENTRALITY}_Iter_{iter}')
            h_same.Add(h_same_3050)
            h_mixed_normalised = h_mixed_normalised_1030.Clone(f'hMixedNormalised_{sign}_{CENTRALITY}_Iter_{iter}')
            h_mixed_normalised.Add(h_mixed_normalised_3050)

            if h_same is None or h_mixed_normalised is None:
                print(f"Error: Could not retrieve histograms for {sign} - Centrality {CENTRALITY} - Iteration {iter}")
                continue

            h_correlation = h_same.Clone(f'hCorrelation_{sign}_{CENTRALITY}_{iter}')
            h_correlation.Divide(h_mixed_normalised)

            outDirSign = outDirIter.mkdir(sign)

            raw_yields[sign] = fitting_routine(outDirSign, workspace,
                                                h_signal=NOMINAL_SIGNAL_HIST,
                                                h_bkg=NOMINAL_BKG_HISTS[sign],
                                                h_data=h_correlation,
                                                h_same_event=h_same,
                                                h_mixed_event=h_mixed_normalised,
                                                sign=sign, centrality=CENTRALITY, use_smoothening=True,
                                                tree_name=f'tree_{sign}_{iter}')

            for obj in (h_same, h_mixed_normalised, h_correlation):
                obj.ResetBit(ROOT.kMustCleanup)
                del obj

        if 'Matter' in raw_yields and 'Antimatter' in raw_yields and raw_yields['Matter'] != 0:
            ratio = raw_yields['Antimatter'] / raw_yields['Matter']
            h_ratios.Fill(ratio)
        else:
            print(f"Warning: skipping ratio for iteration {iter} (missing or zero yield)")

        gc.collect()

    outFile.cd()
    h_ratios.Write()

    outFile.Close()
    infile.Close()
    
if __name__ == "__main__":
    
    silence_roofit()
    ratio_systematic_routine()