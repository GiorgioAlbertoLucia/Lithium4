import numpy as np
from ROOT import TFile
from torchic.core.histogram import load_hist


outfile = TFile.Open('output/PbPb/systematic_uncertainties.root', 'RECREATE')

INPUT_FILES = {
    'Both': {
        '010':  'output/systematics_with_upper_limit_Both_010.root',
        '1030': 'output/hist_systematics_with_upper_limit.root',
        '3050': 'output/hist_systematics_with_upper_limit.root',
        '5080': 'output/hist_systematics_with_upper_limit.root',
        '1050': 'output/systematics_with_upper_limit_Both_1050.root',
        '050':  'output/hist_systematics_with_upper_limit.root',
    },
    'Matter': {
        '010':  'output/hist_systematics_with_upper_limit.root',
        '1030': 'output/hist_systematics_with_upper_limit.root',
        '3050': 'output/hist_systematics_with_upper_limit.root',
        '5080': 'output/hist_systematics_with_upper_limit.root',
        '1050': 'output/hist_systematics_with_upper_limit.root',
        '050':  'output/hist_systematics_with_upper_limit.root',
    },
    'Antimatter': {
        '010':  'output/hist_systematics_with_upper_limit.root',
        '1030': 'output/hist_systematics_with_upper_limit.root',
        '3050': 'output/hist_systematics_with_upper_limit.root',
        '5080': 'output/hist_systematics_with_upper_limit.root',
        '1050': 'output/hist_systematics_with_upper_limit.root',
        '050':  'output/hist_systematics_with_upper_limit.root',
    }
}

for sign in ['Both', 'Matter', 'Antimatter']:
    
    sign_str = '' if sign == 'Both' else sign
    outdir = outfile.mkdir(f'Correlation{sign_str}')
    for centrality in ['010', '1030', '3050', '5080', '1050', '050']:
        
        outdir_centrality = outdir.mkdir(f'{centrality}')
        input_file = INPUT_FILES[sign].get(centrality, None)
        if input_file is None:
            print(f'No input file for sign {sign} and centrality {centrality}. Skipping.')
            continue
        infile = TFile.Open(input_file, 'READ')
        nominal_hist = None

        nominal_hist_name = f'Correlation{sign_str}/Default/hCorrelation{centrality}' if centrality != '1050' else f'Correlation{sign_str}/Default/hCorrelationDirectComputation{centrality}'
        nominal_hist = load_hist('/home/galucia/Lithium4/preparation/output/PbPb/correlation_PbPb_hadronpid.root',
                                 nominal_hist_name)
        nbins = nominal_hist.GetNbinsX()

        points = [[] for _ in range(nbins)]
        iter_hists = []
        for iter in range(100):
            iter_hist = None
            
            if 'hist' not in input_file:
                iter_canvas = infile.Get(f'{sign}/{centrality}/iter_{iter}/nominal/model/data_bkg_comparison')
                list_of_primitives = iter_canvas.GetListOfPrimitives()
                for primitive in list_of_primitives:
                    if f'hCorrelation_{centrality}' in primitive.GetName():
                        iter_hist = primitive
                        break
            elif 'hist' in input_file and centrality != '1050' and centrality != '050':
                iter_hist_se = infile.Get(f'{sign}/{centrality}/iter_{iter}/hSame')
                iter_hist_me = infile.Get(f'{sign}/{centrality}/iter_{iter}/hMixedNormalised')
                iter_hist = iter_hist_se.Clone(f'hCorrelation_{centrality}_iter_{iter}')
                iter_hist.Divide(iter_hist_me)
            elif 'hist' in input_file and centrality == '1050':
                iter_hist_se_1030 = infile.Get(f'{sign}/1030/iter_{iter}/hSame')
                iter_hist_me_1030 = infile.Get(f'{sign}/1030/iter_{iter}/hMixedNormalised')
                iter_hist_se_3050 = infile.Get(f'{sign}/3050/iter_{iter}/hSame')
                iter_hist_me_3050 = infile.Get(f'{sign}/3050/iter_{iter}/hMixedNormalised')
                iter_hist_se = iter_hist_se_1030.Clone(f'hSame_1030_3050_iter_{iter}')
                iter_hist_se.Add(iter_hist_se_3050)
                iter_hist_me = iter_hist_me_1030.Clone(f'hMixedNormalised_1030_3050_iter_{iter}')
                iter_hist_me.Add(iter_hist_me_3050)
                iter_hist = iter_hist_se.Clone(f'hCorrelationDirectComputation_{centrality}_iter_{iter}')
                iter_hist.Divide(iter_hist_me)
            elif 'hist' in input_file and centrality == '050':
                iter_hist_se_010 = infile.Get(f'{sign}/010/iter_{iter}/hSame')
                iter_hist_me_010 = infile.Get(f'{sign}/010/iter_{iter}/hMixedNormalised')
                iter_hist_se_1030 = infile.Get(f'{sign}/1030/iter_{iter}/hSame')
                iter_hist_me_1030 = infile.Get(f'{sign}/1030/iter_{iter}/hMixedNormalised')
                iter_hist_se_3050 = infile.Get(f'{sign}/3050/iter_{iter}/hSame')
                iter_hist_me_3050 = infile.Get(f'{sign}/3050/iter_{iter}/hMixedNormalised')
                iter_hist_se = iter_hist_se_1030.Clone(f'hSame_1030_3050_5080_iter_{iter}')
                iter_hist_se.Add(iter_hist_se_3050)
                iter_hist_se.Add(iter_hist_se_010)
                iter_hist_me = iter_hist_me_1030.Clone(f'hMixedNormalised_1030_3050_5080_iter_{iter}')
                iter_hist_me.Add(iter_hist_me_3050)
                iter_hist_me.Add(iter_hist_me_010)
                iter_hist = iter_hist_se.Clone(f'hCorrelationDirectComputation_{centrality}_iter_{iter}')
                iter_hist.Divide(iter_hist_me)
                
            iter_hists.append(iter_hist)

            for ibin in range(1, nbins + 1):
                value = iter_hist.GetBinContent(ibin)
                points[ibin - 1].append(value)
        
        h_corr_syst_errors = nominal_hist.Clone(f'hCorrelationSyst{centrality}')
        systematic_hist = nominal_hist.Clone(f'systematic_hist_{centrality}')
        for ibin in range(1, nbins + 1):
            bin_residuals = points[ibin - 1]
            std_residual = np.std(bin_residuals)
            systematic_hist.SetBinContent(ibin, std_residual)
            systematic_hist.SetBinError(ibin, 0.0)
            h_corr_syst_errors.SetBinContent(ibin, nominal_hist.GetBinContent(ibin))
            h_corr_syst_errors.SetBinError(ibin, std_residual)

        outdir_centrality.cd()
        nominal_hist.Write()
        systematic_hist.Write()
        h_corr_syst_errors.Write()
        
        #for iter, iter_hist in enumerate(iter_hists):
        #    iter_hist.Write(f'iter_{iter}_hist')