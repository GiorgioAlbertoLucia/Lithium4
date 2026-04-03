
from ROOT import TFile, TCanvas, \
                 RooRealVar, RooDataHist, RooGenericPdf, RooAddPdf, RooFitResult, \
                    kRed, kDashed, kGreen, kOrange

from torchic.core.histogram import load_hist
from torchic.roopdf import RooGausDExp, RooGausExp
from torchic.core.graph import load_graph


# Add this after the imports and before if __name__ == '__main__':

class FitConfig:
    """Configuration for a single fit at a specific pt value"""
    def __init__(self, pt, signal_pars=None, bkg_pars=None, fit_range=(-6, 8), integral_range=(-2, 2)):
        self.pt = pt
        self.fit_range = fit_range
        self.integral_range = integral_range
        
        # Default signal parameters
        self.signal_pars = signal_pars or {
            'mean': (-0.5, 0., 0.5),
            'sigma': (0.5, 1., 2.),
            'rlife0': (-1.2, -0.8, -0.5),
            'rlife1': (0.5, 2., 1.2),
        }
        
        # Default background parameters
        self.bkg_pars = bkg_pars or {
            'alpha': (0.001, 0.14, 0.3),
            'offset': (-20., -5., 0.),
            'mean': (-10, -5., -2),
            'sigma': (0.1, 1.0, 5.),
            'rlife': (0.5, 2., 10.),
        }

def perform_fit(input_th2, config, output_dir):
    """Perform a single fit and return the purity value"""
    pt = config.pt
    
    h_slice = input_th2.ProjectionY(f'purity_slice_{pt}', 
                                     input_th2.GetXaxis().FindBin(pt), 
                                     input_th2.GetXaxis().FindBin(pt))

    nsigma = RooRealVar('nsigma', 'nsigma', *config.fit_range)
    datahist = RooDataHist('datahist', 'datahist', [nsigma], Import=h_slice)

    # Create signal parameters from config (min, init, max)
    signal_pars = {
        'mean': RooRealVar('signal_mean', '#mu', *config.signal_pars['mean']),
        'sigma': RooRealVar('signal_sigma', '#sigma', *config.signal_pars['sigma']),
        'rlife0': RooRealVar('signal_rlife0', 'rlife0', *config.signal_pars['rlife0']),
        'rlife1': RooRealVar('signal_rlife1', 'rlife1', *config.signal_pars['rlife1']),
    }
    signal_pdf = RooGausDExp('signal', 'signal', nsigma, *signal_pars.values())

    # Create background parameters from config
    bkg_pars = {
        'alpha': RooRealVar('alpha', 'alpha', *config.bkg_pars['alpha']),
        'offset': RooRealVar('offset', 'offset', *config.bkg_pars['offset']),
        'mean': RooRealVar('mean_bkg', 'mean_bkg', *config.bkg_pars['mean']),
        'sigma': RooRealVar('sigma_bkg', 'sigma_bkg', *config.bkg_pars['sigma']),
        'rlife': RooRealVar('bkg_rlife', 'rlife', *config.bkg_pars['rlife']),
    }
    
    exp = RooGenericPdf('bkg_exp', 'bkg', 
                        f'exp(-alpha * ({nsigma.GetName()} - offset))', 
                        [nsigma, bkg_pars['alpha'], bkg_pars['offset']])
    gausexp = RooGausExp('bkg_gausexp', 'bkg', nsigma, 
                         bkg_pars['mean'], bkg_pars['sigma'], bkg_pars['rlife'])
    
    norm_signal = RooRealVar('fraction_signal', 'fraction_signal', h_slice.Integral(), 0., h_slice.Integral()*10.)
    norm_exp_bkg = RooRealVar('fraction_exp_bkg', 'fraction_exp_bkg', h_slice.Integral(), 0., h_slice.Integral()*10.)
    norm_gausexp_bkg = RooRealVar('fraction_gausexp_bkg', 'fraction_gausexp_bkg', h_slice.Integral(), 0., h_slice.Integral()*10.)
    model = RooAddPdf('model', 'model', [signal_pdf, exp, gausexp], [norm_signal, norm_exp_bkg, norm_gausexp_bkg])

    fit_result = model.fitTo(datahist, Save=True)

    frame = nsigma.frame()
    datahist.plotOn(frame)
    model.plotOn(frame)
    model.plotOn(frame, Components=signal_pdf, LineColor=kRed, LineStyle=kDashed)
    model.plotOn(frame, Components=exp, LineColor=kGreen, LineStyle=kDashed)
    model.plotOn(frame, Components=gausexp, LineColor=kOrange, LineStyle=kDashed)
    model.paramOn(frame, Layout=(0.6, 0.9, 0.9), Format='NEU', ShowConstants=False)

    canvas = TCanvas(f'canvas_{pt}', f'canvas_{pt}', 800, 600)
    frame.Draw()

    output_dir.cd()
    canvas.Write(f'fit_pt_{pt}')

    nsigma.setRange('integral_range', *config.integral_range)

    int_signal = (signal_pdf.createIntegral(nsigma, nsigma, 'integral_range').getVal() 
                 * norm_signal.getVal())
    
    int_bkg = 0.
    int_bkg += (exp.createIntegral(nsigma, nsigma, 'integral_range').getVal() 
                * norm_exp_bkg.getVal())
    int_bkg += (gausexp.createIntegral(nsigma, nsigma, 'integral_range').getVal() 
                * norm_gausexp_bkg.getVal())

    purity = int_signal / (int_signal + int_bkg) if (int_signal + int_bkg) > 0 else 0.
    
    return purity

if __name__ == '__main__':
    
    input_th2 = load_hist('/home/galucia/Lithium4/calibration/output/purity_23_pass4_p.root',
                          'Had/TOF/h2_nsigmaHadTOF_preselection')
    
    input_purity_graph_matter = load_graph('/home/galucia/Lithium4/calibration/output/purity_23_pass4_p.root',
                                     'Had/TOF/g_purity_matter')
    input_purity_graph_antimatter = load_graph('/home/galucia/Lithium4/calibration/output/purity_23_pass4_p.root',
                                     'Had/TOF/g_purity_antimatter')
    
    output_file = TFile.Open('/home/galucia/Lithium4/calibration/output/purity_23_pass4_p_v1.root', 'recreate')
    output_dir_fits = output_file.mkdir('repeated_fits')

    fit_configs = [
        
        FitConfig(
            pt=2.575,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -2., -0.6),
                'rlife1': (1.5, 0.5, 3.),
            },
            bkg_pars={
                'alpha': (0.001, 0.0001, 0.2),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),
        
        FitConfig(
            pt=3.775,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -2., -0.5),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.14, 0.001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),
        
        FitConfig(
            pt=3.825,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -2., -0.6),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.14, 0.001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),

        FitConfig(
            pt=3.875,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -2., -0.6),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.14, 0.001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),

        FitConfig(
            pt=-3.725,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -3., -0.8),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.001, 0.0001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),
        
        FitConfig(
            pt=-3.775,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -3., -0.8),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.001, 0.0001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),

        FitConfig(
            pt=-3.825,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -3., -0.8),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.08, 0.0001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),
        
        FitConfig(
            pt=-3.875,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -3., -0.8),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.08, 0.0001, 0.3),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),

        FitConfig(
            pt=-3.925,
            signal_pars={
                'mean': (-0.5, 0., 0.5),
                'sigma': (0.5, 1., 2.),
                'rlife0': (-1.2, -3., -0.8),
                'rlife1': (1.2, 0.5, 2.),
            },
            bkg_pars={
                'alpha': (0.08, 0.0001, 0.1),
                'offset': (-5., -20., 0.),
                'mean': (-5., -10, -2),
                'sigma': (1.0, 0.1, 5.),
                'rlife': (2., 0.5, 10.),
            },
            fit_range=(-6, 8),
            integral_range=(-2, 2)
        ),
    ]

    purity_results_matter = {}
    purity_results_antimatter = {}
    
    for config in fit_configs:
        purity = perform_fit(input_th2, config, output_dir_fits)
        
        # Separate results by sign of pt
        if config.pt >= 0:
            purity_results_matter[config.pt] = purity
        else:
            purity_results_antimatter[config.pt] = purity
            
        print(f'Purity at pt={config.pt}: {purity:.4f}')
    
    output_file.cd()
    if purity_results_matter:
        g_purity_matter = input_purity_graph_matter.Clone('purity_graph_matter_updated')
        for pt, purity in purity_results_matter.items():
            for i in range(g_purity_matter.GetN()):
                x = g_purity_matter.GetX()[i]
                if abs(x - pt) < 1e-3:
                    g_purity_matter.SetPoint(i, x, purity)
                    break
        g_purity_matter.Write()
    
    if purity_results_antimatter:
        g_purity_antimatter = input_purity_graph_antimatter.Clone('purity_graph_antimatter_updated')
        for pt, purity in purity_results_antimatter.items():
            pt = abs(pt)
            for i in range(g_purity_antimatter.GetN()):
                x = g_purity_antimatter.GetX()[i]
                if abs(x - pt) < 1e-3:
                    g_purity_antimatter.SetPoint(i, x, purity)
                    break
        g_purity_antimatter.Write()
    
    output_file.Close()