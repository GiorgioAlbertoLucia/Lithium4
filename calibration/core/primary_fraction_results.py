"""
Results drawing and analysis
"""

import numpy as np
import pandas as pd
from ROOT import TFile, TDirectory, TH2F, TGraphErrors, TH2D, TF1

from torchic.core.graph import create_graph
from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object

from core.primary_fraction_config import FitConfig


def accumulate_covariance_matrices(fit_results: pd.DataFrame, 
                                   out_dir: TDirectory,
                                   sign: str) -> None:
    """
    Accumulate covariance matrices from all fits into a single TH2
    
    Args:
        fit_results: DataFrame with fit results including covariance matrices
        out_dir: Output directory
        sign: 'matter' or 'antimatter'
    """
    # Find maximum matrix size
    max_size = 0
    valid_cov_count = 0
    
    for idx, row in fit_results.iterrows():
        if row['covariance'] is not None:
            max_size = max(max_size, row['covariance'].shape[0])
            valid_cov_count += 1
    
    if max_size == 0:
        print(f"No valid covariance matrices found for {sign}")
        return
    
    # Create cumulative histogram
    h_cov_cumulative = TH2D(f'h_covariance_cumulative_{sign}',
                           f'Cumulative Covariance Matrix ({sign});Parameter i;Parameter j',
                           max_size, -0.5, max_size - 0.5,
                           max_size, -0.5, max_size - 0.5)
    
    # Accumulate matrices
    for idx, row in fit_results.iterrows():
        cov = row['covariance']
        if cov is None:
            continue
            
        n_params = cov.shape[0]
        for i in range(n_params):
            for j in range(n_params):
                current_val = h_cov_cumulative.GetBinContent(i + 1, j + 1)
                h_cov_cumulative.SetBinContent(i + 1, j + 1, current_val + cov[i, j])
    
    out_dir.cd()
    h_cov_cumulative.Write()
    
    print(f"Accumulated {valid_cov_count} covariance matrices for {sign}")

def draw_means_sigmas(fit_results: pd.DataFrame,
                      out_dir: TDirectory,
                      sign: str) -> None:
    """
    Draw means and sigmas of the fit components
    Args:
        fit_results: DataFrame with fit results
        out_dir: Output directory
        sign: 'matter' or 'antimatter'
    """

    graph_mean = create_graph(
        fit_results, 'pt', 'mean', 'pt_err', 'mean_error',
        name=f'g_mean_{sign}',
        title=';#it{p}_{T} (GeV/c); #mu'
    )
    graph_sigma = create_graph(
        fit_results, 'pt', 'sigma', 'pt_err', 'sigma_error',
        name=f'g_sigma_{sign}',
        title=';#it{p}_{T} (GeV/c); #sigma'
    )

    graph_log_sigma = graph_sigma.Clone(f'g_log_sigma_{sign}')

    set_root_object(graph_mean, marker_style=20)
    set_root_object(graph_sigma, marker_style=20)
    set_root_object(graph_log_sigma, marker_style=20, title=';#it{p}_{T} (GeV/c); log(#sigma)')

    fit_mean = TF1('fit_mean', 'pol0', fit_results['pt'].min(), fit_results['pt'].max())
    graph_mean.Fit(fit_mean, 'Q')

    for i in range(graph_sigma.GetN()):
        y  = graph_sigma.GetY()[i]
        ey = graph_sigma.GetEY()[i]
        graph_log_sigma.SetPoint(i, graph_sigma.GetX()[i], np.log(y) if y > 0. else 0.)
        graph_log_sigma.SetPointError(i, graph_sigma.GetEX()[i], ey / y if y > 0. else 0.)

    fit_sigma = TF1('fit_sigma', 'log([0] + [1]/pow(x,[2]))',
        fit_results['pt'].min(), fit_results['pt'].max())

    fit_sigma.SetParameters(1.7e-3, 2.5e-3, 1.2)
    fit_sigma.SetParLimits(2, 0.3, 1.6)

    graph_log_sigma.Fit(fit_sigma, 'Q')

    out_dir.cd()
    graph_mean.Write(f'g_mean_{sign}')
    graph_sigma.Write(f'g_sigma_{sign}')
    graph_log_sigma.Write(f'g_log_sigma_{sign}')

        

def draw_results(h2_nsigma: TH2F, fit_results: pd.DataFrame, 
                out_dir: TDirectory, sign: str, 
                fit_config: FitConfig) -> None:
    """
    Draw the results of the fit
    
    Args:
        h2_nsigma: Reference histogram for pt binning
        fit_results: DataFrame with fit results
        out_dir: Output directory
        sign: 'matter' or 'antimatter'
        fit_config: Fit configuration
    """
    fit_results['pt_err'] = h2_nsigma.GetXaxis().GetBinWidth(1) / 2.
    
    # Primary integrals
    graph_primaries = create_graph(
        fit_results, 'pt', 'primaries_integral', 'pt_err', 'primaries_integral_error',
        name=f'g_primaries_{sign}', 
        title=';#it{p}_{T} (GeV/c); #it{N}_{p}'
    )
    graph_primaries_mc = create_graph(
        fit_results, 'pt', 'primaries_integral_mc', 'pt_err', 0., 
        name=f'g_primaries_mc_{sign}', 
        title=';#it{p}_{T} (GeV/c); #it{N}_{p}'
    )
    
    # Total integrals
    graph_total = create_graph(
        fit_results, 'pt', 'total_integral', 'pt_err', 0., 
        name=f'g_total_{sign}', 
        title=';#it{p}_{T} (GeV/c); #it{N}_{total}'
    )
    
    # Primary fraction
    graph_primary_fraction = create_graph(
        fit_results, 'pt', 'primary_fraction', 'pt_err', 'primary_fraction_error',
        name=f'g_primary_fraction_{sign}', 
        title=';#it{p}_{T} (GeV/c); f_{p}'
    )
    graph_primary_fraction_mc = create_graph(
        fit_results, 'pt', 'primary_fraction_mc', 'pt_err', 0., 
        name=f'g_primary_fraction_mc_{sign}', 
        title='MC;#it{p}_{T} (GeV/c); f_{p}'
    )
    graph_primaries_fraction_in_range = create_graph(
        fit_results, 'pt', 'primaries_fraction_in_range', 'pt_err', 0.,
        name=f'g_primaries_fraction_in_range_{sign}', 
        title=';#it{p}_{T} (GeV/c); f_{p} in DCA range'
    )
    
    # Set marker style
    for graph in [graph_primaries, graph_total, graph_primary_fraction, 
                  graph_primary_fraction_mc, graph_primaries_mc, graph_primaries_fraction_in_range]:
        set_root_object(graph, marker_style=20)

    # Write primary graphs
    out_dir.cd()
    graph_primaries.Write(f'g_primaries_integral_{sign}')
    graph_primaries_mc.Write(f'g_primaries_integral_mc_{sign}')
    graph_total.Write(f'g_total_integral_{sign}')
    graph_primary_fraction.Write(f'g_primary_fraction_{sign}')
    graph_primary_fraction_mc.Write(f'g_primary_fraction_mc_{sign}')
    graph_primaries_fraction_in_range.Write(f'g_primaries_fraction_in_range_{sign}')
    
    # Material fraction (if requested)
    if fit_config.store_fractions:
        graph_material_integral = create_graph(
            fit_results, 'pt', 'material_integral', 'pt_err', 'material_integral_error',
            name=f'g_material_{sign}',
            title=';#it{p}_{T} (GeV/c); #it{N}_{material}'
        )
        graph_material_fraction = create_graph(
            fit_results, 'pt', 'material_fraction', 'pt_err', 'material_fraction_error',
            name=f'g_material_fraction_{sign}',
            title=';#it{p}_{T} (GeV/c); f_{material}'
        )
        graph_material_fraction_mc = create_graph(
            fit_results, 'pt', 'material_fraction_mc', 'pt_err', 0.,
            name=f'g_material_fraction_mc_{sign}',
            title='MC;#it{p}_{T} (GeV/c); f_{material}'
        )
        
        set_root_object(graph_material_integral, marker_style=20)
        set_root_object(graph_material_fraction, marker_style=20)
        set_root_object(graph_material_fraction_mc, marker_style=20)
        
        graph_material_integral.Write(f'g_material_integral_{sign}')
        graph_material_fraction.Write(f'g_material_fraction_{sign}')
        graph_material_fraction_mc.Write(f'g_material_fraction_mc_{sign}')
        
        # Weak decay fraction
        graph_weak_decay_integral = create_graph(
            fit_results, 'pt', 'weak_decay_integral', 'pt_err', 'weak_decay_integral_error',
            name=f'g_weak_decay_{sign}',
            title=';#it{p}_{T} (GeV/c); #it{N}_{weak decay}'
        )
        graph_weak_decay_fraction = create_graph(
            fit_results, 'pt', 'weak_decay_fraction', 'pt_err', 'weak_decay_fraction_error',
            name=f'g_weak_decay_fraction_{sign}',
            title=';#it{p}_{T} (GeV/c); f_{weak decay}'
        )
        graph_weak_decay_fraction_mc = create_graph(
            fit_results, 'pt', 'weak_decay_fraction_mc', 'pt_err', 0.,
            name=f'g_weak_decay_fraction_mc_{sign}',
            title='MC;#it{p}_{T} (GeV/c); f_{weak decay}'
        )
        
        set_root_object(graph_weak_decay_integral, marker_style=20)
        set_root_object(graph_weak_decay_fraction, marker_style=20)
        set_root_object(graph_weak_decay_fraction_mc, marker_style=20)
        
        graph_weak_decay_integral.Write(f'g_weak_decay_integral_{sign}')
        graph_weak_decay_fraction.Write(f'g_weak_decay_fraction_{sign}')
        graph_weak_decay_fraction_mc.Write(f'g_weak_decay_fraction_mc_{sign}')
    
    # Store covariance matrices if requested
    if fit_config.store_covariance:
        accumulate_covariance_matrices(fit_results, out_dir, sign)

    draw_means_sigmas(fit_results, out_dir, sign)


def matter_antimatter_ratio(particle: str, outdir: TDirectory, 
                            efficiency_file: str, efficiency_hist_name: str) -> None:
    """
    Calculate and store matter/antimatter ratio
    
    Args:
        particle: Particle type
        outdir: Output directory
        efficiency_file: Path to efficiency file
        efficiency_hist_name: Name of efficiency histogram
    """
    g_matter = outdir.Get('g_primaries_integral_matter')
    g_antimatter = outdir.Get('g_primaries_integral_antimatter')

    g_matter_mc = outdir.Get('g_primaries_integral_mc_matter')
    g_antimatter_mc = outdir.Get('g_primaries_integral_mc_antimatter')

    h_efficiency = load_hist(efficiency_file, efficiency_hist_name)

    npoints = g_matter.GetN()
    g_ratio = TGraphErrors(npoints)
    g_ratio_mc = TGraphErrors(npoints)
    
    for ipoint in range(npoints):
        x_value = g_matter.GetPointX(ipoint)
        efficiency_matter = h_efficiency.GetBinContent(h_efficiency.FindBin(x_value))
        efficiency_antimatter = h_efficiency.GetBinContent(h_efficiency.FindBin(-x_value))

        # Data ratio
        denominator = g_matter.Eval(x_value) * efficiency_antimatter
        if denominator > 0.:
            y_value = (g_antimatter.Eval(x_value) * efficiency_matter / denominator)
            y_error = y_value * np.sqrt(
                1. / g_antimatter.Eval(x_value) + 1. / g_matter.Eval(x_value)
            ) if (g_matter.Eval(x_value) > 0. and g_antimatter.Eval(x_value) > 0.) else 0.
        else:
            y_value, y_error = 0., 0.
            
        g_ratio.SetPoint(ipoint, x_value, y_value)
        g_ratio.SetPointError(ipoint, g_matter.GetErrorX(ipoint), y_error)

        # MC ratio
        denominator_mc = g_matter_mc.Eval(x_value)
        if denominator_mc > 0.:
            y_value_mc = g_antimatter_mc.Eval(x_value) / denominator_mc
            y_error_mc = y_value_mc * np.sqrt(
                1. / g_antimatter_mc.Eval(x_value) + 1. / g_matter_mc.Eval(x_value)
            ) if (g_matter_mc.Eval(x_value) > 0. and g_antimatter_mc.Eval(x_value) > 0.) else 0.
        else:
            y_value_mc, y_error_mc = 0., 0.
            
        g_ratio_mc.SetPoint(ipoint, x_value, y_value_mc)
        g_ratio_mc.SetPointError(ipoint, g_matter_mc.GetErrorX(ipoint), y_error_mc)

    set_root_object(g_ratio, name='g_ratio', 
                   title=';#it{p}_{T} (GeV/#it{c}); Antimatter / Matter', 
                   marker_style=20)
    set_root_object(g_ratio_mc, name='g_ratio_mc', 
                   title='MC;#it{p}_{T} (GeV/#it{c}); Antimatter / Matter', 
                   marker_style=20)
    
    outdir.cd()
    g_ratio.Write()
    g_ratio_mc.Write()