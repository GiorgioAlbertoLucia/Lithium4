"""
Main workflow for primary fraction estimation
"""

import pandas as pd
from ROOT import TFile, TCanvas

from torchic.core.histogram import load_hist
from torchic.utils.terminal_colors import TerminalColors as tc

from core.primary_fraction_config import AnalysisConfig
from core.primary_fraction_models import build_crystal_ball_model, build_core_gaussian
from core.primary_fraction_templates import prepare_hist_material_template
from core.primary_fraction_fitting import fit_slice
from core.primary_fraction_results import draw_results, matter_antimatter_ratio


def template_fitting_routine(h2_data, h2_mc, outdir, outpdf:str, 
                             particle, particle_config, fit_config, 
                             h_fraction_he3=None):
    """
    Main template fitting routine for a single particle
    
    Args:
        h2_data: Data histogram (2D: pt vs DCA)
        h2_mc: Dictionary of MC histograms by flag
        outdir: Output directory
        particle: Particle name ('He' or 'Pr')
        particle_config: Particle configuration
        fit_config: Fit configuration
        h_fraction_he3: He3 from hypertriton fraction (optional)
    """
    
    canvas = TCanvas('c', '')

    pt_min_matter, pt_max_matter = particle_config.pt_range
    pt_min_antimatter, pt_max_antimatter = -pt_max_matter, -pt_min_matter
    dca_min, dca_max = particle_config.dca_range
    
    from ROOT import RooRealVar
    dca = RooRealVar('dca', 'DCA_{xy}', dca_min, dca_max, 'cm')

    # Build core Gaussian for smearing determination
    gaussian_core, gaussian_core_pars = build_core_gaussian(dca)
    dca.setRange('core', *particle_config.dca_core_range)

    for sign in ['matter', 'antimatter']:

        fit_results = None
        canvas.Print(f'{outpdf.split('.')[0]}_{sign}.pdf(')

        if sign == 'antimatter':
            pt_min, pt_max = pt_min_antimatter, pt_max_antimatter
            # Remove material for antimatter
            h2_mc_sign = {k: v for k, v in h2_mc.items() if k != 'material'}
        else:
            pt_min, pt_max = pt_min_matter, pt_max_matter
            h2_mc_sign = h2_mc.copy()

        pt_bin_start = h2_data.GetXaxis().FindBin(pt_min)
        pt_bin_end = h2_data.GetXaxis().FindBin(pt_max) + 1 if pt_max < 0 else h2_data.GetXaxis().FindBin(pt_max)

        for pt_bin in range(pt_bin_start, pt_bin_end):
            
            pt = h2_data.GetXaxis().GetBinCenter(pt_bin)
            rigidity = abs(pt) / 2. if particle == 'He' else abs(pt)
            
            # Build models for this bin
            pdfs, pdf_params = {}, {}
            for flag in h2_mc_sign.keys():
                pdfs[flag], pdf_params[flag] = build_crystal_ball_model(
                    dca, suffix=flag, flag=flag, rigidity=rigidity #, mean=gaussian_core_pars['mean']
                )

            h_dca = h2_data.ProjectionY(f'h_dca_{pt:.2f}', pt_bin, pt_bin, 'e')
            
            if h_dca.GetEntries() < fit_config.min_entries_for_fit:
                continue
            
            bin_outdir = outdir.mkdir(f'pt_{pt:.2f}')

            # Perform fit
            print(tc.RED+f"\nFitting {particle} {sign} in pT bin {pt:.2f} GeV/c..."+tc.RESET)
            dca_frame, ifit_results = fit_slice(
                h2_data, h2_mc_sign, pdfs, pdf_params,
                gaussian_core, gaussian_core_pars,
                dca, pt_bin, particle, bin_outdir, 
                fit_config, h_fraction_he3
            )
            
            if dca_frame is None:
                continue
            
            # Store results
            if fit_results is None:
                fit_results = pd.DataFrame.from_dict([ifit_results])
            else:
                fit_results = pd.concat(
                    [fit_results, pd.DataFrame.from_dict([ifit_results])], 
                    ignore_index=True
                )

            # Draw frame
            dca_frame.Draw()
            canvas.SetLogy()

            bin_outdir.cd()
            canvas.Write(f'dca_frame_{sign}_{pt:.2f}')
            canvas.Print(f'{outpdf.split('.')[0]}_{sign}.pdf')
            canvas.Clear()

            pdfs.clear()
            pdf_params.clear()

        # Draw and save results
        if fit_results is not None:
            draw_results(h2_data, fit_results, outdir, sign, fit_config)
        
        canvas.Print(f'{outpdf.split('.')[0]}_{sign}.pdf)')

    del canvas


def main(config: AnalysisConfig = None):
    """Main analysis function"""
    
    # Load configuration
    if config is None:
        config = AnalysisConfig()
    
    outfile = TFile.Open(config.paths.output_file, 'RECREATE')
    
    # Load He3 fraction from hypertriton
    h_fraction_he3_from_h3l = load_hist(
        config.paths.he3_fraction_file,
        config.paths.he3_fraction_histname
    )
    
    for particle in ['He', 'Pr']:
        outdir_particle = outfile.mkdir(particle)

        for direction in ['xy', 'z']:
            print(f"\n{'='*60}")
            print(f"Processing {particle}")
            print(f"{'='*60}\n")
            
            outdir = outdir_particle.mkdir(f'DCA{direction}')
            particle_config = config.get_particle_config(particle)
            outpdf = config.paths.output_pdf.replace('.pdf', f'_{particle}_{direction}.pdf')
            
            # Load MC templates
            h2_mc = {}
            for flag_name, flag_suffix in particle_config.mc_flags.items():
                    
                if flag_name == 'material' and particle == 'He':
                    if config.fits[particle].use_material_template:
                        
                        # Subtraction method to prepare material template
                        efficiency_path = f"{particle_config.efficiency_hist_path}/{particle_config.efficiency_hist_name}"
                        efficiency_material_path = f"{particle_config.efficiency_hist_path}/{particle_config.efficiency_material_hist_name}"
                        
                        #h2_material = prepare_hist_material_template(
                        #    config.paths.data_input_file,
                        #    f'{particle}/h2DCA{direction}Pt{particle}',
                        #    outdir,
                        #    config.paths.efficiency_file,
                        #    efficiency_path,
                        #    efficiency_material_path,
                        #    particle
                        #)
                        
                        #h2_material = load_hist(config.paths.mc_input_file, f'De_as_He/h2DCA{direction}PtDe_as_He{flag_suffix}')
                        h2_material = load_hist(config.paths.mc_input_file, f'Pr_as_He/h2DCA{direction}PtPr_as_He{flag_suffix}')
                        
                        h2_mc[flag_name] = h2_material

                    elif config.fits[particle].use_material_template_from_mc:
                        hist_path = f'{particle}/h2DCA{direction}Pt{particle}{flag_suffix}'
                        h2_mc[flag_name] = load_hist(config.paths.mc_input_file, hist_path)

                    # If not using template, will be created as Gaussian/pol0 in fitting
                                    
                else:
                    hist_path = f'{particle}/h2DCA{direction}Pt{particle}{flag_suffix}'
                    h2_mc[flag_name] = load_hist(config.paths.mc_input_file, hist_path)

                outdir.cd()
                h2_mc[flag_name].Write(f'h2DCA{direction}Pt{particle}_{flag_name}_mc')

            # Load data
            h2_data = load_hist(
                config.paths.data_input_file, 
                f'{particle}/h2DCA{direction}Pt{particle}'
            )
            outdir.cd()
            h2_data.Write(f'h2DCA{direction}Pt{particle}_data')

            # Run fitting
            h_fraction_he3 = h_fraction_he3_from_h3l if particle == 'He' else None
            template_fitting_routine(
                h2_data, h2_mc, outdir, outpdf, particle,
                particle_config, config.fits[particle], h_fraction_he3
            )
            
            # Calculate matter/antimatter ratio
            efficiency_path = f"{particle_config.efficiency_hist_path}/{particle_config.efficiency_primary_hist_name}"
            matter_antimatter_ratio(
                particle, outdir,
                config.paths.efficiency_file,
                efficiency_path
            )

    outfile.Close()
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)


if __name__ == '__main__':
    
    config = AnalysisConfig()
    #config.fits['He'].use_material_template_from_mc = True
    #config.fits['He'].initialise_from_mc_yields = True
    #config.fits['He'].use_convolutional_smearing = False
    #config.fits['Pr'].initialise_from_mc_yields = True
    #config.fits['Pr'].use_convolutional_smearing = False
    #config.fits['He'].is_mc = True
    #config.fits['Pr'].is_mc = True

    config.fits['He'].add_pol0_to_mc_template_for_material = True
    config.fits['He'].add_pol0_to_mc_template_for_primaries = False
    config.fits['Pr'].add_pol0_to_mc_template_for_material = True

    #config.fits['He'].max_pt_material_template = 2.0

    config.paths.data_input_file = 'output/dca/dca_data_template_smaller_tolerance.root'
    
    config.paths.output_file = 'output/dca/primary_fraction_results_new_smaller_tolerance_pr.root'
    #config.paths.output_file = 'output/primary_fraction_results_mc.root'
    #config.fits['He'].material_background_type = 'pol0'
    config.fits['He'].use_material_template = True

    config.particles['He'].mc_flags = {'primaries': '_IsPhysicalPrimary',
                                       'weak_decay': '_IsSecondaryFromWeakDecay',
                                       'material': '_IsSecondaryFromMaterial',
                                       }
    
    main(config=config)