"""
Configuration module for primary fraction estimation
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional

@dataclass
class ParticleConfig:
    """Configuration for a specific particle type"""
    pt_range: Tuple[float, float]
    dca_range: Tuple[float, float]
    dca_core_range: Tuple[float, float]
    mc_flags: Dict[str, str]
    
    # Efficiency histogram paths
    efficiency_hist_path: str
    efficiency_hist_name: str
    efficiency_material_hist_name: str
    efficiency_primary_hist_name: str

    pt_rebin: int = 1

@dataclass
class FitConfig:
    """Configuration for fitting procedure"""
    # Template options
    use_material_template: bool = True  # If False, use pol0
    use_material_template_from_mc: bool = False
    use_weak_decay_template: bool = True
    use_convolutional_smearing: bool = True
    is_mc : bool = False

    add_pol0_to_mc_template_for_material: bool = False
    add_pol0_to_mc_template_for_primaries: bool = False
    
    # Normalization options
    initialise_from_mc_yields: bool = False
    fix_material_yield: bool = False
    fix_weak_decay_yield: bool = False
    
    # Background model
    material_background_type: str = 'template'  # 'template', 'gaussian', or 'pol0'
    
    # Fitting options
    min_entries_for_fit: int = 100
    min_entries_for_template: int = 100

    max_pt_material_template: float = 999  # If set, only use MC templates up to this pT for material background
    
    # Storage options
    store_covariance: bool = True
    store_fractions: bool = True


@dataclass
class PathConfig:
    """Configuration for input/output paths"""
    data_input_file: str
    data_input_file_he: str
    mc_input_file: str
    output_file: str
    output_pdf: str
    efficiency_file: str
    he3_fraction_file: str
    he3_fraction_histname: str


class AnalysisConfig:
    """Main configuration class"""
    
    def __init__(self):
        self.particles = {
            'He': ParticleConfig(
                pt_rebin=4,
                pt_range=(1., 5.0),
                dca_range=(-0.1, 0.1),
                dca_core_range=(-0.004, 0.006),
                mc_flags={
                    'primaries': '_IsPhysicalPrimary',
                    'weak_decay': '_IsSecondaryFromWeakDecay',
                    'material': '_IsSecondaryFromMaterial',
                },
                efficiency_hist_path='He',
                efficiency_hist_name='efficiency',
                efficiency_material_hist_name='h_efficiency_he3_material',
                efficiency_primary_hist_name='h_efficiency_he3_prim'
            ),
            'Pr': ParticleConfig(
                pt_range=(0.4, 4.0),
                dca_range=(-0.1, 0.1),
                dca_core_range=(-0.004, 0.006),
                mc_flags={
                    'primaries': '_IsPhysicalPrimary',
                    'weak_decay': '_IsFromLambda0',
                    'material': '_IsSecondaryFromMaterial',
                },
                efficiency_hist_path='Pr',  # Note: This might need adjustment
                efficiency_hist_name='h_efficiency_he3',
                efficiency_material_hist_name='h_efficiency_he3_material',
                efficiency_primary_hist_name='h_efficiency_he3_prim'
            )
        }
        
        self.paths = PathConfig(
            data_input_file='output/dca/dca_data_template.root',
            data_input_file_he='output/dca/dca_data_template_nuclei_spectra.root',
            #data_input_file='output/dca/dca_mc_template_check.root',
            mc_input_file='output/dca/dca_mc_template.root',
            output_file='output/dca/primary_fraction_estimation.root',
            output_pdf='output/dca/primary_fraction_estimation.pdf',
            efficiency_file='/home/galucia/Efficiency/nucleiQC/output/efficiency/He/efficiency_LHC25g11.root',
            he3_fraction_file='/home/galucia/Lithium4/phase_space_studies/output/he3_from_hypertriton.root',
            he3_fraction_histname='h_fraction'
        )
        
        self.fits = {'He':  FitConfig(),
                     'Pr': FitConfig()}
    
    def get_particle_config(self, particle: str) -> ParticleConfig:
        """Get configuration for a specific particle"""
        if particle not in self.particles:
            raise ValueError(f"Unknown particle: {particle}")
        return self.particles[particle]