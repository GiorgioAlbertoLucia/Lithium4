"""
PDF models for template fitting
"""

from typing import Dict, Tuple
from ROOT import RooRealVar, RooCrystalBall, RooGaussian, RooPolynomial, RooArgList


def build_crystal_ball_model(x: RooRealVar, suffix: str = '', **kwargs) -> Tuple[RooCrystalBall, Dict[str, RooRealVar]]:
    """
    Build a Crystal Ball PDF model
    
    Args:
        x: Observable variable
        suffix: Suffix for parameter names
        **kwargs: Override default parameter values
    
    Returns:
        Tuple of (pdf, parameters_dict)
    """

    nmax = 3. if (kwargs.get('flag', None) == 'primaries' and kwargs.get('rigidity', 0.) > 0.7) else 10.
    nmax = 30

    pars = {
        'mean': kwargs.get('mean', RooRealVar(f'mean{suffix}', f'#mu^{suffix}', 0., -0.015, 0.015)),
        'sigma': RooRealVar(f'sigma{suffix}', f'#sigma^{suffix}', 0.01, 1.e-3, 0.05),
        'aL': RooRealVar(f'aL{suffix}', f'#alpha_{{L}}^{suffix}', 1.2, 0.5, 10.),
        'nL': RooRealVar(f'nL{suffix}', f'n_{{L}}^{suffix}', 2., 1., nmax),
        'aR': RooRealVar(f'aR{suffix}', f'#alpha_{{R}}^{suffix}', 1.2, 0.5, 10.),
        'nR': RooRealVar(f'nR{suffix}', f'n_{{R}}^{suffix}', 2., 1., nmax)
    }
    
    pdf = RooCrystalBall(f'pdf_{suffix}', f'pdf_{suffix}', x, *pars.values())
    return pdf, pars


def build_gaussian_model(x: RooRealVar, suffix: str = '', **kwargs) -> Tuple[RooGaussian, Dict[str, RooRealVar]]:
    """
    Build a Gaussian PDF model
    
    Args:
        x: Observable variable
        suffix: Suffix for parameter names
        **kwargs: Override default parameter values
    
    Returns:
        Tuple of (pdf, parameters_dict)
    """
    pars = {
        'mean': kwargs.get('mean', RooRealVar(f'mean{suffix}', f'#mu^{suffix}', 0., -0.05, 0.05)),
        'sigma': RooRealVar(f'sigma{suffix}', f'#sigma^{suffix}', 0.01, 0.001, 0.05)
    }
    
    pdf = RooGaussian(f'pdf_{suffix}', f'pdf_{suffix}', x, pars['mean'], pars['sigma'])
    return pdf, pars


def build_pol0_model(x: RooRealVar, suffix: str = '') -> Tuple[RooPolynomial, Dict[str, RooRealVar]]:
    """
    Build a constant (pol0) PDF model
    
    Args:
        x: Observable variable
        suffix: Suffix for parameter names
    
    Returns:
        Tuple of (pdf, parameters_dict)
    """
    pars = {
        'c0': RooRealVar(f'c0{suffix}', f'c_{{0}}^{suffix}', 1.0, 0.0, 10.0)
    }
    
    pdf = RooPolynomial(f'pdf_{suffix}', f'pdf_{suffix}', x, RooArgList(pars['c0']), lowestOrder=0)
    return pdf, pars


def build_smearing_gaussian(x: RooRealVar, mean_val: float = 0.) -> Tuple[RooGaussian, Dict[str, RooRealVar]]:
    """
    Build a Gaussian smearing function
    
    Args:
        x: Observable variable
        mean_val: Initial value for mean (will be fixed)
    
    Returns:
        Tuple of (pdf, parameters_dict)
    """
    pars = {
        'mean': RooRealVar('mean_smearing', '#mu^{smearing}', mean_val, -0.005, 0.005),
        'sigma': RooRealVar('sigma_smearing', '#sigma^{smearing}', 0.003, 1.e-5, 0.005),
    }
    
    pars['mean'].setConstant(True)
    
    pdf = RooGaussian('gaussian_smearing', 'gaussian_smearing', x, *pars.values())
    return pdf, pars


def build_core_gaussian(x: RooRealVar) -> Tuple[RooGaussian, Dict[str, RooRealVar]]:
    """
    Build a core Gaussian for initial smearing determination
    
    Args:
        x: Observable variable
    
    Returns:
        Tuple of (pdf, parameters_dict)
    """
    pars = {
        'mean': RooRealVar('mean_core', '#mu^{core}', 0., -0.05, 0.05),
        'sigma': RooRealVar('sigma_core', '#sigma^{core}', 1.e-3, 1.e-5, 0.1),
    }
    
    pdf = RooGaussian('gaussian_core', 'gaussian_core', x, *pars.values())
    return pdf, pars