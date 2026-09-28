# core/config_loader.py
import yaml
from pathlib import Path
from torchic.core.histogram import HistLoadInfo

def load_yaml(path: str) -> dict:
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def build_hist_load_info_dict(section: dict) -> dict:
    """
    Turns:
      root_file: /path/to/file.root
      centralities: {'010': 'histA', '1030': 'histB', ...}
    into {'010': HistLoadInfo(path, 'histA'), ...}
    """
    root_file = section['root_file']
    return {cent: HistLoadInfo(root_file, name)
            for cent, name in section['centralities'].items()}

def build_hist_load_info_variations_radius(section: dict) -> dict:
    """
    Turns nested {cent: {variation: hist_name}} into
    {cent: {variation: HistLoadInfo(path, hist_name)}}
    """
    root_file = section['root_file']
    return {
        cent: {var: HistLoadInfo(root_file, name) for var, name in variations.items()}
        for cent, variations in section['variations_radius'].items()
    }
    
def build_hist_load_info_variations_model_params(section: dict) -> dict:
    """
    Turns nested {cent: {variation: hist_name}} into
    {cent: {variation: HistLoadInfo(path, hist_name)}}
    """
    root_file = section.get('root_file_variations_model_params', section['root_file'])
    return {
        cent: {var: HistLoadInfo(root_file, name) for var, name in variations.items()}
        for cent, variations in section['variations_model_params'].items()
    }