from typing import List 

from ROOT import TFile

from torchic.core.histogram import load_hist
from torchic.core.graph import load_graph

def compute_average_purity(input_files: List[str], output_file: TFile, particle:str , detector:str,
                           weights: List[float]) -> None:
    """
    Computes the average purity from the input files and saves it to the output file.

    Parameters:
    input_files (list): List of input file paths containing purity histograms.
    output_file (str): Path to the output file where the average purity will be saved.
    weights (list): List of weights for each input file.
    """
    
    outdir = output_file.mkdir(f"{particle}/{detector}")
        
    for charge in ['matter', 'antimatter']:
        graphs = []
        
        graph_name = f"{particle}/{detector}/g_purity_{charge}"
        for file in input_files:
            graph = load_graph(file, graph_name)
            graph.SetName(f"{graph.GetName()}_{file.split('/')[-1].replace('.root', '')}")
            graphs.append(graph)
                
        g_purity_averaged = graphs[0].Clone(f"g_purity_{charge}_averaged")
        for ipoint in range(g_purity_averaged.GetN()):
            y_values = [graph.GetY()[ipoint] for graph in graphs]
            average_y = sum(y * w for y, w in zip(y_values, weights)) / sum(weights)
            g_purity_averaged.SetPoint(ipoint, g_purity_averaged.GetX()[ipoint], average_y)
        
        outdir_charge = outdir.mkdir(f'{charge}_individuals')
        outdir_charge.cd()
        for graph in graphs:
            graph.Write()
        
        outdir.cd()
        g_purity_averaged.Write()
        
        graphs.clear()  # Clear the list for the next charge
        
def compute_average_primary_fraction(input_files: List[str], output_file: TFile, particle:str, weights: List[float],
                                     set_dummy_fraction_to_one: bool = False, set_dummy_fraction_to_zero: bool = False) -> None:
    """
    Computes the average primary fraction from the input files and saves it to the output file.

    Parameters:
    input_files (list): List of input file paths containing primary fraction histograms.
    output_file (str): Path to the output file where the average primary fraction will be saved.
    weights (list): List of weights for each input file.
    """
    
    outdir = output_file.mkdir(f"{particle}")
        
    for charge in ['matter', 'antimatter']:
        graphs = []
        
        graph_name = f"{particle}/DCAxy/g_primary_fraction_{charge}"
        for file in input_files:
            graph = load_graph(file, graph_name)
            graph.SetName(f"{graph.GetName()}_{file.split('/')[-1].replace('.root', '')}")
            graphs.append(graph)
                
        g_primary_fraction_averaged = graphs[0].Clone(f"g_primary_fraction_{charge}_averaged")
        for ipoint in range(g_primary_fraction_averaged.GetN()):
            y_values = [graph.GetY()[ipoint] for graph in graphs]
            average_y = sum(y * w for y, w in zip(y_values, weights)) / sum(weights)
            
            if set_dummy_fraction_to_one and particle == 'He' and charge == 'matter':
                if graphs[0].GetX()[ipoint] < 1.8:
                    average_y = 1.0
            elif set_dummy_fraction_to_zero and particle == 'He' and charge == 'matter':
                if graphs[0].GetX()[ipoint] < 1.8:
                    average_y = 0.0
            g_primary_fraction_averaged.SetPoint(ipoint, g_primary_fraction_averaged.GetX()[ipoint], average_y)
        
        outdir_charge = outdir.mkdir(f'{charge}_individuals')
        outdir_charge.cd()
        for graph in graphs:
            graph.Write()
        
        outdir.cd()
        g_primary_fraction_averaged.Write()
        
        graphs.clear()  # Clear the list for the next charge
        
    
if __name__ == "__main__":
    
    input_files_purity = [
        '/home/galucia/DetectorCalibration/output/purity/LHC23_PbPb_pass5_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC24ar_pass3_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC25_PbPb_pass1_purity.root'
    ]
    input_files_primary_fraction = [
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC23_PbPb_pass5.root',
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC23_PbPb_pass5.root', # dummy for now
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC23_PbPb_pass5.root' # dummy for now
    ]
    output_file = TFile('output/average_purity_and_primary_fraction.root', 'RECREATE')
    
    input_files_weights = [
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC23_PbPb_pass5_hadronpid_same.root',
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC24ar_pass3_hadronpid_same.root',
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC25_PbPb_pass1_hadronpid_same.root',
    ]
    h_weights = [load_hist(file, 'QA/hKstar') for file in input_files_weights]
    weights = [h.Integral(h.FindBin(0), h.FindBin(0.4)) for h in h_weights]
    
    outdir_purity = output_file.mkdir(f"purity")
    outdir_primary_fraction = output_file.mkdir(f"primary_fraction")
        
    for particle in ['He3', 'Had']:
        particle_primary_fraction = 'He' if particle == 'He3' else 'Pr'
        compute_average_primary_fraction(input_files_primary_fraction, outdir_primary_fraction, particle_primary_fraction, weights)
        
        for detector in ['TPC', 'TOF']:
            if particle == 'He3' and detector == 'TOF':
                continue  # Skip He3 TOF as it is not relevant
            
            compute_average_purity(input_files_purity, outdir_purity, particle, detector, weights)