import argparse
import numpy as np
from typing import List 

from ROOT import TFile, TCanvas, TLine

from torchic.core.histogram import load_hist
from torchic.core.graph import load_graph
from torchic.utils.root import set_root_object, set_alice_global_style, init_legend, set_alice_frame_style
from torchic.utils.colors import get_color

INPUT_FILES_PURITY = {
    'He3': [
        '/home/galucia/DetectorCalibration/output/purity/LHC23_PbPb_pass5_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC24ar_pass3_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC25_PbPb_pass1_purity.root'
    ],
    'Had': [
        '/home/galucia/DetectorCalibration/output/purity/LHC23_PbPb_pass5_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC24ar_pass3_purity.root',
        '/home/galucia/DetectorCalibration/output/purity/LHC25_PbPb_pass1_purity.root'
    ]
}

INPUT_FILES_PRIMARY_FRACTION = {
    'He3': [
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC23_PbPb_pass5.root',
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC24ar_pass3.root',
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_LHC25_PbPb_pass1.root'
    ],
    'Had': [
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_trend_LHC23_pass5.root',
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_trend_LHC24ar_pass3.root',
        '/home/galucia/DetectorCalibration/output/dca/primary_fraction_trend_LHC25_pass1.root'
    ]
}

def _weighted_average_by_x(graphs: List, weights: List[float], precision: int = 6) -> List[tuple]:
    """
    Computes a weighted average of a list of TGraphs, matching points by x-value
    (rounded to `precision` decimals) instead of by index. Points missing from a
    given graph (e.g. failed fit) are simply skipped for that graph, and the
    weighted average is renormalized over whichever graphs do have that point.

    Returns a sorted list of (x, average_y) tuples for every x-value present in
    at least one graph.
    """
    graph_dicts = [
        {round(graph.GetX()[i], precision): graph.GetY()[i] for i in range(graph.GetN())}
        for graph in graphs
    ]
    all_x = sorted(set().union(*[d.keys() for d in graph_dicts]))

    averaged_points = []
    for x in all_x:
        y_values, point_weights = [], []
        for d, w in zip(graph_dicts, weights):
            if x in d:
                y_values.append(d[x])
                point_weights.append(w)
        if not y_values:
            continue  # no graph has this point, shouldn't normally happen
        averaged_points.append((x, np.average(y_values, weights=point_weights)))
    return averaged_points

def compute_average_purity(input_files: List[str], output_file: TFile, particle:str , detector:str,
                           weights: List[float], centrality:str=None) -> None:
    """
    Computes the average purity from the input files and saves it to the output file.

    Parameters:
    input_files (list): List of input file paths containing purity histograms.
    output_file (str): Path to the output file where the average purity will be saved.
    weights (list): List of weights for each input file.
    """
    
    outdir = output_file.mkdir(f"{particle}/{detector}")
    centrality_str = f"{centrality}/" if centrality else ""
        
    for charge in ['matter', 'antimatter']:
        graphs = []
        
        graph_name = f"{centrality_str}{particle}/{detector}/g_purity_{charge}"
        for file in input_files:
            print(f"Loading graph {graph_name} from {file}")
            graph = load_graph(file, graph_name)
            graph.SetName(f"{graph.GetName()}_{file.split('/')[-1].replace('.root', '')}")
            graphs.append(graph)
                
        g_purity_averaged = graphs[0].Clone(f"g_purity_{charge}_averaged")
        for ipoint in range(g_purity_averaged.GetN()):
            y_values = [graph.GetY()[ipoint] for graph in graphs]
            average_y = np.average(y_values, weights=weights)
            g_purity_averaged.SetPoint(ipoint, g_purity_averaged.GetX()[ipoint], average_y)
        
        outdir_charge = outdir.mkdir(f'{charge}_individuals')
        outdir_charge.cd()
        for graph in graphs:
            graph.Write()
        
        outdir.cd()
        g_purity_averaged.Write()
        
        graphs.clear()  # Clear the list for the next charge
        
def compute_average_primary_fraction(input_files: List[str], output_file: TFile, particle:str, weights: List[float],
                                     centrality:str,
                                     set_dummy_fraction_to_one: bool = False, set_dummy_fraction_to_zero: bool = False,
                                     set_lower_value: bool = False, set_higher_value: bool = False,
                                     match_ratio: bool = False) -> None:
    """
    Computes the average primary fraction from the input files and saves it to the output file.

    Parameters:
    input_files (list): List of input file paths containing primary fraction histograms.
    output_file (str): Path to the output file where the average primary fraction will be saved.
    weights (list): List of weights for each input file.
    """
    
    outdir = output_file.mkdir(f"{particle}")
        
    for charge in ['matter', 'antimatter']:
        
        graphs, ratio_graphs, weak_decay_graphs = [], [], []
        graph_name = f"{particle}/DCAxy_{centrality}/g_primary_fraction_{charge}"
        graph_name_weak_decay = f"{particle}/DCAxy_{centrality}/g_weak_decay_fraction_{charge}"
        ratio_graph_name = f"{particle}/DCAxy_{centrality}/g_ratio"
        for file in input_files:
            print(f"Loading graph {graph_name} from file {file}")
            graph = load_graph(file, graph_name)
            if graph_name == "He/DCAxy_centrality_0_10/g_primary_fraction_matter" and 'LHC24' in file:
                # for this specific graph, use the antimatter graph instead of the matter graph
                graph = load_graph(file, f"{particle}/DCAxy_{centrality}/g_primary_fraction_antimatter")
            
            graph.SetName(f"{graph_name}_{file.split('/')[-1].replace('.root', '')}")
            graphs.append(graph)

            print(f"Loading weak decay graph {graph_name_weak_decay} from file {file}")
            weak_decay_graph = load_graph(file, graph_name_weak_decay)
            weak_decay_graph.SetName(f"{weak_decay_graph.GetName()}_{file.split('/')[-1].replace('.root', '')}")
            weak_decay_graphs.append(weak_decay_graph)
            
            if particle == 'He' and charge == 'matter':
                print(f"Loading ratio graph {ratio_graph_name} from file {file}")
                ratio_graph = load_graph(file, ratio_graph_name)
                ratio_graph.SetName(f"{ratio_graph.GetName()}_{file.split('/')[-1].replace('.root', '')}")
                ratio_graphs.append(ratio_graph)

        primary_fraction_points = _weighted_average_by_x(graphs, weights)
        weak_decay_points = _weighted_average_by_x(weak_decay_graphs, weights)
        ratio_points = dict(_weighted_average_by_x(ratio_graphs, weights))

        g_primary_fraction_averaged = graphs[0].Clone(f"g_primary_fraction_{charge}_averaged")
        g_primary_fraction_averaged.Set(0)  # clear, refill using the x-matched points below
        for ipoint, (x, average_y) in enumerate(primary_fraction_points):
            average_ratio = ratio_points.get(x)

            if set_dummy_fraction_to_one and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = 1.0
            elif set_dummy_fraction_to_zero and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = 0.0
            elif set_lower_value and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = average_y + 0.3 if average_y < 0.7 else 1.0
            elif set_higher_value and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = average_y - 0.3 if average_y > 0.3 else 0.0
            elif match_ratio and particle == 'He' and charge == 'matter':
                if x < 1.8 and average_ratio is not None:
                    # Set the primary fraction so that Antimatter/Matter = 1
                    if average_ratio > 0:
                        average_y *= average_ratio
                    else:
                        average_y -= 0.3
                    average_y = min(max(average_y, 0.0), 1.0)  # Ensure the value is between 0 and 1
                    
            g_primary_fraction_averaged.SetPoint(ipoint, x, average_y)

        g_weak_decay_fraction_averaged = weak_decay_graphs[0].Clone(f"g_weak_decay_fraction_{charge}_averaged")
        g_weak_decay_fraction_averaged.Set(0)
        for ipoint, (x, average_y) in enumerate(weak_decay_points):
            if set_dummy_fraction_to_one and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = 1.0
            elif set_dummy_fraction_to_zero and particle == 'He' and charge == 'matter':
                if x < 1.8:
                    average_y = 0.0
            g_weak_decay_fraction_averaged.SetPoint(ipoint, x, average_y)
        
        outdir_charge = outdir.mkdir(f'{charge}_individuals')
        outdir_charge.cd()
        for graph in graphs:
            graph.Write()
        for weak_decay_graph in weak_decay_graphs:
            weak_decay_graph.Write()
        
        outdir.cd()
        g_primary_fraction_averaged.Write()
        g_weak_decay_fraction_averaged.Write()
        graphs.clear()  # Clear the list for the next charge
        
def compute_routine(args: argparse.Namespace) -> None:
    
    output_file_path = 'output/average_purity_and_primary_fraction'
    if args.set_dummy_fraction_to_one:
        output_file_path += '_dummy_fraction_to_one'
    elif args.set_dummy_fraction_to_zero:
        output_file_path += '_dummy_fraction_to_zero'
    elif args.set_lower_value:
        output_file_path += '_set_lower_value'
    elif args.set_higher_value:
        output_file_path += '_set_higher_value'
    elif args.match_ratio:
        output_file_path += '_match_ratio'
    output_file_path += '.root'
    output_file = TFile(output_file_path, 'RECREATE')
    
    input_files_weights = [
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC23_PbPb_pass5_hadronpid_same.root',
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC24ar_pass3_hadronpid_same.root',
        '/home/galucia/Lithium4/preparation/output/PbPb/LHC25_PbPb_pass1_hadronpid_same.root',
    ]
    h_weights = [load_hist(file, 'QA/hKstar') for file in input_files_weights]
    weights = [h.Integral(h.FindBin(0), h.FindBin(0.4)) for h in h_weights]
    
    centrality_classes = ['centrality_0_10', 'centrality_10_50']
    for centrality in centrality_classes:
        
        outdir_centrality = output_file.mkdir(f"{centrality}")
    
        outdir_purity = outdir_centrality.mkdir(f"purity")
        outdir_primary_fraction = outdir_centrality.mkdir(f"primary_fraction")
            
        for particle in ['He3', 'Had']:
        #for particle in ['Had']:
            particle_primary_fraction = 'He' if particle == 'He3' else 'Pr'
            compute_average_primary_fraction(INPUT_FILES_PRIMARY_FRACTION[particle], outdir_primary_fraction, particle_primary_fraction, weights,
                                             centrality=centrality,
                                            set_dummy_fraction_to_one=args.set_dummy_fraction_to_one,
                                            set_dummy_fraction_to_zero=args.set_dummy_fraction_to_zero,
                                            set_lower_value=args.set_lower_value,
                                            set_higher_value=args.set_higher_value, 
                                            match_ratio=args.match_ratio)
            
            for detector in ['TPC', 'TOF']:
                if particle == 'He3' and detector == 'TOF':
                    continue  # Skip He3 TOF as it is not relevant
                
                compute_average_purity(INPUT_FILES_PURITY[particle], outdir_purity, particle, detector, weights,
                                       centrality=centrality)
    
def compare_routine(args: argparse.Namespace) -> None:
    
    input_files_lambda = [
        '/home/galucia/Lithium4/calibration/output/lambda_parameters_dummy_fraction_to_one.root', # primary fraction boosted to one
        '/home/galucia/Lithium4/calibration/output/lambda_parameters_dummy_fraction_to_zero.root',  # primary fraction set to zero
        '/home/galucia/Lithium4/calibration/output/lambda_parameters_set_lower_value.root',  # lambda parameters with lower value
        '/home/galucia/Lithium4/calibration/output/lambda_parameters_set_higher_value.root',  # lambda parameters with higher value
        '/home/galucia/Lithium4/calibration/output/lambda_parameters.root'  # nominal lambda parameters
    ]   
    
    outfile = TFile('output/compare_average_purity_and_primary_fraction.root', 'RECREATE')
    h_lambda_fraction_to_one = load_hist(input_files_lambda[0], 'Both/hLambdaParameters')
    h_lambda_fraction_to_zero = load_hist(input_files_lambda[1], 'Both/hLambdaParameters')
    h_lambda_set_lower_value = load_hist(input_files_lambda[2], 'Both/hLambdaParameters')
    h_lambda_set_higher_value = load_hist(input_files_lambda[3], 'Both/hLambdaParameters')
    h_lambda_nominal = load_hist(input_files_lambda[4], 'Both/hLambdaParameters')
    
    c_compare = TCanvas('c_compare', 'c_compare', 800, 600) 
    hframe = c_compare.DrawFrame(0, 0.5, 0.4, 0.9, ';#it{k}* (GeV/#it{c});#it{#lambda}')
    set_root_object(h_lambda_nominal, line_color=get_color(0), line_width=2, title='Nominal #it{#lambda} parameters')
    set_root_object(h_lambda_fraction_to_one, line_color=get_color(1), line_width=2, title='#it{f}_{primary}(^{3}He) = 1, #it{p}_{T} < 1.8 GeV/#it{c}')
    set_root_object(h_lambda_fraction_to_zero, line_color=get_color(2), line_width=2, title='#it{f}_{primary}(^{3}He) = 0, #it{p}_{T} < 1.8 GeV/#it{c}')
    set_root_object(h_lambda_set_lower_value, line_color=get_color(3), line_width=2, title='#it{f}_{primary}(^{3}He) = lower value, #it{p}_{T} < 1.8 GeV/#it{c}')
    set_root_object(h_lambda_set_higher_value, line_color=get_color(4), line_width=2, title='#it{f}_{primary}(^{3}He) = higher value, #it{p}_{T} < 1.8 GeV/#it{c}')
    h_lambda_fraction_to_one.Draw('hist same')
    h_lambda_fraction_to_zero.Draw('hist same')
    h_lambda_nominal.Draw('hist same')
    h_lambda_set_lower_value.Draw('hist same')
    h_lambda_set_higher_value.Draw('hist same')

    lines = [
        TLine(0, 0.7, 0.4, 0.7),
        TLine(0, 0.8, 0.4, 0.8),
    ]
    for line in lines:
        set_root_object(line, line_color=13, line_style=2, line_width=2)
        #line.Draw('same')
    
    legend = init_legend(0.4, 0.22, 0.85, 0.45)
    legend.AddEntry(h_lambda_nominal, h_lambda_nominal.GetTitle(), 'l')
    legend.AddEntry(h_lambda_fraction_to_one, h_lambda_fraction_to_one.GetTitle(), 'l')
    legend.AddEntry(h_lambda_fraction_to_zero, h_lambda_fraction_to_zero.GetTitle(), 'l')
    legend.AddEntry(h_lambda_set_lower_value, h_lambda_set_lower_value.GetTitle(), 'l')
    legend.AddEntry(h_lambda_set_higher_value, h_lambda_set_higher_value.GetTitle(), 'l')
    legend.Draw()
    
    h_ratio_one_to_nominal = h_lambda_fraction_to_one.Clone('h_ratio')
    h_ratio_one_to_nominal.Divide(h_lambda_nominal)
    h_ratio_zero_to_nominal = h_lambda_fraction_to_zero.Clone('h_ratio')
    h_ratio_zero_to_nominal.Divide(h_lambda_nominal)
    h_ratio_set_lower_to_nominal = h_lambda_set_lower_value.Clone('h_ratio')
    h_ratio_set_lower_to_nominal.Divide(h_lambda_nominal)
    h_ratio_set_higher_to_nominal = h_lambda_set_higher_value.Clone('h_ratio')
    h_ratio_set_higher_to_nominal.Divide(h_lambda_nominal)
    set_root_object(h_ratio_one_to_nominal, line_color=get_color(2), line_width=2, title='Ratio: #it{f}_{primary}(^{3}He) = 1 / nominal')
    set_root_object(h_ratio_zero_to_nominal, line_color=get_color(3), line_width=2, title='Ratio: #it{f}_{primary}(^{3}He) = 0 / nominal')
    set_root_object(h_ratio_set_lower_to_nominal, line_color=get_color(4), line_width=2, title='Ratio: #it{f}_{primary}(^{3}He) = lower value / nominal')
    set_root_object(h_ratio_set_higher_to_nominal, line_color=get_color(5), line_width=2, title='Ratio: #it{f}_{primary}(^{3}He) = higher value / nominal')
    
    c_compare_ratio = TCanvas('c_compare_ratio', 'c_compare_ratio', 800, 600)
    hframe_ratio = c_compare_ratio.DrawFrame(0, 0.5, 0.4, 1.3, ';#it{k}* (GeV/#it{c});Ratio')
    line_at_one = TLine(0, 1, 0.4, 1)
    set_root_object(line_at_one, line_color=13, line_style=2, line_width=2)
    line_at_one.Draw('same')
    h_ratio_one_to_nominal.Draw('hist same')
    h_ratio_zero_to_nominal.Draw('hist same')
    h_ratio_set_lower_to_nominal.Draw('hist same')
    h_ratio_set_higher_to_nominal.Draw('hist same')
    
    legend_ratio = init_legend(0.4, 0.2, 0.85, 0.4)
    legend_ratio.AddEntry(h_ratio_one_to_nominal, h_ratio_one_to_nominal.GetTitle(), 'l')
    legend_ratio.AddEntry(h_ratio_zero_to_nominal, h_ratio_zero_to_nominal.GetTitle(), 'l')
    legend_ratio.AddEntry(h_ratio_set_lower_to_nominal, h_ratio_set_lower_to_nominal.GetTitle(), 'l')
    legend_ratio.AddEntry(h_ratio_set_higher_to_nominal, h_ratio_set_higher_to_nominal.GetTitle(), 'l')
    legend_ratio.Draw()
    
    outfile.cd()
    c_compare.Write()
    c_compare_ratio.Write()
    
    outfile.Close()

def compare_routine_match_ratio(args: argparse.Namespace) -> None:
    
    input_files_lambda = [
        '/home/galucia/Lithium4/calibration/output/lambda_parameters_match_ratio.root',  # primary fraction set to match ratio
        '/home/galucia/Lithium4/calibration/output/lambda_parameters.root'  # nominal lambda parameters
    ]   
    
    outfile = TFile('output/compare_average_purity_and_primary_fraction.root', 'RECREATE')
    for centrality in ['centrality_0_10', 'centrality_10_50']:
        
        h_lambda_match_ratio = load_hist(input_files_lambda[0], f'{centrality}/Both/hLambdaParameters')
        h_lambda_nominal = load_hist(input_files_lambda[1], f'{centrality}/Both/hLambdaParameters')
        
        c_compare = TCanvas(f'c_compare_{centrality}', 'c_compare', 800, 600) 
        hframe = c_compare.DrawFrame(0, 0.5, 0.4, 0.9, ';#it{k}* (GeV/#it{c});#it{#lambda}')
        set_root_object(h_lambda_nominal, line_color=get_color(0), line_width=2, title='Nominal #it{#lambda} parameters')
        set_root_object(h_lambda_match_ratio, line_color=get_color(1), line_width=2, title='^{3}He = ^{3}#bar{He}, #it{p}_{T} < 1.8 GeV/#it{c}')
        h_lambda_nominal.Draw('hist same')
        h_lambda_match_ratio.Draw('hist same')

        lines = [
            TLine(0, 0.7, 0.4, 0.7),
            TLine(0, 0.8, 0.4, 0.8),
        ]
        for line in lines:
            set_root_object(line, line_color=13, line_style=2, line_width=2)
            #line.Draw('same')
        
        legend = init_legend(0.4, 0.22, 0.85, 0.45)
        legend.AddEntry(h_lambda_nominal, h_lambda_nominal.GetTitle(), 'l')
        legend.AddEntry(h_lambda_match_ratio, h_lambda_match_ratio.GetTitle(), 'l')
        legend.Draw()
        
        h_ratio_match_to_nominal = h_lambda_match_ratio.Clone('h_ratio')
        h_ratio_match_to_nominal.Divide(h_lambda_nominal)
        set_root_object(h_ratio_match_to_nominal, line_color=get_color(2), line_width=2, title='Ratio: #it{f}_{primary}(^{3}He) = match ratio / nominal')
        
        c_compare_ratio = TCanvas(f'c_compare_ratio_{centrality}', 'c_compare_ratio', 800, 600)
        hframe_ratio = c_compare_ratio.DrawFrame(0, 0.5, 0.4, 1.3, ';#it{k}* (GeV/#it{c});Ratio')
        line_at_one = TLine(0, 1, 0.4, 1)
        set_root_object(line_at_one, line_color=13, line_style=2, line_width=2)
        line_at_one.Draw('same')
        h_ratio_match_to_nominal.Draw('hist same')
        
        legend_ratio = init_legend(0.4, 0.2, 0.85, 0.4)
        legend_ratio.AddEntry(h_ratio_match_to_nominal, h_ratio_match_to_nominal.GetTitle(), 'l')
        legend_ratio.Draw()
        
        outfile.cd()
        c_compare.Write()
        c_compare_ratio.Write()
    
    outfile.Close() 
   
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Compute average purity and primary fraction from multiple ROOT files.')
    parser.add_argument('--mode', type=str, default='compute', choices=['compute', 'compare'], help='Mode of operation: compute or compare.')
    parser.add_argument('--set_dummy_fraction_to_one', action='store_true', help='Set dummy fraction to one for specified particles and charge.')
    parser.add_argument('--set_dummy_fraction_to_zero', action='store_true', help='Set dummy fraction to zero for specified particles and charge.')
    parser.add_argument('--set_lower_value', action='store_true', help='Set lower value for specified particles and charge.')
    parser.add_argument('--set_higher_value', action='store_true', help='Set higher value for specified particles and charge.')
    parser.add_argument('--match_ratio', action='store_true', help='Set primary fraction so that Antimatter/Matter = 1.')
    args = parser.parse_args()
    
    set_alice_global_style()
    
    if args.mode == 'compute':
        compute_routine(args)
    elif args.mode == 'compare':
        if args.match_ratio:
            compare_routine_match_ratio(args)
        else:
            compare_routine(args)
    