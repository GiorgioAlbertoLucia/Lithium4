import yaml

import os
import sys
sys.path.append('..')
from torchic import Plotter
from torchic.utils.root import set_alice_global_style

import argparse

if __name__ == '__main__':
    
    set_alice_global_style()

    #input_file = 'config/lambda_parameters.yml'
    #input_file = 'config/same_mixed_event_comparisons.yml'
    #input_file = 'config/purity.yml' 
    input_file = 'config/correlation_functions.yml' 
    #input_file = 'config/checks_all_vs_cbt.yml'
    #input_file = 'config/check_mixed.yml'
    #input_file = 'config/checks_before_forum.yml'

    parser = argparse.ArgumentParser(description='Plotting script for ROOT files using YAML configuration.')
    parser.add_argument('--config', type=str, default=input_file, help='Path to the YAML configuration file.')

    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    plotter = Plotter(config['outPath'])

    for plot in config['plots']:

        plotter.create_canvas(plot['axisSpecs'], **plot['canvas'])
        plotter.create_multigraph(plot['axisSpecs'])
        if 'legends' in plot.keys():
            for legend in plot['legends']:
                position = [legend['xmin'], legend['ymin'], legend['xmax'], legend['ymax']]
                if legend['bool']:
                    plotter.create_legend(position, **legend['kwargs'])

        if 'graphs' in plot:
            for graph in plot['graphs']:
                if not os.path.exists(graph['inPath']):
                    print(f"Error: File {graph['inPath']} does not exist. Skipping graph {graph['graphName']}.")
                    continue
                plotter.add_graph(graph['inPath'], graph['graphName'], graph['graphLabel'], **graph['kwargs'])

        if 'hists' in plot.keys():
            for hist in plot['hists']:
                if not os.path.exists(hist['inPath']):
                    print(f"Error: File {hist['inPath']} does not exist. Skipping histogram {hist['histName']}.")
                    continue
                plotter.add_hist(hist['inPath'], hist['histName'], hist['histLabel'], **hist['kwargs'])
        
        if len(plotter.graph_dict) == 0 and len(plotter.hist_dict) == 0:
            print(f"Warning: No graphs or histograms were added for plot with output {plot['outPDF']}. Check the configuration and file paths.")
            plotter.reset()
            continue

        if 'ROIs' in plot.keys():
            for roi in plot['ROIs']:
                plotter.add_ROI(roi['lineSpecs'], roi['boxSpecs'], **roi['kwargs'])

        if 'lines' in plot.keys():
            for line in plot['lines']:
                plotter.add_line(line['lineSpecs'], **line['kwargs'])

        if 'funcs' in plot:
            for func in plot['funcs']:
                plotter.add_func(func['inPath'], func['funcName'], func['funcLabel'], **func['kwargs'])
        
        if 'multigraph' in plot:
            plotter.draw_multigraph(**plot['multigraph']['kwargs'])

        if 'texts' in plot:
            for text in plot['texts']:
                plotter.add_text(text['text'], text['position'], **text['kwargs'])
        
        plotter.draw_legend()
        print(f'{plotter._canvas=}')
        plotter.save(plot['outPDF'])
    
    plotter.outfile.Close()