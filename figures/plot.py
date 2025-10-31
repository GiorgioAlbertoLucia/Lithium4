import yaml

import sys
sys.path.append('..')
from torchic import Plotter

if __name__ == '__main__':

    input_file = 'config/checks_all_vs_cbt.yml'

    with open(input_file, 'r') as f:
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
                plotter.add_graph(graph['inPath'], graph['graphName'], graph['graphLabel'], **graph['kwargs'])

        if 'hists' in plot.keys():
            for hist in plot['hists']:
                plotter.add_hist(hist['inPath'], hist['histName'], hist['histLabel'], **hist['kwargs'])

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