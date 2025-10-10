
import numpy as np
from ROOT import TFile, gInterpreter
from torchic import Dataset, AxisSpec
from torchic.core.histogram import build_TH2
from alive_progress import alive_bar
import yaml

import sys
sys.path.append('..')
from preparation.prepare import prepare_selections, prepare_input_tchain, prepare_rdataframe

gInterpreter.ProcessLine(f'#include "EtaPhiUtils.h"')
from ROOT import ComputePhiStar

def run_on_delta_phi_rdf(rdataframe: Dataset, radius:float, iter: int, outfile: TFile):

    h2 = rdataframe.Define('fPhiStarHe3', f'ComputePhiStar(fPhiHe3, fInnerParamTPCHe3 * 2, {radius})') \
                   .Define('fPhiStarHad', f'ComputePhiStar(fPhiHad, fInnerParamTPCHad, {radius})') \
                   .Define('fDeltaPhiStar', 'fPhiStarHe3 - fPhiStarHad') \
                   .Histo2D((f'h2DeltaEtaDeltaPhiStar{iter}',
                             f'#it{{r}}_{{TPC}} = {radius:.2f} m; #Delta#phi* (rad); #Delta#eta;',
                             50, -0.2, 0.2, 50, -0.2, 0.2),  'fDeltaPhiStar', 'fDeltaEta')
    
    outfile.cd()    
    h2.Write()



if __name__ == '__main__':

    config_file = '/home/galucia/Lithium4/preparation/config/config_prepare.yml'
    config = yaml.safe_load(open(config_file, 'r'))

    base_selection, selection = prepare_selections(config)
    chain_data = prepare_input_tchain(config)
    rdf = prepare_rdataframe(chain_data, base_selection, selection)

    outfile = TFile.Open('output/eta_phi_studies.root', 'recreate')
    radii = np.linspace(0.8, 2.5, 18)   # inner radius of the TPC: 80 cm, outer radius: 2.5 m

    with alive_bar(total=len(radii), title='Loop over TPC radii') as bar:
        for iter, radius in enumerate(radii):
            run_on_delta_phi_rdf(rdf, radius, iter, outfile)
            bar()

    outfile.Close()
