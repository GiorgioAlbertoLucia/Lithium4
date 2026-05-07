import numpy as np
import pandas as pd
from ROOT import TFile, RooRealVar, TDirectory, TGraphErrors, TF1, RooDataHist, TCanvas

from torchic import Dataset, AxisSpec
from torchic.roopdf.roopdf_utils import init_roopdf
from torchic.physics.particles import PARTICLES


if __name__ == '__main__':
    
    dataset = Dataset.from_root('/data/galucia/lithium_local/same/LHC23_PbPb_pass5_hadronpid_same.root',
                                tree_name='O2he3hadtable',
                                folder_name='DF*')
    print(f'{dataset.columns=}')
    
    dataset['fPHad'] = abs(dataset['fPtHad']) * np.cosh(dataset['fEtaHad'])
    print(f'{PARTICLES["Pr"].mass=}')
    dataset['fTOFbetagammaHad'] =  dataset['fPHad'] / PARTICLES['Pr'].mass
    
    axis_spec_betagamma = AxisSpec(100, 0, 10)
    axis_spec_dEdx = AxisSpec(2000, 0, 2000)
    h2_dEdx_vs_betagamma = dataset.build_th2('fTOFbetagammaHad', 'fSignalTPCHad', axis_spec_betagamma, axis_spec_dEdx,
                                             name='h2_dEdx_vs_betagamma', title=';#beta#gamma;#LT d#it{E}/d#it{x} #GT (a.u.)')
    
    outfile = TFile('output/dEdx_vs_betagamma.root', 'recreate')
    h2_dEdx_vs_betagamma.Write()
    outfile.Close()