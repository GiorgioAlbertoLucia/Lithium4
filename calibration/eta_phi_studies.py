
import numpy as np
from ROOT import TFile
from torchic import Dataset, AxisSpec

def get_phi_star(phi, pt, radius):
    '''
        Correction for the magnetic field
        radius: radius in the tpc at which the correction is computed (in meters)
    '''

    B = 0.5 # magnetic field in Tesla
    phi_star = (phi + np.arcsin(0.3 * radius * B / (2 * pt)) + np.pi) % 2*np.pi
    phi_star = phi_star - np.pi
    return phi_star

v_get_phi_star = np.vectorize(get_phi_star)

def run_on_delta_phi(dataset: Dataset, radius:float, iter: int, outfile: TFile):

    dataset[f'fPhiStar{iter}He3'] = v_get_phi_star(dataset['fPhiHe3'], dataset['fPtHe3'], radius)
    dataset[f'fPhiStar{iter}Had'] = v_get_phi_star(dataset['fPhiHad'], dataset['fPtHad'], radius)
    dataset[f'fDeltaPhiStar{iter}'] = dataset[f'fPhiStar{iter}He3'] - dataset[f'fPhiStar{iter}Had'] 

    axis_spec_deltaeta = AxisSpec(100, -0.10, 0.10)
    axis_spec_deltaphi = AxisSpec(100, -0.10, 0.10)
    h2_deltaeta_deltaphi = dataset.build_th2(f'fDeltaPhiStar{iter}', 'fDeltaEta', 
                                             axis_spec_deltaphi, axis_spec_deltaeta,
                                             name=f'h2_deltaeta_deltaphi_{iter}',
                                             title=f'radius = {radius} m; #Delta#phi * (radians); #Delta#eta;')

    outfile.cd()    
    h2_deltaeta_deltaphi.Write()

if __name__ == '__main__':

    infiles = ['/data/galucia/lithium_local/same_merged/LHC23_PbPb_pass5_same_small.root',
               '/data/galucia/lithium_local/same_merged/LHC24ar_pass2_same_small.root',]
    columns = ['fPtHe3', 'fEtaHe3', 'fPhiHe3', 'fPtHad', 'fEtaHad', 'fPhiHad']
    dataset = Dataset.from_root(infiles, tree_name='O2he3hadtable', folder_name='DF*', columns=columns)
    dataset['fDeltaEta'] = dataset['fEtaHe3'] - dataset['fEtaHad']

    outfile = TFile.Open('output/eta_phi_studies.root', 'recreate')
    radii = np.linspace(0.8, 2.5, 18)   # inner radius of the TPC: 80 cm, outer radius: 2.5 m

    for iter, radius in enumerate(radii):

        print(f'Running iteration {iter} for radius {radius} m')
        run_on_delta_phi(dataset, radius, iter, outfile)

    outfile.Close()
