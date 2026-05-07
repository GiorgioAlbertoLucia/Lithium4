

from ROOT import TFile
from torchic.core.histogram import load_hist


if __name__ == "__main__":
    
    h_pt_he3 = load_hist('/home/galucia/Lithium4/calibration/output/dca_data_template_smaller_tolerance.root',
                         'He/hPtHe')
    h_efficiency = load_hist('/home/galucia/Lithium4/phase_space_studies/output/single_track_efficiency.root',
                             'He/h_efficiency_he3')
    
    h_pt_he3_matter = h_pt_he3.Clone('h_pt_he3_matter')
    h_pt_he3_antimatter = h_pt_he3.Clone('h_pt_he3_antimatter')

    for ibin in range(1, h_pt_he3.GetNbinsX() + 1):
        pt = h_pt_he3.GetXaxis().GetBinCenter(ibin)
        if pt < 0.:
            continue
        
        h_pt_he3_antimatter.SetBinContent(ibin, h_pt_he3.GetBinContent(h_pt_he3.FindBin(-pt)))

    h_pt_he3_matter_corrected = h_pt_he3_matter.Clone('h_pt_he3_matter_corrected')
    h_pt_he3_antimatter_corrected = h_pt_he3_antimatter.Clone('h_pt_he3_antimatter_corrected')

    for ibin in range(1, h_pt_he3.GetNbinsX() + 1):
        pt = h_pt_he3.GetXaxis().GetBinCenter(ibin)
        if pt < 0.:
            continue
        
        eff_matter = h_efficiency.GetBinContent(h_efficiency.FindBin(pt))
        eff_antimatter = h_efficiency.GetBinContent(h_efficiency.FindBin(-pt))

        if eff_matter > 0:
            h_pt_he3_matter_corrected.SetBinContent(ibin, h_pt_he3_matter.GetBinContent(ibin) / eff_matter)
        
        if eff_antimatter > 0:
            h_pt_he3_antimatter_corrected.SetBinContent(ibin, h_pt_he3_antimatter.GetBinContent(ibin) / eff_antimatter)
    
    h_ratio = h_pt_he3_antimatter.Clone('h_ratio_antimatter_over_matter')
    h_ratio.Divide(h_pt_he3_matter)
    
    h_ratio_corrected = h_pt_he3_antimatter_corrected.Clone('h_ratio_antimatter_over_matter_corrected')
    h_ratio_corrected.Divide(h_pt_he3_matter_corrected)
    
    outfile = TFile('output/he3_matter_antimatter_comparison.root', 'RECREATE')
    h_pt_he3_matter.Write()
    h_pt_he3_antimatter.Write()

    h_pt_he3_matter_corrected.Write()
    h_pt_he3_antimatter_corrected.Write()

    h_ratio.Write()
    h_ratio_corrected.Write()
