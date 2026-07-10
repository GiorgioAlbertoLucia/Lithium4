"""
Template preparation and manipulation
"""

import numpy as np
from ROOT import TFile, TH2F, TCanvas, TLegend, RooDataHist, RooFFTConvPdf, RooNumConvPdf

from torchic.core.histogram import load_hist
from torchic.utils.root import set_root_object


def prepare_hist_material_template(input_data_file: str, input_data_name: str, 
                                   outfile: TFile, efficiency_file: str,
                                   efficiency_name: str, efficiency_material_name: str,
                                   particle: str = 'He') -> TH2F:
    """
    Prepare material template from data using matter/antimatter asymmetry
    
    Args:
        input_data_file: Path to data file
        input_data_name: Name of data histogram
        outfile: Output ROOT file
        efficiency_file: Path to efficiency file
        efficiency_name: Name of efficiency histogram
        efficiency_material_name: Name of material efficiency histogram
        particle: Particle type
    
    Returns:
        Material template histogram
    """
    h2_data_input = load_hist(input_data_file, input_data_name)
    h_efficiency = load_hist(efficiency_file, efficiency_name)
    h_efficiency_material = load_hist(efficiency_file, efficiency_material_name)

    h2_antimatter_template = h2_data_input.Clone(f'h2_antimatter_template_{particle}')
    
    for xbin in range(1, h2_data_input.GetNbinsX() + 1):
        pt = h2_data_input.GetXaxis().GetBinCenter(xbin)
        if pt < 0:
            continue

        pt_bin_antimatter = h2_data_input.GetXaxis().FindBin(-pt)
        efficiency_matter = h_efficiency.GetBinContent(h_efficiency.FindBin(pt))
        efficiency_antimatter = h_efficiency.GetBinContent(h_efficiency.FindBin(-pt))
        
        for ybin in range(1, h2_data_input.GetNbinsY() + 1):
            antimatter_content = h2_data_input.GetBinContent(pt_bin_antimatter, ybin)
            corrected_content = (antimatter_content * efficiency_matter / efficiency_antimatter 
                               if efficiency_antimatter > 0. else 0.)
            h2_antimatter_template.SetBinContent(xbin, ybin, corrected_content)
                
    h2_material_template = h2_data_input.Clone(f'h2_material_template_{particle}')
    h2_material_template.Add(h2_antimatter_template, -1.)

    outfile.cd()
    h2_data_input.Write('h2_data_template')
    h2_material_template.Write('h2_material_before_efficiency_correction')

    for xbin in range(1, h2_material_template.GetNbinsX() + 1):
        pt = h2_material_template.GetXaxis().GetBinCenter(xbin)
        if pt < 0:
            continue

        efficiency = h_efficiency.GetBinContent(h_efficiency.FindBin(pt))
        efficiency_material = h_efficiency_material.GetBinContent(h_efficiency.FindBin(pt))
        
        for ybin in range(1, h2_material_template.GetNbinsY() + 1):
            content = h2_material_template.GetBinContent(xbin, ybin)
            content = content if content > 0. else 0.
            corrected_content = (content * efficiency_material / efficiency 
                               if efficiency > 0. else 0.)
            h2_material_template.SetBinContent(xbin, ybin, corrected_content)

    outfile.cd()
    h2_antimatter_template.Write()
    h2_material_template.Write()

    return h2_material_template


def prepare_convoluted_template(x, pdf, params, h_mc, gaussian_smearing, gaussian_smearing_params,
                                bin_outdir, flag, particle, 
                                use_convolution: bool = True):
    """
    Prepare a convoluted PDF template from MC histogram
    
    Args:
        x: Observable variable
        pdf: Base PDF
        params: PDF parameters
        h_mc: MC histogram
        gaussian_smearing: Smearing Gaussian PDF
        bin_outdir: Output directory
        flag: Template flag (primaries, material, etc.)
        particle: Particle type
        use_convolution: Whether to convolve with smearing
    
    Returns:
        Convoluted PDF
    """

    # use larger range for fitting material templates to better constrain tails
    default_range = (x.getMin(), x.getMax())
    # Save original binning
    original_bins = x.getBins()
    original_min = x.getMin()
    original_max = x.getMax()

    n_uniform = 1000  # >= 1000 as RooFit itself suggests
    x.setBins(n_uniform)

    
    #if flag == 'material' or (flag == 'primaries' and particle == 'He'):
    #    x.setRange(default_range[0] - 0.05, default_range[1] + 0.05)

    dh = RooDataHist('tmp', '', [x], Import=h_mc)
    pdf.fitTo(dh, PrintLevel=1, Save=True)

    for param in params.values():
        param.setConstant(True)
    
    frame = x.frame()
    dh.plotOn(frame)
    pdf.plotOn(frame, Name=pdf.GetName(), LineColor=2)
    
    pdf_convoluted = pdf
    if use_convolution: # and not (flag == 'material' and particle != 'Pr'):
        pdf_convoluted = RooFFTConvPdf(pdf.GetName() + '_convoluted', 
                                      pdf.GetTitle() + '_convoluted', 
                                      x, pdf, gaussian_smearing, ipOrder=10)
        
        #pdf_convoluted = RooNumConvPdf(pdf.GetName() + '_convoluted', 
        #                              pdf.GetTitle() + '_convoluted', 
        #                              x, pdf, gaussian_smearing)
        
    
    #gaussian_smearing_params['mean'].setVal(0)

    pdf_convoluted.plotOn(frame, Name=pdf_convoluted.GetName(), 
                         LineColor=3, LineStyle='--')
    pdf_convoluted.paramOn(frame, ShowConstants=True)

    legend = TLegend(0.11, 0.6, 0.4, 0.89)
    legend.SetBorderSize(0)

    dh.plotOn(frame)
    legend.AddEntry(frame.findObject(pdf.GetName()), pdf.GetTitle(), "l")
    legend.AddEntry(frame.findObject(pdf_convoluted.GetName()), 
                   pdf_convoluted.GetTitle(), "l")

    canvas = TCanvas(pdf.GetName(), '')
    frame.SetMinimum(0.1)
    canvas.SetLogy()
    frame.Draw()
    legend.Draw('same')
    
    bin_outdir.cd()
    canvas.Write()

    x.setRange(*default_range) # reset to default range
    x.setBins(original_bins)

    del dh
    return pdf_convoluted

def prepare_template(x, pdf, params, h_mc, bin_outdir, flag, particle):
    """
    Prepare a convoluted PDF template from MC histogram
    
    Args:
        x: Observable variable
        pdf: Base PDF
        params: PDF parameters
        h_mc: MC histogram
        gaussian_smearing: Smearing Gaussian PDF
        bin_outdir: Output directory
        flag: Template flag (primaries, material, etc.)
        particle: Particle type
        use_convolution: Whether to convolve with smearing
    
    Returns:
        Convoluted PDF
    """

    default_range = (x.getMin(), x.getMax())
    if flag == 'material':
        x.setRange(default_range[0] - 0.05, default_range[1] + 0.05)

    dh = RooDataHist('tmp', '', [x], Import=h_mc)
    pdf.fitTo(dh, PrintLevel=1, Save=True)

    for param in params.values():
        param.setConstant(True)
    
    frame = x.frame()
    dh.plotOn(frame)
    pdf.plotOn(frame, Name=pdf.GetName(), LineColor=2)
    pdf.paramOn(frame, ShowConstants=True)
    
    dh.plotOn(frame)

    canvas = TCanvas(pdf.GetName(), '')
    frame.Draw()
    
    bin_outdir.cd()
    canvas.Write()

    x.setRange(*default_range) # reset to default range

    del dh
    return pdf