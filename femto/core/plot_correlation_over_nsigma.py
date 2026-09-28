from ROOT import TFile, TCanvas, TPad, TH1F, TLine, TLatex, gStyle, TColor, TGraphAsymmErrors, TPaveText
from torchic.utils.root import set_root_object
from torchic.core.histogram import load_hist
from torchic.utils.colors import get_color

SIGNIFICANCE_DICT = {
    'Matter': {
        '010': 0.90,
        '1030': 2.78,
        '3050': 1.45,
        '1050': 3.23
    },
    'Antimatter': {
        '010': 1.04,
        '1030': 1.67,
        '3050': 2.75,
        '1050': 3.12
    },
    'Both': {
        '010': 0.66,
        '1030': 3.23,
        '3050': 3.20,
        '1050': 4.61
    }
}

def compute_bkg_bandwidth(available_bkgs:list, bkg_file_path:str, nominal_bkg:str=None):
    """
    Compute, bin by bin, the fractional up/down deviation of the background
    variations with respect to a nominal background histogram.
    If nominal_bkg is not given, the first entry of available_bkgs is used as nominal
    (by convention this is the 'nominal/..._nominal' histogram).
    Returns lists: x_centers, ex (half bin width), ey_low_rel, ey_high_rel
    """
    h_available_bkgs = [load_hist(bkg_file_path, bkg_name) for bkg_name in available_bkgs]
    h_nominal = load_hist(bkg_file_path, nominal_bkg) if nominal_bkg is not None else h_available_bkgs[0]

    x, ex, ey_low_rel, ey_high_rel = [], [], [], []
    for i in range(1, h_nominal.GetNbinsX() + 1):
        nominal_value = h_nominal.GetBinContent(i)
        x.append(h_nominal.GetBinCenter(i))
        ex.append(h_nominal.GetBinWidth(i) / 2)
        if nominal_value == 0:
            ey_low_rel.append(0.)
            ey_high_rel.append(0.)
            continue
        variations = [h.GetBinContent(i) for h in h_available_bkgs]
        ey_high_rel.append((max(variations) - nominal_value) / nominal_value)
        ey_low_rel.append((nominal_value - min(variations)) / nominal_value)

    return x, ex, ey_low_rel, ey_high_rel


def create_bkg_band_from_model(model_curve, x_values, ex_values, ey_low_rel, ey_high_rel):
    """
    Build a TGraphAsymmErrors whose central values come from model_curve
    (e.g. the fitted background TGraph/RooCurve) and whose up/down widths
    come from the relative bandwidth of the background variations.
    """
    graph = TGraphAsymmErrors(len(x_values))
    for i, (x, ex, dlo, dhi) in enumerate(zip(x_values, ex_values, ey_low_rel, ey_high_rel)):
        y = model_curve.Eval(x)
        graph.SetPoint(i, x, y)
        graph.SetPointError(i, ex, ex, y * dlo, y * dhi)
    return graph

def plot_correlation_over_nsigma(file:TFile, pdf_path:str, x_limits:list, sign:str, centrality:str, 
                                 available_bkgs:list=None, bkg_file_path:str=None, normalisation_value:float=None,
                                 use_systematics:bool=False):

    sign_str = '' if sign == 'Both' else sign
    h_systematics_name = f'Correlation{sign_str}/{centrality}/hCorrelationSyst{centrality}' # if centrality != '1050' else f'{sign}/hCorrelationDirectComputation1050Syst'
    #h_systematics = load_hist('/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root',
    #                          h_systematics_name)
    h_systematics = load_hist('/home/galucia/Lithium4/preparation/output/PbPb/systematic_uncertainties.root',
                              h_systematics_name)
    if h_systematics is not None:
        set_root_object(h_systematics, marker_color=601, fill_color_alpha=(601, 0.3), line_color=601, marker_style=1, marker_size=0)

    y_portion = 0.3
    canvas = TCanvas('cCorrelationOverNsigma', '', 1200, 1400)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    
    canvas.cd()
    upper_pad = TPad('upper_pad', '', 0, y_portion - 0.05, 1, 1.)
    upper_pad.Draw()

    x_limits[1] = x_limits[1] - 0.001

    upper_canvas = file.Get('model/fit_signal')
    upper_canvas.SetTitle('')
    upper_canvas.SetBottomMargin(0.0)
    canvas_primitives = [upper_canvas.FindObject(primitive) for primitive in upper_canvas.GetListOfPrimitives()]
    canvas_primitives_dict = {}
    upper_pad.cd()

    ymax = 0.

    for primitive in canvas_primitives:
        #if 'frame' in primitive.GetName():          
        #    ymax = primitive.GetMaximum()
        
        if 'hCorrelation' in primitive.GetName():
            for ibin in range(1, primitive.GetN()):
                if primitive.GetPointY(ibin) > ymax:
                    ymax = primitive.GetPointY(ibin)

        if primitive.GetName() not in canvas_primitives_dict:
            canvas_primitives_dict[primitive.GetName()] = primitive
        else:
            n_priors = 0
            for names in canvas_primitives_dict.keys():
                if primitive.GetName() in names:
                    n_priors += 1
            canvas_primitives_dict[f'{primitive.GetName()};{n_priors}'] = primitive
        
    hframe_upper = upper_canvas.DrawFrame(x_limits[0], 0., x_limits[1], ymax*1.5, f';;#it{{C}}(#it{{k}}*)')
    hframe_upper.GetYaxis().SetTitleSize(0.05)
    #hframe_upper.GetXaxis().SetLabelSize(0.045)
    hframe_upper.GetYaxis().SetLabelSize(0.045)
    hframe_upper.GetYaxis().SetTitleOffset(0.9)
    hframe_upper.GetXaxis().SetLabelSize(0.)    # <-- add: hide x labels
    hframe_upper.GetXaxis().SetTitleSize(0.)
    hframe_upper.GetYaxis().ChangeLabel(1, -1, 0)  # set first label size to 0 (invisible)
    
    sign_label = 'p#minus^{3}He' if sign == 'Matter' else ('#bar{p}#minus^{3}#bar{He}' if sign == 'Antimatter' 
                                                           else 'p#minus^{3}He #oplus #bar{p}#minus^{3}#bar{He}')
    
    correlation_name = None
    for name, primitive in canvas_primitives_dict.items():
        if 'frame' in name:          
            continue

        if 'hCorrelation' in name:   
            set_root_object(primitive, marker_color=1, line_color=1, marker_style=20, marker_size=1.7)
            primitive.GetXaxis().SetLimits(x_limits[0], x_limits[1])
            correlation_name = name
            continue # save for last

        if 'signal_pdf' in name:
            set_root_object(primitive, line_color=get_color(0))
            #continue

        #if 'model' in name:
        #    continue
        
        if 'TPave' in name:  # legend
            primitive.SetX1(0.32)
            primitive.SetX2(0.75)
            primitive.SetY1(0.15)
            primitive.SetY2(0.35)
            primitive.SetMargin(0.1)
            correlation_names = [_name for _name in canvas_primitives_dict.keys() if 'hCorrelation' in _name]
            #primitive.AddEntry(canvas_primitives_dict[correlation_names[0]], sign_label, 'p')

        primitive.Draw('p same' if 'hCorrelation' in name else 'same')

    if available_bkgs is not None and bkg_file_path is not None:
        available_bkgs = [f'{sign}/{centrality}/{bkg_name}' for bkg_name in available_bkgs]
        x_values, ex_values, ey_low_rel, ey_high_rel = compute_bkg_bandwidth(available_bkgs, bkg_file_path)
        bkg_model_curve = canvas_primitives_dict.get('bkg_pdf')
        bkg_band = create_bkg_band_from_model(bkg_model_curve, x_values, ex_values, ey_low_rel, ey_high_rel)
        set_root_object(bkg_band, fill_color_alpha=(get_color(3), 1), line_color=get_color(3))
        
        model_curve = canvas_primitives_dict.get('model')
        model_band = create_bkg_band_from_model(model_curve, x_values, ex_values, ey_low_rel, ey_high_rel)
        set_root_object(model_band, fill_color_alpha=(get_color(1), 1), line_color=get_color(1))
        
        bkg_band.Draw('e3 same')
        model_band.Draw('e3 same')
    canvas_primitives_dict[correlation_name].Draw('p same')
    
    # Draw systematics as shaded area 
    if h_systematics is not None:
        h_systematics.SetFillColorAlpha(1, 0.3)
        h_systematics.SetLineColor(1)
        h_systematics.SetMarkerColor(1)  
        h_systematics.Draw('e2 same')
        
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.32, 0.8, f'ALICE Preliminary')
    latex.DrawLatex(0.32, 0.8, f'ALICE')
    centrality_label = f'{centrality[:1]}#minus{centrality[1:]}' if len(centrality) < 4 else f'{centrality[:2]}#minus{centrality[2:]}'
    latex.DrawLatex(0.32, 0.74, f'Pb#minusPb #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV')
    #latex.DrawLatex(0.32, 0.26, f'{sign_label}')
    latex.DrawLatex(0.32, 0.68, f'FT0C Centrality: {centrality_label}%')

    canvas.cd()
    lower_pad = TPad('lower_pad', '', 0, 0., 1, y_portion + 0.024)
    lower_pad.SetBottomMargin(0.3)
    lower_pad.SetTopMargin(0.)
    lower_pad.Draw()

    h_nsigma = file.Get('model/nsigma')
    set_root_object(h_nsigma, marker_color=get_color(3), marker_style=20, marker_size=1.7)
    h_nsigma_model = file.Get('model/nsigma_model')
    set_root_object(h_nsigma_model, marker_color=get_color(1), marker_style=20, marker_size=1.7)
    ### h_nsigma_syst = file.Get('model/nsigma_syst')
    ### set_root_object(h_nsigma_syst, marker_color=get_color(3), fill_color_alpha=(get_color(3), 0.3), line_color=get_color(3), marker_style=20, marker_size=1.7)
    ### h_nsigma_model_syst = file.Get('model/nsigma_model_syst')
    ### set_root_object(h_nsigma_model_syst, marker_color=get_color(1), fill_color_alpha=(get_color(1), 0.3), line_color=get_color(1), marker_style=20, marker_size=1.7)

    x_step = h_nsigma.GetBinWidth(1)
    nbins = int((x_limits[1] - x_limits[0])/x_step)
    h_nsigma_canvas = TH1F('h_nsigma_canvas', f';{h_nsigma.GetXaxis().GetTitle()};Pull', nbins, *x_limits)
    for ibin in range(1, nbins+1):
        x_value = h_nsigma_canvas.GetBinCenter(ibin)
        h_nsigma_canvas.SetBinContent(ibin, h_nsigma.GetBinContent(h_nsigma.FindBin(x_value)))
    set_root_object(h_nsigma_canvas, marker_color=797, marker_style=20, x_title_size=0.1, y_title_size=0.1,
                    x_title_offset=0.8, y_title_offset=0.3, x_label_size=0.1, y_label_size=0.1)
    
    lower_pad.cd()
    minimum_nsigma_canvas = h_nsigma_canvas.GetMinimum() * 0.7 if h_nsigma_canvas.GetMinimum() > 0 else h_nsigma_canvas.GetMinimum() * 1.3
    maximum_nsigma_canvas = h_nsigma_canvas.GetMaximum() * 1.3 if h_nsigma_canvas.GetMaximum() > 0 else h_nsigma_canvas.GetMaximum() * 0.7
    hframe = lower_pad.DrawFrame(x_limits[0], minimum_nsigma_canvas, x_limits[1], maximum_nsigma_canvas, 
                                 f';{h_nsigma.GetXaxis().GetTitle()};Pull')
    hframe.GetYaxis().SetNdivisions(5)
    set_root_object(hframe, x_title_size=0.1, y_title_size=0.1,
                    x_title_offset=0.8, y_title_offset=0.3, x_label_size=0.1, y_label_size=0.1)
    line = TLine(x_limits[0], 0., x_limits[1], 0.)
    set_root_object(line, line_style=2, line_color=15, line_width=2)

    #lower_text = TLatex()
    #lower_text.SetNDC()
    #lower_text.SetTextSize(0.07)
    #lower_text.SetTextFont(42)
    #lower_text.DrawLatex(0.5, 0.82, '#it{n}#sigma = #frac{#it{C}_{data}(#it{k}*) #minus #it{C}_{interaction}(#it{k}*)}{#sigma_{stat.}(#it{k}*)}')

    hframe.GetXaxis().SetTitleSize(0.1)
    hframe.GetYaxis().SetTitleSize(0.12)
    hframe.GetXaxis().SetLabelSize(0.1)
    hframe.GetYaxis().SetLabelSize(0.1)
    hframe.GetXaxis().SetTitleOffset(1.)
    hframe.GetYaxis().SetTitleOffset(0.38)
    hframe.GetXaxis().SetLimits(x_limits[0], x_limits[1])


    line.Draw('same')

    ### if use_systematics:
    ###     h_nsigma_syst.Draw('p0 same')
    ###     h_nsigma_model_syst.Draw('p0 same')
    ### else:
    h_nsigma.Draw('p0 same')
    h_nsigma_model.Draw('p0 same')

    canvas.SaveAs(pdf_path)
    

    file.cd()
    canvas.Write()
    
    text = TPaveText(0.24, 0.38, 0.75, 0.48, 'NDC')
    text.SetFillStyle(0)
    text.SetBorderSize(0)
    text.SetTextSize(0.045)
    text.SetTextFont(42)
    
    if centrality in SIGNIFICANCE_DICT[sign].keys():
        #if sign != 'Both':
        #    text.AddText(f'Significance: {SIGNIFICANCE_DICT[sign][centrality]:.2f} #sigma')
        #else:
        #    text.AddText(f'Significance (^{{4}}Li): {SIGNIFICANCE_DICT["Matter"][centrality]:.2f} #sigma')
        #    text.AddText(f'Significance (^{{4}}#bar{{Li}}): {SIGNIFICANCE_DICT["Antimatter"][centrality]:.2f} #sigma')
        text.AddText(f'Significance: {SIGNIFICANCE_DICT[sign][centrality]:.2f} #sigma')
        upper_pad.cd()
        text.Draw('same')
    
    hframe_upper.GetXaxis().SetTitleSize(0.05)
    #hframe_upper.GetXaxis().SetLabelSize(0.045)
    hframe_upper.GetXaxis().SetLabelSize(0.045)
    hframe_upper.GetXaxis().SetTitleOffset(0.9)
    hframe_upper.GetXaxis().SetTitle('#it{k}* (GeV/#it{c})')
    upper_pad.SaveAs(pdf_path.replace('.pdf', '_no_pull.pdf'))
    
    del canvas
