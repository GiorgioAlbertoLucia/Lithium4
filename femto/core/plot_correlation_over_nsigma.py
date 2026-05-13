from ROOT import TFile, TCanvas, TPad, TH1F, TLine, TLatex, gStyle, TColor
from torchic.utils.root import set_root_object
from torchic.core.histogram import load_hist
from torchic.utils.colors import get_color

def plot_correlation_over_nsigma(file:TFile, pdf_path:str, x_limits:list, sign:str, centrality:str):

    h_systematics_name = f'{sign}/hCorrelation{centrality}Syst' if centrality != '1050' else f'{sign}/hCorrelationDirectComputation1050Syst'
    h_systematics = load_hist('/home/galucia/Lithium4/preparation/output/correlation_with_systematics.root',
                              h_systematics_name)
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
        if 'frame' in primitive.GetName():          
            ymax = primitive.GetMaximum()

        if primitive.GetName() not in canvas_primitives_dict:
            canvas_primitives_dict[primitive.GetName()] = primitive
        else:
            n_priors = 0
            for names in canvas_primitives_dict.keys():
                if primitive.GetName() in names:
                    n_priors += 1
            canvas_primitives_dict[f'{primitive.GetName()};{n_priors}'] = primitive
        
    hframe_upper = upper_canvas.DrawFrame(x_limits[0], 0., x_limits[1], ymax*1.3, f';;#it{{C}}(#it{{k}}*)')
    hframe_upper.GetYaxis().SetTitleSize(0.05)
    #hframe_upper.GetXaxis().SetLabelSize(0.045)
    hframe_upper.GetYaxis().SetLabelSize(0.045)
    hframe_upper.GetYaxis().SetTitleOffset(0.9)
    hframe_upper.GetXaxis().SetLabelSize(0.)    # <-- add: hide x labels
    hframe_upper.GetXaxis().SetTitleSize(0.)
    hframe_upper.GetYaxis().ChangeLabel(1, -1, 0)  # set first label size to 0 (invisible)
    
    sign_label = 'p#minus^{3}He' if sign == 'Matter' else '#bar{p}#minus^{3}#bar{He}'

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

    canvas_primitives_dict[correlation_name].Draw('p same')
    
    # Draw systematics as shaded area 
    h_systematics.SetFillColorAlpha(1, 0.3)
    h_systematics.SetLineColor(1)
    h_systematics.SetMarkerColor(1)  
    #h_systematics.Draw('e2 same')
        
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    latex.DrawLatex(0.32, 0.8, f'ALICE Preliminary')
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
    set_root_object(h_nsigma, marker_color=get_color(0), marker_style=20, marker_size=1.7)
    h_nsigma_model = file.Get('model/nsigma_model')
    set_root_object(h_nsigma_model, marker_color=get_color(2), marker_style=20, marker_size=1.7)

    x_step = h_nsigma.GetBinWidth(1)
    nbins = int((x_limits[1] - x_limits[0])/x_step)
    h_nsigma_canvas = TH1F('h_nsigma_canvas', f';{h_nsigma.GetXaxis().GetTitle()};Pull', nbins, *x_limits)
    for ibin in range(1, nbins+1):
        x_value = h_nsigma_canvas.GetBinCenter(ibin)
        h_nsigma_canvas.SetBinContent(ibin, h_nsigma.GetBinContent(h_nsigma.FindBin(x_value)))
    set_root_object(h_nsigma_canvas, marker_color=797, marker_style=20, x_title_size=0.1, y_title_size=0.1,
                    x_title_offset=0.8, y_title_offset=0.3, x_label_size=0.1, y_label_size=0.1)
    
    lower_pad.cd()
    minimum_nsigma_canvas = h_nsigma_canvas.GetMinimum() * 0.9 if h_nsigma_canvas.GetMinimum() > 0 else h_nsigma_canvas.GetMinimum() * 1.1
    maximum_nsigma_canvas = h_nsigma_canvas.GetMaximum() * 1.1 if h_nsigma_canvas.GetMaximum() > 0 else h_nsigma_canvas.GetMaximum() * 0.9
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
    #h_nsigma_canvas.Draw('p0 same')
    h_nsigma.Draw('p0 same')
    h_nsigma_model.Draw('p0 same')

    canvas.SaveAs(pdf_path)

    file.cd()
    canvas.Write()
    del canvas
