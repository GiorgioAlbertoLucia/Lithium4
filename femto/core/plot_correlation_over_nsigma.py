from ROOT import TFile, TCanvas, TPad, TH1F
from torchic.utils.root import set_root_object

def plot_correlation_over_nsigma(file:TFile, pdf_path:str, x_limits:list):

    y_portion = 0.3
    canvas = TCanvas('cCorrelationOverNsigma', '', 800, 800)
    
    canvas.cd()
    upper_pad = TPad('upper_pad', '', 0, y_portion, 1, 1.)
    upper_pad.Draw()

    upper_canvas = file.Get('model/fit_signal')
    upper_canvas.SetTitle('')
    canvas_primitives = [upper_canvas.FindObject(primitive) for primitive in upper_canvas.GetListOfPrimitives()]
    upper_pad.cd()
    for primitive in canvas_primitives:
        if 'frame' in primitive.GetName():          
            primitive.SetBinContent(1, 0)
            primitive.SetTitle(f';{primitive.GetXaxis().GetTitle()};C(#it{{k}}*)')
        if 'hCorrelation' in primitive.GetName():   set_root_object(primitive, marker_color=797, line_color=797)
        primitive.Draw('p same' if 'hCorrelation' in primitive.GetName() else 'same')
    #upper_canvas.Draw()

    canvas.cd()
    lower_pad = TPad('lower_pad', '', 0, 0.05, 1, y_portion)
    lower_pad.SetBottomMargin(0.2)
    lower_pad.Draw()

    h_nsigma = file.Get('model/nsigma')
    x_step = h_nsigma.GetBinWidth(1)
    nbins = int((x_limits[1] - x_limits[0])/x_step)
    h_nsigma_canvas = TH1F('h_nsigma_canvas', f';{h_nsigma.GetXaxis().GetTitle()};#it{{n}}#sigma', nbins, *x_limits)
    for ibin in range(1, nbins+1):
        x_value = h_nsigma_canvas.GetBinCenter(ibin)
        h_nsigma_canvas.SetBinContent(ibin, h_nsigma.GetBinContent(h_nsigma.FindBin(x_value)))
    set_root_object(h_nsigma_canvas, marker_color=797, marker_style=20, x_title_size=0.1, y_title_size=0.1,
                    x_title_offset=0.8, y_title_offset=0.3, x_label_size=0.1, y_label_size=0.1)
    lower_pad.cd()
    h_nsigma_canvas.Draw('p0')

    canvas.SaveAs(pdf_path)
    del canvas
