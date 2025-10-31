from dataclasses import dataclass
from ROOT import TCanvas, TPaveText, TLegend

class PdfHandler:

    def __init__(self, pdf_output_path:str):
        
        self._pdf_path = pdf_output_path
        self._canvas = TCanvas('', '', 900, 600)
        self._canvas.Print(f'{self._pdf_path}(')

    @property
    def canvas(self):
        return self._canvas
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._canvas.Clear()
        self._canvas.Print(f'{self._pdf_path})')

    def cd(self):
        self._canvas.cd()

    def load_canvas(self, canvas):
        self._canvas = canvas
    
    def save(self):
        self._canvas.Print(self._pdf_path)

    def clear(self):
        self._canvas.Clear()

    def save_and_clear(self):
        self._canvas.Print(self._pdf_path)
        self._canvas.Clear()

    def draw_save_and_clear(self, obj, **kwargs):
        self.canvas_settings(**kwargs)
        self._canvas.cd()
        obj.Draw(kwargs.get('draw_option', ''))
        self._canvas.Print(self._pdf_path)
        self._canvas.Clear()

    def draw_and_save(self, obj, **kwargs):
        self.canvas_settings(**kwargs)
        self._canvas.cd()
        obj.Draw(kwargs.get('draw_option', ''))
        self._canvas.Print(self._pdf_path)
    
    def canvas_settings(self, **kwargs):
        self._canvas.SetLogy(kwargs.get('logy', False))

def init_legend(xmin, ymin, xmax, ymax) -> TLegend:
    legend = TLegend(xmin, ymin, xmax, ymax)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    return legend

def write_params_to_text(params: tuple, coordinates:tuple=[0.7, 0.9, 0.7, 0.9]) -> TPaveText:
    '''
        Write the parameters of the fit to the canvas

        Args:
            params: dictionary of the parameters to write (string: RooRealVar)
    '''

    text = TPaveText(coordinates[0], coordinates[1], coordinates[2], coordinates[3], 'NDC')
    text.SetFillColor(0)
    text.SetBorderSize(0)
    for param in params:
        if param.isConstant():
            text.AddText(f'{param.GetTitle()} = {param.getVal():.4f} {param.getUnit()} (fixed)')
        else:
            text.AddText(f'{param.GetTitle()} = ({param.getVal():.4f} #pm {param.getError():.4f}) {param.getUnit()}')
    return text
