from ROOT import TPaveText

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
