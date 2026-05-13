from ROOT import gInterpreter
gInterpreter.ProcessLine(f'#include "../include/Parametrisation.h"')
gInterpreter.ProcessLine(f'#include "../include/Common.h"')

from torchic.utils.terminal_colors import TerminalColors as tc

def load_parametrisation(config: dict):
    """Read parametrisation from config and push values into Common.h via setters."""

    params = config.get('parameterisation', {})
    if not params:
        print(tc.YELLOW + 'No parametrisation found in config, using defaults.' + tc.RESET)
        return

    species_index = {'Pr': 0, 'He': 1, 'Pi': 2}

    def call(setter, *args):
        args_str = ', '.join(str(a) for a in args)
        gInterpreter.ProcessLine(f'parametrisation::{setter}({args_str});')

    for sp, idx in species_index.items():
        if f'kDCAxyResolutionParams_{sp}' in params:
            call('SetDCAxyResolutionParams', idx, *params[f'kDCAxyResolutionParams_{sp}'])
        if f'kDCAzResolutionParams_{sp}' in params:
            call('SetDCAzResolutionParams', idx, *params[f'kDCAzResolutionParams_{sp}'])
        if f'kDCAxyResolutionParamsMC_{sp}' in params:
            call('SetDCAxyResolutionParamsMC', idx, *params[f'kDCAxyResolutionParamsMC_{sp}'])
        if f'kDCAzResolutionParamsMC_{sp}' in params:
            call('SetDCAzResolutionParamsMC', idx, *params[f'kDCAzResolutionParamsMC_{sp}'])
        if f'kITSParams_{sp}' in params:
            call('SetITSParams', idx, *params[f'kITSParams_{sp}'])
        if f'kITSParamsMC_{sp}' in params:
            call('SetITSParamsMC', idx, *params[f'kITSParamsMC_{sp}'])
        if f'kITSResolutionParams_{sp}' in params:
            call('SetITSResolutionParams', idx, *params[f'kITSResolutionParams_{sp}'])
        if f'kITSResolutionParamsMC_{sp}' in params:
            call('SetITSResolutionParamsMC', idx, *params[f'kITSResolutionParamsMC_{sp}'])

    scalar_setters = {
        'kHeTPCParams':         ('SetHeTPCParams',         5),
        'kHeTPCResolution':     ('SetHeTPCResolution',     1),
        'kHeTPCParamsResiduals':('SetHeTPCParamsResiduals',6),
        'kHeTPCParamsMC':       ('SetHeTPCParamsMC',       5),
        'kHeTPCResolutionMC':   ('SetHeTPCResolutionMC',   1),
        'kPrTPCParams':         ('SetPrTPCParams',         5),
        'kPrTPCResolution':     ('SetPrTPCResolution',     1),
        'kPrTOFParams':         ('SetPrTOFParams',         3),
        'kPrTOFResolutionParams':('SetPrTOFResolutionParams',2),
        'kHePidTrkParams':      ('SetHePidTrkParams',      3),
    }
    for key, (setter, _) in scalar_setters.items():
        if key in params:
            val = params[key]
            call(setter, *val) if isinstance(val, list) else call(setter, val)