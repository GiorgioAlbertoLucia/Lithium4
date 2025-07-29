'''
    Upload particle dictionary
'''

ParticleMasses = {
    'Unidentified': 0.,
    'El':   0.0005109989461,
    'Pi':   0.13957018,
    'Ka':   0.493677,
    'Pr':   0.938272081,
    'De':   1.875613,
    'He':   2.80839160743,
}

ParticlePDG = {
    'Unidentified': 0,
    'El': 11, 
    'Pi': 211, 
    'Ka': 321, 
    'Pr': 2212,
    'De': 1000010020,
    'He': 1000020030,}

ParticlePID = {
    'El': 0,
    'Mu': 1,
    'Pi': 2,
    'Ka': 3,
    'Pr': 4,
    'De': 5,
    'Tr': 6,    # H3
    'He': 7,
    'Al': 8,    # He4
    'P0': 9,    # Pi0
    'Ph': 10,   # Photon
    'K0': 11,   # K0
    'La': 12,   # Lambda
    'HT': 13,   # HyperTriton
    'H4L': 14,  # Hyper H4
    'Xi': 15,
    'Om': 16,   # Omega
}

ParticleLabels = {
    'El': 'e',
    'Mu': '#mu',
    'Pi': '#pi',
    'Ka': 'K',
    'Pr': 'p',
    'De': 'd',
    'Tr': '^{3}H',
    'He': '^{3}He',
    'Al': '^{4}He',
    'P0': '#pi^{0}',
    'Ph': '#gamma',
    'K0': 'K^{0}',
    'La': '#Lambda',
    'HT': '^{3}_{#Lambda}H',
    'H4L': '^{4}_{#Lambda}H',
    'Xi': '#Xi^{-}',
    'Om': '#Omega^{-}',
}