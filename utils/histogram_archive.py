'''
    Set of histograms commonly used in the annalysis
'''

from utils.histogram_registry import HistogramRegistry, RegistryEntry
from copy import deepcopy

PT_NBINS, PT_MIN, PT_MAX = 200, -10, 10
KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX = 200, 0, 2.0

class Archive:
    QA_HISTOGRAMS = {
        "h2PtClusterSizeHe": RegistryEntry("h2PtClusterSizeHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});Cluster size", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fClusterSizeCosLamHe3", 90, 0, 15, 'true', 'QA'),
        "h2PtClusterSizePr": RegistryEntry("h2PtClusterSizePr", "p ;#it{p}_{T} (GeV/#it{c});Cluster size", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fClusterSizeCosLamHad", 90, 0, 15, 'true', 'QA'),
        "h2PtNSigmaITSHe": RegistryEntry("h2PtNSigmaITSHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaITSHe3", 100, -4, 4, 'true', 'QA'),
        "h2PtNSigmaTPCHe": RegistryEntry("h2PtNSigmaTPCHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTPCHe3", 100, -4, 4, 'true', 'QA'),
        "h2PpctNSigmaTPCHe": RegistryEntry("h2PpctNSigmaTPCHe", "^{3}He ;#it{p}_{TPC} (GeV/#it{c});n#sigma_{TPC}", "fInnerParamTPCHe3", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTPCHe3", 100, -4, 4, 'true', 'QA'),
        "h2PtMassTOFHe": RegistryEntry("h2PtMassTOFHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});#it{m}_{TOF}", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fMassTOFHe3", 200, 0, 4, 'true', 'QA'),
        "h2PtNSigmaITSPr": RegistryEntry("h2PtNSigmaITSPr", "p ;#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaITSHad", 100, -4, 4, 'true', 'QA'),
        "h2PtNSigmaTPCPr": RegistryEntry("h2PtNSigmaTPCPr", "p ;#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTPCHadPr", 100, -4, 4, 'true', 'QA'),
        "h2PtNSigmaTOFPr": RegistryEntry("h2PtNSigmaTOFPr", "p ;#it{p}_{T} (GeV/#it{c});n#sigma_{TOF}", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTOFHad", 100, -4, 4, 'true', 'QA'),
        "h2DeltaEtaDeltaPhi": RegistryEntry("h2DeltaEtaDeltaPhi", ";#Delta#eta;#Delta#phi (rad)", "fDeltaEta", 100, -0.1, 0.1, "fDeltaPhi", 100, -0.1, 0.1, 'true', 'QA'),

        "hPtHe": RegistryEntry("hPtHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});", "fSignedPtHe3", 200, -10, 10., condition='true', save_directory='QA'),
        "hEtaHe": RegistryEntry("hEtaHe", "^{3}He ;#eta;", "fEtaHe3", 100, -1, 1., condition='true', save_directory='QA'),
        "hPhiHe": RegistryEntry("hPhiHe", "^{3}He ;#phi (rad);", "fPhiHe3", 100, -3.14, 3.14, condition='true', save_directory='QA'),
        "hPtPr": RegistryEntry("hPtPr", "p ;#it{p}_{T} (GeV/#it{c});", "fSignedPtHad", 200, -10, 10., condition='true', save_directory='QA'),
        "hPtLi": RegistryEntry("hPtLi", "^{3}He+p ;#it{p}_{T} (GeV/#it{c});", "fSignedPtLi", 200, -10, 10., condition='true', save_directory='QA'),
        
        "hKstar": RegistryEntry("hKstar", ";#it{k}* (GeV/#it{c});", "fKstar", 100, 0, 0.4, condition='true', save_directory='QA'),
        "hCentrality": RegistryEntry("hCentrality", ";Centrality FT0C (%);", "fCentralityFT0C", 100, 0, 100, condition='true', save_directory='QA'),
        "hInvariantMass": RegistryEntry("hInvariantMass", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition='true', save_directory='QA'),

        "hCentralityKstar": RegistryEntry("hCentralityKstar", ";Centrality FT0C (%);#it{k}* (GeV/#it{c})", "fCentralityFT0C", 100, 0, 100, "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 'true', 'kstar'),
        "hPLiKstar": RegistryEntry("hPLiKstar", ";#it{p} (^{4}Li) (GeV/#it{c});#it{k}* (GeV/#it{c})", "fPLi", 100, 0, 10, "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 'fCentralityFT0C < 50', 'kstar'),
        "hPtHadKstar": RegistryEntry("hPtHadKstar", ";#it{p}_{T} (p) (GeV/#it{c});#it{k}* (GeV/#it{c})", "fPtHad", 100, 0, 10, "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 'fCentralityFT0C < 50', 'kstar'),
        "hPtHe3Kstar": RegistryEntry("hPtHe3Kstar", ";#it{p}_{T} (^{3}He) (GeV/#it{c});#it{k}* (GeV/#it{c})", "fPtHe3", 100, 0, 10, "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 'fCentralityFT0C < 50', 'kstar'),

        "h2PtDCAxyHe3": RegistryEntry("h2PtDCAxyHe3", "^{3}He ;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fDCAxyHe3", 100, -0.05, 0.05, 'true', 'QA'),
        "h2PtDCAzHe3": RegistryEntry("h2PtDCAzHe3", "^{3}He ;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fDCAzHe3", 100, -0.05, 0.05, 'true', 'QA'),
        "h2PtDCAxyHad": RegistryEntry("h2PtDCAxyHad", " p ;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fDCAxyHad", 100, -0.05, 0.05, 'true', 'QA'),
        "h2PtDCAzHad": RegistryEntry("h2PtDCAzHad", "p ;#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fDCAzHad", 100, -0.05, 0.05, 'true', 'QA'),

        "h2PhiChi2He3": RegistryEntry("h2PhiChi2He3", "^{3}He ;#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", "fPhiHe3", 100, -3.14, 3.14, "fChi2TPCHe3", 100, 0, 10, 'true', 'QA'),
        "h2Is23Chi2He3": RegistryEntry("h2Is23Chi2He3", "^{3}He ;Is23 ;#chi^{2}_{TPC} / N_{clusters TPC}", "fIs23", 2, 0, 2, "fChi2TPCHe3", 100, 0, 10, 'true', 'QA'),
        "h2PhiChi2Had": RegistryEntry("h2PhiChi2Had", "p ;#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", "fPhiHad", 100, -3.14, 3.14, "fChi2TPCHad", 100, 0, 10, 'true', 'QA'),
        "h2PtSharedClustersTpcHe": RegistryEntry("h2PtSharedClustersTpcHe", "^{3}He ;#it{p}_{T} (GeV/#it{c});Shared clusters TPC (^{3}He)", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fSharedClustersHe3", 100, 0, 100, 'true', 'QA'),
        "hPidForTrackingHe3": RegistryEntry("hPidForTrackingHe3", ";PID label;", "fPIDtrkHe3", 32, -0.5, 31.5, condition='true', save_directory='QA', 
                                            labels_x=['e', '#mu', '#pi', 'K', 'p', 'd', '^{3}H','^{3}He', '^{4}He', '#pi^{0}', '#gamma','K^{0}', '#Lambda', '^{3}_{#Lambda}H', '^{4}_{#Lambda}H', '#Xi^{-}', '#Omega^{-}', '^{4}_{#Lambda}He', '^{5}_{#Lambda}He'])
    }

    INVMASS_HISTOGRAMS = {}
    for sign, condition in {'Matter': 'fSignedPtHe3 > 0', 'Antimatter': 'fSignedPtHe3 < 0'}.items():
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}"] = RegistryEntry(f"hInvariantMass{sign}", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition}', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLi2"] = RegistryEntry(f"hInvariantMass{sign}PLi2", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi > 2', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLi3"] = RegistryEntry(f"hInvariantMass{sign}PLi3", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi > 3', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLi4"] = RegistryEntry(f"hInvariantMass{sign}PLi4", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi > 4', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLi5"] = RegistryEntry(f"hInvariantMass{sign}PLi5", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi > 5', save_directory='invmass')
        
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLiUnder2"] = RegistryEntry(f"hInvariantMass{sign}PLiUnder2", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi < 2', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLiUnder3"] = RegistryEntry(f"hInvariantMass{sign}PLiUnder3", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi < 3', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLiUnder4"] = RegistryEntry(f"hInvariantMass{sign}PLiUnder4", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi < 4', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PLiUnder5"] = RegistryEntry(f"hInvariantMass{sign}PLiUnder5", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPLi < 5', save_directory='invmass')

        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PtHadronUnder1p0"] = RegistryEntry(f"hInvariantMass{sign}PtHadronUnder1p0", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPtHad < 1.0', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PtHadronUnder1p1"] = RegistryEntry(f"hInvariantMass{sign}PtHadronUnder1p1", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPtHad < 1.1', save_directory='invmass')
        INVMASS_HISTOGRAMS[f"hInvariantMass{sign}PtHadronUnder1p2"] = RegistryEntry(f"hInvariantMass{sign}PtHadronUnder1p2", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition=f'{condition} && fPtHad < 1.2', save_directory='invmass')

    
    CENTRALITY_DICT = { # centrality : centrality_condition
        '010': 'fCentralityFT0C < 10',
        '1030': '10 <= fCentralityFT0C && fCentralityFT0C < 30',
        '3050': '30 <= fCentralityFT0C && fCentralityFT0C < 50',
    }
    KSTAR_HISTOGRAMS = {
    }
    for centrality, centrality_condition in CENTRALITY_DICT.items():
        KSTAR_HISTOGRAMS[f"hKstar{centrality}"] = RegistryEntry(f"hKstar{centrality}", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition}', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKt{centrality}"] = RegistryEntry(f"hKt{centrality}", ";#it{k}_{T} (GeV/#it{c});", "fKt", 500, 0, 5, condition=f'{centrality_condition}', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hMt{centrality}"] = RegistryEntry(f"hMt{centrality}", ";#it{m}_{T} (GeV/#it{c});", "fMt", 500, 0, 5, condition=f'{centrality_condition}', save_directory='kstar')
        
        KSTAR_HISTOGRAMS[f"hKstar{centrality}SharedUnder50"] = RegistryEntry(f"hKstar{centrality}SharedUnder50", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fSharedClustersHe3 < 50', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}SharedUnder40"] = RegistryEntry(f"hKstar{centrality}SharedUnder40", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fSharedClustersHe3 < 40', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}SharedUnder30"] = RegistryEntry(f"hKstar{centrality}SharedUnder30", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fSharedClustersHe3 < 30', save_directory='kstar')
        
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLi2"] = RegistryEntry(f"hKstar{centrality}PLi2", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi > 2', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLi3"] = RegistryEntry(f"hKstar{centrality}PLi3", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi > 3', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLi4"] = RegistryEntry(f"hKstar{centrality}PLi4", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi > 4', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLi5"] = RegistryEntry(f"hKstar{centrality}PLi5", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi > 5', save_directory='kstar')
        
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLiUnder2"] = RegistryEntry(f"hKstar{centrality}PLiUnder2", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi < 2', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLiUnder3"] = RegistryEntry(f"hKstar{centrality}PLiUnder3", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi < 3', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLiUnder4"] = RegistryEntry(f"hKstar{centrality}PLiUnder4", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi < 4', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PLiUnder5"] = RegistryEntry(f"hKstar{centrality}PLiUnder5", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPLi < 5', save_directory='kstar')
        
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PtHadronUnder1p0"] = RegistryEntry(f"hKstar{centrality}PtHadronUnder1p0", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPtHad < 1.0', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PtHadronUnder1p1"] = RegistryEntry(f"hKstar{centrality}PtHadronUnder1p1", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPtHad < 1.1', save_directory='kstar')
        KSTAR_HISTOGRAMS[f"hKstar{centrality}PtHadronUnder1p2"] = RegistryEntry(f"hKstar{centrality}PtHadronUnder1p2", ";#it{k}* (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition=f'{centrality_condition} && fPtHad < 1.2', save_directory='kstar')

def register_qa_histograms(registry: HistogramRegistry):
    """
        Register the QA histograms in the registry.
    """
    for name, entry in Archive.QA_HISTOGRAMS.items():
        registry.register(entry)

def select_qa_histograms(histogram_list: list):
    Archive.QA_HISTOGRAMS = {histogram: Archive.QA_HISTOGRAMS[histogram] for histogram in histogram_list}

def register_invmass_histograms(registry: HistogramRegistry):
    """
        Register the QA histograms in the registry.
    """
    for name, entry in Archive.INVMASS_HISTOGRAMS.items():
        registry.register(entry)

def select_invmass_histograms(histogram_list: list):
    Archive.INVMASS_HISTOGRAMS = {histogram: Archive.INVMASS_HISTOGRAMS[histogram] for histogram in histogram_list}


def register_kstar_histograms(registry: HistogramRegistry):
    """
        Register the kstar histograms in the registry.
    """
    for name, entry in Archive.KSTAR_HISTOGRAMS.items():
        registry.register(entry)


def register_kstar_matter_histograms(registry: HistogramRegistry):
    """
        Register the kstar matter histograms in the registry.
    """
    for name, _entry in Archive.KSTAR_HISTOGRAMS.items():
        entry = deepcopy(_entry)
        entry.name = entry.name + 'Matter'
        entry.save_directory = 'kstarMatter'
        entry.condition += " && fSignedPtHe3 > 0"
        registry.register(entry)

def register_kstar_antimatter_histograms(registry: HistogramRegistry):
    """
        Register the kstar antimatter histograms in the registry.
    """
    for name, _entry in Archive.KSTAR_HISTOGRAMS.items():
        entry = deepcopy(_entry)
        entry.name = entry.name + 'Antimatter'
        entry.save_directory = 'kstarAntimatter'
        entry.condition += " && fSignedPtHe3 < 0"
        registry.register(entry)