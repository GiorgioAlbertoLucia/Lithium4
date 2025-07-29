'''
    Set of histograms commoly used in the annalysis
'''

from utils.histogram_registry import HistogramRegistry, RegistryEntry
from copy import deepcopy

PT_NBINS, PT_MIN, PT_MAX = 200, -10, 10

QA_HISTOGRAMS = {
    "h2PtNSigmaTPCHe": RegistryEntry("h2PtNSigmaTPCHe", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTPCHe3", 100, -4, 4, 'true', 'QA'),
    "h2PtNSigmaTPCPr": RegistryEntry("h2PtNSigmaTPCPr", ";#it{p}_{T} (GeV/#it{c});n#sigma_{TPC}", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaTPCHad", 100, -4, 4, 'true', 'QA'),
    "h2PtClusterSizeHe": RegistryEntry("h2PtClusterSizeHe", ";#it{p}_{T} (GeV/#it{c});Cluster size", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fClusterSizeCosLamHe3", 90, 0, 15, 'true', 'QA'),
    "h2PtExpectedClusterSizeHe": RegistryEntry("h2PtExpectedClusterSizeHe", ";#it{p}_{T} (GeV/#it{c});Expected cluster size", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fExpectedClusterSizeHe3", 90, 0, 15, 'true', 'QA'),
    "h2PtClusterSizePr": RegistryEntry("h2PtClusterSizePr", ";#it{p}_{T} (GeV/#it{c});Cluster size", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fClusterSizeCosLamHad", 90, 0, 15, 'true', 'QA'),
    "h2PtExpectedClusterSizePr": RegistryEntry("h2PtExpectedClusterSizePr", ";#it{p}_{T} (GeV/#it{c});Expected cluster size", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fExpectedClusterSizeHad", 90, 0, 15, 'true', 'QA'),
    "h2PtNSigmaITSHe": RegistryEntry("h2PtNSigmaITSHe", ";#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaITSHe3", 100, -4, 4, 'true', 'QA'),
    "h2PtNSigmaITSPr": RegistryEntry("h2PtNSigmaITSPr", ";#it{p}_{T} (GeV/#it{c});n#sigma_{ITS}", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fNSigmaITSHad", 100, -4, 4, 'true', 'QA'),
    "h2DeltaEtaDeltaPhi": RegistryEntry("h2DeltaEtaDeltaPhi", ";#Delta#eta;#Delta#phi (rad)", "fDeltaEta", 100, -0.1, 0.1, "fDeltaPhi", 100, -0.1, 0.1, 'true', 'QA'),
    "hEta": RegistryEntry("hEtaHe", ";#it{k}^{*} (GeV/#it{c});", "fEtaHe3", 100, -1, 1., condition='true', save_directory='QA'),
    "hKstar": RegistryEntry("hKstar", ";#it{k}^{*} (GeV/#it{c});", "fKstar", 100, 0, 0.4, condition='true', save_directory='QA'),
    "hCentrality": RegistryEntry("hCentrality", ";Centrality (ft0c);", "fCentralityFT0C", 100, 0, 100, condition='true', save_directory='QA'),
    "hInvariantMass": RegistryEntry("hInvariantMass", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition='true', save_directory='QA'),
    "hInvariantMassMatter": RegistryEntry("hInvariantMassMatter", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition='fSignedPtHe3 > 0', save_directory='QA'),
    "hInvariantMassAntimatter": RegistryEntry("hInvariantMassAntimatter", ";Invariant mass (GeV/#it{c}^{2});", "fMassInvLi", 400, 3.747, 3.947, condition='fSignedPtHe3 < 0', save_directory='QA'),
    "hPtDCAxyHe3": RegistryEntry("hPtDCAxyHe3", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fDCAxyHe3", 100, -0.1, 0.1, 'true', 'QA'),
    "hPtDCAzHe3": RegistryEntry("hPtDCAzHe3", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", "fSignedPtHe3", PT_NBINS, PT_MIN, PT_MAX, "fDCAzHe3", 100, -1.0, 1.0, 'true', 'QA'),
    "hPtDCAxyHad": RegistryEntry("hPtDCAxyHad", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fDCAxyHad", 100, -0.1, 0.1, 'true', 'QA'),
    "hPtDCAzHad": RegistryEntry("hPtDCAzHad", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm)", "fSignedPtHad", PT_NBINS, PT_MIN, PT_MAX, "fDCAzHad", 100, -1.0, 1.0, 'true', 'QA'),
    "hPhiChi2He3": RegistryEntry("hPhiChi2He3", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", "fPhiHe3", 100, -3.14, 3.14, "fChi2TPCHe3", 100, 0, 10, 'true', 'QA'),
    "hPhiChi2Had": RegistryEntry("hPhiChi2Had", ";#phi (rad);#chi^{2}_{TPC} / N_{clusters TPC}", "fPhiHad", 100, -3.14, 3.14, "fChi2TPCHad", 100, 0, 10, 'true', 'QA'),
}

def register_qa_histograms(registry: HistogramRegistry):
    """
        Register the QA histograms in the registry.
    """
    for name, entry in QA_HISTOGRAMS.items():
        registry.register(entry)

KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX = 80, 0, 0.8
KSTAR_HISTOGRAMS = {
    "hCentralityKstar": RegistryEntry("hCentralityKstar", ";Centrality (ft0c);#it{k}^{*} (GeV/#it{c})", "fCentralityFT0C", 100, 0, 100, "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, 'true', 'kstar'),
    "hKstar010": RegistryEntry("hKstar010", ";#it{k}^{*} (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition='fCentralityFT0C < 10', save_directory='kstar'),
    "hKt010": RegistryEntry("hKt010", ";#it{k}_{T} (GeV/#it{c});", "fKt", 500, 0, 5, condition='fCentralityFT0C < 10', save_directory='kstar'),
    "hMt010": RegistryEntry("hMt010", ";#it{m}_{T} (GeV/#it{c});", "fMt", 500, 0, 5, condition='fCentralityFT0C < 10', save_directory='kstar'),
    "hKstar1030": RegistryEntry("hKstar1030", ";#it{k}^{*} (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition="10 <= fCentralityFT0C && fCentralityFT0C < 30", save_directory='kstar'),
    "hKt1030": RegistryEntry("hKt1030", ";#it{k}_{T} (GeV/#it{c});", "fKt", 500, 0, 5, condition="10 <= fCentralityFT0C && fCentralityFT0C < 30", save_directory='kstar'),
    "hMt1030": RegistryEntry("hMt1030", ";#it{m}_{T} (GeV/#it{c});", "fMt", 500, 0, 5, condition="10 <= fCentralityFT0C && fCentralityFT0C < 30", save_directory='kstar'),
    "hKstar3050": RegistryEntry("hKstar3050", ";#it{k}^{*} (GeV/#it{c});", "fKstar", KSTAR_NBINS, KSTAR_MIN, KSTAR_MAX, condition="30 <= fCentralityFT0C && fCentralityFT0C < 50", save_directory='kstar'),
    "hKt3050": RegistryEntry("hKt3050", ";#it{k}_{T} (GeV/#it{c});", "fKt", 500, 0, 5, condition="30 <= fCentralityFT0C && fCentralityFT0C < 50", save_directory='kstar'),
    "hMt3050": RegistryEntry("hMt3050", ";#it{m}_{T} (GeV/#it{c});", "fMt", 500, 0, 5, condition="30 <= fCentralityFT0C && fCentralityFT0C < 50", save_directory='kstar'),
}

def register_kstar_histograms(registry: HistogramRegistry):
    """
        Register the kstar histograms in the registry.
    """
    for name, entry in KSTAR_HISTOGRAMS.items():
        registry.register(entry)


def register_kstar_matter_histograms(registry: HistogramRegistry):
    """
        Register the kstar matter histograms in the registry.
    """
    for name, _entry in KSTAR_HISTOGRAMS.items():
        entry = deepcopy(_entry)
        entry.name = entry.name + 'Matter'
        entry.save_directory = 'kstarMatter'
        entry.condition += " && fSignedPtHe3 > 0"
        registry.register(entry)

def register_kstar_antimatter_histograms(registry: HistogramRegistry):
    """
        Register the kstar antimatter histograms in the registry.
    """
    for name, _entry in KSTAR_HISTOGRAMS.items():
        entry = deepcopy(_entry)
        entry.name = entry.name + 'Antimatter'
        entry.save_directory = 'kstarAntimatter'
        entry.condition += " && fSignedPtHe3 < 0"
        registry.register(entry)