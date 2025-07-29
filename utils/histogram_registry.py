from dataclasses import dataclass
from ROOT import RDataFrame, TFile

@dataclass
class RegistryEntry:

    name: str
    title: str
    xvar: str
    nbinsx: int
    xmin: float
    xmax: float
    yvar: str = ''
    nbinsy: int = 0
    ymin: float = 0.0
    ymax: float = 0.0
    condition: str = 'True'
    save_directory: str = ''

class HistogramRegistry:
    
    def __init__(self):
        self._registry = {}
        self._registry_entries = {}

    def register(self, entry: RegistryEntry):
        if entry.name in self._registry_entries:
            raise ValueError(f"Histogram '{entry.name}' is already registered.")
        self._registry_entries[entry.name] = entry

    def draw_histogram(self, rdf: RDataFrame):
        
        for entry in self._registry_entries.values():

            if entry.nbinsy > 0:
                hist = rdf.Filter(entry.condition).Histo2D((entry.name, entry.title, 
                                                            entry.nbinsx, entry.xmin, entry.xmax, 
                                                            entry.nbinsy, entry.ymin, entry.ymax), 
                                                            entry.xvar, entry.yvar)
            else:
                hist = rdf.Filter(entry.condition).Histo1D((entry.name, entry.title, 
                                                            entry.nbinsx, entry.xmin, entry.xmax), 
                                                            entry.xvar)
            self._registry[entry.name] = hist

    def prepare_directories(self, output_file: TFile):
        if not self._registry_entries:
            raise ValueError("No registry entries to prepare directories for.")
        
        seen_directories = set()
        for entry in self._registry_entries.values():
            if entry.save_directory and entry.save_directory not in seen_directories and entry.save_directory != '':
                output_file.mkdir(entry.save_directory)
                seen_directories.add(entry.save_directory)
                print(f"Created directory: {entry.save_directory}")

    def save_histograms(self, output_file: TFile):

        print("\nSaving histograms to output file")
        if not self._registry:
            raise ValueError("No histograms to save.")
        
        entry_names_per_dir = {}
        for entry in self._registry_entries.values():
            idir = entry.save_directory or ''
            entry_names_per_dir.setdefault(idir, []).append(entry)

        for idir, entries in entry_names_per_dir.items():
            if idir:
                output_file.cd(idir)
            else:
                output_file.cd()

            for entry in entries:
                hist = self._registry.get(entry.name)
                if hist:
                    hist.Write(entry.name)
                    print(f"Saved histogram: {idir}:{entry.name}")
                else:
                    print(f"Histogram {entry.name} not found in registry.")

