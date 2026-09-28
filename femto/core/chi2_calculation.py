from typing import List
import numpy as np

from ROOT import TH1F, Math

def integrate_hist_density(hist, x_min, x_max):
    if x_max <= x_min:
        return 0.0

    total = 0.0
    for i_bin in range(1, hist.GetNbinsX() + 1):
        bin_low = hist.GetBinLowEdge(i_bin)
        bin_high = bin_low + hist.GetBinWidth(i_bin)
        overlap = min(bin_high, x_max) - max(bin_low, x_min)
        if overlap > 0.0:
            total += hist.GetBinContent(i_bin) * overlap
    return total / (x_max - x_min)


def bin_averaged_hist_values(hist, reference_hist, scale):
    values = []
    for i_bin in range(1, reference_hist.GetNbinsX() + 1):
        x_min = reference_hist.GetBinLowEdge(i_bin)
        x_max = x_min + reference_hist.GetBinWidth(i_bin)
        values.append(scale * integrate_hist_density(hist, x_min, x_max))
    return np.array(values, dtype=float)


def split_uncertainty_source(unc, corr_fraction):
    if corr_fraction < 0.0 or corr_fraction > 1.0:
        raise ValueError(f"Correlation fraction must be in [0, 1], got {corr_fraction}")
    corr_unc = np.sqrt(corr_fraction) * unc
    uncorr_unc = np.sqrt(1.0 - corr_fraction) * unc
    return corr_unc, uncorr_unc



def covariance_chi2(residuals: np.ndarray, data_unc: np.ndarray, sources: List[np.ndarray], 
                    x_maxs: np.ndarray=None, x_max_chi2: float=None):
    cov = np.diag(data_unc * data_unc)
    for source in sources:
        cov += np.outer(source, source)
    cov_inv = np.linalg.pinv(cov, hermitian=True)
    weighted_residuals = cov_inv @ residuals
    chi2_contrib = residuals * weighted_residuals
    
    chi2 = 0
    if x_maxs is None:
        chi2 = float(np.sum(chi2_contrib))
    else:
        for ix, x_max in enumerate(x_maxs):
            if x_max > x_max_chi2:
                break
            chi2 += chi2_contrib[ix]
    
    return chi2, chi2_contrib, cov

class Chi2Calculation:
    
    
    def __init__(self, hist_stat:TH1F, hist_syst:TH1F, model_nominal:TH1F, 
                 model_variations:List[List[TH1F]], model_signal:TH1F=None):
        
        self.kstar_min = None
        self.kstar_max = None
        
        self.data = None
        self.stat_unc = None
        self.syst_unc = None
        self.data_unc = None
        self.model_signal = None
        
        self.model_nominal = None
        #self.model_up = None
        #self.model_down = None
        self.model_unc = None
        self.model_unc_corr = None
        self.model_unc_uncorr = None
        
        self.chi2_contrib = None # chi2 contribution per bin
        
        self.load_data(hist_stat, hist_syst)
        self.load_model(model_nominal, model_variations, model_signal, hist_stat)
        
    def load_data(self, hist_stat:TH1F, hist_syst:TH1F):
        
        self.kstar_min = np.array(
            [hist_stat.GetBinLowEdge(ibin) for ibin in range(1, hist_stat.GetNbinsX() + 1)],
            dtype=float
        )
        self.kstar_max = np.array(
            [hist_stat.GetBinLowEdge(ibin+1) for ibin in range(1, hist_stat.GetNbinsX() + 1)],
            dtype=float
        )
        self.data = np.array(
            [hist_stat.GetBinContent(ibin) for ibin in range(1, hist_stat.GetNbinsX() + 1)],
            dtype=float
        )
        self.stat_unc = np.array(
            [hist_stat.GetBinError(ibin) for ibin in range(1, hist_stat.GetNbinsX() + 1)],
            dtype=float
        )
        self.syst_unc = np.array(
            [hist_syst.GetBinError(ibin) for ibin in range(1, hist_syst.GetNbinsX() + 1)],
            dtype=float
        ) if hist_syst is not None else np.zeros_like(self.stat_unc)
        self.data_unc = np.sqrt(self.stat_unc**2 + self.syst_unc**2)
    
    def load_model(self, model_nominal:TH1F, model_variations:List[List[TH1F]], model_signal:TH1F, hist_ref:TH1F):
        '''
        
        '''
        
        self.model_up, self.model_down, self.model_unc = [], [], []
        for ivariation, variations in enumerate(model_variations):
            hist_model_up = model_nominal.Clone(f"hist_model_up_{ivariation}")
            hist_model_down = model_nominal.Clone(f"hist_model_down_{ivariation}")
            for ibin in range(1, model_nominal.GetNbinsX() + 1):
                y_value_up = max([variation.GetBinContent(ibin) for variation in variations]) if len(variations) > 0 else model_nominal.GetBinContent(ibin)
                y_value_down = min([variation.GetBinContent(ibin) for variation in variations]) if len(variations) > 0 else model_nominal.GetBinContent(ibin)
                hist_model_up.SetBinContent(ibin, y_value_up)
                hist_model_down.SetBinContent(ibin, y_value_down)
            
            model_up = bin_averaged_hist_values(hist_model_up, hist_ref, scale=1.0)
            model_down = bin_averaged_hist_values(hist_model_down, hist_ref, scale=1.0)
            model_unc = 0.5 * (model_up - model_down)
            #self.model_up.append(model_up)
            #self.model_down.append(model_down)
            self.model_unc.append(model_unc)

        self.model_nominal = bin_averaged_hist_values(model_nominal, hist_ref, scale=1.0)
        if model_signal is not None:
            self.model_signal = bin_averaged_hist_values(model_signal, hist_ref, scale=1.0)
            self.model_nominal += self.model_signal
    
    def calculate_chi2(self, correlation_fractions:List[float]=None, x_max_chi2:float=None):
        
        self.residuals = self.data - self.model_nominal
        if correlation_fractions is None:
            correlation_fractions = [1.0] * len(self.model_unc)
        self.model_unc_corr, self.model_unc_uncorr = [], []
        
        ### print('DEBUG: kstar values:', self.kstar_min)
        ### print('DEBUG: model_nominal:', self.model_nominal)
        
        for i_source, (model_unc, corr_fraction) in enumerate(zip(self.model_unc, correlation_fractions)):
            model_unc_corr, model_unc_uncorr = split_uncertainty_source(model_unc, corr_fraction)
            self.model_unc_corr.append(model_unc_corr)
            self.model_unc_uncorr.append(model_unc_uncorr)

        uncorr_unc = self.data_unc**2
        for model_unc_uncorr in self.model_unc_uncorr:
            uncorr_unc += model_unc_uncorr**2
        uncorr_unc = np.sqrt(uncorr_unc)
        
        x_max_chi2 = x_max_chi2 if x_max_chi2 is not None else self.kstar_max[-1]
        chi2, self.chi2_contrib, _ = covariance_chi2(self.residuals, uncorr_unc, self.model_unc_corr, x_maxs=self.kstar_max, x_max_chi2=x_max_chi2)
        ndf = np.where(self.kstar_max <= x_max_chi2)[0].size - 1
        p_value = float(Math.chisquared_cdf_c(chi2, ndf))

        return chi2, p_value, self.chi2_contrib        
    
    def write_histograms(self, outfile, suffix:str=''):
        nbins = len(self.data)
        kstar_min = self.kstar_min[0]
        kstar_max = self.kstar_max[-1]
        
        h_chi2 = TH1F(f'chi2{suffix}', ';#it{k}* (GeV/#it{c}); #chi^{2}', nbins, kstar_min, kstar_max)
        h_nsigma = TH1F(f'nsigma{suffix}', ';#it{k}* (GeV/#it{c}); #it{n}#sigma', nbins, kstar_min, kstar_max)
        h_chi2_contrib = TH1F(f'chi2_contrib{suffix}', '#chi^{2} bin by bin;#it{k}* (GeV/#it{c}); #chi^{2} contribution', nbins, kstar_min, kstar_max)
        chi2 = 0
        for ibin in range(nbins):
            nsigma = np.sqrt(self.chi2_contrib[ibin]) * np.sign(self.residuals[ibin]) if self.chi2_contrib[ibin] > 0 else 0.0
            h_nsigma.SetBinContent(ibin + 1, nsigma)
            h_chi2_contrib.SetBinContent(ibin + 1, self.chi2_contrib[ibin])
            chi2 += self.chi2_contrib[ibin]
            h_chi2.SetBinContent(ibin + 1, chi2)
            
        outfile.cd()
        for hist in [h_chi2, h_nsigma, h_chi2_contrib]:
            hist.Write()
        
        