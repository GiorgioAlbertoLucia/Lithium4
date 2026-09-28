import numpy as np

from ROOT import TFile, TH1F

from torchic.core.histogram import load_hist, build_efficiency
from torchic.utils.root import set_root_object

def get_weighted_average(hist, weight_hist):
    weighted_sum = 0.0
    total_weight = 0.0
    
    for i in range(1, hist.GetNbinsX() + 1):
        bin_content = hist.GetBinContent(i)
        weight = weight_hist.GetBinContent(i)
        weighted_sum += bin_content * weight
        total_weight += weight
    
    return weighted_sum / total_weight if total_weight != 0 else 0.0

def compute_event_loss(h_gen_events_vs_nch, output_file, h_nch_centralities: list, h_nch_1050: TH1F):
    
    h_gen_events = h_gen_events_vs_nch.ProjectionX("hGenEvents", 1, 1)
    h_rec_at_least_one_and_passed_evsel = h_gen_events_vs_nch.ProjectionX("hRecAtLeastOneAndPassedEvSel", 2, 2)
    
    h_event_loss_vs_nch = h_rec_at_least_one_and_passed_evsel.Clone("hEventLossVsNch")
    h_event_loss_vs_nch.Divide(h_gen_events)
    
    event_loss_010 = get_weighted_average(h_event_loss_vs_nch, h_nch_centralities[0])
    gen_events_010 = get_weighted_average(h_gen_events, h_nch_centralities[0])
    event_loss_error_010 = np.sqrt(event_loss_010 * (1 - event_loss_010) / gen_events_010) if (gen_events_010 != 0 and event_loss_010 < 1) else 0.0
    event_loss_1030 = get_weighted_average(h_event_loss_vs_nch, h_nch_centralities[1])
    gen_events_1030 = get_weighted_average(h_gen_events, h_nch_centralities[1])
    event_loss_error_1030 = np.sqrt(event_loss_1030 * (1 - event_loss_1030) / gen_events_1030) if (gen_events_1030 != 0 and event_loss_1030 < 1) else 0.0
    event_loss_3050 = get_weighted_average(h_event_loss_vs_nch, h_nch_centralities[2])
    gen_events_3050 = get_weighted_average(h_gen_events, h_nch_centralities[2])
    event_loss_error_3050 = np.sqrt(event_loss_3050 * (1 - event_loss_3050) / gen_events_3050) if (gen_events_3050 != 0 and event_loss_3050 < 1) else 0.0
    event_loss_5080 = get_weighted_average(h_event_loss_vs_nch, h_nch_centralities[3])
    gen_events_5080 = get_weighted_average(h_gen_events, h_nch_centralities[3])
    event_loss_error_5080 = np.sqrt(event_loss_5080 * (1 - event_loss_5080) / gen_events_5080) if (gen_events_5080 != 0 and event_loss_5080 < 1) else 0.0
    event_loss_1050 = get_weighted_average(h_event_loss_vs_nch, h_nch_1050)
    gen_events_1050 = get_weighted_average(h_gen_events, h_nch_1050)
    event_loss_error_1050 = np.sqrt(event_loss_1050 * (1 - event_loss_1050) / gen_events_1050) if (gen_events_1050 != 0 and event_loss_1050 < 1) else 0.0
    
    h_event_loss_vs_centrality = TH1F("hEventLossVsCentrality", "Event loss vs Centrality; FT0C Centrality (%); Event loss", 10, 0., 100.)
    for i in range(1, 11):
        h_nch_centrality = h_nch_centralities[i-1]
        event_loss_cent_interval = get_weighted_average(h_event_loss_vs_nch, h_nch_centrality)
        gen_events_cent_interval = get_weighted_average(h_gen_events, h_nch_centrality)
        h_event_loss_vs_centrality.SetBinContent(i, event_loss_cent_interval)
        h_event_loss_vs_centrality.SetBinError(i, np.sqrt(event_loss_cent_interval * (1 - event_loss_cent_interval) / gen_events_cent_interval) if (gen_events_cent_interval != 0 and event_loss_cent_interval < 1) else 0.0)
        #h_event_loss_vs_centrality.SetBinError(i, 1e-12)

    h_event_loss_vs_centrality_analysis_extended = TH1F("hEventLossVsCentralityAnalysisExtended", "Event loss vs Centrality; FT0C Centrality (%); Event loss", 5, -0.5, 6.5)
    h_event_loss_vs_centrality_analysis_extended.SetBinContent(1, event_loss_010)
    h_event_loss_vs_centrality_analysis_extended.SetBinContent(2, event_loss_1030)
    h_event_loss_vs_centrality_analysis_extended.SetBinContent(3, event_loss_3050)
    h_event_loss_vs_centrality_analysis_extended.SetBinContent(4, event_loss_5080)
    h_event_loss_vs_centrality_analysis_extended.SetBinContent(5, event_loss_1050)
    h_event_loss_vs_centrality_analysis_extended.SetBinError(1, event_loss_error_010)
    h_event_loss_vs_centrality_analysis_extended.SetBinError(2, event_loss_error_1030)
    h_event_loss_vs_centrality_analysis_extended.SetBinError(3, event_loss_error_3050)
    h_event_loss_vs_centrality_analysis_extended.SetBinError(4, event_loss_error_5080)
    h_event_loss_vs_centrality_analysis_extended.SetBinError(5, event_loss_error_1050)
    #for i in range(1, 6):
    #    h_event_loss_vs_centrality_analysis_extended.SetBinError(i, 1e-12)
    h_event_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(1, "0-10%")
    h_event_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(2, "10-30%")
    h_event_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(3, "30-50%")
    h_event_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(4, "50-80%")
    h_event_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(5, "10-50%")
    
    h_event_loss_vs_centrality_analysis = TH1F("hEventLossVsCentralityAnalysis", "Event loss vs Centrality; FT0C Centrality (%); Event loss", 2, -0.5, 1.5)
    h_event_loss_vs_centrality_analysis.SetBinContent(1, event_loss_010)
    h_event_loss_vs_centrality_analysis.SetBinContent(2, event_loss_1050)
    h_event_loss_vs_centrality_analysis.SetBinError(1, event_loss_error_010)
    h_event_loss_vs_centrality_analysis.SetBinError(2, event_loss_error_1050)
    #for i in range(1, 3):
    #    h_event_loss_vs_centrality_analysis.SetBinError(i, 1e-12)
    h_event_loss_vs_centrality_analysis.GetXaxis().SetBinLabel(1, "0-10%")
    h_event_loss_vs_centrality_analysis.GetXaxis().SetBinLabel(2, "10-50%")
    
    set_root_object(h_event_loss_vs_nch, title=f'; d#it{{N}}_{{ch}} / d#it{{#eta}} |_{{#it{{#eta}}|<0.5}}; Event loss')
    output_dir = output_file.mkdir('EventLoss')
    output_dir.cd()
    h_gen_events.Write()
    h_rec_at_least_one_and_passed_evsel.Write()
    h_event_loss_vs_nch.Write()
    h_event_loss_vs_centrality.Write()
    h_event_loss_vs_centrality_analysis.Write()
    h_event_loss_vs_centrality_analysis_extended.Write()
                    
    return h_event_loss_vs_centrality_analysis

def compute_inverse_event_splitting(h_rec_cent_vs_nch, h_gen_rec_at_least_one_and_passed_evsel_vs_nch, output_file):
    
    h_rec_cent = h_rec_cent_vs_nch.ProjectionX("hRecCent", 1, h_rec_cent_vs_nch.GetNbinsX()+1)
    h_gen_rec_at_least_one_and_passed_evsel = h_gen_rec_at_least_one_and_passed_evsel_vs_nch.ProjectionX("hGenRecAtLeastOneAndPassedEvSel", 1, h_gen_rec_at_least_one_and_passed_evsel_vs_nch.GetNbinsX()+1)
    
    h_rec_010 = h_rec_cent.Integral(1, 10)
    h_rec_1050 = h_rec_cent.Integral(11, 50)
    
    h_gen_rec_at_least_one_and_passed_evsel_010 = h_gen_rec_at_least_one_and_passed_evsel.Integral(1, 10)
    h_gen_rec_at_least_one_and_passed_evsel_1050 = h_gen_rec_at_least_one_and_passed_evsel.Integral(11, 50)
    
    event_splitting_010 =  h_gen_rec_at_least_one_and_passed_evsel_010 / h_rec_010 if h_rec_010 != 0 else 0
    event_splitting_error_010 = np.sqrt(event_splitting_010 * (1 - event_splitting_010) / h_gen_rec_at_least_one_and_passed_evsel_010) if (h_gen_rec_at_least_one_and_passed_evsel_010 != 0 and event_splitting_010 < 1) else 0.0
    inverse_event_splitting_010 = 1 / event_splitting_010 if event_splitting_010 != 0 else 0
    inverse_event_splitting_error_010 = event_splitting_error_010 / (event_splitting_010 ** 2) if event_splitting_010 != 0 else 0.0
    event_splitting_1030 = h_gen_rec_at_least_one_and_passed_evsel.Integral(11, 30) / h_rec_cent.Integral(11, 30) if h_rec_cent.Integral(11, 30) != 0 else 0
    event_splitting_error_1030 = np.sqrt(event_splitting_1030 * (1 - event_splitting_1030) / h_gen_rec_at_least_one_and_passed_evsel.Integral(11, 30)) if (h_gen_rec_at_least_one_and_passed_evsel.Integral(11, 30) != 0 and event_splitting_1030 < 1) else 0.0
    inverse_event_splitting_1030 = 1 / event_splitting_1030 if event_splitting_1030 != 0 else 0
    inverse_event_splitting_error_1030 = event_splitting_error_1030 / (event_splitting_1030 ** 2) if event_splitting_1030 != 0 else 0.0
    event_splitting_3050 = h_gen_rec_at_least_one_and_passed_evsel.Integral(31, 50) / h_rec_cent.Integral(31, 50) if h_rec_cent.Integral(31, 50) != 0 else 0
    event_splitting_error_3050 = np.sqrt(event_splitting_3050 * (1 - event_splitting_3050) / h_gen_rec_at_least_one_and_passed_evsel.Integral(31, 50)) if (h_gen_rec_at_least_one_and_passed_evsel.Integral(31, 50) != 0 and event_splitting_3050 < 1) else 0.0
    inverse_event_splitting_3050 = 1 / event_splitting_3050 if event_splitting_3050 != 0 else 0
    inverse_event_splitting_error_3050 = event_splitting_error_3050 / (event_splitting_3050 ** 2) if event_splitting_3050 != 0 else 0.0
    event_splitting_5080 = h_gen_rec_at_least_one_and_passed_evsel.Integral(51, 80) / h_rec_cent.Integral(51, 80) if h_rec_cent.Integral(51, 80) != 0 else 0
    event_splitting_error_5080 = np.sqrt(event_splitting_5080 * (1 - event_splitting_5080) / h_gen_rec_at_least_one_and_passed_evsel.Integral(51, 80)) if (h_gen_rec_at_least_one_and_passed_evsel.Integral(51, 80) != 0 and event_splitting_5080 < 1) else 0.0
    inverse_event_splitting_5080 = 1 / event_splitting_5080 if event_splitting_5080 != 0 else 0
    inverse_event_splitting_error_5080 = event_splitting_error_5080 / (event_splitting_5080 ** 2) if event_splitting_5080 != 0 else 0.0
    event_splitting_1050 =  h_gen_rec_at_least_one_and_passed_evsel_1050 / h_rec_1050 if h_rec_1050 != 0 else 0
    event_splitting_error_1050 = np.sqrt(event_splitting_1050 * (1 - event_splitting_1050) / h_gen_rec_at_least_one_and_passed_evsel_1050) if (h_gen_rec_at_least_one_and_passed_evsel_1050 != 0 and event_splitting_1050 < 1) else 0.0
    inverse_event_splitting_1050 = 1 / event_splitting_1050 if event_splitting_1050 != 0 else 0
    inverse_event_splitting_error_1050 = event_splitting_error_1050 / (event_splitting_1050 ** 2) if event_splitting_1050 != 0 else 0.0
    
    h_inverse_event_splitting_vs_centrality = TH1F("hEventSplittingVsCentrality", "Event splitting vs Centrality; FT0C Centrality (%); 1 / Event splitting", 10, 0., 100.)
    for i in range(1, 11):
        h_rec_cent_interval = h_rec_cent.Integral(i, 10*i)
        h_gen_rec_at_least_one_and_passed_evsel_interval = h_gen_rec_at_least_one_and_passed_evsel.Integral(i, 10*i)
        inverse_event_splitting_cent_interval = h_rec_cent_interval / h_gen_rec_at_least_one_and_passed_evsel_interval if h_gen_rec_at_least_one_and_passed_evsel_interval != 0 else 0
        event_splitting_cent_interval = 1 / inverse_event_splitting_cent_interval if inverse_event_splitting_cent_interval != 0 else 0
        event_splitting_error_cent_interval = np.sqrt(event_splitting_cent_interval * (1 - event_splitting_cent_interval) / h_gen_rec_at_least_one_and_passed_evsel_interval) if (h_gen_rec_at_least_one_and_passed_evsel_interval != 0 and event_splitting_cent_interval < 1) else 0.0
        inverse_event_splitting_error_cent_interval = event_splitting_error_cent_interval / (event_splitting_cent_interval ** 2) if inverse_event_splitting_cent_interval != 0 else 0.0
        h_inverse_event_splitting_vs_centrality.SetBinContent(i, inverse_event_splitting_cent_interval)
        h_inverse_event_splitting_vs_centrality.SetBinError(i, inverse_event_splitting_error_cent_interval)
        #h_inverse_event_splitting_vs_centrality.SetBinError(i, 1e-12)

    h_inverse_event_splitting_vs_centrality_analysis_extended = TH1F("hEventSplittingVsCentralityAnalysisExtended", "Event splitting vs Centrality; FT0C Centrality (%); 1 / Event splitting", 5, -0.5, 6.5)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinContent(1, inverse_event_splitting_010)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinContent(2, inverse_event_splitting_1030)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinContent(3, inverse_event_splitting_3050)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinContent(4, inverse_event_splitting_5080)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinContent(5, inverse_event_splitting_1050)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(1, inverse_event_splitting_error_010)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(2, inverse_event_splitting_error_1030)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(3, inverse_event_splitting_error_3050)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(4, inverse_event_splitting_error_5080)
    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(5, inverse_event_splitting_error_1050)
    #for i in range(1, 6):
    #    h_inverse_event_splitting_vs_centrality_analysis_extended.SetBinError(i, 1e-12)
    h_inverse_event_splitting_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(1, "0-10%")
    h_inverse_event_splitting_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(2, "10-30%")
    h_inverse_event_splitting_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(3, "30-50%")
    h_inverse_event_splitting_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(4, "50-80%")
    h_inverse_event_splitting_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(5, "10-50%")

    h_inverse_event_splitting_vs_centrality_analysis = TH1F("hEventSplittingVsCentralityAnalysis", "Event splitting vs Centrality; FT0C Centrality (%); 1 / Event splitting", 2, -0.5, 1.5)
    h_inverse_event_splitting_vs_centrality_analysis.SetBinContent(1, inverse_event_splitting_010)
    h_inverse_event_splitting_vs_centrality_analysis.SetBinContent(2, inverse_event_splitting_1050)
    h_inverse_event_splitting_vs_centrality_analysis.SetBinError(1, inverse_event_splitting_error_010)
    h_inverse_event_splitting_vs_centrality_analysis.SetBinError(2, inverse_event_splitting_error_1050)
    #for i in range(1, 3):
    #    h_inverse_event_splitting_vs_centrality_analysis.SetBinError(i, 1e-12)
    h_inverse_event_splitting_vs_centrality_analysis.GetXaxis().SetBinLabel(1, "0-10%")
    h_inverse_event_splitting_vs_centrality_analysis.GetXaxis().SetBinLabel(2, "10-50%")
                                                           

    output_dir = output_file.mkdir('EventSplitting')
    output_dir.cd()
    h_rec_cent.Write()
    h_gen_rec_at_least_one_and_passed_evsel.Write()
    h_inverse_event_splitting_vs_centrality.Write()
    h_inverse_event_splitting_vs_centrality_analysis.Write()
    h_inverse_event_splitting_vs_centrality_analysis_extended.Write()
                    
    return h_inverse_event_splitting_vs_centrality

def compute_event_correction(h_event_loss_vs_centrality, h_inverse_event_splitting_vs_centrality, output_file):
    
    h_event_correction_vs_centrality_extended = TH1F("hEventCorrectionVsCentralityExtended", "Event correction vs Centrality; FT0C Centrality (%); Event correction", 5, -0.5, 6.5)
    
    for i in range(1, 6):
        event_loss = h_event_loss_vs_centrality.GetBinContent(i)
        event_loss_error = h_event_loss_vs_centrality.GetBinError(i)
        inverse_event_splitting = h_inverse_event_splitting_vs_centrality.GetBinContent(i)
        inverse_event_splitting_error = h_inverse_event_splitting_vs_centrality.GetBinError(i)
        event_correction = event_loss * inverse_event_splitting
        event_correction_error = np.sqrt((inverse_event_splitting * event_loss_error)**2 + (event_loss * inverse_event_splitting_error)**2)
        h_event_correction_vs_centrality_extended.SetBinContent(i, event_correction)
        h_event_correction_vs_centrality_extended.SetBinError(i, event_correction_error)
        #h_event_correction_vs_centrality_extended.SetBinError(i, 1e-12)

    h_event_correction_vs_centrality_extended.GetXaxis().SetBinLabel(1, "0-10%")
    h_event_correction_vs_centrality_extended.GetXaxis().SetBinLabel(2, "10-30%")
    h_event_correction_vs_centrality_extended.GetXaxis().SetBinLabel(3, "30-50%")
    h_event_correction_vs_centrality_extended.GetXaxis().SetBinLabel(4, "50-80%")
    h_event_correction_vs_centrality_extended.GetXaxis().SetBinLabel(5, "10-50%")
    
    h_event_correction_vs_centrality = TH1F("hEventCorrectionVsCentrality", "Event correction vs Centrality; FT0C Centrality (%); Event correction", 2, -0.5, 1.5)
    for i in range(1, 3):
        event_loss = h_event_loss_vs_centrality.GetBinContent(i)
        event_loss_error = h_event_loss_vs_centrality.GetBinError(i)
        inverse_event_splitting = h_inverse_event_splitting_vs_centrality.GetBinContent(i)
        inverse_event_splitting_error = h_inverse_event_splitting_vs_centrality.GetBinError(i)
        event_correction = event_loss * inverse_event_splitting
        event_correction_error = np.sqrt((inverse_event_splitting * event_loss_error)**2 + (event_loss * inverse_event_splitting_error)**2)
        h_event_correction_vs_centrality.SetBinContent(i, event_correction)
        h_event_correction_vs_centrality.SetBinError(i, event_correction_error)
        #h_event_correction_vs_centrality.SetBinError(i, 1e-12)
    
    h_event_correction_vs_centrality.GetXaxis().SetBinLabel(1, "0-10%")
    h_event_correction_vs_centrality.GetXaxis().SetBinLabel(2, "10-50%")

    output_dir = output_file.mkdir('EventCorrection')
    output_dir.cd()
    h_event_correction_vs_centrality_extended.Write()
    h_event_correction_vs_centrality.Write()
    
    return h_event_correction_vs_centrality_extended

def compute_signal_loss(h_gen_signal_passed_ev_sel_pt_vs_nch, h_gen_signal_pt_vs_nch, h_nch_centralities, output_file):
    
    h_gen_signal_passed_ev_sel_nch = h_gen_signal_passed_ev_sel_pt_vs_nch.ProjectionY("hGenSignalPassedEvSelNch")
    h_gen_signal_nch = h_gen_signal_pt_vs_nch.ProjectionY("hGenSignalNch")
    
    gen_signal_passed_ev_sel_010 = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centralities[0])
    gen_signal_010 = get_weighted_average(h_gen_signal_nch, h_nch_centralities[0])
    signal_loss_010 = gen_signal_passed_ev_sel_010 / gen_signal_010 if gen_signal_010 != 0 else 0
    signal_loss_error_010 = np.sqrt(signal_loss_010 * (1 - signal_loss_010) / gen_signal_010) if (gen_signal_010 != 0 and signal_loss_010 < 1) else 0.0
    gen_signal_passed_ev_sel_1030 = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centralities[1])
    gen_signal_1030 = get_weighted_average(h_gen_signal_nch, h_nch_centralities[1])
    signal_loss_1030 = gen_signal_passed_ev_sel_1030 / gen_signal_1030 if gen_signal_1030 != 0 else 0
    signal_loss_error_1030 = np.sqrt(signal_loss_1030 * (1 - signal_loss_1030) / gen_signal_1030) if (gen_signal_1030 != 0 and signal_loss_1030 < 1) else 0.0
    gen_signal_passed_ev_sel_3050 = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centralities[2])
    gen_signal_3050 = get_weighted_average(h_gen_signal_nch, h_nch_centralities[2])
    signal_loss_3050 = gen_signal_passed_ev_sel_3050 / gen_signal_3050 if gen_signal_3050 != 0 else 0
    signal_loss_error_3050 = np.sqrt(signal_loss_3050 * (1 - signal_loss_3050) / gen_signal_3050) if (gen_signal_3050 != 0 and signal_loss_3050 < 1) else 0.0
    gen_signal_passed_ev_sel_5080 = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centralities[3])
    gen_signal_5080 = get_weighted_average(h_gen_signal_nch, h_nch_centralities[3])
    signal_loss_5080 = gen_signal_passed_ev_sel_5080 / gen_signal_5080 if gen_signal_5080 != 0 else 0
    signal_loss_error_5080 = np.sqrt(signal_loss_5080 * (1 - signal_loss_5080) / gen_signal_5080) if (gen_signal_5080 != 0 and signal_loss_5080 < 1) else 0.0
    gen_signal_passed_ev_sel_1050 = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centralities[2])
    gen_signal_1050 = get_weighted_average(h_gen_signal_nch, h_nch_centralities[2])
    signal_loss_1050 = gen_signal_passed_ev_sel_1050 / gen_signal_1050 if gen_signal_1050 != 0 else 0
    signal_loss_error_1050 = np.sqrt(signal_loss_1050 * (1 - signal_loss_1050) / gen_signal_1050) if (gen_signal_1050 != 0 and signal_loss_1050 < 1) else 0.0
    
    h_signal_loss_vs_centrality = TH1F("hSignalLossVsCentrality", "Signal loss vs Centrality; FT0C Centrality (%); Signal loss", 10, 0., 100.)
    for i in range(1, 11):
        h_nch_centrality = h_nch_centralities[i-1]
        signal_loss_cent_interval = get_weighted_average(h_gen_signal_passed_ev_sel_nch, h_nch_centrality) / get_weighted_average(h_gen_signal_nch, h_nch_centrality) if get_weighted_average(h_gen_signal_nch, h_nch_centrality) != 0 else 0
        signal_loss_error_cent_interval = np.sqrt(signal_loss_cent_interval * (1 - signal_loss_cent_interval) / get_weighted_average(h_gen_signal_nch, h_nch_centrality)) if (get_weighted_average(h_gen_signal_nch, h_nch_centrality) != 0 and signal_loss_cent_interval < 1) else 0.0
        h_signal_loss_vs_centrality.SetBinContent(i, signal_loss_cent_interval)
        h_signal_loss_vs_centrality.SetBinError(i, signal_loss_error_cent_interval)
        #h_signal_loss_vs_centrality.SetBinError(i, 1e-12)
        
    h_signal_loss_vs_centrality_analysis_extended = TH1F("hSignalLossVsCentralityAnalysisExtended", "Signal loss vs Centrality; FT0C Centrality (%); Signal loss", 5, -0.5, 6.5)
    h_signal_loss_vs_centrality_analysis_extended.SetBinContent(1, signal_loss_010)
    h_signal_loss_vs_centrality_analysis_extended.SetBinContent(2, signal_loss_1030)
    h_signal_loss_vs_centrality_analysis_extended.SetBinContent(3, signal_loss_3050)
    h_signal_loss_vs_centrality_analysis_extended.SetBinContent(4, signal_loss_5080)
    h_signal_loss_vs_centrality_analysis_extended.SetBinContent(5, signal_loss_1050)
    h_signal_loss_vs_centrality_analysis_extended.SetBinError(1, signal_loss_error_010)
    h_signal_loss_vs_centrality_analysis_extended.SetBinError(2, signal_loss_error_1030)
    h_signal_loss_vs_centrality_analysis_extended.SetBinError(3, signal_loss_error_3050)
    h_signal_loss_vs_centrality_analysis_extended.SetBinError(4, signal_loss_error_5080)
    h_signal_loss_vs_centrality_analysis_extended.SetBinError(5, signal_loss_error_1050)
    #for i in range(1, 6):
    #    h_signal_loss_vs_centrality_analysis_extended.SetBinError(i, 1e-12)
    h_signal_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(1, "0-10%")
    h_signal_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(2, "10-30%")
    h_signal_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(3, "30-50%")
    h_signal_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(4, "50-80%")
    h_signal_loss_vs_centrality_analysis_extended.GetXaxis().SetBinLabel(5, "10-50%")
    
    h_signal_loss_vs_centrality_analysis = TH1F("hSignalLossVsCentralityAnalysis", "Signal loss vs Centrality; FT0C Centrality (%); Signal loss", 2, -0.5, 1.5)
    h_signal_loss_vs_centrality_analysis.SetBinContent(1, signal_loss_010)
    h_signal_loss_vs_centrality_analysis.SetBinContent(2, signal_loss_1050)
    h_signal_loss_vs_centrality_analysis.SetBinError(1, signal_loss_error_010)
    h_signal_loss_vs_centrality_analysis.SetBinError(2, signal_loss_error_1050)
    #for i in range(1, 3):
    #    h_signal_loss_vs_centrality_analysis.SetBinError(i, 1e-12)
    h_signal_loss_vs_centrality_analysis.GetXaxis().SetBinLabel(1, "0-10%")
    h_signal_loss_vs_centrality_analysis.GetXaxis().SetBinLabel(2, "10-50%")
    
    output_dir = output_file.mkdir('SignalLoss')
    output_dir.cd()
    h_gen_signal_passed_ev_sel_nch.Write()
    h_gen_signal_nch.Write()
    h_signal_loss_vs_centrality.Write()
    h_signal_loss_vs_centrality_analysis.Write()
    h_signal_loss_vs_centrality_analysis_extended.Write()
    return h_signal_loss_vs_centrality
    

def convert_nch_to_centrality(h_nch_vs_centrality, output_file):
    
    h_nch_for_centralities = []
    for i in range(1, 11):
        h_nch_for_centrality = h_nch_vs_centrality.ProjectionY(f"hNchForCentrality{i*10}", i, i*10)
        h_nch_for_centralities.append(h_nch_for_centrality)
    
    h_nch_for_centrality_1050 = h_nch_vs_centrality.ProjectionY("hNchForCentrality1050", 11, 50)

    output_dir = output_file.mkdir('NchForCentrality')
    output_dir.cd()
    for ihist, h_nch_for_centrality in enumerate(h_nch_for_centralities):
        h_nch_for_centrality.Write(f'hNchForCentrality{(ihist)*10}{(ihist+1)*10}')
    h_nch_for_centrality_1050.Write()
    
    return h_nch_for_centralities, h_nch_for_centrality_1050

if __name__ == "__main__":
    
    input_files = {
        '2023': '/data/galucia/lithium/event_loss/EventLoss_LHC25g11.root',
        '2024': '/data/galucia/lithium/event_loss/EventLoss_LHC26e5.root',
        '2025': '/data/galucia/lithium/event_loss/EventLoss_LHC26e6.root',
    }
    input_dir = 'he3-hadron-femto/QA/EventLoss'
    
    gen_events_vs_nch = {key: load_hist(value, f"{input_dir}/hGenEventsNchEta05") for key, value in input_files.items()}
    nch_vs_centrality = {key: load_hist(value, f"{input_dir}/hGenCentralityColvsMultiplicityGenEta05") for key, value in input_files.items()}
    rec_cent_vs_nch = {key: load_hist(value, f"{input_dir}/hRecoCentralityColvsMultiplicityRecoEta05") for key, value in input_files.items()}
    gen_rec_at_least_one_and_passed_evsel_vs_nch = {key: load_hist(value, f"{input_dir}/hGenCentralityColvsMultiplicityGenEta05") for key, value in input_files.items()}
    gen_signal_pt_vs_nch = {key: load_hist(value, f"{input_dir}/hGenLi4vsMultiplicityGenEta05BeforeEvtSel") for key, value in input_files.items()}
    gen_signal_passed_ev_sel_pt_vs_nch = {key: load_hist(value, f"{input_dir}/hGenLi4vsMultiplicityGenEta05AfterSel") for key, value in input_files.items()}

    output_file = TFile(f"EventCorrection.root", "RECREATE")

    for year in input_files.keys():
        
        outdir_year = output_file.mkdir(year)
        h_nch_centralities, h_nch_1050 = convert_nch_to_centrality(nch_vs_centrality[year], outdir_year)
        h_event_loss_vs_centrality = compute_event_loss(gen_events_vs_nch[year], outdir_year, h_nch_centralities, h_nch_1050)
        h_inverse_event_splitting_vs_centrality = compute_inverse_event_splitting(rec_cent_vs_nch[year], gen_rec_at_least_one_and_passed_evsel_vs_nch[year], outdir_year)
        h_signal_loss_vs_nch = compute_signal_loss(gen_signal_passed_ev_sel_pt_vs_nch[year], gen_signal_pt_vs_nch[year], h_nch_centralities, outdir_year)
        compute_event_correction(h_event_loss_vs_centrality, h_inverse_event_splitting_vs_centrality, outdir_year)
        
        
    output_file.Close()