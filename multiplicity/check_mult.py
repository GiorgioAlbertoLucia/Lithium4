
from array import array
from math import sqrt
import ROOT

input_centrality = [50, 80]  # Example input centrality range

cent_bins = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
mult_vals = [2047, 1668, 1253, 848, 559, 351, 205, 110, 53, 23.2]
uncertainties = [54, 42, 33, 25, 19, 14, 11, 8, 5, 2.8]

th1_cent_mult = ROOT.TH1D(
	"hCentMult",
	"Centrality multiplicity;Centrality (%);Multiplicity",
	len(mult_vals),
	array("d", cent_bins),
)

for bin_index, (mult_val, unc) in enumerate(zip(mult_vals, uncertainties), start=1):
	th1_cent_mult.SetBinContent(bin_index, mult_val)
	th1_cent_mult.SetBinError(bin_index, unc)

file_centrality_pairs = ROOT.TFile("centrality_pairs.root")
th1_cent = file_centrality_pairs.Get("hCentrality")


bin_low_mult = th1_cent_mult.FindBin(input_centrality[0])
bin_high_mult = th1_cent_mult.FindBin(input_centrality[1])

mult_average = 0
mult_average_error = 0
weight_sum = 0
for bin_index in range(bin_low_mult, bin_high_mult):
	print("------------------------------")
	print("Centrality: ", th1_cent_mult.GetBinLowEdge(bin_index), "-", th1_cent_mult.GetBinLowEdge(bin_index + 1))
	mult_value = th1_cent_mult.GetBinContent(bin_index)
	mult_error = th1_cent_mult.GetBinError(bin_index)
	cent_low = th1_cent_mult.GetBinLowEdge(bin_index)
	cent_high = th1_cent_mult.GetBinLowEdge(bin_index + 1)
	weight = th1_cent.Integral(th1_cent.FindBin(cent_low), th1_cent.FindBin(cent_high - 0.0001))
	print("Weight: ", weight)
	mult_average += mult_value * weight
	mult_average_error += (mult_error * weight) ** 2
	weight_sum += weight

if weight_sum > 0:
	mult_average /= weight_sum
	mult_average_error = sqrt(mult_average_error) / weight_sum
 
print(f"Average multiplicity for centrality range {input_centrality[0]}-{input_centrality[1]}: {mult_average} ± {mult_average_error}")
