# Preparation
Data preparation for the specific studies that follow. 
In this folder the code used to create the interesting output histograms and the QA outputs is stored.

## prepare.py
Reads trees with a RDataFrame and creates QA histograms and histograms of interest used to compute the correlation function in several centrality classes. Requires the input files to contain a single TTree (several trees with different branches cannot be handled). A folder structure with several trees with the same branches is accepted.
Configuration for data selections and input files can be set in a yaml configuration file (e.g. config/config_prepare.yml).

## merge_trees.py
Reads trees with different branches and returns a unique tree that stacks all the original branches horizontally. It can loop over several folders containing the trees to merge. Prepares the input for prepare.py.

## correlation.py
Computes the correlation function from the input kstar distribution prepared in prepare.py.

## systematics.py
Runs the analysis several times creating the kstar distribution and computing the correlation function in order to compute the systematics. It can run either varying all chosen variables or by varying them one at a time.
