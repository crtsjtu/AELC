# Codes for "A randomized controlled trial for response of microbiome network to exercise and diet intervention in patients with nonalcoholic fatty liver disease".
This repository is the collection of codes for our paper: A randomized controlled trial for response of microbiome network to exercise and diet intervention in patients with nonalcoholic fatty liver disease
## Dependencies
In order to run the scripts, you need to install Matlab R2021a+, Cytoscape 3.9.0+, Python 3.0+ and R 4.0+. 
## Construction of co-occurrence network
To reproduce the results in paper, you should download **cag network.xlsx** and **sample annotation.xlsx**. Then, run the script **corr_cag.m**. This will generate correlation network nodes variable and coefficient variable. A sample network file **cag network.xlsx** can be visualized by Cytoscape.
## Construction of single SparCC network
* step1. Generating reference abundance table for each sample. (**process_phyloseq.R**)
* step2. Constructing sample specific network(**sSparsingle19add1_perm.py**)
* step3. Extract graph file (gml) (**sSparsingle19add1_perm.py**)
