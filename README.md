# covid_insertions_paper

This repository contains scripts and data to reproduce the results described in paper "Insertions in SARS-CoV-2 genome caused by template switch and duplications give rise to new variants of potential concern" by Sofya K Garushyants, Igor B Rogozin, and Eugene V. Koonin.

**data** contains all data files:
GISAID_msa_0617_ab_insertions.csv - list of initial inserts and genomes where they were found from the alignment\
GISAID_msa0617_insertions_filtered.csv - output of Inserts_aliJun17_filter_data_and_save.R, this file contains list of filtered inserts that were not verified by alingnment check\
Inserts_msa0617_phylogeny_all.csv - the same data as presented in Supplementary table 2. The final list of verified inserts (8 inserts that were not verified by raw read data analysis are not removed, by marked with 'N' in 'Confirmed_by_reads' column\
external folder contains data from Huston et al 2021  and Kim et al 2020 (see text)


**scripts** contains all scripts utilized for data analysis and construction of Figures 1 and 2, as well as supplementary figures


