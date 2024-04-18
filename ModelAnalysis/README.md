## Model analysis

Requires: Claudia_medium_cons_250723.mat


Model statistics and Fig 4C generated with: driverAnalysis_0723_PAPER.m


Flux Balance Analysis (FBA) per model done by running in MATLAB:

[solf,resEX] = FBA_mediumConc_varO2_v2(modelNr,21.253)

with modelNr from 1 .. 4 (for experimental conditions: 'WT-D30';'WT-D60'; 'PD-D30';'PD-D60'


Fig4DE generated from FBA results shown in SupplementaryTable_XX2 (plotting)


Fig5A-E generated from FBA results saved in FBAmodel[1-4].mat with:

compareFBA_v3_PAPER.m

Respective cell-type specific supplementary figures generated with:

compareFBA_v3_PAPER_perCelltype.m
