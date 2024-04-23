# Phenotypic Rates of Color Pattern Evolution in Ranitomeya Poison Frogs
Code for BAMM analyses from Muell and Brown (2024) for publication in Evolutionary Ecology.

Citation: Muell MR, and Brown JL. 2024. Integrating ecological niche modeling and rates of evolution to model geographic regions of mimetic color pattern selection. Evolutionary Ecology (2024):1-21.

Access the manuscript [here](https://link.springer.com/article/10.1007/s10682-024-10290-8)

Guide to Repo:

"20mil.treefile" is the phylogeny from Muell et al. (2022) used for all analyses.

"rate-weighting.xlsx" contains the tip rates for each phenotype, how they were scaled to relative rates from raw output, and averaging to determine coefficients for GIS layers when building maps. Tab names ending in "GIS" contain final coefficients used in mimicry maps.

"BAMM_output_processing.R" contains BAMMtools code for importing analysis output files, determining rate shift configurations, and plotting output.

BAMM_Analyses folder is subsetted into one folder for each of the four phenotypes. Within each, there are the following:
1) A control file (file extension ".ctl") containing parameters for the BAMM analysis.
2) A trait file (e.g., "fuzzy_stripes.txt") showing the assigned tip states to input into the analysis.
3) Chain_swap file with output chain data.
4) Event_data file containing the inferred rate shift results (import into "BAMM_output_processing.R" to pull results)
5) MCMC_out file with output mcmc information.
6) run info file, output from analysis.

Oh -- and to actually run the program in a Bash terminal, just type:
`bamm -c myControlFile.txt`
