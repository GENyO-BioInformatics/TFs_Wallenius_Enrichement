# TFs Enrichment Analysis Bias Addressing

This github repository contains all the information necessary to perform and exhaustive evaluation of the Singular Enrichment Analysis (SEA) of transcription factors (TFs).

# Repository Structure and Scripts Explained

## Setup  - Install dependencies

```bash
# Python dependencies - requirements.in
python -m venv venv
pip install pip-tools
pip-compile
pip-sync

# R dependencies - R_requirements.in



```


## Data

The data folder contains the base files necessary to carry out this bioinformatics research.

```r
# Run to obtain or update the necessary files
Rscript download_data.R
```

You will download and save:
* The annotation databases:
  * GO_BP.tsv
  * KEGG.tsv
  * Reactome.tsv
  * WikiPathways.tsv
  * annotation_info_table.tsv
* The CollecTRI regulon:
  * collectri_raw.tsv
  * collectri.tsv

## Base Script

Contain the python script (*) to perform the analysis of SEA by the proposed methodologies:
TFs SEA:
* Fisher's Exact Test, Central Hypergeometric Analysis
TFs target genes SEA:
* Fisher's Exact Test, Central Hypergeometric Analysis
* Wallenius' Test, Non-Central Hypergeometric Analysis

#### **TO DO:**
Clean and unify the following scripts in a single one. It should allows a programmatic way to perform the TFs SEA accepting the following parameters:
* TFs (string with comma separated values or a string file path)
* TFs-Targets database
* ?? Organism (not necessary actually, #TODO try using direct ids matching between input and dbs)
* Annotation database
* Methodology
* Output file


Currently we have three scripts:
* analyse_random_lists.py
* random_list_analysis.py
* TFsEnrichment.py
* TFsWalleniusEnrichement.py


## Random Lists Analysis

The folder *random_lists_analysis* contain all the scripts necessary to perform what you intuit :P.

This analysis is divided in the following steps:
1. Generation of universes (*generate_UniversenRandomLists.R*)
2. Creation of random TFs lists (*generate_UniversenRandomLists.R*)
3. Creating bash script to perform the SEA of such random lists.
4. Merging the results and calculating final statistics.

## Use Case

The folder *use_case* contain all the scripts necessary to do the use case SEA

## Figures' data and scripts

The folder *figures_data_n_scripts* exploit the base data and the results from the random lists analysis and the use case. All plots are in .tiff and .png

. figure_collectriTargetsNtfsPerAnnDB.R
. figure_wilcoxontest.R
. figure_probabilities.R
. figure_collectriAnnotationsCoverage.R
. figure_AnotsCurEffortvsProb.R
. randomListsResults_plots.R generate box plots of the annotations significancy frequency:
  .1 RandomTFsHypergeom_terms | RandomTFsHypergeom_ids lanzan
  .2 RandomTFtargetsLists
.


## Other folders

## lib
Contain helper functions of main script

## extra_documentation
Other markdowns used to guide our research by annotating questions and answers.

# TO DOs

1. Be sure to always use annotation id for the unique identifier of the database, and term for the human understandable description of the id. Check scripts.
