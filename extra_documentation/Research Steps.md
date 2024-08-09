## Motivation
The enrichment analysis of TFs via their target genes is biased and functional annotations appear as significant ubiquitously when testing via the central hypergeometric distribution. Using the Wallenius non-central distribution the bias should be pondered. 

## Methodology
#### [Monte Carlo Simulation](https://en.wikipedia.org/wiki/Monte_Carlo_method)
Test random lists of TFs with 3 approaches:
 1. Central hypergeometric distribution + TFs annotations 
 We will recover only functional annotations that describe the very basic function of the TFs
2. Central hypergeometric distribution + TFs-targets annotations
 We will recover as top significant biased functional annotations that are specific of each database, namely, regulation of transcription and cancer
 3.  Wallenius non-central distribution + TFs-targets annotations
The biassed annotations will not appear so up in the rank  and some less represented will be recovered.

#### Databases
1. **Collectri** is the TFs-TargetGene database selected because is a comprehensive and curated effort that includes several other dbs.
2. 4 different functional annotations downloaded from GeneCodis4.
	2.1. Gene Ontology Biological Process
	2.2. Reactome 
	2.3. KEGG
	2.4. Wikipathways



### Steps
###### 1 - Study the context of Collectri in the annotation databases and generate universes
We generate different Venn Diagrams of the TFs, Targets and Interactions of each database but limiting them to the set of TFs which each one of them target at least one gene in the functional annotation database. See script "generate_dbsCommonUniverse.R"




