| Script | Description |
| --- | --- |
| step_01_pulling_gene_data.R | (1) Gathers valid taxa, and generates a .csv file containing valid taxa.<br> (2) Pulls sequence data from NCBI for 5 selected genes, generates a new .csv file containing species hits for each gene. |
| step_02_filtering_datasets.R | (1) Filters each tables to contain only species that had hits. <br> (2) User can then further filter tables to contain only species that have hits in at least **[*user defined #*]** tables. |
| step_03_generate_supermatrix.R | Cleans and trims multiple gene alignments, fills in missing taxa, and outputs FASTA files. Builds a concatenated supermatric for phylogenetic analysis. |
| step_04_filling_nas_alignments.R | For taxa that are not found in all .fasta files, script will fill *_trimmed.fasta files with 'N's. |
| project_functions.R | This script will be loaded when running step_01, in order to pull data properly. |

