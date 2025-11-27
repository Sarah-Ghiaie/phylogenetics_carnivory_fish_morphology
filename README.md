# Cichlid Phylogeny: Mapping feeding strategies across taxa
This repository contains the code, data, and documentation for a final project in BIOL 4300 - 001.<br>
The project involves gathering publicly available genetic data, cleaning and aligning sequences, and generating a Bayesian phylogeny.<br>
This project’s two main objectives are to:<br> 
(1) Reconstruct a phylogeny of the Family Cichlidae<br> 
(2) Map feeding strategies onto the constructed tree to investigate the evolution of these traits across the phylogeny. 

### Genes of Interest:
(1) COI (cytochrome oxidase subunit 1)<br>
(2) CYTB (cytochrome b)<br>
(3) PTR<br>
(4) ND2 (NADH dehydrogenase subunit 2)<br>
(5) ND5 (NADH dehydrogenase subunit 5)<br>

## Contributors:<br>

Sarah-Ghiaie<br>
Bailey-JJ<br>
bacisthebest

## Tools & Software

| Tool / Software | Packages / Description |
|-----------------|------------------------|
| **IQ-TREE / BEAST2** | Used for phylogenetic inference (ML or Bayesian) |
| **FigTree** | Visualization and annotation of phylogenetic trees |
| **R** | **See packages below |

### R Packages Used

| R Package     | Purpose |
|---------------|---------|
| **rfishbase** | Retrieve biological and ecological data for fish taxa |
| **dplyr**     | Data manipulation and cleaning (part of tidyverse) |
| **rentrez**   | Access and download sequence data from NCBI |
| **XML**       | Parsing XML metadata from online databases |
| **httr**      | Handling HTTP requests when retrieving data |
| **seqinr**    | Reading, writing, and manipulating FASTA sequences |
| **ape**       | Core phylogenetics tree handling and analysis |
| **phytools**  | Ancestral state reconstruction and comparative methods |
| **viridis**   | Colorblind-friendly color palettes for figures |
| **tidyverse** | Data cleaning and visualization framework |

## Usage



## Directory Structure

```markdown
├── data
│   ├── aligned
│   │   ├── coi_68_align.fasta
│   │   ├── cytb_68_align.fasta
│   │   ├── nd2_68_align.fasta
│   │   ├── nd5_68_align.fasta
│   │   └── ptr_68_align.fasta
│   ├── raw
│   │   ├── Cichlid_COI.csv
│   │   ├── Cichlid_cytb.csv
│   │   ├── Cichlid_nd2.csv
│   │   ├── Cichlid_nd5.csv
│   │   ├── Cichlid_Ptr.csv
│   │   ├── coi68.csv
│   │   ├── cytb68.csv
│   │   ├── nd268.csv
│   │   ├── nd568.csv
│   │   ├── ptr68.csv
│   │   └── valid_taxa.csv
│   └── trimmed
│       ├── coi_68_trimmed_filled.fasta
│       ├── cytb_68_trimmed_filled.fasta
│       ├── nd2_68_trimmed_filled.fasta
│       ├── nd5_68_trimmed_filled.fasta
│       ├── ptr_68_trimmed_filled.fasta
│       ├── species_supermatrix.fasta
│       ├── species_supermatrix.fasta.bionj
│       ├── species_supermatrix.fasta.ckp.gz
│       ├── species_supermatrix.fasta.contree
│       ├── species_supermatrix.fasta.iqtree
│       ├── species_supermatrix.fasta.log
│       ├── species_supermatrix.fasta.mldist
│       ├── species_supermatrix.fasta.model.gz
│       ├── species_supermatrix.fasta.splits.nex
│       ├── species_supermatrix.fasta.treefile
│       ├── species_supermatrix.fasta.ufboot
│       └── species_supermatrix.partitions.csv
├── dataset_activity
├── draft_trees
├── final_tree
├── project_setup.sh
├── README.md
├── results
│   └── beast_run
└── scripts
    ├── 5_gene_data_script.R
    ├── filling_nas_alignments.R
    ├── filtering_datasets.R
    ├── potential_functions.R
    └── pulling_sequence_data.R
```
