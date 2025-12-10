# Cichlid Phylogeny: Mapping feeding strategies across taxa
This repository contains the code, data, and documentation for a final project in BIOL 4300 - 001.<br>
The project involves gathering publicly available genetic data, cleaning and aligning sequences, and generating a Bayesian phylogeny.<br>
This project’s two main objectives are to:<br> 
(1) Reconstruct a phylogeny of the Family Cichlidae<br> 
(2) Map feeding strategies onto the constructed tree to investigate the evolution of these traits across the phylogeny. 

### Genes of Interest:
(1) COI (cytochrome oxidase subunit 1)<br>
(2) CYTB (cytochrome b)<br>
(3) PTR (Patch-Related)<br>
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
| **R/RStudio** | **See packages below |

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
| **evobiR**    | Collection of tools for evolutionary biology analyses, including sequence and phylogenetic workflows |
| **msaR**      | Interactive visualization of multiple sequence alignments |
| **phangorn**  | Phylogenetic inference and model testing |

## Usage

Files contained in the **scripts/** directory can be run in order. R and RStudio are required to run these scripts.

## Directory Structure

| Sub-Directory | Content Description |
|---------------|---------|
| **data/** | Contains sub-directories: alignments/; raw/; supermatrix/; trimmed/ and all associated files. |
| **scripts/** | Contains all R scripts needed to get project results |
| **Dataset_Activity** | Easy access to files needed for this assignment submission |
| **Draft_Trees** | Easy access to files needed for this assignment submission |
| **Final_Tree** | Easy access to files needed for final assignment submission |
