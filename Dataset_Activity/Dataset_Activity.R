# Install and load required packages
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

library(rentrez)
library(stringr)

# Define queries with synonyms for each gene
gene_queries <- list(
  "cytochrome b" = '(cytochrome b[All Fields] OR cytb[All Fields]) AND txid8113[Organism:exp]',
  "ND2"          = '(ND2[All Fields] OR "NADH dehydrogenase subunit 2"[All Fields]) AND txid8113[Organism:exp]',
  "COI"          = '(COI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit I"[All Fields]) AND txid8113[Organism:exp]',
  "RAG1"         = 'RAG1[All Fields] AND txid8113[Organism:exp]'
)

# Function to fetch up to 50 unique taxa per gene
fetch_unique_taxa <- function(gene_symbol, query, out_file, db = "nuccore", retmax = 500) {
  message("\nSearching for: ", gene_symbol)
  
  # Search NCBI
  search_res <- entrez_search(db = db, term = query, use_history = TRUE, retmax = retmax)
  if (search_res$count == 0) {
    message("No results for ", gene_symbol)
    return(NULL)
  }
  
  # Fetch summaries (accession + taxon)
  summaries <- entrez_summary(db = db, web_history = search_res$web_history, retmax = retmax)
  accession <- sapply(summaries, function(x) x$caption)
  taxon     <- sapply(summaries, function(x) x$organism)
  
  # Deduplicate by taxon
  unique_idx <- !duplicated(taxon)
  accession <- accession[unique_idx]
  taxon <- taxon[unique_idx]
  
  # Limit to 50 unique taxa
  if (length(accession) < 50) {
    warning("Only ", length(accession), " unique taxa found for ", gene_symbol)
  }
  accession <- head(accession, 50)
  taxon <- head(taxon, 50)
  
  # Fetch FASTA for these accession numbers
  fasta <- entrez_fetch(db = db, id = accession, rettype = "fasta", retmode = "text")
  
  # Save FASTA
  write(fasta, file = out_file)
  message("Saved ", length(accession), " unique taxa sequences to ", out_file)
  
  # Print accession and taxon
  cat("\n--- Results for", gene_symbol, "---\n")
  for (i in seq_along(accession)) {
    cat(accession[i], "-", taxon[i], "\n")
  }
  
  return(list(accession = accession, taxon = taxon, fasta = fasta))
}

# Run for all genes
results <- list()
for (gene in names(gene_queries)) {
  out_file <- paste0(gsub(" ", "_", gene), "_cichlidae.fasta")
  results[[gene]] <- fetch_unique_taxa(gene, gene_queries[[gene]], out_file, retmax = 500)
}

message("\nAll sequences downloaded (up to 50 unique taxa per gene).")
