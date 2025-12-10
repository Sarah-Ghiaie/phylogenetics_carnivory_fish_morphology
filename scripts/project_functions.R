#### BIOL 4300 - Project Functions ####

## Valid taxa dataframe function ####
create_ds_optimized <- function(tax_search, list_ids, taxa_target) {
  
  all_ids <- unique(tax_search$ids) # Just double checking there aren't any repeat names
  n_ids <- length(all_ids)
  
  print(paste("Processing", n_ids, "taxonomy IDs..."))
  
  # Batch processing to reduce API calls
  batch_size <- 200
  all_records <- list()
  
  start_time <- Sys.time()
  
  for (i in seq(1, n_ids, by = batch_size)) {
    end_idx <- min(i + batch_size - 1, n_ids)
    batch_ids <- all_ids[i:end_idx]
    
    print(paste("Processing batch", ceiling(i/batch_size), "of", ceiling(n_ids/batch_size)))
    
    # Get summaries for entire batch at once
    batch_records <- entrez_summary(db = "taxonomy", id = batch_ids)
    
    # Handle both single records and lists of records
    if (length(batch_ids) == 1) {
      batch_records <- list(batch_records)
      names(batch_records) <- batch_ids
    }
    
    all_records <- c(all_records, batch_records)
  }
  
  # Extract species names and filter in vectorized operations
  print("Extracting species names and filtering...")
  
  # Vectorized extraction of scientific names
  species_names <- sapply(all_records, function(x) x$scientificname)
  taxids <- names(all_records)
  
  # Vectorized filtering using %in% (much faster than loop)
  valid_mask <- species_names %in% cichlids_set
  
  # Keep only valid species
  valid_taxids <- taxids[valid_mask]
  valid_species <- species_names[valid_mask]
  
  print(paste("Found", sum(valid_mask), "valid species out of", length(species_names), "total"))
  
  # Create final dataframe
  n_valid <- length(valid_taxids)
  tax_df <- data.frame(
    taxid = valid_taxids,
    species_name = valid_species,
    stringsAsFactors = FALSE
  )
  
  # Pre-allocate list columns
  tax_df$accession_num <- vector("list", n_valid)
  tax_df$seq_name <- vector("list", n_valid)
  tax_df$seq_len <- vector("list", n_valid)
  
  end_time <- Sys.time()
  print(paste("Total processing time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
  
  return(tax_df)
}


# DONT USE THE BELOW FUNCTIONS, JUST RUN THE NORMAL SCRIPTS ON OTHER MODULES
## Search for gene hits, generate .cv and .fasta files ####

gene_search <- function(valid_tax, gene, outfile_prefix) {
  
  valid_tax2 <- valid_tax
  gene_sequence <- character()
  seq_count <- 0 # This just counts the number of gene hits
  
  for (i in 1:nrow(valid_tax2)) {
    
    # Uncomment this section, to run tests on function
    if (seq_count >= 5) {
     print("5 taxa hits for COI have been found.")
     break
    }
    
    name <- valid_tax2$species_name[i]
    gene_search <- gene
  
    
    tryCatch({
      
      nuc_search <- entrez_search(db = "nuccore", 
                                  term = paste0('"', name, '"', '[ORGN] AND ', gene_search, '[GENE]'), 
                                  retmax = 5)
      print(paste0('"', name, '"', '[ORGN] AND ', gene_search, '[GENE]'))
      
      
      if(is.null(nuc_search$ids) || length(nuc_search$ids)<1) { # If no hits appear, add NA to the row
        print(paste("No sequence found for: ", name))
        valid_tax2[i,3:5] <- NA
      }
      
      else { # If results appear, add search accession number, gene title, and length.
        
        sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
        valid_tax2[i,3] <- sumout$caption
        valid_tax2[i,4] <- sumout$title
        valid_tax2[i,5] <- sumout$slen
        
        print(paste("Sequence found for", name,))

        #seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
        #coi_sequence <- c(coi_sequence, seq.dat)
        
        #print("Fasta added.")
        seq_count <- seq_count + 1        
      }
      
    }, error = function(e) {
      message("Error processing ", name, ": ", e$message)
      #valid_tax2[i,3:5] <- NA
      
    })
  }
  print(paste("DONE. Number of hits: ", seq_count))

  
  # This table will contain all species that had hits for that gene
  # Convert list columns to character strings
  list_cols <- sapply(valid_tax2, is.list)
  valid_tax2[list_cols] <- lapply(valid_tax2[list_cols], function(x) {
    sapply(x, function(y) paste(y, collapse = "; "))
  })
  
  write.csv(valid_tax2, paste0("./data/raw/metadata/Cichlid_", gene_search, ".csv"), row.names = FALSE)
  
  print(".csv results can be found in the ./data/raw/metadata directory.")
}


## Generate .fasta files for each gene using the master .csv ####

get_fasta <- function(master_csv, gene, outfile_prefix) {
  master_csv2 <- master_csv
  gene_sequence <- character()
  seq_count <- 0 # This just counts the number of gene hits
  
  for (i in 1:nrow(master_csv2)) {
    
    # Uncomment this section, to run tests on function
    # if (seq_count >= 20) {
    #   print("20 taxa hits for COI have been found.")
    #   break
    # }
    
    name <- master_csv2$species_name[i]
    query <- paste0("[ORGN] AND ", gene, "[GENE]")
    
    tryCatch({
      
      nuc_search <- entrez_search(db = "nuccore", 
                                  term = paste0('"', name, '"', query), 
                                  retmax = 5)
      
      if(length(nuc_search$ids) > 0) {
        
        print(paste("Sequence found for", name, "Getting fasta"))
        
        seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
        clean_name <- gsub(" ", "_", name)
        seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
        gene_sequence <- c(gene_sequence, seq.dat2)
        
        print("Fasta added.")
        seq_count <- seq_count + 1
      }
      
      else {
        message("No hits for ", name)
      }
      
    }, error = function(e) {
      message("Error processing ", name, ": ", e$message)
      master_csv2[i,3:5]
      
    })
  }
  print(paste("DONE. Number of hits: ", seq_count))
  
  writeLines(as.character(gene_sequence), paste0("./data/raw/seq_raw", outfile_prefix, ".fasta"))
  
}