#### BIOL 4300 - Project Functions ####

create_ds <- function(tax_search, list_ids, taxa_target) {
  
  tax_stat <- matrix(, nrow = length(unique(tax_search$ids)), ncol = 5) # empty table to hold information
  colnames(tax_stat) <- c("taxid", "species_name", "accession_num", "seq_name", "seq_len")
    
  tax_stat[,1] <- tax_search$ids
  
  tax_stat_df <- data.frame(
    taxid = rep(NA_character_, length(unique(tax_search$ids))),
    species_name = rep(NA_character_, length(unique(tax_search$ids))),
    accession_num = I(vector("list", length(unique(tax_search$ids)))),
    seq_name = I(vector("list", length(unique(tax_search$ids)))),
    seq_len = I(vector("list", length(unique(tax_search$ids)))),
    stringsAsFactors = FALSE
  )
  
  tax_stat_df[,1] <- tax_search$ids
  
  # For-loop that compares species names from NCBI search results with Fishbase.
  # Names from NCBI that are listed as taxonomically valid are kept, the others are not listed in the final data frame.
  
  # NOTE: Loop runs backwards to avoid the error of deleting a row, and then shrinking the number of values to compare.
  
  # Takes roughly 36 minutes, maybe there is a way to cut down time?
  start_time <- Sys.time() # Using to track how long this takes to run
  
  for (i in length(unique(tax_search$ids:1))) {
    print(paste(i, "a", sep = ""))
    
    record <- entrez_summary(db = "taxonomy", id = tax_search$ids[i])
    
    if (record$scientificname %in% cichlids$Species) { # Compares the scientific name associated with a NCBI TAXID to Fishbase.
      tax_stat_df[i,2] <- record$scientificname
      #tax_stat_df[i,3] <- record$
        print("Species added.")
    }
    else { # If name is not found, delete that entire row from the data frame.
      tax_stat_df <- tax_stat_df[-i,]
      print(paste("Deleted row:", i, sep = " "))
    }
  }
  
  end_time <- Sys.time()
  
  print(paste("Total Number of Taxonomically Valid Matches after", as.numeric(round((end_time - start_time), digits = 2), units = "mins"), "minutes:", nrow(tax_stat_df)))
  # Output should look something like: "Total Number of Taxonomically Valid Matches after ~~ minutes: ~~~~"
  
  return(tax_stat_df)
}
