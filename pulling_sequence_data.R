#### Pulling Cichlidae species ####
library(rfishbase)
library(dplyr)
library(rentrez)
library(XML)
library(httr)
library(seqinr)

cich_search <- entrez_search(db = "taxonomy", term = "txid8113[SBTR]", retmax = 5000)
search_ids <- cich_search$ids

cichlids <- load_taxa() %>% dplyr::filter(Family == "Cichlidae") # list of all valid species from fishbase
cich_stat <- matrix(, nrow = length(unique(cich_search$ids)), ncol = 5) # empty table to hold information
colnames(cich_stat) <- c("taxid", "species_name", "accession_num", "seq_name", "seq_len")

cich_stat[,1] <- cich_search$ids



# For-loop that compares species names from NCBI search results with Fishbase.
# Names from NCBI that are listed as taxonomically valid are kept, the others are not listed in the final data frame.

# NOTE: Loop runs backwards to avoid the error of deleting a row, and then shrinking the number of values to compare.

# Takes roughly 36 minutes, maybe there is a way to cut down time?
start_time <- Sys.time() # Using to track how long this takes to run

for (i in length(unique(cich_search$ids)):1) {
  print(paste(i, "a", sep = ""))
  
  record <- entrez_summary(db = "taxonomy", id = cich_search$ids[i])
  
  if (record$scientificname %in% cichlids$Species) { # Compares the scientific name associated with a NCBI TAXID to Fishbase.
    cich_stat[i,2] <- record$scientificname
    print("Species added.")
    }
  else { # If name is not found, delete that entire row from the data frame.
    cich_stat <- cich_stat[-i,]
    print(paste("Deleted row:", i, sep = " "))
    }
}

end_time <- Sys.time()

print(paste("Total Number of Taxonomically Valid Matches after", as.numeric(round((end_time - start_time), digits = 2), units = "mins"), "minutes:", nrow(cich_stat)))
# "Total Number of Taxonomically Valid Matches after 35.69 minutes: 1200"






#seqout <- entrez_search(db = "nuccore", term = paste("txid",taxout$ids[i],"[ORGN] AND MITOCHONDRION[ALL]", sep = ""))


