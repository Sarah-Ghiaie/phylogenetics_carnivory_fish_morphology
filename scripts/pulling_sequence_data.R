#### Pulling Cichlidae species ####
library(rfishbase)
library(dplyr)
library(rentrez)
library(XML)
library(httr)
library(seqinr)
source("./potential_functions.R")
source("./please.R")

#### Defining Variables: ####
set_entrez_key("d9b97feea308d7dbc357193e79a85721f808")


cich_search <- entrez_search(db = "taxonomy", term = "txid8113[SBTR]", retmax = 5000)
search_ids <- cich_search$ids
cichlids <- load_taxa() %>% dplyr::filter(Family == "Cichlidae") # list of all valid species from fishbase

valid_tax <- create_ds_optimized(cich_search, search_ids, cichlids)

# Get taxid for taxa of interest, not all species need to match every genes (we can make a super matrix). I think we should try to have two genes for every species we use, but that might not be possible.


list_cols <- sapply(valid_tax, is.list)
valid_tax[list_cols] <- lapply(valid_tax[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(valid_tax, "./data/raw/valid_taxa.csv", row.names = FALSE)

#### Pulling data and making dataset ####
# Potential genes of interest (only need 3 for the upcoming activity):

#1) COI
#2) Ptr
#3) cytb (cytochrome b gene)


# Finding accession numbers for a gene across at least 50 taxa
seq_count <- 0

for (i in 1:nrow(valid_tax)) {
  
  if (seq_count >= 50) {
    print("50 taxa hits for COI have been found.")
    break
  }
  
  name <- valid_tax$species_name[i]
  
  tryCatch({
    
    nuc_search <- entrez_search(db = "nuccore", 
                                term = paste0('"', name, '"', '[ORGN] AND COI[GENE]'), 
                                retmax = 5)
    
    if(length(nuc_search$ids)<1) {
      valid_tax[i,3:5] <- NA
    }
    
    else {
      sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
      valid_tax[i,3] <- sumout$caption
      valid_tax[i,4] <- sumout$title
      valid_tax[i,5] <- sumout$slen
      
      print(paste("Sequence found for", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      coi_sequence <- c(coi_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
    })
}
print("DONE")

writeLines(as.character(coi_sequence), "./data/raw/metadata/coi_out.fasta")


# This table will contain all species that had hits for COI
# Creating a copy for saving
cich_coi_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_coi_save, is.list)
cich_coi_save[list_cols] <- lapply(cich_coi_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_coi_save, "./data/raw/Cichlid_COI.csv", row.names = FALSE)


cich_coi_save %>% 
  dplyr::filter(accession_num != 'NA') %>% 
  summarise(count = n())



#### Now let's look at the Ptr gene ####

seq_count <- 0

for (i in 1:nrow(valid_tax)) {
  
  if (seq_count >= 50) {
    print("50 taxa hits for Ptr have been found.")
    break
  }
  
  name <- valid_tax$species_name[i]
  
  tryCatch({
    
    nuc_search <- entrez_search(db = "nuccore", 
                                term = paste0('"', name, '"', '[ORGN] AND Ptr[GENE]'), 
                                retmax = 5)
    
    if(length(nuc_search$ids)<1) {
      valid_tax[i,3:5] <- NA
    }
    
    else {
      sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
      valid_tax[i,3] <- sumout$caption
      valid_tax[i,4] <- sumout$title
      valid_tax[i,5] <- sumout$slen
      
      print(paste("Sequence found for", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      coi_sequence <- c(coi_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  })
}
print("DONE")

writeLines(as.character(coi_sequence), "./data/raw/metadata/ptr_out.fasta")


# This table will contain all species that had hits for COI
# Creating a copy for saving
cich_ptr_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_ptr_save, is.list)
cich_ptr_save[list_cols] <- lapply(cich_ptr_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_ptr_save, "./data/raw/Cichlid_Ptr.csv", row.names = FALSE)


cich_ptr_save %>% 
  dplyr::filter(accession_num != 'NA') %>% 
  summarise(count = n())



#### Now let's look at cytb ####

seq_count <- 0

for (i in 1:nrow(valid_tax)) {
  
  if (seq_count >= 50) {
    print("50 taxa hits for cytb have been found.")
    break
  }
  
  name <- valid_tax$species_name[i]
  
  tryCatch({
    
    nuc_search <- entrez_search(db = "nuccore", 
                                term = paste0('"', name, '"', '[ORGN] AND cytb[GENE]'), 
                                retmax = 5)
    
    if(length(nuc_search$ids)<1) {
      valid_tax[i,3:5] <- NA
    }
    
    else {
      sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
      valid_tax[i,3] <- sumout$caption
      valid_tax[i,4] <- sumout$title
      valid_tax[i,5] <- sumout$slen
      
      print(paste("Sequence found for", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      coi_sequence <- c(coi_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  })
}
print("DONE")

writeLines(as.character(coi_sequence), "./data/raw/metadata/cytb_out.fasta")


# This table will contain all species that had hits for COI
# Creating a copy for saving
cich_cytb_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_cytb_save, is.list)
cich_cytb_save[list_cols] <- lapply(cich_cytb_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_cytb_save, "./data/raw/Cichlid_cytb.csv", row.names = FALSE)


cich_cytb_save %>% 
  dplyr::filter(accession_num != 'NA') %>% 
  summarise(count = n())
