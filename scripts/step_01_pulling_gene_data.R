library(rfishbase)
library(dplyr)
library(rentrez)
library(XML)
library(httr)
library(seqinr)
source("./scripts/project_functions.R")

# Defining Variables: ####
set_entrez_key("d9b97feea308d7dbc357193e79a85721f808") # Lets us pull more stuff from NCBI without getting booted.


cich_search <- entrez_search(db = "taxonomy", term = "txid8113[SBTR]", retmax = 5000) # Searching for all taxids under the Cichlid family.
search_ids <- cich_search$ids # List of taxids
cichlids <- load_taxa() %>% dplyr::filter(Family == "Cichlidae") # list of all valid species from fishbase

# Comparing species names found in both NCBI and Fishbase, if there is overlap, those species are kept.
valid_tax <- create_ds_optimized(cich_search, search_ids, cichlids) 


list_cols <- sapply(valid_taxa, is.list)
valid_taxa[list_cols] <- lapply(valid_taxa[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))}) # Making columns easily useable down the line


# Can uncomment this write command, if needed. But valid.taxa.csv is already made.
#write.csv(valid_tax, "./data/raw/valid_taxa.csv", row.names = FALSE)


#### Pulling data and making dataset ####
# Genes of interest
#1) COI
#2) Ptr
#3) cytb (cytochrome b gene)
#4) ND2 (NADH dehydrogenase subunit 2)
#5) ND5 (NADH dehydrogenase subunit 5)


# Creating .csv for all species with a specific gene hit, and then making a .fasta file for matches as well.
valid_tax <- read.csv("./data/raw/valid_taxa.csv")

# Pulling COI data ####
seq_count <- 0
coi_sequence <- character()

for (i in 1:nrow(valid_tax)) {
  
  # if (seq_count >= 50) {
  #   print("50 taxa hits for COI have been found.")
  #   break
  # }  # Can uncomment this to test script runs correctly before searching for a large number of taxa
  
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
      
      print(paste("Sequence found for:", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      coi_sequence <- c(coi_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  }, error = function(e) {
    print(paste("Error processing:", name, "-", e$message))
    valid_tax[i, 3:5] <<- NA
    
  })
}
print("DONE")

writeLines(as.character(coi_sequence), "./data/raw/coi_out.fasta")

# This table will contain all species that had hits for COI
# Creating a copy for saving
cich_coi_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_coi_save, is.list)
cich_coi_save[list_cols] <- lapply(cich_coi_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_coi_save, "./data/raw/Cichlid_COI.csv", row.names = FALSE)








# Pulling Ptr data ####

valid_tax <- read.csv("./data/raw/valid_taxa.csv") # Resetting dataset

seq_count <- 0
ptr_sequence <- character()

for (i in 1:nrow(valid_tax)) {
  
  # if (seq_count >= 50) {
  #   print("50 taxa hits for Ptr have been found.")
  #   break
  # }
  
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
      
      print(paste("Sequence found for:", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      ptr_sequence <- c(ptr_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  }, error = function(e) {
    print(paste("Error processing:", name, "-", e$message))
    valid_tax[i, 3:5] <<- NA
    
  })
}
print("DONE")

writeLines(as.character(ptr_sequence), "./data/raw/ptr_out.fasta")

# This table will contain all species that had hits for Ptr
# Creating a copy for saving
cich_ptr_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_ptr_save, is.list)
cich_ptr_save[list_cols] <- lapply(cich_ptr_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_ptr_save, "./data/raw/Cichlid_Ptr.csv", row.names = FALSE)







# Pulling cytb data ####

valid_tax <- read.csv("./data/raw/valid_taxa.csv") # Resetting dataset

seq_count <- 0
cytb_sequence <- character()

for (i in 1:nrow(valid_tax)) {
  
  # if (seq_count >= 50) {
  #   print("50 taxa hits for cytb have been found.")
  #   break
  # }
  
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
      
      print(paste("Sequence found for:", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      cytb_sequence <- c(cytb_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  }, error = function(e) {
    print(paste("Error processing:", name, "-", e$message))
    valid_tax[i, 3:5] <<- NA
    
  })
}
print("DONE")

writeLines(as.character(cytb_sequence), "./data/raw/cytb_out.fasta")

# This table will contain all species that had hits for Ptr
# Creating a copy for saving
cich_cytb_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_cytb_save, is.list)
cich_cytb_save[list_cols] <- lapply(cich_cytb_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_cytb_save, "./data/raw/Cichlid_cytb.csv", row.names = FALSE)






# Pulling ND2 data ####

valid_tax <- read.csv("./data/raw/valid_taxa.csv") # Resetting dataset

seq_count <- 0
nd2_sequence <- character()

for (i in 1:nrow(valid_tax)) {
  
  # if (seq_count >= 50) {
  #   print("50 taxa hits for cytb have been found.")
  #   break
  # }
  
  name <- valid_tax$species_name[i]
  
  tryCatch({
    
    nuc_search <- entrez_search(db = "nuccore", 
                                term = paste0('"', name, '"', '[ORGN] AND ND2[GENE]'), 
                                retmax = 5)
    
    if(length(nuc_search$ids)<1) {
      valid_tax[i,3:5] <- NA
    }
    
    else {
      sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
      valid_tax[i,3] <- sumout$caption
      valid_tax[i,4] <- sumout$title
      valid_tax[i,5] <- sumout$slen
      
      print(paste("Sequence found for:", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      nd2_sequence <- c(nd2_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  }, error = function(e) {
    print(paste("Error processing:", name, "-", e$message))
    valid_tax[i, 3:5] <<- NA
    
  })
}
print("DONE")

writeLines(as.character(nd2_sequence), "./data/raw/nd2_out.fasta")

# This table will contain all species that had hits for ND2
# Creating a copy for saving
cich_nd2_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_nd2_save, is.list)
cich_nd2_save[list_cols] <- lapply(cich_nd2_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_nd2_save, "./data/raw/Cichlid_nd2.csv", row.names = FALSE)








# Pulling ND5 data ####

valid_tax <- read.csv("./data/raw/valid_taxa.csv") # Resetting dataset

seq_count <- 0
nd5_sequence <- character()

for (i in 1:nrow(valid_tax)) {
  
  # if (seq_count >= 50) {
  #   print("50 taxa hits for cytb have been found.")
  #   break
  # }
  
  name <- valid_tax$species_name[i]
  
  tryCatch({
    
    nuc_search <- entrez_search(db = "nuccore", 
                                term = paste0('"', name, '"', '[ORGN] AND ND5[GENE]'), 
                                retmax = 5)
    
    if(length(nuc_search$ids)<1) {
      valid_tax[i,3:5] <- NA
    }
    
    else {
      sumout <- entrez_summary(db = "nuccore", id = nuc_search$ids[1])
      valid_tax[i,3] <- sumout$caption
      valid_tax[i,4] <- sumout$title
      valid_tax[i,5] <- sumout$slen
      
      print(paste("Sequence found for:", name, "Getting fasta"))
      
      seq.dat <- entrez_fetch(db = "nuccore", id = nuc_search$ids[1], rettype = "fasta")
      nd5_sequence <- c(nd5_sequence, seq.dat)
      
      print("Fasta added.")
      seq_count <- seq_count + 1
    }
    
  }, error = function(e) {
    print(paste("Error processing:", name, "-", e$message))
    valid_tax[i, 3:5] <<- NA
    
  })
}
print("DONE")

writeLines(as.character(nd5_sequence), "./data/raw/nd5_out.fasta")

# This table will contain all species that had hits for ND5
# Creating a copy for saving
cich_nd5_save <- valid_tax

# Convert list columns to character strings
list_cols <- sapply(cich_nd5_save, is.list)
cich_nd5_save[list_cols] <- lapply(cich_nd5_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_nd5_save, "./data/raw/Cichlid_nd5.csv", row.names = FALSE)


