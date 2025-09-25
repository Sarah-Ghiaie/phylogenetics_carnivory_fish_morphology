#### Pulling Cichlidae species ####
library(rfishbase)
library(dplyr)
library(rentrez)
library(XML)
library(httr)
library(seqinr)
source("./potential_functions.R")

#### Defining Variables: ####
cich_search <- entrez_search(db = "taxonomy", term = "txid8113[SBTR]", retmax = 5000)
search_ids <- cich_search$ids
cichlids <- load_taxa() %>% dplyr::filter(Family == "Cichlidae") # list of all valid species from fishbase


#### Pulling data and making dataset ####
# Potential genes of interest (only need 3 for the upcoming activity):

#1) COI
#2) Ptr
#3) cytb (cytochrome b gene)


# Finding accession numbers for a gene across at least 50 taxa

for (i in 1:nrow(cich_stat_safe)) {
  
  name <- cich_stat_safe$species_name[i]
  nuc_search <- entrez_search(db = "nuccore", 
                              term = paste0(name, "[ORGN] AND COI[GENE]"), 
                              retmax = 5)
  
  print(nuc_search$ids)
  cich_stat_safe$accession_num[i] <- nuc_search$ids[1]
}

# This table will contain all species that had hits for COI
cich_coi_df <- cich_stat_safe

# Creating a copy for saving
cich_coi_save <- cich_coi_df

# Convert list columns to character strings
list_cols <- sapply(cich_coi_save, is.list)
cich_coi_save[list_cols] <- lapply(cich_coi_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_coi_save, "./Cichlid_COI.csv", row.names = FALSE)


cich_coi_df %>% 
  filter(sapply(accession_num, length) != 0) %>% 
  count()

      # 400/1200 species had results for the COI gene



# Now using that 400 species, let's do it again but for another gene
# Using the saved table will help use filter even further down to less taxa (We can include more taxa later, I just think working with not 1200 species would be a lot)

# NOTE: I know I wrote above that I was using the saved table, but I forgot to actually create a new table that only includes rows that actually had accession numbers, so we will have to do that at some point. But for now,
        # I really don't care.

        # A solution would be to just compare each table, and remove species that don't have accession numbers in both tables, and then compare that table to the third table. Hopefully that leaves us with more than 50 species
        # so we can do more filter as we add another 2 genes.



#### Now let's look at the Ptr gene? ####

cich_ptr_df <- cich_coi_df

for (i in 1:nrow(cich_ptr_df)) {
  
  name <- cich_ptr_df$species_name[i]
  nuc_search <- entrez_search(db = "nuccore", 
                              term = paste0(name, "[ORGN] AND Ptr gene"), 
                              retmax = 5)
  
  print(nuc_search$ids)
  cich_ptr_df$accession_num[i] <- nuc_search$ids[1]
}

cich_ptr_save <- cich_ptr_df

# Convert list columns to character strings
list_cols <- sapply(cich_ptr_save, is.list)
cich_ptr_save[list_cols] <- lapply(cich_ptr_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_ptr_save, "./Cichlid_Ptr.csv", row.names = FALSE)


cich_ptr_df %>% 
  filter(sapply(accession_num, length) != 0) %>% 
  count()

# 166 species had results for the Ptr gene





#### Now let's look at cytb??

cich_cytb_df <- cich_ptr_df

for (i in 1:nrow(cich_cytb_df)) {
  
  name <- cich_cytb_df$species_name[i]
  nuc_search <- entrez_search(db = "nuccore", 
                              term = paste0(name, "[ORGN] AND cytb gene"), 
                              retmax = 5)
  
  print(nuc_search$ids)
  cich_cytb_df$accession_num[i] <- nuc_search$ids[1]
}

cich_cytb_save <- cich_cytb_df

# Convert list columns to character strings
list_cols <- sapply(cich_cytb_save, is.list)
cich_cytb_save[list_cols] <- lapply(cich_cytb_save[list_cols], function(x) {
  sapply(x, function(y) paste(y, collapse = "; "))
})

write.csv(cich_cytb_save, "./Cichlid_cytb.csv", row.names = FALSE)


cich_cytb_df %>% 
  filter(sapply(accession_num, length) != 0) %>% 
  count()

# 543/1200 species had results for the cytb genes.



# Now we need to compare all of these tables (Or if you guys have time, just rerun the Ptr and cytb script, but just use a new table that doesn't inlcude NULL accession numnber values)
# This will give us a workable table to get our fastas




#### Now let's get the fasta files for all these seuquences!

# We need one file per gene, so using the most up to date table ("./Cichlid_cytb.csv"), we can loop through and make a fasta file for each of the three genes.

