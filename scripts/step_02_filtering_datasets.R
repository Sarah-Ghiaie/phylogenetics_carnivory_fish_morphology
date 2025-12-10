library(tidyverse)
library(rentrez)

# Cleaning the data tables, and creating new versions that only include species that 
# overlap between species found in other tables.

csv_coi <- read.csv("./data/raw/Cichlid_COI.csv")
csv_ptr <- read.csv("./data/raw/Cichlid_Ptr.csv")
csv_cytb <- read.csv("./data/raw/Cichlid_cytb.csv")
csv_nd2 <- read.csv("./data/raw/Cichlid_nd2.csv")
csv_nd5 <- read.csv("./data/raw/Cichlid_nd5.csv")


# Filter out rows that don't have results ####
filtered_coi <- csv_coi %>% 
  filter(!is.na(accession_num))

filtered_ptr <- csv_ptr %>% 
  filter(!is.na(accession_num))

filtered_cytb <- csv_cytb %>% 
  filter(!is.na(accession_num))

filtered_nd2 <- csv_nd2 %>% 
  filter(!is.na(accession_num))

filtered_nd5 <- csv_nd5 %>% 
  filter(!is.na(accession_num))


# Finding common species between 2-5 datasets ####
common_species <- intersect(intersect(filtered_coi$taxid, filtered_ptr$taxid), filtered_cytb$taxid)

common_species_names <- intersect(intersect(filtered_coi$species_name, filtered_ptr$species_name), filtered_cytb$species_name)
write.csv(common_species_names, "~/Desktop/68_diet_data.csv")

length(common_species) # 68 taxa across the first 3 genes

## Testing 2 gene intersections: ####
# 94 across only 2 genes (coi + ptr):
length(intersect(filtered_coi$taxid, filtered_ptr$taxid))

# 275 across only 2 genes (coi + cytb):
length(intersect(filtered_coi$taxid, filtered_cytb$taxid))
 
## Testing 3 gene intersections: ####
# 68 across 3 genes:
length(intersect(intersect(filtered_coi$taxid, filtered_ptr$taxid), filtered_cytb$taxid))
## Testing 4 gene intersections: ####
# 58 across 4 genes:
length(intersect(intersect(intersect(filtered_coi$taxid, filtered_ptr$taxid), filtered_cytb$taxid), filtered_nd2$taxid))

## Testing 5 gene intersections: ####
# 28 across 5 genes:
length(intersect(intersect(intersect(intersect(filtered_coi$taxid, filtered_ptr$taxid), filtered_cytb$taxid), filtered_nd2$taxid), filtered_nd5$taxid))



# Getting .csv for only the 68 taxa ####

coi68 <- csv_coi %>% 
  filter(taxid %in% common_species)

write.csv(coi68, "./data/raw/coi68.csv", row.names = FALSE)

ptr68 <- csv_ptr %>% 
  filter(taxid %in% common_species)

write.csv(ptr68, "./data/raw/ptr68.csv", row.names = FALSE)

cytb68 <- csv_cytb %>% 
  filter(taxid %in% common_species)

write.csv(cytb68, "./data/raw/cytb68.csv", row.names = FALSE)

nd2_filter <- csv_nd2 %>% 
  filter(taxid %in% common_species)

write.csv(nd2_filter, "./data/raw/nd268.csv", row.names = FALSE)

nd5_filter <- csv_nd5 %>% 
  filter(taxid %in% common_species)

write.csv(nd5_filter, "./data/raw/nd568.csv", row.names = FALSE)



# Getting fastas for the 68 ####

### COI ####
set_entrez_key("d9b97feea308d7dbc357193e79a85721f808")
coi68_csvtofasta <- read.csv("./data/raw/coi68.csv")

coi_sequence <- character()

for (i in 1:length(coi68_csvtofasta$accession_num)) {
  
  name <- coi68_csvtofasta$species_name[i]
  
  seq.dat <- entrez_fetch(db = "nuccore", id = coi68_csvtofasta$accession_num[i], rettype = "fasta")
  clean_name <- gsub(" ", "_", name)
  seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
  coi_sequence <- c(coi_sequence, seq.dat2)
  
  print("Fasta added.")
}

writeLines(as.character(coi_sequence), "./data/raw/metadata/coi_68_out.fasta")

### PTR ####
ptr68_csvtofasta <- read.csv("./data/raw/ptr68.csv")

ptr_sequence <- character()

for (i in 1:length(ptr68_csvtofasta$accession_num)) {
  
  name <- ptr68_csvtofasta$species_name[i]
  
  seq.dat <- entrez_fetch(db = "nuccore", id = ptr68_csvtofasta$accession_num[i], rettype = "fasta")
  clean_name <- gsub(" ", "_", name)
  seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
  ptr_sequence <- c(ptr_sequence, seq.dat2)
  
  print("Fasta added.")
}

writeLines(as.character(ptr_sequence), "./data/raw/metadata/ptr_68_out.fasta")



### CYTB ####
cytb68_csvtofasta <- read.csv("./data/raw/cytb68.csv")

cytb_sequence <- character()

for (i in 1:length(cytb68_csvtofasta$accession_num)) {
  
  name <- cytb68_csvtofasta$species_name[i]
  
  seq.dat <- entrez_fetch(db = "nuccore", id = cytb68_csvtofasta$accession_num[i], rettype = "fasta")
  clean_name <- gsub(" ", "_", name)
  seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
  cytb_sequence <- c(cytb_sequence, seq.dat2)
  
  print("Fasta added.")
}

writeLines(as.character(cytb_sequence), "./data/raw/metadata/cytb_68_out.fasta")



### ND2 ####
nd268_csvtofasta <- read.csv("./data/raw/nd268.csv")

nd268_csvtofasta <- nd268_csvtofasta %>% 
  filter(!is.na(accession_num))

nd2_sequence <- character()

for (i in 1:length(nd268_csvtofasta$accession_num)) {
  
  name <- nd268_csvtofasta$species_name[i]
  
  seq.dat <- entrez_fetch(db = "nuccore", id = nd268_csvtofasta$accession_num[i], rettype = "fasta")
  clean_name <- gsub(" ", "_", name)
  seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
  nd2_sequence <- c(nd2_sequence, seq.dat2)
  
  print("Fasta added.")
}

writeLines(as.character(nd2_sequence), "./data/raw/metadata/nd2_68_out.fasta")



### ND5 ####
nd568_csvtofasta <- read.csv("./data/raw/nd568.csv")

nd568_csvtofasta <- nd568_csvtofasta %>% 
  filter(!is.na(accession_num))

nd5_sequence <- character()

for (i in 1:length(nd568_csvtofasta$accession_num)) {
  
  name <- nd568_csvtofasta$species_name[i]
  
  seq.dat <- entrez_fetch(db = "nuccore", id = nd568_csvtofasta$accession_num[i], rettype = "fasta")
  clean_name <- gsub(" ", "_", name)
  seq.dat2 <- sub("^>.*?\n", paste0(">", clean_name, "\n"), seq.dat)
  nd5_sequence <- c(nd5_sequence, seq.dat2)
  
  print("Fasta added.")
}

writeLines(as.character(nd5_sequence), "./data/raw/metadata/nd5_68_out.fasta")
