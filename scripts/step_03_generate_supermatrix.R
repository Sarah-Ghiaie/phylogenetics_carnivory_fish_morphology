library(evobiR)
library(ape)
library(msaR)
library(phangorn)
library(tools)

# Load in the alignment files that were got from online MAFFT
gene <- read.FASTA("./data/aligned/coi_68_align.fasta")
gene <- read.FASTA("./data/aligned/ptr_68_align.fasta")
gene <- read.FASTA("./data/aligned/cytb_68_align.fasta") 
gene <- read.FASTA("./data/aligned/nd2_68_align.fasta")
gene <- read.FASTA("./data/aligned/nd5_68_align.fasta")
length(gene)

msaR(gene)
length(unique(names(gene)))

gene_unique <- gene[!duplicated(names(gene))]

# trimming
gene2 <- phyDat(gene_unique, type = "DNA")

gene3 <- gene2[, colMeans(as.character(gene2) == "-") < 0.2] # Percentage value can change depending on the file being trimmed.
gene4 <- as.DNAbin(gene3)
msaR(gene4)

write.FASTA(gene4, "./data/trimmed/not_filled/coi_68_trimmed.fasta") # Change output file to the correct gene


# Filling in FASTA files ####
# Define trimmed FASTA files
files <- c(
  "./data/trimmed/not_filled/coi_68_trimmed.fasta",
  "./data/trimmed/not_filled/cytb_68_trimmed.fasta",
  "./data/trimmed/not_filled/ptr_68_trimmed.fasta",
  "./data/trimmed/not_filled/nd2_68_trimmed.fasta",
  "./data/trimmed/not_filled/nd5_68_trimmed.fasta"
)

# Function to read a FASTA file into a named vector
read_fasta_base <- function(filepath) {
  lines <- readLines(filepath)
  headers <- grep("^>", lines)
  seqs <- sapply(seq_along(headers), function(i) {
    start <- headers[i] + 1
    end <- ifelse(i < length(headers), headers[i + 1] - 1, length(lines))
    paste(lines[start:end], collapse = "")
  })
  names(seqs) <- gsub("^>", "", lines[headers])
  return(seqs)
}

# Function to write back to FASTA
write_fasta_base <- function(seqs, filepath) {
  out <- character(length(seqs) * 2)
  out[seq(1, length(out), by = 2)] <- paste0(">", names(seqs))
  out[seq(2, length(out), by = 2)] <- seqs
  writeLines(out, filepath)
}

# Read all alignments
alns <- lapply(files, read_fasta_base)
names(alns) <- basename(files)

# Get all taxa across all files
all_taxa <- sort(unique(unlist(lapply(alns, names))))

# Fill missing taxa
filled_alns <- lapply(alns, function(aln) {
  aln_len <- unique(nchar(aln))
  if (length(aln_len) != 1) stop("Alignment has variable lengths!")
  aln_len <- aln_len[1]
  
  missing <- setdiff(all_taxa, names(aln))
  
  if (length(missing) > 0) {
    placeholder <- rep(paste(rep("N", aln_len), collapse = ""), length(missing))
    names(placeholder) <- missing
    aln <- c(aln, placeholder)
  }
  
  aln[all_taxa]  # consistent ordering
})

# Write filled FASTAs
for (i in seq_along(filled_alns)) {
  out_path <- file.path("./data/raw/trimmed/filled/",
                        paste0(file_path_sans_ext(basename(files[i])), ".fasta"))
  write_fasta_base(filled_alns[[i]], out_path)
  cat("Wrote:", out_path, "\n")
}

# Show how many taxa were missing per gene
for (i in seq_along(alns)) {
  missing <- setdiff(all_taxa, names(alns[[i]]))
  cat(basename(files[i]), "was missing", length(missing), "taxa\n")
}


# Generate supermatrix ####
# supermatrix from evobiR; use in the directory with all the aligned_sequences
# prefix = file name prefix
# Change your working directory to/ the "trimmed/filled" directory before running
SuperMatrix(missing = "-", prefix = "species_supermatrix", save = T)
