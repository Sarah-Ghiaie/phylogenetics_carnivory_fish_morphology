library(ape)
library(phytools)
library(viridis)

# Load tree
tree <- read.nexus("~/Desktop/Jannuzzi-4300/learning_beast2/beast/bin/80mil_output.tree")

# Resolve polytomies FIRST (before any rooting)
tree <- multi2di(tree)

# NOW root it
root_tree <- root(tree, outgroup = "Planiliza_haematocheilus", resolve.root=TRUE)

# Load trait data
traits <- read.csv("~/Desktop/cichlid_traits.csv")
trait_vector <- setNames(traits$diet_type, traits$species_name)

# Trim tree to match data
to_remove <- setdiff(root_tree$tip.label, names(trait_vector))
trimmed_tree <- drop.tip(root_tree, to_remove)

# Convert to factor
discrete <- as.factor(trait_vector[trimmed_tree$tip.label])

# GROUP INTO FEEDING STRATEGIES
discrete_grouped <- discrete
levels(discrete_grouped) <- list(
  "Active_Predator" = c("Carnivore", "Carnivore/Piscivore", "Specialized Carnivore"),
  "Micro_Predator" = c("Planktivore", "Micro-predator/Planktivore", "Carnivore/Planktivore", "Omnivore/Micro-predator"),
  "Generalist" = c("Omnivore", "Carnivore/Omnivore", "Omnivore/Carnivore", "Omnivore/Herbivore"),
  "Grazer" = c("Herbivore", "Herbivore/Grazer", "Herbivore/Omnivore", "Herbivore/Planktivore"),
  "Detritivore" = c("Herbivore/Detritivore", "Omnivore/Detritivore")
)

# Check if tree is binary and rooted
is.binary(trimmed_tree)
is.rooted(trimmed_tree)

# Run ancestral state reconstruction
fit <- ace(discrete_grouped, trimmed_tree, type="discrete", model="ER")

# Set up colors
cols <- setNames(viridis(length(levels(discrete_grouped)), end = 0.8),
                 levels(discrete_grouped))

# Plot tree with pie charts
par(mar=c(1,1,1,1))
plotTree(trimmed_tree, ftype="i", fsize=0.7, lwd=2,
         xlim=c(0, max(nodeHeights(trimmed_tree))*1.3))

# Add pie charts at internal nodes
nodelabels(pie=fit$lik.anc, piecol=cols, cex=0.6)

# Add colored circles at tips
tip_colors <- cols[discrete_grouped[trimmed_tree$tip.label]]
tiplabels(pch=21, bg=tip_colors, cex=1, frame="circle")

# Add legend
legend("bottomleft", title = "Feeding Strategy", levels(discrete_grouped), fill=cols, bty="n", cex=0.9)