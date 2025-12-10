library(ape)
library(phytools)
library(viridis)
library(tidyverse)

# Load the ALREADY ROOTED tree from FigTree
fish_tree <- read.tree("./data/rooted_ml_tree.nwk")
fish_tree <- ladderize(fish_tree, right = TRUE) # This will flip the tree so that the outgroup is on the top and not the bottom

# Load trait data
traits <- read.csv("./data/raw/cichlid_traits.csv")
trait_vector <- setNames(traits$diet_type, traits$species_name)
discrete <- as.factor(trait_vector[fish_tree$tip.label])

# Group feeding strategies (I picked these kinda arbitrarily)
discrete_grouped <- discrete
levels(discrete_grouped) <- list(
  "Active Predator" = c("Carnivore", "Carnivore/Piscivore", "Specialized Carnivore"),
  "Micro Predator" = c("Planktivore", "Micro-predator/Planktivore", "Carnivore/Planktivore", "Omnivore/Micro-predator"),
  "Generalist" = c("Omnivore", "Carnivore/Omnivore", "Omnivore/Carnivore", "Omnivore/Herbivore"),
  "Grazer" = c("Herbivore", "Herbivore/Grazer", "Herbivore/Omnivore", "Herbivore/Planktivore"),
  "Detritivore" = c("Herbivore/Detritivore", "Omnivore/Detritivore")
)


# Load in the location data
locations <- read.csv("./data/raw/species_locations.csv")
location_vector <- setNames(locations$location_group, locations$species_name)
location_discrete <- as.factor(location_vector[fish_tree$tip.label])

# Colors for location (colors the species names, maybe change later?)
location_cols <- c(
  "Lake_Tanganyika" = "#E41A1C",
  "Central_American_freshwaters" = "#377EB8",
  "South_American_Malagasy_rivers" = "#4DAF4A",
  "African_rivers_lagoons_widespread" = "#FF7F00",
  "Other_single_lake_endemics" = "#984EA3",
  "Estuarine_marine" = "#000000"
)
tip_location_colors <- location_cols[as.character(location_discrete[fish_tree$tip.label])]


# Ancestral state reconstruction
fit <- ace(discrete_grouped, fish_tree, type = "discrete", model = "ER")

# Color palette
cols <- setNames(viridis(length(levels(discrete_grouped)), end = 1),
                 levels(discrete_grouped))

# Plot
# png("~/Desktop/phylogeny_ASR_powerpoint.png", 
#     width = 1920, height = 1080, res = 150)  # 16:9 aspect ratio, high resolution

pdf("./final_tree/cichlid_ml_phylogeny.pdf", width = 15, height = 10)

par(mar = c(1, 5, 1, 1))
plot.phylo(fish_tree, type = "phylogram", font = 3, use.edge.length = TRUE, show.tip.label = TRUE, tip.color = tip_location_colors,
           label.offset = 0.005, y.lim = c(1, length(fish_tree$tip.label) * 0.999), 
           cex = 0.7, edge.width = 2) # Not completely sure why, but this plot function works better

# Add pie charts
nodelabels(pie = fit$lik.anc, piecol = cols, cex = 0.4)

# Add tip colors
tip_colors <- cols[discrete_grouped[fish_tree$tip.label]]
tiplabels(pch = 21, bg = tip_colors, cex = 1.3, frame = "circle")


## Legends ####

# Legend for fedding strategy
legend(x = 0, y = 15, title = "Feeding Strategy", title.font = 2, title.adj = 0.5, legend = levels(discrete_grouped), 
       fill = cols, bty = "n", cex = 1.1)

# Legend for location
legend(x = 0, y = 27, title = "Location", title.font = 2, title.adj = 0.1,
       legend = c("Lake Tanganyika", "Central America", "S. America/Madagascar", 
                  "African Rivers", "Other Lakes"),
       fill = location_cols, bty = "n", cex = 1.1)
       

dev.off() # this is just saving the png from earlier





# Making some extra figures ####

strat_count <- read.csv("./data/raw/location_feeding_strat.csv")

custom_viridis <- c(
  "active_predator" = viridis(5)[1],
  "micro_predator" = viridis(5)[2],     
  "general" = viridis(5)[3],            
  "grazer" = viridis(5)[4],             
  "detritivore" = viridis(5)[5]         
)

plot1 <- strat_count %>% 
  mutate(feeding_strategy = factor(feeding_strategy, 
                                   levels = c("active_predator", "micro_predator", 
                                              "general", "grazer", "detritivore"))) %>%
  mutate(location = factor(location, levels = c("lake_taganyika", "central_america",
                                                "south_america", "african_rivers", "other"))) %>% 
  ggplot(aes(x = location, y = count, fill = feeding_strategy)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = custom_viridis, labels = c("active_predator" = "Active Predator",
                                                        "micro_predator" = "Micro Predator",
                                                        "general" = "Generalist",
                                                        "grazer" = "Grazer",
                                                        "detritivore" = "Detritivore")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Feeding Strategies by Location",
       x = "Location",
       y = "Count",
       fill = "Feeding Strategy")

plot1
