# Setup -----

# load packages
library(BiocParallel)
library(scran)
library(scater)
library(miloR)
library(tidyverse)
library(patchwork)

# load the SCE object
sce <- readRDS("R_objects/Caron_clustered.PBMMCandETV6RUNX1.rds")

# check the contents of the object
sce

# plot UMAP done on the batch-corrected data
plotReducedDim(sce, dimred = "UMAP_corrected", 
               colour_by = "label", 
               text_by = "label")


# Differential abundance with Milo ----

# create the Milo object can be simply converted from a SCE
milo <- Milo(sce)

milo


# Build KNN graph ----

# add KNN graph to Milo object
milo <- buildGraph(milo, 
                   k = 60, 
                   d = 50, 
                   reduced.dim = "corrected", 
                   BPPARAM = MulticoreParam(7))

# sample index cells to define neighbourhoods
milo <- makeNhoods(milo, 
                   prop = 0.1, 
                   k = 60, 
                   d = 50, 
                   reduced_dims = "corrected")

# check our object again
milo

# distribution of neighbourhood sizes
plotNhoodSizeHist(milo) +
  geom_vline(xintercept = 100, col = "salmon")

# count cells in each neighbourhood
milo <- countCells(milo, 
                   meta.data = colData(milo),
                   samples = "SampleName")

# Milo now has a counts matrix
head(nhoodCounts(milo))


# Run DA analysis ----

# calculate distances between neighbourhoods - for p-value correction
milo <- calcNhoodDistance(milo, d = 50, reduced.dim = "corrected")

# define a table for our model design
sample_info <- unique(colData(milo)[,c("SampleName", "SampleGroup")])
rownames(sample_info) <- sample_info$SampleName

sample_info

# run DA test
da_results <- testNhoods(milo, 
                         design = ~ SampleGroup, 
                         design.df = sample_info, 
                         reduced.dim = "corrected")

# results are returned as a data.frame
da_results %>%
  arrange(SpatialFDR) %>%
  head()


# Visualisations ----

# p-value histogram
ggplot(da_results, aes(PValue)) + 
  geom_histogram(bins = 50)

# volcano plot
# each point in this plot corresponds to a neighbourhood (not a cell)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point(aes(colour = FDR < 0.1)) +
  geom_hline(yintercept = 1) 

# build neighbourhood graph embedding
milo <- buildNhoodGraph(milo)

# our original UMAP with our previously annotated cell labels
umap_plot <- plotReducedDim(milo, 
                            dimred = "UMAP_corrected", 
                            colour_by = "label", 
                            text_by = "label")

# the neighbourhood map adjusted to match UMAP embedding
nh_graph_plot <- plotNhoodGraphDA(milo, 
                                  da_results, 
                                  layout = "UMAP_corrected",
                                  alpha = 0.05)

# the two plots together side-by-side
umap_plot + nh_graph_plot +
  plot_layout(guides="collect")


# Annotate NHoods ----

# annotate our neighbourhood DA results with our cell labels
da_results <- annotateNhoods(milo, da_results, coldata_col = "label")
head(da_results)

# histogram of fraction of cells in the neighbourhood with the same label
ggplot(da_results, aes(label_fraction)) + 
  geom_histogram(bins = 50)

# add "mixed" label to neighbourhoods with less 70% consistency
da_results$label <- ifelse(da_results$label_fraction < 0.7, 
                           "Mixed", 
                           da_results$label)

head(da_results)

# distribution of logFC across neighbourhood labels
plotDAbeeswarm(da_results, group.by = "label")

