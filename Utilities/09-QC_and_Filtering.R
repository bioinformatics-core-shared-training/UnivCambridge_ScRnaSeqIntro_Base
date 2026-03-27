#!/usr/bin/env Rscript
#SBATCH -J QC_and_Filtering
#SBATCH -o QC_and_Filtering.%j.out
#SBATCH -e QC_and_Filtering.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script first runs QC and filtering of genes and cells on all of the
# samples. It then downsamples this data set to just 500 cells per sample.

library(DropletUtils)
library(scater)
library(rtracklayer)
library(BiocParallel)
library(tidyverse)
bp.params <- MulticoreParam(workers = 7)

# load sample sheet
samplesheet <- read_tsv("Data/sample_sheet.tsv")

# create vector of file paths with names
list_of_files <- str_c("Data/Cellranger_Outputs/", 
                       samplesheet$Sample, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samplesheet$Sample

# load data
sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)

# add meta data 
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), samplesheet, by="Sample", sort=FALSE)
rownames(colData(sce)) <- sce$Barcode

# Keep only detected genes
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes,]

# Annotate genes

# or just get the info from the GTF
gtf <- "Data/references/gencode.v41.primary_assembly.annotation.mod_id.gtf.gz"
annot <- readGFF(gtf, 
                 columns="seqid", 
                 tags="gene_id", 
                 filter=list(type="gene")) %>%
    rename(ID = gene_id, Chromosome = seqid)

rowData(sce) <- merge(rowData(sce), annot, by = "ID", sort=FALSE, all.x=TRUE)
rownames(rowData(sce)) <- rowData(sce)$ID


# qc addPerCellQC
is.mito <- rowData(sce)$Chromosome%in%c("chrM", "MT")
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)


# adpative thresholds with Sample as batch
cell_qc_results <- quickPerCellQC(colData(sce),
                                percent_subsets=c("subsets_Mito_percent"),
                                batch=sce$Sample)

sce <- sce[, !cell_qc_results$discard]

# Keep only detected genes
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes,]

# update cell QC metrics
colData(sce) <- colData(sce)[,1:4]
is.mito <- rowData(sce)$Chromosome%in%c("chrM", "MT")
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

saveRDS(sce, "R_objects/Caron_filtered.full.rds")

# Now sub-sample to 500 cells per sample

set.seed(63278)
barcodes <- colData(sce) %>%
	data.frame() %>%
	group_by(SampleName) %>%
	sample_n(500) %>%
    pull(Barcode)	

sce.sub <- sce[,barcodes]

# for each gene in each cell: is it expressed?
detected_genes <- rowSums(counts(sce.sub)) > 0
sce.sub <- sce.sub[detected_genes, ]

# update cell QC metrics
colData(sce.sub) <- colData(sce.sub)[,1:4]
is.mito <- rowData(sce.sub)$Chromosome%in%c("chrM", "MT")
sce.sub <- addPerCellQC(sce.sub, subsets=list(Mito=is.mito), BPPARAM = bp.params)

# Write object to file
saveRDS(sce.sub, "R_objects/Caron_filtered.500.rds")

