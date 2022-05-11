#!/usr/bin/env Rscript

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")
packages <- c("ggplot2", "ggthemes", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser()
parser$add_argument("-o", "--output_folder", type="character",
                    help="Path to output folder for PR-curve pipeline")
parser$add_argument("-e", "--essentials_file", type="character",
                    help="Path to file containing a list of essential genes")
parser$add_argument("-c", "--corum_file", type="character",
                    help="Path to file containing a list of corum complexes")
parser$add_argument("-r", "--run_correlation", action="store_true",
                    help="Whether or not to run correlation analyses")
args <- parser$parse_args()

# Sets important parameters
pr_folder <- args$output_folder
essentials_file <- args$essentials_file
corum_file <- args$corum_file
run_cor <- args$run_correlation

######
# MAIN SCRIPT
######

# Sets important parameters
data_folder <- pr_folder
output_folder <- file.path(data_folder, "output")

# Reads in AUPRC and gene info data
ae_auprc <- read.delim(file.path(data_folder, "ae_complex_AUPRC.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
real_auprc <- read.delim(file.path(data_folder, "real_complex_AUPRC.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
olf_auprc <- read.delim(file.path(data_folder, "olf_complex_AUPRC.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pca_auprc <- read.delim(file.path(data_folder, "pca_complex_AUPRC.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
essentials <- readLines(essentials_file)
corum <- read.delim(corum_file, stringsAsFactors = FALSE, sep = "\t", header = TRUE)

# Makes output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Joins AUPRC files
common_id <- unique(Reduce(intersect, list(ae_auprc$ID, real_auprc$ID, olf_auprc$ID, pca_auprc$ID)))
auprc <- ae_auprc[ae_auprc$ID %in% common_id, 1:4]
auprc <- auprc[order(as.numeric(auprc$ID)),]
real_auprc <- real_auprc[real_auprc$ID %in% common_id,]
olf_auprc <- olf_auprc[olf_auprc$ID %in% common_id,]
pca_auprc <- pca_auprc[pca_auprc$ID %in% common_id,]
colnames(auprc)[4] <- "AE_AUPRC"
auprc$REAL_AUPRC <- real_auprc$AUPRC[order(real_auprc$ID)]
auprc$OLF_AUPRC <- olf_auprc$AUPRC[order(olf_auprc$ID)]
auprc$PCA_AUPRC <- pca_auprc$AUPRC[order(pca_auprc$ID)]
auprc$AE_VS_REAL_DIFF <- auprc$AE_AUPRC - auprc$REAL_AUPRC
auprc <- auprc[order(auprc$AE_VS_REAL_DIFF, decreasing = TRUE),]

# Writes joined AUPRC values to file
write.table(auprc, file.path(output_folder, "all_auprc.tsv"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Runs correlation-based analyses if specified
if (run_cor) {

  # Reads in normalized data
  broad <- read.table(file.path(data_folder, "real_broad.tsv"), sep = "\t", header = TRUE, row.names = 1)
  ae <- read.table(file.path(data_folder, "normalized_ae.tsv"), sep = "\t", header = TRUE, row.names = 1)
  olf <- read.table(file.path(data_folder, "normalized_olf.tsv"), sep = "\t", header = TRUE, row.names = 1)

  # Slightly changes the first value value in real data with all of the same values
  ind <- which(apply(broad, 1, function(x) length(unique(x)) == 1))
  noise <- rnorm(length(ind) * ncol(broad), mean = 0, sd = 0.00001)
  broad[ind,1] <- broad[ind,1] + 0.0001

  # Gets mean PCCs
  fitness <- rowMeans(broad)

  # Gets correlation matrices
  real_cor <- cor(t(broad))
  ae_cor <- cor(t(ae))
  olf_cor <- cor(t(olf))

  # Gets top correlations in normalized data and the corresponding real correlations
  n_top <- 50000
  cor_na <- ae_cor
  cor_na[lower.tri(cor_na, diag = TRUE)] <- NA
  norm_top_cor <- order(cor_na, na.last = TRUE, decreasing = TRUE)[1:n_top]
  all_top_cor <- data.frame(gene1 = rep(NA, length(norm_top_cor)),
                            gene2 = rep(NA, length(norm_top_cor)),
                            gene1_essential = rep(NA, length(norm_top_cor)),
                            gene2_essential = rep(NA, length(norm_top_cor)),
                            ae_cor = rep(NA, length((norm_top_cor))),
                            olf_cor = rep(NA, length((norm_top_cor))),
                            real_cor = rep(NA, length((norm_top_cor))),
                            fitness_effect = rep(NA, length((norm_top_cor))))
  for (i in 1:nrow(all_top_cor)) {
    cor_ind <- arrayInd(norm_top_cor[i], dim(ae_cor))
    genes <- sort(c(rownames(ae_cor)[cor_ind[1,1]], colnames(ae_cor)[cor_ind[1,2]]))
    all_top_cor$gene1[i] <- genes[1]
    all_top_cor$gene2[i] <- genes[2]
    all_top_cor$gene1_essential[i] <- genes[1] %in% essentials
    all_top_cor$gene2_essential[i] <- genes[2] %in% essentials
    all_top_cor$ae_cor[i] <- ae_cor[cor_ind[1,1], cor_ind[1,2]]
    all_top_cor$olf_cor[i] <- olf_cor[cor_ind[1,1], cor_ind[1,2]]
    all_top_cor$real_cor[i] <- real_cor[cor_ind[1,1], cor_ind[1,2]]
    all_top_cor$fitness_effect[i] <- mean(c(fitness[genes[1]], fitness[genes[2]]))
  }

  # Gets all CORUM pairs
  co_complex_pairs <- data.frame(gene1 = rep(NA, 100000),
                                gene2 = rep(NA, 100000))
  counter <- 1
  for (i in 1:nrow(corum)) {
    genes <- toupper(strsplit(corum$Subunits_gene_name[i], ";")[[1]])
    genes <- genes[genes %in% rownames(ae_cor)]
    if (length(genes) < 2) {
      next
    }
    pairs <- combn(genes, 2)
    for (j in 1:ncol(pairs)) {
      co_complex_pairs[counter,] <- sort(pairs[,j])
      counter <- counter + 1
    }
  }
  co_complex_pairs <- co_complex_pairs[complete.cases(co_complex_pairs),]
  co_complex_pairs <- co_complex_pairs[!duplicated(co_complex_pairs),]
  cat(paste("Co-complex size of", nrow(co_complex_pairs), "pairs\n"))
  co_complex_unique <- paste0(co_complex_pairs$gene1, co_complex_pairs$gene2)

  # Marks CORUM co-complex pairs in list of top correlations
  all_top_cor$corum <- FALSE
  all_top_cor_unique <- paste0(all_top_cor$gene1, all_top_cor$gene2)
  all_top_cor$corum[all_top_cor_unique %in% co_complex_unique] <- TRUE

  # Writes data to file
  write.table(all_top_cor, file.path(output_folder, "top_cor.tsv"), sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}
