#!/usr/bin/env Rscript

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")
packages <- c("ggplot2", "ggthemes", "pheatmap", "argparse", "gplots")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", type="character",
                    help="Path to training data")
parser$add_argument("-f", "--olf_file", type="character",
                    help="Path to file containing a list of olfactory receptor genes")
parser$add_argument("-n", "--n_components", type="integer", default=5,
                    help="Number of principal components to remove")
parser$add_argument("-e", "--essentials_file", type="character",
                    help="Path to file containing a list of essential genes")
parser$add_argument("-c", "--corum_file", type="character",
                    help="Path to file containing a list of corum complexes")
parser$add_argument("-o", "--output_folder", type="character", default="output", 
                    help="Path to output folder containing latent space [default %(default)]")
args <- parser$parse_args()

# Sets important parameters
input_file <- args$input_file
olf_file <- args$olf_file
n_components <- args$n_components
essentials_file <- args$essentials_file
corum_file <- args$corum_file
output_folder <- args$output_folder
plot_folder <- file.path(output_folder, "plots")
pr_folder <- file.path(output_folder, "pr")

######
# MAIN SCRIPT
######

# Sets additional parameters
n_profiles <- 250
remove_essential <- FALSE
only_essential <- FALSE
write_cdt <- FALSE
remove_outliers <- TRUE
breaks <- seq(-1, 1, by = (1/150))
pal <- colorRampPalette(c("blue", "black", "yellow"))(n = length(breaks))

if (remove_essential == TRUE & only_essential == TRUE) { 
  stop("Can't keep only essentials and also remove essentials") 
}

# Makes output folder if nonexistent
if (!dir.exists(output_folder)) { dir.create(output_folder) }
if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
if (!dir.exists(pr_folder)) { dir.create(pr_folder) }

# Reads in and formats data
ae_cor <- as.numeric(readLines(file.path(output_folder, "autoencoder_latent_cor.txt")))
ae_latent <- read.table(file.path(output_folder, "autoencoder_latent_mat.tsv"), sep = "\t", stringsAsFactors = FALSE)
ae_profiles <- read.table(file.path(output_folder, "autoencoder_fake_profiles.tsv"), sep = "\t", stringsAsFactors = FALSE)
olf <- read.csv(olf_file, colClasses = c(Approved.symbol ="character"))$Approved.symbol
data <- read.csv(file.path(input_file), sep = "\t", header = TRUE, row.names = 1)
corum <- read.delim(file.path(corum_file), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
essentials <- readLines(file.path(essentials_file))

# Subsets depmap and removes extreme values
ae_profiles <- ae_profiles[,1:ncol(data)]
names(ae_cor) <- rownames(data)
rownames(ae_latent) <- rownames(data)
rownames(ae_profiles) <- rownames(data)
colnames(ae_latent) <- paste("Latent", 1:ncol(ae_latent))

# Removes essential genes if specified
if (remove_essential) {
  keep_ind <- !(rownames(data) %in% essentials)
  data <- data[keep_ind,]
  ae_cor <- ae_cor[keep_ind]
  ae_latent <- ae_latent[keep_ind]
}

# Subsets and formats the rest of the data
olf <- rownames(data)[rownames(data) %in% olf]
broad_olf <- data[olf,]
essentials <- rownames(data)[rownames(data) %in% essentials]
nonessentials <- rownames(data)[!(rownames(data) %in% essentials)]
nonolf <- nonessentials[!(nonessentials %in% olf)]
broad_nonolf <- data[sample(nonolf, nrow(broad_olf)),]

# Makes histogram of latent correlations
temp <- data.frame(cor = as.numeric(ae_cor))
ggplot(temp, aes(cor)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept = -0.25, linetype = 2, color = "blue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray", alpha = 0.7) +
  geom_vline(xintercept = 0.25, linetype = 2, color = "red", alpha = 0.7) +
  xlab("PCC between synthetic data and CERES scores") +
  ylab("Count") +
  theme_tufte(base_size = 18)
ggsave(file.path(plot_folder, "cor_histogram.png"))

# Writes real data to cdt file
if (write_cdt) {
  hc <- hclust(dist(t(data)))
  hr <- hclust(dist(data))
  r2gtr(hc, file=file.path(plot_folder, "real_cluster.atr"))
  r2atr(hr, file=file.path(plot_folder, "real_cluster.gtr"))
  r2cdt(hr, hc, data, file=file.path(plot_folder, "real_cluster.cdt"))
}

# Normalizes data to autoencoder latent space
broad_ae <- data - as.matrix(ae_profiles)
dimnames(broad_ae) <- dimnames(data)

# Writes AE-normalized data to file
if (write_cdt) {
  hc <- hclust(dist(t(ae_profiles)))
  hr <- hclust(dist(ae_profiles))
  r2gtr(hc, file=file.path(plot_folder, "ae_cluster.atr"))
  r2atr(hr, file=file.path(plot_folder, "ae_cluster.gtr"))
  r2cdt(hr, hc, ae_profiles, file=file.path(plot_folder, "ae_cluster.cdt"))
}

# Runs PCA on data as-is
print("Running PCA on real data...")
pca <- prcomp(data, center = TRUE, scale. = FALSE)
real_v <- pca$rotation[,1:n_components]
projected <- as.matrix(data) %*% real_v %*% t(real_v)
broad_pca <- data - projected
dimnames(broad_pca) <- dimnames(data)

# Writes data projected onto real PCs to file
if (write_cdt) {
  hc <- hclust(dist(t(projected)))
  hr <- hclust(dist(projected))
  r2gtr(hc, file=file.path(plot_folder, "pca_cluster.atr"))
  r2atr(hr, file=file.path(plot_folder, "pca_cluster.gtr"))
  r2cdt(hr, hc, projected, file=file.path(plot_folder, "pca_cluster.cdt"))
}

# Runs PCA on data with olfactory receptors as null model
print("Running PCA with Olf null model...")
pca <- prcomp(broad_olf, center = TRUE, scale. = FALSE)
olf_v <- pca$rotation[,1:min(n_components, length(olf) - 1)]
#projected <- scale(data, pca$center, pca$scale) %*% v %*% t(v) 
projected <- as.matrix(data) %*% olf_v %*% t(olf_v)
broad_olf_pca <- data - projected
dimnames(broad_olf_pca) <- dimnames(data)

# Writes data projected onto olf PCs to file
if (write_cdt) {
  hc <- hclust(dist(t(projected)))
  hr <- hclust(dist(projected))
  r2gtr(hc, file=file.path(plot_folder, "olf_cluster.atr"))
  r2atr(hr, file=file.path(plot_folder, "olf_cluster.gtr"))
  r2cdt(hr, hc, projected, file=file.path(plot_folder, "olf_cluster.cdt"))
}

# Runs PCA on data with non-olfactory nonessential genes as null model
print("Running PCA with non-essential, non-Olf null model...")
pca <- prcomp(data[nonolf,], center = FALSE, scale. = FALSE)
nonolf_v <- pca$rotation[,1:n_components]
projected <- as.matrix(data) %*% nonolf_v %*% t(nonolf_v)
broad_nonolf_pca <- data - projected
dimnames(broad_nonolf_pca) <- dimnames(data)

# Writes data projected onto olf PCs to file
if (write_cdt) {
  hc <- hclust(dist(t(projected)))
  hr <- hclust(dist(projected))
  r2gtr(hc, file=file.path(plot_folder, "nonolf_cluster.atr"))
  r2atr(hr, file=file.path(plot_folder, "nonolf_cluster.gtr"))
  r2cdt(hr, hc, projected, file=file.path(plot_folder, "nonolf_cluster.cdt"))
}

# Optionally removes essential genes or keeps only essential genes
if (remove_essential) {
  data <- data[nonessentials,]
  broad_pca <- broad_pca[nonessentials,]
  broad_ae <- broad_ae[nonessentials,]
  broad_olf_pca <- broad_olf_pca[nonessentials,]
} else if (only_essential) {
  data <- data[essentials,]
  broad_pca <- broad_pca[essentials,]
  broad_ae <- broad_ae[essentials,]
  broad_olf_pca <- broad_olf_pca[essentials,]
}

# Writes data to file
write.table(data, file.path(pr_folder, "real_broad.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(broad_pca, file.path(pr_folder, "normalized_pca.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(broad_ae, file.path(pr_folder, "normalized_ae.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(broad_olf_pca, file.path(pr_folder, "normalized_olf.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(ae_latent, file.path(pr_folder, "ae_latent.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

# Makes heatmaps of profiles
if (!only_essential) {
  png(file.path(plot_folder, "ceres_profile_nonessential_heatmap.png"), width = 8, height = 6, units = "in", res = 800)
  heatmap.2(as.matrix(data[sample(nonessentials, n_profiles),]),
            col = pal, trace = "none", key.xlab = NA, key.ylab = NA, key.title = NA,
            labRow = NA, labCol = NA)
  dev.off()
  olf_set <- sample(olf, min(n_profiles, length(olf)))
  png(file.path(plot_folder, "ceres_profile_olf_heatmap.png"), width = 8, height = 6, units = "in", res = 800)
  heatmap.2(as.matrix(data[nonessentials[nonessentials %in% olf_set],]),
            col = pal, trace = "none", key.xlab = NA, key.ylab = NA, key.title = NA,
            labRow = NA, labCol = NA)
  dev.off()
}
if (!remove_essential) {
  png(file.path(plot_folder, "ceres_profile_essential_heatmap.png"), width = 8, height = 6, units = "in", res = 800)
  heatmap.2(as.matrix(data[sample(essentials, n_profiles),]),
            col = pal, trace = "none", key.xlab = NA, key.ylab = NA, key.title = NA,
            labRow = NA, labCol = NA)
  dev.off()
}
png(file.path(plot_folder, "ceres_profile_heatmap.png"), width = 8, height = 6, units = "in", res = 800)
heatmap.2(as.matrix(data[1:n_profiles,]),
          col = pal, trace = "none", key.xlab = NA, key.ylab = NA, key.title = NA,
          labRow = NA, labCol = NA)
dev.off()
png(file.path(plot_folder, "ae_profile_heatmap.png"), width = 8, height = 6, units = "in", res = 800)
heatmap.2(as.matrix(ae_profiles[1:n_profiles,]),
          col = pal, trace = "none", key.xlab = NA, key.ylab = NA, key.title = NA,
          labRow = NA, labCol = NA)
dev.off()


