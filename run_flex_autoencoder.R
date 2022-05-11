#!/usr/bin/env Rscript

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6")
packages <- c("stringi", "devtools", "RColorBrewer", "argparse", "FLEX")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser()
parser$add_argument("-o", "--output_folder", type="character", default="output", 
                    help="Path to output folder containing latent space [default %(default)]")
parser$add_argument("-m", "--mito_removal", action="store_true",
                    help="Whether or not to remove mitochondrial complexes")
args <- parser$parse_args()

# Sets important parameters
output_folder <- args$output_folder
remove_mito <- args$mito_removal
pr_folder <- file.path(output_folder, "pr")

######
# MAIN SCRIPT
######

# Sets paths to depmap files
real_file <- file.path(pr_folder, "real_broad.tsv")
ae_file <- file.path(pr_folder, "normalized_ae.tsv")
olf_file <- file.path(pr_folder, "normalized_olf.tsv")
pca_file <- file.path(pr_folder, "normalized_pca.tsv")
latent_file <- file.path(pr_folder, "ae_latent.tsv")

# Lists mitochondrial complexes
mito_complex_ids <- c("320", "178")

# Reads in FLEX data
data('data_complex', package = 'FLEX')
corum <- MakeCoAnnotationFromGeneSymbols(data_standard = data_complex, overlap_length = 1)

# Reads KEGG co-pathway data
# kegg <- MakeCoAnnotationFromGeneSymbols(data_standard = data_pathway, overlap_length = 1, subset_str = c())

# Gets profile similarities of real data
real_interaction <- GetInteractionData(real_file, delim = "\t")
real_complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, real_interaction)
temp <- as.data.frame(real_complex, stringsAsFactors = FALSE)
if (remove_mito) { temp <- getSubsetOfCoAnnRemoveIDs(data.standard = temp, ids = mito_complex_ids, replace = FALSE) }
pred <- list(real = list(true = temp$true, predicted = temp$predicted))

# Gets profile similarities of GAN-normalized data
ae_interaction <- GetInteractionData(ae_file, delim = "\t")
ae_complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, ae_interaction)
temp <- as.data.frame(ae_complex, stringsAsFactors = FALSE)
if (remove_mito) { temp <- getSubsetOfCoAnnRemoveIDs(data.standard = temp, ids = mito_complex_ids, replace = FALSE) }
pred <- append(pred, list(ae = list(true = temp$true, predicted = temp$predicted)))

# Gets profile similarities of Olf-normalized data
olf_interaction <- GetInteractionData(olf_file, delim = "\t")
olf_complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, olf_interaction)
temp <- as.data.frame(olf_complex, stringsAsFactors = FALSE)
if (remove_mito) { temp <- getSubsetOfCoAnnRemoveIDs(data.standard = temp, ids = mito_complex_ids, replace = FALSE) }
pred <- append(pred, list(olf = list(true = temp$true, predicted = temp$predicted)))

# Gets profile similarities of PCA-normalized data
pca_interaction <- GetInteractionData(pca_file, delim = "\t")
pca_complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, pca_interaction)
temp <- as.data.frame(pca_complex, stringsAsFactors = FALSE)
if (remove_mito) { temp <- getSubsetOfCoAnnRemoveIDs(data.standard = temp, ids = mito_complex_ids, replace = FALSE) }
pred <- append(pred, list(pca = list(true = temp$true, predicted = temp$predicted)))

# Sets certain plotting parameters
curve_names <- c('Real Depmap', 'AE-norm', 'Olf-norm', 'PCA-norm')
cols <- brewer.pal(4, "Set1")

# Gets profile similarities of AE latent space and modifies plotting parameters accordingly
latent_interaction <- GetInteractionData(latent_file, delim = "\t")
if (!is.null(dim(latent_interaction))) {
  latent_complex <- CalculatePredictionAndTrueOnLibraryProfiles(corum, latent_interaction)
  temp <- as.data.frame(latent_complex, stringsAsFactors = FALSE)
  if (remove_mito) { temp <- getSubsetOfCoAnnRemoveIDs(data.standard = temp, ids = mito_complex_ids, replace = FALSE) }
  pred <- append(pred, list(latent = list(true = temp$true, predicted = temp$predicted)))
  curve_names <- c('Real Depmap', 'AE-norm', 'Olf-norm', 'PCA-norm', 'Latent')
  cols <- brewer.pal(5, "Set1")
}

# Plots PR curves with colors from: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=6
title <- 'CORUM profile similarity'
labs <- c('TP', 'Precision')
fname <- file.path(pr_folder, paste0("corum_pr_", basename(output_folder)))
PlotPRSimilarity(pred, fig.title = title, fig.labs = labs, legend.names = curve_names, legend.color = cols, save.figure = TRUE,
                 outfile.name = fname, outfile.type = "png")

# Makes contribution-by-size plot for real data
temp <- as.data.frame(real_complex, stringsAsFactors = FALSE)
size_contrib <- GetAreaUnderPRCurveForEntities(data_complex, corum, temp)
PlotContributionScatter(size_contrib, fig.title = 'Real contribution scatter', fig.labs = c('AUPRC', 'Complex size'), 
                        show.text = FALSE, save.figure = TRUE, outfile.name = file.path(pr_folder, "real_contribution_scatter"),
                        outfile.type = "png")
write.table(size_contrib, file.path(pr_folder, 'real_complex_AUPRC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# GAN-normalized data
temp <- as.data.frame(ae_complex, stringsAsFactors = FALSE)
size_contrib <- GetAreaUnderPRCurveForEntities(data_complex, corum, temp)
PlotContributionScatter(size_contrib, fig.title = 'AE contribution scatter', fig.labs = c('AUPRC', 'Complex size'), 
                        show.text = FALSE, save.figure = TRUE, outfile.name = file.path(pr_folder, "ae_contribution_scatter"),
                        outfile.type = "png")
write.table(size_contrib, file.path(pr_folder, 'ae_complex_AUPRC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# Olf-normalized data
temp <- as.data.frame(olf_complex, stringsAsFactors = FALSE)
size_contrib <- GetAreaUnderPRCurveForEntities(data_complex, corum, temp)
PlotContributionScatter(size_contrib, fig.title = 'Olf contribution scatter', fig.labs = c('AUPRC', 'Complex size'), 
                        show.text = FALSE, save.figure = TRUE, outfile.name = file.path(pr_folder, "olf_contribution_scatter"),
                        outfile.type = "png")
write.table(size_contrib, file.path(pr_folder, 'olf_complex_AUPRC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# PCA-normalized data
temp <- as.data.frame(pca_complex, stringsAsFactors = FALSE)
size_contrib <- GetAreaUnderPRCurveForEntities(data_complex, corum, temp)
PlotContributionScatter(size_contrib, fig.title = 'PCA contribution scatter', fig.labs = c('AUPRC', 'Complex size'), 
                        show.text = FALSE, save.figure = TRUE, outfile.name = file.path(pr_folder, "pca_contribution_scatter"),
                        outfile.type = "png")
write.table(size_contrib, file.path(pr_folder, 'pca_complex_AUPRC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

# Latent space
if (!is.null(dim(latent_interaction))) {
  temp <- as.data.frame(latent_complex, stringsAsFactors = FALSE)
  size_contrib <- GetAreaUnderPRCurveForEntities(data_complex, corum, temp)
  PlotContributionScatter(size_contrib, fig.title = 'Latent contribution scatter', fig.labs = c('AUPRC', 'Complex size'), 
                          show.text = FALSE, save.figure = TRUE, outfile.name = file.path(pr_folder, "pca_contribution_scatter"),
                          outfile.type = "png")
  write.table(size_contrib, file.path(pr_folder, 'latent_complex_AUPRC.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
}

# Makes stepwise contribution plot for real data
pairs_in <- data.frame(true = real_complex$true, 
                       predicted = real_complex$predicted, 
                       ID = real_complex$ID, 
                       stringsAsFactors = FALSE)
stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = real_complex$predicted, 
                                             value.true = real_complex$true, 
                                             x.axis = 'TP', y.axis = 'precision')
precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, 1, 0.025))
precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex)
write.table(stepwise_contrib, file.path(pr_folder, 'real_stepwise_contribution.tsv'), sep="\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, 
                          min.pairs = 10, fig.title = 'Real glacier plot', 
                          save.figure = TRUE, outfile.name = file.path(pr_folder, "real_structure"),
                          outfile.type = "png")

# AE-normalized data
pairs_in <- data.frame(true = ae_complex$true, 
                       predicted = ae_complex$predicted, 
                       ID = ae_complex$ID, 
                       stringsAsFactors = FALSE)
stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = ae_complex$predicted, 
                                             value.true = ae_complex$true, 
                                             x.axis = 'TP', y.axis = 'precision')
precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, , 0.025))
precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex)
write.table(stepwise_contrib, file.path(pr_folder, 'ae_stepwise_contribution.tsv'), sep="\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, 
                          min.pairs = 10, fig.title = 'AE glacier plot', 
                          save.figure = TRUE, outfile.name = file.path(pr_folder, "ae_structure"),
                          outfile.type = "png")

# Olf-normalized data
pairs_in <- data.frame(true = olf_complex$true, 
                       predicted = olf_complex$predicted, 
                       ID = olf_complex$ID, 
                       stringsAsFactors = FALSE)
stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = olf_complex$predicted, 
                                             value.true = olf_complex$true, 
                                             x.axis = 'TP', y.axis = 'precision')
precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, 1, 0.025))
precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex)
write.table(stepwise_contrib, file.path(pr_folder, 'olf_stepwise_contribution.tsv'), sep="\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, 
                          min.pairs = 10, fig.title = 'Olf glacier plot', 
                          save.figure = TRUE, outfile.name = file.path(pr_folder, "olf_structure"),
                          outfile.type = "png")

# PCA-normalized data
pairs_in <- data.frame(true = pca_complex$true, 
                       predicted = pca_complex$predicted, 
                       ID = pca_complex$ID, 
                       stringsAsFactors = FALSE)
stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = pca_complex$predicted, 
                                             value.true = pca_complex$true, 
                                             x.axis = 'TP', y.axis = 'precision')
precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, 1, 0.025))
precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex)
write.table(stepwise_contrib, file.path(pr_folder, 'pca_stepwise_contribution.tsv'), sep="\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, 
                          min.pairs = 10, fig.title = 'PCA glacier plot', 
                          save.figure = TRUE, outfile.name = file.path(pr_folder, "pca_structure"),
                          outfile.type = "png")

# Latent space
if (!is.null(dim(latent_interaction))) {
  pairs_in <- data.frame(true = latent_complex$true, 
                        predicted = latent_complex$predicted, 
                        ID = latent_complex$ID, 
                        stringsAsFactors = FALSE)
  stepwise_contrib <- GenerateDataForPerfCurve(value.predicted = latent_complex$predicted, 
                                              value.true = latent_complex$true, 
                                              x.axis = 'TP', y.axis = 'precision')
  precision_cutoffs <- c(stepwise_contrib$y[length(stepwise_contrib$y)], seq(0.01, 1, 0.025))
  precision_cutoffs[1] <- round(precision_cutoffs[1], 4)
  stepwise_contrib <- GetStepwiseContributionOfEntities(pairs_in, precision_cutoffs, data_complex)
  write.table(stepwise_contrib, file.path(pr_folder, 'pca_stepwise_contribution.tsv'), sep="\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  stepwise_contrib <- stepwise_contrib[!duplicated(stepwise_contrib$Name), ]
  PlotContributionStructure(stepwise_contrib, cutoff.all = precision_cutoffs, 
                            min.pairs = 10, fig.title = 'Latent glacier plot', 
                            save.figure = TRUE, outfile.name = file.path(pr_folder, "latent_structure"),
                            outfile.type = "png")
}