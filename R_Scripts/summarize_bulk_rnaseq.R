############################
# Gopher pipelines version 0.1.0 bulk RNA-seq summary script
# R version 3.5.0
# All relevant R libraries are installed at /panfs/roc/groups/14/msistaff/public/gopher-pipelines/v0/R/
# Usage: Rscript summarize_counts.R <out_dir> <work_dir> <sample sheet path> <merged_counts_path> <min_feature_length> <min_cpm> &> Rout.txt
# See run_summary_stats.pbs for additional run context.
# See https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf for the edgeR user manual.
# Contact help@msi.umn.edu for questions
############################


# Prepend the library path to libPaths
.libPaths(c("/panfs/roc/groups/14/msistaff/public/CHURP_Deps/v0/R/",.libPaths()))

# Load libraries
library('limma')
library('edgeR')
library('ggplot2')
library('reshape2')
library('gplots')
library('gtools')
library('grid')

#grab the working and output directories, as well as the sample sheet, and merged raw counts matrix
args <- commandArgs(trailingOnly = T)
out_dir <- args[1]
work_dir <- args[2]
samp_sheet <- args[3]
fc_mat <- args[4]
min_len <- args[5]
min_cpm <- args[6]

setwd(work_dir)

############################
# Read in data and prep relevant data and set output files
############################

# Get the sample sheet to grab group membership downstream 
sheet <- read.table(samp_sheet, sep = "|", header = F, comment.char = "#")
sheet$V1 <- make.names(sheet$V1)

# Because there may be cases where a subset of individuals in the samplesheet are run. We'll pull in the featureCounts matrix early and grab the relevant IDs
raw_mat <- read.table(fc_mat, header = T, sep = '\t', comment.char = '#')
samp_ids <- names(raw_mat)[-(1:6)]
sheet <- sheet[which(sheet$V1 %in% samp_ids),]

# Now we can proceed with group determinations.
groups <- as.vector(sheet$V2)
uniq_groups <- unique(groups)
true_groups <- groups[which(groups != 'NULL')]
n_true_groups <- length(uniq_groups[which(uniq_groups != 'NULL')])

# Check if the number of groups meets our analysis criteria
if (n_true_groups >= 5){
write("Maximum number of groups exceeded: This experimental design appears to be too complex for this pipeline. If you'd like analysis help, please e-mail help@msi.umn.edu.", stderr())
quit(status = 0, save = "no")
}

# Then check if the number of replicates meets our criteria.
if (length(true_groups[which(table(true_groups) < 3)]) > 0){
write("At least one of the groups has fewer than three replicates, which makes reliable statistical interpretation difficult", stderr())
writeLines(text = paste("This group has fewer than three replicates", unique(groups)[which(table(groups) < 3)], sep=" : "), con = stderr())
quit(status = 0, save = "no")
}

# Set filename variables
mds_plot <- paste(out_dir, "Plots/mds_plot.pdf", sep = "/")
counts_plot <- paste(out_dir, "Plots/cpm_plot.pdf", sep = "/")
counts_list <- paste(out_dir, "Counts/cpm_list.txt", sep = "/")
hmap <- paste(out_dir, "Plots/high_variance_heatmap.pdf", sep = "/")

# Filter the featureCounts matrix (read in above) on variance and feature length.
raw_mat <- raw_mat[-which(apply(raw_mat[,seq(-1,-6)],1,var) < 1),]
raw_mat <- raw_mat[which(raw_mat$Length >= min_len),]

# Convert the raw matrix into a DGE object. Column 1 is Geneid, columns 7+ are the sample counts. The groups list  is generated above. 
edge_mat <- DGEList(counts = raw_mat[,seq(-1,-6)], genes = raw_mat[,1], group = groups)

# Calculate the normalization factors
edge_mat <- calcNormFactors(edge_mat)

############################
# Generate descriptive accounts of the data (MDS, Normalized Counts, Counts Distributions, and a Heatmap) 
############################

# set a diverging color palette. From http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6
#pal <- c('#e41a1c','#377eb8','#4daf4a','#984ea3') Nice colors but they're not colorblind friendly


# This is THE ONLY colorblind acceptable palette with 4 colors from colorbrewer: http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
pal <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c')
col_vec <- pal[match(groups,uniq_groups)]

# legend code adapted from https://support.bioconductor.org/p/101530/
# Set the MDS plot pdf and write the plot
pdf(mds_plot)
opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
plotMDS(edge_mat, cex = 0.75, col = col_vec)
legend(par("usr")[2], mean(par("usr")[3:4]), legend = c('Group', uniq_groups), text.col = c('black', unique(col_vec)), bty = "n")
par(opar)
dev.off()

## set par back to original
#par(opar)

# Set a variable holding the log(1+CPM) counts.
cpm_counts <- cpm(edge_mat, log = T, prior.count = 1)

# Create a dataframe of cpm_counts and gene IDs in 'wide' format. Melt into 'long' format for ggplot
cdf <- data.frame(edge_mat$genes, cpm_counts)
write.table(cdf, file = counts_list, sep = '\t', quote = FALSE, row.names = FALSE)
tidy_cdf <- melt(cdf, id.vars = "genes", variable.name = "sample_id", value.name = "per_feature_count")

# manually add group information to tidy_cdf
tidy_cdf$group <- factor(rep(uniq_groups[match(groups,uniq_groups)], each = nrow(edge_mat$genes)))

# Set the counts plot pdf and write the violin plot of normalized counts per sample
pdf(counts_plot)
p <- ggplot(tidy_cdf, aes(x = sample_id, y = per_feature_count, fill = group)) + geom_violin(trim = F) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7)) + scale_fill_manual("Group", values = pal)
p + labs(x = "Sample ID", y = "Feature count -- log(1+cpm)", fill = "Group")
dev.off()


# Calculate count variance across samples and select the top 500 variance features.
cpm_counts <- cpm(edge_mat, log = T, prior.count = 1)
gene_var <- apply(cpm_counts, 1, var)
select_var <- names(sort(gene_var, decreasing=TRUE))[1:500]
high_var <- cpm_counts[select_var,]

# Set the heatmap pdf and plot the normalized counts heatmap
pdf(hmap)
heatmap.2(high_var,trace="none", main = "Top 500 variance genes", cexCol = 0.75, dendrogram = "column", labRow = "", ColSideColors = col_vec, srtCol = 45, margins = c(8,8))
legend('left', title = 'Group', legend = uniq_groups, fill = unique(col_vec), cex = 0.8, box.lty = 0 )
dev.off()

############################
# Differential expression testing and summaries
############################
n_groups <- length(uniq_groups)

# If there is only one grouping in the data, we don't need to run the subsequent tests.
if (n_groups == 1){
write("Only 1 grouping present, skipping differential expression tests.", stderr())
quit(status = 0, save = "no")
} 

#Check if group IDS are numbers and if so, modify them
if (is.numeric(uniq_groups)){
write("Group IDs are numeric, prepending 'group' to group ID for compatibility with edgeR tools", stderr())
uniq_groups <- sub('^', 'group', uniq_groups)
}

#Subset the data object to get rid of samples with a 'NULL' group
edge_mat <- edge_mat[,which(edge_mat$samples$group != 'NULL')]

#Redefine these variables
uniq_groups <- as.vector(unique(edge_mat$samples$group))
#n_groups = length(uniq_groups)

# Generate the design matrix for GLM fitting and estimate common and tag-wise dispersion in one go.
design <- model.matrix(~0+group, data = edge_mat$samples)
colnames(design) <- uniq_groups
edge_mat <- estimateDisp(edge_mat, design = design)

# Fit the per-feature negative binomial GLM.
fit <- glmQLFit(edge_mat, design)

# Generate the matrix of pairwise combinations
combos <- combinations(n = n_true_groups, r = 2, v = uniq_groups, repeats.allowed = F)

# Now we loop over all unique (non-self) comparisons, performing the DE test via quasi-likelihood F-test (edgeR recommended), and write the significant (at 0.05) DE genes to file.  
for (row in 1:nrow(combos)){
	comp <- paste(combos[row,][2],"-",combos[row][1], sep = "")
	comp_var <- makeContrasts(comp, levels = design)
	qlf <- glmQLFTest(fit, contrast  = comp_var)
	tags = topTags(qlf, n = nrow(qlf$genes))
	de_file <- paste(out_dir, "/DEGs/DE_", comp, "_list.txt", sep = "") 
	write.table(tags$table[which(tags$table$FDR < 0.05),], file = de_file, sep = '\t', quote = FALSE, row.names = FALSE)
}
