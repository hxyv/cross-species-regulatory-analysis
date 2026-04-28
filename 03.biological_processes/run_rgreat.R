args <- commandArgs(trailingOnly = TRUE)
base <- args[1] # project base directory (for outputs and mouse files)
data_root <- args[2] # path to human .bed.gz files

# use a custom R library path via environment variable
user_lib <- Sys.getenv("R_USER_LIB", unset = NA)
if (!is.na(user_lib) && nzchar(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

# load GREAT analysis 
# load tools for genomic ranges, GRanges
library(rGREAT)
library(GenomicRanges)

# output directory for results
outdir <- file.path(base, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# read .bed files
read_bed_simple <- function(file) {
  df <- read.table(file, header = FALSE, sep = "\t")[, 1:3]
  colnames(df) <- c("chr", "start", "end")

  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start + 1, end = df$end)
  )
}

# human using the original hg38 corrdinates 
human_specific <- read_bed_simple(file.path(data_root, "human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed.gz"))
human_shared <- read_bed_simple(file.path(data_root, "human_adrenal_idr_optimal.shared.original_human_coordinates.bed.gz"))
all_human <- c(human_specific, human_shared)

# mouse using mm10 coordinates 
mouse_specific <- read_bed_simple(file.path(base, "mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed"))
shared_mouse <- read_bed_simple(file.path(base, "mouse_adrenal_idr_optimal.shared_with_human_mapped.bed"))
all_mouse <- c(mouse_specific, shared_mouse)

# Run RGreat 

# human using (hg38)
job_human_specific <- submitGreatJob(human_specific, species = "hg38")
job_human_shared <- submitGreatJob(human_shared, species = "hg38")
job_all_human <- submitGreatJob(all_human,species = "hg38")

# mouse using (mm10)
job_mouse_specific <- submitGreatJob(mouse_specific, species = "mm10")
job_shared_mouse <- submitGreatJob(shared_mouse, species = "mm10")
job_all_mouse <- submitGreatJob(all_mouse,species = "mm10")

# Get results
# Human Results 
res_human_specific <- getEnrichmentTables(job_human_specific)
res_human_shared <- getEnrichmentTables(job_human_shared)
res_all_human <- getEnrichmentTables(job_all_human)
# Mouse Results
res_mouse_specific <- getEnrichmentTables(job_mouse_specific)
res_shared_mouse <- getEnrichmentTables(job_shared_mouse)
res_all_mouse <- getEnrichmentTables(job_all_mouse)

# Exrtract the GO BP
get_bp <- function(res) res[["GO Biological Process"]]

# Filter and rank results using fold enrichment 
summarize_go <- function(res, label) {
  bp <- get_bp(res)
  bp_sig <- bp[bp$Binom_Adjp_BH < 0.05, ]
  bp_sig <- bp_sig[order(bp_sig$Binom_Fold_Enrichment, decreasing = TRUE), ]
  bp_sig$set <- rep(label, nrow(bp_sig))

  return(bp_sig)
}

top_human_specific <- summarize_go(res_human_specific, "human_specific")
top_human_shared <- summarize_go(res_human_shared, "human_shared")
top_all_human <- summarize_go(res_all_human, "all_human")

top_mouse_specific <- summarize_go(res_mouse_specific, "mouse_specific")
top_shared_mouse <- summarize_go(res_shared_mouse, "shared_mouse")
top_all_mouse <- summarize_go(res_all_mouse, "all_mouse")

# Save filtered data
write.csv(top_human_specific, file.path(outdir, "human_specific_BP_filtered.csv"), row.names = FALSE)
write.csv(top_human_shared, file.path(outdir, "human_shared_BP_filtered.csv"), row.names = FALSE)
write.csv(top_all_human, file.path(outdir, "all_human_BP_filtered.csv"), row.names = FALSE)

write.csv(top_mouse_specific, file.path(outdir, "mouse_specific_BP_filtered.csv"), row.names = FALSE)
write.csv(top_shared_mouse, file.path(outdir, "shared_mouse_BP_filtered.csv"), row.names = FALSE)
write.csv(top_all_mouse, file.path(outdir, "all_mouse_BP_filtered.csv"), row.names = FALSE)

# Completed
message("Done! Correct species-specific GREAT analysis complete.")
