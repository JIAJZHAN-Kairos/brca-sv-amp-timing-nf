#!/usr/bin/env Rscript
# Run AmplificationTimeR (Wedge-lab, pinned commit) on every event in events.tsv.
# Usage: run_amptimer.R <input_dir> <out_dir> <amptimer_commit>
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]; COMMIT <- args[3]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lib <- file.path(getwd(), "Rlib"); dir.create(lib, showWarnings = FALSE); .libPaths(c(lib, .libPaths()))
options(repos = BiocManager::repositories(), Ncpus = 4)
need <- c("VariantAnnotation", "SomaticSignatures", "BSgenome.Hsapiens.UCSC.hg38", "remotes")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) BiocManager::install(miss, lib = lib, update = FALSE, ask = FALSE)
if (!requireNamespace("AmplificationTimeR", quietly = TRUE))
  remotes::install_github(paste0("Wedge-lab/AmplificationTimeR@", COMMIT), lib = lib, upgrade = "never", dependencies = FALSE)
suppressPackageStartupMessages({ library(AmplificationTimeR); library(BSgenome.Hsapiens.UCSC.hg38) })
ev <- read.delim(file.path(IN, "events.tsv"), stringsAsFactors = FALSE, colClasses = c(seg_chr = "character"))
rd <- function(s, what) read.delim(file.path(IN, paste0(s, ".", what, ".tsv")), stringsAsFactors = FALSE,
                                   colClasses = c(chr = "character"))
res <- list()
for (i in seq_len(nrow(ev))) {
  e <- ev[i, ]
  warn <- character()
  r <- tryCatch(withCallingHandlers({
    if (is.na(e$seg_chr)) stop("no single MutationTimeR segment at gene midpoint")
    out <- time_amplification(cn_data = rd(e$sample, "cn"), multiplicity_data = rd(e$sample, "mult"),
                              mutation_data = rd(e$sample, "muts"), muts_type = "All", sample_id = e$sample,
                              amplification_chrom = e$seg_chr, amplification_start = as.numeric(e$seg_start),
                              amplification_stop = as.numeric(e$seg_end), is_WGD = as.logical(e$is_wgd), genome = "hg38")
    out <- as.data.frame(out)
    cbind(event_id = e$event_id, gene = e$gene, set = e$set, status = "ok", error = "", out)
  }, warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = function(err) data.frame(event_id = e$event_id, gene = e$gene, set = e$set, sample = e$sample,
                                   status = "not_timed", error = conditionMessage(err)))
  r$warnings <- paste(unique(warn), collapse = "; ")
  res[[i]] <- r
  cat(e$event_id, r$status, r$error, "\n")
}
keys <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(x) { for (k in setdiff(keys, names(x))) x[[k]] <- NA; x[keys] }))
write.table(res, file.path(OUT, "amplificationtimer_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(c(capture.output(sessionInfo()), paste("AmplificationTimeR commit:", COMMIT)), file.path(OUT, "sessionInfo.txt"))
print(table(res$set, res$status))
