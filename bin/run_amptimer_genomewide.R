#!/usr/bin/env Rscript
# Genome-wide AmplificationTimeR: time every gained CN segment of the given samples.
# Usage: run_amptimer_genomewide.R <input_dir> <out_tsv> <amptimer_commit> <sample> [<sample> ...]
#   input_dir holds events.tsv and <s>.{cn,mult,muts}.tsv (prepare_amptimer_genomewide.R).
# time_amplification() annotates the trinucleotide context of every mutation it is given on each call,
# but only uses mutations inside the region. To keep runtime linear in the number of events, each call
# gets the mutations inside its region; outputs are identical to passing the whole table. A region
# with < 3 multiplicity rows gets the whole table, so the tool fails with its own message.
# Set AMPT_FULL_MUTS=1 to pass the whole table on every call (used once to check the two agree).
args <- commandArgs(TRUE); IN <- args[1]; OUTF <- args[2]; COMMIT <- args[3]; SAMPLES <- args[-(1:3)]
lib <- file.path(getwd(), "Rlib"); dir.create(lib, showWarnings = FALSE); .libPaths(c(lib, .libPaths()))
if (!requireNamespace("AmplificationTimeR", quietly = TRUE)) {
  options(repos = BiocManager::repositories(), Ncpus = 4)
  need <- c("VariantAnnotation", "SomaticSignatures", "BSgenome.Hsapiens.UCSC.hg38", "remotes")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) BiocManager::install(miss, lib = lib, update = FALSE, ask = FALSE)
  remotes::install_github(paste0("Wedge-lab/AmplificationTimeR@", COMMIT), lib = lib, upgrade = "never", dependencies = FALSE)
}
suppressPackageStartupMessages({ library(AmplificationTimeR); library(BSgenome.Hsapiens.UCSC.hg38) })
ev <- read.delim(file.path(IN, "events.tsv"), stringsAsFactors = FALSE, colClasses = c(chr = "character"))
rd <- function(s, what) read.delim(file.path(IN, paste0(s, ".", what, ".tsv")), stringsAsFactors = FALSE,
                                   colClasses = c(chr = "character"))
res <- list()
for (s in SAMPLES) {
  cn <- rd(s, "cn"); mult <- rd(s, "mult"); muts <- rd(s, "muts")
  e_s <- ev[ev$sample_id == s, ]
  t0 <- Sys.time()
  for (i in seq_len(nrow(e_s))) {
    e <- e_s[i, ]
    inreg <- mult$chr == e$chr & mult$end >= e$startpos & mult$end <= e$endpos
    mu <- if (sum(inreg) >= 3 && Sys.getenv("AMPT_FULL_MUTS") != "1") muts[muts$chr == e$chr & muts$end >= e$startpos & muts$end <= e$endpos, ] else muts
    warn <- character()
    r <- tryCatch(withCallingHandlers({
      out <- time_amplification(cn_data = cn, multiplicity_data = mult, mutation_data = mu, muts_type = "All",
                                sample_id = s, amplification_chrom = e$chr, amplification_start = as.numeric(e$startpos),
                                amplification_stop = as.numeric(e$endpos), is_WGD = as.logical(e$is_wgd), genome = "hg38")
      cbind(event_id = e$event_id, status = "ok", error = "", as.data.frame(out))
    }, warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(err) data.frame(event_id = e$event_id, status = "not_timed", error = conditionMessage(err), sample = s))
    r$warnings <- paste(unique(warn), collapse = "; ")
    res[[length(res) + 1]] <- r
  }
  cat(s, nrow(e_s), "events", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
keys <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(x) { for (k in setdiff(keys, names(x))) x[[k]] <- NA; x[keys] }))
write.table(res, OUTF, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(table(res$status))
