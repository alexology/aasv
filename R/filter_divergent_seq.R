#' Filter divergent sequences from species-level alignments
#'
#' @description
#' Removes sequences that are too genetically divergent from every other
#' sequence in the same species-level alignment.
#'
#' @param paths Character vector of file paths to species-level alignments,
#'   typically all files produced by \code{\link{align_species_seq}} after
#'   \code{\link{improve_alignment}} has cleaned up gap artefacts.
#' @param max_distance Hard distance threshold (proportion of differing
#'   sites, 0-1). A sequence is removed when its distance to its nearest
#'   neighbour within the file exceeds this value. Default \code{0.10}
#'   (10\%).
#' @param min_seqs Minimum number of sequences required in a file before the
#'   check is applied. Files with fewer sequences are left untouched, since
#'   an outlier cannot be reliably identified from too few individuals.
#'   Default \code{4}.
#' @param print_messages Print a progress message for each processed file
#'   reporting how many sequences were retained. Default \code{TRUE}.
#' @param log_dir Optional path to a directory where a removal log will be
#'   written. A file \code{removed_divergent.tsv} is created (or appended
#'   to) with one row per removed sequence recording the source file, taxon
#'   name, sequence ID, its nearest-neighbour distance, and the before/after
#'   counts. If \code{NULL} (the default), no log is written.
#' @param ... Further arguments passed to \code{DistanceMatrix} of the
#'   \code{DECIPHER} package (e.g. \code{processors} to control
#'   parallelism, or \code{verbose = FALSE} to suppress its console output).
#'
#' @details
#' This function is meant to run \emph{after} \code{\link{improve_alignment}}:
#' gap artefacts should already be resolved, so any remaining divergence is
#' due to the sequence content itself (misidentification, contamination,
#' NUMTs) rather than an alignment problem.
#'
#' Distance is computed once per file with a single call to
#' \code{DECIPHER::DistanceMatrix()} (terminal gaps excluded), which is the
#' expensive step. Sequences are then flagged by their \strong{nearest-
#' neighbour distance} (the smallest pairwise distance to any other sequence
#' in the file) rather than their mean/median distance to the group: this
#' catches true singleton outliers while tolerating legitimate multi-modal
#' within-species structure, since a sequence belonging to either subclade
#' still has a close relative in the other. Flagged sequences are all
#' dropped in a single pass - there is no iterative re-computation of the
#' distance matrix as sequences are removed.
#'
#' If \emph{every} sequence in a file is flagged, nothing is removed and a
#' warning is issued instead: this usually signals a threshold or taxonomy
#' problem rather than independent bad sequences, and removing everything
#' would leave no reference for that species.
#'
#' After removal, gap-only columns are trimmed from the remaining
#' sequences; the alignment itself is not recomputed, since
#' \code{\link{improve_alignment}} is expected to have already produced a
#' gap-clean alignment upstream.
#'
#' @return The input \code{paths} vector, returned invisibly.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet width AAStringSet
#' @importFrom DECIPHER DistanceMatrix
#' @importFrom stats setNames
#' @importFrom utils write.table
#'
#' @export

filter_divergent_seq <- function(paths, max_distance = 0.10, min_seqs = 4,
                                 print_messages = TRUE, log_dir = NULL, ...) {

  if (!is.null(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file   <- file.path(log_dir, "removed_divergent.tsv")
    log_header <- !file.exists(log_file) || file.size(log_file) == 0L
  }

  for (ali_file in paths) {

    ali        <- Biostrings::readAAStringSet(ali_file)
    n_start    <- length(ali)
    taxon_name <- gsub("\\.fasta$", "", basename(ali_file))

    if (n_start < max(2L, min_seqs)) next()

    dist_mat      <- DECIPHER::DistanceMatrix(ali, includeTerminalGaps = FALSE, ...)
    diag(dist_mat) <- NA

    nn_dist <- apply(dist_mat, 1L, min, na.rm = TRUE)
    bad     <- nn_dist > max_distance

    if (!any(bad)) {
      if (print_messages) message(taxon_name, ": no divergent sequences")
      next()
    }

    if (all(bad)) {
      warning(taxon_name, ": every sequence exceeds max_distance (",
              max_distance, ") - skipping, check taxonomy or threshold",
              call. = FALSE)
      next()
    }

    ali_out <- ali[names(ali)[!bad]]

    mat       <- do.call(rbind, strsplit(as.character(ali_out), "", fixed = TRUE))
    keep_cols <- colSums(mat != "-") > 0L
    if (!all(keep_cols)) {
      seqs    <- apply(mat[, keep_cols, drop = FALSE], 1L, paste, collapse = "")
      ali_out <- Biostrings::AAStringSet(setNames(seqs, names(ali_out)))
    }

    Biostrings::writeXStringSet(ali_out,
                                filepath = ali_file,
                                width    = max(Biostrings::width(ali_out)) + 1L,
                                format   = "fasta")

    if (!is.null(log_dir)) {
      removed <- names(ali)[bad]
      log_df  <- data.frame(
        file       = basename(ali_file),
        taxon      = taxon_name,
        sequence   = removed,
        distance   = round(nn_dist[bad], 4),
        n_start    = n_start,
        n_retained = length(ali_out),
        stringsAsFactors = FALSE
      )
      write.table(log_df,
                  file      = log_file,
                  sep       = "\t",
                  row.names = FALSE,
                  col.names = log_header,
                  append    = !log_header,
                  quote     = FALSE)
      log_header <- FALSE
    }

    if (print_messages) {
      message(taxon_name, ": ", length(ali_out), "/", n_start,
              " sequences retained (", sum(bad), " divergent removed)")
    }
  }

  invisible(paths)
}
