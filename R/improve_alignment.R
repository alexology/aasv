#' Improve problematic alignments
#'
#' @description
#' Re-aligns sequences flagged by \code{\link{check_alignment}} by iteratively
#' removing gap-causing and self-gapped sequences until the alignment is clean.
#'
#' @param paths Character vector of file paths to problematic alignments,
#'   as returned by \code{\link{check_alignment}}.
#' @param print_messages Print a progress message for each processed file
#'   reporting how many sequences were retained. Default \code{TRUE}.
#' @param ... Further arguments passed to \code{AlignSeqs} and
#'   \code{AlignProfiles} of the \code{DECIPHER} package. Use
#'   \code{verbose = FALSE} here to suppress DECIPHER's console output.
#'
#' @details
#' COI is a protein-coding gene: internal gaps in a species-level alignment
#' are biologically impossible and indicate a bad sequence (chimera,
#' frameshift, low-quality read).
#'
#' \strong{Deduplication.} Before aligning, sequences with identical amino
#' acid content are collapsed to a single representative. All alignment and
#' gap-removal steps operate on the unique representatives only. Duplicates
#' are re-expanded in place once the final alignment is obtained.
#'
#' \strong{Block alignment.} Sequences are grouped by the species prefix in
#' their name (text before \code{"__"}). Within each block \code{AlignSeqs}
#' is used; blocks are then merged with \code{AlignProfiles}. At species
#' level all sequences share the same prefix, so this degenerates to a
#' single \code{AlignSeqs} call. At genus or family level only the affected
#' block is re-aligned after each removal, avoiding a full N-sequence
#' re-alignment from scratch.
#'
#' \strong{Step 1 - gap causers (priority).} For every alignment column, only
#' sequences that \emph{span} that column (non-gap characters on both sides)
#' are considered. If the majority (> 50\%) of spanning sequences carry an
#' internal gap at that column, the column is flagged; any sequence with a
#' residue in a flagged column is the outlier forcing those gaps and is
#' removed first. The spanning-based check avoids penalising sequences that
#' are merely longer or shorter than the majority (terminal gaps, not
#' internal). Step 1 is repeated until no gap-causing sequences remain.
#'
#' \strong{Step 2 - self-gapped sequences.} Once no sequence causes gaps in
#' others, sequences that carry internal gaps within their own span are
#' removed and the alignment is rebuilt.
#'
#' After each removal the function first tries dropping gap-only columns from
#' the remaining aligned sequences (fast path); a full re-alignment is
#' triggered only if problems persist.
#'
#' When all remaining candidates are flagged, the function searches for an
#' alignment file of another species from the same genus in the same
#' directory and uses it as an external reference: each candidate is aligned
#' pairwise against one reference sequence and the one with fewest gap
#' characters is kept. If no sibling species file is found, one is kept at
#' random. Each fixed file overwrites the original so that downstream calls
#' to \code{\link{align_genus_seq}} and \code{\link{align_family_seq}} pick
#' up the corrected alignments.
#'
#' @return The input \code{paths} vector, returned invisibly.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet width AAStringSet
#' @importFrom DECIPHER AlignSeqs AlignProfiles
#' @importFrom stats setNames
#'
#' @export

improve_alignment <- function(paths, print_messages = TRUE, ...) {

  # --- helpers ---------------------------------------------------------------

  strip_gaps <- function(x) {
    Biostrings::AAStringSet(setNames(gsub("-", "", as.character(x)), names(x)))
  }

  build_mat <- function(x) {
    do.call(rbind, strsplit(as.character(x), "", fixed = TRUE))
  }

  trim_cols <- function(x) {
    mat  <- build_mat(x)
    keep <- colSums(mat != "-") > 0L
    if (all(keep)) return(x)
    seqs <- apply(mat[, keep, drop = FALSE], 1L, paste, collapse = "")
    Biostrings::AAStringSet(setNames(seqs, names(x)))
  }

  # Collapse identical sequences to one representative each.
  # Returns: reps   - unique AAStringSet (one per distinct AA string)
  #          rep_nm - named char vector: original name → representative name
  dedup_seqs <- function(x) {
    seqs   <- as.character(x)
    first  <- !duplicated(seqs)
    reps   <- x[first]
    rep_nm <- setNames(names(reps)[match(seqs, seqs[first])], names(x))
    list(reps = reps, rep_nm = rep_nm)
  }

  # Broadcast an aligned set of representatives back to all original sequences.
  expand_to_all <- function(ali_reps, rep_nm, all_names) {
    chars <- setNames(as.character(ali_reps), names(ali_reps))
    Biostrings::AAStringSet(setNames(chars[rep_nm], all_names))
  }

  # Align within species blocks (prefix before "__"), merge with AlignProfiles.
  # Single-block case (species level) degenerates to plain AlignSeqs.
  align_blocks <- function(x, ...) {
    groups   <- sub("__.*$", "", names(x))
    grp_list <- split(x, groups)
    aligned  <- lapply(grp_list, function(blk) {
      if (length(blk) == 1L) blk else DECIPHER::AlignSeqs(blk, ...)
    })
    if (length(aligned) == 1L) return(aligned[[1L]])
    # Merge largest block first for numerical stability
    ord    <- order(lengths(aligned), decreasing = TRUE)
    result <- aligned[[ord[1L]]]
    for (i in ord[-1L]) {
      result <- DECIPHER::AlignProfiles(result, aligned[[i]])
    }
    result
  }

  causes_internal_gaps <- function(mat) {
    nc      <- ncol(mat)
    non_gap <- mat != "-"
    seq_starts <- apply(non_gap, 1L, function(r) {
      w <- which(r)
      if (length(w)) w[1L] else nc + 1L
    })
    seq_ends <- apply(non_gap, 1L, function(r) {
      w <- which(r)
      if (length(w)) w[length(w)] else 0L
    })
    cols     <- seq_len(nc)
    spanning <- outer(seq_starts, cols, "<=") &
                outer(seq_ends,   cols, ">=")
    n_spanning        <- colSums(spanning)
    n_internal_gap    <- colSums(spanning & !non_gap)
    internal_gap_frac <- ifelse(n_spanning > 1L, n_internal_gap / n_spanning, 0)
    apply(non_gap, 1L, function(r) any(r & (internal_gap_frac > 0.5)))
  }

  has_self_gap <- function(mat) {
    apply(mat, 1L, function(row) {
      ng <- which(row != "-")
      length(ng) >= 2L && any(row[ng[1L]:ng[length(ng)]] == "-")
    })
  }

  try_trim <- function(x) {
    ali <- trim_cols(x)
    mat <- build_mat(ali)
    ok  <- !any(causes_internal_gaps(mat)) && !any(has_self_gap(mat))
    list(ali = ali, clean = ok)
  }

  pick_by_reference <- function(candidates, ali_file) {
    genus    <- sub("_.*$", "", gsub("\\.fasta$", "", basename(ali_file)))
    siblings <- setdiff(
      list.files(dirname(ali_file),
                 pattern   = paste0("^", genus, "_.*\\.fasta$"),
                 full.names = TRUE),
      ali_file
    )
    if (length(siblings) == 0L) return(sample.int(length(candidates), 1L))
    ref <- strip_gaps(Biostrings::readAAStringSet(siblings[[1L]])[1L])
    gap_counts <- vapply(seq_along(candidates), function(i) {
      ali <- DECIPHER::AlignSeqs(c(candidates[i], ref), ...)
      sum(build_mat(ali)[1L, ] == "-")
    }, integer(1L))
    which.min(gap_counts)
  }

  # Remove bad representatives, update dd and all_names, try fast trim.
  # Returns list(ali_reps, dd, all_names, done).
  remove_bad <- function(bad, ali_reps, dd, all_names, ali_file) {
    keep_nms <- names(ali_reps)[!bad]
    if (length(keep_nms) == 0L)
      keep_nms <- names(dd$reps)[pick_by_reference(dd$reps, ali_file)]

    keep_orig <- dd$rep_nm %in% keep_nms
    all_names <- all_names[keep_orig]
    dd$rep_nm <- dd$rep_nm[keep_orig]

    if (length(keep_nms) == 1L) {
      ali_reps <- strip_gaps(ali_reps[keep_nms])
      dd$reps  <- ali_reps
      return(list(ali_reps = ali_reps, dd = dd,
                  all_names = all_names, done = TRUE))
    }

    res <- try_trim(ali_reps[keep_nms])
    if (res$clean) {
      return(list(ali_reps = res$ali, dd = dd,
                  all_names = all_names, done = TRUE))
    }
    dd$reps <- strip_gaps(ali_reps[keep_nms])
    list(ali_reps = NULL, dd = dd, all_names = all_names, done = FALSE)
  }

  # --- main loop -------------------------------------------------------------

  for (ali_file in paths) {

    ali        <- Biostrings::readAAStringSet(ali_file)
    taxon_name <- gsub("\\.fasta$", "", basename(ali_file))
    n_start    <- length(ali)

    if (n_start == 0L) next()

    seqs      <- strip_gaps(ali)
    dd        <- dedup_seqs(seqs)
    all_names <- names(seqs)
    t_start   <- proc.time()[["elapsed"]]
    t_last    <- t_start
    iter      <- 0L

    repeat {

      iter <- iter + 1L

      if (length(dd$reps) <= 1L) {
        ali_reps <- dd$reps
        break
      }

      now <- proc.time()[["elapsed"]]
      if (print_messages && (now - t_last) >= 30) {
        message(taxon_name, ": iteration ", iter, " - ",
                length(dd$reps), " unique sequences remaining (",
                round(now - t_start), "s)")
        t_last <- now
      }

      ali_reps <- align_blocks(dd$reps, ...)
      mat      <- build_mat(ali_reps)

      # Step 1 (priority): remove gap causers
      bad <- causes_internal_gaps(mat)
      if (any(bad)) {
        res <- remove_bad(bad, ali_reps, dd, all_names, ali_file)
        dd        <- res$dd
        all_names <- res$all_names
        if (res$done) {
          ali_reps <- res$ali_reps
          break
        }
        next
      }

      # Step 2: remove self-gapped sequences
      bad <- has_self_gap(mat)
      if (!any(bad)) break

      res <- remove_bad(bad, ali_reps, dd, all_names, ali_file)
      dd        <- res$dd
      all_names <- res$all_names
      if (res$done) {
        ali_reps <- res$ali_reps
        break
      }
    }

    ali_out <- expand_to_all(ali_reps, dd$rep_nm, all_names)

    Biostrings::writeXStringSet(ali_out,
                                filepath = ali_file,
                                width    = max(Biostrings::width(ali_out)) + 1L,
                                format   = "fasta")

    if (print_messages) {
      message(taxon_name, ": ", length(ali_out), "/", n_start,
              " sequences retained")
    }
  }

  invisible(paths)
}
