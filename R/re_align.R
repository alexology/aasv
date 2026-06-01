#' Re-align sequences with different alignment parameters
#'
#' @description
#' Re-aligns a set of previously aligned sequences using new parameters passed
#' to \code{AlignSeqs}.
#'
#' @param x Character vector of absolute paths to alignment FASTA files.
#' @param raw_sequences If \code{TRUE}, read the corresponding raw (unaligned)
#'   sequences from the \code{raw_sequences} folder instead of the alignment
#'   file, then save the result back to the species-level alignment path. This
#'   is intended for species-level paths only; genus- and family-level paths
#'   do not map 1-to-1 with individual raw-sequence files.
#' @param verbose Print a progress message for each processed file. Default
#'   \code{TRUE}.
#' @param ... Further arguments passed to \code{AlignSeqs} of the
#'   \code{DECIPHER} package.
#'
#' @details
#' When \code{raw_sequences = FALSE} (the default), existing gap characters are
#' stripped from the alignment before \code{AlignSeqs} is called, so the
#' aligner starts from the raw sequences rather than inheriting the previous
#' gap pattern. Each file is overwritten in place.
#'
#' When \code{raw_sequences = TRUE}, the function derives the raw-sequence path
#' by replacing \code{alignments} with \code{raw_sequences} in each input path
#' and redirecting any genus- or family-level segment to \code{species}. The
#' re-aligned result is saved to the corresponding species-level alignment path.
#'
#' @return The input vector \code{x}, returned invisibly.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet width AAStringSet
#' @importFrom DECIPHER AlignSeqs
#'
#' @export

re_align <- function(x, raw_sequences = FALSE, verbose = TRUE, ...) {

  strip_gaps <- function(ali) {
    Biostrings::AAStringSet(gsub("-", "", as.character(ali)))
  }

  for (path in x) {

    if (raw_sequences) {
      to_read <- gsub("alignments",   "raw_sequences", path)
      to_read <- gsub("family|genus", "species",       to_read)
      to_save <- gsub("family|genus", "species",       path)
    } else {
      to_read <- path
      to_save <- path
    }

    ali <- Biostrings::readAAStringSet(to_read)

    # strip alignment gaps so AlignSeqs treats every sequence as unaligned
    if (!raw_sequences) {
      ali <- strip_gaps(ali)
    }

    new_ali <- DECIPHER::AlignSeqs(ali, ...)

    Biostrings::writeXStringSet(new_ali,
                                filepath = to_save,
                                width    = max(Biostrings::width(new_ali)) + 1,
                                format   = "fasta")

    if (verbose) {
      seq_name <- gsub("\\.fasta$", "", basename(to_read))
      message(gsub("_", " ", seq_name), ": done")
    }
  }

  invisible(x)
}
