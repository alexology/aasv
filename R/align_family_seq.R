#' Align family-level sequences
#'
#' @description
#' Merges genus-level alignments into a single family-level alignment by
#' progressively aligning profiles with \code{AlignProfiles} from the
#' \code{DECIPHER} package.
#'
#' @param taxon Character vector of family names to process. If \code{NULL}
#'   (default), all families found under \code{alignment_dir/alignments/} that
#'   have genus-level files are processed.
#' @param alignment_dir Path to the folder that contains the \code{alignments}
#'   directory produced by \code{\link{align_genus_seq}}.
#' @param overwrite Overwrite existing family-level alignments.
#' @param ... Further arguments passed to \code{AlignProfiles} of the
#'   \code{DECIPHER} package.
#'
#' @details
#' Genus-level alignments are merged progressively: the first profile is used
#' as the anchor and each subsequent profile is added to it with
#' \code{AlignProfiles}. Before merging, files are sorted in descending order
#' of sequence count so that the largest (most informative) profile serves as
#' the anchor. This reduces gap propagation that can occur when a small profile
#' is merged first.
#'
#' @return Called for its side effects (writing FASTA files to disk).
#'   Returns \code{NULL} invisibly.
#'
#' @importFrom Biostrings readAAStringSet writeXStringSet width
#' @importFrom DECIPHER AlignProfiles
#'
#' @examples
#' \dontrun{
#' align_family_seq(alignment_dir = "path/to/results")
#' }
#'
#' @export

align_family_seq <- function(taxon = NULL,
                             alignment_dir = NULL,
                             overwrite = FALSE,
                             ...) {

  if (is.null(alignment_dir)) {
    alignment_dir <- getwd()
  }

  if (is.null(taxon)) {
    all_dirs  <- list.dirs(file.path(alignment_dir, "alignments"), recursive = FALSE)
    has_files <- sapply(file.path(all_dirs, "genus"), function(d) length(list.files(d)) > 0)
    taxon     <- basename(all_dirs[has_files])
  }

  for (fam in taxon) {

    family_dir <- file.path(alignment_dir, "alignments", fam, "family")

    if (dir.exists(family_dir) && !overwrite) {
      message(fam, ": already present")
      next()
    }

    genus_files <- list.files(file.path(alignment_dir, "alignments", fam, "genus"),
                              full.names = TRUE,
                              pattern = "\\.fasta$")

    if (length(genus_files) == 0) {
      message(fam, ": no genus files found, skipping")
      next()
    }

    # sort files largest-first so the most populated profile anchors the
    # progressive merge, reducing gap propagation from small profiles
    seq_counts  <- sapply(genus_files, function(f) length(Biostrings::readAAStringSet(f)))
    genus_files <- genus_files[order(seq_counts, decreasing = TRUE)]

    dir.create(family_dir, showWarnings = FALSE)

    family_path <- file.path(family_dir, paste0(fam, ".fasta"))

    query_fasta <- Biostrings::readAAStringSet(genus_files[1])

    for (i in seq_along(genus_files)[-1]) {
      query_temp  <- Biostrings::readAAStringSet(genus_files[i])
      query_fasta <- DECIPHER::AlignProfiles(query_fasta, query_temp, ...)
    }

    Biostrings::writeXStringSet(query_fasta,
                                filepath = family_path,
                                width = max(Biostrings::width(query_fasta)) + 1,
                                format = "fasta")

    message(fam, ": done")
  }
}
