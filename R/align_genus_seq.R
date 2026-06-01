#' Align genus-level sequences
#'
#' @description
#' Merges species-level alignments into a single genus-level alignment by
#' progressively aligning profiles with \code{AlignProfiles} from the
#' \code{DECIPHER} package.
#'
#' @param taxon Character vector of family names to process. If \code{NULL}
#'   (default), all families found under \code{alignment_dir/alignments/} that
#'   have species-level files are processed.
#' @param alignment_dir Path to the folder that contains the \code{alignments}
#'   directory produced by \code{\link{align_species_seq}}.
#' @param overwrite Overwrite existing genus-level alignments.
#' @param ... Further arguments passed to \code{AlignProfiles} of the
#'   \code{DECIPHER} package.
#'
#' @details
#' Species-level alignments are merged progressively: the first profile is used
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
#' @export
#'
#' @examples
#' \dontrun{
#' align_genus_seq(alignment_dir = "path/to/results")
#' }

align_genus_seq <- function(taxon = NULL,
                            alignment_dir = NULL,
                            overwrite = FALSE,
                            ...) {

  if (is.null(alignment_dir)) {
    alignment_dir <- getwd()
  }

  if (is.null(taxon)) {
    all_dirs  <- list.dirs(file.path(alignment_dir, "alignments"), recursive = FALSE)
    has_files <- sapply(file.path(all_dirs, "species"), function(d) length(list.files(d)) > 0)
    taxon     <- basename(all_dirs[has_files])
  }

  for (fam in taxon) {

    genus_dir <- file.path(alignment_dir, "alignments", fam, "genus")

    if (dir.exists(genus_dir) && !overwrite) {
      message(fam, ": already present")
      next()
    }

    species_path  <- file.path(alignment_dir, "alignments", fam, "species")
    species_files <- list.files(species_path, pattern = "\\.fasta$")

    if (length(species_files) == 0) {
      message(fam, ": no species files found, skipping")
      next()
    }

    species_paths <- file.path(species_path, species_files)

    # extract genus as the first underscore-delimited token of each filename
    genus_unique <- unique(sapply(strsplit(species_files, "_"), `[[`, 1))

    dir.create(genus_dir, showWarnings = FALSE)

    for (genus in genus_unique) {

      query_files <- species_paths[grepl(paste0("^", genus, "_"), species_files)]

      # sort files largest-first so the most populated profile anchors the
      # progressive merge, reducing gap propagation from small profiles
      seq_counts  <- sapply(query_files, function(f) length(Biostrings::readAAStringSet(f)))
      query_files <- query_files[order(seq_counts, decreasing = TRUE)]

      genus_path  <- file.path(genus_dir, paste0(genus, ".fasta"))

      query_fasta <- Biostrings::readAAStringSet(query_files[1])

      for (i in seq_along(query_files)[-1]) {
        query_temp  <- Biostrings::readAAStringSet(query_files[i])
        query_fasta <- DECIPHER::AlignProfiles(query_fasta, query_temp, ...)
      }

      Biostrings::writeXStringSet(query_fasta,
                                  filepath = genus_path,
                                  width = max(Biostrings::width(query_fasta)) + 1,
                                  format = "fasta")

      message(genus, ": done")
    }

    message(fam, ": done")
  }
}
