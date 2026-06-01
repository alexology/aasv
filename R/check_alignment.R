#' Alignment problems
#'
#' @description
#' Check for problems in sequence alignments.
#'
#' @param taxon A family name for which the alignment exists.
#' @param alignment_dir Path to the top-level results folder produced by
#'   \code{\link{align_species_seq}} (the folder that contains the
#'   \code{alignments} subdirectory). Defaults to the working directory.
#' @param type Subfolder to search. Use \code{"alignments"} (default) for
#'   aligned amino acid sequences, or \code{"raw_sequences"} for the
#'   untranslated raw sequences.
#' @param tax_lev Taxonomic level of the alignments. Must be \code{"family"},
#'   \code{"genus"} or \code{"species"}.
#'
#' @details
#' For each alignment file, the function checks whether any sequence contains
#' internal gaps (i.e. gaps that are not at the leading or trailing edge of
#' the sequence). Terminal gaps are expected in a multiple sequence alignment
#' and are considered normal; internal gaps in closely related sequences may
#' indicate an alignment artefact.
#'
#' Each file is reported as either \emph{"Alignment should be fine"} (only
#' terminal gaps found) or \emph{"Potential problems in the alignment"}
#' (internal gaps found). The function returns the paths of problematic files
#' so they can be passed directly to \code{\link{improve_alignment}}.
#'
#' @return A character vector of paths to alignments with potential problems.
#'   An empty vector is returned when no problems are found.
#'
#' @importFrom Biostrings readAAStringSet
#'
#' @export

check_alignment <- function(taxon,
                            alignment_dir = NULL,
                            type    = "alignments",
                            tax_lev = "species") {

  if (is.null(alignment_dir)) {
    alignment_dir <- getwd()
  }

  all_ali <- list.files(file.path(alignment_dir, type, taxon, tax_lev),
                        full.names = TRUE,
                        pattern    = "\\.fasta$")

  if (length(all_ali) == 0) {
    stop("No files in the specified directory.")
  }

  # returns TRUE if a single alignment row contains an internal gap,
  # i.e. a gap that lies between the first and last non-gap character
  has_internal_gap <- function(row) {
    non_gap <- which(row != "-")
    if (length(non_gap) < 2) return(FALSE)
    any(row[non_gap[1]:non_gap[length(non_gap)]] == "-")
  }

  path_store <- character(0)

  for (ali_file in all_ali) {

    ali <- Biostrings::readAAStringSet(ali_file)

    if (length(ali) <= 1) next()

    a_ma <- as.matrix(ali)

    if (!any(a_ma == "-")) next()

    taxon_name <- gsub("\\.fasta$", "", basename(ali_file))
    internal   <- apply(a_ma, 1, has_internal_gap)

    if (any(internal)) {
      path_store <- c(path_store, ali_file)
      message(taxon_name, " - Potential problems in the alignment.")
    } else {
      message(taxon_name, " - Alignment should be fine.")
    }
  }

  if (length(path_store) == 0) {
    message("No problematic alignments found.")
  }

  invisible(path_store)
}
