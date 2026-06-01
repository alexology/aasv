#' ASV Functional Structure Analysis
#'
#' @description
#' Processes ASV taxonomy data, organizes output directories, and performs
#' functional structural analysis by comparing sequences against reference
#' databases at species, genus, and family levels.
#'
#' @param asv_taxonomy A data frame containing ASV sequences and taxonomic
#'   information. Must include columns \code{species}, \code{genus},
#'   \code{family}, and \code{ASV}.
#' @param output_dir Path to the directory where results will be saved.
#' @param alignment_dir Path to the directory containing reference alignments
#'   produced by \code{\link{align_species_seq}}, \code{\link{align_genus_seq}},
#'   and \code{\link{align_family_seq}}.
#' @param vsearch_path Path to the VSEARCH executable. Set to \code{NULL}
#'   (default) to skip read-back mapping; \code{median_quality_prob} will be
#'   \code{NA} in the output.
#' @param trimmed_folder Path to the folder containing trimmed sequences.
#'   Required when \code{vsearch_path} is provided; ignored otherwise.
#' @param ee Expected error threshold (default: 40).
#' @param genetic_code The genetic code for translation (default: \code{"5"},
#'   invertebrate mitochondrial).
#' @param hydro_window Window size for hydrophobicity calculation (default: 7).
#' @param hydro_threshold Threshold for hydrophobicity (default: 0.1).
#' @param ... Further arguments passed to \code{calculate_asv_distance}.
#'
#' @details
#' The function follows a hierarchical workflow for each ASV:
#'
#' \bold{1. ID generation:} \code{asv_taxonomy} is sorted by species and
#' zero-padded identifiers are assigned (e.g. ASV001, ASV100).
#'
#' \bold{2. Reference mapping:} For each ASV the function searches for
#' reference sequences at the species, genus, and family level in that order.
#' Levels for which reference files are missing are skipped; all levels where
#' data are available are processed and their results combined.
#'
#' \bold{3. Functional calculation:} \code{calculate_asv_distance()} is called
#' for each available level. It performs amino acid alignment, translation
#' using \code{genetic_code}, and sliding-window hydrophobicity analysis.
#' When \code{vsearch_path} is provided, reads from \code{trimmed_folder}
#' are also mapped back to each ASV to assess per-codon sequencing quality.
#'
#' \bold{4. Output:} Results are written to \code{output_dir}:
#' \itemize{
#'   \item \code{asv_taxonomy_revised.xlsx} - input taxonomy with ASV IDs.
#'   \item \code{aa_structure_results/<ASV_id>_aa_structure.xlsx} - per-ASV
#'     amino acid structural data (written only when data exist for at least
#'     one level).
#'   \item \code{aa_structure_results/hydrophobicity.xlsx} - summary across
#'     all processed ASVs.
#' }
#'
#' @return The enriched \code{asv_taxonomy} data frame with the added
#'   \code{ASV_id} column, returned invisibly.
#'
#' @importFrom dplyr arrange mutate select where bind_rows
#' @importFrom writexl write_xlsx
#' @importFrom Biostrings readAAStringSet readDNAStringSet
#' @importFrom tibble tibble
#'
#' @export

asv_functional_structure <- function(asv_taxonomy,
                                     output_dir,
                                     alignment_dir,
                                     vsearch_path    = NULL,
                                     trimmed_folder  = NULL,
                                     ee = 40,
                                     genetic_code = "5",
                                     hydro_window = 7,
                                     hydro_threshold = 0.1,
                                     ...) {

  # reads and concatenates one or more DNA FASTA files into a DNAStringSet
  read_dna_files <- function(paths) {
    do.call(c, lapply(paths, Biostrings::readDNAStringSet))
  }

  # sort and assign zero-padded ASV IDs
  asv_taxonomy <- dplyr::arrange(asv_taxonomy, species)
  n_asvs       <- nrow(asv_taxonomy)
  id_fmt       <- paste0("ASV%0", nchar(as.character(n_asvs)), "d")
  asv_taxonomy <- dplyr::mutate(asv_taxonomy,
                                ASV_id = sprintf(id_fmt, seq_len(n_asvs)))

  writexl::write_xlsx(asv_taxonomy,
                      file.path(output_dir, "asv_taxonomy_revised.xlsx"))

  sample_columns <- dplyr::select(asv_taxonomy, dplyr::where(is.numeric))

  # create output dirs and purge any stale files from a previous run
  for (d in c("aa_structure_results", "vsearch_map_back", "query_alignments")) {
    full  <- file.path(output_dir, d)
    dir.create(full, showWarnings = FALSE)
    stale <- list.files(full, full.names = TRUE)
    if (length(stale) > 0) invisible(file.remove(stale))
  }

  # ── colour-depth probe ───────────────────────────────────────────────────
  # 0  = no colour (NO_COLOR set or non-interactive)
  # 8  = basic ANSI 16-colour fallback
  # 256 = 256-colour (RStudio, xterm-256color, …)
  .nc <- if (nchar(Sys.getenv("NO_COLOR")) > 0) {
    0L
  } else if (nchar(Sys.getenv("RSTUDIO")) > 0 ||
             grepl("256", Sys.getenv("TERM"),      fixed = TRUE) ||
             grepl("256", Sys.getenv("COLORTERM"), fixed = TRUE)) {
    256L
  } else {
    8L
  }

  if (.nc >= 256L) {
    # 256-colour: 172 = #d78700 (≈ #D97706), 33 = #0087ff (vivid blue)
    .orange <- "\x1b[48;5;172m\x1b[38;5;0m"
    .blue   <- "\x1b[38;5;33m"
    .off    <- "\x1b[0m"
  } else if (.nc >= 8L) {
    # basic ANSI fallback: yellow bg + black fg / bright blue
    .orange <- "\x1b[43m\x1b[30m"
    .blue   <- "\x1b[94m"
    .off    <- "\x1b[0m"
  } else {
    .orange <- .blue <- .off <- ""
  }

  msg_asv   <- function(...) message(paste0(.orange, ..., .off))
  msg_level <- function(...) message(paste0(.blue,   ..., .off))

  tax_levels <- c("species", "genus", "family")
  nt_base <- function(fam) {
    file.path(alignment_dir, "raw_sequences_nt", fam, "species")
  }
  res_hydro  <- list()

  for (i in seq_len(n_asvs)) {

    asv_id_i    <- asv_taxonomy$ASV_id[i]
    species_i   <- asv_taxonomy$species[i]
    species_u_i <- gsub(" ", "_", species_i)
    genus_i     <- asv_taxonomy$genus[i]
    family_i    <- asv_taxonomy$family[i]
    asv_i       <- asv_taxonomy$ASV[i]

    # sample columns where this ASV is present (abundance > 0)
    asv_samples_i <- colnames(sample_columns)[
      as.logical(sample_columns[i, ] > 0)
    ]

    label_i <- if (!is.na(species_i)) species_i else if (!is.na(genus_i)) genus_i else family_i
    res_aa  <- list()

    for (z in tax_levels) {

      # resolve reference file paths for this taxonomic level
      if (z == "species") {
        aa_path  <- file.path(alignment_dir, "alignments", family_i, "species",
                              paste0(species_u_i, ".fasta"))
        nt_paths <- file.path(nt_base(family_i), paste0(species_u_i, ".fasta"))

      } else if (z == "genus") {
        aa_path  <- file.path(alignment_dir, "alignments", family_i, "genus",
                              paste0(genus_i, ".fasta"))
        nt_paths <- list.files(nt_base(family_i),
                               pattern    = paste0("^", genus_i, "_"),
                               full.names = TRUE)
      } else {
        aa_path  <- file.path(alignment_dir, "alignments", family_i, "family",
                              paste0(family_i, ".fasta"))
        nt_paths <- list.files(nt_base(family_i),
                               pattern    = "\\.fasta$",
                               full.names = TRUE)
      }

      ref_nt_ok <- length(nt_paths) > 0 && any(file.exists(nt_paths))
      if (!file.exists(aa_path) || !ref_nt_ok) next()

      aa_ref_i <- Biostrings::readAAStringSet(aa_path)
      nt_ref_i <- read_dna_files(nt_paths[file.exists(nt_paths)])

      res_z_i <- calculate_asv_distance(
        output_dir      = output_dir,
        vsearch_path    = vsearch_path,
        query_dna       = asv_i,
        query_id        = asv_id_i,
        ref_aa_aligned  = aa_ref_i,
        species_name    = species_u_i,
        family_name     = family_i,
        ref_nt          = nt_ref_i,
        asv_samples     = asv_samples_i,
        genetic_code    = genetic_code,
        hydro_threshold = hydro_threshold,
        hydro_window    = hydro_window,
        trimmed_folder  = trimmed_folder,
        ee              = ee,
        tax_lev         = z,
        ...
      )

      res_hydro <- c(res_hydro, list(tibble::tibble(
        tax_lev = z, species = species_i, ASV_id = asv_id_i, res_z_i$hydro
      )))

      if (length(res_z_i) > 1 && "aa_pos" %in% names(res_z_i$aa)) {
        res_aa[[z]] <- tibble::tibble(
          tax_lev = z, species = species_i, ASV_id = asv_id_i, res_z_i$aa
        )
      }

      msg_level(label_i, " - ", z, ": done")
    }

    if (length(res_aa) > 0) {
      writexl::write_xlsx(
        dplyr::bind_rows(res_aa),
        file.path(output_dir, "aa_structure_results",
                  paste0(asv_id_i, "_aa_structure.xlsx"))
      )
    }

    msg_asv("ASV ", i, "/", n_asvs, " - ", label_i)
  }

  writexl::write_xlsx(
    dplyr::bind_rows(res_hydro),
    file.path(output_dir, "aa_structure_results", "hydrophobicity.xlsx")
  )

  invisible(asv_taxonomy)
}
