#' Download and align species-level sequences
#'
#' @description
#' This function downloads and aligns sequences from species within the same family
#' of the query sequences.
#'
#' @param taxonomy A file with the taxonomy assigned to ASV. It must have at least
#'   3 columns: family, genus and species.
#' @param alignment_dir Path to a folder to store the results.
#' @param overwrite Overwrite existing data.
#' @param timeout Download timeout in seconds passed to \code{options(timeout)}.
#' @param min_length Select only the downloaded sequences longer than this value.
#' @param max_length Select only the downloaded sequences shorter than this value.
#' @param genetic_code Genetic code. Default to invertebrate mitochondrial code.
#' @param ... Further arguments to be passed to \code{AlignSeqs} of the \code{DECIPHER} package.
#'
#' @details
#' This function downloads sequences from BOLD to build a reference database that is used for
#' further calculations. Given a set of query species, \code{align_species_seq} downloads
#' all the available sequences up to family level. Only the sequences with a full
#' species name are retained (e.g. Gorilla sp. is removed). Once the sequences are
#' downloaded and filtered, they are aligned for each species separately with
#' the function \code{AlignSeqs} of the package \code{DECIPHER}. Raw sequences are
#' saved to disk in a folder called \code{raw_sequences}, while the aligned
#' sequences are saved in a folder called \code{alignments}.
#'
#' @return Called for its side effects (downloading raw sequences and writing
#'   aligned FASTA files to disk). Returns \code{NULL} invisibly.
#'
#' @importFrom httr2 request req_url_query req_perform resp_body_json resp_body_string
#' @importFrom jsonlite fromJSON
#' @importFrom purrr pluck map list_rbind
#' @importFrom rlang %||%
#' @importFrom DECIPHER AlignSeqs
#' @importFrom dplyr bind_rows filter pull mutate select inner_join group_by cur_group_id ungroup
#' @importFrom Biostrings DNAStringSet AAStringSet width writeXStringSet
#' @importFrom ShortRead clean
#' @importFrom tibble tibble
#' @importFrom stats na.omit
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' align_species_seq(
#'   data.frame(family  = "Oligoneuridae",
#'              genus   = "Oligoneuriella",
#'              species = "Oligoneuriella rhenana"),
#'   min_length  = 640,
#'   max_length  = 700,
#'   terminalGap = -1
#' )
#' }


align_species_seq <- function(taxonomy = NULL,
                              alignment_dir = NULL,
                              overwrite = FALSE,
                              timeout = 7200,
                              min_length = 500,
                              max_length = 1000,
                              genetic_code = "5",
                              ...) {

  if (is.null(taxonomy) || !all(c("family", "genus", "species") %in% names(taxonomy))) {
    stop("`taxonomy` must be a data.frame with at least columns: family, genus, species.")
  }

  families <- unique(na.omit(taxonomy$family))

  if (is.null(alignment_dir)) {
    alignment_dir <- getwd()
  }

  if (!dir.exists(alignment_dir)) {
    dir.create(alignment_dir)
  }

  # create top-level output folders once, outside the family loop
  for (d in c("raw_sequences", "alignments", "raw_sequences_nt", "aa_tables")) {
    dir.create(file.path(alignment_dir, d), showWarnings = FALSE)
  }

  old_timeout <- options(timeout = timeout)
  on.exit(options(old_timeout), add = TRUE)

  parse_record <- function(x) {
    tryCatch({
      r <- jsonlite::fromJSON(x, simplifyDataFrame = FALSE)
      tibble::tibble(
        processid      = r$processid      %||% NA_character_,
        identification = r$identification %||% NA_character_,
        marker         = sub(".*\\.", "", r$record_id %||% ""),
        sequence       = r$nuc            %||% NA_character_
      )
    }, error = function(e) NULL)
  }

  for (taxon in families) {

    if (dir.exists(file.path(alignment_dir, "raw_sequences", taxon)) && !overwrite) {
      next()
    }

    for (d in c("raw_sequences", "raw_sequences_nt", "alignments", "aa_tables")) {
      dir.create(file.path(alignment_dir, d, taxon), showWarnings = FALSE)
    }
    for (d in c("raw_sequences", "raw_sequences_nt", "alignments")) {
      dir.create(file.path(alignment_dir, d, taxon, "species"), showWarnings = FALSE)
    }

    query_id <- httr2::request(
      "https://portal.boldsystems.org/api/query"
    ) |>
      httr2::req_url_query(
        query = paste0("tax:family:", taxon),
        extent = "full"
      ) |>
      httr2::req_perform() |>
      httr2::resp_body_json() |>
      purrr::pluck("query_id")

    if (is.null(query_id)) {
      message("No query_id returned for family: ", taxon, " - skipping.")
      next()
    }

    doc_url <- paste0("https://portal.boldsystems.org/api/documents/", query_id)

    n_records <- 0L
    for (i in seq_len(12)) {
      n_records <- httr2::request(doc_url) |>
        httr2::req_url_query(length = 0L, start = 0L) |>
        httr2::req_perform() |>
        httr2::resp_body_json() |>
        purrr::pluck("recordsTotal", .default = 0L)
      if (n_records > 0L) break
      Sys.sleep(5)
    }

    if (n_records == 0L) {
      message("No records found for family: ", taxon, " - skipping.")
      next()
    }

    content <- httr2::request(paste0(doc_url, "/download?format=json")) |>
      httr2::req_perform() |>
      httr2::resp_body_string()

    lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(trimws(lines))]

    id_query_info <- purrr::map(lines, parse_record) |>
      purrr::list_rbind() |>
      dplyr::filter(
        marker == "COI-5P",
        !is.na(sequence),
        sequence != ""
      )
    
    
    
    reference_species <- id_query_info |>
      dplyr::filter(
        !grepl("cf.", identification, fixed = TRUE),
        !grepl("sp.", identification, fixed = TRUE),
        !grepl("/",   identification, fixed = TRUE),
        !grepl("[0-9]", identification),
        countSpaces(identification) == 1,
        identification != ""
      )

    species_unique <- unique(reference_species$identification)
    message("Download ", taxon, " sequences...")

    for (sp in species_unique) {

      id_query          <- dplyr::filter(id_query_info, identification == sp)
      id_query_seqs     <- dplyr::pull(id_query, sequence)
      names(id_query_seqs) <- gsub(" ", "_", paste0(id_query$identification, "__", id_query$processid))

      # remove NAs and sequences with non-standard characters
      id_query_seqs <- id_query_seqs[!is.na(id_query_seqs)]
      id_query_seqs <- id_query_seqs[!grepl("[^ACGTN+.-]", id_query_seqs)]

      if (length(id_query_seqs) == 0) {
        message(sp, ": without sequences")
        next()
      }

      id_query_seqs <- Biostrings::DNAStringSet(id_query_seqs)
      id_query_seqs <- ShortRead::clean(id_query_seqs)

      w             <- Biostrings::width(id_query_seqs)
      id_query_seqs <- id_query_seqs[w >= min_length & w <= max_length]

      if (length(id_query_seqs) == 0) {
        message(sp, ": without sequences")
        next()
      }

      id_query_aa <- as.character(id_query_seqs) |>
        sapply(best_translation, genetic_code = genetic_code)

      # drop failed translations and keep nt sequences in sync
      valid         <- !is.na(id_query_aa)
      id_query_aa   <- id_query_aa[valid]
      id_query_seqs <- id_query_seqs[names(id_query_aa)]

      id_query_aa <- Biostrings::AAStringSet(id_query_aa)

      # drop sequences with internal stop codons, keeping nt and aa in sync
      has_stop      <- grepl("*", as.character(id_query_aa), fixed = TRUE)
      id_query_seqs <- id_query_seqs[!has_stop]
      id_query_aa   <- id_query_aa[!has_stop]

      if (length(id_query_aa) == 0) next()

      # build sequence index after all filtering so lengths are consistent
      seq_df <- tibble::tibble(
        processid = names(id_query_aa),
        sequence  = as.character(id_query_aa)
      ) |>
        dplyr::group_by(sequence) |>
        dplyr::mutate(id = paste0("seq_", dplyr::cur_group_id())) |>
        dplyr::ungroup()

      species_name <- gsub(" ", "_", sp)

      paths <- list(
        raw     = file.path(alignment_dir, "raw_sequences",    taxon, "species", paste0(species_name, ".fasta")),
        raw_nt  = file.path(alignment_dir, "raw_sequences_nt", taxon, "species", paste0(species_name, ".fasta")),
        aligned = file.path(alignment_dir, "alignments",       taxon, "species", paste0(species_name, ".fasta")),
        aa_tab  = file.path(alignment_dir, "aa_tables",        taxon,            paste0(species_name, ".txt"))
      )

      Biostrings::writeXStringSet(id_query_aa,   filepath = paths$raw,    width = max_length + 1, format = "fasta")
      Biostrings::writeXStringSet(id_query_seqs, filepath = paths$raw_nt, width = max_length + 1, format = "fasta")
      write.table(aa_table(id_query_seqs, genetic_code = genetic_code), file = paths$aa_tab)

      aa_unique         <- id_query_aa
      names(aa_unique)  <- seq_df$id
      aa_unique         <- unique(aa_unique)

      if (length(aa_unique) == 1) {

        Biostrings::writeXStringSet(id_query_aa, filepath = paths$aligned, width = max_length + 1, format = "fasta")

      } else {

        aligned <- DECIPHER::AlignSeqs(aa_unique, ...)

        aligned_df <- tibble::tibble(
          id       = names(aligned),
          sequence = as.character(aligned)
        ) |>
          dplyr::inner_join(dplyr::select(seq_df, processid, id), by = "id")

        # restore original sequence order
        aligned_df  <- aligned_df[match(names(id_query_aa), aligned_df$processid), ]

        aligned_out <- Biostrings::AAStringSet(aligned_df$sequence)
        names(aligned_out) <- aligned_df$processid

        Biostrings::writeXStringSet(aligned_out, filepath = paths$aligned, width = max_length + 1, format = "fasta")
      }

      message(sp, ": done")
    }
  }
}
