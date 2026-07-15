#' COI Functionality Index for ASVs
#'
#' @description
#' Computes a Functionality Index (FI) for each ASV/haplotype: a continuous
#' 0-1 score quantifying how compatible a query sequence is with the
#' evolutionary constraints expected of a functional mitochondrial COI
#' protein. FI is a compatibility score, not an estimate of the probability
#' that a sequence is mitochondrial vs. nuclear (NUMT) in origin - see
#' \strong{Details}.
#'
#' @param output_dir Folder where the results of functional calculation are stored.
#' @param grantham_max Normalization ceiling for Grantham distance (default
#'   \code{215}, the theoretical maximum of the Grantham scale). Only used
#'   when \code{include_grantham = TRUE}.
#' @param quality_threshold Minimum per-codon sequencing-quality probability
#'   for a mutated position to be treated as reliable evidence. Only used
#'   when \code{include_quality = TRUE}. Mutations with
#'   \code{median_quality_prob} below this threshold are excluded from the FI
#'   calculation entirely (the call cannot be trusted, so it neither helps
#'   nor hurts the score). Mutations where \code{median_quality_prob} is
#'   \code{NA} (i.e. when VSEARCH read-back was not run) are always treated
#'   as reliable.
#' @param include_position If \code{TRUE} (default) the codon-position
#'   penalty (first/second codon position substitutions) contributes to the
#'   total penalty.
#' @param include_grantham If \code{TRUE} (default) the (normalized) Grantham
#'   distance contributes to each mutated site's penalty.
#' @param include_quality If \code{TRUE} (default) mutated positions with
#'   low read-back quality are excluded from the FI calculation.
#' @param include_hydro If \code{TRUE} (default) the sequence's local
#'   hydrophobicity-profile violation rate (against the family reference
#'   envelope) contributes to the total penalty as an independent,
#'   sequence-level term (it is not tied to any individual mutated site).
#' @param include_conservation If \code{TRUE} (default) each mutated site's
#'   penalty is weighted by how conserved that position is across the
#'   family-level reference alignment (an amino-acid change at an invariant
#'   position counts far more than the same change at a naturally variable
#'   position). If \code{FALSE}, every site is treated as maximally
#'   conserved.
#' @param conservation_k Scaling constant (default \code{20}) controlling how
#'   quickly confidence in the conservation estimate grows with the number of
#'   reference sequences behind a mutated position. The raw entropy-based
#'   conservation weight is multiplied by a confidence factor \code{1 -
#'   exp(-n_ref / conservation_k)}, so a position "conserved" across only two
#'   reference sequences counts for much less than the same apparent
#'   conservation backed by hundreds of references. Positions with an unknown
#'   reference count (\code{n_ref_sequences} is \code{NA}) are left
#'   un-adjusted. Only used when \code{include_conservation = TRUE}.
#'
#' @details
#' For each ASV and taxonomic level, the query is first checked for hard
#' incompatibilities with a functional open reading frame: a premature stop
#' codon, a frameshift, a non-triplet indel, or an otherwise invalid/
#' incomplete ORF. In the current implementation these are all detected
#' collectively via the alignment-gap signal already produced upstream by
#' \code{calculate_asv_distance()} (an internal gap in the query's
#' alignment to the reference profile); if present, \code{FI = 0} and the
#' sequence is classified \emph{"Severe incompatibility with functional
#' COI."}, bypassing the soft-penalty calculation below.
#'
#' For sequences with a valid ORF, a soft penalty is computed for each
#' mutated amino-acid position (a position where the query's amino acid
#' does not match any reference at that taxonomic level):
#'
#' \preformatted{site_penalty = conservation_weight * normalized_grantham}
#'
#' \code{normalized_grantham} is the Grantham distance normalized by
#' \code{grantham_max}. \code{conservation_weight} is:
#'
#' \preformatted{conservation_weight = (1 - normalized_entropy) * (1 - exp(-n_ref / conservation_k))}
#'
#' The first factor is one minus the Shannon entropy of the reference amino
#' acids observed at that position (normalized by \code{log(20)}); the
#' second is a confidence factor that discounts conservation estimates built
#' from few reference sequences (\code{n_ref} - see \code{conservation_k}),
#' so a position that merely appears invariant because only two references
#' cover it is not treated as strongly conserved as one backed by hundreds.
#'
#' Per-site penalties are averaged into \code{mean_site_penalty} and
#' combined with a codon-position penalty:
#'
#' \preformatted{codon_penalty = (Npos1 + 2 * Npos2) / (2 * Nmut)}
#'
#' This weights second-codon-position substitutions twice as heavily as
#' first-position ones, and synonymous third-position changes at zero. A
#' sequence-level hydrophobicity term, \code{normalized_hydro} (the query's
#' local violation rate against the family hydrophobicity envelope, computed
#' by \code{calculate_asv_distance()}), is also added. Hydrophobicity is
#' a whole-protein property rather than a per-residue one, so it is kept as
#' an independent, sequence-level component instead of being folded into
#' \code{site_penalty} - assigning the same sequence-level value to every
#' mutated site would otherwise let it count as evidence once per mutation.
#'
#' The total penalty combines all three terms:
#'
#' \preformatted{total_penalty = 0.35 * mean_site_penalty + 0.3 * codon_penalty + 0.35 * normalized_hydro
#' FI = 1 - total_penalty}
#'
#' bounded to \code{[0, 1]}.
#'
#' The Functionality Index estimates compatibility with functional COI
#' evolution based on codon position, amino-acid conservation, Grantham
#' distance, and hydrophobicity shift. It should not be interpreted as a
#' definitive classification of NUMT versus mitochondrial origin.
#'
#' @return A data frame with one row per ASV per taxonomic level - covering
#'   every ASV/level combination processed upstream, not only ones with a
#'   flagged amino-acid substitution (an ASV whose ORF matches the reference
#'   at every position still gets a row, sourced from the hydrophobicity
#'   summary, with \code{n_aa_substitutions = 0}) - with columns
#'   \code{ASV_id}, \code{tax_lev}, \code{FI} (the Functionality Index,
#'   \code{0}-\code{1}), \code{Class} (the interpretation band: "Plausible
#'   functional sequence" (\code{FI >= 0.90}), "Artifact-NUMTs candidate"
#'   (\code{FI < 0.90}), or "Severe incompatibility with functional COI."
#'   for a hard ORF incompatibility), and the submetrics that combine into
#'   \code{FI}:
#'   \itemize{
#'     \item \code{mean_conservation_weight} - average, confidence-adjusted
#'       \code{conservation_weight} across mutated sites (\code{NA} when there
#'       are none)
#'     \item \code{mean_conservation_confidence} - average of the raw
#'       confidence factor \code{1 - exp(-n_ref / conservation_k)} itself,
#'       reported separately from \code{mean_conservation_weight} so it's
#'       clear how much the entropy-based weight was discounted for sparse
#'       reference coverage (\code{NA} when there are no mutated sites)
#'     \item \code{mean_normalized_grantham} - average normalized Grantham
#'       distance across mutated sites (\code{NA} when there are none)
#'     \item \code{normalized_hydro} - the sequence's normalized local
#'       hydrophobicity violation rate, an independent sequence-level metric
#'       (\code{0} when \code{include_hydro = FALSE})
#'     \item \code{mean_site_penalty} - the averaged per-site penalty, based
#'       on conservation and Grantham distance only (\code{0} when there are
#'       no mutated sites)
#'     \item \code{n_pos1}, \code{n_pos2} - number of mutated sites with a
#'       first/second codon-position change (\code{Npos1}/\code{Npos2})
#'     \item \code{codon_penalty} - the codon-position penalty
#'     \item \code{total_penalty} - \code{0.35 * mean_site_penalty + 0.3 *
#'       codon_penalty + 0.35 * normalized_hydro}, i.e. \code{1 - FI} before
#'       the \code{[0, 1]} bound
#'   }
#'   as well as \code{Evidence} (a human-readable summary),
#'   \code{n_aa_substitutions} (number of mutated amino-acid positions found
#'   at that taxonomic level), \code{n_flagged_aa_substitutions} (how many of
#'   those positions contributed a nonzero site penalty), and
#'   \code{n_ref_sequences} (the number of reference sequences the
#'   comparison was based on - the weakest-covered mutated position sets
#'   this value, and it is \code{NA} when no reference-backed mutation data
#'   is available; a low count means the score rests on thin evidence).
#'
#' @export
#'
#' @importFrom readxl read_xlsx
#' @importFrom dplyr mutate group_by summarise case_when bind_rows left_join full_join select filter right_join

classify_asv <- function(output_dir,
                         grantham_max      = 215,
                         quality_threshold  = 0.999,
                         include_position   = TRUE,
                         include_grantham   = TRUE,
                         include_quality    = TRUE,
                         include_hydro      = TRUE,
                         include_conservation = TRUE,
                         conservation_k     = 20) {

  target_dir <- file.path(output_dir, "aa_structure_results")
  all_files  <- list.files(path = target_dir, pattern = "\\.xlsx$", full.names = TRUE)
  hydro_file <- all_files[grepl("hydrophobicity", all_files, ignore.case = TRUE)]
  asv_files  <- all_files[grepl("_aa_structure", all_files)]

  log20 <- log(20)

  # -- helpers -----------------------------------------------------------------
  aa3 <- function(aa) {
    if (is.na(aa)) return(NA_character_)
    lookup <- c(
      A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", E = "Glu",
      Q = "Gln", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
      M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr", W = "Trp",
      Y = "Tyr", V = "Val", X = "Xaa"
    )
    if (aa %in% names(lookup)) lookup[[aa]] else aa
  }

  codon_pos_label <- function(pos_str) {
    if (is.na(pos_str)) return("unknown codon position")
    ordinals <- c("1" = "first", "2" = "second", "3" = "third")
    positions <- trimws(strsplit(as.character(pos_str), ",")[[1]])
    labels <- ordinals[positions]
    labels <- labels[!is.na(labels)]
    if (length(labels) == 0) return(paste("position", pos_str))
    if (length(labels) == 1) return(paste(labels, "codon position"))
    paste(paste(labels[-length(labels)], collapse = ", "), "and",
          labels[length(labels)], "codon positions")
  }

  # -- hydrophobicity ---------------------------------------------------------
  # hydro_data is also the authoritative list of every ASV/tax_lev combination
  # ever processed upstream (calculate_asv_distance always records a
  # hydrophobicity result, even for an ASV whose ORF matches the reference at
  # every position and so never gets an aa_structure_results file) - it is
  # read whenever the file exists, independent of include_hydro, so those
  # ASVs are still reported rather than silently dropped. include_hydro only
  # controls whether local_violation_rate contributes to the FI score.
  hydro_file_present <- length(hydro_file) > 0
  if (include_hydro && !hydro_file_present) stop("Hydrophobicity file not found!")

  hydro_data <- data.frame()
  if (hydro_file_present) {
    hydro_data <- readxl::read_xlsx(hydro_file[1]) %>%
      dplyr::mutate(global_pass = as.logical(global_pass)) %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(
        local_violation_rate = mean(local_violation_rate, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # -- process one aa_structure file -----------------------------------------
  process_asv_file <- function(f) {
    df_asv <- readxl::read_xlsx(f)
    if (!"aa_pos" %in% names(df_asv)) return(NULL)
    if (!"n_sequences_ali" %in% names(df_asv)) df_asv$n_sequences_ali <- NA_integer_

    gap_info <- df_asv %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(has_alignment_gap = any(aa_pos == -1L, na.rm = TRUE),
                       .groups = "drop")

    norm_df <- dplyr::filter(df_asv, aa_pos != -1L)

    if (nrow(norm_df) == 0) {
      summary_df <- gap_info %>%
        dplyr::mutate(n_aa_substitutions         = 0L,
                      n_flagged_aa_substitutions  = 0L,
                      mean_site_penalty           = NA_real_,
                      mean_conservation_weight    = NA_real_,
                      mean_conservation_confidence = NA_real_,
                      mean_normalized_grantham    = NA_real_,
                      n_pos1                      = 0L,
                      n_pos2                      = 0L,
                      n_ref_sequences             = NA_integer_)
      return(list(summary = summary_df, details = data.frame()))
    }

    # a mutated position is reliable evidence only when its read-back quality
    # meets the threshold; unreliable positions are dropped entirely rather
    # than penalised or rewarded, since the call itself cannot be trusted
    reliable <- norm_df %>%
      dplyr::mutate(
        is_reliable = (!include_quality | is.na(median_quality_prob) |
                        median_quality_prob >= quality_threshold)
      ) %>%
      dplyr::filter(is_reliable)

    if (nrow(reliable) == 0) {
      summary_df <- gap_info %>%
        dplyr::mutate(n_aa_substitutions         = 0L,
                      n_flagged_aa_substitutions  = 0L,
                      mean_site_penalty           = NA_real_,
                      mean_conservation_weight    = NA_real_,
                      mean_conservation_confidence = NA_real_,
                      mean_normalized_grantham    = NA_real_,
                      n_pos1                      = 0L,
                      n_pos2                      = 0L,
                      n_ref_sequences             = NA_integer_)
      return(list(summary = summary_df, details = data.frame()))
    }

    # conservation: Shannon entropy of the reference amino acids observed at
    # each mutated position, across ALL matching reference sequences (not
    # just the most charitable one) - normalized by log(20). The raw entropy
    # weight is then discounted by a confidence factor that grows with the
    # number of reference sequences behind the position, so a handful of
    # identical references isn't scored as confidently conserved as hundreds.
    conservation_per_pos <- reliable %>%
      dplyr::group_by(ASV_id, tax_lev, aa_pos) %>%
      dplyr::summarise(
        conservation_weight = {
          tab <- table(ref_aa)
          p   <- as.numeric(tab) / sum(tab)
          norm_entropy <- if (length(p) <= 1) 0 else -sum(p * log(p)) / log20
          1 - min(norm_entropy, 1)
        },
        n_ref_sequences = if (all(is.na(n_sequences_ali))) NA_integer_
                          else as.integer(min(n_sequences_ali, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        conservation_confidence = ifelse(is.na(n_ref_sequences), 1,
                                         1 - exp(-n_ref_sequences / conservation_k)),
        conservation_weight = conservation_weight * conservation_confidence
      )

    # label reference: pick the most charitable (lowest-Grantham) reference
    # per position to report the substitution and its Grantham distance
    label_per_pos <- reliable %>%
      dplyr::group_by(ASV_id, tax_lev, aa_pos) %>%
      dplyr::slice_min(order_by = grantham_dist, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(ASV_id, tax_lev, aa_pos, illegal_aa, ref_aa,
                    triplet_mut_pos, grantham_dist)

    best_per_pos <- dplyr::left_join(conservation_per_pos, label_per_pos,
                                     by = c("ASV_id", "tax_lev", "aa_pos"))

    list(summary = NULL, details = best_per_pos, gap_info = gap_info)
  }

  processed   <- Filter(Negate(is.null), lapply(asv_files, process_asv_file))

  # summaries for gap-only / no-reliable-evidence files, computed directly above
  early_summaries <- dplyr::bind_rows(lapply(processed, function(x) x$summary))

  all_details <- dplyr::bind_rows(lapply(processed, function(x) x$details))
  all_gap_info <- dplyr::bind_rows(lapply(processed, function(x) x$gap_info))
  if (nrow(all_details) == 0) {
    all_details <- data.frame(ASV_id = character(), tax_lev = character())
  }

  # -- per-site penalty: conservation + Grantham only --------------------------
  # hydrophobicity is a whole-sequence property (see below) and is kept out of
  # the per-site penalty so it isn't counted once per mutated position.
  if (nrow(all_details) > 0) {
    all_details <- all_details %>%
      dplyr::mutate(
        normalized_grantham  = if (include_grantham) {
                                 ifelse(is.na(grantham_dist), 0, pmin(grantham_dist / grantham_max, 1))
                               } else 0,
        conservation_weight  = if (include_conservation) conservation_weight else 1,
        site_penalty         = conservation_weight * normalized_grantham,
        is_pos1              = include_position & grepl("1", as.character(triplet_mut_pos)),
        is_pos2              = include_position & grepl("2", as.character(triplet_mut_pos))
      )

    computed_summary <- all_details %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(
        n_aa_substitutions        = dplyr::n(),
        n_flagged_aa_substitutions = sum(site_penalty > 0, na.rm = TRUE),
        mean_site_penalty         = mean(site_penalty, na.rm = TRUE),
        mean_conservation_weight = mean(conservation_weight, na.rm = TRUE),
        mean_conservation_confidence = mean(conservation_confidence, na.rm = TRUE),
        mean_normalized_grantham = mean(normalized_grantham, na.rm = TRUE),
        n_pos1                    = sum(is_pos1, na.rm = TRUE),
        n_pos2                    = sum(is_pos2, na.rm = TRUE),
        n_ref_sequences           = if (all(is.na(n_ref_sequences))) NA_integer_
                                    else min(n_ref_sequences, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::right_join(all_gap_info, by = c("ASV_id", "tax_lev"))
  } else {
    computed_summary <- data.frame()
  }

  all_summary <- dplyr::bind_rows(early_summaries, computed_summary)

  # -- list every processed ASV, not only the flagged ones ---------------------
  # hydro_data covers every ASV/tax_lev combination that was ever compared
  # against a reference (see note above), so full-joining against it brings in
  # ASVs with a clean ORF at this level - these have no mutated sites at all
  # and would otherwise be missing from the output entirely.
  if (nrow(hydro_data) > 0) {
    all_summary <- if (nrow(all_summary) == 0) hydro_data
                   else dplyr::full_join(all_summary, hydro_data, by = c("ASV_id", "tax_lev"))
  } else if (!"local_violation_rate" %in% names(all_summary)) {
    all_summary$local_violation_rate <- NA_real_
  }

  # -- early exit ------------------------------------------------------------
  if (nrow(all_summary) == 0) {
    return(data.frame(ASV_id = character(), tax_lev = character(),
                      FI = numeric(), Class = character(),
                      mean_conservation_weight = numeric(),
                      mean_conservation_confidence = numeric(),
                      mean_normalized_grantham = numeric(),
                      normalized_hydro = numeric(),
                      mean_site_penalty = numeric(),
                      n_pos1 = integer(), n_pos2 = integer(),
                      codon_penalty = numeric(), total_penalty = numeric(),
                      Evidence = character(),
                      n_aa_substitutions = integer(),
                      n_flagged_aa_substitutions = integer(),
                      n_ref_sequences = integer()))
  }

  for (col in c("has_alignment_gap", "n_aa_substitutions", "n_flagged_aa_substitutions",
                "mean_site_penalty", "mean_conservation_weight", "mean_conservation_confidence",
                "mean_normalized_grantham", "n_pos1", "n_pos2", "n_ref_sequences")) {
    if (!col %in% names(all_summary)) all_summary[[col]] <- NA
  }

  # -- Functionality Index -----------------------------------------------------
  verdicts <- all_summary %>%
    dplyr::mutate(
      has_alignment_gap  = ifelse(is.na(has_alignment_gap), FALSE, has_alignment_gap),
      n_aa_substitutions = ifelse(is.na(n_aa_substitutions), 0L, n_aa_substitutions),
      n_flagged_aa_substitutions = ifelse(is.na(n_flagged_aa_substitutions), 0L, n_flagged_aa_substitutions),
      n_pos1             = ifelse(is.na(n_pos1), 0L, n_pos1),
      n_pos2             = ifelse(is.na(n_pos2), 0L, n_pos2),
      mean_site_penalty  = ifelse(n_aa_substitutions == 0 | is.na(mean_site_penalty), 0, mean_site_penalty),
      mean_conservation_weight = ifelse(n_aa_substitutions == 0, NA_real_, mean_conservation_weight),
      mean_conservation_confidence = ifelse(n_aa_substitutions == 0, NA_real_, mean_conservation_confidence),
      mean_normalized_grantham = ifelse(n_aa_substitutions == 0, NA_real_, mean_normalized_grantham),
      normalized_hydro   = if (include_hydro) {
                             ifelse(is.na(local_violation_rate), 0, pmin(local_violation_rate, 1))
                           } else 0,
      codon_penalty      = ifelse(n_aa_substitutions == 0, 0,
                                  (n_pos1 + 2 * n_pos2) / (2 * n_aa_substitutions)),
      total_penalty      = 0.35 * mean_site_penalty + 0.3 * codon_penalty +
                           0.35 * normalized_hydro,
      FI                 = pmax(0, pmin(1, 1 - total_penalty)),
      FI                 = ifelse(has_alignment_gap, 0, FI),
      Class = dplyr::case_when(
        has_alignment_gap ~ "Severe incompatibility with functional COI.",
        FI >= 0.90         ~ "Plausible functional sequence",
        TRUE               ~ "Artifact-NUMTs candidate"
      )
    )

  # -- evidence strings -------------------------------------------------------
  build_evidence <- function(i) {
    row      <- verdicts[i, ]
    asv_id   <- row$ASV_id
    cur_lev  <- row$tax_lev
    is_gap   <- isTRUE(row$has_alignment_gap)

    parts <- character(0)

    if (is_gap) {
      parts <- c(parts, paste0("query has a premature stop codon, frameshift, ",
                               "or non-triplet indel relative to the ", cur_lev, "-level reference"))
    } else {
      muts <- all_details[all_details$ASV_id == asv_id &
                          all_details$tax_lev == cur_lev, ]

      if (nrow(muts) == 0) {
        parts <- c(parts,
                   paste0("no novel amino acid substitutions at ", cur_lev, " level"))
      } else {
        for (k in seq_len(nrow(muts))) {
          m       <- muts[k, ]
          from_aa <- aa3(m$ref_aa)
          to_aa   <- aa3(m$illegal_aa)
          pos_lab <- codon_pos_label(m$triplet_mut_pos)
          gran    <- m$grantham_dist

          desc <- paste0(from_aa, "->", to_aa, " substitution; ", pos_lab)
          if (include_grantham && !is.na(gran))
            desc <- paste0(desc, "; Grantham = ", round(gran))
          if (include_conservation)
            desc <- paste0(desc, "; ",
                           if (m$conservation_weight >= 0.5) "conserved" else "variable",
                           " reference position (conservation weight = ",
                           round(m$conservation_weight, 2), ", confidence = ",
                           round(m$conservation_confidence, 2), ")")
          desc <- paste0(desc, "; site penalty = ", round(m$site_penalty, 2))
          parts <- c(parts, desc)
        }
      }
    }

    parts <- c(parts, paste0("FI = ", round(row$FI, 3), " (", row$Class, ")"))

    # cross-level context using verdicts at all levels for this ASV
    all_v <- verdicts[verdicts$ASV_id == asv_id, ]
    others <- all_v[all_v$tax_lev != cur_lev, ]
    if (nrow(others) > 0) {
      other_desc <- paste0(others$tax_lev, " FI = ", round(others$FI, 2))
      parts <- c(parts, paste0("other levels: ", paste(other_desc, collapse = "; ")))
    }

    paste(parts, collapse = "; ")
  }

  evidences <- vapply(seq_len(nrow(verdicts)), build_evidence, character(1))

  # -- final output -----------------------------------------------------------
  verdicts %>%
    dplyr::mutate(Evidence                 = evidences,
                 FI                        = round(FI, 3),
                 mean_conservation_weight = round(mean_conservation_weight, 3),
                 mean_conservation_confidence = round(mean_conservation_confidence, 3),
                 mean_normalized_grantham = round(mean_normalized_grantham, 3),
                 normalized_hydro         = round(normalized_hydro, 3),
                 mean_site_penalty        = round(mean_site_penalty, 3),
                 codon_penalty            = round(codon_penalty, 3),
                 total_penalty            = round(total_penalty, 3)) %>%
    dplyr::select(ASV_id, tax_lev, FI, Class,
                 mean_conservation_weight, mean_conservation_confidence,
                 mean_normalized_grantham, normalized_hydro,
                 mean_site_penalty, n_pos1, n_pos2, codon_penalty, total_penalty,
                 Evidence, n_aa_substitutions, n_flagged_aa_substitutions, n_ref_sequences)
}
