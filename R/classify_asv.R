#' Filter ASV by functional properties
#'
#' @description
#' This function takes ASVs functional description as input to classify them as 
#' as laboratory artifact or NUMTs.  
#' 
#' @param output_dir Folder where the results of functional calculation are stored.
#' @param grantham_threshold The threshold to retain an amino acid substitution as functional.
#' @param quality_threshold Minimum probability to retain a codon. Ignored for
#'   mutations where \code{median_quality_prob} is \code{NA} (i.e. when VSEARCH
#'   read-back was not run); those mutations are treated as passing quality.
#' @param include_hydro If \code{TRUE} includes the hydrophobicity assessment.
#'
#' @details
#' The function \code{classify_asv} uses hydrophobicity and functional results to
#' to assess if an ASV is likely a laboratory error or a NUMT. For each mutation, this 
#' function select the best result among those available for each taxonomic level separately.
#' When an ASV has more than one mutation, the worst result is choosen to provide
#' the final assessment. Since hydrophobicity is often penalizing and the user can 
#' decide to include it or not in the evaluation. 
#' 
#' 
#' @return A data frame with one row per ASV per taxonomic level, with columns
#'   \code{ASV_id}, \code{tax_lev}, \code{IS_REAL}, \code{REASON},
#'   \code{max_of_best_grantham}, and \code{min_of_best_quality}.
#'
#' @export
#'
#' @importFrom readxl read_xlsx
#' @importFrom dplyr mutate group_by summarise case_when bind_rows left_join full_join select filter

classify_asv <- function(output_dir, 
                         grantham_threshold = 50, 
                         quality_threshold = 0.999,
                         include_hydro = TRUE) {
  
  # set up the path to functional results
  target_dir <- file.path(output_dir, "aa_structure_results")
  
  # list all the files of functional results
  all_files <- list.files(path = target_dir, pattern = "\\.xlsx$", full.names = TRUE)
  
  # identify the hydrophobicity results
  hydro_file <- all_files[grepl("hydrophobicity", all_files, ignore.case = TRUE)]
  
  # functional results based on codon evaluation
  asv_files  <- all_files[grepl("_aa_structure", all_files)]
  
  # process hydrophobicity results
  hydro_data <- data.frame()
  if (include_hydro) {
    if (length(hydro_file) == 0) stop("Hydrophobicity file not found!")
    hydro_data <- readxl::read_xlsx(hydro_file[1]) %>%
      dplyr::mutate(global_pass = as.logical(global_pass)) %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(
        hydro_pass_level = any(global_pass == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # process functional assessment results
  asv_tax_summaries <- lapply(asv_files, function(f) {
    df_asv <- readxl::read_xlsx(f)

    # File was written when all references had gaps at every queried position;
    # no codon-level data is available — skip rather than crashing.
    if (!"aa_pos" %in% names(df_asv)) return(NULL)

    # check for alignment gaps
    gap_info <- df_asv %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(has_alignment_gap = any(aa_pos == -1, na.rm = TRUE), .groups = "drop")
    
    # find the best result for each mutation
    best_per_mutation <- df_asv %>%
      dplyr::filter(aa_pos != -1) %>%
      dplyr::mutate(
        is_valid_ref = dplyr::case_when(
          (!grepl("2", as.character(triplet_mut_pos)) | grantham_dist <= grantham_threshold) &
            (is.na(median_quality_prob) |                  # quality not assessed (no VSEARCH)
               median_quality_prob >= quality_threshold |
               grantham_dist <= grantham_threshold) ~ TRUE,
          TRUE ~ FALSE
        )
      ) %>%
      dplyr::group_by(ASV_id, tax_lev, aa_pos) %>%
      dplyr::summarise(
        mutation_ok = any(is_valid_ref == TRUE, na.rm = TRUE),
        best_grantham_at_pos = if(all(is.na(grantham_dist))) NA_real_ else min(grantham_dist, na.rm = TRUE),
        best_quality_at_pos = if(all(is.na(median_quality_prob))) NA_real_ else max(median_quality_prob, na.rm = TRUE),
        .groups = "drop"
      )
    
    # aggregate at ASV level for each taxonomic level
    asv_summary <- best_per_mutation %>%
      dplyr::group_by(ASV_id, tax_lev) %>%
      dplyr::summarise(
        all_mutations_ok = all(mutation_ok == TRUE, na.rm = TRUE),
        max_of_best_grantham = if(all(is.na(best_grantham_at_pos))) NA_real_ else max(best_grantham_at_pos, na.rm = TRUE),
        min_of_best_quality = if(all(is.na(best_quality_at_pos))) NA_real_ else min(best_quality_at_pos, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::right_join(gap_info, by = c("ASV_id", "tax_lev"))
    
    return(asv_summary)
  })
  
  # merge the results
  combined_asv <- dplyr::bind_rows(asv_tax_summaries)
  
  # FIX: Handle case where no asv_files were found
  if (nrow(combined_asv) == 0) {
    if (include_hydro && nrow(hydro_data) > 0) {
      combined_asv <- hydro_data
    } else {
      # Return an empty dataframe with the expected columns if no data is found at all
      return(data.frame(ASV_id=character(), tax_lev=character(), IS_REAL=logical(), 
                        REASON=character(), max_of_best_grantham=numeric(), 
                        min_of_best_quality=numeric()))
    }
  } else if (include_hydro) {
    # If both exist, merge them to ensure all ASVs from hydro are present
    combined_asv <- combined_asv %>% dplyr::full_join(hydro_data, by = c("ASV_id", "tax_lev"))
  }
  
  # Ensure all necessary mutation columns exist (filling with NA if they were missing due to no asv_files)
  for (col in c("all_mutations_ok", "max_of_best_grantham", "min_of_best_quality", "has_alignment_gap")) {
    if (!col %in% names(combined_asv)) combined_asv[[col]] <- NA
  }
  
  if (!include_hydro) {
    combined_asv$hydro_pass_level <- TRUE
  }
  
  # combine the results
  final_report <- combined_asv %>%
    dplyr::mutate(
      # if a value is NA (e.g. no mutation file found), we assume it's OK/False depending on context
      has_alignment_gap = ifelse(is.na(has_alignment_gap), FALSE, has_alignment_gap),
      all_mutations_ok = ifelse(is.na(all_mutations_ok), TRUE, all_mutations_ok),
      hydro_pass_level = ifelse(is.na(hydro_pass_level), TRUE, hydro_pass_level),
      
      IS_REAL = dplyr::case_when(
        has_alignment_gap == TRUE ~ FALSE,
        all_mutations_ok == TRUE & hydro_pass_level == TRUE ~ TRUE,
        TRUE ~ FALSE
      ),
      REASON = dplyr::case_when(
        has_alignment_gap == TRUE ~ "Excluded: Query causes gaps in alignment (aa_pos = -1)",
        all_mutations_ok == FALSE ~ "At least one mutation position is bio-chemically implausible",
        include_hydro & hydro_pass_level == FALSE ~ "Global hydrophobicity failed at this level",
        TRUE ~ "Verified (Best-of-Worst hierarchical check)"
      )
    ) %>%
    dplyr::select(ASV_id, tax_lev, IS_REAL, REASON, max_of_best_grantham, min_of_best_quality)
  
  return(final_report)
}