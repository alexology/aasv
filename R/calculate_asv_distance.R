#' @importFrom methods as
#' @importFrom stats median quantile
#' @importFrom utils write.table

calculate_asv_distance <- function(output_dir,
                                   vsearch_path,
                                   query_dna, 
                                   query_id,
                                   ref_aa_aligned,   # The aligned AA string from your 'alignments' folder
                                   species_name, 
                                   family_name, 
                                   ref_nt,
                                   asv_samples,
                                   genetic_code,
                                   hydro_threshold,
                                   hydro_window,
                                   trimmed_folder,
                                   ee,
                                   tax_lev,
                                   ...) {
  
  ref_aa_aligned <- unique(ref_aa_aligned)

  # Sentinel result for a query that cannot be reconciled with a functional,
  # indel-free COI ORF (premature stop codon, frameshift, or an alignment-
  # forcing indel found further down). aa_pos = -1 is the signal classify_asv()
  # reads as "Severe incompatibility with functional COI."
  no_orf_result <- function() {
    hydro_results <- data.frame(global_min = numeric(),
                                global_max = numeric(),
                                global_score = numeric(),
                                global_pass = logical(),
                                local_violation_rate = numeric(),
                                structural_verdict = character())

    results_ref <- data.frame(aa_pos = -1,
                              illegal_aa = NA,
                              dna_triplet = NA,
                              nt_coords = NA,
                              n_sequences_ali = NA,
                              n_unique_seq_ali = NA,
                              median_quality_prob = NA,
                              n_sequences_prob = NA,
                              ref_aa = NA,
                              dna_triplet_ref = NA,
                              nt_coords_ref = NA,
                              ref_id = NA,
                              triplet_mut_pos = NA,
                              ham_dist = NA,
                              grantham_dist = NA)

    list(hydro = hydro_results, aa = results_ref)
  }

  # 1. Translate Query
  q_dna_obj <- Biostrings::DNAStringSet(query_dna)

  names(q_dna_obj) <- query_id


  q_dna_start_pos <- best_translation(q_dna_obj,
                                      genetic_code = genetic_code,
                                      res_type = "position")

  q_aa_chr <- best_translation(q_dna_obj, genetic_code = genetic_code)

  # None of the 3 reading frames translates the query end-to-end without a
  # stop codon: a genuine premature stop (e.g. a NUMT or a mutation), not
  # merely a short/incomplete fragment. Report it as a hard ORF
  # incompatibility instead of aligning whatever stop-free tail remains.
  if (is.na(q_aa_chr)) return(no_orf_result())

  q_aa_raw <- Biostrings::AAStringSet(q_aa_chr)

  names(q_aa_raw) <- query_id

  # 2. Align Query AA to the Aligned Reference AA
  # We use global-local to find where the query fits in the pre-aligned ref
  aln <- DECIPHER::AlignProfiles(q_aa_raw,
                                 ref_aa_aligned,
                                 ...)
  
  # create the path to save the alingments
  aln_save <- file.path(output_dir,
                        "query_alignments",
                        paste0(query_id, "_", tax_lev, ".fasta"))
  
  
  # write to disk
  Biostrings::writeXStringSet(aln, aln_save)
  
  # check if the query has inner gaps (a deletion relative to the reference
  # profile) OR forces gap-opening in most reference sequences (an insertion
  # not shared by any reference) - either case is a structural incompatibility
  # with a functional, indel-free COI ORF, so end the function prematurely
  aln_mat       <- as.matrix(aln)
  query_row     <- which(names(aln) == query_id)
  query_has_ig  <- inner_gaps(aln_mat)[query_row] ||
                     causes_internal_gaps(aln_mat)[query_row]

  if(query_has_ig) return(no_orf_result())
  
  
  # get the number of gaps for each sequence
  # do this for very messy alignment
  num_gaps_seq <- Biostrings::letterFrequency(aln, "-")[,1]
  
  #######
  
  # 2. Create the Consensus Matrix
  # This results in a matrix where rows are AAs (A, R, N, D...) and columns are positions
  prof_matrix <- Biostrings::consensusMatrix(aln[which(names(aln) != query_id)])
  
  # 3. Convert counts to proportions (Probabilities)
  # This gives you a profile where each column sums to 1
  prof_prop <- prop.table(prof_matrix, margin = 2)
  
  # 4. Filter for standard Amino Acids (removing gaps '-' and 'X')
  standard_aa <- row.names(prof_prop) %in% Biostrings::AA_STANDARD
  prof_final <- prof_prop[standard_aa, ]
  
  # 3. Get the Offset (How many leading dashes in the reference alignment?)
  # This tells us where the 'biological' index 1 starts
  ref_aln_str <- as.character(aln)
  leading_dashes <- nchar(gsub("^([-]*).*", "\\1", ref_aln_str))
  
  # Assuming 'aln_query' is your ASV aligned to the reference profile
  q_vec <- strsplit(as.character(aln[which(names(aln) %in% query_id)]), "") %>%
    unlist()
  
  # Logic to flag "Forbidden" Amino Acids
  audit_results <- lapply(seq_along(q_vec), function(i) {
    aa_q <- q_vec[i]
    
    if(aa_q == "-" | aa_q == "*") return(NULL) # Skip gaps/stops
    
    # Look up the frequency of this AA at this position in the reference database
    # If the AA is missing from the matrix, frequency is 0
    freq <- ifelse(aa_q %in% rownames(prof_final), prof_final[aa_q, i], 0)
    
    return(data.frame(
      pos = i,
      query_aa = aa_q,
      db_freq = freq,
      is_illegal = freq == 0
    ))
  })
  
  # subset the list
  leading_dashes_query <- leading_dashes[names(leading_dashes) == query_id]  
  leading_dashes_ref <- leading_dashes[names(leading_dashes) != query_id]  
  
  
  ### hydrophobicity module ----------------------------------------------------
  aa_mat <- aln_mat
  
  # find the position of the ending dash in the query sequences
  match_data <- regexpr("-+$", as.character(aln[which(rownames(aa_mat) %in% query_id)]))
  
  # Extract the position
  # match_data returns the starting index of the match
  first_dash_pos <- as.integer(match_data)
  
  # Handle cases with no dashes (regexpr returns -1 if no match). final_pos is
  # used below as an exclusive upper bound (final_pos - 1), so the no-dash case
  # must be ncol(aa_mat) + 1 or the alignment's last column is silently dropped.
  if(first_dash_pos == -1) {
    final_pos <- ncol(aa_mat) + 1
  } else {
    final_pos <- first_dash_pos
  }
  
  
  aa_mat_q_gaps <- aa_mat[, c((leading_dashes_query+1):(final_pos-1)), drop = FALSE]

  aa_query_q_gaps <- aa_mat_q_gaps[rownames(aa_mat_q_gaps) == query_id,  , drop = FALSE]
  aa_mat_q_gaps   <- aa_mat_q_gaps[rownames(aa_mat_q_gaps) != query_id,  , drop = FALSE]

  # Drop columns where the majority of references have a gap: these positions lie
  # beyond the end (or before the start) of most reference sequences and must not
  # inflate or deflate the query's hydrophobicity score simply because the query
  # is longer (or starts earlier) than the typical reference.
  if (nrow(aa_mat_q_gaps) > 0) {
    shared_cols     <- colMeans(aa_mat_q_gaps == "-") < 0.5
    aa_query_q_gaps <- aa_query_q_gaps[, shared_cols, drop = FALSE]
    aa_mat_q_gaps   <- aa_mat_q_gaps[,   shared_cols, drop = FALSE]
  }

  n_gaps_cover  <- apply(aa_mat_q_gaps, 1, function(x) sum(x == "-"))
  aa_mat_q_gaps <- aa_mat_q_gaps[n_gaps_cover == 0, , drop = FALSE]
  
  if(nrow(aa_mat_q_gaps) == 0 || ncol(aa_mat_q_gaps) < hydro_window){
    hydro_results <- data.frame(global_min = numeric(),
                                global_max = numeric(),
                                global_score = numeric(),
                                global_pass = logical(),
                                local_violation_rate = numeric(),
                                structural_verdict = character())
    
  } else {
    # gap_free_cols <- which(colSums(aa_mat == "-") == 0)
    # aa_mat <- aa_mat[, gap_free_cols]
    
    aa_mat_ref <- aa_mat_q_gaps
    hydro_query <- aa_query_q_gaps %>%
      unlist() %>%
      paste0(collapse = "")
    
    aa_clean <- Biostrings::AAStringSet(apply(aa_mat_ref, 1, paste, collapse = ""))
    
    # restore the original names
    names(aa_clean) <- names(aa_mat_ref)
    
    hydro_ref <- aa_clean %>%
      as.character()
    
    # 1. Calculate profiles for all references
    # Assumes ref_aa_set is your AAStringSet of BOLD references
    ref_data <- lapply(hydro_ref,
                       function(x) get_hydro_metrics(x, hydro_window = hydro_window))
    
    # 2. Extract global scores and calculate the 95% interval
    ref_globals <- sapply(ref_data, function(x) x$global)
    global_bounds <- quantile(ref_globals, probs = c(0, 1))
    
    ref_profile_mat <- do.call(cbind, lapply(ref_data, function(x) x$profile))
    lower_envelope <- apply(ref_profile_mat, 1, quantile, 0.025)
    upper_envelope <- apply(ref_profile_mat, 1, quantile, 0.975)
    
    
    hydro_results <- check_asv_structure(hydro_query,
                                         global_bounds,
                                         lower_envelope,
                                         upper_envelope,
                                         violation_threshold = hydro_threshold,
                                         hydro_window = hydro_window)
    
    
    hydro_results <- data.frame(global_min = global_bounds[1],
                                global_max = global_bounds[2],
                                hydro_results)
    
    rownames(hydro_results) <- NULL
  }
  
  
  
  ###---------------------------------------------------------------------------
  
  
  # from list to data.frame, excluding leading dashes  
  audit_df <- do.call(rbind, audit_results[-seq_len(leading_dashes_query)])  
  
  # Filter only the illegal hits
  illegal_hits <- audit_df[audit_df$is_illegal,]
  
  if(nrow(illegal_hits) == 0){
    list(hydro = hydro_results) 
  } else{
    
    results <- lapply(seq_len(nrow(illegal_hits)), function(j) {
      row <- illegal_hits[j, ]
      aa_pos <- row$pos
      
      # Calculate DNA coordinates
      nt_start <- (((aa_pos - leading_dashes_query) * 3) - 2) + (q_dna_start_pos -1)
      nt_end   <- nt_start + 2
      
      # Extract the triplet from the raw DNA sequence
      triplet <- substring(query_dna, first = nt_start, last = nt_end)
      
      check_translation <- best_translation(triplet, genetic_code = genetic_code)
      
      if(illegal_hits$query_aa[j] != check_translation){
        stop("Lost in translation")
      }
      
      data.frame(
        aa_pos = aa_pos,
        illegal_aa = row$query_aa,
        dna_triplet = triplet,
        nt_coords = paste0(nt_start, "-", nt_end)
      )
    })
    
    results <- dplyr::bind_rows(results)  
    
    # number of sequences in th reference alignment
    results$n_sequences_ali <- length(ref_aa_aligned)
    
    # number of sequences in th reference alignment
    results$n_unique_seq_ali <- length(unique(ref_aa_aligned))
    
    # map the results back to reference dataset
    ref_nt_sub <- ref_nt[which(names(ref_nt) %in% names(ref_aa_aligned))]
    starting_postion_ref <- sapply(seq_along(ref_nt_sub),
                                   function(x) best_translation(ref_nt_sub[x],
                                                                 genetic_code = genetic_code,
                                                                 res_type = "position"))
    names(starting_postion_ref) <- names(ref_nt_sub)
    
    results_ref <- list()
    for (j in seq_len(nrow(results))) {
      results_j_pos <- results$aa_pos[j]

      results_ref_j <- lapply(seq_along(ref_nt_sub), function(x) {
        ref_seq_name <- names(ref_nt_sub)[x]

        # Use the reference row from the combined alignment so positions are
        # consistent with results_j_pos (which comes from aln, not ref_aa_aligned)
        ref_aln_seq <- as.character(aln[which(names(aln) == ref_seq_name)])
        aa_ref_x    <- substring(ref_aln_seq, first = results_j_pos, last = results_j_pos)

        # Skip references that have a gap at this position (cross-species insertion
        # column present in the genus alignment but absent in this reference)
        if (aa_ref_x == "-") return(NA)

        # Count non-gap amino acids up to results_j_pos to obtain the codon index.
        # Using nchar(gsub) rather than (pos - leading_dashes) makes the calculation
        # correct even when AlignProfiles introduces internal gap columns.
        k          <- nchar(gsub("-", "", substring(ref_aln_seq, 1, results_j_pos)))
        nt_start_r <- (k * 3 - 2) + (starting_postion_ref[x] - 1)
        nt_end_r   <- nt_start_r + 2

        triplet_ref <- substring(as.character(ref_nt_sub[x]), first = nt_start_r, last = nt_end_r)

        if (nchar(triplet_ref) != 3) {
          NA
        } else {
          check_translation_ref <- aa_ref_x == best_translation(triplet_ref,
                                                                 genetic_code = genetic_code)

          if (!check_translation_ref) {
            stop("Lost in translation.")
          }

          data.frame(
            aa_pos = results_j_pos,
            ref_aa = aa_ref_x,
            dna_triplet_ref = triplet_ref,
            nt_coords_ref = paste0(nt_start_r, "-", nt_end_r),
            ref_id = names(triplet_ref))
        }
      })
      
      results_ref_j <- results_ref_j[sapply(results_ref_j, is.data.frame)]
      results_ref <- c(results_ref, list(dplyr::bind_rows(results_ref_j)))   
    }
    
    results_ref <- dplyr::bind_rows(results_ref)
    
    rownames(results_ref) <- NULL
    
    if(nrow(results_ref) == 0){
      return(list(hydro = hydro_results,
                  aa = results_ref) )
    }
    
    # default: quality columns absent when VSEARCH is not used
    results$median_quality_prob <- NA_real_
    results$n_sequences_prob    <- NA_integer_

    if (!is.null(vsearch_path) && !is.null(trimmed_folder)) {

      fastq_to_check <- list.files(trimmed_folder,
                                   full.names = TRUE,
                                   pattern = paste0(asv_samples, collapse = "|"))

      dir.create(file.path(output_dir, "vsearch_map_back"), showWarnings = FALSE)

      query_fasta     <- file.path(output_dir, "vsearch_map_back",
                                   paste0(query_id, ".fasta"))
      query_fasta_bio <- Biostrings::DNAStringSet(query_dna)
      names(query_fasta_bio) <- query_id
      Biostrings::writeXStringSet(query_fasta_bio, filepath = query_fasta)

      fastq_sub <- character(0)

      for (j in seq_along(fastq_to_check)) {
        output_tab   <- file.path(output_dir, "vsearch_map_back",
                                  paste0(query_id, "_", j, "_mapping_results.tsv"))
        output_fastq <- file.path(output_dir, "vsearch_map_back",
                                  paste0(query_id, "_", j, "_mapping_results.fastq"))
        fastq_sub <- c(fastq_sub, output_fastq)

        system2(vsearch_path,
                args = c("--usearch_global", shQuote(query_fasta),
                         "--db",             shQuote(fastq_to_check[j]),
                         "--id",             "1.0",
                         "--strand",         "plus",
                         "--maxaccepts",     "0",
                         "--maxrejects",     "0",
                         "--userout",        shQuote(output_tab),
                         "--userfields",     "target",
                         "--notrunclabels"),
                stdout = TRUE)

        if (!file.exists(output_tab) || length(readLines(output_tab)) == 0) next()

        read_ids     <- trimws(unique(unlist(read.table(output_tab, sep = "\t",
                                                         stringsAsFactors = FALSE,
                                                         header = FALSE))))
        output_lines <- file.path(output_dir, "vsearch_map_back",
                                  paste0(query_id, "_", j, "_clean_ids.txt"))
        con <- file(output_lines, "wb")
        writeLines(read_ids, con, sep = "\n")
        close(con)

        system2(vsearch_path,
                args = c("--fastx_getseqs", shQuote(fastq_to_check[j]),
                         "--label_words",   shQuote(output_lines),
                         "--fastqout",      shQuote(output_fastq),
                         "--notrunclabels"),
                stdout = FALSE, stderr = FALSE)
      }

      fastq_sub <- fastq_sub[file.exists(fastq_sub)]

      if (length(fastq_sub) > 0) {
        pad_cols <- function(m, k) {
          if (ncol(m) < k) cbind(m, matrix(NA, nrow(m), k - ncol(m))) else m
        }

        fq_list <- lapply(fastq_sub, ShortRead::readFastq) %>%
          lapply(Biostrings::quality) %>%
          lapply(Biostrings::PhredQuality) %>%
          lapply(as, "NumericList") %>%
          lapply(as.matrix)

        max_ncol <- max(sapply(fq_list, ncol))
        fq_list  <- do.call(rbind, lapply(fq_list, pad_cols, k = max_ncol))
        fq_list  <- fq_list[apply(fq_list, 1, sum, na.rm = TRUE) <= ee, , drop = FALSE]

        triplet_to_check      <- do.call(rbind, strsplit(results$nt_coords, "-")) %>%
          as.data.frame()
        triplet_to_check$V1   <- as.numeric(triplet_to_check$V1)
        triplet_to_check$V2   <- as.numeric(triplet_to_check$V2)

        triplet_quality <- apply(triplet_to_check, 1,
                                 function(x) fq_list[, x[1]:x[2], drop = FALSE],
                                 simplify = FALSE)

        results$median_quality_prob <- sapply(
          lapply(triplet_quality, function(x) 1 - rowSums(x)), median)
        results$n_sequences_prob    <- nrow(fq_list)
      }
    }

    # update the results
    results_ref <- tibble::as_tibble(dplyr::inner_join(results, results_ref, by = "aa_pos"))
    
    
    # calculate the mutation position
    results_ref$triplet_mut_pos <- mapply(get_mut_pos,
                                          results_ref$dna_triplet,
                                          results_ref$dna_triplet_ref)
    
    # calculate the hamming distance (stringdist)
    results_ref$ham_dist <- stringdist::stringdist(results_ref$dna_triplet, 
                                                   results_ref$dna_triplet_ref, 
                                                   method = "hamming")
    
    
    results_ref$grantham_dist <- mapply(grantham::grantham_distance,
                                        grantham::as_three_letter(results_ref$illegal_aa), 
                                        grantham::as_three_letter(results_ref$ref_aa)) %>%
      t() %>%
      as.data.frame() %>%
      dplyr::pull("d") %>%
      unlist()
    
    
    list(hydro = hydro_results,
         aa = results_ref) 
    
  }
}

