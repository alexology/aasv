get_hydro_metrics <- function(aa_seq, hydro_window) {
  aa_str <- as.character(aa_seq)
  n <- nchar(aa_str)
  
  # 1. Global Score (The average for the whole 140 AA)
  global_score <- Peptides::hydrophobicity(aa_str, scale = "KyteDoolittle")
  
  # 2. Local Profile (Sliding hydro_window to find 'disasters')
  # We use a hydro_window of 7, typical for alpha-helix spans
  indices <- seq_len(n - hydro_window + 1)
  local_profile <- sapply(indices, function(i) {
    sub_seq <- substr(aa_str, i, i + hydro_window - 1)
    Peptides::hydrophobicity(sub_seq, scale = "KyteDoolittle")
  })
  
  return(list(global = global_score, profile = local_profile))
}


# with the help of google gemini
check_asv_structure <- function(query_aa,
                                global_bounds,
                                lower_env,
                                upper_env,
                                violation_threshold,
                                hydro_window = 7,
                                asv_samples = NULL) {

  q_metrics <- get_hydro_metrics(aa_seq = query_aa,
                                 hydro_window = hydro_window)
  
  # Check A: Is Global score within bounds?
  global_pass <- q_metrics$global >= global_bounds[1] & 
    q_metrics$global <= global_bounds[2]
  
  # Check B: How many hydro_window positions break the local envelope?
  # We look for 'local disasters'
  violations <- which(q_metrics$profile < lower_env | 
                        q_metrics$profile > upper_env)
  
  violation_rate <- length(violations) / length(q_metrics$profile)
  
  return(data.frame(
    global_score = q_metrics$global,
    global_pass = global_pass,
    local_violation_rate = violation_rate,
    # High violation rate (>10%) suggests a non-functional sequence
    structural_verdict = ifelse(global_pass & violation_rate < violation_threshold, "PASS", "FAIL")
  ))
}


# from GEMINI
get_mut_pos <- function(q, r) {
  if (is.na(q) | is.na(r)) return(NA)
  q_v <- strsplit(as.character(q), "")[[1]]
  r_v <- strsplit(as.character(r), "")[[1]]
  
  # Return as a string for easy viewing in a table (e.g., "1" or "1,2")
  return(paste(which(q_v != r_v), collapse = ","))
}


