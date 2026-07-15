# Suppress R CMD check NOTEs for variables used inside dplyr/tidyr NSE calls.
# These are column names referenced without quotes inside mutate(), filter(),
# arrange(), etc.; they are not unbound global variables at runtime.
utils::globalVariables(c(
  # magrittr dot placeholder
  ".",
  # bold / align_species_seq column names
  "marker", "identification", "sequence", "processid",
  # aa_similarity_matrix intermediate names
  "name", "ID", "value", "triplets_1", "triplets_2",
  # asv_functional_structure column names
  "species", "ASV_id",
  # classify_asv column names
  "tax_lev", "aa_pos", "grantham_dist", "median_quality_prob",
  "triplet_mut_pos", "is_valid_ref", "mutation_ok",
  "best_grantham_at_pos", "best_quality_at_pos",
  "all_mutations_ok", "has_alignment_gap",
  "global_pass", "hydro_pass_level",
  "IS_REAL", "REASON",
  "max_of_best_grantham", "min_of_best_quality",
  # classify_asv Functionality Index column names
  "illegal_aa", "ref_aa", "is_reliable", "n_sequences_ali",
  "n_ref_sequences", "conservation_weight", "conservation_confidence",
  "normalized_grantham", "site_penalty", "is_pos1", "is_pos2",
  "local_violation_rate", "mean_site_penalty", "mean_conservation_weight",
  "mean_conservation_confidence", "mean_normalized_grantham",
  "n_pos1", "n_pos2", "codon_penalty", "normalized_hydro",
  "total_penalty", "FI", "Class", "Evidence",
  "n_aa_substitutions", "n_flagged_aa_substitutions",
  # align_species_seq deduplication index column
  "id",
  # classify_asv / grantham result column
  "d"
))
