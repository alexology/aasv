skip_if_not_installed("writexl")
skip_if_not_installed("readxl")

# Helper: create a minimal aa_structure_results directory
make_classify_dir <- function() {
  tmpdir <- tempfile(pattern = "aasv_classify_")
  dir.create(file.path(tmpdir, "aa_structure_results"), recursive = TRUE)
  tmpdir
}

# Helper: write a minimal hydrophobicity Excel file. `local_violation_rate`
# feeds into the FI's independent, sequence-level hydrophobicity term;
# `global_pass` is kept only because classify_asv() still reads the column
# (it is not used for scoring). This file also doubles as the authoritative
# list of every ASV/tax_lev combination that was processed upstream, so
# classify_asv() can report ASVs with no mutated sites at all.
write_hydro <- function(dir, asv_id = "ASV001", tax_lev = "species",
                        global_pass = TRUE, local_violation_rate = 0.0) {
  write_hydro_multi(dir, rows = list(
    list(asv_id = asv_id, tax_lev = tax_lev, global_pass = global_pass,
         local_violation_rate = local_violation_rate)
  ))
}

# Helper: write a hydrophobicity Excel file covering several ASV/tax_lev rows
# at once, to simulate a batch where some ASVs have no mutated sites (and so
# never get an <ASV_id>_aa_structure.xlsx file) alongside ones that do.
write_hydro_multi <- function(dir, rows) {
  df <- dplyr::bind_rows(lapply(rows, function(r) {
    data.frame(
      ASV_id               = r$asv_id,
      tax_lev              = r$tax_lev,
      global_min           = -1,
      global_max           = 1,
      global_score         = 0.5,
      global_pass          = r$global_pass %||% TRUE,
      local_violation_rate = r$local_violation_rate %||% 0.0,
      structural_verdict   = "PASS"
    )
  }))
  writexl::write_xlsx(df, file.path(dir, "aa_structure_results",
                                    "hydrophobicity.xlsx"))
}

# Helper: write an ASV structure Excel file with a given aa_pos
write_struct <- function(dir, asv_id = "ASV001", tax_lev = "species",
                         aa_pos = 5L,
                         grantham_dist = 20, median_quality_prob = 0.9999,
                         triplet_mut_pos = "3",
                         ref_aa = "L", illegal_aa = "F",
                         n_sequences_ali = NA_integer_) {
  df <- data.frame(
    ASV_id               = asv_id,
    tax_lev              = tax_lev,
    aa_pos               = aa_pos,
    illegal_aa           = if (aa_pos == -1L) NA_character_ else illegal_aa,
    dna_triplet          = NA_character_,
    nt_coords            = NA_character_,
    n_sequences_ali      = n_sequences_ali,
    n_unique_seq_ali     = NA_integer_,
    median_quality_prob  = median_quality_prob,
    n_sequences_prob     = NA_integer_,
    ref_aa               = if (aa_pos == -1L) NA_character_ else ref_aa,
    dna_triplet_ref      = NA_character_,
    nt_coords_ref        = NA_character_,
    ref_id               = NA_character_,
    triplet_mut_pos      = triplet_mut_pos,
    ham_dist             = 1L,
    grantham_dist        = grantham_dist
  )
  writexl::write_xlsx(df, file.path(dir, "aa_structure_results",
                                    paste0(asv_id, "_aa_structure.xlsx")))
}

# Helper: write an ASV structure Excel file with several mutation rows at once.
# Each element of `rows` becomes one row; use several rows sharing the same
# `aa_pos` (with different `ref_aa`) to simulate multiple reference sequences
# disagreeing at a position (a "variable" site for the conservation term).
write_struct_multi <- function(dir, asv_id = "ASV001", tax_lev = "species",
                               rows) {
  df <- dplyr::bind_rows(lapply(rows, function(r) {
    data.frame(
      ASV_id               = asv_id,
      tax_lev              = tax_lev,
      aa_pos               = r$aa_pos,
      illegal_aa           = r$illegal_aa,
      dna_triplet          = NA_character_,
      nt_coords            = NA_character_,
      n_sequences_ali      = r$n_sequences_ali %||% NA_integer_,
      n_unique_seq_ali     = NA_integer_,
      median_quality_prob  = r$median_quality_prob %||% 0.9999,
      n_sequences_prob     = NA_integer_,
      ref_aa               = r$ref_aa,
      dna_triplet_ref      = NA_character_,
      nt_coords_ref        = NA_character_,
      ref_id               = NA_character_,
      triplet_mut_pos      = r$triplet_mut_pos %||% "3",
      ham_dist             = 1L,
      grantham_dist        = r$grantham_dist
    )
  }))
  writexl::write_xlsx(df, file.path(dir, "aa_structure_results",
                                    paste0(asv_id, "_aa_structure.xlsx")))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- output structure --------------------------------------------------------

test_that("classify_asv returns a data.frame with required columns", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir)

  result <- classify_asv(tmpdir)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("ASV_id", "tax_lev", "FI", "Class",
                     "mean_conservation_weight", "mean_conservation_confidence",
                     "mean_normalized_grantham",
                     "normalized_hydro", "mean_site_penalty", "n_pos1", "n_pos2",
                     "codon_penalty", "total_penalty", "Evidence",
                     "n_aa_substitutions", "n_flagged_aa_substitutions",
                     "n_ref_sequences") %in% names(result)))
})

# --- reporting every processed ASV, not only the flagged ones ---------------

test_that("an ASV with no mutated sites at any level is still reported (no aa_structure file at all)", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  # ASV001 has a flagged mutation and gets an aa_structure file; ASV002 has a
  # clean ORF at this level, so calculate_asv_distance() upstream never wrote
  # one for it - only the hydrophobicity file records that it was processed
  write_hydro_multi(tmpdir, rows = list(
    list(asv_id = "ASV001", tax_lev = "species", local_violation_rate = 0),
    list(asv_id = "ASV002", tax_lev = "species", local_violation_rate = 0)
  ))
  write_struct(tmpdir, asv_id = "ASV001", aa_pos = 5L, grantham_dist = 20,
               triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)

  expect_true("ASV002" %in% result$ASV_id)
  row2 <- result[result$ASV_id == "ASV002", ]
  expect_equal(row2$n_aa_substitutions, 0L)
  expect_equal(row2$FI, 1)
  expect_equal(row2$Class, "Plausible functional sequence")
})

test_that("clean ASVs are reported per taxonomic level alongside flagged ones", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro_multi(tmpdir, rows = list(
    list(asv_id = "ASV001", tax_lev = "species", local_violation_rate = 0),
    list(asv_id = "ASV001", tax_lev = "genus",   local_violation_rate = 0),
    list(asv_id = "ASV002", tax_lev = "species", local_violation_rate = 0),
    list(asv_id = "ASV002", tax_lev = "genus",   local_violation_rate = 0)
  ))
  write_struct(tmpdir, asv_id = "ASV001", tax_lev = "species",
               aa_pos = 5L, grantham_dist = 200, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)

  # ASV002 never had an aa_structure file at all, but both of its levels
  # (present in the hydrophobicity file) must still be listed
  expect_equal(nrow(result[result$ASV_id == "ASV002", ]), 2L)
  expect_setequal(result$tax_lev[result$ASV_id == "ASV002"], c("species", "genus"))
})

test_that("the reported submetrics reproduce FI by hand", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0.2)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]

  expect_equal(row$mean_normalized_grantham, round(20 / 215, 3), tolerance = 0.001)
  expect_equal(row$normalized_hydro, 0.2, tolerance = 0.001)
  expect_equal(row$mean_conservation_weight, 1, tolerance = 0.001)
  expect_equal(row$n_pos1, 0L)
  expect_equal(row$n_pos2, 1L)
  expect_equal(row$codon_penalty, 1, tolerance = 0.001)

  # recompute from the raw (unrounded) inputs, not from the already-rounded
  # display columns, to avoid double-rounding ambiguity. site_penalty is
  # conservation x Grantham only; hydrophobicity is an independent
  # sequence-level term rather than being folded into the per-site penalty.
  expected_site_penalty <- 1 * (20 / 215)
  expected_total <- 0.35 * expected_site_penalty + 0.3 * 1 + 0.35 * 0.2
  expect_equal(row$mean_site_penalty, expected_site_penalty, tolerance = 0.002)
  expect_equal(row$total_penalty, expected_total, tolerance = 0.002)
  expect_equal(row$FI, 1 - expected_total, tolerance = 0.002)
})

test_that("submetrics are NA/0 when there are no mutated sites, but hydrophobicity still contributes", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0.3)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]

  expect_true(is.na(row$mean_conservation_weight))
  expect_true(is.na(row$mean_normalized_grantham))
  expect_equal(row$normalized_hydro, 0.3, tolerance = 0.001)
  expect_equal(row$mean_site_penalty, 0)
  expect_equal(row$codon_penalty, 0)
  # hydrophobicity is an independent sequence-level term, so it still
  # contributes to the total penalty even with no mutated amino-acid sites
  expect_equal(row$total_penalty, 0.35 * 0.3, tolerance = 0.001)
  expect_equal(row$FI, 1 - 0.35 * 0.3, tolerance = 0.001)
})

test_that("submetrics are all NA for a Severe incompatibility (alignment gap) row", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir, aa_pos = -1L)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]

  expect_true(is.na(row$mean_conservation_weight))
  expect_true(is.na(row$mean_normalized_grantham))
})

test_that("classify_asv errors when hydrophobicity file is missing and include_hydro = TRUE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_error(classify_asv(tmpdir, include_hydro = TRUE),
               "Hydrophobicity file not found")
})

test_that("classify_asv returns empty data.frame when no files are present", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  result <- classify_asv(tmpdir, include_hydro = FALSE)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_true("Class" %in% names(result))
})

# --- hard incompatibility (alignment gap) ------------------------------------

test_that("classify_asv assigns FI = 0 and 'Severe incompatibility' when an alignment gap is present", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir, aa_pos = -1L)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 0)
  expect_equal(row$Class, "Severe incompatibility with functional COI.")
  expect_equal(row$n_aa_substitutions, 0L)
  expect_equal(row$n_flagged_aa_substitutions, 0L)
  expect_true(is.na(row$n_ref_sequences))
})

# --- FI formula: codon position ----------------------------------------------

test_that("a synonymous (3rd position) mutation with low Grantham and no hydro shift scores Plausible", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20,
               median_quality_prob = 0.9999, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 0.967, tolerance = 0.001)
  expect_equal(row$Class, "Plausible functional sequence")
})

test_that("second codon position mutations are penalised twice as heavily as first-position ones", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))
  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "1")
  fi_pos1 <- classify_asv(tmpdir)$FI[1]

  tmpdir2 <- make_classify_dir()
  on.exit(unlink(tmpdir2, recursive = TRUE), add = TRUE)
  write_hydro(tmpdir2, local_violation_rate = 0)
  write_struct(tmpdir2, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "2")
  fi_pos2 <- classify_asv(tmpdir2)$FI[1]

  tmpdir3 <- make_classify_dir()
  on.exit(unlink(tmpdir3, recursive = TRUE), add = TRUE)
  write_hydro(tmpdir3, local_violation_rate = 0)
  write_struct(tmpdir3, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")
  fi_pos3 <- classify_asv(tmpdir3)$FI[1]

  expect_equal(fi_pos1, 0.817, tolerance = 0.001)
  expect_equal(fi_pos2, 0.667, tolerance = 0.001)
  expect_equal(fi_pos3, 0.967, tolerance = 0.001)
  expect_true(fi_pos3 > fi_pos1)
  expect_true(fi_pos1 > fi_pos2)
})

test_that("include_position = FALSE removes the codon-position penalty", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir, include_position = FALSE)
  expect_equal(result$FI[result$ASV_id == "ASV001"], 0.967, tolerance = 0.001)
})

# --- FI formula: Grantham distance -------------------------------------------

test_that("high Grantham distance lowers FI", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 200, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 0.374, tolerance = 0.001)
  expect_equal(row$Class, "Artifact-NUMTs candidate")
})

test_that("include_grantham = FALSE removes the Grantham contribution", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 200, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_grantham = FALSE)
  expect_equal(result$FI[result$ASV_id == "ASV001"], 1, tolerance = 0.001)
  expect_equal(result$Class[result$ASV_id == "ASV001"],
               "Plausible functional sequence")
})

# --- FI formula: hydrophobicity shift -----------------------------------------

test_that("a large local hydrophobicity violation rate lowers FI", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 1.0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 0.617, tolerance = 0.001)
  expect_equal(row$Class, "Artifact-NUMTs candidate")
})

test_that("include_hydro = FALSE removes the hydrophobicity contribution", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")
  # no hydrophobicity file is written at all - include_hydro = FALSE must not
  # require it
  result <- classify_asv(tmpdir, include_hydro = FALSE)
  expect_equal(result$FI[result$ASV_id == "ASV001"], 0.967, tolerance = 0.001)
})

# --- FI formula: conservation weighting --------------------------------------

test_that("a mutation at a variable reference position is penalised less than at a conserved one", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  # four references disagreeing at the same position (L, L, M, V) -> variable
  write_struct_multi(tmpdir, rows = list(
    list(aa_pos = 5L, illegal_aa = "F", ref_aa = "L", grantham_dist = 20, triplet_mut_pos = "3"),
    list(aa_pos = 5L, illegal_aa = "F", ref_aa = "L", grantham_dist = 20, triplet_mut_pos = "3"),
    list(aa_pos = 5L, illegal_aa = "F", ref_aa = "M", grantham_dist = 20, triplet_mut_pos = "3"),
    list(aa_pos = 5L, illegal_aa = "F", ref_aa = "V", grantham_dist = 20, triplet_mut_pos = "3")
  ))

  result_with_cons <- classify_asv(tmpdir, include_conservation = TRUE)
  result_no_cons    <- classify_asv(tmpdir, include_conservation = FALSE)

  fi_with_cons <- result_with_cons$FI[result_with_cons$ASV_id == "ASV001"]
  fi_no_cons   <- result_no_cons$FI[result_no_cons$ASV_id == "ASV001"]

  expect_equal(fi_with_cons, 0.979, tolerance = 0.001)
  expect_equal(fi_no_cons,   0.967, tolerance = 0.001)
  expect_true(fi_with_cons > fi_no_cons)
})

# --- FI formula: conservation confidence (sparse reference sets) ------------

test_that("a conserved position backed by few reference sequences is down-weighted", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3",
               n_sequences_ali = 2L)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]

  # raw entropy-based weight is 1 (single distinct ref aa observed), but the
  # confidence factor 1 - exp(-n_ref / conservation_k) heavily discounts it
  # when only two reference sequences back the call
  expected_weight <- 1 - exp(-2 / 20)
  expect_equal(row$mean_conservation_weight, expected_weight, tolerance = 0.001)
  expect_true(row$mean_conservation_weight < 0.1)
})

test_that("the same apparent conservation is trusted more with many reference sequences", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3",
               n_sequences_ali = 200L)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]

  expected_weight <- 1 - exp(-200 / 20)
  expect_equal(row$mean_conservation_weight, expected_weight, tolerance = 0.001)
  expect_true(row$mean_conservation_weight > 0.99)
})

test_that("conservation confidence down-weighting raises FI for the same mutation backed by fewer references", {
  tmpdir_sparse <- make_classify_dir()
  on.exit(unlink(tmpdir_sparse, recursive = TRUE))
  write_hydro(tmpdir_sparse, local_violation_rate = 0)
  write_struct(tmpdir_sparse, aa_pos = 5L, grantham_dist = 200, triplet_mut_pos = "2",
               n_sequences_ali = 2L)
  fi_sparse <- classify_asv(tmpdir_sparse)$FI[1]

  tmpdir_dense <- make_classify_dir()
  on.exit(unlink(tmpdir_dense, recursive = TRUE), add = TRUE)
  write_hydro(tmpdir_dense, local_violation_rate = 0)
  write_struct(tmpdir_dense, aa_pos = 5L, grantham_dist = 200, triplet_mut_pos = "2",
               n_sequences_ali = 500L)
  fi_dense <- classify_asv(tmpdir_dense)$FI[1]

  # same substitution and Grantham distance, but the position backed by only
  # two references earns a smaller (less confident) conservation weight, so
  # its penalty is smaller and FI higher
  expect_true(fi_sparse > fi_dense)
})

test_that("conservation_k controls how quickly confidence grows with n_ref_sequences", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3",
               n_sequences_ali = 2L)

  weight_k20 <- classify_asv(tmpdir, conservation_k = 20)$mean_conservation_weight[1]
  weight_k2  <- classify_asv(tmpdir, conservation_k = 2)$mean_conservation_weight[1]

  # a smaller k reaches full confidence faster, so n_ref_sequences = 2 is
  # trusted more when conservation_k = 2 than when conservation_k = 20
  expect_true(weight_k2 > weight_k20)
})

test_that("conservation weight is left unadjusted when n_ref_sequences is unknown", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  expect_equal(result$mean_conservation_weight[result$ASV_id == "ASV001"], 1,
               tolerance = 0.001)
})

# --- quality as a reliability pre-filter --------------------------------------

test_that("mutations below the quality threshold are excluded from scoring entirely", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  # severe mutation, but low read-back quality
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 200,
               median_quality_prob = 0.5, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir, include_quality = TRUE)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 1, tolerance = 0.001)
  expect_equal(row$n_aa_substitutions, 0L)
  expect_equal(row$Class, "Plausible functional sequence")
})

test_that("include_quality = FALSE keeps low-quality mutations in the score", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 200,
               median_quality_prob = 0.5, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir, include_quality = FALSE)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$FI, 0.374, tolerance = 0.001)
  expect_equal(row$n_aa_substitutions, 1L)
})

test_that("NA quality (no VSEARCH read-back) is always treated as reliable", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20,
               median_quality_prob = NA_real_, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_quality = TRUE)
  expect_equal(result$n_aa_substitutions[result$ASV_id == "ASV001"], 1L)
  expect_equal(result$FI[result$ASV_id == "ASV001"], 0.967, tolerance = 0.001)
})

# --- Evidence content --------------------------------------------------------

test_that("Evidence contains amino acid substitution arrow notation", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3",
               ref_aa = "L", illegal_aa = "F")

  result <- classify_asv(tmpdir)
  expect_true(grepl("Leu.*Phe", result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence contains codon position label for triplet_mut_pos = '2'", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 30, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)
  expect_true(grepl("second codon position",
                    result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence mentions the hard incompatibility reason when an alignment gap is present", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir, aa_pos = -1L)

  result <- classify_asv(tmpdir)
  expect_true(grepl("premature stop codon, frameshift, or non-triplet indel",
                    result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence contains Grantham distance", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 98, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)
  expect_true(grepl("98", result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence omits Grantham when include_grantham = FALSE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_grantham = FALSE)
  expect_false(grepl("Grantham", result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence mentions conservation when include_conservation = TRUE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_conservation = TRUE)
  expect_true(grepl("conservation weight",
                    result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence omits conservation when include_conservation = FALSE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_conservation = FALSE)
  expect_false(grepl("conservation weight",
                     result$Evidence[result$ASV_id == "ASV001"]))
})

test_that("Evidence reports the FI and Class", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  expect_true(grepl("FI = 0.967", result$Evidence[result$ASV_id == "ASV001"]))
  expect_true(grepl("Plausible functional sequence",
                    result$Evidence[result$ASV_id == "ASV001"]))
})

# --- substitution counts and reference-sequence coverage ---------------------

test_that("classify_asv reports n_aa_substitutions and n_flagged_aa_substitutions", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct_multi(tmpdir, rows = list(
    list(aa_pos = 5L,  illegal_aa = "F", ref_aa = "L",
         grantham_dist = 0,   triplet_mut_pos = "3"),   # zero penalty
    list(aa_pos = 10L, illegal_aa = "W", ref_aa = "G",
         grantham_dist = 200, triplet_mut_pos = "2")    # flagged
  ))

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$n_aa_substitutions, 2L)
  expect_equal(row$n_flagged_aa_substitutions, 1L)
})

test_that("classify_asv reports 0 substitutions when none are found", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)

  result <- classify_asv(tmpdir)
  row <- result[result$ASV_id == "ASV001", ]
  expect_equal(row$n_aa_substitutions, 0L)
  expect_equal(row$n_flagged_aa_substitutions, 0L)
  expect_equal(row$FI, 1)
  expect_equal(row$Class, "Plausible functional sequence")
})

test_that("classify_asv reports n_ref_sequences from n_sequences_ali", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3",
               n_sequences_ali = 42L)

  result <- classify_asv(tmpdir)
  expect_equal(result$n_ref_sequences[result$ASV_id == "ASV001"], 42L)
})

test_that("classify_asv reports the weakest n_ref_sequences across mutation positions", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct_multi(tmpdir, rows = list(
    list(aa_pos = 5L,  illegal_aa = "F", ref_aa = "L",
         grantham_dist = 20, triplet_mut_pos = "3", n_sequences_ali = 50L),
    list(aa_pos = 10L, illegal_aa = "W", ref_aa = "G",
         grantham_dist = 20, triplet_mut_pos = "3", n_sequences_ali = 3L)
  ))

  result <- classify_asv(tmpdir)
  expect_equal(result$n_ref_sequences[result$ASV_id == "ASV001"], 3L)
})

test_that("classify_asv reports NA n_ref_sequences when unavailable", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, local_violation_rate = 0)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  expect_true(is.na(result$n_ref_sequences[result$ASV_id == "ASV001"]))
})
