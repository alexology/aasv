skip_if_not_installed("writexl")
skip_if_not_installed("readxl")

# Helper: create a minimal aa_structure_results directory
make_classify_dir <- function() {
  tmpdir <- tempfile(pattern = "aasv_classify_")
  dir.create(file.path(tmpdir, "aa_structure_results"), recursive = TRUE)
  tmpdir
}

# Helper: write a minimal hydrophobicity Excel file
write_hydro <- function(dir, asv_id = "ASV001", tax_lev = "species",
                        global_pass = TRUE) {
  df <- data.frame(
    ASV_id               = asv_id,
    tax_lev              = tax_lev,
    global_min           = -1,
    global_max           = 1,
    global_score         = 0.5,
    global_pass          = global_pass,
    local_violation_rate = 0.0,
    structural_verdict   = "PASS"
  )
  writexl::write_xlsx(df, file.path(dir, "aa_structure_results", "hydrophobicity.xlsx"))
}

# Helper: write an ASV structure Excel file with a given aa_pos
write_struct <- function(dir, asv_id = "ASV001", tax_lev = "species",
                         aa_pos = 5L,
                         grantham_dist = 20, median_quality_prob = 0.9999,
                         triplet_mut_pos = "3") {
  df <- data.frame(
    ASV_id               = asv_id,
    tax_lev              = tax_lev,
    aa_pos               = aa_pos,
    illegal_aa           = if (aa_pos == -1L) NA_character_ else "X",
    dna_triplet          = NA_character_,
    nt_coords            = NA_character_,
    n_sequences_ali      = NA_integer_,
    n_unique_seq_ali     = NA_integer_,
    median_quality_prob  = median_quality_prob,
    n_sequences_prob     = NA_integer_,
    ref_aa               = NA_character_,
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

# --- output structure --------------------------------------------------------

test_that("classify_asv returns a data.frame with required columns", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir)

  result <- classify_asv(tmpdir)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("ASV_id", "tax_lev", "IS_REAL", "REASON") %in% names(result)))
})

test_that("classify_asv errors when hydrophobicity file is missing and include_hydro = TRUE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_error(classify_asv(tmpdir, include_hydro = TRUE), "Hydrophobicity file not found")
})

# --- IS_REAL classification --------------------------------------------------

test_that("classify_asv marks ASV as FALSE when alignment gap present (aa_pos = -1)", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir)
  write_struct(tmpdir, aa_pos = -1L)

  result <- classify_asv(tmpdir)
  expect_false(result$IS_REAL[result$ASV_id == "ASV001"])
})

test_that("classify_asv marks ASV as TRUE when mutation passes all filters", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, global_pass = TRUE)
  # Wobble position (pos 3), conservative grantham, high quality
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20,
               median_quality_prob = 0.9999, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  expect_true(result$IS_REAL[result$ASV_id == "ASV001"])
})

test_that("classify_asv marks ASV as FALSE when hydrophobicity fails", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, global_pass = FALSE)
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20,
               median_quality_prob = 0.9999, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir)
  expect_false(result$IS_REAL[result$ASV_id == "ASV001"])
})

test_that("classify_asv ignores hydrophobicity when include_hydro = FALSE", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, global_pass = FALSE)  # would normally fail
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 20,
               median_quality_prob = 0.9999, triplet_mut_pos = "3")

  result <- classify_asv(tmpdir, include_hydro = FALSE)
  expect_true(result$IS_REAL[result$ASV_id == "ASV001"])
})

test_that("classify_asv returns empty data.frame when no files are present", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  result <- classify_asv(tmpdir, include_hydro = FALSE)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("classify_asv marks ASV as FALSE when grantham distance is too high", {
  tmpdir <- make_classify_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_hydro(tmpdir, global_pass = TRUE)
  # Position 2 (Forbidden Zone) AND high grantham distance
  write_struct(tmpdir, aa_pos = 5L, grantham_dist = 200,
               median_quality_prob = 0.5, triplet_mut_pos = "2")

  result <- classify_asv(tmpdir)
  expect_false(result$IS_REAL[result$ASV_id == "ASV001"])
})
