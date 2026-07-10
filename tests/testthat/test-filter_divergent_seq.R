run_filter <- function(path, ...) {
  suppressMessages(
    filter_divergent_seq(path, print_messages = FALSE, verbose = FALSE, ...)
  )
}

test_that("filter_divergent_seq handles empty paths without error", {
  expect_no_error(filter_divergent_seq(character(0), print_messages = FALSE))
})

test_that("filter_divergent_seq returns input paths invisibly", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "ACDEFGHIKLMNPQRSTVWQ",
      s4 = "CDEAEACDEFGHIKLMNPQR"),
    ali_path
  )

  result <- run_filter(ali_path)
  expect_equal(result, ali_path)
})

test_that("filter_divergent_seq removes a sequence far from all others", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # s1/s2/s3 are within 5% of each other (1 differing site out of 20);
  # s4 differs from every other sequence at every site (100% divergent).
  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "ACDEFGHIKLMNPQRSTVWQ",
      s4 = "CDEACDEAEACDEFGHIKLM"),
    ali_path
  )

  run_filter(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 3L)
  expect_false("s4" %in% names(result))
})

test_that("filter_divergent_seq keeps all sequences within the threshold", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "ACDEFGHIKLMNPQRSTVWQ",
      s4 = "ACDEFGHIKLMNPQRSTVWQ"),
    ali_path
  )

  run_filter(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 4L)
})

test_that("filter_divergent_seq skips files below min_seqs", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Only 3 sequences (below the default min_seqs = 4), even though s3 is
  # wildly divergent from s1/s2 - the file must be left untouched.
  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "CDEACDEAEACDEFGHIKLM"),
    ali_path
  )

  run_filter(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 3L)
})

test_that("filter_divergent_seq keeps everything and warns when all sequences are divergent", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # A cyclic Latin square of order 4: any two rows differ at every column,
  # so every pairwise distance is 100% and every sequence is flagged.
  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDE", s2 = "CDEA", s3 = "DEAC", s4 = "EACD"),
    ali_path
  )

  expect_warning(
    filter_divergent_seq(ali_path, print_messages = FALSE, verbose = FALSE),
    "every sequence exceeds max_distance"
  )

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 4L)
})

test_that("filter_divergent_seq respects a custom max_distance", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # s3 is the lone sequence differing from s1/s2/s4 at 1 of 20 sites (5%) -
  # kept at the default 10% threshold but removed once tightened to 1%.
  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "ACDEFGHIKLMNPQRSTVWQ",
      s4 = "ACDEFGHIKLMNPQRSTVWY"),
    ali_path
  )

  run_filter(ali_path, max_distance = 0.01)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 3L)
  expect_false("s3" %in% names(result))
})

test_that("filter_divergent_seq writes a removal log", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))
  log_dir <- file.path(tmpdir, "logs")

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKLMNPQRSTVWY",
      s2 = "ACDEFGHIKLMNPQRSTVWY",
      s3 = "ACDEFGHIKLMNPQRSTVWQ",
      s4 = "CDEACDEAEACDEFGHIKLM"),
    ali_path
  )

  suppressMessages(
    filter_divergent_seq(ali_path, print_messages = FALSE, verbose = FALSE,
                         log_dir = log_dir)
  )

  log_file <- file.path(log_dir, "removed_divergent.tsv")
  expect_true(file.exists(log_file))

  log_df <- read.table(log_file, header = TRUE, sep = "\t")
  expect_equal(nrow(log_df), 1L)
  expect_equal(log_df$sequence, "s4")
  expect_equal(log_df$n_start, 4L)
  expect_equal(log_df$n_retained, 3L)
})

test_that("filter_divergent_seq keeps single-sequence files unchanged", {
  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEFGHIKL"), ali_path)

  run_filter(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 1L)
  expect_equal(as.character(result), c(s1 = "ACDEFGHIKL"))
})
