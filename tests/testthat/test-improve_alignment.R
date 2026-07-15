run_improve <- function(path) {
  suppressMessages(
    improve_alignment(path, print_messages = FALSE, verbose = FALSE)
  )
}

test_that("improve_alignment handles empty paths without error", {
  expect_no_error(improve_alignment(character(0), print_messages = FALSE))
})

test_that("improve_alignment returns input paths invisibly", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN", s3 = "ACDEFGHIKLMN"),
    ali_path
  )

  result <- run_improve(ali_path)
  expect_equal(result, ali_path)
})

test_that("improve_alignment overwrites the original file", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN", s3 = "ACDEFGHIKLMN"),
    ali_path
  )

  mtime_before <- file.info(ali_path)$mtime
  Sys.sleep(0.05)
  run_improve(ali_path)

  expect_gt(
    as.numeric(file.info(ali_path)$mtime),
    as.numeric(mtime_before)
  )
})

test_that("improve_alignment produces no messages when silenced", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN", s3 = "ACDEFGHIKLMN"),
    ali_path
  )

  expect_no_message(
    improve_alignment(ali_path, print_messages = FALSE, verbose = FALSE)
  )
})

test_that("improve_alignment retains sequences that re-align without gaps", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # s1 has an apparent internal gap in the input alignment, but after
  # strip_gaps + re-align it becomes a shorter sequence with only trailing
  # gaps, it must NOT be removed.
  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN", s3 = "ACDEFGHIKLMN"),
    ali_path
  )

  run_improve(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 3L)
  expect_true("s1" %in% names(result))
})

test_that("improve_alignment output sequences have equal width", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN", s3 = "ACDEFGHIKLMN"),
    ali_path
  )

  run_improve(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  if (length(result) > 1L) {
    expect_equal(length(unique(Biostrings::width(result))), 1L)
  }
})

test_that("improve_alignment keeps single-sequence files unchanged", {
  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEFGHIKL"), ali_path)

  run_improve(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 1L)
  expect_equal(as.character(result), c(s1 = "ACDEFGHIKL"))
})

test_that("improve_alignment keeps all sequences when none have internal gaps", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(
    c(s1 = "ACDEFGHIKL--", s2 = "ACDEFGHIKLMN", s3 = "--DEFGHIKLMN"),
    ali_path
  )

  run_improve(ali_path)

  result <- Biostrings::readAAStringSet(ali_path)
  expect_equal(length(result), 3L)
})
