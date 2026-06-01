test_that("check_alignment errors on empty directory", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_error(
    check_alignment("Testidae", alignment_dir = tmpdir),
    "No files in the specified directory"
  )
})

test_that("check_alignment returns empty vector when no gaps present", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKLM", s2 = "MGNPQRSTVWY"),   # equal width (11)
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  result <- suppressMessages(check_alignment("Testidae", alignment_dir = tmpdir))
  expect_equal(length(result), 0L)
})

test_that("check_alignment returns empty vector for terminal gaps only", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  # leading and trailing gaps are acceptable — not internal
  write_aa_fasta(
    c(s1 = "ACDEFGHIKL--", s2 = "ACDEFGHIKLMN", s3 = "--DEFGHIKLMN"),
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  result <- suppressMessages(check_alignment("Testidae", alignment_dir = tmpdir))
  expect_equal(length(result), 0L)
})

test_that("check_alignment detects internal gaps and returns the file path", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(
    tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta"
  )
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  result <- suppressMessages(
    check_alignment("Testidae", alignment_dir = tmpdir)
  )
  expect_equal(result, ali_path)
})

test_that("check_alignment skips single-sequence files", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  # single sequence with an internal gap character — must be skipped
  write_aa_fasta(
    c(s1 = "ACDEF--GHIKL"),
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  result <- suppressMessages(check_alignment("Testidae", alignment_dir = tmpdir))
  expect_equal(length(result), 0L)
})

test_that("check_alignment returns only paths with internal gaps (mixed files)", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  sp       <- file.path(tmpdir, "alignments", "Testidae", "species")
  clean    <- file.path(sp, "Testus_alpha.fasta")
  bad      <- file.path(sp, "Testus_beta.fasta")

  write_aa_fasta(c(s1 = "ACDEFGHIKL--", s2 = "ACDEFGHIKLMN"), clean)
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), bad)

  result <- suppressMessages(check_alignment("Testidae", alignment_dir = tmpdir))
  expect_equal(result, bad)
})

test_that("check_alignment respects tax_lev argument", {
  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(
    file.path(tmpdir, "alignments", "Testidae", "genus"), recursive = TRUE
  )
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(
    tmpdir, "alignments", "Testidae", "genus", "Testus.fasta"
  )
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  result <- suppressMessages(
    check_alignment("Testidae", alignment_dir = tmpdir, tax_lev = "genus")
  )
  expect_equal(result, ali_path)
})
