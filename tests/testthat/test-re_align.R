test_that("re_align handles empty paths without error", {
  expect_no_error(re_align(character(0), verbose = FALSE))
})

test_that("re_align returns input paths invisibly", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  result <- suppressMessages(re_align(ali_path, verbose = FALSE))
  expect_equal(result, ali_path)
})

test_that("re_align overwrites the original file", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  mtime_before <- file.info(ali_path)$mtime
  Sys.sleep(0.05)
  suppressMessages(re_align(ali_path, verbose = FALSE))

  expect_gt(as.numeric(file.info(ali_path)$mtime), as.numeric(mtime_before))
})

test_that("re_align with verbose = FALSE produces no messages", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  expect_no_message(re_align(ali_path, verbose = FALSE))
})

test_that("re_align strips gaps before re-aligning (raw_sequences = FALSE)", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(s1 = "ACDEF--GHIKL", s2 = "ACDEFGHIKLMN"), ali_path)

  suppressMessages(re_align(ali_path, verbose = FALSE))

  result <- Biostrings::readAAStringSet(ali_path)
  expect_gte(length(result), 1L)
  # output sequences must all have the same width (valid alignment)
  if (length(result) > 1L) {
    expect_equal(length(unique(Biostrings::width(result))), 1L)
  }
})

test_that("re_align preserves sequence names", {
  skip_if_not_installed("DECIPHER")

  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  ali_path <- file.path(tmpdir, "Testus_alpha.fasta")
  write_aa_fasta(c(SeqA = "ACDEF--GHIKL", SeqB = "ACDEFGHIKLMN"), ali_path)

  suppressMessages(re_align(ali_path, verbose = FALSE))

  result <- Biostrings::readAAStringSet(ali_path)
  expect_setequal(names(result), c("SeqA", "SeqB"))
})
