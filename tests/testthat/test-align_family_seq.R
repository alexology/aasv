make_genus_dir <- function(family = "Testidae") {
  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(file.path(tmpdir, "alignments", family, "genus"), recursive = TRUE)
  tmpdir
}

test_that("align_family_seq skips and messages when no genus files are present", {
  tmpdir <- make_genus_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_message(
    align_family_seq("Testidae", alignment_dir = tmpdir),
    "no genus files found"
  )
  expect_false(dir.exists(file.path(tmpdir, "alignments", "Testidae", "family")))
})

test_that("align_family_seq creates family directory and output file", {
  tmpdir <- make_genus_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKL", s2 = "MGNPQRSTVWY"),
    file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")
  )

  suppressMessages(align_family_seq("Testidae", alignment_dir = tmpdir))

  expect_true(dir.exists(file.path(tmpdir, "alignments", "Testidae", "family")))
  expect_true(file.exists(file.path(tmpdir, "alignments", "Testidae", "family", "Testidae.fasta")))
})

test_that("align_family_seq skips already-processed families when overwrite = FALSE", {
  tmpdir <- make_genus_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKL"),
    file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")
  )

  suppressMessages(align_family_seq("Testidae", alignment_dir = tmpdir))

  family_file <- file.path(tmpdir, "alignments", "Testidae", "family", "Testidae.fasta")
  mtime_first <- file.info(family_file)$mtime
  Sys.sleep(0.05)

  expect_message(
    align_family_seq("Testidae", alignment_dir = tmpdir, overwrite = FALSE),
    "already present"
  )
  expect_equal(file.info(family_file)$mtime, mtime_first)
})

test_that("align_family_seq output is a readable AAStringSet", {
  tmpdir <- make_genus_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKL", s2 = "MGNPQRSTVWY"),
    file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")
  )

  suppressMessages(align_family_seq("Testidae", alignment_dir = tmpdir))

  out <- Biostrings::readAAStringSet(
    file.path(tmpdir, "alignments", "Testidae", "family", "Testidae.fasta")
  )
  expect_s4_class(out, "AAStringSet")
  expect_gte(length(out), 1L)
})
