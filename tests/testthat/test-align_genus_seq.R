test_that("align_genus_seq skips and messages when no species files are present", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_message(
    align_genus_seq("Testidae", alignment_dir = tmpdir),
    "no species files found"
  )
  expect_false(dir.exists(file.path(tmpdir, "alignments", "Testidae", "genus")))
})

test_that("align_genus_seq creates genus directory and output file", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(seq1 = "ACDEFGHIKL", seq2 = "MGNPQRSTVWY"),
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  suppressMessages(align_genus_seq("Testidae", alignment_dir = tmpdir))

  expect_true(dir.exists(file.path(tmpdir, "alignments", "Testidae", "genus")))
  expect_true(file.exists(file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")))
})

test_that("align_genus_seq extracts genus correctly from filename", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  sp <- file.path(tmpdir, "alignments", "Testidae", "species")
  write_aa_fasta(c(s1 = "ACDEFGHIKL"),  file.path(sp, "GenusA_speciesX.fasta"))
  write_aa_fasta(c(s1 = "MGNPQRSTVWY"), file.path(sp, "GenusB_speciesY.fasta"))

  suppressMessages(align_genus_seq("Testidae", alignment_dir = tmpdir))

  genus_dir <- file.path(tmpdir, "alignments", "Testidae", "genus")
  expect_true(file.exists(file.path(genus_dir, "GenusA.fasta")))
  expect_true(file.exists(file.path(genus_dir, "GenusB.fasta")))
})

test_that("align_genus_seq skips already-processed families when overwrite = FALSE", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKL"),
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  suppressMessages(align_genus_seq("Testidae", alignment_dir = tmpdir))

  genus_file  <- file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")
  mtime_first <- file.info(genus_file)$mtime
  Sys.sleep(0.05)

  expect_message(
    align_genus_seq("Testidae", alignment_dir = tmpdir, overwrite = FALSE),
    "already present"
  )
  expect_equal(file.info(genus_file)$mtime, mtime_first)
})

test_that("align_genus_seq output is a readable AAStringSet", {
  tmpdir <- make_alignment_dir()
  on.exit(unlink(tmpdir, recursive = TRUE))

  write_aa_fasta(
    c(s1 = "ACDEFGHIKL", s2 = "MGNPQRSTVWY"),
    file.path(tmpdir, "alignments", "Testidae", "species", "Testus_alpha.fasta")
  )

  suppressMessages(align_genus_seq("Testidae", alignment_dir = tmpdir))

  out <- Biostrings::readAAStringSet(
    file.path(tmpdir, "alignments", "Testidae", "genus", "Testus.fasta")
  )
  expect_s4_class(out, "AAStringSet")
  expect_gte(length(out), 1L)
})
