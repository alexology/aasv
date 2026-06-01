# Shared helpers loaded automatically by testthat before every test file.

# Creates a temp alignment dir with the standard folder structure.
# Returns the root path; caller owns cleanup via on.exit(unlink(..., recursive=TRUE)).
make_alignment_dir <- function(family = "Testidae") {
  tmpdir <- tempfile(pattern = "aasv_test_")
  dir.create(file.path(tmpdir, "alignments", family, "species"), recursive = TRUE)
  tmpdir
}

# Writes a named character vector of AA sequences (gaps allowed) to a FASTA file.
write_aa_fasta <- function(seqs, path) {
  x <- Biostrings::AAStringSet(seqs)
  Biostrings::writeXStringSet(x, filepath = path, format = "fasta")
  invisible(path)
}
