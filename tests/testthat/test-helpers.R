# Internal helpers are not exported; access via :::

# --- countSpaces -------------------------------------------------------------

test_that("countSpaces returns 0 for strings with no spaces", {
  expect_equal(aasv:::countSpaces("helloworld"), 0L)
  expect_equal(aasv:::countSpaces(""),           0L)
})

test_that("countSpaces counts spaces correctly", {
  expect_equal(aasv:::countSpaces("hello world"),   1L)
  expect_equal(aasv:::countSpaces("a b c"),          2L)
  expect_equal(aasv:::countSpaces("Gorilla gorilla"), 1L)
})

test_that("countSpaces is vectorised", {
  result <- aasv:::countSpaces(c("a b", "abc", "a b c"))
  expect_equal(result, c(1L, 0L, 2L))
})

# --- grantham_score ----------------------------------------------------------

test_that("grantham_score assigns correct classes", {
  expect_equal(aasv:::grantham_score(0),   "conservative")
  expect_equal(aasv:::grantham_score(50),  "conservative")
  expect_equal(aasv:::grantham_score(75),  "moderately conservative")
  expect_equal(aasv:::grantham_score(125), "moderately radical")
  expect_equal(aasv:::grantham_score(200), "radical")
})

test_that("grantham_score is vectorised", {
  result <- aasv:::grantham_score(c(25, 75, 125, 200))
  expect_equal(result, c("conservative", "moderately conservative",
                          "moderately radical", "radical"))
})

# --- triplet_similarity ------------------------------------------------------

test_that("triplet_similarity assigns correct classes", {
  expect_equal(aasv:::triplet_similarity(0.05), "radical")
  expect_equal(aasv:::triplet_similarity(0.2),  "moderately radical")
  expect_equal(aasv:::triplet_similarity(0.5),  "moderately conservative")
  expect_equal(aasv:::triplet_similarity(0.8),  "conservative")
})

# --- phylogeny_score ---------------------------------------------------------

test_that("phylogeny_score returns conservative for 1 and radical otherwise", {
  expect_equal(aasv:::phylogeny_score(1),   "conservative")
  expect_equal(aasv:::phylogeny_score(0.5), "radical")
  expect_equal(aasv:::phylogeny_score(0),   "radical")
})

test_that("phylogeny_score is vectorised", {
  result <- aasv:::phylogeny_score(c(1, 0.8, 1))
  expect_equal(result, c("conservative", "radical", "conservative"))
})

# --- inner_gaps --------------------------------------------------------------

test_that("inner_gaps detects internal gaps", {
  m <- matrix(c("A", "C", "-", "E", "F",
                 "A", "C", "D", "E", "F"),
               nrow = 2L, byrow = TRUE)
  result <- aasv:::inner_gaps(m)
  expect_true(result[1L])
  expect_false(result[2L])
})

test_that("inner_gaps ignores leading terminal gaps", {
  m <- matrix(c("-", "-", "D", "E", "F",
                 "A",  "C", "D", "E", "F"),
               nrow = 2L, byrow = TRUE)
  result <- aasv:::inner_gaps(m)
  expect_false(result[1L])
  expect_false(result[2L])
})

test_that("inner_gaps ignores trailing terminal gaps", {
  m <- matrix(c("A", "C", "D", "-", "-",
                 "A", "C", "D", "E", "F"),
               nrow = 2L, byrow = TRUE)
  result <- aasv:::inner_gaps(m)
  expect_false(result[1L])
  expect_false(result[2L])
})

test_that("inner_gaps handles a single-row matrix", {
  m <- matrix(c("A", "-", "D"), nrow = 1L)
  expect_true(aasv:::inner_gaps(m)[1L])

  m_clean <- matrix(c("A", "C", "D"), nrow = 1L)
  expect_false(aasv:::inner_gaps(m_clean)[1L])
})

test_that("inner_gaps returns a logical vector of length equal to nrow", {
  m <- matrix(c("A", "C", "D",
                 "A", "-", "D",
                 "A", "C", "D"),
               nrow = 3L, byrow = TRUE)
  result <- aasv:::inner_gaps(m)
  expect_type(result, "logical")
  expect_length(result, 3L)
})

# --- causes_internal_gaps -----------------------------------------------------

test_that("causes_internal_gaps flags a sequence that inserts residues absent from all others", {
  # row 1 has a residue at column 3 where every other (spanning) row is gapped:
  # row 1 is forcing an internal gap into rows 2-4, not carrying one itself.
  m <- matrix(c("A", "C", "X", "E", "F",
                "A", "C", "-", "E", "F",
                "A", "C", "-", "E", "F",
                "A", "C", "-", "E", "F"),
              nrow = 4L, byrow = TRUE)
  result <- aasv:::causes_internal_gaps(m)
  expect_true(result[1L])
  expect_false(any(result[2:4]))
})

test_that("causes_internal_gaps does not flag a normal alignment", {
  m <- matrix(c("A", "C", "D", "E", "F",
                "A", "C", "D", "E", "F",
                "A", "C", "D", "E", "F"),
              nrow = 3L, byrow = TRUE)
  result <- aasv:::causes_internal_gaps(m)
  expect_false(any(result))
})

test_that("causes_internal_gaps ignores terminal gaps (shorter sequences, not insertions)", {
  m <- matrix(c("A", "C", "D", "E", "F",
                "-", "-", "D", "E", "F",
                "A", "C", "D", "E", "-"),
              nrow = 3L, byrow = TRUE)
  result <- aasv:::causes_internal_gaps(m)
  expect_false(any(result))
})

# --- best_translation --------------------------------------------------------

test_that("best_translation translates a clean sequence", {
  skip_if_not_installed("Biostrings")
  # ATG=M, CCA=P, GCA=A in any standard code
  result <- aasv:::best_translation("ATGCCAGCA", genetic_code = "5")
  expect_type(result, "character")
  expect_false(is.na(result))
  expect_equal(nchar(result), 3L)
})

test_that("best_translation returns start position when res_type = 'position'", {
  skip_if_not_installed("Biostrings")
  pos <- aasv:::best_translation("ATGCCAGCA", genetic_code = "5",
                                  res_type = "position")
  expect_equal(pos, 1L)
})

test_that("best_translation skips a leading stop codon and finds the next frame", {
  skip_if_not_installed("Biostrings")
  # Frame 1: TAA-CCA-GCA -> TAA is stop (code 5) -> skip
  # Frame 2: AAC-CAG-CA  -> N, Q + partial -> no stop -> position 2
  result_pos <- aasv:::best_translation("TAACCAGCA", genetic_code = "5",
                                         res_type = "position")
  expect_equal(result_pos, 2L)
})

test_that("best_translation returns NA for an empty input sequence", {
  skip_if_not_installed("Biostrings")
  # nchar("") == 0 so max_iter == 0; the while-loop never runs; res stays NA
  result <- aasv:::best_translation("", genetic_code = "5")
  expect_true(is.na(result))
})

# --- aa_table ----------------------------------------------------------------

test_that("aa_table returns a data.frame with required columns", {
  skip_if_not_installed("Biostrings")
  x      <- Biostrings::DNAStringSet(c(seq1 = "ATGCCAGCA"))
  result <- aasv:::aa_table(x, genetic_code = "5")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("seq_id", "position", "triplet", "aa") %in% names(result)))
})

test_that("aa_table has one row per codon", {
  skip_if_not_installed("Biostrings")
  x      <- Biostrings::DNAStringSet(c(seq1 = "ATGCCAGCA"))  # 3 codons
  result <- aasv:::aa_table(x, genetic_code = "5")
  expect_equal(nrow(result), 3L)
})

test_that("aa_table preserves sequence names and triplet content", {
  skip_if_not_installed("Biostrings")
  x      <- Biostrings::DNAStringSet(c(myseq = "ATGCCAGCA"))
  result <- aasv:::aa_table(x, genetic_code = "5")
  expect_true(all(result$seq_id == "myseq"))
  expect_equal(result$triplet[1L], "ATG")
  expect_equal(result$position,    c(1L, 2L, 3L))
})

test_that("aa_table handles multiple sequences", {
  skip_if_not_installed("Biostrings")
  x <- Biostrings::DNAStringSet(c(s1 = "ATGCCAGCA", s2 = "ATGAAATTT"))
  result <- aasv:::aa_table(x, genetic_code = "5")
  expect_equal(nrow(result), 6L)
  expect_equal(sort(unique(result$seq_id)), c("s1", "s2"))
})
