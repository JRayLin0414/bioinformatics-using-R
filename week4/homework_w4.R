load("week4/data/U4_seq1 & seqW.RData")

seqW <- tolower(seqW)
unique(strsplit(seqW, split = "")[[1]])
# [1] "a" "g" "c" "t"

nanog_pat_1 <- "taat[gt][gt]"
nanog_pat_2 <- "[cg][ga][cg]c[gc]atta[acgt][cg]"
nanog_pat_name <- "nanog_pat_1"
nanog_pat <- get(nanog_pat_name)

count_hits <- function(seq, pattern) {
  seq <- tolower(seq)
  pattern <- tolower(pattern)
  hit_pos <- gregexpr(pattern, seq, perl = TRUE)[[1]]
  if (hit_pos[1] == -1) {
    return(0L)
  }

  length(hit_pos)
}

mark_hits <- function(seq, pattern) {
  seq <- tolower(seq)
  pattern <- tolower(pattern)
  hit_pos <- gregexpr(pattern, seq, perl = TRUE)
  regmatches(seq, hit_pos) <- lapply(
    regmatches(seq, hit_pos),
    function(x) paste0("<", x, ">")
  )
  seq
}

get_context <- function(seq, pos, motif, flank = 12) {
  left <- max(1, pos - flank)
  right <- min(nchar(seq), pos + nchar(motif) + flank - 1)
  
  local_seq <- substr(seq, left, right)
  mark_hits(local_seq, motif)
}

# find if there is any sequence matching the pattern
hit_pos <- gregexpr(nanog_pat, seqW, perl = TRUE)[[1]]


# if not -1, then there is at least a match
has_hit <- hit_pos[1] != -1


# how many hits?
hit_n <- if (has_hit) length(hit_pos) else 0


# check if the hit motif is the same
hit_motifs <- if (has_hit) {
  regmatches(seqW, gregexpr(nanog_pat, seqW, perl = TRUE))[[1]]
} else {
  character(0)
}

answer_tbl <- if (has_hit) {
  data.frame(
    position = hit_pos,
    motif = hit_motifs,
    marked_context = mapply(
      get_context,
      pos = hit_pos,
      motif = hit_motifs,
      MoreArgs = list(seq = seqW)
    )
  )
} else {
  data.frame(
    position = integer(0),
    motif = character(0),
    marked_context = character(0)
  )
}

same_motif <- if (hit_n == 0L) NA else length(unique(hit_motifs)) == 1L

cat("Q1. Does seqW have NANOG binding site(s) for", nanog_pat_name, "?", has_hit, "\n")
# Q1. Does seqW have NANOG binding site(s) for nanog_pat_1 ? TRUE
# Q1. Does seqW have NANOG binding site(s) for nanog_pat_2 ? TRUE

cat("Q2. How many for", nanog_pat_name, "?", hit_n, "\n")
# Q2. How many for nanog_pat_1 ? 3
# Q2. How many for nanog_pat_2 ? 3

cat("Q3. Where for", nanog_pat_name, "?", if (has_hit) paste(hit_pos, collapse = ", ") else "none", "\n")
# Q3. Where for nanog_pat_1 ? 8, 328, 1237
# Q3. Where for nanog_pat_2 ? 1346, 2001, 2137

cat("Q4. Do these motifs all the same for", nanog_pat_name, "?", same_motif, "\n")
print(answer_tbl)
# Q4. Do these motifs all the same for nanog_pat_1 ? TRUE
#   position  motif                   marked_context
# 1        8 taattt      agcattt<taattt>ttctctttcaag
# 2      328 taattt tatgtcttctat<taattt>caaactgaaatt
# 3     1237 taattt aaatctggcatc<taattt>tttttggtggac

# Q4. Do these motifs all the same for nanog_pat_2 ? FALSE
#   position       motif                        marked_context
# 1     1346 cagccattacg attcttgataga<cagccattacg>aacctactttat
# 2     2001 cggcgattacg tgaggagtgcgg<cggcgattacg>gagcgcgccaga
# 3     2137 gagcgattacg ggattgtaaatg<gagcgattacg>caccaatcagca


set.seed(42)  # for reproducibility
n_sim <- 5000
sim_counts <- integer(n_sim)

for (i in seq_len(n_sim)) {
  sim_seq <- paste(
    sample(unlist(strsplit(seqW, "")), nchar(seqW), replace = FALSE),
    collapse = ""
  )
  sim_counts[i] <- count_hits(sim_seq, nanog_pat)
}

times_gt_q2 <- sum(sim_counts > hit_n)

cat("Q5. Does seqW have higher frequency of NANOG binding site for", nanog_pat_name, "?", times_gt_q2 < 0.05 * n_sim, "\n")
cat("Randomized seq with hit count > Q2:", times_gt_q2, "out of", n_sim, "\n")
print(summary(sim_counts))
# Q5. Does seqW have higher frequency of NANOG binding site for nanog_pat_1 ? FALSE
# Randomized seq with hit count > Q2: 976 out of 5000
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   1.000   2.000   2.291   3.000   9.000

# Q5. Does seqW have higher frequency of NANOG binding site for nanog_pat_2 ? TRUE
# Randomized seq with hit count > Q2: 0 out of 5000
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.0000  0.0000  0.0596  0.0000  2.0000

save.image(file = "week4/week4.RData")