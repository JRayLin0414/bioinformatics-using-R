library(data.table)

# Use fread() to load the data
DNA_seq <- fread("week4/data/DNA_sequences.txt")
#        ID                                         SEQ
#       <int>                                      <char>
#    1:     1 gacaccctcatggtgctgtgACAaatcattgccttcaaatgaa
#    2:     2 tgtgactgtatacctgtcatGATtgtgttatgaaaataataat
#    3:     3 TCATTAGAGATTGGCTGGTGGAAtcatctgataatctatttat
#    4:     4 tactgactcttaaaaaaaaaATAGGTGAAACCTATACGAGACC
#    5:     5 agtagaatgagatttttatgTATaaaaaaacttttaataatta
#   ---
# 1444:  1444 GAGACCCATATTTTACAGGATATCACACAGGAATTACAAGATG
# 1445:  1445 tgtacctactggtttgattaTGCAGGAAAACTTCAAGATTATA
# 1446:  1446 CACCAAGATTTGTGAAGAAGCCTattatggtatgtatccattt
# 1447:  1447 ttttgatggctgaatatagaATTGAAGCTCAGTGACATATCTA
# 1448:  1448 cttttaaaatagcttcttaaCTTtagagtttaaattttgtaca

# convert all char to lowercase
DNA_seq$SEQ <- tolower(DNA_seq$SEQ)

# Check if the $ID is valid
id_ok <- is.integer(DNA_seq$ID) &&
  !anyNA(DNA_seq$ID) &&
  !anyDuplicated(DNA_seq$ID) &&
  identical(DNA_seq$ID, seq_len(nrow(DNA_seq)))

cat("Is $ID valid?", id_ok, "\n")  # TRUE

bad_id <- DNA_seq[
  is.na(ID) | duplicated(ID) | ID != seq_len(.N)
]
bad_id


# Check if the $SEQ is valid
seq_ok <- grepl("^[ACGTacgt]+$", DNA_seq$SEQ)

cat("How many invalid $SEQ?", sum(!seq_ok), "\n") # 7

bad_seq <- DNA_seq[!seq_ok]
bad_seq
#       ID                                         SEQ
#    <int>                                      <char>
# 1:    51 tccggctcaaattcctcaccacagcatgaqcggatgcggaaac
# 2:   118 caccatccttcgtgaagaaggtvagttgtccttgaaaagaaag
# 3:   225 atttttggactagcagtattxattagcacaatatatttgtaca
# 4:   281 ttccgcctaaaattgaagctgtjgaaagttcctgaagaaagga
# 5:   553 aaatataaacatttatcatatktacctctggaacatgtggaag
# 6:  1009                                       #ref!
# 7:  1147 tgtggattcttacetaactttcaagttgtacttgaaaagaaag