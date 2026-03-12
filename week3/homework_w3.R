library(data.table)

file_path <- "week3/data/ENST00000589042_TTN.csv"

# --------------------------------- read.table() version -----------------------------------------

# load the data to the table
TTN.table <- read.table(
  file_path,
  sep = ",",
  header = TRUE,
)

# Create $cls to separate intron, exon and flanking(up & downstream)
TTN.table$cls <- "flanking" # set all records to "flanking" as default for $cls
TTN.table$cls[grepl("^ENSE", TTN.table$`Exon / Intron`)] <- "exon"  # if contain "ENSE" in the name, set $cls to "exon"
TTN.table$cls[grepl("^Intron", TTN.table$`Exon / Intron`)] <- "intron"  # if contain "Intron" in the name, set $cls to "intron"

# Create $len to store processed length
TTN.table$len <- gsub(",", "", TTN.table$Length)  # remvoe the ',' from the $Length
TTN.table$len <- as.integer(TTN.table$len)  # make sure all the value is integer

# Sort the table according to cls and len
TTN.table <- TTN.table[TTN.table$cls %in% c("exon", "intron"), ]  # remove flankings
TTN.table <- TTN.table[order(TTN.table$cls, TTN.table$len), ] # first order by cls and then len

# Get all the introns from the table
TTN.table.intron <- TTN.table[TTN.table$cls == "intron", ]
TTN.table.intron$Sequence <-gsub(" ", "", TTN.table.intron$Sequence)  # Renmove the white space in the $Sequence

TTN.table.intron$seq <- substr(
  TTN.table.intron$Sequence,
  pmax(1, nchar(TTN.table.intron$Sequence) - 60 + 1),
  nchar(TTN.table.intron$Sequence)
) # Get the seq from the tail(60 chars)

# Get the first and last two nucleotides from the $seq
TTN.table.intron$first_two_nt <- substr(TTN.table.intron$seq, 1, 2)
TTN.table.intron$last_two_nt  <- substr(
  TTN.table.intron$seq,
  nchar(TTN.table.intron$seq) - 1,
  nchar(TTN.table.intron$seq)
)

# Just for curiosity
unique(TTN.table.intron$first_two_nt) # [1] "gt" "gc"
unique(TTN.table.intron$last_two_nt) # [1] "ag"


# --------------------------------- fread() / data.table version -----------------------------------------

# Read the data to data.table
TTN.dt <- fread(file_path, check.names = FALSE)

# Create $cls to separate intron, exon and flanking(up & downstream)
TTN.dt[, cls := "flanking"]
TTN.dt[grepl("^ENSE", `Exon / Intron`), cls := "exon"]
TTN.dt[grepl("^Intron", `Exon / Intron`), cls := "intron"]

# Create $len to store processed length
TTN.dt[, len := gsub(",", "", Length)]
TTN.dt[len == "", len := NA_character_]
TTN.dt[, len := as.integer(len)]


TTN.dt <- TTN.dt[cls %in% c("exon", "intron")]
setorder(TTN.dt, cls, len)

# Get all the introns from the table
TTN.dt.intron <- TTN.dt[cls == "intron"]

# Remove whitespace in Sequence
TTN.dt.intron[, Sequence := gsub(" ", "", Sequence)]

# Get the seq from the tail (60 chars)
TTN.dt.intron[, seq := substr(
  Sequence,
  pmax(1, nchar(Sequence) - 60 + 1),
  nchar(Sequence)
)]

# Get the first and last two nucleotides from seq
TTN.dt.intron[, first_two_nt := substr(seq, 1, 2)]
TTN.dt.intron[, last_two_nt := substr(
  seq,
  nchar(seq) - 1,
  nchar(seq)
)]


save.image(file = "week3/week3.RData")
