# getwd()

# Read the csv to R
data <- read.table("./week2/ENST00000589042_TTN.csv", sep = ",", header = TRUE)

# What is the dim of the table?
table_dim <- dim(data)
print(paste("The dimension of the table is:", table_dim[1], "rows and", table_dim[2], "columns"))

# Get the index of the exon and intron
e.idx <- which(substr(data$Exon...Intron, 1, 4) == "ENSE")
i.idx <- which(substr(data$Exon...Intron, 1, 6) == "Intron")

# Get the length of each exon and intron
e.length <- as.numeric(gsub(",", "", data$Length[e.idx]))
i.length <- as.numeric(gsub(",", "", data$Length[i.idx]))

# Check type
# is.numeric(e.length)
# is.numeric(i.length)

# Make sure the amount of the exon is equals to 363 (accroding to the no. in the csv)
# print(length(e.length) == 363)

# Calculate the median and average length of the exons and introns
e.avg <- mean(e.length)
e.med <- median(e.length)

i.avg <- mean(i.length)
i.med <- median(i.length)

print(paste("The median length of the exons is:", e.med, "and the average length of the exons is:", e.avg))
print(paste("The median length of the introns is:", i.med, "and the average length of the introns is:", i.avg))
