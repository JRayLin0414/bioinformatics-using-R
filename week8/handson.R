load('week8/data/U6_ensg2GS2uniprot.Rdata')

setT <- ensg2GS2uniprot
setZ <- ensg2GS2uniprot[
  !is.na(ensg2GS2uniprot$featureNum) & ensg2GS2uniprot$featureNum == 0,
]
setN <- ensg2GS2uniprot[is.nan(ensg2GS2uniprot$featureNum), ]

genesT <- sort(unique(setT$ENSG))
genesZ <- sort(unique(setZ$ENSG))
genesN <- sort(unique(setN$ENSG))

countZN <- length(intersect(genesZ, genesN))
countZOnly <- length(setdiff(genesZ, genesN))
countNOnly <- length(setdiff(genesN, genesZ))
countTOnly <- length(setdiff(genesT, union(genesZ, genesN)))

plotSets <- function() {
  op <- par(mar = c(1, 1, 4, 1))
  on.exit(par(op))

  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 10), asp = 1)

  rect(0.8, 0.8, 9.2, 9.2, border = "grey35", lwd = 2, col = "grey95")
  text(1.2, 9.45, sprintf("setT (total) = %d", length(genesT)), adj = c(0, 0), cex = 1)

  symbols(4.2, 5.0, circles = 2.2, inches = FALSE, add = TRUE,
          bg = adjustcolor("steelblue", alpha.f = 0.35), fg = "steelblue4")
  symbols(5.8, 5.0, circles = 1.9, inches = FALSE, add = TRUE,
          bg = adjustcolor("tomato", alpha.f = 0.35), fg = "tomato4")

  text(3.0, 7.45, sprintf("setZ = %d", length(genesZ)), col = "steelblue4", cex = 1)
  text(7.0, 7.05, sprintf("setN = %d", length(genesN)), col = "tomato4", cex = 1)

  text(3.3, 5.0, countZOnly, cex = 1.2)
  text(6.7, 5.0, countNOnly, cex = 1.2)
  text(5.0, 5.0, countZN, cex = 1.2, font = 2)
  text(5.0, 1.5, sprintf("setT only = %d", countTOnly), cex = 1.1)

  title("Nested set diagram of ENSG sets")
}

if (interactive()) {
  plotSets()
} else {
  png("week8/handson_venn.png", width = 1200, height = 1200, res = 180)
  plotSets()
  dev.off()
}
