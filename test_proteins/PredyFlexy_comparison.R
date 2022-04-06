#!/usr/bin/Rscript
data_flexamino <- read.table("optionals/Q8IU85_headerless.txt", header = F)
data_predyflexy <- read.table("optionals/Q8IU85_headerless.predyflexy", header = F)

data_predyflexy[data_predyflexy[10] == -9.999, 10] = NA

cor.test(data_flexamino[, 3], data_predyflexy[, 10])

pdf("Q8IU85_scores_comparison.png")
plot(data_flexamino[, 3]~data_predyflexy[, 10], ylab = "FlexAmino score",
     xlab = "PredyFlexy score", main = "Comparison of program scores")
dev.off()

