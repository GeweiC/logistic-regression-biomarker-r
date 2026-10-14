# Grouped biomarker disease proportions: exploratory summary.
source(file.path("R", "helpers.R"), local = TRUE)
d <- read_immune_data()
prepare_output_dirs()
summary_table <- d[, c("immune", "disease", "nodisease", "sample_size", "proportion")]
write.csv(summary_table, file.path("results", "grouped_proportions.csv"), row.names = FALSE)

png(file.path("figures", "observed_biomarker_proportions.png"), width=1200, height=750, res=140)
plot(d$immune, 100 * d$proportion, pch=19, col="#2563A6", cex=1.25,
     xlab="Biomarker level", ylab="Observed disease proportion (%)",
     ylim=c(0, 100), main="Observed disease proportion by biomarker level")
grid(col="grey90", lty=3)
points(d$immune, 100 * d$proportion, pch=19, col="#2563A6", cex=1.25)
dev.off()

cat("Exploration complete: ", nrow(d), " groups and ", sum(d$sample_size),
    " people (grouped counts).\n", sep="")
