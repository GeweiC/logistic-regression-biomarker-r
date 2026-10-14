# Polynomial logistic models on the original biomarker scale (degree 1 to 5).
source(file.path("R", "helpers.R"), local = TRUE)
d <- read_immune_data()
prepare_output_dirs()
models <- fit_polynomials(d, "immune")
lrt <- nested_lrt_table(models)
write.csv(lrt, file.path("results", "raw_polynomial_lrt.csv"), row.names=FALSE)
# At alpha=0.05, the second degree improves on the first, but the third does not.
chosen <- models[[2L]]
threshold <- threshold_from_quadratic(chosen, d$immune)
accuracy <- grouped_classification(d, threshold, "Quadratic original scale")
write.csv(accuracy, file.path("results", "raw_quadratic_classification.csv"), row.names=FALSE)

x_grid <- seq(min(d$immune), max(d$immune), length.out=400)
pred <- predict(chosen, newdata=data.frame(immune=x_grid), type="response")
png(file.path("figures", "raw_quadratic_fit.png"), width=1200, height=750, res=140)
plot(d$immune, 100 * d$proportion, pch=19, col="#2563A6", cex=1.2,
     xlab="Biomarker level", ylab="Disease probability (%)", ylim=c(0,100),
     main="Quadratic logistic model: biomarker on the original scale")
lines(x_grid, 100 * pred, col="#E07A38", lwd=3)
abline(v=threshold, h=50, col="grey45", lty=2)
legend("bottomright", legend=c("Observed group proportion", "Fitted quadratic model", "50% probability / threshold"),
       col=c("#2563A6", "#E07A38", "grey45"), pch=c(19,NA,NA),
       lty=c(NA,1,2), lwd=c(NA,3,1), bty="n")
dev.off()

cat("Original-scale polynomial LRT summary:\n")
print(lrt, row.names=FALSE)
cat(sprintf("Quadratic threshold: %.4f; in-sample accuracy: %.2f%% (%d/%d).\n",
            threshold, 100 * accuracy$accuracy, accuracy$correct, accuracy$total))
