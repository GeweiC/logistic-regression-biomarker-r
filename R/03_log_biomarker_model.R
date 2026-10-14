# Log-biomarker logistic regressions and in-sample classification comparison.
source(file.path("R", "helpers.R"), local=TRUE)
d <- read_immune_data()
prepare_output_dirs()
models_log <- fit_polynomials(d, "log_immune")
lrt_log <- nested_lrt_table(models_log)
write.csv(lrt_log, file.path("results", "log_polynomial_lrt.csv"), row.names=FALSE)

# Sequential LRTs with alpha=0.05 favour the parsimonious linear log model.
# Note: AIC of the degree-2 log model is marginally lower; see methodology.
chosen <- models_log[[1L]]
threshold <- threshold_from_log_model(chosen)
acc_log <- grouped_classification(d, threshold, "Linear log-biomarker")
raw_model <- fit_polynomials(d, "immune", degrees=2L)[[1]]
raw_threshold <- threshold_from_quadratic(raw_model, d$immune)
acc_raw <- grouped_classification(d, raw_threshold, "Quadratic original scale")
total <- sum(d$sample_size)
acc_majority <- data.frame(method="Predict disease for all", threshold=NA_real_,
                            true_negative=0, true_positive=sum(d$disease),
                            correct=sum(d$disease), total=total,
                            accuracy=sum(d$disease)/total)
comparison <- rbind(acc_majority, acc_raw, acc_log)
write.csv(comparison, file.path("results", "in_sample_accuracy.csv"), row.names=FALSE)
write.csv(data.frame(term=names(coef(chosen)), coefficient=unname(coef(chosen))),
          file.path("results", "log_model_coefficients.csv"), row.names=FALSE)

x_grid <- seq(min(d$immune), max(d$immune), length.out=400)
pred_log <- predict(chosen, newdata=data.frame(log_immune=log(x_grid)), type="response")
pred_quad <- predict(raw_model, newdata=data.frame(immune=x_grid), type="response")
png(file.path("figures", "log_biomarker_fit.png"), width=1200, height=750, res=140)
plot(d$immune, 100 * d$proportion, pch=19, cex=1.2, col="#2563A6",
     xlab="Biomarker level", ylab="Disease probability (%)", ylim=c(0,100),
     main="Log-biomarker logistic model and 50% classification threshold")
lines(x_grid, 100 * pred_log, col="#E07A38", lwd=3)
abline(v=threshold, h=50, col="grey45", lty=2)
legend("bottomright", legend=c("Observed group proportion", "Fitted log-biomarker model", "50% probability / threshold"),
       col=c("#2563A6", "#E07A38", "grey45"), pch=c(19,NA,NA),
       lty=c(NA,1,2), lwd=c(NA,3,1), bty="n")
dev.off()

png(file.path("figures", "model_curves_comparison.png"), width=1200, height=750, res=140)
plot(d$immune, 100 * d$proportion, pch=19, col="grey30", cex=1.2,
     xlab="Biomarker level", ylab="Disease probability (%)", ylim=c(0,100),
     main="Fitted logistic curves: original versus log-biomarker scale")
lines(x_grid, 100 * pred_quad, col="#2563A6", lwd=3, lty=2)
lines(x_grid, 100 * pred_log, col="#E07A38", lwd=3)
legend("bottomright", legend=c("Observed proportions", "Quadratic / original scale", "Linear / log scale"),
       col=c("grey30", "#2563A6", "#E07A38"), pch=c(19,NA,NA),
       lty=c(NA,2,1), lwd=c(NA,3,3), bty="n")
dev.off()

png(file.path("figures", "in_sample_accuracy.png"), width=1200, height=750, res=140)
bar <- barplot(100 * comparison$accuracy, names.arg=c("Always predict\ndisease", "Quadratic\noriginal scale", "Linear\nlog scale"),
         ylim=c(0,100), col=c("grey65", "#2563A6", "#E07A38"),
         ylab="Correctly classified (%)", main="In-sample classification accuracy (n = 418)")
text(bar, 100 * comparison$accuracy + 3, sprintf("%.2f%%", 100 * comparison$accuracy))
mtext("Calculated on the same data used for fitting; no held-out evaluation", side=1, line=5, cex=0.78)
dev.off()

cat("Log-scale polynomial LRT summary:\n")
print(lrt_log, row.names=FALSE)
cat("Linear log-biomarker coefficients:\n")
print(coef(chosen))
cat(sprintf("Log-model threshold: %.5f; in-sample accuracy: %.2f%% (%d/%d).\n",
            threshold, 100 * acc_log$accuracy, acc_log$correct, acc_log$total))
cat(sprintf("Majority-class baseline: %.2f%%.\n", 100 * acc_majority$accuracy))
