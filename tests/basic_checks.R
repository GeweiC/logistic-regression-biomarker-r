# Basic tests: no external R packages required.
# Run from project root: source("tests/basic_checks.R")
source(file.path("R", "helpers.R"), local=TRUE)
dummy <- data.frame(immune=c(10, 50), disease=c(2, 8), nodisease=c(8, 2))
dummy$sample_size <- dummy$disease + dummy$nodisease
acc <- grouped_classification(dummy, 40, "Synthetic example")
stopifnot(acc$correct == 16, acc$total == 20, abs(acc$accuracy - 0.8) < 1e-12)
stopifnot(abs(plogis(0) - 0.5) < 1e-12)
cat("Basic grouped-classification checks passed.\n")

# When authorised coursework CSV is present, also verify key submitted results.
if (file.exists(file.path("data", "immune.csv"))) {
  d <- read_immune_data()
  stopifnot(nrow(d) == 12L, sum(d$sample_size) == 418L)
  raw <- fit_polynomials(d, "immune")
  log_models <- fit_polynomials(d, "log_immune")
  stopifnot(abs(nested_lrt_table(raw)$LR_statistic[2] - 5.2047072) < 1e-3)
  stopifnot(abs(nested_lrt_table(log_models)$p_value[2] - 0.1434061) < 1e-4)
  t <- threshold_from_log_model(log_models[[1]])
  stopifnot(abs(t - 42.28117) < 0.02)
  r <- grouped_classification(d, t, "Validation")
  stopifnot(r$correct == 322L, r$total == 418L)
  cat("Key coursework model results matched.\n")
} else {
  cat("Source CSV not present: independent formula tests passed; dataset checks skipped.\n")
}
