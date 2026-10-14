# From repository root in RStudio: source("R/run_all.R")
if (!file.exists(file.path("data", "immune.csv"))) {
  stop("Add your authorised immune.csv to data/ before running; see data/README.md")
}
for (script in c("01_explore_data.R", "02_polynomial_logistic.R",
                 "03_log_biomarker_model.R")) {
  cat("\n== Running ", script, " ==\n", sep="")
  source(file.path("R", script), local=new.env(parent=globalenv()))
}
cat("\nDone: inspect results/ and figures/.\n")
