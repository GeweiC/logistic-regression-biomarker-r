# Methodology and provenance

## Source reconstruction

This portfolio is reconstructed from a completed MATH3091 university coursework report, which contains R script excerpts and computed outputs, plus the originally supplied `immune.csv` grouped data. The separate R Markdown file was a largely unfilled **coursework template**, not the source of the full final solution. The organised `.R` scripts in this repository were written later for reproducibility and portfolio use; they are **not represented as the original submitted code files**.

## Grouped-binomial model

The data consist of 12 biomarker levels and counts of diseased/non-diseased individuals. Each row has a binomial response `cbind(disease, nodisease)`, whose modelled probability is:

`logit(p(x)) = beta_0 + beta_1 * x + ... + beta_k * x^k`

For positive biomarker levels, the log-transformed alternative is:

`logit(p(x)) = gamma_0 + gamma_1 * log(x) + ... + gamma_k * log(x)^k`

The original report compared degrees 1 through 5 with **sequential one-degree likelihood-ratio tests** at a 5% significance level.

- On the original biomarker scale, degree 2 improves significantly upon degree 1 (LR ~5.2047, p ~0.0225). Moving from degree 2 to 3 does not (p ~0.1231). The selected original-scale model is quadratic.
- On the log scale, adding the quadratic term to the linear log model does not clear 5% significance (p ~0.1434); the report selected the parsimonious linear log model.
- Selection is **not claimed to minimise AIC globally**. In fact, the log-quadratic model has a marginally lower AIC (approximately 48.762) than the chosen linear log model (48.903). The choice follows the report's sequential LRT and simplicity rationale.

## Classification and limitations

A 0.5 fitted-probability threshold for the selected linear log-biomarker model gives a biomarker cutoff of approximately 42.281. With observations classified using biomarker >= cutoff as positive, the grouped **training-set** accuracy is 322/418 = 77.03%. A rule that always predicts disease correctly classifies 266/418 = 63.64%. A quadratic model on the original scale reaches 321/418 = 76.79%.

These are **apparent (in-sample) accuracies**, not external or held-out test accuracy; model selection and threshold selection both used the same data. Without patient-level information, no individual train/test split has been demonstrated. The percentage difference should not be called statistically significant without a suitable additional test. The fitted association does not establish disease causation or clinical diagnostic validity.

## Verification status

During portfolio preparation, the original grouped counts were independently fitted with Python `statsmodels` binomial GLMs. The coefficient estimates, LRT p-values, threshold, and all in-sample accuracy counts matched the submitted report to displayed precision. The accompanying R code has **not yet been run in this packaging environment** (R is unavailable here). The project owner should run `source("tests/basic_checks.R")` and `source("R/run_all.R")` in RStudio before describing the reorganised R implementation as execution-tested.

## Data licence

The data were supplied for coursework; the right to redistribute the CSV has not been verified. The original assessment document and student-number-bearing report are excluded from the package. Included PNG charts depict model fits and grouped summaries only, not individual records. Publication of any derivative materials should still follow applicable university/data-owner requirements.
