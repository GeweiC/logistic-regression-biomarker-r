# Shared functions for grouped-binomial modelling.
# Run all scripts from the directory containing the top-level README.md.

read_immune_data <- function(path = file.path("data", "immune.csv")) {
  if (!file.exists(path)) {
    stop("Missing data/immune.csv. See data/README.md; keep your coursework CSV local.")
  }
  d <- read.csv(path, check.names = FALSE)
  required <- c("immune", "disease", "nodisease")
  if (!all(required %in% names(d))) {
    stop("CSV must contain: ", paste(required, collapse = ", "))
  }
  d <- d[, required]
  if (anyNA(d) || any(!is.finite(as.matrix(d))) ||
      any(d$immune <= 0) || any(d$disease < 0) || any(d$nodisease < 0) ||
      any(d$disease != floor(d$disease)) ||
      any(d$nodisease != floor(d$nodisease)) ||
      any(d$disease + d$nodisease == 0)) {
    stop("Invalid biomarker or grouped outcome counts in immune.csv")
  }
  if (anyDuplicated(d$immune)) stop("Biomarker levels must be unique")
  d <- d[order(d$immune), ]
  rownames(d) <- NULL
  d$sample_size <- d$disease + d$nodisease
  d$proportion <- d$disease / d$sample_size
  d$log_immune <- log(d$immune)
  d
}

prepare_output_dirs <- function() {
  dir.create("figures", recursive = TRUE, showWarnings = FALSE)
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
}

fit_polynomials <- function(d, predictor = c("immune", "log_immune"), degrees = 1:5) {
  predictor <- match.arg(predictor)
  lapply(degrees, function(k) {
    f <- as.formula(sprintf(
      "cbind(disease, nodisease) ~ poly(%s, %d, raw = TRUE)", predictor, k
    ))
    glm(f, family = binomial(link = "logit"), data = d)
  })
}

nested_lrt_table <- function(models) {
  n <- length(models)
  lr <- rep(NA_real_, n)
  p <- rep(NA_real_, n)
  for (i in seq.int(2L, n)) {
    lr[i] <- 2 * as.numeric(logLik(models[[i]]) - logLik(models[[i - 1L]]))
    df_difference <- attr(logLik(models[[i]]), "df") -
                     attr(logLik(models[[i - 1L]]), "df")
    p[i] <- pchisq(lr[i], df = df_difference, lower.tail = FALSE)
  }
  data.frame(degree = seq_len(n), LR_statistic = lr, p_value = p,
             AIC = vapply(models, AIC, numeric(1)))
}

threshold_from_log_model <- function(model) {
  b <- unname(coef(model))
  if (length(b) != 2L || abs(b[2]) < 1e-12) stop("Expected a linear log-biomarker model")
  exp(-b[1] / b[2])
}

threshold_from_quadratic <- function(model, observed_range) {
  b <- unname(coef(model))
  if (length(b) != 3L) stop("Expected a quadratic model")
  roots <- polyroot(b)
  candidates <- Re(roots)[abs(Im(roots)) < 1e-7 &
                          Re(roots) >= min(observed_range) &
                          Re(roots) <= max(observed_range)]
  if (length(candidates) != 1L) {
    stop("Expected one in-range quadratic threshold; inspect roots")
  }
  as.numeric(candidates[1])
}

# In-sample accuracy on grouped counts. NOT a held-out test-set accuracy.
grouped_classification <- function(d, threshold, label) {
  no_disease_correct <- sum(d$nodisease[d$immune < threshold])
  disease_correct <- sum(d$disease[d$immune >= threshold])
  total <- sum(d$sample_size)
  data.frame(method = label, threshold = threshold,
             true_negative = no_disease_correct,
             true_positive = disease_correct,
             correct = no_disease_correct + disease_correct,
             total = total,
             accuracy = (no_disease_correct + disease_correct) / total)
}
