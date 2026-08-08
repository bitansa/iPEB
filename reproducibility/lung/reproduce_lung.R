# Lung analysis workflow via library(iPEB), on the SYNTHETIC cohort.
# Mirrors the structure of the paper's lung tables: three like-for-like
# comparisons (four-marker, six-marker, selected-four) across the three
# objectives, on BOTH the held-out test centers and the training centers, plus
# the subject-level detection cross-classification (Web Table 6).
#
# The real PLCO data are controlled-access (NCI CDAS) and are NOT distributed.
# To reproduce the paper's actual numbers, point `lung` at your approved extract
# with the same columns; on the shipped synthetic data the numbers are
# illustrative only. First run make_synthetic_lung.R to create synthetic_lung.csv.
# Run:  source("reproduce_lung.R")

source("../simulations/sim_helpers.R")   # eval_score / metrics_subj / roc_auc + library(iPEB)

lung <- read.csv("synthetic_lung.csv", stringsAsFactors = FALSE)
markers6 <- c("CA125", "CEA", "CYFRA21", "PSFTPB", "OPN", "HE4")  # 4MP = first four
markers4 <- markers6[1:4]
covars   <- c("age", "batch", "study_year")
specs    <- c(0.60, 0.95, 0.99)
ALPHA    <- 0.60                          # triage operating point for fitting

train_centers <- c(2, 4, 5, 6, 8, 9)      # six training centers; the rest held out
tr_rows   <- lung$center %in% train_centers
train_ids <- unique(lung$id[tr_rows])
test_rows <- !tr_rows
ctr_train <- lung$id %in% train_ids & lung$case == 0   # training control visits (calibration)

## Evaluate a score on a chosen set of rows; the threshold is always calibrated
## on the training controls (as in the paper).
ev <- function(score, rows) t(sapply(specs, function(sp)
  eval_score(score[rows], lung$case[rows], lung$id[rows], lung$time_to_dx[rows], score[ctr_train], sp)))

## Frozen logistic panel scored through the PEB layer (the comparator).
peb_score <- function(markers) {
  g <- stats::glm(stats::as.formula(paste("case ~", paste(markers, collapse = "+"))),
                  family = binomial, data = lung[tr_rows, ])
  d <- lung; d$panel <- as.numeric(stats::predict(g, newdata = lung))
  ipeb_innovations(d, "panel", reference_ids = train_ids, id = "id", case = "case",
                   time = "time", covariates = covars, slope = "off", innovation = "ar1")[, 1]
}

## iPEB through the package; returns the full per-row score plus a label.
ipeb_score <- function(markers, objective, select = "none", n_markers = NULL) {
  fit <- ipeb(lung[tr_rows, ], markers = markers, id = "id", case = "case", time = "time",
              time_to_dx = "time_to_dx", covariates = covars, objective = objective,
              alpha = ALPHA, slope = "off", innovation = "ar1",
              select = select, n_markers = n_markers)
  list(score = predict(fit, lung), info = paste0(fit$variant, ": ", paste(fit$markers, collapse = "+")))
}

## Print the test-center and training-center tables for one score.
report <- function(score, label, info = NULL) {
  cat("\n---", label, if (!is.null(info)) paste0("[", info, "]") else "", "---\n")
  te <- ev(score, test_rows); tr <- ev(score, tr_rows)
  rownames(te) <- rownames(tr) <- paste0("spec ", specs)
  cat("  TEST centers:\n");  print(round(te, 3))
  cat("  TRAIN centers:\n"); print(round(tr, 3))
}

run_panel <- function(markers, name, sel = "none", nkeep = NULL, sel_markers = markers) {
  cat("\n===", name, "===")
  peb <- peb_score(sel_markers)
  report(peb, "PEB (frozen panel)")
  ip_sens <- ipeb_score(markers, "sensitivity", select = sel, n_markers = nkeep)
  report(ip_sens$score, "iPEB sensitivity", ip_sens$info)
  for (ob in c("leadtime", "combined")) {
    r <- ipeb_score(markers, ob, select = sel, n_markers = nkeep)
    report(r$score, paste("iPEB", ob), r$info)
  }
  list(peb = peb, ipeb = ip_sens$score)          # kept for the paired cross-classification
}

four <- run_panel(markers4, "Four-marker (frozen 4MP vs iPEB)")
six  <- run_panel(markers6, "Six-marker (frozen logistic vs iPEB)")

# selected-four: PEB ranks the six by |z| and keeps the top four; iPEB uses backward selection.
gfull <- stats::glm(stats::as.formula(paste("case ~", paste(markers6, collapse = "+"))),
                    family = binomial, data = lung[tr_rows, ])
sel4 <- names(sort(abs(summary(gfull)$coefficients[markers6, "z value"]), decreasing = TRUE))[1:4]
cat("\n(PEB logistic-selected four:", paste(sel4, collapse = "+"), ")")
sel  <- run_panel(markers6, "Selected-four (each reduces six to four)",
                  sel = "backward", nkeep = 4, sel_markers = sel4)

## ---- subject-level detection cross-classification (Web Table 6) -------------
## Sensitivity-objective iPEB vs the frozen panel: at each specificity, count
## test/train CASES detected by both, by iPEB only, by the panel only, or neither.
detect_cases <- function(score, thr, rows) {
  cids <- unique(lung$id[rows & lung$case == 1])
  vapply(cids, function(id) {
    ix <- which(lung$id == id & rows & lung$time_to_dx > 0)
    any(score[ix] > thr, na.rm = TRUE)
  }, logical(1))
}
paired <- function(cmp, comp_name) {
  for (setname in c("Test", "Train")) {
    rows <- if (setname == "Test") test_rows else tr_rows
    for (sp in specs) {
      thr_i <- stats::quantile(cmp$ipeb[ctr_train], sp, na.rm = TRUE, type = 8)
      thr_p <- stats::quantile(cmp$peb[ctr_train],  sp, na.rm = TRUE, type = 8)
      di <- detect_cases(cmp$ipeb, thr_i, rows); dp <- detect_cases(cmp$peb, thr_p, rows)
      cat(sprintf("  %-14s %-5s spec %.2f | both %3d | iPEB-only %3d | panel-only %3d | neither %3d\n",
                  comp_name, setname, sp, sum(di & dp), sum(di & !dp), sum(!di & dp), sum(!di & !dp)))
    }
  }
}
cat("\n\n=== Subject-level detection cross-classification (iPEB sensitivity vs panel) ===\n")
paired(four, "Four-marker"); paired(six, "Six-marker"); paired(sel, "Selected-four")

cat("\nColumns of the tables: AUC, SensW, LT (median lead, years); one row per specificity.\n")
cat("Note: synthetic data -- illustrative workflow, not the paper's numbers.\n")
