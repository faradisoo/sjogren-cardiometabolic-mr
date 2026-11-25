###############################################################################
# 04. POST HOC POWER - MINIMUM DETECTABLE EFFECTS
#
#   Minimum detectable effect (MDE) at 80% power and the Bonferroni alpha,
#   from the standard error of each primary estimate (non-centrality approach,
#   Brion et al. 2013):   MDE = (z[1 - alpha/2] + z[0.80]) x SE
#
#   Analysis A: per doubling of the genetically predicted odds of pSS (x ln 2),
#               exponentiated (odds ratio) for coronary heart disease  -> Table 5
#   Analysis B: odds ratio for pSS per SD of genetically predicted expression
#
# Input : results/MR/A_results_allSNPs.csv, A_results_noHLA.csv, B_results.csv
#         (written by script 02)
# Output: results/MR/Power_MDE_A.csv, results/MR/Power_MDE_B.csv
###############################################################################

out <- function(f) file.path("results", "MR", f)
multiplier <- function(alpha, power = 0.80) qnorm(1 - alpha / 2) + qnorm(power)

units <- c("ieu-a-299" = "SD", "ieu-a-300" = "SD", "ieu-a-302" = "SD",
           "ieu-b-35" = "ln mg/L", "ieu-b-38" = "mmHg", "ieu-b-39" = "mmHg",
           "ieu-a-7" = "OR")


# =============================================================================
# Analysis A: genetic liability to pSS -> 7 cardiometabolic traits (Table 5)
# =============================================================================
mde_A <- function(file) {
  res   <- read.csv(out(file))
  res   <- res[res$method == "Inverse variance weighted", ]
  alpha <- 0.05 / nrow(res)                    # 7 outcomes -> 0.0071
  m     <- multiplier(alpha)                   # 3.53
  mde   <- m * res$se * log(2)                 # per doubling of the odds of pSS
  binary <- res$id.outcome == "ieu-a-7"        # coronary heart disease
  data.frame(outcome    = sub(" \\|\\|.*$", "", res$outcome),
             unit       = unname(units[res$id.outcome]),
             MDE        = ifelse(binary, exp(mde), mde),
             alpha      = alpha,
             multiplier = m)
}

A_all   <- mde_A("A_results_allSNPs.csv")
A_noHLA <- mde_A("A_results_noHLA.csv")
power_A <- merge(A_all[, c("outcome", "unit", "MDE")],
                 A_noHLA[, c("outcome", "MDE")],
                 by = "outcome", suffixes = c("_allSNPs", "_HLA_excluded"))
print(power_A, digits = 3)
# 19 Nov 2026 results: lipids 0.014-0.033 SD, blood pressure 0.12-0.23 mmHg,
# coronary heart disease OR 1.033 (all SNPs)
write.csv(power_A, out("Power_MDE_A.csv"), row.names = FALSE)


# =============================================================================
# Analysis B: blood expression of B-cell genes -> pSS
# =============================================================================
B <- read.csv(out("B_results.csv"))
B <- B[B$method %in% c("Inverse variance weighted", "Wald ratio"), ]
n_tests <- length(unique(B$exposure)) * length(unique(B$id.outcome))   # 5 genes x 8 outcomes = 40
alpha_B <- 0.05 / n_tests                                             # 0.00125
m_B     <- multiplier(alpha_B)                                        # 4.07
B_pss   <- B[B$id.outcome == "ebi-a-GCST90013879", ]
power_B <- data.frame(gene   = B_pss$exposure,
                      nsnp   = B_pss$nsnp,
                      se     = B_pss$se,
                      MDE_OR = exp(m_B * B_pss$se))
power_B <- power_B[order(power_B$MDE_OR), ]
print(power_B, digits = 3)
# 19 Nov 2026 results: SELL 1.41, CR2 6.5, BANK1 14.8, PTPRC 30.7, MS4A1 79.6
write.csv(power_B, out("Power_MDE_B.csv"), row.names = FALSE)

cat(sprintf("Bonferroni alpha: A = %.4f (multiplier %.2f); B = %.5f (multiplier %.2f)\n",
            A_all$alpha[1], A_all$multiplier[1], alpha_B, m_B))
