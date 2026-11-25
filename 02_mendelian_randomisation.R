###############################################################################
# 02. MENDELIAN RANDOMISATION - SJÖGREN'S SYNDROME AND CARDIOMETABOLIC TRAITS
#
#   ANALYSIS A    : genetic liability to Sjögren's syndrome -> 7 cardiometabolic traits
#   ANALYSIS B    : blood expression of B-cell genes (cis-eQTL, eQTLGen)
#                   -> Sjögren's syndrome and the same 7 traits
#   COLOCALISATION: for every analysis B finding that passes Bonferroni correction
#
# Requires an OpenGWAS token in ~/.Renviron (OPENGWAS_JWT=...), see README.
# Output: results/MR/
# Results reported in the paper were generated on 19 November 2026.
###############################################################################


# =============================================================================
# 0. PACKAGES AND SETTINGS
# =============================================================================
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(org.Hs.eg.db)
  library(coloc)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
# Bioconductor packages mask some dplyr verbs; restore them
select <- dplyr::select
first  <- dplyr::first

out_dir <- "results/MR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(out_dir, f)

ieugwasr::user()   # checks the token; a 401 error means the token is missing or expired

# Weighted-median standard errors are bootstrapped; fix the seed so reruns give identical P values
set.seed(2026)

sjogren_id <- "ebi-a-GCST90013879"   # Sjögren's syndrome (Firth), UK Biobank, Mbatchou 2021

cardio_outcomes <- c(
  LDL = "ieu-a-300",   # GLGC, Willer 2013
  HDL = "ieu-a-299",   # GLGC, Willer 2013
  TG  = "ieu-a-302",   # GLGC, Willer 2013
  CRP = "ieu-b-35",    # CHARGE, Ligthart 2018
  CAD = "ieu-a-7",     # CARDIoGRAMplusC4D, Nikpay 2015
  SBP = "ieu-b-38",    # ICBP, Evangelou 2018 (includes UK Biobank -> sample overlap)
  DBP = "ieu-b-39"     # ICBP, Evangelou 2018 (includes UK Biobank -> sample overlap)
)
binary_outcomes <- c(sjogren_id, "ieu-a-7")

# MR engine: 4 estimators + heterogeneity + pleiotropy + Steiger directionality
run_mr <- function(dat) {
  dat <- dat %>% filter(mr_keep)
  list(
    dat     = dat,
    results = generate_odds_ratios(mr(dat, method_list = c(
                "mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))),
    het     = mr_heterogeneity(dat),
    pleio   = mr_pleiotropy_test(dat),
    steiger = tryCatch(directionality_test(dat), error = function(e) NULL)  # needs 'psych'
  )
}
clean_name <- function(x) sub(" \\|\\|.*", "", x)


###############################################################################
# ANALYSIS A: SJÖGREN'S SYNDROME -> CARDIOMETABOLIC TRAITS
###############################################################################

# A1. Instruments for Sjögren's syndrome
# At P < 5e-8 only 3 SNPs remain, all in the HLA region, so the threshold is
# relaxed to 5e-6 (common for rarer autoimmune diseases; stated in the Methods).
inst_sjo <- extract_instruments(sjogren_id, p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000)
inst_sjo$F_stat <- (inst_sjo$beta.exposure / inst_sjo$se.exposure)^2
inst_sjo <- inst_sjo %>% filter(F_stat > 10)
inst_sjo$HLA <- inst_sjo$chr.exposure == "6" &
                as.numeric(inst_sjo$pos.exposure) > 25e6 &
                as.numeric(inst_sjo$pos.exposure) < 34e6
table(inst_sjo$HLA)   # 19 Nov 2026: 3 HLA, 7 non-HLA

# A2. Outcome data and harmonisation
out_A <- extract_outcome_data(inst_sjo$SNP, cardio_outcomes, proxies = TRUE, rsq = 0.8)
dat_A <- harmonise_data(inst_sjo, out_A, action = 2)

# A3. MR with all SNPs, and excluding HLA variants (sensitivity analysis)
A_all   <- run_mr(dat_A)
A_noHLA <- run_mr(dat_A %>% filter(!SNP %in% inst_sjo$SNP[inst_sjo$HLA]))
threshold_A <- 0.05 / length(cardio_outcomes)   # Bonferroni = 0.0071

A_all$results   %>% filter(method == "Inverse variance weighted") %>% select(outcome, nsnp, b, se, pval)
A_noHLA$results %>% filter(method == "Inverse variance weighted") %>% select(outcome, nsnp, b, se, pval)

# A4. Single-SNP estimates (used to identify the heterogeneity source, rs1811197 in HLA)
A_single <- mr_singlesnp(A_all$dat) %>%
  filter(!grepl("^All", SNP)) %>%
  mutate(HLA = SNP %in% inst_sjo$SNP[inst_sjo$HLA])

# A5. Main table (one row per outcome)
fmt <- function(b, se) sprintf("%.3f (%.3f, %.3f)", b, b - 1.96 * se, b + 1.96 * se)
summarise_ivw <- function(res, label) {
  res$results %>%
    filter(method == "Inverse variance weighted") %>%
    select(id.outcome, outcome, nsnp, b, se, pval) %>%
    left_join(res$het %>% filter(method == "Inverse variance weighted") %>%
                select(id.outcome, Q_pval), by = "id.outcome") %>%
    left_join(res$pleio %>% select(id.outcome, egger_p = pval), by = "id.outcome") %>%
    transmute(Outcome = clean_name(outcome),
              !!paste(label, "nSNP") := nsnp,
              !!paste(label, "Beta (95% CI)") := fmt(b, se),
              !!paste(label, "P") := signif(pval, 2),
              !!paste(label, "P het") := signif(Q_pval, 2),
              !!paste(label, "P Egger int") := signif(egger_p, 2))
}
table_A <- left_join(summarise_ivw(A_all, "All SNPs"),
                     summarise_ivw(A_noHLA, "HLA excluded"), by = "Outcome")
table_A

# A6. Save
write.csv(inst_sjo,          out("A_instruments_Sjogren.csv"),   row.names = FALSE)
write.csv(A_all$results,     out("A_results_allSNPs.csv"),       row.names = FALSE)
write.csv(A_noHLA$results,   out("A_results_noHLA.csv"),         row.names = FALSE)
write.csv(A_all$het,         out("A_heterogeneity_allSNPs.csv"), row.names = FALSE)
write.csv(A_noHLA$het,       out("A_heterogeneity_noHLA.csv"),   row.names = FALSE)
write.csv(A_all$pleio,       out("A_pleiotropy_allSNPs.csv"),    row.names = FALSE)
write.csv(A_noHLA$pleio,     out("A_pleiotropy_noHLA.csv"),      row.names = FALSE)
write.csv(A_single,          out("A_single_SNP.csv"),            row.names = FALSE)
if (!is.null(A_all$steiger))
  write.csv(A_all$steiger,   out("A_steiger.csv"),               row.names = FALSE)
write.csv(table_A,           out("Table_A_main.csv"),            row.names = FALSE)

# Forest plot: run 03_forest_plot_MR_A.R


###############################################################################
# ANALYSIS B: B-CELL GENE EXPRESSION (cis-eQTL) -> SJÖGREN'S SYNDROME & CARDIO
###############################################################################

# B1. B-cell genes up-regulated in pSS (01_transcriptomics.R) -> Ensembl -> eQTLGen ID
mr_genes <- c("IGHM", "MS4A1", "CXCL13", "SELL", "UBD", "CR2", "BANK1", "PAX5", "PTPRC")
ensg <- mapIds(org.Hs.eg.db, keys = mr_genes, keytype = "SYMBOL",
               column = "ENSEMBL", multiVals = "first")
eqtl_ids  <- setNames(paste0("eqtl-a-", ensg), names(ensg))
available <- gwasinfo(eqtl_ids)$id
eqtl_ids  <- eqtl_ids[eqtl_ids %in% available]
eqtl_ids   # CXCL13 and UBD are not available in eQTLGen

# B2. eQTL instruments (P < 5e-8, F > 10)
inst_eqtl <- extract_instruments(eqtl_ids, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
inst_eqtl$F_stat <- (inst_eqtl$beta.exposure / inst_eqtl$se.exposure)^2
inst_eqtl <- inst_eqtl %>% filter(F_stat > 10)
inst_eqtl$gene <- names(eqtl_ids)[match(inst_eqtl$id.exposure, eqtl_ids)]
inst_eqtl$exposure <- inst_eqtl$gene
# IGHM and PAX5 have no SNPs passing the thresholds and drop out

# B3. Keep cis-eQTLs only (within 1 Mb of the gene)
# eQTLGen in OpenGWAS also contains trans-eQTLs (e.g. 12q24 SH2B3/ATXN2, CD40,
# ORMDL3/GSDMB, HLA). These are pleiotropic and produce spurious associations.
gene_loc <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(inst_eqtl$gene),
                                  keytype = "SYMBOL", columns = c("CHRLOC", "CHRLOCEND")) %>%
  filter(CHRLOCCHR %in% c(1:22, "X")) %>%
  group_by(SYMBOL) %>%
  summarise(gene_chr   = first(CHRLOCCHR),
            gene_start = min(abs(CHRLOC)),
            gene_end   = max(abs(CHRLOCEND)), .groups = "drop")
inst_eqtl <- inst_eqtl %>%
  left_join(gene_loc, by = c("gene" = "SYMBOL")) %>%
  mutate(cis = chr.exposure == gene_chr &
               as.numeric(pos.exposure) > gene_start - 1e6 &
               as.numeric(pos.exposure) < gene_end + 1e6)
table(inst_eqtl$gene, inst_eqtl$cis)

# B4. Confirm that cis SNPs within each gene are independent
# Rule: if any SNP is missing from the LD reference panel, or any pair has r2 > 0.01,
#       the gene uses only its strongest eQTL SNP (Wald ratio).
check_ld <- function(snps) {
  if (length(snps) < 2) return(TRUE)
  m <- suppressWarnings(tryCatch(ld_matrix(snps, pop = "EUR"), error = function(e) NULL))
  if (is.null(m) || is.null(dim(m)) || nrow(m) < length(snps)) return(FALSE)
  r2 <- m^2; diag(r2) <- 0
  all(r2 < 0.01)
}
cis_snp <- inst_eqtl %>% filter(cis)
ld_ok <- sapply(split(cis_snp$SNP, cis_snp$gene), check_ld)
ld_ok   # TRUE = independent SNPs; FALSE = lead SNP only

snp_final <- cis_snp %>%
  group_by(gene) %>%
  filter(ld_ok[gene[1]] | pval.exposure == min(pval.exposure)) %>%
  ungroup()
snp_final %>% count(gene)   # 19 Nov 2026: BANK1 1, CR2 1, MS4A1 1, PTPRC 1, SELL 5

# B5. Outcomes: Sjögren's syndrome + 7 cardiometabolic traits
outcome_B <- c(Sjogren = sjogren_id, cardio_outcomes)
out_B <- extract_outcome_data(unique(snp_final$SNP), outcome_B, proxies = TRUE, rsq = 0.8)
dat_B <- harmonise_data(snp_final, out_B, action = 2)
B_res <- run_mr(dat_B)

threshold_B <- 0.05 / (length(unique(snp_final$gene)) * length(outcome_B))   # 0.00125
B_main <- B_res$results %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  mutate(outcome = clean_name(outcome), pass_bonferroni = pval < threshold_B) %>%
  arrange(pval)
B_main %>% select(exposure, outcome, method, nsnp, b, se, pval, pass_bonferroni) %>% head(10)

write.csv(inst_eqtl,       out("B_instruments_cis_and_trans.csv"), row.names = FALSE)
write.csv(snp_final,       out("B_instruments_final_cis.csv"),     row.names = FALSE)
write.csv(B_res$results,   out("B_results.csv"),                   row.names = FALSE)
write.csv(B_res$het,       out("B_heterogeneity.csv"),             row.names = FALSE)


###############################################################################
# COLOCALISATION for analysis B findings passing Bonferroni (continuous outcomes)
###############################################################################
# Tests whether the eQTL and outcome signals share one causal variant (PP.H4)
# or are driven by distinct nearby variants (PP.H3 = linkage confounding).
signif_B <- B_main %>% filter(pass_bonferroni, !id.outcome %in% binary_outcomes)

# Region queries are large; retry if the OpenGWAS server returns a transient error
get_region <- function(region, id, tries = 3) {
  for (k in seq_len(tries)) {
    res <- tryCatch(associations(variants = region, id = id, proxies = 0),
                    error = function(e) e)
    if (!inherits(res, "error")) break
    message("OpenGWAS query failed (", conditionMessage(res), "); retrying in 30 s")
    Sys.sleep(30)
  }
  if (inherits(res, "error")) stop(res)
  res %>%
    filter(!is.na(beta), !is.na(se), se > 0) %>%
    distinct(rsid, .keep_all = TRUE)
}
merge_for_coloc <- function(e, g) {
  inner_join(e, g, by = "rsid", suffix = c(".e", ".g")) %>%
    filter((ea.e == ea.g & nea.e == nea.g) | (ea.e == nea.g & nea.e == ea.g)) %>%
    mutate(beta.e = ifelse(ea.e == ea.g, beta.e, -beta.e),
           eaf    = ifelse(is.na(eaf.g), eaf.e, eaf.g),
           MAF    = pmin(eaf, 1 - eaf)) %>%
    filter(!is.na(MAF), MAF > 0)
}
run_coloc <- function(m, N_eqtl, N_gwas, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {
  res <- suppressWarnings(coloc.abf(
    dataset1 = list(snp = m$rsid, beta = m$beta.e, varbeta = m$se.e^2,
                    type = "quant", N = N_eqtl, MAF = m$MAF),
    dataset2 = list(snp = m$rsid, beta = m$beta.g, varbeta = m$se.g^2,
                    type = "quant", N = N_gwas, MAF = m$MAF),
    p1 = p1, p2 = p2, p12 = p12))
  res$summary
}

# Prior sensitivity: default coloc priors plus three alternatives
priors <- list(
  "default (p1 = p2 = 1e-4, p12 = 1e-5)" = c(1e-4, 1e-4, 1e-5),
  "p12 = 1e-6"                           = c(1e-4, 1e-4, 1e-6),
  "p12 = 1e-4"                           = c(1e-4, 1e-4, 1e-4),
  "p1 = p2 = 1e-3, p12 = 1e-4"           = c(1e-3, 1e-3, 1e-4)
)

N_EQTLGEN <- 31684
coloc_results <- NULL
coloc_priors  <- NULL
for (i in seq_len(nrow(signif_B))) {
  row <- signif_B[i, ]
  lead_snp <- snp_final %>% filter(gene == row$exposure) %>% slice_min(pval.exposure, n = 1)
  pos    <- as.numeric(lead_snp$pos.exposure)
  region <- paste0(lead_snp$chr.exposure, ":", pos - 1e6, "-", pos + 1e6)
  e <- get_region(region, row$id.exposure)
  g <- get_region(region, row$id.outcome)
  m <- merge_for_coloc(e, g)
  N_gwas <- as.numeric(gwasinfo(row$id.outcome)$sample_size)
  lead_outcome <- g %>% slice_min(p, n = 1)

  # LD between the eQTL instrument and the lead outcome SNP in the region
  ld <- tryCatch(ld_matrix(c(lead_snp$SNP, lead_outcome$rsid), pop = "EUR"),
                 error = function(e) NULL)
  r2_lead <- if (!is.null(ld) && all(dim(ld) == 2)) round(ld[1, 2]^2, 3) else NA

  s <- run_coloc(m, N_eqtl = N_EQTLGEN, N_gwas = N_gwas)
  coloc_results <- bind_rows(coloc_results, data.frame(
    gene = row$exposure, outcome = row$outcome, region = region,
    eqtl_snp = lead_snp$SNP, lead_outcome_snp = lead_outcome$rsid,
    p_lead_outcome = lead_outcome$p, r2_eqtl_vs_lead = r2_lead,
    t(round(s, 4)), check.names = FALSE))

  for (pr in names(priors)) {
    s_pr <- run_coloc(m, N_EQTLGEN, N_gwas, priors[[pr]][1], priors[[pr]][2], priors[[pr]][3])
    coloc_priors <- bind_rows(coloc_priors, data.frame(
      gene = row$exposure, outcome = row$outcome, prior = pr,
      t(c(s_pr[1], signif(s_pr[-1], 3))), check.names = FALSE))
  }
}
coloc_results   # 19 Nov 2026: BANK1-DBP and BANK1-SBP PP.H3 = 1.00 (lead SNP rs13107325, SLC39A8)
coloc_priors    # PP.H3 = 1.00 under every prior
if (!is.null(coloc_results)) {
  write.csv(coloc_results, out("B_coloc.csv"),                   row.names = FALSE)
  write.csv(coloc_priors,  out("B_coloc_prior_sensitivity.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
