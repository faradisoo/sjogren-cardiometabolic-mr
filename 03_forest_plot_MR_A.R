###############################################################################
# 03. FOREST PLOT FOR MR ANALYSIS A (supplementary figure)
#
# Input : results/MR/A_results_allSNPs.csv and A_results_noHLA.csv
#         (created by 02_mendelian_randomisation.R)
# Output: results/MR/A_forest.png
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

colours <- c("All SNPs" = "#2a78d6", "HLA excluded" = "#eb6834")
shapes  <- c("All SNPs" = 16,        "HLA excluded" = 15)

labels_out <- c(
  "ieu-a-300" = "LDL cholesterol",
  "ieu-a-299" = "HDL cholesterol",
  "ieu-a-302" = "Triglycerides",
  "ieu-b-35"  = "C-reactive protein",
  "ieu-b-38"  = "Systolic blood pressure",
  "ieu-b-39"  = "Diastolic blood pressure",
  "ieu-a-7"   = "Coronary heart disease"
)
group_out <- c(
  "ieu-a-300" = "Lipids (SD units)", "ieu-a-299" = "Lipids (SD units)",
  "ieu-a-302" = "Lipids (SD units)",
  "ieu-b-35"  = "Inflammation (log CRP)",
  "ieu-b-38"  = "Blood pressure (mmHg)", "ieu-b-39" = "Blood pressure (mmHg)",
  "ieu-a-7"   = "Coronary heart disease (odds ratio)"
)

d <- bind_rows(
  read.csv("results/MR/A_results_allSNPs.csv") %>% mutate(analysis = "All SNPs"),
  read.csv("results/MR/A_results_noHLA.csv")   %>% mutate(analysis = "HLA excluded")
) %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    label    = labels_out[id.outcome],
    group    = group_out[id.outcome],
    binary   = id.outcome == "ieu-a-7",
    est      = ifelse(binary, exp(b), b),
    lo       = ifelse(binary, exp(b - 1.96 * se), b - 1.96 * se),
    hi       = ifelse(binary, exp(b + 1.96 * se), b + 1.96 * se),
    p_txt    = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    txt      = sprintf("  %6.3f (%6.3f, %6.3f)   %s", est, lo, hi, p_txt),
    analysis = factor(analysis, levels = names(colours))
  )

panel <- function(grp, ref = 0, header = FALSE) {
  dd <- d %>% filter(group == grp)
  lv <- rev(unname(labels_out[names(labels_out) %in% dd$id.outcome]))
  dd <- dd %>%
    mutate(y = match(label, lv) + ifelse(analysis == "All SNPs", 0.17, -0.17))

  p <- ggplot(dd, aes(y = y, colour = analysis, shape = analysis)) +
    geom_vline(xintercept = ref, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.7) +
    geom_point(aes(x = est), size = 2.6) +
    geom_text(aes(x = Inf, label = txt), hjust = 0, family = "mono",
              size = 2.9, colour = "grey15", show.legend = FALSE) +
    scale_y_continuous(breaks = seq_along(lv), labels = lv,
                       limits = c(0.5, length(lv) + 0.5), expand = c(0, 0)) +
    scale_colour_manual(values = colours, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = grp, colour = NULL, shape = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor    = element_blank(),
          panel.grid.major.y  = element_blank(),
          axis.text.y         = element_text(colour = "grey15"),
          plot.title          = element_text(face = "bold", size = 10),
          plot.title.position = "plot",
          plot.margin         = margin(4, 190, 4, 4))
  if (ref == 1) p <- p + scale_x_log10()
  if (header) p <- p + annotate("text", x = Inf, y = length(lv) + 0.5,
                                label = "  Estimate (95% CI)        p",
                                hjust = 0, vjust = -0.6, family = "mono",
                                fontface = "bold", size = 2.9, colour = "grey15")
  p
}

p_A <- panel("Lipids (SD units)", header = TRUE) /
  panel("Inflammation (log CRP)") /
  panel("Blood pressure (mmHg)") /
  panel("Coronary heart disease (odds ratio)", ref = 1) +
  plot_layout(heights = c(3, 1, 2, 1), guides = "collect") +
  plot_annotation(
    title    = "Genetic liability to Sjögren's syndrome and cardiometabolic traits",
    subtitle = "Two-sample MR (inverse variance weighted), effect per 1-unit increase in log-odds of Sjögren's syndrome",
    caption  = "Bonferroni threshold p < 0.0071 (0.05/7). HLA excluded: variants on chr6 25–34 Mb removed.",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12.5),
                     plot.subtitle = element_text(size = 9.5, colour = "grey30"),
                     plot.caption  = element_text(size = 8.5, colour = "grey40", hjust = 0))
  ) &
  theme(legend.position = "bottom")

ggsave("results/MR/A_forest.png", p_A, width = 10, height = 7, dpi = 300, bg = "white")
