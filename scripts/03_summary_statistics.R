# =============================================================================
# Script 08: Summary statistics, missing data & figures
# =============================================================================
# Input:
#   chem_avg_flagged.csv  — output of Script 2 (averaged, flagged, all sites)
#   Missing QuEST data CSV — flagged missing samples tracking sheet
#
# All cleaning, averaging, rep variability, and outlier flagging was already
# done in Script 2. This script loads the final averaged dataset and produces
# summary statistics, missing data summaries, and all figures.
#
# Outputs:
#   summary_stats_by_watershed.csv
#   summary_stats_by_site.csv
#   summary_stats_by_sample_type.csv
#   summary_stats_by_date.csv
#   summary_stats_by_event.csv
#   outlier_records.csv
#   rep_variability_summary.csv
#   figures/ (all plots)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(viridis)
library(googledrive)

dir.create("figures", showWarnings = FALSE)

# ── 0. CONFIGURATION ──────────────────────────────────────────────────────────

RPD_THRESHOLD <- 15    # matches Script 2

solutes <- c(
  "NPOC..mg.C.L.", "TDN..mg.N.L.", "NH4..ug.N.L.", "PO4..ug.P.L.",
  "Cl..mg.Cl.L.",  "NO3..mg.N.L.", "SO4..mg.S.L.", "Na..mg.Na.L.",
  "K..mg.K.L.",    "Mg..mg.Mg.L.", "Ca..mg.Ca.L."
)

solute_labels <- c(
  "NPOC..mg.C.L." = "NPOC (mg C/L)", "TDN..mg.N.L."  = "TDN (mg N/L)",
  "NH4..ug.N.L."  = "NH4 (µg N/L)", "PO4..ug.P.L."  = "PO4 (µg P/L)",
  "Cl..mg.Cl.L."  = "Cl (mg Cl/L)", "NO3..mg.N.L."  = "NO3 (mg N/L)",
  "SO4..mg.S.L."  = "SO4 (mg S/L)", "Na..mg.Na.L."  = "Na (mg Na/L)",
  "K..mg.K.L."    = "K (mg K/L)",   "Mg..mg.Mg.L."  = "Mg (mg Mg/L)",
  "Ca..mg.Ca.L."  = "Ca (mg Ca/L)"
)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

# Load chem_avg_flagged from Google Drive
chem_folder <- drive_get("https://drive.google.com/drive/folders/1yNoGf43hkNDC5TX7wMWc27PYvcRD0vul")
chem_csv    <- googledrive::drive_ls(path = chem_folder, type = "csv")

googledrive::drive_download(
  file      = chem_csv$id[chem_csv$name == "chem_avg_flagged.csv"],
  path      = "googledrive/chem_avg_flagged.csv",
  overwrite = TRUE
)

wqual <- read.csv("googledrive/chem_avg_flagged.csv") %>%
  rename(Watershed = Sub_Project) %>%
  mutate(
    Date        = as.Date(Date),
    Watershed   = str_trim(Watershed),
    Site        = str_trim(Site),
    Sample_Type = if_else(is.na(Sample_Type), "stream", str_trim(Sample_Type)),
    across(all_of(solutes), as.numeric)
  )

message(sprintf("✔ Loaded chem_avg_flagged: %d rows, %d watersheds",
                nrow(wqual), n_distinct(wqual$Watershed)))
message("Sample types:")
print(table(wqual$Sample_Type, wqual$Watershed))

# Load missing data tracking sheet
googledrive::drive_download(
  file      = chem_csv$id[chem_csv$name == "Missing QuEST data NM 260514.csv"],
  path      = "googledrive/Missing QuEST data NM 260514.csv",
  overwrite = TRUE
)

missing <- read.csv("googledrive/Missing QuEST data NM 260514.csv") %>%
  rename(Watershed = Sub_Project) %>%
  mutate(
    Date      = as.Date(Collection.Date, format = "%m/%d/%y"),
    Watershed = str_trim(Watershed),
    Site      = str_trim(Site)
  )

# ── 2. SPLIT STREAM VS PRECIPITATION SAMPLES ─────────────────────────────────
# Stream samples are the main dataset for leverage and most analyses.
# Snow/rain are kept separate for precipitation chemistry comparisons.

wqual_stream <- wqual %>% filter(Sample_Type == "stream")
wqual_precip <- wqual %>% filter(Sample_Type %in% c("snow", "rain"))
wqual_blank  <- wqual %>% filter(Sample_Type == "blank")  # kept for reference

message(sprintf("✔ Stream samples: %d | Precipitation: %d | Blanks: %d",
                nrow(wqual_stream), nrow(wqual_precip), nrow(wqual_blank)))

# ── 3. SUMMARY STATISTICS ─────────────────────────────────────────────────────

compact_stats <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      across(
        all_of(solutes),
        list(
          N      = ~ sum(!is.na(.)),
          N_miss = ~ sum(is.na(.)),
          Mean   = ~ round(mean(., na.rm = TRUE), 3),
          SD     = ~ round(sd(., na.rm = TRUE), 3),
          Min    = ~ round(min(., na.rm = TRUE), 3),
          Median = ~ round(median(., na.rm = TRUE), 3),
          Max    = ~ round(max(., na.rm = TRUE), 3)
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(contains("__"), names_to = c("Solute", "Stat"), names_sep = "__") %>%
    mutate(Solute = recode(Solute, !!!solute_labels)) %>%
    pivot_wider(names_from = Stat, values_from = value)
}

# Stream samples — main groupings
write_csv(compact_stats(wqual_stream, "Watershed"),
          "data/summary_stats_by_watershed.csv")
write_csv(compact_stats(wqual_stream, c("Watershed", "Site")),
          "data/summary_stats_by_site.csv")
write_csv(compact_stats(wqual_stream, c("Watershed", "Date")),
          "data/summary_stats_by_date.csv")
write_csv(compact_stats(wqual_stream, c("Watershed", "Site", "event")),
          "data/summary_stats_by_event.csv")

# All sample types together — useful for comparing stream vs snow/rain
write_csv(compact_stats(wqual, c("Watershed", "Sample_Type")),
          "data/summary_stats_by_sample_type.csv")

# Precipitation samples separately
if (nrow(wqual_precip) > 0) {
  write_csv(compact_stats(wqual_precip, c("Watershed", "Sample_Type", "Site")),
            "data/summary_stats_precipitation.csv")
}

message("✔ Summary stats CSVs saved.")

# ── 4. MISSING DATA SUMMARY ───────────────────────────────────────────────────

# % missing per solute × watershed (stream samples only)
miss_pct <- wqual_stream %>%
  group_by(Watershed) %>%
  summarise(
    n_total   = n(),
    across(all_of(solutes), ~ mean(is.na(.)) * 100),
    .groups = "drop"
  ) %>%
  pivot_longer(all_of(solutes), names_to = "Solute", values_to = "pct_missing") %>%
  mutate(Solute = recode(Solute, !!!solute_labels))


# ── 5. REP VARIABILITY SUMMARY ────────────────────────────────────────────────
# rep_variability data comes from Script 2 output.
# Summarise from the flags already in chem_avg_flagged.

rep_summary <- wqual %>%
  filter(n_reps > 1) %>%
  group_by(Watershed, Sample_Type) %>%
  summarise(
    n_multi_rep    = n(),
    n_high_rpd     = sum(high_rpd_any, na.rm = TRUE),
    pct_high_rpd   = round(n_high_rpd / n_multi_rep * 100, 1),
    .groups = "drop"
  )

write_csv(rep_summary, "data/rep_variability_summary.csv")
message("✔ Rep variability summary saved.")

# ── 6. OUTLIER RECORDS ────────────────────────────────────────────────────────
# Pull from the flag columns already in chem_avg_flagged.

out_either_cols <- paste0(solutes, "_out_either")

outlier_records <- wqual %>%
  filter(n_outlier_flags > 0) %>%
  pivot_longer(all_of(solutes), names_to = "Solute", values_to = "Value") %>%
  rowwise() %>%
  mutate(
    out_IQR    = get(paste0(Solute, "_out_IQR")),
    out_zscore = get(paste0(Solute, "_out_zscore")),
    out_either = get(paste0(Solute, "_out_either")),
    Solute     = solute_labels[Solute]
  ) %>%
  ungroup() %>%
  filter(out_either == TRUE, !is.na(Value)) %>%
  select(Watershed, Site, Date, Sample_Type, event, Solute, Value,
         out_IQR, out_zscore, flag_summary) %>%
  arrange(Watershed, desc(Value))

write_csv(outlier_records, "data/outlier_records.csv")
message(sprintf("✔ Outlier records saved: %d records", nrow(outlier_records)))

# ── 7. MISSING DATA HEATMAPS ──────────────────────────────────────────────────

# 7a. By watershed (stream only)
p_heatmap <- ggplot(miss_pct, aes(x = Solute, y = Watershed, fill = pct_missing)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pct_missing, 1)), size = 3) +
  scale_fill_gradient(low = "white", high = "#c0392b", name = "% Missing") +
  labs(title = "Missing Data — % Missing by Watershed & Solute (stream samples)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/01_missing_heatmap.png", p_heatmap, width = 13, height = 5, dpi = 150)

# 7b. By site
miss_site <- wqual_stream %>%
  group_by(Watershed, Site) %>%
  summarise(
    across(all_of(solutes), ~ mean(is.na(.)) * 100),
    .groups = "drop"
  ) %>%
  pivot_longer(all_of(solutes), names_to = "Solute", values_to = "pct_missing") %>%
  mutate(Solute = recode(Solute, !!!solute_labels))

p_heatmap_site <- ggplot(miss_site, aes(x = Solute, y = Site, fill = pct_missing)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct_missing > 0, round(pct_missing, 0), "")), size = 2.2) +
  scale_fill_gradient(low = "white", high = "#c0392b", name = "% Missing") +
  facet_wrap(~ Watershed, scales = "free_y", ncol = 2) +
  labs(title = "Missing Data by Site & Solute", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/02_missing_heatmap_site.png", p_heatmap_site,
       width = 15, height = 14, dpi = 150)

# ── 8. REP VARIABILITY PLOTS ─────────────────────────────────────────────────
# Summarised from flags already in chem_avg_flagged — no raw rep data needed.

# Rep variability bar chart by watershed and sample type
p_rpd_bar <- rep_summary %>%
  ggplot(aes(x = Watershed, y = pct_high_rpd, fill = Sample_Type)) +
  geom_col(position = "dodge", color = "white") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "#c0392b", linewidth = 0.5) +
  scale_fill_viridis_d(name = "Sample type") +
  labs(title = paste0("Rep variability — % of multi-rep rows with RPD > ",
                      RPD_THRESHOLD, "% (above 10× MDL)"),
       x = NULL, y = "% high RPD") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/03_rep_rpd_bar.png", p_rpd_bar, width = 10, height = 5, dpi = 150)
message("✔ Rep variability plot saved.")

# ── 9. BOXPLOTS BY WATERSHED ──────────────────────────────────────────────────
# Stream samples only — main comparison

plots_ws <- map(solutes, function(sol) {
  ggplot(wqual_stream, aes(x = Watershed, y = .data[[sol]], fill = Watershed)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, na.rm = TRUE, alpha = 0.8) +
    scale_fill_viridis_d(guide = "none") +
    labs(title = solute_labels[sol], x = NULL, y = solute_labels[sol]) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
})

ggsave("figures/04_boxplots_by_watershed.png",
       wrap_plots(plots_ws, ncol = 3), width = 16, height = 18, dpi = 150)

# ── 10. BOXPLOTS BY SITE ──────────────────────────────────────────────────────

walk2(solutes, seq_along(solutes), function(sol, i) {
  p <- ggplot(wqual_stream, aes(x = Site, y = .data[[sol]], fill = Site)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1.5, na.rm = TRUE, alpha = 0.8) +
    scale_fill_viridis_d(guide = "none") +
    facet_wrap(~ Watershed, scales = "free", ncol = 2) +
    labs(title = solute_labels[sol], x = NULL, y = solute_labels[sol]) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(sprintf("figures/05_%02d_boxplot_site_%s.png", i,
                 gsub("[^A-Za-z0-9]", "_", sol)), p, width = 13, height = 11, dpi = 150)
})
message("✔ Boxplots saved.")

# ── 11. TIME SERIES BY SITE ───────────────────────────────────────────────────

walk2(solutes, seq_along(solutes), function(sol, i) {
  p <- wqual_stream %>%
    filter(!is.na(.data[[sol]]), !is.na(Date)) %>%
    ggplot(aes(x = Date, y = .data[[sol]], color = Site, group = Site)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = factor(n_reps > 1)), size = 1.8, alpha = 0.9) +
    scale_color_viridis_d() +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                       labels = c("Single sample", "Averaged reps"), name = NULL) +
    scale_x_date(date_labels = "%b %Y") +
    facet_wrap(~ Watershed, scales = "free", ncol = 2) +
    labs(title = solute_labels[sol], x = "Date", y = solute_labels[sol], color = "Site") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right")
  ggsave(sprintf("figures/06_%02d_timeseries_%s.png", i,
                 gsub("[^A-Za-z0-9]", "_", sol)), p, width = 14, height = 10, dpi = 150)
})
message("✔ Time series plots saved.")

# ── 12. PRECIPITATION CHEMISTRY PLOTS ────────────────────────────────────────
# Only produced if there are snow or rain samples

if (nrow(wqual_precip) > 0) {

  # Boxplot comparing snow vs rain vs stream per solute
  p_precip_compare <- map(solutes, function(sol) {
    bind_rows(
      wqual_stream %>% mutate(type_label = "Stream"),
      wqual_precip %>% mutate(type_label = str_to_title(Sample_Type))
    ) %>%
      filter(!is.na(.data[[sol]])) %>%
      ggplot(aes(x = type_label, y = .data[[sol]], fill = type_label)) +
      geom_boxplot(alpha = 0.8, outlier.shape = 21) +
      scale_fill_viridis_d(guide = "none") +
      facet_wrap(~ Watershed, scales = "free_y", ncol = 2) +
      labs(title = solute_labels[sol], x = NULL, y = solute_labels[sol]) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  ggsave("figures/10_precip_vs_stream_boxplots.png",
         wrap_plots(p_precip_compare, ncol = 3),
         width = 16, height = 18, dpi = 150)
  message("✔ Precipitation vs stream comparison plots saved.")
}

# ── 13. OUTLIER STRIP CHARTS ──────────────────────────────────────────────────

outlier_df_long <- wqual_stream %>%
  pivot_longer(all_of(solutes), names_to = "Solute", values_to = "Value") %>%
  rowwise() %>%
  mutate(
    out_IQR    = get(paste0(Solute, "_out_IQR")),
    out_zscore = get(paste0(Solute, "_out_zscore")),
    out_either = get(paste0(Solute, "_out_either"))
  ) %>%
  ungroup() %>%
  filter(!is.na(Value)) %>%
  mutate(Solute = recode(Solute, !!!solute_labels))

walk2(unique(outlier_df_long$Solute), seq_along(unique(outlier_df_long$Solute)),
      function(sol_label, i) {
        df <- outlier_df_long %>%
          filter(Solute == sol_label) %>%
          mutate(flag = case_when(
            out_IQR & out_zscore ~ "Both",
            out_IQR              ~ "IQR only",
            out_zscore           ~ "Z-score only",
            TRUE                 ~ "Normal"
          ))
        p <- ggplot(df, aes(x = Watershed, y = Value,
                            color = flag, shape = flag, size = flag)) +
          geom_jitter(width = 0.25, alpha = 0.75) +
          scale_color_manual(values = c("Both"="#c0392b","IQR only"="#e67e22",
                                        "Z-score only"="#8e44ad","Normal"="grey60")) +
          scale_shape_manual(values = c("Both"=17,"IQR only"=15,
                                        "Z-score only"=18,"Normal"=16)) +
          scale_size_manual(values  = c("Both"=3,"IQR only"=2.5,
                                        "Z-score only"=2.5,"Normal"=1.5)) +
          labs(title = paste("Outliers –", sol_label),
               x = NULL, y = sol_label,
               color = "Type", shape = "Type", size = "Type") +
          theme_minimal(base_size = 11) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1))
        ggsave(sprintf("figures/07_%02d_outliers_%s.png", i,
                       gsub("[^A-Za-z0-9 ]", "_", sol_label)),
               p, width = 10, height = 6, dpi = 150)
      })
message("✔ Outlier plots saved.")

# ── 14. MISSING FLAGS BAR CHART & SAMPLE COUNTS ──────────────────────────────
# 
# p_flags <- missing %>%
#   filter(flags != "no flag") %>%
#   count(Watershed, flags) %>%
#   ggplot(aes(x = Watershed, y = n, fill = flags)) +
#   geom_col(position = "dodge", color = "white") +
#   scale_fill_viridis_d(name = "Flag") +
#   labs(title = "Samples Requiring Resending — by Watershed & Flag Type",
#        x = NULL, y = "Number of Samples") +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1))
# 
# ggsave("figures/08_missing_flags.png", p_flags, width = 10, height = 5, dpi = 150)

# Sample counts — stream sites only
p_obs <- wqual_stream %>%
  count(Watershed, Site) %>%
  ggplot(aes(x = Site, y = n, fill = Watershed)) +
  geom_col(color = "white") +
  facet_wrap(~ Watershed, scales = "free_x", ncol = 2) +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "Number of Averaged Samples by Site & Watershed (stream)",
       x = NULL, y = "n") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/09_sample_counts.png", p_obs, width = 14, height = 8, dpi = 150)

# Sample counts including precipitation
if (nrow(wqual_precip) > 0) {
  p_obs_all <- wqual %>%
    filter(Sample_Type != "blank") %>%
    count(Watershed, Sample_Type, Site) %>%
    ggplot(aes(x = Site, y = n, fill = Sample_Type)) +
    geom_col(color = "white") +
    facet_wrap(~ Watershed, scales = "free_x", ncol = 2) +
    scale_fill_viridis_d(name = "Sample type") +
    labs(title = "Number of Averaged Samples by Site & Watershed (all types)",
         x = NULL, y = "n") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("figures/09b_sample_counts_all_types.png", p_obs_all,
         width = 14, height = 8, dpi = 150)
}

