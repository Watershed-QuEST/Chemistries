# =============================================================================
# Script 2: Rep variability flags, rep averaging, outlier flags, and QC flags
# =============================================================================
# Inputs:
#   chem_data   = output of Script 1 (MDL-adjusted, with _flag columns)
#   valid_sites = site lookup table with Sub_Project and Site columns
#
# This script works on chem_data BEFORE any watershed-specific merging.
# The idea is to flag and average the chemistry data once, cleanly, so that
# each watershed script can simply join the averaged table to its discharge data.
#
# Outputs:
#   chem_flagged  = pre-average table with all flags attached (one row per raw sample)
#   chem_avg      = averaged table (one row per site × date) with all flags
#   outlier_table = long-format table of all outlier records for review
# =============================================================================

library(tidyverse)
library(lubridate)
library(googledrive)

# ── 0. CONFIGURATION ──────────────────────────────────────────────────────────
IQR_MULT      <- 1.5   # IQR multiplier for outlier detection
Z_THRESH      <- 3     # z-score threshold for outlier detection

# Solute columns (must match names in chem_data after Script 1)
solutes <- c(
  "NPOC..mg.C.L.", "TDN..mg.N.L.",  "NH4..ug.N.L.",
  "PO4..ug.P.L.",  "Cl..mg.Cl.L.",  "NO3..mg.N.L.",
  "SO4..mg.S.L.",  "Na..mg.Na.L.",  "K..mg.K.L.",
  "Mg..mg.Mg.L.",  "Ca..mg.Ca.L."
)

# Friendly labels (used in output tables)
solute_labels <- c(
  "NPOC..mg.C.L." = "NPOC (mg C/L)", "TDN..mg.N.L."  = "TDN (mg N/L)",
  "NH4..ug.N.L."  = "NH4 (µg N/L)", "PO4..ug.P.L."  = "PO4 (µg P/L)",
  "Cl..mg.Cl.L."  = "Cl (mg Cl/L)", "NO3..mg.N.L."  = "NO3 (mg N/L)",
  "SO4..mg.S.L."  = "SO4 (mg S/L)", "Na..mg.Na.L."  = "Na (mg Na/L)",
  "K..mg.K.L."    = "K (mg K/L)",   "Mg..mg.Mg.L."  = "Mg (mg Mg/L)",
  "Ca..mg.Ca.L."  = "Ca (mg Ca/L)"
)

# ── 1. CLEAN & PREPARE ────────────────────────────────────────────────────────
# chem_data comes from Script 1 with:
#   - Site column (renamed from Sub_ProjectB)
#   - Collection.Date already parsed
#   - MDL replacement already applied
#   - {solute}_flag columns (TRUE = was below MDL, replaced with MDL/2)

#### load chem data ####
# chem data is for all the sites
chem <- drive_get("https://drive.google.com/drive/folders/1ZCVAoIamyMMtwh-Cy3SpeQx2IWYu6gg2")
3
# list all CSV files in the folder
chem_csv <- googledrive::drive_ls(path = chem, type = "csv")

# call the specific file you want (most recent one)
googledrive::drive_download(file = chem_csv$id[chem_csv$name=="2026-05-15_chem_data.csv"], 
                            path = "googledrive/2026-05-15_chem_data.csv",
                            overwrite = T)
# load it into R
chem_data = read.csv("googledrive/2026-05-15_chem_data.csv")

# Rename Collection.Date → Date to match your watershed scripts
chem_data <- chem_data %>%
  rename(Date = Collection.Date)

# Coerce all solute columns to numeric (safety check)
chem_data <- chem_data %>%
  mutate(across(all_of(solutes), as.numeric))

# Specify Site for NH
chem_data <- chem_data %>%
  mutate(
    Site = if_else(str_trim(Sub_Project) == "QuEST",
                   str_trim(Sample.Name),
                   str_trim(Site))
  )

# Filter to valid sites (uses the same valid_sites lookup from your analysis)
valid_sites <- tribble(
  ~Sub_Project, ~Site,
  # New Mexico — stream sites
  "New Mexico","USF01","New Mexico","USF02","New Mexico","USF03","New Mexico","USF04",
  "New Mexico","USF05","New Mexico","USF06","New Mexico","USF07","New Mexico","USF08",
  "New Mexico","USF09","New Mexico","USF10","New Mexico","USF11","New Mexico","USF12",
  "New Mexico","USF13","New Mexico","USF14","New Mexico","USF15","New Mexico","USF16",
  "New Mexico","USF17","New Mexico","USF18","New Mexico","USF19","New Mexico","USF20",
  "New Mexico","USF21","New Mexico","USF22","New Mexico","USF23","New Mexico","USF24",
  "New Mexico","USF25","New Mexico","USF26","New Mexico","USF27","New Mexico","USF28",
  "New Mexico","USF29","New Mexico","USF30","New Mexico","USF31","New Mexico","USF32",
  "New Mexico","USF33","New Mexico","USF40","New Mexico","USF41",
  # New Mexico — snow and rain (add more as they appear)
  "New Mexico","SNOW1","New Mexico","SNOW2","New Mexico","SNOW3","New Mexico","SNOWMELT",
  "New Mexico","RAIN","New Mexico","RAIN1","New Mexico","RAIN2","New Mexico","RAIN3",
  # Arkansas — stream sites
  "Arkansas","BRM01","Arkansas","BRMQ1","Arkansas","BRAA1","Arkansas","BRA01",
  "Arkansas","BRM02","Arkansas","BRAB1","Arkansas","BRB01","Arkansas","BRM03",
  "Arkansas","BRMQ3","Arkansas","BRM04","Arkansas","BRCD1","Arkansas","BRMQ4",
  "Arkansas","BRD01","Arkansas","BRE01","Arkansas","BRM05","Arkansas","BRF01",
  "Arkansas","BRM06","Arkansas","BRM07","Arkansas","BRC01","Arkansas","BRA02",
  # Nevada — stream sites
  "Nevada","DVO",   "Nevada","DVSB1","Nevada","DVSB2","Nevada","DVMS1","Nevada","DVMS2",
  "Nevada","DVMS3", "Nevada","DVMS4","Nevada","DVMS5","Nevada","DVMS6","Nevada","DVET",
  "Nevada","DVNWT1","Nevada","DVNWT2","Nevada","DVNWT3","Nevada","DVNWT4","Nevada","DVNWT5",
  "Nevada","DVWT1", "Nevada","DVWT2","Nevada","DVWT3","Nevada","DVWT4","Nevada","DVWT5",
  # New Hampshire (Lamprey — site from Sample.Name)
  "New Hampshire","SMB",  "New Hampshire","OMC",  "New Hampshire","PRC",   "New Hampshire","LMP00","New Hampshire","DCR",
  "New Hampshire","CTB",  "New Hampshire","LST01","New Hampshire","LMP01", "New Hampshire","DDB",  "New Hampshire","NCB",
  "New Hampshire","NCB-DOWN","New Hampshire","LMP07","New Hampshire","LMP09","New Hampshire","HRB","New Hampshire","NBR-UP",
  "New Hampshire","NBR",  "New Hampshire","LMP12","New Hampshire","LMP19", "New Hampshire","LMP27","New Hampshire","LMP72",
  # Alabama — stream sites
  "Alabama","SSM01","Alabama","SST02","Alabama","SST03","Alabama","SST04","Alabama","SST05",
  "Alabama","SST06","Alabama","SST07","Alabama","MAYFNEON","Alabama","SST08","Alabama","SST09",
  "Alabama","SSM10","Alabama","SST11","Alabama","SST12","Alabama","SST13","Alabama","SST14",
  "Alabama","SST15","Alabama","SST16","Alabama","SST17","Alabama","SST18","Alabama","SST19",
  "Alabama","SSM20"
)
# ── Fill empty Site from Sample.Name ──────────────────────────────────────────
# Only applies to rows where Site is missing or empty.
# Extraction logic differs by watershed based on Sample.Name patterns.

chem_data <- chem_data %>%
  mutate(
    # ── Step 1: classify sample type ─────────────────────────────────────────
    Sample_Type = case_when(
      str_detect(Sample.Name, regex("snow|snowmelt", ignore_case = TRUE)) ~ "snow",
      str_detect(Sample.Name, regex("rain",          ignore_case = TRUE)) ~ "rain",
      str_detect(Sample.Name, regex("blank",         ignore_case = TRUE)) ~ "blank",
      TRUE ~ "stream"
    ),
    
    # ── Step 2: extract Site ──────────────────────────────────────────────────
    Site = case_when(
      # Already has a usable site and it's a stream — leave it alone
      !is.na(Site) & str_trim(Site) != "" & Sample_Type == "stream" ~ Site,
      
      # Blanks: Label them clearly by watershed so they can pass through
      Sample_Type == "blank" ~ 
        paste0(str_replace_all(str_trim(Sub_Project), " ", ""), "_BLANK"),
      
      # New Mexico Rain/Snow: Extract the exact site identifier (e.g., SNOW1, SNOWMELT, RAIN)
      str_trim(Sub_Project) == "New Mexico" & Sample_Type %in% c("snow", "rain") ~ 
        str_extract(Sample.Name, regex("SNOW\\d*|SNOWMELT|RAIN\\d*", ignore_case = TRUE)) %>% 
        str_to_upper(),
      
      # Nevada: date-prefixed → 4th segment; no date → whole string
      str_trim(Sub_Project) == "Nevada" ~
        if_else(
          str_detect(Sample.Name, "^\\d{4}[_-]\\d{2}[_-]\\d{2}"),
          str_split_i(Sample.Name, "_", 4),
          str_trim(Sample.Name)
        ),
      
      # Arkansas
      str_trim(Sub_Project) == "Arkansas" ~
        str_split_i(Sample.Name, "_", 2),
      
      # Alabama
      str_trim(Sub_Project) == "Alabama" ~
        str_split_i(Sample.Name, "_", 2),
      
      TRUE ~ Site
    ),
    Site = str_trim(Site)
  )

# Check sample types and unified site names
table(chem_data$Sample_Type, chem_data$Sub_Project)
chem_data %>% filter(Sample_Type != "stream") %>% count(Sub_Project, Sample_Type, Site)

# ── Normalize USF site codes ───────────────────────────────────────────────────
# Some sites were entered as USF4 instead of USF04 — standardize to USF + 2-digit number
chem_data <- chem_data %>%
  mutate(
    Site = if_else(
      str_trim(Sub_Project) == "New Mexico" & str_detect(Site, "^USF\\d$"),
      str_replace(Site, "^USF(\\d)$", "USF0\\1"),
      Site
    )
  )

# ── Fix known data entry errors ────────────────────────────────────────────────

# Fix 1: SSM sites that should be SST
# SSM01, SSM10, SSM20 are correct mainstem sites — leave those alone.
# All other SSM entries are typos and should be SST.
chem_data <- chem_data %>%
  mutate(
    Site = if_else(
      str_trim(Sub_Project) == "Alabama" &
        str_detect(Site, "^SSM") &
        !Site %in% c("SSM01", "SSM10", "SSM20"),
      str_replace(Site, "^SSM", "SST"),
      Site
    )
  )

# Fix 2: SSM01 sample that was incorrectly assigned to Arkansas
# Sample 2024-09-16_SSM01_R1 belongs to Alabama
chem_data <- chem_data %>%
  mutate(
    Sub_Project = if_else(
      str_detect(Sample.Name, "2024-09-16_SSM01") &
        str_trim(Sub_Project) == "Arkansas",
      "Alabama",
      Sub_Project
    )
  )

# Fix 3: BRM sites entered with sequential numbers (BRM01–BRM30) instead of
# the correct site code embedded in Sample.Name.
# Extract the correct site from Sample.Name for all affected rows.
# Pattern: "2025-02-17_SITE_Rep N" → extract SITE (2nd segment)
chem_data <- chem_data %>%
  mutate(
    Site = if_else(
      str_trim(Sub_Project) == "Arkansas" &
        str_detect(Sample.Name, "^2025-02-17_BR") &
        str_detect(Site, "^BRM\\d{2}$"),  # catches BRM01 through BRM30
      str_split_i(Sample.Name, "_", 2),    # extract 2nd segment from Sample.Name
      Site
    )
  )

# Fix 4: Nevada site name typos
chem_data <- chem_data %>%
  mutate(
    Site = case_when(
      # DVMT3 → DVWT3 (typo in site name)
      str_trim(Sub_Project) == "Nevada" & Site == "DVMT3" ~ "DVWT3",
      
      # DVSB on 2024-09-26 → DVSB2 (missing rep number)
      str_trim(Sub_Project) == "Nevada" & Site == "DVSB" &
        as.Date(Date) == as.Date("2024-09-26") ~ "DVSB2",
      
      # DVWT4 → DVNWT4 (wrong prefix)
      str_trim(Sub_Project) == "Nevada" & Site == "DVWT4" ~ "DVNWT4",
      
      TRUE ~ Site
    )
  )

# Spot check
message("=== Nevada site fixes check ===")
chem_data %>%
  filter(str_trim(Sub_Project) == "Nevada") %>%
  count(Site) %>%
  filter(str_detect(Site, "^DV")) %>%
  arrange(Site) %>%
  print()

# Fix 5: New Hampshire/Lamprey site name typos
chem_data <- chem_data %>%
  mutate(
    # Fix watershed name: QuEST → New Hampshire
    Sub_Project = if_else(
      str_trim(Sub_Project) == "QuEST",
      "New Hampshire",
      Sub_Project
    ),
    
    # Fix site name typos (QuEST rows — now named Lamprey River Watershed)
    Site = case_when(
      str_trim(Sub_Project) == "New Hampshire" &
        str_to_upper(str_trim(Site)) == "SBM" ~ "SMB",   # normalize any casing variants
      
      str_trim(Sub_Project) == "New Hampshire" &
        Site == "LMP0" ~ "LMP00",
      
      TRUE ~ Site
    )
  )

# Spot check
message("=== Lamprey site fixes check ===")
chem_data %>%
  filter(str_trim(Sub_Project) == "New Hampshire") %>%
  count(Site) %>%
  arrange(Site) %>%
  print()

#-------------------------------------------------------------------------------

# ── SITE SUMMARY (before filtering) ──────────────────────────────────────────
# Pull all unique site names per watershed before dropping non-valid sites.
# Useful for auditing what's in the raw data vs what gets kept.

site_summary_raw <- chem_data %>%
  mutate(Site = str_trim(Site)) %>%
  filter(!is.na(Site) & Site != "") %>%
  group_by(Sub_Project, Sample_Type, Site) %>%
  summarise(
    n_samples = n(),
    date_min  = min(as.Date(Date), na.rm = TRUE),
    date_max  = max(as.Date(Date), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Sub_Project, Sample_Type, Site)

write.csv(site_summary_raw, "data/site_summary_raw.csv", row.names = FALSE)

# Run this before the semi_join line
nrow(chem_data)

# Filter to valid stream/precipitation sites OR keep any blank samples
chem_data <- chem_data %>%
  filter(
    Sample_Type == "blank" | 
      paste0(Sub_Project, "-", Site) %in% paste0(valid_sites$Sub_Project, "-", valid_sites$Site)
  )

# chem_data <- chem_data %>%
#   semi_join(valid_sites, by = c("Sub_Project" = "Sub_Project", "Site" = "Site"))

# Then after
nrow(chem_data)

message(sprintf("✔ Rows after site filtering: %d", nrow(chem_data)))

# ── 2. MDL FLAG SUMMARY ───────────────────────────────────────────────────────
# Script 1 already added {solute}_flag columns (TRUE = below MDL).
# Here we add a convenience column: total number of MDL-flagged solutes per row.

flag_cols <- paste0(solutes, "_flag")
existing_flag_cols <- flag_cols[flag_cols %in% names(chem_data)]

chem_data <- chem_data %>%
  mutate(
    n_below_MDL = rowSums(across(all_of(existing_flag_cols), ~ . == TRUE), na.rm = TRUE)
  )

message(sprintf("✔ MDL flag columns found: %d of %d", length(existing_flag_cols), length(flag_cols)))

"n_below_MDL" %in% names(chem_data)

# ── 3. REP VARIABILITY FLAG (before averaging) ────────────────────────────────
# Uses RPD (relative percent difference) per Jody's lab QC protocol.
# RPD = |A - B| / mean(A, B) × 100 — identical to CV × √2 for 2 reps.
# For >2 reps: mean of all pairwise RPDs.
#
# Jody's rule: only flag if RPD > RPD_THRESHOLD AND both values > 10× MDL.
# Below 10× MDL, high variability near detection is expected and not flagged.

RPD_THRESHOLD <- 15    # % RPD threshold (Jody uses 15%)

# MDL values from Jody's 2025 table
mdl_values <- c(
  "NPOC..mg.C.L." = 0.14,   # mg C/L
  "TDN..mg.N.L."  = 0.03,   # mg N/L
  "NH4..ug.N.L."  = 6,      # ug N/L
  "PO4..ug.P.L."  = 2,      # ug P/L
  "Cl..mg.Cl.L."  = 0.05,   # mg Cl/L
  "NO3..mg.N.L."  = 0.005,  # mg N/L (IC method — all watersheds except Arkansas)
  "SO4..mg.S.L."  = 0.05,   # mg S/L
  "Na..mg.Na.L."  = 0.1,    # mg Na/L
  "K..mg.K.L."    = 0.02,   # mg K/L
  "Mg..mg.Mg.L."  = 0.1,    # mg Mg/L
  "Ca..mg.Ca.L."  = 0.2     # mg Ca/L
)

# 10× MDL thresholds — variability only flagged above these values
mdl_10x <- mdl_values * 10

# Helper: mean pairwise RPD for a vector of values
mean_pairwise_rpd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  pairs <- combn(x, 2)
  rpds  <- abs(pairs[1,] - pairs[2,]) / rowMeans(t(pairs)) * 100
  round(mean(rpds, na.rm = TRUE), 1)
}

rep_rpd <- chem_data %>%
  group_by(Sub_Project, Site, Date) %>%
  mutate(n_reps = n()) %>%
  filter(n_reps > 1) %>%
  summarise(
    n_reps = first(n_reps),
    across(
      all_of(solutes),
      ~ mean_pairwise_rpd(.),
      .names = "rpd_{.col}"
    ),
    # also store mean value per solute to apply 10× MDL rule
    across(
      all_of(solutes),
      ~ mean(., na.rm = TRUE),
      .names = "mean_{.col}"
    ),
    .groups = "drop"
  )

# Long format: one row per site × date × solute
rep_rpd_long <- rep_rpd %>%
  pivot_longer(starts_with("rpd_"), names_to = "Solute", values_to = "RPD") %>%
  mutate(Solute = str_remove(Solute, "^rpd_")) %>%
  left_join(
    rep_rpd %>%
      pivot_longer(starts_with("mean_"), names_to = "Solute", values_to = "mean_val") %>%
      mutate(Solute = str_remove(Solute, "^mean_")) %>%
      select(Sub_Project, Site, Date, Solute, mean_val),
    by = c("Sub_Project", "Site", "Date", "Solute")
  ) %>%
  mutate(
    mdl_10x_val = mdl_10x[Solute],
    above_10x_MDL = !is.na(mean_val) & mean_val > mdl_10x_val,
    # only flag if RPD > threshold AND mean value is above 10× MDL
    high_RPD = !is.na(RPD) & RPD > RPD_THRESHOLD & above_10x_MDL
  )

# Flagged records
rep_rpd_flagged <- rep_rpd_long %>%
  filter(high_RPD) %>%
  arrange(Sub_Project, desc(RPD))

write.csv(rep_rpd_flagged, "data/flags_rep_variability.csv",  row.names = FALSE)
write.csv(rep_rpd,         "data/flags_rep_rpd_wide.csv",     row.names = FALSE)

message(sprintf("✔ Rep variability (RPD): %d site×date×solute combinations flagged (RPD > %d%% AND mean > 10× MDL)",
                nrow(rep_rpd_flagged), RPD_THRESHOLD))

# Join high_rpd_any flag back to raw data
rep_rpd_summary_per_row <- rep_rpd_long %>%
  group_by(Sub_Project, Site, Date) %>%
  summarise(
    n_reps      = first(rep_rpd$n_reps[rep_rpd$Sub_Project == first(Sub_Project) &
                                         rep_rpd$Site == first(Site) &
                                         rep_rpd$Date == first(Date)]),
    high_rpd_any = any(high_RPD, na.rm = TRUE),
    .groups = "drop"
  )

# simpler join approach
rep_rpd_join <- rep_rpd %>%
  select(Sub_Project, Site, Date, n_reps) %>%
  left_join(
    rep_rpd_long %>%
      group_by(Sub_Project, Site, Date) %>%
      summarise(high_rpd_any = any(high_RPD, na.rm = TRUE), .groups = "drop"),
    by = c("Sub_Project", "Site", "Date")
  )

chem_data <- chem_data %>%
  left_join(rep_rpd_join, by = c("Sub_Project", "Site", "Date")) %>%
  mutate(
    n_reps       = replace_na(n_reps, 1),
    high_rpd_any = replace_na(high_rpd_any, FALSE)
  )


# ── 4. AVERAGE REPS (same site × date) ───────────────────────────────────────
# Average solute values and MDL flags across reps.
# For MDL flags: a column is flagged TRUE if ANY rep was below MDL.
# For high_rpd_any: TRUE if ANY rep set for that site×date had high
chem_data$Date <- as.POSIXct(chem_data$Date, format = "%Y-%m-%d")

# Before averaging
chem_data %>% distinct(Sub_Project, Site, Date) %>% nrow()

chem_avg <- chem_data %>%
  group_by(Sub_Project, Site, Date, Sample_Type) %>%
  summarise(
    # Average solute concentrations
    across(
      all_of(solutes),
      ~ if (all(is.na(.))) NA_real_ else mean(., na.rm = TRUE)
    ),
    # MDL flags: TRUE if any rep was below MDL
    across(
      all_of(existing_flag_cols),
      ~ any(. == TRUE, na.rm = TRUE),
      .names = "{.col}"
    ),
    # Rep metadata
    n_reps      = n(),
    high_rpd_any = any(high_rpd_any, na.rm = TRUE),
    n_below_MDL = mean(n_below_MDL, na.rm = TRUE),  # avg MDL flags per rep
    Sample.Name = paste0(first(Site), "_", first(Date), "_Avg"),
    .groups = "drop"
  )

# After averaging — should match the above
nrow(chem_avg)

message(sprintf("✔ Averaged: %d raw rows → %d site×date rows", nrow(chem_data), nrow(chem_avg)))

# ── 5. SYNOPTIC EVENT GROUPING ────────────────────────────────────────────────
# Group dates within 1 day of each other (per watershed + site) into one event.
# This mirrors the cumsum(diff(Date) >= 2) logic in your existing scripts.

chem_avg <- chem_avg %>%
  mutate(Date = as.Date(Date)) %>%          # re-parse after summarise drops the class
  arrange(Sub_Project, Site, Date) %>%
  group_by(Sub_Project, Site) %>%
  mutate(event = cumsum(c(TRUE, diff(Date) >= 2))) %>%
  ungroup()

# Add year and day-of-year (matches your existing scripts)
chem_avg <- chem_avg %>%
  mutate(
    year = year(Date),
    doy  = yday(Date)
  )

# ── 6. OUTLIER FLAGS (on averaged data) ───────────────────────────────────────
# Two methods: IQR (1.5×) and z-score (|z| > 3), computed per watershed × solute.
# Adds {solute}_out_IQR, {solute}_out_zscore, {solute}_out_either columns.

flag_outliers_wide <- function(df) {
  for (sol in solutes) {
    df <- df %>%
      group_by(Sub_Project) %>%
      mutate(
        # IQR method
        .q1  = quantile(.data[[sol]], 0.25, na.rm = TRUE),
        .q3  = quantile(.data[[sol]], 0.75, na.rm = TRUE),
        .iqr = .q3 - .q1,
        !!paste0(sol, "_out_IQR") :=
          !is.na(.data[[sol]]) &
          (.data[[sol]] < (.q1 - IQR_MULT * .iqr) |
           .data[[sol]] > (.q3 + IQR_MULT * .iqr)),
        # Z-score method
        .mn  = mean(.data[[sol]], na.rm = TRUE),
        .sd  = sd(.data[[sol]],   na.rm = TRUE),
        .z   = (.data[[sol]] - .mn) / .sd,
        !!paste0(sol, "_out_zscore") :=
          !is.na(.data[[sol]]) & abs(.z) > Z_THRESH,
        # Either method
        !!paste0(sol, "_out_either") :=
          .data[[paste0(sol, "_out_IQR")]] | .data[[paste0(sol, "_out_zscore")]]
      ) %>%
      ungroup() %>%
      select(-.q1, -.q3, -.iqr, -.mn, -.sd, -.z)
  }
  df
}

chem_avg <- flag_outliers_wide(chem_avg)

# Count total outlier flags per row (for easy filtering)
out_either_cols <- paste0(solutes, "_out_either")
chem_avg <- chem_avg %>%
  mutate(
    n_outlier_flags = rowSums(across(all_of(out_either_cols), ~ . == TRUE), na.rm = TRUE)
  )

message(sprintf("✔ Outlier flags added. Rows with ≥1 outlier flag: %d",
                sum(chem_avg$n_outlier_flags > 0, na.rm = TRUE)))

# ── 7. ADDITIONAL QC FLAGS ────────────────────────────────────────────────────

chem_avg <- chem_avg %>%
  mutate(

    # --- DOC:TDN ratio flag ---
    # Unusually high or low DOC:TDN can signal contamination or analytical issues.
    # Flag if ratio > 50 or < 1 (rough thresholds — adjust for your system).
    doc_tdn_ratio = ifelse(
      !is.na(NPOC..mg.C.L.) & !is.na(TDN..mg.N.L.) & TDN..mg.N.L. > 0,
      NPOC..mg.C.L. / TDN..mg.N.L.,
      NA_real_
    ),
    flag_doc_tdn_ratio = !is.na(doc_tdn_ratio) &
                          (doc_tdn_ratio > 50 | doc_tdn_ratio < 1),

    # --- NH4 > TDN flag ---
    # NH4 (converted to mg N/L) cannot exceed TDN. If it does, something is wrong.
    NH4_as_N_mg = NH4..ug.N.L. / 1000,
    flag_NH4_exceeds_TDN = !is.na(NH4_as_N_mg) & !is.na(TDN..mg.N.L.) &
                            NH4_as_N_mg > TDN..mg.N.L.,

    # --- NO3 > TDN flag ---
    flag_NO3_exceeds_TDN = !is.na(NO3..mg.N.L.) & !is.na(TDN..mg.N.L.) &
                            NO3..mg.N.L. > TDN..mg.N.L.,

    # --- Any nutrient internal consistency flag ---
    flag_N_consistency = flag_NH4_exceeds_TDN | flag_NO3_exceeds_TDN,

    # --- Overall QC summary flag ---
    # TRUE if ANY of the additional QC checks fired (not counting MDL or outliers)
    flag_any_QC =  flag_doc_tdn_ratio | flag_N_consistency
  )

# ── 8. MASTER FLAG SUMMARY COLUMN ────────────────────────────────────────────
# Summary of ALL flags fired for each row.
# Makes it easy to filter or sort by "anything flagged".

chem_avg <- chem_avg %>%
  mutate(
    flag_summary = pmap_chr(
      list(
        n_below_MDL     = n_below_MDL,
        high_rpd_any    = high_rpd_any,
        n_outlier_flags = n_outlier_flags,
        flag_doc_tdn_ratio = flag_doc_tdn_ratio,
        flag_N_consistency = flag_N_consistency
      ),
      function(n_below_MDL, high_rpd_any, n_outlier_flags,
               flag_doc_tdn_ratio, flag_N_consistency) {
        flags <- c()
        if (!is.na(n_below_MDL) && n_below_MDL > 0)
          flags <- c(flags, paste0("MDL(", n_below_MDL, ")"))
        if (isTRUE(high_rpd_any))
          flags <- c(flags, "high_RPD")    # ← changed label too
        if (!is.na(n_outlier_flags) && n_outlier_flags > 0)
          flags <- c(flags, paste0("outlier(", n_outlier_flags, ")"))
        if (isTRUE(flag_doc_tdn_ratio))
          flags <- c(flags, "DOC_TDN_ratio")
        if (isTRUE(flag_N_consistency))
          flags <- c(flags, "N_consistency")
        if (length(flags) == 0) "OK" else paste(flags, collapse = "|")
      }
    )
  )

# ── 9. OUTLIER RECORDS TABLE ─────────────────────────────────────────────────
# Long-format table of all rows flagged by both IQR and z-score, for review.

outlier_table <- chem_avg %>%
  select(Sub_Project, Site, Date, event, Sample.Name, all_of(solutes),
         ends_with("_out_IQR"), ends_with("_out_zscore"), ends_with("_out_either"),
         flag_summary) %>%
  pivot_longer(
    cols      = all_of(solutes),
    names_to  = "Solute",
    values_to = "Value"
  ) %>%
  rowwise() %>%
  mutate(
    out_IQR    = get(paste0(Solute, "_out_IQR")),
    out_zscore = get(paste0(Solute, "_out_zscore")),
    out_either = get(paste0(Solute, "_out_either")),
    Solute = solute_labels[Solute]
  ) %>%
  ungroup() %>%
  filter(out_either == TRUE, !is.na(Value)) %>%
  select(Sub_Project, Site, Date, event, Solute, Value,
         out_IQR, out_zscore, flag_summary) %>%
  arrange(Sub_Project, desc(Value))


nrow(outlier_table)
head(outlier_table)

"flag_summary" %in% names(chem_avg)

write.csv(outlier_table, "data/flags_outliers.csv", row.names = FALSE)

message(sprintf("✔ Outlier table: %d records", nrow(outlier_table)))


# ── 10. SAVE OUTPUTS ──────────────────────────────────────────────────────────

# Save averaged + flagged chemistry table (main output — use this in watershed scripts)
write.csv(chem_avg, "data/chem_avg_flagged.csv", row.names = FALSE)

# Save pre-average table with rep flags (useful for traceability)
write.csv(chem_data, "data/chem_raw_flagged.csv", row.names = FALSE)

# define the target folder ID in Google Drive
# this is the current leverage folder
drive_folder_id <- "1yNoGf43hkNDC5TX7wMWc27PYvcRD0vul"

# upload the file to the specified Google Drive folder
drive_put(media = "data/chem_avg_flagged.csv", path = as_id(drive_folder_id))
drive_put(media = "data/chem_raw_flagged.csv", path = as_id(drive_folder_id))

message("
✔ Script 2 complete. Files saved:
  data/chem_avg_flagged.csv     ← main output for watershed scripts
  data/chem_raw_flagged.csv     ← pre-average with rep flags
  data/site_summary_raw.csv     ← all sites before filtering
  data/flags_rep_variability.csv
  data/flags_rep_rpd_wide.csv
  data/flags_outliers.csv

Flag columns in chem_avg_flagged.csv:
  {solute}_flag        TRUE = value was below MDL (replaced with MDL/2 in Script 1)
  {solute}_out_IQR     TRUE = outlier by IQR method (1.5× IQR)
  {solute}_out_zscore  TRUE = outlier by z-score (|z| > 3)
  {solute}_out_either  TRUE = flagged by either method
  high_rpd_any         TRUE = any solute RPD > 15% AND mean > 10x MDL among reps
  n_below_MDL          number of solutes below MDL in this averaged row
  n_outlier_flags      number of solutes flagged as outliers
  n_reps               number of raw reps averaged into this row
  flag_doc_tdn_ratio   TRUE = DOC:TDN ratio < 1 or > 50
  flag_NH4_exceeds_TDN TRUE = NH4 (as mg N/L) > TDN
  flag_NO3_exceeds_TDN TRUE = NO3 > TDN
  flag_N_consistency   TRUE = any N internal consistency issue
  flag_any_QC          TRUE = any of the above QC flags (not MDL or outlier)
  flag_summary         human-readable pipe-separated list of all flags, or 'OK'
  event                synoptic campaign grouping (dates within 2 days → same event)
  year, doy            year and day-of-year for plotting
")

# ── 11. HOW TO USE IN WATERSHED SCRIPTS ──────────────────────────────────────
# In each watershed script, replace your averaging block with:
#
#   chem_avg <- read.csv("data/chem_avg_flagged.csv") %>%
#     filter(Sub_Project == "New Mexico")   # or "Arkansas", "Nevada", etc.
#
# Then join to your discharge data as before:
#   data_nm <- left_join(discharge_data_nm, chem_avg, by = c("Site", "Date"))
#
# The group_by variables differ by watershed (some have Qs, DateTime, flag, etc.)
# but since averaging is already done here, you no longer need those blocks.
# Just join on Site + Date after loading discharge.


