##===============================================================================================================================
## Project: QuEST
## Script 00: MDL replacement
## Replace chemistry values below detection with 1/2 the MDL value
## MDL values from UNH, 2025 — calculated using USGS long-term approach
## across all 2025 runs (not just QuEST samples), so values are conservative.
##
## IMPORTANT — NO3 split:
##   NO3 + NO2 colorimetric (MDL = 0.010 mg N/L) → Arkansas samples only
##   NO3 by IC              (MDL = 0.005 mg N/L) → all other watersheds
## The column NO3..mg.N.L. in the data represents whichever method was used
## for that watershed, so we apply different MDLs by Sub_Project.
##
## Output: 2026-05-15_chem_data.csv → uploaded to Google Drive edited folder
##         This file is the input for 01_flags_and_averaging.R
##===============================================================================================================================

##################
#### Packages ####
##################
library(googledrive)
library(lubridate)
library(dplyr)

###################################
## Clear folders that we will use ##
###################################
files <- list.files(path = "googledrive", full.names = TRUE)
file.remove(files)

#######################
#### Import & Tidy ####
#######################

#### Load chem data ####
chem <- drive_get("https://drive.google.com/drive/folders/15S3uQuvatW_s61yVBXb0Xv_enkDuJRC8")

# list all CSV files in the folder
chem_csv <- googledrive::drive_ls(path = chem, type = "csv")

googledrive::drive_download(
  file      = chem_csv$id[chem_csv$name == "2026-05-15_chem.csv"],
  path      = "googledrive/2026-05-15_chem.csv",
  overwrite = TRUE
)

wqual <- read.csv("googledrive/2026-05-15_chem.csv")

# rename some columns
wqual <- wqual %>% rename(Site = Sub_ProjectB)

# format date columns
wqual$Collection.Date <- as.Date(wqual$Collection.Date, format = "%d-%b-%y")

####################
#### MDL values ####
####################
# "MDL 2" = second iteration of 2025 MDL calculations using all 2025 runs.
# These are more conservative than MDLs calculated under ideal conditions.
#
# Analyte names must match column names in wqual exactly.
# Note: NO3..mg.N.L. gets two different MDLs depending on watershed — handled below.

mdl_map <- c(
  "PO4..ug.P.L."  = 2,      # Ortho-PO4,  ug P/L
  "NH4..ug.N.L."  = 6,      # NH4,         ug N/L
  "Na..mg.Na.L."  = 0.1,    # Na+,         mg Na/L
  "K..mg.K.L."    = 0.02,   # K+,          mg K/L
  "Mg..mg.Mg.L."  = 0.1,    # Mg2+,        mg Mg/L
  "Ca..mg.Ca.L."  = 0.2,    # Ca2+,        mg Ca/L
  "Cl..mg.Cl.L."  = 0.05,   # Cl-,         mg Cl/L
  "SO4..mg.S.L."  = 0.05,   # SO42-,       mg S/L
  "TDN..mg.N.L."  = 0.03,   # TDN,         mg N/L
  "NPOC..mg.C.L." = 0.14    # DOC/NPOC,    mg C/L
  # NO3 is handled separately below — two MDLs depending on watershed
)

# NO3 MDLs
MDL_NO3_IC          <- 0.005   # NO3 by IC — all watersheds except Arkansas
MDL_NO3_COLORIMETRIC <- 0.010  # NO3 + NO2 colorimetric — Arkansas only

###########################################################
#### Replace values below detection with 1/2 the MDL ####
###########################################################

# Step 1: Apply MDL replacement and flags for all solutes EXCEPT NO3
# (NO3 is handled separately because MDL differs by watershed)

chem_data <- wqual %>%
  mutate(
    # MDL substitution: replace below-MDL values with MDL/2
    across(
      names(mdl_map),
      ~ ifelse(!is.na(.) & . < mdl_map[cur_column()],
               mdl_map[cur_column()] / 2,
               .),
      .names = "{col}"
    ),
    # Flag columns: TRUE = value was below MDL and has been replaced
    across(
      names(mdl_map),
      ~ ifelse(!is.na(.) & . < mdl_map[cur_column()], TRUE, FALSE),
      .names = "{col}_flag"
    )
  )

# Step 2: Apply NO3 MDL replacement — different MDL for Arkansas vs all others
# Arkansas uses NO3+NO2 colorimetric; everyone else uses NO3 by IC
chem_data <- chem_data %>%
  mutate(
    # Choose the correct MDL based on watershed
    .no3_mdl = ifelse(str_trim(Sub_Project) == "Arkansas",
                      MDL_NO3_COLORIMETRIC,
                      MDL_NO3_IC),

    # Replace below-MDL NO3 with MDL/2 using watershed-specific MDL
    NO3..mg.N.L. = ifelse(
      !is.na(NO3..mg.N.L.) & NO3..mg.N.L. < .no3_mdl,
      .no3_mdl / 2,
      NO3..mg.N.L.
    ),

    # Flag: TRUE = was below MDL (using watershed-specific threshold)
    `NO3..mg.N.L._flag` = ifelse(
      !is.na(NO3..mg.N.L.) & NO3..mg.N.L. < .no3_mdl,
      TRUE,
      FALSE
    ),

    # Record which NO3 method was applied — useful for documentation
    NO3_method = ifelse(str_trim(Sub_Project) == "Arkansas",
                        "NO3+NO2_colorimetric",
                        "NO3_IC")
  ) %>%
  select(-.no3_mdl)   # drop the temporary MDL column

# Quick sanity check — count how many values were flagged per analyte
message("=== MDL substitution summary (rows flagged below detection) ===")
flag_cols <- c(paste0(names(mdl_map), "_flag"), "NO3..mg.N.L._flag")
flag_cols_present <- flag_cols[flag_cols %in% names(chem_data)]

chem_data %>%
  summarise(across(all_of(flag_cols_present), ~ sum(. == TRUE, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Analyte", values_to = "n_below_MDL") %>%
  mutate(Analyte = str_remove(Analyte, "_flag$")) %>%
  filter(n_below_MDL > 0) %>%
  arrange(desc(n_below_MDL)) %>%
  print()

# Check NO3 method split
message("=== NO3 method assigned by watershed ===")
print(table(chem_data$Sub_Project, chem_data$NO3_method))

#########################################
#### Save edited chem table to Drive ####
#########################################

write.csv(chem_data, "googledrive/2026-05-15_chem_data.csv", row.names = FALSE)

# Upload (overwrite) to the edited folder in Google Drive
drive_folder_id <- "1ZCVAoIamyMMtwh-Cy3SpeQx2IWYu6gg2"
drive_put(
  media = "googledrive/2026-05-15_chem_data.csv",
  path  = as_id(drive_folder_id)
)
