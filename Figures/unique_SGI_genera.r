library(readxl)
library(dplyr)
library(writexl)

# ==== File path ====
file_path <- "/Users/zisan/Library/CloudStorage/OneDrive-UniversityofManitoba/Animal Science/SGI_Publication/Final_clean_version/Supplementary Tables/PiGICo_taxonomy.xlsx"

# ==== Load data ====
sgis <- read_excel(file_path, sheet = "SGI")
pibacs <- read_excel(file_path, sheet = "PiBac")

# ==== Normalize genus column ====
sgis <- sgis %>% mutate(genus = tolower(trimws(genus)))
pibacs <- pibacs %>% mutate(genus = tolower(trimws(genus)))

# ==== Identify SGI-unique genera ====
unique_sgi <- setdiff(unique(sgis$genus), unique(pibacs$genus))

# ==== Extract full taxonomy rows for unique genera ====
sgi_unique_taxonomy <- sgis %>% filter(genus %in% unique_sgi)

# ==== Save as Excel ====
out_path <- "/Users/zisan/Library/CloudStorage/OneDrive-UniversityofManitoba/Animal Science/SGI_Publication/Final_clean_version/Supplementary Tables/SGI_unique_genera_taxonomy.xlsx"
write_xlsx(sgi_unique_taxonomy, out_path)

cat("✅ Saved:", out_path, "\n")
cat("🧬 Number of unique genera:", length(unique_sgi), "\n")
