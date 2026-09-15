#Loading needed libraries

# install.packages(c("tidyverse", "janitor","dplyr","crayon","writexl","ggrepel","reticulate"))
#py_install("chembl_webresource_client", pip = TRUE)
#install.packages(c("httr", "jsonlite", "purrr","furrr))
#install.packages("rlang")

library("rstudioapi")
library(tidyverse)
library(janitor)
library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(readxl)
library(arrow)
library(rlang)
library(ggrepel)
library(reticulate)
library(crayon)



#Set path to current directory
setwd(dirname(getActiveDocumentContext()$path))

#----------------------Variables to change--------------------------------------

multi_gene_handling_method<-"mean"

ratios_to_use<-c("tissue_ratio", "cell_ratio")

#options - "mean","max"

measures_to_plot_in_model<-c("tissue_ratio", "cell_ratio","europepmc")


#-------------------------------------------------------------------------------


#*******************************************************************************
#******************PREPARATION OF DATABASE**************************************

# ---- 1) Load & clean ----
#************ My data***********************************************************


# -----1.1 Load the dataframe with Tissue, cell and Cells in tissue types enrichment ratios
all_enrichment_ratios<-readr::read_tsv("../extracted_data/combined_enrichment_ratios.tsv", guess_max = 1e5) %>%
  clean_names()





#************UP TO DATE PHASE DATA FROM OPENTARGETS*****************************
# ==============================================================================
# 0. SETUP & DEPENDENCIES
# ==============================================================================
library(tidyverse)
library(arrow)

# Increase timeout for large file streams
options(timeout = 1000)

download_dir <- "./open_targets_data"
if (!dir.exists(download_dir)) dir.create(download_dir)

# ==============================================================================
# 1. DOWNLOAD FROM THE EXACT 26.06 DIRECTORY PATHS
# ==============================================================================
ftp_base_url <- "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/"

# Replaced overall associations with association_by_datatype_direct
datasets_to_download <- c(
  "drug_mechanism_of_action"  = "drug_mechanism_of_action",
  "target"                       = "target",
  "drug_molecule"              = "drug_molecule",
  "drug_warning"  = "drug_warning",
  "clinical_indication" = "clinical_indication",
  "disease" = "disease",
  "evidence_clinical_precedence" = "evidence_clinical_precedence",
  "association_by_datatype_indirect" = "association_by_datatype_indirect",
  "association_by_datatype_direct" = "association_by_datatype_direct",
  "association_by_overall_direct" = "association_overall_direct"
)

cat("--- Downloading datasets via open page directory index ---\n")

for (local_name in names(datasets_to_download)) {
  ftp_folder <- datasets_to_download[local_name]
  local_target_dir <- file.path(download_dir, local_name)
  if (!dir.exists(local_target_dir)) dir.create(local_target_dir, recursive = TRUE)
  
  if (length(list.files(local_target_dir)) == 0) {
    cat("Fetching file listings for:", ftp_folder, "...\n")
    
    folder_url <- paste0(ftp_base_url, ftp_folder, "/")
    temp_html <- tempfile(fileext = ".html")
    
    tryCatch({
      download.file(folder_url, destfile = temp_html, quiet = TRUE, mode = "wb")
      html_lines <- readLines(temp_html, warn = FALSE)
      
      # UPDATED REGEX: Matches both chunked "part-*" files AND single "*.parquet" files
      parquet_files <- stringr::str_extract_all(
        html_lines, 
        "(part-[0-9a-zA-Z\\-]+([^\\\"' ]*)\\.parquet)|([a-zA-Z0-9_\\-]+\\.parquet)"
      ) %>% 
        unlist() %>% 
        unique()
      
      # Clean up any HTML tag remnants that might get caught in the regex
      parquet_files <- parquet_files[!stringr::str_detect(parquet_files, "<|>|=|\"")]
      
      if (length(parquet_files) == 0) {
        cat("⚠️ Warning: No parquet files found in the index for:", ftp_folder, "\n")
        next
      }
      
      cat("  -> Found", length(parquet_files), "file(s). Downloading...\n")
      
      for (file_name in parquet_files) {
        file_url <- paste0(folder_url, file_name)
        dest_file <- file.path(local_target_dir, file_name)
        
        cat("     Downloading:", file_name, "\n")
        download.file(file_url, destfile = dest_file, quiet = TRUE, mode = "wb")
      }
      cat("  ✔️ Finished syncing:", local_name, "\n")
    }, error = function(e) {
      cat("❌ Scraping failed for folder:", ftp_folder, "| Error:", e$message, "\n")
    })
  } else {
    cat("Already exists locally:", local_name, "\n")
  }
}

# ==============================================================================
# Add OPENTARGETS data to my data
# ==============================================================================
cat("\n--- Loading datasets into arrow framework ---\n")
drug_moa            <- arrow::open_dataset(file.path(download_dir, "drug_mechanism_of_action"))
target              <- arrow::open_dataset(file.path(download_dir, "target"))
drug_molecule       <- arrow::open_dataset(file.path(download_dir, "drug_molecule"))
drug_warning        <- arrow::open_dataset(file.path(download_dir, "drug_warning"))
clinical_indication <- arrow::open_dataset(file.path(download_dir, "clinical_indication"))
disease             <- arrow::open_dataset(file.path(download_dir, "disease"))
datatype_indirect <- arrow::open_dataset(file.path(download_dir, "association_by_datatype_indirect"))
datatype_direct     <- arrow::open_dataset(file.path(download_dir, "association_by_datatype_direct"))
overall_direct <- arrow::open_dataset(file.path(download_dir, "association_by_overall_direct"))


cat("--- Running mapping logic & transformations ---\n")

# 2a. Map Drugs to Target Genes using drug_moa directly 
drug_target_table <- drug_moa %>%
  collect() %>%
  # Unnest both the targets (Ensembl IDs) and chemblIds arrays
  tidyr::unnest_longer(targets) %>%
  tidyr::unnest_longer(chemblIds) %>%
  dplyr::transmute(
    gene = targets,
    drug_id = chemblIds
  ) %>%
  filter(!is.na(gene), !is.na(drug_id)) %>%
  distinct() %>%
  # Join your custom ratios
  left_join(all_enrichment_ratios, by = "gene")

# 2b. Add drug names from drug_molecule as a metadata layer
drug_lookup <- drug_molecule %>%
  dplyr::select(drug_id = id, drug_name = name) %>%
  collect()

drug_target_table <- drug_target_table %>%
  left_join(drug_lookup, by = "drug_id")


# ==============================================================================
# 2c. Extract Indication-level Clinical Stages & FILTER to Valid Diseases
# ==============================================================================
cat("--- Processing clinical indications and filtering out non-diseases ---\n")

# 1. Pull the list of valid, canonical disease IDs from your 'disease' metadata
valid_disease_ids <- disease %>%
  dplyr::select(disease_id = id) %>%
  collect() %>%
  distinct()

# 2. Extract and filter clinical indications
clinical_indications_table <- clinical_indication %>%
  dplyr::select(
    drug_id = "drugId",
    disease_id = "diseaseId", 
    indication_max_phase = "maxClinicalStage"
  ) %>%
  collect() %>%
  filter(!is.na(drug_id), !is.na(disease_id)) %>%
  # CRITICAL INNER JOIN: Instantly drops HP, MP, GO, OBA, etc.
  dplyr::inner_join(valid_disease_ids, by = "disease_id") %>%
  distinct()


# ==============================================================================
# Process Drug Warnings Natively 
# ==============================================================================
cat("--- Processing drug warnings ---\n")

raw_drug_warnings <- drug_warning %>%
  collect() %>%
  tidyr::unnest_longer(chemblIds) %>%
  dplyr::select(
    drug_id = chemblIds,
    warning_disease_id = efoId, # Keep it distinct from your analytical disease_id
    warning_type = warningType
  ) %>%
  dplyr::filter(warning_type == "Withdrawn") %>%
  dplyr::distinct()

# 1. Global Withdrawals (No specific disease associated with the withdrawal in ChEMBL)
global_withdrawals <- raw_drug_warnings %>%
  dplyr::filter(is.na(warning_disease_id)) %>%
  dplyr::select(drug_id) %>%
  dplyr::mutate(is_globally_withdrawn = TRUE) %>%
  dplyr::distinct()

# 2. Local/Indication-Specific Withdrawals
# Instead of joining on the raw 'warning_disease_id' (which contains messy HP/GO IDs),
# we flag the drug_ids that have indication-specific withdrawals.
withdrawn_drug_ids <- raw_drug_warnings %>%
  dplyr::filter(!is.na(warning_disease_id)) %>%
  dplyr::select(drug_id) %>%
  dplyr::mutate(has_local_withdrawal = TRUE) %>%
  dplyr::distinct()

# ==============================================================================
# 2d. Fetch Datatype-specific Scores & Filter (INDIRECT VERSION)
# ==============================================================================
cat("--- Fetching and filtering indirect datatype-specific association scores ---\n")

# 1. Define valid human disease prefixes
clinical_disease_prefixes <- c("MONDO", "EFO", "NCIT", "Orphanet", "OTAR")

# 2. Get the disease metadata and restrict strictly to clinical disease prefixes
valid_disease_ids <- disease %>%
  dplyr::select(disease_id = id) %>%
  collect() %>%
  
# Uncomment to use only selected disease IDs
  
#  dplyr::filter(sub("_.*", "", disease_id) %in% clinical_disease_prefixes) %>%
  distinct()

# 3. Load, filter, and pivot the INDIRECT association scores
ot_datatype_scores <- datatype_indirect %>%
  dplyr::select(
    gene = targetId, 
    disease_id = diseaseId, 
    datatype_id = aggregationValue,  # aggregationValue holds the datatype ID!
    score = associationScore
  ) %>%
  collect() %>%
  # Filter out non-disease nodes early
  dplyr::semi_join(valid_disease_ids, by = "disease_id") %>%
  # Exclude known_drug to bypass target variable leakage
  dplyr::filter(datatype_id != "known_drug") %>%
  # Pivot wide
  tidyr::pivot_wider(
    names_from = datatype_id, 
    values_from = score, 
    values_fill = 0 
  )


# ==============================================================================
# 2e. Merge Everything Safely
# ==============================================================================
cat("--- Merging target, indication, datatype scores, and warnings ---\n")

final_analysis_df <- drug_target_table %>%
  # 1. Join clean indications (This establishes your high-quality MONDO/EFO disease keys)
  dplyr::inner_join(clinical_indications_table, by = "drug_id", relationship = "many-to-many") %>%
  
  # 2. Join biological datatype scores using the clean disease keys
  dplyr::left_join(ot_datatype_scores, by = c("gene", "disease_id")) %>%
  
  # 3. Join global withdrawals (drug level)
  dplyr::left_join(global_withdrawals, by = "drug_id") %>%
  
  # 4. Join indication-level withdrawal flag (drug level)
  dplyr::left_join(withdrawn_drug_ids, by = "drug_id") %>%
  
  dplyr::mutate(
    is_globally_withdrawn = tidyr::replace_na(is_globally_withdrawn, FALSE),
    has_local_withdrawal  = tidyr::replace_na(has_local_withdrawal, FALSE),
    
    # A drug is considered withdrawn if it's globally withdrawn OR has local withdrawal flags
    is_withdrawn          = is_globally_withdrawn | has_local_withdrawal,
    
    # Target classification: 1 ONLY if phase is "APPROVAL" AND not withdrawn
    is_approved = ifelse(indication_max_phase == "APPROVAL" & !is_withdrawn, 1, 0)
  )

# Remove all rows without a clinical phase
final_analysis_df <- final_analysis_df %>%
  filter(indication_max_phase != "UNKNOWN")

# Count the percentage of drugs with no association data

total_rows <- nrow(final_analysis_df)
completely_unmatched <- final_analysis_df %>%
  filter(is.na(genetic_association) & is.na(somatic_mutation) & is.na(literature)) %>%
  nrow()

cat("number of rows with no association data", (completely_unmatched ) , "\n")
cat( (completely_unmatched / total_rows) * 100, "%\n")

# ==============================================================================
# 3. CONVERT ASSOCIATION SCORE NAs TO 0s
# ==============================================================================
cat("--- Standardizing remaining NAs to 0 for modeling ---\n")

# Dynamically identify the pivoted datatype columns in your merged dataframe
datatype_cols <- colnames(final_analysis_df)[
  colnames(final_analysis_df) %in% c(
    "genetic_association", "somatic_mutation", "literature","genetic_literature", 
    "known_drug", "affected_pathway", "animal_model", 
    "rna_expression", "somatic_alteration"
  )
]

final_analysis_df <- final_analysis_df %>%
  dplyr::mutate(
    # Ensure drug ratios are numeric and default to 0
    tissue_ratio         = as.numeric(tidyr::replace_na(tissue_ratio, 0)),
    cell_ratio           = as.numeric(tidyr::replace_na(cell_ratio, 0)),
    cell_in_tissue_ratio = as.numeric(tidyr::replace_na(cell_in_tissue_ratio, 0)),
    
    # Handle drug warnings/withdrawal logicals
    is_globally_withdrawn = tidyr::replace_na(is_globally_withdrawn, FALSE),
    has_local_withdrawal  = tidyr::replace_na(has_local_withdrawal, FALSE),
    is_withdrawn          = is_globally_withdrawn | has_local_withdrawal
  ) %>%
  # Explicitly replace all NAs in the datatype score columns with 0
  dplyr::mutate(
    across(
      dplyr::any_of(datatype_cols), 
      ~as.numeric(tidyr::replace_na(., 0))
    )
  )













# ==============================================================================
# Model
# ==============================================================================

#-----------------little data check---------------------------------------------
# ==================================
# 1. Checking correlation
# ==================================


# First I am checking the correlation

cor(
  final_analysis_df %>%
    select(tissue_ratio,
           cell_ratio,
           cell_in_tissue_ratio),
  use = "complete.obs"
)


# Check for skewed data
library(e1071)

data.frame(
  variable = c(
    "tissue_ratio",
    "cell_ratio",
    "cell_in_tissue_ratio"
  ),
  skewness = c(
    skewness(final_analysis_df$tissue_ratio, na.rm = TRUE),
    skewness(final_analysis_df$cell_ratio, na.rm = TRUE),
    skewness(final_analysis_df$cell_in_tissue_ratio, na.rm = TRUE)
  )
)

# plot
library(ggplot2)
library(patchwork)

p1 <- ggplot(final_analysis_df, aes(tissue_ratio)) +
  geom_histogram(bins = 50)

p2 <- ggplot(final_analysis_df, aes(cell_ratio)) +
  geom_histogram(bins = 50)

p3 <- ggplot(final_analysis_df, aes(cell_in_tissue_ratio)) +
  geom_histogram(bins = 50)

p1 + p2 + p3

#box plots

p1 <- ggplot(final_analysis_df, aes(y = tissue_ratio)) +
  geom_boxplot()

p2 <- ggplot(final_analysis_df, aes(y = cell_ratio)) +
  geom_boxplot()

p3 <- ggplot(final_analysis_df, aes(y = cell_in_tissue_ratio)) +
  geom_boxplot()

p1 + p2 + p3

# inspect the tails
final_analysis_df %>%
  summarise(
    tissue_p95 = quantile(tissue_ratio, 0.95),
    tissue_p99 = quantile(tissue_ratio, 0.99),
    tissue_max = max(tissue_ratio),
    
    cell_p95 = quantile(cell_ratio, 0.95),
    cell_p99 = quantile(cell_ratio, 0.99),
    cell_max = max(cell_ratio),
    
    cit_p95 = quantile(cell_in_tissue_ratio, 0.95),
    cit_p99 = quantile(cell_in_tissue_ratio, 0.99),
    cit_max = max(cell_in_tissue_ratio)
  )

# Check for extremes
final_analysis_df %>%
  select(
    gene,
    tissue_ratio,
    cell_ratio,
    cell_in_tissue_ratio
  ) %>%
  distinct() %>%
  arrange(desc(cell_ratio)) %>%
  head(20)

# Calculate how many observations are sitting at the cap:
final_analysis_df %>%
  summarise(
    tissue_cap = mean(tissue_ratio == 11),
    cell_cap = mean(cell_ratio == 11)
  )

# create quantile bins
final_analysis_df <- final_analysis_df %>%
  mutate(
    cell_ratio_q = ntile(cell_ratio, 10)
  )

# summarize
final_analysis_df %>%
  group_by(cell_ratio_q) %>%
  summarise(
    mean_cell_ratio = mean(cell_ratio),
    approval_rate = mean(is_approved)
  )

#-------------------------------------------------------------------------------

# Interpretation:
#   
# tissue_ratio is relatively independent.
# cell_ratio and cell_in_tissue_ratio are strongly correlated (r = 0.80).
# I would not put both cell_ratio and cell_in_tissue_ratio in the same model initially.
# I'd consider:
# 
# Model A: tissue_ratio + cell_ratio
# Model B: tissue_ratio + cell_in_tissue_ratio

# ==================================
# 2. Define question
# ==================================

# My question is 

#*****
#*Are drugs targeting genes with higher cell specificity more likely to succeed?


# But opentargets scores are for gene-disease combinations. Rows roughly represent gene-drug-disease combination
# Many rows belong to the same drug or drug-target pair

#I cannot just get mean values per drug/gene because of this. Therefore I am doing a model comparison 



# ==================================
# 3. Select model
# ==================================

# Opentargets association scores are for gene-disease combinations. My specificity scores are only for genes. 
# rows from the same drug are correlated;
# genes occur repeatedly across many rows.

# Therefore I am checking 

n_distinct(final_analysis_df$drug_id)
n_distinct(final_analysis_df$gene)
n_distinct(final_analysis_df$disease_id)


final_analysis_df %>%
  count(gene) %>%
  summarise(max = max(n), median = median(n))

final_analysis_df %>%
  count(drug_id) %>%
  summarise(max = max(n), median = median(n))

final_analysis_df %>%
  count(disease_id) %>%
  summarise(max = max(n), median = median(n))

# I have clustering at all three levels:
#   
# Drug
# Gene
# Disease
# 
# But the question is not "which levels exist?", it's:
# Which levels need to be accounted for to correctly estimate the effect of my gene-level ratios?

#I have drugs appearing up to:11,271 rows - A model without (1|drug_id) would be hard to defend. -> adding (1 | drug_id)
#and genes appear up to:2436 rows- The same value of cell_ratio is being reused many times. Therefore ading a random effect helps.-> adding (1|gene)
# Opentargets scores are already disease related. Therefore I am not using that


#Therefore I am testing these first




library(lme4)

#************************Models*************************************************

#1. Drug-only model
M_drug <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio +
    (1 | drug_id),
  data = final_analysis_df,
  family = binomial(link = "logit")
)

# 2. Drug + Gene model

M_drug_gene <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio +
    (1 | drug_id) +
    (1 | gene),
  data = final_analysis_df,
  family = binomial(link = "logit")
)

# 3. Drug + Gene + Disease model

M_drug_gene_disease <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio +
    (1 | drug_id) +
    (1 | gene) +
    (1 | disease_id),
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#Compare model fit

AIC(
  M_drug,
  M_drug_gene,
  M_drug_gene_disease
)


BIC(
  M_drug,
  M_drug_gene,
  M_drug_gene_disease
)

# Likelihood ratio tests

anova(
  M_drug,
  M_drug_gene,
  test = "Chisq"
)


anova(
  M_drug_gene,
  M_drug_gene_disease,
  test = "Chisq"
)


# Examine variance explained by each random effect

as.data.frame(
  VarCorr(M_drug_gene_disease)
)

#-------------------------------Conclusions-------------------------------------

#The random effects appear to rank as: Disease >>> Drug >>> Gene

# -----Disease
# M_drug_gene vs M_drug_gene_disease
# AIC:79114 -> 55352
# χ² = 23764
# p < 2e-16
# Knowing which disease we're talking about strongly affects approval odds.

# -----Drug
# variance = 12.4.Also substantial.
#Some drugs repeatedly appear in successful indications.

# -----Gene
# M_drug vs M_drug_gene
# ΔAIC = 79134.7 - 79113.9= 20.8
#
# LRT χ² = 22.8
# p = 1.8e-06
#
# Accounting for repeated observations from the same gene significantly improves model fit.

# -----Random effect variance
# disease_id 27.19
# drug_id 12.35
# gene 0.006
# Different diseases have very different baseline success probabilities.

# -----Effect on hypothesis
# My hypothesis is : Higher cell specificity (cell_ratio) increases the odds of approval.
# Not the same as : Different genes have different approval rates
# Those are different questions.A tiny gene random-effect variance does not mean cell_ratio is unimportant.
# It only means: After accounting for cell_ratio, tissue_ratio, evidence scores, drug, and disease, 
# there is little remaining unexplained variation attributable to gene identity.

# -------------------------------Comment----------------------------------------

# I am going to use all (1|drug_id)+ (1|gene)+ (1|disease_id)

#------------------------Defining model formulas--------------------------------

# Predictors from opentargets

ot_predictors <- c(
  "genetic_association", 
  "somatic_mutation", 
  "affected_pathway", 
  "rna_expression", 
  "animal_model", 
  "literature"
)

# Selected random predictors
random_predictors<-c(
  "(1|drug_id)",
  "(1|gene)",
  "(1|disease_id)"
)

all_predictors<-paste(c(ot_predictors, random_predictors), collapse = " + ")


# Base formula
base_DGD_formula<-as.formula(paste("is_approved ~", all_predictors))

#---------------Running models for comparison-----------------------------------
library(lme4)

# First get an idea about runtime
# system.time({
#   
#   M_drug_gene_disease <- glmer(
#     is_approved ~
#       genetic_association +
#       somatic_mutation +
#       affected_pathway +
#       rna_expression +
#       animal_model +
#       literature +
#       tissue_ratio +
#       cell_ratio +
# #      (1 | drug_id) +
# #      (1 | gene) +
#       (1 | disease_id),
#     data = final_analysis_df,
#     family = binomial(link = "logit"),
#     control = glmerControl(
#       optimizer = "bobyqa",
#       optCtrl = list(maxfun = 1e5)
#     )
#   )
#   
# })


# 0.Baseline model

base_DGD <- glmer(
  base_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(base_DGD,file = "./trained_models/base_DGD.rds")
#read theat
test<-readRDS("./trained_models/base_DGD.rds")

# 1.Baseline with tissue ratio

tissue_DGD_formula <- update(base_DGD_formula,. ~ . + 
                               tissue_ratio)

tissue_DGD <- glmer(
  tissue_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_DGD,file = "./trained_models/tissue_DGD.rds")

# 2. Baseline with cell ratio

cell_DGD_formula <- update(base_DGD_formula,. ~ . +
                            cell_ratio)

cell_DGD <- glmer(
  cell_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_DGD,file = "./trained_models/cell_DGD.rds")

# 3. Baseline with cell in tissue ratio

cell_in_tissue_DGD_formula <- update(base_DGD_formula,. ~ . +
                                       cell_in_tissue_ratio)

cell_in_tissue_DGD <- glmer(
  cell_in_tissue_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_in_tissue_DGD,file = "./trained_models/cell_in_tissue_DGD.rds")

# 4. Baseline with cell ratio and tissue ratio

tissue_cell_DGD_formula <- update(base_DGD_formula,. ~ . + 
                                    tissue_ratio + cell_ratio)

tissue_cell_DGD <- glmer(
  tissue_cell_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cell_DGD,file = "./trained_models/tissue_cell_DGD.rds")

# 5. Baseline with cell in tissue ratio and tissue ratio

tissue_cellintissue_DGD_formula <- update(base_DGD_formula,. ~ . + 
                                            tissue_ratio + cell_in_tissue_ratio)

tissue_cellintissue_DGD <- glmer(
  tissue_cellintissue_DGD_formula,
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cellintissue_DGD,file = "./trained_models/tissue_cellintissue_DGD.rds")

#-----------------------------Compare models------------------------------------

# Step 1: Which individual ratio is strongest?
# Compare AIC:

AIC(
  base_DGD,
  tissue_DGD,
  cell_DGD,
  cell_in_tissue_DGD
)

# Cell DGD looks vastly superior compared to other ratios

# Step 2: Does each ratio add information beyond Open Targets?

anova(base_DGD, tissue_DGD, test = "Chisq")

anova(base_DGD, cell_DGD, test = "Chisq")

anova(base_DGD, cell_in_tissue_DGD, test = "Chisq")

# Step 3: Does tissue add anything beyond cell?

anova(
  cell_DGD,
  tissue_cell_DGD,
  test = "Chisq"
)

# Step 4: Does cell add anything beyond tissue?

anova(
  tissue_DGD,
  tissue_cell_DGD,
  test = "Chisq"
)

# Step 5: Repeat for cell_in_tissue_ratio

# Does tissue add beyond cell_in_tissue?
anova(
  cell_in_tissue_DGD,
  tissue_cellintissue_DGD,
  test = "Chisq"
)

# Does cell_in_tissue add beyond tissue?

anova(
  tissue_DGD,
  tissue_cellintissue_DGD,
  test = "Chisq"
)

# Step 6: Extract odds ratios

exp(
  cbind(
    OR = fixef(tissue_cell_DGD),
    confint(
      tissue_cell_DGD,
      parm = "beta_",
      method = "Wald"
    )
  )
)


# ---------------------Summary of results---------------------------------------

# Create summaries

model_summary <- tibble(
  model = c(
    "Base",
    "Tissue",
    "Cell",
    "Cell-in-tissue"
  ),
  AIC = c(
    AIC(base_DGD),
    AIC(tissue_DGD),
    AIC(cell_DGD),
    AIC(cell_in_tissue_DGD)
  )
)

model_summary %>%
  arrange(AIC)


# Plot odds ratios
library(broom.mixed)

coef_df <- tidy(
  cell_DGD,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

coef_df

# Keep only meaningful variables:
coef_df <- coef_df %>%
  filter(term != "(Intercept)")

# Plot
ggplot(
  coef_df,
  aes(
    x = estimate,
    y = reorder(term, estimate)
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.2
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  scale_x_log10() +
  labs(
    x = "Odds Ratio",
    y = "",
    title = "Effect of predictors on drug approval odds"
  ) +
  theme_bw()


# ΔAIC plot

comparison_df <- tibble(
  model = c(
    "Base",
    "Tissue ratio",
    "Cell ratio",
    "Cell-in-tissue ratio"
  ),
  AIC = c(
    AIC(base_DGD),
    AIC(tissue_DGD),
    AIC(cell_DGD),
    AIC(cell_in_tissue_DGD)
  )
) %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  )

ggplot(
  comparison_df,
  aes(
    x = reorder(model, delta_AIC),
    y = delta_AIC
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "",
    y = "ΔAIC from best model",
    title = "Predictive contribution of custom ratios"
  ) +
  theme_bw()


# A concise summary table for the manuscript

summary_results <- tibble(
  Ratio = c(
    "tissue_ratio",
    "cell_ratio",
    "cell_in_tissue_ratio"
  ),
  Delta_AIC = c(
    AIC(tissue_DGD) - AIC(base_DGD),
    AIC(cell_DGD) - AIC(base_DGD),
    AIC(cell_in_tissue_DGD) - AIC(base_DGD)
  ),
  LRT_p = c(
    0.0004991,
    4.484e-14,
    1.142e-13
  )
)

summary_results


#---------------------------testing linearity-----------------------------------

library(splines)


# Nonlinear spline model
cell_DGD_spline <- glmer(
  is_approved ~
    genetic_association +
    somatic_mutation +
    affected_pathway +
    rna_expression +
    animal_model +
    literature +
    ns(cell_ratio, df = 4) +
    (1 | drug_id) +
    (1 | gene) +
    (1 | disease_id),
  data = final_analysis_df,
  family = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

#save model for later
saveRDS(cell_DGD_spline,file = "./trained_models/cell_DGD_spline.rds")

# Compare the models

anova(
  cell_DGD,
  cell_DGD_spline,
  test = "Chisq"
)

AIC(
  cell_DGD,
  cell_DGD_spline
)

BIC(
  cell_DGD,
  cell_DGD_spline
)


# Visualize the fitted relationship

pred <- ggpredict(
  cell_DGD_spline,
  terms = "cell_ratio [all]"
)

ggplot() +
  geom_histogram(
    data = final_analysis_df,
    aes(x = cell_ratio, y = after_stat(density)),
    bins = 50,
    alpha = 0.3,
    fill = "grey70"
  ) +
  geom_line(
    data = pred,
    aes(x = x, y = predicted),
    colour = "blue",
    linewidth = 1.2
  ) +
  labs(
    x = "Cell ratio",
    y = "Predicted approval probability"
  ) +
  theme_bw()

#*******************************************************************************
#---------------Running models for oncology ------------------------------------
#*******************************************************************************

# I checked parents of cancer diseases in above data 
#MONDO_0004992 is the MONDO term for cancer/neoplasm and should appear in the ancestors field for cancer diseases
#https://www.ebi.ac.uk/ols4/ontologies/mondo/classes/http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMONDO_0004992
oncology_diseases <- disease %>%
  collect() %>%
  filter(
    purrr::map_lgl(
      ancestors,
      ~ "MONDO_0004992" %in% .x
    )
  ) %>%
  dplyr::select(disease_id = id)



oncology_drugs <- clinical_indication %>%
  collect() %>%
  inner_join(
    oncology_diseases,
    by = c("diseaseId" = "disease_id")
  ) %>%
  distinct(drugId) %>%
  mutate(
    is_oncology_drug = TRUE
  )


final_analysis_df <- final_analysis_df %>%
  left_join(
    oncology_drugs,
    by = c("drug_id" = "drugId")
  ) %>%
  mutate(
    is_oncology_drug = dplyr::coalesce(
      is_oncology_drug,
      FALSE
    )
  )

# Separate oncology and non oncology drug data

oncology_drug_data <- final_analysis_df[final_analysis_df$is_oncology_drug == TRUE, ]

non_oncology_drug_data <- final_analysis_df[final_analysis_df$is_oncology_drug == FALSE, ]





#--------------------------------Oncology Drugs---------------------------------

# Base formula
base_DGD_oncology_formula<-as.formula(paste("is_approved ~", all_predictors))


# 0.Baseline model

base_DGD_oncology <- glmer(
  base_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(base_DGD_oncology,file = "./trained_models/base_DGD_oncology.rds")
#read theat
test<-readRDS("./trained_models/base_DGD_oncology.rds")

# 1.Baseline with tissue ratio

tissue_DGD_oncology_formula <- update(base_DGD_oncology_formula,. ~ . + 
                               tissue_ratio)

tissue_DGD_oncology <- glmer(
  tissue_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_DGD_oncology,file = "./trained_models/tissue_DGD_oncology.rds")

# 2. Baseline with cell ratio

cell_DGD_oncology_formula <- update(base_DGD_oncology_formula,. ~ . +
                             cell_ratio)

cell_DGD_oncology <- glmer(
  cell_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_DGD_oncology,file = "./trained_models/cell_DGD_oncology.rds")

# 3. Baseline with cell in tissue ratio

cell_in_tissue_DGD_oncology_formula <- update(base_DGD_oncology_formula,. ~ . +
                                       cell_in_tissue_ratio)

cell_in_tissue_DGD_oncology <- glmer(
  cell_in_tissue_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_in_tissue_DGD_oncology,file = "./trained_models/cell_in_tissue_DGD_oncology.rds")

# 4. Baseline with cell ratio and tissue ratio

tissue_cell_DGD_oncology_formula <- update(base_DGD_oncology_formula,. ~ . + 
                                    tissue_ratio + cell_ratio)

tissue_cell_DGD_oncology <- glmer(
  tissue_cell_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cell_DGD_oncology,file = "./trained_models/tissue_cell_DGD_oncology.rds")

# 5. Baseline with cell in tissue ratio and tissue ratio

tissue_cellintissue_DGD_oncology_formula <- update(base_DGD_oncology_formula,. ~ . + 
                                            tissue_ratio + cell_in_tissue_ratio)

tissue_cellintissue_DGD_oncology <- glmer(
  tissue_cellintissue_DGD_oncology_formula,
  data = oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cellintissue_DGD_oncology,file = "./trained_models/tissue_cellintissue_DGD_oncology.rds")


# -----Plot

# Plot odds ratios
library(broom.mixed)

coef_df_onco <- tidy(
  cell_DGD_oncology,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

coef_df_onco

# Keep only meaningful variables:
coef_df_onco <- coef_df_onco %>%
  filter(term != "(Intercept)")

# Plot
ggplot(
  coef_df_onco,
  aes(
    x = estimate,
    y = reorder(term, estimate)
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.2
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  scale_x_log10() +
  labs(
    x = "Odds Ratio",
    y = "",
    title = "Effect of predictors on drug approval odds-oncology Drugs"
  ) +
  theme_bw()


# ΔAIC plot

comparison_df <- tibble(
  model = c(
    "Base",
    "Tissue ratio",
    "Cell ratio",
    "Cell-in-tissue ratio"
  ),
  AIC = c(
    AIC(base_DGD_oncology),
    AIC(tissue_DGD_oncology),
    AIC(cell_DGD_oncology),
    AIC(cell_in_tissue_DGD_oncology)
  )
) %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  )

ggplot(
  comparison_df,
  aes(
    x = reorder(model, delta_AIC),
    y = delta_AIC
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "",
    y = "ΔAIC from best model",
    title = "Predictive contribution of custom ratios- oncology Drugs"
  ) +
  theme_bw()



#--------------------------------Non Oncology Drugs---------------------------------

# Base formula
base_DGD_non_oncology_formula<-as.formula(paste("is_approved ~", all_predictors))


# 0.Baseline model

base_DGD_non_oncology <- glmer(
  base_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(base_DGD_non_oncology,file = "./trained_models/base_DGD_non_oncology.rds")
#read theat
test<-readRDS("./trained_models/base_DGD_non_oncology.rds")

# 1.Baseline with tissue ratio

tissue_DGD_non_oncology_formula <- update(base_DGD_non_oncology_formula,. ~ . + 
                                        tissue_ratio)

tissue_DGD_non_oncology <- glmer(
  tissue_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_DGD_non_oncology,file = "./trained_models/tissue_DGD_non_oncology.rds")

# 2. Baseline with cell ratio

cell_DGD_non_oncology_formula <- update(base_DGD_non_oncology_formula,. ~ . +
                                      cell_ratio)

cell_DGD_non_oncology <- glmer(
  cell_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_DGD_non_oncology,file = "./trained_models/cell_DGD_non_oncology.rds")

# 3. Baseline with cell in tissue ratio

cell_in_tissue_DGD_non_oncology_formula <- update(base_DGD_non_oncology_formula,. ~ . +
                                                cell_in_tissue_ratio)

cell_in_tissue_DGD_non_oncology <- glmer(
  cell_in_tissue_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(cell_in_tissue_DGD_non_oncology,file = "./trained_models/cell_in_tissue_DGD_non_oncology.rds")

# 4. Baseline with cell ratio and tissue ratio

tissue_cell_DGD_non_oncology_formula <- update(base_DGD_non_oncology_formula,. ~ . + 
                                             tissue_ratio + cell_ratio)

tissue_cell_DGD_non_oncology <- glmer(
  tissue_cell_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cell_DGD_non_oncology,file = "./trained_models/tissue_cell_DGD_non_oncology.rds")

# 5. Baseline with cell in tissue ratio and tissue ratio

tissue_cellintissue_DGD_non_oncology_formula <- update(base_DGD_non_oncology_formula,. ~ . + 
                                                     tissue_ratio + cell_in_tissue_ratio)

tissue_cellintissue_DGD_non_oncology <- glmer(
  tissue_cellintissue_DGD_non_oncology_formula,
  data = non_oncology_drug_data,
  family = binomial(link = "logit")
)

#save model for later
saveRDS(tissue_cellintissue_DGD_non_oncology,file = "./trained_models/tissue_cellintissue_DGD_non_oncology.rds")


# -----Plot

# Plot odds ratios
library(broom.mixed)

coef_df_non_onco <- tidy(
  cell_DGD_non_oncology,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

coef_df_non_onco

# Keep only meaningful variables:
coef_df_non_onco <- coef_df_non_onco %>%
  filter(term != "(Intercept)")

# Plot
ggplot(
  coef_df_non_onco,
  aes(
    x = estimate,
    y = reorder(term, estimate)
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.2
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  scale_x_log10() +
  labs(
    x = "Odds Ratio",
    y = "",
    title = "Effect of predictors on drug approval odds-non_oncology Drugs"
  ) +
  theme_bw()


# ΔAIC plot

comparison_df <- tibble(
  model = c(
    "Base",
    "Tissue ratio",
    "Cell ratio",
    "Cell-in-tissue ratio"
  ),
  AIC = c(
    AIC(base_DGD_non_oncology),
    AIC(tissue_DGD_non_oncology),
    AIC(cell_DGD_non_oncology),
    AIC(cell_in_tissue_DGD_non_oncology)
  )
) %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  )

ggplot(
  comparison_df,
  aes(
    x = reorder(model, delta_AIC),
    y = delta_AIC
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "",
    y = "ΔAIC from best model",
    title = "Predictive contribution of custom ratios- non_oncology Drugs"
  ) +
  theme_bw()


















ggplot(pred_cell,
       aes(x = x, y = predicted)) +
  geom_line(size = 1.2, colour = "blue") +
  geom_ribbon(
    aes(ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  labs(
    x = "Cell specificity ratio",
    y = "Predicted probability of approval",
    title = "Relationship between cell specificity and drug approval"
  ) +
  theme_bw()










# Marginal effects plot
library(ggeffects)

pred <- ggpredict(
  cell_DGD,
  terms = "cell_ratio"
)

plot(pred)


















# Base model: Open Targets evidence only
base_model <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    (1 | drug_id),
  (1 | gene),
  data = final_analysis_df,
  family = binomial(link = "logit"))

# Base model + tissue ratio
base_with_tissue_model <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    (1 | drug_id),
  (1 | gene),
  data = final_analysis_df,
  family = binomial(link = "logit")
)

# Base model + cell ratio

base_with_cell_model <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    cell_ratio +
    (1 | drug_id),
  (1 | gene),
  data = final_analysis_df,
  family = binomial
)

# Base model + cell in tissue ratio

base_with_cell_in_tissue_model <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    cell_in_tissue_ratio +
    (1 | drug_id),
  (1 | gene),
  data = final_analysis_df,
  family = binomial
)


#Base model + cell + tissue ratio

base_with_cell_and_tissue_model <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio +
    (1 | drug_id),
    (1 | gene),
  data = final_analysis_df,
  family = binomial(link = "logit")
)

#*******************************************************************************


# Compare model fit systematically

AIC(
  base_model,
  base_with_tissue_model,
  base_with_cell_model,
  base_with_cell_in_tissue_model
)



BIC(
  base_model,
  base_with_tissue_model,
  base_with_cell_model,
  base_with_cell_in_tissue_model
)

# Does tissue_ratio matter

anova(
  base_model,
  base_with_tissue_model,
  test = "Chisq"
)

# Does cell_ratio matter?

anova(base_model,base_with_cell_model,test = "Chisq")

anova(base_with_cell_model, full_model,test = "Chisq")

#Does cell_ratio still matter after tissue_ratio is included?

anova(
  base_with_tissue_model,
  full_model,
  test = "Chisq"
)


#Put results into one summary table


library(broom.mixed)

tidy(
  full_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

















# compare them using a likelihood ratio test:
anova(base_model, full_model, test = "Chisq")

anova(base_model, base_with_tissue_model, test = "Chisq")



#Get odds ratios
exp(cbind(
  OR = fixef(full_model),
  confint(full_model, parm = "beta_", method = "Wald")
))


AIC(base_model, full_model)


# Calculate marginal and conditional R²

library(performance)

r2(base_model)
r2(full_model)


# Check multicollinearity
library(car)

vif_model <- glm(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio,
  data = final_analysis_df,
  family = binomial)
  

# ==================================
# 3. With both cell and cell in tissue ratios
# ==================================

full_model_cell <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_ratio +
    (1 | drug_id),
  data = final_analysis_df,
  family = binomial
)



full_model_cell_in_tissue <- glmer(
  is_approved ~
    literature +
    genetic_association +
    animal_model +
    somatic_mutation +
    genetic_literature +
    affected_pathway +
    tissue_ratio +
    cell_in_tissue_ratio +
    (1 | drug_id),
  data = final_analysis_df,
  family = binomial
)



AIC(
  base_model,
  base_with_tissue_model,
  full_model_cell,
  full_model_cell_in_tissue
)



anova(base_model,
      base_with_tissue_model,
      base_with_cell_model,
      full_model)



































# Ensure required libraries are loaded
library(dplyr)
library(car)      # For VIF multicollinerity checks
library(pROC)     # For ROC-AUC calculation

# ==============================================================================
# 1. DEFINE PREDICTORS
# ==============================================================================

# Define Open Targets biological predictors (excluding literature-related duplicates 
# like 'genetic_literature' to avoid massive multicollinearity)
ot_predictors <- c(
  
#  "tissue_ratio", 
#  "cell_ratio",
  
  
  "genetic_association", 
  "somatic_mutation", 
  "affected_pathway", 
  "rna_expression", 
  "animal_model", 
  "literature"
)

# Define your custom ratio predictors
custom_ratios <- c(
#  "tissue_ratio", 
  "cell_ratio"
#  "cell_in_tissue_ratio"
)

# Create clean formulas
formula_baseline <- as.formula(paste("is_approved ~", paste(ot_predictors, collapse = " + ")))
formula_full     <- as.formula(paste("is_approved ~", paste(c(ot_predictors, custom_ratios), collapse = " + ")))


# ==============================================================================
# 2. FIT THE MODELS
# ==============================================================================
cat("\n--- Fitting Baseline Model (Open Targets Scores Only) ---\n")
model_baseline <- glm(formula_baseline, data = final_analysis_df, family = binomial(link = "logit"))
summary_baseline <- summary(model_baseline)
print(summary_baseline)

cat("\n--- Fitting Full Model (Open Targets + Custom Ratios) ---\n")
model_full <- glm(formula_full, data = final_analysis_df, family = binomial(link = "logit"))
summary_full <- summary(model_full)
print(summary_full)


# ==============================================================================
# 3. COMPARE MODELS
# ==============================================================================
cat("\n==============================================================================\n")
cat("                       MODEL COMPARISON METRICS\n")
cat("==============================================================================\n")

# 1. Likelihood Ratio Test (LRT / Chi-Square)
# This tests whether the addition of your custom ratios significantly improves model fit
lrt_result <- anova(model_baseline, model_full, test = "Chisq")
cat("\n1. Likelihood Ratio Test (Analysis of Deviance):\n")
print(lrt_result)

# ==============================================================================
# 1. DYNAMICALLY EXTRACT AIC VALUES
# ==============================================================================
aic_baseline <- AIC(model_baseline)
aic_full     <- AIC(model_full)
aic_delta    <- aic_baseline - aic_full

# Print metrics to the console dynamically
cat(sprintf("\n2. Information Criteria:\n"))
cat(sprintf("   - Baseline Model AIC: %.2f\n", aic_baseline))
cat(sprintf("   - Full Model AIC:     %.2f\n", aic_full))
cat(sprintf("   - AIC Improvement:    %.2f (Lower is better)\n", aic_delta))


# ==============================================================================
# 2. GENERATE THE AIC PLOT DYNAMICALLY
# ==============================================================================
# Pass the calculated AIC variables directly to the data frame
aic_df <- data.frame(
  Model = c("Baseline", paste("Full Model\n(With Custom Ratios-",custom_ratios,")")),
  AIC = c(aic_baseline, aic_full)
)

# Dynamically determine the y-axis limits around your calculated values
# This sets the lower limit slightly below the lower AIC and upper limit slightly above the higher AIC
y_min <- min(aic_df$AIC) - (aic_delta * 2)
y_max <- max(aic_df$AIC) + (aic_delta * 2)

aic_plot <- ggplot(aic_df, aes(x = Model, y = AIC, fill = Model)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.85) +
  # Zoom in on the y-axis dynamically based on the model values
  coord_cartesian(ylim = c(y_min, y_max)) + 
  scale_fill_manual(values = c("grey60", "dodgerblue")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Model Fit Comparison (Akaike Information Criterion)",
    subtitle = sprintf("An AIC drop of %.2f indicates a vastly superior model fit", aic_delta),
    x = "",
    y = "AIC Value (Lower is Better)"
  ) +
  # Render the values rounded to 2 decimal places above the bars
  geom_text(aes(label = format(round(AIC, 2), big.mark = ",")), vjust = -0.5, fontface = "bold", size = 4) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(aic_plot)


# ************************  summary *****************************************
# Extract odds ratios and 95% confidence intervals
full_summary <- summary(model_full)

# Combine coefficients, odds ratios, and p-values into a clean summary table
results_table <- data.frame(
  Estimate_LogOdds = full_summary$coefficients[, "Estimate"],
  Std_Error        = full_summary$coefficients[, "Std. Error"],
  P_Value          = full_summary$coefficients[, "Pr(>|z|)"],
  Odds_Ratio       = exp(full_summary$coefficients[, "Estimate"]),
  OR_Lower_95      = exp(confint.default(model_full)[, 1]),
  OR_Upper_95      = exp(confint.default(model_full)[, 2])
)

print(round(results_table, 4))


# ************************ ROC summary *****************************************

library(pROC)

# Generate predicted probabilities
final_analysis_df$pred_prob <- predict(model_full, type = "response")

# Calculate ROC curve and AUC
roc_obj <- pROC::roc(final_analysis_df$is_approved, final_analysis_df$pred_prob)
auc_val <- pROC::auc(roc_obj)

cat(sprintf("Model ROC-AUC: %.4f\n", auc_val))



























# 3. ROC-AUC Performance (Discriminative Ability)
# Generate predictions
pred_baseline <- predict(model_baseline, type = "response")
pred_full     <- predict(model_full, type = "response")

# Calculate AUC using the correct function name: roc()
roc_baseline <- pROC::roc(final_analysis_df$is_approved, pred_baseline, quiet = TRUE)
roc_full     <- pROC::roc(final_analysis_df$is_approved, pred_full, quiet = TRUE)

auc_baseline <- auc(roc_baseline)
auc_full     <- auc(roc_full)

cat(sprintf("\n3. Predictive Power (ROC-AUC):\n"))
cat(sprintf("   - Baseline Model AUC: %.4f\n", auc_baseline))
cat(sprintf("   - Full Model AUC:     %.4f\n", auc_full))
cat(sprintf("   - AUC Delta (Boost):  +%.4f\n", auc_full - auc_baseline))


# ==============================================================================
# 4. SANITY CHECK: MULTICOLLINEARITY (VIF)
# ==============================================================================
# Ensure that your ratio columns and OT columns aren't heavily collinear (VIF > 5-10)
cat("\n4. Multicollinearity Assessment (Variance Inflation Factor):\n")
vif_values <- car::vif(model_full)
print(vif_values)


# Visualize 

# 1. Generate the ROC curves from your earlier calculated roc objects
roc_b <- pROC::roc(final_analysis_df$is_approved, pred_baseline, quiet = TRUE)
roc_f <- pROC::roc(final_analysis_df$is_approved, pred_full, quiet = TRUE)

# 2. Extract coordinates for plotting with ggplot
df_roc_baseline <- data.frame(
  Specificity = roc_b$specificities,
  Sensitivity = roc_b$sensitivities,
  Model = paste0("Baseline (OT Scores Only) - AUC: ", round(pROC::auc(roc_b), 3))
)

df_roc_full <- data.frame(
  Specificity = roc_f$specificities,
  Sensitivity = roc_f$sensitivities,
  Model = paste0("Full Model (OT + Custom Ratios) - AUC: ", round(pROC::auc(roc_f), 3))
)

# Combine datasets
plot_data <- rbind(df_roc_baseline, df_roc_full)

# 3. Create the ggplot
roc_comparison_plot <- ggplot(plot_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  # Draw y = x diagonal reference line (intercept = 0, slope = 1)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_path(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("grey60", "#1f77b4")) + # Clean contrast colors for baseline vs full
  theme_minimal(base_size = 12) +
  labs(
    title = "Receiver Operating Characteristic (ROC) Comparison",
    subtitle = "Visualizing the predictive boost gained by integrating custom cellular/tissue ratios",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(roc_comparison_plot)




















































# 2c. Extract Indication-level Clinical Stages from the clinical_indication dataset
# Using your exact clinical_indication download
clinical_indications_table <- clinical_indication %>%
  # Select using the exact column names present in your dataset
  dplyr::select(
    drug_id = "drugId",
    disease_id = "diseaseId", 
    indication_max_phase = "maxClinicalStage"
  ) %>%
  collect() %>%
  filter(!is.na(drug_id), !is.na(disease_id)) %>%
  distinct()

# 2d. Fetch Datatype-specific Scores & EXPLICITLY Drop the 'known_drug' Layer
ot_datatype_scores <- datatype_direct %>%
  dplyr::select(
    gene = targetId, 
    disease_id = diseaseId, 
    datatype_id = aggregationValue,  # aggregationValue holds the datatype ID!
    score = associationScore
  ) %>%
  collect() %>%
  # CRITICAL FILTER: Exclude known_drug to completely bypass target variable leakage
#  filter(datatype_id != "known_drug") %>%
  # Pivot columns wide so each datatype represents an independent control variable
  tidyr::pivot_wider(
    names_from = datatype_id, 
    values_from = score, 
    values_fill = 0 
  )

# ==============================================================================
# Process Drug Warnings (Identify both Global and Disease-Specific Withdrawals)
# ==============================================================================
cat("--- Processing drug warnings ---\n")

drug_warnings_clean <- drug_warning %>%
  collect() %>%
  # Unnest the chemblIds array to get individual drug IDs
  tidyr::unnest_longer(chemblIds) %>%
  dplyr::select(
    drug_id = chemblIds,
    disease_id = efoId,       # Maps to disease_id
    warning_type = warningType
  ) %>%
  # Focus explicitly on "Withdrawn" events
  dplyr::filter(warning_type == "Withdrawn") %>%
  dplyr::distinct()

# Split into global withdrawals vs. indication-specific withdrawals
global_withdrawals <- drug_warnings_clean %>%
  dplyr::filter(is.na(disease_id)) %>%
  dplyr::select(drug_id) %>%
  dplyr::mutate(is_globally_withdrawn = TRUE) %>%
  dplyr::distinct()

indication_withdrawals <- drug_warnings_clean %>%
  dplyr::filter(!is.na(disease_id)) %>%
  dplyr::select(drug_id, disease_id) %>%
  dplyr::mutate(is_locally_withdrawn = TRUE) %>%
  dplyr::distinct()


# ==============================================================================
# 2e. Merge Everything Together to Drug-Target-Disease Granularity
# ==============================================================================
cat("--- Merging target, indication, datatype scores, and warnings ---\n")

final_analysis_df <- drug_target_table %>%
  # Join indications (Many-to-Many mapping)
  dplyr::inner_join(clinical_indications_table, by = "drug_id", relationship = "many-to-many") %>%
  # Join biological datatype scores
  dplyr::left_join(ot_datatype_scores, by = c("gene", "disease_id")) %>%
  # Join global withdrawals (matching only on drug_id)
  dplyr::left_join(global_withdrawals, by = "drug_id") %>%
  # Join indication-specific withdrawals (matching on both drug_id and disease_id)
  dplyr::left_join(indication_withdrawals, by = c("drug_id", "disease_id")) %>%
  dplyr::mutate(
    # Set withdrawal flags to FALSE if no match was found
    is_globally_withdrawn = tidyr::replace_na(is_globally_withdrawn, FALSE),
    is_locally_withdrawn  = tidyr::replace_na(is_locally_withdrawn, FALSE),
    is_withdrawn          = is_globally_withdrawn | is_locally_withdrawn,
    
    # Target classification: 1 ONLY if phase is "APPROVAL" AND not withdrawn
    is_approved = ifelse(indication_max_phase == "APPROVAL" & !is_withdrawn, 1, 0),
    
    # Ensure your custom ratio columns are numeric and fill missing values with 0
    tissue_ratio         = as.numeric(tissue_ratio),
    tissue_ratio         = ifelse(is.na(tissue_ratio), 0, tissue_ratio),
    
    cell_ratio           = as.numeric(cell_ratio),
    cell_ratio           = ifelse(is.na(cell_ratio), 0, cell_ratio),
    
    cell_in_tissue_ratio = as.numeric(cell_in_tissue_ratio),
    cell_in_tissue_ratio = ifelse(is.na(cell_in_tissue_ratio), 0, cell_in_tissue_ratio)
  ) %>%
  # Safely fill remaining missing values for datatype scores with 0
  dplyr::mutate(
    across(
      where(is.numeric) & !c(tissue_ratio, cell_ratio, cell_in_tissue_ratio, is_approved), 
      ~replace_na(., 0)
    )
  )


















# ==============================================================================
# function for Calculate Harmonic Sum 
# ==============================================================================


calculate_ot_harmonic_sum <- function(score_vector) {
  # 1. Define Open Targets standard datatype weights
  weights <- c(
    
    #************This is drug phase score weight
    "known_drug"          = 1.0,
    #************************************
    
    "genetic_association" = 1.0,
    "somatic_mutation"    = 1.0,
    "affected_pathway"    = 1.0,
    "literature"          = 0.5,
    "rna_expression"      = 0.2,
    "animal_model"        = 0.2  # OT's typical weight for animal models
  )
  
  # Clean up input names to match weight keys
  names(score_vector) <- c(
    "known_drug", "genetic_association", "literature", 
    "somatic_mutation", "rna_expression", "animal_model", "affected_pathway"
  )
  
  # 2. Multiply each score by its corresponding weight
  weighted_scores <- score_vector * weights[names(score_vector)]
  
  # 3. Remove NA, NULL, or 0 values
  active_weighted <- weighted_scores[!is.na(weighted_scores) & weighted_scores > 0]
  
  if (length(active_weighted) == 0) return(0)
  
  # 4. Sort scores in descending order
  sorted_weighted <- sort(active_weighted, decreasing = TRUE)
  
  # 5. Apply the harmonic sum: s_i / i^2
  i <- seq_along(sorted_weighted)
  raw_sum <- sum(sorted_weighted / (i^2))
  
  # 6. Normalize by the theoretical maximum limit (1.644) and cap at 1.0
  normalized_score <- raw_sum / 1.644
  return(min(normalized_score, 1.0))
}














# 2. Apply it to your R DataFrame (excluding clinical precedence / known_drug)
final_analysis_df_with_harmonic_sum <- final_analysis_df %>%
  rowwise() %>%
  mutate(
    # Specify the columns you want to include in the calculation. 
    # Simply omit your clinical precedence column (often named "known_drug" or similar) from this list.
    overall_score_no_clinical = calculate_harmonic_sum(
      c_across(c(
        
        #This layer includes known drug score. keep in commented***************
        known_drug,
        
        genetic_association, 
        literature, 
        somatic_mutation, 
        rna_expression, 
        animal_model, 
        affected_pathway
      ))
    )
  ) %>%
  ungroup()

















# ==============================================================================
# 3. STATISTICAL ANALYSIS OUTPUT
# ==============================================================================
cat("\n--- Running Multi-Variable Logistic Regression ---\n")

# Dynamically construct the formula using all downloaded/pivoted biological datatypes
datatype_columns <- setdiff(names(ot_datatype_scores), c("gene", "disease_id"))
formula_string   <- paste("is_approved ~ custom_ratio +", paste(datatype_columns, collapse = " + "))

cat("Fitting model formula:", formula_string, "\n\n")

logistic_model <- glm(
  as.formula(formula_string), 
  data = final_analysis_df, 
  family = binomial(link = "logit")
)

print(summary(logistic_model))






























#********************ADD Drug-Target data from OPENTARGETS**********************

# ----- Open OPENTARGETS datasets
drug_moa <- open_dataset(
  "../opentargets_data/drug_mechanism_of_action"
)

target <- open_dataset("../opentargets_data/target")

#add drug names
drug_molecule <- open_dataset(
  "../opentargets_data/drug_molecule"
)


drug_warning <- open_dataset(
  "../opentargets_data/drug_warning"
)


clinical_indication <- open_dataset(
  "../opentargets_data/clinical_indication"
)

disease <- open_dataset(
  "../opentargets_data/disease"
)

#**********************************************

# -----
#Join gene metrics to the drug-target mapping

drug_target_table <- drug_moa %>%
  collect() %>%
  tidyr::unnest_longer(targets) %>%
  tidyr::unnest_longer(chemblIds) %>%
  dplyr::transmute(
    gene = targets,
    drug_id = chemblIds
  ) %>%
  distinct()

#Join the enrichment metrics:
drug_target_table <- drug_target_table %>%
  left_join(
    all_enrichment_ratios,
    by = "gene"
  )

#Add drug names
drug_lookup <- drug_molecule %>%
  dplyr::select(
    drug_id = id,
    drug_name = name
  ) %>%
  collect()

drug_target_table <- drug_target_table %>%
  left_join(drug_lookup, by = "drug_id")

# Add max phase data

#Create a lookup table:
phase_lookup <- drug_molecule %>%
  dplyr::select(
    drug_id = id,
    max_phase = maximumClinicalStage
  ) %>%
  collect()


drug_target_table <- drug_target_table %>%
  dplyr::left_join(
    phase_lookup,
    by = "drug_id"
  )

# Add drug warnings
warning_table <- drug_warning %>%
  collect() %>%
  tidyr::unnest_longer(chemblIds) %>%
  dplyr::transmute(
    drug_id = chemblIds,
    warningType,
    toxicityClass,
    description
  ) %>%
  distinct()


warning_summary <- warning_table %>%
  group_by(drug_id) %>%
  summarise(
    warnings = paste(
      unique(warningType),
      collapse = "; "
    ),
    n_warnings = n_distinct(warningType),
    .groups = "drop"
  )


drug_target_table <- drug_target_table %>%
  left_join(
    warning_summary,
    by = "drug_id"
  )


#*******************ADD oncology classification  data***************************

# #inspect the therapeutic area IDs attached to known cancer diseases:
# disease %>%
#   collect() %>%
#   filter(
#     grepl("cancer", name, ignore.case = TRUE)
#   ) %>%
#   dplyr::select(name, therapeuticAreas) %>%
#   head(10)
# 
# #examine one:
# x <- disease %>%
#   collect() %>%
#   filter(name == "respiratory system cancer")
# 
# x$therapeuticAreas[[1]]

# I checked parents of cancer diseases in above data 
#MONDO_0004992 is the MONDO term for cancer/neoplasm and should appear in the ancestors field for cancer diseases
#https://www.ebi.ac.uk/ols4/ontologies/mondo/classes/http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FMONDO_0004992
oncology_diseases <- disease %>%
  collect() %>%
  filter(
    purrr::map_lgl(
      ancestors,
      ~ "MONDO_0004992" %in% .x
    )
  ) %>%
  dplyr::select(disease_id = id)



oncology_drugs <- clinical_indication %>%
  collect() %>%
  inner_join(
    oncology_diseases,
    by = c("diseaseId" = "disease_id")
  ) %>%
  distinct(drugId) %>%
  mutate(
    is_oncology_drug = TRUE
  )


drug_target_table <- drug_target_table %>%
  left_join(
    oncology_drugs,
    by = c("drug_id" = "drugId")
  ) %>%
  mutate(
    is_oncology_drug = dplyr::coalesce(
      is_oncology_drug,
      FALSE
    )
  )


#**********************clean*********************************

#clean the table
drug_target_table_clean <- drug_target_table %>%
  dplyr::mutate(
    # Override phase if withdrawn
    max_phase = dplyr::if_else(
      !is.na(warnings) & grepl("Withdrawn", warnings),
      "Withdrawn",
      max_phase
    )
  ) %>%
  
  # Remove UNKNOWN drugs
  dplyr::filter(max_phase != "UNKNOWN") %>%
  
  # New binary marketed column
  dplyr::mutate(
    marketed = max_phase == "APPROVAL"
  )


#Get single and multi-target gene lists/ onco 
#add target number column
drug_target_table_clean <- drug_target_table_clean %>%
  group_by(drug_id) %>%
  mutate(single_gene_target = n_distinct(na.omit(gene_name)) == 1) %>%
  ungroup()

single_gene_target_drug_data<-drug_target_table_clean%>%
  filter(single_gene_target==TRUE)

multi_gene_target_drug_data<-drug_target_table_clean%>%
  filter(single_gene_target==FALSE)

oncology_drug_data<-drug_target_table_clean%>%
  filter(is_oncology_drug==TRUE)

non_oncology_drug_data<-drug_target_table_clean%>%
  filter(is_oncology_drug==FALSE)

#***************************ADD Association scores******************************
#*******************************************************************************

# Association scores 
association_by_datasource_direct<-open_dataset("../opentargets_data/association_by_datasource_direct/")

# -----------------------------
# Association by dataset scores
# -----------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)

# ==============================================================================
# STEP 1: Process and Pivot Open Targets Baseline Features (Gene Level)
# ==============================================================================

#*****Association scores are explained in
#https://platform-docs.opentargets.org/associations
#*****


baseline_features <- association_by_datasource_direct %>%
  dplyr::select(
    gene = targetId,
    datasourceId = aggregationValue, 
    score = associationScore
  ) %>%
  # Filter out clinical phase (leakage) and built-in expression
  dplyr::filter(!datasourceId %in% c("chembl", "expression_atlas", "progeny")) %>%
  collect() %>%
  # Pivot wide: handling multiple disease scores per gene by taking the max
  pivot_wider(
    id_cols = gene,
    names_from = datasourceId,
    values_from = score,
    values_fn = max, 
    values_fill = 0 
  )

# Extract names of all Open Targets tracks
ot_components <- colnames(baseline_features)[colnames(baseline_features) != "gene"]
all_variables <- c(ot_components, ratios_to_use)

# ==============================================================================
# STEP 2: Merge Layers and Aggregate to the DRUG Level (Fixes Pseudo-replication)
# ==============================================================================

if (multi_gene_handling_method=="max") {
  
drug_level_dataset <- drug_target_table_clean %>%
  # Define binary market success variable
  mutate(
    is_marketed = ifelse(max_phase == "APPROVAL" & !warnings %in% c("Withdrawn", "Withdrawn from market"), 1, 0)
  ) %>%
  # Bring in the Open Targets data at the gene level
  left_join(baseline_features, by = "gene") %>%
  # Replace any post-join NAs in OT tracks with 0
  mutate(across(all_of(ot_components), ~ ifelse(is.na(.), 0, .))) %>%
  # COLLAPSE TO DRUG LEVEL
  group_by(drug_id) %>%
  summarise(
    # Use baseline R indexing to safely capture market status without BiocGenerics conflicts
    marketed = is_marketed[1], 
        
    # Take the MAX score across all targets for your custom ratios
    tissue_ratio = max(tissue_ratio, na.rm = TRUE),
    cell_ratio = max(cell_ratio, na.rm = TRUE),
    cell_in_tissue_ratio = max(cell_in_tissue_ratio, na.rm = TRUE),
    
    # Take the MAX score across all targets for every Open Targets track dynamically
    across(all_of(ot_components), ~ max(.x, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # Clean up any potential -Inf or Inf values caused by max() on empty groups
  mutate(across(everything(), ~ ifelse(is.infinite(.), 0, .)))
}


if (multi_gene_handling_method=="mean") {
  
  drug_level_dataset <- drug_target_table_clean %>%
    # Define binary market success variable
    mutate(
      is_marketed = ifelse(max_phase == "APPROVAL" & !warnings %in% c("Withdrawn", "Withdrawn from market"), 1, 0)
    ) %>%
    # Bring in the Open Targets data at the gene level
    left_join(baseline_features, by = "gene") %>%
    # Replace any post-join NAs in OT tracks with 0
    mutate(across(all_of(ot_components), ~ ifelse(is.na(.), 0, .))) %>%
    # COLLAPSE TO DRUG LEVEL
    group_by(drug_id) %>%
    summarise(
      # Use baseline R indexing to safely capture market status without BiocGenerics conflicts
      marketed = is_marketed[1], 
      
      # Take the MAX score across all targets for my custom ratios
      tissue_ratio = mean(tissue_ratio, na.rm = TRUE),
      cell_ratio = mean(cell_ratio, na.rm = TRUE),
      cell_in_tissue_ratio = mean(cell_in_tissue_ratio, na.rm = TRUE),
      
      # Take the MAX score across all targets for every Open Targets track dynamically
      across(all_of(ot_components), ~ mean(.x, na.rm = TRUE))
    ) %>%
    ungroup() %>%
    # Clean up any potential -Inf or Inf values caused by max() on empty groups
    mutate(across(everything(), ~ ifelse(is.infinite(.), 0, .)))
}

# Filter for complete cases across all variables on our drug-level dataset
final_drug_dataset <- drug_level_dataset %>%
  filter(if_all(all_of(all_variables), ~ !is.na(.)))

# ==============================================================================
# STEP 3: Define Formulas and Run Overall Power Models
# ==============================================================================
formula_baseline <- as.formula(paste("marketed ~", paste(ot_components, collapse = " + ")))
formula_full     <- as.formula(paste("marketed ~", paste(c(ot_components, ratios_to_use), collapse = " + ")))

# Fit the models on true single-row-per-drug data
model_baseline_drug <- glm(formula_baseline, data = final_drug_dataset, family = binomial(link = "logit"))
model_full_drug     <- glm(formula_full,     data = final_drug_dataset, family = binomial(link = "logit"))

# Print overall diagnostic tests to terminal
print(anova(model_baseline_drug, model_full_drug, test = "Chisq"))
print(summary(model_full_drug))

# ==============================================================================
# STEP 4: Stratified Models & Beautiful Forest Plot
# ==============================================================================

# Define your subsets into a named list based on your study groupings
# (Make sure single_gene_target_drug_data, oncology_drug_data, etc. are loaded)
datasets_list <- list(
  "All drugs"           = final_drug_dataset,
  "Single-target drugs" = final_drug_dataset %>% filter(drug_id %in% single_gene_target_drug_data$drug_id),
  "Multi-target drugs"  = final_drug_dataset %>% filter(drug_id %in% multi_gene_target_drug_data$drug_id),
  "Oncology drugs"      = final_drug_dataset %>% filter(drug_id %in% oncology_drug_data$drug_id),
  "Non-Oncology drugs"  = final_drug_dataset %>% filter(drug_id %in% non_oncology_drug_data$drug_id)
)

# Run the full regression across all stratified datasets
all_measures_results <- imap_dfr(datasets_list, function(df, name) {
  if(nrow(df) == 0) return(NULL) # Skip if a subset is empty
  
  fit <- glm(formula_full, data = df, family = binomial(link = "logit"))
  
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in% measures_to_plot_in_model) %>%
    mutate(dataset = name)
})
  
  # Clean and order labels for presentation
  all_measures_results <- all_measures_results %>%
    mutate(
      predictor = case_when(
        term == "tissue_ratio" ~ "Tissue ratio",
        term == "cell_ratio" ~ "Cell ratio",
        term == "cell_in_tissue_ratio" ~ "Cell-in-tissue ratio",
        term == "europepmc" ~ "Europe PMC score"
      ),
      predictor = factor(predictor, levels = c( "Tissue ratio", "Cell ratio",
                                              "Europe PMC score",
                                               "Cell-in-tissue ratio"
                                               )),
      dataset = factor(dataset, levels = c("Single-target drugs", "Multi-target drugs", "Oncology drugs", "Non-Oncology drugs", "All drugs"))
    )
  
  # Generate the Upgraded Forest Plot
  fully_adjusted_plot <- ggplot(all_measures_results, aes(x = estimate, y = predictor)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.15, linewidth = 0.8) +
    geom_point(size = 3, colour = "#4472C4") +
    facet_wrap(~ dataset, ncol = 2) +
    scale_x_log10() +
    labs(
      title = "Independent Predictors of Drug Success",
      subtitle = "Aggregated Per Drug & Fully Adjusted for Continuous Open Targets Tracks",
      x = "Adjusted Odds Ratio (Log Scale)",
      y = NULL
    ) +
    theme_classic(base_size = 12)
  
  # Render plot
  fully_adjusted_plot
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

















#Add gene symbols
# drug_target_table <- drug_target_table %>%
#   dplyr::rename(
#     gene_ID = gene,
#     gene = gene_name
#   )

#***************ADD Mendalian disease classification from OPENTARGETS***********

g2p <- open_dataset("../opentargets_data/evidence_gene2phenotype/")

orphanet <- open_dataset("../opentargets_data/evidence_orphanet/")


# Extract disease annotations
g2p_disease_table <- g2p %>%
  dplyr::select(
    gene_name = targetFromSourceId,
    disease = diseaseFromSource
  ) %>%
  collect()

orpha_disease_table <- orphanet %>%
  dplyr::select(
    gene_name = targetFromSourceId,
    disease = diseaseFromSource
  ) %>%
  collect()

# Combine both datasets
disease_table <- bind_rows(
  g2p_disease_table,
  orpha_disease_table
) %>%
  distinct()

# Collapse diseases per gene into one string
disease_summary <- disease_table %>%
  group_by(gene_name) %>%
  summarise(
    mendelian_diseases = str_c(
      sort(unique(disease)),
      collapse = "; "
    ),
    has_mendelian_evidence = TRUE,
    .groups = "drop"
  )

# Join to main dataframe
drug_target_table <-
  drug_target_table %>%
  left_join(
    disease_summary,
    by = "gene_name"
  ) %>%
  mutate(
    has_mendelian_evidence =
      ifelse(is.na(has_mendelian_evidence), FALSE, TRUE)
  )




#****************ADD oncology classification related data***********************
#*******************************************************************************

# Association scores 
association_by_datasource_direct<-open_dataset("../opentargets_data/association_by_datasource_direct/")



# Cancer data
#Can this gene guide therapy decisions?
evidence_cancer_biomarkers<-open_dataset("../opentargets_data/evidence_cancer_biomarkers/")

#Precision oncology BIOMARKER evidence
#Is this gene involved in cancer?
evidence_cancer_gene_census<-open_dataset("../opentargets_data/evidence_cancer_gene_census/")



# -----------------------------
# Association by dataset scores
# -----------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)

# ==============================================================================
# STEP 1: Process and Pivot Open Targets Baseline Features (Gene Level)
# ==============================================================================
baseline_features <- association_by_datasource_direct %>%
  dplyr::select(
    gene = targetId,
    datasourceId = aggregationValue, 
    score = associationScore
  ) %>%
  # Filter out clinical phase (leakage) and built-in expression
  dplyr::filter(!datasourceId %in% c("chembl", "expression_atlas", "progeny")) %>%
  collect() %>%
  # Pivot wide: handling multiple disease scores per gene by taking the max
  pivot_wider(
    id_cols = gene,
    names_from = datasourceId,
    values_from = score,
    values_fn = max, 
    values_fill = 0 
  )

# Extract names of all Open Targets tracks
ot_components <- colnames(baseline_features)[colnames(baseline_features) != "gene"]
all_variables <- c(ot_components, "tissue_ratio", "cell_ratio", "cell_in_tissue_ratio")

# ==============================================================================
# STEP 2: Merge Layers and Aggregate to the DRUG Level (Fixes Pseudo-replication)
# ==============================================================================
drug_level_dataset <- drug_target_table %>%
  # Define binary market success variable
  mutate(
    is_marketed = ifelse(max_phase == "APPROVAL" & !warnings %in% c("Withdrawn", "Withdrawn from market"), 1, 0)
  ) %>%
  # Bring in the Open Targets data at the gene level
  left_join(baseline_features, by = "gene") %>%
  # Replace any post-join NAs in OT tracks with 0
  mutate(across(all_of(ot_components), ~ ifelse(is.na(.), 0, .))) %>%
  # COLLAPSE TO DRUG LEVEL
  group_by(drug_id) %>%
  summarise(
    # Use baseline R indexing to safely capture market status without BiocGenerics conflicts
    marketed = is_marketed[1], 
    
    # Take the MAX score across all targets for your custom ratios
    tissue_ratio = max(tissue_ratio, na.rm = TRUE),
    cell_ratio = max(cell_ratio, na.rm = TRUE),
    cell_in_tissue_ratio = max(cell_in_tissue_ratio, na.rm = TRUE),
    
    # Take the MAX score across all targets for every Open Targets track dynamically
    across(all_of(ot_components), ~ max(.x, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # Clean up any potential -Inf or Inf values caused by max() on empty groups
  mutate(across(everything(), ~ ifelse(is.infinite(.), 0, .)))

# Filter for complete cases across all variables on our drug-level dataset
final_drug_dataset <- drug_level_dataset %>%
  filter(if_all(all_of(all_variables), ~ !is.na(.)))

# ==============================================================================
# STEP 3: Define Formulas and Run Overall Power Models
# ==============================================================================
formula_baseline <- as.formula(paste("marketed ~", paste(ot_components, collapse = " + ")))
formula_full     <- as.formula(paste("marketed ~", paste(c(ot_components, "tissue_ratio", "cell_ratio", "cell_in_tissue_ratio"), collapse = " + ")))

# Fit the models on true single-row-per-drug data
model_baseline_drug <- glm(formula_baseline, data = final_drug_dataset, family = binomial(link = "logit"))
model_full_drug     <- glm(formula_full,     data = final_drug_dataset, family = binomial(link = "logit"))

# Print overall diagnostic tests to terminal
print(anova(model_baseline_drug, model_full_drug, test = "Chisq"))
print(summary(model_full_drug))

# ==============================================================================
# STEP 4: Stratified Models & Beautiful Forest Plot
# ==============================================================================









# Define your subsets into a named list based on your study groupings
# (Make sure single_gene_target_drug_data, oncology_drug_data, etc. are loaded)
datasets_list <- list(
  "All drugs"           = final_drug_dataset,
  "Single-target drugs" = final_drug_dataset %>% filter(drug_id %in% single_gene_target_drug_data$drug_id),
  "Multi-target drugs"  = final_drug_dataset %>% filter(drug_id %in% multi_gene_target_drug_data$drug_id),
  "Oncology drugs"      = final_drug_dataset %>% filter(drug_id %in% oncology_drug_data$drug_id),
  "Non-Oncology drugs"  = final_drug_dataset %>% filter(drug_id %in% non_oncology_drug_data$drug_id)
)

# Run the full regression across all stratified datasets
all_measures_results <- imap_dfr(datasets_list, function(df, name) {
  if(nrow(df) == 0) return(NULL) # Skip if a subset is empty
  
  fit <- glm(formula_full, data = df, family = binomial(link = "logit"))
  
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term %in= c("tissue_ratio", "cell_ratio", "cell_in_tissue_ratio")) %>%
           mutate(dataset = name)
})
  
  # Clean and order labels for presentation
  all_measures_results <- all_measures_results %>%
    mutate(
      predictor = case_when(
        term == "tissue_ratio" ~ "Tissue ratio",
        term == "cell_ratio" ~ "Cell ratio",
        term == "cell_in_tissue_ratio" ~ "Cell-in-tissue ratio"
      ),
      predictor = factor(predictor, levels = c("Cell-in-tissue ratio", "Tissue ratio", "Cell ratio")),
      dataset = factor(dataset, levels = c("Single-target drugs", "Multi-target drugs", "Oncology drugs", "Non-Oncology drugs", "All drugs"))
    )
  
  # Generate the Upgraded Forest Plot
  fully_adjusted_plot <- ggplot(all_measures_results, aes(x = estimate, y = predictor)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.15, linewidth = 0.8) +
    geom_point(size = 3, colour = "#4472C4") +
    facet_wrap(~ dataset, ncol = 2) +
    scale_x_log10() +
    labs(
      title = "Independent Predictors of Drug Success",
      subtitle = "Aggregated Per Drug & Fully Adjusted for Continuous Open Targets Tracks",
      x = "Adjusted Odds Ratio (Log Scale)",
      y = NULL
    ) +
    theme_classic(base_size = 12)
  
  # Render plot
  fully_adjusted_plot

















# -----------------------------
# Cancer driver genes
# -----------------------------

cancer_driver_table <- evidence_cancer_gene_census %>%
  dplyr::select(
    gene = targetFromSourceId
  ) %>%
  distinct() %>%
  collect() %>%
  mutate(
    has_cancer_driver_evidence = TRUE
  )

# -----------------------------
# Cancer biomarker genes
# -----------------------------

cancer_biomarker_table <- evidence_cancer_biomarkers %>%
  dplyr::select(
    gene = targetId
  ) %>%
  distinct() %>%
  collect() %>%
  mutate(
    has_cancer_biomarker_evidence = TRUE
  )

# -----------------------------
# Join into main dataframe
# -----------------------------
drug_target_table <-
  drug_target_table %>%
  left_join(
    cancer_biomarker_table,
    by = "gene"
  ) %>%
  mutate(
    has_cancer_biomarker_evidence =
      ifelse(
        is.na(has_cancer_biomarker_evidence),
        FALSE,
        TRUE
      )
  )

drug_target_table <-
  drug_target_table %>%
  left_join(
    cancer_driver_table,
    by = "gene"
  ) %>%
  mutate(
    has_cancer_driver_evidence =
      ifelse(
        is.na(has_cancer_driver_evidence),
        FALSE,
        TRUE
      )
  )



#*******************************************************************************
#*#Collapse to one row per drug

drug_level_table <- drug_target_table %>%
  group_by(drug_id, drug_name) %>%
  summarise(
    target_genes = paste(
      sort(unique(gene_name)),
      collapse = "; "
    ),
    n_targets = n_distinct(gene_name),
    
    tissue_ratio = mean(tissue_ratio, na.rm = TRUE),
    cell_ratio = mean(cell_ratio, na.rm = TRUE),
    cell_in_tissue_ratio = mean(cell_in_tissue_ratio, na.rm = TRUE),
    
    max_phase = dplyr::first(max_phase),
    warnings = dplyr::first(warnings),
    n_warnings = dplyr::first(n_warnings),
    
    is_oncology_drug = any(is_oncology_drug),
    has_cancer_biomarker_evidence = any(has_cancer_biomarker_evidence),
    has_cancer_driver_evidence = any(has_cancer_driver_evidence),
    has_mendelian_evidence = any(has_mendelian_evidence),
  
    
    .groups = "drop"
  ) %>%
  mutate(
    single_gene_target = n_targets == 1
  )

#clean the table
drug_level_table_clean <- drug_level_table %>%
  dplyr::mutate(
    # Override phase if withdrawn
    max_phase = dplyr::if_else(
      !is.na(warnings) & grepl("Withdrawn", warnings),
      "Withdrawn",
      max_phase
    )
  ) %>%
  
  # Remove UNKNOWN drugs
  dplyr::filter(max_phase != "UNKNOWN") %>%
  
  # New binary marketed column
  dplyr::mutate(
    marketed = max_phase == "APPROVAL"
  )

#*******************Look at the summary*****************************************


cols_to_summarize <- c(
  "is_oncology_drug",
  "has_cancer_biomarker_evidence",
  "has_cancer_driver_evidence",
  "has_mendelian_evidence",
  "single_gene_target",
  "marketed"
)

summary_table <- purrr::map_dfr(
  cols_to_summarize,
  \(col) {
    drug_level_table_clean %>%
      dplyr::count(value = .data[[col]]) %>%
      dplyr::mutate(
        variable = col,
        percentage = round(100 * n / sum(n), 1)
      )
  }
) %>%
  dplyr::select(variable, value, n, percentage)

cat(yellow("Summary of Data"))

cat(
  green(
    paste(capture.output(print(summary_table)), collapse = "\n")
  )
)




#*******************Functions to calculate ratios*******************************

calculate_gene_percent <- function(data,
                                   gene_col,
                                   ratio_col,
                                   thresholds,
                                   phase_column,
                                   label = "All genes") {
  
  data %>%
    dplyr::select({{ gene_col }}, {{ ratio_col }}) %>%
    distinct({{ gene_col }}, .keep_all = TRUE) %>%
    crossing(threshold = thresholds) %>%
    mutate(exceeds = {{ ratio_col }} > threshold) %>%
    group_by(threshold) %>%
    summarise(
      n_genes = n(),
      n_exceeding = sum(exceeds, na.rm = TRUE),
      percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(!!rlang::sym(phase_column) := label)
}



# ---Function to calculate ratios with phase data

calculate_gene_percentages <- function(data, thresholds, gene_col, phase_col, ratio_col) {
  
  gene_col  <- enquo(gene_col)
  phase_col <- enquo(phase_col)
  ratio_col <- enquo(ratio_col)
  
  gene_percentage_df <- data %>%
    # keep only columns we need
    dplyr::select(!!gene_col, !!phase_col, !!ratio_col) %>%
    
    # ensure one row per gene
    #distinct(!!gene_col, .keep_all = TRUE) %>%   #<=====IMPORTANT**************
    
    # expand thresholds
    crossing(threshold = thresholds) %>%
    
    # compute exceedance
    mutate(exceeds = !!ratio_col > threshold) %>%
    
    # aggregate
    group_by(!!phase_col, threshold) %>%
    summarise(
      n_genes = n(),
      n_exceeding = sum(exceeds, na.rm = TRUE),
      percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(gene_percentage_df)
}

#*******************************************************************************
#**************PLOT WITH DRUG/Gene NUMBERS *************************************
#*******************************************************************************

#--------Set value here---------------------------------------------------------
phase_column<-"max_phase"


analysis_type<-"drug"
#options "gene" "drug"

phases_used<-"marketed_or_other"
#options "all_phases" "marketed_or_other"

# Uncomment to  select threshold values.
threshold_name<-"Range of 1.5 to 5.5 best for Cell specificity"
thresholds <- seq(1.25,5.25,0.5)

# threshold_name<-"Range of 2 to 10 best for Tissue specificity"
# thresholds <- seq(1,10,1)


#Set axis names and plot title
drug_plot_title<-"Percentage of cell specific targets for different specificity levels"
x_axis_name<-"Top VS Runners Up Ratio"
y_axis_name<-"N cell specific/N total %"

#override auto y max with a value
y_max_method<-"manual"
y_max_limit<-50

one_plot_height<-3
one_plot_width<-4
#*********Change level orders according to plot directions
# ----- plot aesthetics
#Plot directions in panel
plots_along_x_axis<-"Ratio_type" #non oncology and oncology
plots_along_y_axis<-"Gene_type" #tissue_ratio cell_ratio z_ratio etc

# Ratios to use
ratios_to_use<-c("tissue_ratio","cell_ratio","cell_in_tissue_ratio")
#options "cell_ratio","cell_in_tissue_ratio","tissue_ratio"

#select datasets to plot

selected_datasets<-c("single_target_drugs","non_mendalian_single_target_drugs","non_CBM_single_target_drugs","non_CDR_single_target_drugs")

#options
# selected_datasets<-c("mean_all_target_drugs","single_target_drugs","mean_multi_target_drugs",
#                      "max_all_target_drugs","max_multi_target_drugs",
#                      "mendalian_single_target_drugs","mean_mendalian_multi_target_drugs",
#                      "non_mendalian_single_target_drugs","mean_non_mendalian_multi_target_drugs",
#                      "CBM_single_target_drugs","mean_CBM_multi_target_drugs",
#                      "non_CBM_single_target_drugs","mean_non_CBM_multi_target_drugs",
#                      "CDR_single_target_drugs","mean_CDR_multi_target_drugs",
#                      "non_CDR_single_target_drugs","mean_non_CDR_multi_target_drugs")

# may have patterns from "single_target_drugs","non_mendalian_single_target_drugs

#only selected levels
x_levels_order = ratios_to_use # uses the same order as above by default
y_levels_order<-selected_datasets # uses the same order as above by default

#-------------------------------------------------------------------------------
#*******************************************************************************
#*******************************************************************************


# First I am loading dataframes

#load all data
full_detail_drug_database<-drug_target_table_clean

paste("OPENTARGETS had info on",length(unique(full_detail_database$drug_name)),"drugs")

# -----check for rows without ratio values
cols_to_check <- c("tissue_ratio","cell_ratio", "cell_in_tissue_ratio","max_phase","drug_name")


rows_with_invalid <- full_detail_drug_database %>%
  filter(if_any(all_of(cols_to_check), ~ is.na(.) | is.nan(.) | is.infinite(.)))


num_rows_with_invalid <- nrow(rows_with_invalid)


invalid_counts <- sapply(
  full_detail_drug_database[cols_to_check],
  function(x) sum(is.na(x) | is.nan(x) | is.infinite(x))
)
paste(num_rows_with_invalid,"rows were found with missing ratios for drug related data and following is a summary")
invalid_counts

#Converting max phase values and ordering them
full_detail_drug_database <- full_detail_drug_database %>%
  mutate(
    max_phase = factor(
      max_phase,
      levels =
        c(
          "Withdrawn",
          "PRECLINICAL",
          "IND",
          "EARLY_PHASE_1",
          "PHASE_1",
          "PHASE_1_2",
          "PHASE_2",
          "PHASE_2_3",
          "PHASE_3",
          "PREAPPROVAL",
          "APPROVAL"
        )
      ,
      ordered = TRUE
    )
  )

#*******************Prepare data for the model**********************************
#load
data_for_model<-full_detail_drug_database

#add target number column
data_for_model <- data_for_model %>%
  group_by(drug_name) %>%
  mutate(single_gene_target = n_distinct(na.omit(gene_name)) == 1) %>%
  ungroup()

#select data
data_for_model_selected<-data_for_model %>%dplyr::select(
  drug_id,
  drug_name,
  gene_name,
  tissue_ratio,
  cell_ratio,
  cell_in_tissue_ratio,
  max_phase,
  single_gene_target,
  has_mendelian_evidence,
  is_oncology_drug,
  has_cancer_biomarker_evidence,
  has_cancer_driver_evidence
)



#----- select single gene target drugs
data_for_model_single_target <- data_for_model_selected %>%
  group_by(drug_name) %>%
  filter(n_distinct(gene_name) == 1) %>%
  ungroup()

#Drop hgnc column to combine with multi target data
data_for_model_single_target<-data_for_model_single_target%>%dplyr::select(-gene_name)

# ----- multi target data
data_for_model_multi_target_data <- data_for_model_selected %>%
  group_by(drug_name,drug_id) %>%
  filter(n_distinct(gene_name) > 1) %>%
  ungroup()

# ----- mean
# ----- Prepare DF with just mean ratio for each drug---------------------------
cols_to_aggregate <- c("tissue_ratio", "cell_ratio", "cell_in_tissue_ratio")

#check for any duplicate values - just for extra safety
data_for_model_multi_target_data %>%
  group_by(drug_name,drug_id) %>%
  summarise(n_phases = n_distinct(max_phase)) %>%
  filter(n_phases > 1)

#calculate mean
# the length should be similar to the number of multi target drugs after mean calculation
data_for_model_multi_target_mean <- data_for_model_multi_target_data %>%
  group_by(drug_name) %>%
  summarise(
    across(
      all_of(cols_to_aggregate),
      ~ {
        x <- suppressWarnings(as.numeric(.x))
        x <- x[!(is.na(x) | is.nan(x) | is.infinite(x))]
        if (length(x) == 0) NA_real_ else mean(x)
      }
    ),
    
    max_phase = max(max_phase, na.rm = TRUE),
    single_gene_target = dplyr::first(single_gene_target),
    
    has_mendelian_evidence = any(has_mendelian_evidence, na.rm = TRUE),
    has_cancer_biomarker_evidence = any(has_cancer_biomarker_evidence, na.rm = TRUE),
    has_cancer_driver_evidence = any(has_cancer_driver_evidence, na.rm = TRUE),
    is_oncology_drug = any(is_oncology_drug, na.rm = TRUE),
    
    .groups = "drop"
  )


#add single target data and multi target mean data together
data_for_model_combined <- bind_rows(data_for_model_single_target,data_for_model_multi_target_mean)

#Change max phase to marketed
data_for_model_combined <- data_for_model_combined %>%mutate(Marketed = max_phase == "APPROVAL")
#data_for_model_combined <-data_for_model_combined %>%dplyr::select(-max_phase)

#save this as a tsv for the model
readr::write_tsv(data_for_model_combined,"data_for_model_combined.tsv")


#*******************GATHER NEEDED DATA TO PLOT**********************************

# ----- Prepare DF with just max ratio for each drug----------------------------


cols_numeric <- c("tissue_ratio", "cell_ratio", "cell_in_tissue_ratio","max_phase")

all_drug_data_drug_max <- full_detail_drug_database %>%
  group_by(drug_name) %>%
  summarise(
    across(
      all_of(cols_numeric),
      ~ {
        x <- .x[!(is.na(.x) | is.nan(.x) | is.infinite(.x))]
        if (length(x) == 0) NA_real_ else max(x)
      }
    ),
    max_phase = {
      x <- max_phase[!is.na(max_phase)]
      if (length(x) == 0) {
        # return NA with same factor structure
        factor(NA, levels = levels(max_phase), ordered = is.ordered(max_phase))
      } else {
        max(x)
      }
    },
    .groups = "drop"
  )



# ----- Prepare DF with just mean ratio for each drug---------------------------
cols_to_aggregate <- c("tissue_ratio", "cell_ratio", "cell_in_tissue_ratio")

all_drug_data_drug_mean <- full_detail_drug_database %>%
  group_by(drug_name) %>%
  summarise(
    
    # ✅ Numeric columns
    across(
      all_of(cols_to_aggregate),
      ~ {
        x <- suppressWarnings(as.numeric(.x))
        x <- x[!(is.na(x) | is.nan(x) | is.infinite(x))]
        if (length(x) == 0) NA_real_ else mean(x)
      }
    ),
    
    # ✅ Robust max_phase handling
    max_phase = max(max_phase),
    
    .groups = "drop"
  )

#*********************Data selection********************************************

# ----- Get drugs that only affect one gene-------------------------------------
single_target_data <- full_detail_drug_database %>%
  group_by(drug_name) %>%
  filter(n_distinct(gene_name) == 1) %>%
  ungroup()

#check for more duplicate drug rows
paste("checking for any left dupicate data")
single_target_data %>%
  group_by(drug_name) %>%
  filter(n() > 1) %>%
  arrange(drug_name)

# ----- Get Drugs targeting more genes------------------------------------------

# ----- just all data
multi_target_data <- full_detail_drug_database %>%
  group_by(drug_name) %>%
  filter(n_distinct(gene_name) > 1) %>%
  ungroup()

# ----- mean

#check for any duplicate values
multi_target_data %>%
  group_by(drug_name) %>%
  summarise(n_phases = n_distinct(max_phase)) %>%
  filter(n_phases > 1)

#calculate mean
multi_target_data_mean <- multi_target_data %>%
  group_by(drug_name) %>%
  summarise(
    across(
      all_of(cols_to_aggregate),
      ~ {
        x <- suppressWarnings(as.numeric(.x))
        x <- x[!(is.na(x) | is.nan(x) | is.infinite(x))]
        if (length(x) == 0) NA_real_ else mean(x)
      }
    ),
    
    max_phase = max(max_phase, na.rm = TRUE),
    
    .groups = "drop"
  )



# -----max value only for multi target drugs
multi_target_data_only_in_max <- all_drug_data_drug_max %>%
  anti_join(
    single_target_data,
    by = "drug_name"
  )

#-------------------Get Drugs With Mendelian evidence---------------------------
#mendelian- considering any drug with a mendelian gene mendelian
mendalian_target_data <- full_detail_drug_database %>%
  filter( has_mendelian_evidence==TRUE)
mendalian_target_drug_list<-unique(mendalian_target_data$drug_name)

mendalian_single_target_data<- single_target_data %>%
  filter(trimws(drug_name) %in% trimws(mendalian_target_drug_list))

mean_mendalian_multi_target_data<- multi_target_data_mean %>%
  filter(trimws(drug_name) %in% trimws(mendalian_target_drug_list))

#non-mendalian

non_mendalian_single_target_data <- single_target_data %>%
  filter(!trimws(drug_name) %in% trimws(mendalian_target_drug_list))

mean_non_mendalian_multi_target_data<- multi_target_data_mean %>%
  filter(!trimws(drug_name) %in% trimws(mendalian_target_drug_list))


#-------------------Get Drugs With Cancer BioMarker evidence--------------------
#cancer biomarker- considering any drug with a cancer biomarker gene

CBM_target_data <- full_detail_drug_database %>%
  filter( has_cancer_biomarker_evidence==TRUE)
CBM_target_drug_list<-unique(CBM_target_data$drug_name)

CBM_single_target_data<- single_target_data %>%
  filter(trimws(drug_name) %in% trimws(CBM_target_drug_list))

mean_CBM_multi_target_data<- multi_target_data_mean %>%
  filter(trimws(drug_name) %in% trimws(CBM_target_drug_list))

#non-CBM
non_CBM_single_target_data <- single_target_data %>%
  filter(!trimws(drug_name) %in% trimws(CBM_target_drug_list))

mean_non_CBM_multi_target_data<- multi_target_data_mean %>%
  filter(!trimws(drug_name) %in% trimws(CBM_target_drug_list))

#-------------------Get Drugs With Cancer Driver evidence--------------------
#cancer Driver- considering any drug with a cancer Driver gene

CDR_target_data <- full_detail_drug_database %>%
  filter( has_cancer_driver_evidence==TRUE)
CDR_target_drug_list<-unique(CDR_target_data$drug_name)

CDR_single_target_data<- single_target_data %>%
  filter(trimws(drug_name) %in% trimws(CDR_target_drug_list))

mean_CDR_multi_target_data<- multi_target_data_mean %>%
  filter(trimws(drug_name) %in% trimws(CDR_target_drug_list))

#non-CDR
non_CDR_single_target_data <- single_target_data %>%
  filter(!trimws(drug_name) %in% trimws(CDR_target_drug_list))

mean_non_CDR_multi_target_data<- multi_target_data_mean %>%
  filter(!trimws(drug_name) %in% trimws(CDR_target_drug_list))

cat("**********************SUMMARY*****************************")
paste("Out of",length(all_drug_data_drug_mean$drug_name),"total drugs with phase data")
paste(length(single_target_data$drug_name),"had single gene targets")
paste(length(multi_target_data_mean$drug_name),"had multiple target genes")
cat("******")
paste("From",length(single_target_data$drug_name),"drugs with single gene targets")
paste(length(mendalian_single_target_data$drug_name),"had mendelian evidence")
paste(length(non_mendalian_single_target_data$drug_name),"did  not have mendelian evidence")

paste("From",length(multi_target_data_mean$drug_name),"drugs with multiple target genes")
paste(length(mean_mendalian_multi_target_data$drug_name),"had mendelian evidence")
paste(length(mean_non_mendalian_multi_target_data$drug_name),"did not have mendelian evidence")

cat("******")
paste("From",length(single_target_data$drug_name),"drugs with single gene targets")
paste(length(CBM_single_target_data$drug_name),"had CBM evidence")
paste(length(non_CBM_single_target_data$drug_name),"did  not have CBM evidence")

paste("From",length(multi_target_data_mean$drug_name),"drugs with multiple target genes")
paste(length(mean_CBM_multi_target_data$drug_name),"had CBM evidence")
paste(length(mean_non_CBM_multi_target_data$drug_name),"did not have CBM evidence")

cat("******")
paste("From",length(single_target_data$drug_name),"drugs with single gene targets")
paste(length(CDR_single_target_data$drug_name),"had CDR evidence")
paste(length(non_CDR_single_target_data$drug_name),"did  not have CDR evidence")

paste("From",length(multi_target_data_mean$drug_name),"drugs with multiple target genes")
paste(length(mean_CDR_multi_target_data$drug_name),"had CDR evidence")
paste(length(mean_non_CDR_multi_target_data$drug_name),"did not have CDR evidence")

#*******************************************************************************
#**********************Ratio Calculations***************************************
#*******************************************************************************

#-----------------------All drug related data------------------------------------

# Apply function to calculate the threshold passes for all drugs- max ratio
for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("max_all_drugs", ratio, sep = "_")
  label_name<-"max_all_drugs"
  
  assign(
    dataframe_name,
    calculate_gene_percent(
      data = all_drug_data_drug_max,
      gene_col = drug_name,
      ratio_col = !!sym(ratio),
      thresholds = thresholds,
      phase_column = phase_column,
      label = label_name
    )
  )
}

# Apply function to calculate the threshold passes for all drugs- MEAN
for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_all_drugs", ratio, sep = "_")
  label_name<-"mean_all_drugs"
  
  assign(
    dataframe_name,
    calculate_gene_percent(
      data = all_drug_data_drug_mean,
      gene_col = drug_name,
      ratio_col = !!sym(ratio),
      thresholds = thresholds,
      phase_column = phase_column,
      label = label_name
    )
  )
}

# Run the second function with phase data for all drugs -MAX
for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("max_all_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = all_drug_data_drug_max,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# Run the second function with phase data for all drugs -MEAN
for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_all_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = all_drug_data_drug_mean,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

#---------------------Only selected data----------------------------------------


# Run the second function with phase data for 
# -----drugs targeting one gene
for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting multiple genes -max

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("max_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = multi_target_data_only_in_max,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting multiple genes -mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = multi_target_data_mean,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Mendalian genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mendalian_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mendalian_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Mendalian genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_mendalian_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_mendalian_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting non-Mendalian genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("non_mendalian_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = non_mendalian_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting non-Mendalian genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_non_mendalian_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_non_mendalian_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting CBM genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("CBM_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = CBM_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting CBM genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_CBM_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_CBM_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Non CBM genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("non_CBM_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = non_CBM_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Non CBM genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_non_CBM_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_non_CBM_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}



# -----drugs targeting CDR genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("CDR_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = CDR_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting CDR genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_CDR_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_CDR_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Non CDR genes -single target data

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("non_CDR_single_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = non_CDR_single_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

# -----drugs targeting Non CDR genes -multi target data mean

for (ratio in ratios_to_use) {
  
  dataframe_name <- paste("mean_non_CDR_multi_target_drugs", ratio, sep = "_")
  
  assign(
    dataframe_name,
    calculate_gene_percentages(
      data        = mean_non_CDR_multi_target_data,
      gene_col    = drug_name,
      phase_col   = !!rlang::sym(phase_column),
      ratio_col   = !!rlang::sym(ratio),  
      thresholds  = thresholds
    )
  )
}

#--------------------combine all these datasets---------------------------------


# Get object names matching pattern
df_names <- ls(pattern = "(target_drugs|all_drugs)")

# Keep only data frames
df_names <- df_names[sapply(df_names, function(x) is.data.frame(get(x)))]


#Remove wrong data frames listed - if there is any
not_proper_raio_dataframes<-c("missing_in_new_data_after_z_ratio")
df_names<-df_names[!df_names %in% not_proper_raio_dataframes]

# ----- Print names and count
cat(yellow("Following is the list of all available datasets\n"))
print(df_names)
cat(green("<======SELET DATASET TYPES FROM HERE*********\n"))
cat("Number of dataframes:", length(df_names), "\n")
cat(yellow("but using following selected datasets as selected by x_levels_order\n"))



#select just selected datasets above
selected_df_names <- unlist(
  lapply(selected_datasets, function(prefix) {
    df_names[startsWith(df_names, prefix)]
  })
)
selected_df_names 

df_list <- mget(df_names)

# Combine into one dataframe
combined_df <- dplyr::bind_rows(df_list, .id = "source_df")

# Separate data info and ratio info into 2

combined_df_separated <- combined_df %>%
  separate(
    col = source_df,
    into = c("Gene_type", "Ratio_type"),
    sep = "(?<=drugs)_"
  )


#---------------Add all drug lines to selected datasets-------------------------

# Add all drugs datasets to it, so it can be used in all the datasets as baseline
all_drugs_keywords<-c("all_drugs","mean_all_drugs","max_all_drugs")

all_drugs_datasets<-unlist(
  lapply(all_drugs_keywords, function(prefix) {
    df_names[startsWith(df_names, prefix)]
  })
)

selected_df_names_and_baseline<-c(selected_df_names,all_drugs_datasets)

# Load into list
df_list <- mget(selected_df_names_and_baseline)

#convert phase into character
df_list <- lapply(df_list, function(df) {
  df$max_phase <- as.character(df$max_phase)
  df
})


# Duplicating all_genes rows for all drug related genes - max value for multi targets
combined_df_separated_with_all_genes_integrated <- combined_df_separated %>%
  bind_rows(
    combined_df_separated %>%
      filter(Gene_type == "max_all_drugs") %>%
      mutate(
        Gene_type = "max_multi_target_drugs",
        # Category  = "all_genes_for_oncology"  # optional: update this too if needed
      )
  )

# Duplicating mean_all_genes rows for all the other selected datasets
for (dataset in selected_datasets) {
  combined_df_separated_with_all_genes_integrated <- combined_df_separated_with_all_genes_integrated %>%
    bind_rows(
      combined_df_separated %>%
        filter(Gene_type == "mean_all_drugs") %>%
        mutate(
          Gene_type = dataset,
          # Category  = "all_genes_for_oncology"  # optional: update this too if needed
        )
    )
}

#*************************
# pass it to the final df

final_df <- combined_df_separated_with_all_genes_integrated%>%
  mutate(
    Gene_type = dplyr::recode(
      Gene_type,
      "max_all_drugs" = "max_all_target_drugs",
      "mean_all_drugs"="mean_all_target_drugs"
    )
  )


# ------------------------- Plot data ------------------------------------------

#set max y according to the selected thresholds first
max_y<- max(final_df$percent_exceeding[final_df$threshold>=min(thresholds)])
max_y<-ceiling(max_y / 5) * 5

#override y max
if (y_max_method=="manual") {
  max_y<-y_max_limit
}

#max_y<-35

#copy final df to plot df
plotting_df<-final_df



#Order facets according to selection above
plotting_df[[plots_along_x_axis]] <- factor(
  plotting_df[[plots_along_x_axis]],
  levels = x_levels_order
)

plotting_df[[plots_along_y_axis]] <- factor(
  plotting_df[[plots_along_y_axis]],
  levels = y_levels_order
)

#Drop rows that are not used in the selected plots- Because they are given NA in previous step
plotting_df <- plotting_df %>%
  filter(!is.na(Gene_type))


gene_plot<-ggplot(plotting_df,
                  aes(x = threshold,
                      y = percent_exceeding,
                      color = !!rlang::sym(phase_column),
                      group = !!rlang::sym(phase_column))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0,max_y)) +
  scale_x_continuous(limits = c(min(thresholds),max(thresholds))) +
  labs(
    title = paste(drug_plot_title,threshold_name,sep = "-"),
    x = x_axis_name,
    y = y_axis_name,
    color = "Max Clinical Phase"
  ) +
  theme_minimal(base_size = 13)+
  scale_color_manual(
    values = c(
      #      "All phases" = "lightblue",
      "max_all_drugs" = "purple",
      "mean_all_drugs" = "black",
      "1" = "pink",
      "2" = "orange",
      "3" = "gold",
      "4" = "forestgreen",
      "Withdrawn"="grey",
      "-1"="white", #genes without proper data are invisible
      #The following are only for old phase data
      "Target of marketed drug"="forestgreen",
      "Target of phase 1 drug"="pink",
      "Target of phase 2 drug"="orange",
      "Target of phase 3 drug"="gold",
      "Target of withdrawn drug"="grey",
      "all genes from paper"="black",
      "paper_all_gene_drugs"="black",
      "Marketed"="forestgreen",
      "Other"="orange",
      "PRECLINICAL"="red",
      "IND"="red",
      "EARLY_PHASE_1"="pink",
      "PHASE_1"="pink",
      "PHASE_1_2"="pink",
      "PHASE_2"="orange",
      "PHASE_2_3"="orange",
      "PHASE_3"="gold",
      "PREAPPROVAL"="gold",
      "APPROVAL"="forestgreen"
    )
  )+
  facet_grid(get(plots_along_y_axis)~get(plots_along_x_axis))

# add number of genes in each point

gene_plot <- gene_plot +
  geom_text_repel(aes(label = n_exceeding), size = 3)



gene_plot

#************************ saving data*******************************************

#create subfolder
gene_subfolder_name <- "drugwise_comparison_plots"

if (!dir.exists(gene_subfolder_name)) {
  dir.create(gene_subfolder_name)
}

#save plot
plot_saving_name <- paste0(paste(na.omit(unique(final_df$Ratio_type)), collapse = "_"),".png")
plot_saving_name<-paste(gene_subfolder_name,"/max_val_",plot_saving_name,sep = "")

#calculating width and height based on the number of plots to keep it consistant and comparable
plot_height<-length(unique(plotting_df[[plots_along_y_axis]])) * one_plot_height
plot_width<-length(unique(plotting_df[[plots_along_x_axis]])) * one_plot_width


#save the plot
ggsave(plot_saving_name,gene_plot, width = plot_width, height = plot_height, dpi = 150)

# save with custom name and size
custom_save_name<-"drug_success_significant_looking_comparisons_zoomed_1_5"
custom_save_name<-paste("./",gene_subfolder_name,"/",custom_save_name,".png",sep = "")
ggsave(custom_save_name,gene_plot, width = plot_width, height = plot_height, dpi = 150,limitsize = FALSE)
