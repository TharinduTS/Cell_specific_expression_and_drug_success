library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

# ============================
# Step 1: Build drug-level data 
# ============================


# install.packages(c("tidyverse", "janitor","dplyr","crayon","writexl"))
library(tidyverse)
library(janitor)
library(dplyr)
library(crayon)
library(writexl)
library(stringr)
library(tidyr)
library(readxl)


# ---- 1) Load & clean ----
# I downloaded this drug and cell type enrichment dataset from chapter 8 of my drug target project

raw_data <- readr::read_tsv("./data_for_stats.tsv", guess_max = 1e5) %>%
  clean_names()

# I am skipping first 2 rows with no data here
tissue_paper_data_drug_related_genes<-read_excel("./t1_genes_drug_targets_z_scores_from_tissue_paper.xls",sheet = "drug targets",skip = 2)

# And read data for summary of all genes
# This is to read values for **Cell + Tissue** based enrichment values *************************************
cell_and_tissue<-readr::read_tsv("./gene_top2_enrichment.tsv", guess_max = 1e5) %>%
  clean_names()

# Read enrichment data based **only on cell type** data
cell_type_only<-readr::read_tsv("./cell_type_only_top2_enrichment.tsv", guess_max = 1e5) %>%
  clean_names()

# And read data from the paper Ryaboshapkina et al table 1
 
tissue_paper_data_all_genes<-read_excel("./t1_genes_drug_targets_z_scores_from_tissue_paper.xls",sheet = "18,377 protein-coding")

# # Prepare and clean dataframe *******************************************
# 
# #rename the messed up column name
names(raw_data)[names(raw_data) == "parent_molecule_ch_embl_id"] <- "parent_molecule_chembl_id"
names(raw_data)[names(raw_data) == "target_ch_embl_id"] <- "target_chembl_id"
names(raw_data)[names(raw_data) == "ch_embl_hgnc"] <- "chembl_hgnc"

#And rename the log2 enrichment value which is labeled as enrichment value to avoid confusions
names(raw_data)[names(raw_data) == "enrichment_value"] <- "Penalized_log2_enrichment_value"
#Drop all rows with unknown drug phase (-1)
raw_data <- raw_data[raw_data$max_phase != -1, ]

# checking for any rows with NA values
incomplete_cases_df <- raw_data[!complete.cases(raw_data), ]

# Getting a summary of data available after filtering 

no_of_unique_parent_molecule_IDs <- dplyr::n_distinct(raw_data$parent_molecule_chembl_id, na.rm = TRUE)
no_of_unique_parent_molecule_names <- dplyr::n_distinct(raw_data$parent_molecule_name, na.rm = TRUE)

paste(no_of_unique_parent_molecule_IDs, "unique parent molecule IDs and", sep = " ", no_of_unique_parent_molecule_names, "unique parent molecule names found")
if (no_of_unique_parent_molecule_IDs == no_of_unique_parent_molecule_names) {
  print("All the unique parent molecule IDs have unique names")
} else {
  cat(red("NOT ALL PARENT MOLECULES HAVE UNIQUE NAMES"), "\n")
}

no_of_unique_target_IDs<-dplyr::n_distinct(raw_data$target_chembl_id, na.rm = TRUE)
no_of_unique_target_names<-dplyr::n_distinct(raw_data$target_name, na.rm = TRUE)
paste(no_of_unique_target_IDs,"unique target IDs and", sep = " ",no_of_unique_target_names, "unique target names found")

if (no_of_unique_target_IDs == no_of_unique_target_names) {
  print("All the unique target IDs have unique names")
} else {
  cat(red("NOT ALL TARGETS HAVE UNIQUE NAMES"), "\n")
}





# Find target names associated with multiple target IDs
targets_with_multiple_IDs <- raw_data %>%
  filter(!is.na(target_chembl_id), !is.na(target_name)) %>%
  distinct(target_chembl_id, target_name) %>%
  group_by(target_name) %>%
  filter(n() > 1) %>%
  arrange(target_name) %>%
  ungroup()

print("Targets with the same name but multiple IDs are -because of different drug mechanisms")
targets_with_multiple_IDs

print("There can be more than one drug mechanisms for some drugs.")


# Find parent molecule names associated with multiple target IDs
parent_molecules_with_multiple_target_IDs <- raw_data %>%
  filter(!is.na(target_chembl_id), !is.na(parent_molecule_name),!is.na(target_name)) %>%
  distinct(target_chembl_id,parent_molecule_name) %>%
  group_by(parent_molecule_name) %>%
  filter(n() > 1) %>%
  arrange(parent_molecule_name) %>%
  ungroup()



# for the downstream analyses

cat(red("Using ",length(unique(raw_data$parent_molecule_name)),"parent molecules with unique names and",length(unique(raw_data$target_chembl_id)),"targets with unique IDs for downstream analyses after dropping target names column"), "\n")





# Check for parent IDs associated with multiple drug success phases - to make this error proof
parent_molecule_IDs_with_multiple_success_phases <- raw_data %>%
  filter(!is.na(max_phase), !is.na(parent_molecule_chembl_id)) %>%
  distinct(max_phase,parent_molecule_chembl_id) %>%
  group_by(parent_molecule_chembl_id) %>%
  filter(n() > 1) %>%
  arrange(parent_molecule_chembl_id) %>%
  ungroup()

parent_molecules_with_errors<-length(parent_molecule_IDs_with_multiple_success_phases)

paste("Found ",parent_molecules_with_errors," drugs with more than one max phase. I found this to be an error in CHeMBL dataset",sep = "")

print("Dropping these from the dataframe")

cleaned_raw_data <- raw_data %>%
  filter(
    !parent_molecule_chembl_id %in% 
      parent_molecule_IDs_with_multiple_success_phases$parent_molecule_chembl_id
  )

paste(length(unique(cleaned_raw_data$parent_molecule_chembl_id))," Unique parent molecules and ", length(unique(cleaned_raw_data$target_name)),"targets(names) are remaining for the analyses after cleaning")


#****************Drug based analysis-My Method******************************************************
#****************Drug based analysis-My Method******************************************************
#****************Drug based analysis-My Method******************************************************


targets_with_multiple_success_phases <- cleaned_raw_data %>%
  group_by(target_name) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  select(target_name)



# Find target names associated with multiple drug success phases
targets_with_multiple_success_phases <- cleaned_raw_data %>%
  filter(!is.na(max_phase), !is.na(target_name)) %>%
  distinct(max_phase,target_name) %>%
  group_by(target_name) %>%
  filter(n() > 1) %>%
  arrange(target_name) %>%
  ungroup()

paste("Out of those ",length(unique(cleaned_raw_data$target_name))," unique drug targets, ",length(unique(targets_with_multiple_success_phases$target_name)),"drug targets have variable drug success for different drugs")

#separate targets with just one max phase
targets_with_one_max_phase <- cleaned_raw_data %>%
  filter(
    !target_name %in% 
      targets_with_multiple_success_phases$target_name
  )

#separate targets with multiple max phases

targets_with_different_max_phases <- cleaned_raw_data %>%
  filter(
    target_name %in% 
      targets_with_multiple_success_phases$target_name
  )

#print summary

paste("All the cleaned data with ",length(unique(cleaned_raw_data$target_name))," targets were saved as cleaned_raw_data dataframe", sep = '')
paste("All the targets with only one max phase with",length(unique(targets_with_one_max_phase$target_name)),"targets were saved as targets_with_one_max_phase dataframe")
paste("All the targets with multiple max phases with",length(unique(targets_with_different_max_phases$target_name)),"targets were saved as targets_with_different_max_phases dataframe")


# Find drugs that may affect multiple genes

drugs_affecting_multiple_genes <- cleaned_raw_data %>%
  filter(!is.na(gene_name), !is.na(parent_molecule_chembl_id)) %>%
  distinct(gene_name,parent_molecule_chembl_id) %>%
  group_by(parent_molecule_chembl_id) %>%
  filter(n() > 1) %>%
  arrange(parent_molecule_chembl_id) %>%
  ungroup()

# ============================
# Averaging expression for drugs affecting multiple genes 
# ============================

print("This has too many data layers. First I am going to get rid of")
print("DIFFERENT GENES WITHIN CELL TYPES BY AVERAGING THE AFFECT OF EACH GENE WITHIN EACH CELL TYPE for drugs affecting nultipple genes")
print("To be used later")

#*******************************************************************************************************
# Function for summarizing various stats

# =============================================================================
# summarize_numeric_by_group():
#   - data:        data frame (e.g., raw_data)
#   - group_cols:  character vector of columns to group by (e.g., c("parent_molecule_name"))
#   - value_col:   single character string of the numeric column to summarize
#                  (e.g., "enrichment_value" OR "log2_enrichment_value")
#   - stats:       which stats to return. Any of:
#                    c("n", "n_non_missing", "mean", "median", "sd", "iqr",
#                      "min", "p10", "p25", "p75", "p90", "max", "topkmean")
#   - top_k:       if "topkmean" requested, take mean of top_k values
#   - na_rm:       remove NAs during numeric summaries
#   - keep_cols:   character vector of metadata columns to carry along;
#                  picks first non-NA and can warn if inconsistent
#   - check_consistency: warn if keep_cols have >1 distinct non-NA value within a group
#   - arrange_by:  optional name of summary column to sort by (desc)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
})


summarize_numeric_by_group <- function(data,
                                       group_cols,
                                       value_col,
                                       stats = c("max", "mean", "p90"),
                                       top_k = 3,
                                       na_rm = TRUE,
                                       keep_cols = NULL,
                                       check_consistency = TRUE,
                                       arrange_by = NULL) {
  
  # check columns exist
  missing_cols <- setdiff(c(group_cols, value_col, keep_cols), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # convert to numeric if needed
  data[[value_col]] <- suppressWarnings(as.numeric(data[[value_col]]))
  
  # consistency helper
  pick_one <- function(vec, nm, key) {
    vals <- unique(vec[!is.na(vec)])
    if (length(vals) > 1 && check_consistency) {
      warning(sprintf(
        "Column '%s' has multiple distinct values within %s: %s",
        nm, key, paste(vals, collapse = ", ")
      ))
    }
    vals[1]
  }
  
  # group and summarize main numeric stats
  out <- data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n = n(),
      n_non_missing = sum(!is.na(.data[[value_col]])),
      mean = if ("mean" %in% stats) mean(.data[[value_col]], na.rm = na_rm) else NULL,
      median = if ("median" %in% stats) median(.data[[value_col]], na.rm = na_rm) else NULL,
      sd = if ("sd" %in% stats) sd(.data[[value_col]], na.rm = na_rm) else NULL,
      iqr = if ("iqr" %in% stats) IQR(.data[[value_col]], na.rm = na_rm) else NULL,
      min = if ("min" %in% stats) min(.data[[value_col]], na.rm = na_rm) else NULL,
      p10 = if ("p10" %in% stats) quantile(.data[[value_col]], 0.1, na.rm = na_rm) else NULL,
      p25 = if ("p25" %in% stats) quantile(.data[[value_col]], 0.25, na.rm = na_rm) else NULL,
      p75 = if ("p75" %in% stats) quantile(.data[[value_col]], 0.75, na.rm = na_rm) else NULL,
      p90 = if ("p90" %in% stats) quantile(.data[[value_col]], 0.9, na.rm = na_rm) else NULL,
      max = if ("max" %in% stats) max(.data[[value_col]], na.rm = na_rm) else NULL,
      topkmean = if ("topkmean" %in% stats) {
        vals <- sort(.data[[value_col]], decreasing = TRUE)
        mean(head(vals, top_k), na.rm = na_rm)
      } else NULL,
      .groups = "drop"
    )
  
  # add keep_cols
  if (!is.null(keep_cols)) {
    meta <- data %>%
      group_by(across(all_of(group_cols))) %>%
      summarise(
        across(all_of(keep_cols),
               ~ pick_one(.x,
                          cur_column(),
                          paste0(paste(group_cols, collapse = ", "), " = ", paste(cur_group(), collapse = " / "))),
               .names = "{.col}"),
        .groups = "drop"
      )
    
    out <- left_join(out, meta, by = group_cols)
  }
  
  if (!is.null(arrange_by) && arrange_by %in% names(out)) {
    out <- out %>% arrange(desc(.data[[arrange_by]]))
  }
  
  out
}


#*****************************************************************************************************

#******************************************************************************************


#selecting drugs affecting more than one gene
drugs_affecting_more_than_one_gene <- cleaned_raw_data %>%
  filter(
    parent_molecule_chembl_id %in% 
      drugs_affecting_multiple_genes$parent_molecule_chembl_id
  )

#*************************************************************************

# Averaging expression for drugs affecting multiple genes

drugs_affecting_more_than_one_gene_averaged <- summarize_numeric_by_group(
  data        = drugs_affecting_more_than_one_gene,
  group_cols  = c("parent_molecule_name","cell_type"),
  value_col   = "Penalized_log2_enrichment_value",
  stats       = c("mean"),
  keep_cols   = c("parent_molecule_chembl_id","parent_molecule_type", "max_phase", "first_approval","warning_type","first_withdrawn_year","efo_id"),
  arrange_by  = "max"  # sort by descending max
)

# renaming this mean enrichment as enrichment to make downstream analysis smooth
drugs_affecting_more_than_one_gene_averaged<-rename(drugs_affecting_more_than_one_gene_averaged,Penalized_log2_enrichment_value=mean)




# Averaging expression for drugs affecting multiple genes - separately for different target IDs/mechanisms

drugs_affecting_more_than_one_gene_averaged_sep_target_IDs <- summarize_numeric_by_group(
  data        = drugs_affecting_more_than_one_gene,
  group_cols  = c("parent_molecule_name","cell_type","target_chembl_id"),
  value_col   = "Penalized_log2_enrichment_value",
  stats       = c("mean"),
  keep_cols   = c("parent_molecule_chembl_id","parent_molecule_type", "max_phase", "first_approval","warning_type","first_withdrawn_year","efo_id"),
  arrange_by  = "max"  # sort by descending max
)

# renaming this mean enrichment as enrichment to make downstream analysis smooth
drugs_affecting_more_than_one_gene_averaged_sep_target_IDs<-rename(drugs_affecting_more_than_one_gene_averaged_sep_target_IDs,Penalized_log2_enrichment_value=mean)




#selecting drugs affecting only one gene

drugs_affecting_just_one_gene <- cleaned_raw_data %>%
  filter(
    !parent_molecule_chembl_id %in% 
      drugs_affecting_multiple_genes$parent_molecule_chembl_id
  )


#*************************************************************************
#preparing cell+tissue data for the analysis - drugs affecting just one gene


drugs_affecting_just_one_gene_cell_and_tissue <- drugs_affecting_just_one_gene %>%
  
  # 1. Decode HTML entity just in case
  mutate(present_tissues = str_replace_all(present_tissues, "&amp;", "&")) %>%
  
  # 2. Split into one row per tissue
  separate_rows(present_tissues, sep = "\\s*&\\s*") %>%
  
  # 3. Extract enrichment value
  mutate(
    enrichment_value = as.numeric(
      str_extract(present_tissues, "^[0-9.]+")
    ),
    
    # 4. Extract tissue name (between ':' and ' with')
    tissue = str_trim(
      str_extract(present_tissues, "(?<=:).*?(?=\\swith)")
    ),
    
    # 5. Create cell_type + tissue combination
    cell_type_tissue = paste(cell_type, tissue, sep = " - ")
  ) %>%
  
  # 6. Keep only relevant columns (adjust as needed)
  select(
#    cell_type,
#    tissue,
#    cell_type_tissue,
#    enrichment_value,
    everything()
  )

#preparing cell+tissue data for the analysis - drugs affecting more than one gene


drugs_affecting_more_than_one_gene_cell_and_tissue <- drugs_affecting_more_than_one_gene %>%
  
  # 1. Decode HTML entity just in case
  mutate(present_tissues = str_replace_all(present_tissues, "&amp;", "&")) %>%
  
  # 2. Split into one row per tissue
  separate_rows(present_tissues, sep = "\\s*&\\s*") %>%
  
  # 3. Extract enrichment value
  mutate(
    enrichment_value = as.numeric(
      str_extract(present_tissues, "^[0-9.]+")
    ),
    
    # 4. Extract tissue name (between ':' and ' with')
    tissue = str_trim(
      str_extract(present_tissues, "(?<=:).*?(?=\\swith)")
    ),
    
    # 5. Create cell_type + tissue combination
    cell_type_tissue = paste(cell_type, tissue, sep = " - ")
  ) %>%
  
  # 6. Keep only relevant columns (adjust as needed)
  select(
    #    cell_type,
    #    tissue,
    #    cell_type_tissue,
    #    enrichment_value,
    everything()
  )
#******************************************************************************************

#*************************************************************************
#preparing cell+tissue data for the analysis - drugs affecting more than one gene *averaged*


drugs_affecting_more_than_one_gene_averaged_cell_and_tissue <-summarize_numeric_by_group(
  data        = drugs_affecting_more_than_one_gene_cell_and_tissue,
  group_cols  = c("parent_molecule_name","cell_type_tissue"),
  value_col   = "Penalized_log2_enrichment_value",
  stats       = c("mean"),
  keep_cols   = c("parent_molecule_chembl_id","parent_molecule_type", "max_phase", "first_approval","warning_type","first_withdrawn_year","efo_id"),
  arrange_by  = "max"  # sort by descending max
)

# renaming this mean enrichment as enrichment to make downstream analysis smooth
drugs_affecting_more_than_one_gene_averaged_cell_and_tissue<-rename(drugs_affecting_more_than_one_gene_averaged_cell_and_tissue,Penalized_log2_enrichment_value=mean)


#******************************************************************************************



print("**************************************************")
paste("All the cleaned data with ",length(unique(cleaned_raw_data$parent_molecule_chembl_id))," drugs(parent molecules) were saved as cleaned_raw_data dataframe", sep = '')


paste("All the cleaned data with ",length(unique(cleaned_raw_data$target_name))," targets were saved as cleaned_raw_data dataframe", sep = '')
paste("All the targets with only one max phase with",length(unique(targets_with_one_max_phase$target_name)),"targets were saved as targets_with_one_max_phase dataframe")
paste("All the targets with multiple max phases with",length(unique(targets_with_different_max_phases$target_name)),"targets were saved as targets_with_different_max_phases dataframe")

paste("out of",length(unique(cleaned_raw_data$parent_molecule_chembl_id)),"total drugs (parent molecule ID)")
paste(length(unique(drugs_affecting_just_one_gene$parent_molecule_chembl_id)), "drugs only affect one gene and were saved in drugs_affecting_just_one_gene" )
paste(length(unique(drugs_affecting_more_than_one_gene$parent_molecule_chembl_id)),"drugs affect more than one gene and were saved in drugs_affecting_more_than_one_gene")


cat(red("TO REDUCE THE COMPLEXITY OF ALL THE DATA LAYERS,\n"))
cat(red("FIRST I AM LOOKING AT DRUGS THAT ONLY AFFECT ONE GENE,\n"))
paste(length(unique(drugs_affecting_just_one_gene$parent_molecule_chembl_id)), "drugs only affect one gene and were saved in drugs_affecting_just_one_gene" )
cat(blue("I am going to calculate the top VS runners up for these first ,\n"))

#*****************************************************************************************************"

# CHAPTER 1 - Drugs affecting only one gene/multiple genes"

#*****************************************************************************************************"

#********************VALUES TO CHANGE*************************************************************************


# select the dataframe here  <<<========
dataframe_name_to_use<-"cleaned_raw_data"

#options 
#For cell type data       - drugs_affecting_just_one_gene   drugs_affecting_more_than_one_gene   cleaned_raw_data <-this is both one gene and multiple genes
#                           drugs_affecting_more_than_one_gene_averaged    drugs_affecting_more_than_one_gene_averaged_sep_target_IDs

#For cell and tissue data - drugs_affecting_just_one_gene_cell_and_tissue   drugs_affecting_more_than_one_gene_cell_and_tissue
#                           drugs_affecting_more_than_one_gene_averaged_cell_and_tissue


#select the dictinct group to consider here
distinct_group<-"cell_type"

#options 
#For cell type data       - "cell_type"  "cell_type_group" "cell_type_class" to be distinct when you select top and runners up values

#For cell and tissue data - "cell_type_tissue"

#select oncology or non oncology genes

selected_gene_type<-"non_oncology"

#options  "oncology" "non_oncology" "both_onco_and_non_onco"

# Set checking threshold values for specificity
thresholds <- seq(1.5,2.5,0.25)

#********************VALUES TO CHANGE ENDS*************************************************************************
#*****************************************************************************************************


dataframe_of_the_selected_dataset<-get(dataframe_name_to_use)

#***************************************************************************
# Trying separating oncology/non oncology genes to see its affect

if (selected_gene_type == "non_oncology") {
  non_onco_genes <- tissue_paper_data_drug_related_genes$gene[
    tissue_paper_data_drug_related_genes$atLeast1DrugForOncologyIndication == 0
  ]
  
  dataframe_of_the_selected_dataset <- dataframe_of_the_selected_dataset[
    dataframe_of_the_selected_dataset$gene_name %in% non_onco_genes, 
  ]
}

if (selected_gene_type == "oncology") {
  non_onco_genes <- tissue_paper_data_drug_related_genes$gene[
    tissue_paper_data_drug_related_genes$atLeast1DrugForOncologyIndication == 1
  ]
  
  dataframe_of_the_selected_dataset <- dataframe_of_the_selected_dataset[
    dataframe_of_the_selected_dataset$gene_name %in% non_onco_genes, 
  ]
}


#***************************************************************************



#using dataframe name for the files
name_used_for_current_analysis<-dataframe_name_to_use

#This is the first/runners up ratio value to be used for the drugs with just one affected cell type.
#This should be higher than or equal to the highest ratio between first/runners up ratio for genes with multiple cell types
#Therefore, I will set up a highest value this ratio can have later in the script using this value
#For now, I am using the max ratio for genes with at least 2 cell types to a lower than the value for genes with just one cell type
#And I keep the x limit to a couple ticks lower than max ratio to keep the proper shape of the plot

custom_value <-10 #max(top_vs_runners_cell_types$Penalized_log2_enrichment_value)   


#Axis and plot names for the plot -optional
plot_title_name<-paste("Cell type sp VS Max phase reached for\n ",name_used_for_current_analysis," by ",as.character(distinct_group),"\n ",selected_gene_type,sep = "")
x_axis_name<-"Top VS Runners Up Ratio"
y_axis_name<-"N cell specific/N total %"


#*************************************************************************************************************



library(dplyr)

##collecting only the two *cell types/groups/ or classes*  with highest enrichment for each drug
top_vs_runners_cell_type_groups<- dataframe_of_the_selected_dataset %>%
  group_by(parent_molecule_name) %>%
  arrange(desc(Penalized_log2_enrichment_value), .by_group = TRUE) %>%
  distinct(get(distinct_group), .keep_all = TRUE) %>%  # keep best cell type only
  slice_head(n = 2) %>%                              # take top 2 cell types
  ungroup()


#collecting only the two *cell types*  with highest enrichment for each drug
#top_vs_runners_cell_types<-dataframe_of_the_selected_dataset %>%
#  group_by(parent_molecule_name) %>%
#  slice_max(Penalized_log2_enrichment_value, n = 2, with_ties = FALSE) %>%
#  ungroup()

# Calculating highest vs runners up ratio

df<-top_vs_runners_cell_type_groups

# conditionally changing where to find proper enrichment value (in enrichment_value if we use both cell and tissue)
#because when I splitted rows into tissues, enrichment_value is the column name I used
value_cols <- c("Penalized_log2_enrichment_value")
if (distinct_group=="cell_type_tissue") {
  value_cols <- c("enrichment_value")
}



df_ratio <- df %>%
  group_by(parent_molecule_name) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~ {
        if (cur_column() %in% value_cols) {
          vals <- sort(.x, decreasing = TRUE)
          
          if (length(vals) >= 2) {
            vals[1] / vals[2]
          } else {
            custom_value
          }
        } else {
          first(.x)
        }
      }
    ),
    .groups = "drop"
  )


#Rename the new ratio column approprately
df_ratio<-rename(df_ratio,Top_VS_runners_up_ratio=Penalized_log2_enrichment_value)

#Use a cutoff value to remove super high ratios
cat(red("Using the custom value to cutoff very large ratios,\n"))
custom_value-1


max_val <- custom_value-1
df_ratio <- df_ratio %>%
  mutate(Top_VS_runners_up_ratio = if_else(Top_VS_runners_up_ratio > max_val, max_val, Top_VS_runners_up_ratio))

cat(red("All the other values in the rows of df_ratio other than the max_vs_runners_up_ratio are for the highest enrichment value,\n"))

# X axis range customization
# higher_x_limit<-max_val-0.5
# intervals_for_plot<-seq(from=1.3,to=higher_x_limit,by=0.25)
# thresholds <- intervals_for_plot

# Outcome as ordered factor (for later steps)
df_ratio <- df_ratio %>%
  mutate(
    max_phase_ord = factor(max_phase, levels = sort(unique(na.omit(max_phase))), ordered = TRUE),
    success = as.integer(max_phase >= 4),  # for target-level summaries later
    success_ord = factor(success, levels = sort(unique(na.omit(success))), ordered = TRUE)
  )

#separate the rows without withdrawn drugs

df_without_withdrawn<-df_ratio %>%
  filter(
    warning_type != "Withdrawn",
    !is.na(warning_type),
  )

#calculating ratios for drugs that were **NOT WITHDRAWN**
percentage_df <- df_without_withdrawn %>%
  # keep only columns we need
  select(parent_molecule_name, max_phase_ord, Top_VS_runners_up_ratio) %>%
  
  # ensure one row per drug
  distinct(parent_molecule_name, .keep_all = TRUE) %>%
  
  # expand thresholds
  crossing(threshold = thresholds) %>%
  
  # compute exceedance
  mutate(exceeds = Top_VS_runners_up_ratio > threshold) %>%
  
  # aggregate
  group_by(max_phase_ord, threshold) %>%
  summarise(
    n_drugs = n(),                                # total drugs in this phase
    n_exceeding = sum(exceeds, na.rm = TRUE),     # NEW: count exceeding threshold
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  )


# with overall drug numbers

overall_df <- df_ratio %>%
  select(parent_molecule_name, Top_VS_runners_up_ratio) %>%
  distinct(parent_molecule_name, .keep_all = TRUE) %>%
  crossing(threshold = thresholds) %>%
  mutate(exceeds = Top_VS_runners_up_ratio > threshold) %>%
  group_by(threshold) %>%
  summarise(
    n_drugs = n(),
    n_exceeding = sum(exceeds, na.rm = TRUE),
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(max_phase_ord = "All phases")

# Add data for overall drug numers to drugs that were not withdrawn
plot_df <- bind_rows(
  percentage_df %>%
    mutate(max_phase_ord = as.character(max_phase_ord)),
  overall_df
)

# with all genes

# all_genes_ratio <- all_genes %>%
#  select(gene_name,ratio) %>%
#   distinct(gene_name, .keep_all = TRUE) %>%
#   crossing(threshold = thresholds) %>%
#   mutate(exceeds = ratio > threshold) %>%
#   group_by(threshold) %>%
#   summarise(
#     n_drugs = n(),
#     n_exceeding = sum(exceeds, na.rm = TRUE),
#     percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(max_phase_ord = "All genes")
# 
# plot_df <- bind_rows(
#   percentage_df %>%
#     mutate(max_phase_ord = as.character(max_phase_ord)),
#   all_genes_ratio
# )

# ****************

# Add Withdrawn numbers

# Selecting all the drugs that reached phase 4 
df_withdrawn<-df_ratio %>%
  filter(
    warning_type == "Withdrawn",
    !is.na(warning_type),
  )

withdrawn__ratio_df <- df_withdrawn %>%
  select(parent_molecule_name, Top_VS_runners_up_ratio) %>%
  distinct(parent_molecule_name, .keep_all = TRUE) %>%
  crossing(threshold = thresholds) %>%
  mutate(exceeds = Top_VS_runners_up_ratio > threshold) %>%
  group_by(threshold) %>%
  summarise(
    n_drugs = n(),
    n_exceeding = sum(exceeds, na.rm = TRUE),
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(max_phase_ord = "Withdrawn")


plot_df <- bind_rows(
  plot_df %>%
    mutate(max_phase_ord = as.character(max_phase_ord)),
  withdrawn__ratio_df
)

max_y<-max(plot_df$percent_exceeding[plot_df$threshold>=min(thresholds)])
#max_y<-max(plot_df$percent_exceeding)
#add a little extra room
round_up_to_5 <- function(x) {
  ceiling(x / 5) * 5
}
max_y<-round_up_to_5(max_y)


current_plot<-ggplot(plot_df,
                                          aes(x = threshold,
                                              y = percent_exceeding,
                                              color = max_phase_ord,
                                              group = max_phase_ord)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, max_y)) +
  #optional x limits
  scale_x_continuous(limits = c(min(thresholds),max(thresholds))) +
  labs(
    title = plot_title_name,
    x = x_axis_name,
    y = y_axis_name,
    color = "Max Clinical Phase"
  ) +
  theme_minimal(base_size = 13)+
  scale_color_manual(
    values = c(
      "All phases" = "black",
      "All genes" = "lightblue",
      "1" = "pink",
      "2" = "orange",
      "3" = "gold",
      "4" = "forestgreen",
      "Withdrawn"="grey"
    )
  )

current_plot



#create a folder for the current combination if it doesn't already exist

folder_name <- paste(name_used_for_current_analysis)

if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

plot_name<-paste(folder_name,"/",name_used_for_current_analysis,"_by_",as.character(distinct_group),"_",selected_gene_type,"_top_vs_runners.png",sep = "")

#save the plot
ggsave(file.path(plot_name),current_plot, width = 10, height = 8, dpi = 150)



#Summarize these data in a nice way
df_final_summary <- top_vs_runners_cell_type_groups %>%
  arrange(parent_molecule_chembl_id,
          desc(Penalized_log2_enrichment_value)) %>%
  group_by(parent_molecule_chembl_id) %>%
  slice_head(n = 2) %>%                 # works for 1 or 2 rows
  mutate(rank = row_number()) %>%       # 1 = highest, 2 = runner-up
  ungroup() %>%
  pivot_wider(
    #unique columns selected
    id_cols = c(
      parent_molecule_chembl_id,
      parent_molecule_name,
      parent_molecule_type,
      max_phase,
      warning_type,
     target_chembl_id,
     target_name,
     action_type,
     chembl_hgnc,
    ),
    #variable columns selected
    names_from  = rank,
    values_from = c(
      Penalized_log2_enrichment_value,
      all_of(distinct_group)
    ),
    names_glue = "{.value}_{ifelse(rank == 1, 'highest', 'runner_up')}"
  )
#add ratio column
df_final_summary$ratio=df_final_summary$Penalized_log2_enrichment_value_highest/df_final_summary$Penalized_log2_enrichment_value_runner_up



#save rough the summary
summary_name<-paste(folder_name,"/",name_used_for_current_analysis,"_by_",as.character(distinct_group),"_top_vs_runners.xlsx",sep = "")

write_xlsx(plot_df,summary_name)

#detailed summary
detailed_summary_name<-paste(folder_name,"/",name_used_for_current_analysis,"_by_",as.character(distinct_group),"_detailed_summary.xlsx",sep = "")

write_xlsx(df_final_summary,detailed_summary_name)



#*****************************************************************************************************

# CHAPTER 2 -  SIDE EFFECTS - NOTE THAT WHATEVER DATA SELECTED IN THE SECTION ABOVE IS USED HERE AS WELL

#*****************************************************************************************************  

# First I have to select only the drugs that reached drug testing phase 4 as 
# other drugs do not have warnings
#Then I should select drugs that were approved before 2010 or so

# I am starting with th df_ratio from the previous step
# MAKE SURE YOU HAVE RUN THE CORRECT df_ratio BUILIDING STEP BEFORE THIS




#**********************WARNINGS*******************************************************************

# Selecting all the drugs that reached phase 4 
df_ratio_phase_4_all<-df_ratio %>%
  filter(
    max_phase_ord >= 4,
    !is.na(warning_type),
  )


#thresholds <- c(1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.254,5)

paste("Using following amount of drugs for the analysis")

paste(sum(df_ratio_phase_4_all$warning_type == "No Warnings found", na.rm = TRUE)," Drugs with no warnings found",sep = "")
paste(sum(df_ratio_phase_4_all$warning_type == "Black Box Warning", na.rm = TRUE)," Drugs with Black Box Warnings found",sep = "")
paste(sum(df_ratio_phase_4_all$warning_type == "Withdrawn", na.rm = TRUE)," Drugs were Withdrawn",sep = "")





percentage_df <- df_ratio_phase_4_all %>%
  # keep only columns we need
  select(parent_molecule_name,warning_type,Top_VS_runners_up_ratio) %>%
  
  # ensure one row per drug
  distinct(parent_molecule_name, .keep_all = TRUE) %>%
  
  # expand thresholds
  crossing(threshold = thresholds) %>%
  
  # compute exceedance
  mutate(exceeds = Top_VS_runners_up_ratio > threshold) %>%
  
  # aggregate
  group_by(warning_type, threshold) %>%
  summarise(
    n_drugs = n(),                                # total drugs in this phase
    n_exceeding = sum(exceeds, na.rm = TRUE),     # NEW: count exceeding threshold
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  )


# with overall drug numbers

overall_df <- df_ratio_phase_4_all %>%
  select(parent_molecule_name, Top_VS_runners_up_ratio) %>%
  distinct(parent_molecule_name, .keep_all = TRUE) %>%
  crossing(threshold = thresholds) %>%
  mutate(exceeds = Top_VS_runners_up_ratio > threshold) %>%
  group_by(threshold) %>%
  summarise(
    n_drugs = n(),
    n_exceeding = sum(exceeds, na.rm = TRUE),
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(warning_type = "All drugs")


plot_df <- bind_rows(
  percentage_df %>%
    mutate(max_phase_ord = as.character(warning_type)),
  overall_df
)

max_y<-max(plot_df$percent_exceeding)
#add a little extra room
round_up_to_5 <- function(x) {
  ceiling(x / 5) * 5
}
max_y<-round_up_to_5(max_y)

warnings_plot<-ggplot(plot_df,
                     aes(x = threshold,
                         y = percent_exceeding,
                         color = warning_type,
                         group = warning_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, max_y)) +
  labs(
    x = x_axis_name,
    y = y_axis_name,
    color = "Warning type"
  ) +
  theme_minimal(base_size = 13)+
  scale_color_manual(
    values = c(
      "Black Box Warning" = "black",
      "Withdrawn" = "red",
      "All drugs" = "grey",
      "No Warnings found" = "forestgreen"
    )
  )
warnings_plot

#make plot name
warnings_plot_name<-paste(folder_name,"/",name_used_for_current_analysis,"_by_",as.character(distinct_group),"_warnings.png",sep = "")

#save the plot
ggsave(file.path(warnings_plot_name),warnings_plot, width = 10, height = 8, dpi = 150)

#save the summary
warnings_summary_name<-paste(folder_name,"/",name_used_for_current_analysis,"_by_",as.character(distinct_group),"_warnings.xlsx",sep = "")

write_xlsx(plot_df,warnings_summary_name)


#*#**********************************Data selected from Ryaboshapkina et al *****************************************************************
#*#**********************************Data selected from Ryaboshapkina et al *****************************************************************
#*#**********************************Data selected from Ryaboshapkina et al *****************************************************************

#********************make the changes here*******************************And then below to use old data

# select data type to plot
selected_data<-"all_final_gene_data"
#options-   "all_final_gene_data"   "oncology_gene_data"    "non_oncology_gene_data"

#Change betweeen Cell+ Tissue type and just Cell type enrichment
#Uncomment following to use just cell type based enrichment values. Otherwise it will use enrichment based on Cell+Tissue type by default
dataset_to_use<-"cell_type_only"

#options cell_and_tissue    cell_type_only    tissue_type_only

# Set checking threshold values for specificity
thresholds <- seq(1.4,4,0.25)
x_axis_name<-"Top VS Runners Up Ratio"
y_axis_name<-"N cell specific/N total %"

#*****UNCOMMENT THE LINES USING OLD PHASE DATA DOWN BELOW (AROUND LINE 1090) TO USE OLD PHASE DATA

#********************



# Read enrichment data based **only on Tissue type** data





# And read data from the paper Ryaboshapkina et al table 1

tissue_paper_data_all_genes<-read_excel("./t1_genes_drug_targets_z_scores_from_tissue_paper.xls",sheet = "18,377 protein-coding")

# I am skipping first 2 rows with no data here
tissue_paper_data_drug_related_genes<-read_excel("./t1_genes_drug_targets_z_scores_from_tissue_paper.xls",sheet = "drug targets",skip = 2)
#*******************





all_genes<-get(dataset_to_use)



#****************************************************************************



# filtering data only for genes present in tissue paper

#all genes are in all_genes df

#Get all the genes that were present in tissue paper
all_genes_in_tissue_paper<-all_genes %>%
  filter(gene %in% tissue_paper_data_all_genes$ensembl_gene_id)

#Fillter all the drug related genes mentioned in the paper from that

drug_related_genes_in_the_paper<-all_genes %>%
  filter(gene %in% tissue_paper_data_drug_related_genes$ensembl_gene_id)

# add maximum phase for each of these genes from my data base *************************


# define phase order
phase_levels <- c(
  "1",
  "2",
  "3",
  "4"
)

# collapse table 2 to highest phase per gene
highest_phase_for_each_drug_related_gene <- cleaned_raw_data %>%
  mutate(phase_rank = match(max_phase, phase_levels)) %>%
  group_by(gene) %>%
  slice_max(phase_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene,warning_type, highest_phase = max_phase)

# join onto table 1
drug_related_genes_in_the_paper_with_max_phase_and_warnings <- drug_related_genes_in_the_paper %>%
  left_join(highest_phase_for_each_drug_related_gene, by = c("gene_name" = "gene"))

# Drop the few genes that were not found in new list of drugs
drug_related_genes_in_the_paper_with_max_phase_and_warnings<-drug_related_genes_in_the_paper_with_max_phase_and_warnings%>%na.omit()

#**********************************
# add maximum phase for each of these genes from OLD data base from paper *************************


# define phase order
phase_levels_paper <- c(
  "Target of phase 1 drug",
  "Target of phase 2 drug",
  "Target of phase 3 drug",
  "Target of marketed drug",
  "Target of withdrawn drug"
)

# collapse table 2 to highest phase per gene
tissue_paper_data_drug_related_genes_highest_phase_only <- tissue_paper_data_drug_related_genes %>%
  mutate(phase_rank = match(Category, phase_levels_paper)) %>%
  group_by(gene) %>%
  slice_max(phase_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene,ensembl_gene_id,atLeast1DrugForOncologyIndication, Category = Category)


#**********************************


#add oncology or non oncology drug data  ** AND OLD CATEGORY JUST TO COMPARE****
drug_related_genes_in_the_paper_with_max_phase_and_warnings<-drug_related_genes_in_the_paper_with_max_phase_and_warnings%>%
  left_join(
    tissue_paper_data_drug_related_genes_highest_phase_only%>%select(ensembl_gene_id,atLeast1DrugForOncologyIndication,Category),
    by=c("gene"="ensembl_gene_id"))

# Outcome as ordered factor (for later steps)
drug_related_genes_in_the_paper_with_max_phase_and_warnings <- drug_related_genes_in_the_paper_with_max_phase_and_warnings %>%
  mutate(
    max_phase_ord = factor(highest_phase, levels = sort(unique(na.omit(highest_phase))), ordered = TRUE)
  )

#start here

#all final gene data
used_phase_data<-"new_phase_data"
all_final_gene_data<-drug_related_genes_in_the_paper_with_max_phase_and_warnings

#*****I am switching to use old max phases the tissue paper had just for comparison here. Keep this commented out to use new data*****
#Fixing the plot name to match data
#first I am dropping max_phase_ord column with new phase data and warning type from new data
#Then renaming "Category" to max_phase_ord
#Dropping all the withdrawn drugs to prevent it from messing up data later- because 

# used_phase_data<-"old_phase_data"
# all_final_gene_data<-all_final_gene_data%>%select(-max_phase_ord,-warning_type)
# all_final_gene_data<-all_final_gene_data%>%rename(max_phase_ord=Category)


#only the oncology genes
oncology_gene_data<-all_final_gene_data%>%filter(atLeast1DrugForOncologyIndication==1)

#only the non oncology genes
non_oncology_gene_data<-all_final_gene_data%>%filter(atLeast1DrugForOncologyIndication==0)


#without withdrawn genes
if (used_phase_data=="new_phase_data") {
  final_gene_data_without_withdrawn<-get(selected_data)%>%filter(warning_type != "Withdrawn")
}
if (used_phase_data=="old_phase_data") {
  final_gene_data_without_withdrawn<-get(selected_data)%>%filter(max_phase_ord != "Target of withdrawn drug")
}


#data for withdrawn
if (used_phase_data=="new_phase_data") {
  final_gene_data_withdrawn_only<-get(selected_data)%>%filter(warning_type == "Withdrawn")
}
if (used_phase_data=="old_phase_data") {
  final_gene_data_withdrawn_only<-get(selected_data)%>%filter(max_phase_ord == "Target of withdrawn drug")
}


#**************************Calculating different ratios*************************************************************




#calculating ratios without withdrawn
gene_percentage_df <- final_gene_data_without_withdrawn %>%
  # keep only columns we need
  select(gene_name, max_phase_ord,ratio) %>%
  
  # ensure one row per drug
  distinct(gene_name, .keep_all = TRUE) %>%
  
  # expand thresholds
  crossing(threshold = thresholds) %>%
  
  # compute exceedance
  mutate(exceeds = ratio > threshold) %>%
  
  # aggregate
  group_by(max_phase_ord, threshold) %>%
  summarise(
    n_genes = n(),                                # total drugs in this phase
    n_exceeding = sum(exceeds, na.rm = TRUE),     # NEW: count exceeding threshold
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  )

#calculating ratios just for withdrawn
withdrawn_percentage_df <- final_gene_data_withdrawn_only %>%
  # keep only columns we need
  select(gene_name, max_phase_ord,ratio) %>%
  
  # ensure one row per drug
  distinct(gene_name, .keep_all = TRUE) %>%
  
  # expand thresholds
  crossing(threshold = thresholds) %>%
  
  # compute exceedance
  mutate(exceeds = ratio > threshold) %>%
  
  # aggregate
  group_by(max_phase_ord, threshold) %>%
  summarise(
    n_genes = n(),                                # total drugs in this phase
    n_exceeding = sum(exceeds, na.rm = TRUE),     # NEW: count exceeding threshold
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  )

#change withdrawn max phase to withdrawn
withdrawn_percentage_df<-withdrawn_percentage_df%>%mutate(max_phase_ord="Withdrawn")

# with all genes in the paper
all_genes_percentage <- all_genes_in_tissue_paper %>%
  select(gene_name,ratio) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  crossing(threshold = thresholds) %>%
  mutate(exceeds = ratio > threshold) %>%
  group_by(threshold) %>%
  summarise(
    n_genes = n(),
    n_exceeding = sum(exceeds, na.rm = TRUE),
    percent_exceeding = 100 * mean(exceeds, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(max_phase_ord = "All genes")

#*******************************************************************************************************************************

#combine all these dataframes

gene_plot_df<-bind_rows(gene_percentage_df,withdrawn_percentage_df,all_genes_percentage)

#****************************************plot this********************************************


max_y<- max(gene_plot_df$percent_exceeding[gene_plot_df$threshold>=min(thresholds)])
max_y<-round_up_to_5(max_y)
# max_y<-max(gene_plot_df$percent_exceeding)
# #add a little extra room
# round_up_to_5 <- function(x) {
#   ceiling(x / 5) * 5
# }
# max_y<-round_up_to_5(max_y)


gene_plot<-ggplot(gene_plot_df,
                     aes(x = threshold,
                         y = percent_exceeding,
                         color = max_phase_ord,
                         group = max_phase_ord)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0,max_y)) +
  scale_x_continuous(limits = c(min(thresholds),max(thresholds))) +
  labs(
    title = paste("Specific expression gene percentage\n",selected_data,used_phase_data,dataset_to_use),
    x = x_axis_name,
    y = y_axis_name,
    color = "Max Clinical Phase"
  ) +
  theme_minimal(base_size = 13)+
  scale_color_manual(
    values = c(
      #      "All phases" = "lightblue",
      "All genes" = "black",
      "1" = "pink",
      "2" = "orange",
      "3" = "gold",
      "4" = "forestgreen",
      "Withdrawn"="grey",
      #The following are only for old phase data
      "Target of marketed drug"="forestgreen",
      "Target of phase 1 drug"="pink",
      "Target of phase 2 drug"="orange",
      "Target of phase 3 drug"="gold",
      "Target of withdrawn drug"="grey"
    )
  )

gene_plot


#create a folder for plots
gene_folder_name <- "gene_plots"


if (!dir.exists(gene_folder_name)) {
  dir.create(gene_folder_name)
}

#save the whole summary here
write_xlsx(all_final_gene_data,"gene_plots/full_summary.xlsx")

#create subfolder
gene_subfolder_name <- paste("gene_plots/",selected_data,used_phase_data,dataset_to_use,sep = "")

if (!dir.exists(gene_subfolder_name)) {
  dir.create(gene_subfolder_name)
}

gene_plot_saving_name<-paste("gene_plot_with",used_phase_data,selected_data,dataset_to_use,sep = "_")
gene_plot_saving_path<-paste(gene_subfolder_name,gene_plot_saving_name,sep = "/")
gene_plot_saving_path<-paste(gene_plot_saving_path,".png",sep = "")

#save the plot
ggsave(gene_plot_saving_path,gene_plot, width = 10, height = 8, dpi = 150)

#save summary
gene_summary_saving_name<-paste("summary",used_phase_data,selected_data,dataset_to_use,sep = "_")
gene_summary_saving_path<-paste(gene_subfolder_name,gene_summary_saving_name,sep = "/")
gene_summary_saving_path<-paste(gene_summary_saving_path,".xlsx",sep = "")

write_xlsx(gene_plot_df,gene_summary_saving_path)

#save detailed summary
gene_detailed_summary_saving_name<-paste("detailed_summary",used_phase_data,selected_data,dataset_to_use,sep = "_")
gene_detailed_summary_saving_path<-paste(gene_subfolder_name,gene_detailed_summary_saving_name,sep = "/")
gene_detailed_summary_saving_path<-paste(gene_detailed_summary_saving_path,".xlsx",sep = "")

gene_detailed_summary_df<-get(selected_data)

write_xlsx(gene_detailed_summary_df,gene_detailed_summary_saving_path)

