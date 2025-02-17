
args = commandArgs(trailingOnly = TRUE)

filename.count= args[1]
filename.output_normalized= args[2]
filename.output_full= args[3]
filename.output_ASE= args[3]

print(args)



####  >>
# loading libraries
    <<  ####
suppressMessages(library(DESeq2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))


###
# read the count-matrix
txt = read.table(filename.count, header=TRUE)

#mother_2 excluding
txt_mother_2 <- txt %>% select(-mother_2)

data <- txt_mother_2 %>% as.data.frame()

## Define SampleTypes
sampletype <- factor(c(rep(".mother",3), rep("child_1", 4), rep("child_2", 4), rep("child_3", 4), rep("child_4", 4), rep("child_5", 4), rep("child_6", 4), rep("child_7", 4), rep("child_8", 4), rep("father", 4) ))

meta_full <- data.frame(sampletype, row.names = colnames(data))


dds_full <- DESeqDataSetFromMatrix(data, colData = meta_full, design = ~ sampletype)

keep <- rowSums(counts(dds_full) >= 10) >= 3

dds_full <- dds_full[keep,]


dds_full <- DESeq(dds_full)

#########

# Defining levels to compare ()
children_levels <- c("child_1", "child_2", "child_3", "child_4", "child_5", "child_6", "child_7", "child_8",".mother")

# Creating an empty list to store results
results_list <- list()

# Looping through each child level and comparing to father
for (child_level in children_levels) {
  # Define the contrast
  contrast <- c("sampletype", child_level, "father")
  
  ### Run DESeq analysis
  child_results <- results(dds_full, contrast = contrast)
  
  ## Store
  results_list[[child_level]] <- child_results
}



convert_to_df <- function(result, child_id) {
  df <- data.frame(
    Gene = rownames(result),
    log2FoldChange = result$log2FoldChange,
    padj = result$padj,
    Child = child_id
  )
  return(df)
}

## Applying the function on the mother and children

child_1_df <- convert_to_df(results_list$child_1, "child_1")
child_2_df <- convert_to_df(results_list$child_2, "child_2")
child_3_df <- convert_to_df(results_list$child_3, "child_3")
child_4_df <- convert_to_df(results_list$child_4, "child_4")
child_5_df <- convert_to_df(results_list$child_5, "child_5")
child_6_df <- convert_to_df(results_list$child_6, "child_6")
child_7_df <- convert_to_df(results_list$child_7, "child_7")
child_8_df <- convert_to_df(results_list$child_8, "child_8")
mother_df <- convert_to_df(results_list$.mother, ".mother")

### Combining all data frames as one data frame
combined_df <- bind_rows(
  child_1_df,
  child_2_df,
  child_3_df,
  child_4_df,
  child_5_df,
  child_6_df,
  child_7_df,
  child_8_df,
  mother_df
)



Log2FC_FULL <- pivot_wider(
  data = combined_df,
  names_from = "Child",
  values_from = c("log2FoldChange","padj")
)

### Combine and rename

Log2FC_FULL <- Log2FC_FULL%>%
  rename(
    .mother_log2fc = log2FoldChange_.mother,
    .mother_padj = padj_.mother,
    child_1_log2fc = log2FoldChange_child_1,
    child_1_padj = padj_child_1,
    child_2_log2fc = log2FoldChange_child_2,
    child_2_padj = padj_child_2,
    child_3_log2fc = log2FoldChange_child_3,
    child_3_padj = padj_child_3,
    child_4_log2fc = log2FoldChange_child_4,
    child_4_padj = padj_child_4,
    child_5_log2fc = log2FoldChange_child_5,
    child_5_padj = padj_child_5,
    child_6_log2fc = log2FoldChange_child_6,
    child_6_padj = padj_child_6,
    child_7_log2fc = log2FoldChange_child_7,
    child_7_padj = padj_child_7,
    child_8_log2fc = log2FoldChange_child_8,
    child_8_padj = padj_child_8
  ) %>%
  select(Gene, .mother_log2fc, .mother_padj, starts_with("child"))

print(Log2FC_FULL)

write.table(Log2FC_FULL, filename.output_full, sep='\t',quote = TRUE)
###############
# Get normalised counts from Deseq2
#############

normalized_counts_full <- counts(dds_full, normalized = TRUE)

write.table(normalized_counts_full, filename.output_normalized, sep='\t',quote = TRUE)


################
# Only rows with at least one individual with log2FC >= 1 or at least one individual <= -1
###############
selected_rows <- Log2FC_FULL %>%
  filter(if_any(ends_with("_log2fc"), ~ abs(.) >= 1))

## Output the potential ASE for downstream analysis
write.table(selected_rows, filename.output_ASE, sep='\t',quote = TRUE)

## ### ##

#  <<< Some Plots >>>


plotCounts(dds_full, gene = GENE, intgroup = "sampletype", returnData = TRUE) %>% ggplot() + aes(sampletype, count) + geom_boxplot(aes(fill=sampletype)) + scale_y_log10() + theme_bw() + 
  ggtitle("Boxplot of Counts for Gene" + GENE)


