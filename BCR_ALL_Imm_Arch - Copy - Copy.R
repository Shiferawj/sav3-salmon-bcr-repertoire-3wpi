
##Data based on 3 wpi HK and Sp samples >> Balanced samples

#install.packages("devtools")
#devtools::install_url("https://github.com/immunomind/immunarch/raw/master/immunarch.tar.gz")
#library(immunarch)
install.packages("immunarch")
install.packages("ggpubr")
library(immunarch)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(tibble)
library(ggrepel)
library(scales)
library(purrr)

##Exploratory data analysis
immdata <- repLoad("C:/Users/sje054/OneDrive - UiT Office 365/BCR/3w_HK_Sp_Balanced") # Loading data


#number of clones/seq-depth
Res = repExplore(immdata$data, "clones")   
vis(Res)

#################################New
##script to modify the title of the above figure
clone_plot <- immdata$data %>%
  repExplore("clones") %>%
  vis()

clone_plot +
  labs(
    title = "Productive IGH clones per sample",
    subtitle = NULL,
    x = "Sample",
    y = "Number of clones"
  )
############# New script to change sample names
clone_plot <- immdata$data %>%
  repExplore("clones") %>%
  vis()

clone_plot +
  scale_x_discrete(
    labels = function(x) sub(
      "^([A-Za-z]{3})_([0-9]+).*$",
      "\\1_\\2",
      x
    )
  ) +
  labs(
    title = "Productive IGH clones per sample",
    subtitle = NULL,
    x = "Sample",
    y = "Number of clones"
  )
#############################End

####No. of clonotypes
repExplore(immdata$data, "volume") %>% vis() 

##############################################New
##Script to change the figure title
clonotype_plot <- immdata$data %>%
  repExplore("volume") %>%
  vis() +
  labs(
    title = "Unique productive IGH clonotypes per sample",
    subtitle = NULL,
    x = "Sample",
    y = "Number of clonotypes"
  )

clonotype_plot

###################New script to change the name of the samples
clonotype_plot <- immdata$data %>%
  repExplore("volume") %>%
  vis() +
  scale_x_discrete(
    labels = function(x) sub(
      "^([A-Za-z]{3})_([0-9]+).*$",
      "\\1_\\2",
      x
    )
  ) +
  labs(
    title = "Unique productive IGH clonotypes per sample",
    subtitle = NULL,
    x = "Sample",
    y = "Number of clonotypes"
  )

clonotype_plot

###################################End

repExplore(immdata$data, "len") %>% vis()  #CDR3 length nucleotide
#repExplore(immdata_filtered$data, "len", .col = "aa") %>% vis() ## CDR3 Length in aa

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

##########New to change sample names
exp_cnt <- repExplore(immdata$data, .method = "count")

count_plot <- vis(exp_cnt) +
  scale_color_discrete(
    name = "Sample",
    labels = function(x) {
      sub("^([A-Za-z]{3})_([0-9]+).*$", "\\1_\\2", x)
    }
  ) +
  labs(
    title = "Distribution of clonotype abundances",
    x = "Abundance",
    y = "Number of clonotypes"
  )

count_plot
###########

#####CDR3 distribution for all samples
# 1. Run the immunarch sequence length exploration
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")

# 2. Filter data and clean up the sample names
filtered_len <- exp_len %>%
  filter(Length >= 6 & Length <= 22) %>%
  mutate(Sample = gsub("\\.clones_IGH", "", Sample)) # Removes the suffix from names

# 3. Create the distribution plot
p <- ggplot(filtered_len, aes(x = Length, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  # Set the axis limits from 6 to 22, breaking at intervals of 4
  scale_x_continuous(
    limits = c(6, 24),
    breaks = seq(6, 24, by = 2)
  ) +
  theme_minimal() +
  labs(
    title = "CDR3 Length Distribution",
    x = "CDR3 Length (Amino Acids)",
    y = "Count / Abundance",
    fill = "Samples" 
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right" 
  )

print(p)

############ CDR3 distribution with tisses merged
# 1. Run the immunarch sequence length exploration
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")

# 2. Filter, clean names, group by prefix, and calculate Mean and SD per group/length
# 2. Filter, clean names, group by sample category, and calculate Mean and SD per group/length
filtered_len <- exp_len %>%
  # Filter length range
  filter(Length >= 6 & Length <= 22) %>%
  
  # Remove file suffix
  mutate(Sample = gsub("\\.clones_IGH", "", Sample)) %>%
  
  # Extract group name (HKC, HKI, SpC, SpI)
  mutate(Group = sub("^((HKC|HKI|SpC|SpI)).*", "\\1", Sample)) %>%
  
  # Assign "Other" if sample name doesn't match expected groups
  mutate(Group = ifelse(Group %in% c("HKC", "HKI", "SpC", "SpI"),
                        Group, "Private")) %>%
  
  # Calculate mean and SD across samples within each group
  group_by(Group, Length) %>%
  summarize(
    Avg_Count = mean(Count, na.rm = TRUE),
    SD_Count  = sd(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # Replace NA SD values (when only one sample contributes)
  mutate(SD_Count = ifelse(is.na(SD_Count), 0, SD_Count))

# 3. Create the aggregated distribution plot with upper error bars
p <- ggplot(filtered_len, aes(x = Length, y = Avg_Count, fill = Group)) +
  # Using position_dodge2 adds clean spacing for error bars automatically
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_errorbar(
    aes(ymin = Avg_Count, ymax = Avg_Count + SD_Count),
    position = position_dodge(width = 0.9),
    width = 0.2,
    color = "black", 
    linewidth = 0.5
  ) +
  # Configure axis parameters
  scale_x_continuous(
    limits = c(6, 24),
    breaks = seq(6, 24, by = 2)
  ) +
  theme_minimal() +
  labs(
    title = "E",
    x = "CDR3 Length (Amino Acids)",
    y = "Average Count / Abundance",
    fill = "Tissue Groups"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

print(p)

############CDR3 length per sample

# 1. Calculate the mean and SD of CDR3 amino acid length per sample
stats_list <- lapply(immdata$data, function(sample_df) {
  lengths <- nchar(sample_df$CDR3.aa)
  data.frame(
    Avg_Length = mean(lengths, na.rm = TRUE),
    SD_Length  = sd(lengths, na.rm = TRUE)
  )
})

# 2. Combine data, clean suffix, and group by sample category
avg_df <- bind_rows(stats_list, .id = "Sample") %>%
  mutate(
    Sample = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample)   # Extract HKC, HKI, SpC, or SpI
  ) %>%
  mutate(
    Group = ifelse(Group %in% c("HKC", "HKI", "SpC", "SpI"),
                   Group, "Private")
  )

# 3. Create the bar plot
p <- ggplot(avg_df, aes(x = Sample, y = Avg_Length, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = Avg_Length, ymax = Avg_Length + SD_Length),
    width = 0.2,
    color = "black", 
    linewidth = 0.7
  ) +
  theme_minimal() +
  labs(
    title = "Average CDR3 Length per Sample",
    y = "Average Length (aa)",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p)


####################Top clones by sample
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
vis(imm_top)

##########################################New
##script to modify theabove figure title
clonal_proportion_plot <- repClonality(
  immdata$data,
  "top",
  .head = c(10, 100, 1000, 3000, 10000)
) %>%
  vis() +
  labs(
    title = "Clonal abundance distribution across IGH repertoires",
    subtitle = NULL,
    x = "Sample",
    y = "Proportion of repertoire"
  )

clonal_proportion_plot

################New script to change the sample name for the above script
# Copy the data and shorten sample names for plotting
immdata_plot <- immdata$data
names(immdata_plot) <- sub(
  "^([A-Za-z]{3})_([0-9]+).*$",
  "\\1_\\2",
  names(immdata_plot)
)

clonal_proportion_plot <- repClonality(
  immdata_plot,
  "top",
  .head = c(10, 100, 1000, 3000, 10000)
) %>%
  vis() +
  labs(
    title = "Clonal abundance distribution across IGH repertoires",
    subtitle = NULL,
    x = "Sample",
    y = "Proportion of repertoire"
  )

clonal_proportion_plot
########################################################end

####After down sampling
#########################################################
#repExplore(immdata$data, "count") %>% vis() #distribution of abundances or structure of clones 

#########Immunarch downsample command with 
##set.seed to ensure reproducibility 
set.seed(123)
imm_ds <- repSample(
  immdata$data,
  "downsample"
)

##scrips to check imm_ds content>>Start
imm_ds_summary <- imap_dfr(imm_ds, function(df, sample_name) {
  tibble(
    Sample = sample_name,
    n_clonotypes = nrow(df),
    total_clones = sum(df$Clones, na.rm = TRUE),
    unique_CDR3_aa = n_distinct(df$CDR3.aa),
    mean_clones_per_clonotype = mean(df$Clones, na.rm = TRUE),
    max_clones_per_clonotype = max(df$Clones, na.rm = TRUE)
  )
})

View(imm_ds_summary)
###END

##Downsampling without reproducibility
##imm_ds <- repSample(immdata$data, "downsample")

repSample(immdata$data, "downsample") %>% repExplore("clones") %>% vis()
repSample(immdata$data, "downsample") %>% 
  repExplore("clones") %>% 
  vis() +
  labs(title = "Number of Clones per Downsampled Sample")

##NEW to determin the downsampling n
###########################################################################
############################################################
## Downsample repertoires and report the downsampling depth
############################################################

## Original repertoire sizes
original_sizes <- tibble(
  Sample = names(immdata$data),
  Original_clonotypes = sapply(
    immdata$data,
    nrow
  ),
  Original_clone_count = sapply(
    immdata$data,
    function(x) sum(x$Clones, na.rm = TRUE)
  )
)

print(original_sizes)

## The smallest total clone count determines the automatic
## downsampling depth used by repSample(..., "downsample")
target_depth <- min(original_sizes$Original_clone_count, na.rm = TRUE)

cat(
  "\nAutomatic downsampling depth:",
  target_depth,
  "clones per sample\n\n"
)

## Sizes after downsampling
downsampled_sizes <- tibble(
  Sample = names(imm_ds),
  Downsampled_clonotypes = sapply(
    imm_ds,
    nrow
  ),
  Downsampled_clone_count = sapply(
    imm_ds,
    function(x) sum(x$Clones, na.rm = TRUE)
  )
)

## Combine original and downsampled summaries
downsampling_summary <- original_sizes %>%
  left_join(
    downsampled_sizes,
    by = "Sample"
  ) %>%
  mutate(
    Target_depth = target_depth,
    Percent_clones_retained =
      100 * Downsampled_clone_count / Original_clone_count,
    Percent_clonotypes_retained =
      100 * Downsampled_clonotypes / Original_clonotypes
  )

print(downsampling_summary)
######################################################################end of the above script

############################################################
## CLEAN DIVERSITY / CLONALITY FIGURE SCRIPT
## Samples: HKC, HKI, SpC, SpI
## Input objects required:
##   imm_ds
##   div_res

library(ggpubr)
library(patchwork)

group_levels <- c("HKC", "HKI", "SpC", "SpI")

get_group <- function(x) {
  x <- gsub("\\.clones_IGH", "", x)
  sub("_.*", "", x)
}

############################################################
## 1. Chao1
############################################################

div_res <- repDiversity(
  imm_ds,
  .method = "chao1"
)

div_df <- data.frame(
  Sample = rownames(div_res),
  Chao1 = div_res[, "Estimator"]
) %>%
  mutate(
    Sample = gsub("\\.clones_IGH", "", Sample),
    Tissue = factor(get_group(Sample), levels = group_levels)
  ) %>%
  filter(!is.na(Tissue))

table(div_df$Tissue)

############################################################
## 2. Hill q = 5
############################################################

div5_df <- repDiversity(
  imm_ds,
  .method = "div",
  .q = 5,
  .do.norm = NA,
  .laplace = 0
) %>%
  mutate(
    Sample = gsub("\\.clones_IGH", "", Sample),
    Tissue = factor(get_group(Sample), levels = group_levels)
  ) %>%
  filter(!is.na(Tissue))

############################################################
## 3. Hill q = 2
############################################################

div2_df <- repDiversity(
  imm_ds,
  .method = "div",
  .q = 2,
  .do.norm = NA,
  .laplace = 0
) %>%
  mutate(
    Sample = gsub("\\.clones_IGH", "", Sample),
    Tissue = factor(get_group(Sample), levels = group_levels)
  ) %>%
  filter(!is.na(Tissue))

############################################################
## 4. Gini
############################################################

gini <- repDiversity(
  imm_ds,
  .method = "gini"
)

gini_df <- data.frame(
  Sample = rownames(gini),
  Gini = gini[, 1]
) %>%
  mutate(
    Sample = gsub("\\.clones_IGH", "", Sample),
    Tissue = factor(get_group(Sample), levels = group_levels)
  ) %>%
  filter(!is.na(Tissue))

############################################################
## 5. Check all sample numbers
############################################################

table(div_df$Tissue)
table(div5_df$Tissue)
table(div2_df$Tissue)
table(gini_df$Tissue)

############################################################
## 6. Statistics
############################################################

kruskal.test(Chao1 ~ Tissue, data = div_df)
kruskal.test(Value ~ Tissue, data = div5_df)
kruskal.test(Value ~ Tissue, data = div2_df)
kruskal.test(Gini ~ Tissue, data = gini_df)

pairwise.wilcox.test(div_df$Chao1, div_df$Tissue, p.adjust.method = "BH")
pairwise.wilcox.test(div5_df$Value, div5_df$Tissue, p.adjust.method = "BH")
pairwise.wilcox.test(div2_df$Value, div2_df$Tissue, p.adjust.method = "BH")
pairwise.wilcox.test(gini_df$Gini, gini_df$Tissue, p.adjust.method = "BH")

############################################################
## 7. Plot theme
############################################################

theme_journal <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

############################################################
## 8. Four-panel figure
############################################################

pA <- ggplot(div_df, aes(Tissue, Chao1, fill = Tissue)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.8) +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "A: Clonotype Richness",
    x = "",
    y = "Chao1"
  ) +
  theme_journal

pB <- ggplot(div5_df, aes(Tissue, Value, fill = Tissue)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.8) +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "B: Dominance-weighted Diversity",
    x = "",
    y = "Hill q = 5"
  ) +
  theme_journal

pC <- ggplot(div2_df, aes(Tissue, Value, fill = Tissue)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.8) +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "C: Effective Diversity",
    x = "",
    y = "Hill q = 2"
  ) +
  theme_journal

pD <- ggplot(gini_df, aes(Tissue, Gini, fill = Tissue)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.8) +
  stat_compare_means(method = "kruskal.test") +
  labs(
    title = "D: Clonal Inequality",
    x = "",
    y = "Gini index"
  ) +
  theme_journal

fig1 <- (pA | pB) / (pC | pD)

fig1


######################################################
#analysis based on genes
##########################################################

###Collapse V genes to family level>>removes quality score and allel level info
## Collapse V-gene names to family labels
clean_v_family <- function(x) {
  x <- sub("\\(.*\\)", "", x)       # remove parenthetical quality score
  x <- sub("\\*.*$", "", x)         # remove allele designation
  x <- sub("-.*$", "", x)           # remove gene-level suffix after hyphen
  x <- sub("_?S[0-9]+$", "", x)     # remove terminal S plus digits
  x
}

##Build V-family abundance table
v_df <- lapply(names(imm_ds), function(samp) {
  
  df <- imm_ds[[samp]]
  
  df$V.clean <- clean_v_family(df$V.name)
  
  df %>%
    group_by(V.clean) %>%
    summarise(
      Proportion = sum(Proportion, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Sample = sub("\\.clones_IGH$", "", samp)
    )
  
}) %>%
  bind_rows()

###########the above script but script to remove the training 'S_'
clean_v_family <- function(x) {
  x <- as.character(x)
  
  # Remove annotations in parentheses and allele information
  x <- sub("\\(.*\\)", "", x)
  x <- sub("\\*.*$", "", x)
  
  # Remove a trailing S-number suffix: IGHV1DS3 -> IGHV1D
  x <- sub("S_?[0-9]+$", "", x, ignore.case = TRUE)
  
  # Collapse gene-level names to family level: IGHV1-100 -> IGHV1
  x <- sub("-.*$", "", x)
  
  x
}

# Build the V-family abundance table
v_df <- lapply(names(imm_ds), function(samp) {
  df <- imm_ds[[samp]]
  
  df$V.clean <- clean_v_family(df$V.name)
  
  df %>%
    group_by(V.clean) %>%
    summarise(
      Proportion = sum(Proportion, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Sample = sub("\\.clones_IGH$", "", samp)
    )
}) %>%
  bind_rows()
#############end of trailing 'S' script

##Select top V genes overall
top_v <- v_df %>%
  group_by(V.clean) %>%
  summarise(Total = sum(Proportion)) %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  pull(V.clean)

##Create heatmap matrix
v_mat <- v_df %>%
  filter(V.clean %in% top_v) %>%
  select(Sample, V.clean, Proportion) %>%
  pivot_wider(
    names_from = Sample,
    values_from = Proportion,
    values_fill = 0
  )

mat <- as.matrix(v_mat[, -1])
rownames(mat) <- v_mat$V.clean

##Z-score normalize rows
mat_z <- t(scale(t(mat)))

##Clustered heatmap
pheatmap(
  mat_z,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  main = "Top 20 V-Gene Usage-Proportion"
)


##Build TRUE public clonotype table >> Extract CDR3s across all samples
## because Immunarch does not compute that via "public-clones".
all_clones <- lapply(names(imm_ds), function(samp) {
  df <- imm_ds[[samp]][, c("CDR3.aa", "Clones")]
  df$Sample <- samp
  df
}) %>% dplyr::bind_rows()

##Find public clones (appear in >1 sample)
public_table <- all_clones %>%
  group_by(CDR3.aa) %>%
  summarise(
    n_samples = n_distinct(Sample),
    total_abundance = sum(Clones)
  ) %>%
  filter(n_samples > 1) %>%
  arrange(desc(n_samples), desc(total_abundance))
##Heatmap of public clone distribution
public_mat <- all_clones %>%
  filter(CDR3.aa %in% public_table$CDR3.aa) %>%
  group_by(CDR3.aa, Sample) %>%
  summarise(Clones = sum(Clones), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Sample,
                     values_from = Clones,
                     values_fill = 0)

mat <- as.matrix(public_mat[, -1])
rownames(mat) <- public_mat$CDR3.aa

##Clustered heatmap
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  angle_col = 45
)

########the above script (from 701-737) but with the appropriate CDR3 QCs
##Productive
##C/N-terminus anchors

# The 12 balanced samples used in the paper
balanced_samples <- c(
  "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
  "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
  "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
  "SpI_1I3w", "SpI_2I3w", "SpI_5I3w"
)

# Check that the required samples are present
available_samples <- sub("\\.clones_IGH$", "", names(imm_ds))
missing_samples <- setdiff(balanced_samples, available_samples)

if (length(missing_samples) > 0) {
  stop(
    "These required samples are missing from imm_ds: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Extract CDR3 amino-acid sequences and MiXCR clone counts
all_clones <- lapply(names(imm_ds), function(samp) {
  df <- imm_ds[[samp]]
  
  required_columns <- c("CDR3.aa", "Clones")
  missing_columns <- setdiff(required_columns, names(df))
  
  if (length(missing_columns) > 0) {
    stop(
      "Sample ", samp, " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data.frame(
    Sample = samp,
    CDR3.aa = as.character(df$CDR3.aa),
    Clones = suppressWarnings(as.numeric(as.character(df$Clones)))
  )
}) %>%
  bind_rows() %>%
  mutate(
    Sample_clean = sub("\\.clones_IGH$", "", Sample),
    Group = sub("_.*$", "", Sample_clean)
  ) %>%
  filter(Sample_clean %in% balanced_samples)

# QC: retain standard amino acids, both CDR3 anchors,
# and positive, finite MiXCR clone counts
all_clones <- all_clones %>%
  filter(
    !is.na(CDR3.aa),
    grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", CDR3.aa),
    grepl("^C.*[WF]$", CDR3.aa),
    !is.na(Clones),
    is.finite(Clones),
    Clones > 0
  )

# Check that valid clonotypes remain for all 12 samples
retained_samples <- unique(all_clones$Sample_clean)
samples_without_valid_clonotypes <- setdiff(
  balanced_samples,
  retained_samples
)

if (length(samples_without_valid_clonotypes) > 0) {
  stop(
    "No clonotypes passed the QC filters for: ",
    paste(samples_without_valid_clonotypes, collapse = ", ")
  )
}

# Identify public CDR3 amino-acid sequences within each group:
# present in at least 2 of the 3 samples in that group
public_table <- all_clones %>%
  group_by(Group, CDR3.aa) %>%
  summarise(
    n_samples = n_distinct(Sample_clean),
    total_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_samples >= 2) %>%
  arrange(Group, desc(n_samples), desc(total_abundance))

# Build a heatmap of public CDR3s across all 12 samples
public_cdr3s <- unique(public_table$CDR3.aa)

public_mat <- all_clones %>%
  filter(CDR3.aa %in% public_cdr3s) %>%
  group_by(CDR3.aa, Sample_clean) %>%
  summarise(
    Clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    CDR3.aa,
    Sample_clean = balanced_samples,
    fill = list(Clones = 0)
  ) %>%
  pivot_wider(
    names_from = Sample_clean,
    values_from = Clones,
    values_fill = 0
  )

# Convert to a matrix and shorten sample labels
mat <- as.matrix(public_mat[, balanced_samples, drop = FALSE])
rownames(mat) <- public_mat$CDR3.aa

colnames(mat) <- sub(
  "^([A-Za-z]{3})_([0-9]+).*$",
  "\\1_\\2",
  colnames(mat)
)

# Draw the heatmap
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  angle_col = 45
)
###############################End of the above script
##
########Identify tissue-specific vs public clones
##Add tissue annotation
meta <- data.frame(
  Sample = names(imm_ds),
  Tissue = c(
    rep("HKC", 3),
    rep("HKI", 3),
    rep("SpC", 3),
    rep("SpI", 3)
  )
)

meta
#View(meta)
str(meta)
##Merge with public clones
public_annot <- public_table %>%
  left_join(all_clones, by = "CDR3.aa") %>%
  left_join(meta, by = "Sample")

#View(public_annot)

##Public clone enrichment by tissue
tissue_public <- public_annot %>%
  group_by(Tissue, CDR3.aa) %>%
  summarise(total = sum(Clones), .groups = "drop")
#view(tissue_public)

####
tissue_public %>%
  arrange(desc(total)) %>%
  head(100) %>%
  print(n = 100)

###
#table(tissue_public$Tissue)

###################New script that replaces 870-903 so that it
##it accommodates the above update
# Add tissue and condition annotations, joining public clonotypes
# to sample-level counts within the same group
public_annot <- public_table %>%
  inner_join(
    all_clones,
    by = c("Group", "CDR3.aa")
  ) %>%
  mutate(
    Tissue = if_else(grepl("^HK", Group), "HK", "Spleen"),
    Condition = if_else(grepl("I$", Group), "Infected", "Control")
  )

# Summarize abundance of public clonotypes in each
# tissue-condition group
tissue_public <- public_annot %>%
  group_by(Tissue, Condition, CDR3.aa) %>%
  summarise(
    total = sum(Clones, na.rm = TRUE),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Tissue, Condition, desc(total))

# Display the 100 most abundant group-public clonotypes
tissue_public %>%
  group_by(Tissue, Condition) %>%
  slice_max(total, n = 100, with_ties = FALSE) %>%
  ungroup() %>%
  print(n = 100)

# Classify each CDR3 according to the tissues in which it
# met the public-clonotype criterion
public_tissue_status <- public_table %>%
  mutate(
    Tissue = if_else(grepl("^HK", Group), "HK", "Spleen"),
    Condition = if_else(grepl("I$", Group), "Infected", "Control")
  ) %>%
  group_by(CDR3.aa) %>%
  summarise(
    public_groups = paste(
      sort(unique(as.character(Group))),
      collapse = "; "
    ),
    public_tissues = paste(
      sort(unique(Tissue)),
      collapse = "; "
    ),
    public_conditions = paste(
      sort(unique(Condition)),
      collapse = "; "
    ),
    public_status = case_when(
      any(Tissue == "HK") & any(Tissue == "Spleen") ~
        "Public in both tissues",
      any(Tissue == "HK") ~
        "Public in HK",
      any(Tissue == "Spleen") ~
        "Public in spleen"
    ),
    .groups = "drop"
  )

print(public_tissue_status)
##################end of the above script
##
############################################################
## Public/private clone analysis for IGH repertoire
## Samples: HKC, HKI, SpC, SpI
## Clone definition: identical CDR3.aa
## Relaxed public clone: CDR3.aa present in >= 2 samples within a group
############################################################
############################################################
## 1. Build metadata
############################################################

meta <- data.frame(
  Sample = names(imm_ds)
) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI"))
  )

meta
table(meta$Group)

############################################################
## 2. Rebuild all_clones with V/J-gene information
############################################################

all_clones <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa,
        CDR3.nt,
        Clones,
        Proportion,
        V.name,
        D.name,
        J.name
      ) %>%
      mutate(Sample = s)
  })
)

all_clones_clean <- all_clones %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    CDR3_length = nchar(CDR3.aa)
  )

#View(all_clones_clean)

############################################################
## 3. Public clones per group: relaxed definition
## Public = same CDR3.aa present in at least 2 samples
############################################################

get_public_clones_relaxed <- function(group_name, min_samples = 2) {
  
  all_clones_clean %>%
    filter(Group == group_name) %>%
    group_by(CDR3.aa) %>%
    summarise(
      n_samples = n_distinct(Sample),
      total_clones = sum(Clones, na.rm = TRUE),
      mean_clones = mean(Clones, na.rm = TRUE),
      total_proportion = sum(Proportion, na.rm = TRUE),
      mean_proportion = mean(Proportion, na.rm = TRUE),
      CDR3_length = first(CDR3_length),
      .groups = "drop"
    ) %>%
    filter(n_samples >= min_samples) %>%
    mutate(Group = group_name)
}

HKC_public <- get_public_clones_relaxed("HKC", min_samples = 2)
HKI_public <- get_public_clones_relaxed("HKI", min_samples = 2)
SpC_public <- get_public_clones_relaxed("SpC", min_samples = 2)
SpI_public <- get_public_clones_relaxed("SpI", min_samples = 2)

public_by_group <- bind_rows(
  HKC_public,
  HKI_public,
  SpC_public,
  SpI_public
)

#View(public_by_group)

############################################################
## 4. Count public clones per group
############################################################

public_counts <- public_by_group %>%
  group_by(Group) %>%
  summarise(
    n_public_clones = n_distinct(CDR3.aa),
    .groups = "drop"
  )

public_counts

ggplot(public_counts, aes(x = Group, y = n_public_clones, fill = Group)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Relaxed Public Clones per Group",
    subtitle = "Public = CDR3.aa shared by at least 2 samples within group",
    x = "Group",
    y = "Number of public clones"
  )

############################################################
## 5. Top public clones by abundance
############################################################

top_public_clones <- public_by_group %>%
  group_by(Group) %>%
  arrange(desc(total_clones), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

#View(top_public_clones)

ggplot(top_public_clones,
       aes(x = reorder(CDR3.aa, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top Relaxed Public Clones per Group",
    x = "CDR3 amino acid sequence",
    y = "Total clone count"
  )

############################################################
## 6. Annotate public clones with V/J genes
############################################################

public_annot <- public_by_group %>%
  left_join(
    all_clones_clean,
    by = c("CDR3.aa", "Group")
  )

#View(public_annot)

############################################################
## 7. V-gene usage among public clones
############################################################

public_v_usage <- public_annot %>%
  group_by(Group, V.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

#View(public_v_usage)

top_public_v_usage <- public_v_usage %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 15, with_ties = FALSE) %>%
  ungroup()

ggplot(top_public_v_usage,
       aes(x = reorder(V.name, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top V-Gene Usage Among Public Clones",
    x = "V gene",
    y = "Total clone count"
  )

############################################################
## 8. J-gene usage among public clones
############################################################

public_j_usage <- public_annot %>%
  group_by(Group, J.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

#View(public_j_usage)

top_public_j_usage <- public_j_usage %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 15, with_ties = FALSE) %>%
  ungroup()

ggplot(top_public_j_usage,
       aes(x = reorder(J.name, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top J-Gene Usage Among Public Clones",
    x = "J gene",
    y = "Total clone count"
  )

############################################################
## 9. CDR3 length distribution of public clones
############################################################

public_lengths <- public_by_group %>%
  mutate(CDR3_length = nchar(CDR3.aa))

ggplot(public_lengths,
       aes(x = CDR3_length, fill = Group)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 Length Distribution of Public Clones",
    x = "CDR3 amino acid length",
    y = "Number of public clones"
  )

ggplot(public_lengths,
       aes(x = Group, y = CDR3_length, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 Length of Public Clones",
    x = "Group",
    y = "CDR3 amino acid length"
  )

############################################################
## 10. Shared public clones between control and infected
############################################################

HK_shared_public <- intersect(HKC_public$CDR3.aa, HKI_public$CDR3.aa)
Sp_shared_public <- intersect(SpC_public$CDR3.aa, SpI_public$CDR3.aa)

HKC_specific_public <- setdiff(HKC_public$CDR3.aa, HKI_public$CDR3.aa)
HKI_specific_public <- setdiff(HKI_public$CDR3.aa, HKC_public$CDR3.aa)

SpC_specific_public <- setdiff(SpC_public$CDR3.aa, SpI_public$CDR3.aa)
SpI_specific_public <- setdiff(SpI_public$CDR3.aa, SpC_public$CDR3.aa)

public_category_summary <- data.frame(
  Category = c(
    "HKC public",
    "HKI public",
    "SpC public",
    "SpI public",
    "HKC-HKI shared public",
    "SpC-SpI shared public",
    "HKC-specific public",
    "HKI-specific public",
    "SpC-specific public",
    "SpI-specific public"
  ),
  N_clones = c(
    n_distinct(HKC_public$CDR3.aa),
    n_distinct(HKI_public$CDR3.aa),
    n_distinct(SpC_public$CDR3.aa),
    n_distinct(SpI_public$CDR3.aa),
    length(HK_shared_public),
    length(Sp_shared_public),
    length(HKC_specific_public),
    length(HKI_specific_public),
    length(SpC_specific_public),
    length(SpI_specific_public)
  )
)

public_category_summary

ggplot(public_category_summary,
       aes(x = reorder(Category, N_clones), y = N_clones)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Summary of Public Clone Categories",
    x = "",
    y = "Number of clones"
  )



############################################################
## 11. Private clones
## Private = CDR3.aa present in only one sample overall
############################################################

clone_occurrence_all <- all_clones_clean %>%
  group_by(CDR3.aa) %>%
  summarise(
    n_samples_total = n_distinct(Sample),
    total_clones = sum(Clones, na.rm = TRUE),
    mean_clones = mean(Clones, na.rm = TRUE),
    .groups = "drop"
  )

private_clones <- clone_occurrence_all %>%
  filter(n_samples_total == 1)

public_overall_clones <- clone_occurrence_all %>%
  filter(n_samples_total >= 2)

nrow(private_clones)
nrow(public_overall_clones)

############################################################
## 12. Private clone counts per group
############################################################

private_annot <- all_clones_clean %>%
  semi_join(private_clones, by = "CDR3.aa")

private_counts <- private_annot %>%
  group_by(Group, Sample) %>%
  summarise(
    n_private_clones = n_distinct(CDR3.aa),
    total_private_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

#View(private_counts)

ggplot(private_counts,
       aes(x = Group, y = n_private_clones, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "Private Clone Counts per Group",
    x = "Group",
    y = "Number of private clones"
  )

####

############################################################
## 13. Public/private proportion per sample
############################################################

sample_public_private <- all_clones_clean %>%
  mutate(
    Clone_type = ifelse(CDR3.aa %in% private_clones$CDR3.aa,
                        "Private", "Shared/Public")
  ) %>%
  group_by(Group, Sample, Clone_type) %>%
  summarise(
    n_clones = n_distinct(CDR3.aa),
    total_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

#View(sample_public_private)

ggplot(sample_public_private,
       aes(x = Sample, y = n_clones, fill = Clone_type)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Private vs Shared/Public Clones per Sample",
    x = "Sample",
    y = "Number of clones"
  )

############################################################
## 16. Clones shared by infected groups but absent from controls
############################################################

infected_public <- union(HKI_public$CDR3.aa, SpI_public$CDR3.aa)
control_public  <- union(HKC_public$CDR3.aa, SpC_public$CDR3.aa)

infected_specific_public <- setdiff(infected_public, control_public)
control_specific_public  <- setdiff(control_public, infected_public)

length(infected_specific_public)
length(control_specific_public)

infected_specific_table <- public_annot %>%
  filter(CDR3.aa %in% infected_specific_public)

#View(infected_specific_table)

############################################################
## Figure: infected-specific vs control-specific public clones
############################################################
## 1. Define clone categories
infected_public <- union(HKI_public$CDR3.aa, SpI_public$CDR3.aa)
control_public  <- union(HKC_public$CDR3.aa, SpC_public$CDR3.aa)

infected_specific_public <- setdiff(infected_public, control_public)
control_specific_public  <- setdiff(control_public, infected_public)
shared_public            <- intersect(infected_public, control_public)

## 2. Summary table
public_specific_summary <- data.frame(
  Category = c(
    "Infected-specific public",
    "Control-specific public",
    "Shared public"
  ),
  N_clones = c(
    length(infected_specific_public),
    length(control_specific_public),
    length(shared_public)
  )
)

public_specific_summary

## 3. Bar plot
p_public_specific <- ggplot(
  public_specific_summary,
  aes(x = Category, y = N_clones, fill = Category)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = N_clones),
    vjust = -0.4,
    size = 5
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Condition-Specific Public Clones",
    #subtitle = "Public clones shared by ???2 samples within infected or control groups",
    x = "",
    y = "Number of public CDR3.aa clones"
  )

p_public_specific


############################################################
## 17. V-gene usage of infected-specific public clones
############################################################

infected_specific_v_usage <- infected_specific_table %>%
  group_by(Group, V.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

#View(infected_specific_v_usage)

ggplot(infected_specific_v_usage,
       aes(x = reorder(V.name, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "V-Gene Usage of Infected-Specific Public Clones",
    x = "V gene",
    y = "Total clone count"
  )

##Top 15 infected
top15_v_usage_pct <- infected_specific_v_usage %>%
  group_by(Group) %>%
  mutate(
    Percent = 100 * total_clones / sum(total_clones)
  ) %>%
  slice_max(
    Percent,
    n = 30,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top15_v_usage_pct,
  aes(
    x = reorder(V.name, Percent),
    y = Percent,
    fill = Group
  )
) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Top 30 V Genes Among Infection-Specific Public Clones",
    x = "V gene",
    y = "Relative abundance (%)"
  )


####Top 15 V genes by total clone abundance
## Top 15 V genes per group
top15_public_v <- public_v_usage %>%
  group_by(Group) %>%
  slice_max(
    order_by = total_clones,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top15_public_v,
  aes(
    x = reorder(V.name, total_clones),
    y = total_clones,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 15 V Genes Among Public Clones",
    x = "V gene",
    y = "Total clone count"
  )

###Top 15 V genes by number of public clones
##Sometimes this is biologically more informative than abundance:
## Top 15 V genes by number of public CDR3s
top15_v_nclones <- public_v_usage %>%
  group_by(Group) %>%
  slice_max(
    order_by = n_public_CDR3,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top15_v_nclones,
  aes(
    x = reorder(V.name, n_public_CDR3),
    y = n_public_CDR3,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 15 V-Gene Usage of Infected-Specific Public Clones",
    x = "V gene",
    y = "Number of public CDR3 clones"
  )

##Even better: compare the same top V genes across groups
#This often makes differences easier to see:
## Top 15 V genes overall
top15_overall_v <- public_v_usage %>%
  group_by(V.name) %>%
  summarise(
    total_clones = sum(total_clones),
    .groups = "drop"
  ) %>%
  slice_max(total_clones, n = 30)

plot_df <- public_v_usage %>%
  filter(V.name %in% top15_overall_v$V.name)

ggplot(
  plot_df,
  aes(
    x = V.name,
    y = total_clones,
    fill = Group
  )
) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 30 Public V Genes Across All Groups",
    x = "V gene",
    y = "Total clone count"
  )

############################################################
###########New script to replace from 972-1580 and compatible with 
##the updated rules
############################################################
## 1. Define the balanced sample set
############################################################

balanced_samples <- c(
  "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
  "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
  "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
  "SpI_1I3w", "SpI_2I3w", "SpI_5I3w"
)

sample_metadata <- data.frame(
  Sample = names(imm_ds),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Sample_clean = sub("\\.clones_IGH$", "", Sample),
    Group = sub("_.*$", "", Sample_clean)
  ) %>%
  filter(Sample_clean %in% balanced_samples)

missing_samples <- setdiff(
  balanced_samples,
  sample_metadata$Sample_clean
)

if (length(missing_samples) > 0) {
  stop(
    "Missing required samples in imm_ds: ",
    paste(missing_samples, collapse = ", ")
  )
}

sample_check <- sample_metadata %>%
  distinct(Sample_clean, Group) %>%
  count(Group, name = "n_samples")

if (
  n_distinct(sample_metadata$Sample_clean) != 12 ||
  any(sample_check$n_samples != 3)
) {
  stop("Expected 12 balanced samples, with 3 samples in each group.")
}

############################################################
## 2. Define V- and J-family cleaning functions
############################################################

clean_v_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\(.*\\)", "", x)                 # Remove parentheses/annotations
  x <- sub("\\*.*$", "", x)                   # Remove allele information
  x <- sub("S_?[0-9]+$", "", x, ignore.case = TRUE)
  x <- sub("-.*$", "", x)                     # Collapse gene to family
  x
}

clean_j_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\(.*\\)", "", x)                 # Remove parentheses/annotations
  x <- sub("\\*.*$", "", x)                   # Remove allele information
  x <- sub("T_?[0-9]+(?=D?$)", "", x, perl = TRUE)
  x <- sub("-.*$", "", x)                     # Collapse gene to family
  x
}

# Examples:
# IGHV1DS3 -> IGHV1D
# IGHV1S1  -> IGHV1
# IGHJ1T2D -> IGHJ1D
# IGHJ2T5  -> IGHJ2

############################################################
## 3. Build one QC-filtered clonotype table
############################################################
## Assumes imm_ds contains productive MiXCR clonotypes.
## CDR3.aa must include the N-terminal C and terminal W/F anchors.

sample_names_to_use <- sample_metadata$Sample

all_clones_clean <- bind_rows(lapply(sample_names_to_use, function(samp) {
  df <- imm_ds[[samp]]
  
  required_columns <- c(
    "CDR3.aa", "CDR3.nt", "Clones", "Proportion",
    "V.name", "J.name"
  )
  
  missing_columns <- setdiff(required_columns, names(df))
  
  if (length(missing_columns) > 0) {
    stop(
      "Sample ", samp, " is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Preserve D.name if available; it is not required for these analyses.
  if (!"D.name" %in% names(df)) {
    df$D.name <- NA_character_
  }
  
  df %>%
    select(
      CDR3.aa, CDR3.nt, Clones, Proportion,
      V.name, D.name, J.name
    ) %>%
    mutate(Sample = samp)
})) %>%
  mutate(
    Sample_clean = sub("\\.clones_IGH$", "", Sample),
    Group = sub("_.*$", "", Sample_clean),
    Tissue = if_else(grepl("^HK", Group), "HK", "Spleen"),
    Condition = if_else(grepl("I$", Group), "Infected", "Control"),
    
    CDR3.aa = as.character(CDR3.aa),
    Clones = suppressWarnings(as.numeric(as.character(Clones))),
    Proportion = suppressWarnings(as.numeric(as.character(Proportion))),
    CDR3_length = nchar(CDR3.aa),
    
    # Preserve raw V.name and J.name; add cleaned family-level columns.
    V_family = clean_v_family(V.name),
    J_family = clean_j_family(J.name)
  ) %>%
  filter(
    Sample_clean %in% balanced_samples,
    !is.na(CDR3.aa),
    grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", CDR3.aa),
    grepl("^C.*[WF]$", CDR3.aa),
    !is.na(Clones),
    is.finite(Clones),
    Clones > 0
  ) %>%
  mutate(
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

# Confirm that QC retained clonotypes from all 12 samples.
retained_samples <- unique(all_clones_clean$Sample_clean)
samples_without_valid_clonotypes <- setdiff(
  balanced_samples,
  retained_samples
)

if (length(samples_without_valid_clonotypes) > 0) {
  stop(
    "No clonotypes passed QC for: ",
    paste(samples_without_valid_clonotypes, collapse = ", ")
  )
}

# Check the cleaned family names.
family_check <- all_clones_clean %>%
  distinct(V.name, V_family, J.name, J_family)

print(family_check, n = 20)

############################################################
## 4. Identify public CDR3 clonotypes within each group
############################################################
## Public = same CDR3.aa detected in at least 2 of 3 samples
## within the same group.

public_by_group <- all_clones_clean %>%
  group_by(Group, CDR3.aa) %>%
  summarise(
    n_samples = n_distinct(Sample_clean),
    total_clones = sum(Clones, na.rm = TRUE),
    mean_clones = mean(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    mean_proportion = mean(Proportion, na.rm = TRUE),
    CDR3_length = dplyr::first(CDR3_length),
    .groups = "drop"
  ) %>%
  filter(n_samples >= 2) %>%
  arrange(Group, desc(total_clones))

# Count public CDR3 clonotypes in each group.
public_counts <- public_by_group %>%
  group_by(Group) %>%
  summarise(
    n_public_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  )

print(public_counts)

ggplot(public_counts, aes(x = Group, y = n_public_CDR3, fill = Group)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Public CDR3 clonotypes per group",
    subtitle = "Public = detected in at least 2 of 3 samples within a group",
    x = "Group",
    y = "Number of public CDR3 clonotypes"
  )

############################################################
## 5. Annotate group-public clonotypes with sample-level data
############################################################
## Join by both Group and CDR3.aa to avoid matching a sequence
## to records from unrelated groups.

public_annot <- public_by_group %>%
  inner_join(
    all_clones_clean,
    by = c("Group", "CDR3.aa")
  )

############################################################
## 6. Top public CDR3 clonotypes by abundance
############################################################

top_public_clones <- public_by_group %>%
  group_by(Group) %>%
  slice_max(
    order_by = total_clones,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top_public_clones,
  aes(
    x = reorder(CDR3.aa, total_clones),
    y = total_clones,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top public CDR3 clonotypes per group",
    x = "CDR3 amino-acid sequence",
    y = "Total clone count"
  )

############################################################
## 7. V-family and J-family usage among public clonotypes
############################################################

public_v_family_usage <- public_annot %>%
  filter(!is.na(V_family), V_family != "") %>%
  group_by(Group, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

public_j_family_usage <- public_annot %>%
  filter(!is.na(J_family), J_family != "") %>%
  group_by(Group, J_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

top_public_v_families <- public_v_family_usage %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 15, with_ties = FALSE) %>%
  ungroup()

ggplot(
  top_public_v_families,
  aes(x = reorder(V_family, total_clones), y = total_clones, fill = Group)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top V-family usage among public clonotypes",
    x = "V family",
    y = "Total clone count"
  )

top_public_j_families <- public_j_family_usage %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 15, with_ties = FALSE) %>%
  ungroup()

ggplot(
  top_public_j_families,
  aes(x = reorder(J_family, total_clones), y = total_clones, fill = Group)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top J-family usage among public clonotypes",
    x = "J family",
    y = "Total clone count"
  )

############################################################
## 8. Retain gene-level usage tables for possible future analyses
############################################################
## These use the original gene calls, including their original
## suffixes and allele-level naming.

public_v_gene_usage <- public_annot %>%
  filter(!is.na(V.name), V.name != "") %>%
  group_by(Group, V.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

public_j_gene_usage <- public_annot %>%
  filter(!is.na(J.name), J.name != "") %>%
  group_by(Group, J.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

############################################################
## 9. CDR3 length distribution of public clonotypes
############################################################
## Apply the 6-22 aa interval specifically to this length
## analysis; it is not used to filter the other analyses.

public_lengths <- public_by_group %>%
  filter(CDR3_length >= 6, CDR3_length <= 22)

ggplot(public_lengths, aes(x = CDR3_length, fill = Group)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 length distribution of public clonotypes",
    x = "CDR3 amino-acid length",
    y = "Number of public clonotypes"
  )

ggplot(public_lengths, aes(x = Group, y = CDR3_length, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 length of public clonotypes",
    x = "Group",
    y = "CDR3 amino-acid length"
  )

############################################################
## 10. Public clonotype representation by tissue
############################################################

public_tissue_status <- public_by_group %>%
  mutate(
    Tissue = if_else(grepl("^HK", Group), "HK", "Spleen"),
    Condition = if_else(grepl("I$", Group), "Infected", "Control")
  ) %>%
  group_by(CDR3.aa) %>%
  summarise(
    public_in_HK = any(Tissue == "HK"),
    public_in_spleen = any(Tissue == "Spleen"),
    public_groups = paste(
      sort(unique(as.character(Group))),
      collapse = "; "
    ),
    tissue_status = case_when(
      public_in_HK & public_in_spleen ~ "Public in both tissues",
      public_in_HK ~ "Public in HK group(s) only",
      public_in_spleen ~ "Public in spleen group(s) only"
    ),
    .groups = "drop"
  )

print(public_tissue_status)

# Sample-level abundance of group-public clonotypes,
# summarized by tissue, condition, and CDR3.
tissue_public <- public_annot %>%
  mutate(
    Tissue = if_else(grepl("^HK", Group), "HK", "Spleen"),
    Condition = if_else(grepl("I$", Group), "Infected", "Control")
  ) %>%
  group_by(Tissue, Condition, CDR3.aa) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Tissue, Condition, desc(total_clones))

############################################################
## 11. Classify public clonotypes by condition
############################################################
## A condition-specific public CDR3 must be public in that
## condition and absent from all samples in the opposite condition.

control_public <- public_by_group %>%
  filter(Group %in% c("HKC", "SpC")) %>%
  pull(CDR3.aa) %>%
  unique()

infected_public <- public_by_group %>%
  filter(Group %in% c("HKI", "SpI")) %>%
  pull(CDR3.aa) %>%
  unique()

all_control_cdr3 <- all_clones_clean %>%
  filter(Condition == "Control") %>%
  pull(CDR3.aa) %>%
  unique()

all_infected_cdr3 <- all_clones_clean %>%
  filter(Condition == "Infected") %>%
  pull(CDR3.aa) %>%
  unique()

infected_specific_public <- setdiff(infected_public, all_control_cdr3)
control_specific_public <- setdiff(control_public, all_infected_cdr3)
shared_public <- intersect(control_public, infected_public)

# Public in one condition but detected, without meeting the public
# criterion, in the other condition.
infected_public_also_detected_in_control <- setdiff(
  setdiff(infected_public, shared_public),
  infected_specific_public
)

control_public_also_detected_in_infected <- setdiff(
  setdiff(control_public, shared_public),
  control_specific_public
)

public_category_summary <- data.frame(
  Category = c(
    "Infection-specific public; absent from controls",
    "Control-specific public; absent from infected fish",
    "Public in both conditions",
    "Infection-public; also detected in controls",
    "Control-public; also detected in infected fish"
  ),
  N_CDR3 = c(
    length(infected_specific_public),
    length(control_specific_public),
    length(shared_public),
    length(infected_public_also_detected_in_control),
    length(control_public_also_detected_in_infected)
  )
)

print(public_category_summary)

ggplot(
  public_category_summary,
  aes(x = reorder(Category, N_CDR3), y = N_CDR3, fill = Category)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N_CDR3), hjust = -0.15) +
  coord_flip() +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Public CDR3 clonotypes by condition",
    x = NULL,
    y = "Number of CDR3 clonotypes"
  )

############################################################
## 12. Private and shared clonotypes
############################################################
## Private = detected in exactly one sample across the balanced
## dataset. This is distinct from the public-within-group rule.

clone_occurrence_all <- all_clones_clean %>%
  group_by(CDR3.aa) %>%
  summarise(
    n_samples_total = n_distinct(Sample_clean),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

private_clones <- clone_occurrence_all %>%
  filter(n_samples_total == 1)

# Identify public CDR3s within each group at the Group-CDR3 level.
group_public_keys <- public_by_group %>%
  select(Group, CDR3.aa) %>%
  distinct()

sample_clone_status <- all_clones_clean %>%
  left_join(clone_occurrence_all, by = "CDR3.aa") %>%
  left_join(
    group_public_keys %>% mutate(is_public_in_group = TRUE),
    by = c("Group", "CDR3.aa")
  ) %>%
  mutate(
    is_public_in_group = coalesce(is_public_in_group, FALSE),
    Clone_type = case_when(
      is_public_in_group ~ "Public within group",
      n_samples_total == 1 ~ "Private to one sample",
      TRUE ~ "Shared across samples; not public within group"
    )
  )

sample_public_private <- sample_clone_status %>%
  group_by(Group, Sample_clean, Clone_type) %>%
  summarise(
    n_CDR3 = n_distinct(CDR3.aa),
    total_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Sample_plot = sub(
      "^([A-Za-z]{3})_([0-9]+).*$",
      "\\1_\\2",
      Sample_clean
    )
  )

ggplot(
  sample_public_private,
  aes(x = Sample_plot, y = n_CDR3, fill = Clone_type)
) +
  geom_col() +
  facet_wrap(~Group, scales = "free_x") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Public, shared, and private CDR3 clonotypes per sample",
    x = "Sample",
    y = "Number of CDR3 clonotypes",
    fill = "Clonotype category"
  )

############################################################
## 13. Family usage among infection-specific public clonotypes
############################################################

infected_specific_table <- all_clones_clean %>%
  filter(
    Condition == "Infected",
    CDR3.aa %in% infected_specific_public
  )

infected_specific_v_family_usage <- infected_specific_table %>%
  filter(!is.na(V_family), V_family != "") %>%
  group_by(Group, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

ggplot(
  infected_specific_v_family_usage,
  aes(x = reorder(V_family, total_clones), y = total_clones, fill = Group)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "V-family usage among infection-specific public clonotypes",
    x = "V family",
    y = "Total clone count"
  )

infected_specific_j_family_usage <- infected_specific_table %>%
  filter(!is.na(J_family), J_family != "") %>%
  group_by(Group, J_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

ggplot(
  infected_specific_j_family_usage,
  aes(x = reorder(J_family, total_clones), y = total_clones, fill = Group)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "J-family usage among infection-specific public clonotypes",
    x = "J family",
    y = "Total clone count"
  )

#############################################
## END of the above script
############################################################

#####################New Script after collapse to family names
############################################################
## IGH PUBLIC/PRIVATE CLONE + V/J FAMILY USAGE ANALYSIS
## Clone definition: unique CDR3.aa
## Relaxed public clone: CDR3.aa shared by >=2 samples in group
############################################################


############################################################
## 1. Metadata
############################################################

meta <- data.frame(Sample = names(imm_ds)) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI"))
  )

table(meta$Group)

############################################################
## 2. Build clone table with V/J information
############################################################

all_clones <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa,
        CDR3.nt,
        Clones,
        Proportion,
        V.name,
        D.name,
        J.name
      ) %>%
      mutate(Sample = s)
  })
)

############################################################
## 3. Clean names and remove allele-level information
############################################################
all_clones_clean <- all_clones %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH$", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    CDR3_length = nchar(CDR3.aa),
    
    # Remove allele designations
    V_gene_no_allele = sub("\\*.*$", "", V.name),
    J_gene_no_allele = sub("\\*.*$", "", J.name),
    
    # Collapse V genes to families, preserving D
    V_family = sub(
      "_?S_?[0-9]+$",
      "",
      sub("-.*$", "", V_gene_no_allele)
    ),
    
    # Collapse J genes to families, preserving a final D
    J_family = sub(
      "-.*$",
      "",
      sub("_?T_?[0-9]+(?=D?$)", "", J_gene_no_allele, perl = TRUE)
    )
  )










#View(all_clones_clean)

############################################################
## 4. Relaxed public clones per group
############################################################

get_public_clones_relaxed <- function(group_name, min_samples = 2) {
  
  all_clones_clean %>%
    filter(Group == group_name) %>%
    group_by(CDR3.aa) %>%
    summarise(
      n_samples = n_distinct(Sample),
      total_clones = sum(Clones, na.rm = TRUE),
      mean_clones = mean(Clones, na.rm = TRUE),
      total_proportion = sum(Proportion, na.rm = TRUE),
      mean_proportion = mean(Proportion, na.rm = TRUE),
      CDR3_length = first(CDR3_length),
      .groups = "drop"
    ) %>%
    filter(n_samples >= min_samples) %>%
    mutate(Group = group_name)
}

HKC_public <- get_public_clones_relaxed("HKC")
HKI_public <- get_public_clones_relaxed("HKI")
SpC_public <- get_public_clones_relaxed("SpC")
SpI_public <- get_public_clones_relaxed("SpI")

public_by_group <- bind_rows(
  HKC_public,
  HKI_public,
  SpC_public,
  SpI_public
)

#View(public_by_group)

############################################################
## 5. Public clone count per group
############################################################

public_counts <- public_by_group %>%
  group_by(Group) %>%
  summarise(
    n_public_clones = n_distinct(CDR3.aa),
    .groups = "drop"
  )

ggplot(public_counts, aes(Group, n_public_clones, fill = Group)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Relaxed Public Clones per Group",
    subtitle = "Public = CDR3.aa shared by at least 2 samples",
    x = "Group",
    y = "Number of public clones"
  )

############################################################
## 6. Annotate public clones with gene/family information
############################################################

public_annot <- public_by_group %>%
  left_join(
    all_clones_clean,
    by = c("CDR3.aa", "Group")
  )

#View(public_annot)

############################################################
## 7. V-family usage among public clones
############################################################

public_v_family <- public_annot %>%
  group_by(Group, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

top15_public_v_family <- public_v_family %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 20, with_ties = FALSE) %>%
  ungroup()

ggplot(top15_public_v_family,
       aes(x = reorder(V_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Top 20 V Families Among Public Clones",
    x = "V family",
    y = "Total clone count"
  )

############################################################
## 8. J-family usage among public clones
############################################################

public_j_family <- public_annot %>%
  group_by(Group, J_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

top15_public_j_family <- public_j_family %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 20, with_ties = FALSE) %>%
  ungroup()

ggplot(top15_public_j_family,
       aes(x = reorder(J_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "J Families Among Public Clones",
    x = "J family",
    y = "Total clone count"
  )

############################################################
## 9. Overall V-family usage across all clones
############################################################

overall_v_family <- all_clones_clean %>%
  group_by(Group, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  )

top15_overall_v_family <- overall_v_family %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 20, with_ties = FALSE) %>%
  ungroup()

ggplot(top15_overall_v_family,
       aes(x = reorder(V_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Top 20 V Families Across All Clones",
    x = "V family",
    y = "Total clone count"
  )

############################################################
## 10. Overall J-family usage across all clones
############################################################

overall_j_family <- all_clones_clean %>%
  group_by(Group, J_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    total_proportion = sum(Proportion, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  )

top15_overall_j_family <- overall_j_family %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 20, with_ties = FALSE) %>%
  ungroup()

ggplot(top15_overall_j_family,
       aes(x = reorder(J_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "J Family usage Across All Clones",
    x = "J family",
    y = "Total clone count"
  )

############################################################
## 11. V-family heatmap for public clones
############################################################

v_family_heatmap_df <- public_v_family %>%
  select(Group, V_family, total_clones) %>%
  pivot_wider(
    names_from = Group,
    values_from = total_clones,
    values_fill = 0
  )

v_family_mat <- v_family_heatmap_df %>%
  column_to_rownames("V_family") %>%
  as.matrix()

pheatmap(
  log10(v_family_mat + 1),
  main = "Public Clone V-Family Usage",
  fontsize_row = 8,
  fontsize_col = 10
)

############################################################
## 12. J-family heatmap for public clones
############################################################

j_family_heatmap_df <- public_j_family %>%
  select(Group, J_family, total_clones) %>%
  pivot_wider(
    names_from = Group,
    values_from = total_clones,
    values_fill = 0
  )

j_family_mat <- j_family_heatmap_df %>%
  column_to_rownames("J_family") %>%
  as.matrix()

pheatmap(
  log10(j_family_mat + 1),
  main = "Public Clone J-Family Usage",
  fontsize_row = 8,
  fontsize_col = 10
)

############################################################
## 13. Public clone CDR3 length
############################################################

public_lengths <- public_by_group %>%
  mutate(CDR3_length = nchar(CDR3.aa))

ggplot(public_lengths,
       aes(x = CDR3_length, fill = Group)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 Length Distribution of Public Clones",
    x = "CDR3 amino acid length",
    y = "Number of public clones"
  )

ggplot(public_lengths,
       aes(x = Group, y = CDR3_length, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CDR3 Length of Public Clones",
    x = "Group",
    y = "CDR3 amino acid length"
  )

############################################################
## 14. Shared and group-specific public clones
############################################################

HK_shared_public <- intersect(HKC_public$CDR3.aa, HKI_public$CDR3.aa)
Sp_shared_public <- intersect(SpC_public$CDR3.aa, SpI_public$CDR3.aa)

HKC_specific_public <- setdiff(HKC_public$CDR3.aa, HKI_public$CDR3.aa)
HKI_specific_public <- setdiff(HKI_public$CDR3.aa, HKC_public$CDR3.aa)

SpC_specific_public <- setdiff(SpC_public$CDR3.aa, SpI_public$CDR3.aa)
SpI_specific_public <- setdiff(SpI_public$CDR3.aa, SpC_public$CDR3.aa)

infected_public <- union(HKI_public$CDR3.aa, SpI_public$CDR3.aa)
control_public  <- union(HKC_public$CDR3.aa, SpC_public$CDR3.aa)

infected_specific_public <- setdiff(infected_public, control_public)
control_specific_public  <- setdiff(control_public, infected_public)

public_category_summary <- data.frame(
  Category = c(
    "HKC public",
    "HKI public",
    "SpC public",
    "SpI public",
    "HKC-HKI shared",
    "SpC-SpI shared",
    "HKC-specific",
    "HKI-specific",
    "SpC-specific",
    "SpI-specific",
    "Infected-specific public",
    "Control-specific public"
  ),
  N_clones = c(
    n_distinct(HKC_public$CDR3.aa),
    n_distinct(HKI_public$CDR3.aa),
    n_distinct(SpC_public$CDR3.aa),
    n_distinct(SpI_public$CDR3.aa),
    length(HK_shared_public),
    length(Sp_shared_public),
    length(HKC_specific_public),
    length(HKI_specific_public),
    length(SpC_specific_public),
    length(SpI_specific_public),
    length(infected_specific_public),
    length(control_specific_public)
  )
)

ggplot(public_category_summary,
       aes(x = reorder(Category, N_clones), y = N_clones)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Public Clone Category Summary",
    x = "",
    y = "Number of clones"
  )


##Color new
############################################################
## Plot public clonotype categories with a complete border
############################################################
public_category_summary <- public_category_summary %>%
  mutate(
    Category = as.character(Category),
    
    ## Classify categories for consistent colouring
    Type = case_when(
      Category %in% c(
        "HKC public",
        "SpC public"
      ) ~ "Control",
      
      Category %in% c(
        "HKI public",
        "SpI public"
      ) ~ "Infected",
      
      Category %in% c(
        "HKC-specific",
        "SpC-specific",
        "Control-specific public"
      ) ~ "Control-specific",
      
      Category %in% c(
        "HKI-specific",
        "SpI-specific",
        "Infected-specific public"
      ) ~ "Infected-specific",
      
      Category %in% c(
        "HKC-HKI shared",
        "SpC-SpI shared"
      ) ~ "Shared",
      
      TRUE ~ "Unclassified"
    ),
    
    ## Factor levels are displayed from bottom to top
    ## after coord_flip()
    Category = factor(
      Category,
      levels = rev(c(
        "Infected-specific public",
        "Control-specific public",
        "SpI public",
        "SpI-specific",
        "SpC public",
        "SpC-specific",
        "HKI public",
        "HKI-specific",
        "HKC public",
        "HKC-specific",
        "SpC-SpI shared",
        "HKC-HKI shared"
      ))
    ),
    
    Type = factor(
      Type,
      levels = c(
        "Control",
        "Control-specific",
        "Infected",
        "Infected-specific",
        "Shared",
        "Unclassified"
      )
    )
  )

## Check category names, counts, and classifications
print(
  public_category_summary %>%
    dplyr::select(Category, N_clones, Type)
)

## Stop if a category name did not match the defined factor levels
if (any(is.na(public_category_summary$Category))) {
  stop(
    "One or more Category names do not match the specified factor levels. ",
    "Check unique(public_category_summary$Category) before factor conversion."
  )
}

############################################################
## Generate figure
############################################################

public_category_plot <- ggplot(
  public_category_summary,
  aes(
    x = Category,
    y = N_clones,
    fill = Type
  )
) +
  geom_col(
    width = 0.85,
    colour = NA
  ) +
  
  ## Add exact values beside bars
  geom_text(
    aes(label = N_clones),
    hjust = -0.15,
    size = 3.8,
    colour = "black"
  ) +
  
  coord_flip(clip = "off") +
  
  scale_fill_manual(
    values = c(
      "Control"           = "#2C6DB2", # dark blue
      "Control-specific"  = "#85B8E3", # light blue
      "Infected"          = "#C43C39", # dark red
      "Infected-specific" = "#F28E85", # light red
      "Shared"            = "#8064A2" # purple
      
    ),
    drop = FALSE
  ) +
  
  ## Leave room for values at the ends of bars
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.14))
  ) +
  
  labs(
    title = "Shared and Group-Specific Public Clonotypes",
    x = NULL,
    y = "Number of public clonotypes",
    fill = "Category type"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 11
    ),
    
    axis.text.x = element_text(
      colour = "black",
      size = 11
    ),
    
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 8)
    ),
    
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    
    ## Box around the plotting panel
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  )

public_category_plot

#end of color


############################################################
## 15. V/J family usage in infected-specific public clones
############################################################

infected_specific_table <- public_annot %>%
  filter(CDR3.aa %in% infected_specific_public)

infected_specific_v_family <- infected_specific_table %>%
  group_by(Group, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

infected_specific_j_family <- infected_specific_table %>%
  group_by(Group, J_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

ggplot(infected_specific_v_family,
       aes(x = reorder(V_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "V-Family Usage of Infected-Specific Public Clones",
    x = "V family",
    y = "Total clone count"
  )

ggplot(infected_specific_j_family,
       aes(x = reorder(J_family, total_clones),
           y = total_clones,
           fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "J-Family Usage of Infected-Specific Public Clones",
    x = "J family",
    y = "Total clone count"
  )

############################################################
## 16. Private clones
## Private = CDR3.aa present in only one sample overall
############################################################

clone_occurrence_all <- all_clones_clean %>%
  group_by(CDR3.aa) %>%
  summarise(
    n_samples_total = n_distinct(Sample),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

private_clones <- clone_occurrence_all %>%
  filter(n_samples_total == 1)

private_annot <- all_clones_clean %>%
  semi_join(private_clones, by = "CDR3.aa")

private_counts <- private_annot %>%
  group_by(Group, Sample) %>%
  summarise(
    n_private_clones = n_distinct(CDR3.aa),
    total_private_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(private_counts,
       aes(x = Group, y = n_private_clones, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "Private Clone Counts per Sample",
    x = "Group",
    y = "Number of private clones"
  )

############################################################
## 17. Public/private summary per sample
############################################################

sample_public_private <- all_clones_clean %>%
  mutate(
    Clone_type = ifelse(
      CDR3.aa %in% private_clones$CDR3.aa,
      "Private",
      "Shared/Public"
    )
  ) %>%
  group_by(Group, Sample, Clone_type) %>%
  summarise(
    n_clones = n_distinct(CDR3.aa),
    total_abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(sample_public_private,
       aes(x = Sample, y = n_clones, fill = Clone_type)) +
  geom_col() +
  facet_wrap(~ Group, scales = "free_x") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(
    title = "Private vs Shared/Public Clones per Sample",
    x = "Sample",
    y = "Number of clones"
  )

################New script replacing 2193-3917 and compatible with the
##new rules and updates

expected_groups <- c("HKC", "HKI", "SpC", "SpI")

# Confirm that the QC-filtered objects from the previous script exist
required_objects <- c("all_clones_clean", "public_by_group", "public_annot")

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Run the previous QC script first. Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

required_columns <- c(
  "Sample_clean", "Group", "Clones", "Proportion",
  "V.name", "J.name", "V_family", "J_family"
)

missing_columns <- setdiff(required_columns, names(all_clones_clean))

if (length(missing_columns) > 0) {
  stop(
    "`all_clones_clean` is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Confirm that the suffixes have been removed from the family columns
unremoved_v_suffix <- grepl(
  "S_?[0-9]+$",
  all_clones_clean$V_family,
  ignore.case = TRUE
)

unremoved_j_suffix <- grepl(
  "T_?[0-9]+D?$",
  all_clones_clean$J_family,
  ignore.case = TRUE
)

if (any(unremoved_v_suffix, na.rm = TRUE)) {
  stop("A V_family still contains a trailing S-number suffix.")
}

if (any(unremoved_j_suffix, na.rm = TRUE)) {
  stop("A J_family still contains a trailing T-number suffix.")
}

# One row per biological sample and group
sample_groups <- all_clones_clean %>%
  distinct(Sample_clean, Group) %>%
  mutate(Group = as.character(Group))

# Ensure every family summary includes zeroes for families not observed
# in a particular sample. This permits fair group-level averages.
summarise_family_per_sample <- function(data, family_column, sample_groups) {
  
  valid_data <- data %>%
    filter(
      !is.na(.data[[family_column]]),
      .data[[family_column]] != "",
      !is.na(Clones),
      is.finite(Clones),
      !is.na(Proportion)
    ) %>%
    transmute(
      Sample_clean,
      Group = as.character(Group),
      Family = .data[[family_column]],
      Clones,
      Proportion
    )
  
  observed <- valid_data %>%
    group_by(Sample_clean, Group, Family) %>%
    summarise(
      total_clones = sum(Clones, na.rm = TRUE),
      repertoire_proportion = sum(Proportion, na.rm = TRUE),
      .groups = "drop"
    )
  
  family_universe <- sort(unique(valid_data$Family))
  
  if (length(family_universe) == 0) {
    stop("No valid family names were found in `", family_column, "`.")
  }
  
  tidyr::crossing(
    sample_groups,
    Family = family_universe
  ) %>%
    left_join(
      observed,
      by = c("Sample_clean", "Group", "Family")
    ) %>%
    mutate(
      total_clones = replace_na(total_clones, 0),
      repertoire_proportion = replace_na(repertoire_proportion, 0)
    )
}

# Summarise the three biological samples within each group
summarise_family_by_group <- function(sample_data) {
  sample_data %>%
    group_by(Group, Family) %>%
    summarise(
      mean_sample_clones = mean(total_clones),
      mean_sample_proportion = mean(repertoire_proportion),
      prevalence = mean(total_clones > 0),
      .groups = "drop"
    )
}

############################################################
## A. Overall V- and J-family usage across all QC-passing clonotypes
############################################################

overall_v_sample <- summarise_family_per_sample(
  all_clones_clean,
  "V_family",
  sample_groups
)

overall_j_sample <- summarise_family_per_sample(
  all_clones_clean,
  "J_family",
  sample_groups
)

overall_v_family <- summarise_family_by_group(overall_v_sample)
overall_j_family <- summarise_family_by_group(overall_j_sample)

# Select the 20 most abundant V families per group by mean sample proportion
top20_overall_v <- overall_v_family %>%
  group_by(Group) %>%
  slice_max(
    order_by = mean_sample_proportion,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top20_overall_v,
  aes(
    x = reorder(Family, mean_sample_proportion),
    y = 100 * mean_sample_proportion,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top V-family usage across all clonotypes",
    x = "V family",
    y = "Mean repertoire proportion across samples (%)"
  )

# Select the 20 most abundant J families per group by mean sample proportion
top20_overall_j <- overall_j_family %>%
  group_by(Group) %>%
  slice_max(
    order_by = mean_sample_proportion,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top20_overall_j,
  aes(
    x = reorder(Family, mean_sample_proportion),
    y = 100 * mean_sample_proportion,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top J-family usage across all clonotypes",
    x = "J family",
    y = "Mean repertoire proportion across samples (%)"
  )

############################################################
## B. V- and J-family usage among public clonotypes
############################################################

public_v_sample <- summarise_family_per_sample(
  public_annot,
  "V_family",
  sample_groups
)

public_j_sample <- summarise_family_per_sample(
  public_annot,
  "J_family",
  sample_groups
)

public_v_family <- summarise_family_by_group(public_v_sample)
public_j_family <- summarise_family_by_group(public_j_sample)

top20_public_v <- public_v_family %>%
  group_by(Group) %>%
  slice_max(
    order_by = mean_sample_proportion,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top20_public_v,
  aes(
    x = reorder(Family, mean_sample_proportion),
    y = 100 * mean_sample_proportion,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top V-family usage among public clonotypes",
    x = "V family",
    y = "Mean repertoire proportion across samples (%)"
  )

top20_public_j <- public_j_family %>%
  group_by(Group) %>%
  slice_max(
    order_by = mean_sample_proportion,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()

ggplot(
  top20_public_j,
  aes(
    x = reorder(Family, mean_sample_proportion),
    y = 100 * mean_sample_proportion,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Top J-family usage among public clonotypes",
    x = "J family",
    y = "Mean repertoire proportion across samples (%)"
  )

############################################################
## C. Family-usage heatmaps
############################################################

plot_family_heatmap <- function(summary_data, plot_title) {
  
  wide_data <- summary_data %>%
    mutate(Group = as.character(Group)) %>%
    select(Family, Group, mean_sample_proportion) %>%
    pivot_wider(
      names_from = Group,
      values_from = mean_sample_proportion,
      values_fill = 0
    )
  
  # Add a zero column if a group is absent from the data
  missing_groups <- setdiff(expected_groups, names(wide_data))
  
  for (group_name in missing_groups) {
    wide_data[[group_name]] <- 0
  }
  
  wide_data <- wide_data %>%
    select(Family, all_of(expected_groups))
  
  heatmap_matrix <- as.matrix(
    wide_data[, expected_groups, drop = FALSE]
  )
  storage.mode(heatmap_matrix) <- "numeric"
  rownames(heatmap_matrix) <- wide_data$Family
  
  pheatmap(
    log10(100 * heatmap_matrix + 1),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = plot_title
  )
}

plot_family_heatmap(
  public_v_family,
  "Public-clonotype V-family usage"
)

plot_family_heatmap(
  public_j_family,
  "Public-clonotype J-family usage"
)

plot_family_heatmap(
  overall_v_family,
  "Overall V-family usage"
)

plot_family_heatmap(
  overall_j_family,
  "Overall J-family usage"
)

############################################################
## D. Preserve raw gene-level usage tables for future analyses
############################################################

public_v_gene_usage <- public_annot %>%
  filter(!is.na(V.name), V.name != "") %>%
  group_by(Group, V.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

public_j_gene_usage <- public_annot %>%
  filter(!is.na(J.name), J.name != "") %>%
  group_by(Group, J.name) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_public_CDR3 = n_distinct(CDR3.aa),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(total_clones))

#####################################end of the script

#############################New analysis V-J pairing
############################################################
## 1. Define and validate the 12 samples
############################################################
shm_dir <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/SHM_mutation_spectrum_HK_Sp_3w"

region_rs_file <- file.path(shm_dir, "region_rs_SAV_3wpi.rds")
sample_region_rs_file <- file.path(shm_dir, "sample_region_rs_SAV_3wpi.csv")
sample_region_class_rs_file <- file.path(
  shm_dir,
  "sample_region_class_rs_SAV_3wpi.csv"
)

required_files <- c(
  region_rs_file,
  sample_region_rs_file,
  sample_region_class_rs_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "These files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

region_rs <- readRDS(region_rs_file)

sample_region_rs <- read.csv(
  sample_region_rs_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_region_class_rs <- read.csv(
  sample_region_class_rs_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check that the loaded data contain the expected groups and regions.
print(table(region_rs$Group, useNA = "ifany"))
print(table(region_rs$Region, useNA = "ifany"))

expected_samples <- c(
  "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
  "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
  "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
  "SpI_1I3w", "SpI_2I3w", "SpI_5I3w"
)

expected_groups <- c("HKC", "HKI", "SpC", "SpI")

if (!exists("region_rs")) {
  stop("The object `region_rs` is not present in the R environment.")
}

required_columns <- c("Sample", "CDR3.aa", "Clones")
missing_columns <- setdiff(required_columns, names(region_rs))

if (length(missing_columns) > 0) {
  stop(
    "region_rs is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!("V.name" %in% names(region_rs) || "V_family" %in% names(region_rs))) {
  stop("region_rs must contain either `V.name` or `V_family`.")
}

if (!("J.name" %in% names(region_rs) || "J_family" %in% names(region_rs))) {
  stop("region_rs must contain either `J.name` or `J_family`.")
}

sample_metadata <- tibble(Sample_clean = expected_samples) %>%
  mutate(
    Group = sub("_.*$", "", Sample_clean),
    Tissue = if_else(startsWith(Group, "HK"), "HK", "Spleen"),
    Condition = if_else(endsWith(Group, "I"), "Infected", "Control"),
    Group = factor(Group, levels = expected_groups),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected")),
    Sample = paste0(Sample_clean, ".clones_IGH")
  )

print(sample_metadata, n = Inf)

if (any(table(sample_metadata$Group) != 3L)) {
  stop("Expected exactly three samples in each group.")
}

############################################################
## 2. Clean V- and J-family names
############################################################

clean_v_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\(.*\\)$", "", x)          # Remove trailing annotation
  x <- sub("\\*.*$", "", x)             # Remove allele
  x <- sub("-.*$", "", x)               # Collapse gene-level suffix
  x <- sub("S_?[0-9]+$", "", x)          # Remove terminal S and number
  trimws(x)
}

clean_j_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\(.*\\)$", "", x)          # Remove trailing annotation
  x <- sub("\\*.*$", "", x)             # Remove allele
  x <- sub("T_?[0-9]+(?=D?$)", "", x, perl = TRUE)  # Remove T and number; keep D
  x <- sub("-.*$", "", x)               # Collapse gene-level suffix
  trimws(x)
}

############################################################
## 3. Standardize sample names and apply clonotype QC
############################################################
vj_input <- region_rs %>%
  mutate(
    Sample_original = as.character(Sample),
    Sample_clean = sub("\\.clones_IGH$", "", Sample_original)
  )

# Report samples present in region_rs but excluded from the balanced analysis.
excluded_samples <- setdiff(
  unique(vj_input$Sample_clean),
  sample_metadata$Sample_clean
)

message(
  "Samples excluded from this balanced analysis: ",
  if (length(excluded_samples) == 0) "none"
  else paste(excluded_samples, collapse = ", ")
)

# Retain only the 12 samples defined in sample_metadata.
vj_input <- vj_input %>%
  filter(Sample_clean %in% sample_metadata$Sample_clean)

# Confirm that all selected samples are present.
missing_samples <- setdiff(
  sample_metadata$Sample_clean,
  unique(vj_input$Sample_clean)
)

if (length(missing_samples) > 0) {
  stop(
    "Selected samples missing from region_rs: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Prefer raw gene names when available; otherwise clean the existing family names.
v_source <- if ("V.name" %in% names(vj_input)) {
  vj_input$V.name
} else {
  vj_input$V_family
}

j_source <- if ("J.name" %in% names(vj_input)) {
  vj_input$J.name
} else {
  vj_input$J_family
}

vj_input$V_family <- clean_v_family(v_source)
vj_input$J_family <- clean_j_family(j_source)

# Replace any existing metadata with metadata derived from the validated sample IDs.
vj_input <- vj_input %>%
  select(-any_of(c("Group", "Tissue", "Condition", "Sample"))) %>%
  left_join(sample_metadata, by = "Sample_clean")

# Convert counts and CDR3 amino-acid sequences to standard formats.
vj_input <- vj_input %>%
  mutate(
    CDR3.aa = toupper(as.character(CDR3.aa)),
    Clones = suppressWarnings(as.numeric(as.character(Clones)))
  )

rows_before_qc <- nrow(vj_input)

# Retain standard amino-acid CDR3 sequences with the productive C and W/F anchors,
# assigned V/J families, and valid positive MiXCR counts.
vj_clonotypes <- vj_input %>%
  filter(
    !is.na(CDR3.aa),
    grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", CDR3.aa),
    grepl("^C[ACDEFGHIKLMNPQRSTVWY]*[WF]$", CDR3.aa),
    !is.na(V_family), V_family != "",
    !is.na(J_family), J_family != "",
    !is.na(Clones), is.finite(Clones), Clones > 0
  )

# Remove repeated rows arising from regional annotations.
# Keep CDR3.nt in the key when it is available.
dedup_columns <- c(
  "Sample_clean", "Group", "Tissue", "Condition",
  "CDR3.aa", intersect("CDR3.nt", names(vj_clonotypes)),
  "V_family", "J_family", "Clones"
)

vj_clonotypes <- vj_clonotypes %>%
  distinct(across(all_of(dedup_columns))) %>%
  mutate(
    VJ_pair = paste(V_family, J_family, sep = "__")
  )

qc_summary <- tibble(
  rows_before_qc = rows_before_qc,
  rows_after_qc_and_deduplication = nrow(vj_clonotypes),
  distinct_samples_retained = n_distinct(vj_clonotypes$Sample_clean),
  distinct_vj_pairs = n_distinct(vj_clonotypes$VJ_pair)
)

print(qc_summary)

if (n_distinct(vj_clonotypes$Sample_clean) != 12) {
  stop("QC did not retain clonotypes from all 12 samples.")
}

# Check that the family names do not retain the suffixes discussed previously.
if (any(grepl("S_?[0-9]+$", vj_clonotypes$V_family))) {
  stop("Some V-family names still end in S followed by a number.")
}

if (any(grepl("T_?[0-9]+(?=D?$)", vj_clonotypes$J_family, perl = TRUE))) {
  stop("Some J-family names still contain T followed by a number.")
}

############################################################
## 4. Calculate sample-level V-J pair proportions
############################################################

# Each row is a sample-level clonotype assignment to a V-J family pair.
# If CDR3.nt is available, distinct nucleotide clonotypes are retained.
vj_observed <- vj_clonotypes %>%
  group_by(
    Sample_clean, Group, Tissue, Condition,
    V_family, J_family, VJ_pair
  ) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

vj_pair_universe <- vj_observed %>%
  distinct(V_family, J_family, VJ_pair)

# Complete the sample-by-pair grid so that absent pairs are represented by zero.
vj_sample <- crossing(
  sample_metadata,
  vj_pair_universe
) %>%
  left_join(
    vj_observed %>%
      select(
        Sample_clean, V_family, J_family, VJ_pair, total_clones
      ),
    by = c("Sample_clean", "V_family", "J_family", "VJ_pair")
  ) %>%
  mutate(total_clones = replace_na(total_clones, 0)) %>%
  group_by(Sample_clean) %>%
  mutate(
    sample_total_assigned_clones = sum(total_clones, na.rm = TRUE),
    VJ_fraction = if_else(
      sample_total_assigned_clones > 0,
      total_clones / sample_total_assigned_clones,
      NA_real_
    )
  ) %>%
  ungroup()

# Confirm that pair proportions sum to 1 within each sample.
sample_fraction_check <- vj_sample %>%
  group_by(Sample_clean, Group) %>%
  summarise(
    total_assigned_clones = first(sample_total_assigned_clones),
    sum_VJ_fraction = sum(VJ_fraction, na.rm = TRUE),
    .groups = "drop"
  )

print(sample_fraction_check, n = Inf)

############################################################
## 5. Summarize V-J pair use by group
############################################################

vj_group <- vj_sample %>%
  group_by(Group, Tissue, Condition, V_family, J_family, VJ_pair) %>%
  summarise(
    mean_fraction = mean(VJ_fraction, na.rm = TRUE),
    sd_fraction = sd(VJ_fraction, na.rm = TRUE),
    median_fraction = median(VJ_fraction, na.rm = TRUE),
    n_samples_present = sum(total_clones > 0),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  )

print(vj_group, n = 30)

# Display the 20 most represented pairs in each group.
top20_vj_by_group <- vj_group %>%
  group_by(Group) %>%
  slice_max(mean_fraction, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Group, desc(mean_fraction))

print(top20_vj_by_group, n = Inf)

############################################################
## 6. Plot the most represented pairs
############################################################

ggplot(
  top20_vj_by_group,
  aes(
    x = reorder(VJ_pair, mean_fraction),
    y = 100 * mean_fraction,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Group, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Most represented V-J family pairs",
    x = "V-J family pair",
    y = "Mean sample proportion of assigned clones (%)"
  )

############################################################
## 7. Display a heatmap for each group
############################################################
library(grid)

draw_vj_heatmap <- function(group_name) {
  heatmap_data <- vj_group %>%
    filter(Group == group_name) %>%
    select(VJ_pair, J_family, mean_fraction) %>%
    pivot_wider(
      names_from = J_family,
      values_from = mean_fraction,
      values_fill = 0
    )
  
  mat <- as.matrix(heatmap_data[, -1, drop = FALSE])
  rownames(mat) <- heatmap_data$VJ_pair
  
  hm <- pheatmap(
    mat,
    cluster_rows = nrow(mat) > 1,
    cluster_cols = ncol(mat) > 1,
    main = paste("V-J family usage:", group_name),
    angle_col = 45,
    silent = TRUE
  )
  
  grid.newpage()
  grid.draw(hm$gtable)
}

for (g in expected_groups) {
  draw_vj_heatmap(g)
}

############################################################
## 8. Combined heatmap of the most represented V-J pairs
############################################################

top_pairs <- vj_group %>%
  group_by(VJ_pair) %>%
  summarise(overall_mean_fraction = mean(mean_fraction), .groups = "drop") %>%
  slice_max(overall_mean_fraction, n = 40, with_ties = FALSE) %>%
  pull(VJ_pair)

combined_heatmap_data <- vj_group %>%
  filter(VJ_pair %in% top_pairs) %>%
  select(VJ_pair, Group, mean_fraction) %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_fraction,
    values_fill = 0
  ) %>%
  select(VJ_pair, all_of(expected_groups))

combined_mat <- as.matrix(combined_heatmap_data[, -1, drop = FALSE])
rownames(combined_mat) <- combined_heatmap_data$VJ_pair

# Row scaling highlights relative group patterns for each pair.
combined_mat_scaled <- t(scale(t(combined_mat)))
combined_mat_scaled[!is.finite(combined_mat_scaled)] <- 0

combined_hm <- pheatmap(
  combined_mat_scaled,
  cluster_rows = nrow(combined_mat_scaled) > 1,
  cluster_cols = FALSE,
  main = "Relative V-J family-pair representation across groups",
  silent = TRUE
)

grid.newpage()
grid.draw(combined_hm$gtable)

############################################################
## 9. Exploratory infected-versus-control comparisons
############################################################

safe_wilcox_p <- function(control_values, infected_values) {
  tryCatch(
    wilcox.test(
      control_values,
      infected_values,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

vj_infection_tests <- vj_sample %>%
  group_by(Tissue, V_family, J_family, VJ_pair) %>%
  summarise(
    pair_total_clones = sum(total_clones, na.rm = TRUE),
    
    n_control_present = sum(
      total_clones[Condition == "Control"] > 0,
      na.rm = TRUE
    ),
    
    n_infected_present = sum(
      total_clones[Condition == "Infected"] > 0,
      na.rm = TRUE
    ),
    
    mean_control_fraction = mean(
      VJ_fraction[Condition == "Control"],
      na.rm = TRUE
    ),
    
    mean_infected_fraction = mean(
      VJ_fraction[Condition == "Infected"],
      na.rm = TRUE
    ),
    
    p_value = if (sum(total_clones, na.rm = TRUE) > 0) {
      safe_wilcox_p(
        VJ_fraction[Condition == "Control"],
        VJ_fraction[Condition == "Infected"]
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  group_by(Tissue) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    log2_fold_change = log2(
      (mean_infected_fraction + 1e-6) /
        (mean_control_fraction + 1e-6)
    )
  ) %>%
  ungroup() %>%
  arrange(p_adj_BH)

print(vj_infection_tests, n = Inf)

###
# Number of pairs tested and significant in each tissue
vj_infection_tests %>%
  group_by(Tissue) %>%
  summarise(
    pairs_in_table = n(),
    pairs_tested = sum(!is.na(p_value)),
    pairs_significant_BH = sum(p_adj_BH < 0.05, na.rm = TRUE),
    smallest_adjusted_p = if (all(is.na(p_adj_BH))) {
      NA_real_
    } else {
      min(p_adj_BH, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  print(n = Inf)

# Display the test results and adjusted p-values for every pair
vj_infection_tests %>%
  select(
    Tissue, VJ_pair,
    n_control_present, n_infected_present,
    mean_control_fraction, mean_infected_fraction,
    p_value, p_adj_BH
  ) %>%
  arrange(Tissue, p_adj_BH) %>%
  print(n = Inf, width = Inf)

# Display pairs significant after adjustment, if any
vj_infection_tests %>%
  filter(!is.na(p_adj_BH), p_adj_BH < 0.05) %>%
  select(
    Tissue, VJ_pair,
    mean_control_fraction, mean_infected_fraction,
    p_value, p_adj_BH
  ) %>%
  print(n = Inf, width = Inf)



# Descriptive candidates: present in at least 2 infected samples and
# with a higher mean fraction in infected than control fish.
# This is a prioritization rule, not evidence of antigen specificity.
infection_enriched_candidates <- vj_infection_tests %>%
  filter(
    n_infected_present >= 2,
    mean_infected_fraction > mean_control_fraction
  ) %>%
  arrange(desc(log2_fold_change))

print(infection_enriched_candidates, n = Inf)

#### figure
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# Define descriptive candidates and labels.
# This rule does not imply statistical significance.
vj_plot_data <- vj_infection_tests %>%
  filter(
    !is.na(mean_control_fraction),
    !is.na(mean_infected_fraction),
    !is.na(log2_fold_change)
  ) %>%
  mutate(
    Candidate = (
      n_infected_present >= 2 &
        mean_infected_fraction > mean_control_fraction
    ),
    Pair_label = paste(V_family, J_family, sep = "_")
  )

# Label the 12 candidates with the highest mean infected proportion
# separately within each tissue.
vj_labels <- vj_plot_data %>%
  filter(Candidate) %>%
  group_by(Tissue) %>%
  slice_max(
    order_by = mean_infected_fraction,
    n = 12,
    with_ties = FALSE
  ) %>%
  ungroup()

# Plot
p_vj_candidates <- ggplot(
  vj_plot_data,
  aes(
    x = log2_fold_change,
    y = mean_infected_fraction,
    color = Candidate
  )
) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_text_repel(
    data = vj_labels,
    aes(label = Pair_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    show.legend = FALSE,
    seed = 1
  ) +
  facet_wrap(~Tissue, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE" = "#E66A13"
    ),
    breaks = c(FALSE, TRUE),
    labels = c("Other V-J pair", "Exploratory candidate"),
    name = NULL
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    title = "Exploratory V-J pairing patterns",
    subtitle = paste(
      "Labels show the 12 candidates with the highest mean infected-sample",
      "proportion in each tissue; no pair was significant after BH correction"
    ),
    x = expression(log[2]~"fold change (infected/control)"),
    y = "Mean frequency in infected samples"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey88"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.3),
    plot.title = element_text(face = "bold")
  )

print(p_vj_candidates)

###############Infection associated pairings
library(dplyr)
library(ggplot2)

# Reapply the agreed family-name cleaning rules to ensure labels are clean.
clean_v_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\*.*$", "", x)                  # Remove allele
  x <- sub("-.*$", "", x)                    # Collapse gene suffix
  x <- sub("S_?[0-9]+$", "", x)              # Remove terminal S and number
  trimws(x)
}

clean_j_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\*.*$", "", x)                  # Remove allele
  x <- sub("T_?[0-9]+(?=D?$)", "", x, perl = TRUE) # Remove T + number; retain D
  x <- sub("-.*$", "", x)                    # Collapse gene suffix
  trimws(x)
}

# Select descriptive infection-enriched candidates.
vj_candidates_plot <- vj_infection_tests %>%
  mutate(
    V_family = clean_v_family(V_family),
    J_family = clean_j_family(J_family),
    VJ_label = paste(V_family, J_family, sep = "_")
  ) %>%
  filter(
    n_infected_present >= 2,
    mean_infected_fraction > mean_control_fraction,
    !is.na(log2_fold_change)
  ) %>%
  group_by(Tissue) %>%
  slice_max(
    order_by = log2_fold_change,
    n = 15,
    with_ties = FALSE
  ) %>%
  arrange(log2_fold_change, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    # Make pair labels unique across facets while retaining their order.
    plot_label = paste(Tissue, VJ_label, sep = "___"),
    plot_label = factor(plot_label, levels = unique(plot_label))
  )

# Display the pairs used in the plot.
vj_candidates_plot %>%
  select(
    Tissue, VJ_label, n_infected_present,
    mean_control_fraction, mean_infected_fraction, log2_fold_change
  ) %>%
  print(n = Inf, width = Inf)

# Plot
p_vj_candidates_bar <- ggplot(
  vj_candidates_plot,
  aes(x = log2_fold_change, y = plot_label, fill = Tissue)
) +
  geom_col(width = 0.82) +
  facet_wrap(~Tissue, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(
    values = c("HK" = "#D95F02", "Spleen" = "#0878B6"),
    guide = "none"
  ) +
  labs(
    title = "V???J pairs meeting descriptive infection-enrichment criteria",
    subtitle = paste(
      "Higher mean proportion in infected fish and detected in at least",
      "two of three infected samples; no pair was significant after BH correction"
    ),
    x = expression(log[2]~"fold change (infected/control)"),
    y = "V???J pair"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey88"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_vj_candidates_bar)

########## V-J-CDR3_length classes
library(dplyr)
library(ggplot2)

# Reapply family-name cleaning to ensure the S-number and T-number suffixes
# are removed while retaining D.
clean_v_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\*.*$", "", x)
  x <- sub("-.*$", "", x)
  x <- sub("S_?[0-9]+$", "", x)
  trimws(x)
}

clean_j_family <- function(x) {
  x <- as.character(x)
  x <- sub("\\*.*$", "", x)
  x <- sub("T_?[0-9]+(?=D?$)", "", x, perl = TRUE)
  x <- sub("-.*$", "", x)
  trimws(x)
}

top_configuration_classes <- vj_clonotypes %>%
  mutate(
    V_family = clean_v_family(V_family),
    J_family = clean_j_family(J_family),
    CDR3_length = nchar(CDR3.aa),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI"))
  ) %>%
  group_by(Group, V_family, J_family, CDR3_length) %>%
  summarise(
    total_clone_count = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  slice_max(
    order_by = total_clone_count,
    n = 20,
    with_ties = FALSE
  ) %>%
  arrange(Group, total_clone_count) %>%
  ungroup() %>%
  mutate(
    configuration = paste(V_family, J_family, CDR3_length, sep = "_"),
    plot_label = paste(Group, configuration, sep = "___"),
    plot_label = factor(plot_label, levels = unique(plot_label))
  )

# Display the values plotted
top_configuration_classes %>%
  select(Group, configuration, total_clone_count) %>%
  print(n = Inf)

# Plot
p_vj_cdr3_classes <- ggplot(
  top_configuration_classes,
  aes(x = total_clone_count, y = plot_label, fill = Group)
) +
  geom_col(width = 0.82) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(
    values = c(
      "HKC" = "#F8766D",
      "HKI" = "#7CAE00",
      "SpC" = "#00BFC4",
      "SpI" = "#C77CFF"
    ),
    guide = "none"
  ) +
  labs(
    title = "Top 20 V-J-CDR3-length configuration classes",
    x = "Total clone count",
    y = "V-J family-CDR3-length class"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_vj_cdr3_classes)

###CDR3-similarity network
library(stringdist)
library(igraph)
library(ggraph)
# One node per unique CDR3 amino-acid sequence.
# Sequences occurring in multiple samples are represented by one node.
node_data <- vj_clonotypes %>%
  filter(
    !is.na(CDR3.aa),
    CDR3.aa != ""
  ) %>%
  group_by(CDR3.aa) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_samples = n_distinct(Sample_clean),
    .groups = "drop"
  ) %>%
  rename(name = CDR3.aa)

seqs <- node_data$name

if (length(seqs) < 2) {
  stop("At least two unique CDR3 amino-acid sequences are needed.")
}

# Calculate pairwise Levenshtein distances and retain pairs
# differing by one or two amino-acid edits.
distance_matrix <- stringdist::stringdistmatrix(
  seqs,
  seqs,
  method = "lv"
)

edge_indices <- which(
  distance_matrix >= 1 &
    distance_matrix <= 2 &
    upper.tri(distance_matrix),
  arr.ind = TRUE
)

edge_data <- data.frame(
  from = seqs[edge_indices[, 1]],
  to = seqs[edge_indices[, 2]],
  distance = distance_matrix[edge_indices]
)

# Include all sequences as nodes, including sequences with no qualifying edges.
cdr3_graph <- igraph::graph_from_data_frame(
  d = edge_data,
  directed = FALSE,
  vertices = node_data
)

network_components <- igraph::components(cdr3_graph)

message("Unique CDR3 sequences: ", igraph::vcount(cdr3_graph))
message("Sequence-similarity edges: ", igraph::ecount(cdr3_graph))
message("Connected components: ", network_components$no)
message(
  "Largest component: ",
  max(network_components$csize),
  " sequences"
)

# Draw the network.
set.seed(123)

p_cdr3_network <- ggraph(cdr3_graph, layout = "fr") +
  geom_edge_link(
    colour = "grey65",
    alpha = 0.55,
    linewidth = 0.35
  ) +
  geom_node_point(
    shape = 21,
    size = 2.8,
    fill = "#E6A400",
    colour = "black",
    stroke = 0.35
  ) +
  theme_void(base_size = 14) +
  labs(title = "CDR3 amino-acid similarity network") +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

print(p_cdr3_network)





























#####################################Nex high level analysis
############################################################
## ADVANCED IGH REPERTOIRE ANALYSIS
## Groups: HKC, HKI, SpC, SpI
## Clone definition: CDR3.aa
############################################################

library(vegan)
library(ggrepel)
library(stringdist)
library(igraph)

############################################################
## 1. Build metadata
############################################################

meta <- data.frame(Sample = names(imm_ds)) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI"))
  )

table(meta$Group)

############################################################
## 2. Build full clone table
############################################################

all_clones <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa, CDR3.nt, Clones, Proportion,
        V.name, D.name, J.name, Sequence,
        CDR1.aa, CDR2.aa, FR1.aa, FR2.aa, FR3.aa, FR4.aa
      ) %>%
      mutate(Sample = s)
  })
)

all_clones <- all_clones %>%
  mutate(
    Sample_clean = sub("\\.clones_IGH$", "", Sample),
    Group = sub("_.*$", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    CDR3_length = nchar(CDR3.aa),
    
    # Remove parenthetical annotations and allele designations
    V_gene_no_allele = sub("\\*.*$", "", sub("\\(.*\\)", "", V.name)),
    J_gene_no_allele = sub("\\*.*$", "", sub("\\(.*\\)", "", J.name)),
    
    # Collapse V genes to families, preserving D
    V_family = sub(
      "_?S_?[0-9]+$",
      "",
      sub("-.*$", "", V_gene_no_allele)
    ),
    
    # Remove terminal T plus digits, preserving a final D; then collapse at hyphen
    J_family = sub(
      "-.*$",
      "",
      sub("_?-?T_?[0-9]+(?=D?$)", "", J_gene_no_allele, perl = TRUE)
    )
  )

############################################################
## 3. Differential clone expansion
############################################################

clone_group_abund <- all_clones %>%
  group_by(Group, CDR3.aa) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    mean_prop = mean(Proportion, na.rm = TRUE),
    V_family = first(V_family),
    J_family = first(J_family),
    CDR3_length = first(CDR3_length),
    .groups = "drop"
  )

clone_wide <- clone_group_abund %>%
  select(Group, CDR3.aa, total_clones) %>%
  pivot_wider(
    names_from = Group,
    values_from = total_clones,
    values_fill = 0
  )

diff_expansion <- clone_wide %>%
  mutate(
    log2FC_HKI_vs_HKC = log2((HKI + 1) / (HKC + 1)),
    log2FC_SpI_vs_SpC = log2((SpI + 1) / (SpC + 1))
  )

top_HK_expanded <- diff_expansion %>%
  arrange(desc(log2FC_HKI_vs_HKC)) %>%
  slice_head(n = 50)

top_Sp_expanded <- diff_expansion %>%
  arrange(desc(log2FC_SpI_vs_SpC)) %>%
  slice_head(n = 50)

#View(top_HK_expanded)
#View(top_Sp_expanded)

############################################################
## 10. CDR3 motif analysis: k-mer enrichment
############################################################

get_kmers <- function(seq, k = 3) {
  if (is.na(seq) || nchar(seq) < k) return(NA)
  sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
}

kmer_table <- all_clones %>%
  filter(!is.na(CDR3.aa)) %>%
  rowwise() %>%
  mutate(Kmer = list(get_kmers(CDR3.aa, k = 3))) %>%
  unnest(Kmer) %>%
  ungroup() %>%
  filter(!is.na(Kmer))

##
kmer_table <- kmer_table %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI"))
  )

colnames(kmer_table)
table(kmer_table$Group)

##
kmer_group <- kmer_table %>%
  group_by(Group, Kmer) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  )


####
## Add Group to your current all_clones object
all_clones <- all_clones %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean)
  )

table(all_clones$Group)


top_kmers <- kmer_group %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 25, with_ties = FALSE) %>%
  ungroup()

ggplot(top_kmers,
       aes(reorder(Kmer, total_clones), total_clones, fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Top 25 CDR3 3-mer Motifs",
    x = "3-mer motif",
    y = "Abundance-weighted motif count"
  )

####New analysis on 3-aa motiffs
########################################################
############################################################
## CDR3 amino-acid motif analysis
## Tests selective remodeling within a preserved motif landscape
############################################################

## Required packages
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(vegan)

set.seed(1234)

############################################################
## 1. Construct clonotype-level dataset
############################################################

motif_input <- immdata$data %>%
  purrr::imap_dfr(
    ~ .x %>%
      dplyr::select(
        CDR3.aa,
        Clones
      ) %>%
      dplyr::mutate(Sample = .y)
  ) %>%
  dplyr::mutate(
    Sample = gsub("\\.clones_IGH$", "", Sample),
    
    Group = sub("_.*", "", Sample),
    
    Tissue = dplyr::case_when(
      grepl("^HK", Group) ~ "HK",
      grepl("^Sp", Group) ~ "Spleen",
      TRUE ~ NA_character_
    ),
    
    Condition = dplyr::case_when(
      Group %in% c("HKC", "SpC") ~ "Control",
      Group %in% c("HKI", "SpI") ~ "Infected",
      TRUE ~ NA_character_
    ),
    
    CDR3.aa = as.character(CDR3.aa),
    Clones = as.numeric(Clones)
  ) %>%
  dplyr::filter(
    !is.na(CDR3.aa),
    CDR3.aa != "",
    !grepl("\\*", CDR3.aa),
    !is.na(Tissue),
    !is.na(Condition),
    !is.na(Clones),
    Clones > 0
  )

## Check sample metadata
sample_metadata <- motif_input %>%
  dplyr::distinct(Sample, Group, Tissue, Condition) %>%
  dplyr::arrange(Tissue, Condition, Sample)

print(sample_metadata)
table(sample_metadata$Tissue, sample_metadata$Condition)

############################################################
## Optional: retain only the definitive balanced 12 samples
############################################################

## Uncomment and edit if motif_input contains additional samples.

# final_samples <- c(
#   "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
#   "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
#   "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
#   "SpI_1I3w", "SpI_2I3w", "SpI_4I3w"
# )
#
# motif_input <- motif_input %>%
#   dplyr::filter(Sample %in% final_samples)
#
# sample_metadata <- sample_metadata %>%
#   dplyr::filter(Sample %in% final_samples)

############################################################
## 2. Collapse duplicated CDR3 amino-acid clonotypes
############################################################

motif_input <- motif_input %>%
  dplyr::group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    CDR3.aa
  ) %>%
  dplyr::summarise(
    Clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 3. Prepare full and terminal-trimmed CDR3 sequences
############################################################

## Trimming two residues from each end reduces the influence
## of recurrent V- and J-associated boundary residues.

terminal_trim <- 2

motif_input <- motif_input %>%
  dplyr::mutate(
    CDR3_full = CDR3.aa,
    CDR3_trimmed = dplyr::if_else(
      nchar(CDR3.aa) > 2 * terminal_trim + 2,
      substr(
        CDR3.aa,
        terminal_trim + 1,
        nchar(CDR3.aa) - terminal_trim
      ),
      NA_character_
    )
  )

############################################################
## 4. Function to extract overlapping amino-acid k-mers
############################################################

extract_kmers <- function(sequence, k = 3) {
  
  if (
    is.na(sequence) ||
    sequence == "" ||
    nchar(sequence) < k
  ) {
    return(character(0))
  }
  
  starts <- seq_len(nchar(sequence) - k + 1)
  
  substring(
    sequence,
    first = starts,
    last = starts + k - 1
  )
}

############################################################
## 5. Generate per-sample motif counts
############################################################
make_motif_table <- function(
    input_data,
    sequence_column,
    k = 3,
    weighting = c("clonotype", "abundance")) {
  
  weighting <- match.arg(weighting)
  sequence_column <- rlang::ensym(sequence_column)
  
  input_data %>%
    dplyr::mutate(
      Sequence_for_motifs = as.character(!!sequence_column),
      
      Motif = purrr::map(
        Sequence_for_motifs,
        extract_kmers,
        k = k
      ),
      
      Motif_weight = if (weighting == "clonotype") {
        rep(1, dplyr::n())
      } else {
        as.numeric(Clones)
      }
    ) %>%
    dplyr::filter(!is.na(Sequence_for_motifs)) %>%
    tidyr::unnest_longer(Motif) %>%
    dplyr::filter(
      !is.na(Motif),
      nchar(Motif) == k
    ) %>%
    dplyr::group_by(
      Sample,
      Group,
      Tissue,
      Condition,
      Motif
    ) %>%
    dplyr::summarise(
      Motif_count = sum(Motif_weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(
      Motif_frequency = Motif_count / sum(Motif_count)
    ) %>%
    dplyr::ungroup()
}

###Then rerun all four tables:
motifs_full_clonotype <- make_motif_table(
  motif_input,
  CDR3_full,
  k = 3,
  weighting = "clonotype"
)

motifs_full_abundance <- make_motif_table(
  motif_input,
  CDR3_full,
  k = 3,
  weighting = "abundance"
)

motifs_trimmed_clonotype <- make_motif_table(
  motif_input,
  CDR3_trimmed,
  k = 3,
  weighting = "clonotype"
)

motifs_trimmed_abundance <- make_motif_table(
  motif_input,
  CDR3_trimmed,
  k = 3,
  weighting = "abundance"
)



## B. Motifs are weighted by clone abundance
motifs_full_abundance <- make_motif_table(
  motif_input,
  CDR3_full,
  k = 3,
  weighting = "abundance"
)

## C. Clonotype-weighted after terminal trimming
motifs_trimmed_clonotype <- make_motif_table(
  motif_input,
  CDR3_trimmed,
  k = 3,
  weighting = "clonotype"
)

## D. Abundance-weighted after terminal trimming
motifs_trimmed_abundance <- make_motif_table(
  motif_input,
  CDR3_trimmed,
  k = 3,
  weighting = "abundance"
)

############################################################
## 7. Filter extremely rare motifs
############################################################

## Retain motifs detected in at least 25% of samples.
## This reduces noise and the multiple-testing burden.

filter_motifs <- function(motif_table, minimum_prevalence = 0.25) {
  
  n_samples <- dplyr::n_distinct(motif_table$Sample)
  
  retained <- motif_table %>%
    dplyr::group_by(Motif) %>%
    dplyr::summarise(
      N_samples = dplyr::n_distinct(Sample),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      N_samples >= ceiling(minimum_prevalence * n_samples)
    ) %>%
    dplyr::pull(Motif)
  
  motif_table %>%
    dplyr::filter(Motif %in% retained)
}

motifs_full_clonotype_f <- filter_motifs(
  motifs_full_clonotype
)

motifs_full_abundance_f <- filter_motifs(
  motifs_full_abundance
)

motifs_trimmed_clonotype_f <- filter_motifs(
  motifs_trimmed_clonotype
)

motifs_trimmed_abundance_f <- filter_motifs(
  motifs_trimmed_abundance
)

############################################################
## 8. Convert motif table into sample-by-motif matrix
############################################################

make_motif_matrix <- function(motif_table) {
  
  motif_table %>%
    dplyr::select(
      Sample,
      Motif,
      Motif_frequency
    ) %>%
    tidyr::pivot_wider(
      names_from = Motif,
      values_from = Motif_frequency,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("Sample") %>%
    as.matrix()
}

############################################################
## 9. Global motif analysis function
############################################################

run_global_motif_analysis <- function(
    motif_table,
    metadata,
    analysis_name) {
  
  motif_matrix <- make_motif_matrix(motif_table)
  
  metadata_used <- metadata %>%
    dplyr::filter(Sample %in% rownames(motif_matrix)) %>%
    dplyr::arrange(match(Sample, rownames(motif_matrix)))
  
  stopifnot(
    identical(
      metadata_used$Sample,
      rownames(motif_matrix)
    )
  )
  
  ##########################################################
  ## Hellinger transformation
  ##########################################################
  
  motif_hellinger <- vegan::decostand(
    motif_matrix,
    method = "hellinger"
  )
  
  motif_distance <- vegan::vegdist(
    motif_hellinger,
    method = "euclidean"
  )
  
  ##########################################################
  ## PERMANOVA: tissue ?? condition
  ##########################################################
  
  permanova_interaction <- vegan::adonis2(
    motif_distance ~ Tissue * Condition,
    data = metadata_used,
    permutations = 9999,
    by = "margin"
  )
  
  ##########################################################
  ## Additive model
  ##########################################################
  
  permanova_additive <- vegan::adonis2(
    motif_distance ~ Tissue + Condition,
    data = metadata_used,
    permutations = 9999,
    by = "margin"
  )
  
  ##########################################################
  ## Dispersion analysis
  ##########################################################
  
  group_factor <- interaction(
    metadata_used$Tissue,
    metadata_used$Condition,
    drop = TRUE
  )
  
  dispersion_model <- vegan::betadisper(
    motif_distance,
    group_factor
  )
  
  dispersion_test <- vegan::permutest(
    dispersion_model,
    permutations = 9999
  )
  
  ##########################################################
  ## PCA
  ##########################################################
  
  pca_model <- prcomp(
    motif_hellinger,
    center = TRUE,
    scale. = FALSE
  )
  
  explained_variance <- 100 *
    pca_model$sdev^2 /
    sum(pca_model$sdev^2)
  
  pca_scores <- as.data.frame(pca_model$x) %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::left_join(
      metadata_used,
      by = "Sample"
    )
  
  pca_plot <- ggplot(
    pca_scores,
    aes(
      x = PC1,
      y = PC2,
      colour = Condition,
      shape = Tissue
    )
  ) +
    geom_point(
      size = 4,
      alpha = 0.9
    ) +
    stat_ellipse(
      aes(group = interaction(Tissue, Condition)),
      type = "norm",
      linewidth = 0.6,
      linetype = 2,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c(
        "Control" = "#2C6DB2",
        "Infected" = "#C43C39"
      )
    ) +
    labs(
      title = analysis_name,
      x = paste0(
        "PC1 (",
        round(explained_variance[1], 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(explained_variance[2], 1),
        "%)"
      ),
      colour = "Condition",
      shape = "Tissue"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 13
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      )
    )
  
  ##########################################################
  ## Return results
  ##########################################################
  
  list(
    matrix = motif_matrix,
    distance = motif_distance,
    metadata = metadata_used,
    permanova_interaction = permanova_interaction,
    permanova_additive = permanova_additive,
    dispersion_test = dispersion_test,
    pca_scores = pca_scores,
    pca_plot = pca_plot
  )
}

############################################################
## 10. Run global analyses
############################################################

global_full_clonotype <- run_global_motif_analysis(
  motifs_full_clonotype_f,
  sample_metadata,
  "Full CDR3: clonotype-weighted motifs"
)

global_full_abundance <- run_global_motif_analysis(
  motifs_full_abundance_f,
  sample_metadata,
  "Full CDR3: abundance-weighted motifs"
)

global_trimmed_clonotype <- run_global_motif_analysis(
  motifs_trimmed_clonotype_f,
  sample_metadata,
  "Trimmed CDR3: clonotype-weighted motifs"
)

global_trimmed_abundance <- run_global_motif_analysis(
  motifs_trimmed_abundance_f,
  sample_metadata,
  "Trimmed CDR3: abundance-weighted motifs"
)

############################################################
## 11. Inspect PERMANOVA and dispersion results
############################################################

global_full_clonotype$permanova_interaction
global_full_clonotype$permanova_additive
global_full_clonotype$dispersion_test

global_full_abundance$permanova_interaction
global_full_abundance$permanova_additive
global_full_abundance$dispersion_test

global_trimmed_clonotype$permanova_interaction
global_trimmed_clonotype$permanova_additive
global_trimmed_clonotype$dispersion_test

global_trimmed_abundance$permanova_interaction
global_trimmed_abundance$permanova_additive
global_trimmed_abundance$dispersion_test

############################################################
## 12. Display PCA figures
############################################################

global_full_clonotype$pca_plot
global_full_abundance$pca_plot
global_trimmed_clonotype$pca_plot
global_trimmed_abundance$pca_plot

############################################################
## 13. Calculate pairwise cosine similarity
############################################################

cosine_similarity <- function(motif_matrix) {
  
  numerator <- motif_matrix %*% t(motif_matrix)
  
  magnitude <- sqrt(
    rowSums(motif_matrix^2)
  )
  
  denominator <- outer(magnitude, magnitude)
  
  numerator / denominator
}

cosine_matrix <- cosine_similarity(
  global_trimmed_clonotype$matrix
)

cosine_long <- as.data.frame(as.table(cosine_matrix)) %>%
  dplyr::rename(
    Sample_1 = Var1,
    Sample_2 = Var2,
    Cosine_similarity = Freq
  ) %>%
  dplyr::filter(
    as.character(Sample_1) <
      as.character(Sample_2)
  ) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(
        Sample,
        Tissue_1 = Tissue,
        Condition_1 = Condition
      ),
    by = c("Sample_1" = "Sample")
  ) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(
        Sample,
        Tissue_2 = Tissue,
        Condition_2 = Condition
      ),
    by = c("Sample_2" = "Sample")
  ) %>%
  dplyr::mutate(
    Comparison = dplyr::case_when(
      Tissue_1 == Tissue_2 &
        Condition_1 == Condition_2 ~ "Within group",
      
      Tissue_1 == Tissue_2 &
        Condition_1 != Condition_2 ~ "Between conditions",
      
      Tissue_1 != Tissue_2 &
        Condition_1 == Condition_2 ~ "Between tissues",
      
      TRUE ~ "Different tissue and condition"
    )
  )

cosine_summary <- cosine_long %>%
  dplyr::group_by(Comparison) %>%
  dplyr::summarise(
    N_pairs = dplyr::n(),
    Median_similarity = median(Cosine_similarity),
    IQR_similarity = IQR(Cosine_similarity),
    Minimum = min(Cosine_similarity),
    Maximum = max(Cosine_similarity),
    .groups = "drop"
  )

print(cosine_summary)

############################################################
## 14. Plot global similarity
############################################################

cosine_plot <- ggplot(
  cosine_long,
  aes(
    x = Comparison,
    y = Cosine_similarity,
    fill = Comparison
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.65
  ) +
  scale_fill_manual(
    values = c(
      "Within group" = "#8064A2",
      "Between conditions" = "#E45756",
      "Between tissues" = "#4C78A8",
      "Different tissue and condition" = "grey65"
    )
  ) +
  labs(
    title = "Similarity of CDR3 amino-acid motif profiles",
    x = NULL,
    y = "Cosine similarity"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    )
  )

cosine_plot

############################################################
## 15. Differential representation of individual motifs
############################################################

## The analysis is performed separately within each tissue.
## With n = 3 per condition, results should be interpreted
## as exploratory and reported with FDR correction.

test_individual_motifs <- function(
    motif_table,
    minimum_samples = 3) {
  
  all_samples <- motif_table %>%
    dplyr::distinct(
      Sample,
      Tissue,
      Condition
    )
  
  motif_complete <- motif_table %>%
    dplyr::select(
      Sample,
      Tissue,
      Condition,
      Motif,
      Motif_frequency
    ) %>%
    tidyr::complete(
      Sample = all_samples$Sample,
      Motif,
      fill = list(Motif_frequency = 0)
    ) %>%
    dplyr::select(
      Sample,
      Motif,
      Motif_frequency
    ) %>%
    dplyr::left_join(
      all_samples,
      by = "Sample"
    )
  
  motif_complete %>%
    dplyr::group_by(Tissue, Motif) %>%
    dplyr::filter(
      sum(Motif_frequency > 0) >= minimum_samples,
      dplyr::n_distinct(Condition) == 2
    ) %>%
    dplyr::summarise(
      Control_median = median(
        Motif_frequency[Condition == "Control"],
        na.rm = TRUE
      ),
      
      Infected_median = median(
        Motif_frequency[Condition == "Infected"],
        na.rm = TRUE
      ),
      
      Log2_fold_change = log2(
        (Infected_median + 1e-7) /
          (Control_median + 1e-7)
      ),
      
      P_value = tryCatch(
        wilcox.test(
          Motif_frequency[Condition == "Infected"],
          Motif_frequency[Condition == "Control"],
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      ),
      
      .groups = "drop"
    ) %>%
    dplyr::group_by(Tissue) %>%
    dplyr::mutate(
      Adjusted_P = p.adjust(
        P_value,
        method = "BH"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Adjusted_P, P_value)
}

differential_full_clonotype <- test_individual_motifs(
  motifs_full_clonotype_f
)

differential_full_abundance <- test_individual_motifs(
  motifs_full_abundance_f
)

differential_trimmed_clonotype <- test_individual_motifs(
  motifs_trimmed_clonotype_f
)

differential_trimmed_abundance <- test_individual_motifs(
  motifs_trimmed_abundance_f
)

## Inspect the strongest motifs
differential_trimmed_clonotype %>%
  dplyr::filter(P_value < 0.05) %>%
  print(n = 50)

differential_trimmed_abundance %>%
  dplyr::filter(P_value < 0.05) %>%
  print(n = 50)

## Motifs significant after FDR correction
differential_trimmed_clonotype %>%
  dplyr::filter(Adjusted_P < 0.05)

differential_trimmed_abundance %>%
  dplyr::filter(Adjusted_P < 0.05)

############################################################
## 16. Compare clonotype- and abundance-weighted results
############################################################

motif_comparison <- differential_trimmed_clonotype %>%
  dplyr::select(
    Tissue,
    Motif,
    Clonotype_log2FC = Log2_fold_change,
    Clonotype_P = P_value,
    Clonotype_FDR = Adjusted_P
  ) %>%
  dplyr::full_join(
    differential_trimmed_abundance %>%
      dplyr::select(
        Tissue,
        Motif,
        Abundance_log2FC = Log2_fold_change,
        Abundance_P = P_value,
        Abundance_FDR = Adjusted_P
      ),
    by = c("Tissue", "Motif")
  ) %>%
  dplyr::mutate(
    Same_direction = sign(Clonotype_log2FC) ==
      sign(Abundance_log2FC),
    
    Interpretation = dplyr::case_when(
      Clonotype_FDR < 0.05 &
        Abundance_FDR < 0.05 &
        Same_direction ~
        "Consistent compositional and expansion effect",
      
      Clonotype_FDR < 0.05 &
        (is.na(Abundance_FDR) | Abundance_FDR >= 0.05) ~
        "Primarily repertoire-composition effect",
      
      Abundance_FDR < 0.05 &
        (is.na(Clonotype_FDR) | Clonotype_FDR >= 0.05) ~
        "Primarily clonal-expansion effect",
      
      TRUE ~ "No FDR-significant effect"
    )
  ) %>%
  dplyr::arrange(
    pmin(
      Clonotype_FDR,
      Abundance_FDR,
      na.rm = TRUE
    )
  )

print(motif_comparison, n = 50)

############################################################
## 17. Plot the strongest altered motifs
############################################################

top_motifs <- differential_trimmed_abundance %>%
  dplyr::filter(!is.na(P_value)) %>%
  dplyr::group_by(Tissue) %>%
  dplyr::slice_min(
    order_by = P_value,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::pull(Motif) %>%
  unique()

top_motif_data <- motifs_trimmed_abundance_f %>%
  dplyr::filter(Motif %in% top_motifs) %>%
  tidyr::complete(
    nesting(
      Sample,
      Group,
      Tissue,
      Condition
    ),
    Motif = top_motifs,
    fill = list(
      Motif_count = 0,
      Motif_frequency = 0
    )
  )

top_motif_plot <- ggplot(
  top_motif_data,
  aes(
    x = Condition,
    y = Motif_frequency,
    colour = Condition
  )
) +
  geom_boxplot(
    aes(fill = Condition),
    width = 0.60,
    alpha = 0.25,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    size = 1.7
  ) +
  facet_grid(
    Tissue ~ Motif,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "Control" = "#2C6DB2",
      "Infected" = "#C43C39"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#2C6DB2",
      "Infected" = "#C43C39"
    )
  ) +
  labs(
    title = "CDR3 motifs showing the strongest infection-associated changes",
    x = NULL,
    y = "Normalized motif frequency"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

top_motif_plot

####End of the above motiff analysis
################################################################


############################################################
## 11. Convergent CDR3 similarity clusters
## Warning: can be slow if many clones; uses top abundant clones
############################################################

top_clones_for_network <- all_clones %>%
  group_by(CDR3.aa) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_clones)) %>%
  slice_head(n = 500)

seqs <- top_clones_for_network$CDR3.aa

dist_mat <- stringdistmatrix(seqs, seqs, method = "lv")

edge_df <- which(dist_mat <= 2 & dist_mat > 0, arr.ind = TRUE) %>%
  as.data.frame() %>%
  filter(row < col) %>%
  mutate(
    from = seqs[row],
    to = seqs[col],
    distance = dist_mat[cbind(row, col)]
  ) %>%
  select(from, to, distance)

cdr3_graph <- graph_from_data_frame(edge_df, directed = FALSE)

plot(
  cdr3_graph,
  vertex.size = 4,
  vertex.label = NA,
  main = "CDR3 amino-acid similarity Network"
)

clusters <- components(cdr3_graph)

cluster_table <- data.frame(
  CDR3.aa = names(clusters$membership),
  Cluster = clusters$membership
) %>%
  left_join(top_clones_for_network, by = "CDR3.aa")

#View(cluster_table)

############################################################
## 12. Approximate lineage analysis
## Lineage definition: same V family + same J family + same CDR3 length
############################################################
all_clones_full <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa,
        CDR3.nt,
        Clones,
        Proportion,
        V.name,
        D.name,
        J.name
      ) %>%
      mutate(Sample = s)
  })
) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    
    V_gene_no_allele = sub("\\*.*", "", V.name),
    J_gene_no_allele = sub("\\*.*", "", J.name),
    
    V_family = sub("-.*", "", V_gene_no_allele),
    J_family = sub("-.*", "", J_gene_no_allele),
    
    CDR3_length = nchar(CDR3.aa)
  )

##
lineage_df <- all_clones_full %>%
  mutate(
    LineageID = paste(
      V_family,
      J_family,
      CDR3_length,
      sep = "_"
    )
  )

lineage_summary <- lineage_df %>%
  group_by(Group, LineageID) %>%
  summarise(
    n_CDR3 = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    mean_CDR3_length = mean(CDR3_length, na.rm = TRUE),
    .groups = "drop"
  )

top_lineages <- lineage_summary %>%
  group_by(Group) %>%
  slice_max(total_clones, n = 20, with_ties = FALSE) %>%
  ungroup()

ggplot(top_lineages,
       aes(reorder(LineageID, total_clones), total_clones, fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Top 20 V-J-CDR3-Length Configration Classes",
    x = "V-J family-CDR3-length classes",
    y = "Total clone count"
  )

############################################################
## 13. Exploratory SHM proxy
## True SHM requires germline reference alignment.
## Here we use sequence diversity within approximate lineages.
############################################################

lineage_diversity_proxy <- lineage_df %>%
  group_by(Group, LineageID) %>%
  summarise(
    n_nt_variants = n_distinct(CDR3.nt),
    n_aa_variants = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(lineage_diversity_proxy,
       aes(x = Group, y = n_nt_variants, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_bw(base_size = 14) +
  labs(
    title = "Exploratory Lineage Sequence Diversity",
    subtitle = "Proxy only; true SHM requires germline alignment",
    x = "Group",
    y = "Number of nucleotide variants per lineage"
  )

###stat
library(ggpubr)

ggplot(
  lineage_diversity_proxy,
  aes(Group, n_nt_variants, fill = Group)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  
  stat_compare_means(
    method = "kruskal.test",
    label.y = max(lineage_diversity_proxy$n_nt_variants)*1.1
  ) +
  
  stat_compare_means(
    comparisons = list(
      c("HKC","HKI"),
      c("SpC","SpI")
    ),
    method = "wilcox.test"
  ) +
  
  theme_bw(base_size = 14)+
  labs(
    title = "Exploratory Lineage Sequence Diversity",
    subtitle = "Proxy only; true SHM requires germline alignment",
    x = "Group",
    y = "Number of nucleotide variants per lineage"
  )




############################################################
## 14. Machine-learning style feature matrix
############################################################

v_features <- all_clones_full %>%
  group_by(Sample, V_family) %>%
  summarise(value = sum(Proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = V_family,
    values_from = value,
    values_fill = 0,
    names_prefix = "Vfam_"
  )

j_features <- all_clones_full %>%
  group_by(Sample, J_family) %>%
  summarise(value = sum(Proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = J_family,
    values_from = value,
    values_fill = 0,
    names_prefix = "Jfam_"
  )

length_features <- all_clones_full %>%
  group_by(Sample) %>%
  summarise(
    mean_CDR3_length = mean(CDR3_length, na.rm = TRUE),
    median_CDR3_length = median(CDR3_length, na.rm = TRUE),
    n_unique_CDR3 = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  )

public_features <- public_by_group %>%
  group_by(Group) %>%
  summarise(n_public_group = n_distinct(CDR3.aa), .groups = "drop")

feature_matrix <- meta %>%
  select(Sample, Group, Tissue, Condition) %>%
  left_join(v_features, by = "Sample") %>%
  left_join(j_features, by = "Sample") %>%
  left_join(length_features, by = "Sample")

#View(feature_matrix)

#########fig
############################################################
## Figure from machine-learning style feature matrix
## PCA + feature heatmap
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)

############################################################
## 1. Replace NA values
############################################################

feature_matrix_clean <- feature_matrix %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

############################################################
## 2. PCA of repertoire features
############################################################

feature_numeric <- feature_matrix_clean %>%
  select(where(is.numeric))

## Remove zero-variance columns
feature_numeric <- feature_numeric[
  ,
  apply(feature_numeric, 2, var, na.rm = TRUE) > 0,
  drop = FALSE
]

pca <- prcomp(
  feature_numeric,
  center = TRUE,
  scale. = TRUE
)

pca_df <- data.frame(
  Sample = feature_matrix_clean$Sample,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(
    feature_matrix_clean %>% select(Sample, Group, Tissue, Condition),
    by = "Sample"
  )

var_exp <- summary(pca)$importance[2, ] * 100

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Group, shape = Tissue)
) +
  geom_point(size = 4) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of Repertoire Feature Matrix",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  )

p_pca

############################################################
## 3. Heatmap of scaled repertoire features
############################################################

feature_heatmap_mat <- feature_numeric %>%
  as.data.frame()

rownames(feature_heatmap_mat) <- feature_matrix_clean$Sample

## Z-score features
feature_heatmap_z <- scale(feature_heatmap_mat)

##
annotation_col <- feature_matrix_clean %>%
  select(Sample, Group, Tissue, Condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample

annotation_col$Sample <- NULL


pheatmap(
  feature_heatmap_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  main = "Repertoire Feature Matrix"
)

############################################################
## 4. Optional: top PCA loadings
############################################################

loadings_df <- data.frame(
  Feature = rownames(pca$rotation),
  PC1 = pca$rotation[, 1],
  PC2 = pca$rotation[, 2]
)

top_loadings <- loadings_df %>%
  mutate(abs_PC1 = abs(PC1)) %>%
  arrange(desc(abs_PC1)) %>%
  slice_head(n = 20)

p_loadings <- ggplot(
  top_loadings,
  aes(x = reorder(Feature, abs_PC1), y = abs_PC1)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top Features Driving PC1",
    x = "Feature",
    y = "Absolute loading"
  )

p_loadings


############################################################
## END
############################################################

############################################################
## TIER 8: Tissue-specific infection response
## Model: Feature ~ Tissue * Condition
## Groups:
## HKC = HK Control
## HKI = HK Infected
## SpC = Spleen Control
## SpI = Spleen Infected
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(ggpubr)

############################################################
## 1. Metadata
############################################################

meta <- data.frame(Sample = names(imm_ds)) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

table(meta$Group)
table(meta$Tissue, meta$Condition)

############################################################
## 2. Build full clone table
############################################################

all_clones_full <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa,
        CDR3.nt,
        Clones,
        Proportion,
        V.name,
        J.name
      ) %>%
      mutate(Sample = s)
  })
) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    
    V_gene_no_allele = sub("\\*.*", "", V.name),
    J_gene_no_allele = sub("\\*.*", "", J.name),
    
    #V_family = sub("-.*", "", V_gene_no_allele),
    #J_family = sub("-.*", "", J_gene_no_allele),
    V_family = sub("-.*$", "", sub("_?S[0-9]+$", "", V_gene_no_allele)),
    J_family = sub("-.*$", "", sub("_?T[0-9]+(?=D?$)", "", J_gene_no_allele, perl = TRUE)),
    
    CDR3_length = nchar(CDR3.aa),
    
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

############################################################
## 3. Build sample-level feature table
############################################################

sample_features_basic <- all_clones_full %>%
  group_by(Sample, Group, Tissue, Condition) %>%
  summarise(
    total_reads = sum(Clones, na.rm = TRUE),
    n_unique_CDR3 = n_distinct(CDR3.aa),
    mean_CDR3_length = mean(CDR3_length, na.rm = TRUE),
    median_CDR3_length = median(CDR3_length, na.rm = TRUE),
    mean_clone_size = mean(Clones, na.rm = TRUE),
    max_clone_size = max(Clones, na.rm = TRUE),
    top_clone_proportion = max(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 4. Add public/private clone features
############################################################

clone_occurrence <- all_clones_full %>%
  group_by(CDR3.aa) %>%
  summarise(
    n_samples_total = n_distinct(Sample),
    .groups = "drop"
  )

private_clones <- clone_occurrence %>%
  filter(n_samples_total == 1) %>%
  pull(CDR3.aa)

shared_clones <- clone_occurrence %>%
  filter(n_samples_total >= 2) %>%
  pull(CDR3.aa)

sample_public_private <- all_clones_full %>%
  mutate(
    Clone_type = ifelse(CDR3.aa %in% private_clones,
                        "Private", "Shared")
  ) %>%
  group_by(Sample, Group, Tissue, Condition, Clone_type) %>%
  summarise(
    n_clones = n_distinct(CDR3.aa),
    abundance = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Clone_type,
    values_from = c(n_clones, abundance),
    values_fill = 0
  )

sample_features <- sample_features_basic %>%
  left_join(
    sample_public_private,
    by = c("Sample", "Group", "Tissue", "Condition")
  ) %>%
  mutate(
    private_fraction =
      n_clones_Private / n_unique_CDR3,
    shared_fraction =
      n_clones_Shared / n_unique_CDR3,
    private_abundance_fraction =
      abundance_Private / total_reads,
    shared_abundance_fraction =
      abundance_Shared / total_reads
  )

############################################################
## 5. Add V-family features
############################################################

v_family_features <- all_clones_full %>%
  group_by(Sample, V_family) %>%
  summarise(
    value = sum(Proportion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = V_family,
    values_from = value,
    values_fill = 0,
    names_prefix = "Vfam_"
  )

############################################################
## 6. Add J-family features
############################################################

j_family_features <- all_clones_full %>%
  group_by(Sample, J_family) %>%
  summarise(
    value = sum(Proportion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = J_family,
    values_from = value,
    values_fill = 0,
    names_prefix = "Jfam_"
  )

############################################################
## 7. Final feature matrix
############################################################

feature_matrix <- sample_features %>%
  left_join(v_family_features, by = "Sample") %>%
  left_join(j_family_features, by = "Sample") %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

#View(feature_matrix)

############################################################
## 8. Function for Tissue * Condition model
############################################################

run_interaction_model <- function(data, feature_name) {
  
  formula_use <- as.formula(
    paste0("`", feature_name, "` ~ Tissue * Condition")
  )
  
  model <- lm(formula_use, data = data)
  
  tidy(model) %>%
    mutate(
      Feature = feature_name
    )
}

numeric_features <- feature_matrix %>%
  select(where(is.numeric)) %>%
  colnames()

interaction_results <- bind_rows(
  lapply(numeric_features, function(f) {
    tryCatch(
      run_interaction_model(feature_matrix, f),
      error = function(e) NULL
    )
  })
)

interaction_results <- interaction_results %>%
  group_by(term) %>%
  mutate(
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()

################################################8888
library(dplyr)
library(broom)

## 1. Keep only numeric features with variation
numeric_features <- feature_matrix %>%
  select(where(is.numeric)) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0)) %>%
  colnames()

numeric_features

## 2. Safer model function
run_interaction_model <- function(data, feature_name) {
  
  model_df <- data %>%
    select(Tissue, Condition, all_of(feature_name)) %>%
    rename(Value = all_of(feature_name)) %>%
    filter(!is.na(Value))
  
  if (nrow(model_df) < 6) return(NULL)
  if (var(model_df$Value, na.rm = TRUE) == 0) return(NULL)
  
  model <- lm(Value ~ Tissue * Condition, data = model_df)
  
  broom::tidy(model) %>%
    mutate(Feature = feature_name)
}

## 3. Run models
interaction_results <- bind_rows(
  lapply(numeric_features, function(f) {
    tryCatch(
      run_interaction_model(feature_matrix, f),
      error = function(e) {
        message("Failed: ", f, " | ", e$message)
        NULL
      }
    )
  })
)

## 4. Check
dim(interaction_results)
head(interaction_results)
colnames(interaction_results)

## 5. Adjust p-values only if results exist
if (nrow(interaction_results) > 0) {
  
  interaction_results <- interaction_results %>%
    group_by(term) %>%
    mutate(
      padj = p.adjust(p.value, method = "BH")
    ) %>%
    ungroup()
  
  interaction_effects <- interaction_results %>%
    filter(term == "TissueSpleen:ConditionInfected") %>%
    arrange(padj)
  
  View(interaction_effects)
}
################################################888





#View(interaction_results)

############################################################
## 9. Extract important effects
############################################################

condition_effects <- interaction_results %>%
  filter(term == "ConditionInfected") %>%
  arrange(padj)

tissue_effects <- interaction_results %>%
  filter(term == "TissueSpleen") %>%
  arrange(padj)

interaction_effects <- interaction_results %>%
  filter(term == "TissueSpleen:ConditionInfected") %>%
  arrange(padj)

#View(condition_effects)
#View(tissue_effects)
#View(interaction_effects)

############################################################
## 10. Plot top interaction effects
############################################################

top_interaction_features <- interaction_effects %>%
  arrange(padj) %>%
  slice_head(n = 12) %>%
  pull(Feature)


length(top_interaction_features)

top_interaction_features
dev.list()
graphics.off()
dev.list()

##
interaction_plot_df <- feature_matrix %>%
  select(Sample, Group, Tissue, Condition, all_of(top_interaction_features)) %>%
  pivot_longer(
    cols = all_of(top_interaction_features),
    names_to = "Feature",
    values_to = "Value"
  )
##

p <- ggplot(
  interaction_plot_df,
  aes(x = Condition, y = Value, fill = Condition)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_grid(Feature ~ Tissue, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Top Tissue-Specific Infection Response Features",
    x = "",
    y = "Feature value"
  )

print(p)

############################################################
## 11. Heatmap of model effects
############################################################

effect_heatmap <- interaction_results %>%
  filter(
    term %in% c(
      "TissueSpleen",
      "ConditionInfected",
      "TissueSpleen:ConditionInfected"
    )
  ) %>%
  select(Feature, term, estimate) %>%
  pivot_wider(
    names_from = term,
    values_from = estimate,
    values_fill = 0
  )

effect_mat <- effect_heatmap %>%
  tibble::column_to_rownames("Feature") %>%
  as.matrix()

effect_mat <- effect_mat[
  apply(effect_mat, 1, function(x) any(abs(x) > 0)),
  ,
  drop = FALSE
]

pheatmap::pheatmap(
  effect_mat,
  scale = "row",
  main = "Tissue, Infection, and Interaction Effects",
  fontsize_row = 6
)

############################################################
## 12. Explicit HKC vs HKI and SpC vs SpI tests
############################################################

pairwise_feature_tests <- bind_rows(
  lapply(numeric_features, function(f) {
    
    hk_df <- feature_matrix %>%
      filter(Group %in% c("HKC", "HKI"))
    
    sp_df <- feature_matrix %>%
      filter(Group %in% c("SpC", "SpI"))
    
    hk_p <- tryCatch(
      wilcox.test(hk_df[[f]] ~ hk_df$Group, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    
    sp_p <- tryCatch(
      wilcox.test(sp_df[[f]] ~ sp_df$Group, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    
    data.frame(
      Feature = f,
      HKC_vs_HKI_p = hk_p,
      SpC_vs_SpI_p = sp_p
    )
  })
) %>%
  mutate(
    HKC_vs_HKI_padj = p.adjust(HKC_vs_HKI_p, method = "BH"),
    SpC_vs_SpI_padj = p.adjust(SpC_vs_SpI_p, method = "BH")
  )

#View(pairwise_feature_tests)

############################################################
## 13. Plot selected biologically important features
############################################################

selected_features <- c(
  "n_unique_CDR3",
  "mean_CDR3_length",
  "top_clone_proportion",
  "private_fraction",
  "shared_fraction",
  "private_abundance_fraction",
  "shared_abundance_fraction"
)

selected_features <- selected_features[
  selected_features %in% colnames(feature_matrix)
]

selected_plot_df <- feature_matrix %>%
  select(Sample, Group, Tissue, Condition, all_of(selected_features)) %>%
  pivot_longer(
    cols = all_of(selected_features),
    names_to = "Feature",
    values_to = "Value"
  )

##########fig
p <- ggplot(
  selected_plot_df,
  aes(x = Condition, y = Value, fill = Condition)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2) +
  facet_grid(
    Feature ~ Tissue,
    scales = "free_y"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text.y = element_text(size = 9),
    strip.text.x = element_text(size = 10, face = "bold"),
    
    ## increase margins around figure
    plot.margin = margin(
      t = 10,
      r = 40,
      b = 10,
      l = 10
    ),
    
    ## allow labels to extend beyond panel
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Tissue-Specific Infection Response",
    x = "",
    y = "Feature value"
  )

print(p)

############################################################
############################################################
## REAL SHM ANALYSIS FOR IGH DATA
## Requires: germline IGHV FASTA file
############################################################
############################################################
## REAL SHM ANALYSIS USING LOCAL IMGT GERMLINE FILES
## Folder:
## C:/Users/sje054/OneDrive - UiT Office 365/AquaTrep/IGH_genes
##
## Uses:
## - IGHV file for real V-region SHM frequency
## - IGHD/IGHJ files are read and checked, but SHM is mainly V-region based
############################################################

library(stringr)
library(Biostrings)
library(ggpubr)

############################################################
## 1. Germline folder
############################################################

GERMLINE_DIR <- "C:/Users/sje054/OneDrive - UiT Office 365/AquaTrep/IGH_genes"

list.files(GERMLINE_DIR, full.names = TRUE)

############################################################
## 2. Find V, D, and J germline files automatically
############################################################

all_gene_files <- list.files(
  GERMLINE_DIR,
  full.names = TRUE
)

IGHV_file <- all_gene_files[grepl("IGHV|V", basename(all_gene_files), ignore.case = TRUE)][1]
IGHD_file <- all_gene_files[grepl("IGHD|D", basename(all_gene_files), ignore.case = TRUE)][1]
IGHJ_file <- all_gene_files[grepl("IGHJ|J", basename(all_gene_files), ignore.case = TRUE)][1]

IGHV_file
IGHD_file
IGHJ_file

############################################################
## 3. Function to read IMGT-style germline files
## Handles dots in aligned sequences
############################################################

read_imgt_germline <- function(file_path, gene_type = c("V", "D", "J")) {
  
  gene_type <- match.arg(gene_type)
  
  raw_lines <- readLines(file_path)
  raw_lines <- raw_lines[nchar(raw_lines) > 0]
  
  header_idx <- grep("^>", raw_lines)
  
  records <- lapply(seq_along(header_idx), function(i) {
    
    start <- header_idx[i]
    end <- ifelse(
      i < length(header_idx),
      header_idx[i + 1] - 1,
      length(raw_lines)
    )
    
    header <- raw_lines[start]
    seq_lines <- raw_lines[(start + 1):end]
    
    seq <- paste(seq_lines, collapse = "")
    
    ## Remove IMGT alignment dots and non-DNA characters
    seq <- gsub("\\.", "", seq)
    seq <- gsub("[^ACGTacgt]", "", seq)
    seq <- toupper(seq)
    
    fields <- str_split_fixed(sub("^>", "", header), "\\|", 8)
    
    data.frame(
      Header = sub("^>", "", header),
      Gene_name = fields[, 2],
      Species = fields[, 3],
      Functionality = fields[, 4],
      Region = fields[, 5],
      Germline_seq = seq,
      Gene_type = gene_type,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(records) %>%
    mutate(
      Gene_no_allele = sub("\\*.*", "", Gene_name),
      Gene_family = sub("-.*", "", Gene_no_allele),
      Germline_length = nchar(Germline_seq)
    )
}

############################################################
## 4. Read germline V, D, J files
############################################################

germline_V <- read_imgt_germline(IGHV_file, gene_type = "V")
germline_D <- read_imgt_germline(IGHD_file, gene_type = "D")
germline_J <- read_imgt_germline(IGHJ_file, gene_type = "J")

head(germline_V)
head(germline_D)
head(germline_J)

table(germline_V$Functionality)
table(germline_D$Functionality)
table(germline_J$Functionality)

############################################################
## 5. Build clone table from imm_ds
############################################################

all_clones_full <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        Clones,
        Proportion,
        CDR3.nt,
        CDR3.aa,
        V.name,
        D.name,
        J.name,
        Sequence,
        V.end,
        CDR1.nt,
        CDR2.nt,
        FR1.nt,
        FR2.nt,
        FR3.nt
      ) %>%
      mutate(Sample = s)
  })
) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected")),
    
    #V_gene_no_allele = sub("\\*.*", "", V.name),
    V_gene_no_allele = sub("S[0-9]+$", "", sub("\\*.*$", "", V.name)),
    D_gene_no_allele = sub("\\*.*", "", D.name),
    J_gene_no_allele = sub("\\*.*", "", J.name),
    
    V_family = sub("-.*", "", V_gene_no_allele),
    D_family = sub("-.*", "", D_gene_no_allele),
    J_family = sub("-.*", "", J_gene_no_allele),
    
    V_region_nt = substr(Sequence, 1, V.end)
  )

############################################################
## 6. Match observed V genes to germline V sequences
############################################################

shm_df <- all_clones_full %>%
  left_join(
    germline_V %>%
      select(
        V_gene_no_allele = Gene_no_allele,
        Germline_V_name = Gene_name,
        Germline_seq
      ),
    by = "V_gene_no_allele"
  ) %>%
  filter(
    !is.na(Germline_seq),
    !is.na(V_region_nt),
    nchar(V_region_nt) > 0
  )

## Check matching success
match_summary <- all_clones_full %>%
  distinct(V.name, V_gene_no_allele) %>%
  left_join(
    germline_V %>%
      distinct(V_gene_no_allele = Gene_no_allele, Germline_seq),
    by = "V_gene_no_allele"
  ) %>%
  mutate(Matched = !is.na(Germline_seq))

table(match_summary$Matched)

match_summary %>%
  filter(!Matched) %>%
  print(n = 50)

############################################################
## 7. Function: mutation frequency by V-region alignment
############################################################
calc_mutation_stats <- function(query, germline) {
  
  if (is.na(query) || is.na(germline)) {
    return(data.frame(
      compared_nt = NA_integer_,
      mutations = NA_integer_,
      mutation_freq = NA_real_
    ))
  }
  
  query <- toupper(query)
  germline <- toupper(germline)
  
  query <- gsub("[^ACGT]", "", query)
  germline <- gsub("[^ACGT]", "", germline)
  
  min_len <- min(nchar(query), nchar(germline))
  
  if (min_len == 0) {
    return(data.frame(
      compared_nt = NA_integer_,
      mutations = NA_integer_,
      mutation_freq = NA_real_
    ))
  }
  
  query <- substr(query, 1, min_len)
  germline <- substr(germline, 1, min_len)
  
  q <- strsplit(query, "")[[1]]
  g <- strsplit(germline, "")[[1]]
  
  mutations <- sum(q != g)
  
  data.frame(
    compared_nt = min_len,
    mutations = mutations,
    mutation_freq = mutations / min_len
  )
}

############################################################
## 8. Run SHM per sequence
############################################################
shm_results <- bind_rows(
  lapply(seq_len(nrow(shm_df)), function(i) {
    
    stats <- calc_mutation_stats(
      query = shm_df$V_region_nt[i],
      germline = shm_df$Germline_seq[i]
    )
    
    bind_cols(
      shm_df[i, ],
      stats
    )
  })
)

#View(shm_results)

############################################################
## 9. SHM per sample
############################################################

sample_shm <- shm_results %>%
  group_by(Sample, Group, Tissue, Condition) %>%
  summarise(
    mean_mutation_freq = weighted.mean(
      mutation_freq,
      w = Clones,
      na.rm = TRUE
    ),
    median_mutation_freq = median(mutation_freq, na.rm = TRUE),
    total_mutations = sum(mutations * Clones, na.rm = TRUE),
    total_compared_nt = sum(compared_nt * Clones, na.rm = TRUE),
    overall_mutation_freq = total_mutations / total_compared_nt,
    .groups = "drop"
  )

#View(sample_shm)

saveRDS(shm_df, "shm_df.rds")
############################################################
## 10. SHM distribution by group
############################################################

p_shm_group <- ggplot(
  sample_shm,
  aes(x = Group, y = overall_mutation_freq, fill = Group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 3) +
  stat_compare_means(method = "kruskal.test") +
  theme_bw(base_size = 14) +
  labs(
    title = "IGH V-region SHM Frequency",
    x = "Group",
    y = "Mutation frequency"
  )

p_shm_group

############################################################
## 11. HKC vs HKI and SpC vs SpI tests
############################################################

HK_shm_test <- wilcox.test(
  overall_mutation_freq ~ Group,
  data = sample_shm %>% filter(Group %in% c("HKC", "HKI")),
  exact = FALSE
)

Sp_shm_test <- wilcox.test(
  overall_mutation_freq ~ Group,
  data = sample_shm %>% filter(Group %in% c("SpC", "SpI")),
  exact = FALSE
)

HK_shm_test
Sp_shm_test

############################################################
## 12. Tissue x infection model
############################################################

shm_lm <- lm(
  overall_mutation_freq ~ Tissue * Condition,
  data = sample_shm
)

summary(shm_lm)

############################################################
## 13. SHM by V family
############################################################

v_family_shm <- shm_results %>%
  group_by(Sample, Group, Tissue, Condition, V_family) %>%
  summarise(
    total_mutations = sum(mutations * Clones, na.rm = TRUE),
    total_compared_nt = sum(compared_nt * Clones, na.rm = TRUE),
    mutation_freq = total_mutations / total_compared_nt,
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_compared_nt > 0)

top_v_families <- v_family_shm %>%
  group_by(V_family) %>%
  summarise(
    total_clones = sum(total_clones),
    .groups = "drop"
  ) %>%
  slice_max(total_clones, n = 15)

v_family_shm_top <- v_family_shm %>%
  filter(V_family %in% top_v_families$V_family)

p_v_family_shm <- ggplot(
  v_family_shm_top,
  aes(x = Group, y = mutation_freq, fill = Group)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~ V_family, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "SHM Frequency by V Family",
    x = "",
    y = "Mutation frequency"
  )

p_v_family_shm

############################################################
## 14. Mutation spectrum
############################################################
get_mutation_spectrum <- function(query, germline) {
  
  if (is.na(query) || is.na(germline)) {
    return(data.frame())
  }
  
  query <- toupper(query)
  germline <- toupper(germline)
  
  query <- gsub("[^ACGT]", "", query)
  germline <- gsub("[^ACGT]", "", germline)
  
  min_len <- min(
    nchar(query),
    nchar(germline)
  )
  
  if (min_len == 0) {
    return(data.frame())
  }
  
  query <- substr(query, 1, min_len)
  germline <- substr(germline, 1, min_len)
  
  q <- strsplit(query, "")[[1]]
  g <- strsplit(germline, "")[[1]]
  
  mut_idx <- which(q != g)
  
  if (length(mut_idx) == 0) {
    return(data.frame())
  }
  
  data.frame(
    Position = mut_idx,
    From = g[mut_idx],
    To = q[mut_idx],
    Mutation = paste0(g[mut_idx], ">", q[mut_idx])
  )
}



#############################################################
#New analysis using the files imported from Saga
##SHM spectrum complete script
############################################################
## Read and analyze Saga SHM mutation-spectrum outputs
############################################################

library(dplyr)
library(ggplot2)
library(tidyr)

############################################################
## 1. Set folder path
############################################################

shm_dir <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/SHM_mutation_spectrum_HK_Sp_3w"

full_file <- file.path(shm_dir, "mutation_spectrum_full.csv")
summary_file <- file.path(shm_dir, "mutation_spectrum_summary.csv")
###check
dir.exists(shm_dir)
list.files("C:/Users/sje054", pattern = "OneDrive", full.names = TRUE)

##Check end

file.exists(full_file)
file.exists(summary_file)

############################################################
## 2. Read files
############################################################

mutation_spectrum <- read.csv(full_file)

mutation_spectrum_summary <- read.csv(summary_file)

head(mutation_spectrum)
head(mutation_spectrum_summary)

str(mutation_spectrum)
str(mutation_spectrum_summary)

############################################################
## 3. Clean group order
############################################################

group_levels <- c("HKC", "HKI", "SpC", "SpI")

mutation_spectrum <- mutation_spectrum %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

mutation_spectrum_summary <- mutation_spectrum_summary %>%
  mutate(
    Group = factor(Group, levels = group_levels)
  )

############################################################
## 4. Relative mutation spectrum by group
############################################################

mut_rel_group <- mutation_spectrum_summary %>%
  group_by(Group) %>%
  mutate(
    total_weighted = sum(weighted_count, na.rm = TRUE),
    relative_freq = weighted_count / total_weighted,
    percent = 100 * relative_freq
  ) %>%
  ungroup()

#View(mut_rel_group)

############################################################
## 5. Plot absolute mutation spectrum
############################################################
p_abs <- ggplot(
  mutation_spectrum_summary,
  aes(x = Mutation, y = weighted_count, fill = Group)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ Group) +
  theme_bw(base_size = 13) +
  labs(
    title = "SHM Mutation Spectrum",
    subtitle = "Clone-weighted mutation counts",
    x = "Mutation type",
    y = "Weighted mutation count"
  )

p_abs

############################################################
## 6. Plot relative mutation spectrum
############################################################

p_rel <- ggplot(
  mut_rel_group,
  aes(x = Mutation, y = percent, fill = Group)
) +
  geom_col(width = 0.8) +
  facet_wrap(~ Group) +
  theme_bw(base_size = 13) +
  labs(
    title = "Relative SHM Mutation Spectrum",
    subtitle = "Percent of clone-weighted mutations within each group",
    x = "Mutation type",
    y = "Relative frequency (%)"
  )

p_rel

############################################################
## 7. Group comparison plot: all groups side-by-side
############################################################

p_compare <- ggplot(
  mut_rel_group,
  aes(x = Mutation, y = percent, fill = Group)
) +
  geom_col(position = "dodge", width = 0.8) +
  theme_bw(base_size = 13) +
  labs(
    title = "Comparison of SHM Mutation Spectra",
    x = "Mutation type",
    y = "Relative frequency (%)"
  )

p_compare

############################################################
## 8. Transition vs transversion summary
############################################################

transition_mutations <- c("A>G", "G>A", "C>T", "T>C")

mut_type_summary <- mutation_spectrum_summary %>%
  mutate(
    Mutation_class = ifelse(
      Mutation %in% transition_mutations,
      "Transition",
      "Transversion"
    )
  ) %>%
  group_by(Group, Mutation_class) %>%
  summarise(
    weighted_count = sum(weighted_count, na.rm = TRUE),
    raw_count = sum(raw_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(
    percent = 100 * weighted_count / sum(weighted_count)
  ) %>%
  ungroup()

#View(mut_type_summary)

p_transition <- ggplot(
  mut_type_summary,
  aes(x = Group, y = percent, fill = Mutation_class)
) +
  geom_col(position = "stack") +
  theme_bw(base_size = 13) +
  labs(
    title = "Transition vs Transversion Mutations",
    x = "",
    y = "Relative frequency (%)",
    fill = "Mutation class"
  )

p_transition

############################################################
## 9. Sample-level mutation spectrum
############################################################

sample_mut_summary <- mutation_spectrum %>%
  group_by(Sample, Group, Tissue, Condition, Mutation) %>%
  summarise(
    weighted_count = sum(weighted_count, na.rm = TRUE),
    raw_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Sample) %>%
  mutate(
    percent = 100 * weighted_count / sum(weighted_count, na.rm = TRUE)
  ) %>%
  ungroup()

#View(sample_mut_summary)

############################################################
## 10. Statistical test per mutation type: HKC vs HKI, SpC vs SpI
############################################################

test_mutation <- function(data, group1, group2) {
  
  data %>%
    filter(Group %in% c(group1, group2)) %>%
    group_by(Mutation) %>%
    summarise(
      pvalue = tryCatch(
        wilcox.test(
          percent ~ Group,
          data = cur_data(),
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      ),
      mean_group1 = mean(percent[Group == group1], na.rm = TRUE),
      mean_group2 = mean(percent[Group == group2], na.rm = TRUE),
      log2FC = log2((mean_group2 + 1e-6) / (mean_group1 + 1e-6)),
      .groups = "drop"
    ) %>%
    mutate(
      Comparison = paste(group2, "vs", group1),
      padj = p.adjust(pvalue, method = "BH")
    )
}

mutation_tests <- bind_rows(
  test_mutation(sample_mut_summary, "HKC", "HKI"),
  test_mutation(sample_mut_summary, "SpC", "SpI")
)

#View(mutation_tests)

############################################################
## 12. V-family mutation spectrum
############################################################

vfamily_mut_summary <- mutation_spectrum %>%
  group_by(Group, V_family, Mutation) %>%
  summarise(
    weighted_count = sum(weighted_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group, V_family) %>%
  mutate(
    percent = 100 * weighted_count / sum(weighted_count)
  ) %>%
  ungroup()

top_vfamilies <- mutation_spectrum %>%
  group_by(V_family) %>%
  summarise(
    total = sum(weighted_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  slice_max(total, n = 12)

vfamily_mut_top <- vfamily_mut_summary %>%
  filter(V_family %in% top_vfamilies$V_family)

p_vfamily <- ggplot(
  vfamily_mut_top,
  aes(x = Mutation, y = percent, fill = Group)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ V_family, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Mutation Spectrum by Top V Families",
    x = "Mutation type",
    y = "Relative frequency (%)"
  )

p_vfamily


###Option 1
###order the V families by abundance before plotting:
v_order <- top_vfamilies %>%
  arrange(desc(total)) %>%
  pull(V_family)

vfamily_mut_top$V_family <- factor(
  vfamily_mut_top$V_family,
  levels = v_order
)


##
p_vfamily <- ggplot(
  vfamily_mut_top,
  aes(x = Mutation, y = percent, fill = Group)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ V_family, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    )
  ) +
  labs(
    title = "Mutation Spectrum by Top Abundance V Families",
    x = "Mutation type",
    y = "Relative frequency (%)"
  )

p_vfamily


##Option three
v_order <- top_vfamilies %>%
  arrange(desc(total)) %>%
  pull(V_family)

vfamily_mut_top$V_family <- factor(
  vfamily_mut_top$V_family,
  levels = v_order
)

###
############################################################
## Same Y-axis across all V families
############################################################

ymax <- max(
  vfamily_mut_top$percent,
  na.rm = TRUE
)

p_vfamily <- ggplot(
  vfamily_mut_top,
  aes(
    x = Mutation,
    y = percent,
    fill = Group
  )
) +
  geom_col(position = "dodge") +
  facet_wrap(
    ~ V_family,
    ncol = 4,
    nrow = 3
  ) +
  coord_cartesian(
    ylim = c(0, ymax * 1.05)
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 7,
      vjust = 0.5,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.position = "bottom"
  ) +
  labs(
    title = "Mutation Spectrum by Top V Families",
    x = "",
    y = "Relative frequency (%)"
  )

p_vfamily


########################################################
#ENd of the SHM script
###########################




############################################################
## Do infection-specific public clones have elevated SHM?
##
## Required objects:
## 1. shm_results
##    Must contain: CDR3.aa, Sample, Group, Tissue, Condition,
##                  Clones, mutation_freq, mutations, compared_nt
##
## 2. HKI_public, SpI_public, HKC_public, SpC_public
##    Must contain: CDR3.aa
############################################################

library(dplyr)
library(ggplot2)
library(ggpubr)

############################################################
## 1. Define infection-specific and control-specific public clones
############################################################

infected_public <- union(
  HKI_public$CDR3.aa,
  SpI_public$CDR3.aa
)

control_public <- union(
  HKC_public$CDR3.aa,
  SpC_public$CDR3.aa
)

infected_specific_public <- setdiff(
  infected_public,
  control_public
)

control_specific_public <- setdiff(
  control_public,
  infected_public
)

shared_public <- intersect(
  infected_public,
  control_public
)

length(infected_specific_public)
length(control_specific_public)
length(shared_public)

############################################################
## 2. Annotate SHM results by clone category
############################################################

shm_clone_category <- shm_results %>%
  mutate(
    Public_category = case_when(
      CDR3.aa %in% infected_specific_public ~ "Infection-specific public",
      CDR3.aa %in% control_specific_public  ~ "Control-specific public",
      CDR3.aa %in% shared_public            ~ "Shared public",
      TRUE                                  ~ "Private"
    )
  ) %>%
  filter(Public_category != "Private")

shm_clone_category$Public_category <- factor(
  shm_clone_category$Public_category,
  levels = c(
    "Control-specific public",
    "Shared public",
    "Infection-specific public"
  )
)

table(shm_clone_category$Public_category)

############################################################
## 3. Per-sequence SHM comparison
############################################################

p_seq <- ggplot(
  shm_clone_category,
  aes(
    x = Public_category,
    y = mutation_freq,
    fill = Public_category
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  stat_compare_means(method = "kruskal.test") +
  stat_compare_means(
    comparisons = list(
      c("Control-specific public", "Infection-specific public"),
      c("Shared public", "Infection-specific public")
    ),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "SHM in Infection-Specific Public Clones",
    x = "",
    y = "V-region mutation frequency"
  )

p_seq

############################################################
## 4. Sample-level SHM comparison
## More conservative: each fish contributes one value
############################################################

sample_public_shm <- shm_clone_category %>%
  group_by(Sample, Group, Tissue, Condition, Public_category) %>%
  summarise(
    total_mutations = sum(mutations * Clones, na.rm = TRUE),
    total_compared_nt = sum(compared_nt * Clones, na.rm = TRUE),
    weighted_mutation_freq = total_mutations / total_compared_nt,
    n_CDR3 = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_compared_nt > 0)

#View(sample_public_shm)

p_sample <- ggplot(
  sample_public_shm,
  aes(
    x = Public_category,
    y = weighted_mutation_freq,
    fill = Public_category
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2.5) +
  stat_compare_means(method = "kruskal.test") +
  stat_compare_means(
    comparisons = list(
      c("Control-specific public", "Infection-specific public"),
      c("Shared public", "Infection-specific public")
    ),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Sample-Level SHM in Public Clone Categories",
    subtitle = "Clone-weighted mutation frequency per fish",
    x = "",
    y = "Weighted mutation frequency"
  )

p_sample

############################################################
## 5. Infection-specific public clones by tissue
############################################################

infected_specific_shm <- shm_results %>%
  filter(CDR3.aa %in% infected_specific_public)

infected_specific_sample <- infected_specific_shm %>%
  group_by(Sample, Group, Tissue, Condition) %>%
  summarise(
    total_mutations = sum(mutations * Clones, na.rm = TRUE),
    total_compared_nt = sum(compared_nt * Clones, na.rm = TRUE),
    weighted_mutation_freq = total_mutations / total_compared_nt,
    n_infection_specific_public = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_compared_nt > 0)

p_infected_tissue <- ggplot(
  infected_specific_sample,
  aes(
    x = Group,
    y = weighted_mutation_freq,
    fill = Group
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2.5) +
  theme_bw(base_size = 13) +
  labs(
    title = "SHM of Infection-Specific Public Clones",
    x = "",
    y = "Weighted mutation frequency"
  )

p_infected_tissue

############################################################
## 6. Statistical tests
############################################################

## Overall public-category comparison
kruskal.test(
  weighted_mutation_freq ~ Public_category,
  data = sample_public_shm
)

pairwise.wilcox.test(
  sample_public_shm$weighted_mutation_freq,
  sample_public_shm$Public_category,
  p.adjust.method = "BH"
)

## Direct test:
## infection-specific public vs control-specific public
wilcox.test(
  weighted_mutation_freq ~ Public_category,
  data = sample_public_shm %>%
    filter(
      Public_category %in% c(
        "Control-specific public",
        "Infection-specific public"
      )
    ),
  exact = FALSE
)

############################################################
## 7. Summary table for reporting
############################################################

public_shm_summary <- sample_public_shm %>%
  group_by(Public_category) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_SHM = mean(weighted_mutation_freq, na.rm = TRUE),
    median_SHM = median(weighted_mutation_freq, na.rm = TRUE),
    SD_SHM = sd(weighted_mutation_freq, na.rm = TRUE),
    mean_n_CDR3 = mean(n_CDR3, na.rm = TRUE),
    total_CDR3 = sum(n_CDR3, na.rm = TRUE),
    .groups = "drop"
  )

public_shm_summary

############################################################
## 8. V-family SHM among infection-specific public clones
############################################################

infected_public_vfamily_shm <- infected_specific_shm %>%
  group_by(Group, V_family) %>%
  summarise(
    total_mutations = sum(mutations * Clones, na.rm = TRUE),
    total_compared_nt = sum(compared_nt * Clones, na.rm = TRUE),
    weighted_mutation_freq = total_mutations / total_compared_nt,
    n_CDR3 = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_compared_nt > 0)

top_vfamilies_infected_public <- infected_public_vfamily_shm %>%
  group_by(V_family) %>%
  summarise(
    total_clones = sum(total_clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  slice_max(total_clones, n = 15)

infected_public_vfamily_top <- infected_public_vfamily_shm %>%
  filter(V_family %in% top_vfamilies_infected_public$V_family)

p_vfamily_shm <- ggplot(
  infected_public_vfamily_top,
  aes(
    x = reorder(V_family, weighted_mutation_freq),
    y = weighted_mutation_freq,
    fill = Group
  )
) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "SHM of Infection-Specific Public Clones by V Family",
    x = "V family",
    y = "Weighted mutation frequency"
  )

p_vfamily_shm

############################################################


############################################################
## TRUE-LIKE BCR CLONAL LINEAGE RECONSTRUCTION
## Definition:
##   Same V gene/family
##   Same J gene/family
##   Similar CDR3.aa
##
## Input required:
##   imm_ds
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringdist)
library(igraph)
library(ggraph)
library(tidygraph)

############################################################
## 1. Build full clone table
############################################################

all_clones_full <- bind_rows(
  lapply(names(imm_ds), function(s) {
    imm_ds[[s]] %>%
      select(
        CDR3.aa,
        CDR3.nt,
        Clones,
        Proportion,
        V.name,
        J.name
      ) %>%
      mutate(Sample = s)
  })
) %>%
  mutate(
    Sample_clean = gsub("\\.clones_IGH", "", Sample),
    Group = sub("_.*", "", Sample_clean),
    Tissue = ifelse(grepl("^HK", Group), "HK", "Spleen"),
    Condition = ifelse(grepl("I$", Group), "Infected", "Control"),
    
    V_gene_no_allele = sub("\\*.*", "", V.name),
    J_gene_no_allele = sub("\\*.*", "", J.name),
    
    #V_family = sub("-.*", "", V_gene_no_allele),
    V_family = sub(
      "-.*$",
      "",
      sub("S[0-9]+$", "", V_gene_no_allele)
    ),
    
    #J_family = sub("-.*", "", J_gene_no_allele),
    J_family = sub(
      "-.*$",
      "",
      sub("T[0-9]+(?=D?$)", "", J_gene_no_allele, perl = TRUE)
    ),
    CDR3_length = nchar(CDR3.aa),
    
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  ) %>%
  filter(
    !is.na(CDR3.aa),
    !is.na(V_family),
    !is.na(J_family)
  )

############################################################
## 2. Collapse identical CDR3s within each sample
############################################################

clone_table <- all_clones_full %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    V_family,
    J_family,
    CDR3.aa,
    CDR3.nt,
    CDR3_length
  ) %>%
  summarise(
    Clones = sum(Clones, na.rm = TRUE),
    Proportion = sum(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 3. Define candidate lineage blocks
## Same V family + J family + similar CDR3 length
############################################################

clone_table <- clone_table %>%
  mutate(
    BlockID = paste(
      V_family,
      J_family,
      CDR3_length,
      sep = "_"
    )
  )

############################################################
## 4. Function to build lineages within each block
############################################################
build_lineages_one_block <- function(df, max_aa_distance = 1) {
  
  seqs <- unique(df$CDR3.aa)
  
  if (length(seqs) == 1) {
    df$LineageID <- paste0(df$BlockID[1], "_L1")
    return(df)
  }
  
  dmat <- stringdist::stringdistmatrix(
    seqs,
    seqs,
    method = "lv"
  )
  
  edges <- which(
    dmat <= max_aa_distance & dmat > 0,
    arr.ind = TRUE
  ) %>%
    as.data.frame() %>%
    filter(row < col)
  
  if (nrow(edges) == 0) {
    
    lineage_map <- data.frame(
      CDR3.aa = seqs,
      Component = seq_along(seqs)
    )
    
  } else {
    
    edge_df <- data.frame(
      from = seqs[edges$row],
      to = seqs[edges$col]
    )
    
    g <- igraph::graph_from_data_frame(
      edge_df,
      vertices = data.frame(name = seqs),
      directed = FALSE
    )
    
    comps <- igraph::components(g)
    
    lineage_map <- data.frame(
      CDR3.aa = names(comps$membership),
      Component = comps$membership
    )
  }
  
  lineage_map <- lineage_map %>%
    mutate(
      LineageID = paste0(df$BlockID[1], "_L", Component)
    ) %>%
    select(CDR3.aa, LineageID)
  
  df %>%
    left_join(lineage_map, by = "CDR3.aa")
}

############################################################
## 5. Build all lineages
## max_aa_distance = 1 is strict
## Try 2 for a more relaxed definition
############################################################

max_aa_distance <- 1

lineage_df <- clone_table %>%
  group_by(BlockID) %>%
  group_modify(~ build_lineages_one_block(.x, max_aa_distance)) %>%
  ungroup()

View(lineage_df)

############################################################
## 6. Lineage-level summary
############################################################
lineage_summary <- lineage_df %>%
  group_by(LineageID) %>%
  summarise(
    V_family = dplyr::first(V_family),
    J_family = dplyr::first(J_family),
    CDR3_length = dplyr::first(CDR3_length),
    
    n_CDR3_aa = n_distinct(CDR3.aa),
    n_CDR3_nt = n_distinct(CDR3.nt),
    n_samples = n_distinct(Sample),
    n_groups = n_distinct(Group),
    
    total_clones = sum(Clones, na.rm = TRUE),
    mean_clone_size = mean(Clones, na.rm = TRUE),
    max_clone_size = max(Clones, na.rm = TRUE),
    
    groups_present = paste(sort(unique(as.character(Group))), collapse = ";"),
    samples_present = paste(sort(unique(Sample)), collapse = ";"),
    
    .groups = "drop"
  ) %>%
  arrange(desc(total_clones))






View(lineage_summary)

############################################################
## 7. Sample-level lineage metrics
############################################################

sample_lineage_metrics <- lineage_df %>%
  group_by(Sample, Group, Tissue, Condition, LineageID) %>%
  summarise(
    lineage_clones = sum(Clones, na.rm = TRUE),
    n_CDR3_aa = n_distinct(CDR3.aa),
    n_CDR3_nt = n_distinct(CDR3.nt),
    .groups = "drop"
  ) %>%
  group_by(Sample, Group, Tissue, Condition) %>%
  summarise(
    n_lineages = n_distinct(LineageID),
    mean_lineage_size = mean(lineage_clones, na.rm = TRUE),
    median_lineage_size = median(lineage_clones, na.rm = TRUE),
    max_lineage_size = max(lineage_clones, na.rm = TRUE),
    
    mean_lineage_aa_diversity = mean(n_CDR3_aa, na.rm = TRUE),
    mean_lineage_nt_diversity = mean(n_CDR3_nt, na.rm = TRUE),
    
    expanded_lineages = sum(lineage_clones >= 10, na.rm = TRUE),
    
    .groups = "drop"
  )

View(sample_lineage_metrics)

############################################################
## 8. Compare lineage metrics across groups
############################################################

metrics_long <- sample_lineage_metrics %>%
  pivot_longer(
    cols = c(
      n_lineages,
      mean_lineage_size,
      median_lineage_size,
      max_lineage_size,
      mean_lineage_aa_diversity,
      mean_lineage_nt_diversity,
      expanded_lineages
    ),
    names_to = "Metric",
    values_to = "Value"
  )

p_lineage_metrics <- ggplot(
  metrics_long,
  aes(x = Group, y = Value, fill = Group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "Clonal Lineage Architecture",
    x = "",
    y = "Value"
  )

p_lineage_metrics


############################################################
## 8. Exact values underlying the lineage-metrics figure
############################################################

library(dplyr)
library(tidyr)

metric_columns <- c(
  "n_lineages",
  "mean_lineage_size",
  "median_lineage_size",
  "max_lineage_size",
  "mean_lineage_aa_diversity",
  "mean_lineage_nt_diversity",
  "expanded_lineages"
)

############################################################
## 8.1 Exact value for every sample
############################################################

lineage_exact_values <- sample_lineage_metrics %>%
  select(
    Sample,
    Group,
    all_of(metric_columns)
  ) %>%
  arrange(Group, Sample)

print(lineage_exact_values, n = Inf, width = Inf)

# Optional: open as a spreadsheet-style table in RStudio
View(lineage_exact_values)





############################################################
## 9. Statistics for lineage metrics
############################################################

lineage_metric_tests <- metrics_long %>%
  group_by(Metric) %>%
  summarise(
    kruskal_p = kruskal.test(Value ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(kruskal_p, method = "BH")
  )

lineage_metric_tests

pairwise_lineage_tests <- bind_rows(
  lapply(unique(metrics_long$Metric), function(m) {
    
    df <- metrics_long %>% filter(Metric == m)
    
    hk_p <- tryCatch(
      wilcox.test(
        Value ~ Group,
        data = df %>% filter(Group %in% c("HKC", "HKI")),
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    
    sp_p <- tryCatch(
      wilcox.test(
        Value ~ Group,
        data = df %>% filter(Group %in% c("SpC", "SpI")),
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    
    data.frame(
      Metric = m,
      HKC_vs_HKI_p = hk_p,
      SpC_vs_SpI_p = sp_p
    )
  })
) %>%
  mutate(
    HKC_vs_HKI_padj = p.adjust(HKC_vs_HKI_p, method = "BH"),
    SpC_vs_SpI_padj = p.adjust(SpC_vs_SpI_p, method = "BH")
  )

pairwise_lineage_tests

############################################################
## 10. Identify infection-associated lineages
############################################################

infected_lineages <- lineage_df %>%
  filter(Condition == "Infected") %>%
  pull(LineageID) %>%
  unique()

control_lineages <- lineage_df %>%
  filter(Condition == "Control") %>%
  pull(LineageID) %>%
  unique()

infection_specific_lineages <- setdiff(
  infected_lineages,
  control_lineages
)

control_specific_lineages <- setdiff(
  control_lineages,
  infected_lineages
)

shared_lineages <- intersect(
  infected_lineages,
  control_lineages
)

lineage_summary <- lineage_summary %>%
  mutate(
    Lineage_category = case_when(
      LineageID %in% infection_specific_lineages ~ "Infection-specific",
      LineageID %in% control_specific_lineages ~ "Control-specific",
      LineageID %in% shared_lineages ~ "Shared",
      TRUE ~ "Private"
    )
  )

lineage_category_counts <- lineage_summary %>%
  count(Lineage_category)

lineage_category_counts

ggplot(
  lineage_category_counts,
  aes(x = Lineage_category, y = n, fill = Lineage_category)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
    ) +
  labs(
    title = "Lineage Categories",
    x = "",
    y = "Number of lineages"
  )

############################################################
## 11. Top expanded lineages
############################################################

top_lineages <- lineage_summary %>%
  slice_max(total_clones, n = 20)

ggplot(
  top_lineages,
  aes(
    x = reorder(LineageID, total_clones),
    y = total_clones,
    fill = Lineage_category
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Top Expanded BCR Lineages",
    x = "Lineage",
    y = "Total clone count"
  )

############################################################
## 12. Build lineage trees / networks
## For the top infection-specific lineages
############################################################

build_lineage_graph <- function(lineage_id) {
  
  df <- lineage_df %>%
    filter(LineageID == lineage_id) %>%
    group_by(CDR3.aa, Group, Tissue, Condition) %>%
    summarise(
      Clones = sum(Clones, na.rm = TRUE),
      .groups = "drop"
    )
  
  seqs <- df$CDR3.aa
  
  if (length(seqs) < 2) return(NULL)
  
  dmat <- stringdistmatrix(seqs, seqs, method = "lv")
  
  g_full <- graph_from_adjacency_matrix(
    as.matrix(dmat),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  mst_g <- mst(g_full)
  
  V(mst_g)$name <- seqs
  
  node_meta <- df %>%
    select(CDR3.aa, Group, Tissue, Condition, Clones)
  
  V(mst_g)$Group <- node_meta$Group[
    match(V(mst_g)$name, node_meta$CDR3.aa)
  ]
  
  V(mst_g)$Condition <- node_meta$Condition[
    match(V(mst_g)$name, node_meta$CDR3.aa)
  ]
  
  V(mst_g)$Clones <- node_meta$Clones[
    match(V(mst_g)$name, node_meta$CDR3.aa)
  ]
  
  mst_g
}

top_infection_lineages <- lineage_summary %>%
  filter(Lineage_category == "Infection-specific") %>%
  slice_max(total_clones, n = 6) %>%
  pull(LineageID)

top_infection_lineages

############################################################
## 13. Plot lineage networks
############################################################
#install.packages(
 # c("ggraph", "tidygraph"),
 # repos = "https://cloud.r-project.org"
#)

library(ggraph)
library(tidygraph)


for (lin in top_infection_lineages) {
  
  g <- build_lineage_graph(lin)
  
  if (is.null(g)) next
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.5) +
    geom_node_point(
      aes(
        color = Group,
        size = Clones
      ),
      alpha = 0.9
    ) +
    theme_void() +
    labs(
      title = paste("Lineage tree:", lin)
    )
  
  print(p)
  
  ggsave(
    paste0("Lineage_tree_", gsub("[^A-Za-z0-9]", "_", lin), ".png"),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

############################################################
## 14. Network topology metrics per lineage
############################################################

lineage_network_metrics <- bind_rows(
  lapply(unique(lineage_df$LineageID), function(lin) {
    
    g <- build_lineage_graph(lin)
    
    if (is.null(g)) {
      return(
        data.frame(
          LineageID = lin,
          n_nodes = 1,
          n_edges = 0,
          mean_degree = 0,
          max_degree = 0,
          diameter = 0
        )
      )
    }
    
    data.frame(
      LineageID = lin,
      n_nodes = vcount(g),
      n_edges = ecount(g),
      mean_degree = mean(degree(g)),
      max_degree = max(degree(g)),
      diameter = diameter(g, weights = NA)
    )
  })
) %>%
  left_join(
    lineage_summary %>%
      select(
        LineageID,
        Lineage_category,
        total_clones,
        n_samples,
        groups_present
      ),
    by = "LineageID"
  )

View(lineage_network_metrics)

############################################################
## 15. Compare topology by lineage category
############################################################

network_long <- lineage_network_metrics %>%
  pivot_longer(
    cols = c(
      n_nodes,
      n_edges,
      mean_degree,
      max_degree,
      diameter,
      total_clones,
      n_samples
    ),
    names_to = "Metric",
    values_to = "Value"
  )

ggplot(
  network_long,
  aes(
    x = Lineage_category,
    y = Value,
    fill = Lineage_category
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Lineage Network Topology",
    x = "",
    y = "Value"
  )

############################################################
##

#######################################################NEW
###New - lineage tree for top expanded 5 BCR
############################################################
## Top 5 most expanded BCR lineages
## Expansion defined by total clone abundance
############################################################

library(dplyr)
library(igraph)
library(stringdist)
library(ggraph)
library(tidygraph)
library(ggplot2)

############################################################
## 1. Build one lineage network
############################################################

build_lineage_graph <- function(lineage_id) {
  
  lineage_raw <- lineage_df %>%
    filter(LineageID == lineage_id) %>%
    filter(
      !is.na(CDR3.aa),
      CDR3.aa != "",
      !is.na(Clones)
    )
  
  if (nrow(lineage_raw) == 0) {
    return(NULL)
  }
  
  ## One node per unique CDR3 sequence
  ## Node colour is assigned from the group contributing
  ## the greatest clone abundance for that CDR3.
  node_meta <- lineage_raw %>%
    group_by(CDR3.aa, Group, Tissue, Condition) %>%
    summarise(
      group_clones = sum(Clones, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(CDR3.aa) %>%
    arrange(desc(group_clones), .by_group = TRUE) %>%
    summarise(
      Clones = sum(group_clones, na.rm = TRUE),
      
      Dominant_group = dplyr::first(Group),
      Dominant_tissue = dplyr::first(Tissue),
      Dominant_condition = dplyr::first(Condition),
      
      Groups_present = paste(
        sort(unique(as.character(Group))),
        collapse = ";"
      ),
      
      n_groups = n_distinct(Group),
      
      .groups = "drop"
    ) %>%
    mutate(
      Node_group = ifelse(
        n_groups > 1,
        "Shared across groups",
        as.character(Dominant_group)
      )
    )
  
  seqs <- node_meta$CDR3.aa
  
  if (length(seqs) < 2) {
    return(NULL)
  }
  
  ##########################################################
  ## Pairwise amino-acid Levenshtein distances
  ##########################################################
  
  dmat <- stringdistmatrix(
    seqs,
    seqs,
    method = "lv"
  )
  
  rownames(dmat) <- seqs
  colnames(dmat) <- seqs
  
  ##########################################################
  ## Complete weighted graph followed by MST
  ##########################################################
  
  g_full <- graph_from_adjacency_matrix(
    as.matrix(dmat),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  mst_g <- mst(
    g_full,
    weights = E(g_full)$weight
  )
  
  ##########################################################
  ## Add node metadata
  ##########################################################
  
  V(mst_g)$name <- seqs
  
  idx <- match(
    V(mst_g)$name,
    node_meta$CDR3.aa
  )
  
  V(mst_g)$Clones <- node_meta$Clones[idx]
  V(mst_g)$Group <- node_meta$Node_group[idx]
  V(mst_g)$Dominant_group <- node_meta$Dominant_group[idx]
  V(mst_g)$Tissue <- node_meta$Dominant_tissue[idx]
  V(mst_g)$Condition <- node_meta$Dominant_condition[idx]
  V(mst_g)$Groups_present <- node_meta$Groups_present[idx]
  
  mst_g
}

############################################################
## 2. Identify the five most expanded lineages
############################################################
## This ranks all lineages, irrespective of whether they are
## shared, infection-specific, or control-specific.

top_expanded_lineages <- lineage_summary %>%
  filter(
    !is.na(LineageID),
    !is.na(total_clones)
  ) %>%
  arrange(desc(total_clones)) %>%
  slice_head(n = 5) %>%
  pull(LineageID)

top_expanded_lineages

## Optional summary table
top_expanded_lineage_summary <- lineage_summary %>%
  filter(LineageID %in% top_expanded_lineages) %>%
  arrange(desc(total_clones))

print(top_expanded_lineage_summary)

write.csv(
  top_expanded_lineage_summary,
  "Top_5_expanded_lineages_summary.csv",
  row.names = FALSE
)

############################################################
## 3. Plot the five lineage networks separately
############################################################

top_lineage_plots <- list()

for (lin in top_expanded_lineages) {
  
  g <- build_lineage_graph(lin)
  
  if (is.null(g)) {
    message(
      "Skipping ", lin,
      ": fewer than two unique CDR3 sequences."
    )
    next
  }
  
  lineage_total <- lineage_summary %>%
    filter(LineageID == lin) %>%
    pull(total_clones)
  
  lineage_category <- lineage_summary %>%
    filter(LineageID == lin) %>%
    pull(Lineage_category)
  
  p <- ggraph(
    g,
    layout = "fr"
  ) +
    geom_edge_link(
      aes(
        width = 1 / weight
      ),
      alpha = 0.45
    ) +
    geom_node_point(
      aes(
        color = Group,
        size = Clones
      ),
      alpha = 0.90
    ) +
    scale_edge_width(
      range = c(0.2, 1.3),
      guide = "none"
    ) +
    scale_size_continuous(
      range = c(3, 14)
    ) +
    theme_void(base_size = 13) +
    labs(
      title = paste(
        "Expanded BCR lineage:",
        lin
      ),
      subtitle = paste0(
        "Total clone abundance = ",
        lineage_total,
        "; category = ",
        lineage_category
      ),
      color = "Group",
      size = "Clone abundance"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  print(p)
  
  top_lineage_plots[[lin]] <- p
  
  ggsave(
    filename = paste0(
      "Top_expanded_lineage_",
      gsub("[^A-Za-z0-9]", "_", lin),
      ".png"
    ),
    plot = p,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}

############################################################
## 4. Optional: combine the five trees into one figure
############################################################

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

library(patchwork)

if (length(top_lineage_plots) > 0) {
  
  p_top5_combined <- wrap_plots(
    top_lineage_plots,
    ncol = 2,
    guides = "collect"
  ) +
    plot_annotation(
      title = "Five most expanded BCR lineages",
      subtitle = paste(
        "Node size represents clone abundance;",
        "node colour represents group origin"
      )
    ) &
    theme(
      legend.position = "bottom"
    )
  
  print(p_top5_combined)
  
  ggsave(
    "Top_5_expanded_BCR_lineages_combined.png",
    p_top5_combined,
    width = 15,
    height = 17,
    dpi = 300
  )
}

############################################################
## 5. Network metrics for only the top five lineages
############################################################

top5_lineage_network_metrics <- bind_rows(
  lapply(
    top_expanded_lineages,
    function(lin) {
      
      g <- build_lineage_graph(lin)
      
      if (is.null(g)) {
        return(
          data.frame(
            LineageID = lin,
            n_nodes = 1,
            n_edges = 0,
            mean_degree = 0,
            max_degree = 0,
            diameter = 0,
            mean_edge_distance = NA_real_
          )
        )
      }
      
      data.frame(
        LineageID = lin,
        n_nodes = vcount(g),
        n_edges = ecount(g),
        mean_degree = mean(degree(g)),
        max_degree = max(degree(g)),
        diameter = diameter(
          g,
          weights = NA
        ),
        mean_edge_distance = mean(
          E(g)$weight,
          na.rm = TRUE
        )
      )
    }
  )
) %>%
  left_join(
    lineage_summary %>%
      select(
        LineageID,
        Lineage_category,
        total_clones,
        n_samples,
        groups_present
      ),
    by = "LineageID"
  ) %>%
  arrange(desc(total_clones))

print(top5_lineage_network_metrics)

write.csv(
  top5_lineage_network_metrics,
  "Top_5_expanded_lineage_network_metrics.csv",
  row.names = FALSE
)




##START NEW
#####Complete REGION analysis 
###Based on the correct samples

###0. Read the Saga outputs
############################################################
## 0. Load SAV SHM region analysis
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

shm_dir <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/SHM_mutation_spectrum_HK_Sp_3w"

region_rs <- readRDS(
  file.path(shm_dir,
            "region_rs_SAV_3wpi.rds")
)

sample_region_rs <-
  read.csv(file.path(
    shm_dir,
    "sample_region_rs_SAV_3wpi.csv"
  ))

sample_region_class_rs <-
  read.csv(file.path(
    shm_dir,
    "sample_region_class_rs_SAV_3wpi.csv"
  ))

table(region_rs$Group)
table(region_rs$Region)

############################################################
## 1. Overall SHM summary
############################################################

region_summary <-
  sample_region_rs %>%
  group_by(Group, Region) %>%
  summarise(
    
    mean_RS = mean(RS_ratio),
    
    sd_RS = sd(RS_ratio),
    
    median_RS = median(RS_ratio),
    
    .groups="drop"
    
  )

region_summary

############################################################
## 2. Region-specific R/S ratio
############################################################

p1 <-
  
  ggplot(
    sample_region_rs,
    
    aes(
      Group,
      RS_ratio,
      fill=Group
    )
  )+
  
  geom_boxplot(
    outlier.shape=NA,
    alpha=.8
  )+
  
  geom_jitter(
    width=.15,
    size=2,
    alpha=.8
  )+
  
  facet_wrap(
    ~Region,
    scales="free_y"
  )+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Replacement/Silent Mutation Ratio across V-gene Regions",
    
    x="",
    
    y="Weighted R/S ratio"
    
  )

p1

############################################################
## 3. CDR vs Framework
############################################################

p2 <-
  
  ggplot(
    
    sample_region_class_rs,
    
    aes(
      
      Region_class,
      
      RS_ratio,
      
      fill=Region_class
      
    )
    
  )+
  
  geom_boxplot(
    
    outlier.shape=NA,
    
    alpha=.8
    
  )+
  
  geom_jitter(
    
    width=.15,
    
    size=2,
    
    alpha=.8
    
  )+
  
  facet_wrap(
    
    ~Group
    
  )+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Replacement/Silent Ratio in CDR and Framework Regions",
    
    x="",
    
    y="Weighted R/S ratio"
    
  )

p2

############################################################
## 4. Infection effect
############################################################

p3 <-
  
  ggplot(
    
    sample_region_class_rs,
    
    aes(
      
      Condition,
      
      RS_ratio,
      
      fill=Condition
      
    )
    
  )+
  
  geom_boxplot(
    
    outlier.shape=NA,
    
    alpha=.8
    
  )+
  
  geom_jitter(
    
    width=.15,
    
    size=2,
    
    alpha=.8
    
  )+
  
  facet_wrap(
    
    ~Tissue
    
  )+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Effect of Infection on SHM Selection",
    
    x="",
    
    y="Weighted R/S ratio"
    
  )

p3

############################################################
## 5. Tissue effect
############################################################

ggplot(
  
  sample_region_class_rs,
  
  aes(
    
    Tissue,
    
    RS_ratio,
    
    fill=Tissue
    
  )
  
)+
  
  geom_boxplot(
    
    outlier.shape=NA,
    
    alpha=.8
    
  )+
  
  geom_jitter(
    
    width=.15,
    
    size=2
    
  )+
  
  facet_wrap(
    
    ~Condition
    
  )+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Tissue-specific SHM",
    
    y="Weighted R/S ratio"
    
  )


############################################################
## 6. Statistics
############################################################

kruskal.test(
  
  RS_ratio~Group,
  
  data=sample_region_class_rs
  
)

pairwise.wilcox.test(
  
  sample_region_class_rs$RS_ratio,
  
  sample_region_class_rs$Group,
  
  p.adjust.method="BH"
  
)

############################################################
## 7. Tissue x Infection
############################################################

lm1 <-
  
  lm(
    
    RS_ratio~
      
      Tissue*
      
      Condition,
    
    data=sample_region_class_rs
    
  )

anova(lm1)

summary(lm1)

############################################################
## 8. V-family SHM
############################################################

vfamily_rs <-
  
  region_rs %>%
  
  group_by(
    
    Group,
    
    V_family
    
  )%>%
  
  summarise(
    
    RS_ratio=
      
      sum(R_weighted)/(sum(S_weighted)+1),
    
    .groups="drop"
    
  )
###Heatmap
vmat<-
  
  vfamily_rs%>%
  
  pivot_wider(
    
    names_from=Group,
    
    values_from=RS_ratio,
    
    values_fill=0
    
  )%>%
  
  column_to_rownames("V_family")

pheatmap(
  
  vmat,
  
  cluster_rows=TRUE,
  
  cluster_cols=TRUE,
  
  color=colorRampPalette(
    
    c("navy","white","firebrick3")
    
  )(100),
  
  main="Selection across V Families"
  
)

###Which regions experience strongest selection?
############################################################
## 9. Compare FR1/CDR1/FR2/CDR2/FR3
############################################################

ggplot(
  
  sample_region_rs,
  
  aes(
    
    Region,
    
    RS_ratio,
    
    fill=Region
    
  )
  
)+
  
  geom_boxplot(
    
    alpha=.8
    
  )+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Selection Strength across V-gene Regions",
    
    y="Weighted R/S ratio"
    
  )

############################################################
## ADVANCED R/S SELECTION ANALYSIS
## SAV 3 wpi HK/Sp BCR dataset
##
## Questions:
## 1. Do infection-specific public clones have higher R/S ratios?
## 2. Which V families show strongest antigenic selection?
## 3. Does spleen accumulate more CDR replacement mutations than HK?
## 4. Which IGHV families show infection-specific positive selection?
## 5. Candidate antiviral IGHV families:
##    V-family enrichment + high CDR R/S ratio
############################################################

############################################################
## 1. Load files
############################################################

shm_dir <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/SHM_mutation_spectrum_HK_Sp_3w"

region_rs <- readRDS(
  file.path(shm_dir, "region_rs_SAV_3wpi.rds")
)

sample_region_rs <- read.csv(
  file.path(shm_dir, "sample_region_rs_SAV_3wpi.csv")
)

sample_region_class_rs <- read.csv(
  file.path(shm_dir, "sample_region_class_rs_SAV_3wpi.csv")
)

region_rs <- region_rs %>%
  mutate(
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected")),
    Region_class = factor(Region_class, levels = c("FR", "CDR"))
  )

table(region_rs$Group)
table(region_rs$Tissue)
table(region_rs$Condition)

############################################################
## 2. OPTIONAL: Load public clone objects
## Required for infection-specific public clone analysis
##
## These objects must already exist:
##   HKC_public, HKI_public, SpC_public, SpI_public
## each containing CDR3.aa
############################################################

if (
  exists("HKC_public") &&
  exists("HKI_public") &&
  exists("SpC_public") &&
  exists("SpI_public")
) {
  
  infected_public <- union(HKI_public$CDR3.aa, SpI_public$CDR3.aa)
  control_public  <- union(HKC_public$CDR3.aa, SpC_public$CDR3.aa)
  
  infected_specific_public <- setdiff(infected_public, control_public)
  control_specific_public  <- setdiff(control_public, infected_public)
  shared_public            <- intersect(infected_public, control_public)
  
} else {
  
  warning(
    "Public clone objects not found. Public/private analysis will be skipped.
     Load HKC_public, HKI_public, SpC_public, and SpI_public first."
  )
}

############################################################
## 3. Do infection-specific public clones have higher R/S?
############################################################

if (exists("infected_specific_public")) {
  
  public_rs <- region_rs %>%
    mutate(
      Public_category = case_when(
        CDR3.aa %in% infected_specific_public ~ "Infection-specific public",
        CDR3.aa %in% control_specific_public ~ "Control-specific public",
        CDR3.aa %in% shared_public ~ "Shared public",
        TRUE ~ "Private"
      )
    )
  
  public_cdr_rs <- public_rs %>%
    filter(Region_class == "CDR") %>%
    group_by(Sample, Group, Tissue, Condition, Public_category) %>%
    summarise(
      total_R = sum(R_weighted, na.rm = TRUE),
      total_S = sum(S_weighted, na.rm = TRUE),
      CDR_RS_ratio = total_R / (total_S + 1),
      R_fraction = total_R / (total_R + total_S + 1),
      n_CDR3 = n_distinct(CDR3.aa),
      total_clones = sum(Clones, na.rm = TRUE),
      .groups = "drop"
    )
  
  public_cdr_rs$Public_category <- factor(
    public_cdr_rs$Public_category,
    levels = c(
      "Private",
      "Control-specific public",
      "Shared public",
      "Infection-specific public"
    )
  )
  
  public_cdr_rs_summary <- public_cdr_rs %>%
    group_by(Public_category) %>%
    summarise(
      n_samples = n_distinct(Sample),
      mean_RS = mean(CDR_RS_ratio, na.rm = TRUE),
      median_RS = median(CDR_RS_ratio, na.rm = TRUE),
      sd_RS = sd(CDR_RS_ratio, na.rm = TRUE),
      mean_R_fraction = mean(R_fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(public_cdr_rs_summary)
  
  public_cdr_rs_tests <- public_cdr_rs %>%
    filter(Public_category != "Private") %>%
    group_by(Public_category) %>%
    summarise(n = n(), .groups = "drop")
  
  print(public_cdr_rs_tests)
  
  p_public_rs <- ggplot(
    public_cdr_rs,
    aes(
      x = Public_category,
      y = CDR_RS_ratio,
      fill = Public_category
    )
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = "CDR R/S Ratio by Public Clone Category",
      x = "",
      y = "CDR R/S ratio"
    )
  
  print(p_public_rs)
  
  public_kruskal <- kruskal.test(
    CDR_RS_ratio ~ Public_category,
    data = public_cdr_rs
  )
  
  print(public_kruskal)
  
  public_pairwise <- pairwise.wilcox.test(
    public_cdr_rs$CDR_RS_ratio,
    public_cdr_rs$Public_category,
    p.adjust.method = "BH"
  )
  
  print(public_pairwise)
  
  write.csv(
    public_cdr_rs,
    file.path(shm_dir, "Selection_public_clone_CDR_RS.csv"),
    row.names = FALSE
  )
  
  write.csv(
    public_cdr_rs_summary,
    file.path(shm_dir, "Selection_public_clone_CDR_RS_summary.csv"),
    row.names = FALSE
  )
}

############################################################
## 4. Which V families are under strongest selection?
############################################################

vfamily_cdr_rs <- region_rs %>%
  filter(Region_class == "CDR") %>%
  group_by(Sample, Group, Tissue, Condition, V_family) %>%
  summarise(
    total_R = sum(R_weighted, na.rm = TRUE),
    total_S = sum(S_weighted, na.rm = TRUE),
    CDR_RS_ratio = total_R / (total_S + 1),
    R_fraction = total_R / (total_R + total_S + 1),
    total_clones = sum(Clones, na.rm = TRUE),
    n_CDR3 = n_distinct(CDR3.aa),
    .groups = "drop"
  )

vfamily_selection_summary <- vfamily_cdr_rs %>%
  group_by(V_family) %>%
  summarise(
    mean_CDR_RS = mean(CDR_RS_ratio, na.rm = TRUE),
    median_CDR_RS = median(CDR_RS_ratio, na.rm = TRUE),
    mean_R_fraction = mean(R_fraction, na.rm = TRUE),
    total_clones = sum(total_clones, na.rm = TRUE),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_CDR_RS))

top_selected_vfamilies <- vfamily_selection_summary %>%
  filter(n_samples >= 3) %>%
  slice_max(mean_CDR_RS, n = 15)

print(top_selected_vfamilies)

p_top_v_selection <- vfamily_cdr_rs %>%
  filter(V_family %in% top_selected_vfamilies$V_family) %>%
  ggplot(
    aes(
      x = reorder(V_family, CDR_RS_ratio),
      y = CDR_RS_ratio,
      fill = Group
    )
  ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(
    aes(color = Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "V Families with Highest CDR R/S Ratios",
    x = "V family",
    y = "CDR R/S ratio"
  )

p_top_v_selection



#############################################new55
#Corrected setup and definitions
############################################################
############################################################
## SHM selection analysis: public clonotypes and IGHV families
## Atlantic salmon SAV experiment ??? balanced 12-sample dataset
##
## Required objects already loaded in the R environment:
##   1. region_rs_raw    regional replacement/silent mutation table
##   2. all_clones_clean complete productive clonotype table
##
## The script explicitly supports sample names ending in:
##   .clones_IGH
############################################################

## 2. User settings --------------------------------------------

balanced_samples <- c(
  "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
  "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
  "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
  "SpI_1I3w", "SpI_2I3w", "SpI_5I3w"
)

expected_groups <- c("HKC", "HKI", "SpC", "SpI")
public_min_samples <- 2L
vfamily_min_samples <- 3L
top_n_vfamilies <- 15L

shm_dir <- file.path(getwd(), "SHM_public_selection_results")
dir.create(shm_dir, recursive = TRUE, showWarnings = FALSE)

## 3. General helper functions --------------------------------

abort_with_values <- function(message, values = character()) {
  if (length(values) > 0L) {
    message <- paste0(message, "\n  ", paste(values, collapse = "\n  "))
  }
  stop(message, call. = FALSE)
}

require_object <- function(object_name) {
  if (!exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
    stop("Required object `", object_name,
         "` is not present in the global R environment.", call. = FALSE)
  }
}

require_columns <- function(df, columns, table_name) {
  missing <- setdiff(columns, names(df))
  if (length(missing) > 0L) {
    stop(table_name, " is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

first_existing_column <- function(df, candidates) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) == 0L) NA_character_ else found[[1L]]
}

safe_numeric <- function(x, column_name, table_name) {
  old <- as.character(x)
  new <- suppressWarnings(as.numeric(old))
  bad <- !is.na(old) & trimws(old) != "" & is.na(new)
  if (any(bad)) {
    examples <- unique(old[bad])[seq_len(min(10L, length(unique(old[bad]))))]
    abort_with_values(
      paste0("Non-numeric values were found in `", column_name,
             "` of ", table_name, ":"), examples
    )
  }
  new
}

safe_ratio <- function(R, S) {
  ## A +1 pseudocount avoids infinite ratios when S = 0.
  R / (S + 1)
}

save_csv <- function(x, filename) {
  readr::write_csv(x, file.path(shm_dir, filename), na = "")
}

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(shm_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(shm_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600)
}

## 4. Sample-name normalization -------------------------------

standardize_sample_names <- function(df, table_name = deparse(substitute(df))) {
  require_columns(df, "Sample", table_name)
  
  original <- trimws(as.character(df$Sample))
  
  ## Primary rule: remove the established MiXCR suffix exactly.
  cleaned <- sub("\\.clones_IGH$", "", original, ignore.case = TRUE)
  
  ## Allow paths while retaining the basename.
  cleaned <- basename(gsub("\\\\", "/", cleaned))
  
  ## Defensive fallback for extra prefixes/suffixes. This extracts only a
  ## canonical SAV sample identifier; it does not invent a sample name.
  extracted <- stringr::str_extract(
    cleaned,
    stringr::regex("(?:HKC|HKI|SpC|SpI)_[0-9]+[CI]3w", ignore_case = TRUE)
  )
  
  canonical <- case_when(
    str_detect(extracted, regex("^HKC_", ignore_case = TRUE)) ~
      str_replace(extracted, regex("^HKC_", ignore_case = TRUE), "HKC_"),
    str_detect(extracted, regex("^HKI_", ignore_case = TRUE)) ~
      str_replace(extracted, regex("^HKI_", ignore_case = TRUE), "HKI_"),
    str_detect(extracted, regex("^SpC_", ignore_case = TRUE)) ~
      str_replace(extracted, regex("^SpC_", ignore_case = TRUE), "SpC_"),
    str_detect(extracted, regex("^SpI_", ignore_case = TRUE)) ~
      str_replace(extracted, regex("^SpI_", ignore_case = TRUE), "SpI_"),
    TRUE ~ NA_character_
  )
  canonical <- sub("W$", "w", canonical)
  
  df$Sample_original <- original
  df$Sample <- canonical
  
  unrecognized <- unique(original[is.na(canonical)])
  if (length(unrecognized) > 0L) {
    abort_with_values(
      paste0("Unrecognized sample names in ", table_name,
             ". Expected names such as `HKC_3C3w.clones_IGH`:"),
      unrecognized
    )
  }
  df
}

derive_sample_metadata <- function(df) {
  df %>%
    mutate(
      Group = str_extract(Sample, "^(HKC|HKI|SpC|SpI)"),
      Tissue = if_else(str_starts(Group, "HK"), "HK", "Spleen"),
      Condition = if_else(str_ends(Group, "I"), "Infected", "Control"),
      Group = factor(Group, levels = expected_groups),
      Tissue = factor(Tissue, levels = c("HK", "Spleen")),
      Condition = factor(Condition, levels = c("Control", "Infected"))
    )
}

validate_balanced_samples <- function(df, table_name) {
  observed <- df %>% distinct(Sample, Group, Tissue, Condition)
  duplicate_metadata <- observed %>% count(Sample) %>% filter(n != 1L)
  if (nrow(duplicate_metadata) > 0L) {
    abort_with_values(paste0("Conflicting metadata in ", table_name, ":"),
                      duplicate_metadata$Sample)
  }
  
  missing <- setdiff(balanced_samples, observed$Sample)
  unexpected <- setdiff(observed$Sample, balanced_samples)
  if (length(missing) > 0L || length(unexpected) > 0L) {
    abort_with_values(
      paste0(table_name, " failed the balanced-sample check. Missing: ",
             paste(missing, collapse = ", "), "; unexpected: ",
             paste(unexpected, collapse = ", "))
    )
  }
  
  counts <- observed %>% count(Group, name = "n_samples") %>%
    complete(Group = factor(expected_groups, levels = expected_groups),
             fill = list(n_samples = 0L))
  if (nrow(observed) != 12L || any(counts$n_samples != 3L)) {
    stop(table_name, " must contain exactly 12 samples and 3 per group.",
         call. = FALSE)
  }
  save_csv(observed %>% arrange(Group, Sample),
           paste0(table_name, "_sample_check.csv"))
  message(table_name, ": verified 12 balanced samples (3 per group).")
  invisible(observed)
}

## 5. Confirm inputs -------------------------------------------
region_rs_raw <- region_rs
require_object("region_rs_raw")
require_object("all_clones_clean")

require_columns(region_rs_raw,
                c("Sample", "CDR3.aa", "Region_class",
                  "R_weighted", "S_weighted", "Clones"),
                "region_rs_raw")
require_columns(all_clones_clean, c("Sample", "CDR3.aa"),
                "all_clones_clean")

## 6. Normalize and inspect all source sample names ------------

region_name_map <- region_rs_raw %>%
  standardize_sample_names("region_rs_raw") %>%
  distinct(Sample_original, Sample) %>% arrange(Sample)

clone_name_map <- all_clones_clean %>%
  standardize_sample_names("all_clones_clean") %>%
  distinct(Sample_original, Sample) %>% arrange(Sample)

save_csv(region_name_map, "region_rs_sample_name_mapping.csv")
save_csv(clone_name_map, "all_clones_sample_name_mapping.csv")

message("Recognized region_rs_raw samples: ",
        paste(region_name_map$Sample, collapse = ", "))

## 7. Prepare regional R/S data --------------------------------
region_rs <- region_rs_raw %>%
  standardize_sample_names("region_rs_raw") %>%
  filter(Sample %in% balanced_samples) %>%
  derive_sample_metadata() %>%
  mutate(
    CDR3.aa = trimws(as.character(CDR3.aa)),
    Region_class = toupper(trimws(as.character(Region_class))),
    R_weighted = safe_numeric(R_weighted, "R_weighted", "region_rs_raw"),
    S_weighted = safe_numeric(S_weighted, "S_weighted", "region_rs_raw"),
    Clones = safe_numeric(Clones, "Clones", "region_rs_raw")
  ) %>%
  filter(!is.na(CDR3.aa), CDR3.aa != "")

if (nrow(region_rs) == 0L) {
  stop("No regional R/S rows remained after balanced-sample filtering.",
       call. = FALSE)
}
if (any(region_rs$R_weighted < 0, na.rm = TRUE) ||
    any(region_rs$S_weighted < 0, na.rm = TRUE) ||
    any(region_rs$Clones < 0, na.rm = TRUE)) {
  stop("R_weighted, S_weighted, and Clones cannot contain negative values.",
       call. = FALSE)
}
validate_balanced_samples(region_rs, "region_rs")

region_levels <- sort(unique(region_rs$Region_class))
message("Region_class values: ", paste(region_levels, collapse = ", "))
if (!"CDR" %in% region_levels) {
  stop("`Region_class` contains no `CDR` rows. Available values: ",
       paste(region_levels, collapse = ", "), call. = FALSE)
}

## 8. Derive IGHV family safely -------------------------------

if (!"V_family" %in% names(region_rs)) {
  v_col <- first_existing_column(region_rs,
                                 c("V.name", "V_gene", "VGene", "bestVHit"))
  if (is.na(v_col)) {
    warning("No V-family or recognized V-gene column was found; IGHV-family ",
            "analysis will be skipped.", call. = FALSE)
  } else {
    region_rs <- region_rs %>%
      mutate(
        V_gene_source = trimws(as.character(.data[[v_col]])),
        V_family = V_gene_source %>%
          sub("\\*.*$", "", .) %>%
          sub("-.*$", "", .) %>%
          trimws()
      )
    message("V_family derived from `", v_col, "`.")
  }
} else {
  region_rs <- region_rs %>%
    mutate(V_family = trimws(as.character(V_family)))
}

## 9. Prepare complete clonotype presence table ----------------

all_clones_balanced <- all_clones_clean %>%
  standardize_sample_names("all_clones_clean") %>%
  filter(Sample %in% balanced_samples) %>%
  derive_sample_metadata() %>%
  mutate(CDR3.aa = trimws(as.character(CDR3.aa))) %>%
  filter(!is.na(CDR3.aa), CDR3.aa != "")

if (nrow(all_clones_balanced) == 0L) {
  stop("No clonotype rows remained after balanced-sample filtering.",
       call. = FALSE)
}
validate_balanced_samples(all_clones_balanced, "all_clones_balanced")

## Use the union of both inputs for presence/absence classification.
## This prevents a clonotype retained in the R/S table from being called absent.
clonotype_presence <- bind_rows(
  all_clones_balanced %>%
    transmute(Sample, Group = as.character(Group),
              Tissue = as.character(Tissue),
              Condition = as.character(Condition), CDR3.aa),
  region_rs %>%
    transmute(Sample, Group = as.character(Group),
              Tissue = as.character(Tissue),
              Condition = as.character(Condition), CDR3.aa)
) %>%
  distinct(Sample, Group, Tissue, Condition, CDR3.aa) %>%
  mutate(
    Group = factor(Group, levels = expected_groups),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

metadata_conflicts <- clonotype_presence %>%
  group_by(Sample) %>%
  summarise(n_groups = n_distinct(Group), n_tissues = n_distinct(Tissue),
            n_conditions = n_distinct(Condition), .groups = "drop") %>%
  filter(n_groups != 1L | n_tissues != 1L | n_conditions != 1L)
if (nrow(metadata_conflicts) > 0L) {
  stop("Conflicting sample metadata were detected between input tables.",
       call. = FALSE)
}

## Record input-table coverage differences for auditability.
rs_pairs <- region_rs %>% distinct(Sample, CDR3.aa)
clone_pairs <- all_clones_balanced %>% distinct(Sample, CDR3.aa)
save_csv(anti_join(rs_pairs, clone_pairs, by = c("Sample", "CDR3.aa")),
         "sample_clonotypes_in_region_RS_only.csv")
save_csv(anti_join(clone_pairs, rs_pairs, by = c("Sample", "CDR3.aa")),
         "sample_clonotypes_in_all_clones_only.csv")

## 10. Define public and condition-specific public clonotypes --

public_by_group <- clonotype_presence %>%
  group_by(Group, CDR3.aa) %>%
  summarise(n_samples_present = n_distinct(Sample), .groups = "drop") %>%
  filter(n_samples_present >= public_min_samples)

public_sets <- setNames(
  lapply(expected_groups, function(g) {
    public_by_group %>% filter(Group == g) %>% pull(CDR3.aa) %>% unique()
  }), expected_groups
)

infected_public_candidates <- union(public_sets$HKI, public_sets$SpI)
control_public_candidates <- union(public_sets$HKC, public_sets$SpC)

all_control_clonotypes <- clonotype_presence %>%
  filter(Condition == "Control") %>% pull(CDR3.aa) %>% unique()
all_infected_clonotypes <- clonotype_presence %>%
  filter(Condition == "Infected") %>% pull(CDR3.aa) %>% unique()

## User-approved definitions:
## infection-specific public = public in HKI and/or SpI AND absent from every
## control sample; control-specific is defined symmetrically.
infected_specific_public <- setdiff(infected_public_candidates,
                                    all_control_clonotypes)
control_specific_public <- setdiff(control_public_candidates,
                                   all_infected_clonotypes)
shared_public <- intersect(infected_public_candidates,
                           control_public_candidates)

## Public in one condition but present non-publicly in the other condition.
## These belong to Other/private in the four-category plot, not Shared public.
cross_condition_nonshared <- setdiff(
  union(intersect(infected_public_candidates, all_control_clonotypes),
        intersect(control_public_candidates, all_infected_clonotypes)),
  shared_public
)

if (length(intersect(infected_specific_public, all_control_clonotypes)) > 0L ||
    length(intersect(control_specific_public, all_infected_clonotypes)) > 0L ||
    length(intersect(infected_specific_public, control_specific_public)) > 0L) {
  stop("Condition-specific public-clonotype validation failed.", call. = FALSE)
}

definition_summary <- tibble(
  Category = c("HKC public", "HKI public", "SpC public", "SpI public",
               "Infected public candidates", "Control public candidates",
               "Infection-specific public", "Control-specific public",
               "Strict shared public", "Cross-condition nonshared"),
  n_clonotypes = c(length(public_sets$HKC), length(public_sets$HKI),
                   length(public_sets$SpC), length(public_sets$SpI),
                   length(infected_public_candidates),
                   length(control_public_candidates),
                   length(infected_specific_public),
                   length(control_specific_public), length(shared_public),
                   length(cross_condition_nonshared))
)
print(definition_summary, n = Inf)
save_csv(definition_summary, "public_clonotype_definition_summary.csv")
save_csv(tibble(CDR3.aa = infected_specific_public),
         "infection_specific_public_clonotypes.csv")
save_csv(tibble(CDR3.aa = control_specific_public),
         "control_specific_public_clonotypes.csv")
save_csv(tibble(CDR3.aa = shared_public),
         "strict_shared_public_clonotypes.csv")
save_csv(tibble(CDR3.aa = cross_condition_nonshared),
         "cross_condition_nonshared_clonotypes.csv")

## 11. Assign categories to regional R/S rows ------------------

category_levels <- c("Private", "Control-specific public",
                     "Shared public", "Infection-specific public")

public_rs <- region_rs %>%
  mutate(
    Public_category = case_when(
      CDR3.aa %in% infected_specific_public ~ "Infection-specific public",
      CDR3.aa %in% control_specific_public ~ "Control-specific public",
      CDR3.aa %in% shared_public ~ "Shared public",
      TRUE ~ "Private"
    ),
    Public_category = factor(Public_category, levels = category_levels)
  )

assignment_conflicts <- public_rs %>%
  filter((Public_category == "Control-specific public" &
            Condition != "Control") |
           (Public_category == "Infection-specific public" &
              Condition != "Infected")) %>%
  distinct(Sample, Group, Tissue, Condition, CDR3.aa, Public_category)
save_csv(assignment_conflicts, "public_category_assignment_conflicts.csv")
if (nrow(assignment_conflicts) > 0L) {
  stop("Invalid condition-specific assignments remain; see conflict CSV.",
       call. = FALSE)
}

## 12. Sample-level CDR R/S summaries -------------------------

public_cdr_rs <- public_rs %>%
  filter(Region_class == "CDR") %>%
  group_by(Sample, Public_category) %>%
  summarise(
    # Explicit namespace prevents BiocGenerics::first() from masking
    # dplyr::first() when Bioconductor packages are attached.
    Group = dplyr::first(as.character(Group)),
    Tissue = dplyr::first(as.character(Tissue)),
    Condition = dplyr::first(as.character(Condition)),
    total_R = sum(R_weighted, na.rm = TRUE),
    total_S = sum(S_weighted, na.rm = TRUE),
    CDR_RS_ratio = safe_ratio(total_R, total_S),
    R_fraction = total_R / (total_R + total_S + 1),
    n_CDR3 = n_distinct(CDR3.aa),
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = expected_groups),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected")),
    Public_category = factor(Public_category, levels = category_levels)
  )

duplicates <- public_cdr_rs %>% count(Sample, Public_category) %>% filter(n > 1L)
if (nrow(duplicates) > 0L) {
  stop("More than one row per Sample ?? Public_category was generated.",
       call. = FALSE)
}

category_condition_check <- public_cdr_rs %>%
  distinct(Public_category, Condition, Sample) %>%
  count(Public_category, Condition, name = "n_samples") %>%
  arrange(Public_category, Condition)
print(category_condition_check, n = Inf)
save_csv(category_condition_check, "public_category_condition_check.csv")

invalid <- category_condition_check %>%
  filter((Public_category == "Control-specific public" &
            Condition != "Control") |
           (Public_category == "Infection-specific public" &
              Condition != "Infected"))
if (nrow(invalid) > 0L) stop("Category-condition validation failed.", call. = FALSE)

sample_category_map <- public_cdr_rs %>%
  distinct(Public_category, Condition, Tissue, Group, Sample) %>%
  arrange(Public_category, Group, Sample)
save_csv(sample_category_map, "public_category_sample_check.csv")

public_cdr_rs_summary <- public_cdr_rs %>%
  group_by(Public_category) %>%
  summarise(
    n_samples = n_distinct(Sample),
    mean_RS = mean(CDR_RS_ratio, na.rm = TRUE),
    median_RS = median(CDR_RS_ratio, na.rm = TRUE),
    sd_RS = if_else(n() > 1L, sd(CDR_RS_ratio, na.rm = TRUE), NA_real_),
    min_RS = min(CDR_RS_ratio, na.rm = TRUE),
    max_RS = max(CDR_RS_ratio, na.rm = TRUE),
    mean_R_fraction = mean(R_fraction, na.rm = TRUE),
    total_CDR3 = sum(n_CDR3),
    total_clone_abundance = sum(total_clones),
    .groups = "drop"
  )
print(public_cdr_rs_summary, n = Inf)
save_csv(public_cdr_rs, "public_clonotype_sample_CDR_RS.csv")
save_csv(public_cdr_rs_summary, "public_clonotype_CDR_RS_summary.csv")


## 13. Public-category plot ------------------------------------

public_colors <- c(
  "Private" = "grey75",
  "Control-specific public" = "#3B82C4",
  "Shared public" = "#8C6BB1",
  "Infection-specific public" = "#E76F51"
)

p_public_rs <- ggplot(public_cdr_rs,
                      aes(Public_category, CDR_RS_ratio,
                          fill = Public_category, color = Public_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
  geom_jitter(width = 0.12, height = 0, size = 2.6, alpha = 0.85) +
  scale_fill_manual(values = public_colors, drop = FALSE) +
  scale_color_manual(values = public_colors, drop = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none", panel.grid.minor = element_blank()) +
  labs(title = "CDR1-CDR2 R/S Ratio by Public Clonotype Category",
        x = NULL, y = "CDR replacement-to-silent mutation ratio")

print(p_public_rs)
save_plot(p_public_rs, "public_clonotype_CDR_RS", 8.5, 6.5)

## 14. Exploratory targeted statistical comparison ------------
##
## An omnibus Kruskal-Wallis test across all four categories is not used:
## categories recur within samples, while condition-specific categories are
## structurally restricted to one condition. Treating all dots as independent
## would violate the test design. The only direct exploratory contrast below
## compares control-specific with infection-specific sample summaries.

specific_test_data <- public_cdr_rs %>%
  filter(Public_category %in% c("Control-specific public",
                                "Infection-specific public")) %>%
  droplevels()

specific_counts <- specific_test_data %>%
  count(Public_category, name = "n_samples")
save_csv(specific_counts, "specific_public_statistical_sample_sizes.csv")

specific_test_result <- tibble(
  comparison = "Control-specific vs infection-specific public",
  test = "Wilcoxon rank-sum (exploratory)",
  statistic = NA_real_, p_value = NA_real_, note = NA_character_
)

if (n_distinct(specific_test_data$Public_category) == 2L &&
    all(specific_counts$n_samples >= 2L)) {
  wt <- wilcox.test(CDR_RS_ratio ~ Public_category,
                    data = specific_test_data, exact = FALSE)
  specific_test_result$statistic <- unname(wt$statistic)
  specific_test_result$p_value <- wt$p.value
  specific_test_result$note <- paste0(
    "Exploratory: category is structurally linked to condition; small n."
  )
} else {
  specific_test_result$note <-
    "Not run: both categories require at least two sample-level observations."
}
print(specific_test_result, n = Inf)
save_csv(specific_test_result, "specific_public_exploratory_Wilcoxon.csv")

## 15. IGHV-family CDR R/S analysis ---------------------------

if ("V_family" %in% names(region_rs)) {
  vfamily_cdr_rs <- region_rs %>%
    filter(Region_class == "CDR", !is.na(V_family), V_family != "") %>%
    group_by(Sample, V_family) %>%
    summarise(
      # Explicit namespace prevents BiocGenerics::first() from masking
      # dplyr::first() when Bioconductor packages are attached.
      Group = dplyr::first(as.character(Group)),
      Tissue = dplyr::first(as.character(Tissue)),
      Condition = dplyr::first(as.character(Condition)),
      total_R = sum(R_weighted, na.rm = TRUE),
      total_S = sum(S_weighted, na.rm = TRUE),
      CDR_RS_ratio = safe_ratio(total_R, total_S),
      R_fraction = total_R / (total_R + total_S + 1),
      total_clones = sum(Clones, na.rm = TRUE),
      n_CDR3 = n_distinct(CDR3.aa),
      .groups = "drop"
    ) %>%
    mutate(
      Group = factor(Group, levels = expected_groups),
      Tissue = factor(Tissue, levels = c("HK", "Spleen")),
      Condition = factor(Condition, levels = c("Control", "Infected"))
    )
  
  vfamily_summary <- vfamily_cdr_rs %>%
    group_by(V_family) %>%
    summarise(
      mean_CDR_RS = mean(CDR_RS_ratio, na.rm = TRUE),
      median_CDR_RS = median(CDR_RS_ratio, na.rm = TRUE),
      sd_CDR_RS = if_else(n() > 1L, sd(CDR_RS_ratio, na.rm = TRUE), NA_real_),
      min_CDR_RS = min(CDR_RS_ratio, na.rm = TRUE),
      max_CDR_RS = max(CDR_RS_ratio, na.rm = TRUE),
      mean_R_fraction = mean(R_fraction, na.rm = TRUE),
      total_clones = sum(total_clones), total_CDR3 = sum(n_CDR3),
      n_samples = n_distinct(Sample), n_groups = n_distinct(Group),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_CDR_RS), desc(n_samples))
  
  eligible_vfamilies <- vfamily_summary %>%
    filter(n_samples >= vfamily_min_samples)
  top_selected_vfamilies <- eligible_vfamilies %>%
    slice_max(mean_CDR_RS, n = top_n_vfamilies, with_ties = FALSE)
  
  save_csv(vfamily_cdr_rs, "Vfamily_sample_CDR_RS.csv")
  save_csv(vfamily_summary, "Vfamily_CDR_RS_summary.csv")
  save_csv(eligible_vfamilies, "eligible_Vfamilies.csv")
  save_csv(top_selected_vfamilies, "top_15_Vfamilies_CDR_RS.csv")
  
  if (nrow(top_selected_vfamilies) > 0L) {
    plot_order <- top_selected_vfamilies %>%
      arrange(mean_CDR_RS) %>% pull(V_family)
    vfamily_plot_data <- vfamily_cdr_rs %>%
      filter(V_family %in% top_selected_vfamilies$V_family) %>%
      mutate(V_family = factor(V_family, levels = plot_order))
    
    group_colors <- c(HKC = "#4C78A8", HKI = "#E45756",
                      SpC = "#72B7B2", SpI = "#F2A541")
    
    ## Points are shown by group. Boxplots are intentionally not dodged because
    ## many V-family ?? group cells may contain fewer than three samples.
    p_top_v <- ggplot(vfamily_plot_data,
                      aes(V_family, CDR_RS_ratio, color = Group)) +
      geom_point(position = position_jitter(width = 0.14, height = 0),
                 size = 2.2, alpha = 0.85) +
      stat_summary(aes(group = V_family), fun = median, geom = "crossbar",
                   width = 0.55, color = "black", linewidth = 0.35) +
      coord_flip() +
      scale_color_manual(values = group_colors, drop = FALSE) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(), legend.position = "right") +
      labs(title = "Top 15 IGHV Families with the Highest CDR1-CDR2 R/S Ratios",
           x = "IGHV family",
           y = "CDR replacement-to-silent mutation ratio", color = "Group")
    
    print(p_top_v)
    save_plot(p_top_v, "top_Vfamilies_CDR_RS", 9, 7.5)
  } else {
    warning("No IGHV family met the minimum sample criterion.", call. = FALSE)
  }
}










#############################################End 55

############################################################
## 5. Does spleen accumulate more CDR replacement mutations?
############################################################

cdr_tissue_rs <- region_rs %>%
  filter(Region_class == "CDR") %>%
  group_by(Sample, Group, Tissue, Condition) %>%
  summarise(
    total_R = sum(R_weighted, na.rm = TRUE),
    total_S = sum(S_weighted, na.rm = TRUE),
    CDR_RS_ratio = total_R / (total_S + 1),
    R_fraction = total_R / (total_R + total_S + 1),
    .groups = "drop"
  )

cdr_tissue_summary <- cdr_tissue_rs %>%
  group_by(Tissue, Condition) %>%
  summarise(
    mean_RS = mean(CDR_RS_ratio, na.rm = TRUE),
    median_RS = median(CDR_RS_ratio, na.rm = TRUE),
    sd_RS = sd(CDR_RS_ratio, na.rm = TRUE),
    mean_R_fraction = mean(R_fraction, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(cdr_tissue_summary)

p_tissue_cdr <- ggplot(
  cdr_tissue_rs,
  aes(
    x = Tissue,
    y = CDR_RS_ratio,
    fill = Tissue
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  facet_wrap(~ Condition) +
  theme_bw(base_size = 13) +
  labs(
    title = "Tissue Effect on CDR R/S Ratio",
    x = "",
    y = "CDR R/S ratio"
  )

p_tissue_cdr

tissue_test <- wilcox.test(
  CDR_RS_ratio ~ Tissue,
  data = cdr_tissue_rs,
  exact = FALSE
)

print(tissue_test)

tissue_condition_model <- lm(
  CDR_RS_ratio ~ Tissue * Condition,
  data = cdr_tissue_rs
)

print(anova(tissue_condition_model))
print(summary(tissue_condition_model))

############################################################
## 6. Which IGHV families show positive selection only after infection?
############################################################

vfamily_infection_tests <- vfamily_cdr_rs %>%
  group_by(V_family, Tissue) %>%
  summarise(
    n_control = sum(Condition == "Control"),
    n_infected = sum(Condition == "Infected"),
    mean_control = mean(CDR_RS_ratio[Condition == "Control"], na.rm = TRUE),
    mean_infected = mean(CDR_RS_ratio[Condition == "Infected"], na.rm = TRUE),
    log2FC_infected_vs_control =
      log2((mean_infected + 1e-6) / (mean_control + 1e-6)),
    pvalue = tryCatch(
      wilcox.test(
        CDR_RS_ratio ~ Condition,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    infection_selected =
      mean_infected > mean_control & log2FC_infected_vs_control > 0.25
  ) %>%
  arrange(padj)

print(vfamily_infection_tests)

p_vfamily_infection <- ggplot(
  vfamily_infection_tests,
  aes(
    x = log2FC_infected_vs_control,
    y = -log10(padj),
    color = infection_selected
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Tissue) +
  theme_bw(base_size = 13) +
  labs(
    title = "Infection-associated Increase in CDR R/S by V Family",
    x = "log2 fold-change infected/control",
    y = "-log10 adjusted P-value",
    color = "Higher in infected"
  )

p_vfamily_infection

############################################################
## 7. Candidate antiviral IGHV families:
##    combine V-family enrichment + high CDR R/S
############################################################
## This section assumes you have all_clones_full or V usage table.
## If not, it will build from region_rs as approximation.

vfamily_abundance <- region_rs %>%
  group_by(Sample, Group, Tissue, Condition, V_family) %>%
  summarise(
    total_clones = sum(Clones, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Sample) %>%
  mutate(
    V_fraction = total_clones / sum(total_clones, na.rm = TRUE)
  ) %>%
  ungroup()

vfamily_enrichment <- vfamily_abundance %>%
  group_by(V_family, Tissue) %>%
  summarise(
    mean_control = mean(V_fraction[Condition == "Control"], na.rm = TRUE),
    mean_infected = mean(V_fraction[Condition == "Infected"], na.rm = TRUE),
    log2FC_abundance =
      log2((mean_infected + 1e-6) / (mean_control + 1e-6)),
    pvalue_abundance = tryCatch(
      wilcox.test(
        V_fraction ~ Condition,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj_abundance = p.adjust(pvalue_abundance, method = "BH")
  )

candidate_vfamilies <- vfamily_enrichment %>%
  left_join(
    vfamily_infection_tests %>%
      select(
        V_family,
        Tissue,
        log2FC_RS = log2FC_infected_vs_control,
        padj_RS = padj,
        mean_control_RS = mean_control,
        mean_infected_RS = mean_infected
      ),
    by = c("V_family", "Tissue")
  ) %>%
  mutate(
    candidate_score =
      scale(log2FC_abundance)[, 1] +
      scale(log2FC_RS)[, 1],
    Candidate = log2FC_abundance > 0 &
      log2FC_RS > 0 &
      mean_infected_RS > median(mean_infected_RS, na.rm = TRUE)
  ) %>%
  arrange(desc(candidate_score))

print(candidate_vfamilies)

top_candidate_vfamilies <- candidate_vfamilies %>%
  filter(Candidate) %>%
  slice_max(candidate_score, n = 15)

print(top_candidate_vfamilies)

p_candidate <- ggplot(
  candidate_vfamilies,
  aes(
    x = log2FC_abundance,
    y = log2FC_RS,
    color = Candidate,
    label = V_family
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ Tissue) +
  theme_bw(base_size = 13) +
  labs(
    title = "Candidate Antiviral V Families",
    subtitle = "Combining infection enrichment and CDR R/S increase",
    x = "V-family abundance log2FC infected/control",
    y = "CDR R/S log2FC infected/control"
  )

p_candidate

############################################################
## 7b. Candidate plot with V-family labels
############################################################

library(ggrepel)

p_candidate_labeled <- ggplot(
  candidate_vfamilies,
  aes(
    x = log2FC_abundance,
    y = log2FC_RS
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    aes(color = Candidate),
    size = 3,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = candidate_vfamilies %>%
      filter(Candidate == TRUE),
    aes(label = V_family),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  facet_wrap(~ Tissue) +
  theme_bw(base_size = 13) +
  labs(
    title = "Candidate Antiviral V Families",
    subtitle = "V families enriched after infection and showing increased CDR R/S ratio",
    x = "V-family abundance log2FC infected/control",
    y = "CDR R/S log2FC infected/control",
    color = "Candidate"
  )

p_candidate_labeled



############################################################
## 8. Heatmap of candidate V-family selection
############################################################

candidate_heatmap_df <- vfamily_cdr_rs %>%
  group_by(Group, V_family) %>%
  summarise(
    mean_RS = mean(CDR_RS_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  semi_join(
    candidate_vfamilies %>%
      slice_max(candidate_score, n = 20),
    by = "V_family"
  )

candidate_mat <- candidate_heatmap_df %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_RS,
    values_fill = 0
  ) %>%
  column_to_rownames("V_family")

pheatmap(
  candidate_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "CDR R/S Ratio of Candidate Antiviral V Families"
)




######################
#Don't run this otherwise you will risk not seeing the figures in the output pane
############3




################################################################################
##New analysis 
############################################################
## Link clonal expansion with SHM and CDR selection
## Question:
## Do expanded clones accumulate more SHM and show higher CDR R/S?
############################################################

library(dplyr)
library(ggplot2)
library(readr)

############################################################
## 1. Set working directory and load data
############################################################

workdir <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/SHM_mutation_spectrum_HK_Sp_3w"
setwd(workdir)

region_rs <- readRDS("region_rs_SAV_3wpi.rds")

############################################################
## 2. Prepare clone-level summary
############################################################

clone_shm <- region_rs %>%
  mutate(
    Tissue = recode(as.character(Tissue),
                    "Sp" = "Spleen",
                    "SP" = "Spleen"),
    Region_class = case_when(
      Region %in% c("CDR1", "CDR2") ~ "CDR",
      Region %in% c("FR1", "FR2", "FR3") ~ "FR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Region_class)) %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    CDR3.aa,
    V_family,
    J_family,
    Clones
  ) %>%
  summarise(
    total_R = sum(R, na.rm = TRUE),
    total_S = sum(S, na.rm = TRUE),
    total_mutations = total_R + total_S,
    
    CDR_R = sum(R[Region_class == "CDR"], na.rm = TRUE),
    CDR_S = sum(S[Region_class == "CDR"], na.rm = TRUE),
    FR_R = sum(R[Region_class == "FR"], na.rm = TRUE),
    FR_S = sum(S[Region_class == "FR"], na.rm = TRUE),
    
    CDR_RS_ratio = CDR_R / (CDR_S + 1),
    FR_RS_ratio = FR_R / (FR_S + 1),
    
    .groups = "drop"
  ) %>%
  mutate(
    log10_clone_size = log10(Clones + 1),
    log10_mutations = log10(total_mutations + 1),
    expanded_clone = Clones >= quantile(Clones, 0.90, na.rm = TRUE),
    highly_mutated = total_mutations >= quantile(total_mutations, 0.75, na.rm = TRUE),
    high_CDR_RS = CDR_RS_ratio >= quantile(CDR_RS_ratio, 0.75, na.rm = TRUE),
    expanded_mutated_selected =
      expanded_clone & highly_mutated & high_CDR_RS
  )

write_csv(clone_shm, "clone_expansion_SHM_CDR_selection_summary.csv")

############################################################
## 3. Correlation: clone size vs SHM burden
############################################################

cor_results <- clone_shm %>%
  group_by(Tissue, Condition) %>%
  summarise(
    n_clones = n(),
    spearman_rho = cor(
      log10_clone_size,
      log10_mutations,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = cor.test(
      log10_clone_size,
      log10_mutations,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH")
  )

write_csv(cor_results, "correlation_clone_size_vs_SHM.csv")
print(cor_results)

############################################################
## 4. Figure A: clone expansion vs SHM burden
############################################################

p1 <- ggplot(
  clone_shm,
  aes(
    x = log10_clone_size,
    y = log10_mutations,
    color = Condition
  )
) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~Tissue) +
  theme_bw(base_size = 14) +
  labs(
    title = "Expanded clones accumulate more SHM",
    x = "log10(Clone abundance + 1)",
    y = "log10(V-region mutations + 1)",
    color = "Condition"
  )

print(p1)

ggsave(
  "Figure_clone_expansion_vs_SHM.png",
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

################STAT Start
library(dplyr)
library(lmerTest)
library(emmeans)

# Keep complete, finite observations and set factor order
clone_shm_test <- clone_shm %>%
  filter(
    is.finite(log10_clone_size),
    is.finite(log10_mutations),
    !is.na(Sample),
    !is.na(Tissue),
    !is.na(Condition)
  ) %>%
  mutate(
    Sample = factor(Sample),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )

# Test whether the clone-size/mutation-burden slope varies by tissue and condition.
# The sample random intercept accounts for clonotypes from the same sample.
shm_model <- lmer(
  log10_mutations ~ log10_clone_size * Tissue * Condition +
    (1 | Sample),
  data = clone_shm_test
)

# Slope within each tissue-condition group; tests whether each slope differs from zero
shm_slopes <- emtrends(
  shm_model,
  ~ Condition | Tissue,
  var = "log10_clone_size"
)

slope_results <- as.data.frame(
  summary(shm_slopes, infer = c(TRUE, TRUE), adjust = "none")
) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

print(slope_results)

# Compare infected and control slopes separately within HK and spleen
slope_comparisons <- as.data.frame(
  pairs(shm_slopes, by = "Tissue", adjust = "none")
) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

print(slope_comparisons)
###STAT End

############################################################
## 5. Bin clones by abundance rank
############################################################

clone_shm_binned <- clone_shm %>%
  group_by(Sample) %>%
  mutate(
    abundance_rank = rank(-Clones, ties.method = "first"),
    abundance_bin = case_when(
      abundance_rank <= 10 ~ "Top 1-10",
      abundance_rank <= 100 ~ "Top 11-100",
      abundance_rank <= 1000 ~ "Top 101-1000",
      TRUE ~ "Rare >1000"
    )
  ) %>%
  ungroup() %>%
  mutate(
    abundance_bin = factor(
      abundance_bin,
      levels = c("Top 1-10", "Top 11-100", "Top 101-1000", "Rare >1000")
    )
  )

############################################################
## 6. Figure B: SHM burden across clone-abundance classes
############################################################

p2 <- ggplot(
  clone_shm_binned,
  aes(
    x = abundance_bin,
    y = total_mutations,
    fill = abundance_bin
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 14) +
  labs(
    title = "SHM burden increases in expanded clone classes",
    x = "Clone abundance class",
    y = "Total V-region mutations",
    fill = "Clone class"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

print(p2)

ggsave(
  "Figure_SHM_by_clone_abundance_class.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)

###########Stat start for abundance clases
library(dplyr)
library(glmmTMB)
library(emmeans)

# Confirm the source data are present
nrow(clone_shm)
names(clone_shm)

# Rebuild the binned data from clone_shm
clone_shm_binned <- clone_shm %>%
  filter(
    !is.na(Clones),
    Clones >= 1,
    !is.na(total_mutations),
    !is.na(Sample),
    !is.na(Tissue),
    !is.na(Condition)
  ) %>%
  mutate(
    abundance_bin = case_when(
      Clones <= 10   ~ "1-10",
      Clones <= 100  ~ "11-100",
      Clones <= 1000 ~ "101-1000",
      Clones > 1000  ~ ">1000"
    ),
    abundance_bin = factor(
      abundance_bin,
      levels = c("1-10", "11-100", "101-1000", ">1000"),
      ordered = TRUE
    ),
    Sample = factor(Sample),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  ) %>%
  droplevels()

# Check that rows and abundance classes are present
nrow(clone_shm_binned)
table(clone_shm_binned$abundance_bin, useNA = "ifany")
with(clone_shm_binned, table(Tissue, Condition, abundance_bin))
clone_shm_binned <- clone_shm_binned %>%
  mutate(
    abundance_bin = factor(
      abundance_bin,
      levels = c("1-10", "11-100")
    )
  ) %>%
  droplevels()

m_full <- glmmTMB(
  total_mutations ~ abundance_bin * Tissue * Condition + (1 | Sample),
  family = nbinom2,
  data = clone_shm_binned
)

m_reduced <- glmmTMB(
  total_mutations ~ Tissue * Condition + (1 | Sample),
  family = nbinom2,
  data = clone_shm_binned
)

# Overall test of whether abundance class contributes to the model
anova(m_reduced, m_full)

# Compare the two classes within each tissue-condition group
bin_emmeans <- emmeans(
  m_full,
  ~ abundance_bin | Tissue * Condition
)

bin_comparisons <- as.data.frame(
  pairs(bin_emmeans, by = c("Tissue", "Condition"), adjust = "none")
) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

print(bin_comparisons)

#######################################################
################
library(dplyr)
library(glmmTMB)
library(emmeans)

# Set the four abundance-rank classes in their intended order
rank_levels <- c("1-10", "11-100", "101-1000", ">1000")

# Rank clonotypes within each sample, from highest to lowest Clones count
clone_shm_binned <- clone_shm %>%
  dplyr::filter(
    !is.na(Sample),
    !is.na(Clones),
    !is.na(total_mutations),
    !is.na(Tissue),
    !is.na(Condition)
  ) %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(
    abundance_rank = rank(-Clones, ties.method = "min"),
    abundance_bin = dplyr::case_when(
      abundance_rank <= 10   ~ "1-10",
      abundance_rank <= 100  ~ "11-100",
      abundance_rank <= 1000 ~ "101-1000",
      TRUE                   ~ ">1000"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    abundance_bin = factor(abundance_bin, levels = rank_levels),
    Sample = factor(Sample),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  ) %>%
  droplevels()

# Check counts in all four classes and every tissue-condition stratum
class_counts <- with(
  clone_shm_binned,
  table(Tissue, Condition, abundance_bin)
)
print(class_counts)

if (!all(rank_levels %in% unique(as.character(clone_shm_binned$abundance_bin)))) {
  stop("At least one abundance-rank class is absent from the dataset.")
}

if (any(class_counts == 0)) {
  stop("At least one abundance-rank class is absent from a tissue-condition group.")
}

# Full model: abundance class, tissue, condition, and their interactions
m_full <- glmmTMB(
  total_mutations ~ abundance_bin * Tissue * Condition + (1 | Sample),
  family = nbinom2,
  data = clone_shm_binned
)

# Reduced model without abundance-class terms
m_reduced <- glmmTMB(
  total_mutations ~ Tissue * Condition + (1 | Sample),
  family = nbinom2,
  data = clone_shm_binned
)

# Overall likelihood-ratio test of abundance class and its interactions
overall_class_test <- anova(m_reduced, m_full)
print(overall_class_test)

# Compare all four classes within each tissue-condition group:
# six pairwise comparisons in each of the four strata
class_means <- emmeans(
  m_full,
  ~ abundance_bin | Tissue * Condition
)

class_comparisons <- as.data.frame(
  pairs(
    class_means,
    by = c("Tissue", "Condition"),
    adjust = "none"
  )
) %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH")
  )

print(class_comparisons, row.names = FALSE)
  


##########Stat end

############################################################
## 7. Figure C: CDR R/S ratio across clone-abundance classes
############################################################

p3 <- ggplot(
  clone_shm_binned,
  aes(
    x = abundance_bin,
    y = CDR_RS_ratio,
    fill = abundance_bin
  )
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 14) +
  labs(
    title = "Expanded clones show evidence of CDR-focused selection",
    x = "Clone abundance class",
    y = "CDR R/S ratio",
    fill = "Clone class"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

print(p3)

ggsave(
  "Figure_CDR_RS_by_clone_abundance_class.png",
  p3,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
## 8. Compare expanded vs non-expanded clones
############################################################

expanded_tests <- clone_shm %>%
  group_by(Tissue, Condition) %>%
  summarise(
    n_expanded = sum(expanded_clone, na.rm = TRUE),
    n_nonexpanded = sum(!expanded_clone, na.rm = TRUE),
    
    p_SHM = wilcox.test(
      total_mutations ~ expanded_clone,
      exact = FALSE
    )$p.value,
    
    p_CDR_RS = wilcox.test(
      CDR_RS_ratio ~ expanded_clone,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(
    padj_SHM = p.adjust(p_SHM, method = "BH"),
    padj_CDR_RS = p.adjust(p_CDR_RS, method = "BH")
  )

write_csv(expanded_tests, "expanded_vs_nonexpanded_SHM_CDR_RS_tests.csv")
print(expanded_tests)

############################################################
## 9. Identify strongest candidate matured antiviral clones
############################################################

candidate_clones <- clone_shm %>%
  filter(expanded_mutated_selected) %>%
  arrange(desc(Clones), desc(total_mutations), desc(CDR_RS_ratio))

write_csv(
  candidate_clones,
  "candidate_expanded_mutated_CDR_selected_clones.csv"
)

print(candidate_clones)

############################################################
## 10. Figure D: integrated clone-level view
############################################################

p4 <- ggplot(
  clone_shm,
  aes(
    x = log10_clone_size,
    y = CDR_RS_ratio
  )
) +
  geom_point(
    aes(
      color = expanded_mutated_selected,
      size = total_mutations
    ),
    alpha = 0.55
  ) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 14) +
  labs(
    title = "Expanded clonotypes with high mutation burden and CDR R/S ratios",
        x = "log10(Clone abundance + 1)",
    y = "CDR R/S ratio",
    color = "Expanded + mutated",
    size = "Total mutations"
  )

print(p4)

ggsave(
  "Figure_integrated_expansion_SHM_selection.png",
  p4,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
## 11. Summary table for manuscript
############################################################

summary_table <- clone_shm %>%
  group_by(Group, Tissue, Condition) %>%
  summarise(
    n_clones = n(),
    median_clone_size = median(Clones, na.rm = TRUE),
    median_total_mutations = median(total_mutations, na.rm = TRUE),
    median_CDR_RS = median(CDR_RS_ratio, na.rm = TRUE),
    n_expanded_mutated_selected =
      sum(expanded_mutated_selected, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_table, "summary_expansion_SHM_selection_by_group.csv")
print(summary_table)

############################################################
## DONE
############################################################

cat("\nAnalysis complete.\n")
cat("Main outputs saved:\n")
cat("- clone_expansion_SHM_CDR_selection_summary.csv\n")
cat("- correlation_clone_size_vs_SHM.csv\n")
cat("- expanded_vs_nonexpanded_SHM_CDR_RS_tests.csv\n")
cat("- candidate_expanded_mutated_CDR_selected_clones.csv\n")
cat("- Figure_clone_expansion_vs_SHM.png\n")
cat("- Figure_SHM_by_clone_abundance_class.png\n")
cat("- Figure_CDR_RS_by_clone_abundance_class.png\n")
cat("- Figure_integrated_expansion_SHM_selection.png\n")









###DO NOT RUN SCRIPTS BELOW THIS<<THEY JUST DILUTE THE ABOVE FINDINGS
##DON'T RUN THE BELOW



############################################################
## Clone-frequency classes linked to SHM and CDR selection
##
## Classes:
## Expanded:     >1% of the sample repertoire
## Intermediate: 0.1???1%
## Small:        0.01???0.1%
## Rare:         <0.01%
############################################################
#install.packages("lmerTest")
#install.packages("emmeans")
library(dplyr)
library(ggplot2)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)

############################################################
## 1. Working directory and input
############################################################

workdir <- paste0(
  "C:/Users/sje054/OneDrive - UiT Office 365/",
  "BCR/SHM_mutation_spectrum_HK_Sp_3w"
)

setwd(workdir)

region_rs <- readRDS("region_rs_SAV_3wpi.rds")

############################################################
## 2. Check essential columns
############################################################

required_columns <- c(
  "Sample", "Group", "Tissue", "Condition",
  "CDR3.aa", "V_family", "J_family",
  "Clones", "Region", "R", "S"
)

missing_columns <- setdiff(required_columns, names(region_rs))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

############################################################
## 3. Clean metadata and define region classes
############################################################

region_rs_clean <- region_rs %>%
  mutate(
    Sample = as.character(Sample),
    Group = as.character(Group),
    Tissue = as.character(Tissue),
    Condition = as.character(Condition),
    CDR3.aa = as.character(CDR3.aa),
    V_family = as.character(V_family),
    J_family = as.character(J_family),
    
    Tissue = case_when(
      Tissue %in% c("Sp", "SP", "Spleen") ~ "Spleen",
      Tissue %in% c("HK", "Head kidney") ~ "HK",
      TRUE ~ Tissue
    ),
    
    Region_class = case_when(
      Region %in% c("CDR1", "CDR2") ~ "CDR",
      Region %in% c("FR1", "FR2", "FR3") ~ "FR",
      TRUE ~ NA_character_
    ),
    
    Clones = as.numeric(Clones),
    R = as.numeric(R),
    S = as.numeric(S)
  ) %>%
  filter(
    !is.na(Sample),
    !is.na(CDR3.aa),
    CDR3.aa != "",
    !is.na(Clones),
    Clones > 0,
    !is.na(Region_class)
  )

############################################################
## 4. Create one unambiguous clone identifier
############################################################
## The same CDR3 can theoretically occur with different V/J
## assignments, so the identifier includes CDR3, V and J.

region_rs_clean <- region_rs_clean %>%
  mutate(
    Clone_ID = paste(
      CDR3.aa,
      V_family,
      J_family,
      sep = "__"
    )
  )

############################################################
## 5. Obtain one abundance record per clone
############################################################
## This step prevents Clones from being counted repeatedly
## across CDR1, CDR2, FR1, FR2 and FR3 rows.

clone_abundance <- region_rs_clean %>%
  distinct(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family,
    Clones
  ) %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family
  ) %>%
  summarise(
    clone_count = max(Clones, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 6. Calculate clone frequency within each sample
############################################################

clone_abundance <- clone_abundance %>%
  group_by(Sample) %>%
  mutate(
    repertoire_count = sum(clone_count, na.rm = TRUE),
    clone_frequency = clone_count / repertoire_count,
    clone_percentage = 100 * clone_frequency
  ) %>%
  ungroup()

############################################################
## 7. Assign biologically interpretable frequency classes
############################################################

clone_abundance <- clone_abundance %>%
  mutate(
    Clone_class = case_when(
      clone_frequency > 0.01   ~ "Expanded (>1%)",
      clone_frequency >= 0.001 ~ "Intermediate (0.1???1%)",
      clone_frequency >= 1e-4  ~ "Small (0.01???0.1%)",
      clone_frequency < 1e-4   ~ "Rare (<0.01%)",
      TRUE ~ NA_character_
    ),
    
    Clone_class = factor(
      Clone_class,
      levels = c(
        "Expanded (>1%)",
        "Intermediate (0.1???1%)",
        "Small (0.01???0.1%)",
        "Rare (<0.01%)"
      )
    )
  )

############################################################
## 8. Calculate clone-level mutation and selection metrics
############################################################
## R and S are used without clone weighting here because each
## row represents a mutation measurement for that clonotype.
## Clone abundance is analyzed separately as clone_frequency.

clone_mutation <- region_rs_clean %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family
  ) %>%
  summarise(
    total_R = sum(R, na.rm = TRUE),
    total_S = sum(S, na.rm = TRUE),
    
    CDR_R = sum(
      R[Region_class == "CDR"],
      na.rm = TRUE
    ),
    
    CDR_S = sum(
      S[Region_class == "CDR"],
      na.rm = TRUE
    ),
    
    FR_R = sum(
      R[Region_class == "FR"],
      na.rm = TRUE
    ),
    
    FR_S = sum(
      S[Region_class == "FR"],
      na.rm = TRUE
    ),
    
    total_mutations = total_R + total_S,
    
    CDR_mutations = CDR_R + CDR_S,
    FR_mutations = FR_R + FR_S,
    
    # The +1 stabilizes ratios when silent mutations are absent
    CDR_RS_ratio = CDR_R / (CDR_S + 1),
    FR_RS_ratio = FR_R / (FR_S + 1),
    
    .groups = "drop"
  )

############################################################
## 9. Combine clone abundance, frequency and SHM metrics
############################################################

clone_frequency_shm <- clone_abundance %>%
  left_join(
    clone_mutation,
    by = c(
      "Sample",
      "Group",
      "Tissue",
      "Condition",
      "Clone_ID",
      "CDR3.aa",
      "V_family",
      "J_family"
    )
  ) %>%
  mutate(
    log10_clone_count = log10(clone_count + 1),
    log10_clone_frequency = log10(clone_frequency + 1e-8),
    log10_total_mutations = log10(total_mutations + 1),
    log1p_CDR_RS = log1p(CDR_RS_ratio),
    
    Tissue = factor(
      Tissue,
      levels = c("HK", "Spleen")
    ),
    
    Condition = factor(
      Condition,
      levels = c("Control", "Infected")
    )
  ) %>%
  filter(!is.na(Clone_class))

write_csv(
  clone_frequency_shm,
  "clone_frequency_classes_SHM_selection.csv"
)

############################################################
## 10. Check the number of clones in each class
############################################################

class_counts <- clone_frequency_shm %>%
  count(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class,
    name = "n_clones"
  )

print(class_counts)

write_csv(
  class_counts,
  "clone_frequency_class_counts_per_sample.csv"
)

############################################################
## 11. Check repertoire occupancy by each clone class
############################################################

class_occupancy <- clone_frequency_shm %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_clones = n(),
    total_clone_count = sum(clone_count, na.rm = TRUE),
    repertoire_fraction = sum(
      clone_frequency,
      na.rm = TRUE
    ),
    repertoire_percentage = 100 * repertoire_fraction,
    .groups = "drop"
  )

print(class_occupancy)

write_csv(
  class_occupancy,
  "clone_class_repertoire_occupancy.csv"
)

############################################################
## 12. Figure A: repertoire space occupied by each class
############################################################

p_class_occupancy <- ggplot(
  class_occupancy,
  aes(
    x = Clone_class,
    y = repertoire_percentage,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.10,
    size = 2,
    alpha = 0.8
  ) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 13) +
  labs(
    title = "Repertoire occupancy by clone-frequency class",
    x = "Clone-frequency class",
    y = "Percentage of repertoire"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_class_occupancy)

ggsave(
  "Figure_clone_class_repertoire_occupancy.png",
  plot = p_class_occupancy,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 13. Figure B: SHM burden by frequency class
############################################################

p_shm_class <- ggplot(
  clone_frequency_shm,
  aes(
    x = Clone_class,
    y = total_mutations,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.14,
    size = 0.55,
    alpha = 0.18
  ) +
  facet_grid(Tissue ~ Condition) +
  coord_cartesian(
    ylim = c(
      0,
      quantile(
        clone_frequency_shm$total_mutations,
        0.99,
        na.rm = TRUE
      )
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "SHM burden across clone-frequency classes",
    subtitle = "Clone classes are based on frequency within each sample",
    x = "Clone-frequency class",
    y = "Total V-region mutations"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_shm_class)

ggsave(
  "Figure_SHM_by_clone_frequency_class.png",
  plot = p_shm_class,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 14. Figure C: CDR R/S ratio by frequency class
############################################################

p_cdr_rs_class <- ggplot(
  clone_frequency_shm,
  aes(
    x = Clone_class,
    y = CDR_RS_ratio,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.14,
    size = 0.55,
    alpha = 0.18
  ) +
  facet_grid(Tissue ~ Condition) +
  coord_cartesian(
    ylim = c(
      0,
      quantile(
        clone_frequency_shm$CDR_RS_ratio,
        0.99,
        na.rm = TRUE
      )
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "CDR selection across clone-frequency classes",
    x = "Clone-frequency class",
    y = "CDR replacement/silent ratio"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_cdr_rs_class)

ggsave(
  "Figure_CDR_RS_by_clone_frequency_class.png",
  plot = p_cdr_rs_class,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 15. Sample-level summaries
############################################################
## These are useful because samples, rather than clones, are
## the independent biological replicates.

sample_class_summary <- clone_frequency_shm %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_clones = n(),
    
    median_clone_frequency = median(
      clone_frequency,
      na.rm = TRUE
    ),
    
    median_mutations = median(
      total_mutations,
      na.rm = TRUE
    ),
    
    mean_mutations = mean(
      total_mutations,
      na.rm = TRUE
    ),
    
    median_CDR_RS = median(
      CDR_RS_ratio,
      na.rm = TRUE
    ),
    
    mean_CDR_RS = mean(
      CDR_RS_ratio,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write_csv(
  sample_class_summary,
  "sample_level_clone_class_SHM_selection_summary.csv"
)

############################################################
## 16. Publication-style sample-level SHM figure
############################################################

p_shm_sample <- ggplot(
  sample_class_summary,
  aes(
    x = Clone_class,
    y = median_mutations,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.10,
    size = 2.3,
    alpha = 0.9
  ) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 13) +
  labs(
    title = "Expanded clone classes carry a higher SHM burden",
    subtitle = "Each point represents the within-sample median",
    x = "Clone-frequency class",
    y = "Median V-region mutations per clone"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_shm_sample)

ggsave(
  "Figure_sample_level_SHM_by_clone_class.png",
  plot = p_shm_sample,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 17. Publication-style sample-level CDR R/S figure
############################################################

p_rs_sample <- ggplot(
  sample_class_summary,
  aes(
    x = Clone_class,
    y = median_CDR_RS,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.10,
    size = 2.3,
    alpha = 0.9
  ) +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 13) +
  labs(
    title = "CDR selection varies with clone-frequency class",
    subtitle = "Each point represents the within-sample median",
    x = "Clone-frequency class",
    y = "Median CDR R/S ratio"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "none"
  )

print(p_rs_sample)

ggsave(
  "Figure_sample_level_CDR_RS_by_clone_class.png",
  plot = p_rs_sample,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
## 18. Mixed model: SHM burden
############################################################
## log1p transformation reduces skew caused by highly
## mutated outlying clones.

model_shm <- lmer(
  log10_total_mutations ~
    Clone_class * Tissue * Condition +
    (1 | Sample),
  data = clone_frequency_shm,
  REML = FALSE
)

print(summary(model_shm))
print(anova(model_shm, type = 3))

capture.output(
  summary(model_shm),
  file = "LMM_SHM_by_clone_frequency_class_summary.txt"
)

write.csv(
  as.data.frame(anova(model_shm, type = 3)),
  "LMM_SHM_by_clone_frequency_class_ANOVA.csv"
)

############################################################
## 19. Holm-corrected SHM comparisons among classes
############################################################

emm_shm <- emmeans(
  model_shm,
  ~ Clone_class | Tissue * Condition
)

pairwise_shm <- pairs(
  emm_shm,
  adjust = "holm"
)

pairwise_shm_df <- as.data.frame(pairwise_shm)

print(pairwise_shm_df)

write_csv(
  pairwise_shm_df,
  "Pairwise_SHM_clone_classes_Holm.csv"
)

############################################################
## 20. Mixed model: CDR R/S ratio
############################################################

model_cdr_rs <- lmer(
  log1p_CDR_RS ~
    Clone_class * Tissue * Condition +
    (1 | Sample),
  data = clone_frequency_shm,
  REML = FALSE
)

print(summary(model_cdr_rs))
print(anova(model_cdr_rs, type = 3))

capture.output(
  summary(model_cdr_rs),
  file = "LMM_CDR_RS_by_clone_frequency_class_summary.txt"
)

write.csv(
  as.data.frame(anova(model_cdr_rs, type = 3)),
  "LMM_CDR_RS_by_clone_frequency_class_ANOVA.csv"
)

############################################################
## 21. Holm-corrected CDR R/S comparisons among classes
############################################################

emm_cdr_rs <- emmeans(
  model_cdr_rs,
  ~ Clone_class | Tissue * Condition
)

pairwise_cdr_rs <- pairs(
  emm_cdr_rs,
  adjust = "holm"
)

pairwise_cdr_rs_df <- as.data.frame(pairwise_cdr_rs)

print(pairwise_cdr_rs_df)

write_csv(
  pairwise_cdr_rs_df,
  "Pairwise_CDR_RS_clone_classes_Holm.csv"
)

############################################################
## 22. Direct expanded-versus-other comparisons
############################################################
## This contrast tests the most important biological question:
## are >1% clones different from all less-expanded classes?

expanded_contrast_shm <- contrast(
  emm_shm,
  method = list(
    "Expanded vs all lower-frequency classes" =
      c(1, -1/3, -1/3, -1/3)
  ),
  adjust = "holm"
)

expanded_contrast_rs <- contrast(
  emm_cdr_rs,
  method = list(
    "Expanded vs all lower-frequency classes" =
      c(1, -1/3, -1/3, -1/3)
  ),
  adjust = "holm"
)

expanded_contrast_shm_df <- as.data.frame(
  expanded_contrast_shm
)

expanded_contrast_rs_df <- as.data.frame(
  expanded_contrast_rs
)

print(expanded_contrast_shm_df)
print(expanded_contrast_rs_df)

write_csv(
  expanded_contrast_shm_df,
  "Contrast_expanded_vs_other_classes_SHM.csv"
)

write_csv(
  expanded_contrast_rs_df,
  "Contrast_expanded_vs_other_classes_CDR_RS.csv"
)

############################################################
## 23. Correlation between continuous clone frequency and SHM
############################################################

correlation_results <- clone_frequency_shm %>%
  group_by(Tissue, Condition) %>%
  group_modify(
    ~ {
      complete_data <- .x %>%
        filter(
          is.finite(log10_clone_frequency),
          is.finite(log10_total_mutations)
        )
      
      if (nrow(complete_data) < 3) {
        return(
          tibble(
            n = nrow(complete_data),
            rho = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      test_result <- cor.test(
        complete_data$log10_clone_frequency,
        complete_data$log10_total_mutations,
        method = "spearman",
        exact = FALSE
      )
      
      tibble(
        n = nrow(complete_data),
        rho = unname(test_result$estimate),
        p_value = test_result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

print(correlation_results)

write_csv(
  correlation_results,
  "Correlation_clone_frequency_vs_SHM.csv"
)

############################################################
## 24. Continuous frequency-versus-SHM figure
############################################################

p_frequency_shm <- ggplot(
  clone_frequency_shm,
  aes(
    x = clone_percentage,
    y = log10_total_mutations,
    color = Condition
  )
) +
  geom_point(
    alpha = 0.25,
    size = 1
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.9
  ) +
  scale_x_log10() +
  facet_wrap(~Tissue) +
  theme_bw(base_size = 13) +
  labs(
    title = "Clone frequency is positively associated with SHM burden",
    x = "Clone frequency within sample (%)",
    y = "log10(V-region mutations + 1)",
    color = "Condition"
  )

print(p_frequency_shm)

ggsave(
  "Figure_clone_frequency_vs_SHM.png",
  plot = p_frequency_shm,
  width = 9,
  height = 5.5,
  dpi = 300
)

############################################################
## 25. Identify expanded, highly mutated, CDR-selected clones
############################################################
## Thresholds for mutation burden and CDR selection are
## calculated within tissue and condition to avoid one group
## defining all candidates.

candidate_clones <- clone_frequency_shm %>%
  group_by(Tissue, Condition) %>%
  mutate(
    high_mutation_burden =
      total_mutations >= quantile(
        total_mutations,
        0.75,
        na.rm = TRUE
      ),
    
    high_CDR_selection =
      CDR_RS_ratio >= quantile(
        CDR_RS_ratio,
        0.75,
        na.rm = TRUE
      ),
    
    affinity_maturation_candidate =
      Clone_class == "Expanded (>1%)" &
      high_mutation_burden &
      high_CDR_selection
  ) %>%
  ungroup() %>%
  arrange(
    desc(affinity_maturation_candidate),
    desc(clone_frequency),
    desc(total_mutations),
    desc(CDR_RS_ratio)
  )

write_csv(
  candidate_clones,
  "candidate_expanded_mutated_CDR_selected_clones_frequency_based.csv"
)

############################################################
## 26. Candidate clone figure
############################################################

p_candidates <- ggplot(
  candidate_clones,
  aes(
    x = clone_percentage,
    y = CDR_RS_ratio
  )
) +
  geom_point(
    aes(
      color = affinity_maturation_candidate,
      size = total_mutations
    ),
    alpha = 0.55
  ) +
  scale_x_log10() +
  facet_grid(Tissue ~ Condition) +
  theme_bw(base_size = 13) +
  labs(
    title = paste(
      "Expanded, mutated and",
      "CDR-selected B-cell clones"
    ),
    subtitle = paste(
      "Candidates combine >1% repertoire frequency,",
      "high SHM burden and high CDR R/S"
    ),
    x = "Clone frequency within sample (%)",
    y = "CDR R/S ratio",
    color = "Candidate",
    size = "Total mutations"
  )

print(p_candidates)

ggsave(
  "Figure_frequency_based_affinity_maturation_candidates.png",
  plot = p_candidates,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
## 27. Final summary
############################################################

final_summary <- clone_frequency_shm %>%
  group_by(
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_clones = n(),
    median_clone_percentage = median(
      clone_percentage,
      na.rm = TRUE
    ),
    median_mutations = median(
      total_mutations,
      na.rm = TRUE
    ),
    median_CDR_RS = median(
      CDR_RS_ratio,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(final_summary)

write_csv(
  final_summary,
  "summary_clone_frequency_classes_SHM_CDR_selection.csv"
)

cat("\nAnalysis complete.\n")
cat("Frequency classes used:\n")
cat("Expanded: >1%\n")
cat("Intermediate: 0.1???1%\n")
cat("Small: 0.01???0.1%\n")
cat("Rare: <0.01%\n")




























###this data is not intersting so I didn't use in the paper

############################################################
## SAMPLE-LEVEL ANALYSIS:
## Clone frequency, SHM burden and CDR-focused selection
##
## Biological unit = individual fish/sample
##
## Clone-frequency classes:
## Expanded:      >1%
## Intermediate:  0.1???1%
## Small:         0.01???0.1%
## Rare:          <0.01%
############################################################


############################################################
## 1. Load packages
############################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


############################################################
## 2. Set working directory and load data
############################################################

workdir <- paste0(
  "C:/Users/sje054/OneDrive - UiT Office 365/",
  "BCR/SHM_mutation_spectrum_HK_Sp_3w"
)

setwd(workdir)

region_rs <- readRDS("region_rs_SAV_3wpi.rds")


############################################################
## 3. Check required columns
############################################################

required_columns <- c(
  "Sample",
  "Group",
  "Tissue",
  "Condition",
  "CDR3.aa",
  "V_family",
  "J_family",
  "Clones",
  "Region",
  "R",
  "S"
)

missing_columns <- setdiff(
  required_columns,
  names(region_rs)
)

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


############################################################
## 4. Clean metadata and define region classes
############################################################

region_rs_clean <- region_rs %>%
  mutate(
    Sample = as.character(Sample),
    Group = as.character(Group),
    Tissue = as.character(Tissue),
    Condition = as.character(Condition),
    CDR3.aa = as.character(CDR3.aa),
    V_family = as.character(V_family),
    J_family = as.character(J_family),
    Region = as.character(Region),
    
    Tissue = case_when(
      Tissue %in% c("Sp", "SP", "Spleen") ~ "Spleen",
      Tissue %in% c("HK", "Head kidney") ~ "HK",
      TRUE ~ Tissue
    ),
    
    Condition = case_when(
      Condition %in% c("C", "Control") ~ "Control",
      Condition %in% c("I", "Infected") ~ "Infected",
      TRUE ~ Condition
    ),
    
    Region_class = case_when(
      Region %in% c("CDR1", "CDR2") ~ "CDR",
      Region %in% c("FR1", "FR2", "FR3") ~ "FR",
      TRUE ~ NA_character_
    ),
    
    Clones = as.numeric(Clones),
    R = as.numeric(R),
    S = as.numeric(S)
  ) %>%
  filter(
    !is.na(Sample),
    !is.na(CDR3.aa),
    CDR3.aa != "",
    !is.na(Clones),
    Clones > 0,
    !is.na(Region_class)
  )


############################################################
## 5. Create a unique clone identifier
############################################################

region_rs_clean <- region_rs_clean %>%
  mutate(
    Clone_ID = paste(
      CDR3.aa,
      V_family,
      J_family,
      sep = "__"
    )
  )


############################################################
## 6. Obtain one abundance record per clonotype
############################################################
## Clones occurs on multiple region rows, so this prevents
## repeated counting of the same clone.

clone_abundance <- region_rs_clean %>%
  distinct(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family,
    Clones
  ) %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family
  ) %>%
  summarise(
    clone_count = max(
      Clones,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


############################################################
## 7. Calculate clone frequency within each sample
############################################################

clone_abundance <- clone_abundance %>%
  group_by(Sample) %>%
  mutate(
    repertoire_count = sum(
      clone_count,
      na.rm = TRUE
    ),
    
    clone_frequency =
      clone_count / repertoire_count,
    
    clone_percentage =
      100 * clone_frequency
  ) %>%
  ungroup()


############################################################
## 8. Assign clone-frequency classes
############################################################

clone_abundance <- clone_abundance %>%
  mutate(
    Clone_class = case_when(
      clone_frequency > 0.01 ~
        "Expanded (>1%)",
      
      clone_frequency >= 0.001 ~
        "Intermediate (0.1???1%)",
      
      clone_frequency >= 0.0001 ~
        "Small (0.01???0.1%)",
      
      clone_frequency < 0.0001 ~
        "Rare (<0.01%)",
      
      TRUE ~ NA_character_
    ),
    
    Clone_class = factor(
      Clone_class,
      levels = c(
        "Expanded (>1%)",
        "Intermediate (0.1???1%)",
        "Small (0.01???0.1%)",
        "Rare (<0.01%)"
      )
    )
  )


############################################################
## 9. Calculate clone-level mutation metrics
############################################################

clone_mutation <- region_rs_clean %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_ID,
    CDR3.aa,
    V_family,
    J_family
  ) %>%
  summarise(
    total_R = sum(
      R,
      na.rm = TRUE
    ),
    
    total_S = sum(
      S,
      na.rm = TRUE
    ),
    
    CDR_R = sum(
      R[Region_class == "CDR"],
      na.rm = TRUE
    ),
    
    CDR_S = sum(
      S[Region_class == "CDR"],
      na.rm = TRUE
    ),
    
    FR_R = sum(
      R[Region_class == "FR"],
      na.rm = TRUE
    ),
    
    FR_S = sum(
      S[Region_class == "FR"],
      na.rm = TRUE
    ),
    
    total_mutations =
      total_R + total_S,
    
    CDR_mutations =
      CDR_R + CDR_S,
    
    FR_mutations =
      FR_R + FR_S,
    
    CDR_RS_ratio =
      CDR_R / (CDR_S + 1),
    
    FR_RS_ratio =
      FR_R / (FR_S + 1),
    
    .groups = "drop"
  )


############################################################
## 10. Combine abundance and mutation data
############################################################

clone_frequency_shm <- clone_abundance %>%
  left_join(
    clone_mutation,
    by = c(
      "Sample",
      "Group",
      "Tissue",
      "Condition",
      "Clone_ID",
      "CDR3.aa",
      "V_family",
      "J_family"
    )
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "HKC",
        "HKI",
        "SpC",
        "SpI"
      )
    ),
    
    Tissue = factor(
      Tissue,
      levels = c(
        "HK",
        "Spleen"
      )
    ),
    
    Condition = factor(
      Condition,
      levels = c(
        "Control",
        "Infected"
      )
    ),
    
    log10_mutations =
      log10(total_mutations + 1),
    
    log1p_CDR_RS =
      log1p(CDR_RS_ratio)
  ) %>%
  filter(
    !is.na(Clone_class),
    !is.na(total_mutations),
    !is.na(CDR_RS_ratio)
  )

write_csv(
  clone_frequency_shm,
  "clone_level_frequency_SHM_selection.csv"
)


############################################################
## 11. Inspect clone-class representation
############################################################

clone_class_counts <- clone_frequency_shm %>%
  count(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class,
    name = "n_clones"
  )

print(clone_class_counts)

write_csv(
  clone_class_counts,
  "clone_class_counts_per_sample.csv"
)


############################################################
## 12. Calculate repertoire occupancy
############################################################

class_occupancy <- clone_frequency_shm %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_clones = n(),
    
    repertoire_percentage =
      100 * sum(
        clone_frequency,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write_csv(
  class_occupancy,
  "clone_class_repertoire_occupancy.csv"
)


############################################################
## 13. Create sample-level summaries
############################################################
## Each row now represents one fish ?? one clone class.

sample_class_summary <- clone_frequency_shm %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_clones = n(),
    
    repertoire_percentage =
      100 * sum(
        clone_frequency,
        na.rm = TRUE
      ),
    
    median_clone_percentage = median(
      clone_percentage,
      na.rm = TRUE
    ),
    
    median_total_mutations = median(
      total_mutations,
      na.rm = TRUE
    ),
    
    mean_total_mutations = mean(
      total_mutations,
      na.rm = TRUE
    ),
    
    median_CDR_mutations = median(
      CDR_mutations,
      na.rm = TRUE
    ),
    
    median_FR_mutations = median(
      FR_mutations,
      na.rm = TRUE
    ),
    
    median_CDR_RS = median(
      CDR_RS_ratio,
      na.rm = TRUE
    ),
    
    mean_CDR_RS = mean(
      CDR_RS_ratio,
      na.rm = TRUE
    ),
    
    median_FR_RS = median(
      FR_RS_ratio,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

write_csv(
  sample_class_summary,
  "sample_level_clone_class_SHM_selection_summary.csv"
)

print(sample_class_summary)


############################################################
## 14. Make sure absent classes are visible
############################################################
## If a sample lacks one class, complete() inserts that
## combination. Statistical tests use only observed values.

sample_class_complete <- sample_class_summary %>%
  complete(
    Sample,
    Clone_class,
    nesting(
      Group,
      Tissue,
      Condition
    )
  )


############################################################
## 15. Figure A: repertoire occupancy
############################################################

p_occupancy <- ggplot(
  sample_class_summary,
  aes(
    x = Clone_class,
    y = repertoire_percentage,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.5,
    alpha = 0.9
  ) +
  facet_grid(
    Tissue ~ Condition
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "Repertoire occupancy by clone-frequency class",
    
    subtitle =
      "Each point represents one fish",
    
    x =
      "Clone-frequency class",
    
    y =
      "Percentage of repertoire"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    
    legend.position = "none"
  )

print(p_occupancy)

ggsave(
  "Figure_sample_level_clone_class_occupancy.png",
  plot = p_occupancy,
  width = 11,
  height = 7,
  dpi = 300
)


############################################################
## 16. Figure B: sample-level SHM burden
############################################################

p_shm_sample <- ggplot(
  sample_class_summary,
  aes(
    x = Clone_class,
    y = median_total_mutations,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.5,
    alpha = 0.9
  ) +
  facet_grid(
    Tissue ~ Condition
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "Expanded clone classes carry a higher SHM burden",
    
    subtitle =
      "Each point represents the within-fish median",
    
    x =
      "Clone-frequency class",
    
    y =
      "Median V-region mutations per clonotype"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    
    legend.position = "none"
  )

print(p_shm_sample)

ggsave(
  "Figure_sample_level_SHM_by_clone_frequency_class.png",
  plot = p_shm_sample,
  width = 11,
  height = 7,
  dpi = 300
)


############################################################
## 17. Figure C: sample-level CDR R/S ratio
############################################################

p_cdr_rs_sample <- ggplot(
  sample_class_summary,
  aes(
    x = Clone_class,
    y = median_CDR_RS,
    fill = Clone_class
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.5,
    alpha = 0.9
  ) +
  facet_grid(
    Tissue ~ Condition
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "CDR-focused selection varies across clone-frequency classes",
    
    subtitle =
      "Each point represents the within-fish median",
    
    x =
      "Clone-frequency class",
    
    y =
      "Median CDR replacement/silent ratio"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    
    legend.position = "none"
  )

print(p_cdr_rs_sample)

ggsave(
  "Figure_sample_level_CDR_RS_by_clone_frequency_class.png",
  plot = p_cdr_rs_sample,
  width = 11,
  height = 7,
  dpi = 300
)


############################################################
## 18. Helper function: Friedman test within each group
############################################################
## Clone class is the repeated factor.
## Sample is the blocking factor.

run_friedman_test <- function(
    data,
    response_variable
) {
  
  complete_data <- data %>%
    select(
      Sample,
      Clone_class,
      all_of(response_variable)
    ) %>%
    pivot_wider(
      names_from = Clone_class,
      values_from = all_of(response_variable)
    ) %>%
    drop_na()
  
  if (
    nrow(complete_data) < 3 ||
    ncol(complete_data) < 3
  ) {
    return(
      tibble(
        n_complete_samples = nrow(complete_data),
        statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  long_complete <- complete_data %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Clone_class",
      values_to = "response"
    )
  
  test_result <- friedman.test(
    response ~ Clone_class | Sample,
    data = long_complete
  )
  
  tibble(
    n_complete_samples = nrow(complete_data),
    statistic = unname(
      test_result$statistic
    ),
    df = unname(
      test_result$parameter
    ),
    p_value = test_result$p.value
  )
}


############################################################
## 19. Friedman tests for SHM across clone classes
############################################################

friedman_shm <- sample_class_summary %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ run_friedman_test(
      .x,
      "median_total_mutations"
    )
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(friedman_shm)

write_csv(
  friedman_shm,
  "Friedman_SHM_across_clone_classes.csv"
)


############################################################
## 20. Friedman tests for CDR R/S across clone classes
############################################################

friedman_cdr_rs <- sample_class_summary %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ run_friedman_test(
      .x,
      "median_CDR_RS"
    )
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(friedman_cdr_rs)

write_csv(
  friedman_cdr_rs,
  "Friedman_CDR_RS_across_clone_classes.csv"
)


############################################################
## 21. Helper function:
## paired Wilcoxon comparisons between clone classes
############################################################

run_paired_wilcox <- function(
    data,
    response_variable
) {
  
  classes <- levels(
    droplevels(
      data$Clone_class
    )
  )
  
  comparisons <- combn(
    classes,
    2,
    simplify = FALSE
  )
  
  results <- lapply(
    comparisons,
    function(pair) {
      
      paired_data <- data %>%
        filter(
          Clone_class %in% pair
        ) %>%
        select(
          Sample,
          Clone_class,
          all_of(response_variable)
        ) %>%
        pivot_wider(
          names_from = Clone_class,
          values_from = all_of(response_variable)
        ) %>%
        drop_na()
      
      if (nrow(paired_data) < 3) {
        return(
          tibble(
            class_1 = pair[1],
            class_2 = pair[2],
            n_pairs = nrow(paired_data),
            statistic = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      test_result <- wilcox.test(
        paired_data[[pair[1]]],
        paired_data[[pair[2]]],
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        class_1 = pair[1],
        class_2 = pair[2],
        n_pairs = nrow(paired_data),
        statistic = unname(
          test_result$statistic
        ),
        p_value = test_result$p.value
      )
    }
  )
  
  bind_rows(results) %>%
    mutate(
      adjusted_p_value = p.adjust(
        p_value,
        method = "holm"
      )
    )
}


############################################################
## 22. Pairwise SHM comparisons among clone classes
############################################################

pairwise_shm_classes <- sample_class_summary %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ run_paired_wilcox(
      .x,
      "median_total_mutations"
    )
  ) %>%
  ungroup()

print(pairwise_shm_classes)

write_csv(
  pairwise_shm_classes,
  "Paired_Wilcoxon_SHM_clone_classes_Holm.csv"
)


############################################################
## 23. Pairwise CDR R/S comparisons among clone classes
############################################################

pairwise_cdr_rs_classes <- sample_class_summary %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ run_paired_wilcox(
      .x,
      "median_CDR_RS"
    )
  ) %>%
  ungroup()

print(pairwise_cdr_rs_classes)

write_csv(
  pairwise_cdr_rs_classes,
  "Paired_Wilcoxon_CDR_RS_clone_classes_Holm.csv"
)


############################################################
## 24. Direct comparison:
## Expanded versus pooled lower-frequency classes
############################################################
## For each fish, compare the expanded-clone median with the
## median across intermediate, small and rare clone classes.

expanded_vs_other <- sample_class_summary %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition
  ) %>%
  summarise(
    expanded_SHM = median_total_mutations[
      Clone_class == "Expanded (>1%)"
    ][1],
    
    lower_frequency_SHM = median(
      median_total_mutations[
        Clone_class != "Expanded (>1%)"
      ],
      na.rm = TRUE
    ),
    
    expanded_CDR_RS = median_CDR_RS[
      Clone_class == "Expanded (>1%)"
    ][1],
    
    lower_frequency_CDR_RS = median(
      median_CDR_RS[
        Clone_class != "Expanded (>1%)"
      ],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


############################################################
## 25. Expanded versus other classes:
## paired Wilcoxon test for SHM
############################################################

expanded_vs_other_shm <- expanded_vs_other %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          is.finite(expanded_SHM),
          is.finite(lower_frequency_SHM)
        )
      
      if (nrow(test_data) < 3) {
        return(
          tibble(
            n_pairs = nrow(test_data),
            statistic = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      result <- wilcox.test(
        test_data$expanded_SHM,
        test_data$lower_frequency_SHM,
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        n_pairs = nrow(test_data),
        statistic = unname(
          result$statistic
        ),
        p_value = result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(expanded_vs_other_shm)

write_csv(
  expanded_vs_other_shm,
  "Expanded_vs_lower_frequency_SHM.csv"
)


############################################################
## 26. Expanded versus other classes:
## paired Wilcoxon test for CDR R/S
############################################################

expanded_vs_other_cdr_rs <- expanded_vs_other %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          is.finite(expanded_CDR_RS),
          is.finite(lower_frequency_CDR_RS)
        )
      
      if (nrow(test_data) < 3) {
        return(
          tibble(
            n_pairs = nrow(test_data),
            statistic = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      result <- wilcox.test(
        test_data$expanded_CDR_RS,
        test_data$lower_frequency_CDR_RS,
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        n_pairs = nrow(test_data),
        statistic = unname(
          result$statistic
        ),
        p_value = result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(expanded_vs_other_cdr_rs)

write_csv(
  expanded_vs_other_cdr_rs,
  "Expanded_vs_lower_frequency_CDR_RS.csv"
)


############################################################
## 27. Compare the four biological groups within each class
############################################################
## HKC, HKI, SpC and SpI are compared separately for each
## clone class using Kruskal-Wallis tests.

group_test_shm <- sample_class_summary %>%
  group_by(Clone_class) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          !is.na(Group),
          is.finite(median_total_mutations)
        )
      
      if (
        n_distinct(test_data$Group) < 2 ||
        nrow(test_data) < 4
      ) {
        return(
          tibble(
            statistic = NA_real_,
            df = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      result <- kruskal.test(
        median_total_mutations ~ Group,
        data = test_data
      )
      
      tibble(
        statistic = unname(
          result$statistic
        ),
        df = unname(
          result$parameter
        ),
        p_value = result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(group_test_shm)

write_csv(
  group_test_shm,
  "Kruskal_Wallis_group_effect_SHM_by_clone_class.csv"
)


############################################################
## 28. Group effects on CDR R/S within each class
############################################################

group_test_cdr_rs <- sample_class_summary %>%
  group_by(Clone_class) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          !is.na(Group),
          is.finite(median_CDR_RS)
        )
      
      if (
        n_distinct(test_data$Group) < 2 ||
        nrow(test_data) < 4
      ) {
        return(
          tibble(
            statistic = NA_real_,
            df = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      result <- kruskal.test(
        median_CDR_RS ~ Group,
        data = test_data
      )
      
      tibble(
        statistic = unname(
          result$statistic
        ),
        df = unname(
          result$parameter
        ),
        p_value = result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "holm"
    )
  )

print(group_test_cdr_rs)

write_csv(
  group_test_cdr_rs,
  "Kruskal_Wallis_group_effect_CDR_RS_by_clone_class.csv"
)


############################################################
## 29. Pairwise group comparisons for SHM
############################################################

pairwise_group_shm <- sample_class_summary %>%
  group_by(Clone_class) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          !is.na(Group),
          is.finite(median_total_mutations)
        )
      
      if (
        n_distinct(test_data$Group) < 2
      ) {
        return(
          tibble()
        )
      }
      
      pw <- pairwise.wilcox.test(
        x = test_data$median_total_mutations,
        g = test_data$Group,
        p.adjust.method = "holm",
        exact = FALSE
      )
      
      as.data.frame(
        as.table(pw$p.value)
      ) %>%
        filter(
          !is.na(Freq)
        ) %>%
        transmute(
          group_1 = Var1,
          group_2 = Var2,
          adjusted_p_value = Freq
        )
    }
  ) %>%
  ungroup()

print(pairwise_group_shm)

write_csv(
  pairwise_group_shm,
  "Pairwise_group_SHM_by_clone_class_Holm.csv"
)


############################################################
## 30. Pairwise group comparisons for CDR R/S
############################################################

pairwise_group_cdr_rs <- sample_class_summary %>%
  group_by(Clone_class) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          !is.na(Group),
          is.finite(median_CDR_RS)
        )
      
      if (
        n_distinct(test_data$Group) < 2
      ) {
        return(
          tibble()
        )
      }
      
      pw <- pairwise.wilcox.test(
        x = test_data$median_CDR_RS,
        g = test_data$Group,
        p.adjust.method = "holm",
        exact = FALSE
      )
      
      as.data.frame(
        as.table(pw$p.value)
      ) %>%
        filter(
          !is.na(Freq)
        ) %>%
        transmute(
          group_1 = Var1,
          group_2 = Var2,
          adjusted_p_value = Freq
        )
    }
  ) %>%
  ungroup()

print(pairwise_group_cdr_rs)

write_csv(
  pairwise_group_cdr_rs,
  "Pairwise_group_CDR_RS_by_clone_class_Holm.csv"
)


############################################################
## 31. Calculate clone-frequency versus SHM correlation
## separately for each fish
############################################################
## This avoids treating all clonotypes from all fish as one
## independent population.

sample_correlations <- clone_frequency_shm %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition
  ) %>%
  group_modify(
    ~ {
      
      test_data <- .x %>%
        filter(
          clone_frequency > 0,
          is.finite(total_mutations)
        )
      
      if (nrow(test_data) < 5) {
        return(
          tibble(
            n_clones = nrow(test_data),
            spearman_rho = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      result <- cor.test(
        log10(
          test_data$clone_frequency
        ),
        log10(
          test_data$total_mutations + 1
        ),
        method = "spearman",
        exact = FALSE
      )
      
      tibble(
        n_clones = nrow(test_data),
        spearman_rho = unname(
          result$estimate
        ),
        p_value = result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    adjusted_p_value = p.adjust(
      p_value,
      method = "BH"
    )
  )

print(sample_correlations)

write_csv(
  sample_correlations,
  "sample_level_clone_frequency_SHM_correlations.csv"
)


############################################################
## 32. Plot sample-level correlation coefficients
############################################################

p_correlations <- ggplot(
  sample_correlations,
  aes(
    x = Group,
    y = spearman_rho,
    fill = Group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    width = 0.60,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.7
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "Association between clone frequency and SHM within individual fish",
    
    subtitle =
      "Each point is a sample-level Spearman correlation",
    
    x =
      "Group",
    
    y =
      "Spearman correlation coefficient"
  ) +
  theme(
    legend.position = "none"
  )

print(p_correlations)

ggsave(
  "Figure_sample_level_clone_frequency_SHM_correlations.png",
  plot = p_correlations,
  width = 7,
  height = 5.5,
  dpi = 300
)


############################################################
## 33. Kruskal-Wallis comparison of sample correlations
############################################################

correlation_group_test <- kruskal.test(
  spearman_rho ~ Group,
  data = sample_correlations
)

correlation_group_result <- tibble(
  statistic = unname(
    correlation_group_test$statistic
  ),
  df = unname(
    correlation_group_test$parameter
  ),
  p_value = correlation_group_test$p.value
)

print(correlation_group_result)

write_csv(
  correlation_group_result,
  "Kruskal_Wallis_sample_correlations_among_groups.csv"
)


############################################################
## 34. Identify candidate expanded, mutated and selected clones
############################################################
## Thresholds are defined within each fish so that unusually
## mutated and selected clones are identified relative to
## that fish's own repertoire.

candidate_clones <- clone_frequency_shm %>%
  group_by(Sample) %>%
  mutate(
    high_mutation_burden =
      total_mutations >= quantile(
        total_mutations,
        0.75,
        na.rm = TRUE
      ),
    
    high_CDR_RS =
      CDR_RS_ratio >= quantile(
        CDR_RS_ratio,
        0.75,
        na.rm = TRUE
      ),
    
    candidate_matured_clone =
      Clone_class == "Expanded (>1%)" &
      high_mutation_burden &
      high_CDR_RS
  ) %>%
  ungroup() %>%
  arrange(
    desc(candidate_matured_clone),
    desc(clone_frequency),
    desc(total_mutations),
    desc(CDR_RS_ratio)
  )

write_csv(
  candidate_clones,
  "candidate_expanded_mutated_CDR_selected_clones.csv"
)


############################################################
## 35. Count candidate clones per fish
############################################################

candidate_counts <- candidate_clones %>%
  group_by(
    Sample,
    Group,
    Tissue,
    Condition
  ) %>%
  summarise(
    total_clones = n(),
    
    candidate_clones = sum(
      candidate_matured_clone,
      na.rm = TRUE
    ),
    
    candidate_percentage =
      100 * candidate_clones / total_clones,
    
    .groups = "drop"
  )

print(candidate_counts)

write_csv(
  candidate_counts,
  "candidate_matured_clone_counts_per_sample.csv"
)


############################################################
## 36. Plot candidate clone counts per fish
############################################################

p_candidate_counts <- ggplot(
  candidate_counts,
  aes(
    x = Group,
    y = candidate_clones,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.60,
    outlier.shape = NA,
    alpha = 0.80
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.7
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "Candidate expanded, mutated and CDR-selected clones",
    
    subtitle =
      "Each point represents one fish",
    
    x =
      "Group",
    
    y =
      "Number of candidate clonotypes"
  ) +
  theme(
    legend.position = "none"
  )

print(p_candidate_counts)

ggsave(
  "Figure_candidate_matured_clone_counts_per_sample.png",
  plot = p_candidate_counts,
  width = 7,
  height = 5.5,
  dpi = 300
)


############################################################
## 37. Integrated clone-level visualization
############################################################
## This plot remains descriptive. Statistical inference is
## based on the sample-level summaries above.

p_integrated <- ggplot(
  candidate_clones,
  aes(
    x = clone_percentage,
    y = CDR_RS_ratio
  )
) +
  geom_point(
    aes(
      color = candidate_matured_clone,
      size = total_mutations
    ),
    alpha = 0.50
  ) +
  scale_x_log10() +
  facet_grid(
    Tissue ~ Condition
  ) +
  theme_bw(
    base_size = 13
  ) +
  labs(
    title =
      "Expanded, mutated and CDR-selected B-cell clones",
    
    subtitle =
      paste(
        "Highlighted candidates occupy >1% of the repertoire",
        "and have high within-sample SHM and CDR R/S values"
      ),
    
    x =
      "Clone frequency within sample (%)",
    
    y =
      "CDR replacement/silent ratio",
    
    color =
      "Candidate clone",
    
    size =
      "Total mutations"
  )

print(p_integrated)

ggsave(
  "Figure_integrated_frequency_SHM_CDR_selection.png",
  plot = p_integrated,
  width = 10,
  height = 7,
  dpi = 300
)


############################################################
## 38. Final manuscript summary table
############################################################

manuscript_summary <- sample_class_summary %>%
  group_by(
    Group,
    Tissue,
    Condition,
    Clone_class
  ) %>%
  summarise(
    n_samples = n_distinct(Sample),
    
    median_repertoire_percentage = median(
      repertoire_percentage,
      na.rm = TRUE
    ),
    
    lower_repertoire_percentage = quantile(
      repertoire_percentage,
      0.25,
      na.rm = TRUE
    ),
    
    upper_repertoire_percentage = quantile(
      repertoire_percentage,
      0.75,
      na.rm = TRUE
    ),
    
    median_SHM = median(
      median_total_mutations,
      na.rm = TRUE
    ),
    
    lower_SHM = quantile(
      median_total_mutations,
      0.25,
      na.rm = TRUE
    ),
    
    upper_SHM = quantile(
      median_total_mutations,
      0.75,
      na.rm = TRUE
    ),
    
    median_CDR_RS = median(
      median_CDR_RS,
      na.rm = TRUE
    ),
    
    lower_CDR_RS = quantile(
      median_CDR_RS,
      0.25,
      na.rm = TRUE
    ),
    
    upper_CDR_RS = quantile(
      median_CDR_RS,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(manuscript_summary)

write_csv(
  manuscript_summary,
  "manuscript_summary_clone_frequency_SHM_selection.csv"
)


############################################################
## 39. Save all important R objects
############################################################

saveRDS(
  clone_frequency_shm,
  "clone_frequency_shm_complete.rds"
)

saveRDS(
  sample_class_summary,
  "sample_class_summary.rds"
)

saveRDS(
  candidate_clones,
  "candidate_matured_clones.rds"
)


############################################################
## 40. Finish
############################################################

cat("\nAnalysis complete.\n\n")

cat("The main statistical outputs are:\n")
cat("- Friedman_SHM_across_clone_classes.csv\n")
cat("- Friedman_CDR_RS_across_clone_classes.csv\n")
cat("- Paired_Wilcoxon_SHM_clone_classes_Holm.csv\n")
cat("- Paired_Wilcoxon_CDR_RS_clone_classes_Holm.csv\n")
cat("- Expanded_vs_lower_frequency_SHM.csv\n")
cat("- Expanded_vs_lower_frequency_CDR_RS.csv\n")
cat("- sample_level_clone_frequency_SHM_correlations.csv\n\n")

cat("The recommended manuscript figures are:\n")
cat("- Figure_sample_level_SHM_by_clone_frequency_class.png\n")
cat("- Figure_sample_level_CDR_RS_by_clone_frequency_class.png\n")
cat("- Figure_sample_level_clone_frequency_SHM_correlations.png\n")
cat("- Figure_candidate_matured_clone_counts_per_sample.png\n")
cat("- Figure_integrated_frequency_SHM_CDR_selection.png\n")






########################################New analysis RS after the reviewer comment
## ---- 1. User settings -------------------------------------------------------

INPUT_DIR  <- "C:/Users/sje054/OneDrive - UiT Office 365/BCR/reviewer_SHM_RS_outputs"
OUTPUT_DIR <- file.path(INPUT_DIR, "integrated_abundance_SHM_analysis")

##Now create the output folder:
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##Confirm both folders are accessible:
dir.exists(INPUT_DIR)
dir.exists(OUTPUT_DIR)

##Now define and verify the three input files:
correlation_file <- file.path(
  INPUT_DIR,
  "sample_level_clone_size_SHM_correlations.csv"
)

clone_file <- file.path(
  INPUT_DIR,
  "clone_level_SHM_summary.csv"
)

bin_file <- file.path(
  INPUT_DIR,
  "sample_abundance_bin_SHM_corrected_RS.csv"
)

file.exists(correlation_file)
file.exists(clone_file)
file.exists(bin_file)


list.files(
  INPUT_DIR,
  full.names = FALSE
)

###Now import the datasets:
library(readr)

cor_df <- read_csv(
  correlation_file,
  show_col_types = FALSE
)

clone_df <- read_csv(
  clone_file,
  show_col_types = FALSE
)

bin_df <- read_csv(
  bin_file,
  show_col_types = FALSE
)

##Then inspect their dimensions and column names:
dim(cor_df)
names(cor_df)

dim(clone_df)
names(clone_df)

dim(bin_df)
names(bin_df)

##Step 2 ??? Prepare the data and verify the sample-level correlations
###############################################################################
# Step 2: Data preparation and verification of sample-level correlations
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)


## 1. Choose the sample set

# Use "all16" for all available samples.
# Use "balanced12" for the balanced manuscript subset.

#ANALYSIS_SET <- "all16"

BALANCED_12 <- c(
  "HKC_3C3w", "HKC_4C3w", "HKC_5C3w",
  "HKI_1I3w", "HKI_4I3w", "HKI_5I3w",
  "SpC_1C3w", "SpC_3C3w", "SpC_5C3w",
  "SpI_1I3w", "SpI_2I3w", "SpI_5I3w"
)


## 2. Standardize metadata columns

cor_df <- cor_df %>%
  mutate(
    Sample_clean = as.character(Sample_clean),
    Group = factor(
      Group,
      levels = c("HKC", "HKI", "SpC", "SpI")
    ),
    Tissue = factor(
      Tissue,
      levels = c("HK", "Spleen")
    ),
    Condition = factor(
      Condition,
      levels = c("Control", "Infected")
    ),
    n_clonotypes = as.numeric(n_clonotypes),
    spearman_rho = as.numeric(spearman_rho)
  )

clone_df <- clone_df %>%
  mutate(
    Sample_clean = as.character(Sample_clean),
    Group = factor(
      Group,
      levels = c("HKC", "HKI", "SpC", "SpI")
    ),
    Tissue = factor(
      Tissue,
      levels = c("HK", "Spleen")
    ),
    Condition = factor(
      Condition,
      levels = c("Control", "Infected")
    ),
    Clones = as.numeric(Clones),
    total_mutations = as.numeric(total_mutations),
    log10_clone_size = as.numeric(log10_clone_size),
    log10_mutations = as.numeric(log10_mutations),
    CDR_RS_corrected = as.numeric(CDR_RS_corrected),
    FR_RS_corrected = as.numeric(FR_RS_corrected)
  )

bin_df <- bin_df %>%
  mutate(
    Sample_clean = as.character(Sample_clean),
    Group = factor(
      Group,
      levels = c("HKC", "HKI", "SpC", "SpI")
    ),
    Tissue = factor(
      Tissue,
      levels = c("HK", "Spleen")
    ),
    Condition = factor(
      Condition,
      levels = c("Control", "Infected")
    ),
    n_clonotypes = as.numeric(n_clonotypes),
    median_mutations = as.numeric(median_mutations),
    mean_mutations = as.numeric(mean_mutations),
    CDR_RS_corrected = as.numeric(CDR_RS_corrected)
  )


## 3. Apply the selected sample set

if (ANALYSIS_SET == "balanced12") {
  
  cor_analysis <- cor_df %>%
    filter(Sample_clean %in% BALANCED_12)
  
  clone_analysis <- clone_df %>%
    filter(Sample_clean %in% BALANCED_12)
  
  bin_analysis <- bin_df %>%
    filter(Sample_clean %in% BALANCED_12)
  
} else if (ANALYSIS_SET == "all16") {
  
  cor_analysis <- cor_df
  clone_analysis <- clone_df
  bin_analysis <- bin_df
  
} else {
  
  stop("ANALYSIS_SET must be either 'all16' or 'balanced12'.")
}


## 4. Confirm the samples included in each table

cat("\nSamples in the correlation table:\n")
print(sort(unique(cor_analysis$Sample_clean)))

cat("\nSamples in the clone-level table:\n")
print(sort(unique(clone_analysis$Sample_clean)))

cat("\nSamples in the abundance-bin table:\n")
print(sort(unique(bin_analysis$Sample_clean)))

cat("\nNumber of samples in each table:\n")

sample_count_check <- data.frame(
  Dataset = c(
    "Sample-level correlations",
    "Clone-level data",
    "Abundance-bin data"
  ),
  Number_of_samples = c(
    n_distinct(cor_analysis$Sample_clean),
    n_distinct(clone_analysis$Sample_clean),
    n_distinct(bin_analysis$Sample_clean)
  )
)

print(sample_count_check)


## 5. Check sample agreement across the three datasets

sample_agreement <- list(
  Missing_from_clone_table = setdiff(
    cor_analysis$Sample_clean,
    clone_analysis$Sample_clean
  ),
  Missing_from_correlation_table = setdiff(
    clone_analysis$Sample_clean,
    cor_analysis$Sample_clean
  ),
  Missing_from_bin_table = setdiff(
    cor_analysis$Sample_clean,
    bin_analysis$Sample_clean
  )
)

print(sample_agreement)


## 6. Check missing values in important clone-level variables

clone_missing_summary <- clone_analysis %>%
  summarise(
    total_rows = n(),
    missing_sample = sum(is.na(Sample_clean)),
    missing_clone_size = sum(is.na(Clones)),
    missing_total_mutations = sum(is.na(total_mutations)),
    missing_log_clone_size = sum(is.na(log10_clone_size)),
    missing_log_mutations = sum(is.na(log10_mutations)),
    missing_CDR_RS = sum(is.na(CDR_RS_corrected)),
    missing_FR_RS = sum(is.na(FR_RS_corrected))
  )

print(clone_missing_summary)


## 7. Check whether clone IDs are duplicated within samples

duplicate_clone_check <- clone_analysis %>%
  count(
    Sample_clean,
    observed_clone_id,
    name = "number_of_rows"
  ) %>%
  filter(number_of_rows > 1)

cat("\nNumber of duplicated sample???clone combinations:\n")
print(nrow(duplicate_clone_check))

if (nrow(duplicate_clone_check) > 0) {
  print(head(duplicate_clone_check, 20))
}


## 8. Create the analysis-ready clone-level table

clone_valid <- clone_analysis %>%
  filter(
    !is.na(Sample_clean),
    !is.na(Clones),
    Clones > 0,
    !is.na(total_mutations)
  )

cat("\nNumber of retained clone-level observations:\n")
print(nrow(clone_valid))


## 9. Recalculate Spearman correlations within each fish

recalculated_correlations <- clone_valid %>%
  group_by(
    Sample_clean,
    Group,
    Tissue,
    Condition
  ) %>%
  summarise(
    n_clonotypes_recalculated = n(),
    
    spearman_rho_recalculated = {
      if (
        n() >= 3 &&
        n_distinct(Clones) > 1 &&
        n_distinct(total_mutations) > 1
      ) {
        cor(
          Clones,
          total_mutations,
          method = "spearman",
          use = "complete.obs"
        )
      } else {
        NA_real_
      }
    },
    
    spearman_p_value = {
      if (
        n() >= 3 &&
        n_distinct(Clones) > 1 &&
        n_distinct(total_mutations) > 1
      ) {
        suppressWarnings(
          cor.test(
            Clones,
            total_mutations,
            method = "spearman",
            exact = FALSE
          )$p.value
        )
      } else {
        NA_real_
      }
    },
    
    .groups = "drop"
  )

print(recalculated_correlations)

## 10. Compare recalculated correlations with exported values

correlation_verification <- cor_analysis %>%
  dplyr::select(
    Sample_clean,
    n_clonotypes_exported = n_clonotypes,
    spearman_rho_exported = spearman_rho
  ) %>%
  left_join(
    recalculated_correlations %>%
      dplyr::select(
        Sample_clean,
        n_clonotypes_recalculated,
        spearman_rho_recalculated,
        spearman_p_value
      ),
    by = "Sample_clean"
  ) %>%
  mutate(
    difference_in_n =
      n_clonotypes_recalculated - n_clonotypes_exported,
    
    difference_in_rho =
      spearman_rho_recalculated - spearman_rho_exported,
    
    correlations_match = case_when(
      is.na(spearman_rho_exported) &
        is.na(spearman_rho_recalculated) ~ TRUE,
      
      !is.na(difference_in_rho) &
        abs(difference_in_rho) < 1e-10 ~ TRUE,
      
      TRUE ~ FALSE
    )
  )

print(
  correlation_verification,
  n = nrow(correlation_verification)
)

## 11. Summarize the verification

verification_summary <- correlation_verification %>%
  summarise(
    number_of_samples = n(),
    samples_with_matching_n = sum(
      difference_in_n == 0,
      na.rm = TRUE
    ),
    samples_with_matching_rho = sum(
      correlations_match,
      na.rm = TRUE
    ),
    maximum_absolute_rho_difference = max(
      abs(difference_in_rho),
      na.rm = TRUE
    )
  )

print(verification_summary)


## 12. Save the verified correlation table

write.csv(
  correlation_verification,
  file.path(
    OUTPUT_DIR,
    paste0(
      "verified_sample_level_correlations_",
      ANALYSIS_SET,
      ".csv"
    )
  ),
  row.names = FALSE
)

## 12. Save the verified correlation table

write.csv(
  correlation_verification,
  file.path(
    OUTPUT_DIR,
    paste0(
      "verified_sample_level_correlations_",
      ANALYSIS_SET,
      ".csv"
    )
  ),
  row.names = FALSE
)


print(verification_summary)

print(
  correlation_verification,
  n = nrow(correlation_verification)
)


##
###############################################################################
# Step 3: Primary sample-level statistical analysis
###############################################################################

library(dplyr)
library(ggplot2)


## 1. Create the verified analysis table

rho_analysis <- recalculated_correlations %>%
  filter(!is.na(spearman_rho_recalculated)) %>%
  mutate(
    Group = factor(Group, levels = c("HKC", "HKI", "SpC", "SpI")),
    Tissue = factor(Tissue, levels = c("HK", "Spleen")),
    Condition = factor(Condition, levels = c("Control", "Infected"))
  )


## 2. Descriptive statistics for individual fish

rho_descriptive <- rho_analysis %>%
  summarise(
    number_of_samples = n(),
    minimum_rho = min(spearman_rho_recalculated),
    maximum_rho = max(spearman_rho_recalculated),
    mean_rho = mean(spearman_rho_recalculated),
    sd_rho = sd(spearman_rho_recalculated),
    median_rho = median(spearman_rho_recalculated),
    Q1_rho = quantile(spearman_rho_recalculated, 0.25),
    Q3_rho = quantile(spearman_rho_recalculated, 0.75)
  )

print(rho_descriptive)


## 3. Group-level descriptive statistics

rho_group_summary <- rho_analysis %>%
  group_by(Tissue, Condition) %>%
  summarise(
    number_of_samples = n(),
    mean_rho = mean(spearman_rho_recalculated),
    sd_rho = sd(spearman_rho_recalculated),
    median_rho = median(spearman_rho_recalculated),
    minimum_rho = min(spearman_rho_recalculated),
    maximum_rho = max(spearman_rho_recalculated),
    .groups = "drop"
  )

print(rho_group_summary)


## 4. Primary test: are fish-level correlations different from zero?

rho_zero_test <- wilcox.test(
  rho_analysis$spearman_rho_recalculated,
  mu = 0,
  alternative = "two.sided",
  exact = FALSE,
  conf.int = TRUE
)

print(rho_zero_test)

rho_zero_test_table <- data.frame(
  test = "One-sample Wilcoxon signed-rank test",
  number_of_samples = nrow(rho_analysis),
  median_rho = median(rho_analysis$spearman_rho_recalculated),
  statistic = unname(rho_zero_test$statistic),
  p_value = rho_zero_test$p.value
)

print(rho_zero_test_table)


## 5. Infection effect within each tissue

rho_condition_by_tissue <- rho_analysis %>%
  group_by(Tissue) %>%
  group_modify(
    ~ {
      test_result <- wilcox.test(
        spearman_rho_recalculated ~ Condition,
        data = .x,
        exact = FALSE
      )
      
      data.frame(
        control_median = median(
          .x$spearman_rho_recalculated[
            .x$Condition == "Control"
          ]
        ),
        infected_median = median(
          .x$spearman_rho_recalculated[
            .x$Condition == "Infected"
          ]
        ),
        statistic = unname(test_result$statistic),
        p_value = test_result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )

print(rho_condition_by_tissue)


## 6. Additive model for tissue and infection effects

# This evaluates whether the fish-level rho values differ by tissue or
# infection condition. Interpret cautiously because there are only 16 fish.

rho_additive_model <- lm(
  spearman_rho_recalculated ~ Tissue + Condition,
  data = rho_analysis
)

summary(rho_additive_model)
anova(rho_additive_model)

rho_model_coefficients <- broom::tidy(
  rho_additive_model,
  conf.int = TRUE
)

print(rho_model_coefficients)


## 7. Bootstrap confidence interval for the overall mean rho

set.seed(20260803)

N_BOOT <- 5000

rho_bootstrap <- replicate(
  N_BOOT,
  mean(
    sample(
      rho_analysis$spearman_rho_recalculated,
      size = nrow(rho_analysis),
      replace = TRUE
    )
  )
)

rho_bootstrap_summary <- data.frame(
  observed_mean_rho =
    mean(rho_analysis$spearman_rho_recalculated),
  
  bootstrap_lower_95 =
    unname(quantile(rho_bootstrap, 0.025)),
  
  bootstrap_upper_95 =
    unname(quantile(rho_bootstrap, 0.975))
)

print(rho_bootstrap_summary)


## 8. Prepare publication-ready figure

rho_plot <- ggplot(
  rho_analysis,
  aes(
    x = Condition,
    y = spearman_rho_recalculated,
    colour = Condition
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey40"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.7,
    alpha = 0.15
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 3,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ Tissue,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      "Control" = "#377EB8",
      "Infected" = "#E41A1C"
    )
  ) +
  scale_y_continuous(
    limits = c(-0.16, 0.06),
    breaks = seq(-0.15, 0.05, 0.05)
  ) +
  labs(
    title = "Within-fish relationship between clonotype abundance and SHM",
    x = NULL,
    y = "Spearman correlation (??)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0.5
    )
  )

print(rho_plot)


## 9. Save the figure

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      "Figure_sample_level_abundance_SHM_correlations_",
      ANALYSIS_SET,
      ".png"
    )
  ),
  plot = rho_plot,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    OUTPUT_DIR,
    paste0(
      "Figure_sample_level_abundance_SHM_correlations_",
      ANALYSIS_SET,
      ".pdf"
    )
  ),
  plot = rho_plot,
  width = 7,
  height = 4.5,
  units = "in",
  device = "pdf"
)











###############################################################################
# POOLED TOP-1,000 CLONOTYPE REGIONAL R/S ANALYSIS
#
# Input:  an R object named `shm_primary`
# Unit:   exact receptor clonotype/variant, pooled within Tissue x Condition
# Scope:  1,000 most abundant clonotypes in each pooled group
# QC:     essential sequence, alignment, length, and codon checks only
# Output: descriptive group summaries and figures (no sample-level objects/tests)
###############################################################################

rm(list = setdiff(ls(), "shm_primary"))

###############################################################################
# 1. Resolve columns once; edit only the candidate lists if necessary
###############################################################################

find_column <- function(data, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) > 0) return(hit[[1]])
  if (required) {
    stop(
      "Could not find the ", label, " column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

col_tissue <- find_column(
  shm_primary, c("Tissue", "tissue"), "tissue"
)
col_condition <- find_column(
  shm_primary, c("Condition", "condition"), "condition"
)
col_abundance <- find_column(
  shm_primary,
  c("Clones", "cloneCount", "Clone.count", "Count", "count", "reads"),
  "clonotype abundance"
)
col_v <- find_column(
  shm_primary,
  c("Germline_V_name", "V.name", "V_gene", "VGene", "allVHitsWithScore"),
  "V-gene"
)
col_j <- find_column(
  shm_primary,
  c("J.name", "J_gene", "JGene", "allJHitsWithScore"),
  "J-gene", required = FALSE
)
col_cdr3nt <- find_column(
  shm_primary,
  c("CDR3.nt", "CDR3_nt", "nSeqCDR3", "cdr3nt"),
  "CDR3 nucleotide", required = FALSE
)
col_cdr3aa <- find_column(
  shm_primary,
  c("CDR3.aa", "CDR3_aa", "aaSeqCDR3", "cdr3aa"),
  "CDR3 amino-acid", required = is.na(col_cdr3nt)
)
col_germline <- find_column(
  shm_primary,
  c("Germline_seq", "Germline.seq", "germline_seq", "germlineSequence"),
  "germline V sequence"
)

region_columns <- c(
  FR1  = find_column(shm_primary, c("FR1.nt",  "FR1_nt"),  "FR1"),
  CDR1 = find_column(shm_primary, c("CDR1.nt", "CDR1_nt"), "CDR1"),
  FR2  = find_column(shm_primary, c("FR2.nt",  "FR2_nt"),  "FR2"),
  CDR2 = find_column(shm_primary, c("CDR2.nt", "CDR2_nt"), "CDR2"),
  FR3  = find_column(shm_primary, c("FR3.nt",  "FR3_nt"),  "FR3")
)

###############################################################################
# 2. Standardize the required data
###############################################################################

clean_nt <- function(x) {
  x <- toupper(gsub("\\s+", "", as.character(x)))
  x[
    is.na(x) |
      x == "" |
      grepl("region_not_covered|not_covered|not covered", x, ignore.case = TRUE)
  ] <- NA_character_
  x
}

standardized <- tibble::tibble(
  Tissue = stringr::str_squish(as.character(shm_primary[[col_tissue]])),
  Condition = stringr::str_squish(as.character(shm_primary[[col_condition]])),
  Abundance = suppressWarnings(as.numeric(shm_primary[[col_abundance]])),
  V_gene = as.character(shm_primary[[col_v]]),
  J_gene = if (is.na(col_j)) NA_character_ else as.character(shm_primary[[col_j]]),
  CDR3_nt = if (is.na(col_cdr3nt)) NA_character_ else clean_nt(shm_primary[[col_cdr3nt]]),
  CDR3_aa = if (is.na(col_cdr3aa)) NA_character_ else as.character(shm_primary[[col_cdr3aa]]),
  Germline_seq = clean_nt(shm_primary[[col_germline]]),
  FR1_nt  = clean_nt(shm_primary[[region_columns[["FR1"]]]]),
  CDR1_nt = clean_nt(shm_primary[[region_columns[["CDR1"]]]]),
  FR2_nt  = clean_nt(shm_primary[[region_columns[["FR2"]]]]),
  CDR2_nt = clean_nt(shm_primary[[region_columns[["CDR2"]]]]),
  FR3_nt  = clean_nt(shm_primary[[region_columns[["FR3"]]]])
) %>%
  dplyr::mutate(
    Tissue = dplyr::case_when(
      grepl("head|kidney|^hk$", Tissue, ignore.case = TRUE) ~ "Head kidney",
      grepl("spleen|^sp$", Tissue, ignore.case = TRUE) ~ "Spleen",
      TRUE ~ Tissue
    ),
    Condition = dplyr::case_when(
      grepl("infect|challenge|^i$", Condition, ignore.case = TRUE) ~ "Infected",
      grepl("control|pbs|mock|^c$", Condition, ignore.case = TRUE) ~ "Control",
      TRUE ~ Condition
    ),
    Group = paste(Tissue, Condition, sep = " | ")
  )

if (anyNA(standardized$Abundance) || any(standardized$Abundance <= 0)) {
  stop("The abundance column contains missing, non-numeric, or non-positive values.")
}

if (anyNA(standardized$Tissue) || anyNA(standardized$Condition)) {
  stop("Tissue or Condition contains missing values.")
}
###############################################################################
# 3. Minimum pre-pooling QC
#
# Exact clonotype variants are defined by V/J, CDR3, germline reference, and
# recovered regional sequences. This avoids merging distinct SHM variants that
# happen to share the same CDR3 amino-acid sequence.
###############################################################################

identity_fields <- c(
  "V_gene", "J_gene",
  if (all(is.na(standardized$CDR3_nt))) "CDR3_aa" else "CDR3_nt",
  "Germline_seq", "FR1_nt", "CDR1_nt", "FR2_nt", "CDR2_nt", "FR3_nt"
)

TOP_N <- 1000L

pooled_clonotypes <- standardized %>%
  dplyr::filter(
    !is.na(V_gene), V_gene != "",
    !is.na(Germline_seq),
    grepl("^[ACGT]+$", Germline_seq)
  ) %>%
  dplyr::group_by(Group, Tissue, Condition, dplyr::across(dplyr::all_of(identity_fields))) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  dplyr::group_by(Group) %>%
  dplyr::arrange(dplyr::desc(Abundance), .by_group = TRUE) %>%
  dplyr::mutate(
    Rank = dplyr::row_number(),
    Group_total_abundance = sum(Abundance),
    Relative_abundance = Abundance / Group_total_abundance
  ) %>%
  dplyr::filter(Rank <= TOP_N) %>%
  dplyr::ungroup()

top1000_counts <- pooled_clonotypes %>%
  dplyr::count(Group, Tissue, Condition, name = "retained_clonotypes")

if (nrow(pooled_clonotypes) == 0) stop("No clonotypes remained after minimum QC.")
print(top1000_counts, n = Inf)

###############################################################################
# 4. Convert the selected clonotypes to region-level records
###############################################################################
regional_long <- pooled_clonotypes %>%
  dplyr::select(
    Group, Tissue, Condition, Rank, Abundance, Relative_abundance,
    V_gene, J_gene, CDR3_nt, Germline_seq,
    FR1_nt, CDR1_nt, FR2_nt, CDR2_nt, FR3_nt
  ) %>%
  tidyr::pivot_longer(
    cols = c(FR1_nt, CDR1_nt, FR2_nt, CDR2_nt, FR3_nt),
    names_to = "Region",
    values_to = "Observed_nt"
  ) %>%
  dplyr::mutate(
    Region = sub("_nt$", "", Region),
    Region = factor(
      Region,
      levels = c("FR1", "CDR1", "FR2", "CDR2", "FR3")
    )
  )


###############################################################################
# 5. Independently align each available observed region to its germline V gene
###############################################################################

align_one_region <- function(observed, germline) {
  empty <- tibble::tibble(
    Germline_region_nt = NA_character_,
    Alignment_identity = NA_real_,
    Alignment_status = NA_character_
  )
  
  if (is.na(observed)) {
    empty$Alignment_status <- "Region unavailable"
    return(empty)
  }
  if (is.na(germline)) {
    empty$Alignment_status <- "Germline unavailable"
    return(empty)
  }
  if (!grepl("^[ACGT]+$", observed) || !grepl("^[ACGT]+$", germline)) {
    empty$Alignment_status <- "Ambiguous nucleotide"
    return(empty)
  }
  if (nchar(observed) < 3L) {
    empty$Alignment_status <- "Region shorter than one codon"
    return(empty)
  }
  
  aln <- tryCatch(
    pwalign::pairwiseAlignment(
      pattern = Biostrings::DNAString(observed),
      subject = Biostrings::DNAString(germline),
      type = "local",
      substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
        match = 2, mismatch = -3, baseOnly = TRUE
      ),
      gapOpening = 8,
      gapExtension = 3
    ),
    error = function(e) NULL
  )
  
  if (is.null(aln)) {
    empty$Alignment_status <- "Alignment failed"
    return(empty)
  }
  
  obs_aln <- as.character(pwalign::alignedPattern(aln))
  germ_aln <- as.character(pwalign::alignedSubject(aln))
  
  if (grepl("-", obs_aln, fixed = TRUE) || grepl("-", germ_aln, fixed = TRUE)) {
    empty$Alignment_status <- "Alignment contains gap"
    return(empty)
  }
  if (nchar(obs_aln) != nchar(observed)) {
    empty$Alignment_status <- "Observed region incompletely aligned"
    return(empty)
  }
  if (nchar(germ_aln) != nchar(observed)) {
    empty$Alignment_status <- "Unequal aligned lengths"
    return(empty)
  }
  
  identity <- as.numeric(pwalign::pid(aln, type = "PID1"))
  if (is.na(identity) || identity < MIN_ALIGNMENT_IDENTITY) {
    empty$Alignment_identity <- identity
    empty$Alignment_status <- "Low alignment identity"
    return(empty)
  }
  
  tibble::tibble(
    Germline_region_nt = germ_aln,
    Alignment_identity = identity,
    Alignment_status = "Included"
  )
}

regional_aligned <- regional_long %>%
  dplyr::mutate(
    Alignment_result = purrr::map2(Observed_nt, Germline_seq, align_one_region)
  ) %>%
  tidyr::unnest(Alignment_result)

##
MIN_ALIGNMENT_IDENTITY <- 80
regional_aligned <- regional_long %>%
  dplyr::mutate(
    Alignment_result = purrr::map2(
      Observed_nt,
      Germline_seq,
      align_one_region
    )
  ) %>%
  tidyr::unnest(Alignment_result)

regional_classified <- regional_aligned %>%
  dplyr::mutate(
    RS_result = purrr::pmap(
      list(
        Observed_nt,
        Germline_region_nt,
        Alignment_status
      ),
      classify_rs
    )
  ) %>%
  tidyr::unnest(RS_result)

included_qc <- regional_classified %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Selected_top_clonotypes = dplyr::n(),
    Region_available = sum(!is.na(Observed_nt)),
    Included_clonotypes =
      sum(Classification_status == "Included", na.rm = TRUE),
    Included_percentage =
      100 * Included_clonotypes / Selected_top_clonotypes,
    Median_alignment_identity = if (
      Included_clonotypes > 0
    ) {
      median(
        Alignment_identity[
          Classification_status == "Included"
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

print(included_qc, n = Inf, width = Inf)

###############################################################################
# 6. Codon-aware replacement/silent classification
###############################################################################

classify_rs <- function(observed, germline, alignment_status) {
  empty <- tibble::tibble(
    Replacement_mutations = NA_integer_,
    Silent_mutations = NA_integer_,
    Total_mutations = NA_integer_,
    Classification_status = NA_character_
  )
  
  if (!identical(alignment_status, "Included")) {
    empty$Classification_status <- alignment_status
    return(empty)
  }
  if (nchar(observed) != nchar(germline)) {
    empty$Classification_status <- "Unequal sequence lengths"
    return(empty)
  }
  if (nchar(observed) %% 3L != 0L) {
    empty$Classification_status <- "Incomplete codon"
    return(empty)
  }
  
  starts <- seq.int(1L, nchar(observed), by = 3L)
  obs_codons <- substring(observed, starts, starts + 2L)
  germ_codons <- substring(germline, starts, starts + 2L)
  differences <- purrr::map2_int(
    obs_codons, germ_codons,
    ~ sum(strsplit(.x, "", fixed = TRUE)[[1]] != strsplit(.y, "", fixed = TRUE)[[1]])
  )
  
  obs_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(obs_codons)))
  germ_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(germ_codons)))
  mutated <- differences > 0L
  
  tibble::tibble(
    Replacement_mutations = sum(differences[mutated & obs_aa != germ_aa]),
    Silent_mutations = sum(differences[mutated & obs_aa == germ_aa]),
    Total_mutations = sum(differences),
    Classification_status = "Included"
  )
}

regional_classified <- regional_aligned %>%
  dplyr::mutate(
    RS_result = purrr::pmap(
      list(Observed_nt, Germline_region_nt, Alignment_status),
      classify_rs
    )
  ) %>%
  tidyr::unnest(RS_result)
###############################################################################
# 7. Minimal QC outputs
###############################################################################

qc_summary <- regional_classified %>%
  dplyr::count(Group, Tissue, Condition, Region, Classification_status, name = "Clonotypes") %>%
  dplyr::group_by(Group, Region) %>%
  dplyr::mutate(Percentage = 100 * Clonotypes / sum(Clonotypes)) %>%
  dplyr::ungroup()

included_qc <- regional_classified %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Selected_top_clonotypes = dplyr::n(),
    Region_available = sum(!is.na(Observed_nt)),
    Included_clonotypes = sum(Classification_status == "Included"),
    Included_percentage = 100 * Included_clonotypes / Selected_top_clonotypes,
    Median_alignment_identity = ifelse(
      Included_clonotypes > 0,
      median(Alignment_identity[Classification_status == "Included"], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

print(included_qc, n = Inf, width = Inf)

library(pwalign)

###############################################################################
# 8. Descriptive pooled-group analysis
#
# FR1 remains in QC/output but is labelled exploratory because its previously
# observed coverage was extremely low. No p-values are calculated: after
# pooling, the four groups are not biological replicates.
###############################################################################

pooled_rs_results <- regional_classified %>%
  dplyr::filter(Classification_status == "Included") %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Analyzed_clonotypes = dplyr::n(),
    Mutated_clonotypes = sum(Total_mutations > 0),
    Total_R = sum(Replacement_mutations),
    Total_S = sum(Silent_mutations),
    Total_mutations = sum(Total_mutations),
    R_fraction = dplyr::if_else(Total_R + Total_S > 0, Total_R / (Total_R + Total_S), NA_real_),
    R_S_ratio = dplyr::if_else(Total_S > 0, Total_R / Total_S, NA_real_),
    Mean_mutations_per_clonotype = mean(Total_mutations),
    Median_mutations_per_clonotype = median(Total_mutations),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Analysis_role = ifelse(Region == "FR1", "Exploratory", "Primary descriptive"))

abundance_weighted_results <- regional_classified %>%
  dplyr::filter(Classification_status == "Included") %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Weighted_R = sum(Replacement_mutations * Abundance),
    Weighted_S = sum(Silent_mutations * Abundance),
    Weighted_total_mutations = sum(Total_mutations * Abundance),
    Weighted_R_fraction = dplyr::if_else(
      Weighted_R + Weighted_S > 0,
      Weighted_R / (Weighted_R + Weighted_S),
      NA_real_
    ),
    Weighted_R_S_ratio = dplyr::if_else(Weighted_S > 0, Weighted_R / Weighted_S, NA_real_),
    .groups = "drop"
  )

###############################################################################
# 8. Descriptive pooled-group analysis
#
# FR1 remains in QC/output but is labelled exploratory because its previously
# observed coverage was extremely low. No p-values are calculated: after
# pooling, the four groups are not biological replicates.
###############################################################################

pooled_rs_results <- regional_classified %>%
  dplyr::filter(Classification_status == "Included") %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Analyzed_clonotypes = dplyr::n(),
    Mutated_clonotypes = sum(Total_mutations > 0),
    Total_R = sum(Replacement_mutations),
    Total_S = sum(Silent_mutations),
    Total_mutations = sum(Total_mutations),
    R_fraction = dplyr::if_else(Total_R + Total_S > 0, Total_R / (Total_R + Total_S), NA_real_),
    R_S_ratio = dplyr::if_else(Total_S > 0, Total_R / Total_S, NA_real_),
    Mean_mutations_per_clonotype = mean(Total_mutations),
    Median_mutations_per_clonotype = median(Total_mutations),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Analysis_role = ifelse(Region == "FR1", "Exploratory", "Primary descriptive"))

abundance_weighted_results <- regional_classified %>%
  dplyr::filter(Classification_status == "Included") %>%
  dplyr::group_by(Group, Tissue, Condition, Region) %>%
  dplyr::summarise(
    Weighted_R = sum(Replacement_mutations * Abundance),
    Weighted_S = sum(Silent_mutations * Abundance),
    Weighted_total_mutations = sum(Total_mutations * Abundance),
    Weighted_R_fraction = dplyr::if_else(
      Weighted_R + Weighted_S > 0,
      Weighted_R / (Weighted_R + Weighted_S),
      NA_real_
    ),
    Weighted_R_S_ratio = dplyr::if_else(Weighted_S > 0, Weighted_R / Weighted_S, NA_real_),
    .groups = "drop"
  )
###############################################################################
# 9. Publication-style descriptive figures
###############################################################################

group_colors <- c(
  "Head kidney | Control" = "#4C78A8",
  "Head kidney | Infected" = "#E45756",
  "Spleen | Control" = "#72B7B2",
  "Spleen | Infected" = "#F58518"
)

plot_theme <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )
#######################
figure_rs_ratio <- pooled_rs_results %>%
  dplyr::filter(Region != "FR1") %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Region,
      y = R_S_ratio,
      fill = Group
    )
  ) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.72
  ) +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(
    title = "Replacement-to-silent mutation ratio across IGHV regions",
    x = NULL,
    y = "Replacement/silent ratio",
    fill = "Group"
  ) +
  plot_theme +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    )
  )

print(figure_rs_ratio)

###########################
figure_r_fraction <- pooled_rs_results %>%
  dplyr::filter(Region != "FR1") %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Region,
      y = R_fraction,
      fill = Group
    )
  ) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.72
  ) +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Replacement mutation fraction across IGHV regions",
    x = NULL,
    y = "Replacement fraction",
    fill = "Group"
  ) +
  plot_theme

print(figure_r_fraction)

##############
figure_coverage <- included_qc %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Region,
      y = Included_percentage,
      fill = Group
    )
  ) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.72
  ) +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(
    title = "Regional sequence coverage after quality control",
    x = NULL,
    y = "Top clonotypes passing regional QC (%)",
    fill = "Group"
  ) +
  plot_theme +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    )
  )

print(figure_coverage)

###############################################################################
# CLONAL EXPANSION AND V-REGION SOMATIC MUTATION
#
# Dataset:
#   shm_primary
#
# Primary question:
#   Does V-region mutation frequency increase with clonotype abundance, and
#   does this relationship differ between control and SAV-infected fish?
#
# Important:
#   - Uses sample-level clonotypes, not pooled_clonotypes.
#   - Recomputes nucleotide substitutions after sequence alignment.
#   - Insertions, deletions, and unobserved germline positions are not counted
#     as nucleotide substitutions.
#   - Primary minimum comparable alignment length = 90 nt.
#   - Sensitivity thresholds = 75, 100, and 150 nt.
###############################################################################


###############################################################################
# 1. Install and load packages
###############################################################################


###############################################################################
# 2. Confirm that the established input exists
###############################################################################

if (!exists("shm_primary")) {
  stop(
    paste0(
      "The object 'shm_primary' is not present in the current R session.\n",
      "Run the completed SHM preprocessing workflow that creates shm_primary, ",
      "then rerun this script."
    )
  )
}

required_columns <- c(
  "Clones",
  "CDR3.aa",
  "V.name",
  "Sample_clean",
  "Group",
  "Tissue",
  "Condition",
  "V_region_nt",
  "Germline_seq",
  "Germline_V_name",
  "number_of_best_matches",
  "germline_match_status"
)

missing_columns <- setdiff(
  required_columns,
  names(shm_primary)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Required columns missing from shm_primary:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


###############################################################################
# 3. Analysis settings
###############################################################################

OUTPUT_DIR <- "Expansion_SHM_final_analysis"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}

PRIMARY_MIN_ALIGNED_NT <- 90
SENSITIVITY_THRESHOLDS <- c(75, 90, 100, 150)

MIN_CLONE_ABUNDANCE <- 2
MIN_ALIGNMENT_IDENTITY <- 80
MIN_ALIGNMENT_COVERAGE <- 0.70

# Final balanced statistical dataset used in the manuscript.
FINAL_SAMPLES <- c(
  "HKC_3C3w",
  "HKC_4C3w",
  "HKC_5C3w",
  "HKI_1I3w",
  "HKI_4I3w",
  "HKI_5I3w",
  "SpC_1C3w",
  "SpC_3C3w",
  "SpC_5C3w",
  "SpI_1I3w",
  "SpI_2I3w",
  "SpI_5I3w"
)


###############################################################################
# 4. Restrict analysis to the final 12 samples
###############################################################################

shm_expansion_input <- shm_primary %>%
  mutate(
    Sample = as.character(Sample_clean),
    Tissue = as.character(Tissue),
    Condition = as.character(Condition),
    Clone_abundance = suppressWarnings(
      as.numeric(Clones)
    ),
    CDR3_aa = as.character(CDR3.aa),
    V_gene = as.character(V.name),
    Germline_V_gene = as.character(Germline_V_name),
    Observed_V = toupper(
      gsub(
        "[[:space:]]",
        "",
        as.character(V_region_nt)
      )
    ),
    Germline_V = toupper(
      gsub(
        "[[:space:]]",
        "",
        as.character(Germline_seq)
      )
    )
  ) %>%
  filter(Sample %in% FINAL_SAMPLES)

missing_final_samples <- setdiff(
  FINAL_SAMPLES,
  unique(shm_expansion_input$Sample)
)

if (length(missing_final_samples) > 0) {
  stop(
    paste0(
      "These final samples were not found in shm_primary:\n",
      paste(missing_final_samples, collapse = ", ")
    )
  )
}

sample_design <- shm_expansion_input %>%
  distinct(
    Sample,
    Tissue,
    Condition
  ) %>%
  arrange(
    Tissue,
    Condition,
    Sample
  )

if (nrow(sample_design) != 12) {
  stop(
    paste(
      "Expected 12 biological samples but found",
      nrow(sample_design)
    )
  )
}

cat("\n===== FINAL SAMPLE DESIGN =====\n")
print(
  sample_design,
  n = Inf,
  width = Inf
)


###############################################################################
# 5. Pre-alignment QC
###############################################################################

shm_pre_alignment <- shm_expansion_input %>%
  mutate(
    Observed_length = nchar(Observed_V),
    Germline_length = nchar(Germline_V),
    
    Pre_alignment_QC = case_when(
      is.na(Sample) |
        Sample == "" ~
        "Missing sample identity",
      
      is.na(Clone_abundance) ~
        "Missing clone abundance",
      
      Clone_abundance < MIN_CLONE_ABUNDANCE ~
        "Clone abundance below 2",
      
      is.na(CDR3_aa) |
        CDR3_aa == "" ~
        "Missing CDR3 amino-acid sequence",
      
      grepl(
        "\\*|stop|region_not_covered",
        CDR3_aa,
        ignore.case = TRUE
      ) ~
        "Invalid or nonproductive CDR3",
      
      is.na(V_gene) |
        V_gene == "" ~
        "Missing IGHV assignment",
      
      is.na(Germline_V_gene) |
        Germline_V_gene == "" ~
        "Missing germline IGHV assignment",
      
      is.na(Observed_V) |
        Observed_V == "" ~
        "Missing observed V-region sequence",
      
      is.na(Germline_V) |
        Germline_V == "" ~
        "Missing germline V-region sequence",
      
      !grepl(
        "^[ACGT]+$",
        Observed_V
      ) ~
        "Ambiguous observed nucleotide",
      
      !grepl(
        "^[ACGT]+$",
        Germline_V
      ) ~
        "Ambiguous germline nucleotide",
      
      is.na(number_of_best_matches) ~
        "Missing germline-match information",
      
      number_of_best_matches != 1 ~
        "Ambiguous best germline match",
      
      TRUE ~
        "Eligible for alignment"
    )
  )

pre_alignment_qc_summary <- shm_pre_alignment %>%
  count(
    Pre_alignment_QC,
    name = "Clonotypes"
  ) %>%
  mutate(
    Percentage = 100 * Clonotypes / sum(Clonotypes)
  ) %>%
  arrange(desc(Clonotypes))

cat("\n===== PRE-ALIGNMENT QC =====\n")
print(
  pre_alignment_qc_summary,
  n = Inf,
  width = Inf
)

alignment_input <- shm_pre_alignment %>%
  filter(
    Pre_alignment_QC == "Eligible for alignment"
  )

if (nrow(alignment_input) == 0) {
  stop("No clonotypes remained after pre-alignment QC.")
}


###############################################################################
# 6. Align observed V-region fragments to assigned germline sequences
###############################################################################

# Local alignment is necessary because observed V-region fragments frequently
# cover only part of the full germline V sequence.
#
# The scoring scheme discourages mismatches and gaps while permitting genuine
# partial-coverage alignment.
#
# Mutation counts are subsequently calculated only where both aligned
# sequences contain A/C/G/T. Alignment gaps are excluded.

ALIGNMENT_CACHE <- file.path(
  OUTPUT_DIR,
  "aligned_V_region_clonotypes.rds"
)

if (file.exists(ALIGNMENT_CACHE)) {
  
  cat("\nLoading previously calculated V-region alignments...\n")
  
  shm_aligned <- readRDS(
    ALIGNMENT_CACHE
  )
  
} else {
  
  cat(
    "\nAligning ",
    nrow(alignment_input),
    " observed V-region sequences to their assigned germlines...\n",
    sep = ""
  )
  
  alignment_result <- pwalign::pairwiseAlignment(
    pattern = Biostrings::DNAStringSet(
      alignment_input$Observed_V
    ),
    subject = Biostrings::DNAStringSet(
      alignment_input$Germline_V
    ),
    type = "local",
    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
      match = 2,
      mismatch = -3,
      baseOnly = TRUE
    ),
    gapOpening = 8,
    gapExtension = 2,
    scoreOnly = FALSE
  )
  
  aligned_observed <- as.character(
    pwalign::alignedPattern(
      alignment_result
    )
  )
  
  aligned_germline <- as.character(
    pwalign::alignedSubject(
      alignment_result
    )
  )
  
  calculate_alignment_metrics <- function(
    observed_alignment,
    germline_alignment,
    observed_full_length,
    germline_full_length
  ) {
    
    observed_vector <- strsplit(
      observed_alignment,
      "",
      fixed = TRUE
    )[[1]]
    
    germline_vector <- strsplit(
      germline_alignment,
      "",
      fixed = TRUE
    )[[1]]
    
    if (length(observed_vector) != length(germline_vector)) {
      return(
        tibble(
          Comparable_nt = NA_integer_,
          Substitutions = NA_integer_,
          Matches = NA_integer_,
          Alignment_identity_percent = NA_real_,
          Observed_coverage = NA_real_,
          Shorter_sequence_coverage = NA_real_
        )
      )
    }
    
    comparable <- observed_vector %in% c("A", "C", "G", "T") &
      germline_vector %in% c("A", "C", "G", "T")
    
    comparable_nt <- sum(comparable)
    
    if (comparable_nt == 0) {
      return(
        tibble(
          Comparable_nt = 0L,
          Substitutions = NA_integer_,
          Matches = NA_integer_,
          Alignment_identity_percent = NA_real_,
          Observed_coverage = 0,
          Shorter_sequence_coverage = 0
        )
      )
    }
    
    substitutions <- sum(
      observed_vector[comparable] !=
        germline_vector[comparable]
    )
    
    matches <- comparable_nt - substitutions
    
    tibble(
      Comparable_nt = as.integer(comparable_nt),
      Substitutions = as.integer(substitutions),
      Matches = as.integer(matches),
      
      Alignment_identity_percent =
        100 * matches / comparable_nt,
      
      Observed_coverage =
        comparable_nt / observed_full_length,
      
      Shorter_sequence_coverage =
        comparable_nt /
        min(
          observed_full_length,
          germline_full_length
        )
    )
  }
  
  alignment_metrics <- purrr::pmap_dfr(
    list(
      aligned_observed,
      aligned_germline,
      alignment_input$Observed_length,
      alignment_input$Germline_length
    ),
    calculate_alignment_metrics
  )
  
  shm_aligned <- bind_cols(
    alignment_input,
    tibble(
      Aligned_observed = aligned_observed,
      Aligned_germline = aligned_germline,
      Alignment_score = as.numeric(
        pwalign::score(
          alignment_result
        )
      )
    ),
    alignment_metrics
  )
  
  saveRDS(
    shm_aligned,
    ALIGNMENT_CACHE
  )
}

cat("\nV-region alignment completed.\n")


###############################################################################
# 7. Alignment QC and mutation-frequency calculation
###############################################################################

shm_alignment_qc <- shm_aligned %>%
  mutate(
    Alignment_QC = case_when(
      is.na(Comparable_nt) |
        Comparable_nt == 0 ~
        "No comparable aligned nucleotides",
      
      Shorter_sequence_coverage <
        MIN_ALIGNMENT_COVERAGE ~
        "Insufficient alignment coverage",
      
      Alignment_identity_percent <
        MIN_ALIGNMENT_IDENTITY ~
        "Alignment identity below 80%",
      
      is.na(Substitutions) |
        Substitutions < 0 ~
        "Invalid substitution count",
      
      Substitutions > Comparable_nt ~
        "Substitutions exceed comparable positions",
      
      TRUE ~
        "Pass"
    ),
    
    Mutation_frequency =
      Substitutions / Comparable_nt,
    
    Mutation_frequency_percent =
      100 * Mutation_frequency,
    
    Log10_abundance =
      log10(Clone_abundance)
  )

alignment_qc_summary <- shm_alignment_qc %>%
  count(
    Alignment_QC,
    name = "Clonotypes"
  ) %>%
  mutate(
    Percentage = 100 * Clonotypes / sum(Clonotypes)
  ) %>%
  arrange(desc(Clonotypes))

cat("\n===== ALIGNMENT QC =====\n")
print(
  alignment_qc_summary,
  n = Inf,
  width = Inf
)


###############################################################################
# 8. Function for constructing threshold-specific datasets
###############################################################################

make_analysis_dataset <- function(min_aligned_nt) {
  
  result <- shm_alignment_qc %>%
    filter(
      Alignment_QC == "Pass",
      Comparable_nt >= min_aligned_nt,
      Clone_abundance >= MIN_CLONE_ABUNDANCE
    ) %>%
    mutate(
      Minimum_alignment_threshold =
        min_aligned_nt,
      
      Expansion_class = case_when(
        Clone_abundance >= 2 &
          Clone_abundance <= 3 ~
          "Small (2???3)",
        
        Clone_abundance >= 4 &
          Clone_abundance <= 10 ~
          "Moderate (4???10)",
        
        Clone_abundance >= 11 &
          Clone_abundance <= 30 ~
          "Expanded (11???30)",
        
        Clone_abundance > 30 ~
          "Highly expanded (>30)"
      ),
      
      Expansion_class = factor(
        Expansion_class,
        levels = c(
          "Small (2???3)",
          "Moderate (4???10)",
          "Expanded (11???30)",
          "Highly expanded (>30)"
        )
      ),
      
      Tissue = case_when(
        Tissue %in% c(
          "HK",
          "Head kidney",
          "Head Kidney"
        ) ~ "Head kidney",
        
        Tissue %in% c(
          "Sp",
          "SP",
          "Spleen"
        ) ~ "Spleen",
        
        TRUE ~ Tissue
      ),
      
      Tissue = factor(
        Tissue,
        levels = c(
          "Head kidney",
          "Spleen"
        )
      ),
      
      Condition = factor(
        Condition,
        levels = c(
          "Control",
          "Infected"
        )
      )
    ) %>%
    filter(
      !is.na(Tissue),
      !is.na(Condition),
      !is.na(Expansion_class),
      is.finite(Mutation_frequency),
      Mutation_frequency >= 0,
      Mutation_frequency <= 1
    )
  
  if (nrow(result) > 0) {
    result <- result %>%
      mutate(
        Centered_log10_abundance =
          Log10_abundance -
          mean(
            Log10_abundance,
            na.rm = TRUE
          )
      )
  }
  
  result
}


###############################################################################
# 9. Primary analysis dataset: minimum 90 comparable nucleotides
###############################################################################

shm_analysis <- make_analysis_dataset(
  PRIMARY_MIN_ALIGNED_NT
)

if (nrow(shm_analysis) == 0) {
  stop(
    paste(
      "No clonotypes passed the primary",
      PRIMARY_MIN_ALIGNED_NT,
      "nt threshold."
    )
  )
}

primary_sample_summary <- shm_analysis %>%
  group_by(
    Sample,
    Tissue,
    Condition
  ) %>%
  summarise(
    Clonotypes = n(),
    Median_comparable_nt = median(
      Comparable_nt
    ),
    Median_clone_abundance = median(
      Clone_abundance
    ),
    Median_mutation_frequency_percent = median(
      Mutation_frequency_percent
    ),
    .groups = "drop"
  )

cat("\n===== PRIMARY DATASET: SAMPLE SUMMARY =====\n")
print(
  primary_sample_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nPrimary clonotypes retained: ",
  nrow(shm_analysis),
  "\nBiological samples retained: ",
  n_distinct(shm_analysis$Sample),
  "\n",
  sep = ""
)


###############################################################################
# 10. Descriptive results by expansion class
###############################################################################

expansion_class_summary <- shm_analysis %>%
  group_by(
    Tissue,
    Condition,
    Expansion_class
  ) %>%
  summarise(
    Biological_samples = n_distinct(Sample),
    Clonotypes = n(),
    
    Median_clone_abundance = median(
      Clone_abundance
    ),
    
    Median_mutation_frequency_percent = median(
      Mutation_frequency_percent
    ),
    
    Q1_mutation_frequency_percent = quantile(
      Mutation_frequency_percent,
      0.25
    ),
    
    Q3_mutation_frequency_percent = quantile(
      Mutation_frequency_percent,
      0.75
    ),
    
    .groups = "drop"
  )

cat("\n===== EXPANSION-CLASS SUMMARY =====\n")
print(
  expansion_class_summary,
  n = Inf,
  width = Inf
)


###############################################################################
# 11. Per-sample Spearman correlations
###############################################################################

sample_correlations <- shm_analysis %>%
  group_by(
    Sample,
    Tissue,
    Condition
  ) %>%
  summarise(
    Clonotypes = n(),
    Distinct_abundances = n_distinct(
      Clone_abundance
    ),
    Distinct_mutation_frequencies = n_distinct(
      Mutation_frequency
    ),
    
    Spearman_rho = {
      valid_test <- n() >= 10 &&
        n_distinct(Clone_abundance) >= 2 &&
        n_distinct(Mutation_frequency) >= 2
      
      if (valid_test) {
        suppressWarnings(
          cor(
            Log10_abundance,
            Mutation_frequency,
            method = "spearman",
            use = "complete.obs"
          )
        )
      } else {
        NA_real_
      }
    },
    
    Spearman_p = {
      valid_test <- n() >= 10 &&
        n_distinct(Clone_abundance) >= 2 &&
        n_distinct(Mutation_frequency) >= 2
      
      if (valid_test) {
        suppressWarnings(
          cor.test(
            Log10_abundance,
            Mutation_frequency,
            method = "spearman",
            exact = FALSE
          )$p.value
        )
      } else {
        NA_real_
      }
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    BH_adjusted_p = p.adjust(
      Spearman_p,
      method = "BH"
    )
  )

cat("\n===== PER-SAMPLE SPEARMAN CORRELATIONS =====\n")
print(
  sample_correlations,
  n = Inf,
  width = Inf
)

sample_correlation_group_summary <- sample_correlations %>%
  group_by(
    Tissue,
    Condition
  ) %>%
  summarise(
    Samples = sum(!is.na(Spearman_rho)),
    Median_rho = median(
      Spearman_rho,
      na.rm = TRUE
    ),
    Minimum_rho = min(
      Spearman_rho,
      na.rm = TRUE
    ),
    Maximum_rho = max(
      Spearman_rho,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


###############################################################################
# 12. Primary beta-binomial mixed model
###############################################################################

# This model uses the underlying nucleotide counts rather than treating
# mutation percentages as continuous normally distributed observations.
#
# Response:
#   substitutions versus matched comparable nucleotides
#
# Random intercept:
#   biological sample
#
# Interaction:
#   tests whether the abundance???mutation relationship differs between
#   control and infected fish.

continuous_interaction_model <- glmmTMB(
  cbind(
    Substitutions,
    Matches
  ) ~
    Centered_log10_abundance * Condition +
    Tissue +
    (1 | Sample),
  family = betabinomial(
    link = "logit"
  ),
  data = shm_analysis
)

continuous_additive_model <- glmmTMB(
  cbind(
    Substitutions,
    Matches
  ) ~
    Centered_log10_abundance +
    Condition +
    Tissue +
    (1 | Sample),
  family = betabinomial(
    link = "logit"
  ),
  data = shm_analysis
)

cat("\n===== CONTINUOUS INTERACTION MODEL =====\n")
print(
  summary(
    continuous_interaction_model
  )
)

cat("\n===== ADDITIVE MODEL =====\n")
print(
  summary(
    continuous_additive_model
  )
)

continuous_model_comparison <- anova(
  continuous_additive_model,
  continuous_interaction_model
)

cat("\n===== ADDITIVE VERSUS INTERACTION MODEL =====\n")
print(
  continuous_model_comparison
)

interaction_p <- continuous_model_comparison$`Pr(>Chisq)`[2]

selected_continuous_model <- if (
  !is.na(interaction_p) &&
  interaction_p < 0.05
) {
  continuous_interaction_model
} else {
  continuous_additive_model
}

selected_model_name <- if (
  identical(
    selected_continuous_model,
    continuous_interaction_model
  )
) {
  "Interaction model"
} else {
  "Additive model"
}

cat(
  "\nSelected primary model: ",
  selected_model_name,
  "\n",
  sep = ""
)

continuous_model_coefficients <- broom.mixed::tidy(
  selected_continuous_model,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

cat("\n===== SELECTED MODEL ODDS RATIOS =====\n")
print(
  continuous_model_coefficients,
  n = Inf,
  width = Inf
)


###############################################################################
# 13. Estimated abundance slopes by condition
###############################################################################

abundance_slopes <- emmeans::emtrends(
  continuous_interaction_model,
  specs = ~ Condition,
  var = "Centered_log10_abundance"
)

abundance_slope_contrast <- pairs(
  abundance_slopes,
  adjust = "none"
)

cat("\n===== ABUNDANCE SLOPES BY CONDITION =====\n")
print(abundance_slopes)

cat("\n===== DIFFERENCE BETWEEN CONDITION-SPECIFIC SLOPES =====\n")
print(abundance_slope_contrast)

abundance_slopes_table <- as.data.frame(
  abundance_slopes
)

abundance_slope_contrast_table <- as.data.frame(
  abundance_slope_contrast
)


###############################################################################
# 14. Expansion-class beta-binomial mixed model
###############################################################################

class_interaction_model <- glmmTMB(
  cbind(
    Substitutions,
    Matches
  ) ~
    Expansion_class * Condition +
    Tissue +
    (1 | Sample),
  family = betabinomial(
    link = "logit"
  ),
  data = shm_analysis
)

class_additive_model <- glmmTMB(
  cbind(
    Substitutions,
    Matches
  ) ~
    Expansion_class +
    Condition +
    Tissue +
    (1 | Sample),
  family = betabinomial(
    link = "logit"
  ),
  data = shm_analysis
)

class_model_comparison <- anova(
  class_additive_model,
  class_interaction_model
)

cat("\n===== EXPANSION-CLASS MODEL COMPARISON =====\n")
print(
  class_model_comparison
)

class_interaction_p <- class_model_comparison$`Pr(>Chisq)`[2]

selected_class_model <- if (
  !is.na(class_interaction_p) &&
  class_interaction_p < 0.05
) {
  class_interaction_model
} else {
  class_additive_model
}

class_emmeans <- emmeans(
  selected_class_model,
  ~ Expansion_class | Condition,
  type = "response"
)

class_pairwise <- pairs(
  class_emmeans,
  adjust = "BH"
)

cat("\n===== EXPANSION-CLASS ESTIMATED MUTATION PROPORTIONS =====\n")
print(class_emmeans)

cat("\n===== EXPANSION-CLASS PAIRWISE COMPARISONS =====\n")
print(class_pairwise)

class_emmeans_table <- as.data.frame(
  class_emmeans
)

class_pairwise_table <- as.data.frame(
  class_pairwise
)


###############################################################################
# 15. Sample-level expansion-class medians
###############################################################################

sample_class_summary <- shm_analysis %>%
  group_by(
    Sample,
    Tissue,
    Condition,
    Expansion_class
  ) %>%
  summarise(
    Clonotypes = n(),
    Median_mutation_frequency_percent = median(
      Mutation_frequency_percent
    ),
    .groups = "drop"
  )


###############################################################################
# 16. Sensitivity analyses at 75, 90, 100, and 150 aligned nucleotides
###############################################################################

fit_sensitivity_model <- function(minimum_nt) {
  
  sensitivity_data <- make_analysis_dataset(
    minimum_nt
  )
  
  represented_samples <- n_distinct(
    sensitivity_data$Sample
  )
  
  if (
    nrow(sensitivity_data) < 100 ||
    represented_samples < 8
  ) {
    return(
      tibble(
        Minimum_aligned_nt = minimum_nt,
        Clonotypes = nrow(sensitivity_data),
        Biological_samples = represented_samples,
        Model_converged = FALSE,
        Term = NA_character_,
        Estimate_log_odds = NA_real_,
        Odds_ratio = NA_real_,
        Standard_error = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  model <- tryCatch(
    glmmTMB(
      cbind(
        Substitutions,
        Matches
      ) ~
        Centered_log10_abundance +
        Condition +
        Tissue +
        (1 | Sample),
      family = betabinomial(
        link = "logit"
      ),
      data = sensitivity_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    return(
      tibble(
        Minimum_aligned_nt = minimum_nt,
        Clonotypes = nrow(sensitivity_data),
        Biological_samples = represented_samples,
        Model_converged = FALSE,
        Term = NA_character_,
        Estimate_log_odds = NA_real_,
        Odds_ratio = NA_real_,
        Standard_error = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  model_results <- broom.mixed::tidy(
    model,
    effects = "fixed"
  ) %>%
    filter(
      term == "Centered_log10_abundance"
    )
  
  tibble(
    Minimum_aligned_nt = minimum_nt,
    Clonotypes = nrow(sensitivity_data),
    Biological_samples = represented_samples,
    Model_converged = isTRUE(
      model$sdr$pdHess
    ),
    Term = model_results$term,
    Estimate_log_odds = model_results$estimate,
    Odds_ratio = exp(
      model_results$estimate
    ),
    Standard_error = model_results$std.error,
    P_value = model_results$p.value
  )
}

sensitivity_results <- purrr::map_dfr(
  SENSITIVITY_THRESHOLDS,
  fit_sensitivity_model
)

cat("\n===== ALIGNMENT-LENGTH SENSITIVITY ANALYSES =====\n")
print(
  sensitivity_results,
  n = Inf,
  width = Inf
)


###############################################################################
# 17. Plot settings
###############################################################################

condition_colors <- c(
  "Control" = "#4D4D4D",
  "Infected" = "#D55E00"
)

expansion_colors <- c(
  "Small (2???3)" = "#D9D9D9",
  "Moderate (4???10)" = "#92C5DE",
  "Expanded (11???30)" = "#4393C3",
  "Highly expanded (>30)" = "#2166AC"
)

publication_theme <- theme_classic(
  base_size = 13
) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "top",
    legend.title = element_text(
      face = "bold"
    )
  )


###############################################################################
# 18. Figure A: continuous abundance???mutation association
###############################################################################

figure_continuous <- ggplot(
  shm_analysis,
  aes(
    x = Log10_abundance,
    y = Mutation_frequency_percent
  )
) +
  geom_point(
    aes(color = Condition),
    alpha = 0.10,
    size = 0.65,
    position = position_jitter(
      width = 0.012,
      height = 0
    )
  ) +
  
  # Thin within-sample trends
  geom_smooth(
    aes(
      group = Sample,
      color = Condition
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.35,
    alpha = 0.35
  ) +
  
  # Thick condition-level descriptive trends
  geom_smooth(
    aes(
      color = Condition
    ),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.15,
    alpha = 0.14
  ) +
  facet_wrap(
    ~ Tissue,
    nrow = 1
  ) +
  scale_color_manual(
    values = condition_colors
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(0.02, 0.07)
    )
  ) +
  labs(
    title = "Association between clonal expansion and V-region mutation frequency",
    subtitle = paste0(
      "Primary analysis: ???",
      PRIMARY_MIN_ALIGNED_NT,
      " comparable V-region nucleotides"
    ),
    x = expression(
      log[10] * "(clone abundance)"
    ),
    y = "V-region mutation frequency (%)",
    color = "Condition"
  ) +
  publication_theme

print(figure_continuous)


###############################################################################
# 19. Figure B: expansion-class analysis using sample-level medians
###############################################################################

figure_expansion_classes <- ggplot(
  sample_class_summary,
  aes(
    x = Expansion_class,
    y = Median_mutation_frequency_percent,
    color = Condition
  )
) +
  geom_boxplot(
    aes(group = interaction(
      Expansion_class,
      Condition
    )),
    width = 0.58,
    outlier.shape = NA,
    color = "black",
    fill = NA,
    linewidth = 0.55
  ) +
  geom_point(
    aes(
      shape = Condition
    ),
    size = 2.7,
    position = position_jitter(
      width = 0.10,
      height = 0
    )
  ) +
  facet_wrap(
    ~ Tissue,
    nrow = 1
  ) +
  scale_color_manual(
    values = condition_colors
  ) +
  scale_shape_manual(
    values = c(
      "Control" = 16,
      "Infected" = 17
    )
  ) +
  labs(
    title = "V-region mutation frequency across clonotype expansion classes",
    subtitle = "Each point represents the median for one biological sample",
    x = "Expansion class",
    y = "Median V-region mutation frequency (%)",
    color = "Condition",
    shape = "Condition"
  ) +
  publication_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

print(figure_expansion_classes)


###############################################################################
# 20. Figure C: within-sample Spearman correlations
###############################################################################

figure_correlations <- sample_correlations %>%
  filter(
    !is.na(Spearman_rho)
  ) %>%
  ggplot(
    aes(
      x = Condition,
      y = Spearman_rho,
      color = Condition
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    color = "black",
    fill = NA,
    linewidth = 0.55
  ) +
  geom_point(
    size = 3,
    position = position_jitter(
      width = 0.07,
      height = 0
    )
  ) +
  facet_wrap(
    ~ Tissue,
    nrow = 1
  ) +
  scale_color_manual(
    values = condition_colors
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(
      -1,
      1,
      0.25
    )
  ) +
  labs(
    title = "Within-sample abundance???mutation associations",
    subtitle = "Each point represents one biological sample",
    x = NULL,
    y = "Spearman correlation coefficient"
  ) +
  publication_theme +
  theme(
    legend.position = "none"
  )

print(figure_correlations)


###############################################################################
# 21. Figure D: sensitivity to minimum aligned length
###############################################################################

figure_sensitivity <- sensitivity_results %>%
  filter(
    !is.na(Odds_ratio)
  ) %>%
  ggplot(
    aes(
      x = factor(
        Minimum_aligned_nt,
        levels = SENSITIVITY_THRESHOLDS
      ),
      y = Odds_ratio
    )
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_line(
    aes(group = 1),
    linewidth = 0.7,
    color = "#2166AC"
  ) +
  geom_point(
    aes(
      size = Clonotypes,
      fill = Model_converged
    ),
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#56B4E9",
      "FALSE" = "#E69F00"
    )
  ) +
  labs(
    title = "Sensitivity to minimum V-region alignment length",
    x = "Minimum comparable aligned nucleotides",
    y = expression(
      "Abundance effect (odds ratio per 1-unit increase in " *
        log[10] * " abundance)"
    ),
    size = "Clonotypes",
    fill = "Converged"
  ) +
  publication_theme

print(figure_sensitivity)


###############################################################################
# 22. Combined publication figure
###############################################################################

combined_figure <- (
  figure_continuous /
    (
      figure_expansion_classes |
        figure_correlations
    )
) +
  plot_annotation(
    title = "Clonal expansion and V-region somatic mutation",
    tag_levels = "A"
  )

print(combined_figure)


###############################################################################
# 23. Save figures
###############################################################################

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_continuous_expansion_vs_SHM.pdf"
  ),
  figure_continuous,
  width = 9.5,
  height = 5.5,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_continuous_expansion_vs_SHM.png"
  ),
  figure_continuous,
  width = 9.5,
  height = 5.5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_expansion_classes_vs_SHM.pdf"
  ),
  figure_expansion_classes,
  width = 10,
  height = 5.8,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_expansion_classes_vs_SHM.png"
  ),
  figure_expansion_classes,
  width = 10,
  height = 5.8,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_sample_Spearman_correlations.pdf"
  ),
  figure_correlations,
  width = 8,
  height = 5.5,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_sample_Spearman_correlations.png"
  ),
  figure_correlations,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_alignment_threshold_sensitivity.pdf"
  ),
  figure_sensitivity,
  width = 8,
  height = 5.5,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_alignment_threshold_sensitivity.png"
  ),
  figure_sensitivity,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_combined_expansion_SHM.pdf"
  ),
  combined_figure,
  width = 12,
  height = 12,
  units = "in"
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_combined_expansion_SHM.png"
  ),
  combined_figure,
  width = 12,
  height = 12,
  units = "in",
  dpi = 600
)


###############################################################################
# 24. Save result tables
###############################################################################

write_csv(
  sample_design,
  file.path(
    OUTPUT_DIR,
    "final_sample_design.csv"
  )
)

write_csv(
  pre_alignment_qc_summary,
  file.path(
    OUTPUT_DIR,
    "pre_alignment_QC_summary.csv"
  )
)

write_csv(
  alignment_qc_summary,
  file.path(
    OUTPUT_DIR,
    "alignment_QC_summary.csv"
  )
)

write_csv(
  primary_sample_summary,
  file.path(
    OUTPUT_DIR,
    "primary_sample_summary.csv"
  )
)

write_csv(
  expansion_class_summary,
  file.path(
    OUTPUT_DIR,
    "expansion_class_summary.csv"
  )
)

write_csv(
  sample_correlations,
  file.path(
    OUTPUT_DIR,
    "per_sample_Spearman_correlations.csv"
  )
)

write_csv(
  sample_correlation_group_summary,
  file.path(
    OUTPUT_DIR,
    "Spearman_summary_by_group.csv"
  )
)

write_csv(
  continuous_model_coefficients,
  file.path(
    OUTPUT_DIR,
    "selected_continuous_model_coefficients.csv"
  )
)

write_csv(
  abundance_slopes_table,
  file.path(
    OUTPUT_DIR,
    "condition_specific_abundance_slopes.csv"
  )
)

write_csv(
  abundance_slope_contrast_table,
  file.path(
    OUTPUT_DIR,
    "abundance_slope_condition_contrast.csv"
  )
)

write_csv(
  class_emmeans_table,
  file.path(
    OUTPUT_DIR,
    "expansion_class_estimated_mutation_proportions.csv"
  )
)

write_csv(
  class_pairwise_table,
  file.path(
    OUTPUT_DIR,
    "expansion_class_pairwise_comparisons.csv"
  )
)

write_csv(
  sensitivity_results,
  file.path(
    OUTPUT_DIR,
    "alignment_length_sensitivity_results.csv"
  )
)

write_csv(
  shm_analysis,
  file.path(
    OUTPUT_DIR,
    "clonotypes_in_primary_SHM_analysis.csv"
  )
)


###############################################################################
# 25. Save models, results, and figures as one R object
###############################################################################

analysis_objects <- list(
  settings = list(
    final_samples = FINAL_SAMPLES,
    primary_minimum_aligned_nt =
      PRIMARY_MIN_ALIGNED_NT,
    sensitivity_thresholds =
      SENSITIVITY_THRESHOLDS,
    minimum_clone_abundance =
      MIN_CLONE_ABUNDANCE,
    minimum_alignment_identity =
      MIN_ALIGNMENT_IDENTITY,
    minimum_alignment_coverage =
      MIN_ALIGNMENT_COVERAGE
  ),
  
  primary_dataset = shm_analysis,
  
  summaries = list(
    pre_alignment_qc =
      pre_alignment_qc_summary,
    alignment_qc =
      alignment_qc_summary,
    primary_sample_summary =
      primary_sample_summary,
    expansion_class_summary =
      expansion_class_summary,
    sample_correlations =
      sample_correlations,
    sensitivity_results =
      sensitivity_results
  ),
  
  models = list(
    continuous_interaction =
      continuous_interaction_model,
    continuous_additive =
      continuous_additive_model,
    selected_continuous =
      selected_continuous_model,
    continuous_model_comparison =
      continuous_model_comparison,
    class_interaction =
      class_interaction_model,
    class_additive =
      class_additive_model,
    selected_class =
      selected_class_model,
    class_model_comparison =
      class_model_comparison
  ),
  
  figures = list(
    continuous =
      figure_continuous,
    expansion_classes =
      figure_expansion_classes,
    sample_correlations =
      figure_correlations,
    sensitivity =
      figure_sensitivity,
    combined =
      combined_figure
  )
)

saveRDS(
  analysis_objects,
  file.path(
    OUTPUT_DIR,
    "complete_expansion_SHM_analysis.rds"
  )
)


###############################################################################
# 26. Save complete statistical output
###############################################################################

sink(
  file.path(
    OUTPUT_DIR,
    "complete_statistical_output.txt"
  )
)

cat("CLONAL EXPANSION AND V-REGION SHM ANALYSIS\n\n")

cat("FINAL SAMPLE DESIGN\n")
print(sample_design)

cat("\nPRE-ALIGNMENT QC\n")
print(pre_alignment_qc_summary)

cat("\nALIGNMENT QC\n")
print(alignment_qc_summary)

cat("\nPRIMARY SAMPLE SUMMARY\n")
print(primary_sample_summary)

cat("\nPER-SAMPLE SPEARMAN CORRELATIONS\n")
print(sample_correlations)

cat("\nCONTINUOUS INTERACTION MODEL\n")
print(summary(continuous_interaction_model))

cat("\nCONTINUOUS ADDITIVE MODEL\n")
print(summary(continuous_additive_model))

cat("\nMODEL COMPARISON\n")
print(continuous_model_comparison)

cat("\nSELECTED PRIMARY MODEL\n")
cat(selected_model_name, "\n")
print(summary(selected_continuous_model))

cat("\nCONDITION-SPECIFIC ABUNDANCE SLOPES\n")
print(abundance_slopes)

cat("\nABUNDANCE-SLOPE CONTRAST\n")
print(abundance_slope_contrast)

cat("\nEXPANSION-CLASS INTERACTION MODEL\n")
print(summary(class_interaction_model))

cat("\nEXPANSION-CLASS ADDITIVE MODEL\n")
print(summary(class_additive_model))

cat("\nEXPANSION-CLASS MODEL COMPARISON\n")
print(class_model_comparison)

cat("\nEXPANSION-CLASS ESTIMATED MARGINAL MEANS\n")
print(class_emmeans)

cat("\nEXPANSION-CLASS PAIRWISE COMPARISONS\n")
print(class_pairwise)

cat("\nALIGNMENT-LENGTH SENSITIVITY ANALYSIS\n")
print(sensitivity_results)

cat("\nSESSION INFORMATION\n")
print(sessionInfo())

sink()


###############################################################################
# 27. Completion message
###############################################################################

cat(
  "\n\n============================================================\n",
  "EXPANSION???SHM ANALYSIS COMPLETED\n",
  "============================================================\n",
  "Primary minimum aligned length: ",
  PRIMARY_MIN_ALIGNED_NT,
  " nt\n",
  "Clonotypes in primary analysis: ",
  nrow(shm_analysis),
  "\n",
  "Biological samples represented: ",
  n_distinct(shm_analysis$Sample),
  "\n",
  "Selected continuous model: ",
  selected_model_name,
  "\n",
  "Results directory: ",
  normalizePath(OUTPUT_DIR),
  "\n",
  "============================================================\n",
  sep = ""
)

print(figure_continuous)
print(figure_expansion_classes)
print(figure_correlations)
print(figure_sensitivity)
print(combined_figure)

