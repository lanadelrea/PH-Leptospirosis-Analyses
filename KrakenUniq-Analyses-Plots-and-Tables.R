library(tidyverse)


# Set folder
path <- "D:/Paper/LEPTO/KrakenUniq/Reports"
path_results <- "D:/Paper/LEPTO/KrakenUniq/Results"

# FILTERING PARAMETERS
# Strict KrakenUniq confidence thresholds
MIN_TAXREADS <- 10
MIN_KMERS <- 1000


# Get KrakenUniq reports
files <- list.files(
  path,
  pattern = "\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No KrakenUniq .tsv reports found in the specified folder.")
}

# Read all KrakenUniq reports
kraken_all <- lapply(files, function(f) {
  
  df <- read.table(
    f,
    header = FALSE,
    sep = "\t",
    skip = 3,
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )
  
  colnames(df) <- c(
    "percent",
    "reads",
    "taxReads",
    "kmers",
    "dup",
    "coverage",
    "taxID",
    "rank",
    "taxName"
  )
  

  # Convert numeric columns
  df <- df %>%
    mutate(
      percent = as.numeric(percent),
      reads = as.numeric(reads),
      taxReads = as.numeric(taxReads),
      kmers = as.numeric(kmers),
      dup = as.numeric(dup),
      coverage = suppressWarnings(as.numeric(coverage)),
      taxID = as.numeric(taxID)
    )

  
  # Clean text
  df$rank <- str_trim(df$rank)
  df$taxName <- str_trim(df$taxName)
  
  
  # Sample name
  df$sample <- tools::file_path_sans_ext(
    basename(f)
  )
  
  df$sample <- gsub(
    "-report-file",
    "",
    df$sample
  )
  
  df
  
}) %>%
  bind_rows()


# Identify host taxa as human
host_taxa <- c(
  "Homo sapiens"
)


# SPECIES-LEVEL DATA
species_data <- kraken_all %>%
  filter(
    rank == "species"
  ) %>%
  mutate(
    is_host = taxName %in% host_taxa,
    is_leptospira = grepl(
      "^Leptospira ",
      taxName
    )
  )


# FILTERING
#
# Three categories:
#
# 1. HIGH_CONFIDENCE:
#    taxReads >= 10 AND kmers >= 1000
#
# 2. DETECTED_LOW_CONFIDENCE:
#    taxReads > 0 AND kmers > 0
#    but does not meet strict threshold
#
# 3. NOT_DETECTED:
#    taxReads == 0 OR kmers == 0
#
# The original KrakenUniq paper used 10 reads + 1000
# unique k-mers as a pathogen-discovery threshold.


lepto_all <- species_data %>%
  filter(is_leptospira) %>%
  mutate(
    
    confidence = case_when(
      
      taxReads >= MIN_TAXREADS &
        kmers >= MIN_KMERS ~
        "HIGH_CONFIDENCE",
      
      taxReads > 0 &
        kmers > 0 ~
        "DETECTED_LOW_CONFIDENCE",
      
      TRUE ~
        "NOT_DETECTED"
    ),
    
    passes_filter =
      taxReads >= MIN_TAXREADS &
      kmers >= MIN_KMERS
  )


# OUTPUT 1: LEPTOSPIRA TABLE
lepto_table <- lepto_all %>%
  select(
    sample,
    taxID,
    taxName,
    reads,
    taxReads,
    kmers,
    dup,
    coverage,
    percent,
    confidence,
    passes_filter
  ) %>%
  arrange(
    sample,
    desc(taxReads)
  )


write.csv(
  lepto_table,
  file.path(
    path_results,
    "Output1_Leptospira_taxReads_filtering.csv"
  ),
  row.names = FALSE
)



# Filtered Leptospira table only
lepto_filtered <- lepto_table %>%
  filter(
    passes_filter
  )


write.csv(
  lepto_filtered,
  file.path(
    path_results,
    "Output1_Leptospira_HIGH_CONFIDENCE.csv"
  ),
  row.names = FALSE
)


# Filtered Leptospira species
lepto_species_filtered <- lepto_all %>%
  filter(
    passes_filter
  ) %>%
  select(
    sample,
    taxName,
    taxReads,
    kmers
  ) %>%
  group_by(
    sample,
    taxName
  ) %>%
  summarise(
    taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    kmers = sum(
      kmers,
      na.rm = TRUE
    ),
    .groups = "drop"
  )



# OUTPUT 2: LEPTOSPIRA + OTHER NON-HOST TAXA
# We use ONLY species-level assignments.
#
# This avoids summing the same reads across:
# species -> genus -> family -> order etc.

nonhost_species <- species_data %>%
  filter(
    !is_host
  )


# Apply strict filter to species-level taxa
nonhost_filtered <- nonhost_species %>%
  filter(
    taxReads >= MIN_TAXREADS,
    kmers >= MIN_KMERS
  )


# Collapse all non-Leptospira taxa into "Other non-host taxa"
figure2_data <- nonhost_filtered %>%
  mutate(
    category = case_when(
      
      is_leptospira ~ taxName,
      
      TRUE ~ "Other non-host taxa"
    )
  ) %>%
  group_by(
    sample,
    category
  ) %>%
  summarise(
    taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Calculate non-host denominator
figure2_data <- figure2_data %>%
  group_by(sample) %>%
  mutate(
    total_nonhost_taxReads =
      sum(taxReads),
    
    Percent =
      100 * taxReads /
      total_nonhost_taxReads
  ) %>%
  ungroup()

# Make sample order
sample_order <- sort(
  unique(
    kraken_all$sample
  )
)

figure2_data$sample <- factor(
  figure2_data$sample,
  levels = sample_order
)

# Plot Output 2
p2 <- ggplot(
  figure2_data,
  aes(
    x = sample,
    y = Percent,
    fill = category
  )
) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Sample",
    y = "Relative abundance of species-level taxonomic assignments (%)",
    fill = "Taxon"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    legend.position =
      "right"
  )


print(p2)

# Save Output 2
ggsave(
  filename = file.path(
    path_results,
    "Output2_Leptospira_other_nonhost.png"
  ),
  plot = p2,
  width = 12,
  height = 7,
  dpi = 300
)


write.csv(
  figure2_data,
  file.path(
    path_results,
    "Output2_Leptospira_other_nonhost_data.csv"
  ),
  row.names = FALSE
)


# OUTPUT 3: NON-HOST TAXA PASSING FILTER
#
# Purpose:
# Show which species-level non-host taxa passed the
# KrakenUniq filtering criteria in each sample.
#
# Leptospira species are shown individually.
# Other non-Leptospira taxa are also retained individually,
# unless minor taxa are intentionally grouped.
#
# Relative abundance denominator:
# Total FILTERED species-level non-host taxReads
# within each sample.
#
# taxReads and kmers are retained in the underlying table.


# 1. PREPARE FILTERED NON-HOST TAXA
figure3_data <- nonhost_filtered %>%
  
  # Only retain taxa with >0 taxReads
  filter(
    taxReads > 0
  ) %>%
  
  # Identify Leptospira
  mutate(
    is_leptospira = grepl(
      "^Leptospira ",
      taxName
    )
  ) %>%
  
  # Combine duplicate sample/taxon entries if present
  group_by(
    sample,
    taxName,
    is_leptospira
  ) %>%
  
  summarise(
    
    # Number of reads uniquely assigned to taxon
    taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    
    # Retain kmer information
    kmers = sum(
      kmers,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


# 2. CALCULATE RELATIVE ABUNDANCE
#
# Denominator:
# total filtered species-level non-host taxReads within each sample.
#
# Host reads are NOT included in the denominator.
#


figure3_data <- figure3_data %>%
  
  group_by(sample) %>%
  
  mutate(
    
    total_nonhost_taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    
    Percent = if_else(
      total_nonhost_taxReads > 0,
      100 * taxReads /
        total_nonhost_taxReads,
      0
    )
    
  ) %>%
  
  ungroup()


# 3. PRESERVE SAMPLE ORDER
figure3_data$sample <- factor(
  figure3_data$sample,
  levels = sample_order
)



# 4. GROUP MINOR NON-LEPTOSPIRA TAXA
#
# TRUE  = group very low-abundance OTHER taxa
# FALSE = show every passing taxon individually
#
# IMPORTANT: Leptospira species are NEVER grouped.
#
# NOTE:
# kmers are retained at the original taxon level.
# We do NOT use grouped kmer values as a biological measure.
#

group_minor_taxa <- TRUE

minor_threshold <- 1
# Taxa contributing <1% of filtered non-host taxReads
# within an individual sample will be grouped.


if (group_minor_taxa) {
  
  figure3_data <- figure3_data %>%
    
    mutate(
      
      plot_taxon = case_when(
        
        # Always retain Leptospira species individually
        is_leptospira ~ taxName,
        
        # Group only minor non-Leptospira taxa
        Percent < minor_threshold ~
          "Other non-Leptospira",
        
        # Keep other sufficiently abundant taxa individually
        TRUE ~ taxName
        
      )
      
    ) %>%
    
    group_by(
      sample,
      plot_taxon
    ) %>%
    
    summarise(
      
      # Sum taxReads for plotting
      taxReads = sum(
        taxReads,
        na.rm = TRUE
      ),
      

      # For individual taxa, kmers are retained.
      kmers = if_else(
        n() == 1,
        first(kmers),
        NA_real_
      ),
      
      .groups = "drop"
      
    )
  
} else {
  
  figure3_data <- figure3_data %>%
    
    mutate(
      plot_taxon = taxName
    )
  
}


# 5. RE-CALCULATE RELATIVE ABUNDANCE AFTER GROUPING
#
# This is important.
#
# If minor taxa were grouped, the grouped values must be recalculated before plotting.
#


figure3_data <- figure3_data %>%
  
  mutate(
    
    # Identify whether the plotted category is Leptospira
    is_leptospira = grepl(
      "^Leptospira ",
      plot_taxon
    )
    
  ) %>%
  
  group_by(sample) %>%
  
  mutate(
    
    total_nonhost_taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    
    Percent = if_else(
      total_nonhost_taxReads > 0,
      100 * taxReads /
        total_nonhost_taxReads,
      0
    )
    
  ) %>%
  
  ungroup()


# 6. DETERMINE TAXON ORDER
#
# Taxa with the greatest total taxReads across all samples appear first in the legend.
#


taxa_order <- figure3_data %>%
  
  group_by(plot_taxon) %>%
  
  summarise(
    
    total_taxReads = sum(
      taxReads,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    desc(total_taxReads)
  ) %>%
  
  pull(plot_taxon)


figure3_data$plot_taxon <- factor(
  figure3_data$plot_taxon,
  levels = taxa_order
)


# 7. CREATE COLOR PALETTE
#
# Leptospira:
#   Gray shades
#
# Other taxa:
#   Distinct colors

# Identify Leptospira categories
lepto_levels <- figure3_data %>%
  
  filter(
    is_leptospira
  ) %>%
  
  distinct(
    plot_taxon
  ) %>%
  
  pull(
    plot_taxon
  ) %>%
  
  as.character()


# Identify other taxa
other_levels <- figure3_data %>%
  
  filter(
    !is_leptospira
  ) %>%
  
  distinct(
    plot_taxon
  ) %>%
  
  pull(
    plot_taxon
  ) %>%
  
  as.character()


# Gray shades for Leptospira
if (length(lepto_levels) > 0) {
  
  lepto_colors <- setNames(
    
    gray.colors(
      n = length(lepto_levels),
      start = 0.25,
      end = 0.70
    ),
    
    lepto_levels
    
  )
  
} else {
  
  lepto_colors <- character(0)
  
}


# Distinct colors for other taxa
if (length(other_levels) > 0) {
  
  other_colors <- setNames(
    
    scales::hue_pal(
      l = 65,
      c = 100
    )(
      length(other_levels)
    ),
    
    other_levels
    
  )
  
} else {
  
  other_colors <- character(0)
  
}


# Combine colors
taxon_colors <- c(
  lepto_colors,
  other_colors
)


# Make sure every plotted taxon has a color
taxon_colors <- taxon_colors[
  levels(figure3_data$plot_taxon)
]

# 8. PLOT OUTPUT 3
p3 <- ggplot(
  
  figure3_data,
  
  aes(
    x = sample,
    y = Percent,
    fill = plot_taxon
  )
  
) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_fill_manual(
    values = taxon_colors,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Sample",
    
    y = paste0(
      "Relative abundance among filtered ",
      "non-host species-level assignments (%)"
    ),
    
    fill = "Taxon"
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    axis.title = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.title = element_text(
      face = "bold"
    )
    
  )


# Display plot
print(p3)


# 9. SAVE OUTPUT 3 PLOT
ggsave(
  
  filename = file.path(
    path_results,
    "Output3_Nonhost_with_Leptospira_highlighted.png"
  ),
  
  plot = p3,
  
  width = 12,
  height = 7,
  dpi = 300
  
)


# 10. SAVE UNDERLYING DATA
#
# Includes:
#   sample
#   taxon
#   taxReads
#   kmers
#   Percent
#   Leptospira status
write.csv(
  
  figure3_data,
  
  file.path(
    path_results,
    "Output3_Nonhost_with_Leptospira_highlighted_data.csv"
  ),
  
  row.names = FALSE
  
)



# 11. SAVE SIMPLE TABLE OF PASSED TAXA
#
# It shows:
#   - sample
#   - taxon
#   - taxReads
#   - kmers
#   - relative percentage
#   - Leptospira status


figure3_passed_taxa <- figure3_data %>%
  
  select(
    sample,
    plot_taxon,
    taxReads,
    kmers,
    Percent,
    is_leptospira
  ) %>%
  
  arrange(
    sample,
    desc(taxReads)
  )


write.csv(
  
  figure3_passed_taxa,
  
  file.path(
    path_results,
    "Output3_Passed_nonhost_taxa.csv"
  ),
  
  row.names = FALSE
  
)



# OUTPUT 4: RELATIVE ABUNDANCE OF LEPTOSPIRA SPECIES
# 
# Denominator:
#
# Total HIGH-CONFIDENCE Leptospira taxReads within each sample.
#
# Therefore:
#
# species taxReads /
# total Leptospira taxReads

figure4_data <- lepto_species_filtered %>%
  
  group_by(sample) %>%
  
  mutate(
    
    total_Leptospira_taxReads =
      sum(
        taxReads,
        na.rm = TRUE
      ),
    
    Percent =
      100 * taxReads /
      total_Leptospira_taxReads
  ) %>%
  
  ungroup()


figure4_data$sample <- factor(
  figure4_data$sample,
  levels = sample_order
)


# Plot Output 4
p4 <- ggplot(
  figure4_data,
  aes(
    x = sample,
    y = Percent,
    fill = taxName
  )
) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Sample",
    y = "Relative abundance within Leptospira (%)",
    fill = "Leptospira species"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


print(p4)


ggsave(
  filename = file.path(
    path_results,
    "Output4_Leptospira_species_relative_abundance.png"
  ),
  plot = p4,
  width = 12,
  height = 7,
  dpi = 300
)


write.csv(
  figure4_data,
  file.path(
    path_results,
    "Output4_Leptospira_species_relative_abundance_data.csv"
  ),
  row.names = FALSE
)