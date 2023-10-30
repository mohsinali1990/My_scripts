# MIT License
# Copyright (c) 2023 Mohsin Ali

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# Load required libraries
library(tidyverse)
library(data.table)

# Clear the environment
rm(list = ls())

# Input file comma-delimited
TASSEL_genoSummary <- file.choose()

# Read the data
tasselSiteSum <- fread(TASSEL_genoSummary, header = TRUE)

# Print data dimensions and structure
print(dim(tasselSiteSum))
glimpse(tasselSiteSum)

# Message for SNP features calculation
message(paste0("SNP features calculation will be running shortly"))

# Calculate SNP features
SNP_Features <- tasselSiteSum %>%
  mutate(
    Total_sites = n_distinct(`Site Name`),
    Minor_Allele_freq_Sqr = (`Minor Allele Frequency` ^ 2),
    Major_Allele_freq_Sqr = (`Major Allele Frequency` ^ 2),
    Major_Allele_freq = `Major Allele Gametes` / (2 * `Number of Taxa`),
    Minor_Allele_freq = `Minor Allele Gametes` / (2 * `Number of Taxa`),
    Obs_het = `Number Heterozygous` / `Number of Taxa`,
    Exp_het = 2 * `Major Allele Frequency` * `Minor Allele Frequency`,
    Nei_Gene_Diversity = 1 - (`Major Allele Frequency` ^ 2) - (`Minor Allele Frequency` ^ 2),
    PIC = 1 - (Major_Allele_freq_Sqr + Minor_Allele_freq_Sqr) - (2 * Major_Allele_freq_Sqr * Minor_Allele_freq_Sqr)
  ) %>%
  select(
    `Site Name`,
    Chromosome,
    Total_sites, 
    Major_Allele_freq, 
    Minor_Allele_freq, 
    Obs_het, 
    Exp_het, 
    Nei_Gene_Diversity, 
    PIC
  )


# Save SNP features to a CSV file
fwrite(SNP_Features, paste0(TASSEL_genoSummary, "_SNP_Features.csv"), row.names = TRUE)

# Message for completion
message(paste0("Writing results done"))

# Display SNP features
glimpse(SNP_Features)

# Calculate overall statistics
overall_stats <- SNP_Features %>%
  summarise(
    Count_site = n(),
    Avg_MinorAllFreq = mean(Minor_Allele_freq),
    Min_MinorAllFreq = min(Minor_Allele_freq),
    Max_MinorAllFreq = max(Minor_Allele_freq),
    Avg_Obs_het = mean(Obs_het),
    Min_Obs_het = min(Obs_het),
    Max_Obs_het = max(Obs_het),
    Avg_Exp_het = mean(Exp_het),
    Min_Exp_het = min(Exp_het),
    Max_Exp_het = max(Exp_het),
    Avg_Nei_Gene_Diversity = mean(Nei_Gene_Diversity),
    Min_Nei_Gene_Diversity = min(Nei_Gene_Diversity),
    Max_Nei_Gene_Diversity = max(Nei_Gene_Diversity),
    Avg_PIC = mean(PIC),
    Min_PIC = min(PIC),
    Max_PIC = max(PIC)
  )

# Print overall statistics
print(overall_stats)
print(paste0("Saving overall SNP features"))

# Save overall statistics to a CSV file
fwrite(overall_stats, paste0(TASSEL_genoSummary, "_SNP_Features_overall.csv"), row.names = TRUE)

# Calculate chromosome-wise statistics
chrwise <- SNP_Features %>%
  group_by(Chromosome) %>%
  summarise(
    Count_site = n(),
    Avg_MinorAllFreq = mean(Minor_Allele_freq),
    Avg_MajorAllFreq = mean(Major_Allele_freq),
    Avg_Obs_het = mean(Obs_het),
    Avg_Exp_het = mean(Exp_het),
    Avg_Nei_Gene_Diversity = mean(Nei_Gene_Diversity),
    Avg_PIC = mean(PIC)
  )

# Print chromosome-wise statistics
print(chrwise)
print(paste0("Saving chromosome-wise SNP features"))

# Save chromosome-wise statistics to a CSV file
fwrite(chrwise, paste0(TASSEL_genoSummary, "_SNP_Features_Chrwise.csv"), row.names = TRUE)

# Message for analysis completion
message(paste0("Analysis Done"))
