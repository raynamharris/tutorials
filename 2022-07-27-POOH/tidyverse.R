# install.packages("dplyr")

library(dplyr)
library(readr)
library(tidyr)

# read df
variants <- read.csv("2022-07-27-POOH/combined_tidy_vcf.csv") 

# make a tibble for easy viewing
variants <- as_tibble(variants)

# view
head(variants)
str(variants)

# select columns
select(variants, sample_id, REF, ALT, DP)

select(variants, -CHROM)

select(variants, ends_with("B"))

select(variants, contains("i"), -Indiv, -FILTER, POS)

# filter rows

filter(variants, sample_id == "SRR2584863")

filter(variants, REF %in% c("T", "G"))

filter(variants, QUAL >= 100)

filter(variants, INDEL)

filter(variants, !is.na(IDV))

filter(variants, sample_id == "SRR2584863", QUAL >= 100)

filter(variants, sample_id == "SRR2584863", (INDEL | QUAL >= 100))

filter(variants, POS >= 1e6 & POS <= 2e6, !INDEL, QUAL > 200)

## pipes

variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)


SRR2584863_variants <- variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)

SRR2584863_variants

SRR2584863_variants %>% head()

variants %>%
  mutate(POLPROB = 1 - (10 ^ -(QUAL/10)))

variants %>%
  mutate(POLPROB = 1 - 10 ^ -(QUAL/10)) %>%
  select(sample_id, POS, QUAL, POLPROB)


variants %>%
  group_by(sample_id) %>%
  summarize(n())

variants %>%
  group_by(ALT) %>%
  count()

variants %>%
  count(ALT)

variants %>%
  count(sample_id)


variants %>%
  group_by(sample_id) %>%
  summarize(
    mean_DP = mean(DP),
    median_DP = median(DP),
    min_DP = min(DP),
    max_DP = max(DP))

# pivot

variants_wide <- variants %>%
  group_by(sample_id, CHROM) %>%
  summarize(mean_DP = mean(DP)) %>%
  pivot_wider(names_from = sample_id, values_from = mean_DP)
variants_wide


variants_wide %>%
  pivot_longer(-CHROM, names_to = "sample_id", values_to = "mean_DP")
