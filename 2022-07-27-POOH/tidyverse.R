# install.packages("tidyverse")

library(readr)
library(dplyr)
library(tidyr)

# read df
variants <- read.csv("2022-07-27-POOH/combined_tidy_vcf.csv") 

head(variants)
glimpse(variants)

select(variants, sample_id, REF, ALT, DP)
select(variants, -CHROM)
select(variants, contains("i"), -Indiv, -FILTER, POS)

filter(variants, sample_id == "SRR2584863")
filter(variants, REF %in% c("T", "G"))
filter(variants, QUAL >= 100)
filter(variants, INDEL)
filter(variants, !is.na(IDV))
filter(variants, sample_id == "SRR2584863", QUAL >= 100)
filter(variants, sample_id == "SRR2584863", (INDEL | QUAL >= 100))
filter(variants, POS >= 1e6 & POS <= 2e6, !INDEL, QUAL > 200)

variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)


SRR2584863_variants <- variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)

SRR2584863_variants %>% slice(1:6)
SRR2584863_variants %>% slice(10:25)

variants %>%
  filter(sample_id == "SRR2584863" & DP >= 10) %>%
  slice(5:11) %>%
  select(REF, ALT, POS)

variants %>%
  mutate(POLPROB = 1 - (10 ^ -(QUAL/10))) 

variants %>%
  mutate(POLPROB = 1 - (10 ^ -(QUAL/10))) %>%
  arrange(POLPROB) %>%
  select(POLPROB, everything())

variants %>%
  mutate(POLPROB = 1 - 10 ^ -(QUAL/10)) %>%
  arrange(POLPROB) %>%
  select(sample_id, POS, QUAL, POLPROB)


variants %>%
  group_by(sample_id) %>%
  summarize(n = n())

variants %>%
  group_by(sample_id)  %>%
  tally()

variants %>%
  group_by(sample_id) %>%
  count()

variants %>%
  group_by(sample_id) %>%
  summarize(
    n = n(),
    mean_DP = mean(DP),
    median_DP = median(DP),
    min_DP = min(DP),
    max_DP = max(DP))

variants_sum <- variants %>%
  group_by(sample_id, CHROM) %>%
  summarize(mean_DP = mean(DP)) 
variants_sum

variants_wide <- variants_sum %>%
  pivot_wider(names_from = sample_id, values_from = mean_DP)
variants_wide


results_wide <- variants %>%
  group_by(sample_id) %>%
  summarize(
    mean_DP = mean(DP),
    median_DP = median(DP),
    min_DP = min(DP),
    max_DP = max(DP))
results_wide 

results_long <- results_wide %>% 
  pivot_longer(-sample_id, names_to = "stat", values_to = "value")
results_long

results_long %>%
  ungroup() %>%
  pivot_wider(names_from = "stat", values_from = "value")

