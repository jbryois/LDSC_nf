#!/usr/bin/env Rscript

list.of.packages <- c("tidyverse", "stringr","crayon","scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(tidyverse)
library(stringr)

files <- list.files(".",pattern=".results",full.names = TRUE)

d <- tibble(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".results","",basename(filename)),
         Trait=sub("_.*","",makenames),
         File=basename(filename)) %>% 
  select(-filename,-makenames) %>%
  unnest(file_contents) %>%
  mutate(Coefficient_P=pnorm(`Coefficient_z-score`,lower.tail=FALSE)) %>%
  select(File,Trait,Category,Prop._SNPs:`Coefficient_z-score`,Coefficient_P) %>%
  mutate(Category=gsub(".bedL2_0|L2_1|L2_0","",Category)) %>%
  arrange(Trait,Enrichment_p)
  
write_tsv(d,path="all.results.txt")