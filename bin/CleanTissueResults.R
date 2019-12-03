#!/usr/bin/env Rscript

library(tidyverse)

files <- list.files(".",pattern=".results",full.names = TRUE)

d <- tibble(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".results","",basename(filename)),
         Trait=sub("_.*","",makenames),
         Tissue=gsub("^_","",str_extract(makenames, "_.*"))) %>% 
  select(-filename,-makenames) %>%
  unnest(file_contents) %>%
  filter(Category=="L2_1") %>% mutate(P=pnorm(`Coefficient_z-score`,lower.tail=FALSE)) %>%
  select(Trait,Tissue,Coefficient,Coefficient_std_error,`Coefficient_z-score`,P) %>%
  arrange(Trait,P)
  
write_tsv(d,path="all.results.txt")