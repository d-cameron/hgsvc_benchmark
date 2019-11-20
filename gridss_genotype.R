
manta_gridss_mapping = inner_join(
  annotate_all_tpfp %>% filter(caller == "manta", !is.na(matchid)),
  annotate_all_tpfp %>% filter(caller == "gridss breakpoints", !is.na(matchid)),
  by="matchid",
  suffix=c(".manta", ".gridss"))

manta_gridss_mapping$mantaGT = geno(manta_vcf)$GT[str_replace(manta_gridss_mapping$id.manta, "_bp[0-9]$", ""), sample]
gt = genotype_simple_non_dup_events(gridss_vcf)
colnames(gt) = colnames(geno(gridss_vcf)$GT)
geno(gridss_vcf)$GT = gt


ggt = gt[,sample]
names(ggt) = names(gridss_vcf)

manta_gridss_mapping$gridssGT=ggt[manta_gridss_mapping$id.gridss]

table(manta_gridss_mapping$mantaGT, manta_gridss_mapping$gridssGT)


