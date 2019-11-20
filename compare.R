library(StructuralVariantAnnotation)
library(tidyverse)
library(stringr)
library(openxlsx)
options(stringsasfactors=FALSE)
theme_set(theme_bw())

maxgap = 50
sample_name = "HG00731"
sample_ordinal = 1

# Issues:
# 1) "SVLEN=." is not valid
#8       127512835       nssv14402396    G       <INS>   .       .       DBVARID;SVTYPE=INS;END=127512836;SVLEN=.;EXPERIMENT=5;SAMPLE=NA19240;REGIONID=nsv3545976
#13      85724540        nssv14389127    C       <INS>   .       .       DBVARID;SVTYPE=INS;END=85724541;SVLEN=.;EXPERIMENT=5;SAMPLE=NA19240;REGIONID=nsv3552124
# 2) <CNV> not supported by StructuralVariantAnnotation)
truth_vcf = readVcf("nstd152.GRCh38.variant_call.vcf.gz")
info(truth_vcf)$SVLEN = as.numeric(info(truth_vcf)$SVLEN)
truth_vcf = truth_vcf[info(truth_vcf)$SVTYPE != "CNV"]
truth_vcf = truth_vcf[info(truth_vcf)$SAMPLE == sample_name]
truth_bpgr = breakpointRanges(truth_vcf)
truth_bpgr$QUAL = 1
seqlevelsStyle(truth_bpgr) = "UCSC"
#truth_bpgr$sample = info(truth_vcf[truth_bpgr$sourceId])$SAMPLE

# Sanity check
truth_csv = read_csv("all_variants_for_nstd152.csv.gz")
names(truth_csv) = str_replace_all(names(truth_csv), " ", "")
truth_csv = truth_csv %>%
  filter(
    AssemblyName == "GRCh38 (hg38)",
    !is.na(truth_csv$VariantSamples),
    truth_csv$VariantSamples == sample_name,
    VariantCalltype != "copy number variation")
length(unique(truth_bpgr$sourceId[!(truth_bpgr$sourceId %in% truth_csv$VariantID)]))
# TODO: what to do with the 6380 regions not in the VCF? Blacklist those regions when comparing?
length(unique(truth_csv$VariantID[!(truth_csv$VariantID %in% truth_bpgr$sourceId)]))

manta_vcf = readVcf("calls/ALL.wgs.UCSD_Manta.20162710.sv.Illumina_high-coverage_PCR-free.genotypes.vcf")
manta_vcf = manta_vcf[geno(manta_vcf)$GT[,sample_name] %in% c("0/1", "1/1")]
manta_bpgr = breakpointRanges(manta_vcf)

lumpy_vcf = readVcf("calls/ALL.lumpy.20160930.genotypes.vcf")
lumpy_vcf= lumpy_vcf[geno(lumpy_vcf)$GT[,sample_name] %in% c("0/1", "1/1")]
lumpy_bpgr = breakpointRanges(lumpy_vcf)

delly_vcf = readVcf("calls/hgsvc.delly.svs.GRCh38.20160931.high_coverage.vcf")
delly_vcf= delly_vcf[geno(delly_vcf)$GT[,sample_name] %in% c("0/1", "1/1")]
delly_bpgr = breakpointRanges(delly_vcf)
delly_bpgr$QUAL = geno(delly_vcf)$RC[delly_bpgr$sourceId, sample_name]

wham_vcf = readVcf("calls/ALL.wham.20160930.genotypes.vcf")
wham_vcf = wham_vcf[geno(wham_vcf)$GT[,sample_name] %in% c("0/1", "1/1")]
wham_bpgr = breakpointRanges(wham_vcf)




# MEI
# melt
# VH
# ForestSV
# wham
# novoBreak

# genomeStrip
# Pindel
# SVelter


# Per trio
#gridss_vcf = readVcf("CHS_2.7.2.vcf")
gridss_vcf = readVcf("PUR_2.7.2.vcf")

gridss_bpgr = breakpointRanges(gridss_vcf)
gridss_begr = breakendRanges(gridss_vcf)
gridss_bpgr$QUAL = geno(gridss_vcf)$QUAL[gridss_bpgr$sourceId, sample_name]
gridss_begr$QUAL = geno(gridss_vcf)$BQ[gridss_begr$sourceId, sample_name]
gridss_begr = gridss_begr[gridss_begr$QUAL >= 50]
gridss_bpgr = gridss_bpgr[gridss_bpgr$QUAL >= 50]
gridss_bpgr = gridss_bpgr[!is.na(gridss_bpgr$svLen) & abs(gridss_bpgr$svLen) + str_length(gridss_bpgr$insSeq) > 50]
gridss_bpgr = gridss_bpgr[gridss_bpgr$partner %in% names(gridss_bpgr)]

gridss1_vcf = readVcf("PUR_1.2.0.vcf.gz")
gridss1_bpgr = breakpointRanges(gridss1_vcf)
gridss1_bpgr$QUAL = gridss1_sample_qual(gridss1_vcf[gridss1_bpgr$sourceId], sample_ordinal)
gridss1_bpgr = gridss1_bpgr[gridss1_bpgr$QUAL >= 50]
gridss1_bpgr = gridss1_bpgr[gridss1_bpgr$partner %in% names(gridss1_bpgr)]

alldata = list(
  bp=list(
    "truth"=truth_bpgr,
    "manta"=manta_bpgr,
    "lumpy"=lumpy_bpgr,
    "gridss2"=gridss_bpgr),
  be=list(
    "gridss2"=gridss_begr,
  ))

filterfunc = NULL
filterfunc = del_only
rocdf = bind_rows(
  rocby(annotate_truth(truth_bpgr, manta_bpgr, NULL, "manta", filterfunc)),
  rocby(annotate_truth(truth_bpgr, lumpy_bpgr, NULL, "lumpy", filterfunc)),
  rocby(annotate_truth(truth_bpgr, delly_bpgr, NULL, "delly", filterfunc)),
  rocby(annotate_truth(truth_bpgr, wham_bpgr, NULL, "wham", filterfunc)),
  rocby(annotate_truth(truth_bpgr, gridss1_bpgr, NULL, "gridss1", filterfunc)),
  rocby(annotate_truth(truth_bpgr, gridss_bpgr, NULL, "gridss2", filterfunc)),
  rocby(annotate_truth(truth_bpgr, gridss_bpgr, gridss_begr, "gridss2 with single breakends", filterfunc))
)
ggplot(rocdf) +
  aes(colour=caller, x=recall, y=precision) +
  geom_point() +
  geom_line() +
  labs(title=paste0(sample_name, " germline calls"), x="Recall", y="Precision")

rtracklayer::export(breakpointgr2pairs(truth_bpgr), con=paste0("processed/", sample_name, ".truth.bedpe"))
rtracklayer::export(breakpointgr2pairs(manta_bpgr), con=paste0("processed/", sample_name, ".manta.bedpe"))
rtracklayer::export(breakpointgr2pairs(gridss_bpgr), con=paste0("processed/", sample_name, ".gridss.bedpe"))
rtracklayer::export(gridss_begr, con=paste0("processed/", sample_name, ".gridss.singlebreakends.bed"))

# TODO: TP/FP fill
annotate_all=bind_rows(
  annotate_truth(truth_bpgr, manta_bpgr, NULL, "manta", NULL),
  annotate_truth(truth_bpgr, lumpy_bpgr, NULL, "lumpy", NULL),
  annotate_truth(truth_bpgr, gridss_bpgr, NULL, "gridss breakpoints", NULL))

annotate_all_tpfp = bind_rows(
    data.frame(caller="truth", eventtype=simpleEventType(truth_bpgr), size=truth_bpgr$svLen, id=names(truth_bpgr)),
    data.frame(caller="manta", eventtype=simpleEventType(manta_bpgr), size=manta_bpgr$svLen, id=names(manta_bpgr)),
    data.frame(caller="lumpy", eventtype=simpleEventType(lumpy_bpgr), size=lumpy_bpgr$svLen, id=names(lumpy_bpgr)),
    data.frame(caller="gridss breakpoints", eventtype=simpleEventType(gridss_bpgr), size=gridss_bpgr$svLen, id=names(gridss_bpgr))) %>%
  left_join(annotate_all, by=c("id", "caller"))

ggplot(annotate_all_tpfp %>% filter(caller!="gridss breakpoints" | qual > 1000)) +
  aes(x=abs(size + 1), colour=!is.na(matchid)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(caller ~ eventtype, scales="free")

gridss_bpgr$tp = FALSE
gridss_bpgr[annotate_all_tpfp %>% filter(caller=="gridss breakpoints", !is.na(matchid)) %>% pull(id)]$tp = TRUE
gridss_bpgr$QUALc = geno(gridss_vcf)$QUAL[gridss_bpgr$sourceId, sample_name]
gridss_bpgr$QUALp1 = geno(gridss_vcf)$QUAL[gridss_bpgr$sourceId, colnames(geno(gridss_vcf)$QUAL)[colnames(geno(gridss_vcf)$QUAL) != sample_name][1]]
gridss_bpgr$QUALp2 = geno(gridss_vcf)$QUAL[gridss_bpgr$sourceId, colnames(geno(gridss_vcf)$QUAL)[colnames(geno(gridss_vcf)$QUAL) != sample_name][2]]
gridss_bpgr = gridss_bpgr[seqnames(gridss_bpgr) %in% seqnames(truth_bpgr)]

# does genotyping fix this issue?
ggplot(as.data.frame(gridss_bpgr) %>% sample_n(25000)) +
  aes(x=QUALp1 + 1, y=QUALp2 + 1, colour=tp) +
  geom_point()
  
ggplot(as.data.frame(gridss_bpgr) %>% sample_n(25000)) +
  aes(y=QUALc / (QUALc + QUALp1 + QUALp2), x=QUALc + 1, colour=tp) +
  geom_point() +
  scale_x_log10()

# random forest:
# - handle imbalanced classes
# - event type
# - genotype likelihoods?
  
# so many supp tables!
xlsx_header = read.xlsx("41467_2018_8148_MOESM3_ESM.xlsx", sheet=1)
illumina_callers = read.xlsx("41467_2018_8148_MOESM3_ESM.xlsx", sheet=9, colNames=TRUE, startRow=2)
write(file="download_illumina_caller_results.sh", paste0("wget -m ", str_replace(illumina_callers$Location, "http://", "ftp://")))


# COLO829 datastore results: gridss 1.2 still running COLO829v1 on datastore

# hgsvc: running CHS_1.2.0 on unix306

# manta-BPI/PURPLE on datastore: not started

# Result: # unexplained copy number transitions in 2018 manta/purple results vs our current gridss/purple result






















