library(StructuralVariantAnnotation)
library(tidyverse)
library(openxlsx)
annotate_truth = function(truth_bpgr, caller_bpgr, caller_begr, caller_name, breakpoint_filter_function=NULL) {
  caller_bpgr = caller_bpgr[seqnames(caller_bpgr) == seqnames(partner(caller_bpgr))]
  if (!is.null(breakpoint_filter_function)) {
    truth_bpgr = breakpoint_filter_function(truth_bpgr)
  }
  .annotate_truth_do_overlap = function(findOverlapFunction, gr) {
    if (is.null(gr)) {
      return(NULL)
    } else {
      gr = gr[seqnames(gr) %in% seqnames(truth_bpgr)]
      if (!is.null(gr$partner)) {
        gr = gr[is.na(gr$partner) | (gr$partner %in% names(gr))]
      }
      if (length(gr) == 0) {
        return(NULL)
      }
      cdf = data.frame(
        callset="caller",
        caller=caller_name,
        id=names(gr),
        qual=gr$QUAL,
        matchid=NA)
      hits = findOverlapFunction(gr, truth_bpgr, ignore.strand=FALSE, maxgap=maxgap) %>%
        as.data.frame() %>%
        mutate(qual=gr$QUAL[queryHits]) %>%
        arrange(qual)
      cdf$qual[hits$queryHits] = hits$qual
      cdf$matchid[hits$queryHits] = names(truth_bpgr)[hits$subjectHits]
      tdf$qual[hits$subjectHits] <<- hits$qual
      tdf$matchid[hits$subjectHits] <<- names(gr)[hits$queryHits]
      return(cdf)
    }
  }
  tdf = data.frame(
    callset="truth",
    caller=caller_name,
    id=names(truth_bpgr),
    qual=rep(NA, length(truth_bpgr)),
    matchid=NA)
  becdf = .annotate_truth_do_overlap(findOverlaps, caller_begr)
  bpcdf = .annotate_truth_do_overlap(findBreakpointOverlaps, caller_bpgr)
  return(bind_rows(tdf, bpcdf, becdf))
}

rocby = function(df, sentinal.qual=-1) {
  ntruth = sum(df$callset=="truth")
  df %>%
    replace_na(list(qual=sentinal.qual)) %>%
    filter(
      # duplicate true positives are ignored
      !(!is.na(matchid) & callset=="caller"),
      # Remove false negatives
      !(callset=="truth" & qual == sentinal.qual)) %>%
    mutate(
      fp=callset=="caller",
      tp=callset=="truth") %>%
    arrange(desc(qual)) %>%
    mutate(tp=cumsum(tp), fp=cumsum(fp)) %>%
    group_by(caller, qual) %>%
    summarise(tp=max(tp), fp=max(fp)) %>%
    ungroup() %>%
    mutate(precision=tp/(tp+fp), recall=tp/ntruth)
}

missed_by_gridss = function(truth_bpgr, caller_bpgr, gridss_bpgr) {
  gdf = bind_rows(
    annotate_truth(truth_bpgr, caller_bpgr, NULL, "caller"),
    annotate_truth(truth_bpgr, gridss_bpgr, NULL, "gridss"))
  gdf %>% filter(callset=="truth" & !is.na(qual)) %>%
    group_by(id) %>%
    filter(n() == 1 & any(caller != "gridss"))
}
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}
filter_breakpoint_to = function(bpgr, eventtype=c("DEL", "DUP", "INV", "INS")) {
  eventtype = match.arg(eventtype)
  bpgr = bpgr[!is.na(bpgr$partner)]
  bpgr$eventtype = simpleEventType(bpgr)
  return(bpgr[bpgr$eventtype %in% eventtype])
}
del_only = function(gr) filter_breakpoint_to(gr, "DEL")
dup_only = function(gr) filter_breakpoint_to(gr, "DUP")
ins_only = function(gr) filter_breakpoint_to(gr, "INS")
inv_only = function(gr) filter_breakpoint_to(gr, "INV")

estimated_event_size = function(vcf) {
  rr = SummarizedExperiment::rowRanges(vcf)
  assertthat::assert_that(all(elementNROWS(rr$ALT) == 1))
  alt = unlist(rr$ALT)
  partner_alt = stringr::str_match(alt, "^[^\\]\\[]*[\\]\\[]([^:]+):([0-9]+)[\\]\\[][^\\]\\[]*$")
  # [,2] partner chr
  # [,3] partner position
  return(ifelse(seqnames(rr) != partner_alt[,2], NA, abs(as.numeric(partner_alt[,3]) - start(rr))))
}

simple_genotype_likelihoods = function(vcf, exclude_refpair_size=500, ref_noise_rate=0.01, alt_noise_rate=ref_noise_rate) {
  vf = geno(vcf)$VF
  bvf = geno(vcf)$BVF
  ref = geno(vcf)$REF
  rp = geno(vcf)$RP
  
  ref[is.na(ref)] = 0
  rp[is.na(rp)] = 0
  vf[is.na(vf)] = 0
  bvf[is.na(bvf)] = 0
  
  size = estimated_event_size(vcf)
  is_small_event = matrix(rep(!is.na(size) & size < exclude_refpair_size, 3), ncol=3)
  # ref fragments
  rf = ref + ifelse(is_small_event, 0, rp)
  # alt fragments
  chralt = as.character(unlist(rowRanges(vcf)$ALT), stringr::fixed("."))
  af = ifelse(startsWith(chralt, ".") | endsWith(chralt, "."), bvf, vf)
  rr = dbinom(af, af+rf, 0/2 + ref_noise_rate, log=TRUE)
  ra = dbinom(af, af+rf, 1/2, log=TRUE)
  rra = dbinom(af, af+rf, 1/3, log=TRUE)
  aa = dbinom(af, af+rf, 2/2 - alt_noise_rate, log=TRUE)
  return(list(af=af, rr=rr, ra=ra, rra=rra, aa=aa))
}
genotype_simple_non_dup_events = function(vcf) {
  gl = simple_genotype_likelihoods(vcf)
  glmax = pmax(gl$rr, gl$ra, gl$aa)
  gt = ifelse(gl$rr == glmax, "0/0", ifelse(gl$ra == glmax, "0/1", "1/1"))
  return(gt)
}

gridss1_sample_qual = function(vcf, ordinal) {
  i = info(vcf)
  sr = as.matrix(i$SR)
  rp = as.matrix(i$RP)
  assr = as.matrix(i$ASSR)
  asrp = as.matrix(i$ASRP)
  sr[is.na(sr)] = 0
  rp[is.na(rp)] = 0
  assr[is.na(assr)] = 0
  asrp[is.na(asrp)] = 0
  q = sr + rp + assr + asrp
  qual = rowRanges(vcf)$QUAL * q[,ordinal] / rowSums(q)
}










