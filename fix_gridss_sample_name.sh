#!/bin/bash
gunzip -c CHS_2.7.2.vcf.gz | sed -s 's/.alt_bwamem_GRCh38DH.20150715.....high_coverage.cram.bam//g' > CHS_2.7.2.vcf
gunzip -c PUR_2.7.2.vcf.gz | sed -s 's/.alt_bwamem_GRCh38DH.20150715.....high_coverage.cram.bam//g' > PUR_2.7.2.vcf
