
########################### 
##
##   Making shell script and config file for snakemake to design control region UCSC guides for 720 CTCF sites (up in mets, H3K27acetylated, found in LNCaP TADs)
##   13 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################


# cleanup
rm(list = ls())
setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")
# get the bed file for which you want to design guides 
bed <- read.delim("input_sites.bed", header = FALSE,sep="\t")
bed$control_begin <- bed$V2-500
bed$control_end <- bed$V2-250
# make a shell script command
bed$f <- paste('/home/t.severson/tools/kent_utils/bigBedToBed -chrom=',bed$V1,' -start=',bed$control_begin,' -end=',bed$control_end,' /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb ',bed$V1,'_',bed$control_begin,'_',bed$control_end,'_control_ucsc_guides.txt', sep="")
bed$g <- paste(bed$V1,'_',bed$control_begin,'_',bed$control_end, sep="")
# write out the shell script for running.
sh <- bed[,6]
write.table(sh, file="control_guides_from_input_sites_bed.sh", quote=FALSE, col.names = FALSE, row.names = FALSE)
conf <- bed[,7]
conf <- as.data.frame(conf)
names(conf) = c('control_site')
write.table(conf, file="control_guides_from_input_sites_bed.tsv", quote=FALSE, col.names = TRUE, row.names = FALSE)