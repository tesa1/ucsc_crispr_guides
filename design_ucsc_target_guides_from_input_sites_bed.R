
########################### 
##
##   Making shell script and config file for snakemake to design UCSC guides for 720 CTCF sites (up in mets, H3K27acetylated, found in LNCaP TADs)
##   8 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################


# cleanup
rm(list = ls())
setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")
# get the bed file for which you want to design guides 
bed <- read.delim("input_sites.bed", header = FALSE,sep="\t")

# make a shell script command
bed$f <- paste('/home/t.severson/tools/kent_utils/bigBedToBed -chrom=',bed$V1,' -start=',bed$V2,' -end=',bed$V3,' /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb ',bed$V1,'_',bed$V2,'_',bed$V3,'_ucsc_guides.txt', sep="")
bed$g <- paste(bed$V1,'_',bed$V2,'_',bed$V3, sep="")
# write out the shell script for running.
sh <- bed[,4]
write.table(sh, file="target_guides_from_input_sites_bed.sh", quote=FALSE, col.names = FALSE, row.names = FALSE)
conf <- bed[,5]
conf <- as.data.frame(conf)
names(conf) = c('site')
write.table(conf, file="target_guides_from_input_sites_bed.tsv", quote=FALSE, col.names = TRUE, row.names = FALSE)