
########################### 
##
##   Making a bed file with homer motif information 
##   8 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################

rm(list = ls())
# cleanup
setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")
met  <- read.delim("target_guides_filtered_homer_motifs.txt", header = FALSE, sep="\t")

names(met) <- c("site_chr","site_start","site_end","motif_chr","motif_start","motif_end","motif_extra","motif_log_odds_score","motif_strand")

#merge homer site info into new column and parse for the highest motif log odds score for each binding site


met$name <- paste(met$site_chr,"_",met$site_start,"_",met$site_end,sep="")
met$rep <- rep("0", length(1:nrow(met)))
met$name_score <- paste(met$name,"_",met$motif_log_odds_score,sep="")

ss <- met[c(4,5,6,12,11,9)]
write.table(ss, file="target_guides_filtered_homer_motifs_all.bed",col.names = F, quote=F, row.names = F, sep="\t")