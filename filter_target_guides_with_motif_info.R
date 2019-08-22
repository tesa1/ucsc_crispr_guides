
########################### 
##
##   Filtering UCSC guides (pre-filtered) with homer motif coverage information 
##   8 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################

rm(list = ls())
setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")
# get motif coverage file
df <- read.delim(file="target_guides_filtered_homer_motifs_all_coverageBed.txt", sep="\t", header=F)

foo <- data.frame(do.call('rbind', strsplit(as.character(df$V4),'_',fixed=TRUE)))
names(foo) <- c("region_chr","region_start", "region_end","region_name", "guide_start","guide_end")
foo$region_name <- paste(foo$region_chr,"_",foo$region_start,"_",foo$region_end,sep="")
foo$guide_name <- paste(foo$region_chr,"_",foo$guide_start,"_",foo$guide_end,sep="")

sdf <- cbind(df, foo$region_name, foo$guide_name)
# get only columns we need
tf <- sdf[c(1,2,3,6,13,14,15)]

names(tf) <- c("guide_chr", "guide_start","guide_end","guide_strand","motif_perc_coverage","binding_region","guide_region")


# get the previous file info file with all ucsc guides filtered (not with homer motif info) for merging together
library(data.table)
gu <- read.csv(file="target_guides_filtered.csv")
gu$guide_region <- paste(gu$chr,"_", gu$start_thick,"_",gu$end_thick,sep="")
# merge the files together
mm <- merge(tf, gu, by='guide_region')
# this file has ALL the regions merged, not just the top 4 by homer motif coverage
write.csv(mm, file="target_guides_filtered_homer_motif_coverage_info.csv")


#select the top 4 guides based on coverage of the homer motif (column V10)
n <-setorder(setDT(tf), -motif_perc_coverage)[, head(.SD, 4), keyby = binding_region]



nn <- merge(n, gu, by='guide_region')
write.csv(nn, file="target_guides_filtered_homer_motif_coverage_info_top4.csv")

# write bed
oo <- as.data.frame(nn)


####
#NOTE, if there aren't 4 guides for the region, this code will repeat the ones that do exist (n <-setorder(setDT(tf), -motif_perc_coverage)[, #head(.SD, 4), keyby = binding_region]).
# so get the guides where there is no duplication (if there are 4 and 3 are duplicated, this will return 1. If there are 4 and 2 are #duplicated, this will return 2.)
# so let's delete the regions that don't have enough guides so we can design for them in a second round.
####

oo$dupes <- duplicated(oo$guide_region)
soo <- subset(oo, dupes=='FALSE')



write.csv(soo, file="target_guides_filtered_homer_motif_coverage_info_top4_filtered.csv")

soo$target_region_guide_region <- paste(soo$target_region_only,"_",soo$start_thick_plus,"_",soo$end_thick,sep="")
soo$rep <-rep("0", length(1:nrow(soo)))
soo$rep_color <- rep("0,165,45",length(1:nrow(soo)))
pp <- soo[c(3,4,5,39,40,14,4,5,41)]
write.table(pp, file="target_guides_filtered_homer_motif_coverage_info_top4_filtered.bed",col.names = FALSE, quote=FALSE, row.names = FALSE, sep="\t")


# so let's delete the regions that don't have enough guides so we can design for them in a second round.
####

x1= xtabs(~ target_region_only, soo)
x1 <- as.data.frame(x1)
write.csv(x1, file="target_guides_filtered_homer_motif_coverage_info_top4_filtered_table.csv")

only_4 <- subset(x1, x1$Freq==4)


final<- subset(soo, (soo$target_region_only %in% only_4$target_region_only))

dim(final)
#[1] 2668   36  
#  = 667 regions with 4 good guides in motifs.

write.csv(final, file="target_guides_filtered_homer_motif_coverage_info_only4_filtered.csv")


#final$target_region_guide_region <- paste(final$target_region_only,"_",final$start_thick_plus,"_",final$end_thick,sep="")
#final$rep <-rep("0", length(1:nrow(final)))
#final$rep_color <- rep("0,165,45",length(1:nrow(final)))
qq<- final[c(3,4,5,39,40,14,4,5,41)]
write.table(qq, file="target_guides_filtered_homer_motif_coverage_info_only4_filtered.bed",col.names = FALSE, quote=FALSE, row.names = FALSE, sep="\t")

### make a figure to show how many regions there are with numbers of guides

require(ggplot2)
pdf(file="target_guides_filtered_homer_motif_coverage_info_top4_filtered.pdf")
ggplot(data=soo, aes(soo$motif_perc_coverage)) + geom_histogram()
dev.off()

