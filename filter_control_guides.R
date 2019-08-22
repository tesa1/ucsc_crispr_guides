########################### 
##
##   R script to concatenate the control guides designed and filter for the unique guides with good specificity and efficiency (mit_specificity_score_only > 70 & $doench2016_percentile_only > 65) (up in mets, H3K27acetylated, found in LNCaP TADs)
##   Connect the final filtered guides with their control guides
##   13 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################
rm(list = ls())
setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")

# loop to read data and add a column with the filename

path<-"/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/"
file.names <- dir(path, pattern ="control_ucsc_guides.txt")
print(file.names)
dataset <- data.frame()
#Read in all files and add name column
for(i in 1:length(file.names)){
  temp_data <- read.delim(file.names[i],header=FALSE, sep="\t", stringsAsFactors = FALSE, colClasses = c("character","numeric","numeric","character","numeric","character","numeric","numeric","character","character","character","character","character","character","character","character","numeric","numeric","character","character","numeric"))
  temp_data$fname <- file.names[i]
  dataset <- rbind(dataset, temp_data)
}

wmet <- dataset
names(wmet) <- c("chr","start","end","name","uint_score","strand","start_thick","end_thick","doench2016_fusi_score_color","moreno_mateos_score_color","mit_spec_score_color","guide_sequence","protospacer_adjascent_motif","mit_spec_score","efficiency_doench2016_score","efficiency_moreno_mateos_score","efficiency_doench2014_score","bae_out_of_frame_score","label_for_mouse_over","unknown","emtpy","target_region")

# filter for guides which are unique in the genome
rm <- subset(wmet, wmet$mit_spec_score!="This guide sequence is not unique in the genome. The specificity scores were not determined.")

# fix the columns so we can use the numbers for filtering
foo <- data.frame(do.call('rbind', strsplit(as.character(rm$label_for_mouse_over),',',fixed=TRUE)))


rm$mit_specificity_score_only <- gsub("MIT Spec. Score: ", "", foo$X1, fixed = TRUE)
rm$doench2016_percentile <- gsub("Doench 2016: ", "", foo$X2, fixed = TRUE)
rm$doench2016_percentile_only <- gsub("%","", rm$doench2016_percentile, fixed = TRUE)
rm$doench2016_percentile_only <- gsub(" ","", rm$doench2016_percentile_only, fixed = TRUE)
rm$moreno_mateos_score_percentile <- gsub("Moreno-Mateos: ", "", foo$X3, fixed = TRUE)
rm$moreno_mateos_score_percentile_only <- gsub("%", "", rm$moreno_mateos_score_percentile, fixed = TRUE)
rm$moreno_mateos_score_percentile_only <- gsub(" ", "", rm$moreno_mateos_score_percentile_only, fixed = TRUE)

rm$moreno_mateos_score_percentile_only <- as.numeric(rm$moreno_mateos_score_percentile_only)
rm$doench2016_percentile_only <- as.numeric(rm$doench2016_percentile_only)
rm$mit_specificity_score_only <-as.numeric(rm$mit_specificity_score_only)
bar <- data.frame(do.call('rbind', strsplit(as.character(rm$target_region),' ',fixed=TRUE)))
names(bar) <- c('X2')
rm$target_region_only <- bar$X2
rm$target_region_only <- gsub("_ucsc_guides.txt","_site",rm$target_region_only,fixed=TRUE)
rm$start_thick_plus1 <- rm$start_thick+1


# Mouse-over a target site to show predicted specificity and efficiency scores:

# The MIT Specificity score summarizes all off-targets into a single number from 0-100. The higher the number, the fewer off-target effects are expected. We recommend guides with an MIT specificity > 50.
# The efficiency score tries to predict if a guide leads to rather strong or weak cleavage. According to (Haeussler et al. 2016), the Doench 2016 Efficiency score should be used to select the guide with the highest cleavage efficiency when expressing guides from RNA PolIII Promoters such as U6. Scores are given as percentiles, e.g. "70%" means that 70% of mammalian guides have a score equal or lower than this guide. The raw score number is also shown in parentheses after the percentile.
# The Moreno-Mateos 2015 Efficiency score should be used instead of the Doench 2016 score when transcribing the guide in vitro with a T7 promoter, e.g. for injections in mouse, zebrafish or Xenopus embryos. The Moreno-Mateos score is given in percentiles and the raw value in parentheses, see the note above.
#colors in the mouse over 
#unable to calculate Doench/Fusi 2016 efficiency score
#	low predicted cleavage: Doench/Fusi 2016 Efficiency percentile <= 30
#	medium predicted cleavage: Doench/Fusi 2016 Efficiency percentile > 30 and < 55
#	high predicted cleavage: Doench/Fusi 2016 Efficiency > 55

# ex (sub1 <- subset(rm, rm$mit_specificity_score_only > 50 & rm$doench2016_percentile_only > 55))



sub2 <- subset(rm, rm$mit_specificity_score_only > 70 & rm$doench2016_percentile_only > 65)
write.csv(sub2, file="control_guides_filtered.csv")


## make a bed file for comparison with CTCF motifs file and add a unique identifier
sub2$target_region_guide_region <- paste(sub2$target_region_only,"_",sub2$start_thick,"_",sub2$end_thick,sep="")
sub2$rep <-rep("0", length(1:nrow(sub2)))
sub2$rep_color <- rep("0,165,45",length(1:nrow(sub2)))
ss <- sub2[c(1,7,8,30,31,6,7,8,32)]

write.table(ss, file="control_guides_filtered.bed",col.names = F, quote=F, row.names = F, sep="\t")


## select 1 control guide per target region based on doench2016_percentile_only score.

library(data.table)
p<-setorder(setDT(sub2), -doench2016_percentile_only)[, head(.SD, 1), keyby = target_region_only]


q <- as.data.frame(p)
write.csv(q, file="control_guides_filtered_top1.csv")
q$target_region_guide_region <- paste(q$target_region_only,"_",q$start_thick_plus,"_",q$end_thick,sep="")
q$rep <-rep("0", length(1:nrow(q)))
q$rep_color <- rep("200,150,50",length(1:nrow(q)))

yy <- q[c(2,8,9,30,31,7,8,9,32)]
write.table(yy, file="control_guides_filtered_top1.bed",col.names = F, quote=F, row.names = F, sep="\t")



## connect the control guides with existing top4_filtered target guides (from design_target_guides_ctcf_motif_only4_filtered.py)

goguides <- read.csv("target_guides_filtered_homer_motif_coverage_info_only4_filtered.csv")
gocontrols <- read.csv("control_guides_filtered_top1.csv")

#remove text in gocontrols and get a number to connect with the good guides (control end + 250bp)
foo <- data.frame(do.call('rbind', strsplit(as.character(gocontrols$target_region_only),'_',fixed=TRUE)))

foo$target_region_start <- as.numeric(as.character(foo$X3))+250
gocontrols$target_chr_target_region_start <- paste(gocontrols$chr,"_",foo$target_region_start, sep="")

# do the same for the good guides so we can merge the two lists
bar <- data.frame(do.call('rbind', strsplit(as.character(goguides$binding_region),'_',fixed=TRUE)))
goguides$target_chr_target_region_start <- paste(goguides$chr,"_",bar$X2, sep="")


# merge the unique good guides regions with the good control guides to get the final list of the control guides for the first round

mcguides <- merge(gocontrols, goguides, by.x='target_chr_target_region_start', by.y='target_chr_target_region_start')

mcguides <- mcguides[c(3:34,71)] 

write.csv(mcguides, file="control_guides_filtered_top1_merged_with_only4_target_guides.csv")

# make a bed file for comparison with homer motifs file and add a unique identifier

mcguides$control_region_guide_region_target_region <- paste(mcguides$target_region_guide_region.x,"_",mcguides$target_region_only.y,sep="")
um <- unique(mcguides)
write.csv(um, file="control_guides_filtered_top1_merged_with_only4_target_guides_unique.csv")
tt <- um[c(2,8,9,34,31,7,8,9,32)]


### get the control guides list and the target guides lists and combine
# get master  list of 644 good target guides with good control guides
controls <-read.csv(file="control_guides_filtered_top1_merged_with_only4_target_guides_unique.csv")
mc <- controls[c(2,3,9,10,8,14,15,25,27,29,35)]
write.csv(mc, file="target_guides_filtered_homer_motif_coverage_info_only4_filtered_final_list_controls.csv")



####################################



# get list of 667 good guides with good homer motif coverage
targets <- read.csv(file="target_guides_filtered_homer_motif_coverage_info_only4_filtered.csv")

# find the final targets list (those final target guides that have control guides)
ftargets <- targets[targets$target_region_only %in% controls$target_region_only.y,]
ft <- ftargets[c(3,4,16,17,15,21,22,32,34,36,40)]
write.csv(ft, file="target_guides_filtered_homer_motif_coverage_info_only4_filtered_final_list_targets.csv")

# find guides from 667 list not in guides from master list =23 regions (4 guides per region = 92)
# these guides now need to have controls guides designed for them.
needs_control <- targets[!targets$target_region_only %in% controls$target_region_only.y,]

write.csv(needs_control, file="target_guides_filtered_homer_motif_coverage_info_only4_filtered_need_controls.csv")
# write a bed file to make a shell script to get guides

fbar <- data.frame(do.call('rbind', strsplit(as.character(needs_control$binding_region),'_',fixed=TRUE)))
ufbar <- unique(fbar)
ufbar$X2 <- levels(droplevels(ufbar$X2))
ufbar$X3 <- levels(droplevels(ufbar$X3))

ufbar$control_begin <- as.numeric(ufbar$X2)-750
ufbar$control_end <- as.numeric(ufbar$X2)-500
# make a shell script command
ufbar$f <- paste('/home/t.severson/tools/kent_utils/bigBedToBed -chrom=',ufbar$X1,' -start=',ufbar$control_begin,' -end=',ufbar$control_end,' /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb ',ufbar$X1,'_',ufbar$control_begin,'_',ufbar$control_end,'_ucsc_control_guides_round3.txt', sep="")
ufbar$g <- paste(ufbar$X1,'_',ufbar$control_begin,'_',ufbar$control_end, sep="")

csh <- ufbar[,6]
write.table(csh, file="target_guides_need_controls_round3.sh", quote=FALSE, col.names = FALSE, row.names = FALSE)
cconf <- ufbar[,7]
cconf <- as.data.frame(cconf)
names(cconf) = c('control_site')
write.table(cconf, file="target_guides_need_controls_round3.tsv", quote=FALSE, col.names = TRUE, row.names = FALSE)



# get the final regions that have no target guides designed and no controls designed and write target file
bed <- read.delim('input_sites.bed', sep="\t", header = FALSE)

bed$target_region_only <- paste(bed$V1,'_',bed$V2,'_',bed$V3,'_site', sep="")

needs_guide <- bed[!bed$target_region_only %in% targets$target_region_only,]

# use bed file to write schell script 
needs_guide$f <- paste('/home/t.severson/tools/kent_utils/bigBedToBed -chrom=',needs_guide$V1,' -start=',needs_guide$V2,' -end=',needs_guide$V3,' /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb ',needs_guide$V1,'_',needs_guide$V2,'_',needs_guide$V3,'_ucsc_guides_round2.txt', sep="")
needs_guide$g <- paste(needs_guide$V1,'_',needs_guide$V2,'_',needs_guide$V3, sep="")
# write out the shell script for running.
sh <- needs_guide[,5]
write.table(sh, file="target_guides_need_guides_round2.sh", quote=FALSE, col.names = FALSE, row.names = FALSE)
conf <- needs_guide[,6]
conf <- as.data.frame(conf)
names(conf) = c('site')
write.table(conf, file="target_guides_need_guides_round2.tsv", quote=FALSE, col.names = TRUE, row.names = FALSE)



# get the final regions that have no target guides designed and no controls designed and make control file
bed <- read.delim('input_sites.bed', sep="\t", header = FALSE)

bed$target_region_only <- paste(bed$V1,'_',bed$V2,'_',bed$V3,'_site', sep="")

needs_guide <- bed[!bed$target_region_only %in% targets$target_region_only,]


needs_guide$control_begin <- as.numeric(needs_guide$V2)-500
needs_guide$control_end <- as.numeric(needs_guide$V2)-250
# make a shell script command

# use bed file to write schell script 
needs_guide$f <- paste('/home/t.severson/tools/kent_utils/bigBedToBed -chrom=',needs_guide$V1,' -start=',needs_guide$control_begin,' -end=',needs_guide$control_end,' /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb ',needs_guide$V1,'_',needs_guide$control_begin,'_',needs_guide$control_end,'_ucsc_control_guides_round2.txt', sep="")
needs_guide$g <- paste(needs_guide$V1,'_',needs_guide$control_begin,'_',needs_guide$control_end, sep="")
# write out the shell script for running.
sh <- needs_guide[,7]
write.table(sh, file="target_guides_need_controls_round2.sh", quote=FALSE, col.names = FALSE, row.names = FALSE)
conf <- needs_guide[,8]
conf <- as.data.frame(conf)
names(conf) = c('control_site')
write.table(conf, file="target_guides_need_controls_round2.tsv", quote=FALSE, col.names = TRUE, row.names = FALSE)



