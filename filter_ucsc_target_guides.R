########################### 
##
##   R script to concatenate the guides designed and filter for the unique guides with good specificity and efficiency (mit_specificity_score_only > 70 & $doench2016_percentile_only > 65) (up in mets, H3K27acetylated, found in LNCaP TADs)
##   8 August 2019, Tesa Severson, Netherlands Cancer Institute (t.severson@nki.nl)
##
###########################
rm(list = ls())

# loop to read data and add a column with the filename

setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")

path<-"/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/"
file.names <- dir(path, pattern ="ucsc_guides.txt")
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

#sub1 <- subset(rm, rm$mit_specificity_score_only > 50 & rm$doench2016_percentile_only > 55)


sub2 <- subset(rm, rm$mit_specificity_score_only > 70 & rm$doench2016_percentile_only > 65)


write.csv(sub2, file="target_guides_filtered.csv")


# make a bed file for comparison with CTCF motifs file and add a unique identifier
sub2$target_region_guide_region <- paste(sub2$target_region_only,"_",sub2$start_thick,"_",sub2$end_thick,sep="")
sub2$rep <-rep("0", length(1:nrow(sub2)))
sub2$rep_color <- rep("0,165,45",length(1:nrow(sub2)))
ss <- sub2[c(1,7,8,30,31,6,7,8,32)]

write.table(ss, file="target_guides_filtered.bed",col.names = F, quote=F, row.names = F, sep="\t")


