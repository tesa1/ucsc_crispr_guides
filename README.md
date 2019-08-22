# Designing guides for a CRISPR-Cas9 experiment to knockout binding sites (including control guides)

This is a manual of how to use publically available tools to obtain previously identified guides and parse for those
that should work. The aim is to design 4 guides per region that are in a CTCF motif and have 1 control guide per region.

Steps for designing guides:
  - Search UCSC database of all known guides 
  - Filter out guides which are not unique in the genome
  - Filter for those that have the best possibility of working
  - Determine coordinates for motif of interest in the regions
  - Select guides which are in the motif of interest
  
Steps for designing control guides:
  - Search UCSC database of all known guides 
  - Filter out guides which are not unique in the genome
  - Filter for those that have the best possibility of working

  
 ## Use of UCSC All predicted Super track ##
Information can be found here:
https://genome-euro.ucsc.edu/cgi-bin/hgTables

It's not necessary but data can be downloaded from https://hgdownload-test.gi.ucsc.edu/gbdb/hg19/crisprAll/crispr.bb.
And it can queried in the UCSC table browser under Genes and Gene Predictions -> CRISPR Targets (https://genome.ucsc.edu/cgi-bin/hgTables) file queried is: https://hgdownload-test.gi.ucsc.edu/gbdb/hg19/crisprAll/crispr.bb

Downloaded data can be extracted with the kent command line tools (bigBedToBed), available:
   http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
   

## Note all commmands expect specific locations of files and path locations. 

These files are dependent on a working directory path. For your own use, you will have to change lines in the python scripts like this:  "WDIR="/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/" and in the Rscripts like this: 
"setwd("/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/")"
to your own path.

The shell scripts produced assume bigBedToBed is here:/home/t.severson/tools/kent_utils/bigBedToBed and the crispr.bb file is here: /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb. You will have to change these locations in the Rscripts to your relevant locations in order to run the shell scripts.

## Query the bed file of target sites and design a shell script to obtain guides in the target regions and the control regions (250-500bp upstream of target site)

 ```bash
snakemake -s make_tsv_and_shell_targets_controls.py 
# note this snakemake requires two R scripts: 1) design_ucsc_target_guides_from_input_sites_bed.R and 
# 2) design_ucsc_control_guides_from_input_sites_bed.R 
```

## Obtain UCSC guides in target regions and and filter for uniqueness, specificity, effeciency and motif coverage 
After initial filtering, this script will identify the top 4 guides per target region based on motif coverage percentage.

homer needs to be installed (eg. /home/t.severson/tools/homer/) and a motif file is necessary with the path defined in design_target_guides_homer_motif_only4_filtered.py

```bash
snakemake -s design_target_guides_homer_motif_only4_filtered.py
# note this snakemake file requires two R scripts 1) filter_ucsc_target_guides.R, 2) make_homer_motif_info_bed.R and
# 3) filter_target_guides_with_motif_info.R.
```

## Obtain UCSC guides in control regions for targets and filter for uniqueness, specificity and effeciency
After initial filtering, this script will identify the top 1 control guides per target region based on motif coverage percentage.


```bash
snakemake -s design_target_guides_homer_motif_only4_filtered.py
# note this snakemake file requires two R scripts 1) filter_ucsc_target_guides.R and 
# 2) filter_target_guides_with_motif_info.R

```
sh ctcf_sites_ucsc_guides.sh 
mkdir guides
cp *.txt guides
cd guides

# add the filename to each file for future use
perl -p -i -e 's/$/ $ARGV/;' *
cat *.txt > ctcf_all_targets_ucsc_guides.txt
```

This will create files for each region with the predicted UCSC guides and the filename for parsing.


## Determining CTCF motifs in the human genome and identifying overlaps with CTCF sites of interest

Steps for determining CTCF motifs in the human genome:
  - Create a motif file from UCSC galaxy human CTCF pssm file
  - Run scanMotifGenomeWide on hg19 with the custom motif file
  - Intersect genomic regions of original CTCF sites of interest with CTCF motif genomic regions

 ```bash
# set path so you can run Homer tools
PATH=$PATH:/home/t.severson/tools/homer/.//bin/
# Run Homer scanMotifGenomeWide to obtain CTCF motif regions
scanMotifGenomeWide.pl /home/t.severson/tools/homer/motifs/consensus_CTCF_galaxy.motif hg19 -bed > galaxy_ctcf_hg19.bed

# intersect with original CTCF binding sites
intersectBed -a ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only.bed -b galaxy_ctcf_hg19.bed -wa -wb > ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only_ctcf_galaxy_motifs.txt
```

## Combining the CTCF motif information with the guides designed in CTCF sites of interest

Steps for final guide selection:
  - Use a custom R script to make ouptut of UCSC file useful.
  - Filter guides based on uniqueness in the genome and (mit_specificity_score > 70 & doench2016_percentile_only > 65 ) with custom R script
  - Get coverage of motif regions with filtered guides
  - When possible, select the top 4 guides based on coverage of the CTCF motif (highest)


```bash
# get coverage of the regions of all CTCF motifs identified with the regions of the filtered CTCF guides
coverageBed -a ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only_ctcf_galaxy_motifs_all.bed -b guides/ctcf_ucsc_guides_filtered.bed > ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only_ctcf_galaxy_motifs_all_guides_coverageBed.txt

```
