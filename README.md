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

All shell scripts produced assume bigBedToBed is here:/home/t.severson/tools/kent_utils/bigBedToBed and the crispr.bb file is here: /home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/crispr.bb. You will have to change these locations in the Rscripts to your relevant locations in order to run the shell scripts.


## Query the bed file of target sites and design a shell script to obtain guides in the target regions and the control regions (250-500bp upstream of target site)

 ```bash
snakemake -s make_tsv_and_shell_targets_controls.py 
# note this snakemake requires two R scripts: 1) design_ucsc_target_guides_from_input_sites_bed.R and 
# 2) design_ucsc_control_guides_from_input_sites_bed.R 
```

## Obtain UCSC guides in target regions and and filter for uniqueness, specificity, effeciency and motif coverage 
After initial filtering, this script will identify the top 4 guides per target region based on motif coverage percentage.

homer genomeWideMotifScan needs to be installed (eg. /home/t.severson/tools/homer/) and a motif file is necessary with the motif path and motif file defined at the top of design_target_guides_homer_motif_only4_filtered.py. Homer home needs to be defined on line 48 of design_target_guides_homer_motif_only4_filtered.py
(http://homer.ucsd.edu/homer/motif/genomeWideMotifScan.html)

```bash
snakemake -s design_target_guides_homer_motif_only4_filtered.py
# note this snakemake file requires three R scripts 1) filter_ucsc_target_guides.R, 2) make_homer_motif_info_bed.R and
# 3) filter_target_guides_with_motif_info.R.
```

## Obtain UCSC guides in control regions for targets and filter for uniqueness, specificity and effeciency. Also determine regions that do not yet have target or control guides and design new shell scripts (for new controls, the region is moved 250bp upstream from original)
After initial filtering, this script will identify the top 1 control guides per target region based on Doench2016 percentile. Provide relevant details for crispr.bb and bigBedToBed.

```bash
snakemake -s design_control_guides_only1_filtered_get_missing.py
# note this snakemake file requires one R script 1) filter_control_guides.R 
```

The files with relevant CRISPR-Cas9 guide sequence information are (targets will have 4 guides and controls 1 per target): 
  target_guides_filtered_homer_motif_coverage_info_only4_filtered_final_list_targets.csv
  target_guides_filtered_homer_motif_coverage_info_only4_filtered_final_list_controls.csv

