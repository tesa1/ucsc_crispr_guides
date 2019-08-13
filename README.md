# Designing guides for a CRISPR-Cas9 experiment to knockout CTCF sites (including control guides)

This is a manual of how to use publically available tools to obtain previously identified guides and parse for those
that should work.

Steps for designing guides:
  - Search UCSC database of all known guides 
  - Filter out guides which are not unique in the genome
  - Filter for those that have the best possibility of working
  - Select guides which are in CTCF motif if possible
  
Steps for designing control guides:
  - Search UCSC database of all known guides 
  - Filter out guides which are not unique in the genome
  - Filter for those that have the best possibility of working

  
 ## Use of UCSC All predicted Super track ##
Information can be found here: # link is dead https://genome-test.gi.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chr2&g=crisprAll 
and here:
https://genome-euro.ucsc.edu/cgi-bin/hgTables

It's not necessary but data can be downloaded from https://hgdownload-test.gi.ucsc.edu/gbdb/hg19/crisprAll/crispr.bb.
And it can queried in the UCSC table browser under Genes and Gene Predictions -> CRISPR Targets (https://genome.ucsc.edu/cgi-bin/hgTables) file queried is:/gbdb/hg19/crisprAll/crispr.bb

Downloaded data can be extracted with the kent command line tools (bigBedToBed), available:
   http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
   
Use the snakemake file to query the bed file of 720 CTCF sites and design a shell script to obtain guides in the target regions and the control regions (250-500bp upstream of target CTCF site)


 ```bash
snakemake -s make_tsv_and_shell_control.py 
# note this snakemake requires two files 1) design_ucsc_guides_from_ctcf_sites_bed.R and 2) design_ucsc_control_guides_from_ucsc_sites_bed.R files
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
