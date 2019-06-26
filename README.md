# Designing guides for a CRISPR-Cas9 experiment to knockout CTCF sites

This is a manual of how to use publically available tools to obtain previously identified guides and parse for those
that should work.

Steps for designing guides:
  - Search UCSC database of all known guides 
  - Filter out guides which are not unique in the genome
  - Filter for those that have the best possibility of working
  - Add CTCF motif information and select guide which is in CTCF site if possible
  
 ## Use of UCSC All predicted Super track ##
Information can be found here: https://genome-test.gi.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chr2&g=crisprAll

Data can be extracted with the kent command line tools (bigBedToBed), available:
   http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
   
Using the bed file of 720 CTCF sites, create a shell script to extract the data. Then run the script to get the CRISPR All guides
in the regions.


 ```bash
 # get the guides
sh ctcf_sites_ucsc_guides.sh 
mkdir guides
cp *.txt guides
cd guides

# add the filename to each file for future use
perl -p -i -e 's/$/ $ARGV/;' *
```

This will create files for each region with the predicted UCSC guides and the filename for parsing.


## Determining CTCF motifs in the human genome

 ```bash
# set path so you can run Homer tools
PATH=$PATH:/home/t.severson/tools/homer/.//bin/
# Run Homer scanMotifGenomeWide to obtain CTCF motif regions
scanMotifGenomeWide.pl /home/t.severson/tools/homer/motifs/consensus_CTCF_galaxy.motif hg19 -bed > galaxy_ctcf_hg19.bed

# intersect with original CTCF binding sites
intersectBed -a ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only.bed -b galaxy_ctcf_hg19.bed -wa -wb > ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only_ctcf_galaxy_motifs.txt

```
