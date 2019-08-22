# From a .tsv input file of regions design guides using UCSC crispr.bb file

import os
import glob
import pandas as pd

WDIR="/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/"
MOTIF_PATH="/home/t.severson/tools/homer/motifs/"
MOTIF_FILE="consensus_CTCF_galaxy.motif"
# mv ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only.bed input_sites.bed
INBED="input_sites.bed"

TARGET_SHSCRIPT="target_guides_from_input_sites_bed.sh"

SAMPLES = pd.read_table(WDIR + "target_guides_from_input_sites_bed.tsv").set_index("site", drop=False)


rule all:
    input: WDIR + "target_guides_filtered.bed", WDIR + "homer_motif_hg19.bed", WDIR + "target_guides_filtered_homer_motifs.txt", WDIR + "target_guides_filtered_homer_motifs_all.bed", WDIR + "target_guides_filtered_homer_motifs_all_coverageBed.txt", WDIR + "target_guides_filtered_homer_motif_coverage_info_top4_filtered.pdf"
 
    
rule run_shell_target_script:
    input:
        sh=WDIR + TARGET_SHSCRIPT
    output: 
        gu=expand(WDIR + "{sample}_ucsc_guides.txt", sample=SAMPLES.index)
    message:
        'Designing Target Guides with UCSC crispr.bb'
    shell:
        "sh {input.sh} > {output.gu}"

rule run_r_filter_script:
    input: 
        ex=expand(WDIR + "{sample}_ucsc_guides.txt", sample=SAMPLES.index),
        sc=WDIR + "filter_ucsc_target_guides.R"
    output: 
        WDIR + "target_guides_filtered.bed"
    shell:
        "Rscript {input.sc}"
        
        
rule homer_motif:
    input:
        mo=MOTIF_PATH + MOTIF_FILE
    output: 
        WDIR + "homer_motif_hg19.bed"
    shell:
        "export PATH=$PATH:/home/t.severson/tools/homer/.//bin/ &&"
        "scanMotifGenomeWide.pl {input.mo} hg19 -bed > {output}"


rule intersect_target_motif:
    input:
        my=WDIR + INBED,
        ctcf=WDIR + "homer_motif_hg19.bed"
    output: 
        WDIR + "target_guides_filtered_homer_motifs.txt"
    shell:
        "intersectBed -a {input.my} -b {input.ctcf} -wa -wb > {output}"

        
rule make_homer_motif_bed:
    input:
        mo=WDIR + "target_guides_filtered_homer_motifs.txt",
        sc=WDIR + "make_homer_motif_info_bed.R"
    output:
        WDIR + "target_guides_filtered_homer_motifs_all.bed"
    shell:
        "Rscript {input.sc}"


rule motif_vs_filter_cov_bed:
    input:
        fi=WDIR + "target_guides_filtered.bed",
        co=WDIR + "target_guides_filtered_homer_motifs_all.bed"
    output:
        WDIR + "target_guides_filtered_homer_motifs_all_coverageBed.txt"
    shell:
        "coverageBed -a {input.fi} -b {input.co} > {output}"


rule final_4_target_guides:
    input:
        fi=WDIR + "target_guides_filtered_homer_motifs_all_coverageBed.txt",
        sc=WDIR + "filter_target_guides_with_motif_info.R"
    output:
        WDIR + "target_guides_filtered_homer_motif_coverage_info_top4_filtered.pdf"
    shell:
        "Rscript {input.sc}"
