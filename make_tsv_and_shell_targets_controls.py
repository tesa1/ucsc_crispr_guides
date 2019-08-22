
# From a bed file, make a shell script and a tsv file to run the guide design 
import os
import glob

WDIR="/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/pipeline/"
# INBED="ctcf_kmeans_co_h3k27ac_up_in_mets_tads_coverage_only.bed" = "input_sites.bed"
INBED="input_sites.bed"

TARGET_SHSCRIPT="target_guides_from_input_sites_bed.sh"
CONTROL_SHSCRIPT="control_guides_from_input_sites_bed.sh"

TARGET_RSCRIPT="design_ucsc_target_guides_from_input_sites_bed.R"
CONTROL_RSCRIPT="design_ucsc_control_guides_from_input_sites_bed.R"

rule all:
    input: WDIR + TARGET_SHSCRIPT, WDIR + CONTROL_SHSCRIPT
    
rule make_shell_script:
    input:
        be=WDIR + INBED,
        sc=WDIR + TARGET_RSCRIPT
    output: 
        WDIR + TARGET_SHSCRIPT
    shell:
        "Rscript {input.sc}"
        
        
rule make_controls_shell_script:
    input:
        be=WDIR + INBED,
        se=WDIR + CONTROL_RSCRIPT
    output: 
        WDIR + CONTROL_SHSCRIPT
    shell:
        "Rscript {input.se}"
        

        
