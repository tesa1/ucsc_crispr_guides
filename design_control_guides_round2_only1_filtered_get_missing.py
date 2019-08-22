
# From a .tsv input file of regions design guides using UCSC crispr.bb file

import os
import glob
import pandas as pd

WDIR="/home/t.severson/zwart/crispr_guides/ctcf_ucsc_guides/test/"

SAMPLES = pd.read_table(WDIR + "target_guides_need_controls_round2.tsv").set_index("control_site", drop=False)

rule all:
    input: WDIR + "target_guides_need_controls_round3.tsv"
    
rule run_control_shell_script:
    input:
        sh=WDIR + "target_guides_need_controls_round2.sh"
    output: 
        gu=expand(WDIR + "{sample}_ucsc_control_guides_round2.txt", sample=SAMPLES.index)
    message:
        'Designing Guides with UCSC crispr.bb'
    shell:
        "sh {input.sh} > {output.gu}"


rule run_r_filter_control_script:
    input:
        i=expand(WDIR + "{sample}_control_ucsc_guides.txt", sample=SAMPLES.index),
        t=WDIR + "filter_control_guides_round2.R"
    output:
        WDIR + "target_guides_need_guides_round2.tsv"
    shell:
        "Rscript {input.t}"

