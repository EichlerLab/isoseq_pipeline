shell.prefix("source config.sh;")

configfile: "config.yaml"

import os

BAX2BAM  = "/net/eichler/vol18/zevk/great_apes/iso_seq/cc2_analysis/pitchfork/deployment/bin/bax2bam"
TOPGROUP = ["adult_brain"]
TYPES    = ["flnc", "nfl"]

NAMES  = []
LOOKUP = {}

for i in config.keys():
    for j in config[i].keys():
        LOOKUP[j] = i
        NAMES.append(j)

if not os.path.exists("log"):
    os.makedirs("log")

def _get_files_by_name(wildcards):
    print(wildcards.names, LOOKUP[wildcards.names], "\t", config[LOOKUP[wildcards.names]][wildcards.names])
    return config[LOOKUP[wildcards.names]][wildcards.names]

localrules: all

rule all  :
     input  : expand("ice_clustering/{top}/{top}.input.fofn", top=TOPGROUP),
              expand("merged_trimmed/{top}/isoseq_{type}.fasta",type=TYPES,top=TOPGROUP),
              expand("merged_trimmed/{top}/isoseq_flnc.fasta.fai",top=TOPGROUP)

rule fofn   :
     input  : expand("ice_clustering/{{top}}/{{top}}_{type}.fastq", type=TYPES)
     output : "ice_clustering/{top}/{top}.input.fofn"
     params : sge_opts="-l mfree=2G -l h_rt=00:20:00 -q eichler-short.q"
     shell  : 'grep {wildcards.top} config.json | perl helper_scripts/fofn-gen.pl  > {output}'

rule faTofq :
     input  : READS="merged_trimmed/{top}/isoseq_{type}.fasta"
     output : "ice_clustering/{top}/{top}_{type}.fastq"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  : 
        """cd ice_clustering/{wildcards.top}
           ln -s ../../{input.READS} {wildcards.top}_{wildcards.type}.fasta
           scripts/fa2fq.py {wildcards.top}_{wildcards.type}.fasta"""

rule catlens :
     input  : STATS=expand("primer_trimmed/{names}/isoseq_flnc.fasta.fai", names=NAMES)
     output : "merged_trimmed/{top}/isoseq_{type}.fasta.fai"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  : "cat {input.STATS} >{output}"

rule catfa  :
     input  : FA=expand("primer_trimmed/{names}/isoseq_{{type}}.fasta", names=NAMES), STATS=expand("primer_trimmed/{names}/isoseq_{{type}}.fasta.len.txt", names=NAMES)
     output : "merged_trimmed/{top}/isoseq_{type}.fasta"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     #shell  : "cat primer_trimmed/{wildcards.top}*/isoseq_{wildcards.type}.fasta > {output}"
     shell : "cat {input.FA} >{output}"

rule lenStat:
     input  : "primer_trimmed/{names}/isoseq_{type}.fasta.fai"
     output : "primer_trimmed/{names}/isoseq_{type}.fasta.len.txt"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  : "awk '{{print $2}}' {input} > {output}"

rule lens   :
     input  : "primer_trimmed/{names}/isoseq_{type}.fasta"
     output : "primer_trimmed/{names}/isoseq_{type}.fasta.fai"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  :  "samtools faidx {input}"

rule trim   :
     input  : FAS="post_pbccs_fasta/{names}.fa", BCC="scripts/do_barcode.sh"
     output : "primer_trimmed/{names}/isoseq_flnc.fasta", "primer_trimmed/{names}/isoseq_nfl.fasta"
     params :  sge_opts="-l mfree=20G -l h_rt=4:00:00 -q eichler-short.q -V -cwd"
     shadow : True
     shell  : 
        "python scripts/fluid_barcode_identification.py --reads_fn {input.FAS} "
                "--primer_fn_forward data/custom_barcode_primers_forward.fa "
                "--primer_fn_reverse custom_barcode_primers_reverse.fa "
                "--ncpus 12 "
                "--flnc_fn_out {output[0]} "
                "--nfl_fn_out {output[1]}"

rule bam2fa :
     input  : BAM="pbccs_results/{names}.pbccs.bam", XML="pbccs_results/{names}.pbccs.consensusreadset.xml"
     output : "post_pbccs_fasta/{names}.fa"
     params :  sge_opts="-l mfree=3G -l h_rt=48:00:00 -q eichler-short.q"
     shell  : "bamtools convert -in {input.BAM} -format fasta > {output}"

rule ccs    :
     input  : BAM="cc2_bams/{names}.subreads.bam", CCS="/net/eichler/vol18/zevk/great_apes/iso_seq/cc2_analysis/pitchfork/deployment/bin/ccs"
     output : "pbccs_results/{names}.pbccs.bam", "pbccs_results/{names}.pbccs.consensusreadset.xml"
     params :  sge_opts="-l mfree=4G -l h_rt=48:00:00 -q eichler-short.q -pe serial 4"
     shell  : "{input.CCS} --numThreads=4 --minLength=200 {input.BAM} {output[0]}"

rule bax2bam:
     input  : BAX2BAM , FL=_get_files_by_name
     output : "cc2_bams/{names}.subreads.bam"
     params : sge_opts="-l mfree=15G -l h_rt=06:00:00 -q eichler-short.q"
     shell  : "{BAX2BAM} -o cc2_bams/{wildcards.names} {input.FL}"
