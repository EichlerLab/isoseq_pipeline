shell.prefix("source env.cfg; ")

configfile: "config.yaml"

import os

BAX2BAM  = "/net/eichler/vol18/zevk/great_apes/iso_seq/cc2_analysis/pitchfork/deployment/bin/bax2bam"
TOPGROUP = ["adult_brain"]
TYPES    = ["flnc", "nfl"]

GMAP_DB = config["gmap_db"]
GMAP_NAME = config["gmap_name"]

NAMES  = []
LOOKUP = {}

for i in TOPGROUP:
    for j in config[i].keys():
        LOOKUP[j] = i
        NAMES.append(j)

if not os.path.exists("log"):
    os.makedirs("log")

def _get_files_by_name(wildcards):
    return config[LOOKUP[wildcards.names]][wildcards.names]

localrules: all

rule all:
    input: expand("{top}/0to100kb_part0/{top}/combined/all_sizes.quivered_hq.fastq", top=TOPGROUP),

rule ice_clustering:
    input: flnc = "merged_trimmed/{top}/isoseq_flnc.fastq",
           nfl = "merged_trimmed/{top}/isoseq_nfl.fastq",
           bas_fofn = "ice_clustering/{top}/{top}.input.fofn",
    output: "{top}/0to100kb_part0/{top}/combined/all_sizes.quivered_hq.fastq"
    params: sge_opts = "-l mfree=4G", max_jobs=20, blasr_procs=4, gcon_procs=2, output_basename="final_consensus_sequence.{top}.fasta"
    shell: # Full path to actual tofu_wrap script: /net/eichler/vol8/home/zevk/projects/VENV_TOFU/lib/python2.7/site-packages/pbtools.pbtranscript-2.2.3-py2.7-linux-x86_64.egg/EGG-INFO/scripts/tofu_wrap.py
        'tofu_wrap.py --nfl_fa {input.nfl} --bas_fofn {input.bas_fofn} -d {wildcards.top} --use_sge \
                      --max_sge_jobs {params.max_jobs} --blasr_nproc {params.blasr_procs} \
                      --gcon_nproc {params.gcon_procs} --quiver_nproc 2 --quiver \
                      --bin_manual "(0,100)" --output_seqid_prefix {wildcards.top} \
                      --sge_env_name serial --gmap_db {GMAP_DB} --gmap_name {GMAP_NAME} \
                      {input.flnc} {params.output_basename}'

rule fofn   :
     input  : expand("merged_trimmed/{{top}}/isoseq_{type}.fastq", type=TYPES)
     output : "ice_clustering/{top}/{top}.input.fofn"
     params : sge_opts="-l mfree=2G -l h_rt=00:20:00 -q eichler-short.q"
     run  : 
        with open(output[0], "w") as out_fn:
            for tg in TOPGROUP:
                for group in config[tg]:
                    print(config[tg][group], file=out_fn)

rule faTofq :
     input  : READS="merged_trimmed/{top}/isoseq_{type}.fasta"
     output : "merged_trimmed/{top}/isoseq_{type}.fastq"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  : 
        "python scripts/fa2fq.py {input.READS}"

rule catlens :
     input  : STATS=expand("primer_trimmed/{names}/isoseq_flnc.fasta.fai", names=NAMES)
     output : "merged_trimmed/{top}/isoseq_{type}.fasta.fai"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
     shell  : "cat {input.STATS} >{output}"

rule catfa  :
     input  : FA=expand("primer_trimmed/{names}/isoseq_{{type}}.fasta", names=NAMES), STATS=expand("primer_trimmed/{names}/isoseq_{{type}}.fasta.len.txt", names=NAMES)
     output : "merged_trimmed/{top}/isoseq_{type}.fasta"
     params :  sge_opts="-l mfree=1G -l h_rt=02:00:00 -q eichler-short.q"
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
     params :  sge_opts="-l mfree=20G -l h_rt=4:00:00 -q eichler-short.q -pe serial 12"
     shadow : True
     shell  : 
        "python scripts/fluid_barcode_identification.py --reads_fn {input.FAS} "
                "--primer_fn_forward data/custom_barcode_primers_forward.fa "
                "--primer_fn_reverse data/custom_barcode_primers_reverse.fa "
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
