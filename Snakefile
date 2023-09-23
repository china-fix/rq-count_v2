import os

relative_dir = os.path.dirname(os.path.abspath(workflow.snakefile))

configfile: relative_dir+"/config/config.yaml"
# print(config["general_config"]["threads_n"])

#SAMPLES, = glob_wildcards("in/{sample}1.(fq|fastq).gz")
#SAMPLES, = glob_wildcards("in/{sample}{read_ext}", read_ext=["1.fq.gz", "1.fastq.gz"])
REFS, = glob_wildcards("in/REF/{ref}.REF")


# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES,EXTENSIONS, = glob_wildcards("in/{sample}1.{read_ext}.gz")
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]



include: "rules/fastp.smk"
include: "rules/mapping.smk"
include: "rules/picard_rmdup.smk"



rule all:
    input:
        #fastp
        # expand("out_clean_read/{sample}1.fq.gz", sample=SAMPLES),
        # expand("out_clean_read/{sample}2.fq.gz", sample=SAMPLES),
        #bwa
        # expand("out_mapping/{sample}_{ref}.sam", sample=SAMPLES, ref=REFS),
        #samtools
        # expand("out_mapping/{sample}_{ref}.bam", sample=SAMPLES, ref=REFS),
        # expand("out_mapping/{sample}_{ref}.bam.sorted", sample=SAMPLES, ref=REFS),
        #bedtools
        #expand("report/{sample}_{ref}.tab", sample=SAMPLES, ref=REFS),
        #expand("report/{sample}_{ref}_rmdup.tab", sample=SAMPLES, ref=REFS),
        expand("report/{sample}_{ref}_depth_rmdup.tab", sample=SAMPLES, ref=REFS),
