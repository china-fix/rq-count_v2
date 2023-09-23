

### apply bwa to map the clean reads to referece seqs ###
rule bwa:
    input:
        read_1 = "out_clean_read/{sample}1.fq.gz",
        read_2 = "out_clean_read/{sample}2.fq.gz",
        ref = "in/REF/{ref}.REF"
    output:
        sam_file = "out_mapping/{sample}_{ref}.sam"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
        # "../envs/mapping.yaml"
    threads: config["general_config"]["threads_n"]
    shell:
        "bwa index {input.ref};"
        "bwa mem -t {threads} {input.ref} {input.read_1} {input.read_2}  > {output.sam_file};"
        # "bwa mem -A 1 -B 8 -w 20 -T 120 -t {threads} {input.ref} {input.read_1} {input.read_2}  > {output.sam_file};"
        #"rm -f {input.ref}.amb {input.ref}.ann {input.ref}.bwt {input.ref}.pac {input.ref}.sa;"

### fitering the bwa alignment according alignment score (AS) ###
#grep -v -P 'AS:i:(1[0-1][0-9]|[1-9][0-9]|[0-9])\t' out_mapping/1_FDSW210197737-1r__REF2.sam > ref2.sam


### change the sam to bam and sorted###
rule samtools:
    input:
        sam_file = "out_mapping/{sample}_{ref}.sam"
    output:
        bam_file = "out_mapping/{sample}_{ref}.bam",
        sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
        # "../envs/mapping.yaml"
    threads: config["general_config"]["threads_n"]
    shell:
        "samtools view -b --threads {threads} {input.sam_file} > {output.bam_file};"
        "samtools sort -o {output.sorted_bam_file} --threads {threads} {output.bam_file} "

# ### get the mapping coverage ###
# rule bedtools:
#     input:
#         sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
#     output:
#         coverage_count = "report/{sample}_{ref}.tab"
#     conda:
#         "../envs/mapping.yaml"
#     shell:
#         "bedtools genomecov -ibam {input.sorted_bam_file} -d > {output.coverage_count};"
#         "sed -i 's/$/\t1/g' {output.coverage_count};"
#         "sed -i '1i reference_seq\t1_index\tcoverage\tcount' {output.coverage_count}"

### get the mapping coverage ###
### get the tab file describe the mapping depth of every base
#rule bedtools:
#    input:
#        sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
#    output:
#        coverage_count = "report/{sample}_{ref}.tab",
#        tab_file = "report/{sample}_{ref}_depth.tab"
#    conda:
#        os.path.join(
#            relative_dir, "envs/mapping.yaml"
#        )
#        # "../envs/mapping.yaml"
#    shell:
#        "bedtools genomecov -ibam {input.sorted_bam_file} -d > {output.coverage_count};"
#        "sed -i 's/$/\t1/g' {output.coverage_count};"
#        "sed -i '1i reference_seq\t1_index\tcoverage\tcount' {output.coverage_count};"
#        "samtools depth -a {input.sorted_bam_file} > {output.tab_file} "

