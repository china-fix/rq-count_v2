
### pre_bwa for indexing the refs ###
rule pre_bwa:
    input:
        ref = "in/REF/{ref}.REF"
    output:
        ref_index = "out_REF_index/{ref}.REF"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
        # "../envs/mapping.yaml"
    threads: config["general_config"]["threads_n"]
    shell:
        "cp {input.ref} {output.ref_index};"
        "bwa index {output.ref_index};"
        # "bwa mem -t {threads} {input.ref} {input.read_1} {input.read_2}  > {output.sam_file};"



### apply bwa to map the clean reads to referece seqs ###
rule bwa:
    input:
        read_1 = "out_clean_read/{sample}1.fq.gz",
        # read_2 = "out_clean_read/{sample}2.fq.gz",
        ref = "out_REF_index/{ref}.REF"
    output:
        #sam_file = "out_mapping/{sample}_{ref}.sam"
        sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
        # "../envs/mapping.yaml"
    threads: config["general_config"]["threads_n"]
    shell:
        # "bwa index {input.ref};"
        "bwa mem -t {threads} {input.ref} {input.read_1} | samtools view -b --threads {threads} | samtools sort -o {output.sorted_bam_file} --threads {threads} ;"
        # "bwa mem -A 1 -B 8 -w 20 -T 120 -t {threads} {input.ref} {input.read_1} {input.read_2}  > {output.sam_file};"
        #"rm -f {input.ref}.amb {input.ref}.ann {input.ref}.bwt {input.ref}.pac {input.ref}.sa;"




