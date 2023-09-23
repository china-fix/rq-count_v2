# ### change the sam to bam and sorted###
# rule samtools:
#     input:
#         sam_file = "out_mapping/{sample}_{ref}.sam"
#     output:
#         bam_file = "out_mapping/{sample}_{ref}.bam",
#         sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
#     conda:
#         os.path.join(
#             relative_dir, "envs/mapping.yaml"
#         )
#         # "../envs/mapping.yaml"
#     threads: config["general_config"]["threads_n"]
#     shell:
#         "samtools view -b --threads {threads} {input.sam_file} > {output.bam_file};"
#         "samtools sort -o {output.sorted_bam_file} --threads {threads} {output.bam_file} "


### input the sorted bam and remove the sequencing dupplicates
rule picard_rm_dup:
    input:
        sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
    output:
        sorted_bam_file_rmdup = "out_mapping/{sample}_{ref}.bam.sorted.rmdup",
        metrict_file = "out_mapping/{sample}_{ref}.metrics"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
    threads: config["general_config"]["threads_n"]
    shell:
        "picard MarkDuplicates -I {input.sorted_bam_file} -O {output.sorted_bam_file_rmdup} -M {output.metrict_file} --REMOVE_SEQUENCING_DUPLICATES true"


rule view_cov:
    input:
        sorted_bam_file_rmdup = "out_mapping/{sample}_{ref}.bam.sorted.rmdup"
    output:
        #coverage_count = "report/{sample}_{ref}_rmdup.tab",
        tab_file = "report/{sample}_{ref}_depth_rmdup.tab"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
    shell:
        #"bedtools genomecov -ibam {input.sorted_bam_file_rmdup} -d > {output.coverage_count};"
        #"sed -i 's/$/\t1/g' {output.coverage_count};"
        #"sed -i '1i reference_seq\t1_index\tcoverage\tcount' {output.coverage_count};"
        "samtools depth -a {input.sorted_bam_file_rmdup} > {output.tab_file} "

