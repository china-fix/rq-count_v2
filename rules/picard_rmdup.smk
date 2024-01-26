
rule view_cov:
    input:
        sorted_bam_file = "out_mapping/{sample}_{ref}.bam.sorted"
    output:
        #coverage_count = "report/{sample}_{ref}_rmdup.tab",
        tab_file = "report/{sample}_{ref}_depth.tab"
    conda:
        os.path.join(
            relative_dir, "envs/mapping.yaml"
        )
    shell:
        #"bedtools genomecov -ibam {input.sorted_bam_file_rmdup} -d > {output.coverage_count};"
        #"sed -i 's/$/\t1/g' {output.coverage_count};"
        #"sed -i '1i reference_seq\t1_index\tcoverage\tcount' {output.coverage_count};"
        "samtools depth -a -d 0 {input.sorted_bam_file} > {output.tab_file} "

