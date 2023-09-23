### all-in-one FASTQ preprocessor ###

rule fastp:
    input:
        #read_1 = "in/{sample}1.(fq|fastq).gz",
        #read_2 = "in/{sample}2.(fq|fastq).gz",
        read_1 = "in/{sample}1."+FQEXTN+".gz",
        read_2 = "in/{sample}2."+FQEXTN+".gz",
    output:
        read_1 = "out_clean_read/{sample}1.fq.gz",
        read_2 = "out_clean_read/{sample}2.fq.gz",
        report_html = "out_clean_read/{sample}.html",
        report_json = "out_clean_read/{sample}.json",
    conda:
        os.path.join(
            relative_dir, "envs/fastp.yaml"
        )
        # "../envs/fastp.yaml"
    threads: config["general_config"]["threads_n"]
    shell:
        "fastp --in1 {input.read_1} --in2 {input.read_2} --out1 {output.read_1} --out2 {output.read_2} \
        -c -h {output.report_html} -j {output.report_json} -w {threads}"

# rule abricate_sum:
#     input:
#         ncbi = "out_abricate/{sample}_ncbi.tab",
#         card = "out_abricate/{sample}_card.tab",
#         vfdb = "out_abricate/{sample}_vfdb.tab",
#         plasmidfinder = "out_abricate/{sample}_plasmidfinder.tab",
#         resfinder = "out_abricate/{sample}_resfinder.tab",
#     output:
#         ncbi = "out_abricate/sum/ncbi/{sample}",
#         card = "out_abricate/sum/card/{sample}",
#         vfdb = "out_abricate/sum/vfdb/{sample}",
#         plasmidfinder = "out_abricate/sum/plasmidfinder/{sample}",
#         resfinder = "out_abricate/sum/resfinder/{sample}",
#     shell:
#         "cp {input.ncbi} {output.ncbi};"
#         "cp {input.card} {output.card};"
#         "cp {input.vfdb} {output.vfdb};"
#         "cp {input.plasmidfinder} {output.plasmidfinder};"
#         "cp {input.resfinder} {output.resfinder}"
#         # in the later summary report phase use the '''abricate --summary results.tab > summary.tab'''