configfile: "config.json"

rule get_snp_verification:
    output:
        touch("get_snp.done")
    params:
        repo = config["SNP"]["REPO"]
    conda:
        config["BUILD"]["ENV"]
    shell:
        "git -C bin/ clone {params.repo};"

rule add_header:
    input:
        OUTDIR + "ResistomeResults/AMR_analytic_matrix.csv"
    output:
        temp(OUTDIR + "ResistomeResults/AMR_analytic_matrix_updated_header.csv")
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/add_header.sh {input} {output}"

rule pass_config_file:
    input:
        OUTDIR + "ResistomeResults/AMR_analytic_matrix_updated_header.csv",
        SAM
    output:
        temp(OUTDIR + "{sample}.config.ini")
    params:
        sample_id = "{sample}"
    run:
        import configparser
        with open(output[0], 'w') as configfile_out:
            config_to_pass = dict(config["SNP"]["CONFIG"])
            config_to_pass["SOURCE_FILES"]["COUNT_MATRIX"] = input[0]
            config_to_pass["SOURCE_FILES"]["SAM_INPUT"] = input[1]
            config_to_pass["OUTPUT_FILES"]["COUNT_MATRIX_FINAL"] = OUTDIR + "SNP/" + params.sample_id + "_AMR_analytic_matrix_with_SNP_confirmation.csv"
            config_to_pass["OUTPUT_FILES"]["OUTPUT_FOLDER"] = OUTDIR + "SNP/" + params.sample_id + "/"
            config_to_pass["TEMP_FILES"]["TEMP_SAM_SORTED"] = "temp/" + params.sample_id + ".sorted.sam"
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)

rule run_snp:
    input:
        "get_snp.done",
        OUTDIR + "ResistomeResults/AMR_analytic_matrix_updated_header.csv",
        sam = SAM,
        config_file = OUTDIR + "{sample}.config.ini"
    output:
        temp(OUTDIR + "SNP/{sample}_AMR_analytic_matrix_with_SNP_confirmation.csv"),
        OUTDIR + "SNP/{sample}/Normal_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Frameshift_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Hypersusceptible_Mutations_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Suppressible_Frameshift_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Intrinsic_Resistance_Genes.csv"
    params:
        outdir = OUTDIR + "SNP",
        config_file = OUTDIR + "{sample}.config_01.ini",
        sorted_sam = "temp/{sample}.sorted.sam"
    conda:
        config["WORKFLOW"]["ENV"]
    envmodules:
        "python/3.8"
    shell:
        "python3 bin/AmrPlusPlus_SNP/SNP_Verification.py "
        "-c {input.config_file} "
        "-a -i {input.sam} "
        "-o {params.outdir}; "
        "rm {params.config_file}; "
        "rm {params.sorted_sam}"

rule merge_snp_matricies:
    input:
        snp_out = expand(OUTDIR + "SNP/{sample}_AMR_analytic_matrix_with_SNP_confirmation.csv", sample = SAMPLES),
        amr_csv = OUTDIR + "ResistomeResults/AMR_analytic_matrix_updated_header.csv"
    output:
        OUTDIR + "ResistomeResults/AMR_analytic_matrix_with_SNP_confirmation.csv"
    shell:
        "rm -rf temp/; "
        "bin/merge_snp_out.py "
        "--infiles {input.snp_out} "
        "--amr_csv {input.amr_csv} "
        "--outfile {output}"


## WRITE RULE THAT COMBINES THE CHANGES MATRACIES