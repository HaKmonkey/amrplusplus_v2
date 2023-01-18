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
        #"chmod -R +x data/AmrPlusPlus_SNP/"


rule pass_config_file:
    input:
        OUTDIR + "ResistomeResults/AMR_analytic_matrix.csv",
        OUTDIR + "{sample}.amr.alignment.sam"
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

            # config_to_pass["OUTPUT_FILES"]["COUNT_MATRIX_FINAL"] = OUTDIR + "SNP/" + params.sample_id + "_AMR_analytic_matrix_with_SNP_confirmation.csv"
            # config_to_pass["OUTPUT_FILES"]["OUTPUT_FOLDER"] = OUTDIR + "SNP/" + params.sample_id + "/"
            # config_to_pass["OUTPUT_FILES"]["NORMAL_TYPE_OUTPUT"] = OUTDIR + "SNP/" + params.sample_id + "_Normal_Type_Genes.csv"
            # config_to_pass["OUTPUT_FILES"]["FRAMESHIFT_TYPE_OUTPUT"] = OUTDIR + "SNP/" + params.sample_id + "_Frameshift_Type_Genes.csv"
            # config_to_pass["OUTPUT_FILES"]["HYPERSUSCEPTIBLE_TYPE_OUTPUT"] = OUTDIR + "SNP/" + params.sample_id + "_Hypersusceptible_Mutations_Type_Genes.csv"
            # config_to_pass["OUTPUT_FILES"]["SUPPRESSIBLE_TYPE_OUTPUT"] = OUTDIR + "SNP/" + params.sample_id + "_Suppressible_Frameshift_Type_Genes.csv"
            # config_to_pass["OUTPUT_FILES"]["INTRINSIC_TYPE_OUTPUT"] = OUTDIR + "SNP/" + params.sample_id + "_Intrinsic_Resistance_Genes.csv"
            # config_to_pass["OUTPUT_FILES"]["DETAILED_FOLDER"] = OUTDIR + "SNP/" + params.sample_id + "_detailed_output"

            # config_to_pass["TEMP_FILES"]["TEMP_SAM_SORTED"] = OUTDIR + "SNP/temp/sorted.sam"

            config_to_pass["OUTPUT_FILES"]["COUNT_MATRIX_FINAL"] = OUTDIR + "SNP/" + params.sample_id + "_AMR_analytic_matrix_with_SNP_confirmation.csv"
            config_to_pass["OUTPUT_FILES"]["OUTPUT_FOLDER"] = OUTDIR + "SNP/" + params.sample_id + "/"
            # config_to_pass["OUTPUT_FILES"]["OUTPUT_FOLDER"] = OUTDIR + "SNP/"

            config_to_pass["TEMP_FILES"]["TEMP_SAM_SORTED"] = "temp/sorted.sam"
            
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)


rule run_snp:
    input:
        "get_snp.done",
        OUTDIR + "ResistomeResults/AMR_analytic_matrix.csv",
        sam = OUTDIR + "{sample}.amr.alignment.sam",
        config_file = OUTDIR + "{sample}.config.ini"
    output:
        # OUTDIR + "SNP/{sample}_AMR_analytic_matrix_with_SNP_confirmation.csv",
        # OUTDIR + "SNP/{sample}.amr.alignment.dedup",
        OUTDIR + "SNP/{sample}/Normal_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Frameshift_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Hypersusceptible_Mutations_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Suppressible_Frameshift_Type_Genes.csv",
        OUTDIR + "SNP/{sample}/Intrinsic_Resistance_Genes.csv"
        #OUTDIR + "SNP/{sample}/detailed_output"
        #OUTDIR + "SNP/temp/sorted.sam"
    params:
        outdir = OUTDIR + "SNP/{sample}",
        config_file = OUTDIR + "{sample}.config_01.ini"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "python3 bin//AmrPlusPlus_SNP/SNP_Verification.py -c {input.config_file}; "
        "rm {params.config_file}; "
        "rm -rf temp/"
        # "-a -i {input.sam} "
        # "-o {params.outdir}"
        # "cut -d ',' -f `awk -v RS=',' "/${sample_id}/{print NR; exit}" ${snp_count_matrix}` ${snp_count_matrix} > ${sample_id}_${prefix}_SNP_count_col"


