configfile: "config.json"

rule get_minikraken:
    output:
        touch("get_minikraken.done")
    shell:
        "bin/download_minikraken.sh"


rule run_kraken:
    input:
        "get_minikraken.done",
        f_reads = OUTDIR + "NonHostReads/{sample}.non.host.R1.fastq.gz",
        r_reads = OUTDIR + "NonHostReads/{sample}.non.host.R2.fastq.gz"
    output:
        report = OUTDIR + "RunKraken/Standard_report/{sample}.kraken.report",
        raw = OUTDIR + "RunKraken/Standard/{sample}.kraken.raw",
        filtered_report = OUTDIR + "RunKraken/Filtered_report/{sample}.kraken.filtered.report",
        filtered_raw = OUTDIR + "RunKraken/Filtered/{sample}.kraken.filtered.raw"
    conda:
        config["KRAKEN"]["ENV"]
    params:
        kraken_db = "--db " + config["KRAKEN"]["DB"],
        confidence = "--confidence 1"
    threads:
        config["KRAKEN"]["THREADS"]
    shell:
        "kraken2 {params.kraken_db} "
        "--paried {input.f_reads} {input.r_reads} "
        "--threads {threads} "
        "--report {output.report} > {output.raw}; "
        "kraken2 {params.kraken_db} "
        "{params.confidence} "
        "--paried {input.f_reads} {input.r_reads} "
        "--threads {threads} "
        "--report {output.filtered_report} > {output.filtered_raw}"


rule kraken_results:
    input:
        expand(OUTDIR + "RunKraken/Standard_report/{sample}.kraken.report", sample = SAMPLES)
    output:
        OUTDIR + "KrakenResults/kraken_analytic_matrix.csv"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/kraken2_long_to_wide.py -i {input} -o {output}"


rule filtered_kraken_results:
    input:
        expand(OUTDIR + "RunKraken/Filtered_report/{sample}.kraken.filtered.report", sample = SAMPLES)
    output:
        OUTDIR + "FilteredKrakenResults/filtered_kraken_analytic_matrix.csv"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/kraken2_long_to_wide.py -i {input} -o {output}"