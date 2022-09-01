configfile: "config.json"

import sys
from os import listdir
from os.path import exists

# Input file setup and check ###################################################

OUTDIR = config["WORKFLOW"]["OUTPUT"]
HOST_FILE = config["BWA"]["HOST"]
AMR_FILE = config["BWA"]["AMR"]
ADAPTERS_FILE = config["TRIMMOMATIC"]["ADAPTERS"]
ANNOTATION_FILE = config["RESISTOME"]["ANNOTATION"]


if not exists(HOST_FILE):
    print(f'HOST failed to open {HOST_FILE} : No such file or directory')
    sys.exit()
if not exists(AMR_FILE):
    print(f'AMR failed to open {AMR_FILE} : No such file or directory')
    sys.exit()
if not exists(ADAPTERS_FILE):
    print(f'ADAPTERS failed to open {ADAPTERS_FILE} : No such file or directory')
    sys.exit()
if not exists(ANNOTATION_FILE):
    print(f'ANNOTATION failed to open {ANNOTATION_FILE} : No such file or directory')
    sys.exit()

# eq to params.reads ###########################################################

def get_samples(reads_source):
    samples = []
    for sample in listdir(reads_source):
        samples.append(sample.split("_")[0])
    return list(set(samples))

SAMPLES = get_samples(config["WORKFLOW"]["READS_SOURCE"])

# set scripts to correct architecture ##########################################

# if config["WORKFLOW"]["ARCHITECTURE"] == "arm64":
#     resistome_script = config["RESISTOME"]["ARM64_SCRIPT"]
#     rarefaction_script = config["RAREFACTION"]["ARM64_SCRIPT"]
# elif config["WORKFLOW"]["ARCHITECTURE"] == "x86":
#     resistome_script = config["RESISTOME"]["X86_SCRIPT"]
#     rarefaction_script = config["RAREFACTION"]["X86_SCRIPT"]

# setting up `all` input so pipeline runs without params #######################

all_input = [
    "build_rarefaction.done",
    "build_resistome.done",
    OUTDIR + "RunQC/trimmomatic.stats",
    OUTDIR + "RemoveHostDNA/HostRemovalStats/host.removal.stats",
    OUTDIR + "ResistomeResults/AMR_analytic_matrix.csv",
    OUTDIR + "SamDedup_ResistomeResults/SamDedup_AMR_analytic_matrix.csv",
    expand(OUTDIR + "RunRarefaction/{sample}.gene.tsv", sample = SAMPLES),
    expand(OUTDIR + "RunRarefaction/{sample}.group.tsv", sample = SAMPLES),
    expand(OUTDIR + "RunRarefaction/{sample}.mechanism.tsv", sample = SAMPLES),
    expand(OUTDIR + "RunRarefaction/{sample}.class.tsv", sample = SAMPLES)
]

if config["KRAKEN"]["INCLUDE"] == "true":
    kraken_results = [
        OUTDIR + "KrakenResults/kraken_analytic_matrix.csv",
        OUTDIR + "FilteredKrakenResults/filtered_kraken_analytic_matrix.csv"
    ]
    all_input.append(kraken_results)

################################################################################


rule all:
    input:
        all_input


rule build_resistome:
    output:
        touch("build_resistome.done")
    conda:
        config["BUILD"]["ENV"]
    shell:
        "bin/build_resistome.sh"


rule build_rarefaction:
    output:
        touch("build_rarefaction.done")
    conda:
        config["BUILD"]["ENV"]
    shell:
        "bin/build_rarefaction.sh"


rule run_qc:
    input:
        f_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_test_R1.fastq.gz",
        r_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_test_R2.fastq.gz"
    output:
        p1 = OUTDIR + "RunQC/Paired/{sample}.1P.fastq.gz",
        p2 = OUTDIR + "RunQC/Paired/{sample}.2P.fastq.gz",
        u1 = OUTDIR + "RunQC/Unpaired/{sample}.1U.fastq.gz",
        u2 = OUTDIR + "RunQC/Unpaired/{sample}.2U.fastq.gz"
    params:
        illumina_clip = "ILLUMINACLIP:" + ADAPTERS_FILE + ":2:30:10:3:TRUE",
        leading = "LEADING:" + config["TRIMMOMATIC"]["LEADING"],
        trailing = "TRAILING:" + config["TRIMMOMATIC"]["TRAILING"],
        sliding_window = "SLIDINGWINDOW:" + config["TRIMMOMATIC"]["SLIDING_WINDOW"],
        minlen = "MINLEN:" + config["TRIMMOMATIC"]["MINLEN"]
    conda:
        config["TRIMMOMATIC"]["ENV"]
    log:
        "logs/run_qc/{sample}.trimmomatic.stats.log"
    threads:
        config["TRIMMOMATIC"]["THREADS"]
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.f_read} {input.r_read} "
        "{output.p1} {output.u1} {output.p2} {output.u2} "
        "{params.illumina_clip} {params.leading} {params.trailing} "
        "{params.sliding_window} {params.minlen} 2> {log}"


rule qc_stats:
    input:
        expand("logs/run_qc/{sample}.trimmomatic.stats.log", sample = SAMPLES)
    output:
        OUTDIR + "RunQC/trimmomatic.stats"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/trimmomatic_stats.py -i {input} -o {output}"


if config["BWA"]["HOST_INDEX"] == "":
    rule build_host_index:
        input:
            HOST_FILE
        output:
            HOST_FILE + ".amb",
            HOST_FILE + ".ann",
            HOST_FILE + ".bwt",
            HOST_FILE + ".pac",
            HOST_FILE + ".sa"
        conda:
            "envs/bwa.yaml"
        shell:
            "bwa index {input}"


rule align_reads_to_host:
    input:
        HOST_FILE + ".amb",
        HOST_FILE + ".ann",
        HOST_FILE + ".bwt",
        HOST_FILE + ".pac",
        HOST_FILE + ".sa",
        host = HOST_FILE,
        fp_reads = OUTDIR + "RunQC/Paired/{sample}.1P.fastq.gz",
        rp_reads = OUTDIR + "RunQC/Paired/{sample}.2P.fastq.gz"
    output:
        temp(OUTDIR + "{sample}.host.sam")
    conda:
        config["BWA"]["ENV"]
    threads:
        config["BWA"]["THREADS"]
    shell:
        "bwa mem -t {threads} {input.host} "
        "{input.fp_reads} {input.rp_reads} > {output}"


rule host_sam_to_bam:
    input:
        OUTDIR + "{sample}.host.sam"
    output:
        OUTDIR + "AlignReadsToHost/{sample}.host.sorted.bam"
    conda:
        config["SAMTOOLS"]["ENV"]
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output}"


rule remove_host_dna:
    input:
        OUTDIR + "AlignReadsToHost/{sample}.host.sorted.bam"
    output:
        idx = OUTDIR + "RemoveHostDNA/{sample}.samtools.idxstats",
        bam = OUTDIR + "RemoveHostDNA/NonHostBAM/{sample}.host.sorted.removed.bam"
    conda:
        config["SAMTOOLS"]["ENV"]
    shell:
        "samtools index {input} && "
        "samtools idxstats {input} > {output.idx}; "
        "samtools view -h -f 4 -b {input} -o {output.bam}"


rule host_removal_stats:
    input:
        expand(OUTDIR + "RemoveHostDNA/{sample}.samtools.idxstats", sample = SAMPLES)
    output:
        OUTDIR + "RemoveHostDNA/HostRemovalStats/host.removal.stats"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/samtools_idxstats.py -i {input} -o {output}"


rule non_host_reads:
    input:
        OUTDIR + "RemoveHostDNA/NonHostBAM/{sample}.host.sorted.removed.bam"
    output:
        fq = OUTDIR + "NonHostReads/{sample}.non.host.R1.fastq.gz",
        fq2 = OUTDIR + "NonHostReads/{sample}.non.host.R2.fastq.gz"
    conda:
        config["BEDTOOLS"]["ENV"]
    shell:
        "bedtools bamtofastq -i {input} -fq {output.fq} -fq2 {output.fq2}"


if config["KRAKEN"]["INCLUDE"] == "true":
    include: "kraken.snakefile"


if config["BWA"]["AMR_INDEX"] == "":
    rule build_amr_index:
        input:
            AMR_FILE
        output:
            AMR_FILE + ".amb",
            AMR_FILE + ".ann",
            AMR_FILE + ".bwt",
            AMR_FILE + ".pac",
            AMR_FILE + ".sa"
        conda:
            "envs/bwa.yaml"
        shell:
            "bwa index {input}"


rule align_to_amr:
    input:
        AMR_FILE + ".amb",
        AMR_FILE + ".ann",
        AMR_FILE + ".bwt",
        AMR_FILE + ".pac",
        AMR_FILE + ".sa",
        amr = AMR_FILE,
        fnh_reads = OUTDIR + "NonHostReads/{sample}.non.host.R1.fastq.gz",
        rnh_reads = OUTDIR + "NonHostReads/{sample}.non.host.R2.fastq.gz"
    output:
        temp(OUTDIR + "{sample}.amr.alignment.sam")
    conda:
        config["BWA"]["ENV"]
    params:
        rg = r"@RG\tID:${sample}\tSM:${sample}"
    threads:
        config["BWA"]["THREADS"]
    shell:
        "bwa mem -t {threads} -R '{params.rg}' {input.amr} "
        "{input.fnh_reads} {input.rnh_reads} > {output}"


rule amr_sam_to_bam:
    input:
        OUTDIR + "{sample}.amr.alignment.sam"
    output:
        sorted_bam = temp(OUTDIR + "AlignToAMR/{sample}.amr.alignment.sorted.bam"),
        sorted_fix_bam = temp(OUTDIR + "AlignToAMR/{sample}.amr.alignment.sorted.fix.bam"),
        sorted2_fix_bam = temp(OUTDIR + "AlignToAMR/{sample}.amr.alignment.sorted.fix.sorted.bam"),
        bam = OUTDIR + "AlignToAMR/{sample}.amr.alignment.dedup.bam",
        sam = OUTDIR + "AlignToAMR/{sample}.amr.alignment.dedup.sam"
    conda:
        config["SAMTOOLS"]["ENV"]
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -n -@ {threads} -o {output.sorted_bam}; "
        "samtools fixmate {output.sorted_bam} {output.sorted_fix_bam}; "
        "samtools sort -@ {threads} {output.sorted_fix_bam} -o {output.sorted2_fix_bam}; "
        "samtools rmdup -S {output.sorted2_fix_bam} {output.bam}; "
        "samtools view -h -o {output.sam} {output.bam}"


# `-type_fp` option not working, according to the github that isn't even an option?
rule run_resistome:
    input:
        "build_resistome.done",
        sam = OUTDIR + "{sample}.amr.alignment.sam",
        amr = AMR_FILE,
        annotation = ANNOTATION_FILE
    output:
        gene_fp = OUTDIR + "RunResistome/{sample}.gene.tsv",
        group_fp = OUTDIR + "RunResistome/{sample}.group.tsv",
        mech_fp = OUTDIR + "RunResistome/{sample}.mechanism.tsv",
        class_fp = OUTDIR + "RunResistome/{sample}.class.tsv"
        #type_fp = OUTDIR + "RunResistome/{sample}.type.tsv"
    conda:
        config["WORKFLOW"]["ENV"]
    params:
        threshold = config["RESISTOME"]["THRESHOLD"]
    shell:
        "bin/resistome "
        "-ref_fp {input.amr} "
        "-annot_fp {input.annotation} "
        "-sam_fp {input.sam} "
        "-gene_fp {output.gene_fp} "
        "-group_fp {output.group_fp} "
        "-mech_fp {output.mech_fp} "
        "-class_fp {output.class_fp} "
        #"-type_fp {output.type_fp} "
        "-t {params.threshold}"


rule resistome_results:
    input:
        expand(OUTDIR + "RunResistome/{sample}.gene.tsv", sample = SAMPLES)
    output:
        OUTDIR + "ResistomeResults/AMR_analytic_matrix.csv"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/amr_long_to_wide.py -i {input} -o {output}"


rule sam_dedup_run_resistome:
    input:
        "build_resistome.done",
        sam = OUTDIR + "AlignToAMR/{sample}.amr.alignment.dedup.sam",
        amr = AMR_FILE,
        annotation = ANNOTATION_FILE
    output:
        gene_fp = OUTDIR + "SamDedupRunResistome/{sample}.gene.tsv",
        group_fp = OUTDIR + "SamDedupRunResistome/{sample}.group.tsv",
        mech_fp = OUTDIR + "SamDedupRunResistome/{sample}.mechanism.tsv",
        class_fp = OUTDIR + "SamDedupRunResistome/{sample}.class.tsv"
    conda:
        config["WORKFLOW"]["ENV"]
    params:
        threshold = config["RESISTOME"]["THRESHOLD"]
    shell:
        "bin/resistome "
        "-ref_fp {input.amr} "
        "-annot_fp {input.annotation} "
        "-sam_fp {input.sam} "
        "-gene_fp {output.gene_fp} "
        "-group_fp {output.group_fp} "
        "-mech_fp {output.mech_fp} "
        "-class_fp {output.class_fp} "
        #"-type_fp {output.type_fp} "
        "-t {params.threshold}"


rule sam_dedup_resistome_results:
    input:
        expand(OUTDIR + "SamDedupRunResistome/{sample}.gene.tsv", sample = SAMPLES)
    output:
        OUTDIR + "SamDedup_ResistomeResults/SamDedup_AMR_analytic_matrix.csv"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/amr_long_to_wide.py -i {input} -o {output}"


rule run_rarefaction:
    input:
        "build_rarefaction.done",
        sam = OUTDIR + "{sample}.amr.alignment.sam",
        amr = AMR_FILE,
        annotation = ANNOTATION_FILE
    output:
        gene_fp = OUTDIR + "RunRarefaction/{sample}.gene.tsv",
        group_fp = OUTDIR + "RunRarefaction/{sample}.group.tsv",
        mech_fp = OUTDIR + "RunRarefaction/{sample}.mechanism.tsv",
        class_fp = OUTDIR + "RunRarefaction/{sample}.class.tsv"
    params:
        rare_min = config["RAREFACTION"]["MIN"],
        rare_max = config["RAREFACTION"]["MAX"],
        rare_skip = config["RAREFACTION"]["SKIP"],
        rare_samples = config["RAREFACTION"]["SAMPLES"],
        threshold = config["RAREFACTION"]["THRESHOLD"]
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/rarefaction "
        "-ref_fp {input.amr} "
        "-sam_fp {input.sam} "
        "-annot_fp {input.annotation} "
        "-gene_fp {output.gene_fp} "
        "-group_fp {output.group_fp} "
        "-mech_fp {output.mech_fp} "
        "-class_fp {output.class_fp} "
        #"-type_fp ${sample_id}.type.tsv "
        "-min {params.rare_min} "
        "-max {params.rare_max} "
        "-skip {params.rare_skip} "
        "-samples {params.rare_samples} "
        "-t {params.threshold}"