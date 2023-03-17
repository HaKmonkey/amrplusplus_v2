configfile: "config.json"

rule make_qiime2_manifest:
    output:
        OUTDIR + "Qiime2Results/manifest.tsv"
    params:
        indir = config["WORKFLOW"]["READS_SOURCE"]
    conda:
        config["BUILD"]["ENV"]
    envmodules:
        "python/3.8"
    shell:
        "bin/qiime2_manifest.py "
        "--indir {params.indir} "
        "--manifest {output}"

rule qiime2_import:
    input:
        OUTDIR + "Qiime2Results/manifest.tsv"
    output:
        OUTDIR + "Qiime2Results/demux.qza"
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33"

rule qiime2_dada2:
    input:
        OUTDIR + "Qiime2Results/demux.qza"
    output:
        dada_table = OUTDIR + "Qiime2Results/dada-table.qza",
        rep_seqs = OUTDIR + "Qiime2Results/rep-seqs.qza",
        denoise_stats = OUTDIR + "Qiime2Results/denoise_stats"
    params:
        trim_left_f = "--p-trim-left-f " + config["QIIME2"]["P_TRIM_LEFT_F"],
        trim_left_r = "--p-trim-left-r " + config["QIIME2"]["P_TRIM_LEFT_R"],
        trunc_len_f = "--p-trunc-len-f " + config["QIIME2"]["P_TRUNC_LEN_F"],
        trunc_len_r = "--p-trunc-len-r " + config["QIIME2"]["P_TRUNC_LEN_R"]
    threads:
        config["QIIME2"]["THREADS"]
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs {input} "
        "--o-table {output.dada_table} "
        "--o-representative-sequences {output.rep_seqs} "
        "--o-denoising-stats {output.denoise_stats} "
        "{params.trim_left_f} {params.trim_left_r} "
        "{params.trunc_len_f} {params.trunc_len_r} "
        "--p-n-threads {threads} "
        "--verbose"

rule qiime2_classify:
    input:
        rep_seqs = OUTDIR + "Qiime2Results/rep-seqs.qza",
        db = config["QIIME2"]["DB"]
    output:
        OUTDIR + "Qiime2Results/taxonomy.qza"
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime feature-classifier classify-sklearn "
        "--i-classifier {input.db} "
        "--i-reads {input.rep_seqs} "
        "--o-classification {output}"

rule qiime2_filter_table:
    input:
        dada_table = OUTDIR + "Qiime2Results/dada-table.qza",
        taxonomy = OUTDIR + "Qiime2Results/taxonomy.qza"
    output:
        OUTDIR + "Qiime2Results/filtered_table.qza",
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime taxa filter-table "
        "--i-table {input.dada_table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-table {output}"

rule qiime2_filter_reads:
    input:
        rep_seqs = OUTDIR + "Qiime2Results/rep-seqs.qza",
        taxonomy = OUTDIR + "Qiime2Results/taxonomy.qza"
    output:
        OUTDIR + "Qiime2Results/filtered_rep-seqs.qza"
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime taxa filter-seqs "
        "--i-sequences {input.rep_seqs} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-sequences {output}"

rule qiime2_tree:
    input:
        OUTDIR + "Qiime2Results/filtered_rep-seqs.qza"
    output:
        OUTDIR + "Qiime2Results/rooted-tree.qza"
    params:
        aligned = OUTDIR + "Qiime2Results/aligned-rep-seqs.qza",
        masked = OUTDIR + "Qiime2Results/masked-aligned-rep-seqs.qza",
        unrooted_tree = OUTDIR + "Qiime2Results/unrooted-tree.qza"
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        "qiime alignment mafft "
        "--i-sequences {input} "
        "--o-alignment {params.aligned}; "
        "qiime alignment mask "
        "--i-alignment {params.aligned} "
        "--o-masked-alignment {params.masked}; "
        "rm {params.aligned}; "
        "qiime phylogeny fasttree "
        "--i-alignment {params.masked} "
        "--o-tree {params.unrooted_tree}; "
        "rm {params.masked}; "
        "qiime phylogeny midpoint-root "
        "--i-tree {params.unrooted_tree} "
        "--o-rooted-tree {output}; "
        "rm {params.unrooted_tree}"

rule qiime2_export: # want output
    input:
        filtered_seqs = OUTDIR + "Qiime2Results/filtered_rep-seqs.qza",
        rooted_tree = OUTDIR + "Qiime2Results/rooted-tree.qza",
        filtered_table = OUTDIR + "Qiime2Results/filtered_table.qza",
        taxonomy = OUTDIR + "Qiime2Results/taxonomy.qza"
    output:
        table_taxa = OUTDIR + "Qiime2Results/table-with-taxonomy.biom",
        tree_nwk = OUTDIR + "Qiime2Results/tree.nwk",
        dna_seqs = OUTDIR + "Qiime2Results/dna-sequences.fasta"
    conda:
        config["QIIME2"]["ENV"]
    envmodules:
        "qiime2/2022.8"
    shell:
        """
        qiime tools export --input-path {input.filtered_seqs} --output-path {output.dna_seqs}
        qiime tools export --input-path {input.taxonomy} --output-path . 
        qiime tools export --input-path {input.rooted_tree} --output-path {output.tree_nwk}
        qiime tools export --input-path {input.filtered_table} --output-path .

        # Change out column headers in taxonomy file
        sed -i '' 's/Feature ID/#OTUID/g' taxonomy.tsv
        sed -i '' 's/Taxon/taxonomy/g' taxonomy.tsv
        sed -i '' 's/Confidence/confidence/g' taxonomy.tsv

        # Add taxonomy file to biom
        biom add-metadata -i feature-table.biom -o {output.table_taxa} --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

        rm feature-table.biom
        """

# rule qiime2_export: # want output
#     input:
#         filtered_seqs = OUTDIR + "Qiime2Results/filtered_rep-seqs.qza",
#         rooted_tree = OUTDIR + "Qiime2Results/rooted-tree.qza",
#         filtered_table = OUTDIR + "Qiime2Results/filtered_table.qza",
#         taxonomy = OUTDIR + "Qiime2Results/taxonomy.qza"
#     output:
#         table_taxa = OUTDIR + "Qiime2Results/table-with-taxonomy.biom",
#         tree_nwk = OUTDIR + "Qiime2Results/tree.nwk",
#         dna_seqs = OUTDIR + "Qiime2Results/dna-sequences.fasta"
#     conda:
#         config["QIIME2"]["ENV"]
#     envmodules:
#         "qiime2/2022.8"
#     shell:
#         """
#         qiime tools export --input-path {input.filtered_seqs} --output-path .
#         qiime tools export --input-path {input.taxonomy} --output-path . 
#         qiime tools export --input-path {input.rooted_tree} --output-path .
#         qiime tools export --input-path {input.filtered_table} --output-path .

#         # Change out column headers in taxonomy file
#         sed -i 's/Feature ID/#OTUID/g' taxonomy.tsv
#         sed -i 's/Taxon/taxonomy/g' taxonomy.tsv
#         sed -i 's/Confidence/confidence/g' taxonomy.tsv

#         # Add taxonomy file to biom
#         biom add-metadata -i feature-table.biom -o table-with-taxonomy.biom --observation-metadata-fp taxonomy.tsv --sc-separated taxonomy

#         rm feature-table.biom
#         """
