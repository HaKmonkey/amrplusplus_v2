Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

## Notes

Qiime2 has been added for feature parity with the NextFlow version of the
pipeline, however I have been unable to successfully test it, so I would
recommend against using it in this version for now. Testing will continue.

## Requirements

Dependancies are managed through [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Install Snakemake and Clone the Repository

Create the environment for amrplusplus

```bash
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n amrplusplus snakemake git
mamba activate amrplusplus
```

Clone the repository.

```bash
mamba activate amrplusplus
git clone https://github.com/HaKmonkey/amrplusplus_v2.git
```

## Usage

The pipeline comes with some test data in `data/raw/`. You can run the pipeline
as is and make sure it is working correctly.

You will want to replace the following with your own files:

### Raw Reads

- `data/raw/S{1,2,3}_test_R{1,2}.fastq.gz`
- Accepts sequence files in standard fastq and gz format
- The pipeline automatically grabs the sample name, so make sure your files end
with `_R1.fastq.gz` or `_R2.fastq.gz`

### Host Genome

- `data/host/chr21.fasta.gz`
- Accepts a fasta formatted host genome
- Adjust `config["BWA"]["HOST"]` in `config.json` to match

### Resistance Database

- `data/amr/megares_modified_database_v2.00.fasta`
- Accepts a fasta formatted resistance database
- Adjust `config["BWA"]["AMR"]` in `config.json` to match

### Annotation Database

- `data/amr/megares_modified_annotations_v2.00.csv`
- Accepts a csv formatted annotation database
- Adjust `config["RESISTOME"]["ANNOTATION"]` in `config.json` to match

### Adapters

- `data/adapters/adapters.fa`
- Accepts a fasta formatted adapter file
- Adjust `config["TRIMMOMATIC"]["ADAPTERS"]` in `config.json` to match

### Run Workflow

```bash
cd amrplusplus_v2
mamba activate amrplusplus

snakemake --use-conda --cores <number of threads available>
```

## Microbial Ecology Group (MEG)
(https://megares.meglab.org/)

Our international multidisciplinary group of scientists and educators is addressing the issues of antimicrobial resistance (AMR) and microbial ecology in agriculture through research, outreach, and education. By characterizing risks related to AMR and microbial ecology, our center will identify agricultural production practices that are harmful and can be avoided, while also identifying and promoting production practices and interventions that are beneficial or do no harm to the ecosystem or public health. This will allow society to realize “sustainable intensification” of agriculture.

## MEGARes and the AMR++ bioinformatic pipeline
(http://megares.meglab.org/amrplusplus/latest/html/v2/)

The MEGARes database contains sequence data for approximately 8,000 hand-curated antimicrobial resistance genes accompanied by an annotation structure that is optimized for use with high throughput sequencing and metagenomic analysis. The acyclical annotation graph of MEGARes allows for accurate, count-based, hierarchical statistical analysis of resistance at the population level, much like microbiome analysis, and is also designed to be used as a training database for the creation of statistical classifiers.

The goal of many metagenomics studies is to characterize the content and relative abundance of sequences of interest from the DNA of a given sample or set of samples. You may want to know what is contained within your sample or how abundant a given sequence is relative to another.

Often, metagenomics is performed when the answer to these questions must be obtained for a large number of targets where techniques like multiplex PCR and other targeted methods would be too cumbersome to perform. AmrPlusPlus can process the raw data from the sequencer, identify the fragments of DNA, and count them. It also provides a count of the polymorphisms that occur in each DNA fragment with respect to the reference database.

Additionally, you may want to know if the depth of your sequencing (how many reads you obtain that are on target) is high enough to identify rare organisms (organisms with low abundance relative to others) in your population. This is referred to as rarefaction and is calculated by randomly subsampling your sequence data at intervals between 0% and 100% in order to determine how many targets are found at each depth.

With AMR++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Output](docs/output.md)
- [FAQs](docs/FAQs.md)
- [Contact](docs/contact.md)
