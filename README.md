# Shovill Nextflow Pipeline

This Nextflow pipeline automates genome assembly using Shovill for paired-end FASTQ reads. It supports both samplesheet input and automatic file pair discovery.

## Features
- Runs Shovill assembly on multiple samples in parallel
- Flexible input: CSV samplesheet or automatic FASTQ pairing
- Configurable assembler (SPAdes, SKESA, etc.)
- Multi-threaded processing
- Publishes assembled contigs as `${sample_id}_contigs.fa`

## Requirements
- Nextflow ≥ 22.10.0
- Shovill (with assembler dependencies)
- Paired-end FASTQ files (.fastq/.fq/.gz)

## Usage

```bash
nextflow run BCCDC-PHL/nf_shovill \
  --outdir ./results \
  --threads 8 \
  [--samplesheet_input samples.csv] \
  [--fastq_input "path/to/fastq/*_{R1,R2}*.fastq.gz"]
```
## Parameters

| Parameter           | Default                    | Description                                |
| ------------------- | -------------------------- | ------------------------------------------ |
| --outdir            | ./results                  | Output directory                           |
| --threads           | 8                          | CPU threads per sample                     |
| --assembler         | spades                     | Shovill assembler (spades, skesa, megahit) |
| --samplesheet_input | NO_FILE                    | CSV samplesheet path                       |
| --fastq_input       | Required if no samplesheet | FASTQ search pattern                       |

## Output

```bash
results/
├── sample1_contigs.fa
├── sample2_contigs.fa
└── ...
```