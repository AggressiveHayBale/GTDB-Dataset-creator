# GTDB Dataset creator
## Overview

The GTDB-Dataset-creator is a command-line tool designed for systematic construction of comprehensive genomic datasets from the Genome Taxonomy Database (GTDB). This tool enables researchers to:

- Precisely select taxonomic subsets from GTDB metadata
- Apply customizable quality filters based on genomic metrics
- Implement advanced sampling strategies for representative dataset construction
- Download FASTA files via the NCBI datasets CLI

The tool's primary purpose is to create well-encompassing datasets specifically designed for benchmarking genomic analysis pipelines or performing bioinformatic analyses on robust, user-specified genomic collections. It ensures dataset quality, diversity, and representativeness according to user-defined parameters.

## Requirements

- R (≥ 4.0.0)
- Required R packages: data.table, R.utils, cluster
- NCBI datasets CLI tool (automatically downloaded if not present)
- Sufficient disk space for genome downloads

## Installation

The tool is self-contained and requires no explicit installation:
```
chmod +x dataset_preparator.R 
```
Dependencies will be automatically installed on first execution.

## Usage
```
./dataset_preparator.R --tax_name "Pseudomonas aeruginosa" --n_samples 20
```
## Parameters

Core Parameters

- --tax_name, -n [required]
Target taxon name (e.g., "Pseudomonas aeruginosa")

- --taxrank, -t [default: species]
Taxonomic rank for selection: species, genus, or family

- --domain, -m [default: bacteria]
Domain to query: bacteria or archaea

- --representative, -r [default: f]
Restrict to representative genomes only (t/f)

Sampling Parameters

- --n_samples, -ns [required]
Total number of genomes to sample

- --sampling_method [default: custom]
Sampling strategy: custom, kmedoids, maxmin, stratified, or hybrid

- --sampling, -ss [default: "A:0.5,O:0.25,R:0.25"]
Proportion schema for custom sampling:
- A: Average/representative genomes
- O: Outlier/divergent genomes
- R: Random genomes

- --quality, -q [default: good]
Quality filter level: good, bad, or mixed

- --min_size, --max_size
Genome size filters (in base pairs)

Operational Parameters

- --output_dir, -o [default: output]
Directory for downloaded genomes

- --database, -db
Path to custom GTDB metadata file

- --dryrun, -d
Perform a dry run without downloading genomes

- --help, -h
Display help information

## Sampling Methods

The tool implements five distinct sampling strategies designed to create robust benchmarking datasets:

1. Custom
Combines representative, outlier, and random selections according to specified proportions to balance diversity and representativeness.

2. k-Medoids
Uses partitioning around medoids clustering to select maximally representative genomes for consistent benchmarking performance.

3. MaxMin
Selects genomes to maximize minimum pairwise distances, ensuring maximum diversity in the dataset.

4. Stratified
Samples proportionally across genome categories (e.g., isolate, metagenome-assembled), preserving dataset composition.

5. Hybrid
Combines stratification with MaxMin sampling within each stratum for optimal diversity within categories.

All distance-based methods utilize Gower's dissimilarity metric across 28 genomic features including:
- Taxonomic classification
- Genome quality metrics (completeness, contamination)
- Assembly statistics (N50, contig count)
- Functional gene counts

## Output

The tool generates:

1. A directory containing downloaded genome FASTA files (*.fna)
2. sample_metadata.csv: Comprehensive metadata for all selected genomes
3. download_accession_list.txt: List of NCBI accessions for downloaded genomes

## Examples

### Basic usage for high-quality Pseudomonas aeruginosa genomes for benchmarking
```
./dataset_preparator.R --tax_name "Pseudomonas aeruginosa" --n_samples 50 --quality good --output_dir paeruginosa_benchmark_set 
```
### Advanced sampling with custom proportions for comparative analysis
```
./dataset_preparator.R --tax_name "Enterobacteriaceae" --taxrank family --n_samples 100 --sampling "A:0.4,O:0.3,R:0.3" --sampling_method custom --domain bacteria
```
### Stratified sampling by genome category for comprehensive benchmarking
```
./dataset_preparator.R \ --tax_name "Mycobacterium tuberculosis"  --n_samples 30 --sampling_method stratified --quality mixed 
```
## Notes

- The tool automatically downloads the latest GTDB metadata (release 226) if not present
- Quality filtering thresholds are dynamically calculated based on quartiles of the dataset
- When insufficient genomes meet quality criteria, thresholds are progressively relaxed
- All distance calculations use Gower's metric to handle mixed data types appropriately
- The resulting datasets are specifically designed to be robust and comprehensive for benchmarking studies

For further information about GTDB, please refer to the official GTDB website.

