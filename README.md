# Micro-C Pipeline

Snakemake pipeline for preprocessing and analysis of Micro-C data. Pipeline follows the [Dovetail Micro-C pipeline](https://micro-c.readthedocs.io/en/latest/index.html) to align reads, build deduplicated pairs files, report library QC, and build `.hic` and `.mcool` contact maps. It also calls A/B compartments, TADs, loops, structural variants, copy number variants, and SV-induced neo-loops/neo-TADs. Runs initial preprocesing on shallowly sequenced libraries to profile complexity before full sequencing. 

## Pipeline steps

### Run for all samples:

| Step | Rule | Tool | Description |
|---|---|---|---|
| 1 | `align` | BWA-MEM | Align each library to the reference |
| 2 | `generate_pairs` | pairtools, samtools | Parse ligation junctions, sort and deduplicate pairs, write pairs file and BAM |
| 3 | `library_qc`, `combine_qc` | Dovetail `get_qc.py` | Library QC metrics, combined across libraries |
| 4 | `merge_pairs` | pairtools | Merge library pairs files per sample |
| 5 | `contact_matrix_hic` | Juicer Tools `pre` | `.hic` contact map |
| 6 | `contact_matrix_cooler` | pairix, cooler | 1 kb `.cool` and balanced multi-resolution `.mcool` |
| 7 | `combine_bam` | samtools | Merge library BAMs per sample |

### Additional steps for shallow sequencing:

| Rule | Tool | Description |
|---|---|---|
| `library_complexity` | preseq | Library complexity curves |

### Additional steps for deep sequencing:

| Rule | Tool | Description |
|---|---|---|
| `convert_hic_h5` | HiCExplorer | 10 kb matrix in HiCExplorer `.h5` format |
| `call_compartment_fanc`, `plot_compartment` | FAN-C | A/B compartments at 100 kb, with per-chromosome plots |
| `call_tads_arrowhead` | Juicer Tools Arrowhead | TADs at 25, 10 and 5 kb |
| `call_tads_hic` | HiCExplorer `hicFindTADs` | TADs at 10 kb |
| `call_loops_mustache` | Mustache | Loops at 10 and 5 kb |
| `hic_breakfinder` | hic_breakfinder | Structural variants |
| `eaglec`, `reformat_eaglec` | EagleC | Structural variants (5 kb) |
| `call_cnv_neoloop` | NeoLoopFinder | Copy number at 25 kb |
| `assemble_sv`, `call_neoloop`, `call_neotad` | NeoLoopFinder | Assemble rearranged regions and call neo-loops and neo-TADs |

## Requirements

- Runs on NIH [Biowulf](https://hpc.nih.gov/) (SLURM and environment modules).
- Must have access to dependencies listed in [config/config.yaml](config/config.yaml).
- Install [cooler](https://cooler.readthedocs.io/en/latest/quickstart.html#installation).

## Installation

```bash
git clone https://github.com/CCRGeneticsBranch/MicroC.git
```

## Quick start

The pipeline is designed to be run in a separate directory (/path/to/my_run).

```bash
# 1. Create a work folder with template config files; choose the depth
bash MicroC/run_microc --runmode=init --workdir=/path/to/my_run --depth=shallow [--samplesheet=/path/to/samples.tsv]

# 2. Edit /path/to/my_run/config/config.yaml and config/samples.tsv (see Inputs)

# 3. Dry run (writes /path/to/my_run/dryrun.log)
bash MicroC/run_microc --runmode=dryrun --workdir=/path/to/my_run

# 4. Submit to SLURM
bash MicroC/run_microc --runmode=run --workdir=/path/to/my_run

# Unlock the work folder after an interrupted run
bash MicroC/run_microc --runmode=unlock --workdir=/path/to/my_run
```

`-d/--depth` (`shallow`, the default, or `deep`) is only used at `init`. It sets `depth` in `config.yaml` and copies the matching cluster resources: [resources/cluster_shallow.yaml](resources/cluster_shallow.yaml) or [resources/cluster.yaml](resources/cluster.yaml), to `WORKDIR/config/cluster.yaml`. `-h` shows all options; `-v` prints the pipeline and launcher versions.

## Inputs

### Sample sheet (`config/samples.tsv`)

One row per library. Libraries with the same `Sample` are merged.

| Column | Description | Example |
|---|---|---|
| `Library` | Library name | `library1` |
| `Sample` | Sample name | `sample1` |
| `R1` | Absolute path to read 1 FASTQ | `/path/to/library1_R1.fastq.gz` |
| `R2` | Absolute path to read 2 FASTQ | `/path/to/library1_R2.fastq.gz` |

### Configuration (`config/config.yaml`)

| Key | Default | Description (source) | Edit? |
|---|---|---|---|
| `depth` | `shallow` | `shallow` or `deep` | Set by `--depth` at `init` |
| `genome` | `hg38` | Only `hg38` is configured |  |
| `hg38: ref`, `bwa_ref` | `/data/khanlab3/gb_omics/resources/genomes/hg38/hg38_main.fa` | Main-chromosome hg38 FASTA and BWA index ([UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)) | |
| `hg38: chrom_sizes` | `/data/khanlab3/gb_omics/resources/genomes/hg38/hg38_main.genome` | Chromosome sizes derived from the reference (`samtools faidx`) |  |
| `hg38: exp_file_inter`, `exp_file_intra` | `/data/khanlab3/gb_omics/resources/hic/hic_breakfinder/...` | Expected inter- and intrachromosomal contacts (see [hic_breakfinder](https://github.com/dixonlab/hic_breakfinder)) | |
| `ploidy` | `2` | Ploidy used by NeoLoopFinder | As needed |
| `microcdir` | `/data/khanlab3/gb_omics/software/Micro-C/` | Dovetail Micro-C repository ([source](https://github.com/dovetail-genomics/Micro-C)) |  |
| `pairixdir` | `/data/khanlab3/gb_omics/software/pairix/bin/` | pairix repository ([source](https://github.com/4dn-dcic/pairix)) |  |
| `workdir`, `manifest`, `tools`, `cluster` | `WORKDIR/...`, `PIPELINE_HOME/...` | Work folder and launcher/repository paths | Filled by `init` |

## Outputs

All paths are relative to the work folder.

| Path | Depth | Description |
|---|---|---|
| `results/bam/<library>.final.bam`, `.pairs.gz`, `.dedup_stats.txt` | all | Per-library alignments, pairs, and deduplication statistics |
| `results/qc/<library>.qc.txt`, `results/qc/combined_qc.txt` | all | Library QC |
| `results/merged_pairs/<sample>.pairs.gz` | all | Merged pairs per sample |
| `results/bam/<sample>.combined.bam` | all | Merged BAM per sample |
| `results/hic/<sample>.hic` | all | Juicer contact map |
| `results/cool/<sample>.1kb.cool`, `.1kb.mcool` | all | Cooler contact maps |
| `results/qc/<library>.preseq`, `results/qc/multiqc_report.html` | shallow | Library complexity |
| `results/h5/<sample>.10kb.h5` | deep | HiCExplorer matrix |
| `results/topology/compartment/fanc/` | deep | A/B compartments (`_100kb_AB.bed`, eigenvector, plots) |
| `results/topology/tads/arrowhead/`, `results/topology/tads/hicexplorer/` | deep | TADs |
| `results/topology/loops/mustache/` | deep | Loops |
| `results/topology/sv/hic_breakfinder/`, `results/topology/sv/eaglec/` | deep | Structural variants |
| `results/topology/cnv/` | deep | Copy number profile and segments (25 kb) |
| `results/topology/sv/neoloopfinder/` | deep | SV assemblies, neo-loops, and neo-TADs |
| `results/pipeline_version.txt` | all | Pipeline version used for the run |
| `runinfo.yaml`, `snakemake.log`, `logs/` | all | Run information, Snakemake log and per-rule SLURM logs |

## Tools and references

Modules are listed in [resources/tools.yaml](resources/tools.yaml).

| Tool | Module / source | Used in | Reference |
|---|---|---|---|
| Snakemake | `snakemake/7.19.1` | workflow | Mölder et al. 2021, *F1000Research* [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| Micro-C | — | method | Hsieh et al. 2020, *Mol Cell* [doi:10.1016/j.molcel.2020.03.002](https://doi.org/10.1016/j.molcel.2020.03.002); [Dovetail Micro-C docs](https://micro-c.readthedocs.io/) |
| BWA-MEM | `bwa` | `align` | Li 2013, [arXiv:1303.3997](https://arxiv.org/abs/1303.3997) |
| pairtools | `pairtools` | `generate_pairs`, `merge_pairs` | Open2C et al. 2024, *PLoS Comput Biol* [doi:10.1371/journal.pcbi.1012164](https://doi.org/10.1371/journal.pcbi.1012164) |
| SAMtools | `samtools` | `generate_pairs`, `combine_bam` | Danecek et al. 2021, *GigaScience* [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| preseq | `preseq` | `library_complexity` | Daley & Smith 2013, *Nat Methods* [doi:10.1038/nmeth.2375](https://doi.org/10.1038/nmeth.2375) |
| MultiQC | `multiqc` | `multiqc` | Ewels et al. 2016, *Bioinformatics* [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |
| Juicer Tools | `juicer` | `contact_matrix_hic`, `call_tads_arrowhead` | Durand et al. 2016, *Cell Syst* [doi:10.1016/j.cels.2016.07.002](https://doi.org/10.1016/j.cels.2016.07.002) |
| pairix | `pairixdir` | `contact_matrix_cooler` | [4dn-dcic/pairix](https://github.com/4dn-dcic/pairix) |
| cooler | `PATH` | `contact_matrix_cooler` | Abdennur & Mirny 2020, *Bioinformatics* [doi:10.1093/bioinformatics/btz540](https://doi.org/10.1093/bioinformatics/btz540) |
| HiCExplorer | `hicexplorer` | `convert_hic_h5`, `call_tads_hic` | Wolff et al. 2020, *Nucleic Acids Res* [doi:10.1093/nar/gkaa220](https://doi.org/10.1093/nar/gkaa220) |
| FAN-C | `fanc` | `call_compartment_fanc`, `plot_compartment` | Kruse et al. 2020, *Genome Biol* [doi:10.1186/s13059-020-02215-9](https://doi.org/10.1186/s13059-020-02215-9) |
| Mustache | `mustache` | `call_loops_mustache` | Roayaei Ardakany et al. 2020, *Genome Biol* [doi:10.1186/s13059-020-02167-0](https://doi.org/10.1186/s13059-020-02167-0) |
| hic_breakfinder | `hic_breakfinder` | `hic_breakfinder` | Dixon et al. 2018, *Nat Genet* [doi:10.1038/s41588-018-0195-8](https://doi.org/10.1038/s41588-018-0195-8) |
| EagleC | `eaglec` | `eaglec` | Wang et al. 2022, *Sci Adv* [doi:10.1126/sciadv.abn9215](https://doi.org/10.1126/sciadv.abn9215) |
| NeoLoopFinder | `neoloopfinder` | `call_cnv_neoloop`, `assemble_sv`, `call_neoloop`, `call_neotad` | Wang et al. 2021, *Nat Methods* [doi:10.1038/s41592-021-01164-w](https://doi.org/10.1038/s41592-021-01164-w) |

## Versioning

The pipeline version is in [VERSION](VERSION) and changes are listed in [CHANGELOG.md](CHANGELOG.md). Each run records the pipeline version, git commit and launcher version in `runinfo.yaml` and `results/pipeline_version.txt`. The launch script is a copy of [biowulf-snakemake-launcher](https://github.com/CCRGeneticsBranch/biowulf-snakemake-launcher) (`bin/launcher.sh`).

## Acknowledgments

Thank you to Hsien-Chao Chou for sharing library merging and topology scripts.

## Contact

Erica Pehrsson, Genetics Branch, Center for Cancer Research, National Cancer Institute. Please report problems through GitHub issues.
