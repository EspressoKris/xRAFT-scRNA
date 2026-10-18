# xRAFT-scRNA: 10x-express RAw FASTQ To single-cell RNA-seq objects

# Project information:

![Static Badge](https://img.shields.io/badge/Nextflow-25.04-23aa62?style=flat&logo=nextflow) ![Static Badge](https://img.shields.io/badge/Cell_Ranger-9.0-0b5fa5?style=flat) ![Static Badge](https://img.shields.io/badge/R-4.4-276DC3?style=flat&logo=r) ![Static Badge](https://img.shields.io/badge/Slurm-executor-5c5c5c?style=flat) ![Static Badge](https://img.shields.io/badge/Apptainer-containers-e8834a?style=flat)

---

## General Description

**Author:** Kristian Gurashi (kristian.gurashi@gmail.com)

**Pipeline:** xRAFT-scRNA: 10x-express RAw FASTQ To single-cell RNA-seq objects

**Scope:** 10x Genomics single-cell RNA-seq, from raw sequencing archives to analysis-ready Seurat (R) and AnnData (Scanpy, Python) objects

Nextflow DSL2 pipeline that takes scRNA-seq data from raw sequencing archives (for example as delivered by Novogene) or existing FASTQs to analysis-ready single-cell objects, in one unattended run on a Slurm cluster. It covers FASTQ extraction, Cell Ranger count, Souporcell genetic demultiplexing for pooled libraries, CellBender ambient-RNA removal, and object creation for either ecosystem: a Seurat object (`.rds`) for R or an AnnData object (`.h5ad`) for Scanpy, or both. Steps whose outputs already exist are skipped, so a failed or extended run resumes where it stopped.

**Why?** Every 10x project begins with the same manual stage: extract and verify archives, run Cell Ranger, remove ambient RNA, demultiplex pooled donors, build an object. Done by hand it is slow, inconsistent and not reproducible. nf-core/scrnaseq skips archive extraction and genetic demultiplexing; xRAFT-scRNA runs exactly these steps from one command and one sample sheet, and returns the same object every time, in R or Python.

## Graphical Workflow

<p align="center">
  <img src="docs/images/pipeline_map.svg" width="100%" alt="Graphical workflow">
</p>

Four lines show the combinations the sample sheet allows: raw archives or existing FASTQs, pooled or single-donor libraries. Filled stations run by default; hollow stations (Extract, Souporcell) run only when the sample sheet asks for them. The terminal station writes a Seurat `.rds`, an AnnData `.h5ad`, or both, according to `params.output_format`.

## Content

One subworkflow per step in `subworkflows/`:

| Step | Subworkflow | Runs when |
| --- | --- | --- |
| 1. Extract FASTQs from tar archives (e.g. Novogene), with MD5 verification | `extraction.nf` | `Raw_Data` is set |
| 2. Build the Cell Ranger libraries CSV (auto-detects FASTQ prefixes) | `cellranger_libraries.nf` | always |
| 3. Cell Ranger count | `cellranger_count.nf` | always |
| 4. Souporcell genetic demultiplexing | `souporcell.nf` | `souporcell_k > 1` |
| 5. CellBender ambient-RNA removal (GPU) | `cellbender.nf` | always |
| 6. Seurat object (RNA + RNA_cellbender assays, Souporcell metadata) | `create_seurat_object.nf` | `output_format` is `seurat` or `both` |
| 7. AnnData object (raw and CellBender layers, Souporcell metadata) | `create_anndata_object.nf` | `output_format` is `scanpy` or `both` |

## Run time

Measured with Slurm accounting over 20 runs on the lab cluster (16 CPUs, 96 GB for Cell Ranger and Souporcell; 4 CPUs, 24 GB, one GPU for CellBender). Queue time excluded.

<p align="center">
  <img src="docs/images/runtime.png" width="100%" alt="Run time per step and end to end">
</p>

**End to end:** a single-donor library in about 4.5 h, a pooled library in about 8.5 h (range 5 to 12.5 h, n = 10). Libraries run in parallel: eight pooled libraries finished overnight, in 13 h.

## How to run it

You need three things: your sequencing data on the cluster, a sample sheet, and one command.

**1. Put your data where the pipeline can see it.** Copy the tar archives from the sequencing provider (for example Novogene) into a folder on the cluster, or point to a folder that already holds FASTQ files.

**2. Fill in the sample sheet.** Open `data/sample_sheet.csv` in Excel or a text editor. One row per library, four columns:

| Column | What to write | Example |
| --- | --- | --- |
| `library_id` | A name for the library (the same name as its folder inside the archive) | `EGO_12` |
| `fastq_dir` | Where the FASTQs should go, or already are | `/data/fastq/EGO_12` |
| `Raw_Data` | The archive file(s) to extract, separated by `;`. Leave empty if your FASTQs already exist | `/data/raw/EGO_12.tar` |
| `souporcell_k` | How many donors were pooled in this library. Write `1` or leave empty if it is a single donor | `4` |

Any extra columns you add (for example `experiment` or `batch_10x`) are kept and copied into the final object.

**3. Choose your output.** In `nextflow.config`, set `params.output_format` to `seurat` (R), `scanpy` (Python) or `both`. Everything else in that file is set up for the lab cluster; you only need to touch it if you want results in a different folder.

**4. Run.** From the pipeline folder:

```bash
sbatch start_pipeline.sh
```

Then log off. The pipeline queues every step on Slurm, runs the libraries in parallel and finishes on its own; a pooled library takes around 8 hours, so a batch is usually ready the next morning.

**5. Collect your objects.** Each library gets its own folder with the finished object, ready to load:

- Seurat: `<seurat_results>/<library_id>/<library_id>_seurat.rds`
- AnnData: `<anndata_results>/<library_id>/<library_id>.h5ad`

The objects hold raw and CellBender-corrected counts, with the donor assignment for pooled libraries, so analysis can start straight away.

**If something goes wrong,** submit the same command again. The pipeline sees which steps already finished and only reruns what is missing. The per-step logs are in the `logs/` folder, and `reports/pipeline_report.html` shows what ran and how long it took.

## Dependencies

| Component | Version | Used for |
| --- | --- | --- |
| Nextflow | 25.04 or later, Slurm executor | Workflow engine |
| Cell Ranger | 9.0 or later | Alignment and counting |
| Apptainer or Singularity | any recent | Containers for Souporcell and CellBender |
| CellBender | container, GPU | Ambient-RNA removal |
| Souporcell | container | Genetic demultiplexing of pooled libraries |
| R | 4.4 or later, with Seurat and scCustomize | Seurat object output |
| Python | 3.10 or later, with Scanpy and anndata | AnnData object output |
| References | Cell Ranger transcriptome, reference FASTA, common-variants VCF | Counting and Souporcell |

Module names are set in `start_pipeline.sh` and `nextflow.config`. Result directories and references live in [nextflow.config](nextflow.config) (`params.samplesheet`, `params.cellranger_results`, `params.souporcell_results`, `params.cellbender_results`, `params.seurat_results`, `params.anndata_results`, `params.transcriptome`, `params.ref_fasta`, `params.common_variants`, `params.logs_dir`); map the Slurm labels `basic`, `build`, `long` and `gpu` to your partitions. The wrapper runs Nextflow from a job-specific scratch directory with `-resume` and copies `.nextflow.log` back when the run finishes; `nextflow run main.nf -resume` works too. Use `clean_nf_logs.sh` to remove old `.nextflow.log*` files.

---

# Networking:

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/kristiangurashi/) [![GitHub](https://img.shields.io/badge/GitHub-181717?style=flat&logo=github&logoColor=white)](https://github.com/EspressoKris)