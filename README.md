# Killifish Atlas Reproducibility Repository

## Description

This repository hosts R and Python code supporting the publication of the **Killifish Aging Atlas**.

The repository provides code and configs for reproducing analyses and figures described in the manuscript.


## Citation

**Manuscript:** Costa, E.K., Chen, J., et al. Multi-tissue transcriptomic aging atlas reveals predictive aging biomarkers in the killifish. Nat Aging 6, 636–664 (2026). https://doi.org/10.1038/s43587-026-01074-6

**Code:** https://doi.org/10.5281/zenodo.16416327

**Data:** 
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1274512
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308970

## Repository Structure
- `clock_code_for_reuse/`: more user-friendly jupyter notebooks, python code, and bash scripts for aging clocks (empty in this version)
- `common_robjects/`: naming of R objects generated that are reused across R code
- `configs/`: YAML/TOML config files for reproducible runs (empty in this version)
- `data/`: raw, interim, processed data (raw empty; interim gitignored; proc contains external dataset zip files & empty dir for files provided to reviewers)
- `environment/`: lockfiles (R requirements), YAML files (Python requirements)
- `main_manuscript_workflow_R/`: R pipelines (QC, DESeq2, Correlation analysis, figures)
- `notebooks_manuscript_python/`: exploratory notebooks (ipynb)
- `other_manuscript_scripts/`: runnable scripts for clocks in python rather than jupyter notebooks, and R helper functions
- `results/`: output tables and figures (mostly gitignored)
- `tests/`: smoke/unit tests (empty in this version)

## System requirements
### Hardware requirements
- **Memory/disk:** ~16 GB RAM recommended

### Software requirements
- **OS:** Linux, macOS, or Windows with R ≥ 4.2 and Python ≥ 3.9  
- **Dependencies:**
  - R (managed with `renv`; `renv.lock` provided)  
  - Python (managed with `.yml` files; check the scripts and notebooks themselves as well) 

### Example compute environment
- **Sequencing QC and read alignment:** Executed on Stanford's [Sherlock HPC cluster](https://www.sherlock.stanford.edu/docs/tech/#resources) using SLURM scheduler. >2 TB disk space needed for raw FASTQs and intermediate files.
  - Genome: from NCBI genome assembly Nfu_20140520; RefSeq assembly [GCF_001465895.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001465895.1/)
  - Dependencies:  
    - Trim-galore (0.5.0)
    - STAR aligner (2.7.10b)  
    - Samtools (1.16.1)
    - Subread/featureCounts (2.0.6)
- **Primary transcriptomic analysis:** Executed on 2019 MacBook Pro (Intel)
	- Platform: x86_64-apple
	- Running under: macOS Sequoia 15.6.1
	- R version 4.3.3 (2024-03-01)
	- RStudio version 2024.09.1+394 (2024.09.1+394) 
- **Aging clocks:** Executed on Stanford's [SCG HPC cluster](https://med.stanford.edu/gbsc/services/scg-cluster/scg-cluster.html) using SLURM scheduler and Jupyter Notebook Server.
	- Jupyter module: Python 3.9, 4 cores, 8 GB memory, 1 GPU


## Software installation guide

- If **starting from FASTQs**, it is recommended to conduct count matrix generation on a HPC cluster. Often institutions have pre-installed bioinformatics tools on their clusters. Check your cluster's documentation for more information. [Example from Stanford's Sherlock HPC cluster](https://www.sherlock.stanford.edu/docs/software/list/).

- Clone the **repository**:

```bash
git clone https://github.com/emkcosta/KillifishAtlas.git
cd KillifishAtlas
```

- Set up **R environment**: 
	Run the following and it will reinstall the packages listed in `./environment/renv.lock`.

```bash
R -e "install.packages('renv'); renv::restore()"
```

- Set up **Python environment**:
YAML files are in `./environment`

**Elastic Net**
```bash
conda env create -f elasticnet_env.yml -n elasticnet
conda activate elasticnet
```

**BayesAge 2.0**
```bash
conda env create -f bayes_env.yml -n bayes
conda activate bayes
```


## Instructions for use

### Step 1: Prepare inputs
- Download FASTQs from SRA (study accession provided in the manuscript).  
- Align with STAR to produce BAMs.  
- Run featureCounts with consistent parameters (see manuscript) to generate per-sample `.featureCounts` files.  
- List them in a metadata file: `workflows/01_GetCounts/Input/ExperimentDesign_batch_Z01-F001_F002_combined_update.csv`.

### Step 2: Run R workflows
- **01_GetCounts**: Combine `.featureCounts` files → expression matrix  
- **02_QC**: Quality control of expression data  
- **03_DESeq2**: Differential expression analysis across age bins and sexes  
  - Data matrices (raw and DESeq2 normalized) are provided to reviewers and under GEO accession [GSE308970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308970). 
  - You can also start with this directory using raw counts.
  - Note: Scripts associated with DESeq2 pipeline are among the longest-running scripts for the R workflows (several hours-long run time as the for loop runs).
- **04_Correlation**: Find age-correlated genes
- **05_HClustering_GeneTrajectories**: Find clusters of gene trajectories over age 
- **06_GSEA_alltissue**: Perform GSEA for highly correlated genes
- **07_Conservation**: Run conservation analysis
- **08_scDeconvolution**: Run single cell deconvolution for kidney
- **09_HClustering_evaluation**: Evaluate the quality of clusters of gene trajectories
- **010_Batch_effects**: Evaluate batch effects
- **011_HyperGO_alltissue**: Perform Hypergeometric GO enrichment for highly correlated genes
- **012_Additional_Figures**: Additional figures

### Step 3: Run Python notebooks
##### Recapitulating Atlas Aging Clock Estimates
- Note: scripts associated with aging clocks will take several hours to run. Run time for PC-R < BayesAge 2.0 < Elastic net.
  - It is recommended to use a compute cluster. 
- Each subfolder contains a `README.md` describing inputs, outputs, and commands.

##### Running the Aging Clocks on Your Own Data
- Clone the github repository to your directory of interest.
- Note: scripts associated with aging clocks will take several hours to run. Run time for PC-R < BayesAge 2.0 < Elastic net.
  - It is recommended to use a compute cluster for this analysis.
  - Example bash script and scripts exported from Jupyter notebook to python script are provided in `other_manuscript_scripts/`. Improved SLURM-job-compatible scripts are forthcoming.
- See directory-specific README for more information.

### Reproducibility
 <TODO>
- Use the central config (`configs/default.yaml`) for consistent paths and parameters.
- Example commands are provided in each workflow README.  

## License

This repository is released under the [MIT License](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository).