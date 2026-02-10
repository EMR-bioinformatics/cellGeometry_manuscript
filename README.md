# cellGeometry\_manuscript

This repository contains scripts for the cellGeometry manuscript



## Repository Structure

```
.
├── 01_Workflow_visualisations
│   ├── celltypist_workflow.R
│   └── plotly sphere demo.R
├── 02_Simulated_benchmarking
│   ├── No_noise
│   │   ├── AMP_dirichlet_plot.R
│   │   ├── AMP_dirichlet.R
│   │   ├── AMP_NMF.R
│   │   ├── AMP_samplesize.R
│   │   ├── Brain_dirichlet_plot.R
│   │   ├── Brain_dirichlet.R
│   │   ├── Brain_NMF.R
│   │   ├── celltypist_dirichlet_plot.R
│   │   ├── celltypist_dirichlet.R
│   │   ├── celltypist_NMF.R
│   │   ├── combined_metric_plot.R
│   │   ├── combined_NMF_metric_plot.R
│   │   ├── specificity_plot.R
│   │   ├── Tabula_cs.R
│   │   ├── Tabula_dirichlet_plot.R
│   │   ├── Tabula_dirichlet.R
│   │   ├── Tabula_NMF.R
│   │   └── Tabula_similarity.R
│   └── With_noise
│       ├── AMP_dirichlet_noise_lin_Music.R
│       ├── AMP_dirichlet_noise_plot.R
│       ├── AMP_dirichlet_noise.R
│       ├── AMP_dirichlet_noise_sqrt.R
│       ├── Brain SE noise.R
│       └── Tabula SE noise.R
├── 03_Real_world_blood
│   ├── Blood_markers.R
│   ├── CellTypist PEAC SE check plot.R
│   ├── CellTypist PEAC SE check.R
│   ├── Music_tabula.R
│   ├── PEAC_bld_plot.R
│   └── PEAC_bld.R
└── 04_Real_world_synovium
    ├── DWLS_reprex
    │   ├── DWLS_reprex.R
    │   ├── R4RA_AMP1_DWLS.rds
    │   ├── R4RA_counts_anon.rds
    │   └── Signature.rdata
    ├── Pathotype_redef.R
    ├── PEAC_AMP.R
    ├── R4RA_AMP_benchmark_plot.R
    ├── R4RA_AMP_benchmark.R
    └── STRAP_AMP.R

```

## Data availability 

* Single-cell datasets for AMP phase 1 synovium are available from ImmPort accession code [SDY998](https://www.immport.org/shared/study/SDY998). 

* Single-cell datasets for Cell Typist are available from  [https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3](https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3)

* Single-cell datasets for Tabula Sapiens are available from [https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) 

* Single-cell datasets for the Human Brain Cell Atlas is available from [https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443](https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443)

* The PEAC RA blood and synovium bulk RNA-Seq datasets, R4RA and STRAP bulk RA synovium RNA-Seq datasets are available from ArrayExpress under accession codes [E-MTAB-6141](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-6141), [E-MTAB-11611](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-11611) and [E-MTAB-13733](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-13733).
