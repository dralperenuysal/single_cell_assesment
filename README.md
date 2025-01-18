# Single Cell Assesment Repository
Single cell assesment for a serous carcinoma transcriptomic data.

A single-cell dataset, "Chemotherapy induces myeloid-driven spatial
T-cell exhaustion in ovarian cancer", was analyzed in this assesment.

## INDEX
  - [output](output/) :Storing the output files and plots.
  - [scripts](scripts/) :Directory, Storing the codes and notebooks.
    - [Preprocesing Task Notebook](scripts/preprocessed.md)
    - [Dimensionality Reduction Task Notebook](scripts/dimension_reduction.md)
    - [Cell Type Annotation Task Notebook](scripts/annotation.md)
    - [Trajectory Inference and Pseudotime Analysis](scripts/trajectory.md)

## Overview

This repository provides a comprehensive analysis of a single-cell RNA
sequencing (scRNA-seq) dataset, focusing on trajectory inference and pseudotime
analysis. The primary aim is to explore the tumor microenvironment (TME)
dynamics in high-grade serous ovarian cancer (HGSC) through advanced
computational tools. This analysis leverages a publicly available dataset,
GSE266577, to examine the effects of chemotherapy on the spatial organization
and immune dynamics of the TME.

## Dataset Information
Dataset: GSE266577

Title: Chemotherapy induces myeloid-driven spatial T-cell exhaustion in ovarian
cancer

Organism: Homo sapiens

Experiment Type: Expression profiling by high-throughput sequencing

Summary:
This dataset characterizes the spatial TME of HGSC at the single-cell level
using scRNA-seq from 48 tumor or ascites samples of 29 HGSC patients. The data
includes paired scRNA-seq samples from chemo-naïve and post-neoadjuvant
chemotherapy (IDS) conditions, focusing on spatial immune dynamics and T-cell
exhaustion mechanisms.

Study Design:

Single-cell gene expression profiles were generated from tumor and ascites
samples.
Data includes chemo-naïve and IDS conditions for 22 patients.
Raw data is available in the European Genome-phenome Archive (EGA) due to
privacy concerns.

Publication:

Citation: Launonen IM, Niemiec I, Hincapié-Otero M, Erkan EP et al.
Chemotherapy induces myeloid-driven spatially confined T cell exhaustion in
ovarian cancer. Cancer Cell 2024 Dec 9;42(12):2045-2063.e10. PMID: 39658541
DOI: https://doi.org/10.1016/j.ccell.2024.11.005