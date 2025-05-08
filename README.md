# Predicting peptides from chelicerate proteins associated with host detection suppression 

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

This repository predicts peptides from chelicerate species that are associated with host detection suppression.
In this context, we refer to suppression of the triad of inflammation, pain, and itch as "host detection suppression." 

## Installation and Setup

This repository uses Snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n ticktides --file envs/dev.yml
conda activate ticktides
```

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the [envs/](./envs/) directory. Snakemake itself is installed in the main development conda environment as specified in the [dev.yml](./envs/dev.yml) file.

To start the pipeline, run:

```{bash}
snakemake -s analyze-peptigate-outputs.snakefile --software-deployment-method conda -j 8
```

## Data

This repository analyzes the outputs of four previous analyses to predict peptide sequences from proteins that are associated with host detection suppression.
These analyses were previous pilots or work done by others at Arcadia.
The four upstream analyses are:
1. [**protein-data-curation**](https://github.com/Arcadia-Science/protein-data-curation/): This repository provides a workflow to download, process, and annotate genomic and transcriptomic data that can then be fed to NovelTree (see next point). The input sequences are all from chelicerate species and are documented [here](https://zenodo.org/records/14113178/files/chelicerate-samples.tsv?download=1). This repository also annotated these sequences, and that metadata is now incorporated into the file [2024-06-26-top-positive-significant-clusters-orthogroups-annotations.tsv.gz](inputs/2024-06-26-top-positive-significant-clusters-orthogroups-annotations.tsv.gz). 
2. **Chelicerate NovelTree run**: [NovelTree](https://github.com/Arcadia-Science/noveltree) applied to the curated chelicerate dataset. The input dataset to NovelTree is protein sequences (either genome gene annotations or open reading frames predicted from transcriptomes). NovelTree applies evolutionary analyses to identify proteins that are "novel" under models of speciation, loss, or transfer. 
3. [**Chelicerate trait mapping**](https://github.com/Arcadia-Science/2024-chelicerate-phylogenomics): trait mapping applied to the results of the chelicerate NovelTree run to identify proteins (from clusters of orthologous groups) that are associated with the host detection-suppression trait. The protein sequences associated with host detection suppression are recorded in the file [`2024-06-26-top-positive-significant-clusters-orthogroups-proteins.fasta.gz`](./inputs/2024-06-26-top-positive-significant-clusters-orthogroups-proteins.fasta.gz).
4. [**Tick salivary gland transcriptome peptides**](https://github.com/Arcadia-Science/2024-tick-sg-peptides-tsa/): peptide sequences predicted from publicly available tick salivary gland transcriptomes. The results from this analysis are in the file [tsa_sg_peptides.faa.gz](inputs/tsa_sg_peptides.faa.gz). 

## Overview

Using the sequences of proteins predicted to be associated with the host detection-suppression trait, we predicted peptide sequences from these proteins.
We define peptides as any protein sequence of less than or equal to 100 amino acids in length.

The peptigate pipeline predicted peptides from transcriptome assemblies.
It had two prediction modules that targeted small open reading frames (sORFs) and cleavage peptides, respectively.
Because we had access to protein sequences but not to transcriptome assemblies, we modified our peptide prediction approach to retrieve as many peptides as possible given our input data.
We:
1. **sORF prediction**: filtered the host detection-suppression-associated proteins to those that were 100 amino acids or less. This was where we suffered the most potential loss in our predictive capabilities -- we think we'd have many more potential sORF sequences if we had transcriptomes as input instead of protein sequences. However, these would require a new trait mapping analysis.
2. **Cleavage peptide prediction**: ran the cleavage peptide prediction portion of the peptigate pipeline on protein sequences longer than 100 amino acids.
3. **Peptide annotation**: annotated the predicted peptides with the annotation and analysis portion of the peptigate pipeline.

We ran peptigate from commit `148823239aad41a8f03da37f5499b00c8a79de40` with the following command:

```
snakemake -s protein_as_input.snakefile --software-deployment-method conda -j 1 -k --configfile input_configs/tot_protein_peptigate_config.yml
```

where the `tot_protein_peptigate_config.yml` looked like:
```
input_dir: "inputs/"
output_dir: "outputs/ToT_20240626"
orfs_amino_acids: "input_data/ToT_20240626/2024-06-26-top-positive-significant-clusters-orthogroups-proteins.fasta.gz"
```

With these peptide predictions, we then proceeded with the rest of the analysis.

We ran additional annotation and comparison analysis with the snakefile in this repository:

```
snakemake -s analyze-peptigate-outputs.snakefile --software-deployment-method conda -j 8
```

This snakefile:
* Compares (BLASTp) peptide predictions against peptides in the Human Peptide Atlas. The idea is that if a peptide resembles a human peptide, it may have evolved to mimic a human process or to interact with human molecules.
* Compares (BLASTp) peptide predictions against peptides predicted in tick salivary glands. We want to target peptides in tick salivary glands because this is where we expect to see most molecules that have evolved to influence human biology.
* Compares (BLASTp) peptide predictions against peptides with known anti-pruritic effects. We only have 5 peptides (see the [anti-pruritic peptides folder](./inputs/antipruritic_peptides)) that are known to have anti-host detection activity, but we compare against them.
* Clusters (MMseqs2, 80% identity) peptide predictions to determine how similar different peptides are to each other. 80% is somewhat arbitrary, but it's commonly used in machine learning algorithms as a cutoff for shared information between proteins.
* Predicts anti-inflammatory bioactivity (AutoPeptideML) for peptide predictions. The model that does this has 70% accuracy. We expect there to be some overlap between inflammation suppression and host detection suppression, so we include this information.


Last, we ran the following [notebooks](./notebooks) to combine and filter the results to produce the most promising candidates for peptides that may suppress host detection. 
We ran the notebooks using the `envs/tidyjupyter.yml` environment.

* [20240626-00-explore-peptide-results.ipynb](./notebooks/20240626-00-explore-peptide-results.ipynb): Initial exploration of results and metadata. This notebook is messy and was more of a scratch space, but I included it because it has some useful visualizations.
* [20240626-01-combine-traitmapping-info-with-peptide-predictions.ipynb](./notebooks/20240626-01-combine-traitmapping-info-with-peptide-predictions.ipynb): combines information from trait mapping analysis and peptide predictions and summarizes information about peptide predictions per orthogroup.
* [20240626-02-combine-peptide-predictions-with-metadata.ipynb](./notebooks/20240626-02-combine-peptide-predictions-with-metadata.ipynb): combines peptide predictions with peptide metadata.
* [20240626-03-peptides-into-pools.ipynb](./notebooks/20240626-03-peptides-into-pools.ipynb): documents the rationale for selecting peptides for synthesis.
* [20240628-04-pub-numbers-and-figures.ipynb](./notebooks/20240628-04-pub-numbers-and-figures.ipynb): additional analysis and summaries for results presented in the pub about this project.

## Results 

We started with 3,534 input protein sequences from 78 orthogroups.
We initially predicted 631 peptides (602 distinct sequences) in 37 orthogroups.
We applied the following initial filters (see [this notebook](./notebooks/20240626-02-combine-peptide-predictions-with-metadata.ipynb)):
* Filtered (removed) propeptides predicted by DeepPeptide. DeepPeptide uses the [UniProt definition of a propeptide](https://www.uniprot.org/help/propep), a part of a protein that's cleaved during maturation or activation. Once cleaved, a propeptide generally has no independent biological function. This reduced the number of peptide predictions to 308 (305 distinct sequences) in 32 orthogroups.
* Filtered (removed) peptides in orthogroups where no peptide had a hit to a predicted peptide from a [tick salivary gland transcriptome](https://github.com/Arcadia-Science/2024-tick-sg-peptides-tsa/). We want to target things expressed in the tick salivary gland because they're more likely to be biologically active in host detection suppression. However, we don’t know if the tick salivary gland transcriptomes we worked with are complete (many are heavily filtered) so we relaxed this filter to function at the orthogroup level. This filter reduced the number of peptides down to 281 (278 distinct peptide sequences). These peptides were from 13 orthogroups and 206 clusters (MMseqs2 80% identity).

Using these sequences, we then further filtered to select the sequences most likely to suppress host detection and easiest to work with for experimental validation (see [this notebook](./notebooks/20240626-03-peptides-into-pools.ipynb)):
* Filtered to orthogroups where the majority of proteins had a predicted peptide. We were most interested in orthogroups where the majority of predicted peptides were of the same type (sORF or cleavage), although we did not use this as a strict filtering class.
* Filtered to orthogroups where the majority of peptides (sORF) or parent proteins (cleavage) had signal peptides. We only kept peptides with signal peptides.

These filters gave us a set of 88 peptides from 3 orthogroups.
We then selected peptides within each of these orthogroups to synthesize.
We made selections based on [ease of synthesis & solubility](https://www.genscript.com/tools/peptide%2danalyzing%2dtool), similar peptide expressed in tick salivary gland transcriptome, and similarity it other peptides in the group.
We ended up with 12 peptides (5 from OG0008102, 3 from OG0001774, 4 from OG0000880).

### Compute Specifications

* peptigate pipeline: Ran on an AWS EC2 instance type `g4dn.2xlarge` running AMI Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 20.04) 20240122 (AMI ID ami-07eb000b3340966b0). Note the pipeline runs many tools that use GPUs.
* [`analyze-peptigate-outputs.snakefile`](./analyze-peptigate-outputs.snakefile): Ran on an AWS EC2 instance type `g4dn.2xlarge` running AMI Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 20.04) 20240122 (AMI ID ami-07eb000b3340966b0). One tool in the pipeline, AutoPeptideML, uses a GPU.
* Notebooks: Ran on a MacBook Pro (2021, M1) using operating system Ventura 13.4 and with 64GB of RAM.   

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
