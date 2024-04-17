import pandas as pd

PEPTIGATE_OUTPUT="outputs/ToT/predictions/peptides.faa"

#########################################################
## Cluster peptide sequences
#########################################################
"""
Clustering groups peptides with similar sequences.
This allows us to approximate whether the same (or similar) peptides were predicted in many
species. The ortholog group information from the ToT run also provides some of this evidence,
but different cleavage peptides could be predicted from different regions of sequence, so we
re-perform clustering here.

We selected 80% as a percent identity cutoff for clustering.
This is the percent identity cutoff frequently used prior to training a machine learning classifier
to make sure there is no shared sequences/signal between training and testing data sets.
We assume that 80% similarlity is a nice balance between clustering peptides that have shared
sequence and function and retaining diversity/heterogeneity of peptides in the data set.
"""


rule cluster_peptigate_protein_peptide_sequences:
    input:
        faa=PEPTIGATE_OUTPUT
    output:
        faa="outputs/analysis/clustering/all_peptides_0.8_rep_seq.fasta",
        tsv="outputs/analysis/clustering/all_peptides_0.8_cluster.tsv",
    params:
        out_prefix="outputs/analysis/clustering/all_peptides_0.8",
    conda:
        "envs/mmseqs2.yml"
    shell:
        """
        mkdir -p tmp
        mmseqs easy-cluster {input} {params.out_prefix} tmp --min-seq-id 0.8 
        """


#########################################################
## Predict anti-inflammatory bioactivity
#########################################################


rule unzip_autopeptideml_antiinflammatory_model:
    """
    The anti-inflammatory model is part of this github repository.
    Issue [#2](https://github.com/Arcadia-Science/2024-tick-sg-peptides-tsa/issues/2) documents how
    we built it.
    """
    input:
        "inputs/autopeptideml_antiinflammatory/apml_antiinflammatory_length15.zip",
    output:
        json="inputs/autopeptideml_antiinflammatory/apml_antiinflammatory_length15/config.json",
    params:
        outdir="inputs/autopeptideml_antiinflammatory/apml_antiinflammatory_length15/",
    shell:
        """
        mkdir -p {params.outdir}
        unzip {input} -d {params.outdir}
        """


rule download_autopeptideml_run_script_from_peptigate:
    """
        Note this won't work until the peptigate repo is public.
        I did this step by hand but I'm adding the rule as a placeholder.
        I think this is preferable over checking in the script here as well so that it doesn't become
        duplicated and need to be updated as the peptigate repo changes.
        """
    output:
        script="scripts/run_autopeptideml.py",
    shell:
        """
            curl -JLo {output} https://raw.githubusercontent.com/Arcadia-Science/peptigate/52a93c07cb46b950d9b379a8e3812d57c41b800a/scripts/run_autopeptideml.py
            """


rule predict_antiinflammatory_bioactivity_with_autopeptideml:
    input:
        script=rules.download_autopeptideml_run_script_from_peptigate.output.script,
        model=rules.unzip_autopeptideml_antiinflammatory_model.output.json,
        faa=PEPTIGATE_OUTPUT
    output:
        tsv="outputs/analysis/predict_antiinflammatory/autopeptideml_antiinflammatory_predictions.tsv",
    params:
        model_dir="inputs/autopeptideml_antiinflammatory/apml_antiinflammatory_length15/ensemble",
    conda:
        "envs/autopeptideml.yml"
    shell:
        """
        python scripts/run_autopeptideml.py \
            --input_fasta  {input.faa} \
            --model_folder {params.model_dir} \
            --model_name antiinflammatory \
            --output_tsv {output.tsv}
        """


#########################################################
## Compare against human peptides
#########################################################
"""
With this analysis, we hope to identify tick peptides that potentially mimic human peptides.
These peptides may have evolved to interact with humans and to control specific aspects of human
phsyiology (coagulation, itch, inflammation, etc.). 

Note that peptipedia does contain peptides from the human peptide atlas so some hits will be
redundant with those reported by peptigate. It is difficult to extract taxonomy information from
peptipedia and BLASTp is relatively inexpensive to run, so we implemented this specific comparison
step.

The human peptide atlas contains tryptic peptides (degradation products from enzymatic digestion.).
While we are not interested in these, we are relying on our tools to limit hits to these peptides
because we anticipate that peptigate will only predict bioactive peptides.
"""


rule download_human_peptide_atlas_peptide_sequences:
    output:
        faa="inputs/databases/humanpeptideatlas/APD_Hs_all.fasta",
    shell:
        """
        curl -JLo {output.faa} https://peptideatlas.org/builds/human/202401/APD_Hs_all.fasta
        """


rule make_diamond_blastdb_for_human_peptide_atlas:
    input:
        rules.download_human_peptide_atlas_peptide_sequences.output.faa,
    output:
        dmnd="inputs/databases/humanpeptideatlas/APD_Hs_all.dmnd",
    conda:
        "envs/diamond.yml"
    params:
        dbprefix="inputs/databases/humanpeptideatlas/APD_Hs_all",
    shell:
        """
        diamond makedb --in {input} -d {params.dbprefix}
        """


rule blastp_peptide_predictions_against_human_peptide_atlas:
    input:
        db=rules.make_diamond_blastdb_for_human_peptide_atlas.output.dmnd,
        faa=PEPTIGATE_OUTPUT
    output:
        tsv="outputs/analysis/compare_human/humanpeptideatlas_blastp_matches.tsv",
    params:
        dbprefix="inputs/databases/humanpeptideatlas/APD_Hs_all",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.faa} -o {output.tsv} --header simple \
         --outfmt 6 qseqid sseqid full_sseq pident length qlen slen qcovhsp scovhsp mismatch gapopen qstart qend sstart send evalue bitscore
        """


#########################################################
## Compare agaist known anti-pruritic peptides
#########################################################
"""
Unlike with anti-inflammatory peptides, there are relatively few examples of anti-pruritic peptides
in the literature. This means that are not enough examples of anti-pruritic peptides to train a
machine learning model that can detect this bioactivity. This set of rules uses BLASTp to compare
peptigate predicted peptides against anti-pruritic peptides.
"""


rule combine_antipruritic_peptide_sequences:
    input:
        "inputs/antipruritic_peptides/calcitonin_gene-related_peptide.faa.gz",
        "inputs/antipruritic_peptides/dynorphin.faa.gz",
        "inputs/antipruritic_peptides/tachykinin-4.faa.gz",
        "inputs/antipruritic_peptides/votucalis.faa.gz",
        "inputs/antipruritic_peptides/ziconotide.faa.gz",
    output:
        faa="inputs/antipruritic_peptides/antipruritic_peptides.faa",
    shell:
        """
        cat {input} | gunzip > {output}
        """


rule make_diamond_blastdb_for_antipruritic_peptides:
    input:
        rules.combine_antipruritic_peptide_sequences.output.faa,
    output:
        dmnd="inputs/databases/antipruritic_peptides/antipruritic_peptides.dmnd",
    conda:
        "envs/diamond.yml"
    params:
        dbprefix="inputs/databases/antipruritic_peptides/antipruritic_peptides",
    shell:
        """
        diamond makedb --in {input} -d {params.dbprefix}
        """


rule blastp_peptide_predictions_against_antipruritic_peptides:
    input:
        db=rules.make_diamond_blastdb_for_antipruritic_peptides.output.dmnd,
        faa=PEPTIGATE_OUTPUT
    output:
        tsv="outputs/analysis/predict_antipruritic.tsv",
    params:
        dbprefix="inputs/databases/antipruritic_peptides/antipruritic_peptides",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.faa} -o {output.tsv} --header simple \
         --outfmt 6 qseqid sseqid full_sseq pident length qlen slen qcovhsp scovhsp mismatch gapopen qstart qend sstart send evalue bitscore
        """


#########################################################
## Compare against tick salivary gland predicted peptides
#########################################################
"""
This set of rules compares peptides predicted from the ToT itch-associated proteins to peptides
predicted from tick salivary gland transcriptomes.
We expect peptides that have anti-itch functions to be expressed in tick salivary glands as this
would allow the peptides to reach host skin where itch suppression is active.
While we can confirm that peptides are expressed in the salivary gland, this analysis does not:
1. Test whether the peptide is expressed in other tissues
2. Definitively state that the peptide is not expressed in the salivary gland; the transcriptomes
   we predicted peptides from are incomplete, so just because we don't see something doesn't mean
   it isn't there. 
"""


rule make_diamond_blastdb_for_tsa_sg_peptide_predictions:
    input: "inputs/tsa_sg_peptides.faa.gz"
    output:
        dmnd="inputs/databases/tsa_sg_peptides/tsa_sg_peptides.dmnd",
    conda:
        "envs/diamond.yml"
    params:
        dbprefix="inputs/databases/tsa_sg_peptides/tsa_sg_peptides",
    shell:
        """
        diamond makedb --in {input} -d {params.dbprefix}
        """


rule blastp_peptide_predictions_against_tsa_sg_peptide_predictions:
    input:
        db=rules.make_diamond_blastdb_for_tsa_sg_peptide_predictions.output.dmnd,
        faa=PEPTIGATE_OUTPUT
    output:
        tsv="outputs/analysis/compare_tsa_sg/tsa_sg_peptides_blastp_matches.tsv",
    params:
        dbprefix="inputs/databases/tsa_sg_peptides/tsa_sg_peptides",
    conda:
        "envs/diamond.yml"
    shell:
        """
        diamond blastp -d {params.dbprefix} -q {input.faa} -o {output.tsv} --header simple \
         --outfmt 6 qseqid sseqid full_sseq pident length qlen slen qcovhsp scovhsp mismatch gapopen qstart qend sstart send evalue bitscore
        """

#########################################################
## Collect outputs
#########################################################


rule all:
    default_target: True
    input:
        rules.blastp_peptide_predictions_against_human_peptide_atlas.output.tsv,
        rules.blastp_peptide_predictions_against_antipruritic_peptides.output.tsv,
        rules.blastp_peptide_predictions_against_tsa_sg_peptide_predictions.output.tsv,
        rules.predict_antiinflammatory_bioactivity_with_autopeptideml.output.tsv,
        rules.cluster_peptigate_protein_peptide_sequences.output.tsv,
