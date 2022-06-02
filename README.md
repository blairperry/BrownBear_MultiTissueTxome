# BrownBear_MultiTissueTxome
This repo contains code for data processing and analysis used in our publication of a multi-tissue transcriptomic dataset for hibernating brown bears.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
### Newly generated data:
- RNA-seq from multiple tissues from two hibernating brown bears
	- Raw data available at NCBI Bioproject [PRJNA835146](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA835146)
		- Tissues: lung (left and right), heart ventricle (left and right), heart atrium (left and right), small intestine, kidney (medulla and cortex), gall bladder, subcutaneous adipose, stomach, gastrocnemius muscle, liver, skin, and spleen. 
	- The following count tables are available in this repo:
		- [Raw gene-level counts](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/data/raw_counts/a_geneLevelQuants/bearNecropsy_geneLevelCounts_10.12.21_bwp.cleaned.txt)
		- [TPM normalized gene-level counts](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/rnaseq_normalization/bear_TPMnorm_geneLevelCounts.csv)
		- [TPM normalized transcript-level counts](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/isoformExpression/Transcripts_DTUScaledTPMCounts.csv.gz)
		- [Normalized orthogroup-level counts for bear and human](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/rnaseq_normalization/combined_oglevel_tmmCounts.txt.gz)

### Previously generated data:
- Brown bear mRNA-seq data (from [Jansen et al. 2019](https://www.nature.com/articles/s42003-019-0574-4))
  - Raw data available at NCBI BioProject: [PRJNA413091](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA413091)
    - Three physiological stages: active, hyperphagia, hibernation
    - Three metabolic tissues: liver, muscle, adipose
- Iso-seq Transcriptome Assembly (from [Tseng et al. 2021](https://academic.oup.com/g3journal/article/12/3/jkab422/6472356))
  - Raw data available at NCBI BioProject: [PRJNA727613](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA727613)
  - Assembly and annotation available at: https://github.com/jokelley/brownbear-isoseq-act-hib

## Software
The following software was used for data processing and analysis:
- AGAT v0.8.0 (extracting and translating longest CDS per gene)
- OrthoFinder v2.5.4 (identifying orthogroups between human and brown bear)
- Trim Galore v0.4.2 (quality trimming raw RNA-seq reads)
- STAR v2.7.6a (mapping RNA-seq reads to bear genome assembly)
- MultiQC v1.11 (generating reports and summarizing analysis results)
- featureCounts v1.6.3 (quuantifying gene-level expression)
- Kallisto v0.46.2 (mapping and quantifying RNA-seq data against bear transcriptome assembly)
- R packages
	- edgeR v3.36.0 (for TMM normalization of count data)
	- ggplot2 v3.3.5 (plotting analysis results)
	- pheatmap v1.0.12 (plotting heatmaps)
	- TissueEnrich v1.14.0 (identifying genes with tissue-specific expression)
	- DTUrtle v1.0.2 (analyses of differential transcript usage)


## Contents
1. [Orthogroup identification](#1-orthogroup-identification)
2. [Gene-level expression quantification](#2-gene-level-expression-quantification)
3. [Transcript-level expression quantification](#3-transcript-level-expression-quantification)
4. [Incorporation of existing human gene expression data from GTEX](#4-incorporation-of-existing-human-gene-expression-data-from-gtex)
5. [Orthogroup-level normalization of gene expression data](#5-orthogroup-level-normalization-of-gene-expression-data)
6. [Evaluation of orthogroup-level normalization with known housekeeping genes](#6-evaluation-of-orthogroup-level-normalization-with-known-housekeeping-genes)
7. [Analysis of differential transcript usage](#7-analysis-of-differential-transcript-usage)
8. [Tissue-specific expression analysis](#8-tissue-specific-expression-analysis)
9. [Analyses of candidate genes](#9-analyses-of-candidate-genes)


---
### 1. Orthogroup identification
Note: These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.
- Assemblies used:
	- Brown bear: [NCBI Assembly ASM358476v1](https://www.ncbi.nlm.nih.gov/assembly/GCF_003584765.1/)
	- Human: NCBI Assembly [GRCh38.p10](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.36/)

#### Extract longest CDS per gene using AGAT
```bash
agat_sp_keep_longest_isoform.pl -gff GCF_003584765.1_ASM358476v1_modified.gtf -o Uarctos.longestCDS.gtf
```

#### Translate longest CDS using AGAT
```bash
agat_sp_extract_sequences.pl -g Uarctos.longestCDS.gtf -f GCF_003584765.1_ASM358476v1_genomic.fna -p -o Uarctos.longestCDS.prots.faa
```

#### Run OrthoFinder with default parameters
```bash
# Run OrthoFinder
orthofinder -t 12 -a 12 -f [directory/with/longest/protein/fasatas] -o ./orthofinderResults

# Parse OrthoFinder output using custom python script - generates table used to annotate GFF with orthogroup info
python2 ./utility_scripts/parseOrthogroups.py ./Orthogroups.tsv homSap 3 # output file: Orthogroups.homSap.tsv
python2 ./utility_scripts/parseOrthogroups.py ./Orthogroups.tsv ursArc 5 # output file: Orthogroups.ursArc.tsv 

# Add orthogroup information to gff file for downstream analyses using custom python script 
python2 ./utility_scripts/annotateGFF.py ./input/longest_cds/homSap.longCDS.gff Orthogroups.homSap.tsv ./orthoAnnotatedGffs/homSap.longCDS.orthogroups.gff

python2 ./utility_scripts/annotateGFF.py ./input/longest_cds/ursArc.longCDS.gff Orthogroups.ursArc.tsv ./orthoAnnotatedGffs/ursArc.longCDS.orthogroups.gff
```
Utility scripts used in above code block can be found [here](https://github.com/blairperry/BrownBear_MultiTissueTxome/tree/master/utility_scripts).

The following script was then used to parse orthoFinder results to generate lookup tables linking bear and human orthologs.
- Link to Rscript: [idConversionTables_10.20.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/orthoFinder/id_conversion/idConversionTables_10.20.21.R)

---
### 2. Gene-level expression quantification
Note: These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.

#### Quality trimming of raw reads using Trim Galore
```bash
trim_galore --paired -q 20 --fastqc --fastqc_args "--noextract --nogroup --outdir naseq_qtrim_10.07.21/fastqc" --stringency 5 --illumina --length 50 -o rnaseq_qtrim_10.07.21 --clip_R1 12 --clip_R2 12 [read1] [read2]
```

#### Mapping trimmed reads to bear genome with STAR and quantifying counts with featureCounts
```bash
# Run STAR
STAR --genomeDir star_reference/ --runThreadN 8 --outFilterMultimapNmax 1 --twopassMode Basic --sjdbGTFfile GCF_003584765.1_ASM358476v1_modified.gtf --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STARmapped_10.07.21/[file_name] --readFilesIn [read1] [read2]

# Summarize STAR results w/ multiQC
multiqc STARmapped_10.07.21/*final.out -i bearMultiTissueRNA_STARoutput

# Convert gff to gtf w/ AGAT
agat_convert_sp_gff2gtf.pl --gff ursArc.longCDS.orthogroups.gff -o ursArc.longCDS.orthogroups.gtf

# Convert gtf to SAF using custom python script 
python2 ./utility_scripts/gtf_to_saf.allGeneInfo.py ursArc.longCDS.orthogroups.gtf exon
# Output file is titled ursArc.longCDS.orthogroups.allGeneInfo.saf

# Quantify gene-level counts using featureCounts
featureCounts -p -F 'SAF' -T 8 -t exon -g gene_id -a ursArc.longCDS.orthogroups.allGeneInfo.saf -o ./a_geneLevelQuants/bearNecropsy_geneLevelCounts_10.12.21_bwp.txt /scratch/user/blair.perry/20211004_112224/STARmapped_10.07.21/*.sortedByCoord.out.bam

```
Utility scripts used in above code block can be found [here](https://github.com/blairperry/BrownBear_MultiTissueTxome/tree/master/utility_scripts)

---
### 3. Transcript-level expression quantification
Note: These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.

```bash
# Index transcriptome with Kallisto
kallisto index -i transcripts.idx ./new_merge.fa

# Quantify transcript-level counts with Kalliso
kallisto quant -i transcripts.idx --rf-stranded -o kallisto_quants_Feb2022/[file_name] -b 100 -t 5 [read1] [read2]
```

---
### 4. Incorporation of existing human gene expression data from GTEX
Downloading raw gene-level expression from GTEX v8.
```bash
# Download count table
wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz

# Download sample info 
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt 
wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
```

GTEX data was then filtered to relevant tissues, and 20 samples per tissue were randomly selected (to reduce computational demands in downstream analyses). 
	- Subsampling of human gene expression data was conducted using this Rscript: [gtex_sampleParsing_12.01.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/human_rnaseq/gtex_sampleParsing_12.01.21.R)
	
---
### 5. Orthogroup-level normalization of gene expression data
In order to compare bear and human gene expression, we normalized expression data at the level of orthogroup rather than at the level of genes. In theory, this will account for the fact that 1) different species may have different numbers of orthologous genes (i.e., 1 ortholog in humans but 2 in bear) and 2) that gene lengths of orthologous genes may differ between species. Please see the manuscript text for additional details on the motivation and limitations of this approach. 

This process involves the following steps:
1. Remove genes with no expression (i.e., raw count == 0) in all human and bear samples
2. Separately for each species, sum raw gene-level counts for genes belonging to the same orthogroup. This produces a set of raw orthogroup-level counts for each species. 
3. Separately for each species, sum the length of all exons for all genes belonging to the same orthogroup. This produces a set of orthogroup "lengths" for each species.
4. Separately for each species, normalize raw orthogrou-level counts to produce transcripts per million (TPM) counts.
5. Merge TPM count tables for both species into a single table.
6. Using merged table with TPM counts for both species, normalize again using TMM method in edgeR. 

Commands used for the above steps and downstream evaluation and analysis can be found in the following script.
- Link to Rscript: [_OrthoLevelCountsAndNorms_01.13.22.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/rnaseq_normalization/_OrthoLevelCountsAndNorms_01.13.22.R)

---
### 6. Evaluation of orthogroup-level normalization with known housekeeping genes

The following script was used to assess orthogroup-level normalized counts using known housekeeping genes (from [Eisenberg and Levanon 2013](https://pubmed.ncbi.nlm.nih.gov/23810203/))
- Link to Rscript: [HousekeepingOGgroups_11.29.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/candidate_genes/housekeeping/HousekeepingOGgroups_11.29.21.R)

---
### 7. Analysis of differential transcript usage

The following scripts were used to quantify and plot differential transcript usage (DTU) between tissues.
- [DTUrtle_FilePrep_01.04.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/isoformExpression/DTUrtle_FilePrep_01.04.21.R) - prepare files and data for pairwise DTU analyses
- [AllTissues_dTurtle_02.04.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/isoformExpression/7_DTUanalyses/AllTissues_dTurtle_02.04.21.R) - run all pairwise DTU analyses (this was run on Kamiak)
- [DTurtleResultParser.py](https://github.com/blairperry/BrownBear_MultiTissueTxome/tree/master/analyses/isoformExpression#:~:text=DTurtleResultParser.py) - used to parwse result files from above script into simpler format
- [dturtle_MultipleTissueResults_01.06.21.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/isoformExpression/dturtle_MultipleTissueResults_01.06.21.R) - parsing and plotting DTU analysis results

---
### 8. Tissue-specific expression analysis

The following script was used to quantify and plot genes and transcripts with tissue-specific expression.
- Link to Rscript: [TissueSpecificTranscripts_01.10.22.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/tissueSpecificExpression/TissueSpecificTranscripts_01.10.22.R)

---
### 9. Analyses of candidate genes

The following script was used to analyze and plot expression of candidate genes across bear tissues.
- Link to Rscript: [Metabolism_CandGenes_04.25.22.R](https://github.com/blairperry/BrownBear_MultiTissueTxome/blob/master/analyses/candidate_genes/CommBio_MetabolismGenes_Fig4and5/Metabolism_CandGenes_04.25.22.R)