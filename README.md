# BrownBear_MultiTissueTxome
This repo contains code for data processing and analysis used in our publication of a multi-tissue transcriptomic dataset for hibernating brown bears.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
### Newly generated data:
- RNA-seq from multiple tissues from two hibernating brown bears
	- Raw data available at NCBI Bioproject [PRJNA835146](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA835146)
		- Tissues: lung (left and right), heart ventricle (left and right), heart atrium (left and right), small intestine, kidney (medulla and cortex), gall bladder, subcutaneous adipose, stomach, gastrocnemius muscle, liver, skin, and spleen. 
	- Count tables are available in the count_tables directory of this repo
		- ADD FILENAME - raw gene-level counts
		- ADD FILENAME - raw transcript-level counts
		- ADD FILENAME - normalized gene-level counts
		- ADD FILENAME - normalized transcript-level counts
		- ADD FILENAME - normalized orthogroup-level counts for bear and human

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
1. Orthogroup identification

---
### 1. Orthogroup identification
These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.
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

---
### 2. Gene-level expression quantification
These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak.

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

