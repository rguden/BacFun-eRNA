# BacFun-RNA

### Author: Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
**Contact:** rodgeemae.guden@ugent.be  

---

## Overview

**BacFun-RNA** is a modular Bash workflow for functional metatranscriptomics analysis of microbial communities. This pipeline processes RNA-Seq data to generate clean mRNA reads, performs taxonomic profiling with MetaPhlAn, functional profiling with HUMAnN3, and downstream gene and pathway analyses.

The workflow is optimized for environmental or host-associated microbial communities and provides bacterial and functional gene abundance tables ready for downstream statistical analysis.

---

## Key Features

✔️ Quality control with FastQC & MultiQC  
✔️ Adapter trimming with Fastp  
✔️ rRNA removal with SortMeRNA  
✔️ MetaPhlAn taxonomic profiling (rRNA-depleted reads)  
✔️ HUMAnN3 functional profiling (genes, pathways)  
✔️ Feature regrouping to KEGG, EggNOG, Pfam, EC, GO  
✔️ Normalization and stratification of abundance tables  
✔️ KEGG pathway enrichment and DESeq2-based differential abundance analysis  

---

## Dependencies

| Tool         | Purpose                               | Installation link                           |
|--------------|---------------------------------------|----------------------------------------------|
| `fastqc`     | Read quality control                 | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| `multiqc`    | Aggregate QC reports                 | [MultiQC](https://multiqc.info/)            |
| `fastp`      | Adapter trimming                     | [fastp](https://github.com/OpenGene/fastp)  |
| `sortmerna`  | rRNA filtering                       | [SortMeRNA](https://github.com/biocore/sortmerna) |
| `metaphlan`  | Taxonomic profiling                  | [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) |
| `humann`     | Functional profiling (HUMAnN3)       | [HUMAnN3](https://github.com/biobakery/humann) |
| `R`          | Downstream OTU processing & DESeq2   | [R Project](https://www.r-project.org/) |

---

## Usage

### 1. Make the script executable

```bash
chmod +x metatranscriptomics_mrna_pipeline.sh
```

### 2. Run the workflow

Specify steps to run as a comma-separated list:

```bash
bash metatranscriptomics_mrna_pipeline.sh -s qc,trim_fastp,postqc_trim_fastp,download_db,index_db,remove_rrna_unmerged,postqc_sortmerna,summary_sortmerna,download_humann3_db,merge_mrna_reads,run_metaphlan,run_humann3,join_humann3_tables,normalize_humann3,split_humann3,download_utility_mapping,regroup_genefamilies_raw,rename_features_raw,split_renamed_raw,regroup_genefamilies,rename_features,analyze_deseq2_kegg_enrichment
```

You can run steps individually or in groups.

---

## Output

The workflow generates:

- Quality control reports before and after trimming  
- Trimmed RNA-Seq reads  
- rRNA-depleted mRNA reads  
- Taxonomic profiles from MetaPhlAn  
- HUMAnN3 gene family and pathway abundance tables  
- KEGG, EggNOG, EC, Pfam, GO regrouped abundance tables  
- Normalized and unnormalized feature abundance tables  
- KEGG enrichment and differential abundance results (DESeq2)  

---

## Notes

- Paths to input, output, and database directories are set at the top of the script. Adjust as needed.  
- Assumes paired-end RNA-Seq FASTQ files with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes.  
- MetaPhlAn and HUMAnN3 databases must be pre-downloaded.  
- Optional email notifications are implemented.  

---

## Citation & License

© 2025 Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
All rights reserved.

For permission to reuse or adapt this workflow, contact the author.

---

## Acknowledgments

Developed for microbial functional profiling within the ORCHESTRA project.

---

## Questions & Issues

For questions or issues, contact the author or open an Issue in the repository.
