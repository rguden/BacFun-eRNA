    #!/bin/bash
set -e  # Stop on error

# =====================================================
# BacFun-RNA workflow: analysing mRNA reads for functional analyses of bacterial communities
# Author: Rodgee Mae Guden

# Copyright © 2025 Rodgee Mae Guden. Marine Biology Research Group, Ghent University
# All rights reserved. This script and its associated workflow are the intellectual property of the author. No part of this code may be reproduced, distributed, or transmitted in any form or by any means, including photocopying, recording, or other electronic or mechanical methods, without the prior written permission of the author, except in the case of brief quotations used in academic or research contexts with proper attribution.

# For permissions, contact: rodgeemae.guden@ugent.be

# =====================================================

# Parameters
THREADS=10
EMAIL="rodgeemae.guden@ugent.be"
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# Input and output directories
INPUT_DIR="/home/rguden/Documents/ORCHESTRA/sequences/RNA_raw_data"
OUTPUT_DIR="/home/rguden/Documents/ORCHESTRA/RNA_output"
LOG_DIR="$OUTPUT_DIR/logs"
LOG_FILE="$LOG_DIR/mrna_humann.log"

# SortMeRNA database directory
SORTMERNA_DB_DIR="/home/rguden/Documents/databases/sortmerna"

# Ensure directories exist
mkdir -p "$OUTPUT_DIR/quality_control"
mkdir -p "$OUTPUT_DIR/trimming/fastp"
mkdir -p "$OUTPUT_DIR/trimming/fastqc/fastp"
mkdir -p "$OUTPUT_DIR/sortmerna"
mkdir -p "$LOG_DIR"
mkdir -p "$SORTMERNA_DB_DIR"


# Function to send email notifications
send_email() {
    local status="$1"
    local message="$2"
    echo -e "Subject: Metatranscriptomics Pipeline - $status\n\n$message" | sendmail "$EMAIL"
}

# Function to log messages
echo_log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
}

# Error Handling
trap 'send_email "FAILED" "Metatranscriptomics Pipeline failed at: $(date)"; exit 1' ERR
set -e  # Exit on error

# Start logging
echo_log "Pipeline started at $START_TIME"

# Parse command-line arguments (allow multiple steps)
while getopts "s:" opt; do
    case ${opt} in
        s)
            IFS=',' read -r -a STEPS <<< "$OPTARG"  # Split steps by commas
            ;;
        *)
            echo_log "Usage: $0 -s <step1,step2,...>" | tee -a $LOG_FILE
            exit 1
            ;;
    esac
done

# Loop through each specified step
for STEP in "${STEPS[@]}"; do
    case "$STEP" in

    qc)
        echo_log "Step 2: Quality Control (FastQC)" | tee -a $LOG_FILE
        fastqc --noextract -t 8 -o $OUTPUT_DIR/quality_control $INPUT_DIR/*.fastq.gz &>> $LOG_FILE
        multiqc $OUTPUT_DIR/quality_control -o $OUTPUT_DIR/quality_control/multiqc_report &>> $LOG_FILE
        echo_log "Step 2: Quality Control (FastQC) completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 2: Quality Control (FastQC) completed successfully!"
        ;;

    trim_fastp)
        echo_log "Step 3: Trimming adpaters with fastp..." | tee -a $LOG_FILE
        for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do
            if [[ -f "$R1" ]]; then
                R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
                SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)

                FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"

                TRIMMED_R1="$FASTP_OUTPUT/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
                TRIMMED_R2="$FASTP_OUTPUT/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
                JSON_REPORT="$FASTP_OUTPUT/${SAMPLE_NAME}_fastp.json"
                HTML_REPORT="$FASTP_OUTPUT/${SAMPLE_NAME}_fastp.html"

                echo_log "Running fastp for $SAMPLE_NAME..." | tee -a $LOG_FILE
                fastp --in1 "$R1" --in2 "$R2" \
                      --out1 "$TRIMMED_R1" --out2 "$TRIMMED_R2" \
                      --detect_adapter_for_pe \
                      --length_required 30 \
                      --thread $THREADS \
                      --json "$JSON_REPORT" \
                      --html "$HTML_REPORT" &>> "$LOG_FILE"
                echo_log "fastp trimming completed for $SAMPLE_NAME!" | tee -a $LOG_FILE
            fi
        done
        echo_log "Step 3: Trimming with fastp completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 3: Trimming with fastp completed successfully!"
        ;;

    postqc_trim_fastp)
        echo_log "Step 4: Running FastQC and MultiQC after fastp trimming..." | tee -a $LOG_FILE

        FASTQC_FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastqc/fastp"

        fastqc --noextract -t $THREADS -o "$FASTQC_FASTP_OUTPUT" "$FASTP_OUTPUT"/*.fastq.gz &>> "$LOG_FILE"
        multiqc "$FASTQC_FASTP_OUTPUT" -o "$FASTQC_FASTP_OUTPUT/multiqc_report" &>> "$LOG_FILE"

        echo_log "Step 4: Post-fastp Quality Control completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 4: Post-fastp Quality Control completed successfully!"
        ;;

    download_db)
        echo_log "Step 7: Downloading SortMeRNA database..." | tee -a $LOG_FILE
        wget -O "$SORTMERNA_DB_DIR/database.tar.gz" "https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz"
        tar -xzf "$SORTMERNA_DB_DIR/database.tar.gz" -C "$SORTMERNA_DB_DIR"
        rm "$SORTMERNA_DB_DIR/database.tar.gz"
        echo_log "SortMeRNA database downloaded!" | tee -a $LOG_FILE
        ;;

    index_db)
        echo_log "Step 8: Indexing SortMeRNA database..." | tee -a $LOG_FILE
        indexdb_rna --ref "$SORTMERNA_DB_DIR/smr_v4.3_default_db.fasta","$SORTMERNA_DB_DIR/smr_v4.3_default_db.idx"
        echo_log "Indexing completed!" | tee -a $LOG_FILE
        ;;

    remove_rrna_unmerged)
        echo_log "Step 10: Removing rRNA for unmerged fast-p reads..." | tee -a $LOG_FILE


        FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"
        SORTMERNA_OUTPUT="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp"

        IDX_DIR="$SORTMERNA_OUTPUT/idx"
        KVDB_DIR="$SORTMERNA_OUTPUT/kvdb"
        READB_DIR="$SORTMERNA_OUTPUT/readb"
        OUT_ALIGNED="$SORTMERNA_OUTPUT/out/aligned"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"

        mkdir -p "$SORTMERNA_OUTPUT"
        mkdir -p "$IDX_DIR"
        mkdir -p "$KVDB_DIR"
        mkdir -p "$READB_DIR"
        mkdir -p "$OUT_ALIGNED"
        mkdir -p "$OUT_MRNA"

        # Loop through paired-end files
        for forward in "$FASTP_OUTPUT/"*_R1.fastq.gz; do
            if [[ -f "$forward" ]]; then
                SAMPLE_NAME=$(basename "$forward" _R1.fastq.gz)
                reverse="$FASTP_OUTPUT/${SAMPLE_NAME}_R2.fastq.gz"

                if [[ ! -f "$reverse" ]]; then
                    echo_log "Skipping $SAMPLE_NAME: Missing reverse read ($reverse)!" | tee -a "$LOG_FILE"
                    continue
                fi

                # Create per-sample directories
                SAMPLE_KVDB="$KVDB_DIR/$SAMPLE_NAME"
                SAMPLE_READB="$READB_DIR/$SAMPLE_NAME"
                mkdir -p "$SAMPLE_KVDB" "$SAMPLE_READB"


        echo_log "Processing $SAMPLE_NAME for rRNA removal..." | tee -a "$LOG_FILE"

        # Run SortMeRNA
        sortmerna \
            --ref "$SORTMERNA_DB_DIR/silva-euk-28s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-euk-18s-id95.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-bac-23s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-bac-16s-id90.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-arc-23s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-arc-16s-id95.fasta" \
            --ref "$SORTMERNA_DB_DIR/rfam-5s-database-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/rfam-5.8s-database-id98.fasta" \
            --reads "$forward" --reads "$reverse" \
            --out2 \
            --paired_out \
            --fastx \
            --blast "1 cigar qcov" \
            --aligned "$OUT_ALIGNED/${SAMPLE_NAME}_aligned" \
            --other "$OUT_MRNA/${SAMPLE_NAME}_mrna" \
            --kvdb "$SAMPLE_KVDB" \
            --readb "$SAMPLE_READB" \
            --v \
            --threads "$THREADS" --index 2 \
            --workdir "$SORTMERNA_OUTPUT" \
            --idx-dir "$IDX_DIR" &>> "$LOG_FILE"

        # Capture exit status
        SORTMERNA_EXIT=$?
        if [[ $SORTMERNA_EXIT -ne 0 ]]; then
            echo_log "Error: SortMeRNA failed for $SAMPLE_NAME with exit code $SORTMERNA_EXIT" | tee -a "$LOG_FILE"
            send_email "FAILED" "Step 10: SortMeRNA failed for $SAMPLE_NAME (Exit Code: $SORTMERNA_EXIT)"
            exit 1
        fi

            echo_log "rRNA removal completed for $SAMPLE_NAME!" | tee -a $LOG_FILE
            else
            echo_log "Skipping $SAMPLE_NAME: File does not exist!" | tee -a $LOG_FILE
            fi
        done
        ;;

    postqc_sortmerna)
        echo_log "Step: Running FastQC and MultiQC after removing rRNA..." | tee -a $LOG_FILE

        SORTMERNA_OUTPUT="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"
        SORTMERNA_QC="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/qc"

        fastqc --noextract -t $THREADS -o "$SORTMERNA_QC" "$OUT_MRNA"/*.fq.gz &>> "$LOG_FILE"
        multiqc "$SORTMERNA_QC" -o "$SORTMERNA_QC/multiqc_report" &>> "$LOG_FILE"

        echo_log "Step: Post-sortmerna Quality Control completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step: Post-sorting Quality Control completed successfully!"
        ;;

    summary_sortmerna)
        echo_log "Step: Summarizing sortmerna results..." | tee -a $LOG_FILE

        # Directories
        FASTP_OUTPUT="/home/rguden/Documents/ORCHESTRA/RNA_output/trimming/fastp"
        SORTMERNA_OUTPUT="/home/rguden/Documents/ORCHESTRA/RNA_output/sortmerna/sortmerna_unmerged_fastp"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"
        OUT_RRNA="$SORTMERNA_OUTPUT/out/aligned"

        # Sample metadata file (make sure this is tab-separated!)
        METADATA_FILE="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"

        # Output summary files
        SUMMARY_FILE="$SORTMERNA_OUTPUT/sortmerna_summary.tsv"
        SUMMARY_COMBINED="$SORTMERNA_OUTPUT/sortmerna_summary_combined.tsv"

        # Create TSV headers (now include condition)
        echo -e "Sample\tCondition\tTotal_Reads_Trimmed_R1\tmRNA_Reads_R1\trRNA_Reads_R1\tmRNA_Percentage_R1\trRNA_Percentage_R1\tTotal_Reads_Trimmed_R2\tmRNA_Reads_R2\trRNA_Reads_R2\tmRNA_Percentage_R2\trRNA_Percentage_R2" > "$SUMMARY_FILE"
        echo -e "Sample\tCondition\tTotal_Reads_Trimmed\tmRNA_Reads\trRNA_Reads\tmRNA_Percentage\trRNA_Percentage" > "$SUMMARY_COMBINED"

        # Read metadata into associative arrays (requires Bash >= 4)
        declare -A SAMPLE_TO_CONDITION

        while IFS=$'\t' read -r sample condition replicate; do
            [[ $sample == "sample" ]] && continue  # skip header
            SAMPLE_TO_CONDITION["$sample"]="$condition"
        done < "$METADATA_FILE"

        # Loop through forward reads (_R1)
        for forward in "$FASTP_OUTPUT/"*_R1.fastq.gz; do
            if [[ -f "$forward" ]]; then
                # Get the true file base (with _trimmed) for file lookup
                SAMPLE_FILEBASE=$(basename "$forward" _R1.fastq.gz)
                # Get the sample name (strip _trimmed) for metadata lookup
                SAMPLE_NAME=${SAMPLE_FILEBASE%_trimmed}

                # Now find the correct reverse filename (needs _trimmed)
                reverse="$FASTP_OUTPUT/${SAMPLE_FILEBASE}_R2.fastq.gz"

                # Metadata: use SAMPLE_NAME
                CONDITION="${SAMPLE_TO_CONDITION[$SAMPLE_NAME]}"
                [[ -z "$CONDITION" ]] && CONDITION="UNKNOWN"

                # Ensure reverse file exists
                if [[ ! -f "$reverse" ]]; then
                    echo "Skipping $SAMPLE_NAME: Missing reverse read!" | tee -a "$LOG_FILE"
                    continue
                fi

                # Read counts
                TOTAL_R1=$(zcat "$forward" | awk '{s++} END {print s/4}')
                TOTAL_R2=$(zcat "$reverse" | awk '{s++} END {print s/4}')

                # mRNA
                if [[ -f "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_fwd.fq.gz" ]]; then
                    MRNA_R1=$(zcat "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_fwd.fq.gz" | awk '{s++} END {print s/4}')
                else
                    MRNA_R1=0
                fi
                if [[ -f "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_rev.fq.gz" ]]; then
                    MRNA_R2=$(zcat "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_rev.fq.gz" | awk '{s++} END {print s/4}')
                else
                    MRNA_R2=0
                fi

                # rRNA
                if [[ -f "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_fwd.fq.gz" ]]; then
                    RRNA_R1=$(zcat "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_fwd.fq.gz" | awk '{s++} END {print s/4}')
                else
                    RRNA_R1=0
                fi
                if [[ -f "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_rev.fq.gz" ]]; then
                    RRNA_R2=$(zcat "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_rev.fq.gz" | awk '{s++} END {print s/4}')
                else
                    RRNA_R2=0
                fi

                # Percentages
                if [[ $TOTAL_R1 -gt 0 ]]; then
                    MRNA_PERC_R1=$(awk "BEGIN {print ($MRNA_R1/$TOTAL_R1)*100}")
                    RRNA_PERC_R1=$(awk "BEGIN {print ($RRNA_R1/$TOTAL_R1)*100}")
                else
                    MRNA_PERC_R1=0
                    RRNA_PERC_R1=0
                fi
                if [[ $TOTAL_R2 -gt 0 ]]; then
                    MRNA_PERC_R2=$(awk "BEGIN {print ($MRNA_R2/$TOTAL_R2)*100}")
                    RRNA_PERC_R2=$(awk "BEGIN {print ($RRNA_R2/$TOTAL_R2)*100}")
                else
                    MRNA_PERC_R2=0
                    RRNA_PERC_R2=0
                fi

                # Output detailed per-direction summary
                echo -e "$SAMPLE_NAME\t$CONDITION\t$TOTAL_R1\t$MRNA_R1\t$RRNA_R1\t$MRNA_PERC_R1\t$RRNA_PERC_R1\t$TOTAL_R2\t$MRNA_R2\t$RRNA_R2\t$MRNA_PERC_R2\t$RRNA_PERC_R2" >> "$SUMMARY_FILE"

                # Combined counts
                TOTAL_COMBINED=$((TOTAL_R1 + TOTAL_R2))
                MRNA_COMBINED=$((MRNA_R1 + MRNA_R2))
                RRNA_COMBINED=$((RRNA_R1 + RRNA_R2))

                if [[ $TOTAL_COMBINED -gt 0 ]]; then
                    MRNA_PERC_COMBINED=$(awk "BEGIN {print ($MRNA_COMBINED/$TOTAL_COMBINED)*100}")
                    RRNA_PERC_COMBINED=$(awk "BEGIN {print ($RRNA_COMBINED/$TOTAL_COMBINED)*100}")
                else
                    MRNA_PERC_COMBINED=0
                    RRNA_PERC_COMBINED=0
                fi

                echo -e "$SAMPLE_NAME\t$CONDITION\t$TOTAL_COMBINED\t$MRNA_COMBINED\t$RRNA_COMBINED\t$MRNA_PERC_COMBINED\t$RRNA_PERC_COMBINED" >> "$SUMMARY_COMBINED"
            fi
        done

        echo "SortMeRNA summary table generated: $SUMMARY_FILE"
        echo_log "Summarizing sortmerna results completed successfully" | tee -a $LOG_FILE
        ;;

    # conda install mamba -n base -c conda-forge
    # mamba create --name humann3_env python=3.9
    # conda activate humann3_env
    # mamba install humann -c biobakery -c defaults -c conda-forge -c bioconda # Keep all channels
    # metaphlan install: to download the metaphlan databases
    download_humann3_db)
        echo_log "Step: Downloading HUMAnN3 databases..." | tee -a $LOG_FILE

        DB_DIR="/home/rguden/Documents/databases/humann3"


        echo_log "[INFO] Downloading ChocoPhlAn database..." | tee -a $LOG_FILE
        humann_databases --download chocophlan full "$DB_DIR" --update-config yes &>> "$LOG_FILE"

        echo_log "[INFO] Downloading UniRef90 database (diamond)..." | tee -a $LOG_FILE
        humann_databases --download uniref uniref90_diamond "$DB_DIR" --update-config yes &>> "$LOG_FILE"

        echo_log "[INFO] HUMAnN3 database downloads completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step: HUMAnN3 database downloads completed."
        ;;


   merge_mrna_reads)
        echo_log "Step: Merging paired-end mRNA reads for HUMAnN3..." | tee -a $LOG_FILE

        MRNA_DIR="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/out/mrna"
        MERGED_DIR="$MRNA_DIR/merged"
        mkdir -p "$MERGED_DIR"

        for R1 in "$MRNA_DIR"/*_fwd.fq.gz; do
            if [[ -f "$R1" ]]; then
                SAMPLE=$(basename "$R1" _fwd.fq.gz)
                R2="$MRNA_DIR/${SAMPLE}_rev.fq.gz"
                if [[ -f "$R2" ]]; then
                    echo_log "Merging $SAMPLE R1 and R2 reads..." | tee -a $LOG_FILE
                    cat "$R1" "$R2" > "$MERGED_DIR/${SAMPLE}_merged.fq.gz"
                else
                    echo_log "[WARNING] Reverse read for $SAMPLE not found, skipping." | tee -a $LOG_FILE
                fi
            fi
        done

        echo_log "Step: Merging mRNA paired-end reads completed." | tee -a $LOG_FILE
        send_email "COMPLETED" "Step: Merging mRNA paired-end reads completed."
        ;;

    # run metaphlan separately to avoid database compatability issues (metaphlan v 3.0; database v31 installed using metaphlan --install)
    run_metaphlan)
        echo_log "Step: Running MetaPhlAn on rRNA-depleted fastq files..." | tee -a "$LOG_FILE"

        DB_BOWTIE="/home/rguden/miniconda3/envs/humann3_env/lib/python3.9/site-packages/metaphlan/metaphlan_databases"
        INDEX_NAME="mpa_v31_CHOCOPhlAn_201901"
        MERGED_DIR="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/out/mrna/merged"
        METAPHLAN_OUT_DIR="$OUTPUT_DIR/metaphlan"

        # Verify that the index exists
        if [[ ! -f "$DB_BOWTIE/${INDEX_NAME}.pkl" ]]; then
            echo_log "[ERROR] MetaPhlAn index file $INDEX_NAME.pkl not found in $DB_BOWTIE" | tee -a "$LOG_FILE"
            send_email "FAILED" "MetaPhlAn index file not found: $INDEX_NAME"
            exit 1
        fi

        mkdir -p "$METAPHLAN_OUT_DIR"

        for fq in "$MERGED_DIR"/*.fq.gz; do
            if [[ -f "$fq" ]]; then
                SAMPLE_NAME=$(basename "$fq" _merged.fq.gz)
                echo_log "Running MetaPhlAn for $SAMPLE_NAME..." | tee -a "$LOG_FILE"

                metaphlan "$fq" \
                    --input_type fastq \
                    --nproc "$THREADS" \
                    --bowtie2db "$DB_BOWTIE" \
                    --index "$INDEX_NAME" \
                    --bowtie2out "$METAPHLAN_OUT_DIR/${SAMPLE_NAME}_metaphlan.bowtie2.bz2" \
                    -o "$METAPHLAN_OUT_DIR/${SAMPLE_NAME}_profile.tsv"

                if [[ $? -ne 0 ]]; then
                    echo_log "[ERROR] MetaPhlAn failed for $SAMPLE_NAME" | tee -a "$LOG_FILE"
                    send_email "FAILED" "MetaPhlAn failed for $SAMPLE_NAME"
                    exit 1
                else
                    echo_log "MetaPhlAn completed for $SAMPLE_NAME!" | tee -a "$LOG_FILE"
                fi
            fi
        done

        echo_log "Step: MetaPhlAn profiling completed successfully!" | tee -a "$LOG_FILE"
        send_email "COMPLETED" "Step: MetaPhlAn profiling completed."
        ;;

    run_humann3)
        echo_log "Step: Running HUMAnN3 using external MetaPhlAn profiles..." | tee -a "$LOG_FILE"

        DB_DIR="/home/rguden/Documents/databases/humann3"
        MERGED_DIR="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/out/mrna/merged"
        METAPHLAN_OUT_DIR="$OUTPUT_DIR/metaphlan"
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"

        mkdir -p "$HUMANN3_OUT_DIR/raw"

        for fq in "$MERGED_DIR"/*.fq.gz; do
            if [[ -f "$fq" ]]; then
                SAMPLE_NAME=$(basename "$fq" _merged.fq.gz)
                TAX_PROFILE="$METAPHLAN_OUT_DIR/${SAMPLE_NAME}_profile.tsv"

                echo_log "Running HUMAnN3 for $SAMPLE_NAME using MetaPhlAn output..." | tee -a "$LOG_FILE"

                humann \
                    --input "$fq" \
                    --output "$HUMANN3_OUT_DIR/raw" \
                    --threads "$THREADS" \
                    --nucleotide-database "$DB_DIR/chocophlan" \
                    --protein-database "$DB_DIR/uniref" \
                    --remove-temp-output \
                    --taxonomic-profile "$TAX_PROFILE"

                if [[ $? -ne 0 ]]; then
                    echo_log "[ERROR] HUMAnN3 failed for $SAMPLE_NAME" | tee -a "$LOG_FILE"
                    send_email "FAILED" "HUMAnN3 failed for $SAMPLE_NAME"
                    exit 1
                fi
                echo_log "HUMAnN3 completed for $SAMPLE_NAME!" | tee -a "$LOG_FILE"
            fi
        done

        echo_log "Step: HUMAnN3 functional profiling completed successfully!" | tee -a "$LOG_FILE"
        send_email "COMPLETED" "HUMAnN3 functional profiling completed."
        ;;

    join_humann3_tables)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        mkdir -p "$HUMANN3_OUT_DIR/join"
        echo_log "Step: Joining HUMAnN3 tables across samples..." | tee -a $LOG_FILE
        for table_type in genefamilies pathabundance pathcoverage; do
            humann_join_tables \
                --input "$HUMANN3_OUT_DIR/raw" \
                --file_name "${table_type}.tsv" \
                --output "$HUMANN3_OUT_DIR/join/joined_${table_type}.tsv"
        done
        echo_log "Step: HUMAnN3 table joining completed successfully!" | tee -a $LOG_FILE
        ;;

    normalize_humann3)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        mkdir -p "$HUMANN3_OUT_DIR/normalized"
        echo_log "Step: Normalizing HUMAnN3 tables (cpm and relab)..." | tee -a $LOG_FILE
        for method in cpm relab; do
            for table_type in genefamilies pathabundance; do
                humann_renorm_table \
                    --input "$HUMANN3_OUT_DIR/join/joined_${table_type}.tsv" \
                    --output "$HUMANN3_OUT_DIR/normalized/joined_${table_type}_${method}.tsv" \
                    --units $method \
                    --update-snames
            done
        done
        echo_log "Step: HUMAnN3 normalization completed successfully!" | tee -a $LOG_FILE
        ;;

    split_humann3)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        mkdir -p "$HUMANN3_OUT_DIR/split"
        echo_log "Step: Splitting stratified and unstratified HUMAnN3 tables..." | tee -a $LOG_FILE
        for method in cpm relab; do
            for table_type in genefamilies pathabundance; do
                humann_split_stratified_table \
                    --input "$HUMANN3_OUT_DIR/normalized/joined_${table_type}_${method}.tsv" \
                    --output "$HUMANN3_OUT_DIR/split/${table_type}_${method}"
            done
        done
        echo_log "Step: HUMAnN3 split completed successfully!" | tee -a $LOG_FILE
        ;;

    download_utility_mapping)
        DB_DIR="/home/rguden/Documents/databases/humann3"
        echo_log "Step: Downloading HUMAnN3 utility mapping files..." | tee -a $LOG_FILE
        humann_databases --download utility_mapping full "$DB_DIR" --update-config yes
        echo_log "Step: Utility mapping download completed!" | tee -a $LOG_FILE
        ;;

    regroup_genefamilies_raw)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        REGROUP_DIR="$HUMANN3_OUT_DIR/regroup"

        mkdir -p "$REGROUP_DIR"

        # Define the path to your utility mapping database
        UTILITY_DB_DIR="/home/rguden/Documents/databases/humann3/utility_mapping"

        echo_log "Step: Regrouping HUMAnN3 gene families frow raw count tables..." | tee -a "$LOG_FILE"

        # Value: "group_map_file|name_map_file" (separated by '|')
        declare -A MAP_CONFIGS
        MAP_CONFIGS=(
            [ko]="map_ko_uniref90.txt.gz|map_ko_name.txt.gz"
            [eggnog]="map_eggnog_uniref90.txt.gz|map_eggnog_name.txt.gz"
            [pfam]="map_pfam_uniref90.txt.gz|map_pfam_name.txt.gz"
            [ec]="map_level4ec_uniref90.txt.gz|map_ec_name.txt.gz"
            [go]="map_go_uniref90.txt.gz|map_go_name.txt.gz"
        )

            for map_label in "${!MAP_CONFIGS[@]}"; do
                # Split the combined string into group_map_file and name_map_file
                IFS='|' read -r GROUP_MAP_FILENAME NAME_MAP_FILENAME <<< "${MAP_CONFIGS[$map_label]}"

                GROUP_MAP_FILE="$UTILITY_DB_DIR/$GROUP_MAP_FILENAME"
                NAME_MAP_FILE="$UTILITY_DB_DIR/$NAME_MAP_FILENAME" # This will be the name map for --function-description

                # Determine the correct --groups parameter type based on map_label
                # For example, for 'ko', --groups is 'KO'. For 'ec', --groups is 'EC'.

                case "$map_label" in
                    ko) GROUP_TYPE="uniref90_ko" ;; # Changed from "KO" for consistency and correctness
                    eggnog) GROUP_TYPE="uniref90_eggnog" ;; # THIS IS THE MAIN FIX
                    pfam) GROUP_TYPE="uniref90_pfam" ;; # Changed from "Pfam"
                    ec) GROUP_TYPE="uniref90_level4ec" ;; # Changed from "EC"
                    go) GROUP_TYPE="uniref90_go" ;; # Changed from "GO"
                    *) echo_log "[ERROR] Unknown map_label: $map_label" | tee -a "$LOG_FILE"; continue ;;
                esac

                # Regroup non-normalized data (raw counts)
                RAW_FILE="$HUMANN3_OUT_DIR/join/joined_genefamilies.tsv"
                OUT_RAW="$REGROUP_DIR/non_normalized"
                mkdir -p "$OUT_RAW"

                OUT_RAW="$REGROUP_DIR/non_normalized/genefamilies_${map_label}_joined.tsv"
                echo_log "Regrouping raw table to $map_label ($GROUP_TYPE)..." | tee -a "$LOG_FILE"

                humann_regroup_table \
                --input "$RAW_FILE" \
                --groups "$GROUP_TYPE" \
                --output "$OUT_RAW"

            done

        echo_log "Step: Regrouping for all table types and mapping systems completed." | tee -a "$LOG_FILE"
        send_email "COMPLETED" "Step: HUMAnN3 gene family regrouping (all tables) completed."
        ;;

    regroup_genefamilies)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        REGROUP_DIR="$HUMANN3_OUT_DIR/regroup"

        mkdir -p "$REGROUP_DIR"

        # Define the path to your utility mapping database
        UTILITY_DB_DIR="/home/rguden/Documents/databases/humann3/utility_mapping"

        echo_log "Step: Regrouping HUMAnN3 gene families (EC, KO, GO, Pfam, EggNOG) from both CPM and RelAb tables (joined, stratified, unstratified)..." | tee -a "$LOG_FILE"

        # Value: "group_map_file|name_map_file" (separated by '|')
        declare -A MAP_CONFIGS
        MAP_CONFIGS=(
            [ko]="map_ko_uniref90.txt.gz|map_ko_name.txt.gz"
            [eggnog]="map_eggnog_uniref90.txt.gz|map_eggnog_name.txt.gz"
            [pfam]="map_pfam_uniref90.txt.gz|map_pfam_name.txt.gz"
            [ec]="map_level4ec_uniref90.txt.gz|map_ec_name.txt.gz"
            [go]="map_go_uniref90.txt.gz|map_go_name.txt.gz"
        )

        for norm in relab cpm; do
            # Define output directory for each normalization type
            mkdir -p "$REGROUP_DIR/$norm"

            for map_label in "${!MAP_CONFIGS[@]}"; do
                # Split the combined string into group_map_file and name_map_file
                IFS='|' read -r GROUP_MAP_FILENAME NAME_MAP_FILENAME <<< "${MAP_CONFIGS[$map_label]}"

                GROUP_MAP_FILE="$UTILITY_DB_DIR/$GROUP_MAP_FILENAME"
                NAME_MAP_FILE="$UTILITY_DB_DIR/$NAME_MAP_FILENAME" # This will be the name map for --function-description

                # Determine the correct --groups parameter type based on map_label
                # For example, for 'ko', --groups is 'KO'. For 'ec', --groups is 'EC'.

                case "$map_label" in
                    ko) GROUP_TYPE="uniref90_ko" ;; # Changed from "KO" for consistency and correctness
                    eggnog) GROUP_TYPE="uniref90_eggnog" ;; # THIS IS THE MAIN FIX
                    pfam) GROUP_TYPE="uniref90_pfam" ;; # Changed from "Pfam"
                    ec) GROUP_TYPE="uniref90_level4ec" ;; # Changed from "EC"
                    go) GROUP_TYPE="uniref90_go" ;; # Changed from "GO"
                    *) echo_log "[ERROR] Unknown map_label: $map_label" | tee -a "$LOG_FILE"; continue ;;
                esac

                # Regroup non-normalized data (raw counts)
                RAW_FILE="$HUMANN3_OUT_DIR/joined/joined_genefamilies_${norm}.tsv"
                OUT_RAW="$REGROUP_DIR/non_normalized"
                mkdir -p "$OUT_RAW"

                OUT_RAW="$REGROUP_DIR/non_normalized/genefamilies_${map_label}_joined.tsv"
                echo_log "Regrouping joined $norm to $map_label ($GROUP_TYPE)..." | tee -a "$LOG_FILE"

                humann_regroup_table \
                --input "$RAW_FILE" \
                --groups "$GROUP_TYPE" \
                --output "$OUT_RAW"

                # Regroup joined table
                JOINED_FILE="$HUMANN3_OUT_DIR/normalized/joined_genefamilies_${norm}.tsv"
                OUT_JOIN="$REGROUP_DIR/$norm/joined"
                mkdir -p "$OUT_JOIN"

                OUT_FILE="$REGROUP_DIR/$norm/joined/genefamilies_${map_label}_joined.tsv"
                echo_log "Regrouping joined $norm to $map_label ($GROUP_TYPE)..." | tee -a "$LOG_FILE"

                humann_regroup_table \
                --input "$JOINED_FILE" \
                --groups "$GROUP_TYPE" \
                --output "$OUT_FILE"

                # Regroup stratified table
                STRAT_FILE="$HUMANN3_OUT_DIR/split/genefamilies_${norm}/joined_genefamilies_${norm}_stratified.tsv"
                OUT_STRAT="$REGROUP_DIR/$norm/strat"
                mkdir -p "$OUT_STRAT"

                OUT_FILE="$REGROUP_DIR/$norm/strat/joined_genefamilies_${map_label}_stratified.tsv"
                echo_log "Regrouping stratified $norm to $map_label ($GROUP_TYPE)..." | tee -a "$LOG_FILE"

                humann_regroup_table \
                --input "$STRAT_FILE" \
                --groups "$GROUP_TYPE" \
                --output "$OUT_FILE"

                # Regroup unstratified table
                UNSTRAT_FILE="$HUMANN3_OUT_DIR/split/genefamilies_${norm}/joined_genefamilies_${norm}_stratified.tsv"
                OUT_UNSTRAT="$REGROUP_DIR/$norm/unstrat"
                mkdir -p "$OUT_UNSTRAT"

                OUT_FILE="$REGROUP_DIR/$norm/unstrat/genefamilies_${map_label}_unstratified.tsv"
                echo_log "Regrouping unstratified $norm to $map_label ($GROUP_TYPE)..." | tee -a "$LOG_FILE"

                humann_regroup_table \
                    --input "$UNSTRAT_FILE" \
                    --groups "$GROUP_TYPE" \
                    --output "$OUT_FILE"
            done
        done

        echo_log "Step: Regrouping for all table types and mapping systems completed." | tee -a "$LOG_FILE"
        send_email "COMPLETED" "Step: HUMAnN3 gene family regrouping (all tables) completed."
        ;;


    rename_features_raw)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        REGROUP_DIR="$HUMANN3_OUT_DIR/regroup"
        RENAME_DIR="$HUMANN3_OUT_DIR/renamed"
        RENAME_DIR_RAW="$RENAME_DIR/non_normalized"

        mkdir -p "$RENAME_DIR_RAW"

        echo_log "Step: Renaming HUMAnN3 feature IDs for all regrouped and normalized outputs..." | tee -a "$LOG_FILE"

        # Mapping label to rename type expected by humann_rename_table
        declare -A RENAME_TYPES_MAP=(
            [ko]="kegg-orthology"
            [eggnog]="eggnog"
            [pfam]="pfam"
            [ec]="ec"
            [go]="go"
        )

            ######################
            # Gene Families
            ######################

            for map_label in "${!RENAME_TYPES_MAP[@]}"; do
                RENAME_TYPE="${RENAME_TYPES_MAP[$map_label]}"

                # JOINED
                INPUT_FILE_JOINED="$REGROUP_DIR/non_normalized/genefamilies_${map_label}_joined.tsv"
                OUTPUT_FILE_JOINED="$RENAME_DIR/non_normalized/genefamilies_${map_label}_joined_renamed.tsv"
                if [[ -f "$INPUT_FILE_JOINED" ]]; then
                    echo_log "Renaming $map_label features in $norm-normalized joined genefamilies..." | tee -a "$LOG_FILE"
                    humann_rename_table \
                        --input "$INPUT_FILE_JOINED" \
                        --output "$OUTPUT_FILE_JOINED" \
                        --names "$RENAME_TYPE"
                else
                    echo_log "[WARNING] Missing file: $INPUT_FILE_JOINED" | tee -a "$LOG_FILE"
                fi
            done


            ######################
            # Path Abundance
            ######################

            # JOINED
            INPUT_PATH_JOINED="$HUMANN3_OUT_DIR/join/joined_pathabundance.tsv"
            OUTPUT_PATH_JOINED="$RENAME_DIR/non_normalized/pathabundance_joined_renamed.tsv"
            if [[ -f "$INPUT_PATH_JOINED" ]]; then
                echo_log "Renaming pathway features in non-normalized joined pathabundance..." | tee -a "$LOG_FILE"
                humann_rename_table \
                    --input "$INPUT_PATH_JOINED" \
                    --output "$OUTPUT_PATH_JOINED" \
                    --names "metacyc-pwy"
            else
                echo_log "[WARNING] Missing file: $INPUT_PATH_JOINED" | tee -a "$LOG_FILE"
            fi

        echo_log "Step: HUMAnN3 renaming completed for all gene families and pathway abundance tables." | tee -a "$LOG_FILE"

        ;;
        
   split_renamed_raw)
      echo_log "Step: Splitting stratified and unstratified HUMAnN3 tables (non-normalized, renamed)..." | tee -a "$LOG_FILE"
      HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
      REGROUP_DIR="$HUMANN3_OUT_DIR/regroup"
      RENAME_DIR="$HUMANN3_OUT_DIR/renamed"
      RENAME_DIR_RAW="$RENAME_DIR/non_normalized"
      
      # Mapping label to rename type expected by humann_rename_table
        declare -A RENAME_TYPES_MAP=(
            [ko]="kegg-orthology"
            [eggnog]="eggnog"
            [pfam]="pfam"
            [ec]="ec"
            [go]="go"
        )
  
      for map_label in "${!RENAME_TYPES_MAP[@]}"; do
          INPUT_JOINED="$RENAME_DIR_RAW/genefamilies_${map_label}_joined_renamed.tsv"
          OUTPUT_SPLIT_DIR="$RENAME_DIR_RAW/split_${map_label}"
  
          if [[ -f "$INPUT_JOINED" ]]; then
              mkdir -p "$OUTPUT_SPLIT_DIR"
              humann_split_stratified_table \
                  --input "$INPUT_JOINED" \
                  --output "$OUTPUT_SPLIT_DIR"
              echo_log "Split completed for $map_label. Output in $OUTPUT_SPLIT_DIR." | tee -a "$LOG_FILE"
          else
              echo_log "[WARNING] Missing file for splitting: $INPUT_JOINED" | tee -a "$LOG_FILE"
          fi
      done
  
      send_email "COMPLETED" "Step: HUMAnN3 splitting of raw renamed tables completed."
      ;;


    rename_features)
        HUMANN3_OUT_DIR="$OUTPUT_DIR/humann3"
        REGROUP_DIR="$HUMANN3_OUT_DIR/regroup"
        RENAME_DIR="$HUMANN3_OUT_DIR/renamed"

        mkdir -p "$RENAME_DIR"

        echo_log "Step: Renaming HUMAnN3 feature IDs for all regrouped and normalized outputs..." | tee -a "$LOG_FILE"

        # Mapping label to rename type expected by humann_rename_table
        declare -A RENAME_TYPES_MAP=(
            [ko]="kegg-orthology"
            [eggnog]="eggnog"
            [pfam]="pfam"
            [ec]="ec"
            [go]="go"
        )

        for norm in relab cpm; do
            mkdir -p "$RENAME_DIR/$norm"

            ######################
            # Gene Families
            ######################

            for map_label in "${!RENAME_TYPES_MAP[@]}"; do
                RENAME_TYPE="${RENAME_TYPES_MAP[$map_label]}"

                # JOINED
                INPUT_FILE_JOINED="$REGROUP_DIR/$norm/joined/genefamilies_${map_label}_joined.tsv"
                OUTPUT_FILE_JOINED="$RENAME_DIR/$norm/genefamilies_${map_label}_joined_renamed.tsv"
                if [[ -f "$INPUT_FILE_JOINED" ]]; then
                    echo_log "Renaming $map_label features in $norm-normalized joined genefamilies..." | tee -a "$LOG_FILE"
                    humann_rename_table \
                        --input "$INPUT_FILE_JOINED" \
                        --output "$OUTPUT_FILE_JOINED" \
                        --names "$RENAME_TYPE"
                else
                    echo_log "[WARNING] Missing file: $INPUT_FILE_JOINED" | tee -a "$LOG_FILE"
                fi

                # STRATIFIED
                INPUT_FILE_STRAT="$REGROUP_DIR/$norm/strat/joined_genefamilies_${map_label}_stratified.tsv"
                OUTPUT_FILE_STRAT="$RENAME_DIR/$norm/genefamilies_${map_label}_stratified_renamed.tsv"
                if [[ -f "$INPUT_FILE_STRAT" ]]; then
                    echo_log "Renaming $map_label features in $norm-normalized stratified genefamilies..." | tee -a "$LOG_FILE"
                    humann_rename_table \
                        --input "$INPUT_FILE_STRAT" \
                        --output "$OUTPUT_FILE_STRAT" \
                        --names "$RENAME_TYPE"
                else
                    echo_log "[WARNING] Missing file: $INPUT_FILE_STRAT" | tee -a "$LOG_FILE"
                fi

                # UNSTRATIFIED
                INPUT_FILE_UNSTRAT="$REGROUP_DIR/$norm/unstrat/genefamilies_${map_label}_unstratified.tsv"
                OUTPUT_FILE_UNSTRAT="$RENAME_DIR/$norm/genefamilies_${map_label}_unstratified_renamed.tsv"
                if [[ -f "$INPUT_FILE_UNSTRAT" ]]; then
                    echo_log "Renaming $map_label features in $norm-normalized unstratified genefamilies..." | tee -a "$LOG_FILE"
                    humann_rename_table \
                        --input "$INPUT_FILE_UNSTRAT" \
                        --output "$OUTPUT_FILE_UNSTRAT" \
                        --names "$RENAME_TYPE"
                else
                    echo_log "[WARNING] Missing file: $INPUT_FILE_UNSTRAT" | tee -a "$LOG_FILE"
                fi
            done

            ######################
            # Path Abundance
            ######################

            # JOINED
            INPUT_PATH_JOINED="$HUMANN3_OUT_DIR/normalized/joined_pathabundance_${norm}.tsv"
            OUTPUT_PATH_JOINED="$RENAME_DIR/$norm/pathabundance_joined_renamed.tsv"
            if [[ -f "$INPUT_PATH_JOINED" ]]; then
                echo_log "Renaming pathway features in $norm-normalized joined pathabundance..." | tee -a "$LOG_FILE"
                humann_rename_table \
                    --input "$INPUT_PATH_JOINED" \
                    --output "$OUTPUT_PATH_JOINED" \
                    --names "metacyc-pwy"
            else
                echo_log "[WARNING] Missing file: $INPUT_PATH_JOINED" | tee -a "$LOG_FILE"
            fi

            # STRATIFIED
            INPUT_PATH_STRAT="$HUMANN3_OUT_DIR/split/pathabundance_${norm}/joined_pathabundance_${norm}_stratified.tsv"
            OUTPUT_PATH_STRAT="$RENAME_DIR/$norm/pathabundance_stratified_renamed.tsv"
            if [[ -f "$INPUT_PATH_STRAT" ]]; then
                echo_log "Renaming pathway features in $norm-normalized stratified pathabundance..." | tee -a "$LOG_FILE"
                humann_rename_table \
                    --input "$INPUT_PATH_STRAT" \
                    --output "$OUTPUT_PATH_STRAT" \
                    --names "metacyc-pwy"
            else
                echo_log "[WARNING] Missing file: $INPUT_PATH_STRAT" | tee -a "$LOG_FILE"
            fi

            # UNSTRATIFIED
            INPUT_PATH_UNSTRAT="$HUMANN3_OUT_DIR/split/pathabundance_${norm}/joined_pathabundance_${norm}_unstratified.tsv"
            OUTPUT_PATH_UNSTRAT="$RENAME_DIR/$norm/pathabundance_unstratified_renamed.tsv"
            if [[ -f "$INPUT_PATH_UNSTRAT" ]]; then
                echo_log "Renaming pathway features in $norm-normalized unstratified pathabundance..." | tee -a "$LOG_FILE"
                humann_rename_table \
                    --input "$INPUT_PATH_UNSTRAT" \
                    --output "$OUTPUT_PATH_UNSTRAT" \
                    --names "metacyc-pwy"
            else
                echo_log "[WARNING] Missing file: $INPUT_PATH_UNSTRAT" | tee -a "$LOG_FILE"
            fi

        done

        echo_log "Step: HUMAnN3 renaming completed for all gene families and pathway abundance tables." | tee -a "$LOG_FILE"
        send_email "COMPLETED" "Step: HUMAnN3 renaming for gene families and pathways completed."

        ;;

    analyze_deseq2_kegg_enrichment)

        echo_log "Step: Running DESeq2 analysis and KEGG enrichment for KO features..." | tee -a "$LOG_FILE"

        # Define input/output
        RENAME_DIR="$OUTPUT_DIR/humann3/renamed"
        COMPARISON_OUT="$OUTPUT_DIR/humann3/functional_kegg_deseq2"
        METADATA_FILE="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"
        R_SCRIPT="/home/rguden/Documents/scripts/deseq2_kegg_pathway.R"

        # Create output directory
        mkdir -p "$COMPARISON_OUT"

        # Use KO gene family joined + renamed table for analysis (cpm normalization)
        INPUT_TSV="$RENAME_DIR/non_normalized/genefamilies_ko_joined_renamed.tsv"

        if [[ ! -f "$INPUT_TSV" ]]; then
            echo_log "[ERROR] KO input table not found at: $INPUT_TSV" | tee -a "$LOG_FILE"
            send_email "FAILED" "KO DESeq2 input missing. Check $INPUT_TSV"
            exit 1
        fi

        if [[ ! -f "$R_SCRIPT" ]]; then
            echo_log "[ERROR] R script not found: $R_SCRIPT" | tee -a "$LOG_FILE"
            send_email "FAILED" "R script for KO DESeq2 not found."
            exit 1
        fi

        # Run R script
        Rscript "$R_SCRIPT" \
            "$INPUT_TSV" \
            "$METADATA_FILE" \
            "$COMPARISON_OUT"

        if [[ $? -eq 0 ]]; then
            echo_log "DESeq2 and KEGG enrichment analysis completed." | tee -a "$LOG_FILE"
            send_email "COMPLETED" "Functional analysis (DESeq2 + KEGG) finished successfully."
        else
            echo_log "[ERROR] R analysis failed. Check logs and input files." | tee -a "$LOG_FILE"
            send_email "FAILED" "R analysis for DESeq2 + KEGG enrichment failed."
            exit 1
        fi
        ;;

    *)
        echo_log "Invalid step: $STEP" | tee -a "$LOG_FILE"
        exit 1
        ;;

    esac
done

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo_log "Pipeline completed at $END_TIME! Logs stored at $LOG_FILE" | tee -a $LOG_FILE
send_email "COMPLETED" "Metatranscriptomics Pipeline finished successfully at $END_TIME"

