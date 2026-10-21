process EXTRACT_RAW_DATA {
    tag "${library_id}"
    label 'extract'

    input:
    tuple val(library_id), val(tar_files), val(fastq_dir)

    output:
    tuple val(library_id), val(fastq_dir), emit: fastqs
    path("${library_id}_extract.log"), emit: log

    script:
    """
    set -o pipefail
    mkdir -p ${fastq_dir}

    {
    echo "============================================"
    echo "  EXTRACT_RAW_DATA: ${library_id}"
    echo "  \$(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================"

    IFS=';' read -ra TARS <<< "${tar_files}"
    TOTAL=\${#TARS[@]}
    CURRENT=0
    TOTAL_EXTRACTED=0
    TOTAL_PASSED=0
    TOTAL_FAILED=0
    ALL_OK=true

    echo "[INFO] ${library_id} has \${TOTAL} tar file(s) to process"
    echo ""

    for TAR in "\${TARS[@]}"; do
        TAR=\$(echo "\$TAR" | xargs)
        [ -z "\$TAR" ] && continue
        CURRENT=\$((CURRENT + 1))

        echo "--------------------------------------------"
        echo "[TAR \${CURRENT}/\${TOTAL}] \$TAR"

        if [ ! -f "\$TAR" ]; then
            echo "[ERROR] File not found: \$TAR"
            ALL_OK=false
            continue
        fi

        TMPDIR_EXTRACT=\$(mktemp -d)
        tar -xf \$TAR -C \${TMPDIR_EXTRACT}

        SAMPLE_DIR=\$(find \${TMPDIR_EXTRACT} -type d -name "${library_id}" | head -1)

        if [ -z "\$SAMPLE_DIR" ]; then
            echo "[WARNING] No folder named '${library_id}' found in this tar"
            rm -rf \${TMPDIR_EXTRACT}
            continue
        fi

        MD5_FILE=""
        if [ -f "\${SAMPLE_DIR}/MD5.txt" ]; then
            MD5_FILE="\${SAMPLE_DIR}/MD5.txt"
        else
            SEARCH_DIR=\$(dirname \$SAMPLE_DIR)
            while [ "\$SEARCH_DIR" != "\${TMPDIR_EXTRACT}" ] && [ "\$SEARCH_DIR" != "/" ]; do
                if [ -f "\${SEARCH_DIR}/MD5.txt" ]; then
                    MD5_FILE="\${SEARCH_DIR}/MD5.txt"
                    break
                fi
                SEARCH_DIR=\$(dirname \$SEARCH_DIR)
            done
            if [ -z "\$MD5_FILE" ] && [ -f "\${TMPDIR_EXTRACT}/MD5.txt" ]; then
                MD5_FILE="\${TMPDIR_EXTRACT}/MD5.txt"
            fi
        fi

        if [ -n "\$MD5_FILE" ]; then
            echo "[INFO] Found MD5 file: \$MD5_FILE"
        else
            echo "[WARNING] No MD5.txt found — skipping checksum verification"
        fi

        FASTQS_IN_TAR=\$(find \$SAMPLE_DIR -type f -name "*.fastq.gz")
        NEW=0

        for FQ in \$FASTQS_IN_TAR; do
            BASENAME=\$(basename \$FQ)
            mv \$FQ ${fastq_dir}/
            NEW=\$((NEW + 1))

            if [ -n "\$MD5_FILE" ]; then
                EXPECTED_MD5=\$(grep "\$BASENAME" "\$MD5_FILE" | awk '{print \$1}' | head -1)
                if [ -n "\$EXPECTED_MD5" ]; then
                    ACTUAL_MD5=\$(md5sum ${fastq_dir}/\$BASENAME | awk '{print \$1}')
                    if [ "\$EXPECTED_MD5" = "\$ACTUAL_MD5" ]; then
                        echo "  [OK]   \$BASENAME (md5: \$ACTUAL_MD5)"
                        TOTAL_PASSED=\$((TOTAL_PASSED + 1))
                    else
                        echo "  [FAIL] \$BASENAME"
                        echo "         Expected: \$EXPECTED_MD5"
                        echo "         Got:      \$ACTUAL_MD5"
                        TOTAL_FAILED=\$((TOTAL_FAILED + 1))
                        ALL_OK=false
                    fi
                else
                    echo "  [WARN] \$BASENAME — no MD5 entry found in MD5.txt"
                fi
            else
                echo "  [NEW]  \$BASENAME (no MD5 check)"
            fi
        done

        echo "[TAR \${CURRENT}/\${TOTAL}] Done: \$NEW files extracted"
        TOTAL_EXTRACTED=\$((TOTAL_EXTRACTED + NEW))

        rm -rf \${TMPDIR_EXTRACT}
    done

    echo ""
    echo "============================================"
    echo "  SUMMARY: ${library_id}"
    echo "============================================"
    FINAL_COUNT=\$(find ${fastq_dir} -name "*.fastq.gz" | wc -l)
    echo "  Tar files processed  : \${TOTAL}"
    echo "  Files extracted      : \${TOTAL_EXTRACTED}"
    echo "  MD5 passed           : \${TOTAL_PASSED}"
    echo "  MD5 failed           : \${TOTAL_FAILED}"
    echo "  Total fastqs in dir  : \${FINAL_COUNT}"
    echo "  Output directory     : ${fastq_dir}"
    echo ""
    echo "Final contents:"
    ls -lh ${fastq_dir}/*.fastq.gz 2>/dev/null || echo "  (no fastq.gz files)"

    if [ "\$ALL_OK" = true ] && [ "\$TOTAL_EXTRACTED" -gt 0 ]; then
        echo ""
        echo "  *** ALL MD5 CHECKS PASSED ***"
        echo "============================================"
    else
        echo ""
        [ "\$TOTAL_FAILED" -gt 0 ] && echo "  *** MD5 VERIFICATION FAILED ***"
        [ "\$TOTAL_EXTRACTED" -eq 0 ] && echo "  *** NO FASTQ FILES EXTRACTED ***"
        echo "============================================"
        exit 1
    fi

    } 2>&1 | tee ${library_id}_extract.log
    """
}

process MERGE_EXTRACT_LOGS {
    label 'basic'

    input:
    path(logs)
    val(logs_dir)

    output:
    path("01_Extract_raw_data.log")

    publishDir "${logs_dir}", mode: 'copy', overwrite: true

    script:
    """
    echo "================================================================" > 01_Extract_raw_data.log
    echo "  EXTRACTION REPORT" >> 01_Extract_raw_data.log
    echo "  Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> 01_Extract_raw_data.log
    echo "  Samples processed: \$(ls -1 ${logs} | wc -l)" >> 01_Extract_raw_data.log
    echo "================================================================" >> 01_Extract_raw_data.log
    echo "" >> 01_Extract_raw_data.log

    for LOG in \$(ls -1 ${logs} | sort); do
        cat \$LOG >> 01_Extract_raw_data.log
        echo "" >> 01_Extract_raw_data.log
        echo "" >> 01_Extract_raw_data.log
    done
    """
}
