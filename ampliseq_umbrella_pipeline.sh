for input_dir in /g/korbel/weber/TMP/VINCENT_GROUP/Samples/112548392SC2; do
    echo "$input_dir"
    nextflow run main.nf \
        -profile mamba \
        --outdir "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus" \
        --input "$input_dir/Sequencing" \
        --extension "/*.{1,2}.fastq.gz" \
        --dada_ref_taxonomy "silva=138" \
        --ignore_failed_trimming \
        --ignore_empty_input_files \
        --FW_primer "GTGYCAGCMGCCGCGGTAA" \
        --RV_primer "CCGYCAATTYMTTTRAGTTT" \
        --concatenate_reads "consensus" \
        --skip_dada_addspecies \
        -resume


    mv "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus/qiime2" \
        "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus/qiime2_silva"

    nextflow run main.nf \
        -profile mamba \
        --outdir "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus" \
        --input "$input_dir/Sequencing" \
        --extension "/*.{1,2}.fastq.gz" \
        --dada_ref_taxonomy pr2 \
        --ignore_failed_trimming \
        --ignore_empty_input_files \
        --FW_primer "GTGYCAGCMGCCGCGGTAA" \
        --RV_primer "CCGYCAATTYMTTTRAGTTT" \
        --concatenate_reads consensus \
        --skip_dada_addspecies \
        -resume

    mv "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus/qiime2" \
        "$input_dir/nfcore-ampliseq-results/nf-core-ampliseq-v261-2annot_pr2_silva-mergePairs_consensus/qiime2_pr2"

done
