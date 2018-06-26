#!/bin/bash

for i in "$@"
do
case $i in
    -s=*|--symbol=*)
    SYMBOL="${i#*=}"
    shift # past argument=value
    ;;
    -r=*|--region=*)
    TARGET_REGION="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
esac
done

export BASE=/opt/BaseSpace/Projects/WGS-TNv2/AppResults
export FN_SA=/metadata/human_sequence_prostate_WGS/metadata/20180123_sample_attributes.txt
export SAMTOOLS=/opt/samtools/samtools
export FN_SAMPLE_IDS=/metadata/human_sequence_prostate_WGS/metadata/sample_ids.txt
export DIR_OUT=/metadata/human_sequence_prostate_WGS/results/pindel
export EXTRACT=/metadata/human_sequence_prostate_WGS/scripts/extract_region_from_illumina_bam.py
export GATK_BIN=/opt/GenomeAnalysisTK.jar
export REFERENCE_FA=/reference/UCSC/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa
export N_CORES=4
export PINDEL=/opt/pindel/pindel
export PINDEL2VCF=/opt/pindel/pindel2vcf
export SNPSIFT_BIN=/opt/snpEff/SnpSift.jar
export SNPEFF_BIN=/opt/snpEff/snpEff.jar
export SIFT_FILTER="( GEN[NORMAL_001].GT = '0/0' ) & (GEN[TUMOR_001].GT = '0/1' ) & (GEN[NORMAL_001].AD[1] = 0 ) & (GEN[TUMOR_001].AD[1] > 5 )"

#while read SAMPLE_ID; do
#    echo $SAMPLE_ID
#    FN_N_OUT=${DIR_OUT}/${SAMPLE_ID}_${SYMBOL}_normal.bam
#    FN_T_OUT=${DIR_OUT}/${SAMPLE_ID}_${SYMBOL}_tumor.bam
#    python ${EXTRACT} -f ${FN_SA} -t normal -s ${SAMPLE_ID} -d ${BASE} -o ${FN_N_OUT} -r ${TARGET_REGION}
#    python ${EXTRACT} -f ${FN_SA} -t tumor -s ${SAMPLE_ID} -d ${BASE} -o ${FN_T_OUT} -r ${TARGET_REGION}
#    ${SAMTOOLS} index ${FN_N_OUT} ${FN_N_OUT}.bai
#    ${SAMTOOLS} index ${FN_T_OUT} ${FN_T_OUT}.bai
#done <$FN_SAMPLE_IDS



# Call pindel
# Pindel version 0.2.5b9, 20160729.
while read SAMPLE_ID; do
    FN_CONFIG=${DIR_OUT}/pindel_config_${SAMPLE_ID}_${SYMBOL}
    FN_ROOT=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}
    FN_DEL=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_D
    FN_INS=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_SI
    FN_DEL_VCF=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_D.vcf
    FN_INS_VCF=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_SI.vcf
    FN_COMBINED_VCF=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_combined.vcf
    FN_FILTERED_VCF=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_combined_filtered.vcf
    FN_SNPEFF_VCF=${DIR_OUT}/pindel_${SAMPLE_ID}_${SYMBOL}_combined_filtered_snpeff.vcf

    # make config
    FN_N_OUT=${DIR_OUT}/${SAMPLE_ID}_${SYMBOL}_normal.bam
    FN_T_OUT=${DIR_OUT}/${SAMPLE_ID}_${SYMBOL}_tumor.bam
    echo -e "${FN_T_OUT} 160 TUMOR_001\\n${FN_N_OUT} 160 NORMAL_001" > $FN_CONFIG

    ${PINDEL} -f ${REFERENCE_FA} \
      --config-file ${FN_CONFIG} \
      --chromosome ${TARGET_REGION} \
      --output-prefix ${FN_ROOT} \
      --number_of_threads ${N_CORES}

    # Call pindel2vcf for deletions
    ${PINDEL2VCF} \
        -r ${REFERENCE_FA} -R GRCh38Decoy -d 20170101 \
        -p ${FN_DEL} \
        --min_supporting_reads 3 \
        --max_internal_repeats 2 \
        -G

    ${PINDEL2VCF} \
        -r ${REFERENCE_FA} -R GRCh38Decoy -d 20170101 \
        -p ${FN_INS} \
        --min_supporting_reads 3 \
        --max_internal_repeats 2 \
        -G

    # combine deletion and small insertion VCF file
    # filter out events with any reads in NORMAL_001
    # call snpEff on the combined file
    
    java -cp ${GATK_BIN} org.broadinstitute.gatk.tools.CatVariants \
      -R ${REFERENCE_FA} \
      -V ${FN_DEL_VCF} \
      -V ${FN_INS_VCF} \
      -out ${FN_COMBINED_VCF} \
      -assumeSorted
    
    java -jar ${SNPSIFT_BIN} filter "${SIFT_FILTER}" ${FN_COMBINED_VCF} > ${FN_FILTERED_VCF}
    
    java -jar ${SNPEFF_BIN} eff -noStats \
        -c /opt/snpEff/snpEff.config \
        -v -canon \
        -no-intergenic \
        -no-intron \
        -no SYNONYMOUS_CODING \
        -no INTRAGENIC \
        GRCh38.86 \
        ${FN_FILTERED_VCF} > ${FN_SNPEFF_VCF}
        
done <$FN_SAMPLE_IDS  

# Report out snpEff results for all samples

