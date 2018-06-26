#!/bin/bash

for i in "$@"
do
case $i in
    -i=*|--fn_sample_attributes=*)
    FN_SA="${i#*=}"
    shift # past argument=value
    ;;
    -d=*|--dir_in=*)
    DIR_IN="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--dir_out=*)
    DIR_OUT="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--fn_ids=*)
    FN_IDS="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--suffix=*)
    SUFFIX="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--build_id=*)
    SUFFIX="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--snpsift=*)
    CALL_SNPSIFT="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--chromosome_names=*)
    CHROMOSOME_NAMES="${i#*=}"
    shift # past argument=value
    ;;    
    -k=*|--snpsift=*)
    STRELKA_MUTATIONSIGS="${i#*=}"
    shift # past argument=value
    ;;
    -v=*|--snp_counting=*)
    SNP_COUNTING="${i#*=}"
    shift # past argument=value
    ;;    
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
esac
done

echo "DIR_IN: ${DIR_IN}" 
grep -v incomplete ${FN_SA} | awk 'NR>1' | cut -f 1 > ${FN_IDS}
while read SAMPLE_ID; do
    echo "Somatic filtering for ${SAMPLE_ID}"
    FN=${DIR_IN}/${SAMPLE_ID}${SUFFIX}.vcf 
    FN_OUT=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_pass.vcf
    FN_TEMP=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_temp.vcf
    FN_OUT_SNP=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_pass_SNP.vcf
    FN_INACTIVE=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_pass_inactivating.vcf
    FN_MISSENSE=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_pass_missense.vcf
    FN_OUT_CPRA=${DIR_OUT}/${SAMPLE_ID}${SUFFIX}_pass_SNP_CPRA.txt # for mutsig
    
    echo "  Calling SNPSift"
    ${CALL_SNPSIFT} filter \
      -s ${CHROMOSOME_NAMES} \
      -f ${FN} "(FILTER = 'PASS') & (CHROM in SET[0] )"> ${FN_OUT}      
    grep "^#" ${FN_OUT} > ${FN_INACTIVE}  # create header
    grep "^#" ${FN_OUT} > ${FN_MISSENSE}  # create header
    grep -e frameshift -e splice_acceptor -e splice_donor -e stop_gain ${FN_OUT} >> ${FN_INACTIVE}
    grep -e missense ${FN_OUT} >> ${FN_MISSENSE}
    
    if [ "$SNP_COUNTING" == "1" ]; then
        ${CALL_SNPSIFT} varType ${FN_OUT} > ${FN_TEMP}
        ${CALL_SNPSIFT} filter -f ${FN_TEMP} "(VARTYPE='SNP')" > ${FN_OUT_SNP}
        rm ${FN_TEMP}
        ${CALL_SNPSIFT} extractFields ${FN_OUT_SNP} CHROM POS REF ALT > ${FN_OUT_CPRA}
    
        if [ ! -f ${DIR_OUT}/${SAMPLE_ID}_strelka_mutation_percentages.txt ]; then       
            echo "  Calling mutation signatures"
            Rscript ${STRELKA_MUTATIONSIGS} \
                --sample_id ${SAMPLE_ID} \
                --fn_CPRA ${FN_OUT_CPRA} \
                --out ${DIR_OUT}/${SAMPLE_ID}_strelka_mutation_percentages.txt
        fi
    fi 
done <$FN_IDS