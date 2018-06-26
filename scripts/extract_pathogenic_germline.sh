#!/bin/bash
# for each sample ID defined in -i, write sample ids to -t,
# unzip each vcf file in -d with name {sample_id}{suffix}.vcf, extract
# clinvar candidates, write a properly formatted VCF to folder -o with name
# {sample_id}{suffix}_pathogenic.vcf

#extract_pathogenic_germline.sh \
#-i=/out/human_sequence_prostate_WGS/reproduce/metadata/20171205_sample_attributes.txt \
#-t=/out/human_sequence_prostate_WGS/reproduce/metadata/20171205_sample_attributes_ID_WGS.txt \
#-d=/temp/illumina \
#-o=/temp/illumina \
#-s=_illumina_germline

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
    -p=*|--snpsift=*)
    CALL_SNPSIFT="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--chromosome_names=*)
    CHROMOSOME_NAMES="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
esac
done

if [ ! -f ${CHROMOSOME_NAMES} ]; then
    echo "Could not find chromosome names file ${CHROMOSOME_NAMES}"
    exit 1
fi

if [ ! -f ${FN_SA} ]; then
    echo "Could not find sample attributes file ${FN_SA}"
    exit 1
fi

grep -v incomplete ${FN_SA} | awk 'NR>1' ${FN_SA} | cut -f 1 > ${FN_IDS}
while read SAMPLE_ID; do
    FN=${DIR_IN}/${SAMPLE_ID}${SUFFIX}.vcf
    FN_PATH_FILTERED=${DIR_IN}/${SAMPLE_ID}${SUFFIX}_filtered.vcf
    FN_PATH_CAND=${DIR_IN}/${SAMPLE_ID}${SUFFIX}_pathogenic_candidate.vcf
    FN_PATH=${DIR_IN}/${SAMPLE_ID}${SUFFIX}_pathogenic.vcf
    echo ${SAMPLE_ID}
    if [ ! -f $FN_PATH ]; then
        echo "Did not find existing $FN_PATH, extracting header..."
        head -n 3000 ${FN_PATH_FILTERED} | grep "#" > $FN_PATH_CAND
        #if [ ! -f $FN_PATH_FILTERED ]; then
        #    head -n 3000 ${FN} | grep "#" > $FN_PATH_FILTERED
        #    echo "Did not find existing ${FN_PATH_FILTERED}, filtering to this sample..."
        #    grep -v "0/0" ${FN} | grep -v $'PL\t0:' | grep PASS >> $FN_PATH_FILTERED
        #fi
        echo "filtering pathogenic candidates to ${FN_PATH_CAND}"
        grep -e "pathogenic,1" -e "splice_donor" -e "splice_acceptor" -e "stop_gain" -e "frameshift" ${FN_PATH_FILTERED} >> ${FN_PATH_CAND}
        echo "sifting candidate to produce ${FN_PATH}"
        ${CALL_SNPSIFT} filter \
          -s ${CHROMOSOME_NAMES} \
          -f ${FN_PATH_CAND} "(CHROM in SET[0]) & (( ! exists AF1000G ) | (AF1000G < 0.01)) & (( ! exists EVS ) | (EVS[0] < 0.01))"> ${FN_PATH}
        rm ${FN_PATH_CAND}    
    fi
done <$FN_IDS