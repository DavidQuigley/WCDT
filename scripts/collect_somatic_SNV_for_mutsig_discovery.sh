for i in "$@"
do
case $i in
    -d=*|--dir_in=*)
    DIR_IN="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--dir_out=*)
    DIR_OUT="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--fn_sa=*)
    SA="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--bgzip=*)
    BGZIP="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--bcftools=*)
    BCFTOOLS="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--tabix=*)
    TABIX="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
esac
done

while read ID_WGS ID_RNA_T ID_RNA_T2 ID_patient biopsy_site biopsy_group clinic ID_illumina_report illumina_report_status ID_illumina_normal ID_illumina_tumor; do
  if [ ! $ID_WGS == "ID_WGS" ]; then
    echo $ID_WGS
    grep -v -e JTF -e KN7 -e chrUn -e chrEBV -e _KI -e _GL ${DIR_IN}/${ID_WGS}_somatic.vcf | grep -e "#" -e "PASS" | cut -f 1,2,3,4,5,6,7,8,9,11 > ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf
    ${BGZIP} ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf
    ${TABIX} -p vcf ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf.gz
    $BCFTOOLS view -m2 -M2 -v snps ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf.gz > ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean_SNP.vcf
    rm ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf.gz
    rm ${DIR_OUT}/${ID_WGS}_somatic_PASS_clean.vcf.gz.tbi
  fi
done <${SA}
