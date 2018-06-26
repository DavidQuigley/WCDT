# pull illumina analysis files from BaseSpace
# files to attempt are defined in sample attributes, which must have headers:
# ID_illumina_normal, ID_illumina_tumor, ID_WGS
# by default does not overwrite files that already exist at target

# FASTQ are at
#/opt/BaseSpace/Projects/WGS-TN/AppResults/LP6008367-DNA-F05_WFR2381991-LP6008328-DNA-G12_WFR2343604/Properties/Input.AppResults/0/Properties/Input.Samples/0/Files
#/opt/BaseSpace/Projects/WGS-TN/AppResults/LP6008367-DNA-F05_WFR2381991-LP6008328-DNA-G12_WFR2343604/Properties/Input.AppResults/1/Properties/Input.Samples/0/Files
import os
import shutil 
from optparse import OptionParser
import glob
import sys


def get_file( distant, local, force_download):
    if force_download or not os.path.isfile( local ):
        print( "cp " + distant + "  " + local )
        try:
            shutil.copyfile( distant, local )
        except IOError:
            print( "Error copying " + distant + " to " + local )
            f_fail.write( distant + '\n')

parser = OptionParser()

parser.add_option("-i", "--fn_sa", dest="fn_sa", \
                  help="sample attribute file", \
                  default="/out/human_sequence_prostate_WGS/reproduce/metadata/20171205_sample_attributes.txt")
parser.add_option("-d", "--dir_root", dest="dir_root", \
                  help="base folder for illumina basespace mount", \
              default="/opt/BaseSpace/Projects/WGS-TN/AppResults")
parser.add_option("-r", "--dir_root_rna", dest="dir_root_RNA", \
                  help="base folder for illumina RNA basespace mount", \
                  default="")                  
parser.add_option("-o", "--dir_out", dest="dir_out", \
                  help="directory to which files will be copied", \
                  default="/temp/illumina")
parser.add_option("-f", "--force", dest="force", \
                  help="force copy of file even if already exists", \
                  action="store_true", default=False)

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)                  
fn_sa = options.fn_sa
dir_root = options.dir_root
dir_root_rna = options.dir_root_RNA
force_download = options.force
dir_out = options.dir_out

f_fail = open( dir_out + 'failed_copy_attempts.txt', 'w' )
f = open( fn_sa ) 
header = f.readline().rstrip('\r\n').split('\t')
idx_n_nor = header.index("ID_illumina_normal")
idx_n_tum = header.index("ID_illumina_tumor")
idx_wgs = header.index("ID_WGS")
idx_rna_t = header.index("ID_RNA_T")
for line in f:
    a = line.rstrip('\r\n').split('\t')
    n_nor = a[idx_n_nor]
    n_tum = a[idx_n_tum]
    id_wgs = a[idx_wgs]
    id_rna_t = a[idx_rna_t]
    print( id_wgs )
    if( n_nor == "NA" ):
        continue
    local_header = dir_out + "/" + id_wgs
    fn_genome_illumina= dir_root + "/" + n_nor + "-" + n_tum + "/Properties/Input.AppResults/0/Files/" + n_nor + "_S1.genome.vcf.gz"
    fn_manta= dir_root + '/' + n_nor + "-" + n_tum + '/Files/' + n_nor + "_" + n_tum + "_G1_P1.somatic.SV.vcf.gz"
    fn_somatic = dir_root + '/' + n_nor + "-" + n_tum + '/Files/' + n_nor + "_" + n_tum + "_G1_P1.somatic.vcf.gz"
    
    fn_genome_local = local_header + "_illumina_germline.vcf.gz"
    fn_manta_local = local_header + "_manta.vcf.gz"
    fn_somatic_local = local_header + "_somatic.vcf.gz"
    
    #if force_download or (not os.path.isfile( fn_genome_local) and not os.path.isfile( fn_genome_local.replace(".gz", "") ) ):
    #    print( "creating " + fn_genome_local )
    #    try:
    #        shutil.copyfile( fn_genome_illumina, fn_genome_local )
    #    except IOError:
    #        print( "Error copying " + fn_genome_illumina + " to " + fn_genome_local )
    #    try:
    #        shutil.copyfile( fn_genome_illumina + ".tbi", fn_genome_local + ".tbi")
    #    except IOError:
    #        print( "Error copying " + fn_genome_illumina + ".tbi" + " to " + fn_genome_local + ".tbi" )            
    if force_download or (not os.path.isfile( fn_manta_local) and not os.path.isfile( fn_manta_local.replace(".gz", "") ) ):
        print( "creating " + fn_manta_local )
        try:
            shutil.copyfile( fn_manta, fn_manta_local )
        except IOError:
            print( "Error copying " + fn_manta + " to " + fn_manta_local )
    if force_download or (not os.path.isfile( fn_somatic_local) and not os.path.isfile( fn_somatic_local.replace(".gz", "") ) ):
        print( "creating " + fn_somatic_local )
        try:
            shutil.copyfile( fn_somatic, fn_somatic_local)
        except IOError:
            print( "Error copying " + fn_somatic + " to " + fn_somatic_local )
    
    fn_bam_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/alignments/" + id_rna_t + '.alignments.bam'
    fn_bai_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/alignments/" + id_rna_t + '.alignments.bam.bai'
    fn_bam_rna_t_local = local_header + "_RNA_aligned_tumor.bam"
    fn_bai_rna_t_local = local_header + "_RNA_aligned_tumor.bam.bai"
    
    fn_metrics_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/metrics/metrics.json"
    fn_counts_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/counts/" + id_rna_t + '.counts.genes'
    fn_fusions_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/fusions/" + id_rna_t + '.fusions.csv'
    fn_mantafusions_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/MantaFusions/" + id_rna_t + '.fusions.csv'
    fn_variants_rna_t = dir_root_rna + "/" + id_rna_t + "/Files/variants/" + id_rna_t + '.vcf.gz'
    
    fn_metrics_rna_t_local = local_header + "_RNA_metrics_tumor.json"
    fn_counts_rna_t_local = local_header + "_RNA_counts_tumor.txt"
    fn_fusions_rna_t_local = local_header + "_RNA_fusions_tumor.csv"
    fn_mantafusions_rna_t_local = local_header + "_RNA_mantafusions_tumor.csv"
    fn_variants_rna_t_local = local_header + "_RNA_variants_tumor.vcf.gz"
#    if id_rna_t != "NA":
#        get_file( fn_bam_rna_t, fn_bam_rna_t_local, force_download )
#        get_file( fn_bai_rna_t, fn_bai_rna_t_local, force_download )
#        get_file( fn_metrics_rna_t, fn_metrics_rna_t_local, force_download )
#        get_file( fn_counts_rna_t, fn_counts_rna_t_local, force_download )
#        get_file( fn_fusions_rna_t, fn_fusions_rna_t_local, force_download )
#        get_file( fn_mantafusions_rna_t, fn_mantafusions_rna_t_local, force_download )
#        get_file( fn_variants_rna_t, fn_variants_rna_t_local, force_download )

f_fail.close()