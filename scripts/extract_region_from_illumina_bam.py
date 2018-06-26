# given sample attributes -f and sample id -s, extract bam of type -t 
# {normal,tumor} from basespace at root directory -d using samtools at -a 
# in the region -r

import os
from optparse import OptionParser
import sys
import subprocess
parser = OptionParser()

parser.add_option("-f", "--fn_sa", dest="fn_sa", \
                  help="path to sample attributes file", \
                  default="")
parser.add_option("-t", "--sample_type", dest="sample_type", \
                  help="path to sample attributes file", \
                  default="")                  
parser.add_option("-s", "--sample_id", dest="sample_id", \
                  help="sample id (ID_wgs)", \
                  default="")
parser.add_option("-d", "--dir_root", dest="dir_root", \
                  help="base folder for illumina basespace mount", \
                  default="/opt/BaseSpace/Projects/WGS-TN-NS6/AppResults")
parser.add_option("-o", "--fn_out", dest="fn_out", \
                  help="path indicating the new file to create", \
                  default="/temp/illumina")
parser.add_option("-a", "--samtools", dest="fn_samtools", \
                  help="path to samtools", \
                  default="/opt/samtools/samtools")
parser.add_option("-r", "--region", dest="region", \
                  help="region to extract", \
                  default="chr22")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)                  
fn_sa = options.fn_sa
sample_type = options.sample_type
sample_id = options.sample_id
dir_root = options.dir_root
fn_out = options.fn_out
fn_samtools = options.fn_samtools
region = options.region
if sample_type != "normal" and sample_type != "tumor":
    print "ERROR: sample_type must be one of {normal,tumor}"
    sys.exit(1)

f = open( fn_sa )
header = f.readline().rstrip('\r\n').split('\t')
idx_report = header.index("ID_illumina_report")
idx_n_tum = header.index("ID_illumina_tumor")
idx_n_nor = header.index("ID_illumina_normal")
idx_wgs = header.index("ID_WGS")
report = ""
for line in f:
    a = line.rstrip('\r\n').split('\t')
    if a[idx_wgs]==sample_id:
        report = a[idx_report]
        if sample_type=="normal":
            sample_string = a[idx_n_nor]
        else:
            sample_string = a[idx_n_tum]

f.close()

if report == "":
    print "ERROR: report not found for " + sample_id + " in sample attributes file " + fn_sa
    sys.exit(1)

if sample_type=="normal":
    dir_no="0"
else:
    dir_no="1"
fn_bam = dir_root + "/" + report + "/Properties/Input.AppResults/" + dir_no + "/Files/" + sample_string + "_S1.bam"

if os.path.exists( fn_bam ):     
    print( " ".join( [fn_samtools, "view", "-b", "-h", fn_bam, region, "-o" + fn_out ] ) )
    subprocess.call( [fn_samtools, "view", "-b", "-h", fn_bam, region, "-o" + fn_out ] )
    print( " ".join( [fn_samtools, "index", fn_out, fn_out + ".bai" ] ) )
    subprocess.call( [fn_samtools, "index", fn_out, fn_out + ".bai" ] )
else:
    print( "ERROR: could not find file to read: " + fn_bam )