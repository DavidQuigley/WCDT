# read all vcf files with given --suffix in --dir_in
# write tab-delimited summary file to --out
from optparse import OptionParser
import glob
import sys
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-o", "--out", dest="fn_out", \
                  help="path indicating file to write", \
                  default="illumina_germline_summary.txt")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_illumina_germline_pathogenic.vcf")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
fn_out = options.fn_out
suffix = options.suffix

s2vals = {}
signature_names = []
is_first=True
for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    s2vals[sample_id] = []
    header = f.readline()
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        if is_first:
            signature_names.append(a[0].replace("Signature.", "S"))
        s2vals[sample_id].append(a[1])
    is_first=False
    f.close()

samples = sorted(s2vals.keys())
fo = open(fn_out, "w")
fo.write( "sample_id\t" + '\t'.join(signature_names) + '\n')
for sample_id in sorted(s2vals.keys()):
    out = [sample_id]
    for val in s2vals[sample_id]:
        out.append(val)
    fo.write( '\t'.join( out ) + '\n' )

fo.close()

