# read all vcf files with given --suffix in --dir_in
# write tab-delimited summary file to --out
from optparse import OptionParser
import glob
import sys
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_somatic.vcf")                  
parser.add_option("-l", "--out_list", dest="fn_out_list", \
                  help="path indicating list file to write", \
                  default="list_somatic_PASS_mutations.txt")
parser.add_option("-c", "--known_chromosomes", dest="fn_chrom", \
                  help="path to allowed chromosomes", \
                  default="")
parser.add_option("-o", "--counts", dest="fn_out_counts", \
                  help="path indicating mutation count file to write", \
                  default="matrix_mutation_count_summary.txt")                    
(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
fn_out_list = options.fn_out_list
suffix = options.suffix
fn_chrom = options.fn_chrom
fn_count = options.fn_out_counts

chroms = []
f = open(fn_chrom)
for line in f:
    chrom=line.rstrip('\r\n')
    chroms.append( chrom )

f.close()

s2count = {}

print("Collecting from directory " + dir_in + " files with suffix " + suffix )
fo = open(fn_out_list, 'w')
fo.write( "chrom\tpos\tsample_id\tchrom_order\n")

for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    n_muts = 0
    print("Reading from " + fn)
    for line in f:
        if "PASS" in line:
            a = line.rstrip('\r\n').split('\t')
            chrom, pos = a[0], a[1]
            chrom_order = chrom.replace("chr", "")
            if chrom_order=="X":
                chrom_order="23"
            if chrom_order=="Y":
                chrom_order="24"
            if chrom_order=="M":
                chrom_order="25"   
            if chrom in chroms:
                fo.write( chrom + '\t' + pos + '\t' + sample_id + '\t' + chrom_order + '\n')
                n_muts = n_muts+1
    
    f.close()
    s2count[sample_id] = n_muts

fo.close()

fo_counts = open(fn_count, 'w')
fo_counts.write( "sample_id\tmutation_count\n")
for sample_id in sorted( s2count.keys() ):
    fo_counts.write( sample_id + '\t' + str(s2count[sample_id]) + '\n' )

fo_counts.close()

print("Wrote list of mutations to " + fn_out_list )
print("Wrote mutation counts to " + fn_count )