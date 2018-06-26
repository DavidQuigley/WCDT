from optparse import OptionParser
import glob
import sys
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-o", "--out", dest="fn_out_matrix", \
                  help="path indicating matrix to write", \
                  default="rna_metrics.txt")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_RNA_metrics_tumor.json")
                  
(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
fn_out_matrix = options.fn_out_matrix
suffix = options.suffix
sample_ids = []
sample2symbol2value = {}
symbols = {}
print("Collecting from directory " + dir_in + " files with suffix " + suffix )
for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    print("Reading from " + fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    sample_ids.append(sample_id)
    sample2symbol2value[ sample_id ] = {}
    for line in f:
        symbol, counts = line.rstrip('\r\n').split('\t')
        symbols[symbol]=1
        sample2symbol2value[ sample_id ][symbol] = counts

symbols = sorted(symbols.keys())
sample_ids = sorted(sample2symbol2value.keys())
fo = open(fn_out_matrix, "w")
fo.write( "symbol\t" + '\t'.join( sample_ids ) + '\n')
for symbol in sorted(symbols):
    out = [symbol]
    for sample_id in sample_ids:
        if symbol in sample2symbol2value[sample_id]:
            out.append(sample2symbol2value[sample_id][symbol])
        else:
            out.append('.')
    fo.write( '\t'.join(out) + '\n' )

fo.close()
print("Wrote to " + fn_out_matrix)