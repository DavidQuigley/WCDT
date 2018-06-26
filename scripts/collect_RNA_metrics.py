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
sample2metric2value = {}
metrics = {}
print("Collecting from directory " + dir_in + " files with suffix " + suffix )
for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    print("Reading from " + fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    sample_ids.append(sample_id)
    sample2metric2value[ sample_id ] = {}
    for line in f:
        for token in line.rstrip('\r\n').replace('"', '').rstrip("}").lstrip("{").split(","):
            key, value = token.split(":")
            key = key.replace(" ", "_")
            metrics[key] = 1
            sample2metric2value[ sample_id ][key] = value

metrics = sorted(metrics.keys())
sample_ids = sorted(sample2metric2value.keys())
fo = open(fn_out_matrix, "w")
fo.write( "sample_id\t" + '\t'.join( metrics ) + '\n')
for sample_id in sorted(sample_ids):
    out = [sample_id]
    for metric in metrics:
        if metric in sample2metric2value[sample_id]:
            out.append(sample2metric2value[sample_id][metric])
        else:
            out.append('.')
    fo.write( '\t'.join(out) + '\n' )

fo.close()
print("Wrote to " + fn_out_matrix)