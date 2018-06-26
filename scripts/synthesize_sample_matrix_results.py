from optparse import OptionParser
import os
import glob
import sys
parser = OptionParser()

parser.add_option("-i", "--files", dest="fn_in", \
                  help="comma-delimited list of matrix files to aggregate", \
                  default="")
parser.add_option("-o", "--out", dest="fn_out", \
                  help="aggregated matrix file to write", \
                  default="sample_summary.txt")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

fn_in = options.fn_in
fn_out = options.fn_out

# confirm all files exist and collect sample names, values, headers, and blanks
# blanks would be used if we are missing a type of data for a given sample
files = fn_in.split(',')
sample_ids = {}
f2sample2vals = {}
f2blank = {}
f2header = {}
for file in files:
    print( "Reading " + file )
    f2sample2vals[file] = {}
    if not os.path.exists(file):
        print("ERROR: cannot find file: " + file)
        sys.exit(1)
    else:
        f = open(file)
        header='\t'.join( f.readline().rstrip('\r\n').split('\t')[1:] )
        f2header[file] = header
        for line in f:
            a=line.rstrip('\r\n').split('\t')
            sample_ids[ a[0] ] = 1
            vals = '\t'.join( a[1:] )
            f2blank[ file ] = '\t'.join( ['NA'] * (len(a)-1) )
            f2sample2vals[ file ][ a[0] ] = vals
        
        f.close()
        
sample_ids = sorted( sample_ids.keys() )
print("Writing summary to " + fn_out )
fo = open( fn_out, 'w' )
fo.write( "IDENTIFIER" )
for file in files:
    fo.write( "\t" + f2header[file] )
fo.write("\n")
for sample_id in sample_ids:
    fo.write( sample_id )
    for file in files:
        if sample_id in f2sample2vals[ file ]:
            fo.write( "\t" + f2sample2vals[file][ sample_id ] )
        else:
            fo.write( "\t" + f2blank[file] )
    fo.write( "\n" )
