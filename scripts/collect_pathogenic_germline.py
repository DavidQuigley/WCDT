# read all vcf files with given --suffix in --dir_in
# write tab-delimited summary file to --out
from optparse import OptionParser
import glob
import sys
parser = OptionParser()

#default_key="BRCA1,BRCA2,ATM,CDK12,CHEK2,PALB2,RAD51D,ATR,MSH2,MSH6,NBN,GEN1,PMS2,BRIP1,MRE11A,FAM175A"
parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-o", "--out_matrix", dest="fn_out_matrix", \
                  help="path indicating matrix to write key genes only", \
                  default="illumina_germline_summary.txt")
parser.add_option("-l", "--out_list", dest="fn_out_list", \
                  help="path indicating list file to write", \
                  default="illumina_germline_summary.txt")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_illumina_germline_pathogenic.vcf")
parser.add_option("-k", "--key_genes", dest="key_genes", \
                  help="comma-delimited list of key genes for matrix, default to all", \
                  default="")
                  
(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
fn_out_list = options.fn_out_list
fn_out_matrix = options.fn_out_matrix
suffix = options.suffix
key_genes = options.key_genes.split(',')
gene2idx = {}
sample_ids = []
symbols_seen = {}
sample2symbol2seen = {}
fo = open(fn_out_list, "w")
fo.write( "sample_id\tsymbol\trefseq\tvariant_type\tchrom\tpos\tref\talt\tcosmic\tdepth_ref\tdepth_alt\n")
print("Collecting from directory " + dir_in + " files with suffix " + suffix )
for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    print("Reading from " + fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    sample_ids.append(sample_id)
    sample2symbol2seen[ sample_id ] = {}
    for line in f:
        if line[0] != "#":
            a = line.rstrip('\r\n').split('\t')
            chrom, pos, ref, alt, info = a[0], a[1], a[3], a[4], a[7]
            cosmic = []
            symbol="NA"
            refseq="NA"
            vartype="NA"
            depth_ref="NA"
            depth_alt="NA"
            GENO = {}
            for k,v in zip( a[8].split(":"), a[9].split(":") ):
                GENO[k]=v
            depths = GENO["AD"].split(",")
            if len( depths ) == 2:
                depth_ref, depth_alt = depths
            else:
                depth_ref = depths[0]
                depth_alt = str( sum( [int(x) for x in depths[1:]] ) )
            for token in info.split(";"):
                name, val = token.split("=")
                if name=="cosmic":
                    for item in val.split(","):
                        cosmic.append( item.split("|")[1] )
                
                elif name=="CSQT":
                    for item in val.split(","):
                        if "NM_" in item:
                            num, symbol, refseq, vartype = item.split("|")
                        elif "ENST" in item and symbol == "NA":
                            num, symbol, refseq, vartype = item.split("|")            
            cosmic = ','.join(cosmic)
            if cosmic=="":
                cosmic = "NA"
            sample2symbol2seen[sample_id][symbol]=chrom+":"+pos+"_"+ref+">"+alt+"_"+depth_ref+"_"+depth_alt
            symbols_seen[symbol] = 1
            out = [sample_id, symbol, refseq, vartype, chrom, pos, ref, alt, cosmic, depth_ref, depth_alt]
            fo.write( '\t'.join(out) + '\n' )

fo.close()
print( "Wrote list of candidates to " + fn_out_list )
symbols = sorted(symbols_seen.keys())

fo = open(fn_out_matrix, "w")
fo.write( "sample_id\t" + '\t'.join( symbols ) + '\n')
for sample_id in sorted(sample_ids):
    out = [sample_id]
    for symbol in symbols:
        if symbol in sample2symbol2seen[sample_id]:
            out.append(sample2symbol2seen[sample_id][symbol])
        else:
            out.append('.')
    fo.write( '\t'.join(out) + '\n' )

fo.close()
print("Wrote to " + fn_out_matrix)