from optparse import OptionParser
import glob
import sys
import subprocess
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-o", "--dir_out", dest="dir_out", \
                  help="directory to which to write files", \
                  default="/temp/illumina")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_copycat.bed")
parser.add_option("-u", "--fn_summary_base", dest="fn_summary_base", \
                  help="base of file to write refseq and symbol summaries", \
                  default="/out/human_sequence_prostate_WGS/reproduce/results/illumina_CNA_summary")
parser.add_option("-b", "--bedtools", dest="bedtools", \
                  help="directory containing bedtools binaries", \
                  default="/opt/bedtools2/bin")
parser.add_option("-g", "--gene_locs", dest="gene_loc_bed", \
                  help="bed file of gene locations to intersect with CNA", \
                  default="/out/human_sequence_prostate_WGS/reproduce/metadata/GRCh38Decoy_refseq_genelocs_from_refFlat.bed")
  
(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
dir_out = options.dir_out
suffix = options.suffix
dir_bedtools = options.bedtools
gene_loc_bed = options.gene_loc_bed
fn_summary_base = options.fn_summary_base
gene_loc_AR_enhancer = dir_out + "/AR_enhancer.bed"
s2g2val = {}
s2g2cn = {}
s2rs2val = {}
s2rs2cn = {}
sample_ids = []
print( "Reading from " + dir_in )
# intersectBed to get copy number calls for each gene locus
for fn in (glob.glob(dir_in  + "/*" + suffix )):
    f = open(fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    sample_ids.append( sample_id )
    print( "Opened " + fn + " for sample " + sample_id )
    s2g2val[ sample_id ] = {}
    s2rs2val[ sample_id ] = {}
    s2g2cn[ sample_id ] = {}
    s2g2cn[ sample_id ] = {}
    fn_out_genes = dir_out + '/' + sample_id + '_copycat_gene_summary.bed'
    f = open(fn_out_genes, "w")
    subprocess.call( [dir_bedtools+"/intersectBed", "-a", fn, "-b", gene_loc_bed, "-wa", "-wb"], stdout=f )
    f.close()
    
    # seperate analysis for AR enhancer region
    fo = open(gene_loc_AR_enhancer, 'w')
    fo.write( "chrX\t66886717\t66927988\tAR_enhancer~AR_enhancer\t0\t+\n" )
    fo.close()
    fn_out_AR_enhancer = dir_out + '/' + sample_id + '_copycat_AR_enhancer.bed'
    f = open(fn_out_AR_enhancer, "w")
    subprocess.call( [dir_bedtools+"/intersectBed", "-a", fn, "-b", gene_loc_AR_enhancer, "-wa", "-wb"], stdout=f )
    f.close()

# extract AR enhancer calls into single table
s2ar = {}
for fn in (glob.glob(dir_out  + "/*_copycat_AR_enhancer.bed" )):
    sample_id = fn.replace("_copycat_AR_enhancer.bed","").replace(dir_out + "/", "")
    f = open( fn )
    for line in f:
        a = line.rstrip('\n\r').split('\t')
        val = float( a[4] )
        if sample_id in s2ar:
            if val > s2ar[ sample_id ]:
                s2ar[ sample_id ] = val
        else:
            s2ar[ sample_id ] = val
    f.close()

fo = open( fn_summary_base + '_AR_enhancer_CN.txt', 'w')
fo.write( "IDENTIFIER\tcopies\n" )
for sample_id in sorted(sample_ids):
    fo.write( sample_id + '\t' + str(s2ar[sample_id]) + '\n' )

fo.close()
print( "wrote AR enhancer calls to " + fn_summary_base + '_AR_enhancer_CN.txt' )

# extract gene calls from intersectBed into hashes
for fn in (glob.glob(dir_out  + "/*_copycat_gene_summary.bed" )):
    sample_id = fn.replace("_copycat_gene_summary.bed","").replace(dir_out + "/", "")
    f = open( fn )
    for line in f:
        a = line.rstrip('\n\r').split('\t')
        symbol, refseq = a[9].split('~')
        if not symbol in s2g2val[sample_id]:
            s2g2val[sample_id][symbol] = a[3]
            s2g2cn[sample_id][symbol] = a[4]
        s2rs2val[sample_id][refseq] = a[3]
        s2g2cn[sample_id][refseq] = a[4]
    
    f.close()

sample_ids = sorted( s2g2val.keys() )
symbols = sorted( s2g2val[sample_ids[0]].keys() )
refseqs = sorted( s2rs2val[sample_ids[0]].keys() )

# Write matrix files
fo = open( fn_summary_base + '_CNA_symbol_copycat.txt', 'w' )
fo.write( "IDENTIFIER\t" + '\t'.join( sample_ids ) + '\n' )
for symbol in symbols:
    out = [symbol]
    for sample_id in sample_ids:
        try:
            out.append( s2g2val[sample_id][symbol] )
        except KeyError:
            out.append( "NA" )
    fo.write( '\t'.join( out ) + '\n' )
fo.close()
print( "Wrote " + fn_summary_base + '_CNA_symbol_copycat.txt' )

fo = open( fn_summary_base + '_CN_integer_symbol_copycat.txt', 'w' )
fo.write( "IDENTIFIER\t" + '\t'.join( sample_ids ) + '\n' )
for symbol in symbols:
    out = [symbol]
    for sample_id in sample_ids:
        try:
            out.append( s2g2cn[sample_id][symbol] )
        except KeyError:
            out.append( "NA" )
    fo.write( '\t'.join( out ) + '\n' )
fo.close()
print( "Wrote " + fn_summary_base + '_CN_integer_symbol_copycat.txt' )
