#read *_manta.vcf
#1) convert all canvas entries into BED file _canvas.bed
#2) convert all manta entries into CIRCOS-ready format
#3) report summary data for sample: sample_id, ploidy, coverage, percent CNA, 
#               number of INS, DEL, INV, BND, bases_CNA_ref, bases_cna_notref
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
                  default="_manta.vcf")
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

s2g2val = {} # sample to gene to value
s2rs2val = {} # sample to refseq to value
s2g2cn = {} # sample to gene to copy number
s2rs2cn = {} # sample to refseq to copy number
sample2chrom2dels = {}
sample2chrom2trans = {}

valid_chroms = []
for i in range(1,23):
    valid_chroms.append("chr" + str(i) )
valid_chroms.append("chrX")
valid_chroms.append("chrY")
valid_chroms_sorted = valid_chroms
valid_chroms = set(valid_chroms)

sample_ids = []

print( "Reading from " + dir_in )
for fn in (glob.glob(dir_in  + "/*" + suffix )):
    f = open(fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")
    sample_ids.append( sample_id )
    print( "Opened " + fn + " for sample " + sample_id )
    s2g2val[ sample_id ] = {}
    s2rs2val[ sample_id ] = {}
    s2g2cn[ sample_id ] = {}
    s2g2cn[ sample_id ] = {}
    sample2chrom2dels[ sample_id ] = {}
    sample2chrom2trans[ sample_id ] = {}
    for chrom in valid_chroms:
        sample2chrom2dels[ sample_id ][ chrom ] = 0
        sample2chrom2trans[ sample_id ][ chrom ] = 0
    
    fn_out_bed = dir_out + '/' + sample_id + '_canvas.bed'
    fn_out_circos = dir_out + '/' + sample_id + '_circos.txt'
    fn_out_summary = dir_out + '/' + sample_id + '_CNA_summary.txt'
    fn_out_genes = dir_out + '/' + sample_id + '_gene_summary.bed'
    fo_bed = open(fn_out_bed, "w")
    fo_circos = open(fn_out_circos, "w")
    fo_circos.write( "Chromosome\tchromStart\tchromEnd\tChromosome.1\tchromStart.1\tchromEnd.1\n")
    n_bnd=0
    n_del=0
    n_inv=0
    n_ins=0
    sum_CNA_ref = 0.0
    sum_CNA_notref = 0.0
    ploidy = "NA"
    coverage = "NA"
    
    for line in f:
        if line[0] == '#':
            if "OverallPloidy" in line:
                ploidy = line.rstrip('\r\n').split("=")[1]
            elif "DiploidCoverage" in line:
                coverage = line.rstrip('\r\n').split("=")[1]
        else:
            a = line.rstrip('\r\n').split('\t')
            if "Canvas" in a[2]:
                canvas, ctype, chrom, positions = a[2].split(":")
                RC, BC, CN, MCC = a[10].split(":")
                pos_start, pos_end = positions.split("-") 
                fo_bed.write( chrom + '\t' + "\t".join( [pos_start, pos_end] )  + '\t' + ctype + '\t' + CN + '\t' + '+\n' )
                if ctype == "REF":
                    sum_CNA_ref += float( pos_end ) - float( pos_start )
                else:
                    sum_CNA_notref += float( pos_end ) - float( pos_start )
            elif "Manta" in a[2]:
                if "PASS" in line: 
                    if "IMPRECISE" not in line and "MantaBND" in line:
                        chr_1=a[0]        
                        start_1 = a[1]
                        end_1 = str(int(a[1])+1)
                        muddle = a[4]
                        if not "chr" in muddle or not "chr" in chr_1: 
                            continue
                        n_pr, n_sr = a[9].split(':')
                        t_pr, t_sr = a[10].split(':')
                        n_pr = n_pr.split(",")[1]
                        n_sr = n_sr.split(",")[1]
                        t_pr = t_pr.split(",")[1]
                        t_sr = t_sr.split(",")[1]
                        chr_2 = "chr" + muddle.split(":")[0].split("chr")[1]
                        start_2 = muddle.split(":")[1].split("[")[0].split("]")[0]
                        end_2 = str(int(start_2)+1)
                        if int(n_pr)==0 and int(n_sr)==0 and int(t_sr)>10:
                            fo_circos.write( '\t'.join( [chr_1, start_1, end_1, chr_2, start_2, end_2] ) + '\n')
                            if chr_1 in valid_chroms:
                                sample2chrom2trans[sample_id][ chr_1 ] = sample2chrom2trans[sample_id][ chr_1 ] + 1 
                            if chr_2 in valid_chroms:
                                sample2chrom2trans[sample_id][ chr_2 ] = sample2chrom2trans[sample_id][ chr_2 ] + 1 
                            
                    if "MantaINS" in a[2]:
                        n_ins += 1
                    if "MantaDEL" in a[2]:
                        n_del += 1
                        if a[0] in valid_chroms:
                            sample2chrom2dels[sample_id][ a[0] ] = sample2chrom2dels[sample_id][ a[0] ] + 1 
                    
                    if "MantaINV" in a[2]:
                        n_inv += 1
                    if "MantaBND" in a[2]:
                        n_bnd += 1
    
    fo_summary = open(fn_out_summary, "w")
    fo_summary.write( '\t'.join(["sample_id", "ploidy", "coverage", "bases_CNA_ref", "bases_cna_notref", "percent_CNA_ref"]) + '\n')
    percent_CNA_ref = round( sum_CNA_ref / (sum_CNA_ref + sum_CNA_notref), 5)
    fo_summary.write( '\t'.join([str(x) for x in [sample_id, ploidy, coverage, sum_CNA_ref, sum_CNA_notref, percent_CNA_ref]]) + '\n')
    fo_summary.close()
    fo_circos.close()
    fo_bed.close()
    
    f = open(fn_out_genes, "w")
    subprocess.call( [dir_bedtools+"/intersectBed", "-a", fn_out_bed, "-b", gene_loc_bed, "-wa", "-wb"], stdout=f )
    f.close()

for fn in (glob.glob(dir_out  + "/*_gene_summary.bed" )):
    sample_id = fn.replace("_gene_summary.bed","").replace(dir_out + "/", "")
    if not "copycat" in sample_id:
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
fo = open( fn_summary_base + '_CNA_symbol.txt', 'w' )
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
print( "Wrote " + fn_summary_base + '_CNA_symbol.txt' )

fo = open( fn_summary_base + '_CNA_refseq.txt', 'w' )
fo.write( "IDENTIFIER\t" + '\t'.join( sample_ids ) + '\n' )
for refseq in refseqs:
    out = [refseq]
    for sample_id in sample_ids:
        try:
            out.append( s2rs2val[sample_id][refseq] )
        except KeyError:
            out.append( "NA" )
    fo.write( '\t'.join( out ) + '\n' )
fo.close()
print( "Wrote " + fn_summary_base + '_CNA_refseq.txt' )


fo = open( fn_summary_base + '_CN_integer_symbol.txt', 'w' )
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
print( "Wrote " + fn_summary_base + '_CN_integer_symbol.txt' )

fo = open( fn_summary_base + '_CN_integer_refseq.txt', 'w' )
fo.write( "IDENTIFIER\t" + '\t'.join( sample_ids ) + '\n' )
for refseq in refseqs:
    out = [refseq]
    for sample_id in sample_ids:
        try:
            out.append( s2g2cn[sample_id][refseq] )
        except KeyError:
            out.append( "NA" )
    fo.write( '\t'.join( out ) + '\n' )
fo.close()

print( "Wrote " + fn_summary_base + '_CN_integer_refseq.txt' )


fo = open( fn_summary_base + '_CNA_summary_statistics.txt', 'w' )
wrote_header=False
for fn in sorted( (glob.glob(dir_out  + "/*_CNA_summary.txt" )) ):
    f = open( fn )
    print( fn )
    if not wrote_header:
        fo.write( f.readline() )
        fo.write( f.readline() )
        wrote_header=True
    else:
        header = f.readline() 
        fo.write( f.readline() )
    
    f.close()

fo.close()

print( "Wrote " + fn_summary_base + '_CNA_summary_statistics.txt' )

fo = open( fn_summary_base + '_dels_by_chrom.txt', 'w' )
fo.write( "sample_id" + '\t' + '\t'.join( valid_chroms_sorted ) + '\n' )
for sample_id in sorted( sample_ids ):
    out = [ sample_id ]
    for chrom in valid_chroms_sorted:
        out.append( sample2chrom2dels[sample_id][chrom] )
    
    fo.write( '\t'.join( [str(x) for x in out] ) + '\n' )

fo.close()

fo = open( fn_summary_base + '_trans_by_chrom.txt', 'w' )
fo.write( "sample_id" + '\t' + '\t'.join( valid_chroms_sorted ) + '\n' )
for sample_id in sorted( sample_ids ):
    out = [ sample_id ]
    for chrom in valid_chroms_sorted:
        out.append( sample2chrom2trans[sample_id][chrom] )
    
    fo.write( '\t'.join( [str(x) for x in out] ) + '\n' )

fo.close()

