from pyfaidx import Fasta
from optparse import OptionParser
import glob
import sys

def identify_microhomology( chrom, pos_start, pos_end, n_shoulder, ref, debug=False  ):    
    # inner5: ACAGAGCGA
    #               \ / pos_start                       \ / pos_end
    # ----sh5--------|-----inner5-----   -----inner3-----|------sh3------------
    #    TCCAGCCTGGGCGACAGAGCGA                 AACCTGGGCAACAGAGCAAGGCTCTGTCTC
    #       ^^^^^^^^^ ~~~~~~~~~                  ^^^^^^^^ ~~~~~~~~
    n_microhomology = 0
    if pos_start > pos_end:
        tmp = pos_start
        pos_start = pos_end
        pos_end = tmp    
    if pos_start - 1 - n_shoulder < 0:
        return 0
    sh5 = str(ref[chrom][ (pos_start-n_shoulder-1) : (pos_start-1)])
    sh3 = str(ref[chrom][ (pos_end) : (pos_end+n_shoulder) ])
    inner5 = str(ref[chrom][ (pos_start) : (pos_start+n_shoulder)])
    inner3 = str(ref[chrom][ (pos_end)-n_shoulder : (pos_end)])
    n_hom_inner5 = 0
    n_hom_inner3 = 0
    for idx5 in range(0, n_shoulder):
        if inner5[idx5] == sh3[idx5]:
            n_hom_inner5 += 1
        else:
            break
    
    for idx3 in range(1, n_shoulder):
        inv_idx = idx3*-1
        if inner3[inv_idx] == sh5[inv_idx]:
            n_hom_inner3 += 1
        else:
            break
    n_microhomology = max( n_hom_inner5, n_hom_inner3 )
            
    return n_microhomology


def extract_PR_SR( a ):
    PR_t=0; SR_t=0; PR_n=0; SR_n=0

    if a[8]=="GT:SU:PE:SR":
        # lumpy
        info_n, info_t = a[10], a[9] # order reversed between manta, lumpy
        x,y,PR_n, SR_n=info_n.split(":")
        x,y,PR_t, SR_t=info_t.split(":")
    else:
        # manta  
        info_n, info_t = a[9], a[10]
        if a[8] == "PR":
            PR_t = info_t.split(",")[1]
            PR_n = info_n.split(",")[1]
        elif a[8]=="PR:SR":
            PR_t, SR_t = info_t.split(":")
            PR_n, SR_n = info_n.split(":")
            PR_t = PR_t.split(",")[1]
            SR_t = SR_t.split(",")[1]
            PR_n = PR_n.split(",")[1]
            SR_n = SR_n.split(",")[1]
    
    return [PR_t, SR_t, PR_n, SR_n]
    

def extract_info( token_str ):
    info = {}
    for token in token_str.split(";"):
        kv = token.split("=")
        if len(kv)==2:
            info[kv[0]] = kv[1]
        else:
            info[kv[0]] = 1
    
    return info

parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/temp/illumina")
parser.add_option("-o", "--out", dest="fn_out_table", \
                  help="path indicating table to write", \
                  default="table_manta_SV.txt")
parser.add_option("-x", "--fn_matrix", dest="fn_del_mh", \
                  help="path indicating summary matrix to write", \
                  default="matrix_manta_SV.txt")
parser.add_option("-f", "--fusion", dest="fn_out_fusion", \
                  help="path indicating fusion table to write", \
                  default="table_manta_fusion.txt")
parser.add_option("-r", "--fn_chromoplexy", dest="fn_out_chromoplexy", \
                  help="path indicating BND table to write for chromoplexy analysis", \
                  default="table_BND.txt")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_manta.vcf")
parser.add_option("-g", "--fn_genome_fa", dest="fn_genome_fa", \
                  help="genome fasta file with .fai index in same directory", \
                  default="/reference/UCSC/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa")
parser.add_option("-m", "--min_match", dest="min_match", \
                  help="minimum number of microhomology matches to declare a hit, default 5", \
                  default=5)
parser.add_option("-l", "--len_shoulder", dest="n_shoulder", \
                  help="number of nucleotides for microhomology search", \
                  default=20)  
parser.add_option("-c", "--chroms", dest="chroms", \
                  help="file containing permitted chromosomes for VCFs", \
                  default="/notebook/human_sequence_prostate_WCDT/WCDT/metadata/chromosome_names.txt")

parser.add_option("-t", '--min_tumor_reads', dest="min_reads", help='minimum tumor read support, default 10', default="10")
parser.add_option("-n", '--max_normal_reads', dest="max_reads", help='maximum normal read support, default 0', default="0")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
fn_out_table = options.fn_out_table
suffix = options.suffix
fn_out_fusion = options.fn_out_fusion
fn_out_chromoplexy = options.fn_out_chromoplexy
min_tumor_reads = int( options.min_reads )
max_normal_reads = int( options.max_reads )
n_shoulder = int(options.n_shoulder)
min_match = int(options.min_match)
fn_genome_fa = options.fn_genome_fa
fn_chroms = options.chroms
fn_del_mh = options.fn_del_mh

ref = Fasta(fn_genome_fa)

chrom_allowed = {}
f = open(fn_chroms)
for line in f:
    chrom_allowed[ line.rstrip('\r\n') ] = 1
f.close()

s2ndel = {}
s2ndup = {}
s2nins = {}
s2ninv = {}
s2ndelmh = {}

print("Collecting from directory " + dir_in + " files with suffix " + suffix )
fo = open(fn_out_table, "w")
fo_fusion = open(fn_out_fusion, "w")
fo_fusion.write( "sample_id\tsymbol_1\tsymbol_2\n" )
fo.write( "sample_id\tsvtype\tsvsize\tnt_of_microhomology\tchrom_start\tpos_start\tchrom_end\tpos_end\n" )

foc = open(fn_out_chromoplexy, "w")
foc.write( "sample\tnum\tchr1\tpos1\tstr1\tchr2\tpos2\tstr2\tsite1\tsite2\n")
bnd_number=1 # to uniquely identify each rearrangement
bnd_seen = {}

for fn in (glob.glob(dir_in  + "/*" + suffix)):
    f = open(fn)
    print("Reading from " + fn)
    sample_id = fn.replace(suffix,"").replace(dir_in + "/", "")    
    s2ndel[sample_id] = 0 
    s2ndup[sample_id] = 0
    s2nins[sample_id] = 0
    s2ninv[sample_id] = 0
    s2ndelmh[sample_id]=0
    for line in f:
        if line[0] != "#":
            a = line.rstrip('\r\n').split('\t')
            chrom, pos_start = a[0], a[1]
            if not chrom in chrom_allowed:
                continue
            
            info = extract_info( a[7] )
            n_microhomology = "NA" # only set for DEL or MantaDEL
            try:
                SVLEN=str(abs(int(info["SVLEN"])))
            except KeyError:
                SVLEN="NA" # not set for BND
            
            if ( ("MGE10kb" in line or "PASS" in line) and "SVTYPE" in info and info["SVTYPE"]=="DUP" and not "MinSomaticScore" in line) or "<DUP>" in line:
                PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                chrom_end = chrom
                pos_end = info["END"]
                if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                    s2ndup[sample_id] = s2ndup[sample_id] + 1
                    fo.write( sample_id + "\tTANDEM\t" + SVLEN + '\t' + n_microhomology + '\t' + chrom + '\t' + pos_start + '\t' + chrom_end + '\t' + pos_end + '\n')
            
            if ("PASS" in line and "SVTYPE" in info and info["SVTYPE"]=="INS") or "<INS>" in line:
                PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                chrom_end = chrom
                pos_end = info["END"]
                if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                    s2nins[sample_id] = s2nins[sample_id] + 1
                    fo.write( sample_id + "\tINS\t" + SVLEN + '\t' + n_microhomology + '\t' + chrom + '\t' + pos_start + '\t' + chrom_end + '\t' + pos_end + '\n')

            if ("PASS" in line and "SVTYPE" in info and info["SVTYPE"]=="INV") or "<INV>" in line:
                PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                chrom_end = chrom
                pos_end = info["END"]
                if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                    s2ninv[sample_id] = s2ninv[sample_id] + 1
                    fo.write( sample_id + "\tINV\t" + SVLEN + '\t' + n_microhomology + '\t' + chrom + '\t' + pos_start + '\t' + chrom_end + '\t' + pos_end + '\n')
            
            if ("PASS" in line and "SVTYPE" in info and info["SVTYPE"]=="DEL" ) or "<DEL>" in line:
                PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                chrom_end = chrom
                pos_end = info["END"]
                debug=False
                n_microhomology = identify_microhomology( chrom, int(pos_start), int(pos_end), n_shoulder, ref, debug )
                if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                    s2ndel[sample_id] = s2ndel[sample_id]+1
                    if n_microhomology>2:
                        s2ndelmh[sample_id] = s2ndelmh[sample_id]+1
                    fo.write( sample_id + '\t' + "DEL" + '\t' + SVLEN + '\t' + str(n_microhomology) + '\t' + chrom + '\t' + pos_start + '\t' + chrom_end + '\t' + pos_end + '\n')
            
            if not "Canvas" in line:
                if info["SVTYPE"]=="BND":
                    PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                    muddle = a[4]
                    strand_start = "0" # forward strand, as understood by ChainFinder
                    if "[" in muddle:
                    	strand_start = "1" 
					strand_end = strand_start # always same as strand_start in Manta
                    try:
                        chrom_end = "chr" + muddle.split(":")[0].split("chr")[1]
                        pos_end = muddle.split(":")[1].split("[")[0].split("]")[0]
                    except IndexError:
                        chrom_end = "NA"
                        pos_end = "NA"
                    if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                        fo.write( sample_id + '\t' + "BND" + '\t' + SVLEN + '\t' + n_microhomology + '\t' + chrom + '\t' + pos_start + '\t' + chrom_end + '\t' + pos_end + '\n')
                    	
                    	# only write once for each translocation, slug is sorted by chromosome alphanumeric order
                    	if chrom > chrom_end: 
                    		bnd_slug = chrom+pos_start+chrom_end+pos_end
                    	else:
                    		bnd_slug = chrom_end+pos_end+chrom+pos_start
                    	if not bnd_slug in bnd_seen and not "random" in bnd_slug and not "Un" in bnd_slug:
	                    	out = [sample_id, str(bnd_number), chrom, pos_start, strand_start, chrom_end, pos_end, strand_end, chrom+":"+pos_start, chrom_end+":"+pos_end]
    	                	foc.write( '\t'.join(out) + '\n')
        	            	bnd_number += 1
        	            	bnd_seen[bnd_slug]=1
            
            if "gene_fusion" in line:
                # fusions will be written to a separate file, fo_fusion
                fusion_partners = []
                PR_t, SR_t, PR_n, SR_n = extract_PR_SR(a)
                if int(PR_t) + int(SR_t) >= min_tumor_reads and int(PR_n) + int(SR_n) <= max_normal_reads:
                    for member in info["CSQT"].split(","):
                        if "NM_" in member:
                            fusion_partners.append( member.split("|")[1])
                    fusion_partners = sorted(fusion_partners)
                    if len(fusion_partners)>1:
                        fo_fusion.write( sample_id + '\t' + fusion_partners[0] + '\t' + fusion_partners[1] + '\n')

fo.close()
fo_fusion.close()

print("Wrote " + fn_out_fusion )
print("Wrote " + fn_out_table )

if suffix=="_manta.vcf":
    fo = open( fn_del_mh, 'w' )
    fo.write( "sample_id\tn_del_manta\tn_del_MHgt2_manta\tn_ins_manta\tn_inv_manta\tn_dup_manta\n")
    for sample_id in sorted( s2ndel.keys() ):
        fo.write( sample_id + '\t' + str(s2ndel[sample_id]) + '\t' + str(s2ndelmh[sample_id]) )
        fo.write( '\t' + str(s2nins[sample_id]) )
        fo.write( '\t' + str(s2ninv[sample_id]) )
        fo.write( '\t' + str(s2ndup[sample_id]) + '\n' )

    fo.close()
    print( "Wrote " + fn_del_mh )

