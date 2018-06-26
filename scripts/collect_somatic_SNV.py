from optparse import OptionParser
import glob
import sys
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/out/human_sequence_prostate_WGS/reproduce/results")
parser.add_option("-o", "--dir_out", dest="dir_out", \
                  help="directory to which to write files", \
                  default="/out/human_sequence_prostate_WGS/reproduce/results")
parser.add_option("-s", "--suffix_inactivating", dest="suffix_inactivating", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_pass_inactivating.vcf")
parser.add_option("-m", "--suffix_missense", dest="suffix_missense", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_pass_missense.vcf")
parser.add_option("-t", "--tool", dest="mutation_tool", \
                  help="tool used for mutation analyis, one of {strelka,mutect}", \
                  default="_pass_missense.vcf")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
dir_out = options.dir_out
mutation_tool = options.mutation_tool
if mutation_tool != "strelka" and mutation_tool != "mutect":
    print("parameter mutation_tool must be one of strelka,mutect")
    sys.exit(1)
    
fn_out_inactive = dir_out + '_matrix_somatic_inactivating_' + mutation_tool + '.txt'
fn_out_missense = dir_out + '_matrix_somatic_missense_' + mutation_tool + '.txt'
fn_out_somatic_list = dir_out + '_list_somatic_' + mutation_tool + '.txt'
suffix_missense = options.suffix_missense
suffix_inactivating = options.suffix_inactivating

    
fo_list = open(fn_out_somatic_list, 'w')
fo_list.write( "sample_id\tsymbol\tvariant_type\tchrom\tpos\treference\talt\tdepth_ref\tdepth_alt\n")
samp2symbol2muts = {}
symbols = {}
for fn in (glob.glob(dir_in  + "/*" + suffix_inactivating )):
    f = open(fn)
    sample_id = fn.replace(dir_in + '/',"").split("_")[0]
    samp2symbol2muts[sample_id] = {}
    if mutation_tool=="strelka" and "mutect" in fn:
        continue
    if mutation_tool=="mutect" and not "mutect" in fn:
        continue            
    print( "Read inactivating sample " + sample_id + " in file " + fn )
    for line in f:
        if line[0]=="#":
            continue
        
        a = line.rstrip('\n\r').split('\t')
        info = a[7].split(";")
        chrom, pos,ref,alt = a[0], a[1],a[3], a[4]
        depth_ref, depth_alt = "NA", "NA"
        GN, GT = {}, {}
        for k,v in zip( a[8].split(":"), a[9].split(":") ):
            GN[k]=v
        for k,v in zip( a[8].split(":"), a[10].split(":") ):
            GT[k]=v            
        depth_ref, depth_alt = GN["DP"], GT["DP"]
        symbol="NA"
        consequence="NA"
        for token in info:
            try:
                k,v = token.split("=")
                if k=="CSQT":
                    for blurb in v.split(","):
                        n, symbol, accession, consequence = blurb.split("|")
                        if "NM_" in accession:
                            blub = chrom+":"+pos+"_"+ref+">"+alt+"_"+consequence+"_"+depth_ref+"_"+depth_alt
                            if symbol in samp2symbol2muts[sample_id]:
                                samp2symbol2muts[sample_id][symbol].append( blub )
                            else:
                                samp2symbol2muts[sample_id][symbol] = [ blub ]
                            symbols[symbol] = 1
                            out = [sample_id, symbol, consequence, chrom, pos, ref, alt, depth_ref, depth_alt]
                            fo_list.write( '\t'.join(out) + '\n' )
                
                elif k=="ANN": # mutect/snpEff pipeline
                    #G|missense_variant|MODERATE|NUTM2B|ENSG00000188199|transcript|ENST00000429828.5|protein_coding|7/7|c.1993C>G|p.Gln665Glu|2376/3292|1993/2637|665/878||
                    blurbs = v.split("|")
                    symbol, accession, consequence =blurbs[3], blurbs[4], blurbs[1]
                    blub = chrom+":"+pos+"_"+ref+">"+alt+"_"+consequence+"_"+depth_ref+"_"+depth_alt
                    if symbol in samp2symbol2muts[sample_id]:
                        samp2symbol2muts[sample_id][symbol].append( blub )
                    else:
	                    samp2symbol2muts[sample_id][symbol] = [blub]
                    out = [sample_id, symbol, consequence, chrom, pos, ref, alt, depth_ref, depth_alt]
                    fo_list.write( '\t'.join(out) + '\n' )
                    symbols[symbol] = 1
            
            except ValueError:
                pass
    
    f.close()

symbols = sorted(symbols.keys())
fo = open(fn_out_inactive, 'w')
fo.write( "IDENTIFIER\t" + '\t'.join(symbols) + '\n' )
for sample_id in sorted(samp2symbol2muts.keys()):
    vals = []
    for symbol in symbols:
        if symbol in samp2symbol2muts[sample_id]:
            vals.append( '|'.join( samp2symbol2muts[sample_id][symbol] ) )
        else:
            vals.append('.')
    
    fo.write( sample_id + '\t' + '\t'.join(vals) + '\n')

fo.close()
print("Wrote " + fn_out_inactive )

samp2symbol2muts = {}
symbols ={}
for fn in (glob.glob(dir_in  + "/*" + suffix_missense )):
    f = open(fn)
    sample_id = fn.replace(dir_in + '/',"").split("_")[0]
    if mutation_tool=="strelka" and "mutect" in fn:
        continue
    if mutation_tool=="mutect" and not "mutect" in fn:
        continue    
    print( "Read missense sample " + sample_id + " in file " + fn )
    samp2symbol2muts[sample_id] = {}
    for line in f:
        if line[0]=="#":
            continue
        
        a = line.rstrip('\n\r').split('\t')
        info = a[7].split(";")
        chrom, pos,ref,alt = a[0], a[1],a[3], a[4]
        depth_ref, depth_alt = "NA", "NA"
        GN, GT = {}, {}
        for k,v in zip( a[8].split(":"), a[9].split(":") ):
            GN[k]=v
        for k,v in zip( a[8].split(":"), a[10].split(":") ):
            GT[k]=v            
        depth_ref, depth_alt = GN["DP"], GT["DP"]
        
        symbol="NA"
        consequence="NA"
        for token in info:
            try:
                k,v = token.split("=")
                if k=="CSQT": # strelka from illumina pipeline
                    for blurb in v.split(","):
                        n, symbol, accession, consequence = blurb.split("|")
                        if "NM_" in accession:
                            blub = chrom+":"+pos+"_"+ref+">"+alt+"_"+consequence+"_"+depth_ref+"_"+depth_alt
                            if symbol in samp2symbol2muts[sample_id]:
                                samp2symbol2muts[sample_id][symbol].append( blub )
                            else:
	                            samp2symbol2muts[sample_id][symbol] = [blub]
                            out = [sample_id, symbol, consequence, chrom, pos, ref, alt, depth_ref, depth_alt]
                            fo_list.write( '\t'.join(out) + '\n' )
                            symbols[symbol] = 1
                elif k=="ANN": # mutect/snpEff pipeline
                    #G|missense_variant|MODERATE|NUTM2B|ENSG00000188199|transcript|ENST00000429828.5|protein_coding|7/7|c.1993C>G|p.Gln665Glu|2376/3292|1993/2637|665/878||
                    blurbs = v.split("|")
                    symbol, accession, consequence =blurbs[3], blurbs[4], blurbs[1]
                    blub = chrom+":"+pos+"_"+ref+">"+alt+"_"+consequence+"_"+depth_ref+"_"+depth_alt
                    if symbol in samp2symbol2muts[sample_id]:
                        samp2symbol2muts[sample_id][symbol].append( blub )
                    else:
	                    samp2symbol2muts[sample_id][symbol] = [blub]
                    out = [sample_id, symbol, consequence, chrom, pos, ref, alt, depth_ref, depth_alt]
                    fo_list.write( '\t'.join(out) + '\n' )
                    symbols[symbol] = 1
            
            except ValueError:
                pass
    f.close()

symbols = sorted(symbols.keys())
fo = open(fn_out_missense, 'w')
fo.write( "IDENTIFIER\t" + '\t'.join(symbols) + '\n' )
for sample_id in sorted(samp2symbol2muts.keys()):
    vals = []
    for symbol in symbols:
        if symbol in samp2symbol2muts[sample_id]:
            vals.append( '|'.join( samp2symbol2muts[sample_id][symbol] ) )
        else:
             vals.append('.')
    
    fo.write( sample_id + '\t' + '\t'.join(vals) + '\n')

fo.close()

fo_list.close()
print("Wrote " + fn_out_missense )
print("Wrote " + fn_out_somatic_list )

