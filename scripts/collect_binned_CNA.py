# binned weighted summary

from optparse import OptionParser
import glob
import sys
import subprocess
import math
parser = OptionParser()

parser.add_option("-d", "--dir_in", dest="dir_in", \
                  help="directory from which to read files", \
                  default="/out/human_sequence_prostate_WGS/reproduce/results")
parser.add_option("-s", "--suffix", dest="suffix", \
                  help="suffix defining files to read, {SAMPLE_ID}suffix", \
                  default="_canvas.bed")
parser.add_option("-o", "--dir_out", dest="dir_out", \
                  help="directory to which to write files", \
                  default="/out/human_sequence_prostate_WGS/reproduce/results")
parser.add_option("-c", "--chrom_lengths", dest="fn_chrom", \
                  help="file with chromosome lengths", \
                  default="/out/human_sequence_prostate_WGS/reproduce/metadata/HG38_chromosome_lengths.txt")
parser.add_option("-b", "--bedtools", dest="bedtools", \
                  help="directory containing bedtools binaries", \
                  default="/opt/bedtools2/bin")
parser.add_option("-w", "--window_width", dest="window_width", \
                  help="directory containing bedtools binaries", \
                  default="3000000")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

dir_in = options.dir_in
dir_out = options.dir_out
suffix = options.suffix
fn_chrom = options.fn_chrom
dir_bedtools = options.bedtools
window_width = int(options.window_width)
if "titan" in suffix:
    tool = "TITAN"
elif "canvas" in suffix:
    tool = "CANVAS"
else:
    tool = "COPYCAT"
print( "Autodetected tool: " + tool)

# create windowed bed file to intersect with segments
f = open(fn_chrom)
fn_chrom_bed = dir_out + '_HG38_3mb_windows.bed'
print("Creating windowed bed file " + fn_chrom_bed )
fo = open(fn_chrom_bed, 'w')
for line in f:
    chrom, length = line.rstrip('\r\n').split('\t')
    length = int(length)
    begin = 1
    end = begin + window_width
    while end < length:
        fo.write( chrom + "\t" + str(begin) + '\t' + str(end) + '\n' )
        begin += window_width
        end += window_width-1
    
    fo.write( chrom + "\t" + str(begin) + '\t' + str(end) + '\n' )

f.close()
fo.close()

print( "Reading from " + dir_in )
is_first=True
sample2segment2vals = {}
sample2segment2sum = {}
for fn in (glob.glob(dir_in  + "/*" + suffix )):
    sample_id = fn.replace(dir_in + '/','').replace(suffix, "")
    sample2segment2vals[ sample_id ] = {}
    sample2segment2sum[ sample_id ] = {}
    if tool=="CANVAS":
        fn_out = fn.replace( "_canvas.bed", "_canvas_window.bed")
    elif tool == "TITAN":
        fn_out = fn.replace( "_titan_segments.bed", "_titan_segments_window.bed")
    else:
        fn_out = fn.replace( "_copycat.bed", "_copycat_segments_window.bed")
    print( "Bedtools on sample " + sample_id + " from file " + fn )
    cmd = [dir_bedtools + "/intersectBed", "-loj", "-a", fn_chrom_bed, "-b", fn]
    fo = open(fn_out, "w")
    subprocess.call( cmd, stdout=fo)
    fo.close()
    f = open( fn_out )
    if is_first:
        segments = []
    
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        segment = a[0] + ":" + a[1]
        if not segment in sample2segment2vals[ sample_id ]:
            sample2segment2vals[ sample_id ][segment] = []
            sample2segment2sum[ sample_id ][segment] = 0
            if is_first:   
                segments.append( segment )
        
        if int(a[4])<int(a[1]):
            left = int(a[1])
        
        else:
            left = int(a[4])
        
        if int(a[5]) < int(a[2]):
            right = int(a[5])
        
        else:
            right = int(a[2])
        
        span = right-left
        if tool=="CANVAS" or tool=="COPYCAT":
            copy_number = float(a[7])
        else:
            try:
                logr = float(a[6].split('~')[2])
                copy_number = round( 2*math.pow(2, logr), 3)
            except IndexError:
                copy_number=-1

        if copy_number == -1:
            sample2segment2vals[sample_id][ segment ]="NA"
        
        else:
            sample2segment2vals[ sample_id ][ segment ].append( copy_number * span )
            sample2segment2sum[ sample_id ][ segment ] += (right-left)
    
    if is_first:
        is_first = False
    
    f.close()

samples = sorted( sample2segment2vals.keys() )
if tool=="CANVAS":
    fn_out = dir_out + '_matrix_binned_weighted_CN_canvas.txt'
elif tool=="TITAN":
    fn_out = dir_out + '_matrix_binned_weighted_CN_titan.txt'
else:
    fn_out = dir_out + '_matrix_binned_weighted_CN_copycat.txt'
print("Writing binned weight matrix " + fn_out)

fo = open(fn_out, 'w')
fo.write( "chrom\tbin_start\t" + '\t'.join( samples ) + '\n')
for segment in segments:
    fo.write( '\t'.join( segment.split(":"))  )
    for sample_id in samples:
        if sample2segment2vals[sample_id][segment] == "NA":
            fo.write( '\tNA')
    
        else:
            weighted_sum = sum(sample2segment2vals[sample_id][segment])
            total_span = sample2segment2sum[sample_id][segment]
            fo.write( '\t' + str( round( weighted_sum / total_span, 1) ))
    
    fo.write('\n')

fo.close()

