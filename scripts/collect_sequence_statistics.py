from optparse import OptionParser

parser.add_option("-s", "--sample_attributes", dest="fn_sa", \
                  help="path to sample attributes file", \
                  default="")

parser.add_option("-o", "--dir_out", dest="dir_out", \
                  help="folder in which to write output", \
                  default="")

(options, args) = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)
fn_sa = options.fn_sa
dir_out = options.dir_out

BASEDIR="/opt/BaseSpace/Projects/WGS-ISAS/AppResults"
f = open(fn_sa)
h = f.readline().strip('\r\n').split('\t')
idx_id_n = h.index('ID_illumina_normal')
idx_id_t = h.index('ID_illumina_tumor')
s2illumina = {}
for line in f:
    a = line.strip('\r\n').split('\t')
    sample_id = a[0]
    s2illumina[ sample_id ] = [ {}, {} ]
    id_n = a[idx_id_n]
    id_t = a[idx_id_t]
    print(sample_id)
    fn_csv_n = BASEDIR + '/' + id_n + '/Files/' + id_n + '_S1.summary.csv'
    fn_csv_t = BASEDIR + '/' + id_t + '/Files/' + id_t + '_S1.summary.csv'
    print(fn_csv_n)
    f_n=open( fn_csv_n )
    for fline in f_n:
        k,v = fline.rstrip('\r\n').split(',')
        s2illumina[ sample_id ][0][k] = v
    
    f_n.close()
    print(fn_csv_t)
    f_t=open( fn_csv_t )
    for fline in f_t:
        k,v = fline.rstrip('\r\n').split(',')
        s2illumina[ sample_id ][1][k] = v
    
    f_t.close()

f.close()

sample_ids = sorted( s2illumina.keys())
params = s2illumina[ sample_ids[0] ][ 0 ].keys()
for param in s2illumina[ sample_ids[0] ][ 1 ].keys():
    if not param in params:
        params.append(param)

params = sorted(params)

fo = open( dir_out + '/matrix_alignment_summary_normal.txt', 'w')
fo.write( 'parameter\t' + '\t'.join(sample_ids)  + '\n')
for param in params:
    fo.write( param.replace(' ', '_' ) )
    for sample_id in sample_ids:
        try:
            fo.write( '\t' + s2illumina[ sample_id ][0][param] )
        except KeyError:
            fo.write( '\tNA')
    
    fo.write('\n')

fo.close()

fo = open( dir_out+ '/matrix_alignment_summary_tumor.txt', 'w')
fo.write( 'parameter\t' + '\t'.join(sample_ids)  + '\n')
for param in params:
    fo.write( param.replace(' ', '_' ) )
    for sample_id in sample_ids:
        try:
            fo.write( '\t' + s2illumina[ sample_id ][1][param] )
        except KeyError:
            fo.write( '\tNA')
    
    fo.write('\n')

fo.close()

