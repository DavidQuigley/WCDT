suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)
dir_meta='/metadata/human_sequence_prostate_WGS/metadata'
dir_build='/metadata/human_sequence_prostate_WGS/results/build_2018_04_15'
chrom_lengths = read.table(paste(dir_meta,'HG38_chromosome_lengths.txt',sep='/'),
                           row.names=1, stringsAsFactors=FALSE)
chrom_names = rownames(chrom_lengths)
centromeres = read.table( paste(dir_meta,'HG38_centromere_loci.txt',sep='/'), 
                          header=TRUE, stringsAsFactors = FALSE)

chromo_score = function( sample_id="", chrom="", 
                         window_density=10000, window_width=2500000 ){
    
    
    
    total_length=chrom_lengths$V2[ which(rownames(chrom_lengths)==chrom) ]
    half_window=(window_width/2)
    window_starts = seq(from=1, to=total_length+window_width-1, by=window_density)
    interval_start = window_starts - half_window
    interval_end = window_starts+half_window
       
    n_inv = rep(0, length(window_starts))
    n_del = rep(0, length(window_starts))
    n_cna_switch = rep(0, length(window_starts))
    med_CNA = rep(0, length(window_starts)) 
    
    list_sv_m_cur = list_sv_m[list_sv_m$sample_id==sample_id ,]
    bed = get_bed( sample_id, dir_build, allow_small = TRUE )
    bed_chrom = bed[bed$chrom==chrom,]
    
    for( i in 1:length(window_starts) ){
        n_inv[i] = sum( list_sv_m_cur$svtype=="INV" & 
                        list_sv_m_cur$chrom_start==chrom &
                        (
                          (list_sv_m_cur$pos_start>interval_start[i] &
                           list_sv_m_cur$pos_start<interval_end[i]) | 
                          (list_sv_m_cur$pos_end>interval_start[i] &
                           list_sv_m_cur$pos_end<interval_end[i])
                        )
                       )
        n_del[i] = sum( list_sv_m_cur$svtype=="DEL" & 
                        list_sv_m_cur$chrom_start==chrom &
                            (
                                (list_sv_m_cur$pos_start>interval_start[i] &
                                     list_sv_m_cur$pos_start<interval_end[i]) | 
                                    (list_sv_m_cur$pos_end>interval_start[i] &
                                         list_sv_m_cur$pos_end<interval_end[i])
                            )
        )
        n_bounces=0
        bed_idx = which( bed_chrom$start>interval_start[i] &
                             bed_chrom$start<interval_end[i] )
        if( length(bed_idx) > 1 ){
            cnas = bed_chrom$cn.int[bed_idx]
            cur_dir = "greater"
            for(j in 1:(length(cnas)-1)){
                if( cur_dir=="greater" ){
                    if( cnas[j+1] > cnas[j] ){
                        n_cna_switch[i] = n_cna_switch[i]+1
                        cur_dir = "lesser"            
                    }
                }else{
                    if( cnas[j+1] < cnas[j] ){
                        n_cna_switch[i] = n_cna_switch[i]+1
                        cur_dir = "greater"            
                    }
                }
            }
        }
        #med_CNA[i] = median( (bed_chrom$cn.int[bed_idx] * (bed_chrom$stop[bed_idx]-bed_chrom$start[bed_idx]))/window_width, na.rm=TRUE )
    }
    #plot(n_inv, pch=19, cex=0.25, ylim=c(0,40), col="#0000cc33")
    #points(n_del, pch=19, cex=0.25, col="#00cc0033")   
    #points(n_cna_switch, pch=19, cex=0.25, col="#cc000033")   
    data.frame(sample_id=rep(sample_id, length(window_starts)),
               chrom=rep(chrom, length(window_starts)),
               window_starts, 
               n_inv, 
               n_del, 
               n_cna_switch,
               med_CNA)
}

get_bed = function( sample_id, dir='/notebook/human_sequence_prostate_WCDT/WCDT/results/build_2018_04_15', allow_small=FALSE ){
    dir_bed=dir
    bed = read.table(paste(dir_bed,"/",sample_id,'_copycat.bed',sep=''),stringsAsFactors=FALSE)
    names(bed)=c("chrom", "start", "stop", "call", "cn.int", "strand")
    if(!allow_small){
        bed = bed[bed$stop-bed$start > 1000,]
    }
    bed
}

option_list <- list(
	make_option(c("-s", "--sample_id"),     type = "character", help = "Sample identifier. [Required]"),
	make_option(c("-w", "--window_width"),  type = "character", help = "window width. [Required]"),
	make_option(c("-o", "--dir_out"),       type = "character", help = "file name. [Required]")
)
parseobj = OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt = parse_args(parseobj)

sample_id = opt$sample_id 
dir_out = opt$dir_out
window_width= as.numeric(opt$window_width)

list_sv_m = read.table( '/metadata/human_sequence_prostate_WGS/results/2018_04_27_list_sv_m.txt', header=TRUE, stringsAsFactors=FALSE)
fn_pdf = paste(dir_out, '/', sample_id, '_chromothripsis_by_chrom.pdf', sep='')
pdf(fn_pdf, height=8, width=12, paper="USr") 
layout(matrix(1:24,4,6,byrow = TRUE))
for(i in 1:24){
    chrom = chrom_names[i]
    cs = chromo_score( sample_id = sample_id, chrom=chrom, window_width=window_width )
    ymax = max( c( cs$n_inv[cs$chrom==chrom], cs$n_del[cs$chrom==chrom], cs$n_cna_switch[cs$chrom==chrom]))
    if(ymax<10){
        ymax=10
    }
    par(mar=c(4,4,2,0.4))
    plot(cs$n_inv[cs$chrom==chrom], pch=19, cex=0.25, ylim=c(0,ymax), col="#0000cc33", main=paste(sample_id, chrom) )
    points(cs$n_del[cs$chrom==chrom], pch=19, cex=0.25, col="#00cc0033")   
    points(cs$n_cna_switch[cs$chrom==chrom], pch=19, cex=0.25, col="#cc000033")      
    points(cs$med_CNA[cs$chrom==chrom], cex=0.25, col="lightgrey", pch=19)   
    if( i==1 ){
        chromo_scores = cs
    }else{
        chromo_scores = rbind(chromo_scores, cs)
    }
}
dev.off()
chromo_scores$med_CNA = round( chromo_scores$med_CNA, 3 )

fn=paste(dir_out,'/',sample_id,'_chromo_scores_', window_width, '.txt', sep='')
write.table( chromo_scores, fn, sep='\t', quote=FALSE)
