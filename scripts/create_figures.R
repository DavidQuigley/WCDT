library(SomaticSignatures)
library(ggplot2)
library(plyr)
library(deconstructSigs)
source('/notebook/code/src/R/quantitative_genetics.R')
source('/notebook/human_sequence_prostate_WCDT/WCDT/figures/prepare_data_for_secondary_analysis.R')
dir_out = '/notebook/human_sequence_prostate_WCDT/WCDT/figures'

################################################################################
#                                                                              #
# STATISTICAL COMPARISONS                                                      #
#                                                                              #
################################################################################

# ploidy vs. SV
cor.test( log10( matrix_samples$mutation_count / 3000 ), matrix_samples$ploidy, 
          method="spearman" )
#S = 123370, p-value = 0.004347
cor.test( matrix_sv_m$n_bnd , matrix_samples$ploidy, method="spearman" )  #S = 126150, p-value = 0.007342
cor.test( matrix_sv_m$n_dels , matrix_samples$ploidy, method="spearman" ) #S = 141370, p-value = 0.07725
cor.test( matrix_sv_m$n_ins , matrix_samples$ploidy, method="spearman" )  #S = 127710, p-value = 0.009707
cor.test( matrix_sv_m$n_inv , matrix_samples$ploidy, method="spearman" )  #S = 141310, p-value = 0.07665

# SV vs BRCA2, CDK12, TP53
#-------------------------------------------------------------------------------
n_dup = rep(NA, dim(matrix_somatic)[2])
dup = table( list_sv_m$sample_id[list_sv_m$svtype=="TANDEM"] )
m = match.idx( sample_ids, names(dup))
n_dup[m$idx.A] = as.numeric( dup[m$idx.B] )
rm(dup)
n_del = rep(NA, dim(matrix_somatic)[2])
del = table( list_sv_m$sample_id[list_sv_m$svtype=="DEL"] )
m = match.idx( sample_ids, names(del))
n_del[m$idx.A] = as.numeric( del[m$idx.B] )
rm(del)
n_bnd = rep(NA, dim(matrix_somatic)[2])
bnd = table( list_sv_m$sample_id[list_sv_m$svtype=="BND"] )
m = match.idx( sample_ids, names(bnd))
n_bnd[m$idx.A] = as.numeric( bnd[m$idx.B] )
rm(bnd)

n_inv = rep(NA, dim(matrix_somatic)[2])
inv = table( list_sv_m$sample_id[list_sv_m$svtype=="INV"] )
m = match.idx( sample_ids, names(inv))
n_inv[m$idx.A] = as.numeric( inv[m$idx.B] )
rm(inv)

p.inv = rep(NA, dim(matrix_somatic)[1])
p.del = rep(NA, dim(matrix_somatic)[1])
p.dup = rep(NA, dim(matrix_somatic)[1])
p.bnd = rep(NA, dim(matrix_somatic)[1])

MS = data.matrix( matrix_CNA_int_ploidy < LOSS_DOUBLE_NONSEX )
MS[is.na(MS)] = FALSE
MS = MS | matrix_somatic

for(i in 1:dim(MS)[1]){
    if( sum(MS[i,]) > 2 ){
        p.del[i] = wilcox.test( n_del ~ MS[i,] )$p.value
        p.dup[i] = wilcox.test( n_dup ~ MS[i,] )$p.value
        p.bnd[i] = wilcox.test( n_bnd ~ MS[i,] )$p.value
        p.inv[i] = wilcox.test( n_inv ~ MS[i,] )$p.value
    }
}
idx.remove = match.idx( dimnames(MS)[[1]], FLAGS)$idx.A
res = data.frame( p.del=signif(p.del,2), 
                  p.dup=signif(p.dup,2), 
                  p.bnd=signif(p.bnd,2), 
                  p.inv=signif(p.inv,2),
                  row.names=dimnames(MS)[[1]],
                  stringsAsFactors=FALSE)
res = res[ setdiff( 1: (dim(res)[1]), idx.remove ),]

association.del = res[ which( res$p.del<0.01),]
association.bnd = res[ which( res$p.bnd<0.01),]
association.dup = res[ which( res$p.dup<0.01),]
association.inv = res[ which( res$p.inv<0.01),]

# Association between BRCA2, CDK12, TP53 and SV types
#-------------------------------------------------------------------------------

BRCA2 = allele_effect("BRCA2")$alleles$bi
CDK12 = allele_effect("CDK12")$alleles$bi
TP53 = allele_effect("TP53")$alleles$bi

wilcox.test( matrix_samples$n_del_manta ~ BRCA2 ) #W = 2, p-value = 3.365e-06
wilcox.test( matrix_samples$n_dup_manta ~ CDK12 ) #W = 0, p-value = 0.003368
wilcox.test( matrix_samples$n_inv_manta ~ TP53 )  #W = 768.5, p-value = 0.0006627

# Association between TP53 and chromothripsis
#-------------------------------------------------------------------------------
fisher.test( table( matrix_samples$has_chromothripsis, TP53 ) )
# p-value = 9.524e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     2.433795 36.892488
# sample estimates:
#     odds ratio 
# 8.298174 

# translocations vs. chromothripsis
wilcox.test( matrix_sv_m$n_bnd ~ matrix_samples$has_chromothripsis)  # W = 744, p-value = 0.2168

# translocations vs. copy number alterations
cor.test( matrix_sv_m$n_bnd, matrix_samples$bases_cna_notref, method="spearman")
# S = 117520, p-value = 0.001308
#     rho 
# 0.3155405 

# translocations vs. tp53
wilcox.test( matrix_sv_m$n_bnd ~ TP53 )  # W = 946.5, p-value = 0.02833

# Association BRCA2 loss and mutation count
#-------------------------------------------------------------------------------
n_muts = log10( matrix_samples$mutation_count/3000 )
ae_brca2 = allele_effect( "BRCA2", do_plot=FALSE, axis_override=NULL )
no_BRCA2_no_chromo = !ae_brca2$alleles$bi & !matrix_samples$has_chromothripsis
no_BRCA2_chromo = !ae_brca2$alleles$bi & matrix_samples$has_chromothripsis
BRCA2_no_chromo = ae_brca2$alleles$bi & !matrix_samples$has_chromothripsis

wilcox.test( n_muts[ no_BRCA2_no_chromo ], 
             n_muts[ BRCA2_no_chromo ])
# data:  n_muts[no_BRCA2_no_chromo] and n_muts[BRCA2_no_chromo]
# W = 64, p-value = 0.0003864

wilcox.test( n_muts[ no_BRCA2_chromo ], 
             n_muts[ BRCA2_no_chromo ])
# data:  n_muts[no_BRCA2_chromo] and n_muts[BRCA2_no_chromo]
# W = 16, p-value = 0.0002016


# TP53 and inversion counts
#-------------------------------------------------------------------------------

wilcox.test( matrix_samples$n_inv_manta ~ TP53 ) # W = 768.5, p-value = 0.0006627
wilcox.test( matrix_samples$n_inv_manta[!matrix_samples$has_chromothripsis] ~ 
                 TP53[!matrix_samples$has_chromothripsis] ) # W = 529.5, p-value = 0.07655

median( matrix_samples$n_inv_manta[ !TP53 ] )  # [1] 59.5
median( matrix_samples$n_inv_manta[ TP53 ] )  # [1] 81
median( matrix_samples$n_inv_manta[ !matrix_samples$has_chromothripsis & !TP53 ] ) # [1] 56.5
median( matrix_samples$n_inv_manta[ !matrix_samples$has_chromothripsis &  TP53 ] ) # [1] 70


# BRCA2 and mutation counts
#-------------------------------------------------------------------------------

median( matrix_samples$mutation_count[ BRCA2 ]/3000 )  # [1] 6.959667
median( matrix_samples$mutation_count[ !BRCA2 ]/3000 ) # [1] 4.025667

# BRCA1/2 compound hets and mutation counts
#-------------------------------------------------------------------------------

ae_brca2 = allele_effect( "BRCA2", do_plot=FALSE, axis_override=NULL )
ae_brca1 = allele_effect( "BRCA1", do_plot=FALSE, axis_override=NULL )
only_bi = ae_brca2$alleles$bi | ae_brca1$alleles$bi
compound_het = !only_bi & ( ae_brca1$alleles$mono & ae_brca2$alleles$mono )

wilcox.test( matrix_mutsig_cosmic[only_bi,3], matrix_mutsig_cosmic[compound_het,3]   )  # W = 37, p-value = 0.1079
wilcox.test( matrix_mutsig_cosmic[!only_bi & !compound_het,3], matrix_mutsig_cosmic[compound_het,3]   )
# W = 162.5, p-value = 0.1248

wilcox.test( nmf[only_bi,8], nmf[compound_het,8]   )             # W = 42, p-value = 0.01998
wilcox.test( nmf[only_bi,8], nmf[!only_bi & !compound_het,8]   ) # W = 691, p-value = 4.433e-06

mc=matrix_samples$mutation_count
wilcox.test( log10( mc[!only_bi & !compound_het]/3000),log10( mc[!only_bi & compound_het]/3000))
# data:  log10(mc[!only_bi & !compound_het]/3000) and log10(mc[!only_bi & compound_het]/3000)
# W = 99, p-value = 0.01155



# test for association between gain at enhancer and second-line therapy with 
# abi or enza

clinical = read.table(
    '/notebook/human_sequence_prostate_WCDT/WCDT/metadata/second_line_therapy.txt',
    stringsAsFactors = FALSE, sep='\t', header=TRUE)

had_enza = clinical$enza == "Resistant"
had_abi = clinical$abi == "Resistant"
secondgen = rep("none", 101)
secondgen[had_enza & !had_abi] = "only_enza"
secondgen[!had_enza & had_abi] = "only_abi"
secondgen[had_enza & had_abi] = "both"
secondgen[had_enza | had_abi] = "either"
table( gain=curated_ar$EN_gain, secondgen )

table(  had_enza, gain=curated_ar$EN_gain)
table( had_abi, gain=curated_ar$EN_gain)
table( had_enza | had_abi, gain=curated_ar$EN_gain)


fisher.test( table( curated_ar$EN_gain, clinical$enza=="Resistant" & clinical$enza=="Resistant") )
fisher.test( table( curated_ar$EN_gain, clinical$abi=="Resistant") )
fisher.test( table( curated_ar$EN_gain, clinical$enza=="Resistant"|clinical$abi=="Resistant") )

fisher.test( table( curated_ar$AR_gain, clinical$enza=="Resistant") )
fisher.test( table( curated_ar$AR_gain, clinical$abi=="Resistant") )
fisher.test( table( curated_ar$AR_gain, clinical$enza=="Resistant"|clinical$abi=="Resistant") )


################################################################################
#                                                                              #
# MAIN TEXT FIGURES                                                            #
#                                                                              #
################################################################################


################################################################################
# FIGURE 1A: genome-wide plot of frequent SV, CNA
################################################################################

pdf( paste(dir_out, 'fig_sv_freq_genomewide.pdf', sep='/'), height=6, width=12)
layout(matrix(1,1,1))
par(mar=c(2,4,1,1))
chr_width=1
plot( -100, -100, xlim=c(1, 25), ylim=c(0, 250 ), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
for( i in 1:24 ){
    if( i %%2==0 ){
        col.alt="white"
    }else{
        col.alt="white"
    }
    rect( i, 0, i+chr_width, 1 + (chrom_lengths$V2[i] / 1000000), col=col.alt )   
    
}
axis(2, seq(from=0, to=250, by=20),las=1)
axis(1, at=(2:25)-0.5, labels = get.split.col(rownames(chrom_lengths),"chr",last=TRUE),las=1)

color_palatte=c("white", "#bdbdbd","black")
cmap = grDevices::colorRampPalette( colors=color_palatte )(50)
MIN_SV_BOUND=0.050
MAX_SV_BOUND=0.40
for(i in 1:24){
    sv_bins_colors = sv_bins_fig1a$sv[sv_bins_fig1a$chrom==rownames(chrom_lengths)[i]]
    cna_bins_amp = sv_bins_fig1a$amp[sv_bins_fig1a$chrom==rownames(chrom_lengths)[i]]
    cna_bins_del = sv_bins_fig1a$del[sv_bins_fig1a$chrom==rownames(chrom_lengths)[i]]
    sv_bins_colors[sv_bins_colors > MAX_SV_BOUND] = MAX_SV_BOUND-0.001
    colors_sv = color_scale( sv_bins_colors, cmap, color_bounds = c(MIN_SV_BOUND,MAX_SV_BOUND) )
    for(j in 1:length(sv_bins_colors)-1){
        rect( i+0.05, j, i + 0.95, j+1, col=colors_sv[j], border=NA)
    }  
    for(j in 1:(length(sv_bins_colors)-1)){
        if( cna_bins_del[j] < 0 ){
            rect( i+0.25, j, i + 0.75 , j+1, col="white", border=NA)
        }else{
            rect( i+0.5, j, i + 0.5 + cna_bins_amp[j] , j+1, col="#2ca25fee", border=NA)
            rect( i+0.5, j, i + 0.5 + cna_bins_del[j] , j+1, col="#0000ffee", border=NA)
        }
    }
    bins_centromere = floor( round( centromeres$chromStart[i] / 1000000, 0 ) ) :
        ceiling( round( centromeres$chromEnd[i] / 1000000, 0 ) ) 
    centro_start = centromeres$chromStart[i]
    centro_end = centromeres$chromEnd[i]
    
    top_cent = max(bins_centromere)
    bottom_cent = min(bins_centromere)
    polygon( c( i, i+1, mean(c(i, i+1))),
             c( top_cent, top_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
    polygon( c( i, i+1, mean(c(i, i+1))),
             c( bottom_cent, bottom_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
    lines( c(i,i), c(bottom_cent, top_cent), col="white")
    lines( c(i+1,i+1), c(bottom_cent, top_cent), col="white")
}
dev.off()

#sv_bins_fig1a$amp = sv_bins_fig1a$amp * -200
#sv_bins_fig1a$del = sv_bins_fig1a$del * 200
#sv_bins_fig1a$sv = sv_bins_fig1a$sv * 100

# Plot SV frequency, showing loci >3SD from mean
#-------------------------------------------------------------------------------

sv_bins_fig1a$sv[ which(sv_bins_fig1a$sv<0) ] = NA

mean( sv_bins_fig1a$sv, na.rm=TRUE)
sd( sv_bins_fig1a$sv, na.rm=TRUE)

sv_bins_fig1a$zscore = (sv_bins_fig1a$sv-mean(sv_bins_fig1a$sv,na.rm=TRUE))/sd(sv_bins_fig1a$sv,na.rm=TRUE)
sv_bins_fig1a$pval = signif( 1-pnorm( sv_bins_fig1a$sv, 
                                     mean=mean(sv_bins_fig1a$sv,na.rm=TRUE), 
                                     sd=sd(sv_bins_fig1a$sv,na.rm=TRUE)), 3)
sv_bins_fig1a$logp = -1*log10(sv_bins_fig1a$pval)
colors = rep("#0000ff33", dim(sv_bins_fig1a)[1])
colors_dark = rep("#0000ff", dim(sv_bins_fig1a)[1])

even_chroms = c(paste( "chr",seq(from=2, to=22, by=2), sep='' ), 'chrY')
for( i in 1:length(even_chroms)){
    colors[sv_bins_fig1a$chrom== even_chroms[i]] = "#ff000033"
    colors_dark[sv_bins_fig1a$chrom== even_chroms[i]] = "#ff0000"
}
sv_bins_fig1a$colors = colors
plot( sv_bins_fig1a$sv, pch=19, las=1, ylab="SV frequency", 
      xaxs="i", cex=0.75, col=sv_bins_fig1a$colors)
three_sd = mean(sv_bins_fig1a$sv, na.rm=TRUE) + (3*sd(sv_bins_fig1a$sv, na.rm=TRUE))
two_sd = mean(sv_bins_fig1a$sv, na.rm=TRUE) + (2*sd(sv_bins_fig1a$sv, na.rm=TRUE))
abline( three_sd, 0 )
text(10,25,"mean + 3sd", adj=0)
idx_high=which(sv_bins_fig1a$sv>=three_sd)
points( idx_high, sv_bins_fig1a$sv[idx_high], pch=1, cex=0.75, col=colors_dark[idx_high] )


################################################################################
# Genome wide sliding window analysis for peaks for TD, DELBND
################################################################################
# BIN_WIDTH   = 50000
# BIN_SPACING = 10000
if( FALSE ){
#for(cc in 1:24){
    chrom = rownames(chrom_lengths)[cc]
    chrom_start = 0
    chrom_end = round( chrom_lengths$V2[cc] )
    bin_span = chrom_end-chrom_start
    n_bins = round(bin_span/BIN_SPACING)
    sv_bins_td = rep(0, n_bins )
    sv_bins_delbndinv = rep(0, n_bins )
    
    bin_start = 0
    bin_end = bin_start + BIN_WIDTH
    bin_starts = rep(0, length(sv_bins_td))
    centro_start = centromeres$chromStart[cc]
    centro_end = centromeres$chromEnd[cc]
    for(j in 1:length(sv_bins_td)){
        if( j %% 100 == 0 ){
            print(paste(chrom, j, length(sv_bins_td)))
        }
        bin_starts[j] = bin_start
        if( (bin_start > centro_start & bin_start < centro_end ) |
            (bin_end > centro_start & bin_end < centro_end)){
            sv_bins_td[j] = NA
            sv_bins_delbndinv[j] = NA
        }else{
            idx_delbndinv = which( 
                   (list_sv_m$svtype=="DEL"| list_sv_m$svtype=="INV" | list_sv_m$svtype=="BND" )&
                    (
                       (list_sv_m$chrom_start==chrom & 
                         list_sv_m$pos_start>bin_start & list_sv_m$pos_start <= bin_end )|
                        (list_sv_m$chrom_end==chrom & 
                         list_sv_m$pos_end>bin_start & list_sv_m$pos_end <= bin_end )
                     )
                   )
            
            idx_td = list_sv_m$svtype=="TANDEM" & list_sv_m$chrom_start==chrom & 
                    (
                        ( list_sv_m$pos_start<bin_start & list_sv_m$pos_end > bin_start ) |
                        ( list_sv_m$pos_start<bin_end & list_sv_m$pos_end > bin_end ) |
                        ( list_sv_m$pos_start>bin_start & list_sv_m$pos_end < bin_end )
                    )

            sv_bins_delbndinv[j] = length( unique( list_sv_m$sample_id[idx_delbndinv] ) )
            sv_bins_td[j] = length( unique( list_sv_m$sample_id[idx_td] ) )
        }
        bin_start = bin_start + BIN_SPACING
        bin_end = bin_end + BIN_SPACING
    }
    
    if(cc==1){
    sv_bins_all = data.frame( chrom=rep(chrom,length(sv_bins_delbndinv)),
                              delbndinv=sv_bins_delbndinv,
                              td=sv_bins_td,
                              start=bin_starts, end=bin_starts+BIN_SPACING,
                              stringsAsFactors=FALSE)
    }else{
    sv_bins_all = rbind(sv_bins_all,
                        data.frame( chrom=rep(chrom,length(sv_bins_delbndinv)),
                                    delbndinv=sv_bins_delbndinv,
                                    td=sv_bins_td,
                                    start=bin_starts, end=bin_starts+BIN_SPACING,
                                    stringsAsFactors=FALSE))
    }
}
#write.table( sv_bins_all, '/notebook/human_sequence_prostate_WCDT/WCDT/results/genome_wide_SV_scan.txt',
#             sep='\t', quote=FALSE)

sv_bins_all = read.table( fn_genomewide_sliding_scan, 
                          header=TRUE, sep='\t' )
head(sv_bins_all)

pdf("/notebook/human_sequence_prostate_WCDT/WCDT/results/genome_wide_SV_scan.pdf",
    height=6, width=8)
for(cc in 1:24){
    chrom = rownames(chrom_lengths)[cc]
    layout(matrix(1:2,2,1))
    par(mar=c(3,3,0.5,0.5))
    ymax=max(sv_bins_all$td[sv_bins_all$chrom==chrom],na.rm=TRUE) + 1
    plot(sv_bins_all$td[sv_bins_all$chrom==chrom], pch=19, cex=0.5, 
         axes=FALSE, ylim=c(0,ymax) )
    axis(2, seq(from=0, to=ymax,by=5), las=1)
    idx_x=seq(from=0, to=sum( sv_bins_all$chrom==chrom ), by=1000)
    axis(1, at=idx_x, labels=c(0, round( sv_bins_all$start[idx_x] )/1000000))
    box()
    
    ymax=max(sv_bins_all$delbndinv[sv_bins_all$chrom==chrom],na.rm=TRUE) + 1
    plot(sv_bins_all$delbndinv[sv_bins_all$chrom==chrom], pch=19, cex=0.5, 
         axes=FALSE, ylim=c(0,ymax))
    axis(2, seq(from=0, to=ymax,by=5), las=1)
    idx_x=seq(from=0, to=sum( sv_bins_all$chrom==chrom ), by=1000)
    axis(1, at=idx_x, labels=c(0, round( sv_bins_all$start[idx_x] )/1000000))
    box()
}
dev.off()

s3=sv_bins_fig1a
s3 = s3[order(s3$sv, decreasing = TRUE),]
s3=s3[,c(1,5,6,2,3,4)]
s3$amp=s3$amp*-200
s3$del=s3$del*100

################################################################################
# FIGURE 1B: PTEN, CDH1  SV and CNA consequences on expression
################################################################################
sv_out = boxplot_inactivations( "PTEN", plot = FALSE, window_size=10000 )
pten.order = order(sv_out$xp)
sv_out = boxplot_inactivations( "TP53", plot = FALSE, window_size=10000 )
tp53.order = order(sv_out$xp)
sv_out = boxplot_inactivations( "RB1", plot = FALSE, window_size=200000 )
rb1.order = order(sv_out$xp)
sv_out = boxplot_inactivations( "CDKN1B", plot = FALSE, window_size=10000 )
cdkn1b.order = order(sv_out$xp)
sv_out = boxplot_inactivations( "CHD1", plot = FALSE, window_size=200000 )
chd1.order = order(sv_out$xp)

postscript( paste(dir_out, 'fig_SV_tumor_suppressor_RNA_barplots.ps', sep='/'), height=3.5, width=4.5)
layout( matrix(1:10,10,1), heights = rep(c(1,1), 5))
barplot_inactivations( "PTEN", window_size=1000, new.order=pten.order, show_events="break" )
barplot_inactivations( "TP53", window_size=1000, new.order=tp53.order, show_events="break" )
barplot_inactivations( "RB1", window_size=2000, new.order=rb1.order, show_events="break" )
barplot_inactivations( "CDKN1B", window_size=1000, new.order=cdkn1b.order, show_events="break" )
barplot_inactivations( "CHD1", window_size=2000, new.order=chd1.order, show_events="break" )
dev.off()

postscript( paste(dir_out, 'fig_SV_tumor_suppressor_RNA_allele_boxplots.ps', sep='/'), height=3.5, width=1)
layout(matrix(1:5,5,1))
par(mar=c(1,0.5,0.5,3))
ae=allele_effect("PTEN", do_plot=TRUE, axis_override=c(0,3,6))
ae=allele_effect("TP53", do_plot=TRUE, axis_override=c(0,1,2,3,4,5))
ae=allele_effect("RB1", do_plot=TRUE, axis_override=c(2,4,6))
ae=allele_effect("CDKN1B", do_plot=TRUE, axis_override=c(2,4,6))
ae=allele_effect("CHD1", do_plot=TRUE, axis_override=c(0,3,6))
dev.off()


################################################################################
# FIGURE 2 FOXA1 TANDEM DUPLICATION
################################################################################

# Calulate counts_cna and counts_sv
#-------------------------------------------------------------------------------
chrom="chr14"
n_counts=ceiling( chrom_lengths$V2[rownames(chrom_lengths)==chrom]/10000)
counts_cna = rep(0,n_counts)
counts_sv = rep(0,n_counts)
lup_means = rep(0, n_counts)
has_amp_foxa1 = rep(FALSE, 101)
has_td_foxa1 = rep(FALSE,101)
for(i in 1:101){
    sample_id = sample_ids[i]
    local_counts = rep(0, n_counts)
    local_counts_sv = rep(0, n_counts)
    bed = get_bed( sample_id, dir_build, allow_small = TRUE )
    bed=bed[bed$chrom==chrom,]
    for(j in 1:dim(bed)[1]){
        if( bed$cn.int[j]>=GAIN_NONSEX ){
            x1 = round( bed$start[j]/10000, 0 )
            x2 = round( bed$stop[j]/10000, 0 )
            local_counts[x1:x2] = local_counts[x1:x2] + 1
        }
    }
    if( length(local_counts)>length(counts_cna)){
        local_counts = local_counts[1:length(counts_cna)]   
    }
    local_counts[local_counts>1]=1
    counts_cna = counts_cna + local_counts
    
    sv = list_sv_m[list_sv_m$chrom_start==chrom & 
                       list_sv_m$sample_id==sample_id & 
                       list_sv_m$svtype=="TANDEM",]
    if( dim(sv)[1] ){
        for(j in 1:dim(sv)[1]){
            x1 = round( sv$pos_start[j] / 10000, 0 )
            x2 = round( sv$pos_end[j] / 10000, 0 )
            local_counts_sv[x1:x2] = local_counts_sv[x1:x2] + 1
        }
    }
    has_td_foxa1[i] = local_counts_sv[3754]>0
    has_amp_foxa1[i] = local_counts[3759]>0
    
    local_counts_sv[local_counts_sv>1]=1
    counts_sv = counts_sv + local_counts_sv
}
which(has_td_foxa1)
which(allele_effect("FOXA1")$alleles$nonsense)

# Plot all individual tandem duplications in a stack sorted by pos_start
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_FOXA1_TD_stack_closeup.pdf", 
     height=1.5, width=3.5)
par(mar=c(0,4,1,1))
x1 = 3650
x2 = 3800
sv_fox=list_sv_m[ (list_sv_m$svtype=="TANDEM" ) & list_sv_m$chrom_start=="chr14" &
                      list_sv_m$pos_start < 37490000 & list_sv_m$pos_end > 37490000, ]
sv_fox = sv_fox[order(sv_fox$pos_start),]
layout(matrix(1,1,1))
plot(-1,-1, ylim=c(0,5), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
cur_y=1
for(i in 1:dim(sv_fox)[1]){
    lines( c(sv_fox$pos_start[i]/10000, sv_fox$pos_end[i]/10000), c(cur_y,cur_y), lwd=3, col="cornflowerblue" )
    cur_y=cur_y+0.25
}
box()
show_gene( "MIPOL1", height=0.5, height_text = -1000) # show locus but not name
show_gene( "FOXA1", height=0.5, height_text = -1000)
dev.off()

# Plot ChIP-seq hotspots
#-------------------------------------------------------------------------------

# load lupien data and bin into 10000 nt
lup_fox = read.table( '/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_04_19_Lupien_FOXA1.bg',
                      sep='\t', stringsAsFactors=FALSE)
names(lup_fox) = c("chrom", "start" , "stop", "score")

for( i in x1:x2){
    bound_start = i * 10000
    bound_end = (i+1) * 10000
    lup_means[i] = mean( lup_fox$score[
        (lup_fox$stop > bound_start & lup_fox$start < bound_end) |
            (lup_fox$start < bound_end & lup_fox$stop > bound_start) ])
}
lup_means[is.nan(lup_means)] = 0
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_FOXA1_lupien.pdf", 
     height=1.5, width=3.5)
layout(matrix(1,1,1))
par(mar=c(3,4,0.5,1))

plot(-1,-1, ylim=c(0,15), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in x1:x2){
    lines( c(i,i), c(0, lup_means[i] ),lwd=3, col="black" )
}
axis(2, seq(from=0, to=15,by=5), las=1, cex=0.75 )
box()
axis(1, at=seq(from=x1, to=x2, by=50),
     labels =seq(from=x1/100, to=x2/100, by=0.5), las=1 )
dev.off()

# plot  DNA amplifications and tandem duplications
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_FOXA1_TD_closeup.pdf", 
     height=3, width=3.5)
layout(matrix(1:2,2,1), heights = c(1,3) )
par(mar=c(0,4,1,1))
x1 = 3650
x2 = 3800
plot(-1,-1, ylim=c(0,40), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=40,by=5) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=40,by=20), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, counts_cna[i] ),lwd=6, col="#00000033" )
}
box()

par(mar=c(4,4,1,1))
plot(-1,-1, ylim=c(0,20), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=20,by=5) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=20,by=10), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, counts_sv[i] ),lwd=6, col="#0000cc33" )
}
box()
dev.off()


# plot whole chromosome
#-------------------------------------------------------------------------------
x1 = 10
x2 = n_counts
par(mar=c(0,4,2,1))
plot(-1,-1, ylim=c(0,40), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="", main=chrom)
for(i in seq(from=10, to=80,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=40,by=10), las=1, cex=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, counts_cna[i] ),lwd=2, col="#00000033" )
}
text( 35, 35, "amp frequency", adj=0, font=2)
box()
par(mar=c(4,4,1,1))
plot(-1,-1, ylim=c(0,20), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=70,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=20,by=10), las=1, cex=0.75 )
axis(1, at=seq(from=x1, to=x2, by=50),
     labels =seq(from=x1/100, to=x2/100, by=0.5), las=1 )
text( 50, 15, "TD frequency", adj=0, font=2)
for(i in x1:x2){
    lines( c(i,i), c(0, counts_sv[i]),lwd=2, col="#0000cc33" )
}
axis(2, seq(from=0, to=90,by=10), las=1, cex=0.75 )
axis(1, at=seq(from=x1, to=x2, by=50),
     labels =seq(from=x1/100, to=x2/100, by=0.5), las=1 )
box()


################################################################################
# MYC TANDEM DUPLICATION
################################################################################

# Calulate chr_counts and sv_counts
#-------------------------------------------------------------------------------
chrom="chr8"
n_counts=ceiling( chrom_lengths$V2[rownames(chrom_lengths)==chrom]/10000)
chr_counts = rep(0,n_counts)
sv_counts = rep(0, n_counts)
has_amp_myc = rep(FALSE, 101)
has_amp_pcat = rep(FALSE,101)
has_td_pcat = rep(FALSE,101)
has_td_myc = rep(FALSE,101)
has_td_pou = rep(FALSE,101)
has_td_myc = rep(FALSE,101)
has_amp_pou = rep(FALSE,101)
for(i in 1:101){
    sample_id = sample_ids[i]
    local_counts = rep(0, n_counts)
    local_counts_sv = rep(0, n_counts)
    bed = get_bed( sample_id, dir_build, allow_small = TRUE )
    bed=bed[bed$chrom==chrom,]
    for(j in 1:dim(bed)[1]){
        if( bed$cn.int[j]>=GAIN_NONSEX ){
            x1 = round( bed$start[j]/10000, 0 )
            x2 = round( bed$stop[j]/10000, 0 )
            local_counts[x1:x2] = local_counts[x1:x2] + 1
        }
    }
    if( length(local_counts)>length(chr_counts)){
        local_counts = local_counts[1:length(chr_counts)]   
    }
    local_counts[local_counts>1]=1
    chr_counts = chr_counts + local_counts
    
    sv = list_sv_m[list_sv_m$chrom_start==chrom & 
                       list_sv_m$sample_id==sample_id & 
                       list_sv_m$svtype=="TANDEM",]
    if( dim(sv)[1] ){
        for(j in 1:dim(sv)[1]){
            x1 = round( sv$pos_start[j] / 10000, 0 )
            x2 = round( sv$pos_end[j] / 10000, 0 )
            local_counts_sv[x1:x2] = local_counts_sv[x1:x2] + 1
        }
    }
    
    has_td_pcat[i] = local_counts_sv[12701]>0
    has_amp_pcat[i] = local_counts[12701]>0
    has_td_myc[i] = local_counts_sv[12773]>0
    has_amp_myc[i] = local_counts[12773]>0
    has_td_pou[i] = local_counts_sv[12741]>0
    has_amp_pou[i] = local_counts[12741]>0
    has_td_myc[i] = local_counts_sv[12773]>0
    
    local_counts_sv[local_counts_sv>1]=1
    sv_counts = sv_counts + local_counts_sv
}

sum(has_amp_myc)
sum(has_td_myc)
sum(has_amp_pcat)
sum(has_td_pcat)
sum(has_amp_pou)
sum(has_td_pou)

MYC = data.frame(
    has_amp_myc,
    has_td_myc,
    has_amp_pcat,
    has_td_pcat,
    has_amp_pou,
    has_td_pou,
    has_nothing=!has_amp_myc & !has_td_myc & !has_amp_pcat & !has_td_pcat,
    xp=as.numeric(matrix_tpm["MYC",]),
    logxp=log( as.numeric(matrix_tpm["MYC",]) ),
    xp.pou=as.numeric(matrix_tpm["POU5F1B",]),
    logxp.pou=log( as.numeric(matrix_tpm["POU5F1B",]) ),
    stringsAsFactors = FALSE)
boxplot( 
    MYC$logxp[MYC$has_nothing],
    MYC$logxp[MYC$has_amp_myc],
    MYC$logxp[MYC$has_td_myc],
    MYC$logxp[MYC$has_td_pou])

sample_ids[ which(MYC$has_td_pou) ]

cor.test( as.numeric(matrix_tpm["MYC",]), as.numeric(matrix_CNA_int_ploidy["MYC",]), method="spearman")

wilcox.test( MYC$logxp[MYC$has_nothing], MYC$logxp[MYC$has_td_myc] )
wilcox.test( MYC$logxp[MYC$has_nothing], MYC$logxp[MYC$has_td_pou] )
median(MYC$xp[MYC$has_nothing])
median( MYC$xp[MYC$has_td_pou], na.rm=TRUE )

# Plot all individual tandem duplications in a stack sorted by pos_start
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_MYC_TD_stack_closeup.pdf", 
     height=1.5, width=3.5)
par(mar=c(0,4,1,1))
x1 = 12650
x2 = 12800

midpoint1 = 127416611
midpoint2 = 127000000
sv_myc=list_sv_m[ list_sv_m$svtype=="TANDEM" & 
                      list_sv_m$chrom_start=="chr8" &
                      (
                          (list_sv_m$pos_start < midpoint1 & list_sv_m$pos_end > midpoint1) |
                              (list_sv_m$pos_start < midpoint2 & list_sv_m$pos_end > midpoint2) 
                      ), ]
sv_myc = sv_myc[order(sv_myc$pos_start),]
layout(matrix(1,1,1))
plot(-1,-1, ylim=c(0,40), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
cur_y=5
for(i in 1:dim(sv_myc)[1]){
    lines( c(sv_myc$pos_start[i]/10000, sv_myc$pos_end[i]/10000), c(cur_y,cur_y), 
           col="cornflowerblue", lwd=2 )
    cur_y=cur_y+0.65
}
box()

lines( c( sv_myc$pos_start[12]/10000, sv_myc$pos_start[12]/10000),
       c(0, 100), col="#33333366")
lines( c( sv_myc$pos_end[10]/10000, sv_myc$pos_end[10]/10000),
       c(0, 100), col="#33333366")
show_gene( "PRNCR1", height=1, height_text = -1000)
show_gene( "MYC", height=1, height_text = -1000)
show_gene( "PCAT1", height=1, height_text = -1000)
show_gene( "POU5F1B", height=1, height_text = -1000)
show_gene("CCAT2", height=3, height_text=-1000)
show_gene("CASC8", height=3, height_text=-1000)
lines( c(127072694/10000, 127082221/10000), c(3,3), lwd=6) # PCAT2
dev.off()


# plot whole chromosome pileup of DNA amplifications and tandem duplications
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_MYC_TD.pdf", height=4, width=3)
layout(matrix(1:2,2,1), heights = c(0.5,1.5) )
par(mar=c(0,4,2,1))

x1 = 1
x2 = n_counts
plot(-1,-1, ylim=c(0,80), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=80,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=80,by=40), las=1, cex=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, chr_counts[i] ),lwd=1, col="#999999" )
}
rect( centromeres[8,2]/10000, 0, centromeres[8,3]/10000, 65, col="white", border="white" )
lines( c(60,60), c(0, chr_counts[i] ),lwd=1, col="#999999" )

box()
par(mar=c(4,4,1,1))
plot(-1,-1, ylim=c(0,30), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=70,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=30,by=10), las=1, cex=0.75 )
axis(1, at=seq(from=0, to=14000, by=7000),
     labels = seq(from=0, to=140, by=70), las=1 )
for(i in x1:x2){
    lines( c(i,i), c(0, sv_counts[i]),lwd=1, col="#999999" )
}
rect( centromeres[8,2]/10000, 0, centromeres[8,3]/10000, 15, col="white", border="white" )
box()
dev.off()


# plot MYC region pileup of DNA amplifications and tandem duplications
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_MYC_TD_closeup.pdf", 
     height=3, width=3.5)
layout(matrix(1:2,2,1), heights = c(1,3) )
par(mar=c(0,4,1,1))
x1 = 12650
x2 = 12800
plot(-1,-1, ylim=c(0,100), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=100,by=20) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=100,by=20), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, chr_counts[i] ),lwd=4, col="#00000033" )
}
box()
par(mar=c(4,4,1,1))
plot(-1,-1, ylim=c(0,30), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=30,by=5) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=30,by=10), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, sv_counts[i] ),lwd=4, col="#0000cc33" )
}
box()
dev.off()

# Plot ChIP-seq hotspots
#-------------------------------------------------------------------------------

# load lupien data and bin into 10000 nt
lup_myc = read.table( '/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_04_19_Lupien_MYC.bg',
                      sep='\t', stringsAsFactors=FALSE)
names(lup_myc) = c("chrom", "start" , "stop", "score")

x1=12650
x2=12800
for( i in x1:x2){
    bound_start = i * 10000
    bound_end = (i+1) * 10000
    lup_means[i] = mean( lup_myc$score[
        (lup_myc$stop > bound_start & lup_myc$start < bound_end) |
            (lup_myc$start < bound_end & lup_myc$stop > bound_start) ])
}
lup_means[is.nan(lup_means)] = 0
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_MYC_lupien.pdf", 
     height=1.5, width=3.5)
layout(matrix(1,1,1))
par(mar=c(3,4,0.5,1))

plot(-1,-1, ylim=c(0,15), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in x1:x2){
    lines( c(i,i), c(0, lup_means[i] ),lwd=1, col="black" )
}
axis(2, seq(from=0, to=15,by=5), las=1, cex=0.75 )
box()
axis(1, at=seq(from=x1, to=x2, by=50),
     labels =seq(from=x1/100, to=x2/100, by=0.5), las=1 )
dev.off()



pcat1 = read.table( '~/Documents/notebook/human_sequence_prostate_WCDT/WCDT/results/PCAT1.abundances.txt', sep='\t', header=TRUE, check.names = FALSE)
m=match.idx( sample_ids, dimnames(pcat1)[[2]] )
plot( log(1+as.numeric( matrix_tpm["PCAT1",m$idx.A] )),
      log(1+as.numeric(pcat1[,m$idx.B])))

plot_cna_expr = function(symbol){
    plot( as.numeric(matrix_tpm[symbol,]),
          as.numeric(matrix_CNA_int_ploidy[symbol,])
    )
}




################################################################################
# FIGURE 3A: SV frequency across samples
################################################################################

fn=paste(dir_fig, 'fig_relative_absolute_SV.ps', sep='/')
postscript(fn, height=2, width=8 )
layout(matrix(1,1,1))
par(mar=c(0.5,3,0.5,1))
data = matrix(0,ncol=length(sample_ids),nrow=5)
dimnames(data)[[2]] = sample_ids
x =table(list_sv_m$sample_id[list_sv_m$svtype=="DEL"])
m = match.idx( sample_ids, names(x) )
data[1,m$idx.A] = as.numeric(x)[m$idx.B]
x =table(list_sv_m$sample_id[list_sv_m$svtype=="BND"])
m = match.idx( sample_ids, names(x) )
data[2,m$idx.A] = as.numeric(x)[m$idx.B]
x =table(list_sv_m$sample_id[list_sv_m$svtype=="INV"])
m = match.idx( sample_ids, names(x) )
data[3,m$idx.A] = as.numeric(x)[m$idx.B]
x =table(list_sv_m$sample_id[list_sv_m$svtype=="TANDEM"])
m = match.idx( sample_ids, names(x) )
data[4,m$idx.A] = as.numeric(x)[m$idx.B]
x =table(list_sv_m$sample_id[list_sv_m$svtype=="INS"])
m = match.idx( sample_ids, names(x) )
data[5,m$idx.A] = as.numeric(x)[m$idx.B]

new.order = order(data[1,])
colors = c("#2c7bb6","#fdae61","#abd9e9","#d7191c","#ffffbf")
barplot( data[,new.order], col=colors, 
         border="black", las=2, cex.names =0.5, 
         names.arg = rep("", dim(data)[2]),
         axes=FALSE, ylim=c(0, 1200), xaxs="i")
axis(2, at=c(0, 600, 1200), labels = c(0,600,1200), cex.axis=1, las=1)
legend(0, 1200, c("deletions","translocations","inversions","tandem duplications","insertions"), 
       fill=colors,
       bty="n", cex=1)
box()

dev.off()

fn=paste(dir_fig, 'fig_brca2_cdk12_tp53_SV.ps', sep='/')
postscript(fn, height=1, width=8 )
col_sv = matrix("white", nrow=4, ncol=101)
col_sv[1,INACTIVE["BRCA2",]] = "black"
col_sv[2,matrix_samples$has_chromothripsis] = "black"
col_sv[3,INACTIVE["CDK12",]] = "black"
col_sv[4,INACTIVE["TP53",]] = "black"
col_sv[1,which(sample_ids=="DTB-165-PRO")]="white"
dimnames(col_sv)[[1]] = c("","","","")
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot.color.grid(col_sv[,new.order],block.height = 4,block.width=2, space.X = 1, 
                space.Y=1,show.axis.Y = FALSE, show.axis.X=FALSE, cex.y=1)
dev.off()


################################################################################
# FIGURE 3B: EXAMPLE CIRCOS PLOTS
################################################################################
par(mar=c(4,4,4,4))
fn_cyto='/notebook/human_sequence_prostate_WCDT/WCDT/metadata/HG38_cytoband_std_chrom.txt'
cyto_info = read.table(fn_cyto, stringsAsFactors=FALSE,header=TRUE, sep='\t')

sample_id="DTB-097-PRO"
layout(matrix(1,1,1))
plot_circos( cyto_info=cyto_info, 
             cna=get_bed( sample_id, dir_build ), 
             sv=get_sv_by_sample( sample_id ), 
             main="", 
             gene_info = NULL,
             fn_out = paste( fn=paste('/notebook/human_sequence_prostate_WCDT/WCDT/figures/circos_', sample_id, '.pdf', sep='' ) ) )
sample_id="DTB-104-BL"
plot_circos( cyto_info=cyto_info, 
             cna=get_bed( sample_id, dir_build ), 
             sv=get_sv_by_sample( sample_id ), 
             main="", 
             gene_info = NULL,
             fn_out = paste( fn=paste('/notebook/human_sequence_prostate_WCDT/WCDT/figures/circos_', sample_id, '.pdf', sep='' ) ) )
sample_id="DTB-063-BL"
plot_circos( cyto_info=cyto_info, 
             cna=get_bed( sample_id, dir_build ), 
             sv=get_sv_by_sample( sample_id ), 
             main="", 
             gene_info = NULL,
             fn_out = paste( fn=paste('/notebook/human_sequence_prostate_WCDT/WCDT/figures/circos_', sample_id, '.pdf', sep='' ) ) )


#for(i in 1:length(sample_ids)){
#    sample_id=sample_ids[i]
#    print(i)
#    plot_circos( cyto_info=cyto_info, 
#                 cna=get_bed( sample_id, dir_build ), 
#                 sv=get_sv_by_sample( sample_id ), 
#                 main="", 
#                 gene_info = NULL,
#                 fn_out = paste( fn=paste('/notebook/human_sequence_prostate_WCDT/WCDT/figures/circos_', sample_id, '.pdf', sep='' ) ) )
#}

################################################################################
# FIGURE 3C: SV boxplots of association between SV in CDK12, BRCA2 and TP53 
################################################################################

fn=paste(dir_fig, 'fig_sv_CDK12_BRCA2_TP53_combined.pdf', sep='/')
pdf(fn, height=4, width=4 )

layout(matrix(1:3,3,1))
par(mar=c(1,3,0.5,1))
colors = c("#2c7bb6","#2c7bb6","#d7191c","#d7191c","#abd9e9","#abd9e9")
pch_style=1
pch_color="#333333ee"
jit_width=1
BRCA2 = allele_effect("BRCA2")$alleles$bi
CDK12 = allele_effect("CDK12")$alleles$bi
TP53 = allele_effect("TP53")$alleles$bi

bp=boxplot( 
    matrix_samples$n_del_manta[ !BRCA2 ],
    matrix_samples$n_del_manta[ BRCA2 ],
    matrix_samples$n_del_manta[ !CDK12 ],
    matrix_samples$n_del_manta[ CDK12 ],
    matrix_samples$n_del_manta[ !TP53 ],
    matrix_samples$n_del_manta[ TP53 ],
    lwd=1, borders=colors,
    boxwex=0.5, outline=FALSE, ylim=c(0,400), axes=FALSE, col=colors, xaxs="i")
axis(2, c(0,200,400), las=1)
box()
par(mar=c(1,3,1,1))
bp=boxplot( 
    matrix_samples$n_dup_manta[ !BRCA2 ],
    matrix_samples$n_dup_manta[ BRCA2 ],
    matrix_samples$n_dup_manta[ !CDK12 ],
    matrix_samples$n_dup_manta[ CDK12 ],
    matrix_samples$n_dup_manta[ !TP53 ],
    matrix_samples$n_dup_manta[ TP53 ],
    lwd=1, borders=colors,
    boxwex=0.5, outline=FALSE, ylim=c(0,500), axes=FALSE, col=colors, xaxs="i")
axis(2, c(0,250,500), las=1)
box()

par(mar=c(1,3,1,1))
bp=boxplot( 
    matrix_samples$n_inv_manta[ !BRCA2 ],
    matrix_samples$n_inv_manta[ BRCA2 ],
    matrix_samples$n_inv_manta[ !CDK12 ],
    matrix_samples$n_inv_manta[ CDK12 ],
    matrix_samples$n_inv_manta[ !TP53 ],
    matrix_samples$n_inv_manta[ TP53 ],
    lwd=1, borders=colors,
    boxwex=0.5, outline=FALSE, ylim=c(0,200), axes=FALSE, col=colors, xaxs="i")
axis(2, c(0,100,200), las=1)
box()
dev.off()



################################################################################
# FIGURE 3D: chromothripsy and TP53
################################################################################

postscript("/notebook/human_sequence_prostate_WCDT/WCDT/figures/TP53_inversion_vs_deletions.ps", height=4, width=4)
layout(matrix(1,1,1 ))
par(mar=c(4,4,1,1))
colors = rep("black", 101)
colors[BRCA2] = "blue"
colors[matrix_samples$has_chromothripsis] = "darkorange"
plot(matrix_sv_m$n_dels, matrix_sv_m$n_inv, pch=19, col=colors, xlab="", ylab="", 
     las=1, xlim=c(0,500), ylim=c(0,500), xaxs="i", yaxs="i") 
legend(0, 500, c("BRCA2 alteration", "chromothripsis"), fill=c("blue", "darkorange"),  bty="n", cex=1)
dev.off()

################################################################################
# FIGURE 3E mutation count boxplot, chromothripsis & BRCA2
################################################################################
n_muts = log10( matrix_samples$mutation_count/3000 )
ae_brca2 = allele_effect( "BRCA2", do_plot=FALSE, axis_override=NULL )
no_BRCA2_no_chromo = !ae_brca2$alleles$bi & !matrix_samples$has_chromothripsis
no_BRCA2_chromo = !ae_brca2$alleles$bi & matrix_samples$has_chromothripsis
BRCA2_no_chromo = ae_brca2$alleles$bi & !matrix_samples$has_chromothripsis

postscript( paste(dir_out, 'fig_mutation_frequency_boxplot.ps', sep='/'), height=3, width=4)
layout(matrix(1,1,1))
par(mar=c(0.5,0.5,1,1))
bp=boxplot( n_muts[ no_BRCA2_no_chromo ], 
            n_muts[no_BRCA2_chromo ], 
            n_muts[BRCA2_no_chromo ], 
            names=c("","",""), axes=FALSE, las=1, ylim=c(0, 2.1),
            box.wex=0.5, col=c("white", "orange", "cornflowerblue"), yaxs="i")
axis( 2, at=c(0, 0.39794, 0.69897, 1, 1.5, 2), 
      labels= c(0, round(10^0.39794,1), 5, 10, round(10^1.5,0), 10^2 ) )
box()
dev.off()


################################################################################
# AR tandem duplication
################################################################################

# AR enhancer copy number
fn_curated='/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_05_07_curated_AR_SV.txt'
fn_cna_en='/notebook/human_sequence_prostate_WCDT/WCDT/results/build_2018_04_15/2018_04_15_matrix_AR_enhancer_CN.txt'
ARCN = read.table(fn_cna_en,
           header=TRUE, stringsAsFactors=FALSE, sep='\t')
curated_ar = read.table(fn_curated, header=TRUE, stringsAsFactors=FALSE, 
                        sep='\t', row.names=1)
curated_ar$logxp = log(1+curated_ar$tpm)
curated_ar$AR_copies = as.numeric(matrix_CNA_int_ploidy["AR",])
curated_ar$EN_copies = ARCN$copies
curated_ar$AR_gain = curated_ar$AR_copies >= GAIN_SEX
curated_ar$EN_gain = curated_ar$EN_copies >= GAIN_SEX

curated_ar$AR_td = rep(FALSE, 101)
curated_ar$EN_td = rep(FALSE, 101)
curated_ar$EN_td[curated_ar$code=="EATD" | 
                 curated_ar$code=="EATD-EAG" |     
                 curated_ar$code=="ETD" | 
                     curated_ar$code=="ETD-AG" |
                     curated_ar$code=="ETD-EAG" | 
                     curated_ar$code=="ETD-EATD" | 
                     curated_ar$code=="ETD-EG-EAG"] = TRUE

curated_ar$AR_td[ curated_ar$code=="ATD"|      
                  curated_ar$code=="EATD" |    
                  curated_ar$code=="EATD-EAG"|
                  curated_ar$code=="ETD-EATD"] = TRUE
curated_ar = curated_ar[order( curated_ar$AR_copies ),]
postscript( '/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_AR_enhancer_cna.ps',
     width=3.5, height=1.5)
layout(matrix(1,1,1))
par(mar=c(0.5, 4,0.5, 0.5))
AE00 = !curated_ar$AR_gain & !curated_ar$EN_gain
AE10 = curated_ar$AR_gain & !curated_ar$EN_gain
AE01 = !curated_ar$AR_gain & curated_ar$EN_gain
AE11 = curated_ar$AR_gain & curated_ar$EN_gain
boxplot( curated_ar$logxp[ AE00 ],
         curated_ar$logxp[ AE10 ],
         curated_ar$logxp[ AE01 ],
         curated_ar$logxp[ AE11 ], 
         las=1, ylim=c(0,12),
         col="lightgrey", lwd=0.5, cex=0.5)
dev.off()




postscript( '/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_AR_enhancer_barplot.ps',
            width=3.5, height=1.5)
colors = rep("grey", 101)
colors[curated_ar$EN_td & !curated_ar$AR_td ] = "#d7191c" #red
colors[!curated_ar$EN_td & curated_ar$AR_td] = "#fdae61" # orange
colors[curated_ar$EN_td & curated_ar$AR_td] = "#2b83ba" #blue
layout(matrix(1:2,2,1), heights=c(1,2))
par(mar=c(1, 3, 0.5, 0.5))
cut_off_copies = curated_ar$AR_copies
cut_off_copies[cut_off_copies<7]=0
barplot( cut_off_copies, xaxs="i", col=colors, border=NA, las=1, ylim=c(8,70), axes=FALSE, width=2 )
axis(2, seq(from=10, to=70, by=30), las=1, cex.axis=0.5)
par(mar=c(0.5, 3, 0, 0.5))
cut_off_copies = curated_ar$AR_copies
cut_off_copies[cut_off_copies>6] = 6
barplot( cut_off_copies, xaxs="i", col=colors, border=NA, las=1, ylim=c(0,6), axes=FALSE, width=2 )
axis(2, seq(from=0, to=8, by=2), las=1, cex.axis=0.75)
dev.off()

sum(!curated_ar$EN_td & !curated_ar$AR_td )
sum(curated_ar$EN_td & !curated_ar$AR_td )
sum(!curated_ar$EN_td & curated_ar$AR_td )
sum(curated_ar$EN_td & curated_ar$AR_td )

sum(!curated_ar$EN_gain & curated_ar$EN_td  ) # TD without gain doesn't happen
sum(!curated_ar$EN_gain & !curated_ar$AR_gain ) # no gain anywhere = 18
sum(curated_ar$EN_gain & !curated_ar$AR_gain )

sum(curated_ar$EN_td & !curated_ar$AR_td & !curated_ar$AR_gain )
sum(curated_ar$EN_gain & !curated_ar$AR_gain )


colors=matrix("white", nrow=2, ncol=101)
independent_TD = curated_ar$code=="ETD-EATD"
en_gain_no_TD = curated_ar$code=="EG"
colors[1,independent_TD] = "black"
colors[2,en_gain_no_TD] = "black"

postscript( '/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_AR_enhancer_subannotations.ps',
            width=3.5, height=1.5)
layout(matrix(1,1,1))
par(mar=c(1, 3, 0.5, 0.5))
plot.color.grid( colors, block.height = 5, block.width = 5, 
                 space.X = 0, space.Y=0, show.axis.X = FALSE, show.axis.Y = FALSE,
                 border=FALSE)
dev.off()

AE00 = !curated_ar$AR_gain & !curated_ar$EN_gain
AE10 = curated_ar$AR_gain & !curated_ar$EN_gain
AE01 = !curated_ar$AR_gain & curated_ar$EN_gain
AE11 = curated_ar$AR_gain & curated_ar$EN_gain

layout(matrix(1:4,1,4),heights=4,1)
par(mar=c(3,3,0.5,0.5))
plot( curated_ar$AR_copies[AE00],curated_ar$EN_copies[AE00], xlim=c(0,2), ylim=c(0,2))
abline(0,1)
plot( curated_ar$AR_copies[AE01],curated_ar$EN_copies[AE01], xlim=c(0,3), ylim=c(0,3))
abline(0,1)
plot( curated_ar$AR_copies[AE10],curated_ar$EN_copies[AE10], xlim=c(0,4), ylim=c(0,4))
abline(0,1)
plot( curated_ar$AR_copies[AE11],curated_ar$EN_copies[AE11], xlim=c(0,70), ylim=c(0,70))
abline(0,1)

colors=rep("lightgrey", 101)
colors[AE10] = "yellow"
colors[AE01] = "orange"
colors[AE11] = "red"
layout(matrix(1:2,2,1),heights=4,1)
barplot( curated_ar$AR_copies, xaxs="i", col=colors, border="darkgrey", las=1,ylim=c(0,10) )
box()
colors = matrix("white", nrow=2, ncol=101)
colors[1, curated_ar$EN_td ] = "black"
colors[2, curated_ar$AR_td ] = "black"
plot.color.grid( colors, block.height = 5, block.width = 5, 
                 space.X = 0, space.Y=0, show.axis.X = FALSE)

colors = matrix("white", nrow=4, ncol=101)

colors[1, !curated_ar$AR_gain & !curated_ar$EN_gain] = "black"
colors[2, curated_ar$AR_gain & !curated_ar$EN_gain] = "black"
colors[3, !curated_ar$AR_gain & curated_ar$EN_gain] = "black"
colors[4, curated_ar$AR_gain & curated_ar$EN_gain] = "black"
plot.color.grid( colors, block.height = 5, block.width = 5, 
                 space.X = 0, space.Y=0, show.axis.X = FALSE)



sum(!curated_ar$AR_td & curated_ar$EN_td)

colors=rep("black",101)
colors[curated_ar$code2=="0TD"] = "blue"
layout(matrix(1:2,1,2))
par(mar=c(5,5,4,1))
cna_ar = as.numeric(matrix_CNA_int_ploidy["AR",])
cna_en = ARCN$copies
xp_ar = log(1+as.numeric(matrix_tpm["AR",]))
plot( cna_ar, cna_en, las=1,
      xlab="CNA AR gene", ylab="CNA AR enhancer", main="All", col=colors)
abline(0,1)
legend( 0, 70, c("no TD", "TD"), fill=c("blue", "black"), bty = "n")

plot( cna_ar, cna_en, xlim=c(0,7), ylim=c(0,7), las=1, col=colors,
      xlab="CNA AR gene", ylab="CNA AR enhancer", main="Zoomed in")
legend( 0, 7, c("no TD", "TD"), fill=c("blue", "black"), bty = "n")
abline(0,1)
sum( cna_en>cna_ar)

has_TD = curated_ar$code2!="0TD"
anova( lm( curated_ar$tpm ~ curated_ar$AR_gain + curated_ar$EN_gain) )

has_nothing = curated_ar$cn<GAIN_SEX & cna_en<GAIN_SEX
has_AR_amp = cna_ar>=GAIN_SEX
has_en_amp = cna_en>=GAIN_SEX
has_only_TD = has_TD & cna_ar<GAIN_SEX
sum(has_nothing)
sum(has_AR_amp)
sum(has_en_amp)
sum(has_only_TD)

layout(matrix(1,1,1))
b=boxplot( 
    xp_ar[has_nothing], 
    xp_ar[has_AR_amp],
    xp_ar[has_TD],
    xp_ar[has_only_TD], 
    boxwex=0.5, outline=FALSE,
    col="lightgrey",
    las=1, ylim=c(0,16))


points( jitter( rep( 1, sum( has_nothing )), 2), 
        xp_ar[has_nothing], pch=1, col="#000000cc", cex=0.75 )  #66666cc6

points( jitter( rep( 2, sum( has_AR_amp )), 0.5), 
        xp_ar[has_AR_amp], pch=1, col="#000000cc", cex=0.75 )

points( jitter( rep( 3, sum( has_TD )), 2), 
        xp_ar[has_TD], pch=1, col="#000000cc", cex=0.75 )

points( jitter( rep( 4, sum( has_only_TD )), 2), 
        xp_ar[has_only_TD], pch=1, col="#000000cc", cex=0.75 )

table(curated_ar$code)
AR = data.frame( has_TD, has_only_TD, xp_ar, cna_ar, cna_en, td_ar, td_en )
AR = AR[order(cna_ar),]

layout(matrix(1:2,2,1))
par(mar=c(0.5,4,0.5, 0.5))
colors = rep("lightgrey", 101)
#colors[AR$has_TD] = "red"
barplot( AR$cna_ar, col=colors, ylim=c(0,10), xaxs="i" )
abline( GAIN_SEX, 0)
colors = matrix("white", nrow=2, ncol=101)
colors[1,AR$has_only_TD] = "black"
colors[2,AR$has_TD] = "black"
dimnames(colors)[[1]] = c("only TD", "has TD")
plot.color.grid( colors, block.height = 5, block.width = 5, 
                 space.X = 0, space.Y=0, show.axis.X = FALSE)

fisher.test( table( AR$has_only_TD, AR$cna_ar<GAIN_SEX ) )

# calculate chrX_counts and td_amp_counts
#-------------------------------------------------------------------------------

chrX_counts = rep(0,15601)
td_noamp_counts = rep(0,15601)
td_amp_counts = rep(0,15601)
has_TD_locus = rep(FALSE, 101)
has_AR_amp = rep(FALSE, 101)
has_TD_AR = rep(FALSE, 101)
has_locus_amp = rep(FALSE, 101)
for(i in 1:101){
    sample_id = sample_ids[i]
    local_counts = rep(0, 15601)
    bed = get_bed( sample_id, dir_build, allow_small = TRUE )
    bed=bed[bed$chrom=="chrX",]
    for(j in 1:dim(bed)[1]){
        if( bed$cn.int[j]>=GAIN_SEX ){
            x1 = round( bed$start[j]/10000, 0 )
            x2 = round( bed$stop[j]/10000, 0 )
            local_counts[x1:x2] = local_counts[x1:x2] + 1
        }
    }
    if( length(local_counts)>length(chrX_counts)){
        local_counts = local_counts[1:length(chrX_counts)]   
    }
    local_counts[local_counts>1]=1
    chrX_counts = chrX_counts + local_counts
    has_AR_amp[i] = local_counts[ 6762 ] > 0
    has_locus_amp[i] = local_counts[ 6692 ] > 0
    
    local_counts_td_noamp = rep(0, 15601)
    local_counts_td_amp = rep(0, 15601)
    sv = list_sv_m[list_sv_m$chrom_start=="chrX" & 
                       list_sv_m$sample_id==sample_id & list_sv_m$svtype=="TANDEM",]
    if( dim(sv)[1] ){
        for(j in 1:dim(sv)[1]){
            x1 = round( sv$pos_start[j] / 10000, 0 )
            x2 = round( sv$pos_end[j] / 10000, 0 )
            if( has_AR_amp[i] ){
                local_counts_td_amp[x1:x2] = local_counts_td_amp[x1:x2] + 1
            }else{
                local_counts_td_noamp[x1:x2] = local_counts_td_noamp[x1:x2] + 1
            }
            
        }
    }
    local_counts_td_amp[local_counts_td_amp>1]=1
    local_counts_td_noamp[local_counts_td_noamp>1]=1
    has_TD_locus[i]=local_counts_td_amp[6692]==1 | local_counts_td_noamp[6692]==1
    has_TD_AR[i] = local_counts_td_noamp[6762]==1 | local_counts_td_amp[6762]==1
    td_noamp_counts = td_noamp_counts + local_counts_td_noamp
    td_amp_counts = td_amp_counts + local_counts_td_amp
}



# Plot all individual tandem duplications in a stack sorted by pos_start
#-------------------------------------------------------------------------------
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_AR_TD_stack_closeup.pdf", 
     height=1.5, width=3.5)
par(mar=c(0,4,1,1))
x1 = 6600
x2 = 6850
sv_ar=list_sv_m[ (list_sv_m$svtype=="TANDEM" ) & list_sv_m$chrom_start=="chrX" &
                     list_sv_m$pos_start < 66900000 & list_sv_m$pos_end > 66900000, ]
sv_ar = sv_ar[order(sv_ar$pos_start),]
layout(matrix(1,1,1))
plot(-1,-1, ylim=c(0,25), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
cur_y=1
for(i in 1:dim(sv_ar)[1]){
    cur_col = "#cc00cc"
    sample_id = sv_ar$sample_id[i]
    if( ! has_AR_amp[ which(sample_ids==sample_id)]){
        cur_col="#000000" 
    }
    lines( c(sv_ar$pos_start[i]/10000, sv_ar$pos_end[i]/10000), c(cur_y,cur_y), lwd=1, col=cur_col )
    cur_y=cur_y+0.25
}
box()
show_gene( "AR", height=0.5, height_text = -1000) # show locus but not name
show_gene( "EDA2R", height=0.5, height_text = -1000) # show locus but not name
show_gene( "OPHN1", height=0.5, height_text = -1000) # show locus but not name
dev.off()


# plot copy number and SV
#-------------------------------------------------------------------------------

pdf( paste(dir_out, 'fig_ar_CN_SV.pdf', sep='/'), height=3, width=3.5)

layout(matrix(1:2,2,1), heights = c(0.5,1.5) )
par(mar=c(0,4,1,1))
x1 = 6600
x2 = 6850
plot(-1,-1, ylim=c(0,100), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=100,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=100,by=20), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, chrX_counts[i] ),lwd=4, col="#00000033" )
}
box()
par(mar=c(4,4,0.5,1))
plot(-1,-1, ylim=c(0,70), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in seq(from=10, to=70,by=10) ){
    abline(i,0,col="lightgrey")
}
axis(2, seq(from=0, to=60,by=20), las=1, cex.axis=0.75 )
for(i in x1:x2){
    lines( c(i,i), c(0, td_noamp_counts[i]+ td_amp_counts[i] ),lwd=2, col="#0000cc66" )
    #lines( c(i,i), c(0, td_noamp_counts[i]+ td_amp_counts[i] ),lwd=1, col="#d95f0e" )
}

for(i in x1:x2){
    lines( c(i,i), c(0, td_amp_counts[i] ),lwd=2, col="#cc00cc66" )
    #lines( c(i,i), c(0, td_amp_counts[i] ),lwd=1, col="#fec44f" )
}
for(i in x1:x2){
    lines( c(i,i), c(0, td_noamp_counts[i] ),lwd=2, col="#000000ee" )
    #lines( c(i,i), c(0, td_noamp_counts[i] ),lwd=1, col="#fff7bc" )
}
box()

dev.off()

# Plot AR expression boxplot
#-------------------------------------------------------------------------------

AR=data.frame( xp=as.numeric(matrix_tpm["AR",] ), 
               has_AR_amp, 
               has_TD_locus, 
               row.names=sample_ids )
AR$logxp = log2( AR$xp+1 )

write.table( AR, '/notebook/human_sequence_prostate_WCDT/WCDT/results/2018_04_24_table_AR_enhancer.txt',
             sep='\t', quote=FALSE)
anova( lm( AR$logxp ~ AR$has_AR_amp + AR$has_TD_locus ) )


write.table( AR, '/notebook/human_sequence_prostate_WCDT/WCDT/results/2018_04_19_AR_TD_amp_table.txt',
             quote=FALSE, sep='\t')
has_nothing =   !AR$has_TD_locus & !AR$has_AR_amp
has_only_amp =  !AR$has_TD_locus & AR$has_AR_amp
has_AR_amp = AR$has_AR_amp
has_only_amp =  !AR$has_TD_locus & AR$has_AR_amp
has_TD =  AR$has_TD_locus
has_only_TD =  AR$has_TD_locus & !AR$has_AR_amp
has_either = AR$has_TD_locus | AR$has_AR_amp

table(curated_ar$code)
rownames(curated_ar)[curated_ar$code=="ETD"]
sample_ids[has_only_TD]
# ha reassigned 

setdiff( sample_ids[has_nothing], rownames(curated_ar)[curated_ar$code=="NEU"])
setdiff( rownames(curated_ar)[curated_ar$code=="NEU"],  sample_ids[has_nothing])
# DTB-042-BL DTB-188-BL from has_only_TD to NEU (nothing)
# "DTB-021-BL" "DTB-083-BL" from nothing to EG (enhancer gain)

curated_has_nothing = curated_ar$code=="NEU"
curated_has_AR_amp = curated_ar$code=="EAG" | 
                     curated_ar$code=="AG" | 
                     curated_ar$code=="EG-AG" | 
                     curated_ar$code=="ETD-EAG" | 
                     curated_ar$code=="EATD-EAG" | 
                     curated_ar$code=="ETD-AG" 
sum( curated_has_AR_amp ) # 58
sum( as.numeric(matrix_CNA_int_ploidy["AR",] ) >= GAIN_SEX) # 70
sum(has_AR_amp) # 71


me_not_ha = setdiff( sample_ids[has_AR_amp] , sample_ids[curated_has_AR_amp] )
matrix_CNA_int_ploidy["AR",me_not_ha] 


#-------------------------------------------------------------------------------

pdf( paste(dir_out, 'fig_ar_CN_SV_boxplot.pdf', sep='/'), height=6, width=4)
layout(matrix(1,1,1))
b=boxplot( 
    AR$logxp[has_nothing], 
    AR$logxp[has_AR_amp],
    AR$logxp[has_TD],
    AR$logxp[has_only_TD], 
    boxwex=0.5, outline=FALSE,
    col=c("white", "white","white","white"),
    las=1, ylim=c(0,16))

points( jitter( rep( 1, sum(has_nothing)), 2), 
        AR$logxp[has_nothing], pch=1, col="#000000cc", cex=0.75 )  #66666cc6

points( jitter( rep( 2, sum(has_AR_amp)), 0.5), 
        AR$logxp[has_AR_amp], pch=1, col="#000000cc", cex=0.75 )

points( jitter( rep( 3, sum(has_TD)), 2), 
        AR$logxp[has_TD], pch=1, col="#000000cc", cex=0.75 )

points( jitter( rep( 4, sum( has_only_TD)), 2), 
        AR$logxp[has_only_TD], pch=1, col="#000000cc", cex=0.75 )
dev.off()

wilcox.test(AR$logxp[has_nothing], AR$logxp[has_AR_amp] ) 
wilcox.test(AR$logxp[has_nothing], AR$logxp[has_TD] ) 
wilcox.test(AR$logxp[has_nothing], AR$logxp[has_only_TD] ) 



# Plot ChIP-seq hotspots
#-------------------------------------------------------------------------------

# load lupien data and bin into 10000 nt
lup_AR = read.table( '/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_04_19_Lupien_AR.bg',
                     sep='\t', stringsAsFactors=FALSE)
names(lup_AR) = c("chrom", "start" , "stop", "score")

x1 = 6600
x2 = 6850
lup_means = rep(NA, x2-x1)
for( i in x1:x2){
    bound_start = i * 10000
    bound_end = (i+1) * 10000
    lup_means[i] = mean( lup_AR$score[
        (lup_AR$stop > bound_start & lup_AR$start < bound_end) |
            (lup_AR$start < bound_end & lup_AR$stop > bound_start) ])
}
lup_means[is.nan(lup_means)] = 0
lup_AR[!is.na(lup_means) & lup_means>2,]

pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_AR_lupien.pdf", 
     height=1.5, width=3.5)
layout(matrix(1,1,1))
par(mar=c(3,4,0.5,1))

plot(-1,-1, ylim=c(0,5), xlim=c(x1, x2), pch='.', las=1, 
     xaxs="i", yaxs="i", xlab="", axes=FALSE, ylab="")
for(i in x1:x2){
    lines( c(i,i), c(0, lup_means[i] ),lwd=1, col="black" )
}
axis(2, seq(from=0, to=15,by=5), las=1, cex.axis=0.75 )
box()
axis(1, at=seq(from=x1, to=x2, by=50),
     labels =seq(from=x1/100, to=x2/100, by=0.5), las=1 )
dev.off()

################################################################################
# FIGURE 4: Mutation signatures
################################################################################

sc = matrix_mutsig_cosmic
nmf = data.matrix(samples(sigs_nmf))
nmf = nmf[match.idx( paste("DTB-",get.split.col(sample_ids, "-", col=2 ),sep=''), dimnames(nmf)[[1]] )$idx.B,]
order.mutsig = order(matrix_samples$n_del_MHgt2_manta, decreasing = TRUE)

pdf( paste(dir_out, 'fig_signatures_HRD_barplot.pdf', sep='/'), height=1.5, width=5)
layout(matrix(1,1,1))
par(mar=c(0.1,0.5,0.5,1))
barplot(matrix_samples$n_del_MHgt2_manta[order.mutsig],
        names.arg = rep("", length(sample_ids)), xaxs="i", lwd=2 ,
        col="black", border="white", las=1, ylim=c(0,200),
        ylab="")
box()
dev.off()

postscript( paste(dir_out, 'fig_signatures_HRD_heatmap.ps', sep='/'), height=2, width=5)
layout(matrix(1,1,1))
par(mar=c(0.1,0.5,0.5,1))
cmap = grDevices::colorRampPalette( colors=color_palatte )(10)
sig_to_show=c("S3","S8")
scores_cosmic=t(matrix_mutsig_cosmic[,sig_to_show])
colors = color_scale( scores_cosmic, cmap, color_bounds = c(0, 0.5) )
dimnames(colors)[[1]] = paste("COSMIC",c(3,8))

color_palatte=c("white", "grey","black")
cmap = grDevices::colorRampPalette( colors=color_palatte )(10)
nmf_transformed = t(nmf)
colors_denovo = color_scale( nmf_transformed, cmap, color_bounds = c(0, 0.07) )
colors_denovo = colors_denovo[c(8,1),]
dimnames(colors_denovo)[[1]] = paste( "de novo", c(8,1))

colors_denovo = rbind( colors, colors_denovo[1,])
dimnames(colors_denovo)[[1]][1:3] = c("","","")
plot.color.grid( colors_denovo[,order.mutsig], block.height = 6, block.width=6, 
                 border = TRUE, cex.y=1, space.X=0 , space.Y=0, show.axis.X=FALSE,
                 show.axis.Y=TRUE)

dev.off()

postscript( paste(dir_out, 'fig_SV_BRCA_chromothripsis.ps', sep='/'), height=2, width=5)
layout(matrix(1,1,1))
par(mar=c(0.1,0.5,0.5,1))
colors = rep("white", 101)
colors[matrix_samples$has_chromothripsis] = "black"
colors = colors[order.mutsig]
colors = matrix(colors, nrow=1, ncol=101)
plot.color.grid( colors, block.height = 6, block.width=6, 
                 border = TRUE, cex.y=1, space.X=0 , space.Y=0, show.axis.X=FALSE,
                 show.axis.Y=TRUE)
dev.off()

postscript( paste(dir_out, 'fig_SV_BRCA_RNA_allele_boxplots.ps', sep='/'), height=2, width=5)
layout(matrix(1,1,1))
par(mar=c(0.1,0.5,0.5,0.5))
o=plot_landscape_segment( c("BRCA2", "BRCA1", "ATM", "CDK12"),
                          order_precomputed = order.mutsig,
                          show_histogram = FALSE, show_events = "break")
dev.off()

cr = matrix_samples$has_chromothripsis
summary( glm(cr ~ ALLELES_INACTIVATED["TP53",],
             family=binomial(link='logit')) )
ord.chr = order( ALLELES_INACTIVATED["TP53",], decreasing=TRUE )

layout(matrix(1:2,2,1), heights=c(1,1))
par(mar=c(0.1,0.5,0.5,0.5))
o=plot_landscape_segment( c("TP53"),
                          order_active=FALSE, order_inactive=FALSE,
                          order_precomputed = ord.chr,
                          show_histogram = FALSE, show_events = "break")
colors = rep("white", 101)
colors[matrix_samples$has_chromothripsis] = "black"
colors = colors[ord.chr]
colors = matrix(colors, nrow=1, ncol=101)
par(mar=c(0.1,0.5,0,0.5))
plot.color.grid( colors, block.height = 6, block.width=6, 
                 border = TRUE, cex.y=1, space.X=0 , space.Y=0, show.axis.X=FALSE,
                 show.axis.Y=TRUE)


ae_brca2 = allele_effect( "BRCA2", do_plot=FALSE, axis_override=NULL )
ae_brca1 = allele_effect( "BRCA1", do_plot=FALSE, axis_override=NULL )
only_bi = ae_brca2$alleles$bi | ae_brca1$alleles$bi
compound_het = !only_bi & ( ae_brca1$alleles$mono & ae_brca2$alleles$mono )
anova( lm( matrix_mutsig_cosmic[,3] ~ only_bi ) )
wilcox.test( matrix_mutsig_cosmic[only_bi,3], matrix_mutsig_cosmic[compound_het,3]   )
wilcox.test( matrix_mutsig_cosmic[!only_bi & !compound_het,3], matrix_mutsig_cosmic[compound_het,3]   )

postscript( paste(dir_out, 'fig_SV_BRCA_RNA_COSMIC3_vs_bi.ps', sep='/'), height=4, width=6)
layout(matrix(1,1,1))
par(mar=c(1,0.5,0.5,0.5))
boxplot( matrix_mutsig_cosmic[!only_bi & !compound_het,3],
         matrix_mutsig_cosmic[!only_bi & compound_het,3],
         matrix_mutsig_cosmic[only_bi & !compound_het,3], lwd=0.5, box.wex=0.5,
         col=c("white", "lightgrey", "darkgrey"), axes=FALSE, ylim=c(0,0.6))
axis(2, at=c(0, 0.2, 0.4), labels = c("","",""), las=1)
box()
dev.off()


postscript( paste(dir_out, 'fig_SV_BRCA1_expr_vs_alleles.ps', sep='/'), height=4, width=6)
layout(matrix(1,1,1))
xp = log( as.numeric(matrix_tpm["BRCA1",]+1) )
par(mar=c(1,0.5,0.5,3))
boxplot( xp[ae_brca1$alleles$bi] ,
         xp[ae_brca1$alleles$mono],
         xp[!ae_brca1$alleles$bi & ! ae_brca1$alleles$mono], lwd=0.5, 
         col=c("white", "lightgrey", "darkgrey"),axes=FALSE, ylim=c(0,5))
axis(2, at=0:4, labels=c("","","","",""))
box()
dev.off()

mc=matrix_samples$mutation_count
postscript( paste(dir_out, 'fig_SV_BRCA2_compound_mutation_freq.ps', sep='/'), height=4, width=6)
layout(matrix(1,1,1))
par(mar=c(3,3,1,1))
boxplot( log10( mc[!only_bi & !compound_het]/3000),
         log10( mc[!only_bi & compound_het]/3000),
                log10( mc[only_bi & !compound_het]/3000), lwd=0.5, box.wex=0.5,
         col=c("white", "lightgrey", "darkgrey"), axes=FALSE)
axis( 2, at=c(0.69897, 1, 1.69897, 2), labels = c("","","",""))
box()
dev.off()

logmc = log10( mc /3000)
ord = order(logmc, decreasing = TRUE)
compound_het[ord]

# There is a modest correlation between larger CNA freq and mutation rate
plot( 1-matrix_samples$percent_CNA_ref, logmc)
cor.test( 1-matrix_samples$percent_CNA_ref, logmc, method="spearman")

# Do compound het samples simply draw from those with higher overall CNA and 
# therefore by chance they have higher mutational freq?
# ANSWER: no group has higher overall percent CNA
boxplot( 
    1-matrix_samples$percent_CNA_ref[only_bi],
    1-matrix_samples$percent_CNA_ref[compound_het],
    1-matrix_samples$percent_CNA_ref[!only_bi & !compound_het])
    

################################################################################
# FIGURE 5A: LANDSCAPE 
################################################################################

x = matrix(FALSE, nrow=10, ncol=101)
x[1,] = allele_effect("ERG")$alleles$activating_sv
x[2,] = allele_effect("ETV1")$alleles$activating_sv
x[3,] = allele_effect("ETV4")$alleles$activating_sv
x[4,] = allele_effect("ETV5")$alleles$activating_sv
x[4,] = x[4,] | allele_effect("ETV5")$alleles$activating_cna
x[5,] = allele_effect("BRAF")$alleles$activating_missense
x[5,] = x[5,] | allele_effect("BRAF")$alleles$activating_sv
x[5,] = x[5,] | allele_effect("BRAF")$alleles$activating_cna
x[6,] = allele_effect("HRAS")$alleles$activating_missense
x[7,] = allele_effect("MYC")$alleles$activating_sv
x[8,] = allele_effect("CHD1")$alleles$n_alleles_inactivated==2
x[9,] = allele_effect("SPOP")$alleles$n_alleles_inactivated>0
x[10,] = allele_effect("IDH1")$alleles$inactivating_missense
ord = mutually_exclusive_order(x)

x_TP53 = allele_effect("TP53")$alleles$n_alleles_inactivated==2
x_sub = x_TP53[ ord[1:43] ]
ord[1:43] = ord[1:43][ order(x_sub, decreasing=TRUE) ]

pdf( paste(dir_out, 'fig_landscape_single_block.pdf', sep='/'), height=6, width=7.5)
layout(matrix(1:2,1,2), widths=c(7,1))
par(mar=c(1,5,0,0.5))

s5 = read.table('/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_05_09_symbols_figure_5.txt',
                header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)

o=plot_landscape_segment( rownames(s5),
                          cex.y=0.75, order_precomputed = ord,
                          show_nonfunctional=FALSE, show_events = s5$action)
dev.off()

o$n_samples_altered

pdf( paste(dir_out, 'fig_mutation_frequency.pdf', sep='/'), height=2, width=7.5)
layout(matrix(1:2,1,2), widths=c(8,1))
par(mar=c(1,5,0,0.5))
barplot( log10(matrix_mutcount$mutation_count[o]/3000), col='cornflowerblue',axes=FALSE,xaxs="i")
axis(2, labels = c("","",""), at=c(0, 1, 2), las=1)
dev.off()

bi=cbind( gene_locs[ which(INACTIVE[,1]>0),], INACTIVE[which(INACTIVE[,1]>0),] ) 
bi = bi[order(bi$chrom, bi$start),]



################################################################################
#                                                                              #
# SUPPLEMENTARY FIGURES                                                        #
#                                                                              #
################################################################################


################################################################################
# FIGURE: sequencing depth
################################################################################

pdf( paste(dir_out, 'fig_supplementary_sequencing_depth.pdf', sep='/'), height=6, width=6)
layout(matrix(1:4,2,2))
par(mar=c(4,6,2,1))
COL_N="#00ff0033"
COL_T="#0000ff33"
b=boxplot( as.numeric( matrix_summary_normal["Mean_coverage",] ), 
           as.numeric( matrix_summary_tumor["Mean_coverage",] ),
           names=c("normal", "tumor"), ylab="mean coverage",
           col=c( COL_N, COL_T ), las=1, box.wex=0.5, main="Mean coverage", ylim=c(0,130))

boxplot( as.numeric( matrix_summary_normal["Total_aligned_reads",] ), 
         as.numeric( matrix_summary_tumor["Total_aligned_reads",] ),
         names=c("normal", "tumor"),
         col=c( COL_N, COL_T ), las=1, box.wex=0.5, main="Total aligned reads")

boxplot( as.numeric( matrix_summary_normal["Percent_Q30_bases",] ), 
         as.numeric( matrix_summary_tumor["Percent_Q30_bases",] ),
         names=c("normal", "tumor"), ylab = "% bases with quality at Q30",
         col=c( COL_N, COL_T ), las=1, box.wex=0.5, main="Percent Q30 bases")

fl_n = as.numeric( matrix_summary_normal["Fragment_length_median",] )
fl_t =as.numeric( matrix_summary_tumor["Fragment_length_median",] )
boxplot( fl_n, fl_t, names=c("normal", "tumor"), main="Median insert length", 
         ylab="median insert length",
         col=c( COL_N, COL_T ), las=1, box.wex=0.5)
dev.off()


################################################################################
# FIGURE: Mean copy number
################################################################################

pdf("/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_mean_cna.pdf",
    height=4, width=10)
layout(matrix(1:1,1,1))
par(mar=c(1,3,2,1))
ys=rowMeans( data.matrix(cnbin[,3:103]), na.rm=TRUE )-2
xs = 1:length(ys)         
plot( -1, xlim=c(1,length(xs)),
      ylim=c(-1,3), xaxs="i", col=cols_cnbin, pch='.',
      main="mean CNA", axes=FALSE, ylab="")
for(i in 1:length(xs)){
    lines( c(xs[i], xs[i]), c(0, ys[i]), col=cols_cnbin[i] )
}
axis(2, at=-1:3, labels = 1:5, las=1)
for(i in seq(from=-1,to=3, by=1)){
    abline(i, 0, col="lightgrey")
}
box()
dev.off()

################################################################################
# FIGURE: ploidy
################################################################################

pdf("/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_ploidy.pdf",
    height=5, width=8)
par(mar=c(0.5,5,1,1))
ord = order(matrix_samples$ploidy,decreasing=TRUE)
colors=rep("grey", 101)
colors[ matrix_samples$ploidy[ord]<=2.5] = "#ffff0066"
colors[ matrix_samples$ploidy[ord]>=3] = "#00ffff66"
barplot(matrix_samples$ploidy[ord], xlab="Ploidy", las=1, xaxs="i",
        ylim=c(0,4), main="Computational estimate of mean ploidy",col=colors)
box()
legend( 70,3.75,
        c("diploid", "triploid", "undetermined"), 
        fill = c("#ffff0066", "#00ffff66", "grey"), 
        bty="n", cex=1.5)
dev.off()

################################################################################
# FIGURE: ploidy vs. mutations/translocations
################################################################################

pdf("/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_ploidy_mutcount_bnd.pdf",
    height=7, width=5)
layout(matrix(1:2,2,1))
par(mar=c(4,5,2,2))
d2 = density( log10( matrix_samples$mutation_count[ matrix_samples$ploidy<2.5]/3000 ) )
d3 = density( log10( matrix_samples$mutation_count[ matrix_samples$ploidy>3]/3000 ) )
plot(d2, ylim=c(0, 3), main="# mutations/Mb", 
     xaxs="i", yaxs="i", axes=FALSE, lwd=1, xlab="")
box()

polygon(d2, col="#ffff0066", border="blue", lwd=2) 
polygon(d3, col="#00ffff66", border="blue", lwd=2) 
axis(1, at=seq(from=1, to=3, by=1), labels = c(10, 100, 1000), las=1)
legend( 1.25,2.5,
        c("diploid", "triploid"), fill = c("#ffff0066", "#00ffff66"), 
        bty="n", cex=1.5)

par(mar=c(4,5,1,1))
d2 = density( matrix_sv_m$n_bnd[ matrix_samples$ploidy<2.5] )
d3 = density( matrix_sv_m$n_bnd[ matrix_samples$ploidy>3] )

plot(d2, xlim=c(0,400), ylim=c(0,0.01), main="# translocations", 
     xaxs="i", yaxs="i", axes=FALSE, lwd=1, xlab="megabases")
box()
polygon(d2, col="#ffff0066", border="blue", lwd=2) 
polygon(d3, col="#00ffff66", border="blue", lwd=2) 
axis(1,seq(from=0, to=400, by=50), las=1)
legend( 250,0.009,
        c("diploid", "triploid"), fill = c("#ffff0066", "#00ffff66"), 
        bty="n", cex=1.5)
dev.off()


################################################################################
#  FIGURE: PTEN SV detail diagram
################################################################################

# Top: locations of SV predicted to affect PTEN
curated_pten = curated_sv[curated_sv$threeprime=="PTEN" & curated_sv$reported=="yes",]

xlim_start =        87850000
pten_coding_start = 87859153
pten_coding_end =   87968775
xlim_end =          88000000

postscript( paste(dir_out, 'fig_pten_SV.ps', sep='/'), height=4, width=12)
par(mar=c(3,0.5, 0.5, 0.5))
plot( -100,-100, xlim=c(xlim_start,xlim_end), ylim=c(0,0.56), axes=FALSE, ylab="", xlab="")
box()
axis(1, at=seq(from=87850000, to=88000000, by=10000), labels=seq(from=87850000, to=88000000, by=10000)/1000000)
baseline = 0.05
lines(c(pten_coding_start,pten_coding_end), c(baseline, baseline) )

exon_starts = c(87864469,87894024,87925512,87931045,87933012,87952117,87957852,87960893,87965286)
exon_ends = c(87864548,87894109,87925557,87931089,87933251,87952259,87958019,87961118,87965473)
rect( 87863437, baseline-0.02, 87864469, baseline+0.02, col="black", border="black")
rect( 87965473, baseline-0.02, 87971930, baseline+0.02, col="black", border="black")
for(i in 1:length(exon_starts)){
    rect( exon_starts[i], baseline+0.05, exon_ends[i], baseline-0.05, col="black", border="black")
}
idx_breaks = which(curated_pten$mechanism=="breakend")
for(i in 1:length(idx_breaks)){
    if(curated_pten$pos3[idx_breaks[i]]==87868175 | curated_pten$pos3[idx_breaks[i]]==87870175){
        offset=0.05
    }else{
        offset=0
    }
    points( curated_pten$pos3[idx_breaks[i]], baseline+0.1-offset, pch="x", col="black", cex=2)
}

cur_y = baseline+0.15
idx_dels = which(curated_pten$mechanism=="deletion")
for(i in 1:length(idx_dels)){
    del_start = curated_pten$pos5[ idx_dels[i] ]
    del_end = curated_pten$pos3[ idx_dels[i] ]
    lines( c( del_start, del_end), c(cur_y, cur_y), col="red", lwd=1 )
    cur_y = cur_y + 0.01
}
cur_y = cur_y + 0.02
idx_inv = which(curated_pten$mechanism=="inversion")
for(i in 1:length(idx_inv)){
    inv_start = curated_pten$pos5[ idx_inv[i] ]
    inv_end = curated_pten$pos3[ idx_inv[i] ]
    lines( c( inv_start, inv_end), c(cur_y, cur_y), col="darkgreen", lwd=1 )
    cur_y = cur_y + 0.01
}
dev.off()

# Bottom: Lollypop plot of PTEN mutation loci
#-------------------------------------------------------------------------------
domain_starts_PTEN = c(1,   7, 186, 352,400)
domain_ends_PTEN =   c(7, 186, 352, 400,403)
domain_names_PTEN = c("PBD", "Phosphotase", "C2", "C-terminal", "PDZ")
domain_colors = c("lightblue", "cornflowerblue", "grey", "darkblue", "red")
loci=c(135,130,98,13)
locus_labels = c("p.I135R","p.R130Q","p.L98R","p.K13E")
loci_nonsense = c(32, 317, 258, 130, 233, 85 )
locus_labels_nonsense = c("p.I32fs", "p.V317fs", "p.F258fs", "p.R130X", "p.R233X", "p.V85X" )

postscript( paste(dir_out, 'fig_pten_mutations.ps', sep='/'), height=4, width=12)
par(mar=c(3,3,1, 1))
plot_lollyplot( loci=c(loci, loci_nonsense), 
                locus_labels=c(locus_labels, locus_labels_nonsense),
                locus_colors = c( rep("lightgreen",4), rep("darkgreen",6) ),
                domain_starts = domain_starts_PTEN,
                domain_ends = domain_ends_PTEN,
                domain_names = domain_names_PTEN,
                domain_colors= domain_colors )
dev.off()


################################################################################
# FIGURE: activating fusions
################################################################################

axl.order = order( as.numeric( matrix_tpm["AXL",] ) )
braf.order = order( as.numeric( matrix_tpm["BRAF",] ) )
myc.order = order( as.numeric( matrix_tpm["MYC",] ) )
pik3ca.order = order( as.numeric( matrix_tpm["PIK3CA",] ) )
map2k4.order = order( as.numeric( matrix_tpm["MAP2K4",] ) )

postscript( paste(dir_out, 'fig_SV_oncogene_RNA_barplots.ps', sep='/'), height=3.5, width=6)
layout( matrix(1:10,10,1), heights = rep(c(1,1), 5))
barplot_inactivations( "BRAF", window_size=1000, new.order=braf.order, show_events="activate" )
barplot_inactivations( "MYC", window_size=1000, new.order=myc.order, show_events="activate" )
barplot_inactivations( "AXL", window_size=1000, new.order=axl.order, show_events="activate" )
barplot_inactivations( "PIK3CA", window_size=1000, new.order=pik3ca.order, show_events="activate" )
barplot_inactivations( "MAP2K4", window_size=1000, new.order=map2k4.order, show_events="activate" )
dev.off()


################################################################################
# FIGURE: CDK12 Tandem duplication length plot
################################################################################

cdk12_bi = sample_ids[ allele_effect("CDK12")$alleles$bi ]
pdf( "/notebook/human_sequence_prostate_WCDT/WCDT/figures/fig_TD_size.pdf",
     height=5, width=8)

td_cdk12 = list_sv_m$svsize[list_sv_m$svtype=="TANDEM" & list_sv_m$sample_id %in% cdk12_bi]/1000000 
td_nocdk12 = list_sv_m$svsize[list_sv_m$svtype=="TANDEM" & !list_sv_m$sample_id %in% cdk12_bi]/1000000 
layout(matrix(1,1,1))
par(mar=c(4,4,2,1))
d=density(td_cdk12)
dno = density(td_nocdk12)
plot(d, xlim=c(0,10), main="Tandem duplication length", 
     xaxs="i", yaxs="i", axes=FALSE, lwd=1, xlab="megabases")
polygon(d, col="#ffff0066", border="blue", lwd=2) 
polygon(dno, col="#00ffff66", border="blue", lwd=2) 
box()
axis(1, seq(from=0, to=10, by=0.5), las=1)
legend( 6,0.7,
        c("CDK12 mutant", "CDK12 WT"), fill = c("#ffff0066", "#00ffff66"), 
        bty="n", cex=1.5)
dev.off()

################################################################################
# FIGURE: Mutation trinucleotide 
################################################################################

pdf( paste(dir_out, 'fig_supplementary_signature_profiles.pdf', sep='/'), 
     height=6, width=8)
layout(matrix(1,1,1))
SomaticSignatures::plotSignatures(sigs_nmf, normalize=TRUE)
dev.off()



################################################################################
# FIGURE: correlation of COSMIC And de novo signature scores and cluster
################################################################################

sc = matrix_mutsig_cosmic
nmf = data.matrix(samples(sigs_nmf))
nmf = nmf[match.idx( paste("DTB-",get.split.col(sample_ids, "-", col=2 ),sep=''), dimnames(nmf)[[1]] )$idx.B,]
nmf_no2 = nmf[,c(1,3,4,5,6,7,8)]
h = hclust(dist(t(nmf_no2)))

pdf( paste(dir_out,'/fig_denovo_cluster.pdf',sep='/'), height=1, width=4)
par(mar=c(0.5, 0.5, 0.5, 0.5))
plot(h, lwd=2, hang=1, axes=FALSE, xlab="", ylab="", main="", sub="",
     labels=rep("", 7))
dev.off()

nmf.cluster = nmf_no2[, h$order]
cors_fit = matrix( NA, nrow=dim(sc)[2], ncol=dim(nmf_no2)[2])
for(rr in 1:dim(cors_fit)[1]){
    for(cc in 1:dim(cors_fit)[2]){
        cors_fit[rr,cc] = cor.test( sc[,rr], nmf.cluster[,cc], method="spearman")$estimate
    }
}
dimnames(cors_fit)[[2]] = paste("S",c(1,3,4,5,6,7,8)[h$order])
dimnames(cors_fit)[[1]] = paste("COSMIC",1:dim(sc)[2] )
color_palatte=c("white", "slateblue1","darkblue")
cmap = grDevices::colorRampPalette( colors=color_palatte )(10)
cors.filter = cors_fit
colors = color_scale( cors.filter, cmap, color_bounds = c(0.2, 1) )
dimnames(colors)[[1]] = paste("COSMIC",1:30)
dimnames(colors)[[2]] = paste("de novo",c(1,3,4,5,6,7,8)[h$order])

pdf( paste(dir_out,'/SomSig_cor_denovo_cosmic.pdf',sep='/'), height=4, width=4)
layout(matrix(1,1,1))
par(mar=c(5,6,1,1))
show = rowSums( cors_fit>0.2, na.rm=TRUE ) > 0
plot.color.grid( colors[show,], block.height = 6, block.width=6, 
                 border = TRUE, cex.y=1, space.X=0 , space.Y=0, show.axis.X=TRUE,
                 show.axis.Y=TRUE)
dev.off()


################################################################################
#                                                                              #
# SUPPLEMENTARY TABLE DATA                                                     #
#                                                                              #
################################################################################


################################################################################
# Table S3
################################################################################
sv3sd = sv_bins_fig1a[which(sv_bins_fig1a$sv>=three_sd),]
sv3sd=sv3sd[order(sv3sd$sv, decreasing = TRUE),]

sv3sd$assignment = rep("none", dim(sv3sd)[1])
for(i in 1:dim(sv3sd)[1]){
    idx = which( assignments$chrom==sv3sd$chrom[i] & assignments$start==sv3sd$start[i])
    if(length(idx)==1){
        sv3sd$assignment[i] = assignments$assignment[idx]   
    }
}
sv3sd=sv3sd[,c("amp", "del", "sv", "chrom", "start", "end", "assignment", "zscore", "pval")]
write.table( sv3sd, '/notebook/human_sequence_prostate_WCDT/WCDT/results/sv3sd.txt',
             quote=FALSE, sep='\t', row.names=FALSE)


HC = matrix_samples$has_chromothripsis

pvals = rep(NA, dim(ALLELES_INACTIVATED)[1])
for(i in 1:dim(ALLELES_INACTIVATED)[1]){
    ae = ALLELES_INACTIVATED[i,] != 2
    if( sum( ae )>2 ){
        ss = summary( glm(cr ~ TP53 + ae, family=binomial(link='logit')) )
        if( dim(ss$coefficients)[1] == 3 ){
            pvals[i] = ss$coefficients[3,4]
        }
    }
}


write.table( list_sv_m, '/notebook/human_sequence_prostate_WCDT/drafts/Cell/resubmission/table_SV.txt',
             sep='\t', quote=FALSE, row.names=FALSE)

# Chromoplexy analysis
##########################

# 11 samples with zero chains
sum(matrix_chain$Maximum_number_of_chromosomes_involved_in_one_chain==0, na.rm=TRUE)
# 50 with chains involving 3 or more chromosomes;
sum(matrix_chain$Maximum_number_of_chromosomes_involved_in_one_chain>2, na.rm=TRUE)


#noncoding

nc = read.table('/Users/david/Documents/notebook/human_sequence_prostate_WCDT/noncoding/2018_06_09_recur_merge.csv',
                sep=',', header=TRUE, stringsAsFactors=FALSE)
nc$loc = paste(nc$Chr,nc$Pos,sep=':')
ncp = nc[nc$nctype=="Promoter(1kb)",]
ncu = nc[nc$nctype=="UTR",]

cosmic = read.table('/notebook/human_sequence_prostate_WCDT/WCDT/metadata/2018_04_15_COSMIC_tier1_cancer_genes.txt',
                    stringsAsFactors=FALSE)$V1

head(cosmic)
m = match.idx( cosmic, nc$Gene, allow.multiple.B = TRUE)
nc_cosmic = nc[m$idx.B,]
write.table( nc_cosmic, '/notebook/human_sequence_prostate_WCDT/WCDT/results/secondary_analysis/noncoding_cosmic_mutations.txt',
             sep='\t', row.names=FALSE, quote=FALSE)

freq.locs = names( table(ncp$loc)[ which( as.numeric( table( ncp$loc ) )>3 ) ] )
nc[nc$loc %in% freq.locs, ]

freq.locs = names( table(ncu$loc)[ which( as.numeric( table( ncu$loc ) )>2 ) ] )

tail(nc[nc$loc %in% freq.locs, ],50)


#write.table( 
#    data.frame(
#        CDK12=allele_effect("CDK12")$alleles$bi,
#        BRCA2=allele_effect("BRCA2")$alleles$bi,
#        row.names=sample_ids),
#    '/notebook/human_sequence_prostate_WCDT/prepub_WCDT/results/WCDT_BRCA2_CDK12.txt', sep='\t',quote=FALSE, row.names=TRUE)
