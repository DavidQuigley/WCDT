library(SomaticSignatures)
library(RCircos)
library(trackViewer)
source('/notebook/code/src/R/quantitative_genetics.R')

N_SAMPLES = 101
GAIN_NONSEX = 3
LOSS_SINGLE_NONSEX = 1.65
LOSS_DOUBLE_NONSEX = 0.6
GAIN_SEX = 1.4
LOSS_SEX = 0.6

build_id = '2018_04_15'
dir_root = "/notebook/human_sequence_prostate_WCDT/WCDT"
#dir_results = paste( '/Volumes/datasets_1/human_sequence_prostate_WCDT/build_' , build_id, sep='')
dir_results = paste( dir_root, "/results/build_", build_id, sep='')
dir_meta = paste( dir_root, "/metadata", sep='')

fn_curated = paste( dir_meta,'/2018_05_07_curated_AR_SV.txt',sep='')
fn_sa = paste(dir_meta, '20180123_sample_attributes.txt', sep='/')
fn_clinical = paste(dir_meta, '2018_01_10_clinical_attributes.txt', sep='/')
fn_cyto = paste( dir_meta, '/HG38_cytoband_std_chrom.txt', sep='')
fn_gene_locs = paste(dir_meta, 'GRCh38Decoy_refseq_genelocs_from_refFlat.bed',sep='/')
fn_chrom_lengths = paste(dir_meta, 'HG38_chromosome_lengths.txt',sep='/')
fn_centromeres = paste(dir_meta, 'HG38_centromere_loci.txt',sep='/')
fn_curated_sv = paste(dir_meta, '/2018_03_15_curated_SV_table.txt', sep='')
fn_curated_missense = paste(dir_meta, '/2018_05_14_curated_missense_table.txt', sep='')
fn_curated_fs = paste(dir_meta, '/2018_05_14_curated_frameshift_table.txt', sep='')
fn_symbols_fig5 = paste(dir_meta, '/2018_05_09_symbols_figure_5.txt', sep='')

fn_mutsig=paste(dir_results, '/', build_id, '_matrix_mutation_signature_summary.txt', sep='')
fn_cna_en = paste( dir_results, '/', build_id, '_matrix_AR_enhancer_CN.txt', sep='')
fn_sample_summary = paste(dir_results, '/', build_id, '_matrix_sample_summary.txt', sep='')
fn_binned_CN = paste(dir_results, '/', build_id, '_matrix_binned_weighted_CN_copycat.txt', sep='')
fn_cna = paste(dir_results, '/', build_id, '_matrix_CNA_symbol_copycat.txt', sep='')
fn_cna_int = paste(dir_results, '/', build_id, '_matrix_CN_integer_symbol_copycat.txt', sep='')
fn_germline = paste(dir_results, '/', build_id, '_matrix_germline_summary.txt', sep='')
dir_scores_chromothripsis = paste( dir_results, "/../secondary_analysis/chromothripsis", sep='')
fn_somatic = paste(dir_results, '/', build_id, '_list_somatic_strelka.txt', sep='')
fn_somatic_mutect = paste(dir_results, '/', build_id, '_list_somatic_mutect.txt', sep='')
fn_missense = paste(dir_results, '/', build_id, '_matrix_somatic_missense_strelka.txt', sep='')
fn_inactivating = paste(dir_results, '/', build_id, '_matrix_somatic_inactivating_strelka.txt', sep='')
fn_mutcount = paste(dir_results, '/', build_id, '_matrix_mutation_count_summary.txt', sep='')
fn_sigs_nmf = paste(dir_results, '/../secondary_analysis/signature_analysis/sigs_nmf.Rdata', sep='')
fn_gof_nmf = paste(dir_results, '/../secondary_analysis/signature_analysis/gof_nmf.Rdata', sep='')
fn_sigs_mm = paste(dir_results, '/../secondary_analysis/signature_analysis/sigs_mm.Rdata', sep='')
fn_manta_sv = paste(dir_results, '/', build_id, '_list_manta_SV.txt', sep='')
fn_tpm = paste(dir_results, '/', build_id, '_matrix_rna_tpm.txt', sep='')
fn_alignment_summary_t = paste(dir_results, '/', build_id, '_matrix_alignment_summary_tumor.txt', sep='')
fn_alignment_summary_n = paste(dir_results, '/', build_id, '_matrix_alignment_summary_normal.txt', sep='')
fn_genomewide_sliding_scan = paste(dir_results, '/genome_wide_SV_scan.txt', sep='')
fn_assignment = paste( dir_meta, '/sv_window_assignment.txt', sep='')
fn_chain =  paste(dir_results, '/', build_id, '_matrix_chainfinder_del278_fdr5.txt', sep='')

# manually curated assignments for table S3
assignments = read.table(fn_assignment, sep='\t', stringsAsFactors=FALSE, 
                         header=TRUE)

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
FLAGS = c('TTN','MUC16','OBSCN','AHNAK2','SYNE1','FLG','MUC5B','DNAH17','PLEC',
          'DST','SYNE2','NEB','HSPG2','LAMA5','AHNAK','HMCN1','USH2A','DNAH11','MACF1',
          'MUC17','DNAH5','GPR98','FAT1','PKD1','MDN1','RNF213','RYR1','DNAH2','DNAH3',
          'DNAH8','DNAH1','DNAH9','ABCA13','APOB','SRRM2','CUBN','SPTBN5','PKHD1','LRP2',
          'FBN3','CDH23','DNAH10','FAT4','RYR3','PKHD1L1','FAT2','CSMD1','PCNT','COL6A3',
          'FRAS1','FCGBP','DNAH7','RP1L1','PCLO','ZFHX3','COL7A1','LRP1B','FAT3','EPPK1',
          'VPS13C','HRNR','MKI67','MYO15A','STAB1','ZAN','UBR4','VPS13B','LAMA1','XIRP2',
          'BSN','KMT2C','ALMS1','CELSR1','TG','LAMA3','DYNC2H1','KMT2D','BRCA2','CMYA5',
          'SACS','STAB2','AKAP13','UTRN','VWF','VPS13D','ANK3','FREM2','PKD1L1','LAMA2',
          'ABCA7','LRP1','ASPM','MYOM2','PDE4DIP','TACC2','MUC2','TEP1','HELZ2','HERC2',
          'ABCA4')

# manually added FLAGS
FLAGS.add = c("PCDHA1","PCDHA2","PCDHA3","PCDHA4","PCDHA5","PCDHA6","PCDHGA1",
              "PCDHGA2","PCDHGA3","PCDHGA4","PCDHGA5","PCDHGB1","PCDHGB2",
              "PCDHA10","PCDHA7","PCDHA8","PCDHA9","PCDHGA6","PCDHGA7","PCDHGA8",
              "PCDHGB3","PCDHGB4","PCDHGB5","PCDHGA9","PCDHGB6",
              "FLG2", "TCHH", "IVL", "RPTN", "DSPP",
              "MUC7","MUC17", "MUC19", "MUC21", "MUC22", "MUC4", "MUC5AC",
              "ZFHX4","ZNF283","ZNF428","ZNF430","ZNF512","ZNF519","ZNF544","ZNF57","ZNF676",
              "ZNF678","ZNF681","ZNF700","ZNF709","ZNF717","ZNF813","ZNF90","ZNF92","ZNF93","ZNF99",
              "ZNF208","ZNF257","ZNF43","ZNF845","ZNF850","ZNF91")
FLAGS.remove = c( "KMT2C", "KTM2D", "BRCA2")
FLAGS = sort( setdiff( c(FLAGS, FLAGS.add), FLAGS.remove ) )


assess_tandem_dup = function( chrom, locus, list_sv, chrom_lengths, samples_ids ){
    # return vector of which samples have a TD intersecting a given locus
    has_TD_locus = rep(FALSE, N_SAMPLES)
    for(i in 1:N_SAMPLES){
        sample_id = sample_ids[i]
        n_bins=ceiling( chrom_lengths$V2[ which(rownames(chrom_lengths)==chrom) ] /10000 )
        local_counts_td = rep(0,n_bins)
        sv = list_sv[list_sv$chrom_start==chrom & 
                     list_sv$sample_id==sample_id & 
                     list_sv$svtype=="TANDEM",]
        if( dim(sv)[1] ){
            for(j in 1:dim(sv)[1]){
                x1 = round( sv$pos_start[j] / 10000, 0 )
                x2 = round( sv$pos_end[j] / 10000, 0 )
                local_counts_td[x1:x2] = local_counts_td[x1:x2] + 1
            }
        }
        has_TD_locus[i]=local_counts_td[locus]>=1 
    }
    has_TD_locus
}



allele_effect=function( symbol, do_plot=FALSE, axis_override=NULL ){
    # missense mutations all come from curated data 
    # external references:
    #sample_ids gene_locs 
    #matrix_SV matrix_tpm matrix_germline matrix_CNA_int_ploidy
    #curated_fs curated_missense
    
    has_somatic_inactivation = rep(FALSE,N_SAMPLES)
    has_activating_missense = rep(FALSE, length(sample_ids))
    has_inactivating_missense = rep(FALSE, length(sample_ids))
    has_activating_sv = rep(FALSE, length(sample_ids))
    has_inactivating_sv = rep(FALSE, length(sample_ids))
    has_inactive = matrix_inactive[symbol,] != '.'
    has_germline = matrix_germline[symbol,] != '.'
    has_loh = as.character(matrix_CNA[symbol,]) == "LOH"
    
    
    chrom = gene_locs$chrom[ which( rownames(gene_locs)==symbol ) ]
    
    # copy number
    cna = as.numeric( matrix_CNA_int_ploidy[symbol,] )
    threshold_single =  LOSS_SINGLE_NONSEX   
    threshold_double =  LOSS_DOUBLE_NONSEX   
    threshold_amp = GAIN_NONSEX   
    cn.hi.threshold = GAIN_NONSEX
    if( chrom=="chrX" | chrom=="chrY" ){
        threshold_single =  LOSS_SEX
        threshold_double =  LOSS_SEX   
        threshold_amp = GAIN_SEX
        cn.hi.threshold = GAIN_SEX
    }
    
    # expression
    xp =  as.numeric( matrix_tpm[symbol,] )
    val_mad_all=mad(xp, na.rm=TRUE)
    val_median_all = median(xp, na.rm=TRUE)
    
    # mutations
    if( symbol %in% unique(curated_fs$symbol) ){
        samples_inactive = curated_fs$sample_id[curated_fs$gene==symbol]
        has_somatic_inactivation[ match.idx( sample_ids, samples_inactive)$idx.A ] = TRUE
    }
    
    if( symbol %in% unique( curated_missense$gene ) ){
        samples_active = curated_missense$sample_id[curated_missense$gene==symbol & curated_missense$consequence=="activate"] 
        has_activating_missense[ match.idx( sample_ids, samples_active)$idx.A ] = TRUE
        samples_inactive = curated_missense$sample_id[curated_missense$gene==symbol & curated_missense$consequence!="activate"] 
        has_inactivating_missense[ match.idx( sample_ids, samples_inactive)$idx.A ] = TRUE
        matrix_missense[symbol, !has_activating_missense & !has_inactivating_missense] = '.'
    }
    
    # sv
    idx = which( curated_sv$threeprime==symbol & curated_sv$consequence=="break" ) 
    if(length(idx)>0){
        for(i in 1:length(idx)){
            sample_id = curated_sv$sample[idx][i]
            idx_sample_id = which(sample_ids==sample_id)
            has_inactivating_sv[idx_sample_id] = TRUE
            matrix_SV[ symbol, idx_sample_id ] = curated_sv$mechanism[idx[i]]
        }
    }
    idx = which( curated_sv$threeprime==symbol & curated_sv$consequence=="activate" ) 
    if(length(idx)>0){
        for(i in 1:length(idx)){
            sample_id = curated_sv$sample[idx][i]
            idx_sample_id = which(sample_ids==sample_id)
            has_activating_sv[idx_sample_id] = TRUE
            matrix_SV[ symbol, idx_sample_id ] = curated_sv$mechanism[idx[i]]
        }
    }
    
    df = data.frame(
        activating_missense = has_activating_missense,
        inactivating_missense = has_inactivating_missense,
        nonsense = has_inactive,
        activating_sv = has_activating_sv,
        inactivating_sv = has_inactivating_sv,
        inactivating_germline = has_germline,
        CNA_1 = cna <= threshold_single,
        CNA_2 = cna <= threshold_double,
        CNA_amp = cna >= threshold_amp,
        LOH = has_loh,
        xp=xp,
        activating_cna = cna >= cn.hi.threshold & !is.na(xp) & xp >= 100
    )
    df$activated = df$activating_missense | df$activating_cna | df$activating_sv
    
    alleles_gone = rowSums( df[,c("inactivating_missense", "LOH", "nonsense", 
                                  "inactivating_sv", "inactivating_germline", 
                                  "CNA_1")])
    has_biallelic =    df$CNA_2 | alleles_gone > 1
    has_monoallelic = !df$CNA_2 & alleles_gone == 1
    
    # special case for CDK12 with two somatic mutations
    if( symbol=="CDK12" ){
        has_biallelic[ which(sample_ids=="DTB-214-BL")]=TRUE
        has_monoallelic[ which(sample_ids=="DTB-214-BL")]=FALSE
    }
    
    n_mono = sum( has_monoallelic )
    n_bi = sum( has_biallelic )
    n_amp = sum( df$CNA_amp )
    biallelic_because_sv = has_biallelic & 
        !df$CNA_2 & df$inactivating_sv  & 
        rowSums( df[,c("inactivating_missense", "LOH", "nonsense", 
                       "inactivating_germline", "CNA_1")] ) == 1
    df$bi = has_biallelic
    df$mono = has_monoallelic
    df$bi_from_sv = biallelic_because_sv
    n_bi_from_sv = sum( biallelic_because_sv )
    
    n_alleles = rep(0, N_SAMPLES)
    n_alleles[has_biallelic] = 2
    n_alleles[has_monoallelic] = 1
    logxp = log(xp+1)
    if(do_plot){
        if( is.null(axis_override)){
            boxplot( logxp[n_alleles==2], logxp[n_alleles==1], logxp[n_alleles==0], 
                     las=1, box.wex=0.25, 
                     cex.axis=0.75, lwd=0.5,
                     col=c("darkgrey", "lightgrey", "white"),
                     names=c("","",""))
        }else{
            boxplot( logxp[n_alleles==2], logxp[n_alleles==1], logxp[n_alleles==0], 
                     las=1, box.wex=0.25, 
                     cex.axis=0.75, lwd=0.5, axes=FALSE,
                     col=c("darkgrey", "lightgrey", "white"),
                     names=c("","",""), ylim=c( min(axis_override), max(axis_override) ))
            axis( 4, axis_override, las=1, cex.axis=0.75)
            axis( 1, c(1,2,3), labels = c("","","") )
            box()
        }
    }
    df$n_alleles_inactivated = n_alleles
    cor_estimate = NA
    cor_pval = NA
    if( sum( !is.na(xp) ) > 50 & ( n_bi > 2 | n_mono > 2) ){
        cc = cor.test( xp, n_alleles )
        cor_estimate = cc$estimate
        cor_pval = cc$p.value
    }
    list( 
        n_bi = n_bi, n_mono = n_mono, n_bi_from_sv = n_bi_from_sv, n_amp=n_amp,
        cor_alleles_xp = signif( as.numeric(cor_estimate), 3),
        pval_alleles_xp = signif( cor_pval, 3 ),
        alleles = df
    )
}


allele_effect_sv = function( sv_name, fn_curated, fn_cna_en ){
    # specialized version of allele_effect for AR enhancer locus that uses
    # curated calls
    all_FALSE = rep(FALSE, N_SAMPLES)
    df=data.frame( 
        activating_missense = all_FALSE,
        inactivating_missense = all_FALSE,
        nonsense = all_FALSE,
        activating_sv = all_FALSE,
        inactivating_sv = all_FALSE,
        inactivating_germline = all_FALSE,
        CNA_1 = all_FALSE,
        CNA_2 = all_FALSE,
        CNA_amp = all_FALSE,
        LOH = all_FALSE,
        xp = rep(NA, N_SAMPLES),
        activating_cna = all_FALSE
    )
    if( sv_name=="AR enhancer"){
        ARCN = read.table(fn_cna_en,
                          header=TRUE, stringsAsFactors=FALSE, sep='\t')
        curated_ar = read.table(fn_curated, header=TRUE, stringsAsFactors=FALSE, 
                                sep='\t', row.names=1)
        df$CNA_amp = ARCN$copies >= GAIN_SEX
        df$activating_cna = df$CNA_amp
    }
    df
}


plot_lollyplot = function(loci,
                          locus_labels,
                          xaxis_positions=NULL,
                          locus_colors="cornflowerblue",
                          domain_starts,
                          domain_ends,
                          domain_names,
                          domain_colors){
    unique_labels = unique(locus_labels)
    unique_loci = rep(0, length(unique_labels))
    unique_colors = rep("", length(unique_labels))
    n_seen = rep(1, length(unique_labels))
    for(i in 1:length(unique_labels)){
        unique_colors[i] = locus_colors[ locus_labels==unique_labels[i] ][1]
        n_seen[i] = sum( locus_labels == unique_labels[i] )
        unique_loci[i] = loci[ locus_labels==unique_labels[i] ][1]
    }
    if( is.null( xaxis_positions ) ){
        xaxis_positions = xaxis_positions   
    }
    mutsites <- GRanges("chr1", IRanges(unique_loci, width=1, names=unique_labels))
    mutsites$score = n_seen
    mutsites$color = unique_colors
    domain_lengths = domain_ends - domain_starts
    domains = GRanges("chr1", IRanges(domain_starts, 
                                      width=domain_lengths,
                                      names=domain_names))
    domains$fill = domain_colors
    domains$label = domain_names
    domains$height = 0.05
    layout(matrix(1,1,1))
    xaxis = sort( c( domain_starts, domain_ends ) )
    par(mar=c(5,4,2,1))
    mutsites$label.parameter.rot = 90
    mutsites$label.parameter.cex = 0.5
    lolliplot(mutsites, features=domains, yaxis=FALSE, xaxis=xaxis)
}


get_bed = function( sample_id, dir='/notebook/human_sequence_prostate_WCDT/WCDT/results/build_2018_04_15', allow_small=FALSE ){
    dir_bed=dir
    bed = read.table(paste(dir_bed,"/copycat/",sample_id,'_copycat.bed',sep=''),stringsAsFactors=FALSE)
    names(bed)=c("chrom", "start", "stop", "call", "cn.int", "strand")
    if(!allow_small){
        bed = bed[bed$stop-bed$start > 1000,]
    }
    bed
}

get_sv_by_sample = function( sample_id ){
    sv=list_sv_m[list_sv_m$sample_id==sample_id, ]
    sv = data.frame( "Chromosome" = sv$chrom_start,
                     "chromStart" = sv$pos_start,
                     "chromEnd" = sv$pos_start+1,
                     "Chromosome.1" = sv$chrom_end,
                     "chromStart.1" = sv$pos_end,
                     "chromEnd.1" = sv$pos_end+1,
                     "svtype" = sv$svtype,
                     "sample_id" = sv$sample_id,
                     stringsAsFactors=FALSE
    )
    m = match.idx( rownames(chrom_lengths), sv$Chromosome, allow.multiple.B = TRUE)$idx.B
    m.1 = match.idx( rownames(chrom_lengths), sv$Chromosome.1, allow.multiple.B = TRUE)$idx.B
    sv = sv[ intersect( m, m.1) , ]
    sv$PlotColor=rep("#2c7bb6", dim(sv)[1])
    sv$PlotColor[sv$svtype=="DEL"] = "#2C7BB6"
    sv$PlotColor[sv$svtype=="BND"] = "#FDAE61"
    sv$PlotColor[sv$svtype=="INV"] =  "#ABD9E9"
    sv$PlotColor[sv$svtype=="TANDEM"] =  "#D7191C"
    sv$PlotColor[sv$svtype=="INS"] = "#FFFFBF"    
    sv
}

plot_circos = function( cyto_info, cna, sv, gene_info=NULL, fn_out="", main="",
                        chr.include=NULL){
    
    if( is.null(cna) ){
        track_sv = 1   
    }else{
        track_cna = 1
        track_sv = 2
    }
    if(is.null( chr.include )){
        cna_show=cna
        sv_show=sv
        chr.exclude=NULL
    }else{
        if(!is.null(cna)){
            m = match.idx( chr.include, cna$chrom, allow.multiple.B = TRUE)
            cna_show = cna[m$idx.B,]
        }
        m1 = match.idx( chr.include, sv$Chromosome, allow.multiple.B = TRUE)
        m2 = match.idx( chr.include, sv$Chromosome.1, allow.multiple.B = TRUE)
        sv_show = sv[intersect( m1$idx.B, m2$idx.B ), ]
        chr.exclude = setdiff( rownames(chrom_lengths), chr.include )
    }
    
    cytoBandData = RCircos.Validate.Cyto.Info( cyto_info, chr.exclude )
    RCircos.Initialize.Plot.Parameters( tracks.inside=track_sv, tracks.outside=0 )
    RCircos.Set.Cytoband.Data( cytoBandData )
    RCircos.Set.Base.Plot.Positions()
    RCircos.Set.Plot.Area();
    rcircos.params <- RCircos.Get.Plot.Parameters();
    rcircos.params$text.size = .75;
    RCircos.Reset.Plot.Parameters(rcircos.params);
    
    plot.window(c(-1.3,1.3), c(-1.3, 1.3))
    if( fn_out != "" ){
        pdf(file=fn_out, height=10, width=10)
        par(mar=c(2,2,2,2))
        plot.new()
        plot.window(c(-1.5,1.5), c(-1.5, 1.5))
    }
    RCircos.Chromosome.Ideogram.Plot( )
    
    if( !is.null(gene_info) ){
        if( !is.null(chr.include)){
            m1 = match.idx( chr.include, gene_info$chrom, allow.multiple.B = TRUE)
            gene_info = gene_info[m1$idx.B,]
        }
        RCircos.Gene.Name.Plot(gene_info, "symbol", track_sv, "in");
    }
    if( !is.null( cna ) ){
        RCircos.Line.Plot(cna_show, data.col=5, track.num=track_cna, side="in")
    }
    if( sum( sv_show$svtype=="INV")>0 ){
        RCircos.Link.Plot(
            link.data=sv_show[sv_show$svtype=="INV",],
            track.num=track_sv, 
            by.chromosome=FALSE)
    }
    if( sum( sv_show$svtype=="TANDEM")>0 ){
        RCircos.Link.Plot(
            link.data=sv_show[sv_show$svtype=="TANDEM",],
            track.num=track_sv, 
            by.chromosome=FALSE)    
    }
    if( sum( sv_show$svtype=="BND")>0 ){
        RCircos.Link.Plot(
            link.data=sv_show[sv_show$svtype=="BND",],
            track.num=track_sv, 
            by.chromosome=FALSE)
    }
    if( sum( sv_show$svtype=="DEL")>0 ){
        RCircos.Link.Plot(
            link.data=sv_show[sv_show$svtype=="DEL",],
            track.num=track_sv, 
            by.chromosome=FALSE)
    }
    if( main != "" ){
        text(-1.1,1.25, main)   
    }
    if( fn_out != "" ){
        dev.off()
    }
}



mutually_exclusive_order = function(x){
    # order(...) for arbitrary number of columns   
    cols = 1:dim(x)[2]
    final_order = cols
    last_in_col=1
    for( rr in 1:dim(x)[1]){
        x_cur = x[,cols]
        if(!is.vector(x_cur)){
            if( sum( x_cur[rr,] ) > 0 ){
                ord = order(x_cur[rr,], decreasing=TRUE)
                x[,cols] = x_cur[,ord]
                final_order[cols] = final_order[cols][ord]
                cur_last_in_col = tail( which(x[rr,] > 0 ), n=1 )
                if( length(cur_last_in_col)>0 ){
                    if( cur_last_in_col>last_in_col){
                        last_in_col=cur_last_in_col
                    }
                    if( last_in_col < dim(x)[2] ){
                        cols = (last_in_col+1) : (dim(x)[2])
                    }
                }else{
                    break   
                }
            }
        }
    }
    final_order
}

plot_landscape_segment = function( fig_symbols, 
                                   order_precomputed=NULL,
                                   cex.y=1, 
                                   show_histogram=TRUE,
                                   show_nonfunctional=TRUE,
                                   show_events="all"){
    if( is.null(order_precomputed) ){
        order.FV = 1:N_SAMPLES
    }else{
        order.FV = order_precomputed
    }
    n.alt = plot.landscape.v4(
        symbol_list=fig_symbols, sample_order=order.FV,
        block.height=2, block.width=3, show_x_axis = FALSE,
        cex.y=cex.y, show_nonfunctional=show_nonfunctional, 
        show_events=show_events) 
    par(mar=c(1,0,0,0.5))
    if( show_histogram ){
        plot.landscape.freq( n.alt, block.height=2)
    }
    n.alt
}


plot.landscape.freq = function(n.alt, block.height){
    total.height = ( block.height * n.alt$n_symbols )
    lwd=2
    plot( -1, -1, xlim=c(0, 100), ylim=c(0, total.height), 
          axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
    rect(0,0,100, total.height, col="#eeeeee")
    box()
    cur.y = total.height
    for(i in seq(from=10, to=100, by=10)){
        lines(c(i,i), c(0, total.height), col="white")
    }
    # group frequencies of each combination of cn, sv, and mut together
    for( i in 1 : n.alt$n_symbols ){
        cur_x=1
        x = matrix(FALSE, nrow=7, ncol=N_SAMPLES)
        x[1,] = !n.alt$cn[i,] & !n.alt$sv[i,] & n.alt$mut[i,]
        x[2,] = !n.alt$cn[i,] & n.alt$sv[i,] & n.alt$mut[i,]
        x[3,] = n.alt$cn[i,] & n.alt$sv[i,] & n.alt$mut[i,]
        x[4,] = !n.alt$cn[i,] & n.alt$sv[i,] & !n.alt$mut[i,]
        x[5,] = n.alt$cn[i,] & n.alt$sv[i,] & !n.alt$mut[i,]
        x[6,] = n.alt$cn[i,] & !n.alt$sv[i,] & !n.alt$mut[i,]
        x[7,] = n.alt$cn[i,] & !n.alt$sv[i,] & n.alt$mut[i,]
        ord = mutually_exclusive_order(x)
        for( oo in 1 : length(ord) ){
            j = ord[oo]
            upper_third = cur.y - (block.height/3)
            lower_third = cur.y - block.height + (block.height/3)
            mid = cur.y - (block.height/2)
            colors = c()
            
            if( n.alt$cn[i,j] ){ colors = c(colors, "#d7191c")} # red
            if( n.alt$sv[i,j] ){ colors = c(colors, "#fdae61")} # orange
            if( n.alt$mut[i,j] ){ colors = c(colors, "chartreuse4")} # green
            if( length(colors)==1 ){
                lines( c(cur_x,cur_x), c(cur.y, cur.y-block.height), col=colors[1], lwd=lwd)
            }else if(length(colors)==2){
                lines( c(cur_x,cur_x), c(cur.y, mid), col=colors[1], lwd=lwd)
                lines( c(cur_x,cur_x), c(mid, cur.y-block.height), col=colors[2], lwd=lwd)
            }else if(length(colors)==3){
                lines( c(cur_x,cur_x), c(cur.y, upper_third), col=colors[1], lwd=lwd)
                lines( c(cur_x,cur_x), c(upper_third, lower_third), col=colors[2], lwd=lwd)
                lines( c(cur_x,cur_x), c(lower_third, cur.y-block.height), col=colors[3], lwd=lwd)
            }
            if( length(colors>0)){
                cur_x=cur_x+1
            }
        }
        lines(c(0,100), c(cur.y, cur.y), col="black")
        lines(c(0,100), c(cur.y-block.height, cur.y-block.height), col="black")
        cur.y = cur.y - block.height
    }
}

plot.landscape.v4=function( symbol_list, sample_order=1:N_SAMPLES,
                            block.height=6, block.width=10, 
                            cex.y=1, show_x_axis=TRUE, 
                            show.y.labels=TRUE, show_nonfunctional=TRUE,
                            show_events="all"){
    
    # show_nonfunctional means show monoallelic inactivating events
    show_mono = show_nonfunctional
    col.noevent="#eeeeee"
    col.CN.amp="red"
    col.CN.single = "cornflowerblue"
    col.CN.double = "darkblue"
    col.CN.LOH = "darkgrey"
    col.fusion = "purple"
    col.missense="chartreuse4"
    col.missense.active="yellow"
    col.germline="black"
    col.svbroken="orange"
    col.svactivated="hotpink"
    col.inactive="saddlebrown"
    spacer_labels=0
    
    if( length( symbol_list )==1 ){
        # R is obnoxious about 1 row matrixes
        symbol_list = c(symbol_list, "KRT5" )
        spacer_labels = 1
    }
    n.rows = length(symbol_list) - spacer_labels
    n.cols = length(sample_ids)
    total.width =  ( block.width * n.cols) 
    total.height = ( block.height * n.rows )
    has_alt_sv = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    has_alt_cna = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    has_alt_mut = matrix(FALSE, nrow=length(symbol_list), ncol=n.cols)
    n.alterations = rep(0, length(symbol_list) )
    
    plot(0,0,col=col.noevent, xlim=c(0,total.width), ylim=c(0,total.height), 
         axes=F, xlab="", ylab="", bg=col.noevent, xaxs="i", yaxs="i")
    rect(0, 0, total.width, total.height, col=col.noevent)
    cur.y = total.height
    symbol_ctr=1
    xlab_locs=rep(0, n.cols)
    ylab_locs=rep(0, length(symbol_list) )
    
    if( length(show_events)==1){
        show_events = rep(show_events, n.rows)   
    }
    for(rr in 1:n.rows){
        this.x.left = 0
        if( symbol_list[rr]=="AR enhancer" ){
            alleles = allele_effect_sv( "AR enhancer", fn_curated, fn_cna_en )
        }else{
            ae=allele_effect(symbol_list[rr])
            alleles = ae$alleles
        }
            
        for(i in 1:n.cols){
            cc = sample_order[i]
            this.x.right = this.x.left+block.width
            this.x.right.halfway = this.x.left+(block.width/2)
            this.y.top = cur.y 
            this.y.mid = cur.y - ( (block.height)/2 )
            this.y.bot = cur.y - block.height
            this.y.upperthird = cur.y - ( (block.height)/3 )
            this.y.lowerthird = cur.y - block.height + ( (block.height)/3 )
            
            rect( this.x.left, cur.y - block.height, 
                  this.x.right, cur.y, col=col.noevent, border=NA)
            
            col.top = col.noevent
            col.bot = col.noevent
            show_amp = FALSE
            show_cn1 = FALSE
            show_cn2 = FALSE
            show_loh = FALSE
            show_SV_ACT = FALSE
            show_SV_INACT = FALSE
            show_mis_ACT = FALSE
            show_mis_INACT = FALSE
            show_non = FALSE
            show_germ = FALSE
            
            colors = c()
            xlab_locs[cc] = mean( c(this.x.left, this.x.right) )
            show_amp = alleles$activating_cna[cc]
            show_SV_ACT = alleles$activating_sv[cc]
            show_mis_ACT = alleles$activating_missense[cc]
            show_cn2 =  alleles$CNA_2[cc]
            show_cn1 =  alleles$CNA_1[cc]
            show_loh =  alleles$LOH[cc]
            show_SV_INACT = alleles$inactivating_sv[cc]
            show_mis_INACT =  alleles$inactivating_missense[ cc ]
            show_non = alleles$nonsense[ cc ]
            show_germ = alleles$inactivating_germline[ cc ]
            
            if( show_SV_ACT & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="activate")){
                colors = c(colors, col.svactivated); 
                has_alt_sv[rr,cc]=1
            }
            if( show_SV_INACT & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.svbroken)
                has_alt_sv[rr,cc]=TRUE
            }
            if( show_germ & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.germline)
                has_alt_mut[rr,cc]=TRUE
            }
            if( show_non & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.inactive)
                has_alt_mut[rr,cc]=TRUE
            }
            if( show_mis_INACT & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.missense)
                has_alt_mut[rr,cc]=TRUE
                # show two missense for CDK12 DTB-214-BL
                if( symbol_list[rr]=="CDK12" & dimnames(matrix_CNA_int_ploidy)[[2]][cc]=="DTB-214-BL"){
                    colors = c(colors, col.missense)
                }
            }
            if( show_loh & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.CN.LOH)
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_amp & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="activate")){
                colors = c(colors, col.CN.amp)
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_mis_ACT & !show_mis_INACT & !show_cn2 & (show_events[rr]=="all" | show_events[rr]=="activate")){
                colors = c(colors, col.missense.active)
                has_alt_mut[rr,cc]=TRUE
            }
            if( show_cn2 & (show_events[rr]=="all" | show_events[rr]=="break")){
                colors = c(colors, col.CN.double)
                has_alt_cna[rr,cc]=TRUE
            }
            if( show_cn1 & !show_cn2  & (show_events[rr]=="all" | show_events[rr]=="break")){
                if( show_nonfunctional | length(colors)>0 ){
                    colors = c(colors, col.CN.single)
                    has_alt_cna[rr,cc]=TRUE
                }
            }
            if( length(colors)>0 ){
                n.alterations[rr]=n.alterations[rr]+1
            }
            if( length(colors)==0 ){
                colors = c(col.noevent, col.noevent)
            }else if( length(colors)==1){
                if( colors[1]==col.CN.double ){
                    colors = c(col.CN.double,col.CN.double)
                }else{
                    colors = c(colors, col.noevent)
                }
            }
            if( length( colors) == 2 ){
                rect( this.x.left, this.y.top, this.x.right, this.y.mid, 
                      col=colors[1], border=NA)
                rect( this.x.left, this.y.mid, this.x.right, this.y.bot, 
                      col=colors[2], border=NA)
            }else{
                rect( this.x.left, this.y.top, this.x.right, this.y.upperthird, 
                      col=colors[1], border=NA)
                rect( this.x.left, this.y.upperthird, this.x.right, this.y.lowerthird, 
                      col=colors[2], border=NA)
                rect( this.x.left, this.y.lowerthird, this.x.right, this.y.bot, 
                      col=colors[3], border=NA)
            }
            this.x.left = this.x.left+block.width
        }
        cur.y = cur.y - block.height 
        symbol_ctr=symbol_ctr+1
    }
    for( cc in 1:n.cols ){
        line_loc = (cc-1)*block.width
        lines( c(line_loc,line_loc  ), c(0, total.height), col="white")    
    }
    if( spacer_labels==0 ){
        for( rr in seq(from=1, to=n.rows, by=1 ) ){
            line_loc = (rr-1)*block.height
            abline( line_loc, 0, col="black")
        }
    }
    if( spacer_labels==1 ){
        symbol_list = symbol_list[1]
        n.alterations = n.alterations[1]
    }
    results = list( cna=has_alt_cna, 
                    mut=has_alt_mut, 
                    sv=has_alt_sv, 
                    n_symbols=length(symbol_list),
                    n_samples_altered = n.alterations)
    ylab_locs = 0:(length(symbol_list)-1) * block.height + (block.height/2) 
    ylab_locs = ylab_locs[length(ylab_locs):1]
    if( !show.y.labels ){
        symbol_list=rep("", length(symbol_list))
    }
    axis(2, at=ylab_locs, labels=symbol_list, las=2, cex.axis=cex.y, 
         tick=FALSE, hadj=1, line=-0.5, font=3)
    if( show_x_axis ){
        axis(1, at=xlab_locs, labels = sample_ids, cex.axis=0.75, las=2 )
    }
    results
}

barplot_inactivations = function( symbol, window_size=1000, new.order=NULL, show_events="all" ){
    xp = log( as.numeric( matrix_tpm[symbol,] ) + 1 )
    if( is.null( new.order ) ){
        new.order = order(xp)
    }else{
        xp = xp[new.order]
    }
    ylim = c(0, ceiling( max ( xp, na.rm=TRUE) ) )
    par(mar=c(0,5,0.5,1))
    barplot( xp, las=1, xaxs="i", axes=FALSE, col="black", border=NA, ylim=ylim )
    box()
    axis(2, at=ylim, labels=ylim, las=1 )
    par(mar=c(1,5,0.5,1))
    o=plot_landscape_segment( symbol, order_precomputed=new.order, 
                              show_histogram = FALSE, show_events=show_events)
}

sv_targets = function( SV, symbol, window ){  
    idx = which(rownames(gene_locs)==symbol)
    chrom = gene_locs$chrom[ idx ]
    pos_5 = gene_locs$start[ idx ] - window
    pos_3 = gene_locs$end[ idx ] + window
    res = SV[ ( SV$chrom_start==chrom & 
                     SV$pos_start >= pos_5 & SV$pos_start <= pos_3 ),]
    sv2 = SV[ ( SV$chrom_end==chrom & 
                    SV$pos_end >= pos_5 & SV$pos_end <= pos_3 ),]
    names(sv2) = names(res)
    idx=match.idx( rownames(sv2), setdiff( rownames( sv2 ), rownames(res )) )$idx.A
    res=rbind( res, sv2[idx,] )
    
    res = res[order(res$chrom_end, res$pos_end ),]
    res
}

color_scale = function( V, color_map, color_bounds=NA, color_NA=NA ){
    if( is.na( color_bounds[1] ) ){
        Vmax = max(V, na.rm=TRUE)
        Vmin = min(V, na.rm=TRUE)
    }else{
        Vmin = color_bounds[1]
        Vmax = color_bounds[2]
    }
    increment = (Vmax-Vmin) / (length(color_map)-1)
    lookup = seq(from=Vmin, to=Vmax, by=increment )
    if( is.vector(V) ){
        out = rep("", length(V))
        for(i in 1:length(V)){
            out[i] = color_map[ which( lookup>V[i] )[1] ]
        }
    }else if( is.matrix(V) ){
        out = matrix("", nrow=dim(V)[1], ncol=dim(V)[2], 
                     dimnames=list( dimnames(V)[[1]], dimnames(V)[[2]]))
        for(rr in 1:dim(V)[1]){
            for(cc in 1:dim(V)[2]){
                col=color_map[ which( lookup>=V[rr,cc] )[1]]
                out[rr,cc] = col 
            }
        }
    }
    if( !is.na( color_NA) )
        out[is.na(out)] = color_NA
    out
}


################################################################################
# Genome data
################################################################################
cyto_info = read.table(fn_cyto, stringsAsFactors=FALSE,header=TRUE, sep='\t')
chrom_lengths = read.table(fn_chrom_lengths, row.names=1, stringsAsFactors=FALSE)
centromeres = read.table( fn_centromeres, header=TRUE, stringsAsFactors = FALSE)
gene_locs = read.table( fn_gene_locs, stringsAsFactors = FALSE)
names(gene_locs) = c("chrom", "start", "end", "symbol", "score", "strand")
gene_locs$symbol = get.split.col(gene_locs$symbol, "~", first=TRUE)
seen = hsh_new()
chrom_valid = hsh_from_vectors(rownames(chrom_lengths), 1:dim(chrom_lengths)[1])
keep = rep(FALSE, dim(gene_locs)[1])
for(i in 1:dim(gene_locs)[1]){
    keep[i] = hsh_in( chrom_valid, gene_locs$chrom[i]) &
              !hsh_in( seen, gene_locs$symbol[i] )
    hsh_set( seen, gene_locs$symbol[i], 1 )
}
gene_locs = gene_locs[keep,]
rownames(gene_locs) = gene_locs$symbol
chrom_idx = get.split.col( gene_locs$chrom, "chr", last=TRUE)
chrom_idx[chrom_idx=="X"]=23
chrom_idx[chrom_idx=="Y"]=24
gene_locs$chrom_idx = as.numeric(chrom_idx)
rm(keep)
rm(seen)
rm(chrom_idx)


################################################################################
# Sample matrix file
################################################################################

matrix_samples = read.table(fn_sample_summary, header=TRUE, sep='\t',
                            stringsAsFactors=FALSE, check.names = FALSE, 
                            row.names = 1)

#canonical list and order set here
sample_ids = rownames(matrix_samples)  
n_samples = length(sample_ids)
sa = load.matrix( fn_sa )


################################################################################
# CNA: binned in 3mb windows across genome
################################################################################

cnbin = read.table(fn_binned_CN,
                   header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                   check.names = FALSE)

cols_cnbin = rep("blue", dim(cnbin)[1])
cols_cnbin[ cnbin$chrom %in% c("chr2","chr4","chr6","chr8","chr10","chr12",
                               "chr14","chr16","chr18","chr20","chr22",
                               "chrY")] = "coral"

matrix_CNA = read.table(fn_cna,
                        header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                        check.names = FALSE, row.names = 1)
matrix_CNA_int = read.table(fn_cna_int,
                            header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                            check.names = FALSE, row.names = 1)

cnamp = rowSums( cnbin[,3:dim(cnbin)[2]] >=3 ) / dim(cnbin)[2]
cndel = rowSums( cnbin[,3:dim(cnbin)[2]] <=1 ) / dim(cnbin)[2]
sexchrom = cnbin$chrom=="chrX" | cnbin$chrom=="chrY"
cnamp[sexchrom] = rowSums( cnbin[sexchrom,3:dim(cnbin)[2]] >=1.5 ) / dim(cnbin)[2]
cndel[sexchrom] = rowSums( cnbin[sexchrom,3:dim(cnbin)[2]] <=0.75 ) / dim(cnbin)[2]
cndel = -1 * cndel

# make sure genes in matrix_CNA, matirx_CNA_int match gene_locs 
m = match.idx( dimnames(matrix_CNA)[[1]], rownames(gene_locs))
matrix_CNA = matrix_CNA[m$idx.A,]
matrix_CNA_int = matrix_CNA_int[m$idx.A,]
gene_locs = gene_locs[m$idx.B,]

matrix_CNA_int_ploidy=matrix_CNA_int

# canonical symbol list and order set here
#-------------------------------------------------------------------------------
symbols = dimnames(matrix_CNA)[[1]]
n_symbols = length(symbols)


# recalculate percent altered 
is_sex = cnbin$chrom=="chrX" | cnbin$chrom=="chrY"
n_alt_notsex = colSums( cnbin[!is_sex,3:103] >= GAIN_NONSEX, na.rm=TRUE ) +
               colSums( cnbin[!is_sex,3:103] <= LOSS_SINGLE_NONSEX, na.rm=TRUE )
n_alt_sex = colSums( cnbin[is_sex,3:103] >= GAIN_SEX, na.rm=TRUE ) +
            colSums( cnbin[is_sex,3:103] <= LOSS_SEX, na.rm=TRUE )

percent_alt = (n_alt_notsex + n_alt_sex) / dim(cnbin)[1]
matrix_samples$percent_CNA_ref = 1-percent_alt


################################################################################
# Germline mutations that potentially are pathogenic
################################################################################

germ = read.table(fn_germline, header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                  check.names = FALSE, row.names = 1)
germ = germ[ match.idx( sample_ids, rownames(germ))$idx.B,]
germ = t(germ)
matrix_germline = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))

m_symbol = match.idx( rownames(matrix_germline), rownames(germ))
m_sample = match.idx( dimnames(matrix_germline)[[2]], dimnames(germ)[[2]] )
matrix_germline[ m_symbol$idx.A, m_sample$idx.B] = germ[m_symbol$idx.B, m_sample$idx.B]


######################################################################
# Calculate chromothripsis
######################################################################

chrom_maxima = matrix(0, N_SAMPLES,24)
for(i in 1:length(sample_ids)){
    sc = read.table( paste(dir_scores_chromothripsis,'/',sample_ids[i],
                           '_chromo_scores_2e+07.txt',sep=''), sep='\t', 
                     header=TRUE)
    print( paste("Loading chromothripsis scores:", sample_ids[i] ))
    for(j in 1:24){
        chrom=rownames(chrom_lengths)[j]   
        sc_chrom = sc[sc$chrom==chrom,]
        for(x in 1:50){
            if( sum( sc_chrom$n_inv>=x & sc_chrom$n_del>=10 & 
                     sc_chrom$n_cna_switch>=x )>0){
                chrom_maxima[i,j]=x   
            }else{
                break
            }
        }
    }
}

chrom_maxima_filtered = chrom_maxima
chrom_maxima_filtered[chrom_maxima_filtered<15]=0
list_chromo=data.frame( which(chrom_maxima_filtered>0, arr.ind = TRUE) )
list_chromo[,1] = sample_ids[list_chromo[,1]]
list_chromo$col = paste("chr", list_chromo$col, sep='')
list_chromo = list_chromo[order(list_chromo[,1]),]
names(list_chromo) = c("sample_id", "chrom")
list_chromo$chrom[list_chromo$chrom=="chr23"] = "chrX"
has_chromothripsis = rep(FALSE, N_SAMPLES)
has_chromothripsis[match.idx( sample_ids, list_chromo$sample_id)$idx.A] = TRUE
matrix_samples$has_chromothripsis = has_chromothripsis

# review question: how stable is the TP53 result to varying definitions of chromothripsis?
chrom_maxima_filtered = chrom_maxima
chrom_maxima_filtered[chrom_maxima_filtered<20]=0
list_chromo=data.frame( which(chrom_maxima_filtered>0, arr.ind = TRUE) )
list_chromo[,1] = sample_ids[list_chromo[,1]]
list_chromo$col = paste("chr", list_chromo$col, sep='')
list_chromo = list_chromo[order(list_chromo[,1]),]
names(list_chromo) = c("sample_id", "chrom")
list_chromo$chrom[list_chromo$chrom=="chr23"] = "chrX"
has_chromothripsis = rep(FALSE, N_SAMPLES)
has_chromothripsis[match.idx( sample_ids, list_chromo$sample_id)$idx.A] = TRUE

###############################################
# Somatic mutation data
###############################################

list_somatic_data = read.table(fn_somatic,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                          check.names = FALSE)

# load inactivating somatic mutations; row is sample, col is symbol
som_inactive = read.table(fn_inactivating, header=TRUE, sep='\t', 
                          stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
m = match.idx( sample_ids, rownames(som_inactive) )
som_inactive = som_inactive[m$idx.B,]
matrix_inactive = matrix('.', nrow=n_symbols, ncol=n_samples, 
                  dimnames=list( y=symbols, x=sample_ids))
som_inactive=t(som_inactive)
m = match.idx( rownames(matrix_inactive), rownames(som_inactive))
matrix_inactive[ m$idx.A,] = som_inactive[m$idx.B,]
rm(som_inactive)

# load missense somatic mutations; row is sample, col is symbol
som_missense = read.table(fn_missense,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
matrix_missense = matrix('.', nrow=n_symbols, ncol=n_samples, 
                  dimnames=list( y=symbols, x=sample_ids))
m = match.idx( sample_ids, rownames(som_missense) )
som_missense = som_missense[m$idx.B,]
som_missense=t(som_missense)
m = match.idx( rownames(matrix_missense), rownames(som_missense))
matrix_missense[ m$idx.A,] = som_missense[m$idx.B,]
rm(som_missense)

matrix_mutcount  = read.table(fn_mutcount, header=TRUE, sep='\t', 
                              stringsAsFactors=FALSE, check.names = FALSE,
                            row.names = 1)
matrix_mutcount = cbind(matrix_mutcount, dummy=rep('.', dim(matrix_mutcount)[1]), stringsAsFactors=FALSE)
m = match.idx( sample_ids, rownames(matrix_mutcount))
matrix_mutcount = matrix_mutcount[m$idx.A,]

list_somatic_data_mutect = read.table(fn_somatic_mutect,
                               header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                               check.names = FALSE)

# load inactivating somatic mutations; row is sample, col is symbol
som_inactive_mutect = read.table(fn_inactivating, header=TRUE, sep='\t', 
                          stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
m = match.idx( sample_ids, rownames(som_inactive_mutect) )
som_inactive_mutect = som_inactive_mutect[m$idx.B,]
matrix_inactive_mutect = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))
som_inactive_mutect=t(som_inactive_mutect)
m = match.idx( rownames(matrix_inactive_mutect), rownames(som_inactive_mutect))
matrix_inactive_mutect[ m$idx.A,] = som_inactive_mutect[m$idx.B,]
rm(som_inactive_mutect)

# load missense somatic mutations; row is sample, col is symbol
som_missense_mutect = read.table(fn_missense,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
matrix_missense_mutect = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))
m = match.idx( sample_ids, rownames(som_missense_mutect) )
som_missense_mutect = som_missense_mutect[m$idx.B,]
som_missense_mutect=t(som_missense_mutect)
m = match.idx( rownames(matrix_missense_mutect), rownames(som_missense_mutect))
matrix_missense_mutect[ m$idx.A,] = som_missense_mutect[m$idx.B,]
rm(som_missense_mutect)

matrix_mutcount  = read.table(fn_mutcount, header=TRUE, sep='\t', 
                              stringsAsFactors=FALSE, check.names = FALSE,
                              row.names = 1)
matrix_mutcount = cbind(matrix_mutcount, dummy=rep('.', dim(matrix_mutcount)[1]), stringsAsFactors=FALSE)
m = match.idx( sample_ids, rownames(matrix_mutcount))
matrix_mutcount = matrix_mutcount[m$idx.A,]


# compare mutect and strelka
mu = paste(list_somatic_data_mutect$sample_id,list_somatic_data_mutect$symbol, sep='|')
st = paste(list_somatic_data$sample_id,list_somatic_data$symbol, sep='|')

################################################################################
# Mutation signatures
################################################################################

matrix_mutsig_cosmic = read.table(fn_mutsig, header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                    check.names = FALSE, row.names = 1)
m = match.idx( sample_ids, rownames(matrix_mutsig_cosmic))
matrix_mutsig_cosmic = matrix_mutsig_cosmic[m$idx.A,]

# De novo
load(fn_sigs_nmf)
#load(fn_gof_nmf)
load(fn_sigs_mm)
nmf=samples(sigs_nmf)
m=match.idx( dimnames(nmf)[[1]], 
           paste( "DTB", get.split.col(sample_ids, "-", col=2), sep='-' ))
nmf = nmf[m$idx.A,]

################################################################################
# Structural variants from manta (list_sv_m) 
################################################################################
list_sv_m = read.table(fn_manta_sv,
                  header=TRUE, sep='\t', 
                  stringsAsFactors=FALSE, check.names = FALSE)
list_sv_m = list_sv_m[list_sv_m$chrom_start != "chrM",]

# eliminate artifact on chromosome 22
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3996660/
list_sv_m= list_sv_m[!(list_sv_m$chrom_start=="chr22" & list_sv_m$svtype=="BND" & list_sv_m$pos_start>=28597176 & list_sv_m$pos_start<=29037622),]

chrom_idx = get.split.col( list_sv_m$chrom_start, "chr", last=TRUE)
chrom_idx[chrom_idx=="X"] = 23
chrom_idx[chrom_idx=="Y"] = 24
chrom_idx = as.numeric(chrom_idx)
list_sv_m = cbind( list_sv_m, chrom_idx, stringsAsFactors=FALSE)

# remove duplicate BND events
keep = rep(TRUE, dim(list_sv_m)[1])
idx_BND = which(list_sv_m$svtype=="BND")
seen = hsh_new()

as_written = paste(list_sv_m$chrom_start, list_sv_m$pos_start, list_sv_m$chrom_end, list_sv_m$pos_end)
reversed = paste(list_sv_m$chrom_end, list_sv_m$pos_end, list_sv_m$chrom_start, list_sv_m$pos_start)

for(i in 1:length(idx_BND)){
    if( !hsh_in( seen, as_written[idx_BND[i]] ) & 
        !hsh_in( seen, reversed[idx_BND[i]] ) ){
        hsh_set( seen, as_written[idx_BND[i]], 1)
        keep[idx_BND[i]] = TRUE
    }else{
        keep[idx_BND[i]] = FALSE
    }
}
list_sv_m = list_sv_m[keep,]

write.table( list_sv_m, 
             '/notebook/human_sequence_prostate_WCDT/WCDT/results/build_2018_04_15/2018_04_27_list_sv_m.txt',
              quote=FALSE, sep='\t')

##############################################
# Summarize counts of SV events
##############################################
matrix_sv_m = data.frame( n_dels=matrix_samples$n_del_manta,
                          n_tandems=matrix_samples$n_dup_manta,
                          n_mh = matrix_samples$n_del_MHgt2_manta,
                          n_ins=matrix_samples$n_ins_manta,
                          n_inv=matrix_samples$n_inv_manta,
                          per=round(matrix_samples$n_del_MHgt2_manta/matrix_samples$n_del_manta,2), 
                          row.names=dimnames(matrix_samples)[[1]], 
                          stringsAsFactors=FALSE)
matrix_sv_m$n_bnd = as.numeric(table(list_sv_m$sample_id[list_sv_m$svtype=="BND"]))



################################################################################
# SV bins (for figure 1a, others)
################################################################################

for(i in 1:24){
  chrom = rownames(chrom_lengths)[i]
  n_bins = ceiling( round( chrom_lengths$V2[i] / 1000000 ) )
  sv_bins = rep(0, n_bins )
  cna_bins_amp = rep(0, n_bins )
  cna_bins_del = rep(0, n_bins )
  bin_start = 0
  bin_end = 1000000
  bins_centromere = floor( round( centromeres$chromStart[i] / 1000000, 0 ) ) :
    ceiling( round( centromeres$chromEnd[i] / 1000000, 0 ) ) 
  bin_starts = rep(0, length(sv_bins))
  centro_start = centromeres$chromStart[i]
  centro_end = centromeres$chromEnd[i]
  for(j in 1:length(sv_bins)){
    bin_starts[j] = bin_start
    idx = which( list_sv_m$chrom_start==chrom & 
                   ((list_sv_m$pos_start>bin_start & list_sv_m$pos_start <= bin_end ) |
                      (list_sv_m$pos_end>bin_start & list_sv_m$pos_end <= bin_end ) |
                      (list_sv_m$svtype=="TANDEM" & list_sv_m$pos_start>bin_start & list_sv_m$pos_end < bin_end )))
    sv_bins[j] = length( unique( list_sv_m$sample_id[idx] ) )
    if( chrom=="chr6" & bin_start==32000001){
      sv_bins[j] = 0
    }
    if( (bin_start > (centro_start-500000) & bin_start < (centro_end+500000) ) |
        (bin_end > (centro_start-500000) & bin_end < (centro_end+500000))){
      sv_bins[j] = 0
    }
    idx = which( cnbin$chrom==chrom & 
                   (bin_start+1)>=cnbin$bin_start & 
                   bin_end<(cnbin$bin_start+3000000))
    if( chrom=="chrX" | chrom=="chrY"){
      cna_bins_amp[j] = sum( cnbin[idx,3:103] > GAIN_SEX, na.rm=TRUE )
      cna_bins_del[j] = sum( cnbin[idx,3:103] < LOSS_SEX, na.rm=TRUE )
    }else{
      cna_bins_amp[j] = sum( cnbin[idx,3:103] > GAIN_NONSEX, na.rm=TRUE )
      cna_bins_del[j] = sum( cnbin[idx,3:103] < LOSS_SINGLE_NONSEX, na.rm=TRUE )
    }
    bin_start = bin_end+1
    bin_end = bin_end + 1000000
  }
  sv_bins[ bins_centromere ] = -1
  cna_bins_del[ bins_centromere ] = -1
  cna_bins_amp = cna_bins_amp/200 * -1
  cna_bins_del = cna_bins_del/200 
  cna_bins_amp[is.na(cna_bins_amp)]=0
  cna_bins_del[is.na(cna_bins_del)]=0
  sv_bins = sv_bins/100
  sv_bins[1] = 0 # telomeres
  if(i==1){
    sv_bins_fig1a = data.frame( chrom=rep(chrom,length(cna_bins_del)),
                                amp=cna_bins_amp,
                                del=cna_bins_del,
                                sv=sv_bins,
                                start=bin_starts, end=bin_starts+1000000,
                                stringsAsFactors=FALSE)
  }else{
    sv_bins_fig1a = rbind(sv_bins_fig1a,
                          data.frame( chrom=rep(chrom,length(cna_bins_del)),
                                      amp=cna_bins_amp,
                                      del=cna_bins_del,
                                      sv=sv_bins,
                                      start=bin_starts, end=bin_starts+1000000,
                                      stringsAsFactors=FALSE))
  }
}

# write out sv_bins_fig1a
SVB = sv_bins_fig1a
SVB$amp = SVB$amp * -200
SVB$del = SVB$del * 200
SVB = SVB[,c(1,5,6,2,3,4)]
SVB$Z = round( (SVB$sv - mean(SVB$sv)) / sd(SVB$sv), 3)
SVB = SVB[order(SVB$sv, decreasing=TRUE),]
fn_data_fig1a='/notebook/human_sequence_prostate_WCDT/WCDT/results/secondary_analysis/data_supporting_figure_1a.txt'
write.table(SVB, fn_data_fig1a, quote=FALSE, sep='\t', row.names=FALSE)

################################################################################
# RNA
################################################################################

tpm = read.table(fn_tpm,
                 header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                 row.names = 1)
m=match.idx( sample_ids, names(tpm) )
tpm = data.matrix(tpm)
matrix_tpm = matrix( NA, nrow=dim(tpm)[1], ncol=length(sample_ids) )
matrix_tpm[,m$idx.A] = tpm[,m$idx.B ]
matrix_tpm = data.frame(matrix_tpm)
names(matrix_tpm) = sample_ids
rownames(matrix_tpm) = dimnames(tpm)[[1]]

################################################################################
# sequencing summary data
################################################################################

matrix_summary_normal = read.table( fn_alignment_summary_t,
                                    row.names=1, stringsAsFactors=FALSE, sep='\t', 
                                    header=TRUE)
matrix_summary_tumor = read.table( fn_alignment_summary_n,
                                   row.names=1, stringsAsFactors=FALSE, sep='\t', 
                                   header=TRUE)
coverage = data.frame(rbind( matrix_summary_normal["Mean_coverage",], 
                             matrix_summary_tumor["Mean_coverage",] ))
rownames(coverage) = c("normal","tumor")
write.table( coverage, '/notebook/human_sequence_prostate_WCDT/prepub_WCDT/results/coverage_summary.txt',
             sep='\t', quote=FALSE)

################################################################################
# Gene activation/ inactivation
################################################################################

matrix_somatic = matrix_inactive != '.' | matrix_missense != '.' 
sums_inactive = rowSums( matrix_inactive != '.' )
sums_missense = rowSums( matrix_missense != '.' )
sums_somatic = rowSums( matrix_missense != '.' | matrix_inactive != '.')

genes = rownames(matrix_CNA_int_ploidy)
n_genes = length(genes)
sexchrom = which( gene_locs$chrom == "chrX" | gene_locs$chrom=="chrY" )

curated_fs = read.table(fn_curated_fs,header=TRUE,sep='\t',stringsAsFactors=FALSE)
symbols_curated_fs = unique( curated_fs$symbol )
curated_sv = read.table(fn_curated_sv,header=TRUE,sep='\t',stringsAsFactors=FALSE)
curated_sv = curated_sv[curated_sv$reported=="yes",]
whitelist_SV = c( unique(curated_sv$threeprime), "AR" )
curated_missense = read.table(fn_curated_missense,header=TRUE,sep='\t',stringsAsFactors=FALSE)
symbols_curated_missense = unique( curated_missense$gene )
curated_missense = curated_missense[curated_missense$reported_pathogenic!="no",]


INACTIVE = matrix(0, nrow=length(genes), ncol=3)
ALLELES_INACTIVATED = matrix(0, nrow=length(genes), ncol=length(sample_ids))
matrix_SV = matrix('.', nrow=length(genes), ncol=length(sample_ids))
dimnames(ALLELES_INACTIVATED)[[1]] = genes
dimnames(ALLELES_INACTIVATED)[[2]] = sample_ids
dimnames(INACTIVE)[[1]] = genes
dimnames(INACTIVE)[[2]] = c("bi","mono", "none")
dimnames(matrix_SV)[[1]] = genes
dimnames(matrix_SV)[[2]] = sample_ids

# samples with tandem duplication hotspots identified in figure 2
matrix_SV["AR", assess_tandem_dup( chrom="chrX", locus=6692, 
                                   list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'
matrix_SV["FOXA1", assess_tandem_dup( chrom="chr14", locus=3754, 
                                      list_sv_m, chrom_lengths, sample_ids) ] = 'tandem'
matrix_SV["MYC", assess_tandem_dup( chrom="chr8", locus=12701, 
                                    list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'
matrix_SV["MYC", assess_tandem_dup( chrom="chr8", locus=12741, 
                                    list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'

# samples with curated SV
for(idx in 1:dim(curated_sv)[1]){
    sample_id = curated_sv$sample[ idx ]
    symbol = curated_sv$threeprime[ idx ]
    idx_in_samples = which( sample_ids==sample_id ) 
    consequence = curated_sv$consequence[ idx ]
    mechanism = curated_sv$mechanism[ idx ]
    matrix_SV[ symbol, sample_id ] = mechanism
}

matrix_CNA["CDK12", "DTB-183-BL"] = "LOH"

# Allele effects on expression

cors = rep(NA, dim(gene_locs)[1])
pvals = rep(NA, dim(gene_locs)[1])
n_0 = rep(NA, dim(gene_locs)[1])
n_1 = rep(NA, dim(gene_locs)[1])
n_2 = rep(NA, dim(gene_locs)[1])

s5 = read.table(fn_symbols_fig5,
                header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)
symbols = rownames(s5)
for(i in 1:length(symbols)){
    if( i%%100==0){
        print(i)
    }
    symbol=symbols[i]
    if( symbol != "AR enhancer"){
        ae=allele_effect( symbol, do_plot=FALSE, axis_override=NULL )
        INACTIVE[symbol, 1] = ae$n_bi
        INACTIVE[symbol, 2] = ae$n_mono
        INACTIVE[symbol, 3] = N_SAMPLES-ae$n_bi-ae$n_mono
        ALLELES_INACTIVATED[symbol,] = ae$alleles$n_alleles_inactivated
        idx = which(rownames(gene_locs)==symbols[i])
        cors[idx] = ae$cor_alleles_xp
        pvals[idx] = ae$pval_alleles_xp
        n_0[idx] = sum(ae$alleles$n_alleles_inactivated==2)
        n_1[idx] = sum(ae$alleles$n_alleles_inactivated==1)
        n_2[idx] = sum(ae$alleles$n_alleles_inactivated==0)
    }
}

matrix_samples$has_HRD = allele_effect("BRCA2")$alleles$n_alleles_inactivated==2
matrix_samples$has_ETS = allele_effect("ERG")$alleles$activating_sv |
                         allele_effect("ETV1")$alleles$activating_sv |
                         allele_effect("ETV4")$alleles$activating_sv |
                         allele_effect("ETV5")$alleles$activating_sv 


non_sex = cnbin$chrom!="chrX" & cnbin$chrom!="chrY"
percent_CNA = round( 
    (colSums(cnbin[non_sex,]>=GAIN_NONSEX|cnbin[non_sex,]<=LOSS_SINGLE_NONSEX, na.rm=TRUE) + 
    colSums(cnbin[!non_sex,]>=GAIN_SEX|cnbin[!non_sex,]<=GAIN_SEX,na.rm=TRUE) ) / 1042 , 2)[3:103]
matrix_samples$percent_CNA = percent_CNA

write.table( matrix_samples, '/notebook/human_sequence_prostate_WCDT/drafts/tables/S2_sample_information.txt', sep='\t', quote=FALSE)

################################################################################
# Chromoplexy Analysis
################################################################################

matrix_chain = read.table(fn_chain, header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)
