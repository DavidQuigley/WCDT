match.idx=function(A, B, allow.multiple.B=F){
    # return dataframe of indices into A and B restricted to perfect matches
    # between A and B, where idx.A[i] == idx.B[i] for each i in matched pairs
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}
suppressPackageStartupMessages(
    require(deconstructSigs, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(SomaticSignatures, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(VariantAnnotation, quietly=TRUE, warn.conflicts=FALSE)
)
option_list <- list(
	make_option(c("-a", "--sample_attributes"), type = "character", help = "File containing sample attributes [Required]"),
	make_option(c("-d", "--dir_in"),   type = "character", help = "directory containing VCFs [Required]"),	
	make_option(c("-s", "--suffix"), type = "character", help = "suffix for VCF files to load [Required]"),
	make_option(c("-c", "--suffix_CPRA"), type = "character", help = "suffix for CPRA files to load [Required]"),
	make_option(c("-o", "--dir_out"),       type = "character", help = "directory to which files will be written. [Required]")
)

genome = BSgenome.Hsapiens.UCSC.hg38

parseobj = OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt = parse_args(parseobj)

fn_sa = opt$sample_attributes
dir_in = opt$dir_in
suffix = opt$suffix
dir_out = opt$dir_out
suffix_CPRA = opt$suffix_CPRA

sa = read.table(fn_sa, stringsAsFactors=FALSE, header=TRUE, sep='\t', row.names=1)
mm = NA
for( i in 1:dim(sa)[1]){
    fn=paste(rownames(sa)[i], suffix, sep="")
    print( rownames(sa)[i] )
    vcf = readVcf( paste(dir_in, fn, sep='/'), "hg38" )
    vr = VRanges(
        seqnames = seqnames(vcf),
        ranges = ranges(vcf),
        ref=as.character(mcols(vcf)$REF),
        alt=as.character(unlist(mcols(vcf)$ALT)),
        sampleNames = as.character(samples(header(vcf))),
        seqinfo = seqinfo(vcf)
    )
    motifs = mutationContext(vr, ref=genome)
    mm_cur = motifMatrix(motifs, group = "sampleNames", normalize = TRUE) 
    if( is.na( mm ) ){
        mm = mm_cur   
    }else{
        mm = cbind( mm, mm_cur )
    }
}
m = match.idx( dimnames(mm)[[2]], sa$ID_illumina_tumor)
dimnames(mm)[[2]] = paste( "DTB-", get.split.col( rownames(sa)[m$idx.B], "-", col=2), sep="" )

gof_nmf = assessNumberSignatures(mm, 2:16, nReplicates = 5)
sigs_nmf = identifySignatures(mm, 14, nmfDecomposition)
sig_denovo=signatures(sigs_nmf)

fn_out_sig = paste( dir_out, 'matrix_de_novo_mutsig.txt', sep='/')
write.table( round(data.frame(t(sig_denovo)), 3), fn_out_sig,
             quote = FALSE, row.names = TRUE, col.names = TRUE)
fn_out_sig = paste( dir_out, 'matrix_de_novo_mutsigs_by_samples.txt', sep='/')
write.table( signif( samples(sigs_nmf),3), fn_out_sig,
             quote = FALSE, row.names = TRUE, col.names = TRUE)

sig.for.de = data.frame(t(sig_denovo))
names(sig.for.de)= names(signatures.nature2013)
for( i in 1:dim(sa)[1]){
    fn_CPRA=paste(dir_in, '/', rownames(sa)[i], suffix_CPRA, sep="")
    muts=read.table(fn_CPRA, stringsAsFactors=FALSE, header=TRUE)
    muts=cbind(muts, sample=rep(sample_id, dim(muts)[1], stringsAsFactors=FALSE))
    bad_contigs=setdiff( unique(muts$CHROM), names(BSgenome.Hsapiens.UCSC.hg38))
    muts=muts[! muts$CHROM %in% bad_contigs,]
    sig_input = mut.to.sigs.input( mut.ref=muts, sample.id="sample", chr="CHROM", 
                               pos="POS", ref="REF", alt="ALT", 
                               bsg=BSgenome.Hsapiens.UCSC.hg38)

    sig=whichSignatures(tumor.ref=sig_input, sample.id=sample_id,
      signatures.ref = sig.for.de, associated = c(),
      signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
      tri.counts.method = "default")
    
    new_weights = data.frame(
                        signature_name=names(sig$weights),
                        weight=round(as.numeric(sig$weights),3),
                        stringsAsFactors=FALSE)
    if(i==1){
        weights = new_weights
    }else{
        weights = cbind( weights, new_weights, stringsAsFactors=FALSE)
    }
    rownames(weights) = weights[,1]
    weights = weights[,seq(from=2, to=208,by=2)]
    names(weights) = rownames(sa)
}
fn_out_applied = paste( dir_out, 'matrix_de_novo_weights_applied_to_deconstructSigs.txt', sep='/')
write.table( weights, fn_out_applied,  sep="\t", quote=FALSE, row.names=FALSE)

