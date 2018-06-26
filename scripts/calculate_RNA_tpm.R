library(matrixStats)

write.matrix = function(o, fn_out, na.string="NA"){
    # Write data frame in standard format used by load.matrix
    rownames = append( 'IDENTIFIER', names(o) )
    write( rownames, fn_out, sep='\t', ncolumns=length(rownames) )
    write.table(o, fn_out, quote=F, sep='\t', row.names=T, na=na.string, col.names=F, append=T)
}

match.idx = function(A, B, allow.multiple.B=F){
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

get.split.col=function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")
        
    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}

# Code from https://gist.github.com/slowkow/c6ab0348747f86e2748b
# 
#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
    
    # Ensure valid arguments.
    stopifnot(length(featureLength) == nrow(counts))
    stopifnot(length(meanFragmentLength) == ncol(counts))
    
    # Compute effective lengths of features in each library.
    effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        featureLength - meanFragmentLength[i] + 1
    }))
    
    # Exclude genes with length less than the mean fragment length.
    idx <- apply(effLen, 1, function(x) min(x) > 1)
    counts <- counts[idx,]
    effLen <- effLen[idx,]
    featureLength <- featureLength[idx]
    
    # Process one column at a time.
    tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
        rate = log(counts[,i]) - log(effLen[,i])
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }))
    
    # Copy the row and column names from the original matrix.
    colnames(tpm) <- colnames(counts)
    rownames(tpm) <- rownames(counts)
    return(tpm)
}



option_list = list(
)
parser = optparse::OptionParser(
  "Rscript calculate_RNA_tpm.R' [options] fn_refflat fn_metrics fn_counts fn_tpm",
  description=c("Extract SV calls.\n"),
  option_list=option_list
)
opt = optparse::parse_args(parser, positional_arguments=TRUE)
if (length(opt$args) < 4) {
    optparse::print_help(parser)
    write("Some of fn_refflat fn_metrics fn_counts fn_tpm are missing.\n", stderr())
    quit("no", 1)
}

#fn_refflat = '/notebook/human_sequence_prostate_WCDT/WCDT/metadata/GRCh38Decoy_refseq_genelocs_from_refFlat.bed'
#fn_metrics = paste(dir_results, '2018_01_02_matrix_rna_metrics.txt', sep='/')
#fn_counts = paste(dir_results, '2018_01_02_matrix_rna_counts.txt', sep='/')
#fn_tpm = paste( dir_results, '2018_01_02_matrix_rna_tpm.txt', sep='/')

#fn_refflat = opt$args[1]
fn_genelength = opt$args[1]
fn_metrics = opt$args[2]
fn_counts = opt$args[3]
fn_tpm = opt$args[4]

#print(paste("Reading refflat",fn_refflat))
print(paste("Reading gene length",fn_genelength))
print(paste("Reading metrics",fn_metrics))
print(paste("Reading counts",fn_counts))

metrics = read.table(fn_metrics,
                    row.names=1, header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                    check.names = FALSE)
insert_length = round( metrics$combined_MeanInsertLength, 1)

RNAcounts = read.table(fn_counts,
                     header=TRUE, sep='\t', 
                     stringsAsFactors=FALSE, check.names = FALSE)

# determine median gene length
symbols = sort(unique(RNAcounts$symbol))
feature_length = rep(0, length(symbols))
#refflat = read.table(fn_refflat, header=FALSE, sep='\t', stringsAsFactors=FALSE, 
#                     check.names = FALSE)
#names(refflat) = c("chrom", "start", "end", "name", "score", "strand")
#refflat = cbind(refflat, 
#                symbol=get.split.col( refflat$name, "~", first=TRUE),
#                length=refflat$end - refflat$start, 
#                stringsAsFactors=FALSE)
#gene_length = read.table( fn_genelength
#for(i in 1:length(symbols)){
#    idx = which( refflat$symbol==symbols[i] )
#    if( length( idx ) > 0 ){
#        feature_length[i] = median( refflat$length[ idx ] )     
#    }
#}
gene_length = read.table(fn_gene_length,
                     header=FALSE, sep='\t', 
                     stringsAsFactors=FALSE, check.names = FALSE)
gene_length = gene_length[,2:3]
names(gene_length) = c("symbol", "length")

for(i in 1:length(symbols)){
    idx = which( gene_length$symbol==symbols[i] )
    if( length( idx ) > 0 ){
        feature_length[i] = gene_length$length[idx]
    }
}

n_samples = dim(RNAcounts)[2]
samples = dimnames(RNAcounts)[[2]][2:n_samples]
m = match.idx( symbols, RNAcounts$symbol)
RNAcounts = RNAcounts[m$idx.B,]
counts = data.matrix( RNAcounts[, 2: dim(RNAcounts)[2] ] )
dimnames(counts)[[1]] = symbols
dimnames(counts)[[2]] = samples
mu.counts = rowMeans(counts)
max.counts = rowMaxs(counts)
keep = max.counts>100 | mu.counts > 100
#keep = rep(TRUE, dim(counts)[1])
tpm = counts_to_tpm(counts[keep,], feature_length[keep], insert_length) 

tpm_df = data.frame(round(tpm,3))
names(tpm_df) = dimnames(tpm)[[2]]
print(paste("Writing tpm to",fn_tpm))
write.matrix( tpm_df, fn_tpm )