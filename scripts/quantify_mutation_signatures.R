suppressPackageStartupMessages(
    require(deconstructSigs, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)
option_list <- list(
	make_option(c("-s", "--sample_id"), type = "character", help = "Sample identifier. [Required]"),
	make_option(c("-i", "--fn_CPRA"),   type = "character", help = "File containing CHROM,POS,REF,ALT. [Required]"),
	make_option(c("-o", "--out"),       type = "character", help = "file name with signature weights. [Required]")
)
parseobj = OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt = parse_args(parseobj)

sample_id = opt$sample_id 
fn_CPRA = opt$fn_CPRA 
fn_out = opt$out
message(paste("MESSAGE: parsing mutations in",fn_CPRA))
message(paste("MESSAGE: writing percentages to ", fn_out))
muts=read.table(fn_CPRA, stringsAsFactors=FALSE, header=TRUE)
muts=cbind(muts, sample=rep(sample_id, dim(muts)[1], stringsAsFactors=FALSE))
bad_contigs=setdiff( unique(muts$CHROM), names(BSgenome.Hsapiens.UCSC.hg38))
muts=muts[! muts$CHROM %in% bad_contigs,]


sig_input = mut.to.sigs.input( mut.ref=muts, sample.id="sample", chr="CHROM", 
                               pos="POS", ref="REF", alt="ALT", 
                               bsg=BSgenome.Hsapiens.UCSC.hg38)
samples = dimnames(sig_input)[[1]]
sig=whichSignatures(tumor.ref=sig_input, sample.id=samples[i],
  signatures.ref = signatures.cosmic, associated = c(),
  signatures.limit = NA, signature.cutoff = 0.06, contexts.needed = TRUE,
  tri.counts.method = "default")

write.table( 
  data.frame(
    signature_name=names(sig$weights),
    weight=round(as.numeric(sig$weights),3),
    stringsAsFactors=FALSE), fn_out,
  sep="\t", quote=FALSE, row.names=FALSE)

message("completed mutation signature analysis")