suppressPackageStartupMessages(
    require(copyCat, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)
option_list <- list(
	make_option(c("-s", "--sample_id"), type = "character", help = "sample_id. [Required]"),
	make_option(c("-n", "--fn_normal"), type = "character", help = "normal BAM window. [Required]"),
	make_option(c("-t", "--fn_tumor"), type = "character", help = "tumor BAM window. [Required]"),
	make_option(c("-o", "--dir_out"), type = "character", help = "output folder. [Required]")
)
parseobj = OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt = parse_args(parseobj)
sample_id=opt$sample_id
dir_annot='/reference/UCSC/GRCh38Decoy/Annotation/copycat'
dir_output= paste(opt$dir_out, '/', sample_id, sep='')
path_window_n=opt$fn_normal
path_window_t=opt$fn_tumor
runPairedSampleAnalysis(annotationDirectory=dir_annot,
                    outputDirectory=dir_output,
                    normal=path_window_n,
                    tumor=path_window_t,
                    inputType="bins",
                    maxCores=2,
                    binSize=0, #infer automatically from bam-window output
                    readLength=150,
                    perLibrary=FALSE,
                    perReadLength=FALSE, #correct each read-length independently
                    verbose=TRUE,
                    minWidth=3, #minimum number of consecutive winds need to call CN
                    minMapability=0.6, #a good default
                    dumpBins=TRUE,
                    doGcCorrection=TRUE,
                    samtoolsFileFormat=NULL, #will infer automatically - mpileup 10col or VCF
                    purity=1,
                    normalSamtoolsFile=NULL,
                    tumorSamtoolsFile=NULL)  #uses the VAFs of mpileup SNPs to infer copy-neutral regions
