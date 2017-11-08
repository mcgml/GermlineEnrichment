#Description: Rscript to run ExomeDepth on a set of BAM files
#Author: Matt Lyon
#Mode: BY_COHORT

#load modules
require(ExomeDepth)
require(optparse)
require(Rsamtools)

#parse command line args
script.help <- "Args -b <bamlist> -r <bed> -f <fasta>"
option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', type='character', help="Path to list of BAMs"),
  make_option(c("-r", "--bed"), action="store", default='', type='character', help="Path to BED file"),
  make_option(c("-f", "--fasta"), action="store", default='', type='character', help="Path to FASTA")
)
opt = parse_args(OptionParser(option_list=option_list))

#check args are provided
if (opt$bamlist == "" || opt$bed == "" || opt$fasta == ""){
  print(script.help)
  stop()
}

# get bamlist
bams = read.table(opt$bamlist,header=FALSE,stringsAsFactors=FALSE)$V1

#load FASTA index
fasta.file <- FaFile(opt$fasta)

# get bam counts
ExomeCount = getBamCounts(bed.file = opt$bed,bam.files = bams,min.mapq = 20,include.chr = FALSE,referenceFasta = opt$fasta)

# prepare the main matrix of read count data
ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
nsamples <- ncol(ExomeCount.mat)

# start looping over each sample
for (i in 1:nsamples) {

  #extract run info from filename
  originalFilename <- gsub("^X", "", gsub("-bam", "\\.bam", gsub("\\.", "-", colnames(ExomeCount.mat)[i])))
  sampleName <- gsub("\\.bam", "", gsub(".*_", "", originalFilename))
  outputPrefix <- gsub("\\.bam", "", originalFilename)
  
  print(paste("Processing sample", i, ":", sampleName, sep = " ", collapse = NULL))
  
  # Create the aggregate reference set for this sample
  my.choice <- select.reference.set(test.counts = ExomeCount.mat[,i],reference.counts = ExomeCount.mat[,-i],bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,n.bins.reduced = 0)

  #print reference sample(s)
  print(paste("Reference sample(s):", my.choice$reference.choice, sep = " ", collapse = NULL))

  #create reference set
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],MAR = 1,FUN = sum)
  
  # Now creating the ExomeDepth object
  all.exons <- new('ExomeDepth',test = ExomeCount.mat[,i], reference = my.reference.selected,formula = 'cbind(test, reference) ~ 1')
  
  # Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,transition.probability = 0.05,chromosome = ExomeCount.dafr$space,start = ExomeCount.dafr$start,end = ExomeCount.dafr$end,name = ExomeCount.dafr$names)

  #write header
  vcf.filename <- paste(c(outputPrefix, "_cnv.vcf"), collapse="")
  write(
    paste(
        "##fileformat=VCFv4.1\n",
        paste("##fileDate=", format(Sys.Date(), format="%Y%m%d"), "\n", collapse = NULL, sep = ""),
        paste("##source=ExomeDepth", packageDescription("ExomeDepth")$Version, "\n", collapse = NULL, sep = ""),
        paste("##reference=file://", opt$fasta, "\n", collapse = NULL, sep=""),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##FORMAT=<ID=BF,Number=1,Type=Float,Description=\"Bayes factor\">\n",
        "##FORMAT=<ID=RE,Number=1,Type=Integer,Description=\"Reads expected\">\n",
        "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reads observed\">\n",
        "##FORMAT=<ID=RR,Number=1,Type=Float,Description=\"Reads ratio\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n",
        "##INFO=<ID=Regions,Number=1,Type=Integer,Description=\"Number of affected regions of interest\">\n",
        "##INFO=<ID=StartPosition,Number=1,Type=Integer,Description=\"Position of first affected region\">\n",
        "##INFO=<ID=EndPosition,Number=1,Type=Integer,Description=\"Position of last affected region\">\n",
        "##ALT=<ID=DEL,Description=\"Deletion\">\n",
        "##ALT=<ID=DUP,Description=\"Duplication\">\n",
        paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sampleName, sep="\t", collapse=NULL),
      collapse = NULL, sep = ""), vcf.filename)

  print(paste(length(all.exons@CNV.calls), "calls identified for", i, ":", sampleName, sep = " ", collapse = NULL))

  #check calls are available for printing
  if (length(all.exons@CNV.calls) > 0){
    
    #add reference base allele to each CNV call for VCF output later
    fasta.ranges <- GRanges(seqnames=all.exons@CNV.calls$chromosome, ranges=IRanges(all.exons@CNV.calls$start, all.exons@CNV.calls$start))
    fasta.sequence <- getSeq(fasta.file, fasta.ranges)
    all.exons@CNV.calls$ref.allele <- as.character(fasta.sequence)

    #genotype using reads.ratio
    all.exons@CNV.calls$genotypes <- "0/1"
    all.exons@CNV.calls$genotypes[which(all.exons@CNV.calls$reads.ratio < 0.25)] <- "1/1"

    #write results to VCF file
    write(
      paste(
        all.exons@CNV.calls$chromosome,all.exons@CNV.calls$start,".",all.exons@CNV.calls$ref.allele,gsub("duplication","<DUP>",gsub("deletion", "<DEL>", all.exons@CNV.calls$type)),".","PASS",
        paste("END=",all.exons@CNV.calls$end,";","Regions=",all.exons@CNV.calls$nexons,";","StartPosition=",all.exons@CNV.calls$start.p,";","EndPosition=", all.exons@CNV.calls$end.p, sep="", collapse=NULL), 
        "GT:BF:RE:RO:RR",
        paste(all.exons@CNV.calls$genotypes, all.exons@CNV.calls$BF, all.exons@CNV.calls$reads.expected, all.exons@CNV.calls$reads.observed, all.exons@CNV.calls$reads.ratio, collapse=NULL, sep=":"),
        sep="\t"), vcf.filename, append=TRUE
    )
  }
    
}

#print session info for logging purposes
sessionInfo()

#print warnings
if (!is.null(warnings())){
  warnings()
}