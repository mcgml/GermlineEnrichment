#Description: Rscript to run ExomeDepth on a set of BAM files
#Author: Matt Lyon
#Status: Development
#Mode: BY_COHORT

#load modules
require(ExomeDepth)
require(optparse)

#parse command line args
option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', type='character', help="Path to list of BAMs"),
  make_option(c("-r", "--bed"), action="store", default='', type='character', help="Path to BED file"),
  make_option(c("-f", "--fasta"), action="store", default='', type='character', help="Path to FASTA")
)
opt = parse_args(OptionParser(option_list=option_list))

# get bamlist
bams = read.table(opt$bamlist,
header=FALSE,
stringsAsFactors=FALSE)$V1

# get bam counts
ExomeCount = getBamCounts(bed.file = opt$bed,
bam.files = bams,
min.mapq = 20,
include.chr = FALSE,
referenceFasta = opt$fasta)

# prepare the main matrix of read count data
ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
nsamples <- ncol(ExomeCount.mat)

# start looping over each sample
for (i in 1:nsamples) {
  print(paste("Processing sample", i, colnames(ExomeCount.mat)[i], sep = " ", collapse = NULL))
  
  # Create the aggregate reference set for this sample
  my.choice <- select.reference.set(test.counts = ExomeCount.mat[,i],
  reference.counts = ExomeCount.mat[,-i],
  bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
  n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
  MAR = 1,
  FUN = sum)
  
  # Now creating the ExomeDepth object
  all.exons <- new('ExomeDepth',
  test = ExomeCount.mat[,i], 
  reference = my.reference.selected,
  formula = 'cbind(test, reference) ~ 1')
  
  # Now call the CNVs
  all.exons <- CallCNVs(x = all.exons,
  transition.probability = 10^-2,
  chromosome = ExomeCount.dafr$space,
  start = ExomeCount.dafr$start,
  end = ExomeCount.dafr$end,
  name = ExomeCount.dafr$names)
  
  #write to file
  write.table(file = paste(colnames(ExomeCount.mat)[i], 'txt', sep = '.'), x = all.exons@CNV.calls, row.names = FALSE, quote = FALSE, sep = "\t")
  
}

#print session info for logging
sessionInfo()

#print warnings
warnings()
