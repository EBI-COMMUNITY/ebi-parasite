#!/usr/bin/env Rscript

library(vcfR)
#library(poppr)
library(ape)
#library(RColorBrewer)
library(optparse)

# Function to plot distribution of SNP 
# On crypto chromosomes.

getArgs <- function() {
    option_list = list(
        make_option(c("-s", "--seqid"), type="character", default=NULL, 
              help="seqid", metavar="character"),
        make_option(c("-v", "--vcf_file"), type="character", default=NULL, 
              help="vcf file path", metavar="character"),
        make_option(c("-f", "--fasta_file"), type="character", default=NULL,
              help="fasta file path", metavar="character"),
        make_option(c("-g", "--gff_file"), type="character", default=NULL,
              help="gff3 file path", metavar="character")
    ); 

    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$seqid) | is.null(opt$vcf_file) | is.null(opt$fasta_file) | is.null(opt$gff_file)){
         print_help(opt_parser)
         stop("All arguements need to be supplied.\n", call.=FALSE)
    }
    assign("seqid", opt$seqid, envir = .GlobalEnv)
    assign("vcf_file", opt$vcf_file, envir = .GlobalEnv)
    assign("fasta_file", opt$fasta_file, envir = .GlobalEnv)
    assign("gff_file", opt$gff_file, envir = .GlobalEnv)
    print ("here0")
}

distPlotter <- function(){
  vcf       <- read.vcfR(vcf_file, verbose = TRUE )
  fasta_dna <- ape::read.dna(fasta_file, format = "fasta")
  gff       <- read.table(gff_file, sep="\t", quote="")
  # Plot chromosomes
  chrome <- create.chromR(name=paste(seqid), vcf=vcf, seq=fasta_dna, ann=gff)
  chrome.nomasker <- proc.chromR(chrome, verbose=TRUE)
  plot(chrome.nomasker,main=seqid)
  return(chrome.nomasker)
}

getArgs()
chromoqc(distPlotter(), dp.alpha=40)


