#!/usr/bin/env Rscript

library(ape)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(gplots)
library(optparse)
library(poppr)
library(UpSetR)
library(vcfR)
library("tibble")
 
getArgs <- function() {
    option_list = list(
        make_option(c("-t", "--genotype_file"), default=NULL, help="genotype file path"),
        make_option(c("-v", "--vcf_file"), default=NULL, help="vcf zip file path"),
        make_option(c("-f", "--fasta_file"), default=NULL, help="fasta file path"),
        make_option(c("-g", "--gff_file"), default=NULL, help="gff3 file path"),
        make_option("--mUpset", default=NULL, metavar="upset_mat_file",
              help="mat file path_for_upset"),
        make_option("--mInter", default=NULL, metavar="intersect_mat_file",
              help="intersect mat file path_for_heatMap"),
        make_option("--mJacc", default=NULL, metavar="jacc_mat_file",
              help="jaccard mat file path_for_heatMap")
    );

    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);
    if ( is.null(opt$genotype_file) | is.null(opt$vcf_file) | is.null(opt$fasta_file) | is.null(opt$gff_file) | is.null(opt$mUpset) | is.null(opt$mInter) |is.null(opt$mJacc) 
       )
    {
         print_help(opt_parser)
         stop("All arguements need to be supplied.\n", call.=FALSE)
    }
    assign("genotype_file", opt$genotype_file, envir = .GlobalEnv)
    assign("vcf_file", opt$vcf_file, envir = .GlobalEnv)
    assign("fasta_file", opt$fasta_file, envir = .GlobalEnv)
    assign("gff_file", opt$gff_file, envir = .GlobalEnv)
    assign("upset_mat_file", opt$mUpset, envir = .GlobalEnv)
    assign("intersect_mat_file", opt$mInter, envir = .GlobalEnv)
    assign("jaccard_mat_file", opt$mJacc, envir = .GlobalEnv)
}

build_genlight <- function() {
    sample.data <- read.table(genotype_file, sep = "\t", header = FALSE, stringsAsFactors=F)
    vcf         <- read.vcfR(vcf_file, verbose = FALSE)
    dna         <- ape::read.dna(fasta_file, format = "fasta")
    gff         <- read.table(gff_file, sep="\t", quote="")
    colnames(vcf@gt) <- gsub('_','',colnames(vcf@gt), perl=TRUE)
    colnames(vcf@gt)[-1] <-  sample.data$V2 # add column name
    # Subset run_id with variants call sets
    gl.bowtie2 <<- vcfR2genlight(vcf)
    pop(gl.bowtie2) <- sample.data$V3
    gl.bowtie2 <<- gl.bowtie2
    ploidy(gl.bowtie2) <- 1
    if (nPop(gl.bowtie2)>2) {
        cols <<- brewer.pal(n = nPop(gl.bowtie2), name = "Paired")
    }
    else {
        cols <<- c('blue','red')
    }
}

build_phylogenetic_tree <- function() {
    print("building phylogenetic tree")
    tree.bowtie2 <- aboot(gl.bowtie2, tree = "upgma", distance = "bitwise.dist", sample = 100, showtree = F, cutoff = 0, quiet = FALSE, root = T)
    print(tree.bowtie2)
    plot.phylo(tree.bowtie2, cex = 0.8, font = 2, adj = 0, tip.color = cols[pop(gl.bowtie2)])
    nodelabels(tree.bowtie2$node.label, adj = c(1.3, -0.5), cex = 0.5, frame = "n", font = 2, xpd = TRUE)
    axis(side = 1)
    title(xlab = "Genetic distance (proportion of loci that are different)")
    title(main = "SNPs-based Isolates phylogeny")
} 

build_PCA <- function() {
    print("building PCA")
    bowtie2.pca <- glPca(gl.bowtie2, nf = 3)
    barplot(100*bowtie2.pca$eig/sum(bowtie2.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
    title(ylab="Percent of variance\nexplained", line = 2)
    title(xlab="Eigenvalues", line = 1)
    bowtie2.pca.scores <- as.data.frame(bowtie2.pca$scores)
    set.seed(9)
    p <- ggplot(bowtie2.pca.scores, aes(x=PC1, y=PC2, colour = cols[pop(gl.bowtie2)]))
    p <- p + geom_point(size=2) + geom_label_repel(aes(label = rownames(bowtie2.pca.scores)))
    p <- p + scale_color_manual(values = cols) 
    p <- p + geom_hline(yintercept = 0) 
    p <- p + geom_vline(xintercept = 0) 
    p <- p + ggtitle("PCA of observed SNPs") +theme_bw() + 
         theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    p
}

build_heatMap <- function() {
    print("building heatMap")
    input_table            <- read.table(intersect_mat_file, header=TRUE)
    row.names(input_table) <- input_table$name
    input_table            <- input_table[, -1]
    input_matrix           <- as.matrix(input_table)
    heatmap.2(input_matrix,
          col=brewer.pal(9,"Oranges"),
          margins = c(14, 14),
          density.info = "none",
          lhei = c(2, 8),
          trace="none",
          srtCol=90)
}

build_UpSetPlot <- function() {
    print("building UpSet plot")
    snp           <- read.delim(file = upset_mat_file)
    sample_name   <- colnames(snp)[-c(1:5)]
    colnames(snp) <- gsub('.PASS.vcf','',colnames(snp))
    sample_name   <- gsub('.PASS.vcf','',sample_name)
    snp$chrom     <- paste(snp$chrom,snp$start,sep="-")
    upset(snp,
          sets=sample_name, 
          mainbar.y.label="Intersection size",
          sets.x.label="set size",
          #number.angles=15,
          #text.scale = 0.7, 
          mb.ratio=c(0.55, 0.45),
          order.by = "freq",
          #keep.order = T
          sets.bar.color = "#56B4E9",
          empty.intersections = NULL, 
          nintersects = 50, 
          number.angles = 0, 
          set_size.angles=0,
          query.legend = 'top', 
          text.scale = .7
         )
}
#    upset(snp,sets=sample_name, mainbar.y.label="Intersection size",
#          sets.x.label="set size",
#          number.angles=15, text.scale = 0.7, mb.ratio=c(0.55, 0.45),
#          order.by="freq", keep.order=TRUE)

#upset(snps, sets =isolate , sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = NULL, nintersects = 50, number.angles = 0, set_size.angles=0,  query.legend = 'top', set_size.show = T, text.scale = .75) 

getArgs()
build_genlight()
build_phylogenetic_tree()
build_PCA()
build_heatMap()
build_UpSetPlot()
