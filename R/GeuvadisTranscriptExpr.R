################################################################################
### Data documentation
################################################################################



#' Sample data for sQTL analysis
#' 
#' Subsets of raw data available in this package and saved as Rdata objects for
#' faster loading.
#' 
#' @format \code{counts} is a data frame with subset of counts from
#' TrQuantCount_CEU_chr19.tsv
#' 
#' \code{gene_ranges} is a GRanges object containing subset of gene coordinates
#' from genes_chr19.bed
#' 
#' \code{genotypes} is a data frame with subset of genotypes from
#' genotypes_CEU_chr19.tsv
#' 
#' \code{snp_ranges} is a Granges object containing subset of SNP coordinates
#' from genotypes_CEU_chr19.tsv
#' 
#' For all the details on how these data sets were produced, see examples.
#' 
#' @source Lappalainen T, Sammeth M, Friedlander MR, et al. Transcriptome and
#' genome sequencing uncovers functional variation in humans. Nature.
#' 2013;501(7468):506-11
#' 
#' @return \code{counts}, \code{gene_ranges}, \code{genotypes},
#' \code{snp_ranges}
#' 
#' @examples
#' library(rtracklayer)
#' data_dir <- system.file("extdata", package = "GeuvadisTranscriptExpr")
#' 
#' gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
#' snp_id_subset <- readLines(file.path(data_dir, "snp_id_subset.txt"))
#' 
#' # Load gene ranges with names!
#' gene_ranges <- import(file.path(data_dir, "genes_chr19.bed"))
#' names(gene_ranges) <- mcols(gene_ranges)$name
#' 
#' gene_ranges <- gene_ranges[gene_id_subset, ]
#' 
#' # Load transcript counts
#' counts <- read.table(file.path(data_dir, "TrQuantCount_CEU_chr19.tsv"),
#'                      header = TRUE, sep = "\t", as.is = TRUE)
#' 
#' counts <- counts[counts$Gene_Symbol %in% gene_id_subset, ]
#' 
#' # Load genotypes
#' genotypes <- read.table(file.path(data_dir, "genotypes_CEU_chr19.tsv"),
#'                         header = TRUE, sep = "\t", as.is = TRUE)
#' 
#' genotypes <- genotypes[genotypes$snpId %in% snp_id_subset, ]
#' 
#' # Create SNP ranges with names!
#' snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start,
#'                                                   genotypes$end))
#' names(snp_ranges) <- genotypes$snpId
#' 
"counts"



#' @rdname counts
"gene_ranges"

#' @rdname counts
"genotypes"

#' @rdname counts
"snp_ranges"















