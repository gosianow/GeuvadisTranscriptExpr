
setwd("/home/gosia/R/multinomial_project/package_devel/GeuvadisTranscriptExpr/data-raw/")


########################################################
# process data downloaded from the GEUVADIS website
########################################################

library(GenomicRanges)
library(rtracklayer)
library(limma)

### annotation

gtf_dir <- "gencode.v12.annotation.gtf"

gtf <- import(gtf_dir)

# keep protein coding genes
keep_index <- mcols(gtf)$gene_type == "protein_coding" & mcols(gtf)$type == "gene"
gtf <- gtf[keep_index]
# remove 'chr'
seqlevels(gtf) <- gsub(pattern = "chr", replacement = "", x = seqlevels(gtf))

genes_bed <- data.frame(chr = seqnames(gtf), start =  start(gtf), end = end(gtf), geneId = mcols(gtf)$gene_id, stringsAsFactors = FALSE)


for(i in as.character(1:22)){

  genes_bed_sub <- genes_bed[genes_bed$chr == i, ]

  write.table(genes_bed_sub, "genes_chr", i ,".bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  }



### samples

metadata_dir <- "E-GEUV-1.sdrf.txt"

samples <- read.table(metadata_dir, header = TRUE, sep="\t", as.is=TRUE)

samples <- unique(samples[c("Assay.Name", "Characteristics.population.")])
colnames(samples) <- c("sample_id", "population")

samples$sample_id_short <- strsplit2(samples$sample_id, "\\.")[,1]


### expression in counts

expr_dir <- "GD660.TrQuantCount.txt"

expr_all <- read.table(expr_dir, header = TRUE, sep="\t", as.is = TRUE)

expr_all <- expr_all[, c("TargetID", "Gene_Symbol", "Chr", samples$sample_id)]
colnames(expr_all) <- c("TargetID", "Gene_Symbol", "Chr", samples$sample_id_short)


for(j in "CEU"){

  for(i in 1:22){

    expr <- expr_all[expr_all$Chr == i, c("TargetID", "Gene_Symbol", samples$sample_id_short[samples$population == j])]

    write.table(expr, paste0("TrQuantCount_", j, "_chr", i, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  }

}




### genotypes
library(Rsamtools)
library(VariantAnnotation)
library(tools)


files <- list.files(path = ".", pattern = "genotypes.vcf.gz", full.names = TRUE, include.dirs = FALSE)

## bigzip and index the vcf files
for(i in 1:length(files)){

  zipped <- bgzip(files[i])
  idx <- indexTabix(zipped, format = "vcf")

}

## extended gene ranges
window <- 5000
gene_ranges <- resize(gtf, GenomicRanges::width(gtf) + 2 * window, fix = "center")

chr <- gsub("chr", "", strsplit2(files, split = "\\.")[, 2])


for(j in "CEU"){

  for(i in 1:length(files)){

    cat(j, chr[i], fill = TRUE)

    zipped <- paste0(file_path_sans_ext(files[i]), ".bgz")
    idx <- paste0(file_path_sans_ext(files[i]), ".bgz.tbi")
    tab <- TabixFile(zipped, idx)

    ## Explore the file header with scanVcfHeader
    hdr <- scanVcfHeader(tab)
    print(all(samples$sample_id_short %in% samples(hdr)))

    ## Read VCF file
    gene_ranges_tmp <- gene_ranges[seqnames(gene_ranges) == chr[i]]

    param <- ScanVcfParam(which = gene_ranges_tmp, samples = samples$sample_id_short[samples$population == j])
    vcf <- readVcf(tab, "hg19", param)


    ## Keep only the bi-allelic SNPs

    # width of ref seq
    rw <- width(ref(vcf))
    # width of first alt seq
    aw <- unlist(lapply(alt(vcf), function(x) {width(x[1])}))
    # number of alternate genotypes
    nalt <- elementLengths(alt(vcf))
    # select only bi-allelic SNPs (monomorphic OK, so aw can be 0 or 1)
    snp <- rw == 1 & aw <= 1 & nalt == 1
    # subset vcf
    vcfbi <- vcf[snp,]

    rowdata <- rowData(vcfbi)

    ## Convert genotype into number of alleles different from reference
    geno <- geno(vcfbi)$GT
    geno01 <- geno
    geno01[,] <- -1
    geno01[geno %in% c("0/0", "0|0")] <- 0 # REF/REF
    geno01[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # REF/ALT
    geno01[geno %in% c("1/1", "1|1")] <- 2 # ALT/ALT
    # geno01 should be integer, not character
    mode(geno01) <- "integer"

    genotypes <- unique(data.frame(chr = seqnames(rowdata), start = start(rowdata), end = end(rowdata), snpId = rownames(geno01), geno01, stringsAsFactors = FALSE))

    ### sorting
    genotypes <- genotypes[order(genotypes[ ,2]), ]

    write.table(genotypes, file = paste0("genotypes_", j, "_chr", chr[i], ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


  }
}





########################################################
# Save data as Rdata files
########################################################

setwd("/home/gosia/R/multinomial_project/package_devel/GeuvadisTranscriptExpr")

library(devtools)


data_dir <- system.file("extdata", package = "GeuvadisTranscriptExpr")


gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
snp_id_subset <- readLines(file.path(data_dir, "snp_id_subset.txt"))



# Load gene ranges with names!
gene_ranges <- import(file.path(data_dir, "genes_chr19.bed"))
names(gene_ranges) <- mcols(gene_ranges)$name

gene_ranges <- gene_ranges[gene_id_subset, ]

use_data(gene_ranges, overwrite = TRUE)


# Load transcript counts
counts <- read.table(file.path(data_dir, "TrQuantCount_CEU_chr19.tsv"),
                     header = TRUE, sep = "\t", as.is = TRUE)

counts <- counts[counts$Gene_Symbol %in% gene_id_subset, ]

use_data(counts, overwrite = TRUE)


# Load genotypes
genotypes <- read.table(file.path(data_dir, "genotypes_CEU_chr19.tsv"),
                        header = TRUE, sep = "\t", as.is = TRUE)

genotypes <- genotypes[genotypes$snpId %in% snp_id_subset, ]

use_data(genotypes, overwrite = TRUE)


# Create SNP ranges with names!
snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start,
                                                  genotypes$end))
names(snp_ranges) <- genotypes$snpId


use_data(snp_ranges, overwrite = TRUE)




########################################################
# Choose a subset of genes and SNPs
########################################################

library(DRIMSeq)

setwd("/home/gosia/R/multinomial_project/package_devel/GeuvadisTranscriptExpr/data-raw/")

data_dir <- system.file("extdata", package = "GeuvadisTranscriptExpr")


### Create dmSQTLdata object


# gene_ranges with names!
gene_ranges <- import(file.path(data_dir, "genes_chr19.bed"))
names(gene_ranges) <- mcols(gene_ranges)$name

counts <- read.table(file.path(data_dir, "TrQuantCount_CEU_chr19.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

genotypes <- read.table(file.path(data_dir, "genotypes_CEU_chr19.tsv"), header = TRUE, sep = "\t", as.is = TRUE)

# snp_ranges with names!
snp_ranges <- GRanges(Rle(genotypes$chr), IRanges(genotypes$start, genotypes$end))
names(snp_ranges) <- genotypes$snpId

## Check if samples in count and genotypes are in the same order
all(colnames(counts[, -(1:2)]) == colnames(genotypes[, -(1:4)]))
sample_id <- colnames(counts[, -(1:2)])


d <- d_org <- dmSQTLdataFromRanges(counts = counts[, -(1:2)], gene_id = counts$Gene_Symbol, feature_id = counts$TargetID, gene_ranges = gene_ranges, genotypes = genotypes[, -(1:4)], snp_id = genotypes$snpId, snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = 2))

plotData(d, out_dir = "./")



### Filtering
d <- dmFilter(d_org, min_samps_gene_expr = 70, min_samps_feature_expr = 5, min_samps_feature_prop = 5, minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0.01, max_features = Inf, BPPARAM = BiocParallel::MulticoreParam(workers = 2))


plotData(d, out_dir = "./filtering_")


# set.seed(123)
# genes_subset <- names(d)[sample(length(d), size = 50, replace = FALSE)]


oo <- order(width(d@genotypes), decreasing = FALSE)
genes_subset <- names(d)[oo][1:50]



d_sub1 <- d_org[genes_subset, ]
snps_subset1 <- unique(d_sub1@blocks@unlistData[, "snp_id"])


d_sub2 <- d[genes_subset, ]
snps_subset2 <- unique(d_sub2@blocks@unlistData[, "snp_id"])


snps_extra <- setdiff(snps_subset1, snps_subset2)
snps_extra <- snps_extra[sample(length(snps_extra), size = 500, replace = FALSE)]


write.table(genes_subset, "gene_id_subset.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(c(snps_subset2, snps_extra), "snp_id_subset.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



























