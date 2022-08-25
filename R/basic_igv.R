#demo_available_syms = c("C20orf96", "ZCCHC3", "SOX12", "NRSN2", "NRSN2-AS1", "TRIB3", 
#"RBCK1", "TBC1D20", "CSNK2A1", "TCF15", "SRXN1", "SLC52A3", "FAM110A", 
#"RPS10L", "ANGPT4")
#
#unsuff = function(x) gsub("\\..*", "", x)
#
#data("gtex_b38_lung_chr20_exc", package="tnt4dn8")
#
#dn8like_ens_pheno_2_sym = function(dn8like_ens) {
#  rownames(dn8like_ens) = unsuff(dn8like_ens$pheno)
#}

#' filter a table with dn8like columns
#' @param dn8like data.frame
#' @param chr character(1)
#' @param start numeric(1)
#' @param end numeric(1)
#' @param posvbl character(1) column named for position
#' @examples
#' data("gtex_b38_lung_chr20_exc", package="tnt4dn8")
#' filter_dn8like( gtex_b38_lung_chr20_exc, 20, 60000, 62000, "pos_b38")
#' @export 
filter_dn8like = function(dn8like, chr, start, end, posvbl) {
 if (posvbl != "pos_b38") stop("code not using non-b38 positions at this time")
 dn8like$posvbl = dn8like[[posvbl]] # not good
 dn8like |> dplyr::filter(chr==chr & posvbl >= start & posvbl <= end)
}

#' produce a list of chrom, start, end, given gene symbol
#' @param sym character(1)
#' @param build character(1)
#' @examples
#' sym2region()
#' @export
sym2region = function(sym="SOX12", build="hg38") {
 if (build=="hg38") gdat = genes(EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79)
 else if (build=="hg19") gdat = genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
 else stop("only hg19 or hg38 available at this time")
 GenomeInfoDb::seqlevelsStyle(gdat) = "UCSC"
 rec = gdat[which(gdat$gene_name == sym)]
 if (length(rec)!=1) stop("symbol request does not yield single record.")
 list(chrom=seqnames(rec)[[1]], start=GenomicRanges::start(rec)[[1]],
     end=GenomicRanges::end(rec)[[1]])
}
 
#' run igvR to plot SNP p-values over gene models
#' @importFrom igvR igvR showGenomicRegion setGenome getGenomicRegion DataFrameQuantitativeTrack displayTrack
#' @importFrom BrowserViz setBrowserWindowTitle
#' @param chr character(1)
#' @param dn8like data.frame
#' @param start numeric(1)
#' @param end numeric(1)
#' @param posvbl character(1) column named for position
#' @param use character(1) defaults to "mlog10p"
#' @param posvbl character(1) what column in dn8like table holds genomic position
#' @param pvbl character(1) what column in dn8like table holds p-value
#' @examples
#' data("gtex_b38_lung_chr20_exc", package="tnt4dn8")
#' runBasic(dn8like=gtex_b38_lung_chr20_exc)
#' @export
runBasic = function(chr="20", start=60000, end=1.1e6, dn8like, use="mlog10p", posvbl="pos_b38", pvbl="p") {
 igv <- igvR()
 setBrowserWindowTitle(igv, "simple igvR demo")
 setGenome(igv, "hg38")
 showGenomicRegion(igv, list(chrom=chr, start=start, end=end))
 loc <- getGenomicRegion(igv)
 dn8like = filter_dn8like(dn8like, chr, start, end, posvbl=posvbl)
 starts = dn8like[[posvbl]]
 ends = starts+1
 values <- -log10(dn8like[[pvbl]]) 
 seqnames = dn8like[["chr"]]
 if (is.numeric(seqnames) | length(grep("^chr", seqnames[1]))==0) seqnames=paste0("chr", seqnames)
 tbl.bedGraph <- data.frame(chrom=seqnames, start=starts, end=ends,
                           value=-log10(dn8like[[pvbl]]), stringsAsFactors=FALSE)
 track <- DataFrameQuantitativeTrack("bedGraph", tbl.bedGraph, color="red", autoscale=TRUE)
 displayTrack(igv, track)
}
