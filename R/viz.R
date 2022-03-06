#' simple plotter
#' @param tab table
#' @export
simpleviz = function(tab)
    ggplot2::ggplot(tab, aes(x=pos, y=-log10(p),
           text=paste("snp", "<br>", rsid, sep=""))) + 
           ggplot2::geom_point()

#' simple converter
#' @param tab table
#' @export
tab2grngs = function(tab) {
    if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this")
# avoid GenomeInfoDb
    t1 = tab$chr[1]
    pref = ""
    if (length(grep("chr", t1)) == 0) pref="chr"
    GenomicRanges::GRanges(paste0(pref, tab$chr), IRanges::IRanges(tab$pos, width=1), value=-log10(tab$p),
      id=tab$rsid)
}
  
#' TnT approach
#' @param tab table
#' @param snpcolor character(1)
#' @param genecolor character(1)
#' @param gt defaults to NULL, otherwise a GeneTrackFromTxDb-like object from TnT
#' @param maxp numeric(1) if non-NULL loci with p-values greater than this are excluded
#' @examples
#' if (requireNamespace("TnT") & requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")) {
#' data(gtex_b38_lung_chr20_exc)
#' cands = avail_syms_gtex()
#' cands
#' chk2 = gtex_b38_lung_chr20_exc |> filter_sym(cands[3], radius=5e5) |> as.data.frame()
#' print(tntplot(chk2))
#' }
#' @export
tntplot = function(tab, snpcolor="lightblue", genecolor="gold", gt=NULL, maxp = .1) {
  if (!requireNamespace("TnT")) stop("install TnT to use this")
  if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this")
  if (is.null(gt)) gt = TnT::GeneTrackFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
      height=100, color=genecolor)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=gt@Data$id, keytype="ENTREZID", column="SYMBOL")
  })
  gt@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(gt@Data))
  if (!is.null(maxp)) {
   todrop = which(tab$p > maxp)
   if (length(todrop)>0) tab = tab[-todrop,]
   }
  t2g = tab2grngs(tab)
  tab$value = t2g$value
  pt = TnT::PinTrack( t2g, height=400, tooltip = as.data.frame(tab), color=snpcolor )
  TnT::TnTGenome(list(pt, gt), view.range=(range(t2g)+100000), coord.range=GenomicRanges::ranges(range(t2g)+5e6)[1])
}
    
#' detailed plot
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' @param tab data.frame
#' @param snpcolor character(1)
#' @param genecolor character(1)
#' @param txcolor character(1)
#' @param GT defaults to NULL, otherwise a GeneTrackFromTxDb-like object from TnT
#' @param tt defaults to NULL, otherwise a TxTrackFromTxDb-like object from TnT
#' @param maxp numeric(1) if non-NULL loci with p-values greater than this are excluded
#' @param trunc_mlp numeric(1) defaults to 300, -log10p greater than this are reset to this value
#' @examples
#' if (requireNamespace("TnT") & requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")) {
#' data(gtex_b38_lung_chr20_exc)
#' cands = avail_syms_gtex()
#' cands
#' chk2 = gtex_b38_lung_chr20_exc |> filter_sym(cands[3], radius=5e4) |> as.data.frame()
#' print(tntplot2(chk2))
#' }
#' @export
tntplot2 = function(tab, snpcolor="lightblue", genecolor="darkgreen", txcolor="darkred", GT=NULL,
       tt = NULL, maxp = .1, trunc_mlp = 300) {
  if (!requireNamespace("TnT")) stop("install TnT to use this")
  if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this")
  uchr = as.character(unique(tab$chr))
  stopifnot(length(uchr)==1)
  seqlevelsStyle(uchr) = "UCSC"
  txdb38 = keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
  odb = org.Hs.eg.db::org.Hs.eg.db
# GENE
  if (is.null(GT)) GT = TnT::GeneTrackFromTxDb(txdb38,
      height=100, color=genecolor, seqlevel=uchr)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = AnnotationDbi::mapIds(odb, keys=GT@Data$id, keytype="ENTREZID", column="SYMBOL")
  })
  GT@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(GT@Data))
# TRANSCRIPT
  if (is.null(tt)) tt = TnT::TxTrackFromTxDb(txdb38,
      height=400, color=txcolor, seqlevel=uchr)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = tt@Data$tooltip$tx_name 
  })
  tt@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(tt@Data))
  if (!is.null(maxp)) {
   todrop = which(tab$p > maxp)
   if (length(todrop)>0) tab = tab[-todrop,]
   }
  t2g = tab2grngs(tab)
  val = t2g$value
  biginds = which(val > trunc_mlp)
  if (length(biginds)>0) {
    t2g$value[biginds] = trunc_mlp
    }
  tab$value = t2g$value
  pt = TnT::PinTrack( t2g, height=200, tooltip = as.data.frame(tab), color=snpcolor )
  TnT::TnTGenome(list(pt, GT, tt), view.range=(range(t2g)+10000), coord.range=GenomicRanges::ranges(range(t2g)+1e5)[1])
}
    
