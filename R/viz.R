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
#' if (requireNamespace("TnT") & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
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
  if (is.null(gt)) gt = TnT::GeneTrackFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
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
#' @param tab data.frame
#' @param snpcolor character(1)
#' @param genecolor character(1)
#' @param txcolor character(1)
#' @param gt defaults to NULL, otherwise a GeneTrackFromTxDb-like object from TnT
#' @param maxp numeric(1) if non-NULL loci with p-values greater than this are excluded
#' @examples
#' if (requireNamespace("TnT") & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
#' data(gtex_b38_lung_chr20_exc)
#' cands = avail_syms_gtex()
#' cands
#' chk2 = gtex_b38_lung_chr20_exc |> filter_sym(cands[3], radius=5e4) |> as.data.frame()
#' print(tntplot2(chk2))
#' }
#' @export
tntplot2 = function(tab, snpcolor="lightblue", genecolor="darkgreen", txcolor="darkred", gt=NULL, maxp = .1) {
  if (!requireNamespace("TnT")) stop("install TnT to use this")
  if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this")
  txdb19 = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  odb = org.Hs.eg.db::org.Hs.eg.db
# GENE
  GT = TnT::GeneTrackFromTxDb(txdb19,
      height=100, color=genecolor)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = AnnotationDbi::mapIds(odb, keys=GT@Data$id, keytype="ENTREZID", column="SYMBOL")
  })
  GT@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(GT@Data))
# TRANSCRIPT
  if (is.null(gt)) gt = TnT::TxTrackFromTxDb(txdb19,
      height=400, color=txcolor)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = gt@Data$tooltip$tx_name 
  })
  gt@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(gt@Data))
  if (!is.null(maxp)) {
   todrop = which(tab$p > maxp)
   if (length(todrop)>0) tab = tab[-todrop,]
   }
  t2g = tab2grngs(tab)
  tab$value = t2g$value
  pt = TnT::PinTrack( t2g, height=200, tooltip = as.data.frame(tab), color=snpcolor )
  TnT::TnTGenome(list(pt, GT, gt), view.range=(range(t2g)+10000), coord.range=GenomicRanges::ranges(range(t2g)+1e5)[1])
}
    
