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
#' @examples
#' if (requireNamespace("TnT") & requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
#' h = dn8_hub()
#' chk2 = h |> tbl("leic_dn8") |> filter_sym("HHIP", radius=1e6) |> as.data.frame()
#' dn8_hub_close(h)
#' print(tntplot(chk2))
#' }
#' @export
tntplot = function(tab, snpcolor="lightblue", genecolor="gold", gt=NULL) {
  if (!requireNamespace("TnT")) stop("install TnT to use this")
  if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this")
  if (is.null(gt)) gt = TnT::GeneTrackFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
      height=100, color=genecolor)  # consider making this optionally passed as a fixed object
  suppressMessages({
    syms = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=gt@Data$id, keytype="ENTREZID", column="SYMBOL")
  })
  gt@Data$display_label = TnT::strandlabel(syms, GenomicRanges::strand(gt@Data))
  t2g = tab2grngs(tab)
  tab$value = t2g$value
  pt = TnT::PinTrack( t2g, height=400, tooltip = as.data.frame(tab), color=snpcolor )
  TnT::TnTGenome(list(pt, gt), view.range=(range(t2g)+100000), coord.range=GenomicRanges::ranges(range(t2g)+5e6)[1])
}
    
