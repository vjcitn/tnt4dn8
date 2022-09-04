#' import signif_variant_gene_pairs from GTEx to data.frame, adding symbol and other metadata about genes
#' @importFrom data.table fread
#' @importFrom dplyr left_join mutate
#' @param pairs_path character(1), field "gene_id" will be join key
#' @param gencode_data.frame data.frame instance with additional metadata, field "Tx" will be join key
#' @note the variant_id field will be split by underscore and first 2 items are chr and variant position
#' @return data.frame with pairs records joined with gencode metadata and chr.variant and chr.pos
#' @examples
#' pairs_path = system.file("gtex_pairs/lung_demo.pairs.txt.gz", package="tnt4dn8")
#' data(gencode_26_df)
#' lit = gtex_pairs_to_data.frame(pairs_path, gencode_26_df)
#' head(lit)
#' @export
gtex_pairs_to_data.frame = function(pairs_path, gencode_data.frame) {
  inp = fread(pairs_path)
  vid = inp$variant_id
  svid = strsplit(vid, "_")
  chr = sapply(svid, "[", 1)
  pos = sapply(svid, "[", 2)
  ini = as.data.frame(left_join(inp, mutate(gencode_data.frame, gene_id=Tx), by="gene_id"))
  ini$chr.variant = chr
  ini$pos.variant = pos
  ini
}
