

resources_in_sqlite_lap = function() list(
  copdg_nhw_eqtl = Sys.getenv("COPDG_NHW_EQTL_SQLITE_PATH"),
  ipf_gwas = Sys.getenv("IPF_SQLITE_PATH"),
  ltcopd_mqtl = Sys.getenv("LTCOPD_MQTL_SQLITE_PATH"),
  leic_copd_gwas = Sys.getenv("LEIC_SQLITE_PATH"))

rsq_connect = function(x) dbConnect(SQLite(), x, flags=SQLITE_RO)

#' provisional constructor of a hub to collect all resources and manage connections
#' @import RSQLite
#' @import methods
#' @importFrom dplyr tbl
#' @importFrom tibble as_tibble
#' @note The structure is a flat list with components connections, tables (just
#' a list of tablenames per connection), and status (defaults to "open")
#' only for use at channing
dn8_hub = function() {
  res = resources_in_sqlite_lap()
  conlist =  lapply(res, rsq_connect)
  names(conlist) = names(res)
  tablenames = lapply(conlist, dbListTables)
  if (any(duplicated(unlist(tablenames)))) stop("can't have duplicate tablenames across databases")
  names(tablenames) = names(res)
  ans = list(connections=conlist, tables=tablenames, status="open")
  class(ans) = "dn8_hub"
  ans
}


#' view hub
#' @param x instance of dn8_hub
#' @param \dots not used
#' @export
print.dn8_hub = function(x, ...) {
 cat("PDD data hub.\n")
 ncons = length(x$connections)
 ntabs = length(unlist(x$tables))
 cat(" ", ncons, " databases, ", ntabs, " tables.\n", sep="")
 cat(" tablenames: ", paste(select_some(tables(x)), collapse=", "), "\n")
 cat(" use [hub] |> tbl([resource]) |> filter_sym(sym, radius) to grab a table for a gene.\n")
}


#' close hub
#' @import DBI
#' @param x instance of dn8_hub
#' @note will run DBI::dbDisconnect
#' @export
dn8_hub_close = function(x) {
 stopifnot(inherits(x, "dn8_hub"))
 w = lapply(x$connections, function(z) DBI::dbDisconnect(z))
 x$status = "closed"
 invisible(w)
}

#' list the 'tables' available in a hub (really the database connection names)
#' @param x instance of dn8_hub
#' @export
tables = function(x) UseMethod("tables")
#' @export
tables.dn8_hub = function(x) unlist(x$tables)

#' select a tbl from a hub
#' @param src a dn8_hub instance
#' @param \dots should be character(1)
#' @export
tbl.dn8_hub = function(src, ...) {
 tbname = list(...)[[1]]
 avail = unlist(src$tables)
 stopifnot(tbname %in% avail)
 dbind = which(vapply(seq_len(length(src$tables)), function(x) tbname %in% src$tables[[x]],
      logical(1)))
 src$connections[[dbind]] |> dplyr::tbl(tbname)
}


#' filter a symbol from a table meeting dn8 naming conditions
#' @import rtracklayer
#' @param .data a tbl
#' @param sym a character(1) gene symbol
#' @param radius numeric(1) include `radius` bp upstream and downstream of gene body limits
#' @note if `pos_b38` is found in the dn8 table, gene addresses from EnsDb.Hsapiens.v79 are used to get gene
#' addresses in hg38 coordinates.  Otherwise, EnsDb.Hsapiens.v75 is used for gene coordinates.
#' @examples
#' if (interactive()) {
#'  data(gtex_b38_lung_chr20_exc)
#'  chk =  gtex_b38_lung_chr20_exc |> filter_sym("ANGPT4", radius=5e5) |> as.data.frame()
#'  print(names(chk))
#'  print(dim(chk))
#'  print(summary(-log10(chk$p)))
#'  chk2 = gtex_b38_lung_chr20_exc |> filter_sym("FAM110A", radius=1e6) |> as.data.frame()
#'  print(dim(chk2))
#'  print(summary(-log10(chk2$p)))
#' }
#' @export
filter_sym = function (.data, sym, radius=0 )
{
    if (!requireNamespace("GenomicRanges")) stop("install GenomicRanges to use this function")
    addr = get_addr_v75(sym)
    if ("pos_b38" %in% colnames(.data)) {
        addr = get_addr_v79(sym)
        }
    if (nrow(addr)==0) stop(paste("symbol", sym, "not found"))
    if (nrow(addr)>1) message("multiple addresses for", sym, "found, using first.")
    mychr = addr[1, "seqnames"]
    mys = addr[1, "start"]
    mye = addr[1, "end"]
    .data = harmonize_cols_for_filter(.data)  # ensure pos and chr are present for eqtl/mqtl sources
    if (!inherits(.data, "data.frame")) {
        dbname = .data$src$con@dbname
        if (length(grep("mQTLcis.sqlite", dbname))>0) mychr = paste0("chr", mychr)
        if (length(grep("subsetted_eqtl_cis_nhw.sqlite", dbname))>0) mychr = paste0("chr", mychr)
        }
    if ("pos_b38" %in% colnames(.data)) {
        if (!requireNamespace("rtracklayer")) stop("install rtracklayer to use this application")
        tmp = .data |> filter(chr == mychr & pos_b38 >= (mys-radius) & pos_b38 <= (mye+radius)) |> as.data.frame()
        tmpg = GenomicRanges::GRanges(tmp$chr, IRanges::IRanges(tmp$pos_b38, width=1))
        GenomicRanges::mcols(tmpg) = tmp
        testchr = as.character(GenomeInfoDb::seqnames(tmpg)[1])
        lchr = length(grep("^chr", testchr))
        if (lchr == 0) GenomeInfoDb::seqlevels(tmpg) = paste0("chr", GenomeInfoDb::seqlevels(tmpg))
        data(hg38ToHg19)
        remapped = rtracklayer::liftOver(tmpg, hg38ToHg19)
        GenomeInfoDb::seqlevelsStyle(remapped) = "UCSC"
        ans = unlist(remapped)
        ans = tibble::as_tibble(ans)  # for compatibility, don't return GRanges
        ans$pos = ans$start
# at this point, liftOver may have introduced another chromosome in the table!
        ans$seqnames = as.character(ans$seqnames)
        chrtab = table(ans$seqnames)
        nchr = length(chrtab)
        if (nchr > 1) {
           if (nchr > 2) stop("liftOver leads to multiple chromosomes")
           minor = chrtab[which.min(chrtab)]
           bad = which(ans$seqnames == names(minor))
           ans = ans[-bad,]
           }
        }       
    else ans = .data |> filter(chr == mychr & pos >= (mys-radius) & pos <= (mye+radius))
    ans
}


harmonize_cols_for_filter = function(tab) {
  dbname = try(tab$src$con@dbname, silent=TRUE)
  if (inherits(dbname, "try-error")) return(tab)
  if (length(grep("mQTLcis.sqlite", dbname))>0) 
     return(mutate(tab, chr=SNP.chr, pos=SNP.pos, p=pvalue, rsid=snps))
  if (length(grep("mQTLcis.all.sqlite", dbname))>0)  # FIXME 'snps' is 'marker format'
     return(mutate(tab, chr=site.chr, pos=SNP.pos, p=pvalue, rsid=snps)) # CIS so site.chr OK
  if (length(grep("subsetted_eqtl_cis_nhw.sqlite", dbname))>0) 
     return(mutate(tab, chr=chrom, pos=start, p=pvalue, rsid=rsnumber))
  return(tab)
  }


get_addr_v75 = function(sym) {
   data(genes_ensv75)
   genes_ensv75 |> 
   dplyr::select(seqnames, start, end, gene_name) |> dplyr::filter(gene_name == sym)
}

get_addr_v79 = function(sym) {
   data(genes_ensv79)
   genes_ensv79 |> 
   dplyr::select(seqnames, start, end, gene_name) |> dplyr::filter(gene_name == sym)
}

select_some = function (obj, maxToShow = 5) 
{
# from Biobase::selectSome
    len <- length(obj)
    if (maxToShow < 3) 
        maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        bot <- ceiling(maxToShow/2)
        top <- len - (maxToShow - bot - 1)
        nms <- obj[c(1:bot, top:len)]
        c(as.character(nms[1:bot]), "...", as.character(nms[-c(1:bot)]))
    }
    else if (is.factor(obj)) 
        as.character(obj)
    else obj
}

tbl.dn8_hub_old = function(src, ...) { x = list(...)[[1]]; toget = paste0(x, "_tbl"); src[[toget]] }
