
suppressPackageStartupMessages({
library(shiny)
library(tnt4dn8)
data(gtex_b38_lung_chr20_exc)
})

#' data(gtex_b38_lung_chr20_exc)
#' cands = avail_syms_gtex()
#' cands
#' chk2 = gtex_b38_lung_chr20_exc |> filter_sym(cands[3], radius=5e5) |> as.data.frame()

server = function(input, output, session) {
  output$gname = renderPrint({
    cat(input$gene)
  })
  output$tntmanh = renderTnT({
   sym = input$gene
   rad = input$radius
   tab = gtex_b38_lung_chr20_exc |> filter_sym(sym, radius=rad) |> as.data.frame()
   tntplot(tab)
  })
  output$desc = renderPrint({
   packageDescription("tnt4dn8")
  })
}
  
