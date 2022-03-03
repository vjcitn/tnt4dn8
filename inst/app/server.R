
suppressPackageStartupMessages({
library(shiny)
library(tnt4dn8)  # need to speed up
library(shinytoastr)
data(gtex_b38_lung_chr20_exc)
})

#' data(gtex_b38_lung_chr20_exc)
#' cands = avail_syms_gtex()
#' cands
#' chk2 = gtex_b38_lung_chr20_exc |> filter_sym(cands[3], radius=5e4) |> as.data.frame()

server = function(input, output, session) {
  output$gname = renderPrint({
    cat(paste(input$gene, "in GTEx Lung."))
  })
  output$tntmanh = renderTnT({
   toastr_info("filtering annotations...", position="top-center")
   sym = input$gene
   rad = input$radius
   tab = gtex_b38_lung_chr20_exc |> filter_sym(sym, radius=rad) |> as.data.frame()
   tntplot2(tab)
  })
  output$desc = renderPrint({
   packageDescription("tnt4dn8")
  })
}
  
