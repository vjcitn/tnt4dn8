library(shiny)
library(tnt4dn8)

#' demo server
#' @importFrom shinytoastr toastr_info
#' @importFrom TnT renderTnT
server = function(input, output, session) {
  build_gt = reactive({
    shinytoastr::toastr_info("building global gene models...")
    txdb19 = TxDb.Hsapiens.UCSC.hg19.knownGene
    tt = TnT::TxTrackFromTxDb(txdb19, height=400, color="darkred")
    GT = TnT::GeneTrackFromTxDb(txdb19, height=100, color="darkgreen")
    list(tt=tt, GT=GT)
    })
  output$gname = renderPrint({
    cat(paste(input$gene, "in GTEx Lung."))
  })
  output$tntmanh = TnT::renderTnT({
   shinytoastr::toastr_info("filtering annotations...", position="top-center")
   sym = input$gene
   rad = input$radius
   data(gtex_b38_lung_chr20_exc)
   tab = gtex_b38_lung_chr20_exc |> filter_sym(sym, radius=rad) |> as.data.frame()
   bgt = build_gt()
   tntplot2(tab, tt=bgt$tt, GT=bgt$GT)
  })
  output$desc = renderPrint({
   packageDescription("tnt4dn8")
  })
  observeEvent(input$stopBtn, {
    stopApp(NULL)
    })
}
