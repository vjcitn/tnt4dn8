

#' demo ui
#' @import shiny
#' @importFrom shinytoastr useToastr
#' @importFrom TnT TnTOutput
make_ui = function() fluidPage(
 useToastr(),
 sidebarLayout(
  sidebarPanel(
   helpText("tnt4dn8 demo"),
   selectInput("gene", "gene", choices = sort(avail_syms_gtex()), selected="RBCK1"),
   numericInput("radius", "radius", min=0, max=1e6, step=1e4, value=5e4),
   actionButton("stopBtn", "stop app"),
   width=2
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("manh", verbatimTextOutput("gname"), TnTOutput("tntmanh")),
    tabPanel("about", helpText("simple demonstration, small set of genes on chr20, GTEx lung excerpt"),
      verbatimTextOutput("desc"))
   )
  )
 )
)



#' demo server
#' @importFrom shinytoastr toastr_info
#' @importFrom TnT renderTnT
server = function(input, output, session) {
  build_gt = reactive({
    toastr_info("building global gene models...")
    txdb38 = TxDb.Hsapiens.UCSC.hg38.knownGene
    tt = TnT::TxTrackFromTxDb(txdb38, height=400, color="darkred")
    GT = TnT::GeneTrackFromTxDb(txdb38, height=100, color="darkgreen")
    list(tt=tt, GT=GT)
    })
  output$gname = renderPrint({
    cat(paste(input$gene, "in GTEx Lung."))
  })
  output$tntmanh = renderTnT({
   toastr_info("filtering annotations...", position="top-center")
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
  
#' app code
#' @export
tntapp = function() {
 ui = make_ui()
 shiny::runApp( list(ui=ui, server=server) )
}
