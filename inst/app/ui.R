

library(shiny)
fluidPage(
 shinytoastr::useToastr(),
 sidebarLayout(
  sidebarPanel(
   helpText("tnt4dn8 demo"),
   helpText("Use clickwheel to zoom, pointer to pan/drag, click on glyphs for details"),
   selectInput("gene", "gene", choices = sort(tnt4dn8::avail_syms_gtex()), selected="RBCK1"),
   numericInput("radius", "radius", min=0, max=1e6, step=1e4, value=5e4),
   actionButton("stopBtn", "stop app"),
   width=2
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("manh", verbatimTextOutput("gname"), TnT::TnTOutput("tntmanh")),
    tabPanel("about", helpText("simple demonstration, small set of genes on chr20, GTEx lung excerpt"),
      helpText("second track demo: EBI GWAS catalog; note all -log10 p values greater than 90 are truncated to 90; get actual p-value by clicking on glyph"),
      verbatimTextOutput("desc"))
   )
  )
 )
)
