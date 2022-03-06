

library(shiny)
fluidPage(
 shinytoastr::useToastr(),
 sidebarLayout(
  sidebarPanel(
   helpText("tnt4dn8 demo"),
   selectInput("gene", "gene", choices = sort(tnt4dn8::avail_syms_gtex()), selected="RBCK1"),
   numericInput("radius", "radius", min=0, max=1e6, step=1e4, value=5e4),
   actionButton("stopBtn", "stop app"),
   width=2
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("manh", verbatimTextOutput("gname"), TnT::TnTOutput("tntmanh")),
    tabPanel("about", helpText("simple demonstration, small set of genes on chr20, GTEx lung excerpt"),
      verbatimTextOutput("desc"))
   )
  )
 )
)
