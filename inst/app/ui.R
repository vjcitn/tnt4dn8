suppressPackageStartupMessages({
library(shiny)
library(tnt4dn8)
})

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("tnt4dn8 demo"),
   selectInput("gene", "gene", choices = sort(avail_syms_gtex())),
   numericInput("radius", "radius", min=0, max=5e6, step=5e5, value=5e5),
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
