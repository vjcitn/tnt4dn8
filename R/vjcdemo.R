#   GWASTrack(
#       trackName,
#       table,
#       chrom.col,
#       pos.col,
#       pval.col,
#       color = "darkBlue",
#       trackHeight = 50,
#       visibilityWindow = 1e+05
#     )

#stvjc@stvjc-XPS-13-9300:~/tnt4dn8$ head Pax*
#CHR	SNP	BP	A1	A2	FRQ	INFO	BETA	SE	P
#19	19:57859145:G:T	57859145	G	T	0.9259	0.8377	-0.0423	0.0306	0.1679

#' test single GWAS track
#' @param table data.frame
#' @param name character(1) used as trackName
#' @param chrom.col numeric(1)
#' @param pos.col numeric(1)
#' @param pval.col numeric(1)
#' @param \dots passed to GWASTrack
#' @examples
#' data(pax_demo)
#' vjcdemo(pax_demo)
#' @export
vjcdemo = function(table, name="gwasdemo",
   chrom.col=1, pos.col=3, pval.col=10, ... ) {
library(igvR)
library(igvShiny)
printf <- function(...) print(noquote(sprintf(...)))

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
ui = shinyUI(

   fluidPage(
    sidebarLayout(
     sidebarPanel(
      helpText("Enter symbol"),
      textInput("genesym", "sym", value="HHIP"),
      actionButton("gotSym", "GO"),
      width=2
      ),
     mainPanel(
       igvShinyOutput('igvShiny_0'),
       width=8
      )
     )
    )

) # ui
#----------------------------------------------------------------------------------------------------
server = function(input, output, session)
{
   gwtr = GWASTrack(table, trackName=name, chrom.col=chrom.col, pos.col=pos.col, pval.col=pval.col, ...)
   output$igvShiny_0 <- renderIgvShiny({
     cat("--- starting renderIgvShiny\n");
     # for a single track application
     genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38",  initialLocus="TERT")
     x <- igvShiny(genomeOptions,
                   displayMode="SQUISHED",
                   tracks=list(gwtr)
                   )
     displayTrack(x, gwtr)
     cat("--- ending renderIgvShiny\n");
     return(x)
     })
   observeEvent(input$gotSym, {
      printf("--- symbol")
      if(nchar(input$genesym) > 0)
        igvShiny::showGenomicRegion(session, id="igvShiny_0", input$genesym)
      })
} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui=ui, server=server))
}


