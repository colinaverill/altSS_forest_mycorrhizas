#pGenerating shiny visualization.
rm(list=ls())
source('paths.r')
library(shiny)
library(tidyverse)

#load data.
data <- readRDS(shiny_viz_data.path)
#grab variable names.
var_names <- names(data)

#Setup variable drop down menu.
ui <- fluidPage(
  fluidRow(column(4, selectInput(inputId = "var", label = "Choose a variable:", choices = var_names)),
           column(8, plotOutput("plot"))
  )
)

#Setup server connection.
server <- function(input, output) {
  
  selectedData <- reactive({
    #grab the right data frame.
    data.frame(data[input$var]) 
  })
  #format variable name.
  varName <- reactive({str_replace_all(input$var, "_", " ")})
  
  #plot the relationship.
  output$plot <- renderPlot({
    df <- selectedData()
    names(df) <- c('x', 'y', 'lo95', 'hi95')
    #plot.
    limy <- c(min(df$lo95)*0.95, max(df$hi95)*1.05)
    plot(y ~ x, data = df, bty='l', cex = 0, ylim=limy,
         ylab = 'Relative Abundance EM Trees',
         xlab = varName())
    #trend line and uncertainty band.
    lines(smooth.spline(df$y ~ df$x), lwd = 2)
    polygon(c(df$x, rev(df$x)),c(df$hi95, rev(df$lo95)), col=adjustcolor('black', 0.3), lty=0)
  })
}

#Send the app to the server.
shinyApp(ui, server)
