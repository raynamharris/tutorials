#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

shinyServer(function(input, output) {
   
  output$scatterplot <- renderPlot({
  
   mpg %>%
      ggplot(aes_string(x = input$xvar)) + 
            geom_point(aes(y = hwy, 
                       color = cyl))
  })
  
})







