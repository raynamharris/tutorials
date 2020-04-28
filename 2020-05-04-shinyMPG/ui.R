#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "xvar",
                  label = "x-axis variable ",
                  choices = c("displ", "cty"))
    ),
    
    mainPanel(
       plotOutput("scatterplot")
    )
  )
))




