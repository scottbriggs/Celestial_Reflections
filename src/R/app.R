library(shiny)
library(bslib)
library(here)
library(arrow)
library(stringr)

source(here("src", "R", "scripts.R"))

options(digits=12)

ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "minty"),
  
  titlePanel("Ephemeris"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("year", "Year", value = 2000, step = 1, min = -13200, max = 17191),
      numericInput("month", "Month", value = 1, step= 1, min = 1, max = 12),
      numericInput("day", "Day", value = 1, step = 1, min = 1, max = 31),
      textOutput("jdnText"),
      selectInput(
        "selectBody",
        "Solar System Body", choices = 
          c("Mercury" = "mer",
            "Venus" = "ven",
            "Earth-Moon Barycenter" = "emb",
            "Earth" = "earth",
            "Mars" = "mar",
            "Jupiter" = "jup",
            "Saturn" = "sat",
            "Uranus" = "ura",
            "Neptune" = "nep",
            "Pluto" = "plu",
            "Moon" = "moon",
            "Sun" = "sun",
            "Nutations" = "nut"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("DE441 Raw", tableOutput("de441Table")),
        tabPanel("Planetary Ephemeris")
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
    jd <- reactive({julianDayNumber(input$year, input$month,input$day)})
    
    output$jdnText <- renderText({str_c("Julian Day Number = ", jd())})
    
    body <- reactive({
      if (input$selectBody == "mer"){
        positionMercurySSB(jd())
      } else if (input$selectBody == "ven"){
        positionVenusSSB(jd())
      } else if (input$selectBody == "emb"){
        positionEMBSSB(jd())
      } else if (input$selectBody == "earth"){
        positionEarthSSB(jd())
      } else if (input$selectBody == "mar"){
        positionMarsSSB(jd())
      } else if (input$selectBody == "jup"){
        positionJupiterSSB(jd())
      } else if (input$selectBody == "sat"){
        positionSaturnSSB(jd())
      } else if (input$selectBody == "ura"){
        positionUranusSSB(jd())
      } else if (input$selectBody == "nep"){
        positionNeptuneSSB(jd())
      } else if (input$selectBody == "plu"){
        positionPlutoSSB(jd())
      } else if (input$selectBody == "moon"){
        positionMoonGEO(jd())
      } else if (input$selectBody == "sun"){
        positionSunSSB(jd())
      } else {
        nutationAngles(jd())
      }
    })
    
    output$de441Table <- renderTable(body(), digits = 12)
    
}

shinyApp(ui, server)