#########################################
#### PRO-simat Bioinformatics WepApp ####
#### Written by Rana Salihoglu       ####
#### Update:15_sep_2022              ####
#### Ⓒ                              ####
#########################################
source('build-jimena-app.R')

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

library(reshape2)                      
library(shiny)
library(htmlwidgets)
library(shinythemes)
library(shinydashboard)
library(shinyjs)
library(DT)
library(fresh)
library(graph)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(png)
library(knitr)
library(visNetwork)
library(plotly)
library(SQUADD)
library(igraph)
library(RMySQL)
library(shinyWidgets)
library(DOSE)
library(GO.db)
library(topGO)
library(tidyr)
library(gridSVG)
library(svglite)
library(svgPanZoom)
library(ggvis)
library(magicfor)
library(shinycssloaders)
library(shinyalert)
library(AnnotationDbi)
library(ReactomePA)
library(magrittr)
library(grid)
library(RColorBrewer)
library(rintrojs)
library(shinypanel)
library(spsComps)
library(stringr)
library(shinyBS)


#library(stringdist)
#library(tidyverse)




# Create the theme
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#2e2d2f"
  ),
  adminlte_sidebar(
    width = "300px",
    dark_bg = "#2e2d2f",
    dark_hover_bg = "#2c3e50"
  ),
  adminlte_global(
    content_bg = "#fdfefe"
  )
)


#----------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------


shinyUI(tagList(
  rintrojs::introjsUI(),
  #spsDepend("shinyCatch"),
  spsDepend("toastr"),
  
  
  dashboardPage(
    
    dashboardHeader(title = span("PRO ", 
                                 span("-Simat", 
                                      style = "color:#69b9da; font-size: 28px")),
                    tags$li(
                      class = "dropdown",
                      tags$style(".main-header {max-height: 75px}"),
                      tags$style(".main-header .logo {height:75px;}"),
                      tags$style(".navbar {min-height:75px !important}"),
                      tags$a(sidebarMenu(id = "tab2",
                                         menuItem(text = span("Home", style = "font-size: 20px"), icon = icon("glyphicon glyphicon-home",lib = "glyphicon"),  tabName = "Category_0"),
                                         menuItem(text = span("NETWORK", style = "font-size: 20px"),tabName = "Category_1"),
                                         
                                         menuItem(text = span("GO ENRICHMENT", style = "font-size: 20px"),tabName = "Category_2"),
                                         menuItem(text = span("KEGG PATHWAY", style = "font-size: 20px"),tabName = "Category_3"),
                                         menuItem(text = span("JIMENA", style = "font-size: 20px"),tabName = "Category_4"),
                                         menuItem(text = span("SQUADD", style = "font-size: 20px"),tabName = "Category_5")),
                             menuItem(uiOutput("tab")),
                             style = "padding-top: 0px;
                                 padding-right: 0px;
                                 padding-bottom: 0px;
                                 padding-left: 0px;"
                      ))
    ),
    sidebar = dashboardSidebar(
      hr(),
      width = 300,
      collapsed = TRUE,
      sidebarMenu(id="tabs",
                  hr(),
                  fluidRow(column(1),
                           column(10,
                                  
                                  menuItem(text = span("Home", style = "font-size: 16px"),icon = tags$i(br(),class = "fas fa-home", style="icon-align:center;font-size: 15px; color:#5a5861"), tabName = "Category_0"),
                                  hr(),
                                  menuItem(text = span("Protein Network", style = "font-size: 16px"),icon = tags$i(br(),class = "fas fa-network-wired", style="icon-align:center;font-size: 15px; color:#5a5861"),tabName = "Category_1"),
                                  hr(),
                                  menuItem(text = span("GO Enrichment", style = "font-size: 16px"), icon = tags$i(br(),class = "fas fa-database", style="icon-align:center;font-size: 15px; color:#5a5861"), tabName = "Category_2"),
                                  hr(),
                                  menuItem(text = span("KEGG Pathway", style = "font-size: 16px"), icon = tags$i(br(),class = "fas fa-bezier-curve", style="icon-align:center;font-size: 15px; color:#5a5861"),tabName = "Category_3"),
                                  hr(),
                                  menuItem(text = span("JIMENA Simulation", style = "font-size: 16px"),icon = tags$i(br(),class = "far fa-chart-bar", style="icon-align:center;font-size: 15px; color:#5a5861"),tabName = "Category_4"),
                                  hr(),
                                  menuItem(text = span("SQUADD", style = "font-size: 16px"), icon = tags$i(br(),class = "far fa-chart-bar", style="icon-align:center;font-size: 15px; color:#5a5861"),tabName = "Category_5"),
                                  hr()),
                           
                           column(1),
                           column(10,
                                  actionBttn(
                                    inputId = "reFresh",
                                    label = "Refresh",
                                    style = "jelly", 
                                    color = "primary",
                                    icon = tags$i(class = "fa fa-refresh")
                                  ),
                                  
                           )
                  )
      ),
      
      conditionalPanel("input.tabs == 'Category_1'",
                       fluidRow(column(1),
                                column(10,
                                       #hr(),
                                       introjsUI(),
                                       actionBttn(inputId="help",label="Help!", style = "jelly", 
                                                  color = "success", icon = icon("info-circle",verify_fa = FALSE)),
                                       
                                       hr(),
                                       tags$head(tags$style(HTML("#hostID + div>.selectize-input{
                                       border-color: red !important;
                                       font-style: italic;
                                       //color:#ff0000;}
                                       #pathogenID + div>.selectize-input{
                                       border-color: red !important;
                                       font-style: italic;
                                       //color:#ff0000;}
                                                                 "))),
                                       selectInput(inputId = "hostID",h4("Organism 1"), 
                                                   choices = c("",
                                                               "Arabidopsis thaliana" = '3702',
                                                               "Bacillus subtilis" = '1423,224308',
                                                               "Bacillus anthracis" = '1392',
                                                               "Bos taurus" = '9913',
                                                               "Caenorhabditis elegans" = '6239',
                                                               #"Cavia porcellus" = "10141",
                                                               "Danio rerio" = "7955",
                                                               "Drosophila melanogaster" = '7227',
                                                               "Escherichia coli" = '562,83333,83334,574521,155864', 
                                                               
                                                               "Gallus gallus" = '9031',
                                                               "Homo sapiens" = '9606',
                                                               "Helicobacter pylori" = "85962",
                                                               "Listeria monocytogenes" = "1639,169963,1906951,1334565,1230340,265669",
                                                               #"Listeria monocytogenes serotype 1/2a (strain EGD / Mackaness)"="1334565",
                                                               #"Listeria monocytogenes serotype 4b str. LL195"="1230340",
                                                               #"Listeria monocytogenes serotype 4b (strain F2365)"="265669",
                                                               "Mus musculus" = "10090",
                                                               #"Oryzias latipes" = "8090",
                                                               #"Ovis aries" = "9540",
                                                               "Plasmodium falciparum" ='5833,36329',
                                                               "Rattus norvegicus" = '10116',
                                                               "Saccharomyces cerevisiae" = '4932,559292',
                                                               #"Saccharomyces cerevisiae (strain ATCC 204508 / S288c)" = "559292",
                                                               "Schizosaccharomyces pombe (strain 972 / ATCC 24843)" = "284812",
                                                               "Staphylococcus aureus" = '1280,93061,158879',
                                                               "Treponema pallidum" = "160,243276",
                                                               "Vaccinia virus" = '10245',
                                                               "Vaccinia virus Copenhagen" = '10249',
                                                               "Vaccinia virus WR" = '10254',
                                                               "Vaccinia virus L-IPV" = '31531',
                                                               "Vaccinia virus Ankara" = '126794',
                                                               "All Vaccinia virus strain"='10245,10249,10254,31531,126794',
                                                               "Xenopus laevis" = "8355"
                                                               
                                                   ), selected = 1),
                                       
                                       
                                       selectInput("pathogenID",h4("Organism 2"), 
                                                   choices = c("",
                                                               "Arabidopsis thaliana" = '3702',
                                                               "Bacillus subtilis" = '1423,224308',
                                                               "Bacillus anthracis" = '1392',
                                                               "Bos taurus" = '9913',
                                                               "Caenorhabditis elegans" = '6239',
                                                               "Cavia porcellus" = "10141",
                                                               "Danio rerio" = "7955",
                                                               "Drosophila melanogaster" = '7227',
                                                               "Escherichia coli" = '562,83333,83334,574521,155864',
                                                               "Gallus gallus" = '9031',
                                                               "Homo sapiens" = '9606',
                                                               "Helicobacter pylori" = "85962",
                                                               "Listeria monocytogenes" = "1639,169963,1906951,1334565,1230340,265669",
                                                               #"Listeria monocytogenes serotype 1/2a (strain EGD / Mackaness)"="1334565",
                                                               #"Listeria monocytogenes serotype 4b str. LL195"="1230340",
                                                               #"Listeria monocytogenes serotype 4b (strain F2365)"="265669",
                                                               "Mus musculus" = "10090",
                                                               "Oryzias latipes" = "8090",
                                                               "Ovis aries" = "9540",
                                                               "Plasmodium falciparum" ='5833,36329',
                                                               "Rattus norvegicus" = '10116',
                                                               "Saccharomyces cerevisiae" = '4932,559292',
                                                               #"Saccharomyces cerevisiae (strain ATCC 204508 / S288c)" = "559292",
                                                               "Schizosaccharomyces pombe (strain 972 / ATCC 24843)" = "284812",
                                                               "Staphylococcus aureus" = '1280,93061,158879',
                                                               "Treponema pallidum" = "160,243276",
                                                               "Vaccinia virus" = '10245',
                                                               "Vaccinia virus Copenhagen" = '10249',
                                                               "Vaccinia virus WR" = '10254',
                                                               "Vaccinia virus L-IPV" = '31531',
                                                               "Vaccinia virus Ankara" = '126794',
                                                               "All Vaccinia virus strain"='10245,10249,10254,31531,126794',
                                                               "Xenopus laevis" = "8355",
                                                               "Vaccinia virus GLV-1h68" = "0",
                                                               "Vaccinia virus VVΔTKΔN1L" ="10"), selected = 1),
                                       
                                       selectInput("nodeShape", h4("Node Shape"), choices= c("",
                                                                                             "Circle" = "circle",
                                                                                             "Database" ="database",
                                                                                             "Diamond" = "diamond",
                                                                                             "Dot" = "dot",
                                                                                             "Star" = "star",
                                                                                             "Ellipse" = "ellipse", 
                                                                                             "Circle" = "circle",
                                                                                             "Box" = "box", 
                                                                                             "Text" = "text",
                                                                                             "Triangle" = "triangle",
                                                                                             "Triangle Down" = "triangleDown",
                                                                                             "Square" = "square", 
                                                                                             "Hexagon" = "hexagon",
                                                                                             "Icon" = "icon"),selected = 1),
                                       
                                       selectInput("layout", h4("Layout"),
                                                   choices=c("",
                                                             "Sugiyama" = "layout_with_sugiyama",
                                                             "Circle" = "layout_in_circle",
                                                             "Nicely" = "layout_nicely",
                                                             "Grid" = "layout_on_grid",
                                                             "Sphere" = "layout_on_sphere",
                                                             "Randomly" = "layout_randomly",
                                                             "Dh" = "layout_with_dh",
                                                             "Fr" = "layout_with_fr",
                                                             "Gem" = "layout_with_gem",
                                                             "Graphopt" = "layout_with_graphopt",
                                                             "KK" = "layout_with_kk",
                                                             "LGL" = "layout_with_lgl",
                                                             "MDS" = "layout_with_mds"),selected = 1),
                                       
                                       
                                       textAreaInput("multiProtein", h4("Multiple Protein")),
                                       actionBttn(
                                         inputId = "runProtein",
                                         label = "Submit",
                                         style = "gradient", 
                                         color = "primary",
                                         icon = icon("glyphicon glyphicon-play",lib = "glyphicon")),
                                       
                                )),
      ),
      
      
      conditionalPanel("input.tabs == 'Category_2'",
                       fluidRow(column(1),
                                column(10,
                                       #hr(),
                                       introjsUI(),
                                       actionBttn(inputId="help2",label="Help!", style = "jelly", 
                                                  color = "success", icon = icon("info-circle", verify_fa = FALSE)),
                                       hr(),
                                       
                                       prettyRadioButtons(
                                         inputId = "chosen_data",
                                         label = "Choose:", 
                                         choices = c("Continue with PPI data" = "PPI_data","Use example data" = "example_csv", "Upload your data"="user_csv"),
                                         icon = icon("check"), 
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly"
                                       ),
                                       
                                       
                                       conditionalPanel(condition = "input.chosen_data=='user_csv'",
                                                        introBox(
                                                          data.step = 1,
                                                          data.intro ="Upload your Differential Gene Expression data in .csv format.",fileInput("userFile", "Choose CSV File",
                                                                                                                                                accept = c(
                                                                                                                                                  "text/csv",
                                                                                                                                                  "text/comma-separated-values,text/plain",
                                                                                                                                                  ".csv"),
                                                          )),
                                                        checkboxInput("header", "Header", TRUE),
                                                        introBox(
                                                          data.step = 2,
                                                          data.intro ="If you upload wrong data you can remove it.",actionLink('reset', 'Delete csv file'),
                                                          tags$style("#reset {
                                                      background:#2d2a32  ;
                                                      color:#ab1a30;
                                                      font-size: 18px;
                                                      font-style: italic;}
                                                        #reset:hover {
                                                        background-color:#ab1a30;
                                                        color:#2d2a32;}"))
                                                        
                                       )
                                )
                       )
      ),
      
      conditionalPanel("input.tabs == 'Category_3'",
                       fluidRow(column(1),
                                column(10,
                                       #hr(),
                                       introjsUI(),
                                       actionBttn(inputId="help3",label="Help!", style = "jelly", 
                                                  color = "success", icon = icon("info-circle", verify_fa = FALSE)),
                                       hr(),
                                       sidebarSearchForm(textId = "searchText", buttonId = "searchButton", 
                                                         label = "KEGG Browser", icon = shiny::icon("search", verify_fa = FALSE))
                                ))),
      
      conditionalPanel("input.tabs == 'Category_4'",
                       fluidRow(column(1),
                                column(10,
                                       #hr(),
                                       introjsUI(),
                                       actionBttn(inputId="help4",label="Help!", style = "jelly", 
                                                  color = "success", icon = icon("info-circle", verify_fa = FALSE)),
                                       hr(),
                                       prettyRadioButtons(
                                         inputId = "choseJimena",
                                         label = "Choose:", 
                                         choices = c("Use example data" = "exampleJim", "Upload your data"="userJim"),
                                         icon = icon("check"), 
                                         bigger = TRUE,
                                         status = "info",
                                         animation = "jelly"
                                       ),
                                       conditionalPanel(condition = "input.choseJimena=='userJim'",
                                                        
                                                        fileInput(inputId = "file1", h4("Choose Txt File"), 
                                                                  accept = c("text/csv","text/comma-separated-values,text/plain", ".csv"))),))),
      
      conditionalPanel("input.tabs == 'Category_5'",
                       fluidRow(column(1),
                                column(10,
                                       #hr(),
                                       introjsUI(),
                                       actionBttn(inputId="help5",label="Help!", style = "jelly", 
                                                  color = "success", icon = icon("info-circle", verify_fa = FALSE)),
                                       hr())))
      
      
      
      
      
    ),
    body = dashboardBody(
      shinyjs::useShinyjs(),
      use_theme(mytheme),
      
      
      tabItems(
        ################################  H o m e // UI // #######################################################################
        tabItem(tabName = "Category_0",
                fluidPage(tags$head(tags$style(HTML(".small-box.bg-yellow {cursor: pointer; top: 0px;background-image: linear-gradient(to left top,#373f51,#1b1b1e,#373f51) !important;
          border-radius: 25px; background-size: auto !important; left: 9%;margin-left: 45px; height: 250px}"))), 
                          tags$head(tags$style(HTML(".small-box.bg-maroon {cursor: pointer;top: 0px;background-image: linear-gradient(to left top,#b8b8ff,#3d348b,#b8b8ff) !important;
          border-radius: 25px;background-size: auto !important; left: 22%;;margin-left: 45px;height: 250px}"))),
                          tags$head(tags$style(HTML(".small-box.bg-purple {cursor: pointer;background-image: linear-gradient(to left top,#b8b8ff,#3d348b,#b8b8ff) !important;
          border-radius: 25px;background-size: auto !important;left: 22%;;margin-right: 45px;height: 250px}"))),
                          tags$head(tags$style(HTML(".small-box.bg-olive {cursor: pointer;background-image: linear-gradient(to left top,#0f0c29,#302b63,#24243e) !important;
          border-radius: 25px;background-size: auto !important;left: 22%;;margin-left: 45px; height: 250px}}"))),
                          tags$head(tags$style(HTML(".small-box.bg-green {cursor: pointer;background-image: linear-gradient(to left top,#0f0c29,#302b63,#24243e) !important;
          border-radius: 25px;background-size: auto !important;left: 22%;;margin-right: 45px; height: 250px}}"))),
                          
                          fluidRow(
                            div(id='clickdiv',valueBox(br(),subtitle = tags$p("NETWORK ANALYSIS", 
                                                                              style = "text-align:center; font-size: 50px;"),
                                                       icon = tags$i(class = "fas fa-dna", style="icon-align:center;font-size: 64px; color:#000000"),  color = "yellow", width = 10)
                            ),
                            
                            div(id='clickdiv2',valueBox(br(), subtitle = tags$p("GO ENRICHMENT ANALYSIS", style = "text-align:center;font-size: 50px",),
                                                        icon = tags$i(class = "fas fa-database", style="font-size: 64px; color:#3d348b"),color = "maroon", width = 5),
                            ),
                            
                            div(id='clickdiv3',valueBox(br(), subtitle = tags$p("KEGG PATHWAY ANALYSIS", style = "text-align:center; font-size: 50px"),
                                                        icon = tags$i(class = "fas fa-project-diagram", style="font-size: 64px; color:#3d348b"), color = "purple", width = 5),
                            ),
                            
                            div(id='clickdiv4',valueBox(br(), subtitle = tags$p("SIMULATION WITH JIMENA", style = "text-align:center; font-size: 50px"),
                                                        icon = tags$i(class = "far fa-chart-bar", style="font-size: 64px; color:#302b63"), color = "olive", width = 5),
                            ),
                            div(id='clickdiv5',valueBox(br(), subtitle = tags$p("VISUAL WITH SQUADD", style = "text-align:center; font-size: 50px"),
                                                        icon = tags$i(class = "far fa-chart-bar", style="font-size: 64px; color: #302b63"), color = "green", width = 5),
                            ),
                            
                          ),
                          tags$footer(HTML("
                    <!-- Footer -->
                           <footer class='page-footer font-large indigo'>
                            <footer style='padding:8px;background:#322F34;color: white;'>
                           <!-- Copyright -->
                           <div class='footer-copyright text-center py-3'>© 2022 Copyright: Dandekar's Group
                           <p></p>
                           <a href='https://www.biozentrum.uni-wuerzburg.de/bioinfo/research/groups/funct-genomics-systems-biology/'> www.biozentrum.uni-wuerzburg.de</a>
                           </div>
                           <!-- Copyright -->

                           </footer>
                           <!-- Footer -->"))
                )
        ),
        ################################### P P I -  N e t w o r k // UI // ##################################################################
        
        tabItem(tabName = "Category_1",
                fluidPage(
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  headerPanel("Protein-Protein Interaction",
                  ),
                  
                  helpText(h4(em("Please select the organisms you want or add your protein list in as UniProtKB to the multiple protein box and choose the layout for the network to be created. If you want to add and delete Node and Edge, click the edit button. You can add nodes anywhere in the network and change the label. Click the Add Edge button to add Edge. Click on a node and drag the edge to another node to connect them. "))),
                  #uiOutput("visNodeL"),
                  selectizeInputWithButtons(
                    
                    inputId = 'SearchVisNode',
                    label = 'Search Gene',
                    label_title = 'Enter the UniprotKB id of the gene you want to search in the protein-protein network you created.',
                    
                    actionButton(
                      style = "color: white; 
                     background-color: #009BFF; 
                     position: relative; 
                     left: 3%;
                     height: 35px;
                     width: 50px;
                     text-align:center;
                     text-indent: -2px;
                     border-radius: 0px;
                     border-width: 0px",
                      'btn1',
                      '',
                      icon('glyphicon glyphicon-search', verify_fa = FALSE, lib = "glyphicon"),
                      title = 'Click to search'
                    ),
                    
                    options = list(create = TRUE, placeholder = 'Enter UniprotKB')
                  ),
                  visNetworkOutput("editable_network", width = "100%", height = "700px"),
                  
                  ##fluidRow(
                  column(1,
                         actionBttn(
                           inputId = "MoReBtn",
                           label = "More",
                           style = "gradient", 
                           color = "primary",
                           icon = icon("glyphicon glyphicon-plus", lib = "glyphicon")
                         )),
                  column(1,
                         actionBttn(
                           inputId = "LessBtn",
                           label = "Less",
                           style = "gradient", 
                           color = "primary",
                           icon = icon("glyphicon glyphicon-minus",lib = "glyphicon")
                         )),
                  column(
                    7,
                    actionBttn(inputId='ab1', style = "gradient",
                              color = "primary",
                              div(HTML(
                                       "<a href='https://string-db.org/' target='_blank'>Learn More</a>")))
                   
                     
                    #actionBttn(inputId='ab1', label="Learn More", style = "gradient",
                     #          color = "primary",
                               #icon = icon("th"), 
                      #         onclick ="window.open('https://string-db.org/', '_blank')")
                    ),
                  column(3,
                         downloadLink('downloadNetwork', 'Download network')),
                  
                  
                  #),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  #mainPanel(  
                  #)
                  
                  box( status = "primary", id = "interaction_table",
                       hr(),DT::dataTableOutput("table_int"),width = "600%",
                       ),data.step = 1
                )
        ),
        #####################################  G O - E n r i c h m e n t // UI // ##################################################
        tabItem(tabName = "Category_2",
                fluidPage(
                  
                  #spsDepend("shinyCatch"),
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  headerPanel("GO Enrichment Analysis"),
                  mainPanel(width = 12,
                            
                            tags$head(tags$style(HTML(".box.box-solid.box-primary {top: 0px;background-image: linear-gradient(to right top,  #2471a3 , #7fb3d5 ) !important;
          background-size: auto !important}"))),
                            br(),
                            helpText(h4(em("Load your differential gene expression data and select the columns containing the Gene names and log2FC values. Choose the keytype of your genes.
Set your parameters for GO Enrichment and KEGG pathway analysis. Then click the Run Analysis button. "))),
                            #),
                            br(),
                            #conditionalPanel(condition = "input.chosen_data=='example_csv'",
                            #conditionalPanel(condition = "input.chosen_data=='user_csv'",
                            
                            shinyjs::useShinyjs(),
                            
                            fluidRow(
                              
                              column(width = 5,
                                     shinyjs::hidden(
                                       div(id = "advanced",
                                           box(width = NULL,height = '400px',status = "primary", solidHeader = TRUE,collapsed = FALSE, collapsible = TRUE, id="prePar1",
                                               
                                               title = HTML("Select Data Columns", "<font size='5'>",
                                                            as.character(actionLink(inputId = "info1", 
                                                                                    label = "", 
                                                                                    icon = icon("info-circle", verify_fa = FALSE),)), "</font>"),
                                               wellPanel(
                                                 
                                                 column(
                                                   4,
                                                   introBox(
                                                     data.step = 4,
                                                     data.intro ="Select the column with Gene in the data you have loaded.",selectInput('geneColumn',h4('Gene column'), ""))),
                                                 
                                                 
                                                 
                                                 column(
                                                   4,
                                                   introBox(
                                                     data.step = 5,
                                                     data.intro ="Select the column with Log2FC values in the data you have loaded.",selectInput('log2FColumn',h4('log2FC column'), ""))),
                                                 
                                                 column(
                                                   4,
                                                   introBox(
                                                     data.step = 6,
                                                     data.intro ="Select the column with q-value in the data you have loaded.",selectInput('pAdColumn',h4('pAdjust column'), ""))),
                                                 
                                                 
                                                 introBox(
                                                   data.step = 7,
                                                   data.intro ="Set treshold for log2FC",numericInput("logFCut", h4("min log2FC"), 1,
                                                                                                      min = -25, max = 25)),
                                                 
                                                 
                                                 introBox(
                                                   data.step = 8,
                                                   data.intro ="Set treshold for q-value",numericInput("pAdjCut", h4("cut pAdjust"), 0.05,
                                                                                                       min = 0, max = 1, step =0.01))
                                                 
                                                 
                                                 
                                               )))),
                                     
                                     box(width = NULL,height = '320px',status = "primary", solidHeader = TRUE, collapsible = TRUE,
                                         title = HTML("Parameters for GO and KEGG analysis", "<font size='5'>",
                                                      as.character(actionLink(inputId = "info2", 
                                                                              label = "", 
                                                                              icon = icon("info-circle", verify_fa = FALSE),)), "</font>"),
                                         
                                         column(
                                           4,
                                           introBox(
                                             data.step = 9,
                                             data.intro ="The q-value refers to the ratio of the number of differential expression genes to the number of total annotated genes in a certain pathway. The q-value is a multiple hypothesis-corrected p-value and consists of values between 0 and 1; values closer to 0 indicate a more significant enrichment.
                                      Set q-value treshold for GO and KEGG analysis (between 0-1)",numericInput("qvalCut", h4("qValue"), 0.05,
                                                                                                                min = 0, max =1, step =0.01))),
                                         column(
                                           4,
                                           introBox(
                                             data.step = 10,
                                             data.intro ="P-value represents the probability or chance of observing at least x number of genes out of the total n genes in the list annotated to a particular GO term, given the proportion of genes in the whole genome that are annotated to that GO Term. The closer the p-value is to zero, the more significant the particular GO term associated with the group of genes is.
                                      Set p-value treshold for GO and KEGG analysis (between 0-1)",numericInput("pvalCut", h4("pValue"), 0.05,
                                                                                                                min = 0, max = 1, step =0.01))),
                                         
                                         column(
                                           4,
                                           
                                           introBox(
                                             data.step = 11,
                                             data.intro ="Select the organism to which your data belongs",selectInput("GO_species",h4("Species for GO"), 
                                                                                                                      choices = c(""), selected = 1))),
                                         
                                         column(
                                           4,
                                           
                                           introBox(
                                             data.step = 12,
                                             data.intro ="Select which keytpe your genes are in.",selectInput('keytype',h4('Select keytype'), 
                                                                                                              choices = c(""), selected = NULL))),
                                         
                                         column(
                                           4,
                                           introBox(
                                             data.step = 13,
                                             data.intro ="Select the organism to which your data belongs for KEGG analysis. 
                                      Go to the KEGG tab to see the KEGG analysis results.",selectInput("KEGG_species",h4("Species for KEGG"), 
                                                                                                        choices = c("",
                                                                                                                    "7165-Anopheles gambiae" = "aga",
                                                                                                                    "3702-Arabidopsis thaliana" = "ath",
                                                                                                                    #"224308-Bacillus subtilis" ="bsu",
                                                                                                                    #"198094-Bacillus anthracis" ="ban",
                                                                                                                    "9913-Bos taurus" = "bta",
                                                                                                                    #"5476-Candida albicans" = "cal",
                                                                                                                    "9612-Canis familiaris" = "cfa",
                                                                                                                    "6239-Caenorhabditis elegans" = "cel",
                                                                                                                    "7955-Danio rerio" = "dre",
                                                                                                                    "7227-Drosophila melanogaster" = "dme",
                                                                                                                    "562-Escherishia coli" = "eco",
                                                                                                                    "386585-Escherichia coli O157:H7 Sakai" = "ecs",
                                                                                                                    "155864-Escherichia coli O157:H7 EDL933" ="ece",
                                                                                                                    "9031-Gallus gallus" = "gga",
                                                                                                                    "9606-Homo sapiens" = "hsa",
                                                                                                                    "1639-Listeria monocytogenes" = "lmo",
                                                                                                                    "10090-Mus musculus" = "mmu",
                                                                                                                    "9544-Macaca mulatta" = "mcc",
                                                                                                                    "9598-Pan troglodytes" = "ptr",
                                                                                                                    "Myxococcus xanthus" ="mxa",
                                                                                                                    #"5833-Plasmodium falciparum" = "pfa",
                                                                                                                    "10116-Rattus norvegicus" = "rno",
                                                                                                                    #"243276-Treponema pallidum" = "tpa",
                                                                                                                    #"1280-Staphylococcus aureus" = "sau",
                                                                                                                    #"4896-Schizosaccharomyces pombe" = "spo",
                                                                                                                    "4932-Saccharomyces cerevisiae" = "sce",
                                                                                                                    "9823-Sus scrofa" = "ssc",
                                                                                                                    "8355-Xenopus laevis" ="xla"
                                                                                                        )))),
                                         
                                         column(
                                           4,
                                           
                                           introBox(
                                             data.step = 14,
                                             data.intro ="Select the method you want GO enrichment analysis to be done.",selectInput('pAdMeth',h4('pAdjust method'), 
                                                                                                                                     choices = c("","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")))),
                                         
                                         introBox(
                                           data.step = 15,
                                           data.intro ="You can now start the GO and KEGG analyzes by clicking the run Analysis button.",actionBttn(
                                             inputId = "submitButton",
                                             label = "Run Analysis",
                                             style = "gradient", 
                                             color = "success",
                                             icon = icon("glyphicon glyphicon-play",lib = "glyphicon")
                                           )))
                                     
                                     
                              ),
                              
                              box(width = 7, #height = '740px',
                                  title = "Gene Expression Data",status = "primary", collapsed = FALSE,collapsible = TRUE, id="upCsv",
                                  introBox(
                                    data.step = 3,
                                    data.intro ="You can see your uploaded data."),
                                  
                                  conditionalPanel(condition = "input.chosen_data == input$chosen_data",
                                                   DT::dataTableOutput(outputId = "whichosen"))
                              )
                              
                            ),     
                            
                            #)),
                            
                            
                            
                            
                            fluidRow(
                              introBox(
                                data.step = 16,
                                data.intro ="Genes and special sub gene ontology terms are divided into 5 clusters according to their functionality.",
                                box(status = "primary", collapsible = TRUE, width = 12,
                                    
                                    title = HTML("Cnetplot", "<font size='5'>",
                                                 as.character(actionLink(inputId = "info12", 
                                                                         label = "", 
                                                                         icon = icon("info-circle",verify_fa = FALSE),)), "</font>"), 
                                    
                                    selectInput('selectCluster',h4('Select Cluster'), 
                                                choices = c("","Apoptotic Process","Negative regulation","Positive regulation","Proliferation","Kinase","Other")),
                                    visNetworkOutput(outputId = "fig3" ,width = "100%", height = "500px")))),
                            
                            fluidRow(
                              
                              box(
                                status = "primary", collapsible = TRUE, width = 12,
                                tabsetPanel(
                                  tabPanel(
                                    title = HTML("Volcano Plot", "<font size='5'>",
                                                 as.character(actionLink(inputId = "info3", 
                                                                         label = "", 
                                                                         icon = icon("info-circle", verify_fa = FALSE),)), "</font>"), 
                                    introBox(
                                      data.step = 17,
                                      data.intro ="An interactive volcanoplot is created in line with the tresholds you have determined.",plotlyOutput("nemrut",width = "100%", height = "500px"))),
                                  
                                  tabPanel(
                                    title = HTML("Gene Ontology Graphs", "<font size='5'>",
                                                 as.character(actionLink(inputId = "info4", 
                                                                         label = "", 
                                                                         icon = icon("info-circle", verify_fa = FALSE),)), "</font>"), 
                                    column(6,
                                           introBox(
                                             data.step = 18,
                                             data.intro ="Choose the type of chart you want to see",selectInput("graphType",h4("Graph Type"), 
                                                                                                                choices = c("",
                                                                                                                            "EGO ALL Dotplot" = "fig4",
                                                                                                                            "EGO ALL Treeplot" = "fig5",
                                                                                                                            "EGO ALL Heatmap" = "heatplt",
                                                                                                                            "EGO ALL Cnetplot" ="plot4",
                                                                                                                            "Publications trend" = "pmc1",
                                                                                                                            "EGO CC Barplot" = "plot1",
                                                                                                                            "EGO CC Dotplot" = "plot2",
                                                                                                                            "EGO CC Cnetplot" ="plot3",
                                                                                                                            "EGO CC GO Graph" = "plotGO1",
                                                                                                                            "EGO MF Barplot" ="plot5",
                                                                                                                            "EGO MF Dotplot" ="plot6",
                                                                                                                            "EGO MF Cnetplot" ="plot7",
                                                                                                                            "EGO MF GO Graph" = "plotGO2",
                                                                                                                            "EGO BP Barplot" ="plot9",
                                                                                                                            "EGO BP Dotplot" ="plot10",
                                                                                                                            "EGO BP Cnetplot" ="plot11",
                                                                                                                            "EGO BP GO Graph" = "plotGO3"
                                                                                                                ), selected = 1))),
                                    
                                    column(6,
                                           introBox(
                                             data.step = 19,
                                             data.intro ="You can see the updated chart by changing the category size.",numericInput("categorySize", h4("Show Category"), 0,
                                                                                                                                     min = 0, max = 25))),
                                    
                                    
                                    conditionalPanel(condition = "input.graphType == input$graphType",
                                                     
                                                     svgPanZoomOutput(outputId = "whichplot2", width = "100%", height = "500px"),
                                                     br(),
                                                     
                                                     introBox(
                                                       data.step = 20,
                                                       data.intro ="You can download the graphic of your choice in png format",downloadButton(outputId = "downloadPlot", label = "Download"),
                                                       helpText(h4(em("Rectangle nodes in the GO Graphs are significant nodes")))),  
                                                     br(),
                                                     
                                    ),
                                    
                                    br()
                                  ),
                                ))
                            ), 
                            
                            
                            
                            box(title = "Table",status = "primary", collapsible = TRUE,width = "100%",
                                introBox(
                                  data.step = 21,
                                  data.intro ="You can see the tables created as a result of the GO enrichment analysis by selecting them here.",selectInput("selecTable",h4("Select Table"), 
                                                                                                                                                             choices = c("",
                                                                                                                                                                         "Biological Process" = "table1",
                                                                                                                                                                         "Molecular Function" = "table2",
                                                                                                                                                                         "Cellular Component" = "table3",
                                                                                                                                                                         "All Categories" = "table4"
                                                                                                                                                             ), selected = 1)),
                                
                                conditionalPanel(condition = "input.selecTable == input$selecTable",
                                                 hr(),DT::dataTableOutput(outputId = "whichtable")), 
                                
                                
                                
                                
                            )
                            
                            
                  )
                )
        ),
        ################################   K e g g  // UI // ###################################################################
        tabItem(tabName = "Category_3",
                fluidPage(
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  headerPanel("KEGG Pathway Analysis"),
                  mainPanel(width = 12, 
                            br(),
                            helpText(h4(em("To see the pathway colored according to the logFC values you have given, paste an ID of your interest from the data obtained as a result of the KEGG Pathway analysis into the input pathway box. For example: hsa05212 "))),
                            br(),
                            tabsetPanel(id = "keggtabs",
                                        tabPanel(
                                          "Output Table",
                                          br(),
                                          box( 
                                            title = HTML("KEGG Enrich Output Data", "<font size='5'>",
                                                         as.character(actionLink(inputId = "info5", 
                                                                                 label = "", 
                                                                                 icon = icon("info-circle", verify_fa = FALSE),)), "</font>"),status = "danger", id = "kegg1",
                                            hr(),DT::dataTableOutput("kk_df"),width = "600%"),data.step = 1,
                                          data.intro = "Took a look at this histogram"
                                          
                                          
                                        ),
                                        
                                        tabPanel("Plot",
                                                 
                                                 br(),
                                                 box(
                                                   title = HTML("Emaplot", "<font size='5'>",
                                                                as.character(actionLink(inputId = "info6", 
                                                                                        label = "", 
                                                                                        icon = icon("info-circle", verify_fa = FALSE),)), "</font>"), 
                                                   status = "primary", collapsible = TRUE, width = 12, id = "kegg2",
                                                   numericInput("categorySizeEma", h4("Show Category"), 0,
                                                                min = 0, max = 25),
                                                   svgPanZoomOutput(outputId = "plot13"),
                                                   downloadButton(outputId = "downloadEma", label = "Download"),
                                                   br()
                                                 )
                                        ),
                                        tabPanel("Pathway",
                                                 
                                                 br(),
                                                 box(
                                                   title = HTML("Pathview", "<font size='5'>",
                                                                as.character(actionLink(inputId = "info7", 
                                                                                        label = "", 
                                                                                        icon = icon("info-circle", verify_fa = FALSE),)), "</font>"),
                                                   status = "primary", collapsible = TRUE, width = 12, id = "pathview1",
                                                   searchInput(
                                                     inputId = "PathviEw",
                                                     placeholder = "Input pathway",
                                                     btnSearch = icon("search", verify_fa = FALSE), 
                                                     btnReset = icon("trash"),
                                                     width = "100%",
                                                     
                                                   ),
                                                   #textInput("PathviEw", h4("Input pathway")),
                                                   tags$div(class = "clearBoth"),
                                                   
                                                   imageOutput(outputId = "pathway1", inline = T),
                                                   br(),
                                                   conditionalPanel("output.pathviewPlotsAvailable",
                                                                    
                                                                    downloadBttn('downloadPathviewPng','Download .png',style = "bordered",
                                                                                 color = "primary"
                                                                    ),
                                                                    
                                                                    downloadBttn('downloadPathviewPdf','Download .pdf',style = "bordered",
                                                                                 color = "primary"
                                                                    ),
                                                   )
                                                 ),
                                                 
                                        ))
                  )
                )
        ),
        #################################  J i m e n a // UI // ############################################################### 
        tabItem(tabName = "Category_4",
                fluidPage(
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  headerPanel("JIMENA"),
                  mainPanel(
                    width = 12,
                    fluidRow(
                      
                      br(),
                      column(width = 6,
                             box(status = "primary", collapsible = TRUE, width = NULL, id = "bxjim",
                                 title = HTML("",
                                              as.character(actionLink(inputId = "info8", 
                                                                      label = "Step 1",
                                                                      icon = icon("info-circle", verify_fa = FALSE, "</font size='5'>"),))),
                                 helpText(em("Upload your signaling network data.
                                    Your (.txt) file should contain activation/inhibition information, node1, label and node2 as column headers.
                                    Then please click the (convert graphml) button.
                                    Now, you can start the jimena simulation. ")),hr(),
                                 #conditionalPanel(condition = "input.choseJimena == input$choseJimena",
                                 DT::dataTableOutput(outputId = "whichjim")
                             ),
                             
                             
                             box(status = "primary", collapsible = TRUE, width = NULL, id = "bxjim2",
                                 actionBttn(
                                   inputId = "runGraphml",
                                   label = "Convert Graphml",
                                   style = "fill", 
                                   color = "warning"
                                 ),
                                 br(),
                                 hr(),
                                 actionBttn(
                                   inputId = "runJimena",
                                   label = "Run Jimena",
                                   style = "gradient", 
                                   color = "primary",
                                   icon = icon("glyphicon glyphicon-play",lib = "glyphicon")
                                 ))),
                      
                      box(width = 6,
                          status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = TRUE,
                          title = HTML("Perturbation setting", 
                                       as.character(actionLink(inputId = "info9", 
                                                               label = "Step 2", 
                                                               icon = icon("info-circle", verify_fa = FALSE,"<font size='5'>"),))), 
                          
                          uiOutput('inVar2'),
                          
                          numericInput("strt", "Start", 0,
                                       min = 0, max = 1000),
                          numericInput("end", "End", 1000,
                                       min = 1, max = 1000),
                          numericInput("val", "Value", 1,
                                       min = 0, max = 1, step =0.01),
                          
                          
                          actionBttn(inputId = "add_pert",
                                     label = "Add New Perturbation",
                                     color = "success",
                                     style = "simple"
                          ),
                          
                          actionBttn(inputId = "remove_per",
                                     label = "Remove Perturbation",
                                     color = "danger",
                                     style = "simple"),
                          
                      ),
                      
                      helpText(em("Select the node you want to add and determine the time interval when the node is activated or inhibited in the simulation and click the add perturbation button.
(start : start time, end: end time, value:1-0 (activation/inhibition)).
You can also give values like 0.5 so that it can be partially active. Then run the simulation again by clicking the Run Jimena button. ")),hr(),
                    ),
                    
                    
                    DT::dataTableOutput("shiny_table"),
                    fluidRow(
                      
                      br(),
                      box(status = "primary", collapsible = TRUE, width = 12,
                          
                          plotlyOutput("jimenaResult")%>% withSpinner(type = 6, color = "#6a57ad", size = 1.5),
                          br()),
                    ),
                    
                    fluidRow(
                      box( status = "primary",collapsible = TRUE, id = "jimena_table",width = 12,
                           hr(),DT::dataTableOutput("jimenaoutput_tb")),data.step = 1)
                    
                  )
                )),
        ##############################  S Q  U A D D // UI // ###############################################################################             
        
        tabItem(tabName = "Category_5",
                fluidPage(
                  tags$style(type="text/css",
                             ".shiny-output-error { visibility: hidden; }",
                             ".shiny-output-error:before { visibility: hidden; }"
                  ),
                  headerPanel("SQUADD"),
                  mainPanel("Simulation with SQUADD Package", id="mainsq",
                            br(),
                            helpText(em("Please run Jimena before using the SQUADD part and identify the perturbation nodes you want to add.
                                          The data obtained as a result of Add Perturbation are used. 
                                          In the tab below, the prediction heatmap and PCA circle graphics that belong to the nodes you selected are created.")),hr(),
                            
                            width = 12,
                            fluidRow(column(2,
                                            
                                            br()),
                                     box(width = 12,
                                         uiOutput("squad_node"),
                                         status = "info", solidHeader = TRUE, collapsible = TRUE, closable = TRUE,
                                         title = HTML("Simulation matrix for the selected nodes", 
                                                      as.character(actionLink(inputId = "info11", 
                                                                              label = "", 
                                                                              icon = icon("info-circle", verify_fa = FALSE,"<font size='5'>"),))), 
                                         
                                         
                                         actionBttn(inputId = "squaddButton",
                                                    label = "Run",
                                                    color = "primary",
                                                    style = "gradient",
                                                    icon = icon("glyphicon glyphicon-play",lib = "glyphicon")
                                         ),
                                         br(),
                                         textOutput("text1"),
                                         tags$head(tags$style("#text1{color:#34495E;
                                 font-size: 20px;
                                 font-style: normal;
                                 position:static;
                                 font-weight:bold;
                                 text-align:center;
                                 
                                 }"
                                         )
                                         ),
                                         
                                         tableOutput("text2"),
                                         tags$head(tags$style(type="text/css","#text2{color:#34495E;
                                         white-space: break-spaces;
                                 font-size: 20px;
                                 font-style: normal;
                                 position:static;
                                 font-weight:bold;
                                 text-align:left;
                                 }"
                                         )
                                         ),
                                         
                                         br(),
                                         hr(), helpText(""), plotOutput('squadd1',width = "100%",height = "400px")),
                                     box(width = 6, id= "heatsquad",
                                         title = "Prediction heatmap", status = "success", solidHeader = TRUE, collapsible = TRUE, closable = TRUE,
                                         
                                         hr(), helpText(""), plotOutput('squadd2',width = "100%",height = "500px")),
                                     
                                     box(width = 6, id = "CCsquad",
                                         title = "Correlation circle", status = "success", solidHeader = TRUE, collapsible = TRUE, closable = TRUE,
                                         
                                         hr(), helpText(""), plotOutput('squadd3',width = "100%",height = "500px"))
                            ))))            
        #################################################################################################################################               
        
      ),
    ),
    tagList(
      tags$head(
        tags$style(
          ".main-header .navbar-custom-menu {
                                        float: left;
                                      }
                                   .sidebar-menu {
                                        display: flex;
                                   }"
        )
      )
    ),
    
    
    
  ))
  
)
