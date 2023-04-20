#########################################
#### PRO-simat Bioinformatics WebApp ####
#### Written by Rana Salihoglu       ####
#### Update:23_Oct_2022              ####
#### Ⓒ                              ####
#########################################
############################################## S  E  R  V  E  R  // P A R T // #######################################################
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#===   S  E  R  V  E  R   ===#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-




shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  
  #------------------------------------------------------------------------- G O _ A N A L Y S I S
  url <- a("eGOSynthetic", href="https://gosyn.bioapps.biozentrum.uni-wuerzburg.de/index.php",target="_blank")
  
  
  
  output$tab <- renderUI({
    tagList(url)
  })
  
  output$tab_string <- renderUI({
    tagList(url2)
  })
  
  observeEvent(input$help, {
    if (input$tabs == "Category_1") {
      rintrojs::introjs(session, options = list(
        steps = data.frame(element = c("#hostID","#pathogenID","#nodeShape","#layout","#runProtein","#MoReBtn","#LessBtn"),
                           intro = c(
                             "Select first Organism for which you want to see protein-protein interactions.",
                             "Select second Organism for which you want to see protein-protein interactions.",
                             "Select the node shape where your proteins will be displayed.",
                             "Select one of the layouts from the select input box to see the network.",
                             "Now you can run analysis with Submit button.",
                             "You can see more interaction.",
                             "You can see less interaction."
                           ))
      ))
    }
  })
  
  
  
  observeEvent(input$help2, {
    if (input$tabs == "Category_2") {
      introjs(session, options = list("showBullets"="true", "showProgress"="true", "showStepNumbers"="false","nextLabel"="Next","prevLabel"="Prev","skipLabel"="Skip"),events = list(onbeforechange = readCallback("switchTabs")))
      
      #rintrojs::introjs(session, options = list("nextLabel" = "Continue",
      #                                         "prevLabel" = "Back",
      #steps = data.frame(element = c("#userFile","#reset"),
      #                  intro = c(
      #c("#userFile","#geneColumn","#logFColumn","#pAdColumn","#submitButton"),
      #"Upload your Differential Gene Expression data in .csv format.",
      #"Select the columns with Gene in the data you have loaded.",
      #"Select the columns withLogFC values in the data you have loaded.",
      #"Select the columns with pAdjust values in the data you have loaded.",
      #                   "Run analysis",
      #                  "delete"
      #               ))
      #))
    }
  })
  
  observeEvent(input$help3, {
    if (input$tabs == "Category_3") {
      rintrojs::introjs(session, options = list(
        steps = data.frame(element = c("#kegg1","#kegg2","#PathviEw","#downloadPathviewPng"),
                           intro = c(
                             "The data obtained as a result of the KEGG analysis are shown here in tabular form.",
                             "You can access the Enrichment Map for enrichment result of overrepresentation analysis (ORA) in the KEGG enrichment emapplot tab.",
                             "According to the data obtained as a result of ORA, when you write the pathway ID to the input pathway text box, an enriched KEGG pathway will be created for your data.",
                             "You can download the colorized version of the pathway you have selected by pathview."
                             
                           ))
      ),events = list(onbeforechange = readCallback("switchTabs")))
      
    }
  })
  
  
  observeEvent(input$help4, {
    if (input$tabs == "Category_4") {
      rintrojs::introjs(session, options = list(
        steps = data.frame(element = c("#bxjim", "#runGraphml","#runJimena","#strt","#end","#val","#add_pert","#runJimena"),
                           intro = c(
                             "Your signaling network data must contain node1, label, node2 columns and must be in txt format. Upload the data you have prepared in the desired format.",
                             "Click the Convert Graphml button. This ensures that your data is optimized for the jimena dynamic simulation.",
                             "Then start the analysis by clicking the Run Jime.",
                             "A reactive node list is created and select the node you want to add from there and set start time",
                             "Set end time",
                             "Set the activation value.When you give a value of 1 to the node you want to add, it will be simulated as an overexpression, while that node will be knocked out with a value of 0.",
                             "Click the Add perturbation button to add the nodes you have determined to your perturbation list.",
                             "Click the run analysis button again and you can see the modified graph."
                             
                             
                           ))
      ))
      
    }
  })
  
  
  observeEvent(input$help5, {
    if (input$tabs == "Category_5") {
      rintrojs::introjs(session, options = list(
        steps = data.frame(element = c("#mainsq","#squad_node", "#squaddButton","#heatsquad","#CCsquad"),
                           intro = c(
                             "Please run Jimena before using the SQUADD part and identify the perturbation nodes you want to add. The data obtained as a result of Add Perturbation are used in here.",
                             "Select the nodes you want from the Select nodes box",
                             "Click the Run button.",
                             "The darker the color, the higher the increase of the node activation.",
                             "The correlation circle was used to assess the consistency of a model. If the angle between each vector is greater than 90 degrees, we can say that there is no correlation."
                           ))
      ))
      
    }
  })
  
  source('connection-mysql.R')
  
  observeEvent(input$reFresh, {
    refresh()
  })
  
  #output$keepAlive <- renderText({
  #req(input$count)
  #paste("keep alive ", input$count)
  #})
  
  
  ################################# QUERY FOR DATABASE // SR // ###################################################################
  v <- reactiveValues(data=NULL,
                      edgevis=NULL,
                      for_symbol=NULL
                      )
 
  
  ################################# Select Organisms >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  observeEvent(input$runProtein,{
    
    ##################### Input Multiple Protein PART #####################
    
    source("multipleInput.R",local = TRUE)
    
    ###################### v i s - n e t w o r k // SR // #####################
   
    observe({
      edges <-e$data
      
      if(nrow(edges) < 1) {
        Sys.sleep(5)
        ask_confirmation(
          "form3",
          title = " ",
          type = "info",
          btn_labels = c("Cancel", "OK"),
          allowEscapeKey = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          text =  
            div(HTML("Oops!", "No protein-protein interaction was found between the organisms you selected. Please choose other organisms. "
            ))
        )
        
      }
      
      
      else {
        
        colnames(edges) <- c("from","Gene_1", "link","tax_from","organism_1","to","Gene_2", "tax_to","organism_2")
        
        nodevistb <- edges
        
        colnames(nodevistb) <- c("From_A","Gene_A", "Interaction_Type","Taxonomy_A","Organism_A","To_B","Gene_B", "Taxonomy_B","Organism_B")
        nodevistb <- nodevistb[,c(1,6,3,2,7,4,8,5,9)]
        v$edgevis <- as.data.frame(nodevistb)
        
        Edges <- as.data.frame(edges$from)
        Edges$to <- edges$to
        Edges$link <- edges$link
        colnames(Edges) <- c("from","to","link")
        
        
        nodes1 <- as.data.frame(edges$from)
        nodes1$tax <- edges$tax_from
        nodes1$group <- edges$organism_1
        nodes1$genes <- edges$Gene_1
        colnames(nodes1) <- c("id","tax","group","genes")
        
        nodes2 <- as.data.frame(Edges$to)
        nodes2$tax <- edges$tax_to
        nodes2$group <- edges$organism_2
        nodes2$genes <- edges$Gene_2
        colnames(nodes2) <- c("id","tax","group","genes")
        
        Nodess <- rbind(nodes1,nodes2)
        Nodess <- distinct(Nodess)
        
        Nodess$id <- Nodess$id
        Nodess$label <- Nodess$id
        Nodess <- Nodess%>% mutate(font.size = ifelse(id=="GL068-",34,
                                                      ifelse(id=="GL121-",34,
                                                             ifelse(id=="GL245-",34,
                                                                    ifelse(id=="gusA",27,
                                                                           ifelse(id=="lacZ",27,
                                                                                  ifelse(id=="RUC-GFP",27,
                                                                                         16)))))))
        
        for_symbol <- NULL
        for_symbol$id <- as.data.frame(Nodess$id)
        for_symbol$gene <- as.data.frame(Nodess$genes)
        for_symbol <- as.data.frame(for_symbol)
        colnames(for_symbol) <- c("id","symbol")
        v$for_symbol <- for_symbol
        
        Nodess$title <- paste0("Gene :",Nodess$genes,
                               "<br>UniprotKB :",Nodess$id,
                               "<br>Species :",Nodess$group,
                               "<br>Taxonomy :",Nodess$tax)
        
        
        
        nodes <- as.data.frame(Nodess)
        #v$data<- as.data.frame(nodes)
        
        nodes <- subset(nodes, select = -c(genes))
        nodes <- nodes %>% group_by(id) %>% filter (! duplicated(id)) 
        
        edges <- as.data.frame(Edges)
        
        edges<-edges %>% mutate(color = ifelse(link=="activation","#5ebf87",
                                               ifelse(link=="inhibition","#e0392b",
                                                      ifelse(link=="MI:0194","#d0d676",
                                                             ifelse(link=="MI:0203","#d6b076",
                                                                    ifelse(link=="MI:0217","#8ec4b7",
                                                                           ifelse(link=="MI:0220","#8eb5c4",
                                                                                  ifelse(link=="MI:0403","#998ec4",
                                                                                         ifelse(link=="MI:0407","#bd8ec4",
                                                                                                ifelse(link=="MI:0414","#c48ea3",
                                                                                                       ifelse(link=="MI:0570","#d3815f",
                                                                                                              ifelse(link=="MI:0844","#bfd35f",
                                                                                                                     ifelse(link=="MI:0882","#76d35f",
                                                                                                                            ifelse(link=="MI:0914","#a994bf",
                                                                                                                                   ifelse(link=="MI:0915","#797e85",
                                                                                                                                          ifelse(link=="MI:0945","#5fcfd3",
                                                                                                                                                 ifelse(link=="MI:1126","#5f94d3",
                                                                                                                                                        ifelse(link=="MI:2364","#d35fcc",
                                                                                                                                                               "#9c95a2"))))))))))))))))))
        
        edges<-edges %>% mutate(arrows.to.type = ifelse(link=="activation","arrow",
                                                        ifelse(link=="inhibition","bar", "image")))
        edges<-edges %>% mutate(title = edges$link)
        edges <- distinct(edges)
        edges$id <- row.names(edges)
        
        
        ledges <-edges[, c("color", "title","arrows.to.type")]
        ledges <- distinct(ledges)
        
        colnames(ledges) <- c("color", "label","arrows")
        
        gene_lst <- reactiveVal(nodes)
        graph_data$nodes <- nodes
        graph_data$edges <- edges
        graph_data$ledges <-ledges
        v$data <- nodes
       
      } 
      
    })
    
    
    graph_data = reactiveValues(
      nodes = NULL,
      edges =NULL,
      ledges=NULL
    )
    if (input$layout == "")
      return(NULL) 
    output$editable_network <- renderVisNetwork({
      
      PPI_network <- reactive({
        visNetwork(graph_data$nodes, graph_data$edges,height = "800px", width = "120%") %>%visIgraphLayout(
          layout = input$layout)%>% 
          visGroups(groupname = "Arabidopsis thaliana", color =list(highlight ="#C0392B",background = "#E57373",border ="#E57373")) %>%
          visGroups(groupname = "Bacillus subtilis (strain 168)",color =list(highlight ="#C0392B",background = "#F06292",border ="#F06292")) %>%
          visGroups(groupname = "Bacillus subtilis",color =list(highlight ="#C0392B",background = "#D77DC1",border ="#D77DC1")) %>%
          visGroups(groupname = "Bos taurus",color =list(highlight ="#C0392B",background = "#BF7785",border ="#BF7785")) %>%
          visGroups(groupname = "Bacillus anthracis",color =list(highlight ="#C0392B",background = "#D77DC1",border ="#D77DC1")) %>%
          visGroups(groupname = "Caenorhabditis elegans", color =list(highlight ="#C0392B",background = "#BA68C8",border ="#BA68C8")) %>%
          visGroups(groupname = "Cavia porcellus", color =list(highlight ="#C0392B",background = "#9575CD",border ="#9575CD")) %>%
          visGroups(groupname = "Danio rerio", color =list(highlight ="#C0392B", background ="#64B5F6",border ="#64B5F6")) %>%
          visGroups(groupname = "Drosophila melanogaster", color =list(highlight ="#C0392B",background = "#4FC3F7",border ="#4FC3F7")) %>%
          visGroups(groupname = "Escherichia coli", color =list(highlight ="#C0392B",background = "#4DD0E1",border ="#4DD0E1")) %>%
          visGroups(groupname = "Escherichia coli O157:H7", color =list(highlight ="#C0392B",background = "#5FE3AD",border ="#5FE3AD")) %>%
          visGroups(groupname = "Escherichia coli O127:H6 (strain E2348/69 / EPEC)", color =list(highlight ="#C0392B",background = "#5FE3AD",border ="#5FE3AD")) %>%
          visGroups(groupname = "Escherichia coli (strain K12)", color =list(highlight ="#C0392B",background = "#47CFC3",border ="#47CFC3")) %>%
          visGroups(groupname = "Escherichia coli O157:H7 str. EDL933", color =list(highlight ="#C0392B",background = "#47CFC3",border ="#47CFC3")) %>%
          visGroups(groupname = "Escherichia coli O55:H7 (strain CB9615 / EPEC)", color =list(highlight ="#C0392B",background = "#47CFC3",border ="#47CFC3")) %>%
          visGroups(groupname = "Gallus gallus", color =list(highlight ="#C0392B",background = "#4DB6AC",border ="#4DB6AC")) %>%
          visGroups(groupname = "Homo sapiens", color =list(highlight ="#C0392B",background = "#7986CB",border ="#7986CB")) %>%
          visGroups(groupname = "Helicobacter pylori (strain ATCC 700392 / 26695)", color =list(highlight ="#C0392B",background = "#81C784",border ="#81C784")) %>%
          visGroups(groupname = "Helicobacter pylori", color =list(highlight ="#C0392B",background = "#81C784",border ="#81C784")) %>%
          visGroups(groupname = "Listeria monocytogenes serovar 1/2a (strain ATCC BAA-679 / EGD-e)", color =list(highlight ="#C0392B",background = "#DCE775",border ="#DCE775")) %>%
          visGroups(groupname = "Listeria monocytogenes serotype 1/2a", color =list(highlight ="#C0392B",background = "#6DAF63",border ="#6DAF63")) %>%
          visGroups(groupname = "Listeria monocytogenes serotype 4b str. LL195", color =list(highlight ="#C0392B",background = "#63AF8E",border ="#63AF8E")) %>%
          visGroups(groupname = "Listeria monocytogenes serotype 4b (strain F2365)", color =list(highlight ="#C0392B",background = "#47A08C",border ="#47A08C")) %>%
          visGroups(groupname = "Listeria monocytogenes", color =list(highlight ="#C0392B",background = "#DCE775",border ="#DCE775")) %>%
          visGroups(groupname = "Mus musculus", color =list(highlight ="#C0392B",background = "#FFF176",border ="#FFF176")) %>%
          visGroups(groupname = "Oryzias latipes", color =list(highlight ="#C0392B",background = "#FFD54F",border ="#FFD54F")) %>%
          visGroups(groupname = "Ovis aries", color =list(highlight ="#C0392B",background = "#FFA726",border ="#FFA726")) %>%
          visGroups(groupname = "Plasmodium falciparum (isolate 3D7)", color =list(highlight ="#C0392B",background = "#FF8A65",border ="#FF8A65")) %>%
          visGroups(groupname = "Plasmodium falciparum", color =list(highlight ="#C0392B",background = "#A3E4D7",border ="#A3E4D7")) %>%
          visGroups(groupname = "Rattus norvegicus", color =list(highlight ="#C0392B",background = "#D32F2F",border ="#D32F2F")) %>%
          visGroups(groupname = "Saccharomyces cerevisiae", color =list(highlight ="#C0392B",background = "#C2185B",border ="#C2185B")) %>%
          visGroups(groupname = "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", color =list(highlight ="#C0392B",background = "#C87CBC",border ="#C87CBC")) %>%
          visGroups(groupname = "Schizosaccharomyces pombe (strain 972 / ATCC 24843)", color =list(highlight ="#C0392B",background = "#9765AD",border ="#9765AD")) %>%
          visGroups(groupname = "Staphylococcus aureus", color =list(highlight ="#C0392B",background = "#7878D0",border ="#7878D0")) %>%
          visGroups(groupname = "Staphylococcus aureus (strain NCTC 8325 / PS 47)", color =list(highlight ="#C0392B",background = "#6565AD",border ="#6565AD")) %>%
          visGroups(groupname = "Staphylococcus aureus (strain N315)", color =list(highlight ="#C0392B",background = "#8484b5",border ="#8484b5")) %>%
          visGroups(groupname = "Treponema pallidum", color =list(highlight ="#C0392B",background = "#1976D2",border ="#1976D2")) %>%
          visGroups(groupname = "Treponema pallidum (strain Nichols)", color =list(highlight ="#C0392B",background = "#78A6D0",border ="#78A6D0")) %>%
          visGroups(groupname = "Vaccinia virus", color =list(highlight ="#C0392B",background = "#0097A7",border ="#0097A7")) %>%
          visGroups(groupname = "Vaccinia virus Copenhagen", color =list(highlight ="#C0392B",background = "#00695C",border ="#00695C")) %>%
          visGroups(groupname = "Vaccinia virus WR", color =list(highlight ="#C0392B",background = "#AFB42B",border ="#AFB42B")) %>%
          visGroups(groupname = "Vaccinia virus L-IPV", color =list(highlight ="#C0392B",background = "#01579B", border ="#01579B")) %>%
          visGroups(groupname = "Vaccinia virus Ankara", color =list(highlight ="#C0392B",background = "#E65100", border ="#E65100")) %>%
          visGroups(groupname = "Vaccinia virus GLV-1h68", color =list(highlight ="#C0392B",background = "#0066FF", border ="#0066FF")) %>%
          visGroups(groupname = "Xenopus laevis", color =list(highlight ="#C0392B",background = "#0066FF", border ="#0066FF")) %>%
          visGroups(groupname = "insertion", color = "#db3058", shape = "triangle") %>%
          visGroups(groupname = "Deletion", color = "#db3058", shape = "triangleDown") %>%
          visNodes(shape = input$nodeShape, size=45, shadow = TRUE)%>%
          visEdges(arrows = "to",label = graph_data$edges$link, width =4, shadow = TRUE, smooth = list(enabled = TRUE, roundness= 0.1)) %>%
          visInteraction(multiselect = TRUE, navigationButtons = TRUE, hover = TRUE)%>%
          visExport(type = "pdf")%>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_selection', nodes.nodes);
                ;}")%>%
          
          visLegend(addEdges = graph_data$ledges, main = list(text = "Legend",
                                                              style = "font-family:Comic Sans MS;color: #212f3c;font-size:20px;text-align:center;"),
                    
                    width = 0.20, stepY= 53, position = "right")%>%
          visOptions(manipulation = T)
      })
      
      output$downloadNetwork <- downloadHandler(
        filename = function() {
          paste('network-', Sys.Date(), '.html', sep='')
        },
        content = function(con) {
          PPI_network() %>% visSave(con)
        })
      
      PPI_network() 
    })
    
    
    ###########
    observeEvent(input$editable_network_graphChange, {
      
      # If the user added a node, add it to the data frame of nodes.
      if(input$editable_network_graphChange$cmd == "addNode") {
        if(input$hostID == '562,83333,83334,574521,155864' || input$hostID == '3702' ||input$pathogenID == '562,83333,83334,574521,155864' || input$pathogenID == '3702' ){
          ask_confirmation(
            "form3",
            title = " ",
            type = "info",
            btn_labels = c("Cancel", "OK"),
            allowEscapeKey = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            text =  
              div(HTML("Enter the Symbol of your protein in the box to write Node id and click save. 
                       Please make sure you type the protein ids as specified on the <a href='https://www.uniprot.org/' target='_blank'>Uniprot web page.</a>
                       If you wrote it correctly, ignore the warning."
              ),
              ),
            
          )
        } else {
          ask_confirmation(
            "form3",
            title = " ",
            type = "info",
            btn_labels = c("Cancel", "OK"),
            allowEscapeKey = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            text =  
              div(HTML("Enter the UniprotID of your protein in the box to write Node id and click save. 
                       Please make sure you type the protein ids as specified on the <a href='https://www.uniprot.org/' target='_blank'>Uniprot web page.</a>
                         If you wrote it correctly, ignore the warning."
              ),
              ),
            
          )
        }
        
        
        temp = bind_rows(
          graph_data$nodes,
          data.frame(id = input$editable_network_graphChange$id,
                     label = input$editable_network_graphChange$label,
                     stringsAsFactors = F)
        )
        graph_data$nodes = temp
        
       
      }
      # If the user added an edge, add it to the data frame of edges.
      else if(input$editable_network_graphChange$cmd == "addEdge") {
        temp = bind_rows(
          graph_data$edges,
          data.frame(id = input$editable_network_graphChange$id,
                     from = input$editable_network_graphChange$from,
                     to = input$editable_network_graphChange$to,
                     stringsAsFactors = F)
        )
        graph_data$edges = temp
      }
      # If the user edited a node, update that record.
      else if(input$editable_network_graphChange$cmd == "editNode") {
        temp = graph_data$nodes
        temp$label[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$label
        graph_data$nodes = temp
      }
      # If the user edited an edge, update that record.
      else if(input$editable_network_graphChange$cmd == "editEdge") {
        temp = graph_data$edges
        temp$from[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$from
        temp$to[temp$id == input$editable_network_graphChange$id] = input$editable_network_graphChange$to
        graph_data$edges = temp
      }
      # If the user deleted something, remove those records.
      else if(input$editable_network_graphChange$cmd == "deleteElements") {
        for(node.id in input$editable_network_graphChange$nodes) {
          temp = graph_data$nodes
          temp = temp[temp$id != node.id,]
          graph_data$nodes = temp
        }
        for(edge.id in input$editable_network_graphChange$edges) {
          temp = graph_data$edges
          temp = temp[temp$id != edge.id,]
          graph_data$edges = temp
        }
      }
    })
    
    #############
    
      vis_cho <- reactive({
        
        vis_node_list <- as.data.frame(graph_data$nodes)
        vis_node <- as.data.frame(vis_node_list$id)
        #colnames(vis_node) <- "Node List"
        
        v$updateDat <- vis_node
        return(vis_node)
      })
      
    
    
    observe({
      updateSelectizeInput(session, 'SearchVisNode', choices = vis_cho())
    })
    
    observeEvent(input$btn1,{
      if(input$SearchVisNode > 0){
        isolate({
          nodes_visSearch <- as.data.frame(v$data)
          current_vis_node <- nodes_visSearch[grep(input$SearchVisNode, nodes_visSearch$label), "id"]
          visNetworkProxy("editable_network") %>% visSelectNodes(id  = current_vis_node)
        })
      }
    })
    
    #, smooth = list(enabled = TRUE, roundness= 0.1)
    # display go info by selected nodes
    observeEvent(input$current_node_selection, {
      isolate({
        
        query<-paste0("
              SELECT DISTINCT
              IGOSynthetic.UniProtKB2GO.GeneInfo_ID  AS Uniprotkb,
              IGOSynthetic.UniProtKB2GO.GOInfo_ID AS GO_id,
              IGOSynthetic.GOInfo.Definition AS Definition
              FROM IGOSynthetic.UniProtKB2GO
              INNER JOIN IGOSynthetic.GOInfo ON IGOSynthetic.UniProtKB2GO.GOInfo_ID = IGOSynthetic.GOInfo.GO_term
              WHERE IGOSynthetic.UniProtKB2GO.GeneInfo_ID = '",input$current_node_selection,"';")
        
        
        qmain <- dbSendQuery(con(), query)
        
        go_info <- print(dbFetch(qmain, n=-1))
        
        showModal(modalDialog(
          output$goInfo <- renderDT(
            go_info
          ),
          easyClose = TRUE
        ))
       
      })
      
      
      # get position info
      observeEvent(input$store_position, {
        visNetworkProxy("editable_network") %>% visGetPositions()
      })
      
      # format positions
      nodes_positions <- reactive({
        positions <- input$network_positions
        if(!is.null(positions)){
          nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
          nodes_positions$id <- names(positions)
          nodes_positions
        } else {
          NULL
        }
      })
      
    })
    
  })
  
  
  ###### Vis Network Data Table ##########################
  
  output$table_int <- DT::renderDataTable(server = FALSE,{
    DT::datatable(v$edgevis,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  options = list(
                    dom = 'Bfrtip',
                    lengthChange = FALSE,
                    deferRender = TRUE,
                    ordering = TRUE,
                    scrollY = 400,
                    scroller = TRUE,
                    buttons = 
                      list( list(
                        extend = 'collection',
                        
                        buttons = list(list(extend='csv',
                                            filename = 'PRO-simat_interaction_table'),
                                       list(extend='excel',
                                            filename = 'PRO-simat_interaction_table'),
                                       list(extend='pdf',
                                            filename= 'PRO-simat_interaction_table')),
                        text = '<span class="glyphicon glyphicon-download-alt"></span>'
                      )),
                    
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#0080F1', 'color': '#fff'});",
                      "}")))})
  
  ######### GO and KEGG analysis PART ##################################
  ############### GO analysis info part 
  
  observeEvent(input$info1, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("For GO and KEGG analysis, you can upload your own DEGs file in csv format or select sample LUAD DEGs data.
Specify the desired columns of the csv file you have uploaded/selected from the select input box. 
                 For example, select the column with gen symbols from the 'select gene column' tab. 
                 Repeat this procedure for all desired columns.
                 Your genes are filtered according to the Log2FC and padjust treshold you will determine and 
                 the filtered genes are used in the analysis."
        ))
    )
  })
  
  observeEvent(input$info2, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("In order to start the GO and KEGG analyzes, fill in the desired parameters in accordance with your data. 
                 Find and select the organism to which your genes belong and the id type to which the genes are given from the 'keytpe'."
        ))
    )
  })
  
  observeEvent(input$info12, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("The result of GO enrichment analysis is presented as cnetplot.
Results are shown for genes associated with apoptosis, proliferation, negative and positive regulation, and other functionalities."
        ))
    )
  })
  
  observeEvent(input$info3, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("If you have performed GO and KEGG enrichment analysis by uploading your DEGs data, you can see the volcanoplot created for your data according to the log2FC and padjust treshold you have determined.
However, if you have performed the analysis using the genes in the Prosimat Network analysis tab, the volcano plot will not be created."
        ))
    )
  })
  
  observeEvent(input$info4, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("Examine the graphics created as a result of the GO enrichment analysis. You can see more GO terms on the chart by increasing the category size."
        ))
    )
  })
  
  observeEvent(input$info5, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("The KEGG analysis result is displayed in tabular form. You can see the colored pathway using the 'ID' found here in the 'Pathway' tab."
        ))
    )
  })
  
  observeEvent(input$info6, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("Emapplot arranges the enriched terms in a network with edges connecting overlapping gene clusters.
You can see more or less pathways by changing the category size."
        ))
    )
  })
  
  observeEvent(input$info7, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("
Start the GO and KEGG analysis using your DEGs data to obtain pathway images colored using the Log2FC' value.
PRO-simat only provides log2FC values for human with PDAAD DEGs data obtained from Gepia2 database.
                 To see the colored pathway you should use the 'ID' in the 'Output Table' where the result of the Kegg analysis is presented. No pathway can be displayed when you use an id that is not in the table."
        ))
    )
  })
  
  observeEvent(input$info8, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("
Upload your signaling network file in txt format containing node1,label,node2 column names.
Click the 'Convert Graphml' button to convert it to Graphml format.
Press the Run Jimena Button to start the Jimena simulation analysis."
        ))
    )
  })
  
  observeEvent(input$info9, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("
Select the node you want to add and determine the time interval when the node is activated or inhibited in the simulation and click the 'Add perturbation' button.
(start : start time, end: end time, value:1-0 (activation/inhibition)).
You can also give values like 0.5 so that it can be partially active. Then run the simulation again by clicking the 'Run Jimena' button. "
        ))
    )
  })
  
  observeEvent(input$info10, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("
Start the GO and KEGG analysis using your DEGs data to obtain pathway images colored using the Log2FC' value.
PRO-simat only provides log2FC values for human with PDAAD DEGs data obtained from Gepia2 database.
                 To see the colored pathway you should use the 'ID' in the 'Output Table' where the result of the Kegg analysis is presented. No pathway can be displayed when you use an id that is not in the table."
        ))
    )
  })
  
  
  observeEvent(input$info11, {
    ask_confirmation(
      "form3",
      title = " ",
      type = "info",
      btn_labels = c("Cancel", "OK"),
      allowEscapeKey = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      text =  
        div(HTML("
For each added perturbation you use in the Jimena tab, your outputs are kept in the program and used in the Squadd visualization. Graphs are created for the genes you will select for each of your analysis results indexed below.
Select the genes of interest from the Select node box.
Thus, it is graphed how the genes you selected vary in different perturbation results."
        ))
    )
  })
  
  
  ################# select data form (example or own csv)
  
  #observeEvent(input$chosen_data,{
  # if("example_csv" == input$chosen_data){
  #  shinyjs::show("collapse1")
  # shinyjs::hide("collapse2")
  #}
  #})
  observeEvent(input$chosen_data,{
    if("example_csv" == input$chosen_data){
      #shinyjs::toggle(id = "advanced", anim = TRUE)
      shinyjs::show(id = "advanced", anim = TRUE)
    }
    if("user_csv" == input$chosen_data){
      #shinyjs::toggle(id = "advanced", anim = TRUE)
      shinyjs::show(id = "advanced", anim = TRUE)
    }
    if("PPI_data" == input$chosen_data){
      #shinyjs::toggle(id = "advanced", anim = TRUE)
      shinyjs::hide(id = "advanced", anim = TRUE)
    }
  })
  ############# If Input table from user ##############    
  
  #Reactive to store loaded data
  reactives <- reactiveValues(
    
    userData = NULL
    
  )
  
  #Observe file being selected
  
  observeEvent(input$userFile, {
    
    #Store loaded data in reactive
    reactives$userData <- read.csv(file = input$userFile$datapath)
    
    #Update select input
    updateSelectInput(session, inputId = 'geneColumn', label = 'Select Gene Column', choices  = colnames(reactives$userData))
    updateSelectInput(session, inputId = 'log2FColumn', label = 'Select log2FC Column', choices  = colnames(reactives$userData))
    updateSelectInput(session, inputId = 'pAdColumn', label = 'Select pAdj Column', choices  = colnames(reactives$userData))
    
  })
  
  #### LUAD example files ######################################################
  
  
  #####################
  
  #Data table
  
  output$whichosen <-  DT::renderDataTable(server=FALSE,{
    if (input$chosen_data == "user_csv" && !is.null(reactives$userData)) {
      
      DT::datatable(reactives$userData,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
              
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 535,
                      scroller = TRUE
                     ))
    }
    else if(input$chosen_data == "example_csv" && is.null(reactives$userData)) {
      reactives$exampData <- read.csv("LUAD_gepia2.csv")
      
      #Update select input
      updateSelectInput(session, inputId = 'geneColumn', label = 'Select Gene Column', choices  = colnames(reactives$exampData))
      updateSelectInput(session, inputId = 'log2FColumn', label = 'Select log2FC Column', choices  = colnames(reactives$exampData))
      updateSelectInput(session, inputId = 'pAdColumn', label = 'Select pAdj Column', choices  = colnames(reactives$exampData))
      
      DT::datatable(reactives$exampData,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
                      
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 535,
                      scroller = TRUE
                      ))
    }
    
  })
  
  
  observeEvent(input$reset, {
    reactives$userData <- NULL
  })
  ########################___________________________________________________________________________________________________________________________
  prosim <- reactiveValues()
  
  observe({
    GoKegg()
  })
  
  GO_speciesDb <- c("",
                    "Arabidopsis thaliana" = "org.At.tair.db",
                    "Anopheles gambiae" = "org.Ag.eg.db",
                    "Bos taurus" = "org.Bt.eg.db",
                    "Canis familiaris" = "org.Cf.eg.db",
                    "Chaenorhabditilis elegans" = "org.Ce.eg.db",
                    "Danio rerio" = "org.Dr.eg.db",
                    "Drosophla melanogaster" = "org.Dm.eg.db",
                    "Escherishia coli" = "org.EcK12.eg.db",
                    "Escherishia coli strain Sakai" = "org.EcSakai.eg.db",
                    "Gallus gallus" = "org.Gg.eg.db",
                    "Mus musculus" = "org.Mm.eg.db",
                    "Rattus norvegicus" = "org.Rn.eg.db",
                    "Homo sapiens" = "org.Hs.eg.db",
                    "Macaca mulatta" = "org.Mmu.eg.db",
                    "Myxococcus xanthus" = "org.Mxanthus.db",
                    "Pan troglodytes" = "org.Pt.eg.db",
                    #"Plasmodium falciparum" = "org.Pf.plasmo.db",
                    "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
                    "Sus scrofa" = "org.Ss.eg.db",
                    "Xenopus laevis" = "org.Xl.eg.db"
  )
  
  
  updateSelectInput(session, "GO_species", choices = GO_speciesDb)
  observeEvent(input$GO_species,{
    if(input$GO_species == "")
      return(NULL)
    
    library(input$GO_species, character.only = T)
    
    annDb = eval(parse(text = input$GO_species))
    keytypes = keytypes(annDb)
    updateSelectInput(session, "keytype", choices = keytypes, selected = "GENENAME")
    
  })
  
  
  userSigDF <- reactiveValues()
  GoKegg <- eventReactive(input$submitButton,{
    Sys.sleep(5)
    if(!is.null(reactives$userData) && input$logFCut == FALSE && input$pAdjCut == FALSE){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please set the pAdjust and min log2FC values."
          ))
      )
      
      
    } else if(input$GO_species <0 || input$keytype <0 || input$KEGG_species <0 || input$pAdMeth <0){
      ask_confirmation(
        "form1",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please make sure you choose the parameters correctly. Remember to select the organism for GO and KEGG analyses. The keytype must be the same as the id type of your genes."
          ))
      )
      
      
    }
    
    else {
      
      withProgress(message = "Running ...",{
        
        if (input$chosen_data =="user_csv"  && input$chosen_data !="PPI_data" && input$chosen_data !="example_csv"){
          
          user_df <- as.data.frame(reactives$userData)
          
          logFc <- (user_df[input$log2FColumn])
          GeneCol <- (user_df[input$geneColumn])
          pAdjCol <- (user_df[input$pAdColumn])
          dForConv <-data.frame(GeneCol,logFc, pAdjCol)
          colnames(dForConv) <- c("geneSymbol", "Log2FC","pAdj")
          
          dForConv = dForConv[!duplicated(dForConv$geneSymbol),]
          userSigDF$dForConv <- dForConv
          
          shiny::reactiveConsole(TRUE)
          spsComps::shinyCatch({
            if (nrow(suppressWarnings(ids <- bitr(dForConv$geneSymbol, fromType=input$keytype, toType=c("ENTREZID"), OrgDb=input$GO_species)))<1){
              stop("Select correct keytpe!")
              
            }else{
              message("Everything is CORRECT :)")
              
            }
            
          }, blocking_level = "error")
          
          
          ids <- as.data.frame(ids)
          
          colnames(ids) <- c("geneSymbol","entrezid")
          UserGenes <- merge(ids, dForConv, by="geneSymbol")
          
          UserGenes <- UserGenes[!duplicated(UserGenes$geneSymbol),]
          
          UserGenes <- subset(UserGenes, select = -c(geneSymbol))
          
          userGeneList <- UserGenes[["Log2FC"]] 
          names(userGeneList) <- UserGenes[["entrezid"]]
          GeneList<-na.omit(userGeneList)
          GeneList <- sort(GeneList, decreasing = TRUE)
          
          sig_genes_df = subset(UserGenes, pAdj < input$pAdjCut)
          sig_genes_df = na.omit(sig_genes_df)
          
          geneSig <- sig_genes_df$Log2FC
          
          names(geneSig) <- sig_genes_df$entrezid
          
          geneSig <- na.omit(geneSig)
          #geneSig <- names(geneSig)[abs(geneSig) > input$logFCut]
          geneSig <- geneSig[abs(geneSig) > input$logFCut]
          
          userSigDF$SiGenesList <- names(geneSig)
          userSigDF$SiGeneListFC <- geneSig
         
          #userSigDF$GenesList <-GeneList
          
          if (input$GO_species == "org.Sc.sgd.db") {
            readableStat <- "FALSE"
           
          }
          
          else{
            readableStat <- "TRUE"
            
          }
        }
        
        ############ example data 
        else if (is.null(reactives$userData) && input$chosen_data !="PPI_data" && input$chosen_data =="example_csv"){
          
          user_df <- as.data.frame(reactives$exampData)
          
          logFc <- (user_df[input$log2FColumn])
          GeneCol <- (user_df[input$geneColumn])
          pAdjCol <- (user_df[input$pAdColumn])
          dForConv <-data.frame(GeneCol,logFc, pAdjCol)
          colnames(dForConv) <- c("geneSymbol", "Log2FC","pAdj")
          
          
          dForConv = dForConv[!duplicated(dForConv$geneSymbol),]
          userSigDF$dForConv <- dForConv
          
          shiny::reactiveConsole(TRUE)
          spsComps::shinyCatch({
            if (nrow(suppressWarnings(ids <- bitr(dForConv$geneSymbol, fromType=input$keytype, toType=c("ENTREZID"), OrgDb=input$GO_species)))<1){
              stop("Select correct keytpe!")
              
            }else{
              message("Everything is CORRECT :)")
              
            }
            
          }, blocking_level = "error")
          
          
          ids <- as.data.frame(ids)
          
          colnames(ids) <- c("geneSymbol","entrezid")
          UserGenes <- merge(ids, dForConv, by="geneSymbol")
          
          UserGenes <- UserGenes[!duplicated(UserGenes$geneSymbol),]
          
          UserGenes <- subset(UserGenes, select = -c(geneSymbol))
          
          userGeneList <- UserGenes[["Log2FC"]] 
          names(userGeneList) <- UserGenes[["entrezid"]]
          GeneList<-na.omit(userGeneList)
          GeneList <- sort(GeneList, decreasing = TRUE)
          
          sig_genes_df = subset(UserGenes, pAdj < input$pAdjCut)
          sig_genes_df = na.omit(sig_genes_df)
          
          geneSig <- sig_genes_df$Log2FC
          
          names(geneSig) <- sig_genes_df$entrezid
          
          geneSig <- na.omit(geneSig)
          #geneSig <- names(geneSig)[abs(geneSig) > input$logFCut]
          geneSig <- geneSig[abs(geneSig) > input$logFCut]
          userSigDF$SiGenesList <- names(geneSig)
          userSigDF$SiGeneListFC <- geneSig
          
          
          #userSigDF$GenesList <-GeneList
          
          if (input$GO_species == "org.Sc.sgd.db") {
            readableStat <- "FALSE"
            
          }
          
          else{
            readableStat <- "TRUE"
           
          }
        }
################################################################################
        ##### data coming PPI
        
        else{
          
          clean_data <- v$for_symbol
          genes_uniprot <- v$updateDat
          colnames(genes_uniprot) <- "id"
          
          
          
          
# ------------------------------------------------------------------------------
          
          
####################### PPI Human Data for analysis ############################
          
          if(input$chosen_data =="PPI_data" && input$GO_species == "org.Hs.eg.db" && input$chosen_data !="user_csv" && input$chosen_data !="example_csv"){
           
            ids <- bitr(genes_uniprot$id, fromType="UNIPROT", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
            
            readableStat <- "TRUE"
            
            PDAAD <- read.csv("PDAAD_Gepia2.csv")
            selected_PDAAD <- merge(PDAAD,ids,by="SYMBOL")
            
            if(nrow(selected_PDAAD) <1){
              
              geneSig <-ids$ENTREZID
              names(geneSig) <- ids$ENTREZID
              
              geneSig <- na.omit(geneSig)
              userSigDF$SiGeneListFC <- geneSig
              userSigDF$SiGeneList <- names(geneSig)
              
            }
            else{
             
              geneSig <- selected_PDAAD$Log2FC
              names(geneSig) <- selected_PDAAD$ENTREZID
              
              geneSig <- na.omit(geneSig)
              userSigDF$SiGeneListFC <- geneSig
              userSigDF$SiGeneList <- names(geneSig)
            }
          }
          
          ################### PPI  Data for Analysis
          if (input$GO_species == "org.At.tair.db" || input$GO_species =="org.EcK12.eg.db" ||input$GO_species == "org.EcSakai.eg.db") {
            uniques <- genes_uniprot[!genes_uniprot$id %in% clean_data$id,]
            uniques <- as.data.frame(uniques)
            colnames(uniques) <- "symbol"
            symbols <- clean_data$symbol
            symbols <- as.data.frame(symbols)
            colnames(symbols) <- "symbol"
            genes_symbol <-rbind(symbols, uniques)
            genes_symbol<- genes_symbol[!apply(genes_symbol == "", 1, all),]
            
            ids <- bitr(genes_symbol, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb= input$GO_species)
            readableStat <- "TRUE"
            
            
          }
          
          if(input$GO_species != "org.Hs.eg.db" && input$GO_species != "org.At.tair.db" && input$GO_species != "org.EcK12.eg.db" && input$GO_species != "org.EcSakai.eg.db" ){
            ids <- bitr(genes_uniprot, fromType="UNIPROT", toType=c("ENTREZID"), OrgDb=input$GO_species)
            readableStat <- "TRUE"
           
            
          }
          
          genes <- ids[[2]]
          genes <- as.data.frame(genes)
          colnames(genes) <- c("gene_name")
          
          row.names(genes) <- NULL
          
          userSigDF$SiGenesList <- genes$gene_name
          userSigDF$SiGenesListFC <- genes$gene_name
         
        }
        
        Ex <- readableStat
        
        
        ego_all <- enrichGO(gene = userSigDF$SiGenesList,
                            OrgDb         = input$GO_species,
                            ont           = "ALL",
                            maxGSSize     = 500,
                            minGSSize = 3,
                            pAdjustMethod = input$pAdMeth,
                            pvalueCutoff  = input$pvalCut,
                            qvalueCutoff  = input$qvalCut,
                            readable = Ex
        )
        
        ego_all_df <- as.data.frame(ego_all)
        row.names(ego_all_df) <- NULL
        
        if(nrow(ego_all_df) < 1) {
          
          Sys.sleep(5)
          ask_confirmation(
            "form3",
            title = " ",
            type = "info",
            btn_labels = c("Cancel", "OK"),
            allowEscapeKey = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            text =  
              div(HTML("Oops!", "You have determined a threshold that is not suitable for your data. Please enter different threshold values and run the analysis again."
              ))
          )
          
        }
        
        shinyCatch({
          if (nrow(go2 <- pairwise_termsim(ego_all))>1){message("Correct")}
          else {stop("Select correct keytpe!")}
        }, blocking_level = "error")
        
        
        
        ego_cc <- enrichGO(gene          = userSigDF$SiGenesList,
                           OrgDb         = input$GO_species,
                           ont           = "CC",
                           maxGSSize     = 500,
                           minGSSize = 3,
                           pAdjustMethod = input$pAdMeth,
                           pvalueCutoff  = input$pvalCut,
                           qvalueCutoff  = input$qvalCut,
                           readable      = Ex)
        
        ego_cc_df <- data.frame(ego_cc)
        row.names(ego_cc_df) <- NULL
        
        ego_mf <- enrichGO(gene          = userSigDF$SiGenesList,
                           OrgDb         = input$GO_species,
                           ont           = "MF",
                           maxGSSize     = 500,
                           minGSSize = 3,
                           pAdjustMethod = input$pAdMeth,
                           pvalueCutoff  = input$pvalCut,
                           qvalueCutoff  = input$qvalCut,
                           readable      = Ex)
        ego_mf_df <- data.frame(ego_mf)
        row.names(ego_mf_df) <- NULL
        
        
        ego_bp <- enrichGO(gene          = userSigDF$SiGenesList,
                           OrgDb         = input$GO_species,
                           ont           = "BP",
                           maxGSSize     = 500,
                           minGSSize = 3,
                           pAdjustMethod = input$pAdMeth,
                           pvalueCutoff  = input$pvalCut,
                           qvalueCutoff  = input$qvalCut,
                           readable      = Ex)
        
        
        ego_bp_df <- data.frame(ego_bp)
        row.names(ego_bp_df) <- NULL
        
      })
      
      
      
      
      output$whichtable <-  DT::renderDataTable(server=FALSE,{
        if (input$selecTable == "table1") { 
          DT::datatable(ego_bp_df,
                        extensions = c('Buttons','Responsive','Scroller'), 
                        options = list(
                          dom = 'Bfrtip',
                          lengthChange = FALSE,
                          deferRender = TRUE,
                          ordering = TRUE,
                          scrollY = 400,
                          scroller = TRUE,
                          buttons = 
                            list( list(
                              extend = 'collection',
                              
                              buttons = list(list(extend='csv',
                                                  filename = 'PRO-simat_GO_BP_table'),
                                             list(extend='excel',
                                                  filename = 'PRO-simat_GO_BP_table'),
                                             list(extend='pdf',
                                                  filename= 'PRO-simat_GO_BP_table')),
                              text = '<span class="glyphicon glyphicon-download-alt"></span>'
                            )),
                          
                          initComplete = JS(
                            "function(settings, json) {",
                            "$(this.api().table().header()).css({'background-color': '#3498DB', 'color': '#fff'});",
                            "}")))
          
          
        } else if  (input$selecTable == "table2") {
          
          DT::datatable(ego_mf_df, 
                        extensions = c('Buttons','Responsive','Scroller'), 
                        options = list(
                          dom = 'Bfrtip',
                          lengthChange = FALSE,
                          deferRender = TRUE,
                          ordering = TRUE,
                          scrollY = 400,
                          scroller = TRUE,
                          buttons = 
                            list( list(
                              extend = 'collection',
                              
                              buttons = list(list(extend='csv',
                                                  filename = 'PRO-simat_GO_MF_table'),
                                             list(extend='excel',
                                                  filename = 'PRO-simat_GO_MF_table'),
                                             list(extend='pdf',
                                                  filename= 'PRO-simat_GO_MF_table')),
                              text = '<span class="glyphicon glyphicon-download-alt"></span>'
                            )),
                          
                          initComplete = JS(
                            "function(settings, json) {",
                            "$(this.api().table().header()).css({'background-color': '#1ABC9C', 'color': '#fff'});",
                            "}")))
          
        } else if  (input$selecTable == "table3") {
          DT::datatable(ego_cc_df, 
                        extensions = c('Buttons','Responsive','Scroller'), 
                        options = list(
                          dom = 'Bfrtip',
                          lengthChange = FALSE,
                          deferRender = TRUE,
                          ordering = TRUE,
                          scrollY = 400,
                          scroller = TRUE,
                          buttons = 
                            list( list(
                              extend = 'collection',
                              
                              buttons = list(list(extend='csv',
                                                  filename = 'PRO-simat_GO_CC_table'),
                                             list(extend='excel',
                                                  filename = 'PRO-simat_GO_CC_table'),
                                             list(extend='pdf',
                                                  filename= 'PRO-simat_GO_CC_table')),
                              text = '<span class="glyphicon glyphicon-download-alt"></span>'
                            )),
                          
                          initComplete = JS(
                            "function(settings, json) {",
                            "$(this.api().table().header()).css({'background-color': '#34495E', 'color': '#fff'});",
                            "}")))
          
          
        } else if  (input$selecTable == "table4") {
          
          DT::datatable(ego_all_df, 
                        extensions = c('Buttons','Responsive','Scroller'), 
                        options = list(
                          dom = 'Bfrtip',
                          lengthChange = FALSE,
                          deferRender = TRUE,
                          ordering = TRUE,
                          scrollY = 400,
                          scroller = TRUE,
                          buttons = 
                            list( list(
                              extend = 'collection',
                              
                              buttons = list(list(extend='csv',
                                                  filename = 'PRO-simat_GO_all_table'),
                                             list(extend='excel',
                                                  filename = 'PRO-simat_GO_all_table'),
                                             list(extend='pdf',
                                                  filename= 'PRO-simat_GO_all_table')),
                              text = '<span class="glyphicon glyphicon-download-alt"></span>'
                            )),
                          
                          initComplete = JS(
                            "function(settings, json) {",
                            "$(this.api().table().header()).css({'background-color': '#E74C3C', 'color': '#fff'});",
                            "}")))
          
          
          
        }
        
      })
      
      
      if(nrow(ego_all_df) < 1) {
        ask_confirmation(
          "form2",
          title = " ",
          type = "info",
          btn_labels = c("Cancel", "OK"),
          allowEscapeKey = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          text =  
            div(HTML("Oops!", "Please change your parameters! "
            ))
        )
        
      }
      
      else {
        download_plots <- reactiveValues()
        output$whichplot2 <- renderSvgPanZoom({
          
          if (input$graphType == "plot1") { 
            withProgress(message = "Plotting CC barplot ...",{
              
              bar1 <- barplot(ego_cc, showCategory=input$categorySize)
              download_plots$images <- bar1
              svgPanZoom(svglite:::stringSVG(show(bar1), width = 19,
                                             height = 6,pointsize = 16
              ), controlIconsEnabled = T)
              
            })
            
          } else if (input$graphType =="plot2") {
            withProgress(message = "Plotting CC dotplot ...",{
              dot1 <- dotplot(ego_cc, showCategory = input$categorySize)
              download_plots$images <- dot1
              svgPanZoom(svglite:::inlineSVG(
                print(dot1),width = 19,
                height = 6,pointsize = 16,
                standalone = TRUE), controlIconsEnabled = T,panEnabled = TRUE,
                dblClickZoomEnabled = TRUE)})
            
          } else if  (input$graphType =="plot5") {
            withProgress(message = "Plotting MF barplot ...",{
              
              bar2 <- barplot(ego_mf, showCategory=input$categorySize)
              download_plots$images <- bar2
              svgPanZoom(svglite:::inlineSVG(print(bar2), width = 19,
                                             height = 6,pointsize = 16,
                                             standalone = TRUE),
                         controlIconsEnabled = T,dblClickZoomEnabled = TRUE)
            })
            
          } else if  (input$graphType =="plot6") {
            withProgress(message = "Plotting MF dotplot ...",{
              dot2 <- dotplot(ego_mf, showCategory =input$categorySize)
              download_plots$images <- dot2
              svgPanZoom(
                svglite::stringSVG(show(dot2),width = 19,
                                   height = 6,pointsize = 16,standalone = F),
                controlIconsEnabled = T, viewBox = FALSE)
              
              
            })
            
          } else if  (input$graphType =="plot9") {
            withProgress(message = "Plotting BP barplot ...",{
              bar3 <- barplot(ego_bp, showCategory=input$categorySize)
              download_plots$images <- bar3
              svgPanZoom(svglite::stringSVG(print(bar3), width = 19,
                                            height = 6,pointsize = 16,standalone = F), controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plot10") {
            withProgress(message = "Plotting BP dotplot ...",{
              r <- scale_colour_viridis_c()
              s <- labs(colour = "pAdjust")
              dot3 <- dotplot(ego_bp,showCategory = input$categorySize)+r+s
              download_plots$images <- dot3
              svgPanZoom(svglite::stringSVG(show(dot3),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
          } else if  (input$graphType =="fig4") {
            withProgress(message = "Plotting ALL categories dotplot ...",{
              dot4 <-  dotplot(ego_all, split="ONTOLOGY", showCategory =input$categorySize) + facet_grid(ONTOLOGY~., scale="free")
              download_plots$images <- dot4
              svgPanZoom(svglite::stringSVG(show(dot4),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
            
          } else if (input$graphType == "fig5") { 
            withProgress(message = "Plotting treeplot ...",{
              outfile7 <- treeplot(go2)
              download_plots$images <- outfile7
              svgPanZoom(svglite::stringSVG(show(outfile7),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
            
          } else if (input$graphType == "heatplt") { 
            withProgress(message = "Plotting Heatplot ...",{
              
              heaTp <- heatplot(ego_all, foldChange=userSigDF$SiGeneListFC, showCategory =input$categorySize)
              download_plots$images <- heaTp
              svgPanZoom(svglite::stringSVG(show(heaTp),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
            
          } else if  (input$graphType =="plot7") {
            withProgress(message = "Plotting MF cnetplot ...",{
              
              
              r <- scale_colour_viridis_c()
              s <- labs(colour = "log2FC")
              
              cnetMf <-  cnetplot(ego_mf, categorySize="pvalue", foldChange=userSigDF$SiGeneListFC, showCategory = input$categorySize, colorEdge = TRUE)+r +s
              download_plots$images <- cnetMf
              svgPanZoom(
                svglite:::inlineSVG(
                  print(cnetMf),width = 19,
                  height = 6,pointsize = 16),
                controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plot3") {
            withProgress(message = "Plotting CC cnetplot ...",{
              r <- scale_colour_viridis_c()
              s <- labs(colour = "log2FC")
              cn2 <- cnetplot(ego_cc, categorySize="pvalue", foldChange=userSigDF$SiGeneListFC, showCategory = input$categorySize, colorEdge = TRUE)+r +s
              download_plots$images <- cn2
              svgPanZoom(svglite::stringSVG(show(cn2),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plot4") {
            withProgress(message = "Plotting ALL cnetplot ...",{
              
              r <- scale_colour_viridis_c()
              s <- labs(colour = "log2FC")
              cn2 <- cnetplot(ego_all, categorySize="pvalue", foldChange=userSigDF$SiGeneListFC, showCategory = input$categorySize, circular = TRUE, colorEdge = TRUE)+r +s
              download_plots$images <- cn2
              svgPanZoom(svglite::stringSVG(show(cn2),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plot11") {
            withProgress(message = "Plotting BP cnetplot ...",{
              r <- scale_colour_viridis_c()
              s <- labs(colour = "log2FC")
              cn3 <- cnetplot(ego_bp, categorySize="pvalue", foldChange=userSigDF$SiGeneListFC, showCategory = input$categorySize, colorEdge = TRUE)+r +s
              download_plots$images <- cn3
              svgPanZoom(svglite::stringSVG(show(cn3),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plotGO1") { 
            withProgress(message = "Plotting CC GOplot ...",{
              outfile1 <- tempfile(fileext='.svg')
              svg(outfile1)
              plotGOgraph(ego_cc,sigForAll = TRUE)
              dev.off()
              download_plots$images <- outfile1
              svgPanZoom((outfile1), controlIconsEnabled = T)
            })
          } else if  (input$graphType =="plotGO2") { 
            withProgress(message = "Plotting MF GOplot ...",{
              outfile2 <- tempfile(fileext='.svg')
              svg(outfile2)
              plotGOgraph(ego_mf,sigForAll = TRUE)
              dev.off()
              download_plots$images <- outfile2
              svgPanZoom(outfile2, controlIconsEnabled = T)})
            
          } else if  (input$graphType =="plotGO3") { 
            withProgress(message = "Plotting BP GOplot ...",{
              outfile3 <- tempfile(fileext='.svg')
              svg(outfile3)
              plotGOgraph(ego_bp, firstSigNodes = 10, useInfo = "all", sigForAll = TRUE,
                          useFullNames = TRUE)
              dev.off() 
              download_plots$images <- outfile3
              svgPanZoom((outfile3), zoomEnabled= TRUE,
                         controlIconsEnabled=TRUE,
                         fit=TRUE,
                         center=TRUE,
                         minZoom= 0.1, viewBox = FALSE)
              
              
            })
            
          } else if  (input$graphType =="pmc1") { 
            withProgress(message = "Plotting publication trend ...",{
              terms <- ego_all$Description[1:5]
              pmc1 <- pmcplot(terms, 2010:2022)
              
              download_plots$images <- pmc1
              svgPanZoom(svglite::stringSVG(show(pmc1),width = 19,
                                            height = 6,pointsize = 16), controlIconsEnabled = T)})
            
            
          }
          
        })
        
        
        output$nemrut <- renderPlotly({
          volcano_data <- userSigDF$dForConv
          
          # add a column of NAs
          volcano_data$diffexpressed <- "NO"
          # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
          volcano_data$diffexpressed[volcano_data$Log2FC > input$logFCut & volcano_data$pAdj < input$pAdjCut] <- "UP"
          # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
          volcano_data$diffexpressed[volcano_data$Log2FC < -input$logFCut & volcano_data$pAdj < input$pAdjCut] <- "DOWN"
          
          
          
          # significant is coded as TRUE, not sig as FALSE
          #volcano_data$sig <- as.factor(abs(volcano_data$LogFC) > input$logFCut & volcano_data$pAdj < input$pAdjCut)
          #convert FDR to -log(FDR)
          volcano_data$negLogpAdj <- -log10(volcano_data$pAdj)
          
          nemrut <- ggplot(volcano_data,aes(label = geneSymbol,  x=Log2FC, y=negLogpAdj, color=diffexpressed)) +
            geom_point() +
            coord_cartesian() +
            scale_color_manual(values=c("blue", "black", "red")) +
            ylab("-log10 FDR") +
            xlab("log2 fold change")+
            geom_vline(xintercept=c(-1, 1), linetype='dashed',col="gray") +
            geom_hline(yintercept= 0, linetype='dashed',col="gray")+
            theme_bw()
          ggplotly(nemrut)
          
          
        })
        
        
        
        # Downloadable csv of selected dataset ----
        output$downloadPlot <- downloadHandler(
          filename =  function() {
            paste(input$graphType, "svg",sep=".")
          },
          # content is a function with argument file. content writes the plot to the device
          content = function(file){
            
            ggsave(file,width=19, height=7, device = "svg")
            
            print(download_plots$images)
            dev.off()
          } 
        )
        
        
        
        
      }
      
      ######### G O _ A N A L Y S I S //SR//  #######################################################  
      
      ########### G O enrichment- Cnetplot ################
      observeEvent(input$selectCluster >0,{
        output$fig3 <-renderVisNetwork({
          
          #######################################################################################
          
          # edge manipulation
          ego_all_df <- as.data.frame(ego_all)
          
          row.names(ego_all_df) <- NULL
          #ego_all_sorted <- ego_all_df %>% arrange(pvalue)
          edges2 <-as.data.frame(ego_all_df[, c("Description","geneID","Description")])
          row.names(edges2) <- NULL
          colnames(edges2) <- c("from","to","link")
          edges2 <- separate_rows(edges2,to, sep = "/")
          edges <- distinct(edges2)
          
          
          # node manipulation
          node2 <- as.data.frame(ego_all_df[, c("geneID")])
          
          node2$group  <-c("Genes")
          node2$size <-c("25")
          node2$cluster <-c("Gene")
          colnames(node2) <- c("id", "group","size","cluster")
          node2 <- separate_rows(node2, id, sep = "/")
          node2 <- distinct(node2)
          node3 <-as.data.frame(ego_all_df[, c("Description", "ONTOLOGY")])
          
          row.names(node3) <- NULL
          
          node3$size <-c("40")
          node3$cluster <-c("0")
          colnames(node3) <- c("id", "group","size","cluster")
          node3 <- distinct(node3)
          
          ### clustering for cnet
          #clusters <-map_dfr(node3$id, ~ {
          # i <- which(stringdist(., node3$id, "jw") < 0.43)
          #tibble(index = i, id = node3$id[i])
          #}, .id = "goClust") %>%
          #distinct(index, .keep_all = T) %>% 
          #mutate(goClust = as.integer(goClust))
          
          #df_merge <- merge(node3, clusters, by = "id",
          #all.x = TRUE)
          
          #df_merge <- subset(df_merge, select=c("id", "group","size","goClust"))
          #colnames(df_merge) <- c("id", "group","size","cluster")
          
          
          regulation <- node3[grep("regulat", node3$id), ]
          positive <- node3[grep("positi", node3$id), ]
          negative <- node3[grep("negati", node3$id), ]
          kinases <- node3[grep("kinase", node3$id), ]
          proliferation <- node3[grep("proliferation", node3$id), ]
          apoptotic <- node3[grep("apopto", node3$id), ]
          node3$cluster <-c("Other")
          node3 <- rows_update(node3, tibble(id = regulation$id , cluster = "Regulation"))
          node3 <- rows_update(node3, tibble(id = apoptotic$id , cluster = "Apoptotic Process"))
          node3 <- rows_update(node3, tibble(id = proliferation$id , cluster = "Proliferation"))
          node3 <- rows_update(node3, tibble(id = positive$id , cluster = "Positive regulation"))
          node3 <- rows_update(node3, tibble(id = negative$id , cluster = "Negative regulation"))
          node3 <- rows_update(node3, tibble(id = kinases$id , cluster = "Kinase"))
          ###########################
          nodes_all <- rbind(node2,node3)
          nodes_all <- distinct(nodes_all)
          #nodes_all <- head(nodes_all,300)
          
          
          nodes_all$label <- nodes_all$id
          nodes_all$title <- nodes_all$id
          nodes <-distinct(nodes_all)
          
          nodes_selected_cluster <-nodes%>%filter(grepl(input$selectCluster,cluster))
          selected_edges <-edges %>%
            filter(from %in% nodes_selected_cluster$id)
          modified_nodes <- nodes%>%
            filter(id %in% selected_edges$to)
          
          
          nodes <-rbind(nodes_selected_cluster,modified_nodes)
          edges <- as.data.frame(selected_edges)
          
          
          #######################################################################################
          visNetwork(nodes, edges,width = "100%", height = "500px") %>%
            visLayout(randomSeed = 123)%>%
            visGroups(groupname = "MF",color = list(background = " #bf6ca9 ", border = " #bf6ca9 ", 
                                                    highlight = " #ca97bd"),shape = "dot") %>%
            visGroups(groupname = "BP", color = list(background = "#8dcbb0", border = "#8dcbb0", 
                                                     highlight = "#b7d1c6"),shape = "dot") %>%
            visGroups(groupname = "CC",  color = list(background = "#d6de77", border = "#d6de77", 
                                                      highlight = "#d2d5aa"),shape = "dot") %>%
            visGroups(groupname = "Genes",  color = list(background = "#77c2de", border = "#77c2de", 
                                                         highlight = "#aed1de"), shape = "triangle") %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1), 
                       nodesIdSelection = TRUE)%>%
            #visLayout(randomSeed = 123)%>%
            #visPhysics(stabilization = TRUE)%>%
            visPhysics(solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -10))%>%
            visLegend(main = list(text = "Legend",
                                  style = "font-family:Comic Sans MS;color: #212f3c;font-size:12px;text-align:center;"),
                      width = 0.05, position = "right")
        })
        
      })
      
      
      
      #### K E G G _ A N A L Y S I S //SR// #################################################################################
      
      
      if(input$GO_species=="org.Ce.eg.db" ){
        
        kkt <- enrichPathway(gene = userSigDF$SiGenesList, organism = "celegans", pvalueCutoff = input$pvalCut, readable = TRUE)
        kk_df <- data.frame(kkt)
        
        
      }
      
      else{
        
        id_cnvrt<-bitr_kegg(userSigDF$SiGenesList, fromType='ncbi-geneid', toType='kegg', organism = input$KEGG_species)
        userSigDF$SiGenesList <- id_cnvrt$kegg
       
        ekk <- enrichKEGG(gene         = userSigDF$SiGenesList,
                          organism     = input$KEGG_species,
                          minGSSize = 1,
                          maxGSSize = 500,
                          pAdjustMethod = input$pAdMeth,
                          pvalueCutoff  = input$pvalCut,
                          qvalueCutoff  = input$qvalCut
                          
        )
        
        kk_df <- data.frame(ekk)
      }
      
      row.names(kk_df) <- NULL
      
      
      output$kk_df <- DT::renderDataTable(kk_df,
                                          extensions = c('Buttons','Responsive','Scroller'), 
                                          options = list(
                                            dom = 'Bfrtip',
                                            lengthChange = FALSE,
                                            deferRender = TRUE,
                                            #ordering = TRUE,
                                            scrollY = 400,
                                            scroller = TRUE,
                                            buttons = 
                                              list( list(
                                                extend = 'collection',
                                                
                                                buttons = list(list(extend='csv',
                                                                    filename = 'PRO-simat_KEGG_table'),
                                                               list(extend='excel',
                                                                    filename = 'PRO-simat_KEGG_table'),
                                                               list(extend='pdf',
                                                                    filename= 'PRO-simat_KEGG_table')),
                                                text = '<span class="glyphicon glyphicon-download-alt"></span>'
                                              )),
                                            
                                            initComplete = JS(
                                              "function(settings, json) {",
                                              "$(this.api().table().header()).css({'background-color': '#8E44AD', 'color': '#fff'});",
                                              "}")))
      
      
      
      gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
      c5 <- read.gmt(gmtfile)
      
      egmt <- enricher(userSigDF$SiGenesList, TERM2GENE=c5)
      
      if(nrow(kk_df) < 1) {
        ask_confirmation(
          "form3",
          title = " ",
          type = "info",
          btn_labels = c("Cancel", "OK"),
          allowEscapeKey = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          text =  
            div(HTML("Oops!", "KEGG Pathway analysis could not be performed for the species you selected."
            ))
        ) 
      }
      
      
      
      
      #______________________________________________________Search Box for KEGG______________________________________________ 
      observeEvent(input$searchButton > 0,{
        if(input$searchButton > 0){
          isolate({
            
            browseKEGG(ekk, input$searchText)
          })
        }
      })
      
      output$plot13 <-renderSvgPanZoom({
        x2 <- pairwise_termsim(ekk)
        ema1 <- emapplot(x2,showCategory =input$categorySizeEma, layout="kk")
        download_plots$images <- ema1
        svgPanZoom(svglite::stringSVG(show(ema1),width = 19,
                                      height = 6,pointsize = 16), controlIconsEnabled = T)
        
      })
      
      ############### emaplot download ###########################
      
      output$downloadEma <- downloadHandler(
        filename =  function() {
          paste("plot13", "png",sep=".")
        },
        content = function(file){
          
          png(file, units="in", width=19, height=7, res=300)
          
          
          dev.off()
        } 
      )
      
      
      prosim <- reactiveValues()
      pathviewReactive <- eventReactive(input$PathviEw,{
        withProgress(message = 'Plotting Pathview Pathway...', {
          
          isolate({
            setProgress(value = 0.7, detail = paste0("Pathview ID ",input$PathviEw," ..."))
            dme <- pathview(gene.data  =userSigDF$SiGeneListFC,
                            pathway.id = input$PathviEw,
                            species    = input$KEGG_species)
            file.copy(paste0(input$PathviEw,".pathview.png"),paste0(input$PathviEw))
            
            #setProgress(value = 0.7, detail = paste0("Pathview ID ",input$PathviEw," generating pdf ..."))
            dmePdf <- pathview(gene.data  =userSigDF$SiGeneListFC,
                               pathway.id = input$PathviEw,
                               species    = input$KEGG_species, kegg.native = F)
            
            prosim$imagePath = paste0(input$PathviEw,".pathview.")
            return(list(
              src = paste0(input$PathviEw),
              filetype = tempfile("image/png"),
              alt = "pathview image"
            ))
          })
          
        })
      })
      
      output$pathway1  = renderImage({
        
        return(pathviewReactive())
        
      })
      
      
      
      output$downloadPathviewPng <- downloadHandler(
        filename = function()  {paste0(prosim$imagePath,"png")},
        content = function(file) {
          file.copy(paste0(getwd(),'/',prosim$imagePath,"png"), file)
        }
      )
      
      output$downloadPathviewPdf <- downloadHandler(
        filename = function()  {paste0(prosim$imagePath,"pdf")},
        content = function(file) {
          file.copy(paste0(getwd(),'/',prosim$imagePath,"pdf"), file)
        }
      )
      
    } 
  })#analysis submit button
  # }) 
  # })# multiprotein input part
  
  
  #------------------------------------------------------------------------- K E G G _ A N A L Y S I S
  
  
  
  
  #____________________________________________________________________________________________________________________ 
  
  ## Jimena example txt file 
  observeEvent(input$choseJimena,{
    if("userJim" == input$choseJimena){
      #shinyjs::toggle(id = "advanced", anim = TRUE)
      shinyjs::show(id = "bxjim", anim = TRUE)
    }
    
    if("exampleJim" == input$choseJimena){
      #shinyjs::toggle(id = "advanced", anim = TRUE)
      shinyjs::show(id = "bxjim", anim = TRUE)
    }
  })
  
  
  
  #############################
  #----------------------------------------------------------- Upload csv and convert graphml file
  
  
  #### uploaded txt file in box
  
  output$whichjim <-  DT::renderDataTable(server=FALSE,{
    if (input$choseJimena == "userJim" && input$file1>0) {
      
      #observeEvent(input$file1,{
      withProgress(message = "Importing your file...",{
        isolate({
          inFile <- input$file1
          if (is.null(inFile)) return(NULL)   
          dff <- read.table(inFile$datapath,header = TRUE)
          
          write.table(dff, "input.txt",row.names = FALSE,quote=FALSE)
          
        }) 
      })
      
      #}) 
      
      user_txt <- read.table("input.txt",header = TRUE)
      DT::datatable(user_txt,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 200,
                      scroller = TRUE))
    }
    else if(input$choseJimena == "exampleJim") {
      test_file <- read.table("testJimena.txt",header = TRUE)
      
      write.table(test_file, "input.txt",row.names = FALSE,quote=FALSE)
      DT::datatable(test_file,
                    extensions = c('Responsive','Scroller'), 
                    options = list(
                      lengthChange = FALSE,
                      deferRender = TRUE,
                      ordering = TRUE,
                      scrollY = 200,
                      scroller = TRUE
                    ))}
  })
  
  
  #### run button for jimena 
  observeEvent(input$runGraphml,{
    withProgress(message = "Converting graphml...",{
      system(paste("perl txt_to_graphml.pl"))
    })
  })
  
  
  ################### JIMENA //SR// #####################################################################################
  
  
  
  
  daf <- reactiveVal()
  add_per <- reactiveValues(data=NULL)
  result_val <- reactiveVal()
  add_per <- reactiveValues(NodeInfo=NULL)
  add_per <- reactiveValues(NodeIndex=NULL)
  add_per <- reactiveValues(remained_nodes=NULL)
  remain_per<- reactiveValues(NodeIndex=NULL)
  remain_per<- reactiveValues(data_name=NULL)
  add_parameter <- reactiveValues()
  remain_parameter <- reactiveValues()
  squad_nodes <- reactiveValues(data=NULL)
  
  jimenaout <- reactiveValues(data = NULL)
  observeEvent(input$runJimena,{
    withProgress(message = 'Please wait for Jimena calculation...', {
      req(input$runGraphml)
      GRAPHML_PATH <- "input.graphml"
      
      
      system("rm jimena_time_series_data.csv")
      system(paste("java -classpath jimena-app/jimena.jar:jimena-app/  App ", GRAPHML_PATH, sep=" "))
      jimena_output <- read.csv("jimena_time_series_data.csv")
      jimena_output_long <- melt(jimena_output, id.vars = "time")   
     
      jimenaout$data <-jimena_output
      
      accumulate_by <- function(dat, var) {
        var <- lazyeval::f_eval(var, dat)
        lvls <- plotly:::getLevels(var)
        dats <- lapply(seq_along(lvls), function(x) {
          cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
        })
        dplyr::bind_rows(dats)
      }
      
      fig <- jimena_output_long %>% accumulate_by(~time)
     
      
      jimena_time_series <- 
        ##################### works for without animation ######################################
      fig%>%plot_ly()%>%
        add_lines(x =~time,  y = ~value, split = ~variable)
      jimena_time_series <- jimena_time_series %>% layout(
        xaxis = list(
          title = "Time",
          zeroline = F
        ),
        yaxis = list(
          title = "Activation",
          zeroline = F
        )
      )
      #######################################
      
      ###########################################################
      
      #fig %>%
      #plot_ly(
      # x = ~time, 
      #y = ~value,
      #split = ~variable,
      #frame = ~frame, 
      #type = 'scatter',
      #mode = 'lines', 
      #line = list(simplyfy = F)
      #)
      #jimena_time_series <- jimena_time_series %>% layout(
      # xaxis = list(
      #  title = "Time",
      # zeroline = F
      #),
      #yaxis = list(
      #  title = "Activation",
      #  zeroline = F
      #)
      #)
      
      #jimena_time_series <- jimena_time_series %>% animation_opts(
      # frame = 100, 
      #transition = 0, 
      #redraw = TRUE
      #)
      #jimena_time_series <- jimena_time_series %>% animation_slider(
      # hide = T
      #)
      #jimena_time_series <- jimena_time_series %>% animation_button(
      # x = 1, xanchor = "right", y = 0, yanchor = "bottom"
      #)
      
      result_val(jimena_time_series)
      
      add_per$data <- data.frame(jimena_output)
      Snode <- as.data.frame(colnames(add_per$data))
      squad_nodes$data <- Snode[-1,]
      
      
      write.table(add_per$data, file = paste0("for_SQUADD/Jimena_output/Index_", remain_per$data_name,".txt"), sep="\t", row.names=FALSE,col.names=FALSE)
      
      outVar <- reactive({
        vars <- as.data.frame(colnames(add_per$data))
        vars <- vars[-1,]
        print(vars)
        
        return(vars)
        
      })
      
      output$inVar2 <- renderUI({
        selectInput(inputId = "inVar2", label = h4("Add Perturbation"), choices =  outVar())
        
        
      })
      
    })
    
    
  }) 
  
  output$jimenaResult <-renderPlotly(
    result_val())
  
  
  ##### JIMENA ADD PERTURBATION PART ####################################################################################
  
  
  observeEvent(input$add_pert, {
    t = rbind(data.frame(Nodes = input$inVar2,
                         Start = input$strt,
                         End = input$end,
                         Value = input$val), daf())
    t <- distinct(t)
    daf(t)
    all_node_ <- as.data.frame(colnames(add_per$data))
    all_node_ <- all_node_[-1,]
    all_node_ <- as.data.frame(all_node_)
    node_length <- nrow(all_node_)
    all_node_$index <- rep(0:(node_length-1))
    colnames(all_node_) <- c("Nodes", "NodeIndex")
    
    
    
    add_per$NodeInfo <-all_node_
    r <- all_node_[which(all_node_$Nodes==input$inVar2),]
    remain_parameter<- r
    remain_per$data_name <- remain_parameter$Nodes
   
    add_per$NodeIndex <- r$NodeIndex
    add_parameter <- rbind(data.frame(    Index = 0,
                                          Nodes = add_per$NodeIndex,
                                          Start = input$strt,
                                          End = input$end,
                                          Value = input$val))
    
    add_parameter <- distinct(add_parameter)
    csv_fname = "jimena-app/ParameterInputs.csv"
    write.table(add_parameter, file = csv_fname, sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
  })
  
  
  
  observeEvent(input$remove_per, {
    t = daf()
   
    
    if (!is.null(input$shiny_table_rows_selected)) {
      t <- t[-as.numeric(input$shiny_table_rows_selected),]
      add_per$remained_nodes <- t
      
    }
    daf(t)
    
    remain <- as.data.frame(add_per$remained_nodes)
    
    node_info <- add_per$NodeInfo
    
    
    last_remain_info <- subset(node_info, (Nodes %in% remain$Nodes))
    
    if (nrow(last_remain_info)>0){
      
     
      merged_part <- merge(last_remain_info, remain, by= "Nodes")
      
      
      remain_parameter <- rbind(data.frame(    Index = 0,
                                               Nodes = merged_part$NodeIndex,
                                               Start = merged_part$Start,
                                               End = merged_part$End,
                                               Value = merged_part$Value,
                                               Fix = NA))
      
      remain_parameter <- distinct(remain_parameter)
      colnames(remain_parameter) <- c(0,0.05,10,1,0.1,0.5)
      
     
      csv_fname = "jimena-app/ParameterInputs.csv"
      write.table(remain_parameter, file = csv_fname, sep = ",",
                  append = FALSE, quote = FALSE,
                  col.names = TRUE, row.names = FALSE)
    }
    
    if (nrow(last_remain_info)==0){
      
      remain_parameter <- rbind(data.frame(    A = 0,
                                               B = 0.05,
                                               C = 10,
                                               D = 1,
                                               E = 0.1,
                                               G = 0.5))
      
      
      csv_fname = "jimena-app/ParameterInputs.csv"
      write.table(remain_parameter, file = csv_fname, sep = ",",
                  append = FALSE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
    
  })
  
  output$shiny_table <- renderDataTable(
    daf(), selection = 'single', extensions = c('Buttons','Responsive','Scroller'), 
    options = list(
      dom = 'Bfrtip',
      lengthChange = FALSE,
      deferRender = TRUE,
      ordering = TRUE,
      scrollY = 100,
      scroller = TRUE,
      buttons = 
        list( list(
          extend = 'collection',
          
          buttons = list(list(extend='csv',
                              filename = 'PRO-simat_perturbation_table'),
                         list(extend='excel',
                              filename = 'PRO-simat_perturbation_table'),
                         list(extend='pdf',
                              filename= 'PRO-simat_perturbation_table')),
          text = '<span class="glyphicon glyphicon-download-alt"></span>'
        )),
      
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#0080F1', 'color': '#fff'});",
        "}")))
  
  
  
  ############################ Jimena output Data Table ####################
  
  output$jimenaoutput_tb <- DT::renderDataTable(server = FALSE,{
    DT::datatable(jimenaout$data,
                  extensions = c('Buttons','Responsive','Scroller'), 
                  options = list(
                    dom = 'Bfrtip',
                    lengthChange = FALSE,
                    deferRender = TRUE,
                    ordering = TRUE,
                    scrollY = 400,
                    scroller = TRUE,
                    buttons = 
                      list( list(
                        extend = 'collection',
                        
                        buttons = list(list(extend='csv',
                                            filename = 'PRO-simat_jimenaout_table'),
                                       list(extend='excel',
                                            filename = 'PRO-simat_jimenaout_table'),
                                       list(extend='pdf',
                                            filename= 'PRO-simat_jimenaout_table')),
                        text = '<span class="glyphicon glyphicon-download-alt"></span>'
                      )),
                    
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#0080F1', 'color': '#fff'});",
                      "}")))})
  
  
  
  
  
  #____________________________________________________________________________________________________________________ 
  ####### SQUADD PACKAGE //SR// ####
  
  #if(nrow(squad_nodes) < 1) {
  # shinyalert("Oops!", "First you have to run Jimena and add nodes using add perturbation option. ", type = "info",
  #           time = 4000, showConfirmButton= FALSE)
  
  #}
  
  squad_legend <- reactive({
    squadd_node_list <- as.data.frame(colnames(add_per$data))
    squadd_node_list <- squadd_node_list[-1,]
   
    return(squadd_node_list)
  })
  output$squad_node <- renderUI({
    selectInput(inputId = "squad_node", label = h4("Select Nodes"), choices = squad_legend(), multiple = TRUE)
    
  })
  
  
  
  observeEvent(input$squaddButton, {
    
    clean_parameter <- rbind(data.frame(    A = 0,
                                            B = 0.05,
                                            C = 10,
                                            D = 1,
                                            E = 0.1,
                                            G = 0.5))
    
   
    csv_fname = "jimena-app/ParameterInputs.csv"
    write.table(clean_parameter, file = csv_fname, sep = ",",
                append = FALSE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
    req(input$add_pert)
    
    if(length(input$squad_node) == 0) {
      ask_confirmation(
        "form3",
        title = " ",
        type = "info",
        btn_labels = c("Cancel", "OK"),
        allowEscapeKey = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        text =  
          div(HTML("Oops!", "Please select the nodes you want to use! "
          ))
      ) 
      
    } else {
      
      
      #folder <- file.path(".","for_SQUADD","Jimena_output")
      
      sim <- simResService (folder= "for_SQUADD/Jimena_output",
                            time= 10,
                            ncolor=6,
                            legend= squad_legend(),
                            indexDeno=1,
                            method="lowess")
      
      list_squadd2 <- input$squad_node
      list_squadd2 <- list_squadd2[-1]
      
      
      sim["selectNode"] <- input$squad_node
      
      
      output$squadd1 <- renderPlot({plotSimMatrix(sim)
      })
      
      output$text1 <- renderText({ 
        xaxis <-  input$squad_node
        #xaxis<- data.frame(xaxis)
        
        #HTML(paste("Labels from left to right for each chart:"),xaxis)
        #paste("Labels from top to bottom:", xaxis)
      })
      
      output$text2 <- renderTable({
        
        yaxis <- list.files(path = "for_SQUADD/Jimena_output/")
        
        a <- str_extract(yaxis, "Index*(\\w+)")
        a <- data.frame(a)
        colnames(a) <- "Data used from top to bottom"
        a
        
        
        #cat(unlist(a), sep="\n")
        
        #HTML(paste("top down label list:"),a)
        
        #paste("Labels from left to right for each chart:", a)
      })
      
      
      tab <- getFittedTable(sim)
      output$squadd2 <- renderPlot({plotPredMap(sim)
      })
      
      
      sim["conditionList"] <- list_squadd2
      output$squadd3 <- renderPlot({plotCC(sim)
        
      })
    }
  })
  
  ###############################################################################################################################
  ## SQUADD PACKAGE ##
  ###############################################################################################################################
  
  
  
  observeEvent(input$savePNGbutton, ignoreInit=TRUE, {
    file.name <- tempfile(fileext=".png")
    savePNGtoFile(session, file.name)
  })
  
  
  observeEvent(input$pngData, ignoreInit=TRUE, {
   
    png.parsed <- fromJSON(input$pngData)
    substr(png.parsed, 1, 30) 
    nchar(png.parsed)  
    png.parsed.headless <- substr(png.parsed, 23, nchar(png.parsed)) 
    png.parsed.binary <- base64decode(png.parsed.headless)
   
    conn <- file("network.png", "wb")
    writeBin(png.parsed.binary, conn)
    close(conn)
    
  })
  #_____________________________________________________________________________________________________________________  
  
  text_reactive <- eventReactive( input$submitButton, {
    
  })
  
  # text output
  output$text <- renderText({
    text_reactive()
  })
  
  #_________________________________________________________________________________________________________________   
  
  ########################### ValueBox of Homepage ########################################################
  
  
  onclick('clickdiv', updateTabsetPanel(session, inputId = "tab2", selected = "Category_1")
  )
  
  onclick('clickdiv2', updateTabsetPanel(session, inputId = "tab2", selected = "Category_2")
  )
  
  onclick('clickdiv3', updateTabsetPanel(session, inputId = "tab2", selected = "Category_3")
  )
  
  onclick('clickdiv4', updateTabsetPanel(session, inputId = "tab2", selected = "Category_4")
  )
  
  onclick('clickdiv5', updateTabsetPanel(session, inputId = "tab2", selected = "Category_5")
  )
  
  ######### Kill DB Connections________________________________________________
  #track_usage(storage_mode = store_json(path = "logs/"))
  #con <- file("shiny.log")
  #sink(con, append=TRUE)
  #sink(con, append=TRUE, type="message")
  #message(paste0("Logged in at ", Sys.time()))
  
  #___________________________________________________________________
  FILES <- list.files(pattern = "hsa*")
  FILES2 <- list.files(pattern = "for_SQUADD/Jimena_output/*.txt")
  
  session$onSessionEnded(function() {
    if (!is.null(FILES) || !is.null(FILES2)) {
      system("rm hsa*")
      system("rm for_SQUADD/Jimena_output/*.txt")
      system("rm input.txt")
      system("rm input.graphml")
      
    }
    #system("rm for_SQUADD/Jimena_output/*.txt")
    #system("rm hsa*")
    system("rm jimena_time_series_data.csv")
    source("killConnection.R", local = TRUE)
    killDbConnections()
    clean_parameter <- rbind(data.frame(    A = 0,
                                            B = 0.05,
                                            C = 10,
                                            D = 1,
                                            E = 0.1,
                                            G = 0.5))
    
    
    csv_fname = "jimena-app/ParameterInputs.csv"
    write.table(clean_parameter, file = csv_fname, sep = ",",
                append = FALSE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    #stopApp()
    
  })
  #system("rm shiny.log")
  
}
)
# server
#----------------------------------------------------------------------------------------------------
