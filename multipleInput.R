e<-reactiveValues(data=NULL)


observe({
  multipleInput()
})




multipleInput<- eventReactive(input$multiProtein,{
  x = reactiveVal(10)
  dbGetQuery(con(), "USE IGOSynthetic;")
  
  dbGetQuery(con(), "SET SESSION max_statement_time = 900;")
  
  if(input$multiProtein > 0){
    
    isolate({
      
      print(input$multiProtein)
      multi = data_frame(input$multiProtein)
      colnames(multi) <- c("X")
      multi <- separate_rows(multi,X, sep = "\n")
      
      multi <- distinct(multi)
      
      
      edges <- data_frame()
      withProgress(message = 'Please wait for protein interaction...', {
        
        qmain3 <- reactiveValues(data=NULL)
        searchList <- data_frame(NULL)
        for (i in multi$X){
          
          query3<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",i,"'or GeneInfoIDto='",i,"')limit ",x(),";")
          qmain3 <-dbGetQuery(con(), query3)
          print(qmain3)
          
          searchList <- rbind(searchList, qmain3)
          searchall <-  data.frame(GeneInfoIDto = c(searchList[,"GeneInfoIDto"], searchList[,"GeneInfoIDfrom"]))
          searchall <- distinct(searchall)
         
          
        }
        #____________________________________________________________________________
        
        observeEvent(input$MoReBtn,{
          withProgress(message = 'More protein interaction...', {
            
            x(x()+4)
            searchList <- data_frame(NULL)
            
            for (i in multi$X){
              query3<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",i,"' or GeneInfoIDto='",i,"')limit ",x(),";")
              qmain3 <-dbGetQuery(con(), query3)
              searchList <- rbind(searchList, qmain3)
              searchall <-  data.frame(GeneInfoIDto = c(searchList[,"GeneInfoIDto"], searchList[,"GeneInfoIDfrom"]))
              searchall <- distinct(searchall)
              
            }
            
            inter_ <- data_frame()
            for (j in searchall$GeneInfoIDto){
              for (k in searchall$GeneInfoIDto){
                if(j==k){next
                }
                query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')")
                qmain4 <-dbGetQuery(con(), query4)
                
                inter_ <- rbind(inter_, qmain4)
              }
            }
            
            if(nrow(inter_) > 0){
              for (r in 1:nrow(inter_)){
                dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
                query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        limit ",x(),";")
                #and GeneInfo.Taxonomy_ID IN('",input$hostID,"')
                
                dbGetQuery(con(), query5)
                
                query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                         ;")
                #WHERE 
                #GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
                
                qmain2 <- dbGetQuery(con(),query2)
                
                edges <- rbind(edges, qmain2)
                edges <- distinct(edges)
              }
              
              e$data <- as.data.frame(edges)
            }
            
            else if (nrow(inter_) < 1){
              
              dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
              query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",i,"'
                        or
                        Gene2Gene.GeneInfoIDto = '",i,"'
                        limit ",x(),"
                        ;")
              dbGetQuery(con(), query5)
              
              query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                         ;")
              
              qmain2 <- dbGetQuery(con(),query2)
              
              edges <- rbind(edges, qmain2)
              edges <- distinct(edges)
            }
            
            e$data <- as.data.frame(edges)
            
          })
        }) #more
        
        #____________________________________________________________________________
        #____________________________________________________________________________
        
        observeEvent(input$LessBtn,{
          withProgress(message = 'Less protein interaction...', {
            x(abs(x()-4))
            
            searchList <- data_frame(NULL)
            
            for (i in multi$X){
              query3<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",i,"' or GeneInfoIDto='",i,"')limit ",x(),";")
              qmain3 <-dbGetQuery(con(), query3)
              searchList <- rbind(searchList, qmain3)
              searchall <-  data.frame(GeneInfoIDto = c(searchList[,"GeneInfoIDto"], searchList[,"GeneInfoIDfrom"]))
              searchall <- distinct(searchall)
            }
            
            inter_ <- data_frame()
            for (j in searchall$GeneInfoIDto){
              for (k in searchall$GeneInfoIDto){
                if(j==k){next
                }
                query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')")
                qmain4 <-dbGetQuery(con(), query4)
                
                inter_ <- rbind(inter_, qmain4)
              }
            }
            
            if(nrow(inter_) > 0){
              for (r in 1:nrow(inter_)){
                dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
                query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        
                        limit ",x(),";")
                
                #and GeneInfo.Taxonomy_ID IN('",input$hostID,"')
                dbGetQuery(con(), query5)
                
                query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                         ;")
                #WHERE 
                #GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
                
                qmain2 <- dbGetQuery(con(),query2)
                
                edges <- rbind(edges, qmain2)
                edges <- distinct(edges)
              }
              
              e$data <- as.data.frame(edges)
            }
            else if (nrow(inter_) < 1){
              
              dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
              query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",i,"'
                        or
                        Gene2Gene.GeneInfoIDto = '",i,"'
                        limit ",x(),"
                        ;")
              dbGetQuery(con(), query5)
              
              query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                         ;")
              
              qmain2 <- dbGetQuery(con(),query2)
              
              edges <- rbind(edges, qmain2)
              edges <- distinct(edges)
            }
            
            e$data <- as.data.frame(edges)
            
          })
        }) #less
        
        #____________________________________________________________________________
        #____________________________________________________________________________
        
        
        inter_ <- data_frame()
        for (j in searchall$GeneInfoIDto){
          for (k in searchall$GeneInfoIDto){
            if(j==k){next
            }
            query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')
                         ;")
            qmain4 <-dbGetQuery(con(), query4)
            
            inter_ <- rbind(inter_, qmain4)
          }
        }
        if(nrow(inter_) > 0){
          
          for (r in 1:nrow(inter_)){
            dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
            
            query1<- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        
                        limit ",x(),";")
            
            #and GeneInfo.Taxonomy_ID IN('",input$hostID,"')
            dbGetQuery(con(), query1)
            
            query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        
                         ;")
            #WHERE 
            #GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
            
            qmain2 <- dbGetQuery(con(),query2)
            
            edges <- rbind(edges, qmain2)
            edges <- distinct(edges)
          }
          
          e$data <- as.data.frame(edges)
        } else if (nrow(inter_) < 1){
          
          
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",i,"'
                        or
                        Gene2Gene.GeneInfoIDto = '",i,"'
                        limit ",x(),"
                        ;")
          dbGetQuery(con(), query5)
          
          query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                         ;")
          
          qmain2 <- dbGetQuery(con(),query2)
          
          edges <- rbind(edges, qmain2)
          edges <- distinct(edges)
        }
        
        e$data <- as.data.frame(edges)
        
        
      })
      
    }) #isolate
    
  }  #if multiple protein
  
  ################################################################################################################  
  else if (input$pathogenID == "0"){
    isolate({
      edges <- read.csv("GLV-1h68.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
  ################################################################################################################  
  else if (input$pathogenID == "10"){
    isolate({
      edges <- read.csv("VVTKN1L.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
  
  else if (!is.null(input$multiProtein) && input$hostID == "9606" && input$pathogenID == "9606"){
    
    edges <- data_frame()
    withProgress(message = 'Please wait for protein interaction...', {
      
      
      query3<-paste0("
            SELECT DISTINCT
            Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
            Gene2Gene.GeneInfoIDto AS GeneInfoIDto
            FROM 
            Gene2Gene 
            INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
            INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
            WHERE 
            GeneInfo.Taxonomy_ID IN(",input$hostID,")
            limit ",x(),";")
      
      qmain3 <-dbGetQuery(con(), query3)
      
      inter_ <- data_frame()
      for (j in qmain3$GeneInfoIDto){
        for (k in qmain3$GeneInfoIDto){
          if(j==k){next
          }
          query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')")
          qmain4 <-dbGetQuery(con(), query4)
          
          inter_ <- rbind(inter_, qmain4)
        }
      }
      
      for (r in 1:nrow(inter_)){
        dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
        query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        and GeneInfo.Taxonomy_ID IN(",input$hostID,")
                        limit ",x(),";")
        dbGetQuery(con(), query5)
        
        query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
                         ;")
        
        qmain2 <- dbGetQuery(con(),query2)
        
        edges <- rbind(edges, qmain2)
        edges <- distinct(edges)
        
      }
      
      e$data <- as.data.frame(edges)
      
      
      
      
    })
    
    
    observeEvent(input$MoReBtn,{
      withProgress(message = 'More protein interaction...', {
        x(x()+4)
        
        searchList <- data_frame(NULL)   
        
        query3<-paste0("
            SELECT DISTINCT
            Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
            Gene2Gene.GeneInfoIDto AS GeneInfoIDto
            FROM 
            Gene2Gene 
            INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
            INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
            WHERE 
            GeneInfo.Taxonomy_ID IN (",input$hostID,")
            limit ",x(),";")
        
        qmain3 <-dbGetQuery(con(), query3)
        searchList <- rbind(searchList, qmain3)
        
        
        inter_ <- data_frame()
        for (j in searchList$GeneInfoIDto){
          for (k in searchList$GeneInfoIDto){
            if(j==k){next
            }
            query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')")
            qmain4 <-dbGetQuery(con(), query4)
            
            inter_ <- rbind(inter_, qmain4)
          }
        }
        
        
        
        for (r in 1:nrow(inter_)){
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        and GeneInfo.Taxonomy_ID IN(",input$hostID,")
                        limit ",x(),";")
          dbGetQuery(con(), query5)
          
          query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
                         ;")
          
          qmain2 <- dbGetQuery(con(),query2)
          
          edges <- rbind(edges, qmain2)
          edges <- distinct(edges)
        }
        
        e$data <- as.data.frame(edges)
        
        
      })
    }) #more
    
    #___________________________________
    
    
    observeEvent(input$LessBtn,{
      withProgress(message = 'Less protein interaction...', {
        x(abs(x()-4))
        
        
        query3<-paste0("
            SELECT DISTINCT
            Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
            Gene2Gene.GeneInfoIDto AS GeneInfoIDto
            FROM 
            Gene2Gene 
            INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
            INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
            WHERE 
            GeneInfo.Taxonomy_ID IN(",input$hostID,")
            limit ",x(),";")
        
        qmain3 <-dbGetQuery(con(), query3)
        
        
        inter_ <- data_frame()
        for (j in qmain3$GeneInfoIDto){
          for (k in qmain3$GeneInfoIDto){
            if(j==k){next
            }
            query4<-paste0("SELECT * FROM Gene2Gene WHERE (GeneInfoIDfrom='",j,"') AND (GeneInfoIDto='",k,"')")
            qmain4 <-dbGetQuery(con(), query4)
            
            inter_ <- rbind(inter_, qmain4)
          }
        }
        
        
        
        for (r in 1:nrow(inter_)){
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query5 <- paste0("CREATE VIEW temptable AS 
                        SELECT DISTINCT
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        Gene2Gene.GeneInfoIDfrom = '",inter_[r,1],"'
                        and
                        Gene2Gene.GeneInfoIDto = '",inter_[r,2],"'
                        and GeneInfo.Taxonomy_ID IN(",input$hostID,")
                        limit ",x(),";")
          dbGetQuery(con(), query5)
          
          query2<-paste0("SELECT DISTINCT
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN(",input$pathogenID,")
                         ;")
          
          qmain2 <- dbGetQuery(con(),query2)
          
          edges <- rbind(edges, qmain2)
          edges <- distinct(edges)
        }
        
        e$data <- as.data.frame(edges)
        
      })
    }) #less
    
    #___________________________________
    
    
    
  }
  
  ####################################################################################################
  
  else if(!is.null(input$multiProtein) && input$hostID != "9606" && input$pathogenID == "9606"){
    
    x <- reactiveVal(100)
    withProgress(message = 'Please wait for protein interaction...', {
      dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
      query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
      
      
      dbGetQuery(con(), query1)
      
      
      dbGetQuery(con(), "DROP TABLE IF exists set2;")
      query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ;")
      
      dbGetQuery(con(), query2)                
      dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
      query3<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
      
      dbGetQuery(con(), query3)                
      dbGetQuery(con(), "DROP TABLE IF exists set3;")
      query4<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ;")
      
      dbGetQuery(con(), query4)                
      query5<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set3
                        LIMIT ",x(),";")
      
      
      qmain1 <- dbGetQuery(con(),query5)
      #edges <- distinct(qmain2)
      
      e$data <- as.data.frame(qmain1) 
      
    })
    
    
    observeEvent(input$MoReBtn,{
      withProgress(message = 'More protein interaction...', {
        x(x()+150)
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
          
          dbGetQuery(con(), query1)
          
          dbGetQuery(con(), "DROP TABLE IF exists set2;")
          query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$pathogenID,") DESC
                        ;")
          
          
          dbGetQuery(con(), query2)
          
          dbGetQuery(con(), "DROP TABLE IF exists set3;")
          query3<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        LIMIT ",x(),"
                        ;")
          dbGetQuery(con(), query3) 
          
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
          query4<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
          
          
          dbGetQuery(con(), query4)                
          dbGetQuery(con(), "DROP TABLE IF exists set4;")
          query5<-paste0("create temporary table set4
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$pathogenID,") DESC
                        ;")
          
          dbGetQuery(con(), query5)                
          query6<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set4
                         UNION
                          SELECT * FROM set3
                         LIMIT ",x(),";")
          
          
          qmain1 <- dbGetQuery(con(),query6)
          #edges <- distinct(qmain2)
          
          
          edges <- distinct(qmain1)
          e$data <- as.data.frame(qmain1)
          
          
        })
      })
    })
    
    
    observeEvent(input$LessBtn,{
      withProgress(message = 'Less protein interaction...', {
        x(abs(x()-140))
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
          
          dbGetQuery(con(), query1)
          
          dbGetQuery(con(), "DROP TABLE IF exists set2;")
          query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$pathogenID,") DESC
                        ;")
          
          
          dbGetQuery(con(), query2)
          
          dbGetQuery(con(), "DROP TABLE IF exists set3;")
          query3<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        LIMIT ",x(),"
                        ;")
          dbGetQuery(con(), query3) 
          
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
          query4<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$hostID,");")
          
          
          dbGetQuery(con(), query4)                
          dbGetQuery(con(), "DROP TABLE IF exists set4;")
          query5<-paste0("create temporary table set4
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$pathogenID,") DESC
                        ;")
          
          dbGetQuery(con(), query5)                
          query6<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set4
                         UNION
                          SELECT * FROM set3
                         LIMIT ",x(),";")
          
          
          qmain1 <- dbGetQuery(con(),query6)
          #edges <- distinct(qmain2)
          
          
          edges <- distinct(qmain1)
          
          e$data <- as.data.frame(qmain1)
         
          
        })
      })
    })
    
    ############################___________________________________________________    
    
  }
  #####################################################################################
  
  else if(!is.null(input$multiProtein) && input$hostID == input$pathogenID){
    x <- reactiveVal(100)
    withProgress(message = 'Please wait for protein interaction...', {
      dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
      query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
      
      
      dbGetQuery(con(), query1)
      
      
      query2<-paste0("
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        LIMIT ",x(),"
                        ;")
      
      
      qmain1 <- dbGetQuery(con(),query2)
      #edges <- distinct(qmain2)
      
      e$data <- as.data.frame(qmain1)
    })
    
    observeEvent(input$MoReBtn,{
      
      withProgress(message = 'More protein interaction...', {
        x(x()+150)
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
          
          dbGetQuery(con(), query1)
          
          
          query2<-paste0("
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        LIMIT ",x(),"
                        ;")
          
          
          qmain1 <- dbGetQuery(con(),query2)
          #edges <- distinct(qmain2)
          
          e$data <- as.data.frame(qmain1)
        })
      })
    })
    
    observeEvent(input$LessBtn,{
      withProgress(message = 'Less protein interaction...', {
        x(abs(x()-140))
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
          
          dbGetQuery(con(), query1)
          
          
          query2<-paste0("
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        LIMIT ",x(),"
                        ;")
          
          
          qmain1 <- dbGetQuery(con(),query2)
          #edges <- distinct(qmain2)
          
          e$data <- as.data.frame(qmain1)
          
        })
      })
    })
  }  
  #####################################################################################
  else{
    x <- reactiveVal(100)
    
    withProgress(message = 'Please wait for protein interaction...', {
      dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
      query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
      
      
      dbGetQuery(con(), query1)
      
      dbGetQuery(con(), "DROP TABLE IF exists set2;")
      query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
      
      dbGetQuery(con(), query2)                
      dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
      query3<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
      
      dbGetQuery(con(), query3)                
      dbGetQuery(con(), "DROP TABLE IF exists set3;")
      query4<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
      
      dbGetQuery(con(), query4)                
      query5<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set3
                        LIMIT ",x(),";")
      
      
      qmain1 <- dbGetQuery(con(),query5)
      #edges <- distinct(qmain2)
      
      e$data <- as.data.frame(qmain1)
      
    })
    
    
    observeEvent(input$MoReBtn,{
      withProgress(message = 'More protein interaction...', {
        x(x()+150)
        
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
          dbGetQuery(con(), query1)
          
          dbGetQuery(con(), "DROP TABLE IF exists set2;")
          query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
          
        
          dbGetQuery(con(), query2)
          
          dbGetQuery(con(), "DROP TABLE IF exists set3;")
          query3<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        LIMIT ",x(),"
                        ;")
          dbGetQuery(con(), query3) 
          
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
          query4<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
                  
          dbGetQuery(con(), query4)                
          dbGetQuery(con(), "DROP TABLE IF exists set4;")
          query5<-paste0("create temporary table set4
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
          
          dbGetQuery(con(), query5)                
          query6<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set4
                         UNION
                          SELECT * FROM set3
                         LIMIT ",x(),";")
          
          
          qmain1 <- dbGetQuery(con(),query6)
          #edges <- distinct(qmain2)
          
          
          edges <- distinct(qmain1)
          
          e$data <- as.data.frame(qmain1)
          
        })
      })
    })
    
    
    observeEvent(input$LessBtn,{
      withProgress(message = 'Less protein interaction...', {
        x(abs(x()-140))
        
        reactiveVal({
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")
          query1<-paste0("
         CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
          dbGetQuery(con(), query1)
          
          dbGetQuery(con(), "DROP TABLE IF exists set2;")
          query2<-paste0("create temporary table set2
            SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
          
          
          dbGetQuery(con(), query2)
          
          dbGetQuery(con(), "DROP TABLE IF exists set3;")
          query3<-paste0("create temporary table set3
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        temptable.Interaction_ID,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1,
                        temptable.GeneInfoIDto,
                        temptable.Gene_2,
                        temptable.Tax_2,
                        temptable.Organism_Name_2
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,")
                        LIMIT ",x(),"
                        ;")
          dbGetQuery(con(), query3) 
          
          dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")                  
          query4<-paste0("CREATE VIEW temptable AS
            SELECT 
                        Gene2Gene.GeneInfoIDfrom AS GeneInfoIDfrom,
                        GeneInfo.Gene_Symbol AS Gene_1,
                        Gene2Gene.Interaction_ID AS Interaction_ID, 
                        Gene2Gene.GeneInfoIDto AS GeneInfoIDto,
                        GeneInfo.Taxonomy_ID AS Tax_1,
                        Taxonomy.Organism_Name AS Organism_Name_1
                        FROM 
                        Gene2Gene 
                        INNER JOIN GeneInfo ON Gene2Gene.GeneInfoIDfrom = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE  
                        GeneInfo.Taxonomy_ID IN (",input$pathogenID,");")
          
          
          dbGetQuery(con(), query4)                
          dbGetQuery(con(), "DROP TABLE IF exists set4;")
          query5<-paste0("create temporary table set4
                  SELECT 
                        temptable.GeneInfoIDfrom,
                        temptable.Gene_1,
                        temptable.Interaction_ID,
                        temptable.Tax_1,
                        temptable.Organism_Name_1,
                        temptable.GeneInfoIDto,
                        GeneInfo.Gene_Symbol AS Gene_2,
                        GeneInfo.Taxonomy_ID AS Tax_2,
                        Taxonomy.Organism_Name AS Organism_Name_2
                        
                        FROM 
                        temptable 
                        INNER JOIN GeneInfo ON temptable.GeneInfoIDto = GeneInfo.UniProtKB
                        INNER JOIN Taxonomy ON GeneInfo.Taxonomy_ID = Taxonomy.ID
                        WHERE 
                        GeneInfo.Taxonomy_ID IN (",input$hostID,")
                        ORDER BY GeneInfo.Taxonomy_ID IN (",input$hostID,") DESC
                        ;")
          
          dbGetQuery(con(), query5)                
          query6<-paste0("SELECT * FROM set2 UNION
                          SELECT * FROM set4
                         UNION
                          SELECT * FROM set3
                         LIMIT ",x(),";")
          
          
          qmain1 <- dbGetQuery(con(),query6)
          #edges <- distinct(qmain2)
          
          
          edges <- distinct(qmain1)
         
          e$data <- as.data.frame(qmain1)
          
        })
      })
    })
  }
  
  
})
