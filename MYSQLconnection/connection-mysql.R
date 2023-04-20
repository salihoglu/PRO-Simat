library(RMySQL)

con <-reactive({
  dbConnect(RMySQL::MySQL(),
                 user = 'root',
                 host = 'localhost',
                 port = 3306,
                 dbname='*****',
                 password = "***")



}) 
#on.exit(dbDisconnect(con), add = TRUE)

