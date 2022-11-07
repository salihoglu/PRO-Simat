# Create empty csv to avoid null exception
system('touch jimena_time_series_data.csv')
# Build jimena application with external Jimena library
system("javac -classpath jimena-app/jimena.jar:jimena-app jimena-app/App.java")