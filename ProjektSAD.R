library(rstudioapi) 
library(dplyr)
library(tidyr)
library(psych)
library(factoextra)
library(cluster)
library(clusterSim)
library(ggplot2)
library(GGally)
library(eurostat)
library(sf)
library(writexl)

#Porządkowanie liniowe metodą hellwiga
hellwig <- function(data, decyzje, wagi){
  data[,decyzje == '-'] <- -data[,decyzje == '-']      #zamiana destymulant na stymulanty
  data <- scale(data)                                  #normalizacja danych
  wzorzec <- apply(data, 2, max)                       #wyznaczanie wzorca
  if(missing(wagi)){                                   #funckja licząca odległosci uwzgledniając wagi
    odleglosc <- function(x, y) 
                    {sqrt(sum((x - y)^2))}
  } else {
    odleglosc <- function(x,y) 
                    {sqrt(sum(diag(wagi)%*%(x-y)^2))}    
  }
  di <- apply(data, 1, odleglosc, wzorzec)            #odległości od wzorca
  d0 <- mean(di) + 2*sd(di)                           #odległosć "mozliwie daleka"
  si <- as.vector(1 - di/d0)                          #finalny wskaźnik
  ranking <- data.frame(row.names = rownames(data), rank = rank(-si), score = si)
  return(arrange(ranking, desc(score)))
}
#Porządkowanie liniowe metoda TOPSIS


topsis <- function(data, decyzje, wagi) {
  data[,decyzje == '-'] <- -data[,decyzje == '-']     #zamiana destymulant na stymulanty
  temp <- rownames(data) 
  data <-sapply(data, function(x) x / sqrt(sum(x^2))) #standaryzacja zmiennych
  rownames(data) <- temp
  wzorzec <- apply(data, 2, max)                      #wyznaczenie wzorca
  antywzorzec <- apply(data, 2, min)                  #wyznaczenie antywzorca
  if(missing(wagi)){                                  #funkcja obliczająca odleglości uzwględniając wagi
    odleglosc <- function(x, y) 
                    {sqrt(sum((x - y)^2))}
  } else {
    odleglosc <- function(x,y) 
                    {sqrt(sum(diag(wagi)%*%(x-y)^2))}    
  }           
  di_max <- apply(data, 1, odleglosc, wzorzec)        #odległości od wzorca
  di_min <- apply(data, 1, odleglosc, antywzorzec)    #odległości od antywzorca
  si <- as.vector(di_min/(di_min + di_max))           #finalny wskaźnik
  ranking <- data.frame(row.names = rownames(data), rank = rank(-si), score = si)
  return(arrange(ranking, desc(score)))
}

#PROJEKT
######WCZYTYWANIE DANYCH############
setwd(dirname(getActiveDocumentContext()$path))
path = getwd()

files <- list.files(path = path, pattern = "*.csv") #Wczytanie danych do listy “data” 

data <- list()
for (file in files){
  temp_name <- substring(file,1,nchar(file)-4) #usuwanie .csv (4 ostatnich liter) z nazwy pliku, będzie to nazwa kolumny
  file <- tibble(read.csv(paste(path,"/", file, sep="")))
  file <- file %>%
    dplyr::select(geo,OBS_VALUE) %>%
    rename(!!temp_name := OBS_VALUE)
  data <- append(data,list(file))
}
#Łączenie tabel (bierzemy część wspólną krajów oraz dodajemy do nich kolumny z danymi poszczególnych wskaźników)
data <- Reduce(
  function(x, y) inner_join(x, y, by = c("geo")), data
)
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- dplyr::select(data,-c("hicp_prc","satisfaction","geo")) #usuniecie zmiennych z cv < 10%

#wykres korelacji
par(mfcol=c(1,1))
ggpairs(data,
        axisLabels = "none",
        upper = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"),
        lower = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"),
        diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"))

#boxplot
ggplot(gather(as.data.frame(scale(data))),aes(x=key,y=value), axisLabels = "none") +
  geom_boxplot() +
  xlab(NULL) +
  theme_light()


#############Porzadkowanie liniowe############
decyzje = c("+","+","-","-","-")
r_hellwig<-hellwig(data,decyzje)
r_topsis<-topsis(data,decyzje)
#zapisanie wynikow do xlsx
#r_hellwig["country"] <- row.names(r_hellwig)
#r_topsis["country"] <- row.names(r_topsis)
#path <- getwd()
#write_xlsx(r_hellwig,"r_hellwig.xlsx")
#write_xlsx(r_topsis,"r_topsis.xlsx")


###########Analiza skupien############
#podziałowe
par(mfcol=c(1,2))
fviz_nbclust(data, pam, method = "silhouette")
fviz_nbclust(data, pam, method = "wss")
gr_medoids<-pam(scale(data),3)
fviz_cluster(gr_medoids, data) 

#hierarchiczne
distances <- dist(scale(data))
gr1<-hclust(distances,method = "ward.D")
plot(gr1)

indexes <- setNames(data.frame(matrix(ncol = 4, nrow = 1)), c("G1","G2","G3","S"))
for(i in 2:6){
  gr1cut<-cutree(gr1,i)
  indexes[i-1,"G1"] <- index.G1(scale(data),gr1cut)
  indexes[i-1,"G2"] <- index.G2(distances,gr1cut)
  indexes[i-1,"G3"] <- index.G3(distances,gr1cut)
  indexes[i-1,"S"] <- index.S(distances,gr1cut)
}
par(mfcol = c(2,2))
for (i in 1:4){
  plot(2:(length(indexes$G1)+1),indexes[,i],'b', xlab="Number of groups",
       ylab=NA, main = paste(colnames(indexes)[i],"index"))
}

gr_ward<-cutree(gr1,3)

##########wizualizacja grup na mapie##########
EU27 <- eu_countries %>% 
  filter(code != 'UK') %>% 
  dplyr::select(geo = code, name)

geometry <- get_eurostat_geospatial(resolution = 60, 
                                    nuts_level = 0, 
                                    year = 2021)
EU27 <- geometry %>% 
  dplyr::select(geo = NUTS_ID, geometry) %>% 
  inner_join(EU27, by = "geo") %>% 
  arrange(geo) %>% 
  st_as_sf()

#k-medoids
EU27kmedoids<- EU27 %>%
  mutate(group = as.factor(gr_medoids$clustering))

EU27kmedoids %>%
  ggplot(aes(fill = group)) +
  geom_sf() +
  geom_sf_text(aes(label = name), size = 3) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 65)) + 
  theme_void() 
#ward
EU27ward<- EU27 %>%
  mutate(group = as.factor(gr_ward))

EU27ward %>%
  ggplot(aes(fill = group)) +
  geom_sf() +
  geom_sf_text(aes(label = name), size = 3) +
  scale_x_continuous(limits = c(-10, 35)) +
  scale_y_continuous(limits = c(35, 65)) + 
  theme_void() 
