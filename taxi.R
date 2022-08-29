library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(igraph)
library(chron)
source('src/gnr.R')
source('src/lnr.R')
source("src/kerFctn.R")

startDate <- '2020-04-11'
endDate <- '2020-10-01'
covidNYC <- read.csv("data/NYC_cases_by_day.csv")
covidNYC$date_of_interest <- as.Date(covidNYC$date_of_interest, format = "%m/%d/%Y")
covidNYC <- covidNYC[covidNYC$date_of_interest > startDate & covidNYC$date_of_interest < endDate, c('date_of_interest', 'MN_CASE_COUNT_7DAY_AVG')]
covidNYC$weekendHoliday <- as.numeric(is.weekend(covidNYC$date_of_interest))
covidNYC$weekendHoliday[covidNYC$date_of_interest%in%as.Date(c('2020-05-25', '2020-07-03', '2020-09-07'))] <- 1
pred <- covidNYC[, c('MN_CASE_COUNT_7DAY_AVG', 'weekendHoliday')]

neighborID <- c(4, 12, 13, 24, 41, 42, 43, 45, 48, 50, 
                68, 74, 75, 79, 87, 88, 90, 100, 103, 104, 
                105, 107, 113, 114, 116, 120, 125, 127, 128, 137, 
                140, 141, 142, 143, 144, 148, 151, 152, 153, 158, 
                161, 162, 163, 164, 166, 170, 186, 194, 202, 209, 
                211, 224, 229, 230, 231, 232, 233, 234, 236, 237, 
                238, 239, 243, 244, 246, 249, 261, 262, 263)
neighborName <- c("Alphabet City", "Battery Park", "Battery Park City", 
                  "Bloomingdale", "Central Harlem", "Central Harlem North", 
                  "Central Park", "Chinatown", "Clinton East", 
                  "Clinton West", "East Chelsea", "East Harlem North", 
                  "East Harlem South", "East Village", "Financial District North", 
                  "Financial District South", "Flatiron", "Garment District", 
                  "Liberty Island", "Ellis Island", "Governor's Island", 
                  "Gramercy", "Greenwich Village North", "Greenwich Village South", 
                  "Hamilton Heights", "Highbridge Park", "Hudson Sq", 
                  "Inwood", "Inwood Hill Park", "Kips Bay", 
                  "Lenox Hill East", "Lenox Hill West", "Lincoln Square East", 
                  "Lincoln Square West", "Little Italy/NoLiTa", "Lower East Side", 
                  "Manhattan Valley", "Manhattanville", "Marble Hill", 
                  "Meatpacking/West Village West", "Midtown Center", "Midtown East", 
                  "Midtown North", "Midtown South", "Morningside Heights", 
                  "Murray Hill", "Penn Station/Madison Sq West", "Randalls Island", 
                  "Roosevelt Island", "Seaport", "SoHo", 
                  "Stuy Town/Peter Cooper Village", "Sutton Place/Turtle Bay North", "Times Sq/Theatre District", 
                  "TriBeCa/Civic Center", "Two Bridges/Seward Park", "UN/Turtle Bay South", 
                  "Union Sq", "Upper East Side North", "Upper East Side South", 
                  "Upper West Side North", "Upper West Side South", "Washington Heights North", 
                  "Washington Heights South", "West Chelsea/Hudson Yards", "West Village", 
                  "World Trade Center", "Yorkville East", "Yorkville West")
names(neighborID) <- neighborName
neighborID <- neighborID[!is.element(neighborID, 103:105)]# exclude the islands
# 66 neighborhoods in Manhattan Borough
# Financial District: 12, 13, 87, 88, 209, 261
# Penn Station: 186
# Grand Central Station: 170
zoneID <- list(c(12, 13, 87, 88, 209, 231, 261), 
               c(113, 114, 125, 144, 158, 211, 249), 
               c(4, 45, 79, 148, 232), 
               c(48, 50, 68, 90, 246), 
               c(100, 161, 163, 164, 186, 230, 234),
               c(107, 137, 162, 170, 224, 229, 233), 
               c(24, 142, 143, 151, 238, 239), 
               c(140, 141, 202, 236, 237, 262, 263), 
               c(116, 152, 166), 
               c(41, 42), 
               c(74, 75, 194), 
               c(120, 127, 128, 153, 243, 244), 
               c(43))# based on community districts
m <- length(zoneID)
names(zoneID) <- c(101:112, 164)

taxi <- list()
for(i in 4:9){
  taxi[[i]] <- read.csv(paste0("data/yellowTripdata2020//yellow_tripdata_2020", '-', paste0('0', i), '.csv'))
}
for(i in 4:9){
  taxi[[i]] <- taxi[[i]] %>% 
    drop_na(c(tpep_pickup_datetime, PULocationID, DOLocationID, passenger_count)) %>% # remove rows with missing values
    filter((PULocationID %in% neighborID) & (DOLocationID %in% neighborID)) %>% 
    mutate(date = substr(tpep_pickup_datetime, start = 1, stop = 10)) %>%
    filter((date > startDate) & (date < endDate)) %>%
    select(date, PULocationID, DOLocationID, weight = passenger_count)
}
taxi <- bind_rows(taxi[[4]], taxi[[5]], taxi[[6]], taxi[[7]], taxi[[8]], taxi[[9]])
# 01-02, 01-04, 01-06~01-17
taxi$PULocationID <- apply(sapply(zoneID, function(zoneIDi) is.element(taxi$PULocationID, zoneIDi)), 1, which)
taxi$DOLocationID <- apply(sapply(zoneID, function(zoneIDi) is.element(taxi$DOLocationID, zoneIDi)), 1, which)
taxi <- taxi %>% filter(date > startDate & date < endDate)

taxigl <- taxi %>% dlply(.(date), function(taxi) {
  g <- graph_from_data_frame(d = taxi %>% dplyr::select(PULocationID, DOLocationID, weight),
                             vertices = 1:m, directed = FALSE)
  g <- simplify(g)# returns a simple graph (remove self-loops and combine multi-edges)
  gl <- graph.laplacian(g, sparse = FALSE)
  dimnames(gl) <- list(names(zoneID), names(zoneID))
  gl
})

res <- list()
# effect of COVID 19 new cases and weekend
# fix time as the start or the end
xOut <- cbind(MN_CASE_COUNT_7DAY_AVG = rep(c(50, 200, 400), 2),
              weekend = rep(0:1, each = 3))
res[[1]] <- gnr(taxigl, pred, xOut)
res[[2]] <- gnr(taxigl, pred, xOut, optns = list(metric = 'power', alpha = 0.5))# better
summary(res[[1]]$residuals)
summary(res[[2]]$residuals)

### Mean Square Prediction Error (MSPE) ###
library(doParallel)
NUM_CORES <- 6
cl <- makeCluster(NUM_CORES)
registerDoParallel(cl)

set.seed(123)
nSim <- 100
nFolds <- 10
n <- nrow(pred)
group <- sort(c(rep.int(1:nFolds, n%/%nFolds), rep.int(nFolds, n%%nFolds)))
mspe <- foreach(icount(nSim), .combine = 'rbind') %dopar% {
  f <- sample(group)
  predi <- split(pred, f)
  gl <- split(taxigl, f)
  spe1 <- 0
  spe2 <- 0
  for(j in 1:nFolds){
    testPred <- predi[[j]]
    trainPred <- Reduce(rbind, predi[-j])
    testGl <- gl[[j]]
    trainGl <- Reduce(append, gl[-j])
    fitGl1 <- gnr(trainGl, trainPred, testPred, 
                      optns = list(metric = 'power', alpha = 0.5))$predict
    fitGl2 <- gnr(trainGl, trainPred, testPred)$predict
    spe1 <- spe1 + sum(sapply(1:sum(f==j), function(k) sum((testGl[[k]]-fitGl1[[k]])^2)))
    spe2 <- spe2 + sum(sapply(1:sum(f==j), function(k) sum((testGl[[k]]-fitGl2[[k]])^2)))
  }
  # print(c(spe1, spe2)/n)
  c(spe1, spe2)/n
}
colMeans(mspe)
# 48898842/50750379 = 96.4%