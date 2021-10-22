# Opioid Overdose Prevention Research Project
# Jacob Englert
# Dept. of Biostatistics and Bioinformatics, Emory University
# Purpose: Conduct a power analysis of various clustering/smoothing schemes

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(broom)    # tidy
library(rgdal)    # readOGR
library(spdep)    # poly2nb, mat2listw
library(rgeos)    # gCentroid
library(smerc)    # cepp.test, bn.test, scan.test
library(INLA) #install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Manually Edit smerc::stat.poisson ---------------------------------------
# tall[good][which(tall[good] < 0)] = 0 <- copy this at the end of the function
# trace(stat.poisson, edit = TRUE)

# Retrieve FindClusters Function from Cluster.R ---------------------------
source("Cluster.R")

# Load GA County Shapefile and create Neighborhood Objects ----------------
GA <- readOGR(dsn = 'Data/GA', verbose = F)
GA.nb <- poly2nb(GA, queen = F) # neighborhood object
W <- nb2mat(GA.nb, style = "B") # matrix form (needed for CARBayes/INLA)
coords <- as.matrix(gCentroid(GA, byid = T)@coords) # centroid coordinates
d <- as.matrix(dist(coords)) # distance matrix

# Create copy of neighborhood object with central region included
nb <- list()
for(i in 1:length(GA)) nb[[i]] <- c(i, GA.nb[[i]])

# Simulate Clustering Scenarios -------------------------------------------

# Where are observed SMR worst?
p1 <- GA %>% st_as_sf() %>% ggplot() + geom_sf(aes(fill = cut(SMR,5))) +
  scale_fill_brewer(type = 'div', palette = "PiYG") +
  labs(fill = 'SMR')

# Distribution of population?
p2 <- GA %>% st_as_sf() %>% ggplot() + geom_sf(aes(fill = cut(log(GA$pop),5))) +
  scale_fill_brewer(type = 'div', palette = "PiYG") +
  labs(fill = 'log(Pop)')

gridExtra::grid.arrange(p1, p2, nrow = 1)

# Specify center of circular cluster(s)
#cen.counties <- c('Gilmer','Columbia','Charlton','Butts') # 2017
cen.counties <- c('DeKalb','Early','Gilmer') # 2020
cen.sizes <- c(6, 6, 6)

# Get k-nearest neighbors for each region
knnb <- knn(coords, longlat = TRUE, k = 10, d = d)
nclus <- length(cen.counties)
center <- GA$OBJECTID[GA$NAME10 %in% cen.counties][order(na.omit(match(GA$NAME10, cen.counties)))]

# Create vector of k nearest neighbors for each center specified
cluster <- list()
for(i in 1:nclus) cluster[[i]] <- knnb[[center[i]]][1:cen.sizes[i]]

# Detailed cluster label vector
cluster.det <- factor(x = rep('None', 159),
                      levels = c(sapply(cluster, function(x) GA$NAME10[x][1]), 'None'))
for(i in 1:nclus) cluster.det[cluster[[i]]] <- cen.counties[i]

# Set delta and epsilon (Waller, Hill, Rudd 2006)
D <- GA$OBJECTID %in% unlist(cluster) # 1 if in cluster, 0 otherwise
E <- c(0.5,1,1.5,2) # increase in risk for regions within the cluster

# Confirm names and locations of counties included in clusters
sapply(cluster, function(x) GA$NAME10[x])
GA %>% st_as_sf() %>% ggplot() + geom_sf(aes(fill = cluster.det))

# Simulate datasets, keeping the total number of cases fixed
sims <- 100
set.seed(4303)
sim.data <- sapply(E, function(x) rmultinom(n = sims, size = sum(GA$mort),
                      prob = GA$pop*(1+D*x)/sum(GA$pop*(1+D*x))), simplify = 'array')

# Confirm regions in the cluster are showing higher SMR
sim.data %>%
  matrix(ncol = prod(dim(sim.data)[2:3])) %>%
  as.data.frame() %>%
  mutate(county = GA$NAME10) %>%
  pivot_longer(cols = -county) %>%
  mutate(e = factor(rep(rep(E, each = sims), times = length(GA))),
         exp = rep(sum(GA$mort)*(GA$pop)/sum(GA$pop), each = sims*length(E))) %>%
  mutate(SMR = value/exp,
         cluster = rep(cluster.det, each = sims*length(E))) %>%
  group_by(e, county) %>%
  filter(median(SMR) > 1) %>%
  ggplot(aes(x = SMR, y = reorder(county, SMR, FUN = median), fill = cluster)) +
  geom_boxplot() +
  facet_wrap(~e, scales = 'free')

# Check size of cluster to help determine parameters

# Print observed (true) cluster members, total populations and cases
sapply(cluster, function(x) GA$NAME10[x])
sapply(cluster, function(x) sum(GA$pop[x]))
sapply(cluster, function(x) sum(GA$mort[x]))

# Print mean total observed (simulated) cases for each value of E
apply(sim.data[unlist(cluster),,], 3, function(x) sum(rowMeans(x)))

# Print mean cluster observed (simulated) cases for each value of E
sapply(cluster, function(x) apply(sim.data[x,,], 3, function(y) sum(rowMeans(y))))

# Run Simulation ----------------------------------------------------------

# Create storage
temp <- array(dim = c(length(E), sims, 4, 3, 159), 
              dimnames = list(E = E,
                              sim = 1:sims,
                              smooth = list('Obs','Adj','ID','BYM'),
                              test = list('TB','BN','KD'),
                              county = GA$NAME10))

# Specify Cluster Test Parameters
# Current Method:
#  - n*: the maximum aggregated population among all clusters
#  - c*: the minimum aggregated cases among all clusters for each E (vector)
#  - min.cases: the minimum aggregated cases among all clusters for each E (vector)
#  - kd.ub: n* / total population
tb.nstar <- max(sapply(cluster, function(x) sum(GA$pop[x])))
bn.cstar <- apply(sapply(cluster, function(x) apply(sim.data[x,,], 3, function(y) sum(rowMeans(y)))), 1, min)
#min.cases <- floor(min(apply(cluster, 2, function(x) sum(GA$pop[x])))/sum(GA$pop)*sum(GA$mort)) # Maybe change this depending on E?

# Run Simulation
set.seed(1919)
for(e in 1:length(E)){
  
  for(s in 1:sims){

    if(s == 1) start <- Sys.time()
    s.start <- Sys.time()
  
    # 0) Get new "observed" dataset and update constant risk expectation
    obMort <- sim.data[,s,e]
    GA$obMort <- obMort
    GA$eMort <- GA$pop * (sum(obMort)/sum(GA$pop))
    
    # 1) Get smoothed rates for each county
    asMort <- numeric(length(GA))
    dsMort <- numeric(length(GA))
    for(j in 1:length(GA)){
    
      # Adjacency Smoothing
      asMort[j] <- mean(obMort[nb[[j]]] / GA$pop[nb[[j]]]) * GA$pop[j]
    
      # Centroid Distance Smoothing
      w.j <- d[nb[[j]], j]
      w.j[1] <- min(w.j[w.j > 0]) # Weight the region of interest
      dsMort[j] <- sum((1/w.j)*(obMort[nb[[j]]] / GA$pop[nb[[j]]])) / sum(1/w.j) * GA$pop[j]
    }
  
    # Bayesian Smoothing (INLA)
    bmod <- inla(obMort ~ 1 + f(OBJECTID, model = "besag", graph = W),
                 data = as.data.frame(GA), E = GA$eMort, family = "poisson",
                 control.predictor = list(compute = TRUE),
                 control.compute = list(dic = FALSE, waic = FALSE))
  
    bsMort <- bmod$summary.fitted.values[, '0.5quant']*GA$eMort
  
    # 2) Apply Cluster Detection Tests
    c1 <- quiet(FindClusters(o = obMort, p = GA$pop, tb.nstar = tb.nstar, kd.min = 2, bn.cstar = bn.cstar[e]))
    c2 <- quiet(FindClusters(o = asMort, p = GA$pop, tb.nstar = tb.nstar, kd.min = 2, bn.cstar = bn.cstar[e]))
    c3 <- quiet(FindClusters(o = dsMort, p = GA$pop, tb.nstar = tb.nstar, kd.min = 2, bn.cstar = bn.cstar[e]))
    c4 <- quiet(FindClusters(o = bsMort, p = GA$pop, tb.nstar = tb.nstar, kd.min = 2, bn.cstar = bn.cstar[e]))

    # 3) Store Results
    temp[e,s,1,1,] <- GA$OBJECTID %in% c1$tb$id
    temp[e,s,1,2,] <- GA$OBJECTID %in% c1$bn$id
    temp[e,s,1,3,] <- GA$OBJECTID %in% c1$kd$id
  
    temp[e,s,2,1,] <- GA$OBJECTID %in% c2$tb$id
    temp[e,s,2,2,] <- GA$OBJECTID %in% c2$bn$id
    temp[e,s,2,3,] <- GA$OBJECTID %in% c2$kd$id
    
    temp[e,s,3,1,] <- GA$OBJECTID %in% c3$tb$id
    temp[e,s,3,2,] <- GA$OBJECTID %in% c3$bn$id
    temp[e,s,3,3,] <- GA$OBJECTID %in% c3$kd$id
  
    temp[e,s,4,1,] <- GA$OBJECTID %in% c4$tb$id
    temp[e,s,4,2,] <- GA$OBJECTID %in% c4$bn$id
    temp[e,s,4,3,] <- GA$OBJECTID %in% c4$kd$id

    # Print update every 10 datasets
    #if(s %% 10 == 0) print(paste(rep('=',10), s, 'datasets processed', rep('=',10)))
    print(paste('E =', E[e], ', Dataset', s, 'finished processing.'))
    print(paste0('Dataset time: ', round(Sys.time() - s.start, 2), ' minutes.',
                 ' Total time: ', round(Sys.time() - start, 2), ' minutes.'))
    print("==================================================================")
  }
  print(paste("E =", E[e], "completed."))
}

# Calculate Accuracy Measures ---------------------------------------------

# E, sim, smooth, test, county

# Mean Overall Sensitivity (ability to detect positives as positive)
sens <- apply(temp[,,,,D], c(1,3,4), function(x) mean(rowSums(x)/sum(D))) # Verified

# Mean Specificity (ability to detect negatives as negative)
spec <- apply(!temp[,,,,!D], c(1,3,4), function(x) mean(rowSums(x)/sum(!D))) # Verified

# Mean PPV (probability a positive test is actually positive)
#P.C <- sum(D)/length(D)               # Probability of cluster
#P.D <- apply(temp[,,,,], c(1,3,4), function(x) mean(rowSums(x)/length(D))) # Probability of detecting a cluster
#PPV <- sens*P.C/P.D

PPV <- apply(temp, c(1,3,4), function(x) mean(apply(x, 1, function(y) mean(which(y) %in% which(D))))) # Verified?

# Mean NPV (probability a negative test is actually negative)
#NPV <- spec*(1 - P.C)/(1 - P.D)

NPV <- apply(temp, c(1,3,4), function(x) mean(apply(x, 1, function(y) mean(which(!y) %in% which(!D))))) # Verified?

# Exact Power (Schuldeln et al. 2021)
EP <- apply(temp, c(1,3,4), function(x) mean(apply(x[,D], 1, all) & apply(x[,!D] == FALSE, 1, all))) # Verified

# Minimum Power (Schuldeln et al. 2021)
MP <- apply(temp[,,,,D], c(1,3,4), function(x) mean(rowSums(x) > 0)) # Verified, only for one cluster

# Correct Classification
CC <- apply(temp, c(1,3,4), function(x) mean(colMeans(apply(x, 1, function(y) y == D)))) # Verified


# Visual Comparison -------------------------------------------------------

# Create data frame of results
res <- data.frame(rbind(matrix(sens, nrow = prod(dim(sens)[1:2])),
                        matrix(spec, nrow = prod(dim(spec)[1:2])),
                        matrix(PPV, nrow = prod(dim(PPV)[1:2])),
                        matrix(NPV, nrow = prod(dim(NPV)[1:2])),
                        matrix(EP, nrow = prod(dim(EP)[1:2])),
                        matrix(MP, nrow = prod(dim(MP)[1:2])),
                        matrix(CC, nrow = prod(dim(CC)[1:2])))) %>%
  mutate(Smoothing = rep(c('Obs','Adj','ID','BYM'), each = 4, times = 7), # 4 smooth x 7 measures
         Measure = rep(c('Sens','Spec','PPV','NPV','EP','MP','CC'), each = 4*length(E)),
         E = rep(E, times = 4*7)) %>% # 4 smooth x 7 measures
  mutate(Smoothing = factor(Smoothing, levels = c('Obs','Adj','ID','BYM'))) %>%
  rename(TB = X1, BN = X2, KD = X3) %>%
  pivot_longer(cols = c(TB, BN, KD), names_to = 'Test') %>%
  mutate(Test = case_when(Test == 'TB' ~ 'Turnbull (1990)',
                          Test == 'KD' ~ 'Kulldorff (1997)',
                          Test == 'BN' ~ 'Besag-Newell (1991)'))

# Save Results
write_csv(res, file = paste0('Simulation Results/', 
                             paste0(sort(substr(cen.counties, 1, 3)), lengths(cluster), collapse = '_'), 
                             '_sim', sims, '.csv'))

# Plot Accuracy Results
res %>%
  #filter(!(Test == 'Turnbull (1990)')) %>%
  ggplot(aes(x = E, y = value, color = Smoothing)) +
  geom_point() +
  geom_line() +
  facet_grid(Measure~Test, scales = 'free')
  theme_minimal()


  

# Extras ------------------------------------------------------------------

# Plot frequency captured
test <- as.data.frame(ftable(apply(temp, c(1,3,4,5), mean)))
test2 <- tidy(GA, region = "NAME10") %>%
  left_join(test, by = c("id" = "county")) %>%
  rename(county = id)

test2 %>% 
  filter(test == 'KD') %>%
  ggplot(aes(x = long, y = lat, group = county)) +
  geom_polygon(aes(fill = Freq), color = 'black') +
  coord_map("tetra") +
  #scale_fill_brewer(type = 'div', palette = "PiYG") +
  facet_grid(E~smooth)
  
# Mean Cluster Sensitivity
c.sens <- matrix(nrow = 4*length(center), ncol = 3)
for(i in 1:length(center)){
  c.D <- GA$OBJECTID %in% cluster[,i]
  c.sens[((i*4)-3):(i*4), 1:3] <- apply(temp[,,,c.D], c(2,3), function(x) mean(rowSums(x)/sum(c.D)))
}
c.sens
  
# Bayesian smoothing with CARBayes
# library(CARBayes)
# cbmod <- quiet(S.CARleroux(formula = obMort ~ offset(log(pop)), data = GA, 
#                            family = "poisson", W = W,
#                            burnin = 20000, n.sample = 100000, thin = 40))
#                           #burnin = 100000, n.sample = 300000, thin = 100))
# bsMort <- cbmod$fitted.values
  
# 2) Apply cluster detection tests
# c1 <- quiet(suppressMessages(FindClusters(o = obMort, tb.nstar = 4000000)))
# c2 <- quiet(suppressMessages(FindClusters(o = asMort, tb.nstar = 4000000)))
# c3 <- quiet(suppressMessages(FindClusters(o = dsMort, tb.nstar = 4000000)))
# c4 <- quiet(suppressMessages(FindClusters(o = bsMort, tb.nstar = 4000000)))