# Opioid Overdose Prevention Research Project
# Jacob Englert
# Dept. of Biostatistics and Bioinformatics, Emory University
# Purpose: Identify clusters using various techniques

# Load Packages -----------------------------------------------------------
# library(tidyverse)
# library(rgdal)      # readOGR
# library(spdep)      # poly2nb, mat2listw, moran, geary
# library(rgeos)      # gCentroid
# library(broom)      # tidy
# library(smerc)      # cepp.test, bn.test, scan.test
# library(SpatialEpi) # (besag_newell, kulldorff - not using)

# Load Estimates obtained from Smoothing.R --------------------------------
#GA <- readOGR(dsn = 'GA', verbose = F)

# Ingredients
# coords <- as.matrix(gCentroid(GA, byid = T)@coords)
# pop <- GA$pop
# obs <- GA$mort
# exp <- GA$eMort

# Create function to apply all three methods
FindClusters <- function(cds = coords, p = pop, o = obs, e = sum(o)/sum(p)*p,
                         tb.nstar = 100000, 
                         bn.cstar = tb.nstar/sum(p)*sum(o), 
                         kd.ub = tb.nstar/sum(GA$pop), kd.min = 2,
                         sig.level = 0.1, nsims = 999, lonlat = TRUE){
  
  # Turnbull (1990)
  tb.test <- cepp.test(coords = cds, cases = o, pop = p, ex = e,
                       nstar = tb.nstar, alpha = sig.level, nsim = nsims,
                       longlat = lonlat, simdist = "poisson")
  
  # Besag and Newell (1991)
  bn.test <- bn.test(coords = cds, cases = o, pop = p, ex = e,
                     cstar = bn.cstar, alpha = sig.level, 
                     longlat = lonlat)
  
  # Kulldorff (1997)
  kd.test <- quiet(scan.test(coords = cds, cases = o, pop = p, ex = e,
                       nsim = nsims, alpha = sig.level, ubpop = kd.ub, 
                       longlat = lonlat, min.cases = kd.min, simdist = "poisson"))
  
  # Store region IDs, cluster association, and associated cluster p-values
  for(test in c("tb", "bn", "kd")){
    
    # Retrieve test object
    t <- get(paste0(test,'.test'))
    
    # Vector(s) containing cluster sizes
    if(length(t$clusters) == 1) l <- length(t$clusters[[1]][["locids"]])
    if(length(t$clusters) > 1) l <- lengths(sapply(t$clusters, "[[", "locids", simplify = FALSE))
    
    # Define output vector names
    t.ids <- paste0(test, '.ids')
    t.clus <- paste0(test, '.clus')
    t.pv <- paste0(test, '.pv')
    
    # Instantiate output vectors
    assign(t.ids, numeric(0))
    assign(t.clus, numeric(0))
    assign(t.pv, numeric(0))
    
    # Assign *.ids, *.clus, and *.pv vectors for each test
    for(clus in 1:length(t$clusters)){
      
      # Only add non-overlapping clusters
      if(clus > 1 & length(intersect(get(t.ids), sapply(t$clusters, "[[", "locids", simplify = FALSE)[[clus]])) > 0) break
      
      assign(t.ids, append(get(t.ids), t$clusters[[clus]][["locids"]]))
      assign(t.clus, append(get(t.clus), rep(clus, times = l[clus])))
      assign(t.pv, append(get(t.pv), rep(sapply(t$clusters, "[[", "pvalue")[[clus]], times = l[clus])))
    }           
          
    # Now remove clusters that are not significant (by default, most significant is returned and we DON'T want that)
    assign(t.ids, get(t.ids)[which(get(t.pv) < sig.level)])
    assign(t.clus, get(t.clus)[which(get(t.pv) < sig.level)])
    assign(t.pv, get(t.pv)[which(get(t.pv) < sig.level)])
    
  }

  return(list(tb = list(id = tb.ids, clus.id = tb.clus, pvalue = tb.pv), 
              bn = list(id = bn.ids, clus.id = bn.clus, pvalue = bn.pv),
              kd = list(id = kd.ids, clus.id = kd.clus, pvalue = kd.pv)))
}

# # Run tests under each smoothing scenario
# set.seed(1)
# c1 <- FindClusters()
# c2 <- FindClusters(o = GA$asMort)
# c3 <- FindClusters(o = GA$dsMort)
# c4 <- FindClusters(o = GA$bsMort)
# 
#  
# # # Merge variables into GA SpatialPolygonsDataFrame
# GA$tbclus1 <- (GA$OBJECTID %in% c1$tb$id)
# GA$tbclus2 <- (GA$OBJECTID %in% c2$tb$id)
# GA$tbclus3 <- (GA$OBJECTID %in% c3$tb$id)
# GA$tbclus4 <- (GA$OBJECTID %in% c4$tb$id)
# 
# GA$bnclus1 <- (GA$OBJECTID %in% c1$bn$id)
# GA$bnclus2 <- (GA$OBJECTID %in% c2$bn$id)
# GA$bnclus3 <- (GA$OBJECTID %in% c3$bn$id)
# GA$bnclus4 <- (GA$OBJECTID %in% c4$bn$id)
# 
# GA$kdclus1 <- (GA$OBJECTID %in% c1$kd$id)
# GA$kdclus2 <- (GA$OBJECTID %in% c2$kd$id)
# GA$kdclus3 <- (GA$OBJECTID %in% c3$kd$id)
# GA$kdclus4 <- (GA$OBJECTID %in% c4$kd$id)
# 
# # Put data into long format for lattice plotting
# GA.clus <- tidy(GA, region = "NAME10") %>%
#   left_join(GA@data, by = c("id" = "NAME10")) %>%
#   rename(county = id) %>%
#   pivot_longer(cols = contains('clus'), names_to = c("C.Method","S.Method"),
#                names_sep = 'clus') %>%
#   mutate(MR = case_when(S.Method == 1 ~ MR,
#                         S.Method == 2 ~ asMR,
#                         S.Method == 3 ~ dsMR,
#                         S.Method == 4 ~ bsMR),
#          S.Method = factor(S.Method, levels = 1:4,
#                            labels = c('Observed', 'Adjacency Smoothing',
#                                       'Inverse-Distance Smoothing',
#                                       'Bayesian Smoothing')),
#          C.Method = factor(C.Method, levels = c('tb','bn','kd'),
#                           labels = c('Turnbull', 'Besag-Newell', "Kuldorff"))) %>%
#   mutate(SMR = MR/1000*pop/eMort) %>%
#   mutate(cMR = cut(MR, breaks = seq(0, 0.6, 0.1), include.lowest = TRUE),
#          cSMR = cut(SMR, breaks = quantile(SMR, seq(0,1,.1)), include.lowest = TRUE))
# 
# 
# 
# 
# 
# # Plot lattice
# latPlot <- GA.clus %>%
#   mutate(cluster = ifelse(value, 'Cluster', 'Non-Cluster')) %>%
#   ggplot(aes(x = long, y = lat, group = county)) +
#   geom_polygon(aes(fill = cluster, alpha = cSMR), color = 'black') +
#   coord_map("tetra") +
#   scale_fill_manual(values = c('yellow','navy'), breaks = c('Cluster')) +
#   facet_grid(C.Method ~ S.Method) +
#   labs(title = 'Opioid Overdose Mortality Clusters in Georgia (2017)',
#        fill = '',
#        alpha = 'Standardized Morbitidy Ratio') +
#   theme_bw() +
#   theme(axis.text = element_blank(), axis.ticks = element_blank(),
#         axis.title = element_blank())
# 
# ggsave(plot = latPlot, filename = "Figures/GA_Mortality_Clusters_Lattice_Plot.png",
#        width = 10, height = 6.5)


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# stat.poisson = function(yin, yout, ein, eout, a = 0, shape = 1) {
#   # determine if there will be any problematic statistics
#   good = which(yin > 0)
#   # create vector for storage
#   tall = numeric(length(yin))
#   # log ratio observed/expected (in and out) for good locations
#   lrin =  log(yin[good]) - log(ein[good])
#   lrout = log(yout[good]) - log(eout[good])
#   # compute statistics for good locations
#   tall[good] = yin[good] * lrin + yout[good] * lrout
#   # if indicator not satisfied, set to 0
#   tall[good][lrin < lrout] = 0
#   # Adjust tiny values close to 0 that R calls negative to be 0 to avoid errors
#   tall[good][tall[good] < 0] = 0
#   if (a > 0) {
#     i = which(shape > 1)
#     tall[i] = tall[i] * (4 * shape[i] / (shape[i] + 1) ^ 2) ^ a
#   }
#   return(tall)
# }


# # Global Indexes of Spatial Autocorrelation -------------------------------
# 
# # Moran's I
# 
# # Construct spatial weight matrix
# nb <- poly2nb(GA, queen = F)
# #d <- as.matrix(dist(gCentroid(ests, byid = T)@coords))
# n <- length(GA)
# W <- matrix(0, n, n)
# for(i in 1:n) W[i, nb[[i]]] <- 1#/(length(nb[[i]]))
# 
# set.seed(2187)
# o <- rbind(matrix(rpois(n = 999*n, lambda = GA$eMort), 999, byrow = TRUE), obs)
# getMoran <- function(x, w){
#   moran(x = x, listw = mat2listw(w), n = length(x), S0 = sum(w))$I
# }
# 
# I <- apply(o, 1, getMoran, w = W)
# hist(I)
# abline(v = I[1000], col = 'red')
# mean(I >= I[1000]) # P-value
# 
# 
# # Geary's c
# getGeary <- function(x, w){
#   geary(x = x, listw = mat2listw(w), n = length(x), n1 = length(x) - 1, S0 = sum(w))$C
# }
# C <- apply(o, 1, getGeary, w = W)
# hist(C)
# abline(v = C[1000], col = 'red')
# mean(C <= C[1000]) # P-value
# 
# 
# 
# # Local Indicators of Spatial Association ---------------------------------
# 
# # Local Moran's I
# 
# lisa <- localmoran_perm(x = GA$MR, listw = mat2listw(W), nsim = 999)
# 
# GA@data$lmipv <- lisa[,5]