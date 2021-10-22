# Opioid Overdose Prevention Research Project
# Jacob Englert
# Dept. of Biostatistics and Bioinformatics, Emory University
# Purpose: Obtained smoothed estimates

# Load Packages -----------------------------------------------------------
library(tidyverse) 
library(readxl)     
library(rgdal)      # readOGR, writeOGR
library(spdep)      # poly2nb
library(rgeos)      # gCentroid
library(INLA)
library(broom)      # tidy

# Georgia County Shapefile (GARC) -----------------------------------------
GA <- readOGR(dsn = 'Data/Counties_Georgia', verbose = F)

# Mortality and Population Data (OASIS) -----------------------------------
mort <- read_xlsx('Data/mort2020.xlsx', skip = 1) %>%
  as.data.frame() %>%
  filter(!is.na(Geography)) %>%
  slice(1:(n() - 4)) %>%
  pivot_longer(cols = '2010':'2020', names_to = 'year', values_to = 'mort') %>%
  mutate(year = as.integer(year), mort = as.integer(mort)) %>%
  select(county = Geography, year, mort) %>%
  arrange(year, match(county, GA$NAME10))

pop <- read_xlsx('Data/pop2020.xlsx', skip = 1) %>%
  as.data.frame() %>%
  filter(!is.na(Geography)) %>%
  slice(1:(n() - 4)) %>%
  pivot_longer(cols = '2010':'2020', names_to = 'year', values_to = 'pop') %>%
  mutate(year = as.integer(year), pop = parse_number(pop)) %>%
  select(county = Geography, year, pop) %>%
  arrange(year, match(county, GA$NAME10))

# Check to see county order agrees between GARC and OASIS
table(GA$NAME10 == unique(mort$county))

# Calculate observed and expected mortality rates
GA@data <- left_join(GA@data,
  filter(left_join(mort, pop, by = c("year", "county")), year == 2020),
  by = c("NAME10" = "county")) %>%
  mutate(MR = mort/pop*1000,
         eMR = sum(mort)/sum(pop)*1000,
         eMort = sum(mort)/sum(pop)*pop) %>%
  mutate(SMR = mort/eMort)


# Calculate Smoothed Mortality Rates --------------------------------------

# Create neighborhood object
GA.nb <- poly2nb(GA, queen = F)
W <- nb2mat(GA.nb, style = "B") # In matrix form

# Create distance matrix based on county centroids
d <- as.matrix(dist(gCentroid(GA, byid = T)@coords))

# Get smoothed rates
GA$asMR <- numeric(length(GA))
GA$dsMR <- numeric(length(GA))
GA$asMort <- numeric(length(GA))
GA$dsMort <- numeric(length(GA))

for(i in 1:length(GA)){
  
  # Create vector of current county and its neighbors
  nb.i <- c(i, GA.nb[[i]])
  
  # Adjacency Smoothing (mean of nearby observed mortality rates per 1,000)
  GA$asMR[i] <- mean(GA$mort[nb.i] / GA$pop[nb.i]) * 1000
  GA$asMort[i] <- GA$asMR[i] / 1000 * GA$pop[i]
  
  # Centroid Distance Smoothing
  w.i <- d[nb.i, i]
  w.i[1] <- min(w.i[w.i > 0]) # Weight the region of interest
  GA$dsMR[i] <- sum((1/w.i)*(GA$mort[nb.i] / GA$pop[nb.i]) * 1000) / sum(1/w.i)
  GA$dsMort[i] <- GA$dsMR[i] / 1000 * GA$pop[i]
}

# Bayesian Smoothing with INLA
bmod <- inla(mort ~ 1 + f(OBJECTID, model = "besag", graph = W),
             data = as.data.frame(GA), E = GA$eMort, family = "poisson",
             control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, waic = TRUE))

GA$bsMort <- bmod$summary.fitted.values[, '0.5quant']*GA$eMort
GA$bsMR <- GA$bsMort/GA$pop * 1000

# Bayesian smoothing with CARBayes
# library(CARBayes)
# cbmod <- S.CARbym(formula = mort ~ offset(log(pop)), data = GA, 
#                   family = "poisson", W = W,
#                   burnin = 100000, n.sample = 300000, thin = 100)
# 
# GA$bsMort <- cbmod$fitted.values
# GA$bsMR <- GA$bsMort/GA$pop * 1000

# Export GA shapefile with all new variables
writeOGR(obj = GA, dsn = 'Data/GA', layer = 'counties', 
         driver = 'ESRI Shapefile', overwrite_layer = TRUE)


# Plotting ----------------------------------------------------------------

# Tidy GA shapefile and rejoin variables
GA.tidy <- tidy(GA, region = "NAME10") %>%
  left_join(GA@data, by = c("id" = "NAME10")) %>%
  rename(county = id)

# Plotting with GGPlot2
# GA %>% spplot(zcol = c("MR","Sm.Adj.MR"))
# p1 <- GA.tidy %>%
#   ggplot(aes(x = long, y = lat, group = county)) +
#   geom_polygon(aes(fill = MR), color = 'black') +
#   coord_map("tetra") +
#   labs(title = 'Mortality Rate', 
#        fill = '')
# 
# p2 <- GA.tidy %>%
#   ggplot(aes(x = long, y = lat, group = county)) +
#   geom_polygon(aes(fill = asMR), color = 'black') +
#   coord_map("tetra") +
#   labs(title = 'Adjacency Smoothed Mortality Rate',
#        fill = '')
# 
# p3 <- GA.tidy %>%
#   ggplot(aes(x = long, y = lat, group = county)) +
#   geom_polygon(aes(fill = dsMR), color = 'black') +
#   coord_map("tetra") +
#   labs(title = 'Inverse-Distance Smoothed Mortality Rate',
#        fill = '')
# 
# gridExtra::grid.arrange(p1, p2, p3, ncol = 2)

# Using Faceting (preferred, scales need reworking)
GA.all <- GA.tidy %>%
  pivot_longer(cols = c('MR','asMR','dsMR', 'bsMR'), 
               names_to = 'Method', values_to = 'MR') %>%
  mutate(SMR = case_when(Method == 'MR' ~ mort/eMort,
                         Method == 'asMR' ~ asMort/eMort,
                         Method == 'bsMR' ~ bsMort/eMort,
                         Method == 'dsMR' ~ dsMort/eMort),
         Method = factor(Method, levels = c('MR','asMR','dsMR','bsMR'),
                         labels = c('Observed', 'Adjacency Smoothing',
                         'Inverse-Distance Smoothing', 'Bayesian Smoothing'))) %>%
  mutate(cMR = cut(MR, breaks = seq(0, 0.6, 0.1), include.lowest = TRUE),
         cSMR = cut(SMR, breaks = quantile(SMR, seq(0,1,.1)), include.lowest = TRUE))
  
  
# Plot Mortality Rates
GA.all %>%
  ggplot(aes(x = long, y = lat, group = county)) +
  geom_polygon(aes(fill = cMR), color = 'black') +
  coord_map("tetra") +
  scale_fill_brewer(type = 'div', palette = "PiYG") +
  facet_wrap(~Method, nrow = 2) +
  labs(title = 'Opioid Overdose Mortality Rates in Georgia (2020)',
       fill = 'MR') +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())
ggsave(filename = 'Figures/GA_MR_Lattice_Plot_2020.png')

# Plot Standardized Mortality Ratio (SMR)
GA.all %>%
  ggplot(aes(x = long, y = lat, group = county)) +
  geom_polygon(aes(fill = cSMR), color = 'black') +
  coord_map("tetra") +
  scale_fill_brewer(type = 'div', palette = "PiYG", direction = -1) +
  facet_wrap(~Method, nrow = 2) +
  labs(title = 'Opioid Overdose Standardized Mortality Ratios in Georgia (2020)',
       fill = 'SMR') +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())
ggsave(filename = 'Figures/GA_SMR_Lattice_Plot_2020.png')

# Plotting with Leaflet
library(leaflet)
colors <- colorNumeric(palette = "PiYG", domain = GA$bsMR)
leaflet(data = GA) %>%
  addTiles() %>%
  addPolygons(fillColor = ~ colors(bsMR), color = "", weight = 1,
              fillOpacity = 0.7) %>%
  addLegend(pal = colors, values = GA$bsMR, opacity = 1, title = "bsMR") %>%
  addScaleBar(position = "bottomleft")

# Can also plot simple features
st_as_sf(GA) %>%
  ggplot() +
  geom_sf(aes(fill = SMR))
