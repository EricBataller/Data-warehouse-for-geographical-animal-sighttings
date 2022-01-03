library(sf)
library(sp)
library(raster)
library(rgbif)
library(protolite)
library(tidyverse)
library(bayestestR)
library(rstanarm)
library(mombf)


## We add the bioclim names to make it easier to follow
clim_names <- c("Annual_Mean_Temperature",
                "Mean_Diurnal_Range",
                "Isothermality",
                "Temperature_Seasonality",
                "Max_Temperature_of_Warmest_Month",
                "Min_Temperature_of_Coldest_Month",
                "Temperature_Annual_Range",
                "Mean_Temperature_of_Wettest_Quarter",
                "Mean_Temperature_of_Driest_Quarter",
                "Mean_Temperature_of_Warmest_Quarter",
                "Mean_Temperature_of_Coldest_Quarter",
                "Annual_Precipitation",
                "Precipitation_of_Wettest_Month",
                "Precipitation_of_Driest_Month",
                "Precipitation_Seasonality",
                "Precipitation_of_Wettest_Quarter",
                "Precipitation_of_Driest_Quarter",
                "Precipitation_of_Warmest_Quarter",
                "Precipitation_of_Coldest_Quarter")

## Load the climate data  and elevation from worldclim
clim_data <- raster::getData("worldclim",var="bio",res=10)
names(clim_data) <- clim_names

altitude <- raster::getData("worldclim", var="alt", res= 10)

##Chose the fileds we want 
fields <- c("key",
            "decimalLatitude",
            "decimalLongitude",
            "family",
            "species")

## Data set key corresponds to IFN3
quercus <- occ_search(country = "ES",
                      datasetKey = "fab4c599-802a-4bfc-8a59-fc7515001bfa",
                      familyKey = 4689,
                      fields = fields,
                      limit  = 80000)
quercus_data <- quercus$data

##Data set key corresponds to Jardin botánico, limit added because error was
## produced for bigger samples
pastos <- occ_search(country = "ES",
                     datasetKey = "4cf3eec1-b902-40c9-b15b-05c5fe5928b6",
                     familyKey = 3073,
                     fields = fields,
                     limit = 100000)

pastos_data <- pastos$data


##make the points spatial/raterize them
quercus_points <- SpatialPoints(quercus_data[,4:5], proj4string = clim_data@crs)
pastos_points <- SpatialPoints(pastos_data[,4:5], proj4string = clim_data@crs)

## get the info from the raster layer of clim_data
quercus_clim <- raster::extract(clim_data, quercus_points)
altitude_quercus <- raster::extract(altitude, quercus_points)



#Get a quercus database with climate variables and unique rows
quercus_alt <- cbind.data.frame(coordinates(quercus_points), altitude_quercus)
quercus_final <- cbind.data.frame(coordinates(quercus_points), quercus_clim)  %>% 
 inner_join(quercus_data) %>%
  distinct(key, .keep_all = TRUE) %>%
  inner_join(quercus_alt) %>%
  distinct(key, .keep_all = TRUE)

names(quercus_final)[names(quercus_final) == 'altitude_quercus'] <- 'elevation'

#Get the same info but for grasses (poaceae)
pastos_clim <- raster::extract(clim_data, pastos_points)
altitude_pastos <- raster::extract(altitude, pastos_points)

pastos_alt <- cbind.data.frame(coordinates(pastos_points), altitude_pastos)
pastos_final <- cbind.data.frame(coordinates(pastos_points), pastos_clim) %>% 
  inner_join(pastos_data) %>%
  distinct(key, .keep_all = TRUE) %>%
  inner_join(pastos_alt) %>%
  distinct(key, .keep_all = TRUE)
names(pastos_final)[names(pastos_final) == 'altitude_pastos'] <- 'elevation'

#Combine the two db's

data_full <- bind_rows(quercus_final, pastos_final)

#Create the sections of Spain, assuming that each 13x13 section is unique
unique_bioclim <- distinct(data_full,
                           Annual_Mean_Temperature,
                           Mean_Diurnal_Range,
                           Isothermality,
                           Temperature_Seasonality,
                           Max_Temperature_of_Warmest_Month,
                           Min_Temperature_of_Coldest_Month,
                           Temperature_Annual_Range,
                           Mean_Temperature_of_Wettest_Quarter,
                           Mean_Temperature_of_Driest_Quarter,
                           Mean_Temperature_of_Warmest_Quarter,
                           Mean_Temperature_of_Coldest_Quarter,
                           Annual_Precipitation,
                           Precipitation_of_Wettest_Month,
                           Precipitation_of_Driest_Month,
                           Precipitation_Seasonality,
                           Precipitation_of_Warmest_Quarter,
                           Precipitation_of_Wettest_Quarter,
                           Precipitation_of_Driest_Quarter,
                           Precipitation_of_Coldest_Quarter,
                           elevation)


group_number <- 1:nrow(unique_bioclim)
unique_bioclim_group <- bind_cols(unique_bioclim, group_number)
names(unique_bioclim_group)[names(unique_bioclim_group) == '...21'] <- 'group'


#Group our data by these section and get a count just for preliminatory analysis
data_full_group <- left_join(data_full, unique_bioclim_group) %>% 
  group_by(group, species) %>%
  add_count(name = "Count") %>%
  distinct(species, .keep_all = TRUE) %>%
  select(-c(key, decimalLatitude, decimalLongitude))


## Created a presence list for each species
data_list <- list()

for (i in group_number) {
  test <- data_full_group %>%
    ungroup() %>%
    filter(group == i) %>%
    select(c(group, species))
  
  dat <- right_join(test, unique_tree_species) %>%
    mutate(new_group = i) %>%
    mutate(is_present = case_when(group >= 1 ~ 1,
                                  TRUE ~ 0)) %>%
    select(-c(group))
  names(dat)[names(dat) == 'new_group'] <- 'group'
  data_list[[i]] <- dat # add it to your list
}

##Get the species for the final list
get_species_name <- function(df) {
  return(df[1,1])
}

species_name_vector <- map(final_list, get_species_name) %>%
  unlist()

## Forward fill na's for the bioclim variables
data_precence_group <- bind_rows(data_list) %>%
  full_join(data_full_group) %>%
  fill(Annual_Mean_Temperature,
       Mean_Diurnal_Range,
       Isothermality,
       Temperature_Seasonality,
       Max_Temperature_of_Warmest_Month,
       Min_Temperature_of_Coldest_Month,
       Temperature_Annual_Range,
       Mean_Temperature_of_Wettest_Quarter,
       Mean_Temperature_of_Driest_Quarter,
       Mean_Temperature_of_Warmest_Quarter,
       Mean_Temperature_of_Coldest_Quarter,
       Annual_Precipitation,
       Precipitation_of_Wettest_Month,
       Precipitation_of_Driest_Month,
       Precipitation_Seasonality,
       Precipitation_of_Warmest_Quarter,
       Precipitation_of_Wettest_Quarter,
       Precipitation_of_Driest_Quarter,
       Precipitation_of_Coldest_Quarter,
       elevation) %>%
  select(-c(Count, group))

#Get the full list without na's for the bioclim variables

final_list <- data_precence_group %>%
  group_by(species) %>%
  filter(sum(is_present) > 399)%>%
  group_split(.keep = TRUE)

#Function to calculate the model selection
matrix_transform <- function(df){
  X <- df %>%
    select(-c(species, is_present, family, Isothermality, Temperature_Annual_Range)) %>%
    as.matrix()
  y <- df$is_present
  pc= bicprior()  
  pm= modelbinomprior()  
  
  ms = modelSelection(y, X, priorCoef = pc, priorDelta = pm, family = 'binomial', method='ALA')
  
  return(ms)
}

#Matrix of model selection
final_matrix <-map(final_list, matrix_transform)

#Calculate the betas from the matrix given by the model selection
beta_extraction <- function(mat){
  beta <- coefByModel(mat, length(postProb(mat)[,3]))
  beta_mat <- beta$postmean[, -1]*(postProb(mat)[,3])
  final_beta <- colSums(beta_mat)
  return(final_beta)
}

#Calculate the lower bound
lower_bound <- function(mat) {
  beta <- coefByModel(mat, length(postProb(mat)[,3]))
  beta_mat <- beta$ci.low[, -1]*(postProb(mat)[,3])
  final_beta <- colSums(beta_mat)
  return(final_beta)
}

##Claculate the upper bound
upper_bound <- function(mat) {
  beta <- coefByModel(mat, length(postProb(mat)[,3]))
  beta_mat <- beta$ci.up[, -1]*(postProb(mat)[,3])
  final_beta <- colSums(beta_mat)
  return(final_beta)
}

## Create the lists
beta_list <- map(final_matrix, beta_extraction)
names(beta_list) <- species_name_vector
ci_low_list <- map(final_matrix, lower_bound)
names(ci_low_list) <- species_name_vector
ci_up_list <- map(final_matrix, upper_bound)
names(ci_up_list) <- species_name_vector

#Turn them into long df's for ggplot and joining them
beta_df <- bind_rows(beta_list)

ci_low <- bind_rows(ci_low_list) %>%
  gather("Variables", "Lower") %>%
  bind_cols(rep(species_name_vector, 18))
names(ci_low) <- c("Variables", "Lower", "Species")

ci_up <- bind_rows(ci_up_list) %>%
  gather("Variables", "Upper") %>%
  bind_cols(rep(species_name_vector, 18))
names(ci_low) <- c("Variables", "Upper", "Species")

new_beta <- beta_df %>%
  gather("Variables", "Values") %>%
  bind_cols(rep(species_name_vector, 18))
names(new_beta) <- c("Variables", "Values", "Species")

final_beta <- left_join(new_beta, ci_low) %>%
  bind_cols(ci_up$Upper)
names(final_beta) <- c("Variables", "Beta", "Species", "Lower", "Upper")


##Facet_warp for 4 species at a time
ggplot(final_beta %>%
         filter(Species %in% species_name_vector[1:4]), aes(x = Variables, y = Beta)) +
  geom_point(aes(col = Species)) +
  geom_errorbar(aes(ymin= Lower, ymax=Upper, col = Species), width=.2,
                position=position_dodge(0.05)) +
  facet_wrap(~ Species) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=5),
        axis.text.x = element_blank())
ggsave("bma_grass.png", dpi = 500)

ggplot(final_beta %>%
         filter(Species %in% species_name_vector[5:8]), aes(x = Variables, y = Beta)) +
  geom_point(aes(col = Species)) +
  geom_errorbar(aes(ymin= Lower, ymax=Upper, col = Species), width=.2,
                position=position_dodge(0.05)) +
  facet_wrap(~ Species) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=5),
        axis.text.x = element_blank())
ggsave("bma_grass_1.png", dpi = 500)


ggplot(final_beta %>%
         filter(Species %in% species_name_vector[9:12]), aes(x = Variables, y = Beta)) +
  geom_point(aes(col = Species)) +
  geom_errorbar(aes(ymin= Lower, ymax=Upper, col = Species), width=.2,
                position=position_dodge(0.05)) +
  facet_wrap(~ Species) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=8), 
        legend.text=element_text(size=5),
        axis.text.x = element_blank())

ggsave("bma_tree.png", dpi = 500)



#Creation of Model size plots

for (i in 1:length(species_name_vector)) {
  nvars= rowSums(final_matrix[[i]]$postSample)
  par(mar=c(5,4,1,1), cex.lab=1.3, cex.axis=1.3)
  plot(nvars, type='l', xlab='Gibbs iteration', ylab='Model size', main = species_name_vector[i])
  
}

##creation of trace plot for log-posterior model probabilities
for (i in 1:length(species_name_vector)) {
  par(mar=c(4,5,1,1), cex.lab=0.75, cex.axis=1)
  plot(final_matrix[[i]]$postProb,
       type='l',
       xlab='Gibbs iteration',
       ylab='log p(y | gamma) + log p(gamma)',
       main= species_name_vector[i])
}

