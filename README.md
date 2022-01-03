# Geospatial analysis and Data warehouse creation for geographical plants/animals sighttings

This github repository contains two projects with very different goals. The only reason they were put together is because they work with the same type of data and have most of the data sources in common.

Part 1 - Geospatial analysis and Bayesian Model Selection.

This is perhaps the most interesting part of the repository. Our goal is to understand more about the species of plants present in a Spanish Dehesa (mainly Oaks and Pasture species). We perform Latent Dirichlet Allocation (LDA) on occurrances data in order to understand the biogeographical patterns of the different groups of species across Spain's municipalities. Lastly, we perform Bayesian Model Selection (BMS) to further understand the true predictors of such occurrences.

Part 2 - Creation of a Data warehouse for animal sighttings.

In this part of the repository we go over the whole process of creation of a Data warehouse (DW). Our objective is to build a minimum viable product capable of storing not only occurrance data of each species, but also their main carachteristics along with the climatic data of the regions they have been observed in. This effectively unifies all the important information required for most of the analyis on animal sighttings. Only geographical occurrences and biodiversity data on carnivore mammals is used, but the DW can be easily extended. 
Note: This project was done using PostgreSQL and Pentaho software.
