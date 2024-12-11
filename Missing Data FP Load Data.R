library(tidyverse)
library(dplyr)
library(readr)
library(mice)

#### Load Data and Metadata:

# Archetype Data for three different k values
load(file = paste0("Data/ArchetypeFit10_24_k", 3,"rep1.RData"))
archetypes3 = fit
load(file = paste0("Data/ArchetypeFit10_24_k", 6,"rep1.RData"))
archetypes6 = fit
load(file = paste0("Data/ArchetypeFit10_24_k", 12,"rep1.RData"))
archetypes12 = fit

## getting archetype weights
W3=archetypes3$alphas
W3=archetypes3$alphas / rowSums(archetypes3$alphas)
W6=archetypes6$alphas
W6=archetypes6$alphas / rowSums(archetypes6$alphas)
W12=archetypes12$alphas
W12=archetypes12$alphas / rowSums(archetypes12$alphas)

#Pull Metadata (It is the same for each archetype count)
# Type = Level of Pertubation for the Subject: Small, Large, Random, None
# Trial = Subject #, gaps come from excluded subjects
# Duration = Duration of the Stride in seconds.

load("Data/wide10_24.RData")
MetaData <- wide.dat[-which(wide.dat[,708] > 2000),] %>% select(Type, Trial, Duration)

make_mis_MCAR = function(d, 
                         Columns_with_missingness = c(2,3),
                         missingness_pct = .1, 
                         seed = 1) {
  set.seed(seed)
  n = nrow(d)
  mis = rbinom(n, 1, missingness_pct)
  d[mis == 1, Columns_with_missingness] = NA
  d
}

####Generate partial Matrices
#arch = number of initial archetypes
#mis = number of columns with partial missingness
#pct = pct of rows with missingness
arch3mis2pct10 = make_mis_MCAR(d = W3, Columns_with_missingness = c(2,3), missingness_pct = .1)
arch3mis2pct30 = make_mis_MCAR(d = W3, Columns_with_missingness = c(2,3), missingness_pct = .3)
arch3mis2pct90 = make_mis_MCAR(d = W3, Columns_with_missingness = c(2,3), missingness_pct = .9)

arch6mis2pct10 = make_mis_MCAR(d = W6, Columns_with_missingness = c(5,6), missingness_pct = .1)
arch6mis2pct30 = make_mis_MCAR(d = W6, Columns_with_missingness = c(5,6), missingness_pct = .3)
arch6mis2pct90 = make_mis_MCAR(d = W6, Columns_with_missingness = c(5,6), missingness_pct = .9)

arch6mis4pct10 = make_mis_MCAR(d = W6, Columns_with_missingness = c(3,4,5,6), missingness_pct = .1)
arch6mis4pct30 = make_mis_MCAR(d = W6, Columns_with_missingness = c(3,4,5,6), missingness_pct = .3)
arch6mis4pct90 = make_mis_MCAR(d = W6, Columns_with_missingness = c(3,4,5,6), missingness_pct = .9)

arch12mis2pct10 = make_mis_MCAR(d = W12, Columns_with_missingness = c(11,12), missingness_pct = .1)
arch12mis2pct30 = make_mis_MCAR(d = W12, Columns_with_missingness = c(11,12), missingness_pct = .3)
arch12mis2pct90 = make_mis_MCAR(d = W12, Columns_with_missingness = c(11,12), missingness_pct = .9)

arch12mis4pct10 = make_mis_MCAR(d = W12, Columns_with_missingness = c(9,10,11,12), missingness_pct = .1)
arch12mis4pct30 = make_mis_MCAR(d = W12, Columns_with_missingness = c(9,10,11,12), missingness_pct = .3)
arch12mis4pct90 = make_mis_MCAR(d = W12, Columns_with_missingness = c(9,10,11,12), missingness_pct = .9)

arch12mis11pct10 = make_mis_MCAR(d = W12, Columns_with_missingness = c(2,3,4,5,6,7,8,9,10,11,12), missingness_pct = .1)
arch12mis11pct30 = make_mis_MCAR(d = W12, Columns_with_missingness = c(2,3,4,5,6,7,8,9,10,11,12), missingness_pct = .3)
arch12mis11pct90 = make_mis_MCAR(d = W12, Columns_with_missingness = c(2,3,4,5,6,7,8,9,10,11,12), missingness_pct = .9)

#### Store Generated Data
dataStore <- list(MetaData,
arch3mis2pct10, arch3mis2pct30, arch3mis2pct90,
arch6mis2pct10,arch6mis2pct30, arch6mis2pct90,
arch6mis4pct10,arch6mis4pct30, arch6mis4pct90,
arch12mis2pct10,arch12mis2pct30, arch12mis2pct90,
arch12mis4pct10, arch12mis4pct30, arch12mis4pct90,
arch12mis11pct10, arch12mis11pct30, arch12mis11pct90)

save(dataStore, file = "Data/MetaData_and_Missing_Data")

#### Load Saved Data
load("Data/MetaData_and_Missing_Data")

MetaData = dataStore[[1]]
arch3mis2pct10 = dataStore[[2]]
arch3mis2pct30 = dataStore[[3]]
arch3mis2pct90 = dataStore[[4]]

arch6mis2pct10 = dataStore[[5]]
arch6mis2pct30 = dataStore[[6]]
arch6mis2pct90 = dataStore[[7]]

arch6mis4pct10 = dataStore[[8]]
arch6mis4pct30 = dataStore[[9]]
arch6mis4pct90 = dataStore[[10]]

arch12mis2pct10 = dataStore[[11]]
arch12mis2pct30 = dataStore[[12]]
arch12mis2pct90 = dataStore[[13]]

arch12mis4pct10 = dataStore[[14]]
arch12mis4pct30 = dataStore[[15]]
arch12mis4pct90 = dataStore[[16]]

arch12mis11pct10 = dataStore[[17]]
arch12mis11pct30 = dataStore[[18]]
arch12mis11pct90 = dataStore[[19]]
