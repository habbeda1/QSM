##### Yeast data - QSM Example #####
#### Read Data ####
library(readr)
library(magrittr)
library(ggplot2)
library(circlize)
library(nnet)

## Download data from https://www.openml.org/search?type=data&sort=runs&id=181&status=active

yeast <- read_csv("dataset_185_yeast.csv")

class(yeast) <- "data.frame"

yeast$class_protein_localization <- as.factor(yeast$class_protein_localization)

#### Import QSM ####
source("qsm.R")

#### Models ####
set.seed(11082021)
ml_multinom <- multinom(class_protein_localization ~ ., yeast)

#### QSM ####

##
## multinom - mit
##
set.seed(31082021)
nCO_randomly <- neighborChordsObj(ml_multinom,
                                  yeast,
                                  features = c("mit"),
                                  method = "randomly", 
                                  nrandom = 10,
                                  q = list(75/(nrow(yeast)+1))
)
nCO_all <- neighborChordsObj(ml_multinom,
                             yeast,
                             features = c("mit"),
                             method = "all", 
                             nrandom = 1,
                             q = list(75/(nrow(yeast)+1))
)
nCO_all$migration_matrix
diag(nCO_all$migration_matrix) <- 0

nCO_randomly$migration_matrix

#### Plotting ####
plot(nCO_all,
     directional = 1,
     direction.type = c("arrows", "diffHeight"), 
     diffHeight  = -0.04,
     link.arr.type = "big.arrow",
     col = c("darkolivegreen3", "darkslategray3", "coral2", "burlywood3", "tomato3", "darkorange2", "deeppink3", "deepskyblue3", "bisque4", "chartreuse4"),
     grid.col = c("darkolivegreen3", "darkslategray3", "coral2", "burlywood3", "tomato3", "darkorange2", "deeppink3", "deepskyblue3", "bisque4", "chartreuse4"),
     annotationTrack = c("grid"#, "axis"
     ))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), paste0(si), sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.7, adj = c(0.5,0.4))
  circos.xaxis(sector.index = si, major.at = seq(0,55,2), labels.cex = 0.7)
}

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), paste0(si), sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.7, adj = c(0.5,0.4))
  circos.xaxis(sector.index = si, labels.cex = 0.7)
}

plot(nCO_randomly,
     directional = 1,
     direction.type = c("arrows", "diffHeight"), 
     diffHeight  = -0.04,
     link.arr.type = "big.arrow",
     col = c("darkolivegreen3", "darkslategray3", "coral2", "burlywood3", "tomato3", "darkorange2", "deeppink3", "deepskyblue3", "bisque4", "chartreuse4"),
     grid.col = c("darkolivegreen3", "darkslategray3", "coral2", "burlywood3", "tomato3", "darkorange2", "deeppink3", "deepskyblue3", "bisque4", "chartreuse4"),
     annotationTrack = c("grid", "axis"))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), paste0(si), sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.7, adj = c(0.5,0.4))
}


#### Tie-Method comparison ####
test_all <- vapply(1:nrow(yeast), 
                   FUN = function(X){
                     
                     temp <- neighborChordsObj(ml_multinom,
                                               yeast,
                                               #                                               predictfn = function(object, data, ...){predict(object, data, type = "response", ...)$prediction},
                                               features = c("mit"),
                                               method = "all", 
                                               nrandom = 1,
                                               q = list(X/(nrow(yeast)+1))
                     )
                     
                     nrow(yeast)-temp$migration_matrix %>% diag %>% sum
                   },
                   FUN.VALUE = 0
)

plot(1:nrow(yeast),test_all,type = "l")


test_randomly <- vector("numeric", nrow(yeast))
set.seed(11082021)
for(i in 1:nrow(yeast)){
  temp <- neighborChordsObj(ml_multinom,
                            yeast,
                            features = c("mit"),
                            method = "randomly",
                            nrandom = 10,
                            q = list(i/(nrow(yeast)+1))
  )
  
  test_randomly[i] <- nrow(yeast)*10-temp$migration_matrix %>% diag %>% sum
  print(i)
}
plot(1:nrow(yeast),test_randomly,type = "l", col = "blue", lty = 2, ylab = "Changed Predictions", xlab = "v")
lines(test_all*10, col = "red", lty = 2)

## ggplot
df <- data.frame(value = c(0, test_all, 0, test_randomly/10), method = c(rep(c("all", "randomly"), each=nrow(yeast)+1)), q = (0:nrow(yeast))/(nrow(yeast)+1))

ggplot(df[df$q < 0.15,], aes(x = q, y = value, color = method)) + 
  geom_line(size = 0.8) +
  theme_minimal() +
  xlab(expression('q'[mit])) + 
  ylab("Switched Predictions")

#### Found neighborhoods with different q ####

test_all <- lapply(1:nrow(yeast), 
                   FUN = function(X){
                     
                     temp <- neighborChordsObj(ml_multinom,
                                               yeast,
                                               features = c("mit"),
                                               method = "all", 
                                               nrandom = 1,
                                               q = list(X/(nrow(yeast)+1))
                     )
                     
                     temp$migration_matrix %>% reshape2:::melt()
                   }
)

dt_final <- data.table()
for(i in 1:length(test_all)){
  
  dt_final <- rbind(dt_final, 
                    data.table(test_all[[i]], q = i/(nrow(yeast)+1)))
  
}

dt_final[, direction := paste0(Var1, " -> ", Var2)]

found <- dt_final[, .(sum(value)), .(direction)]

dt_final_flt <- dt_final[direction %in% found[V1 > 0, direction] & !Var1 == Var2]

dt_final_found <- dt_final_flt[Var1 == "ME3", .(found = sum(value != 0)), .(q)]
dt_final_found <- dt_final_flt[, .(found = sum(value != 0)), .(q)]

library(ggplot2)
ggplot(dt_final_flt[Var1 == "ME3"], aes(x = q, y = value, group = direction, color = direction)) +
  geom_line(linewidth = 2) +
  ylim(c(0,5))

ggplot(dt_final_found, aes(x = q, y = found)) +
  geom_line() +
  ylim(c(0,21)) +
  ylab("Found Neighbors") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "grey85"))

#### Boxplot with the mit values per class ####

library(readr)
library(magrittr)
library(ggplot2)
library(circlize)
library(nnet)
library(ranger)
library(data.table)

yeast <- read_csv("Methoden_Paper/data/dataset_185_yeast.csv")

class(yeast) <- "data.frame"

yeast$class_protein_localization <- as.factor(yeast$class_protein_localization)

yeast <- as.data.table(yeast)

ggplot(yeast, aes(x = class_protein_localization, y = mit)) +
  geom_boxplot(fill = "#FFD900") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "grey85")) +
  ylab("mit") +
  xlab("Target Class")
