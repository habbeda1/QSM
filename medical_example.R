##### Data Example #####

source("qsm.R")

set.seed(24082019)
x1 <- runif(2500, 0, 10)
x2 <- runif(2500, 0, 5)
plot(x1,x2)

y <- ifelse(x1 < 3, 
            "medium pain", 
            ifelse(x2 > 3.5, 
                   "high pain",
                   ifelse(x2 > 1.5, 
                          "medium pain",
                          "no pain")
            )
)

y2 <- ifelse(x1 < 3, 
             5,
             x2*2)

df <- data.frame(x1, x2, health_discrete = as.factor(y))
df2 <- data.frame(x1, x2, health_continuous = y2)

ggplot(df2, aes(x = x1, y = x2, color = health_continuous)) +
  geom_point(size = 3) +
  scale_color_gradient2( low = "yellow", mid = "red", high = "blue", midpoint = 5)


ggplot(df, aes(x = x1, y = x2, group = health_discrete, color = health_discrete)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("coral2", "darkslategray3", "darkolivegreen3"), name = "Health Discrete") +
  ggtitle("Categorical Pain modeled with two features") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        panel.background = element_rect(fill = "grey85")) +
  xlab(expression('x'[1])) +
  ylab(expression('x'[2])) 


#### Data Manipulation ####

library(rpart)

rpart_model <- rpart(health_discrete ~ x1 + x2, data = df)

medical_x1_nCO <- neighborChordsObj(rpart_model,
                                    df,
                                    predictfn = function(object, data, ...){predict(object, data, type = "class", ...)},
                                    features = c("x1"),
                                    method = "all", 
                                    nrandom = 1,
                                    q = list(250/2501))

medical_x1_nCO$migration_matrix <- medical_x1_nCO$migration_matrix[c(1,3,2),]

plot(medical_x1_nCO,
     col = c("coral2", "darkslategray3", "darkolivegreen3"),
     grid.col = c("coral2", "darkslategray3", "darkolivegreen3"),
     transparency = 0.2,
     directional = 1,
     direction.type = c("arrows", "diffHeight"), 
     diffHeight  = -0.04,
     link.arr.type = "big.arrow")

medical_x1_x2_nCO <- neighborChordsObj(rpart_model,
                                       df,
                                       predictfn = function(object, data, ...){predict(object, data, type = "class", ...)},
                                       features = c("x1", "x2"),
                                       method = "all", 
                                       nrandom = 1,
                                       q = list(750/2501, -750/2501))

plot(medical_x1_x2_nCO,
     col = c("coral2", "darkslategray3", "darkolivegreen3"),
     grid.col = c("coral2", "darkslategray3", "darkolivegreen3"),
     transparency = 0.2,
     directional = 1,
     direction.type = c("arrows", "diffHeight"), 
     diffHeight  = -0.04,
     
     link.arr.type = "big.arrow")
