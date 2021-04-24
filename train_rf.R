require(caret)
# require(plotly)
require(magrittr)
# install.packages("doMC")
# library(doMC)

## DATA
params <- read.csv("parameters.csv", header = FALSE) # features
colnames(params) <- c("surface_area", "seq_length", "surface_area/seq_length")
solu_train <- read.csv("data/training/crystal_structs/solubility_values.csv") # solubility values from training data
comb_train <- cbind(params, train) # joined features + solubility

# means vs. sds plot with plotly
# moments <- data.frame(means = colMeans(params),
#                       sds = apply(params, 2, sd))
# moments <- tibble::rownames_to_column(moments, "rowname")
# 
# plot_ly(moments, x = moments$means, y = moments$sds, mode = "markers", text = moments$rowname)

## PLOTS
# Plot features vs solubility (scatter + paired)

features <- 1:3

featurePlot(x = comb_train[, features],
            y = comb_train$solubility,
            plot ="pairs",
            type = c("p", "smooth")
            )


# 10-fold cross validation 
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3)

mtry <- sqrt(ncol(joined[, features]))
tunegrid <- expand.grid(.mtry = mtry)
rf <- train(solubility ~ .,
            data = joined,
            method = "rf",
            num.trees = 500,
            #importance = "permutation",
            #metric = "accuracy",
            tuneGrid = tunegrid,
            trControl = ctrl)




