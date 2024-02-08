##### QSM - Main function #####

require(data.table)
require(magrittr)

neighborChordsObj <- function(model, ## Model to use for the evaluation
                              data, ## Data to use for manipulation
                              predictfn = NULL,  ## predictfunction to get predicted classes [ARGUMENTS : OBJECT, DATA]
                              features = c(), ## Featurenames to be manipulated
                              q = list(), ## Manipulation for the features (as list),
                              method = c("all", "fractional", "randomly"), ## Method for tied observations [FRACTIONAL NOT IMPLEMENTED RIGHT NOW!]
                              data_emp = NULL, ## Original data for estimating the empirical distribution, needs same column names as data
                              nrandom = NULL, ## How many repeats, when method == "randomly"
                              ...){
  
  ## Some input checks
  stopifnot(nrow(data) > 0,
            length(features) == length(q),
            length(features) > 0,
            is.null(data_emp) || (!is.null(data_emp) & all.equal(names(data_emp), names(data))),
            method != "randomly" | !is.null(nrandom)
  )
  
  ## if the method is not repeated random jittering, fix at 1
  if(method != "randomly"){nrandom <- 1}
  
  ## add identifier ID
  data <- data %>% 
    setDT %>% 
    .[, id := 1:nrow(data)]
  
  ## Calculate predictions before manipulation
  if(is.null(predictfn)){
    
    predictions_stored <- predict(model, data, ...)
    
  }else{
    
    predictions_stored <- predictfn(object = model, data = data, ...)
    
  }
  
  predictions_stored <- data.frame(predictions_stored, data[, id])
  
  ## Check for additional empirical data
  if(!is.null(data_emp)){
    
    data_emp <- data_emp[, id := 0]
    
    data <- rbind(data, data_emp)
    
  }
  
  ## Save original data
  data_store <- copy(data)
  result <- NULL
  
  ## do multiple if method = "randomly"
  for(j in 1:nrandom){
    
    ## Manipulating feature by feature
    for(i in 1:length(features)){
      
      ## for numeric feature
      if(data[, eval(as.name(features[i]))] %>% is.numeric){
        
        ## just 1 Argument for q recquired
        stopifnot({q[[i]] %>% length} == 1)
        
        if(method == "all"){
          
          ## Calculate new q-values
          if(q[[i]] >= 0){
            
            ## Estimate ecdf
            ecdffn <- ecdf(data[, eval(as.name(features[i]))])
            
            qtemp <- data[, 
                          .(q = ecdffn(eval(as.name(features[i]))) + q[[i]])] %>% 
              .[q > 1, q := 1]
            
            ## Take new q-values to calculate new feature-values
            data[, eval(features[i]) := quantile(data[, eval(as.name(features[i]))], 
                                                 qtemp %>% unlist,
                                                 type = 1)
            ]
          }else{
            
            ## for negative q -> flip distribution with *(-1)
            
            ## Estimate ecdf
            ecdffn <- ecdf(-data[, eval(as.name(features[i]))])
            
            qtemp <- data[, 
                          .(q = ecdffn(-eval(as.name(features[i]))) - q[[i]])] %>% 
              .[q >= 1, q := 1]
            
            ## Take new q-values to calculate new feature-values
            data[, eval(features[i]) := -quantile(data[, -eval(as.name(features[i]))], 
                                                 qtemp %>% unlist,
                                                 type = 1)
            ]
            
          }
          
          qtemp <- NULL
          
          ## end Method "all"
        }else{
          
          if(method == "fractional"){
            
            ## THIS IS CURRENTLY NOT IN THE PAPER, BUT A TALKED ABOUT OPTION
            ## PRODUCES VALUES, WHICH DID NOT EXIST BEFORE!!!
            
            # ## Estimate ecdf
            # ecdffn <- ecdf(data[, eval(as.name(features[i]))])
            # 
            # if(q[[i]] >= 0){
            #   ## Calculate new q-values
            #   qtemp <- data[, 
            #                 .(q = ecdffn(eval(as.name(features[i]))) + q[[i]])] %>% 
            #     .[q > 1, q := 1] %>% 
            #     .[q < 0, q := 0]
            #   
            #   ## Count number of times a value exists
            #   Nvaltemp <- data[, .(Nval = sum(data[, eval(as.name(features[i]))] == eval(as.name(features[i])))),
            #                    by = row.names(data)]
            #   
            #   ## Take new q-values to calculate new feature-values by using fractional shifting in ties
            #   data[, eval(features[i]) := eval(as.name(features[i])) + (
            #     quantile(data[, eval(as.name(features[i]))], qtemp %>% unlist, type = 1) - eval(as.name(features[i]))
            #   ) / Nvaltemp[, Nval]
            #   ]
            # }else{
            #   
            #   orig <- data[,
            #                .(eval(as.name(features[i])))]
            #   
            #   qtemp <- data[,
            #                 .(q = ecdffn(eval(as.name(features[i]))) + abs(q[[i]]))] %>%
            #     .[q > 1, q := 1] %>%
            #     .[q < 0, q := 0]
            #   
            #   new <- data[, .(quantile(data[, eval(as.name(features[i]))], 
            #                            qtemp %>% unlist,
            #                            type = 1)
            #   )]
            #   
            #   orig_min <- min(orig)
            #   new_min <- min(new)
            #   
            #   ## Look up new value for full shift
            #   xtilde <- vapply(X = orig %>% unlist,
            #                    FUN = function(X){
            #                      
            #                      ## If X is bigger than smallest mapped value, then take the
            #                      ## minimal original value that was mapped to all values
            #                      ## biggerer or equal X;
            #                      ## otherwise just take the minimum
            #                      ifelse(X >= new_min, 
            #                             min(orig[new[, V1] >= X]), 
            #                             orig_min)
            #                      
            #                    },
            #                    FUN.VALUE = 0)
            #   
            #   ## Count number of times a value exists
            #   Nvaltemp <- data[, .(Nval = sum(data[, eval(as.name(features[i]))] == eval(as.name(features[i])))),
            #                    by = row.names(data)]
            #   
            #   data[, eval(features[i]) := eval(as.name(features[i])) + (
            #     xtilde - eval(as.name(features[i]))
            #   ) / Nvaltemp[, Nval]
            #   ]
            #   
            #   xtilde <- NULL
            #   orig <- NULL
            #   new <- NULL
            #   orig_min <- NULL
            #   new_min <- NULL
            #   
            # }
            # 
            # qtemp <- NULL
            # Nvaltemp <- NULL
            # 
            # ## end method "fractional"
            warning("Not implemented right now! Please use one of the other two methods!")
            
          }else{
            
            if(method == "randomly"){
              
              ## Get back original unmanipulated data
              data <- copy(data_store)
              
              ## Count number of times a value exists
              Nvaltemp <- data[, .(Nval = sum(data[, eval(as.name(features[i]))] == eval(as.name(features[i])))),
                               by = row.names(data)]
              v <- var(data[, eval(as.name(features[i]))])
              
              ## Jitter feature slightly and store manipulation
              manip <- Nvaltemp[, original := data[, eval(as.name(features[i]))]] %>% 
                .[Nval > 1, new := original + runif(sum(Nvaltemp[, Nval] > 1), -v*10^-3, +v*10^-3)] %>% 
                .[Nval == 1, new := original]
              
              
              
              
              if(q[[i]] >= 0){
                
                ## Estimate ecdf
                ecdffn <- ecdf(manip[, new])
                
                ## Calculate new q-values
                qtemp <- manip[, 
                               .(q = ecdffn(new) + q[[i]])] %>% 
                  .[q > 1, q := 1]
                
                
                ## Take new q-values to calculate new feature-values by using fractional shifting in ties
                data[, eval(features[i]) := quantile(manip[, new], 
                                                     qtemp %>% unlist,
                                                     type = 1)
                ]
              }else{
                ## for negative q -> flip distribution with *(-1)
                
                ## Estimate ecdf
                ecdffn <- ecdf(-manip[, new])
                
                ## Calculate new q-values
                qtemp <- manip[, 
                               .(q = ecdffn(-new) - q[[i]])] %>% 
                  .[q >= 1, q := 1]
                
                
                ## Take new q-values to calculate new feature-values by using fractional shifting in ties
                data[, eval(features[i]) := -quantile(manip[, -new], 
                                                     qtemp %>% unlist,
                                                     type = 1)
                ]
                
              }
              
              
              ## Transform jittered back to original value
              if(max(Nvaltemp[, Nval]) > 1){
                
                for(m in 1:nrow(manip)){
                  data[manip[m, new] == eval(as.name(features[i])), eval(features[i]) := manip[m, original]]
                }
                
              }
              
              qtemp <- NULL
              Nvaltemp <- NULL
              manip <- NULL
              v <- NULL
            }
            
          }
        } 
        
      }else{
        
        ## If not numeric, then 2 Arguments for q are required
        stopifnot({q[[i]] %>% length} == 2)
        
        data <- data[eval(as.name(features[i])) == q[[i]][[1]], 
                     eval(features[i]) := q[[i]][[2]]]
        
      }
      
    }
    
    ## Drop additional rows from data_emp
    
    data <- data[id != 0]
    
    ## Order data again (for models like kknn)
    
    data <- data[order(id, decreasing = FALSE)]
    
    ## Calculate new predictions
    
    if(is.null(predictfn)){
      
      new_predictions <- predict(model, data, ...)
      
    }else{
      
      new_predictions <- predictfn(object = model, data = data, ...)
      
    }
    
    ## Give resulting list
    if(is.null(result)){
      result <- list(manip_data = data[order(id, decreasing = FALSE)], ## Give changed dataset
                     predictions_before = predictions_stored[,1],
                     predictions_after = new_predictions[order(data[, id], decreasing = FALSE)],
                     migration_matrix = table(predictions_stored[,1], new_predictions[order(data[, id], decreasing = FALSE)]),
                     manipulation = data.frame(features, q = q %>% unlist)
      )
    }else{
      
      result <- list(manip_data = data[order(id, decreasing = FALSE)], ## Give LAST changed dataset
                     predictions_before = c(result$predictions_before, predictions_stored[,1]),
                     predictions_after = c(result$predictions_after, new_predictions[order(data[, id], decreasing = FALSE)]),
                     migration_matrix = result$migration_matrix + table(predictions_stored[,1], new_predictions[order(data[, id], decreasing = FALSE)]),
                     manipulation = data.frame(features, q = q %>% unlist)
      )
      
    }
  }
  
  class(result) <- "nCO"
  
  result
  
}

##### QSM - Plot-Function #####
plot.nCO <- function(object, ## nCO-Objekt
                     ...){
  
  stopifnot(class(object) == "nCO")
  
  ## Melt Migration Matrix
  df <- reshape2::melt(object$migration_matrix)
  
  ## Draw plot
  chordDiagram(df, ...)
  
}