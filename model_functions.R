library(caret)

fit_knn <- function(train_response_factor, 
                    train_var_mx, 
                    trcontrol) {
  # We set a random seed to ensure reproducibility of results
  set.seed(123)
  # We split the training data in training and validation sets
  # to check for accuracy metrics and evaluate the best k for the model
  knn_fit <- train(V1~ .,
                   method     = "knn",
                   tuneGrid   = expand.grid(k = 1:20),
                   trControl  = trcontrol,
                   metric     = "Accuracy",
                   data       = cbind(V1 = train_response_factor, 
                                      train_var_mx))
  # We extract the best k based on accuracy
  knn_best_k <- knn_fit$bestTune$k
  # We fit a kNN model with the best k in the training set
  knn_fit <- knn3(x=train_var_mx, # training set
                  y=train_response_factor, # training set class labels
                  k=knn_best_k)
  return(knn_fit)
}

confusion_matrix <- function(var_mx, 
                             response_factor, 
                             model){
  # predictions on the training set
  Pred <- 
    tryCatch(
      predict(model,
              var_mx,
              type="raw")
      ,
      error = function(e) 
        predict(model,
                var_mx,
                type="class")
      
      )
  #Probs <- predict(model,
   #                var_mx,
    #               type="prob")
  # compare the predicted labels to real labels
  # get different performance metrics
  conf_mx <-
    confusionMatrix(data = response_factor,
                    reference = Pred)
  # make confusion matrix a dataframe
  cmx <- 
    data.frame(
      Accuracy = conf_mx$overall["Accuracy"],
      CI95 = str_c("(",conf_mx$overall["AccuracyLower"],", ",
                   conf_mx$overall["AccuracyUpper"], ")"),
      No_inf_rate = conf_mx$overall["AccuracyNull"],
      Acc_pval = conf_mx$overall["AccuracyPValue"],
      Sensitivity = conf_mx$byClass["Sensitivity"],
      Specificity = conf_mx$byClass["Specificity"],
      PPV = conf_mx$byClass["Pos Pred Value"],
      NPV = conf_mx$byClass["Neg Pred Value"],
      Precision = conf_mx$byClass["Precision"],
      Recall = conf_mx$byClass["Recall"],
      Prevalence = conf_mx$byClass["Prevalence"],
      Detection_rate = conf_mx$byClass["Detection Rate"],
      Detection_prevalence = conf_mx$byClass["Detection Prevalence"],
      Balanced_Accuracy = conf_mx$byClass["Balanced Accuracy"]
    )
  rownames(cmx) <- NULL
  return(cmx)
}

tt_confusion_matrix <- function(train_var_mx,
                                train_response_factor,
                                test_var_mx,
                                test_response_factor,
                                model){
  # We generate a dataframe extracting the confusion matrices of 
  # the model with the train and test sets
  ttcmx <- 
      bind_rows(
        # We generate the conf matrix of train set
        train = confusion_matrix(var_mx = train_var_mx, 
                                 response_factor = train_response_factor,
                                 model = model) %>% 
          mutate(size = nrow(train_var_mx))
        , 
        # Conf matrix of test set
        test = confusion_matrix(var_mx = test_var_mx, 
                                response_factor = test_response_factor,
                                model = model)%>% 
          mutate(size = nrow(test_var_mx))
        ,
        # We add a column to identify train and test set rows
        .id = "cohort"
      )
}

roc_curve <- function(response_mx,
                      response_factor,
                      model) {
  
  Probs <- predict(model,
                   response_mx,
                   type="prob")
  
  roccurve <- pROC::roc(response = response_factor,
                                  predictor = Probs[,1],
                                  ## This function assumes that the second
                                  ## class is the class of interest, so we
                                  ## reverse the labels.
                                  levels = rev(levels(response_factor)))
  # We need to make a datafarme with the model data to plot it together with
  # the rest of models for comparison
  rocdf <- 
    data.frame(
      Sensitivity = roccurve$sensitivities,
      OmSpecificity = 1 - roccurve$specificities,
      Thresholds = roccurve$thresholds,
      AUC = rep(auc(roccurve), length(roccurve$thresholds))
    ) %>%
    mutate(order = factor(seq(nrow(.),1,-1), levels=seq(nrow(.),1,-1))) %>%
    arrange(desc(order))
  
  return(rocdf)
  
}

tt_roc_curve <- function(train_var_mx,
                         train_response_factor,
                         test_var_mx,
                         test_response_factor,
                         model,
                         name){
  
  # We generate a dataframe extracting the confusion matrices of 
  # the model with the train and test sets
  ttcmx <- 
    bind_rows(
      # We generate the conf matrix of train set
      train = roc_curve(response_mx = train_var_mx,
                        response_factor = train_response_factor,
                        model = model) %>% 
        mutate(size = nrow(train_var_mx),
               model = name)
      , 
      # Conf matrix of test set
      test = roc_curve(response_mx = test_var_mx,
                       response_factor = test_response_factor,
                       model = model)%>% 
        mutate(size = nrow(test_var_mx),
               model = name)
      ,
      # We add a column to identify train and test set rows
      .id = "cohort"
    ) %>%
    mutate(test = str_c(model, ": ", cohort, " (AUC=",round(AUC,3),")"))
}

fit_random_forest <- function(train_response_factor, 
                              train_var_mx, 
                              trcontrol) {
  
  # We set a random seed to ensure reproducibility of results
  set.seed(123)
  # We split the training data in training and validation sets
  # to check for accuracy metrics and evaluate the best k for the model
  rf_m <- train(V1~., 
                data        = cbind(V1 = train_response_factor, 
                                    train_var_mx), 
                method      = "ranger",
                trControl   = trcontrol,
                metric      = "Accuracy",
                importance  = "permutation" #, # calculate importance
                #tuneGrid    = data.frame( mtry=mtry,
                 #                         min.node.size = min.node.size,
                  #                        splitrule="gini")
                )
  return(rf_m)
}


fit_log_reg <- function(train_response_factor, 
                        train_var_mx, 
                        trcontrol) {
  
  # We split the training data in training and validation sets
  # to check for accuracy metrics and evaluate the best k for the model
  logregFit <- train(V1 ~ ., 
                     data = cbind(V1 = train_response_factor, 
                                  train_var_mx), 
                     method = "glm",
                     family = "binomial",
                     trControl=    trcontrol,
                     metric      = "Accuracy"
  )
  return(logregFit)
}


fit_elastic_net <- function(train_response_factor, 
                            train_var_mx, 
                            trcontrol) {
  
  # We split the training data in training and validation sets
  # to check for accuracy metrics and evaluate the best k for the model
  elanetFit <- train(V1 ~ ., 
                     data = cbind(V1 = train_response_factor, 
                                  train_var_mx), 
                     method = "glmnet",
                     trControl=    trcontrol,
                     metric      = "Accuracy"
  )
  return(elanetFit)
}

save(fit_knn, confusion_matrix, tt_confusion_matrix, roc_curve,
     tt_roc_curve, fit_random_forest, fit_log_reg, fit_elastic_net,
     file = "./model_functions.RData")
