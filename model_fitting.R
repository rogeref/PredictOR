#### Packages ####
if (!require('tidyverse')) install.packages('tidyverse')
library(tidyverse)
if (!require('caret')) install.packages('caret')
library(caret)
if (!require('cowplot')) install.packages('cowplot')
library(cowplot)
if (!require('corrplot')) install.packages('corrplot')
library(corrplot)
if (!require('pROC')) install.packages('pROC')
library(pROC)
if (!require('devtools')) install.packages('devtools')
library(devtools)
if (!require('ggbiplot')) install_github("vqv/ggbiplot")
library(ggbiplot)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)
if (!require('ComplexHeatmap')) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
if (!require('ggrepel')) install.packages('ggrepel')
library(ggrepel)
if (!require('arules')) install.packages('arules')
library(arules)
if (!require('RColorBrewer')) install.packages('RColorBrewer')
library(RColorBrewer)
if (!require('ggridges')) install.packages('ggridges')
library(ggridges)
load("./model_functions.RData")

#### Data import - Training ####

print("Importing data...")

# We import the labels to predict (response)
train_response <- 
  read.delim("./in/train_response.txt") %>%
  filter(!is.na(response))
# We import the expression matrix of the genes of interest and
# filter out of the expression matrix the samples that have no
# response data available
train_gexp <- 
  read.delim("./in/train_gene_matrix.txt") %>%
  as.data.frame() %>%
  (function(df){ names(df)[1] <- "sample" ; return(df) }) %>%
  inner_join(train_response, by="sample") %>%
  select(-response) %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

#### Data import - Testing ####

# We import the labels to predict (response)
test_response <- 
  read.delim("./in/test_response.txt") %>%
  filter(!is.na(response))
# We import the expression matrix of the genes of interest and
# filter out of the expression matrix the samples that have no
# response data available
test_gexp <- 
  read.delim("./in/test_gene_matrix.txt") %>%
  as.data.frame() %>%
  (function(df){ names(df)[1] <- "sample" ; return(df) }) %>%
  inner_join(test_response, by="sample") %>%
  select(-response) %>%
  column_to_rownames(var="sample") %>%
  as.matrix()


#### Data transformations - Training AND Testing ####

print("Data preprocessing: log-tansforming, filtering and scaling...")

# If the range of values is in the positive numbers (>0),
# we smoothen the extreme expression values by taking a log10
if (min(train_gexp)<0){
  # Expression values remain the same without log transform (although name changes)
  train_log_gexp <- train_gexp
} else {
  train_log_gexp <- log10(train_gexp + 1)
}
if (min(test_gexp)<0){
  # Expression values remain the same without log transform (although name changes)
  test_log_gexp <- test_gexp
} else {
  test_log_gexp <- log10(test_gexp + 1)
}

#### Data filtering and scaling - Training ####

# remove near zero variation for the columns at least
# 85% of the values are the same
# this function creates the filter but doesn't apply it yet
train_nzv = preProcess(train_log_gexp, method="nzv", uniqueCut = 15)

# apply the filter using "predict" function
# return the filtered dataset and assign it to nzv_tgexp
# variable
train_nzv_log_gexp = predict(train_nzv, train_log_gexp)

# center the data so that predictor variables will have zero means
train_processCenter = preProcess(train_nzv_log_gexp, method = c("center"))
train_processScale = preProcess(train_nzv_log_gexp, method = c("center", "scale"))
train_cent_nzv_log_gexp = predict(train_processCenter, train_nzv_log_gexp)
train_sc_nzv_log_gexp = predict(train_processScale, train_nzv_log_gexp)

# remove genes with NA values
if (anyNA(train_sc_nzv_log_gexp)){
  train_sc_nzv_log_gexp = train_sc_nzv_log_gexp[ , colSums(is.na(train_sc_nzv_log_gexp)) == 0]
}

# filter out the variables that are highly correlated
# create a filter for removing highly correlated variables
# if two variables are highly correlated only one of them
# is removed
# The absolute values of pair-wise correlations are considered. 
# If two variables have a high correlation, the function looks at 
# the mean absolute correlation of each variable and removes the 
# variable with the largest mean absolute correlation.
train_corrFilt=preProcess(train_sc_nzv_log_gexp, method = "corr", cutoff = 0.9)
train_filt_gexp = predict(train_corrFilt, train_sc_nzv_log_gexp)

#### Data filtering and scaling - Testing ####

# Get the same filtered matrix of genes than the training set
test_log_gexp = test_log_gexp[,colnames(train_filt_gexp)]
test_processCenter = preProcess(test_log_gexp, method = c("center"))
test_processScale = preProcess(test_log_gexp, method = c("center", "scale"))
test_cent_nzv_log_gexp = predict(test_processCenter, test_log_gexp)
test_sc_nzv_log_gexp = predict(test_processScale, test_log_gexp)
# remove genes with NA values
if (anyNA(test_sc_nzv_log_gexp)){
  print("THERE ARE MISSING VALUES IN THE TESTING SET!!!")
}
test_filt_gexp = test_sc_nzv_log_gexp


#### Data representation for sanity check ####

print("Representing the data transformations as plots...")

# We plot how data transformations impact gene expression distribution per sample
pdf("./out/data_transformations.pdf", width = 24, height=18)
cowplot::plot_grid(ncol=2, nrow=4, align = "hv",
                   t(train_gexp) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(train_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="Expression", 
                          title="Expr. by sample")
                   ,
                   log10(t(train_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(train_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10")
                   ,
                   log10(t(train_cent_nzv_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(train_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10 + CENTERING")
                   ,
                   log10(t(train_sc_nzv_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(train_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10 + SCALING")
                   
                   ,
                   
                   t(test_gexp) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(test_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="Expression", 
                          title="Expr. by sample")
                   ,
                   log10(t(test_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(test_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10")
                   ,
                   log10(t(test_cent_nzv_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(test_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10 + CENTERING")
                   ,
                   log10(t(test_sc_nzv_log_gexp)+1) %>%
                     as.data.frame() %>%
                     rownames_to_column(var="gene") %>%
                     gather(sample, expr, -gene) %>%
                     inner_join(test_response, by="sample") %>%
                     ggplot(aes(x=sample, y=expr, fill=response)) +
                     geom_boxplot() +
                     theme_classic() +
                     theme(axis.text.x = element_text(angle=45, hjust=1)) +
                     labs(y="log10(expression)", 
                          title="Expr. by sample after LOG10 + SCALING")
)
dev.off()

# Plot correlation matrices
pdf("./out/correlations.pdf", 
    width = if (ncol(train_sc_nzv_log_gexp) > 100) {16
      } else if (ncol(train_sc_nzv_log_gexp) < 10) {8
      } else {12}, 
    height= if (ncol(train_sc_nzv_log_gexp) > 100) {16
    } else if (ncol(train_sc_nzv_log_gexp) < 10) {8
    } else {12}
    )

train_M = cor(train_sc_nzv_log_gexp)
corrplot(train_M)

test_M = cor(test_sc_nzv_log_gexp)
corrplot(test_M)

cowplot::plot_grid( align="hv", ncol=1, nrow=2,
                    ggplot(data=NULL, aes(x=as.numeric(train_M))) + 
                      geom_density() + 
                      theme_classic() +
                      labs(title="Regression coefs. in Training")
                    ,
                    ggplot(data=NULL, aes(x=as.numeric(test_M))) + 
                      geom_density() + 
                      theme_classic() +
                      labs(title="Regression coefs. in Test")
)

dev.off()


#### Get the objects for model testing ####

print("Getting the objects for model testing...")

train_response_mx <-
  train_response %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

train_response_expr <- 
  merge(train_response_mx, train_filt_gexp, by="row.names") %>%
  column_to_rownames(var="Row.names")

train_var_mx <- train_response_expr[,-1]

train_response_factor <- factor(train_response_expr[,1], levels=c("OR","NR"))

test_response_mx <-
  test_response %>%
  column_to_rownames(var="sample") %>%
  as.matrix()

test_response_expr <- 
  merge(test_response_mx, test_filt_gexp, by="row.names") %>%
  column_to_rownames(var="Row.names")

test_var_mx <- test_response_expr[,-1]

test_response_factor <- factor(test_response_expr[,1], levels=c("OR","NR"))

rm(train_gexp, train_log_gexp, train_M, train_nzv, train_corrFilt,
   train_nzv_log_gexp, train_processCenter, train_processScale, train_response,
   train_cent_nzv_log_gexp, train_sc_nzv_log_gexp,
   train_response_mx,
   test_gexp, test_log_gexp, test_M,test_processCenter, test_processScale, 
   test_response, test_cent_nzv_log_gexp, test_sc_nzv_log_gexp,
   test_response_mx)

#### PCA ####

print("Principal Component Analysis...")

train_pca <- prcomp(train_var_mx)
test_pca <- prcomp(test_var_mx)

pdf("./out/pca.pdf", width = 15, height=8)

cowplot::plot_grid( align="hv",
                    
                    ggbiplot(train_pca,
                             groups = train_response_factor,
                             ellipse=TRUE) +
                      theme_bw() +
                      labs(title="PCA with training cohort")
                    ,
                    ggbiplot(test_pca,
                             groups = test_response_factor,
                             ellipse=TRUE) +
                      theme_bw() +
                      labs(title="PCA with test cohort")
)

dev.off()

#### Ridge plot #####

print("Ridge Plot...")

pdf("./out/ridges.pdf", 
    width = if (ncol(train_response_expr) > 100) {20
    } else if (ncol(train_response_expr) < 10) {10
    } else {15}, 
    height= if (ncol(train_response_expr) > 100) {20
    } else if (ncol(train_response_expr) < 10) {7
    } else {15}
    )

cowplot::plot_grid( ncol=2, nrow=1, align = "hv",
                    train_response_expr %>% 
                      rownames_to_column(var="sample") %>%
                      gather(gene, expr, -sample, -response) %>%
                      ggplot(aes(x=expr, y=gene, fill=response)) +
                      geom_density_ridges(alpha=0.6) +
                      theme_classic() +
                      labs(title="Gene expression in Training cohort")
                    ,
                    test_response_expr %>% 
                      rownames_to_column(var="sample") %>%
                      gather(gene, expr, -sample, -response) %>%
                      ggplot(aes(x=expr, y=gene, fill=response)) +
                      geom_density_ridges(alpha=0.6) +
                      theme_classic()+
                      labs(title="Gene expression in Test cohort")
)

dev.off()


#### Methods of train-control set separation ####
trctrl <-
  list(
    LOOCV = trainControl(method  = "LOOCV")
    ,
    CV10 = trainControl(method = "cv",
                        number=10)
    ,
    Bootstrap = trainControl( method = "boot",
                              number=100,
                              returnResamp="all")
  )

#### k-nearest neighbors prediction####

print("Fitting a k-NN models...")

# Fit a knn model to predict responders based on each trainControl method
knn_models <- 
  lapply(trctrl,
         function(x){
           fit_knn(train_response_factor,
                   train_var_mx, 
                   trcontrol = x)
         }) %>%
  (function(lt){
    lt_names <- str_c( map_chr(lt, 
                               function(md){return(md$k)}), 
                       "-nn_", 
                       names(lt))
    names(lt) <- lt_names
    return(lt)
  })


#### Random Forest prediction####

print("Fitting Random Forest models...")

# Fit random forest models
rfo_models <- lapply(trctrl[3],
                     function(x){
                       fit_random_forest( train_response_factor,
                                          train_var_mx, 
                                          trcontrol = x)
                     }) %>%
  (function(lt){
    names(lt) <- str_c("RandomForest_", names(lt))
    return(lt)
  })

# print Out Of the Box error
#   print(paste("Random forest prediction error:",
#   rfo_models$RandomForest_Bootstrap$finalModel$prediction.error))

# get the list of variable importance
rf_var_importance <- varImp(rfo_models$RandomForest_Bootstrap)$importance
#   plot(varImp(rfo_models$RandomForest_Bootstrap),top=25)


#### Logistic regression ####

print("Fitting a logistic regression model...")

# fit logistic regression model
# method and family defines the type of regression
# in this case these arguments mean that we are doing logistic
# regression
glm_models <- lapply(trctrl[1],
                     function(x){
                       fit_log_reg( train_response_factor,
                                    train_var_mx, 
                                    trcontrol = x)
                     }) %>%
  (function(lt){
    names(lt) <- str_c("GLM-", ncol(train_var_mx),"var_", names(lt))
    return(lt)
  })


#### Elastic net ####

print("Fitting an elastic net model...")

# we will now train elastic net model
# it will try
enet_models <- 
  lapply(trctrl,
         function(x){
           fit_elastic_net( train_response_factor,
                            train_var_mx,
                            trcontrol = x)
         }) %>%
  (function(lt){
    names(lt) <- str_c("ElasticNet-", ncol(train_var_mx),"var_", names(lt))
    return(lt)
  })


#### Variable importance filtering with Random Forest ####

RF_vi <- varImp(rfo_models $RandomForest_Bootstrap)$importance %>%
  rownames_to_column(var="gene") %>%
  dplyr::rename(imp = Overall) %>%
  mutate(gene_lab = gene) %>%
  arrange(desc(imp)) 

# Save variable importance matrix
write.table(RF_vi, "./out/Variable_importance.txt",
            sep="\t", col.names=T, row.names=F, quote=F)


# We do it if the number of variables is greater than 10
if (ncol(train_var_mx) > 30){

  print("Variable importance filtering...")
  
  RF_vi <- 
    RF_vi %>%
    (function(df){
      # If gene_lab is longer than 10 elements, keep label of only first 10
      df$gene_lab <- c(df$gene_lab[1:10], 
                       rep(NA,length(11:length(df$gene_lab))))
      return(df)}) %>%
    # We discretize the importance metric into high and low using k-means
    mutate(importance_cat = arules::discretize(imp, 
                                               method="cluster",
                                               breaks=4, 
                                               labels=c("VLow","Low","Mod","Imp"))
    )
}
  
pdf("./out/variable_importance.pdf", width = 6, height=6)

color <- if (ncol(train_var_mx) > 30) { as.character(RF_vi$importance_cat) 
  } else {  RF_vi$imp}

cowplot::plot_grid(nrow=1, ncol=2, align = "hv",
                   ggplot(RF_vi, 
                          aes(x=0, y=imp)) + 
                     geom_boxplot() + 
                     geom_point(aes(color=color)) +
                     geom_text_repel(aes(label=gene_lab),
                                     size = 3,
                                     max.overlaps = 20) +
                     theme_classic() +
                     theme(axis.line.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.title.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           legend.position = "none") +
                     labs(y = "importance",
                          title= "Variable Importance")
                   ,
                   ggplot(RF_vi,
                          aes(x=imp)) +
                     geom_density(fill="gray88") +
                     coord_flip() +
                     theme_classic() +
                     theme(axis.line.y = element_blank(),
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank())
)

dev.off()

#### New data with variable importance ####
if (ncol(train_var_mx) > 30){
  
  train_var_mx_imp <- 
    train_var_mx[,(RF_vi %>% filter(importance_cat == "Imp"))$gene]
  
  test_var_mx_imp <-
    test_var_mx[,(RF_vi %>% filter(importance_cat == "Imp"))$gene]
  
  
  #### Logistic regression with important genes ####
  
  print("Fitting a logistic regression model with top important genes...")
  
  glm_models_imp <- 
    lapply(trctrl[3],
           function(x){
             fit_log_reg( train_response_factor,
                          train_var_mx_imp,
                          trcontrol = x)
           }) %>%
    (function(lt){
      names(lt) <- str_c("GLM-", ncol(train_var_mx_imp),"var_", names(lt))
      return(lt) })
  
  # we will now train elastic net model
  # it will try
  enet_models_imp <- 
    lapply(trctrl,
           function(x){
             fit_elastic_net(train_response_factor,
                             train_var_mx_imp,
                             trcontrol=x)
             }) %>%
    (function(lt){
      names(lt) <- str_c("ElasticNet-", ncol(train_var_mx_imp),"var_", names(lt))
      return(lt)  })

}

#### Gradient boosting ####



#### Support vector machine ####





#### Merge all models in the same list ####

all_models <- 
  append(knn_models, 
         rfo_models
  ) %>% 
  append(glm_models) %>%
  append(enet_models)

if (ncol(train_var_mx) > 30){
  
all_models_imp <-
  append(enet_models_imp,
         glm_models_imp)

}

#### Generate the confusion matrices ####

print("Creating confusion matrix...")

# with the performance metrics of each model
all_confusion_mx <- 
  lapply(all_models,
           function(x){
             tt_confusion_matrix(train_var_mx,
                                 train_response_factor,
                                 test_var_mx,
                                 test_response_factor,
                                 x)
           }) %>%
      bind_rows(.id="idcol") %>%
      separate(idcol, into=c("method", "split_method"), sep="_") %>%
      mutate(method = 
               ifelse(str_detect(method, "ElasticNet"),
                      str_c(method, "_", split_method),
                      method))
  
    
if (ncol(train_var_mx) > 30){
  all_confusion_mx <-
    bind_rows(
      all_confusion_mx
      ,
      lapply(all_models_imp,
             function(x){
               tt_confusion_matrix(train_var_mx_imp,
                                   train_response_factor,
                                   test_var_mx_imp,
                                   test_response_factor,
                                   x)
             }) %>%
        bind_rows(.id="idcol") %>%
        separate(idcol, into=c("method", "split_method"), sep="_") %>%
        mutate(method = 
                 ifelse(str_detect(method, "ElasticNet"),
                        str_c(method, "_", split_method),
                        method))
    )
}

write.table(all_confusion_mx, "./out/Confusion_matrix.txt",
            sep="\t", col.names=T, row.names=F, quote=F)

#### ROC Curves ####

print("Potting ROC curves and AUC...")

roc_curves <- 
    # Models with original matrix of expression
    lapply(names(all_models)[!str_detect(names(all_models), "RandomForest")],
           function(x){
             tt_roc_curve(train_var_mx,
                          train_response_factor,
                          test_var_mx,
                          test_response_factor,
                          all_models[[x]],
                          x)
           }) %>%
      bind_rows(.id="split_method") %>%
      mutate(cohort = str_to_sentence(cohort)) %>%
      mutate(cohort = factor(cohort, levels=c("Train", "Test")))
    
if (ncol(train_var_mx) > 30){
  roc_curves <- 
    bind_rows(
      roc_curves
      ,
        # Models with importance-filtered matrix of expression
        lapply(names(all_models_imp)[!str_detect(names(all_models_imp), 
                                                 "RandomForest")],
               function(x){
                 tt_roc_curve(train_var_mx_imp,
                              train_response_factor,
                              test_var_mx_imp,
                              test_response_factor,
                              all_models_imp[[x]],
                              x)
               }) %>%
          bind_rows(.id="split_method") %>%
          mutate(cohort = str_to_sentence(cohort)) %>%
          mutate(cohort = factor(cohort, levels=c("Train", "Test")))
      )
}

pdf("./out/ROC_curves.pdf", 
    width = if (length(unique(roc_curves$model)) > 15) {14
    }  else {12}, 
    height=4.5
    )

ggplot(roc_curves,
       aes(x=OmSpecificity,
           y=Sensitivity,
           group=test,
           color=test)) +
  geom_line(size=0.8) +
  theme_bw() +
  facet_wrap(~cohort)  +
  labs(x = "1 - Specificity")

dev.off()


#### Save all models ####

save(all_models, file = "models.RData")


#### Finalize ####

print("End of run!!")