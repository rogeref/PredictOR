# PredictOR
**Train models for binary category prediction (Objective Response to therapy)**

## Summary

PredictOR allows the user to fit different predictive models using a training and validation cohort of choice. In short, the tool uses the training set to fit logistic regression (with and without elastic net regularization), k-nearest neighbours and random forest prediction models with cross-validation techniques (10-fold, leave-one-out and bootstrapping). Furthermore, it provides information on predictive performance in train/validation sets and returns the metrics of variable importance for each model (see the repository for more information). The initial version of this code aims to fit models to discriminate between responders (patients with Objective Response) and non-responders to atezolizumab + bevacizumab therapy based on bulk RNA-seq transcriptomic data from hepatocellular carcinoma tumors. 

## How to run

The code in its current form can be executed from the terminal, inside the PredictOR directory, with the following command:

```
Rscript model_functions.R
Rscript model_fitting.R
```

Make sure to have the following files witin the ./in/ directory: 
* The expression matrices of the train and test cohorts (with the names *train_gene_matrix.txt* and *test_gene_matrix.txt*). 
* The classes to predict (with the names *train_response.txt* and *test_response.txt*).

An example fo the format that these files need to have can be found in the files ./in/example_gene_matrix.txt and ./in/example_response.txt.

## What you will get

The tool first cleans the data removing variables under the following criteria: variables with any missing value, with near-zero variation (where 85% of observations have the same values, found with `preProcess` in the R package `caret` v.6.0-93), with a high correlation with other variables (cutoff of 0.9 in caret’s `findCorrelation`) and variables found only in one of the train or test datasets. Subsequently, data is centered and scaled (Z-score). PredictOR then generates an exploratory data analysis including correlation plots, PCA and the assessment of value distributions (with `corrplot` v0.92, `ggbiplot` v0.55 and `ggridges` v0.5.4). Subsequently, it uses the training dataset to fit logistic regression (with and without elastic net regularization), k-nearest neighbors and random forest prediction models with cross-validation techniques (10-fold, leave-one-out and bootstrapping) leveraging the functionalities of the caret R package (predict function with `method`=”glmnet”, “glm”, “ranger” or “knn” for elastic net, logistic regression, random forest and k-nearest neighbors, respectively). Furthermore, PredictOR provides information on predictive performance in the train and validation datasets, including accuracy, sensitivity, specificity, positive and negative predictive values and precision among other metrics (from caret’s `confusionMatrix`), it returns the metrics of variable importance (from caret’s `varImp`) and it also provides ROC curves and AUC values for each model (from `pROC` v1.18.0).
