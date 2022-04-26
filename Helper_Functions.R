#Title:	Helper functions for preprocessing data and training models
#Project Title:	Rapid determination of microbial drug resistance using FTIR spectroscopy
#Author: Hewa Godage Sithija Wijesinghe
#Author contact: hewa.wijesinghe@uqconnect.edu.au
#Script info: Contains the helper functions to analyse the Microbial spectra
#Usage:	source("./Helper_Functions.R")
#Copyright statement: TBA.

#---------------------------------------------------------------------------------------------------------------
#set the library location for R
#The libraries I used in the HPCs are stored locally and the directory needs to be set in each run
set_local_libPath <- function(path){
  mypaths <- .libPaths()
  mypaths <- c(path, mypaths)
  .libPaths(mypaths)
}

set_local_libPath('/home/s4187725/R/x86_64-pc-linux-gnu-library/3.5')

###########################################
#read list of files from the working directory

#add labels tag to files using bash
#$ for x in *_HOPS*.spc; do mv $x ${x%.spc}_HOPS.spc;done
#$ for x in *_MER*.spc; do mv $x ${x%.spc}_MER.spc;done
#$ for x in *_ATCC*.spc; do mv $x ${x%.spc}_QC.spc;done

############################################
#required libraries
library(hyperSpec)
library(signal)
library(tidyr)
library(dplyr)
library(caret)
library(mvoutlier)
library(pROC)
library(ggplot2)


#############################################
#generate metadata matrix based on the filenames

generate_metadata <- function (files, format, substitute){
  metadata = data.frame(filename=files, file_processed=sub(substitute, "\\1", files)) %>% 
    separate(file_processed,
             into = format,
             sep = "_"
    ) %>% 
    group_by(isolate, date, initials) %>% 
    mutate(replicate_no=1:n())
  
  metadata$group <- group_indices(metadata)
  
  return(metadata)
}

#-----------------------------------------------------------------------------------------------
#reads individual spc files into a hyperSpec object
#contains a wrapper for binning and filtering functions if needed
process_spc <- function(file) {
  raw <- read.spc(file)
  #processed <- spc.bin(raw, 1, na.rm=TRUE) #binning factor of 1 = no binning
  #processed <- apply(raw, 1, sgolayfilt)
  #raw <- processed
  return(raw)
}


#reads files off the working directory, processes one by one and appends to large HyperSpec object
read_spectra <- function (metadata){
  
  spectra <- NULL
  
  for (filename in metadata$filename){
    raw <- process_spc(filename)
    
    if (!is.null(spectra)) {
      spectra <- rbind(spectra, raw)
    } else {
      spectra <- raw
    }
  }
  
  return(spectra)
  
}

normalise_spectra <- function(spc, method="amideI"){
  #' Normalises the hyperspec matrix. 
  #' 
  #' Methods available: 
  #'     Normalisation by Amide I bond: amideI
  #'     Normalisation by area: area
  #'     Min/Max scaling (0-1): minmax
  #' 
  if (method == "amideI"){
    #normalise by Amide I band. Butler, 2018. page 5.
    print('normalise by Amide I bond')
    factors <- 1 / apply (spc [, , 1600 ~ 1700], 1, mean)
    spectra.normalised <- sweep (spc, 1, factors, "*")
    return(spectra.normalised)
  } else if (method == 'area'){
    print('normalise by area')
    spectra.normalised <- sweep(spc, 1, mean, "/")
    return(spectra.normalised)
  } else if (method == 'minmax') {
    print('min-max scaling')
    #perform baseline correction
    bl <- spc.fit.poly.below(spc)
    base.corrected <- spc - bl
    #scaling
    spectra.normalised <- sweep(spc, 1, max, "/")
    return(spectra.normalised)
  } else {
    print("Invalid method")
    return(spc)
  }
}


filter_spectra <- function(spc, method="savgol", derivative=1, window=5){
  #' Filters the the hyperspec matrix. 
  #' 
  #' Methods available: 
  #'     "savgol": Savitzskly-Golay Filtering.
  #'     m = derivative, n = filter length
  #' 
  if (method=="savgol"){
    filtered<- apply(spc, 1, sgolayfilt, m=derivative, n=window)
  } else {
    print("Invalid method")
    return(spc)    
  }
}


#baseline correction (polynomial correction. Butler, 2018 page 5)
#bl <- spc.fit.poly.below(spectra.normalised)
#spectra.normalised.baseline <- spectra.normalised - bl


##############################
#train the model
###############################

#takes a string as a param and trains the corresponding method through caret::train()
process_model <- function(method, trainData, labels, metric) {
  #' Metric can be 'Accuracy' or 'ROC'
  set.seed(seed)
  if (metric == 'Accuracy'){
    file <- train(
    trainData,
    labels,
    method = method,
    preProc = c("center", "scale"),
    metric = 'Accuracy'
  )
  } else if (metric == 'ROC'){
    file <- train(
    trainData,
    labels,
    method = method,
    preProc = c("center", "scale"),
    metric = 'ROC',
    trControl=trainControl(method="repeatedcv", number=10, repeats=5, classProbs = T, summaryFunction = twoClassSummary)
  )
  } else {
    file <- NULL
  }
  
  return(file)
}



#installs dependencies for each caret model
install_dependencies <- function(method) {
  
  mod <- getModelInfo(method)
  for (i in 1:length(mod)){
    #print(c('method:  ', method))
    for (j in 1:length(mod[[i]]$library)){
      pkg <- mod[[i]]$library[[j]]
      #print(c(i, j, pkg))
      if (length(pkg) != 0){
        if (!(pkg %in% installed.packages())){
          try(install.packages(pkg))
        }
      }
    }
  }
}

attempt_process <- function(method, trainData, labels, metric) {
  install_dependencies(method)
  try(process_model(method, trainData, labels, metric))
}


####trying to pull accuracy element from a trained caret model object
grab_accuracy <- function(trainobj){
  try(max.accuracy<- max(trainobj[[4]][[4]]))
  try(return(max.accuracy))
}

#training function for hyperparameter optimization
controlled_train <- function(trainData, trainLabels, model, tune = 10, cvrepeat=20){
  #' trains a model under varying hyperparameter combinations and returns the trained model
  control.random <- trainControl(method="repeatedcv", repeats=cvrepeat, classProbs = T, savePredictions = TRUE,
                                 summaryFunction = twoClassSummary, search = 'random')
  control.grid <- trainControl(method="repeatedcv", repeats=cvrepeat, classProbs = T, savePredictions = TRUE,
                               summaryFunction = twoClassSummary, search = 'grid')
  if (is.data.frame(tune)){
    trained.model <- train(
      trainData,
      trainLabels,
      method = model,
      preProc = c("center", "scale"),
      metric = 'ROC',
      tuneGrid =  tune,
      trControl=control.grid)
  } else if (is.numeric(tune)) {
    trained.model <- train(
      trainData,
      trainLabels,
      method = model,
      preProc = c("center", "scale"),
      metric = 'ROC',
      tuneLength =  tune,
      trControl=control.random)
  } else {
    trained.model <- train(
      trainData,
      trainLabels,
      method = model,
      preProc = c("center", "scale"),
      metric = 'ROC',
      trControl=control.random)
  }
  return(trained.model)
}

######################################
#running exploraratory test to find the optimal model
#train.prelim <- sapply(caret.models, attempt_process)
#train.accuracy <- sapply(train.accuracy, grab.accuracy)


attempt_predict <- function(model, testing, type){
  predictions <-NULL
  #try(predictions <- predict(model, newdata = testing))
  try(predictions <- predict(model, newdata = testing, type=type))
  return(predictions)
}

generate_conf <- function(prediction, labels, pos){
  conf.mat <- NULL
  try(conf.mat<- confusionMatrix(table(prediction, as.factor(labels)), positive=pos))#, dnn = c("Prediction", "Reference"), norm='none', mode = "sens_spec")))
  return(conf.mat)
}


summary.stats.colnames <- c("Method", "Metric", "CVMetric", "HoldOutMetric",
                            "Kappa",  "P value", "McNemarsPValue", "Sensitivity", "Specificity", "Precision",
                            "TrueNegatives", "FalseNegatives", "FalsePositives", "TruePositives")

extract_summary_stats <- function(trainobj, confus.mat, summaryColnames){
  #' takes a Caret train object and a confusionMatrix object as input
  #' extracts summary statistics from the objects and returns a vector
  #' currently extracted stats:
  #' "Method", "Metric", "Cross validation Accuracy/ROC", "Hold out set Accuracy", "Kappa",  "P value",
  #' "McNemar's Test P-Value", "Sensitivity", "Specificity", "Precision",
  #' "True Negatives", "False Negatives", "False Positives", "True Positives"
  
  stats <- list()
  
  stats[[1]] <- trainobj[[1]] #name
  stats[[2]] <- trainobj[[9]] #metric
  stats[[3]] <- max(trainobj[[4]][[2]]) #Cross validation metric value
  stats[[4]] <- confus.mat$overall[1] #hold out set accuracy
  stats[[5]] <- confus.mat$overall[2] #kappa
  stats[[6]] <- confus.mat$overall[6] #Pvalue
  stats[[7]] <- confus.mat$overall[2] #McNemar
  stats[[8]] <- confus.mat$byClass[1] #Sensitivity
  stats[[9]] <- confus.mat$byClass[2] #Specificity
  stats[[10]] <- confus.mat$byClass[5] #Precision
  stats[[11]] <- as.numeric(confus.mat$table)[1] # TN
  stats[[12]] <- as.numeric(confus.mat$table)[2] # FN
  stats[[13]] <- as.numeric(confus.mat$table)[3] # FP
  stats[[14]] <- as.numeric(confus.mat$table)[4] # TP
  
  stats <- as.data.frame(stats)
  names(stats) <- summaryColnames
  row.names(stats) <- c(trainobj[[1]])
  
  return(stats)
}


return_summary_table <- function(trainobj, confus.mat, summaryColnames){
  #' If a summary table cannot be extracted from the model, returns a blank one.
  tryCatch(
    {stats <- extract_summary_stats(trainobj, confus.mat, summaryColnames)
    return(stats)
  }, error = function(cond) {
    name <- trainobj[[1]]
    stats <- data.frame(matrix(c(name, rep(NA, 13)),nrow=1,ncol=14), stringsAsFactors=FALSE)
    colnames(stats) <- summaryColnames
    row.names(stats) <- c(name)
    return(stats)
    
  }
  )
}


summary_with_model <- function(trainobj, confus.mat, summaryColnames){
  #' function still under construction
  #' The idea is to append the model object itself as an element and return a tibble object
  summary.table <- extract_summary_stats(trainobj, confus.mat, summaryColnames)
}





#-----------------------------------------------------------------------------------------------------
#Functions for the the main method call


model_training <- function(name, trainData, trainLabels, metric, testData, testLabels){
  #' trains models in caret, tests and exports the trained model and confusion matrix as .RDS files
  #' 
  #' Args:
  #' name: name of the model
  #' trainData: numeric matrix for train data
  #' trainLabels: corresponding labels vector
  #' metric: cross validation metric for training the model
  #' testData: numeric matrix for test data
  #' testLabels: corresponding labels vector
  #' 
  tryCatch(
    {
      #train the model
      model.name <- paste("model", name, sep=".")
      print(paste("currently training : ", name))
      assign(model.name, attempt_process(name, trainData, trainLabels, metric))
      trained.model <- eval(parse(text=paste(model.name)))
      
      #made predictions with the test data
      predictions <- paste("prediction", name, sep = '.')
      assign(predictions, attempt_predict(trained.model, testData, type='prob'))
      
      #generate confusion matrix
      confusion.mat <- paste("confusionMatrix", name, sep = '.')
      assign(confusion.mat, generate_conf(eval(parse(text=paste(predictions))), as.factor(testLabels), 'MER'))
      confus.mat <- eval(parse(text=paste(confusion.mat)))
      
      #save outputs
      filename.model <- paste("train_", name, ".RDS", sep='')
      filename.conf <- paste("conf_", name, ".RDS", sep='')
      saveRDS(trained.model, file=filename.model)
      saveRDS(confus.mat, file=filename.conf)
      

    }, error = function(cond) {
      message(paste("Error returned in training method", model.name))
      message(cond)
      
    }, warning = function(cond) {
      message(paste("Warning returned in training method", model.name))
      message(cond)
      
    }, finally={
      #print('done')
      
    }
  )
}



#######################################################################################
#functions to wrap preprocessing 

preprocess <- function(spc, m=1, w=5){
  spectra <- spc
  spectra <- filter_spectra(spectra, method="savgol", derivative=m, window=w)
  spectra <- normalise_spectra(spectra)
  return(spectra)
}
  
#remove outliers
remove_outliers <- function(spectra) {
  outliers <- pcout(spectra$spc, makeplot = F, explvar = 0.99, crit.M1 = 1/20, crit.c1 = 10,crit.M2 = 1/20, crit.c2 = 0.99, cs = 0.25, outbound = 0.05)
  return (which(outliers$wfinal01==1))
}


dataset_split <- function(spectra, file_df) {
  
  #grouping isolates
  isolates <- file_df[, c("label", "isolate")]
  isolates.uniq <- unique(isolates)
  #isolates <- isolates[-outlier.indices]
  
  
  #extract the data matrix from the hyperspec object
  spectra.data <- spectra$spc
  colnames(spectra.data) <-  spectra@wavelength
  row.names(spectra.data) <- spectra$filename
  
  #generate the data split per isolate
  set.seed(seed)
  inTrain <- createDataPartition(
    y=isolates.uniq$label,
    p=0.7,
    list = FALSE
  )
  
  #split data by isolate
  train.isolates <- isolates.uniq$isolate[inTrain]
  test.isolates <- isolates.uniq$isolate[-inTrain]
  
  split <- list()
  split[[1]] <- spectra.data[file_df$isolate %in% train.isolates,]
  split[[2]]  <- spectra.data[file_df$isolate %in% test.isolates,]
  split[[3]] <- file_df$label[file_df$isolate %in% train.isolates]
  split[[4]]  <- file_df$label[file_df$isolate %in% test.isolates]
  names(split) <- c("training", "testing", "train.labels", "test.labels")
  return(split)
}



##########################
  
generate_filename <- function(name, derivative, window, wl.range)  {
  outfile <- paste(name, wl.range[1], wl.range[2], derivative, window, sep='_')
  return(outfile)
}

model_training_tune <- function(name, split, outfile){
  #' same as model_training, but takes output filename as a parameter instead of automatically generating it
  #' 
  #' Args:
  #' name: name of the model
  #' trainData: numeric matrix for train data
  #' trainLabels: corresponding labels vector
  #' metric: cross validation metric for training the model
  #' testData: numeric matrix for test data
  #' testLabels: corresponding labels vector
  #' 
  tryCatch(
    {
      #extract the train/test datasets and labels
      trainData <- split$training
      testData <- split$testing
      trainLabels <- split$train.labels
      testLabels <- split$test.labels
      
      #train the model
      model.name <- paste("model", name, sep=".")
      print(paste("currently training : ", name))
      assign(model.name, attempt_process(name, trainData, trainLabels, "ROC"))
      trained.model <- eval(parse(text=paste(model.name)))
      
      #made predictions with the test data
      predictions <- paste("prediction", name, sep = '.')
      assign(predictions, attempt_predict(trained.model, testData, type='prob'))
      
      #generate confusion matrix
      confusion.mat <- paste("confusionMatrix", name, sep = '.')
      assign(confusion.mat, generate_conf(eval(parse(text=paste(predictions))), as.factor(testLabels), 'MER'))
      confus.mat <- eval(parse(text=paste(confusion.mat)))
      
      #save outputs
      filename.model <- paste("train_", outfile, ".RDS", sep='')
      filename.conf <- paste("conf_", outfile, ".RDS", sep='')
      saveRDS(trained.model, file=filename.model)
      saveRDS(confus.mat, file=filename.conf)
      
      
    }, error = function(cond) {
      message(paste("Error returned in training method", model.name))
      message(cond)
      
    }, warning = function(cond) {
      message(paste("Warning returned in training method", model.name))
      message(cond)
      
    }, finally={
      #print('done')
      
    }
  )
}


#---------------------------------------------------------------------------------------------------------------------
#plotting functions used 

draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Susceptible', cex=1.4)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Resistant', cex=1.4)
  text(125, 370, 'Predicted', cex=1.5, srt=90, font=2)
  text(245, 450, 'True', cex=1.5, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Susceptible', cex=1.4, srt=90)
  text(140, 335, 'Resistant', cex=1.4, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.5, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.4)
  text(30, 85, names(cm$byClass[2]), cex=1.5, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.4)
  text(50, 85, names(cm$byClass[5]), cex=1.5, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.4)
  text(70, 85, names(cm$byClass[6]), cex=1.5, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.4)
  text(90, 85, names(cm$byClass[7]), cex=1.5, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.4)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  


draw_hyperparameter_plots <- function(confusion.mat, prediction.roc, trained.model, tuneparameter){
  #' draws the important plots in the hyperparameter optimisation step and saves them in the R workspace
  #' 
  #' 
  
  #confusion matrix
  draw_confusion_matrix(confusion.mat)
  #Plot the ROC curve. will also return the plot object
  aucval <- paste0("Area under the curve: ", round(prediction.roc$auc, digits=3))
  grob <- grobTree(textGrob(aucval, x=0.4,  y=0.20, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  ggroc(prediction.roc,
                legacy.axes = FALSE)+ theme_minimal() +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") + 
    annotation_custom(grob)  + theme(text = element_text(size=18)) + xlab('Specificity') + ylab('Sensitivity')
  
  #plot the hyperparameter distribution
  trellis.par.set(caretTheme())
  plot(trained.model, metric='ROC' , plotType = "level",
       scales = list(x = list(rot = 90)))
  if (is.data.frame(tuneparameter)){
    ggplot(trained.model)+ scale_x_continuous(trans='log10') 
  } else {
    ggplot(trained.model)
  }
  return(plot)
}
