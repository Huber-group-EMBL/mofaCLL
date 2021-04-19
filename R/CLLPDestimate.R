#codes for CLLPDestimate function

processTrainingData <- function(exprObj, identifier) {
  if (identifier == "ensembl_gene_id") {
    X <- assay(exprObj)
    y <- exprObj$y
  } else if (identifier == "gene_symbol") {
    subObj <- exprObj[!rowData(exprObj)$symbol %in% c("",NA)]
    subObj <- subObj[!duplicated(rowData(exprObj)$symbol),]
    X <- assay(subObj)
    rownames(X) <- rowData(subObj)$symbol
    y <- subObj$y
  }
  return(list(X=X,y=y))
}

#' Estimate CLL-PD in user-specified expression dataset
#'
#' \code{CLLPDestimate} returns the estimated CLL-PD value in the user-specified cohort.
#'
#' This function takes an gene expression dataset (RNAseq or microarray) of a external CLL cohort (user-specified),
#' build a regularized linear model using the expression values of overlapped features in the built-in training cohort and
#' use the selected model to esimate CLL-PD in the external cohort.
#'
#' @param exprMatrix A numeric matrix that contains the expression level of genes or probes for samples in a cohort. Rownames are gene identifiers and column names are sample identifiers.
#' @param identifier A charachater variable that specifies the type of gene identifier. Currently only "ensembl_gene_id" and "gene_symbol" are allowed.
#' @param topVariant If specified, it should be a numeric value indicating the number of most variant features (genes or probes) in the user-specificed data. The default value is NULL, which means all rows in the input matrix will be used.
#' @param normalize A boolean value indicating wether the user-specified expression matrix should be centered by mean and scaled by standard deviation. The default value is TRUE.
#' @param repeats A numeric variable specifying the number of repeats for cross-validation to select prediction model. The default value is 20.
#'
#' @return A list containing three objects 1) estimated_CLLPD: A numeric vector of the estimated CLL-PD values in the user-specified cohort;
#' 2) A dataframe of the features with non-zero coefficients and their coefficients in the selected model (model with highest R2 value);
#' 3) A numeric vector of variance explained (R2) values for CLL-PD of the built-in cohort along the repeated cross-validation runs.
#'
#' @import glmnet DESeq2 SummarizedExperiment
#' @export
#'
CLLPDestimate <- function(exprMatrix, identifier = "ensembl_gene_id", topVariant = NULL, normalize = TRUE, repeats = 20) {

  cat("Preparing model. \n")
  data("rnaTrain")

  #Pre-process the training data, change row names to specified gene identifier
  trainData <- processTrainingData(rnaTrain, identifier)

  #normalize the training data
  trainX <- mscale(trainData$X)
  trainY <- trainData$y

  #process input matrix
  testX <- as.matrix(exprMatrix)
  if (!is.null(topVariant)) {
    if (topVariant > nrow(testX)) {
      message("The number specified is larget than the rows of input matrix. All data will be used.")
    } else {
      sds <- apply(testX, 1, sd)
      testX <- testX[order(sds,decreasing = T)[seq(topVariant)],]
    }
  }
  if (normalize) testX <- mscale(testX)

  #subset for common genes
  commonGene <- intersect(rownames(trainX), rownames(testX))
  trainX <- t(trainX[commonGene,])
  trainX <- removeCorrelated(trainX, cutoff = 0.9, record = FALSE)$reduced

  #train model
  cat("Running repeated cross-validation to select features.")
  lassoRes <- runGlm(trainX,trainY, method = "lasso", repeats=20, folds=5, lambda ="lambda.1se")
  useModel <- lassoRes$modelList[[which.max(lassoRes$r2Train)]]
  cat("\nDone.")
  #predict
  testY <- glmnet:::predict.cv.glmnet(useModel, newx = t(testX[colnames(trainX),]))[,1]

  #process output
  ##feature coefficient
  coefMat <- as.matrix(coef(useModel)[-1,])
  coefMat <- data.frame(coefMat[coefMat[,1]!=0,,drop=FALSE])
  colnames(coefMat) <- "coefficient"

  return(list(estimated_CLLPD = testY,
              featureCoefficient = coefMat,
              trainingR2 = lassoRes$r2Train))
}
