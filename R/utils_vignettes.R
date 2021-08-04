
########################################
###Functions for running the analysis###
########################################

# Function for scale a data matrix
#' @export
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y))/(1.4826*mad(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y))/(sd(y)))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/(1.4826*mad(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) y/(sd(y)))
    }
  } else {
    x.scale = x
  }

  if (!is.null(censor)) {
    x.scaled[x.scaled > censor] <- censor
    x.scaled[x.scaled < -censor] <- -censor
  }
  return(t(as.matrix(x.scaled)))
}

#modified version of plotDataOver view from MOFA package, with black borders in the heatmap
#' @export
plotDataOverview.m <- function (object, colors = NULL)
{
  if (!is(object, "MOFAmodel"))
    stop("'object' has to be an instance of MOFAmodel")
  TrainData <- getTrainData(object)
  M <- getDimensions(object)[["M"]]
  N <- getDimensions(object)[["N"]]
  if (is.null(colors)) {
    palette <- c("#D95F02", "#377EB8", "#E6AB02", "#31A354",
                 "#7570B3", "#E7298A", "#66A61E", "#A6761D", "#666666",
                 "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                 "#A65628", "#F781BF", "#1B9E77")
    if (M < 17)
      colors <- palette[seq_len(M)]
    else colors <- rainbow(M)
  }
  if (length(colors) != M)
    stop("Length of 'colors' does not match the number of views")
  names(colors) <- sort(viewNames(object))
  ovw <- vapply(TrainData, function(dat) apply(dat, 2, function(s) !all(is.na(s))),
                logical(N))
  ovw <- ovw[apply(ovw, 1, any), , drop = FALSE]
  #molten_ovw <- melt(ovw, varnames = c("sample", "view"))
  molten_ovw <- data.frame(ovw) %>% rownames_to_column("sample") %>%
    gather(key = "view", value = "value", -sample) %>% data.frame(stringsAsFactors = FALSE)
  molten_ovw$sample <- factor(molten_ovw$sample, levels = rownames(ovw)[order(rowSums(ovw),
                                                                              decreasing = TRUE)])
  n <- length(unique(molten_ovw$sample))
  molten_ovw$combi <- ifelse(molten_ovw$value, as.character(molten_ovw$view),
                             "missing")
  molten_ovw$ntotal <- paste("n=", colSums(ovw)[as.character(molten_ovw$view)],
                             sep = "")
  molten_ovw$ptotal <- paste("d=", vapply(TrainData, nrow,
                                          numeric(1))[as.character(molten_ovw$view)], sep = "")
  molten_ovw$view_label = paste(molten_ovw$view, molten_ovw$ptotal,
                                sep = "\n")
  molten_ovw$label_pos <- levels(molten_ovw$sample)[n/2]
  p <- ggplot(molten_ovw, aes_string(x = "sample", y = "view_label",
                                     fill = "combi")) + geom_tile(color = "grey50") + geom_text(data = dplyr::filter(molten_ovw,
                                                                                                sample == levels(molten_ovw$sample)[1]), aes_string(x = "label_pos",
                                                                                                                                                    label = "ntotal"), size = 6) + scale_fill_manual(values = c(missing = "grey",
                                                                                                                                                                                                                colors)) + xlab(paste0("Samples (n=", n, ")")) + ylab("") +
    guides(fill = FALSE) + theme(panel.background = element_rect(fill = "white"),
                                 text = element_text(size = 16), axis.ticks = element_blank(),
                                 axis.text.x = element_blank(), axis.text.y = element_text(color = "black"),
                                 panel.grid = element_blank(), plot.margin = unit(c(5.5,
                                                                                    2, 5.5, 5.5), "pt"))
  return(p)
}


#function for cox regression
#' @export
comSurv <- function(response, time, endpoint, scale =FALSE) {

  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)


  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}

#Function to format floats
#' @export
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}


# Function for Kaplan-Meier plot
#' @export
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (years)",
               table_ratio = c(0.7,0.3), yLabelAdjust = 0) {
  #function for km plot
  survS <- tibble(time = time,
                  endpoint = endpoint)

  if (!is.null(maxTime))
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))

  if (stat == "maxstat") {
    ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                       data = survS,
                       smethod = "LogRank",
                       minprop = 0.2,
                       maxprop = 0.8,
                       alpha = NULL)

    survS$group <- factor(ifelse(response >= ms$estimate, "high", "low"))
    p <- comSurv(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "median") {
    med <- median(response, na.rm = TRUE)
    survS$group <- factor(ifelse(response >= med, "high", "low"))
    p <- comSurv(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "binary") {
    survS$group <- factor(response)
    if (nlevels(survS$group) > 2) {
      sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
      p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    } else {
      p <- comSurv(survS$group, survS$time, survS$endpoint)$p
    }
  }

  if (is.null(pval)) {
    if(p< 1e-16) {
      pAnno <- bquote(italic("P")~"< 1e-16")
    } else {
      pval <- formatNum(p, digits = 1)
      pAnno <- bquote(italic("P")~"="~.(pval))
    }

  } else {
     pval <- formatNum(pval, digits = 1)
     pAnno <- bquote(italic("P")~"="~.(pval))
  }

  if (!showP) pAnno <- ""

  colorPal <- colList[1:length(unique(survS$group))]
  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE, palette = colorPal,
                  legend = ifelse(showTable, "none","top"),
                  ylab = "Fraction", xlab = "Time (years)", title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, legend.labs = sort(unique(survS$group)),
                  ggtheme = theme_half + theme(plot.title = element_text(hjust =0.5),
                                               panel.border = element_blank(),
                                               axis.title.y = element_text(vjust =yLabelAdjust)))
  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5)
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size=5)
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}

#function for italicize gene names
#' @export
italicizeGene <- function(featureNames) {
  geneNameList <- c("SF3B1","NOTCH1","TP53")
  specialCase <- "TP53/del17p"
  formatedList <- sapply(featureNames, function(n) {
    if (n %in% geneNameList) {
      sprintf("italic('%s')",n)
    } else if (n == specialCase) {
      sprintf("italic('TP53')~'/del17p'")
    } else if (n == "CLLPD") {
      "CLL-PD"
    }
    else {
      sprintf("'%s'",n)
    }
  })
  return(formatedList)
}

#function for plot hazard ratio
#' @export
plotHazard <- function(survRes, title = "") {
  sumTab <- summary(survRes)$coefficients
  confTab <- summary(survRes)$conf.int
  #correct feature name
  nameOri <- rownames(sumTab)
  nameMod <- substr(nameOri, 1, nchar(nameOri) -1)
  plotTab <- tibble(feature = rownames(sumTab),
                    nameMod = substr(nameOri, 1, nchar(nameOri) -1),
                    HR = sumTab[,2],
                    p = sumTab[,5],
                    Upper = confTab[,4],
                    Lower = confTab[,3]) %>%
    mutate(feature = ifelse(nameMod %in% names(survRes$xlevels), nameMod, feature)) %>%
    mutate(feature = str_replace(feature, "[.]","/")) %>%
    mutate(feature = str_replace(feature, "[_]","-")) %>%
    arrange(desc(abs(p))) %>%
    mutate(feature = italicizeGene(feature)) %>%
    mutate(feature = factor(feature, levels = feature)) %>%
    mutate(type = ifelse(HR >1 ,"up","down")) %>%
    mutate(Upper = ifelse(Upper > 10, 10, Upper))

  ggplot(plotTab, aes(x=feature, y = HR, color = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    geom_point(position = position_dodge(width=0.8), size=3, color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1,color = "grey20") +
    geom_text(position = position_nudge(x = 0.3),
              aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                           formatNum(p, digits = 1))),
              color = "black", size =4.5, parse = TRUE) +
    expand_limits(y=c(-0.5,0))+
    #scale_color_manual(values = c(up = colList[1], down = colList[2])) +
    ggtitle(title) + scale_y_log10() +
    scale_x_discrete(labels = parse(text = levels(plotTab$feature))) +
    ylab("Hazard ratio") +
    coord_flip() +
    theme_full +
    theme(legend.position = "none", axis.title.y = element_blank())
}

#Function to run multivariate Cox model on test table
#' @export
runCox <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    dplyr::filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Modified version of plot variance explained, smaller title
#' @export
plotVarianceExplained.m <- function (object, cluster = TRUE, censor = NULL, ...) {
  R2_list <- calculateVarianceExplained(object, ...)
  fvar_m <- R2_list$R2Total
  fvar_mk <- R2_list$R2PerFactor
  fvar_mk_df <- reshape2::melt(fvar_mk, varnames = c("factor",
                                                     "view"))
  if (!is.null(censor)) {
    fvar_mk_df$value <- ifelse(fvar_mk_df$value > censor, censor, fvar_mk_df$value)
  }
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  if (cluster & ncol(fvar_mk) > 1) {
    hc <- hclust(dist(t(fvar_mk)))
    fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  #change factor name from LF* to F*
  fvar_mk_df$factor <- gsub("LF"," F",fvar_mk_df$factor)

  hm <- ggplot(fvar_mk_df, aes_string(x = "view", y = "factor")) +
    geom_tile(aes_string(fill = "value"), color = "black") +
    #guides(fill = guide_colorbar("R2")) +
    scale_fill_gradientn(colors = c("white","darkblue"), guide = "colorbar") + ylab("Latent factor") +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 13), axis.line = element_blank(),
          axis.ticks = element_blank(), panel.background = element_blank())
  hm <- hm + ggtitle("Variance explained per factor") + guides(fill = guide_colorbar(title = bquote(R^2)))
  fvar_m_df <- data.frame(view = factor(names(fvar_m), levels = names(fvar_m)),
                          R2 = fvar_m)
  if (cluster == TRUE & ncol(fvar_mk) > 1) {
    fvar_m_df$view <- factor(fvar_m_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  bplt <- ggplot(fvar_m_df, aes_string(x = "view", y = "R2")) +
    ggtitle("Total variance explained per view") + geom_bar(stat = "identity",
                                                            fill = "deepskyblue4", width = 0.9) + xlab("") + ylab(bquote(R^2)) +
    scale_y_continuous(expand = c(0.01, 0.01)) + theme(plot.margin = margin(0,2.5, 0, 0, "cm"),
                                                       panel.background = element_blank(),
                                                       plot.title = element_text(size = 13, hjust = 0.5),
                                                       axis.ticks.x = element_blank(),
                                                       axis.text.x = element_blank(),
                                                       axis.text.y = element_text(size = 12, color = "black"),
                                                       axis.title.y = element_text(size = 13,color = "black"),
                                                       axis.line = element_line(size = rel(1),color = "black"))
  p <- plot_grid(bplt, hm, align = "v", nrow = 2, rel_heights = c(0.3,0.7), axis = "l")
  return(p)
}

#Function to generate pretty scientific notation format for plot label
#' @export
sciPretty <- function(n, digits = 2, bold = FALSE) {
  nForm <- strsplit(format(n, digits = digits, scientific = TRUE),split = "e")
  b <- nForm[[1]][1]
  i <- as.integer(nForm[[1]][2])
  #bquote(.(b)%*%10^.(i))
  if(bold) {
    sprintf("bold(%s%%*%%10^%s)",b,i)
  } else sprintf("%s%%*%%10^%s",b,i)
}

#function to remove highly correlated values and keep track of the removing
#' @export
removeCorrelated <- function(x, cutoff = 0.6, method = "pearson", keep = NULL, record = TRUE) {

  if (!is.null(keep)) {
    #if specified some feature to keep, then reorder the matrix, to make sure they are not collapsed
    posNow <- grepl(paste(keep, collapse = "|"), colnames(x))
    posNew <- rev(c(colnames(x)[posNow],colnames(x)[!posNow]))
    x <- x[,posNew]
  }
  #input is a feature matrix, features are in columns
  if (method == "binary") {
    #use binary similarity if input is a binary matrix,
    #maybe also usefull is the input is a sparse matrix
    simiMat <- 1 - as.matrix(dist(t(x), method = "binary"))
  } else if (method == "pearson") {
    #otherwise, using pearson correlation
    simiMat <- cor(x)
  } else if (method == "euclidean") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "euclidean"))
  } else if (method == "cosine") {
    # cosine similarity maybe prefered for sparse matrix
    cosineSimi <- function(x){
      x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    simiMat <- cosineSimi(t(x))
  } else if (method == "canberra") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "canberra"))/nrow(x)
  }

  #generate reduced matrix
  simiMat.ori <- simiMat
  simiMat[upper.tri(simiMat)] <- 0
  diag(simiMat) <- 0
  x.re <- x[,!apply(simiMat, 2, function(n) any(abs(n) >= cutoff))]

  if (record) {
    #a matrix keeping track of the removed features
    mapReduce <- simiMat.ori
    diag(mapReduce) <- 0
    mapList <- lapply(colnames(x.re), function(i) colnames(mapReduce)[mapReduce[i,]>=cutoff])
    names(mapList) <- colnames(x.re)
  } else mapList = NULL

  return(list(reduced = x.re,
              mapReduce = mapList))
}

# Run glmnet with gaussian familiy outcome
#' @export
runGlm <- function(X, y, method = "lasso", repeats=20, folds = 3, testRatio = NULL, lambda = "lambda.1se") {
  modelList <- list()
  lambdaList <- c()
  r2Train <- c()
  r2Test <- c()
  coefMat  <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)

  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }

  for (i in seq(repeats)) {
    if (!is.null(testRatio)) {
      testIdx <- sample(seq_along(y), length(y)*testRatio)
      X.test <- X[testIdx,]
      X.train <- X[-testIdx,]
      y.test <- y[testIdx]
      y.train <- y[-testIdx]
    } else {
      X.train <- X
      y.train <- y
    }


    #train model
    res <- cv.glmnet(X.train,y.train, type.measure = "mse",
                     nfolds = folds, alpha = alpha, standardize = FALSE,
                     intercept = TRUE, family = "gaussian")
    lambdaList <- c(lambdaList, res[[lambda]])

    #calculate variance explained for training model
    y.pred <- glmnet:::predict.cv.glmnet(res, s = lambda, newx = X.train)[,1]
    varExp <- cor(y.train,y.pred)^2
    r2Train <- c(r2Train, varExp)

    modelList[[i]] <- res
    coefMat[,i] <- coef(res, s = lambda)[-1]

    #test model if testRatio is speficied
    if(!is.null(testRatio)) {
      y.pred <- glmnet:::predict.cv.glmnet(res, s = lambda, newx = X.test)
      varExp <- cor(as.vector(y.test),as.vector(y.pred))^2
      r2Test <- c(r2Test, varExp)
    }
  }
  list(modelList = modelList, lambdaList = lambdaList, r2Train = r2Train, coefMat = coefMat,
       r2Test = r2Test)
}

#Get number of selected features from a LASSO model
#' @export
nFeature <- function(lassoRes) {
  coefMat <- lassoRes$coefMat
  return(colSums(coefMat != 0))
}

#modified version of plotTopWeights from MOFA package, take table as input instead of MOFA object
#' @export
plotTopWeights.m <- function (W, view, factor, nfeatures = 10, abs = TRUE, scale = TRUE,
          sign = "both") {
  facName <- factor
  viewName <- view
  W <- dplyr::filter(W, factor == facName, view == viewName)
  if (scale)
    W$value <- W$value/max(abs(W$value))
  W <- W[W$value != 0, ]
  W$sign <- ifelse(W$value > 0, "+", "-")
  if (sign == "positive") {
    W <- W[W$value > 0, ]
  }
  else if (sign == "negative") {
    W <- W[W$value < 0, ]
  }
  if (abs)
    W$value <- abs(W$value)
  W <- W[with(W, order(-abs(value))), ]
  if (nfeatures > 0)
    features <- head(W$feature, nfeatures)
  W <- W[W$feature %in% features, ]
  W <- W[with(W, order(-value, decreasing = TRUE)), ]

  #iatalicize genes
  W$feature <- italicizeGene(W$feature)
  W$feature <- factor(W$feature, levels = W$feature)
  W$start <- 0
  p <- ggplot(W, aes_string(x = "feature", y = "value")) +
    geom_point(size = 2) + geom_segment(aes_string(yend = "start", xend = "feature"), size = 0.75) +
    scale_colour_gradient(low = "grey", high = "black") +
    scale_x_discrete(labels = parse(text = levels(W$feature))) +
    coord_flip() + theme_full +
    theme(axis.title.x = element_text(size = rel(1.0), color = "black"), axis.title.y = element_blank(),
          axis.text.y = element_text(size = rel(1.0), hjust = 1, color = "black"),
          axis.text.x = element_text(size = rel(1.0), color = "black"),
          axis.ticks.y = element_blank(), axis.ticks.x = element_line(),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = rel(1.0), color = "black"),
          legend.key = element_rect(fill = "transparent"), panel.background = element_blank())
  if (sign == "negative")
    p <- p + scale_x_discrete(position = "top")
  if (abs)
    p <- p + ylim(0, max(W$value) + 0.1) + geom_text(label = W$sign,
                                                     y = max(W$value) + 0.1, size = 8)
  if (abs & scale)
    p <- p + ylab(paste("Absolute loading on ", factor))
  else if (abs & !scale)
    p <- p + ylab(paste("Absolute loading on ", factor))
  else if (!abs & scale)
    p <- p + ylab(paste("Loading on factor", factor))
  else p <- p + ylab(paste("Loading on factor", factor))
  return(p)
}

#Function for volcano plot
#' @export
plotVolcano <- function(pTab, fdrCut = 0.05, posCol = "red", negCol = "blue",
                        x_lab = "dm", plotTitle = "",ifLabel = FALSE) {
  plotTab <- pTab %>% mutate(ifSig = ifelse(p.adj > fdrCut, "n.s.",
                                            ifelse(df > 0, "up","down"))) %>%
    mutate(ifSig = factor(ifSig, levels = c("up","down","n.s."))) %>%
    mutate(gene = italicizeGene(gene))
  pCut <- -log10((dplyr::filter(plotTab, ifSig != "n.s.") %>% arrange(desc(p)))$p[1])
  g <- ggplot(plotTab, aes(x=df, y=-log10(p))) +
    geom_point(shape = 21, aes(fill = ifSig),size=3) +
    geom_hline(yintercept = pCut, linetype = "dashed") +
    annotate("text", x = -Inf, y = pCut, label = paste0(fdrCut*100,"% FDR"),
             size = 5, vjust = -1.2, hjust=-0.1) +
    scale_fill_manual(values = c(n.s. = "grey70",
                                  up = posCol, down = negCol)) +
    theme_full +
    theme(legend.title = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 15)) +
    ylab(expression(-log[10]*'('*italic(P)~value*')')) +
    xlab(x_lab) + ggtitle(plotTitle)

  if (ifLabel)
    g <- g + ggrepel::geom_text_repel(data = dplyr::filter(plotTab, ifSig != "n.s."),
                                      aes(x= df, y = -log10(p), label = gene),
                                      size=5, force = 2, parse = TRUE)

  return(g)
}

#Function to load gene signature collections
#' @export
loadGSC <- function(setDir) {
  allLine <- readLines(setDir)
  gsc <- list()
  lapply(allLine, function(eachLine) {
    sp <- strsplit(eachLine, "\t")[[1]]
    gsc[[sp[1]]] <<- sp[3:length(sp)]
  }) %>% invisible()
  return(list(gsc=gsc))
}


#Wrapper function to run gene enrichment analysis using Camera from Limma package
#' @export
runCamera <- function(exprMat, design, gmtFile, id = NULL,
                      contrast = ncol(design),  method = "camera", pCut = 0.05, direction = "both",
                      ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE,
                      setToHighlight = c(), setMap = NULL) {
  scaleFUN <- function(x) sprintf("%.1f", x)

  #prepare indices
  if (is.null(id)) id <- rownames(exprMat)

  if (is.character(gmtFile)) {
    idx <- limma::ids2indices(loadGSC(gmtFile)$gsc, id)
  } else {
    idx <- limma::ids2indices(gmtFile,id)
  }

  #run camera for fry
  if (method == "camera") {
    res <- limma::camera(exprMat, idx, design, contrast)
  } else if (method == "fry") {
    res <- limma::fry(exprMat, idx, design, contrast)
  }

  #plot enrichment results as bar plot

  plotTab <- res %>% rownames_to_column("Name")

  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
  if(!is.null(setMap)) {
    plotTab <- mutate(plotTab, newName = setMap[match(Name, setMap$setName),]$pathwayName) %>%
      mutate(Name = ifelse(is.na(newName), Name, newName))
  }

  plotTab <- plotTab %>%
    mutate(Direction= factor(Direction, levels =c("Down","Up"))) %>%
    arrange(desc(Direction),desc(PValue)) %>%
    mutate(Name = factor(Name, levels = Name))

  if (direction != "both") {
    plotTab <- filter(plotTab, Direction == direction)
  }

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  colAxis <- ifelse(plotTab$Name %in% setToHighlight, "red", "black")

  #for higlighting sets
  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = res, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue), fill = Direction)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5) +
      scale_fill_manual(values = c(Up = colList[1], Down = colList[2])) +
      coord_flip() + xlab("") +
      scale_y_continuous(labels=scaleFUN) +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_full +
      theme(axis.text.y = element_text(size= 12, color = colAxis),
            axis.text.x = element_text(size=12),
            plot.title = element_text(size=14, face = "bold"))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.1))
    } else {
      p <- p + theme(legend.position = "right")
    }

    return(list(enrichTab = res, enrichPlot = p))
  }
}


#function to plot heatmap for each gene set
#' @export
plotSetHeatmap <- function(geneSigTab, geneSet, setName, exprMat, colAnno, scale = TRUE, italicize = TRUE) {

  if (length(geneSet) == 1) {#whether it's a list of direction of genes
    geneList <- loadGSC(geneSet)[["gsc"]][[setName]]
  } else {
    geneList <- geneSet
  }

  sigGene <- dplyr::filter(geneSigTab, symbol %in% geneList) %>%
    arrange(desc(coef))
  colAnno <- colAnno[order(colAnno[,1]),,drop = FALSE]
  plotMat <- exprMat[sigGene[["id"]],rownames(colAnno)]

  if (scale) {
    #calculate z-score and sensor
    plotMat <- t(scale(t(plotMat)))
    plotMat[plotMat >= 4] <- 4
    plotMat[plotMat <= -4] <- -4
  }

  rowLabs <- sigGene$symbol

  if (italicize) {
    rowLabs <- sapply(rowLabs, function(n) {
      sprintf("italic('%s')",n)
    })
  }

  pheatmap(plotMat, color = colorRampPalette(c(colList[2],"white",colList[1]))(100),
           cluster_cols = FALSE, cluster_rows = FALSE,
           annotation_col = colAnno, labels_row = parse(text = rowLabs),
           show_colnames = FALSE, fontsize_row = 8, breaks = seq(-5,5, length.out = 101), treeheight_row = 0,
           border_color = NA, main = setName)

}

#function to plot heatmap for each gene set (using complex heatmap)
#' @export
plotSetHeatmapComplex <- function(geneSigTab, geneSet, setName, exprMat, colAnno, scale = TRUE,
                           rowAnno = NULL, annoCol = NULL, highLight = NULL, italicize = TRUE, plotTitle = NULL) {

  if (length(geneSet) == 1) {#whether it's a list of direction of genes
    geneList <- loadGSC(geneSet)[["gsc"]][[setName]]
  } else {
    geneList <- geneSet
  }
   
  if (is.null(plotTitle)) {
    plotTitle <- setName
  }
  
  sigGene <- dplyr::filter(geneSigTab, symbol %in% geneList) %>%
    arrange(desc(coef))

  colAnno <- colAnno[order(colAnno[,1]),,drop = FALSE]
  #colAnno <- colAnno[,rev(colnames(colAnno)),drop=FALSE]
  plotMat <- exprMat[sigGene[["id"]],rownames(colAnno)]

  if (scale) {
    #calculate z-score and sensor
    plotMat <- t(scale(t(plotMat)))
    plotMat[plotMat >= 4] <- 4
    plotMat[plotMat <= -4] <- -4
  }

  rowLabs <- sigGene$symbol

  if (italicize) {
    labFont <- 3
  } else labFont <- 1

  haCol <- ComplexHeatmap::HeatmapAnnotation(df = colAnno, col=annoCol, which = "column",annotation_name_gp = gpar(fontface = "bold"),
                                             simple_anno_size_adjust = TRUE)

  if (!is.null(rowAnno)) {
    haRow <- ComplexHeatmap::HeatmapAnnotation(df = rowAnno[rownames(plotMat),,drop=FALSE], col=annoCol, which = "row", annotation_name_gp = gpar(fontface = "bold"))
  } else haRow <- NULL

  labelCol <- rep("black",nrow(plotMat))

  if (!is.null(highLight)) {
    labelCol[rownames(plotMat) %in% highLight] <-"red"
  }

  ComplexHeatmap::Heatmap(plotMat, col = colorRampPalette(c(colList[2],"white",colList[1]))(100),name = "z-score",
                          top_annotation = haCol, left_annotation = haRow, show_column_names = FALSE,
                          cluster_columns  = FALSE, cluster_rows = FALSE,
                          row_names_gp = gpar(col = labelCol, fontface = labFont),
                          row_labels = rowLabs,
                          column_title = plotTitle, column_title_gp = gpar(cex= 1.5, fontface = "bold")
  )
}

#function to run FGSA
#' @export
runFGSEA <- function(limmaRes, signatures, stat = "t", name = "symbol",minSize=15,
                     maxSize =1000,nperm = 10000, pcut = 0.1, ifFDR = TRUE, showFDR = FALSE) {
  geneRanks <- limmaRes[[stat]]
  names(geneRanks) <- limmaRes[[name]]
  fgseaRes <- fgsea::fgsea(pathways = sigList,
                    stats = geneRanks,
                    minSize=minSize,
                    maxSize=maxSize,
                    nperm=nperm)
  if(ifFDR) {
    plotSets <- dplyr::filter(fgseaRes, padj <= pcut)$pathway
  } else {
    plotSets <- dplyr::filter(fgseaRes, pval <= pcut)$pathway
  }

  #plot gsea plots
  plots <- lapply(plotSets, function(setName) {
    pval <- formatNum(dplyr::filter(fgseaRes, pathway == setName)$pval, digits = 1)
    FDR <- formatNum(dplyr::filter(fgseaRes, pathway == setName)$padj, digits = 1)
    nes <- dplyr::filter(fgseaRes, pathway == setName)$NES
    NES <- format(nes, digits = 1,nsmall = 1)
    #showValue <- ifelse(showFDR, sprintf("italic(P)~'value = %s\nFDR = %s\nNormalized ES=%s'", pval, FDR, NES),
    #                    sprintf("italic(P)~'value = %s\nNormalized ES=%s'", pval, NES))
    pp <- fgsea::plotEnrichment(signatures[[setName]],
                         geneRanks)
    pp <- pp + ggtitle(setName) +
      ylab("Enrichment score (ES)") + xlab("Rank") +
      theme(plot.title = element_text(face = "bold", size = 15),
            plot.margin = margin(0.2,0.2,0.2,0.2,"cm"),
            axis.text = element_text(size=12),
            axis.title = element_text(size=14)
            )
    if (nes >0) {
      showValue <- sprintf("atop('          '~italic(P)~'=%s','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = Inf, y = Inf,
                          label = showValue, parse= TRUE,
                          hjust=1, vjust =1.3, size = 5)
    } else {
      showValue <- sprintf("atop(italic(P)~'=%s            ','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = 0, y = -Inf,
                          label = showValue, parse=TRUE,
                          hjust=0, vjust =-0.5, size = 5)
    }
    pp
  })
  names(plots) <- plotSets
  return(list(table = fgseaRes, plots = plots))
}


#Function to get reduced dimension table
#' @export
getReducedDimTab <- function(sceObj, dimName = "UMAP", annoMarker = NULL) {
  dim1 <- reducedDims(sceObj)[[dimName]][,1]
  sceRe <- sceObj[,!is.na(dim1)]
  reTab <- reducedDims(sceRe)[[dimName]] %>% data.frame()
  colnames(reTab) <- c("x","y")

  if (!is.null(annoMarker)) {
    markerExpr <- assay(sceObj)[annoMarker,!is.na(dim1)]
    markerExpr <- scaleExprs(markerExpr)
    reTab[[annoMarker]] <- markerExpr
  }

  reTab <- reTab %>% cbind(colData(sceRe) %>% data.frame())
}


#Function for subsetting single cell dataset
#' @export
subsetSCE <- function(sceObj, condiList = NULL, patList = NULL, cluster_id = NULL, k=NULL) {
  #subset object itself
  if (is.null(condiList)) {
    condiList <- unique(sceObj$condition)
  }

  if (is.null(patList)) {
    patList <- unique(sceObj$patient_id)
  }

  if (!is.null(condiList) | !is.null(patList)) {
    sceSub <- sceObj[,sceObj$condition %in% condiList & sceObj$patient_id %in% patList]

    #change experimental info
    expInfo <- metadata(sceSub)$experiment_info %>%
      dplyr::filter(condition %in% condiList, patient_id %in% patList) %>%
      mutate_if(is.factor, droplevels) %>%
      mutate(condition = factor(condition, levels= condiList))
    metadata(sceSub)$experiment_info <- expInfo
  }


  # subset for a certain cluster
  if (!is.null(cluster_id) & !is.null(k)) {
    clusterTab <- metadata(sceSub)$cluster_codes
    clusterTab <- clusterTab[clusterTab[[k]] %in% cluster_id,]
    clusterTab <- clusterTab %>% mutate_if(is.factor, droplevels)
    sceSub <- sceSub[,sceSub$cluster_id %in% clusterTab$som100]
    metadata(sceSub)$cluster_codes <- clusterTab
  }

  #droplevels of columns
  for (allCol in colnames(colData(sceSub))) {
    if (is.factor(sceSub[[allCol]])) {
      sceSub[[allCol]] <- droplevels(sceSub[[allCol]])
    }
  }

  return(sceSub)
}

# Function for scaling marker intensity from 0 to 1
#' @export
scaleExprs <- function (x, margin = 2, q = 0.01)
{
  if (!is(x, "matrix"))
    x <- as.matrix(x)
  qs <- c(rowQuantiles, colQuantiles)[[margin]]
  qs <- qs(x, probs = c(q, 1 - q))
  qs <- matrix(qs, ncol = 2)
  x <- switch(margin, `1` = (x - qs[, 1])/(qs[, 2] - qs[, 1]),
              `2` = t((t(x) - qs[, 1])/(qs[, 2] - qs[, 1])))
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1
  return(x)
}

# Function to get median marker expression
#' @export
getMedian <- function(sceObj, useCluster) {
  clustTab <- metadata(sceObj)$cluster_code
  allCluster <- clustTab[match(sceObj$cluster_id, clustTab$som100),][[useCluster]]
  exprMat <- assays(sceObj)[["exprs"]]

  sumTab <- lapply(seq(nrow(exprMat)), function(i) {
    data.frame(expr = exprMat[i,],
               sample_id = sceObj$sample_id,
               cluster = allCluster) %>%
      group_by(sample_id, cluster) %>%
      summarise(medVal = median(expr, na.rm=TRUE),
                .groups = "drop") %>%
      mutate(antigen = rownames(exprMat)[i])
  }) %>% bind_rows()
  return(sumTab)
}



# A function to make a list of plots share legends in cowplot
#' @export

shareLegend <- function(plotList, position = "left", ratio = c(1,0.1), ncol=2,
                        rel_widths = NULL) {
  #get legend
  legend <- get_legend(plotList[[1]])

  if (is.null(rel_widths)) rel_widths <- rep(1, length(plotList))

  #main plot
  plotList <- lapply(plotList, function(p) p+theme(legend.position = "none"))
  mainPlot <- plot_grid(plotlist = plotList, ncol=ncol,
                        rel_widths = rel_widths)

  #plot the whole plot
  if (position == "bottom") {
    plot_grid(mainPlot, legend, ncol=1, rel_heights  = ratio)
  } else if (position == "left") {
    plot_grid(mainPlot, legend, ncol=2, rel_widths  = ratio)

  }
}


# Function to calculate cell population from an SCE object
#' @export

getPop <- function(sceObj, k) {

  clusterTab <- metadata(sceObj)$cluster_codes

  #get sample annotation an cluster id
  sampleTab <- colData(sceObj)[,c("sample_id","cluster_id")] %>% data.frame() %>%
    mutate(cluster = clusterTab[match(cluster_id, clusterTab$som100),][[k]])

  #calculate cell numbers
  cellTab <- group_by(sampleTab, sample_id) %>%
    summarise(n_cell = length(sample_id)) %>% ungroup()

  popTab <- sampleTab %>%
    group_by(sample_id, cluster) %>% summarise(nClust = sum(length(sample_id))) %>%
    ungroup() %>%
    pivot_wider(names_from = "cluster", values_from = "nClust") %>%
    pivot_longer(-sample_id, names_to = "cluster", values_to = "nClust") %>%
    mutate(nClust = replace_na(nClust,0)) %>%
    left_join(cellTab, by = "sample_id") %>%
    mutate(popSize = 100*nClust/n_cell)

  return(popTab)
}

#plot reduced dimensions
#' @export

plotDM <- function(dmTab, sceObj, colorBy = "cluster", facetBy = NULL, facetCol = 4,  x="DC1", y="DC2",
                   fast = FALSE, colPal = "Zissou1", xLab = "", yLab = "", plotTitle = TRUE,pointsize=2) {
  a_pal = hcl.colors(10, colPal)
  plotTab <- dmTab

  if (!colorBy %in% colnames(plotTab)) {
    plotTab[["colorBy"]] <- scaleExprs(assay(sceObj)[colorBy,],2)
    colorScale <- scale_colour_gradientn("scaled intensity", colors = a_pal)
  } else {
    plotTab[["colorBy"]] <- dmTab[[colorBy]]
    colorScale <- scale_color_manual(values = colList, name = "")
  }

  if (is.null(facetBy)) {
    fct <- NULL
  } else if (facetBy == "patient_id") {
    fct <- facet_wrap(~ CLLPD + patient_id, ncol = facetCol)
  } else if (facetBy == "CLL-PD") {
    plotTab <- mutate(plotTab, CLLPD = paste0(CLLPD, " CLL-PD"))
    fct <- facet_wrap(~ CLLPD, ncol = facetCol)
  }

  if (fast) {
    plotPoint <- geom_scattermore(pointsize = 1) #faster scatter plot using scatter more
  } else {
    plotPoint <- geom_scattermore(pointsize = pointsize, alpha=0.8, pixels = c(1000,1000)) #looks better
  }
  p <- ggplot(plotTab, aes_string(x, y, col = "colorBy")) +
    plotPoint +
    colorScale +
    guides(col = guide_legend(override.aes = list(alpha = 1, size = 5)))+
    ggtitle(ifelse(plotTitle, colorBy, "")) +
    theme_void() + xlab(xlab) + ylab(ylab) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size=22, face = "bold"),
          strip.background = element_rect(fill = "grey80"),
          strip.text = element_text(face = "plain", size=20),
          legend.position = "bottom", legend.text = element_text(size=18), legend.title = element_text(size=18)) +
    fct
  return(p)
}



# Function for performing differential marker expression test
#' @export

diffCondiDE <- function(sceTest, compareList, clusterUse) {

  resDE_condition <-lapply(compareList, function(eachPair) {
    #subset
    usePat <- unique(sceTest[,sceTest$condition == eachPair[2]]$patient_id) #only use patients with the treatment
    sceSub <- subsetSCE(sceTest, condiList = eachPair, patList = usePat)

    ei <- metadata(sceSub)$experiment_info

    designMat <- createDesignMatrix(ei, cols_design = c("condition", "patient_id"))
    contrast <- createContrast(c(0, 1, rep(0, ncol(designMat)-2)))

    ds_res <- diffcyt(sceSub,
                      design = designMat, contrast = contrast,
                      analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                      clustering_to_use = clusterUse, verbose = FALSE,
                      transform = FALSE)


    #get results
    res <- rowData(ds_res$res) %>%
      data.frame() %>% mutate(control = eachPair[1],treat = eachPair[2],
                              p_val = ifelse(is.na(p_val),1, p_val),
                              p_adj = ifelse(is.na(p_adj), 1, p_adj))

    res
  }) %>% bind_rows()

  return(resDE_condition)
}

#Function for CyTOF volcano plots
#' @export
plotCyToVolcano <- function(pTab, fdrCut = 0.10, labSize =5, x_lab = NULL, plotTitle = "",
                            useFdr = TRUE, globalPcut = NULL) {

  pathMap <- c("c-Myc" = "MYC", "GLUT1" = "MYC", "CDK4" = "MYC",
               "P-S6K" = "mTOR","P-4E-BP1" = "mTOR", "P-S6" = "mTOR",
               "P-PLC-gamma 2" = "BTK", "P-ZAP70/Syk" = "BTK" , "P-BTK" = "BTK")

  pathColor <- c("MYC" = colList[3], "mTOR" = colList[5], "BTK" = colList[4], "other" = "black")


  if (!useFdr) pTab <- mutate(pTab, p_adj = p_val)
  if (is.null(x_lab))  x_lab <- expression(log[2]*'(fold change)')

  plotTab <- pTab %>% mutate(ifSig = ifelse(p_adj > fdrCut, "n.s.",
                                            ifelse(logFC > 0, "up","down"))) %>%
    mutate(ifSig = factor(ifSig, levels = c("up","down","n.s."))) %>%
    mutate(pathway = pathMap[as.character(marker_id)]) %>%
    mutate(pathway = ifelse(is.na(pathway),"other",pathway))

  if (is.null(globalPcut)) {
    pCut <- -log10((dplyr::filter(plotTab, ifSig != "n.s.") %>% arrange(desc(p_val)))$p_val[1])
  } else {
    pCut <- -log10(globalPcut)
  }

  g <- ggplot(plotTab, aes(x=logFC, y=-log10(p_val))) +
    geom_point(shape = 21, aes(fill = ifSig),size=3) +
    geom_hline(yintercept = pCut, linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha=0.5) +
    scale_fill_manual(values = c(n.s. = "grey80",
                                 up = colList[1], down = colList[2]),name = "") +
    ggrepel::geom_label_repel(aes(label = marker_id, col = pathway), size=labSize, force = 5, max.overlaps = 100) +
    scale_color_manual(values = pathColor, guide = FALSE) +
    theme_full +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 20),
          legend.margin=margin(0,0,0,0),
          axis.title = element_text(size=20),
          axis.text = element_text(size=18),
          plot.title = element_text(size=20, face= "bold")) +
    ylab(expression(-log[10]*'('*italic(P)~value*')')) +
    xlab(x_lab) + ggtitle(plotTitle)

  if (useFdr) {
    g <- g +  annotate("text", x = -Inf, y = pCut, label = paste0(fdrCut*100,"% FDR"),
                       size = 6, vjust = -1.2, hjust=-0.1)
  } else {
    g <- g +  annotate("text", x = -Inf, y = pCut, label = paste("italic(P)~'=",formatNum(fdrCut, digits = 1),"'"),
                       size = 6, vjust = -1.2, hjust=-0.1, parse = TRUE)
  }
  return(g)
}



# Funciton for differential marker expression related to CLL-PD
#' @export
diffPdDE <- function(sceTest, clusterUse) {

  condiList <- levels(sceTest$condition)

  resDE <- lapply(condiList, function(eachCondi) {
    #subset
    sceSub <- subsetSCE(sceTest, eachCondi)

    #define experiment design
    ei <- metadata(sceSub)$experiment_info
    ei$CLLPD <- relevel(ei$CLLPD, ref = "low")
    designMat <- createDesignMatrix(ei, cols_design = c("CLLPD"))
    contrast <- createContrast(c(0, 1))

    #test using diffcyte
    ds_res <- diffcyt(sceSub,
                      design = designMat, contrast = contrast,
                      analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                      clustering_to_use = clusterUse, verbose = FALSE,
                      transform = FALSE)

    #get results
    res <- rowData(ds_res$res) %>%
      data.frame() %>% mutate(condition = eachCondi)
    res
  }) %>% bind_rows()

  return(resDE)
}

# Funciton for plotting differential  marker expression
#' @export
plotDiffExpr <- function(exprTab, markers, useCondi, resTab,
                         clusterPlot = "CLL",colorBy = "CLLPD", ncol=3, labelSample = NA, scale="free_y") {

  plotTab <- exprTab %>% dplyr::filter(condition == useCondi, cluster == clusterPlot, antigen %in% markers) %>%
    mutate(condition = as.character(condition)) %>%
    mutate(antigen = factor(antigen, levels = markers)) %>%
    select(patient_id, IGHV, CLLPD, condition, exprVal, antigen) %>%
    mutate(patient_id = as.character(patient_id)) %>%
    mutate(sampleLabel = ifelse(patient_id %in% labelSample, patient_id, ""),
           CLLPD = paste0(CLLPD," CLL-PD"))

  rangeTab <- group_by(plotTab, antigen) %>% summarise(yMax = max(exprVal), yMin = min(exprVal)) %>%
    mutate(yRange = yMax-yMin)

  pTab <- dplyr::filter(resTab, marker_id %in% markers, condition == useCondi) %>%
    mutate(x1 = 1, x2 = 2, xp = 1.5) %>%
    left_join(rangeTab, by = c(marker_id="antigen")) %>%
    mutate(annoP = ifelse(p_val < 1e-12, paste("italic(P)~'<",formatNum(1e-12, digits = 1),"'"),
                          paste("italic(P)~'=",formatNum(p_val, digits = 1),"'"))) %>%
    mutate(antigen = marker_id)

  scaleFUN <- function(x) sprintf("%.1f", x)

  p<-ggplot(plotTab, aes(x=CLLPD, y=exprVal)) +
    geom_boxplot(width=0.4, outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(aes_string(fill = colorBy), shape =21, cex=3, width = 0.1) +
    #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
    #             geom = "crossbar", width = 0.5, alpha=0.5) +
    scale_y_continuous(labels=scaleFUN)+
    ggrepel::geom_text_repel(aes(label = sampleLabel), size=4) +
    geom_blank(data = pTab, aes(x=1, y=yMax + yRange*0.6)) +
    geom_segment(data = pTab, aes(x=x1, xend =x2, y = yMax+yRange*0.3, yend = yMax+yRange*0.3)) +
    geom_text(data = pTab, aes(label = annoP, x=xp, y=yMax+yRange*0.3), vjust=-1, parse=TRUE) +
    scale_fill_manual(values = colorMap, name =colorBy) +
    ylab(paste0("Median intensity")) + xlab("")  +
    facet_wrap(~antigen, scale = scale) +
    theme_half  + theme(legend.position = "bottom",
                        strip.text = element_text(size=15, face= "bold"),
                        strip.background = element_rect(fill = "white",color = "white"),
                        axis.text.x = element_text(size=14),
                        legend.text = element_text(size=13),
                        legend.title = element_text(size=18),
                        legend.box.margin=margin(-20,0,0,0),
                        legend.margin=margin(0,0,0,0))

  if (colorBy == "CLLPD") {p <- p + theme(legend.position = "none")}

  return(p)
}
