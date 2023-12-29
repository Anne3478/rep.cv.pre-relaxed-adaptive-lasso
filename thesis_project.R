########################################################################
##################### Thesis analysis ##################################
########################################################################

#setwd("~/Master Statistics and Data Science/Thesis project/Data_Analysis")

# packages for reading and cleaning datasts:
require(readr)
require(haven)
require(readxl)
require(dplyr)

# packages used in rep_cv_pre function
require(Formula)
require(pROC)
library(pre)
require(glmnet)

# packages used for model comparison
require(colorspace)
require(ggplot2)
require(lmerTest)

# packages to print df in LaTeX format
require(knitr)
require(kableExtra)

#########################################################
############## Function for Latex formatting ############
#########################################################
footnote = "\\footnotesize \\\\ \\textit{Note.} Significance level: < 0.001 ‘***’ >0.001 ‘**’ >0.01 ‘*’ >0.05 ‘.’ >0.1 ‘ ’"

latex <- function(df, label =NULL, caption = NULL, footnote = NA, model = NA) {
  if(!is.na(model)) {
    for(i in 1:nrow(df)) {
      if(round(as.numeric(df$`Pr(>|t|)`[i]), 3) == 0 ) {
        df$`Pr(>|t|)`[i] = ">0.001 ***"
      } else if (df$`Pr(>|t|)`[i] < 0.01) {
        df$`Pr(>|t|)`[i] = paste(format(round(as.numeric(df$`Pr(>|t|)`[i]), 3),
                                        nsmall = 3), "**")
      } else if(df$`Pr(>|t|)`[i] < 0.05) {
        df$`Pr(>|t|)`[i] = paste(format(round(as.numeric(df$`Pr(>|t|)`[i]), 3),
                                        nsmall = 3), "*")
      } else if(df$`Pr(>|t|)`[i] < 0.1) {
        df$`Pr(>|t|)`[i] = paste(format(round(as.numeric(df$`Pr(>|t|)`[i]), 3),
                                        nsmall = 3), "$.$")
      } else {
        df$`Pr(>|t|)`[i] = format(round(as.numeric(df$`Pr(>|t|)`[i]), 3),
                                  nsmall = 3)
      }
    }
    if(model == "lm") {
      names(df) <- c("{Estimate}",
                     "{Std. Error}",
                     "{t-value}",
                     "{$P(> |t|)$}")
    } else {
      names(df) <- c("{Estimate}",
                     "{Std. Error}",
                     "{df}",
                     "{t-value}",
                     "{$P(> |t|)$}")
    }
    
    latex_table <- kable_styling(add_footnote(kable(df, format = "latex", 
                                                    booktabs = TRUE, label = label, 
                                                    align = "S", escape=FALSE,
                                                    caption = caption, row.names = T, 
                                                    digits = 3L),
                                              footnote, threeparttable = T,
                                              escape = F, notation = "none"),
                                 latex_options = "hold_position")
  } else {
    latex_table <- kable(df, format = "latex", 
                         booktabs = TRUE, label = label, 
                         align = "S", escape=FALSE,
                         caption = caption, row.names = T, 
                         digits = 3L)
  }
  return(latex_table)
}

##########################################################
######## repeated 10-fold cross-validation ###############
##########################################################

rep_cv_pre <- function(formula, data, k = 10, r = 10, CI_alpha = 0.05, verbose = F,
                       foldids = NULL, family = "gaussian", pclass = 0.5, 
                       alasso = F,...) {
  if (!(length(k) == 1L && k == as.integer(k))) {
    stop("Argument k should be a single positive integer.")
  }
  if (!(length(r) == 1L && r == as.integer(r))){
    stop("Argument r should be a single positive integer.")
  }
  if(!(length(CI_alpha) == 1L && is.numeric(CI_alpha))){
    stop("Argument CI_alpha should be a single positive numeric value.")
  } else if ((CI_alpha <= 0 || CI_alpha >= 1)) {
    stop("Argument CI_alpha should lie between 0 and 1.")
  }
  if(!is.null(foldids) && nrow(foldids) != nrow(data)) {
    stop("Foldids should have the same number of rows as the data")
  } else if(!is.null(foldids) && ncol(foldids) != r){
    stop("The number of columns of foldids should be equal to the number of repeats r")
  }
  
  if(is.null(foldids)) {
    foldids <- replicate(r,
                         sample(rep(1:k, length = nrow(data)),
                                replace = FALSE))
  }
  y_names <- attr(terms(Formula(formula), rhs = 0, data = data), 
                  "variables")
  y_names <- as.character(y_names[-1])
  x_names <- attr(terms(Formula(formula), lhs = 0, data = data), 
                  "term.labels")
  
  y_obs <- data[[y_names]]
  n <- nrow(data)
  x_ncol <- length(x_names)
  y_ncol <- ifelse(family != "multinomial", length(y_names), 
                   nlevels(data[[y_names]]) )
  var_y <- as.numeric(var(data[,y_names]))
  
  
  # Initialize variables  
  cvpreds.1se <- replicate(n = r * y_ncol, rep(NA, times = nrow(data)))
  cvpreds.min <- replicate(n = r * y_ncol, rep(NA, times = nrow(data)))
  if (y_ncol == 1){
    colnames(cvpreds.1se) <- paste("repeat", seq(1,r))
    colnames(cvpreds.min) <- paste("repeat", seq(1,r))
  } else {
    col_names <- paste("category", rep(levels(data[[y_names]]), r), "repeat", 
                       rep(1:r, each = y_ncol))
    colnames(cvpreds.1se) = col_names
    colnames(cvpreds.min) = col_names
  }
  
  if (family %in% c("gaussian", "poisson")) {
    accuracy.1se <- data.frame('repeat.' = 1:r,
                               'MSE' = rep(NA,r),
                               'MSE_se' = rep(NA, r),
                               'MAE' = rep(NA, r),
                               'MAE_se' = rep(NA, r),
                               'VAF' = rep(NA, r),
                               'VAF_adj' = rep(NA, r))
    
    accuracy.min <- data.frame('repeat.' = 1:r,
                               'MSE' = rep(NA,r),
                               'MSE_se' = rep(NA, r),
                               'MAE' = rep(NA, r),
                               'MAE_se' = rep(NA, r),
                               'VAF' = rep(NA, r),
                               'VAF_adj' = rep(NA, r))
    
  } else if (family == "binomial") {
    accuracy.1se <- data.frame('repeat.' = 1:r,
                               'SEL' = rep(NA,r),
                               'SEL_se' = rep(NA, r),
                               'AEL' = rep(NA, r),
                               'AEL_se' = rep(NA, r),
                               'VAF' = rep(NA, r),
                               'VAF_adj' = rep(NA, r),
                               'MCR' = rep(NA, r),
                               'PAC' = rep(NA, r),
                               'AUC' = rep(NA, r),
                               'Sensitivity' = rep(NA, r),
                               'Specificity' = rep(NA, r))
    
    accuracy.min <- data.frame('repeat.' = 1:r,
                               'SEL' = rep(NA,r),
                               'SEL_se' = rep(NA, r),
                               'AEL' = rep(NA, r),
                               'AEL_se' = rep(NA, r),
                               'VAF' = rep(NA, r),
                               'VAF_adj' = rep(NA, r),
                               'MCR' = rep(NA,r),
                               'PAC' = rep(NA, r),
                               'AUC' = rep(NA, r),
                               'Sensitivity' = rep(NA, r),
                               'Specificity' = rep(NA, r))
    
    classification_table.1se = list()
    classification_table.min = list()
    predicted.1se = predicted.min = list()
  } else if (family == "multinomial") {
    accuracy.1se <- data.frame('repeat.' = 1:r,
                               'MCR' = rep(NA,r),
                               'PAC' = rep(NA, r),
                               'AUC' = rep(NA, r))
    accuracy.min <- data.frame('repeat.' = 1:r,
                               'MCR' = rep(NA,r),
                               'PAC' = rep(NA, r),
                               'AUC' = rep(NA, r))
    accuracy_list.1se = accuracy_list.min = list()
    classification_table.1se = list()
    classification_table.min = list()
    predicted.1se = predicted.min = list()
  }
  
  
  lambda <- data.frame('fold' = rep(seq(from = 1, to = k), r),
                       'repeat.' = rep(1:r, each = k),
                       'lambda.1se' = rep(NA, k*r),
                       'lambda.min' = rep(NA, k*r))
  
  # Importance measures
  ## 1 se 
  df_imps.1se <- matrix(rep(NA, (x_ncol + 2) *(r*k)), 
                        ncol = x_ncol + 2, 
                        nrow = r*k)
  
  colnames(df_imps.1se) <- c("repeat.", "fold", x_names)
  df_imps.1se <- data.frame(df_imps.1se)
  df_imps.1se$repeat. <- rep(1:r, each = k)
  df_imps.1se$fold <- rep(seq(from = 1, to = k), r)
  
  imps.1se <-list()
  
  ## min rule
  
  df_imps.min <- matrix(rep(NA, (x_ncol + 2) *(r*k)), 
                        ncol = x_ncol + 2, 
                        nrow = r*k)
  
  colnames(df_imps.min) <- c("repeat.", "fold", x_names)
  df_imps.min <- data.frame(df_imps.min)
  df_imps.min$repeat. <- rep(1:r, each = k)
  df_imps.min$fold <- rep(seq(from = 1, to = k), r)
  
  imps.min <-list()
  
  
  
  sparsity <- data.frame('count_predictors.1se' = rep(NA, k*r), 
                         'count_predictors.min' = rep(NA, k*r), 
                         'count_rules.1se' = rep(NA, k*r),
                         'count_rules.min' = rep(NA, k*r),
                         'fold' = rep(seq(from = 1, to = k), r), 
                         'repeat.' = rep(1:r, each = k))
  predictors.1se <- c()
  predictors.min <- c()
  
  ## fit models
  for (rep_ in 1:r)  {
    
    for (fold in 1:k) {
      mod <- pre(formula = formula, data = data[foldids[, rep_] != fold, ],
                 family = family, ...)
      
      if(alasso == TRUE){
        x <- mod$modmat
        y <- mod$data[ , mod$y_names]
        ridge1_cv <- cv.glmnet(x, y,
                               type.measure = "mse",
                               nfold = 10,
                               alpha = 0)
        best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
        mod$glmnet.fit <- cv.glmnet(x, y, alpha = 1,
                                    penalty.factor = 1 / abs(best_ridge_coef))
      }
      
      lambda[lambda$fold == fold & lambda$repeat. == rep_,]$lambda.1se <- mod$glmnet.fit$lambda.1se
      lambda[lambda$fold == fold & lambda$repeat. == rep_,]$lambda.min <- mod$glmnet.fit$lambda.min
      
      imps.1se[[fold]] <- importance(mod, plot = FALSE)
      imps.min[[fold]] <- importance(mod, plot = FALSE,
                                     penalty.par.val = "lambda.min")
      
      
      if(family %in% c("gaussian", "poisson", "binomial")) {
        cvpreds.1se[foldids[, rep_] == fold, rep_] <- predict(mod, newdata = data[foldids[, rep_] == fold, ],
                                                              type = "response")
        cvpreds.min[foldids[, rep_] == fold, rep_] <- predict(mod, newdata = data[foldids[, rep_] == fold, ],
                                                              type = "response",
                                                              penalty.par.val = "lambda.min")
        df_imps.1se[df_imps.1se$fold == fold & df_imps.1se$repeat. == rep_, imps.1se[[fold]]$varimps$varname] <- imps.1se[[fold]]$varimps$imp
        df_imps.min[df_imps.1se$fold == fold & df_imps.min$repeat. == rep_, imps.min[[fold]]$varimps$varname] <- imps.min[[fold]]$varimps$imp
        
      } else if(family == "multinomial") {
        cvpreds.1se[foldids[, rep_] == fold, (rep_*y_ncol-y_ncol + 1):(rep_*y_ncol)] <- predict(mod, newdata = data[foldids[, rep_] == fold, ],
                                                                                                type = "response")
        cvpreds.min[foldids[, rep_] == fold, (rep_*y_ncol-y_ncol + 1):(rep_*y_ncol)] <- predict(mod, newdata = data[foldids[, rep_] == fold, ],
                                                                                                type = "response",
                                                                                                penalty.par.val = "lambda.min")
      }
      
      
      sparsity[sparsity$fold == fold & sparsity$repeat. == rep_, ]$count_predictors.1se <- nrow(imps.1se[[fold]]$varimps)
      sparsity[sparsity$fold == fold & sparsity$repeat. == rep_, ]$count_rules.1se <- nrow(imps.1se[[fold]]$baseimps)
      
      sparsity[sparsity$fold == fold & sparsity$repeat. == rep_, ]$count_predictors.min <- nrow(imps.min[[fold]]$varimps)
      sparsity[sparsity$fold == fold & sparsity$repeat. == rep_, ]$count_rules.min <- nrow(imps.min[[fold]]$baseimps)
      
      predictors.1se <- c(predictors.1se, imps.1se[[fold]]$varimps$varname)
      predictors.min <- c(predictors.min, imps.min[[fold]]$varimps$varname)
      
    }
    p.1se <- mean(sparsity[sparsity$repeat. == rep_,]$count_rules.1se)
    p.min <- mean(sparsity[sparsity$repeat. == rep_,]$count_rules.min)
    
    if (family %in% c("gaussian", "poisson")) {
      accuracy.1se[rep_, ]$repeat. <- rep_
      accuracy.1se[rep_, ]$MSE <- mean((y_obs - cvpreds.1se[,rep_])^2, na.rm = TRUE)
      accuracy.1se[rep_, ]$MSE_se <- sd((y_obs - cvpreds.1se[,rep_])^2, na.rm = TRUE)/sqrt(n)
      accuracy.1se[rep_, ]$MAE <- mean(abs(y_obs - cvpreds.1se[,rep_] ), na.rm = TRUE)
      accuracy.1se[rep_, ]$MAE_se <- sd(abs(y_obs - cvpreds.1se[,rep_]), na.rm = TRUE)/sqrt(n)
      R2 <- 1 - (accuracy.1se[rep_, ]$MSE/var_y)
      accuracy.1se[rep_, ]$VAF <- R2
      accuracy.1se[rep_, ]$VAF_adj <- 1 - ( (1-R2)*(n-1) / (n-p.1se-1) )
      
      accuracy.min[rep_, ]$repeat. <- rep_
      accuracy.min[rep_, ]$MSE <- mean((y_obs - cvpreds.min[,rep_])^2, na.rm = TRUE)
      accuracy.min[rep_, ]$MSE_se <- sd((y_obs - cvpreds.min[,rep_])^2, na.rm = TRUE)/sqrt(n)
      accuracy.min[rep_, ]$MAE <- mean(abs(y_obs - cvpreds.min[,rep_] ), na.rm = TRUE)
      accuracy.min[rep_, ]$MAE_se <- sd(abs(y_obs - cvpreds.min[,rep_]), na.rm = TRUE)/sqrt(n)
      R2 <- 1 - (accuracy.min[rep_, ]$MSE/var_y)
      accuracy.min[rep_, ]$VAF <- R2
      accuracy.min[rep_, ]$VAF_adj <- 1 - ( (1-R2)*(n-1) / (n-p.min-1) )
      
      
    } else if (family == "binomial") {
      observed <- data[[y_names]]
      y_obs <- as.numeric(observed) - 1
      accuracy.1se[rep_, ]$SEL <- mean((y_obs - cvpreds.1se[, rep_])^2, na.rm = TRUE)
      accuracy.1se[rep_, ]$SEL_se <- sd((y_obs - cvpreds.1se[, rep_])^2, na.rm = TRUE)/sqrt(n)
      accuracy.1se[rep_, ]$AEL <- mean(abs(y_obs - cvpreds.1se[, rep_] ), na.rm = TRUE)
      accuracy.1se[rep_, ]$AEL_se <- sd(abs(y_obs - cvpreds.1se[, rep_]), na.rm = TRUE)/sqrt(n)
      R2 <- 1 - (accuracy.1se[rep_, ]$SEL/var_y)
      accuracy.1se[rep_, ]$VAF <- R2
      accuracy.1se[rep_, ]$VAF_adj <- 1 - ( (1-R2)*(n-1) / (n-p.1se-1) )
      
      accuracy.min[rep_, ]$SEL <- mean((y_obs - cvpreds.min[, rep_])^2, na.rm = TRUE)
      accuracy.min[rep_, ]$SEL_se <- sd((y_obs - cvpreds.min[, rep_])^2, na.rm = TRUE)/sqrt(n)
      accuracy.min[rep_, ]$AEL <- mean(abs(y_obs - cvpreds.min[, rep_] ), na.rm = TRUE)
      accuracy.min[rep_, ]$AEL_se <- sd(abs(y_obs - cvpreds.min[, rep_]), na.rm = TRUE)/sqrt(n)
      R2 <- 1 - (accuracy.min[rep_, ]$SEL/var_y)
      accuracy.min[rep_, ]$VAF <- R2
      accuracy.min[rep_, ]$VAF_adj <- 1 - ( (1-R2)*(n-1) / (n-p.min-1) )
      
      predicted.1se[[rep_]] <- factor(cvpreds.1se[, rep_] > pclass)
      levels(predicted.1se[[rep_]]) <- levels(observed)
      classification_table.1se[[rep_]] <- prop.table(table(predicted.1se[[rep_]], observed))
      accuracy.1se[rep_, ]$MCR <- round((1 - sum(diag(classification_table.1se[[rep_]]))),4)
      accuracy.1se[rep_, ]$PAC <- 100 * (1 - accuracy.1se[rep_,]$MCR)
      auc <- suppressMessages(roc(as.numeric(predicted.1se[[rep_]]), as.numeric(observed)))
      accuracy.1se[rep_, ]$AUC <- auc$auc[1]
      accuracy.1se[rep_, ]$Sensitivity <- classification_table.1se[[rep_]][1] / sum(classification_table.1se[[rep_]][,1])
      accuracy.1se[rep_, ]$Specificity <- classification_table.1se[[rep_]][4] / sum(classification_table.1se[[rep_]][,2])
      
      
      predicted.min[[rep_]] <- factor(cvpreds.min[, rep_] > pclass)
      levels(predicted.min[[rep_]]) <- levels(observed)
      classification_table.min[[rep_]] <- prop.table(table(predicted.min[[rep_]], observed))
      accuracy.min[rep_, ]$MCR <- round((1 - sum(diag(classification_table.min[[rep_]]))),4)
      accuracy.min[rep_, ]$PAC <- 100 * (1 - accuracy.min[rep_,]$MCR)
      auc <- suppressMessages(roc(as.numeric(predicted.min[[rep_]]), as.numeric(observed)))
      accuracy.min[rep_, ]$AUC <- auc$auc[1]
      accuracy.min[rep_, ]$Sensitivity <- classification_table.min[[rep_]][1] / sum(classification_table.min[[rep_]][,1])
      accuracy.min[rep_, ]$Specificity <- classification_table.min[[rep_]][4] / sum(classification_table.min[[rep_]][,2])
      
    } else if (family == "multinomial") {
      observed <- data[[y_names]]
      y_obs <- model.matrix(~observed + 0)
      colnames(y_obs) <- levels(observed)
      
      accuracy_list.1se$SEL[[rep_]] <- data.frame(SEL = colMeans((y_obs - cvpreds.1se[,rep_])^2,
                                                                 na.rm = TRUE), 
                                                  se = apply((y_obs - cvpreds.1se[,rep_])^2, 2, sd, na.rm = TRUE)/sqrt(n))
      accuracy_list.1se$AEL <- data.frame(AEL = colMeans(abs(y_obs - cvpreds.1se[,rep_]), na.rm = TRUE),
                                          se = apply(abs(y_obs - cvpreds.1se[,rep_]), 2, sd, na.rm = TRUE)/sqrt(n))
      
      accuracy_list.min$SEL[[rep_]] <- data.frame(SEL = colMeans((y_obs - cvpreds.min[,rep_])^2,
                                                                 na.rm = TRUE), 
                                                  se = apply((y_obs - cvpreds.min[,rep_])^2, 2, sd, 
                                                             na.rm = TRUE)/sqrt(n))
      accuracy_list.min$AEL <- data.frame(AEL = colMeans(abs(y_obs - cvpreds.min[,rep_]), 
                                                         na.rm = TRUE),
                                          se = apply(abs(y_obs - cvpreds.min[,rep_]), 2, sd, 
                                                     na.rm = TRUE)/sqrt(n))
      
      predicted.1se[[rep_]] <- factor(apply(cvpreds.1se[,(rep_*y_ncol-y_ncol + 1):(rep_*y_ncol)], 1, 
                                            function(x) which(x == 
                                                                max(x))), levels = 1:y_ncol)
      levels(predicted.1se[[rep_]]) <- levels(observed)
      classification_table.1se[[rep_]] <- prop.table(table(predicted.1se[[rep_]], observed))
      accuracy.1se[rep_, ]$MCR <- round((1 - sum(diag(classification_table.1se[[rep_]]))),4)
      accuracy.1se[rep_, ]$PAC <- 100 * (1 - accuracy.1se[rep_,]$MCR)
      auc <- suppressMessages(multiclass.roc(as.numeric(predicted.1se[[rep_]]), 
                                             as.numeric(observed)))
      accuracy.1se[rep_, ]$AUC <- round(auc$auc[1], 4)
      
      predicted.min[[rep_]] <- factor(apply(cvpreds.min[,(rep_*y_ncol-y_ncol + 1):(rep_*y_ncol)], 1, 
                                            function(x) which(x == 
                                                                max(x))), levels = 1:y_ncol)
      levels(predicted.min[[rep_]]) <- levels(observed)
      classification_table.min[[rep_]] <- prop.table(table(predicted.min[[rep_]], observed))
      accuracy.min[rep_, ]$MCR <- round((1 - sum(diag(classification_table.min[[rep_]]))),4)
      accuracy.min[rep_, ]$PAC <- 100 * (1 - accuracy.min[rep_,]$MCR)
      auc <- suppressMessages(multiclass.roc(as.numeric(predicted.min[[rep_]]), 
                                             as.numeric(observed)))
      accuracy.min[rep_, ]$AUC <- round(auc$auc[1], 4)
    }
    if (verbose == TRUE) {
      print(paste("Repetition", rep_, "of", r))
      if(rep_ == r) {
        "Done! \n"
      }
    }
  }
  
  df_imps.1se <- replace(df_imps.1se, is.na(df_imps.1se), 0) # replace NA values with 0
  df_imps.min <- replace(df_imps.min, is.na(df_imps.min), 0) # replace NA values with 0
  
  nonzero_preds.1se <- replace(df_imps.1se[,-c(1:2)], df_imps.1se[,-(1:2)] > 0, 1)
  nonzero_preds.min <- replace(df_imps.min[,-c(1:2)], df_imps.min[,-c(1:2)] > 0, 1)
  
  # Accuracy: VAF and adjusted VAF
  coef_determination.1se <- round(mean(accuracy.1se$VAF), 4)
  coef_determination.min <- round(mean(accuracy.min$VAF), 4)
  
  VAF_adj.1se <- round(mean(accuracy.1se$VAF_adj), 4)
  VAF_adj.min <- round(mean(accuracy.min$VAF_adj), 4)
  
  if(family %in% c("gaussian", "poisson")) {
    # Accuracy: 
    ## confidence interval
    
    MSE_confint.1se <- quantile(accuracy.1se$MSE, c(CI_alpha / 2, 1 - CI_alpha / 2))
    MSE_confint.min <- quantile(accuracy.min$MSE, c(CI_alpha / 2, 1 - CI_alpha / 2))
    
    ### summarize results
    accuracy_summary <- data.frame("avg_MSE" = c(round(mean(accuracy.1se$MSE), 3),
                                                 round(mean(accuracy.min$MSE), 3)),
                                   "sd_MSE"   = c(round(sd(accuracy.1se$MSE), 2),
                                                  round(sd(accuracy.min$MSE), 2)),
                                   "var_MSE" = c(round(var(accuracy.1se$MSE), 2),
                                                 round(var(accuracy.min$MSE), 2)),
                                   "Lower CI" = c(round(MSE_confint.1se[[1]], 2),
                                                  round(MSE_confint.min[[1]], 2)),
                                   "Upper CI" = c(round(MSE_confint.1se[[2]], 2),
                                                  round(MSE_confint.min[[2]], 2)),
                                   "VAF" = c(coef_determination.1se, 
                                             coef_determination.min),
                                   "VAF_adj" = c(round(VAF_adj.1se, 3),
                                                 round(VAF_adj.min, 3)))
    accuracy_sd <- c(accuracy_summary$sd_MSE[1],
                     accuracy_summary$sd_MSE[2])
    
    classification_table.1se = classification_table.min = NULL
    
  } else if (family == "binomial") {
    # Accuracy: 
    ## coef_determination and confidence interval
    coef_determination.1se <- round((1 - (mean(accuracy.1se$SEL)/var(y_obs))), 4)
    coef_determination.min <- round((1 - (mean(accuracy.min$SEL)/var(y_obs))), 4)
    
    SEL_confint.1se <- quantile(accuracy.1se$SEL, c(CI_alpha / 2, 1 - CI_alpha / 2))
    SEL_confint.min <- quantile(accuracy.min$SEL, c(CI_alpha / 2, 1 - CI_alpha / 2))
    
    ### summarize results
    accuracy_summary <- data.frame("avg_SEL" = c(round(mean(accuracy.1se$SEL), 3),
                                                 round(mean(accuracy.min$SEL), 3)),
                                   "sd_SEL"   = c(round(sd(accuracy.1se$SEL), 2),
                                                  round(sd(accuracy.min$SEL), 2)),
                                   "var_SEL" = c(round(var(accuracy.1se$SEL), 2),
                                                 round(var(accuracy.min$SEL), 2)),
                                   "Lower CI" = c(round(SEL_confint.1se[[1]], 2),
                                                  round(SEL_confint.min[[1]], 2)),
                                   "Upper CI" = c(round(SEL_confint.1se[[2]], 2),
                                                  round(SEL_confint.min[[2]], 2)),
                                   "VAF" = c(coef_determination.1se, 
                                             coef_determination.min),
                                   "VAF_adj" = c(round(VAF_adj.1se, 3),
                                                 round(VAF_adj.min, 3)),
                                   "avg_MCR" = c(round(mean(accuracy.1se$MCR), 2),
                                                 round(mean(accuracy.min$MCR), 2)),
                                   "avg_AUC" = c(round(mean(accuracy.1se$AUC), 2),
                                                 round(mean(accuracy.min$AUC), 2)))
    accuracy_sd <- c(accuracy_summary$sd_SEL[1],
                     accuracy_summary$sd_SEL[2])
  } else if (family == "multinomial") {
    accuracy_summary <- data.frame("Avg_PAC" = c(round(mean(accuracy.1se$PAC), 2),
                                                 round(mean(accuracy.min$PAC), 2)),
                                   "avg_AUC" = c(round(mean(accuracy.1se$AUC), 2),
                                                 round(mean(accuracy.min$AUC), 2)))
  }
  
  rownames(accuracy_summary) <- c("1se", "min")
  
  # Sparsity
  avg_num_predictors.1se <- round(mean(sparsity$count_predictors.1se), 2)
  sd_num_predictors.1se <- round(sd(sparsity$count_predictors.1se), 2)
  num_predictors_summary.1se <- summary(as.factor(sparsity$count_predictors.1se))
  num_predictors_confint.1se <- quantile(sparsity$count_predictors.1se, c(CI_alpha /2 ,
                                                                          1-CI_alpha/2))
  
  avg_num_predictors.min <- round(mean(sparsity$count_predictors.min), 2)
  sd_num_predictors.min <- round(sd(sparsity$count_predictors.min), 2)
  num_predictors_summary.min <- summary(as.factor(sparsity$count_predictors.min))
  num_predictors_confint.min <- quantile(sparsity$count_predictors.min, c(CI_alpha /2 , 
                                                                          1-CI_alpha/2))
  frequency_predictors.1se = table(predictors.1se)
  frequency_predictors.min = table(predictors.min)
  
  predictors_summary <- data.frame('average' = c(avg_num_predictors.1se,
                                                 avg_num_predictors.min),
                                   'sd' = c(sd_num_predictors.1se,
                                            sd_num_predictors.min),
                                   'lower_bound' = c(round(num_predictors_confint.1se[[1]], 2),
                                                     round(num_predictors_confint.min[[1]], 2)),
                                   'upper_bound' = c(round(num_predictors_confint.1se[[2]], 2),
                                                     round(num_predictors_confint.min[[2]], 2)))
  
  rownames(predictors_summary) <- c("1se", "min")
  
  # Rules
  avg_num_rules.1se <- round(mean(sparsity$count_rules.1se), 2)
  sd_num_rules.1se <- round(sd(sparsity$count_rules.1se), 2)
  num_rules_summary.1se <- summary(as.factor(sparsity$count_rules.1se))
  num_rules_confint.1se <- quantile(sparsity$count_rules.1se, c(CI_alpha /2 , 
                                                                1-CI_alpha/2)) 
  
  avg_num_rules.min <- round(mean(sparsity$count_rules.min), 2)
  sd_num_rules.min <- round(sd(sparsity$count_rules.min), 2)
  num_rules_summary.min <- summary(as.factor(sparsity$count_rules.min))
  num_rules_confint.min <- quantile(sparsity$count_rules.min, c(CI_alpha /2 , 
                                                                1-CI_alpha/2))
  
  rules_summary <- data.frame('average' = c(avg_num_rules.1se,
                                            avg_num_rules.min),
                              'sd' = c(sd_num_rules.1se,
                                       sd_num_rules.min),
                              'lower_bound' = c(round(num_rules_confint.1se[[1]]),
                                                round(num_rules_confint.min[[1]])),
                              'upper_bound' = c(round(num_rules_confint.1se[[2]]),
                                                round(num_rules_confint.min[[2]])))
  rownames(rules_summary) <- c("1se", "min")
  
  if(family %in% c("gaussian", "poisson", "binomial")) {
    # Stability
    ## of Importance measures
    imps_sd.1se <- apply(df_imps.1se[-c(1,2)], 2, function(x) sd(x, na.rm = T)) # variance
    imps_sd.min <- apply(df_imps.min[-c(1,2)], 2, function(x) sd(x, na.rm = T)) # variance
    imps_dist.1se <- as.vector(dist(df_imps.1se[-c(1,2)])) # distance
    imps_dist.min <- as.vector(dist(df_imps.min[-c(1,2)])) # distance
    
    imps_dist_avg.1se <- mean(imps_dist.1se, na.rm = T)
    imps_dist_avg.min <- mean(imps_dist.min, na.rm = T)
    
    
    ## of selection of predictors 
    nonzero_sd.1se <- sapply(nonzero_preds.1se, sd)
    nonzero_sd.min <- sapply(nonzero_preds.min, sd)
    
    ### stability measure proposed by Nogueira et al., 2018:
    d <- length(x_names)
    kbar.1se <- mean(sparsity$count_predictors.1se)
    kbar.min <- mean(sparsity$count_predictors.min)
    v_rand.1se <- (kbar.1se/d*(1- kbar.1se/d))
    v_rand.min <- (kbar.min/d*(1- kbar.min/d))
    
    stability_features.1se <- 1 - ((1/d * sum(nonzero_sd.1se^2))/ v_rand.1se)
    stability_features.min <- 1 - ((1/d * sum(nonzero_sd.min^2))/ v_rand.min)
    
    M <- nrow(nonzero_preds.1se)
    ki.1se <-rowSums(nonzero_preds.1se)
    hatPF.1se <- colMeans(nonzero_preds.1se)
    
    ki.min <-rowSums(nonzero_preds.min)
    hatPF.min <- colMeans(nonzero_preds.min)
    
    phi_i.1se = phi_i.min = rep(0,M)
    for(i in 1:M){ 
      phi_i.1se[i]<-(1/ v_rand.1se)*((1/d)*sum(nonzero_preds.1se[i,]*hatPF.1se)-(ki.1se[i]*kbar.1se)/d^2-(stability_features.1se/2)*((2*kbar.1se*   ki.1se[i])/d^2-ki.1se[i]/d-kbar.1se/d+1))
      phi_i.min[i]<-(1/ v_rand.min)*((1/d)*sum(nonzero_preds.min[i,]*hatPF.min)-(ki.min[i]*kbar.min)/d^2-(stability_features.min/2)*((2*kbar.min*   ki.min[i])/d^2-ki.min[i]/d-kbar.min/d+1))
    }
    
    phi_bar.1se = mean(phi_i.1se)
    var_stab.1se = (4/M^2)*sum((phi_i.1se-phi_bar.1se)^2)
    phi_bar.min = mean(phi_i.min)
    var_stab.min = (4/M^2)*sum((phi_i.min-phi_bar.min)^2)
    
    # CI of stability of features
    z <- qnorm(1-CI_alpha/2)
    stab_upper.1se <- stability_features.1se + z*sqrt(var_stab.1se)
    stab_lower.1se <- stability_features.1se - z*sqrt(var_stab.1se)
    stab_upper.min <- stability_features.min + z*sqrt(var_stab.min)
    stab_lower.min <- stability_features.min - z*sqrt(var_stab.min)
    
    stability_features <- data.frame("stability" = c(round(stability_features.1se, 2),
                                                     round(stability_features.min, 2)),
                                     "var_stab" = c(var_stab.1se,
                                                    var_stab.min),
                                     "lower" = c(stab_lower.1se,
                                                 stab_lower.min),
                                     "upper" = c(stab_upper.1se,
                                                 stab_upper.min))
    rownames(stability_features) = c("1se", "min")
    
    nonzero_preds_dist.1se <- as.vector(dist(nonzero_preds.1se, method = "binary"))
    nonzero_preds_dist.min <- as.vector(dist(nonzero_preds.min, method = "binary"))
    
    nonzero_preds_dist_avg.1se <- mean(nonzero_preds_dist.1se, na.rm = T)
    nonzero_preds_dist_avg.min <- mean(nonzero_preds_dist.min, na.rm = T)
    
    ## stability of the predictions
    preds_dist.1se <- as.vector(dist(t(cvpreds.1se)))
    preds_dist.min <- as.vector(dist(t(cvpreds.min)))
    
    preds_dist_avg.1se <- mean(preds_dist.1se, na.rm = T)
    preds_dist_avg.min <- mean(preds_dist.min, na.rm = T)
    
    stability <- data.frame("dist_imps" = c(round(imps_dist_avg.1se, 2),
                                            round(imps_dist_avg.min, 2)),
                            "dist_preds" = c(round(preds_dist_avg.1se, 2),
                                             round(preds_dist_avg.min, 2)),
                            "dist_preds_selection" = c(round(nonzero_preds_dist_avg.1se, 2), 
                                                       round(nonzero_preds_dist_avg.min, 2)),
                            "stability_features" = c(stability_features.1se,
                                                     stability_features.min),
                            "sd_MSE" = accuracy_sd,
                            "sd_imps" = c(round(mean(imps_sd.1se), 2),
                                          round(mean(imps_sd.min), 2)),
                            "sd_SELected" = c(mean(nonzero_sd.1se),
                                              mean(nonzero_sd.min)))
    
    rownames(stability) <- c("1se", "min")
    
    stability_predictors <- rbind(imps_sd.1se, imps_sd.min, 
                                  nonzero_sd.1se, nonzero_sd.min)
  }  
  
  
  ## return results
  return(invisible(list(accuracy.1se = accuracy.1se,
                        accuracy.min = accuracy.min,
                        cvpreds.1se = cvpreds.1se,
                        cvpreds.min = cvpreds.min,
                        accuracy_summary = accuracy_summary,
                        classification_table.1se = classification_table.1se,
                        classification_table.min = classification_table.min,
                        predictors_summary = predictors_summary,
                        rules_summary = rules_summary,
                        sparsity = sparsity,
                        df_imps.1se = df_imps.1se,
                        df_imps.min = df_imps.min,
                        stability_features = stability_features,
                        stability = stability,
                        stability_predictors = stability_predictors,
                        nonzero_preds.1se = nonzero_preds.1se,
                        nonzero_preds.min = nonzero_preds.min,
                        lambda = lambda,
                        imps_sd.1se = imps_sd.1se,
                        imps_sd.min = imps_sd.min,
                        imps_dist.1se = imps_dist.1se,
                        imps_dist.min = imps_dist.min,
                        nonzero_preds_dist.1se = nonzero_preds_dist.1se,
                        nonzero_preds_dist.min = nonzero_preds_dist.min,
                        preds_dist.1se = preds_dist.1se,
                        preds_dist.min = preds_dist.min,
                        num_predictors_summary.1se = num_predictors_summary.1se,
                        num_predictors_summary.min = num_predictors_summary.min,
                        frequency_predictors.1se = frequency_predictors.1se,
                        frequency_predictors.min = frequency_predictors.min,
                        imps_dist.1se = imps_dist.1se,
                        imps_dist.min = imps_dist.min,
                        preds_dist.1se = preds_dist.1se,
                        preds_dist.min = preds_dist.min,
                        y_names = y_names,
                        x_names = x_names,
                        n_obs = n,
                        y_ncol = y_ncol)))
  
}


##################################
# combine results function   #####
##################################

labels = c("Standard", "Relaxed", "Adaptive", "Relaxed Adaptive")

pre_combine <- function(models, labels, r = 10, k = 10, 
                        alpha = 0.05, family = "gaussian"){
  
  n_models <- length(models)
  MSE = SEL = mcr = auc = VAF_adj = count_terms = count_pred = c()
  
  
  imps_dist = nonzero_dist = preds_dist = c()
  
  stability_features <- data.frame()
  report_accuracy <- data.frame()
  report_features <- data.frame()
  report_terms <- data.frame()
  report_stability <- data.frame()
  
  for(model in models) {
    
    VAF_adj <- c(VAF_adj, 
                 model$accuracy.1se$VAF_adj,
                 model$accuracy.min$VAF_adj)
    
    count_terms <- c(count_terms, 
                     model$sparsity$count_rules.1se,
                     model$sparsity$count_rules.min)
    
    count_pred <- c(count_pred, 
                    model$sparsity$count_predictors.1se,
                    model$sparsity$count_predictors.min)
    
    stability_features <- rbind(stability_features, model$stability_features)
    
    
    
    imps_dist <- c(imps_dist, model$imps_dist.1se, 
                   model$imps_dist.min)
    nonzero_dist <- c(nonzero_dist, model$nonzero_preds_dist.1se, 
                      model$nonzero_preds_dist.min)
    preds_dist <- c(preds_dist, model$preds_dist.1se, 
                    model$preds_dist.min)
    
    report_features <- rbind(report_features, cbind(
      model$predictors_summary[,1],
      model$predictors_summary[,2]))
    
    report_terms <- rbind(report_terms, cbind(
      model$rules_summary[,1], 
      model$rules_summary[,2]))
    
    report_stability <- rbind(report_stability, cbind(
      model$stability_features$stability,
      model$stability$dist_preds,
      model$stability$dist_imps,
      model$stability$sd_imps))
    
    if(family %in% c("gaussian", "poisson")) {
      MSE <- c(MSE, model$accuracy.1se$MSE,
               model$accuracy.min$MSE)
      
      report_accuracy <- rbind(report_accuracy,cbind(
        model$accuracy_summary$avg_MSE, 
        model$accuracy_summary$sd_MSE, 
        model$accuracy_summary$VAF_adj))
      
    } else if(family == "binomial") {
      SEL <- c(SEL, model$accuracy.1se$SEL,
               model$accuracy.min$SEL)
      auc <- c(auc, model$accuracy.1se$AUC,
               model$accuracy.min$AUC)
      mcr <- c(mcr, model$accuracy.1se$MCR,
               model$accuracy.min$MCR)
      
      report_accuracy <- rbind(report_accuracy, cbind(
        model$accuracy_summary$avg_SEL,
        model$accuracy_summary$sd_SEL,
        model$accuracy_summary$VAF_adj,
        model$accuracy_summary$avg_AUC,
        model$accuracy_summary$avg_MCR))
      
      
    } 
  }
  
  if(family %in% c("gaussian", "poisson")) {
    accuracy <- data.frame("MSE" = MSE,
                           "VAF_adj" = VAF_adj,
                           "Criterion" = rep(c("1se", "min"), each = r, 
                                             times = n_models),
                           "Method" = rep(labels, each = r*2),
                           "Repeat" = rep(1:r, times = 2*n_models))
    accuracy2 <- accuracy[, 1:2]
    boxplot_accuracy <- ggplot(data = accuracy, aes(x = Criterion, y = MSE, 
                                                    fill = Method)) +
      geom_boxplot() +
      labs(x = "Criterion", y = "MSE") +
      scale_fill_manual(values = rainbow_hcl(n = 4)) +
      theme_minimal()  +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 13),
            strip.text = element_text(size = 13),
            legend.position = "none")
    
    names(report_accuracy) <- c("& $M$", "& $SD$", "& $R_{adj}^2$")
    
  } else if(family == "binomial") {
    accuracy <- data.frame("SEL" = SEL,
                           "VAF_adj" = VAF_adj,
                           "auc" = auc,
                           "mcr" = mcr,
                           "Criterion" = rep(c("1se", "min"), each = r, 
                                             times = n_models),
                           "Method" = rep(labels, each = r*2),
                           "Repeat" = rep(1:r, times = 2*n_models))
    accuracy2 <- accuracy[, 1:4]
    boxplot_accuracy <- ggplot(data = accuracy, aes(x = Criterion, y = SEL, 
                                                    fill = Method)) +
      geom_boxplot() +
      labs(x = "Criterion", y = "SEL") +
      scale_fill_manual(values = rainbow_hcl(n = 4)) +
      theme_minimal()  +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 13),
            strip.text = element_text(size = 13),
            legend.position = "none")
    
    
    names(report_accuracy) <- c("& $M$", "& $SD$", "& $R_{adj}^2$", "& $AUC$", "& $mcr$")
    
  }
  
  sparsity <- data.frame("count_terms" = count_terms,
                         "count_pred" = count_pred,
                         "Criterion" = rep(c("1se", "min"), each = r*k, 
                                           times = n_models),
                         "Method" = rep(labels, each = r*k*2),
                         "Repeat" = rep(1:r, each = k, times = 2*n_models))
  
  distances <- data.frame("imps_dist" = imps_dist,
                          "nonzero_dist" = nonzero_dist,
                          "criterion" = rep(c("1se", "min"), 
                                            each = choose(k*r, 2), 
                                            times = n_models),
                          "method" = rep(labels, each = choose(k*r, 2) * 2, 
                                         times = 1))
  distances_preds <- data.frame("preds_dist" = preds_dist,
                                "criterion" = rep(c("1se", "min"), 
                                                  each = choose(r, 2), 
                                                  times = n_models),
                                "method" = rep(labels, each = choose(r, 2) * 2))
  
  stability_features$method <- rep(labels, each = 2)
  stability_features$criterion <- rep(c("1se", "min"), times = n_models)
  
  
  names(report_stability) <- c("& $\\hat{\\phi}(Z)$", "& d(predictions)",
                               "& d(importances)", "& $SD$ of importances")
  names(report_features) = names(report_terms) = c("& $M$", "& $SD$")
  
  
  # Creating data frames with dummy coding
  accuracy2$One_se <- ifelse(accuracy$Criterion == "1se", 1, 0)
  accuracy2$Relaxed <- ifelse(accuracy$Method %in% 
                                c("Relaxed", "Relaxed Adaptive"),
                              1, 0)
  accuracy2$Adaptive <- ifelse(accuracy$Method %in% 
                                 c("Adaptive", "Relaxed Adaptive"),
                               1, 0)
  accuracy2$Repeat <- accuracy$Repeat
  
  sparsity2 <- data.frame("count_terms" = count_terms,
                          "count_pred" = count_pred,
                          "One_se" = ifelse(sparsity$Criterion == "1se", 1, 0),
                          "Relaxed" = ifelse(sparsity$Method %in% 
                                               c("Relaxed", "Relaxed Adaptive"),
                                             1, 0),
                          "Adaptive" = ifelse(sparsity$Method %in% 
                                                c("Adaptive", "Relaxed Adaptive"),
                                              1, 0),
                          "Repeat" = rep(1:r, each = k, times = 2*n_models))
  
  distances2 <- data.frame("imps_dist" = imps_dist,
                           "nonzero_dist" = nonzero_dist,
                           "One_se" = ifelse(distances$criterion == "1se", 1, 0),
                           "Relaxed" = ifelse(distances$method %in% 
                                                c("Relaxed", "Relaxed Adaptive"),
                                              1, 0),
                           "Adaptive" = ifelse(distances$method %in% 
                                                 c("Adaptive", "Relaxed Adaptive"),
                                               1, 0))
  distances_preds2 <- data.frame("preds_dist" = preds_dist,
                                 "One_se" = ifelse(distances_preds$criterion == "1se", 
                                                   1, 0),
                                 "Relaxed" = ifelse(distances_preds$method %in% 
                                                      c("Relaxed", "Relaxed Adaptive"),
                                                    1, 0),
                                 "Adaptive" = ifelse(distances_preds$method %in% 
                                                       c("Adaptive", "Relaxed Adaptive"),
                                                     1, 0))
  
  
  # Model Significance Test
  ## accuracy
  if(family  %in% c("gaussian", "poisson")) {
    model_MSE <- lmer(MSE ~ One_se*Relaxed*Adaptive + (1|Repeat),
                      data = accuracy2)
    ranova_MSE <- ranova(model_MSE)
    
    if(ranova_MSE$`Pr(>Chisq)`[2] > alpha) {
      model_MSE <- lm(MSE ~ One_se*Relaxed*Adaptive,
                      data = accuracy2)
    } 
    model_SEL <- NULL
    df <-as.data.frame(summary(model_MSE)$coefficients)
    
    model_MSE_latex <- latex(df, 
                             caption = "Mixed model of MSE by method",
                             model = class(model_MSE),
                             footnote = footnote)
    model_SEL_latex <- NULL
  } else if(family == "binomial") {
    model_SEL <- lmer(SEL ~ One_se*Relaxed*Adaptive + (1|Repeat),
                      data = accuracy2)
    ranova_SEL <- ranova(model_SEL)
    
    if(ranova_SEL$`Pr(>Chisq)`[2] > alpha) {
      model_SEL <- lm(SEL ~ One_se*Relaxed*Adaptive,
                      data = accuracy2)
    } 
    model_MSE <- NULL
    df <-as.data.frame(summary(model_SEL)$coefficients)
    
    model_SEL_latex <- latex(df, 
                             caption = "Mixed model of SEL by method",
                             model = class(model_SEL),
                             footnote = footnote)
    model_MSE_latex <- NULL
  }
  
  
  
  model_VAF <- lmer(VAF_adj ~ One_se*Relaxed*Adaptive + (1|Repeat),
                    data = accuracy2)
  ranova_VAF <- ranova(model_VAF)
  
  if(ranova_VAF$`Pr(>Chisq)`[2] > alpha) {
    model_VAF <- lm(VAF_adj ~ One_se*Relaxed*Adaptive,
                    data = accuracy2)
  }
  
  df <-as.data.frame(summary(model_VAF)$coefficients)
  
  model_VAF_latex <- latex(df, 
                           caption = "Mixed model of adjusted VAF by method",
                           model = class(model_VAF),
                           footnote = footnote)
  
  ## sparsity
  model_sparsity <- lmer(count_pred ~ One_se*Relaxed*Adaptive + (1|Repeat),
                         data = sparsity2)
  ranova_sparsity <- ranova(model_sparsity)
  
  if (ranova_sparsity$`Pr(>Chisq)`[2] > alpha) {
    model_sparsity <- lm(count_pred ~ One_se*Relaxed*Adaptive,
                         data = sparsity2)
  }
  
  df <-as.data.frame(summary(model_sparsity)$coefficients)
  
  model_sparsity_latex <- latex(df, 
                                caption = "Mixed model of adjusted VAF by method",
                                model = class(model_sparsity),
                                footnote = footnote)
  
  ## stability
  
  model_dist_imps <- lm(imps_dist ~ One_se*Relaxed*Adaptive,
                        data = distances2)
  model_dist_preds <- lm(preds_dist ~ One_se*Relaxed*Adaptive,
                         data = distances_preds2)
  
  df <-as.data.frame(summary(model_dist_imps)$coefficients)
  
  model_dist_imps_latex <- latex(df, 
                                 caption = "Linear model of distances between importance measures by method",
                                 model = "lm",
                                 footnote = footnote)
  
  df <-as.data.frame(summary(model_dist_preds)$coefficients)
  
  model_dist_preds_latex <- latex(df, 
                                  caption = "Linear model of distances between predictions by method",
                                  model = "lm",
                                  footnote = footnote)
  
  # visualizations
  boxplot_VAF <- ggplot(data = accuracy, 
                        aes(x = Criterion, y = VAF_adj, fill = Method)) +
    geom_boxplot() +
    labs(x = "Criterion", y = "adjusted VAF") +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  boxplot_features <- ggplot(sparsity, 
                             aes(x = Method, y = count_pred, fill = Method)) +
    geom_boxplot(position = "identity") + 
    facet_wrap( ~ Criterion) +
    labs(y = "Number of Predictors") +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    scale_x_discrete(breaks = c("", "", "", "")) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  boxplot_terms <- ggplot(sparsity,
                          aes(x = Method, y = count_terms, fill = Method)) +
    geom_boxplot(position = "identity") + 
    facet_wrap( ~ Criterion) +
    labs(y = "Number of Terms") +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    scale_x_discrete(breaks = c("", "", "", "")) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  plot_stability <- ggplot(stability_features,
                           aes(x = method, y = stability, color = method)) +
    geom_segment(aes(xend = method, yend = upper)) +
    geom_segment(aes(xend = method, yend = lower)) +
    geom_point(aes(y = stability, color = method)) +
    facet_wrap( ~ criterion) +
    labs(x = "Method", y = "Stability", color = "Method") +
    scale_x_discrete(breaks = c("", "", "", "")) +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  density_imps_dist <- ggplot(distances, 
                              aes(x = imps_dist, color = method)) +
    geom_density() +
    facet_wrap(~ criterion) +
    labs(x = "Distances between Importances", 
         y = "Density", color = "Method") +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  boxplot_imps_dist <- ggplot(distances, 
                              aes(y = imps_dist, fill = method)) +
    geom_boxplot() +
    facet_wrap(~ criterion) +
    labs(y = "Distances between Importances", 
         fill = "Method") +
    scale_x_discrete(breaks = c("", "", "", "")) +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  density_preds_dist <- ggplot(distances_preds, 
                               aes(x = preds_dist, color = method)) +
    geom_density() +
    facet_wrap(~ criterion) +
    labs(x = "Distances between Predictions", 
         y = "Density", color = "Method") +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal()  +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  boxplot_preds_dist <- ggplot(distances_preds, 
                               aes(y = preds_dist, fill = method)) +
    geom_boxplot() +
    facet_wrap(~ criterion) +
    labs(y = "Distances between Predictions",
         x = "Method",
         fill = "Method") +
    scale_x_discrete(breaks = c("", "", "", "")) +
    scale_fill_manual(values = rainbow_hcl(n = 4)) +
    theme_minimal() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 13),
          legend.position = "none")
  
  
  return(invisible(list("report_accuracy" = report_accuracy,
                        "report_features" = report_features,
                        "report_terms" = report_terms,
                        "report_stability" = report_stability,
                        "accuracy" = accuracy,
                        "sparsity" = sparsity,
                        "stability_features" = stability_features,
                        "distances" = distances,
                        "distances_preds" = distances_preds,
                        "accuracy2" = accuracy2,
                        "sparsity2" = sparsity2,
                        "distances2" = distances2,
                        "distances_preds2" = distances_preds2,
                        "boxplot_accuracy" = boxplot_accuracy,
                        "boxplot_VAF" = boxplot_VAF,
                        "boxplot_features" = boxplot_features,
                        "boxplot_terms" = boxplot_terms,
                        "plot_stability" = plot_stability,
                        "density_imps_dist" = density_imps_dist,
                        "boxplot_imps_dist" = boxplot_imps_dist,
                        "density_preds_dist" = density_preds_dist,
                        "boxplot_preds_dist" = boxplot_preds_dist,
                        "model_sparsity" = model_sparsity,
                        "model_MSE" = model_MSE,
                        "model_SEL" = model_SEL,
                        "model_VAF" = model_VAF,
                        "model_dist_imps" = model_dist_imps,
                        "model_dist_preds" = model_dist_preds,
                        "model_MSE_latex" = model_MSE_latex,
                        "model_SEL_latex" = model_SEL_latex,
                        "model_VAF_latex" = model_VAF_latex,
                        "model_sparsity_latex" = model_sparsity_latex,
                        "model_dist_imps_latex" = model_dist_imps_latex,
                        "model_dist_preds_latex" = model_dist_preds_latex)))
  
}





################################################################################
########################## Dataset 1: Highschool Grades ########################
################################################################################

data1 <- read_delim("1_student-por.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

data1 <- data1[,c(1:30, 33)]

data1 <- data1 %>%
  mutate_if(is.character, as.factor)


set.seed(1978594)
mod1 <- rep_cv_pre(G3 ~ ., data1,
                   verbose = T)

set.seed(1978594)
mod1.r <- rep_cv_pre(G3 ~ ., data1,
                     verbose = T,
                     relax = T)
set.seed(1978594)
mod1.a <- rep_cv_pre(G3 ~ ., data1,
                     verbose = T,
                     ad.alpha = 0)
set.seed(1978594)
mod1.ra <- rep_cv_pre(G3 ~ ., data1,
                      verbose = T,
                      relax = T,
                      ad.alpha = 0)

result1 <- pre_combine(list(mod1, mod1.r, mod1.a, mod1.ra), 
                       labels)

# print latex tables
print(latex(t(result1$report_accuracy)))
print(latex(t(result1$report_features)))
print(latex(t(result1$report_terms)))
print(latex(t(result1$report_stability)))

ggsave("MSE1.png", plot =result1$boxplot_accuracy, 
       dpi = 300, bg = "white")

ggsave("VAF1.png", plot =result1$boxplot_VAF, 
       dpi = 300, bg = "white")

ggsave("boxplot_features1.png", plot =result1$boxplot_features,
       dpi = 300, bg = "white")

ggsave("boxplot_terms1.png", plot =result1$boxplot_terms,
       dpi = 300, bg = "white")

ggsave("boxplot_preds_dist1.png", plot =result1$boxplot_preds_dist,
       dpi = 300, bg = "white")

ggsave("boxplot_imps_dist1.png", plot =result1$boxplot_imps_dist,
       dpi = 300, bg = "white")

ggsave("density_preds_dist1.png", plot =result1$density_preds_dist,
       dpi = 300, bg = "white")

ggsave("density_imps_dist1.png", plot =result1$density_imps_dist,
       dpi = 300, bg = "white")

ggsave("stability1.png", plot =result1$plot_stability,
       dpi = 300, bg = "white")

print(summary(result1$model_sparsity))
print(summary(result1$model_MSE))
print(summary(result1$model_VAF))
print(summary(result1$model_dist_imps))
print(summary(result1$model_dist_preds))

print(result1$model_sparsity_latex)
print(result1$model_MSE_latex)
print(result1$model_VAF_latex)
print(result1$model_dist_imps_latex)
print(result1$model_dist_preds_latex)



################################################################################
############### Dataset 2: Sensation Seeking and Delinquency ###################
################################################################################

data2 <- read_sav("2_sensation seeking.sav")

data2 <- data2 %>%
  mutate_if(is.labelled, as.numeric)

data2 <- data2[complete.cases(data2),]

# Outcome: delinquent behaviour
data2$delinquent <- data2$YSR40 + data2$YSR72 + data2$YSR81 + data2$YSR82 + data2$YSR101

data2 <- data2[, c(2:26, 32)]

set.seed(1978594)
mod2 <- rep_cv_pre(delinquent ~ ., data2,
                   verbose = T)

set.seed(1978594)
mod2.r <- rep_cv_pre(delinquent ~ ., data2,
                     verbose = T,
                     relax = T)
set.seed(1978594)
mod2.a <- rep_cv_pre(delinquent ~ ., data2,
                     verbose = T,
                     ad.alpha = 0)
set.seed(1978594)
mod2.ra <- rep_cv_pre(delinquent ~ ., data2,
                      verbose = T,
                      relax = T,
                      ad.alpha = 0)

result2 <- pre_combine(list(mod2, mod2.r, mod2.a, mod2.ra), labels)

# print latex tables
print(latex(t(result2$report_accuracy)))
print(latex(t(result2$report_features)))
print(latex(t(result2$report_terms)))
print(latex(t(result2$report_stability)))

ggsave("MSE2.png", plot =result2$boxplot_accuracy, 
       dpi = 300, bg = "white")

ggsave("VAF2.png", plot =result2$boxplot_VAF, 
       dpi = 300, bg = "white")

ggsave("boxplot_features2.png", plot =result2$boxplot_features,
       dpi = 300, bg = "white")

ggsave("boxplot_terms2.png", plot =result2$boxplot_terms,
       dpi = 300, bg = "white")

ggsave("boxplot_preds_dist2.png", plot =result2$boxplot_preds_dist,
       dpi = 300, bg = "white")

ggsave("boxplot_imps_dist2.png", plot =result2$boxplot_imps_dist,
       dpi = 300, bg = "white")

ggsave("density_preds_dist2.png", plot =result2$density_preds_dist,
       dpi = 300, bg = "white")

ggsave("density_imps_dist2.png", plot =result2$density_imps_dist,
       dpi = 300, bg = "white")

ggsave("stability2.png", plot =result2$plot_stability,
       dpi = 300, bg = "white")

print(summary(result2$model_sparsity))
print(summary(result2$model_MSE))
print(summary(result2$model_VAF))
print(summary(result2$model_dist_imps))
print(summary(result2$model_dist_preds))


print(result2$model_sparsity_latex)
print(result2$model_MSE_latex)
print(result2$model_VAF_latex)
print(result2$model_dist_imps_latex)
print(result2$model_dist_preds_latex)



