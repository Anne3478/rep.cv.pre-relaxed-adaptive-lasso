#' repcvpre
#'
#' @description
#' Conducts repeated k-fold cross-validation of prediction rule ensembles with
#' function \code{\link{pre}} (Fokkema, 2020). In each fold and repeat the lambda values,
#' importance measures, accuracy on training data, variable selection, and base
#' learner selection are recorded. In addition, summaries for each can be returned.
#' The function can be used to assess model performance by providing measures of
#' accuracy, sparsity and stability.
#'
#' @title Repeated k-fold Cross-validation of Prediction Rule Ensembles
#'
#'
#'
#'
#' @param formula provides a symbolic description of the model to be fitted
#' @param data specifies a data frame containing the variables in the model
#' @param k integer. Number of folds in the k-fold cross-validation
#' @param r integer. Number of repeats in the repeated k-fold cross-validation
#' @param CI_alpha numeric. Significance level used to calculate confidence intervals
#' @param verbose boolean. Specifies whether information on progress of repeats
#' should be printed to the command line
#' @param foldids numeric matrix with number of rows equal to number of rows in
#' the data and number of columns equal to number of repeats r. Each column
#' contains the group number between 1 to k to which the observation belongs
#' in the cross-validation, for a given repeat. Defaults to NULL, resulting in
#' the original training observations being randomly assigned to one of the k
#' folds in each repeat.
#' @param family character. Specifies a GLM-family object or a character string.
#' By default, family = gaussian, which assumes a single continuous response
#' variable. Alternatively, if the response is a single binary factor "binomial",
#' for a count response "poisson", and for a factor with > 2 levels "multinomial"
#' @param pclass numeric. Only used for binary classification. Cut-off value for
#' the predicted probabilities that should be used to classify observations to
#' the second class
#' @param ... specifies additional arguments that will be passed on to function
#' 'pre'
#'
#'
#' @return A list containing dataframes with the results from the repeated
#' k-fold cross-validation.
#'
#'
#' \describe{
#'   \item{\code{accuracy.min} and \code{accuracy.1se}}{Dataframe containing the
#'   Mean Squared Error (MSE), Mean Absolute Error (MAE), variance accounted for (VAF),
#'   and adjusted VAF for each repeat, either based on the min or 1se criterion
#'   respectively. For binary outcomes, the Squared Error Loss (SEL) and Absolute
#'   Error Loss (AEL) instead of the MSE and MAE are reported. In addition, for
#'   binary outcomes the dataframe contains the AUC, misclassification rate (mcr),
#'   percentage accurately classified (PAC), the sensitivity, and specificity.}
#'   \item{\code{accuracy_summary}}{Dataframe containing a summary of the accuracy
#'   measures for both the 1se and min criterion.}
#'   \item{\code{predictors_summary}}{Dataframe containing a summary of the number of
#'   predictors selected, including mean, standard deviation, and confidence interval.}
#'   \item{\code{rules_summary}}{Dataframe containing a summary of the number of
#'   base learners selected; including mean, standard deviation, and confidence interval.}
#'   \item{\code{sparsity}}{Dataframe containing the number of predictors and
#'   base learners selected in each fold and repeat.}
#'   \item{\code{df_imps.1se} and \code{df_imps.min}}{Dataframe containing
#'   importance measures for each predictor for each fold and repeat.}
#'   \item{\code{stability_features}}{Stability of number of predictors selected;
#'   calculated based on a formula by Nogueira et al. (2018). Contains the stability,
#'   its variance, and confidence interval.}
#'   \item{\code{stability}}{Dataframe containing a summary of all stability
#'   measures. This includes euclidean distances between importances, and
#'   between predictions, the standard deviation of the MSE/SEL, of the importances,
#'   and of the variable selection, and stability of predictors selected.}
#'   \item{\code{nonzero_preds.1se} and \code{nonzero_preds.min}}{Dataframe
#'   with a column for each predictor and a row for each fold and repeat, containing
#'   information if the predictor has been selected (1) or not (0).}
#'   \item{\code{lambda}}{Dataframe containing lambda values for each fold and
#'   repeat.}
#'   \item{\code{cvpreds.1se} and \code{cvpreds.min}}{Dataframe containing the
#'   predictions made on the training data in the k-fold cross-validation for
#'   each repeat.}
#'   \item{\code{classification_table.1se} and \code{classification_table.min}}{
#'   Confusion matrix with predictions. Only for binary and multinomial outcomes.}
#'   \item{\code{stability_predictors}}{Dataframe containing for each predictor,
#'   the standard deviation of importance measures and the standard deviation of
#'   times the importance measure was nonzero (predictor selected).}
#'   \item{\code{imps_dist.1se} and \code{imps_dist.min}}{Numeric vector containing
#'   the euclidean distances between the importance measures of each predictor
#'   across the folds and repeats. Used to assess the stability of important measures.}
#'   \item{\code{preds_dist.1se} and \code{preds_dist.min}}{Numeric vector containing
#'   the euclidean distances between the predictions made for the same observations
#'   across different repeats. Used to assess the stability of the predictions across
#'   various test and training data assignments.}
#'   \item{\code{num_predictors_summary.1se} and \code{num_predictors_summary.min}}{
#'   Table containing the frequency a number of predictors selection has occurred.}
#'   \item{\code{frequency_predictors.1se} and \code{frequency_predictors.min}}{
#'   Frequency table of how many times each predictors has been selected across
#'   the folds and repeats.}
#'   \item{\code{n_obs}}{Number of observations in dataset.}
#'}
#'
#'
#' @seealso \code{\link{pre}}, \code{\link{cvpre}}, \code{\link{importance}}
#'
#' @export
#'
#' @importFrom Formula Formula
#' @importFrom pROC roc
#' @importFrom pROC multiclass.roc
#' @importFrom pre pre
#' @importFrom pre importance
#'
#'
#' @importFrom stats dist model.matrix predict qnorm quantile sd terms var
#'
#' @details
#' In each cross-validation, the results were recorded for the selection of the lambda
#' parameter both based on the 1 standard deviation and minimum criterion, denoted
#' by 1se and min respectively. According to the min criterion, lambda resulting
#' in the lowest MSE is selected while according to the 1se criterion, lambda
#' resulting in the MSE one standard deviation above the minimum is selected. The
#' summary output dataframes contain a row for the results based on the 1se criterion
#' and min criterion. For the other dataframes, a dataframe for each of the two
#' criteria exists, denoted in the ending '.1se' or '.min'.
#'
#' @references
#' Fokkema, M. (2020). Fitting prediction rule ensembles with R package pre.
#' \emph{Journal of Statistical Software, 92}(12), 1–30.
#' https://doi.org/10.18637/jss.v092.i12
#'
#' Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of feature
#' selection algorithms. \emph{The Journal of Machine Learning Research, 18}(1), 6345–6398.
#'
#'
#'
#' @author Anne Hilbert
#'
#'
#' @examples
#' ens <- repcvpre(Ozone ~ .,  data = airquality[complete.cases(airquality),])
#' print(ens$accuracy_summary)
#' print(ens$predictors_summary)
#' print(ens$rules_summary)
#' print(ens$stability)
#' print(ens$frequency_predictors.1se)
#'
repcvpre <- function(formula, data, k = 10, r = 10, CI_alpha = 0.05, verbose = F,
                     foldids = NULL, family = "gaussian", pclass = 0.5, ...) {
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
