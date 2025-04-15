
## Create function: Enet Regression + ROC curve + KM plot

run_enet_train_validate <- function(train_df, test_list, gene_list, alpha = 0.1) {
  library(survival)
  library(glmnet)
  library(timeROC)
  library(survminer)
  library(ggplot2)
  
  # ==== 1. TCGA: Train Ridge Model ====
  message("Training on TCGA...")
  
  train_df <- train_df[, c("ID", "OS.time", "OS", gene_list)]
  rownames(train_df) <- train_df$ID
  train_df <- train_df[, -1]
  
  surv_train <- Surv(train_df$OS.time, train_df$OS)
  X_train <- as.matrix(train_df[, !(colnames(train_df) %in% c("OS.time", "OS"))])
  
  ridge_model <- cv.glmnet(X_train, surv_train, alpha = alpha, family = "cox")
  best_lambda <- ridge_model$lambda.min
  message("[TCGA] Best lambda: ", round(best_lambda, 5))
  
  # ==== 1.1 Diagnostic plots ====
  par(mfrow = c(1, 2))
  plot(ridge_model$glmnet.fit, xvar = "lambda", label = TRUE, main = "C. Coefficient Path")
  plot(ridge_model, xvar = "lambda", label = TRUE, main = "D. CV Deviance")
  abline(v = log(best_lambda), col = "red", lty = 2, lwd = 2)
  text(x = log(best_lambda), y = min(ridge_model$cvm), 
       labels = paste("λ.min =", round(log(best_lambda), 2)),
       col = "red", pos = 4)
  
  # ==== 1.2 Select final genes ====
  lasso_coef <- as.matrix(coef(ridge_model, s = best_lambda))
  selected_genes <- rownames(lasso_coef)[lasso_coef != 0]
  selected_genes <- setdiff(selected_genes, "(Intercept)")
  message("Selected genes (", length(selected_genes), "): ", paste(selected_genes, collapse = ", "))
  
  # ==== 2. Validation cohorts ====  
  results <- list()
  colors <- c("#d42921", "#124f7b", "#fabf3d")  # consistent ROC colors

  # Risk Score = exp(-0.0136531033 × ESR1 - 0.0057401447 × APOA1 - 0.0127939803 × IGF1 - 0.0601279664 × PPARGC1A + 0.0273195307 × SERPINE1 + 0.0476269799 × HMOX1 - 0.0001355368 × APCS - 0.0629751524 × ACADS - 0.0080130999 × SLC10A1 - 0.0245156728 × SLC2A2 - 0.0469738128 × ACAT1 - 0.0089196823 × C1S - 0.0848350866 × LCAT)
  
  for (name in names(test_list)) {
    message("Validating: ", name)
    
    test_df <- test_list[[name]]
    test_df <- test_df[, c("ID", "OS.time", "OS", selected_genes)]
    rownames(test_df) <- test_df$ID
    test_df <- test_df[, -1]
    
    surv_test <- Surv(test_df$OS.time, test_df$OS)
    X_test <- as.matrix(test_df[, !(colnames(test_df) %in% c("OS.time", "OS"))])
    
    predictions <- predict(ridge_model, newx = X_test, s = best_lambda)
    test_df$PFASHR.Sig <- as.numeric(predictions)
    test_df$group <- ifelse(test_df$PFASHR.Sig > median(test_df$PFASHR.Sig),
                            "High PFASHR.Sig", "Low PFASHR.Sig")
    test_df$group <- factor(test_df$group, levels = c("High PFASHR.Sig", "Low PFASHR.Sig"))
    
    ROC <- timeROC(T = test_df$OS.time, delta = test_df$OS,
                   marker = test_df$PFASHR.Sig,
                   cause = 1, weighting = "marginal",
                   times = c(12, 36, 50), iid = TRUE)
    
    # ==== ROC Plot ====
    plot(ROC, time = 12, col = colors[1], lwd = 2, title = paste(name, "- ROC"))
    plot(ROC, time = 36, col = colors[2], add = TRUE, lwd = 2)
    plot(ROC, time = 50, col = colors[3], add = TRUE, lwd = 2)
    legend("bottomright",
           c(paste0("AUC@1yr: ", round(ROC$AUC[1], 3)),
             paste0("AUC@3yr: ", round(ROC$AUC[2], 3)),
             paste0("AUC@5yr: ", round(ROC$AUC[3], 3))),
           col = colors, lty = 1, lwd = 2, bty = "n")
    
    # ==== KM Plot ====
    km_fit <- survfit(surv_test ~ group, data = test_df)
    print(ggsurvplot(km_fit, data = test_df,
                     conf.int = TRUE, pval = TRUE, fun = "pct",
                     size = 1, linetype = "strata",
                     palette = c("#e8acb3", "#46a5c7"),
                     legend = c(0.8, 0.85),
                     legend.title = name,
                     legend.labs = c("High PFASHRSig", "Low PFASHRSig")))
    
    results[[name]] <- list(roc = ROC, km = km_fit, data = test_df)
  }
  
  return(list(results = results,
              selected_genes = selected_genes,
              ridge_model = ridge_model,
              best_lambda = best_lambda))
}

## Function application

train_df <- list_train_vali_Data$TCGA
test_list <- list_train_vali_Data[names(list_train_vali_Data) != "TCGA"]
res <- run_enet_train_validate(train_df, test_list, final_sig)
