
library(Mime1)

# Load HCC survival data from multiple cohorts
HCC_cohorts_survival <- readRDS("HCC_cohorts_survival.rds")

# Format input data: standardize column names and order
HCC_cohorts_survival_101ML <- list()
for(i in 1:length(HCC_cohorts_survival)) {
  df <- HCC_cohorts_survival[[i]]
  colnames(df)[2] <- "OS"       # Overall Survival status
  colnames(df)[3] <- "OS.time"  # Overall Survival time
  df <- df[, c(1, 3, 2, 4:ncol(df))]  # Reorder columns
  df_name <- names(HCC_cohorts_survival)[i]
  HCC_cohorts_survival_101ML[[df_name]] <- df
}

# Replace OS.time = 0 with small value (0.003) to avoid errors in survival models
for(i in 1:length(HCC_cohorts_survival_101ML)) {
  df <- HCC_cohorts_survival_101ML[[i]]
  colnames(df)[1] <- "ID"
  if("OS.time" %in% colnames(df)) {
    df$OS.time[df$OS.time == 0] <- 0.003
  }
  HCC_cohorts_survival_101ML[[i]] <- df
}

list_train_vali_Data <- HCC_cohorts_survival_101ML

# Load candidate feature gene list
HCC_PFAS <- read.csv("HCC_PFAS_genes.txt", sep="")
genelist <- HCC_PFAS$Gene_HCC_PFAS_83_0.4

# Construct machine learning prognostic model using 101 model combinations
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize = 5,
                       seed = 5201314)

# Evaluate model C-index across cohorts
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-4],
               order = names(list_train_vali_Data),
               width = 0.35)

save(res, file = "res.Rdata")