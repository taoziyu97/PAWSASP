library("survival")
library("survminer")
load("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/signature_by_source/TCGA_sig_cox.rds")

# TCGA KM analysis
TCGA_sig_list$`cut_0.00` <- lapply(TCGA_sig_list$data, function(x){km_cut_list(x, 0.00)})
TCGA_sig_list$`cut_0.25` <- lapply(TCGA_sig_list$data, function(x){km_cut_list(x, 0.25)})
TCGA_sig_list$`cut_0.50` <- lapply(TCGA_sig_list$data, function(x){km_cut_list(x, 0.50)})

km_cut_list <- function(sig_data, cut){
  # sig_data = test_data
  sig_data <- as.data.frame(sig_data)
  sig_data[,5:86] <- apply(sig_data[,5:86], 2, function(x){cut_off_set(x, cut)})
  return(sig_data)
}

cut_off_set <- function(sig_data, cut){
  quan <- quantile(sig_data, probs = cut)
  cut_0 <- ifelse(sig_data < quan, 1, 2)
}
TCGA_PAN_km <- matrix(data = NA, nrow = 64, ncol = 6, 
                      dimnames = list(c(sig_name),
                                      c("source","scale",
                                        "0%","25%","50%","var"))) %>%
  as.data.frame()

# TCGA_sig_list$km <- lapply(TCGA_sig_list$data, TCGA_km)

TCGA_km <- function(data){
  for (i in 1:length(sig_name)){
    f <- as.formula(paste("Surv(OS.time, OS) ~ ", paste(sig_name[i])))
    fit_0.00 <- survfit(f, data = as.data.frame(TCGA_sig_list$cut_0.00[i]))
    fit_0.25 <- survfit(f, data = as.data.frame(TCGA_sig_list$cut_0.25[i]))
    fit_0.50 <- survfit(f, data = as.data.frame(TCGA_sig_list$cut_0.50[i]))
    summary_km_0 <- summary(fit_0.00)
    summary_km_0.25 <- summary(fit_0.25)
    summary_km_0.50 <- summary(fit_0.50)
    TCGA_PAN_km[i, 1:2] <- c( "TCGA", "WES")
    TCGA_PAN_km[[3]][[i]] <- summary_km_0[17]
    TCGA_PAN_km <- TCGA_PAN_km %>% na.omit()
  }
  return(TCGA_PAN_km)
}

# TCGA_sig_list$km <- lapply(TCGA_sig_list$data, TCGA_km)
