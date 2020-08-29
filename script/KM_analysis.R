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
TCGA_PAN_km <- matrix(data = NA, nrow = 64, ncol = 4, 
                      dimnames = list(c(sig_name),
                                      c("source","scale",
                                        "km list","var"))) %>%
  as.data.frame()

TCGA_sig_list$`km_0` <- lapply(TCGA_sig_list$cut_0.00, function(x){TCGA_km(x)})
TCGA_sig_list$`km_0.25` <- lapply(TCGA_sig_list$cut_0.25, function(x){TCGA_km(x)})
TCGA_sig_list$`km_0.50` <- lapply(TCGA_sig_list$cut_0.50, function(x){TCGA_km(x)})

TCGA_km <- function(data){
  # data = TCGA_sig_list$cut_0.25[1] %>% as.data.frame()
  for (i in 1:length(sig_name)){
    f <- as.formula(paste("Surv(OS.time, OS) ~ ", paste(sig_name[i])))
    fit <- survfit(f, data = data)
    summary_km <- summary(fit)
    TCGA_PAN_km[i, 1:2] <- c( "TCGA", "WES")
    TCGA_PAN_km[i, 4] <- sig_name[i]
    # TCGA_PAN_km[[4]][[i]] <- summary_km_0[17]
    TCGA_PAN_km[[3]][[i]] <- TCGA_PAN_km[[3]][[i]] %>% list()
    TCGA_PAN_km[[3]][[i]] <- summary_km$table %>% as.data.frame()
  }
  return(TCGA_PAN_km)
}

# TCGA_sig_list$km <- lapply(TCGA_sig_list$data, TCGA_km)
