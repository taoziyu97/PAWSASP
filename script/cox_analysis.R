library("survival")
library("survminer")
clinical <- readRDS("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/TCGA_cilinical_info.rds")
OS_time <- clinical[,c("bcr_patient_barcode","OS" ,"OS.time"  )]

# 先针对TCGA SBS数据
TCGA_WES_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
TCGA_WES_ID <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv")
TCGA_WES_sig <- dplyr::full_join(TCGA_WES_SBS, TCGA_WES_ID[,c(2,4:20)], by = "Sample.Names")
TCGA_WES_sig[is.na(TCGA_WES_sig)] <- 0

TCGA_WES_sig$Sample.Names <- substr(TCGA_WES_sig$Sample.Names, 1 ,12)
names(TCGA_WES_sig)[2] <- "bcr_patient_barcode"
TCGA_cox_data <- dplyr::inner_join(OS_time,TCGA_WES_sig, by = "bcr_patient_barcode")

## 单个signature
sig_name <- names(TCGA_WES_sig)[4:85]
## 去除人工引入的signature
load("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/artefacts.Rdata")
sig_name <- subset(sig_name, !(sig_name %in% c(artefacts)))

## 从做cox分析开始，把代码分成数据预处理和cox分析

TCGA_PAN_cox <- matrix(data = NA, nrow = 64, ncol = 7, 
                       dimnames = list(c(sig_name),
                                       c("source","scale",
                                         "HR","P value","CI upper", 
                                         "CI lower", "var"))) %>%
  as.data.frame()

# 针对所有癌症类型
TCGA_sig_list <- TCGA_cox_data %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()

single_cox <- function(data){
  for (i in sig_name){
    f <- as.formula(paste("Surv(OS.time, OS) ~ ", paste(i)))
    res.cox <- coxph(f, data = as.data.frame(data))
    summary_cox <- summary(res.cox)
    TCGA_PAN_cox[i, ] <- c( "TCGA", "WES",
                           summary_cox$conf.int[1],
                           summary_cox$coefficients[5],
                           summary_cox$conf.int[3],
                           summary_cox$conf.int[4],
                           i)
    TCGA_PAN_cox <- TCGA_PAN_cox %>% na.omit()
  }
  return(TCGA_PAN_cox)
}

TCGA_sig_list$cox <- lapply(TCGA_sig_list$data, single_cox)


# 多因素cox分析有一些问题，暂时先不考虑
# (f <- as.formula(paste("Surv(OS.time, OS) ~ ", paste(test$var, collapse = "+"))))
# test_data <- TCGA_sig_list$data[6] %>% as.data.frame()
# res.cox_test <- coxph(Surv(OS.time, OS) ~ SBS1 + SBS2 + SBS3 + SBS5 + SBS9 + SBS10a + 
#                         SBS10b + SBS13 + SBS15 + SBS17a , data = test_data)
