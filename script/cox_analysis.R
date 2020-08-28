library("survival")
library("survminer")
clinical <- readRDS("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/TCGA_cilinical_info.rds")
OS_time <- clinical[,c("bcr_patient_barcode","OS" ,"OS.time"  )]

# 先针对TCGA SBS数据
TCGA_WES_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
TCGA_WES_ID <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv")
TCGA_WES_sig <- dplyr::full_join(TCGA_WES_SBS, TCGA_WES_ID[,c(2,4:20)], by = "Sample.Names")
TCGA_WES_sig[is.na(TCGA_WES_sig)] <- 0
save(TCGA_WES_sig, file = "TCGA_WES_sig.rds")

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

# 针对PCAWG
PCAWG_clinical <- readxl::read_xlsx("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/PCAWG-cli/pcawg_donor_clinical_August2016_v9.xlsx")
clinical_TCGA_name <- PCAWG_clinical$submitted_donor_id[which(substr(PCAWG_clinical$submitted_donor_id, 1,4) == "TCGA")]
PCAWG_clinical2 <- PCAWG_clinical[which(!(PCAWG_clinical$submitted_donor_id %in% clinical_TCGA_name)),] %>% subset(donor_survival_time!="NA")
PCAWG_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWG_DBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
PCAWG_ID <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_SigProfiler_ID_signatures_in_samples.csv")
PCAWG_sig <- dplyr::full_join(PCAWG_SBS, PCAWG_DBS[,c(2,4:14)], by = "Sample.Names")
PCAWG_sig <- dplyr::full_join(PCAWG_sig, PCAWG_ID[,c(2,4:20)], by = "Sample.Names")
PCAWG_sig[is.na(PCAWG_sig)] <- 0
save(PCAWG_sig, file = "PCAWG_sig.rds")

PCAWG_CT <- PCAWG_sig$Cancer.Types[!duplicated(PCAWG_sig$Cancer.Types)]

# 样本名不匹配，匹配一下
PCAWG_sample_sheet <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/PCAWG-cli/pcawg_sample_sheet.tsv",
                               sep = "\t")
match_sample <- PCAWG_sample_sheet[ ,c("icgc_donor_id", "icgc_specimen_id")]

PCAWG_OS_time <- PCAWG_clinical2[, c("icgc_donor_id","donor_vital_status",
                                     "donor_survival_time")]
PCAWG_OS_time <- dplyr::inner_join(match_sample, PCAWG_OS_time, by = "icgc_donor_id") %>%
  na.omit()
names(PCAWG_OS_time)[2] = "Sample.Names"

PCAWG_OS_time$donor_vital_status <- ifelse(PCAWG_OS_time$donor_vital_status == "alive", 0,1)

PCAWG_cox_data <- dplyr::inner_join(PCAWG_OS_time[2:4],PCAWG_sig, by = "Sample.Names")

# 单因素cox分析
PCAWG_sig_list <- PCAWG_cox_data %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()

PCAWG_single_cox <- function(data){
  for (i in sig_name){
    f <- as.formula(paste("Surv(donor_survival_time, donor_vital_status) ~ ", paste(i)))
    res.cox <- coxph(f, data = as.data.frame(data))
    summary_cox <- summary(res.cox)
    TCGA_PAN_cox[i, ] <- c( "PCAWG", "WGS",
                            summary_cox$conf.int[1],
                            summary_cox$coefficients[5],
                            summary_cox$conf.int[3],
                            summary_cox$conf.int[4],
                            i)
    TCGA_PAN_cox <- TCGA_PAN_cox %>% na.omit()
  }
  return(TCGA_PAN_cox)
}

PCAWG_sig_list$cox <- lapply(PCAWG_sig_list$data, PCAWG_single_cox)

# 多因素cox分析
# res.cox_test <- coxph(f, data = test_data, control = coxph.control(iter.max = 50))


# 针对nonPCAWG
nonPCAWG_WES_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WES_SBS$scale <- "WES"
nonPCAWG_WGS_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WGS_SBS$scale <- "WGS"
nonPCAWG_sig <- rbind(nonPCAWG_WES_SBS,nonPCAWG_WGS_SBS)
save(nonPCAWG_sig, file = "nonPCAWG_sig.rds")


# PCAWG_SBS <- read.csv("C:/Users/lenovo/Documents/GitHub/PAWSASP/data/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")



# 多因素cox分析有一些问题，暂时先不考虑
# (f <- as.formula(paste("Surv(OS.time, OS) ~ ", paste(test$var, collapse = "+"))))
# test_data <- TCGA_sig_list$data[6] %>% as.data.frame()
# res.cox_test <- coxph(Surv(OS.time, OS) ~ SBS1 + SBS2 + SBS3 + SBS5 + SBS9 + SBS10a + 
#                         SBS10b + SBS13 + SBS15 + SBS17a , data = test_data)
