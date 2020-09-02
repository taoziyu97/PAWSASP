## COSMIC2 signature 
cancer_type <- c('Adrenocortical carcinoma','Acute lymphoblastic leukemia', 'Acute myelogenous leukemia', 'Bladder', 'Breast', 
                 'Cervix','Chondrosarcoma', 'B-cell chronic lymphocytic leukemia', 'Colorectum', 'Glioblastoma', 
                 'Glioma Low grade','Head and neck', 'Kidney Chromophobe', 'Kidney Clear Cell', 'Kidney Papillary', 
                 'Liver','Lung Adeno', 'Lung Small Cell', 'Lung Squamous', 'Lymphoma B-cell', 
                 'Lymphoma Hodgkin','Medulloblastoma', 'Melanoma', 'Masopharyngeal Carcinoma', 'Masopharyngeal Carcinoma',
                 'Neuroblastoma','Oesophagus', 'Oral gingivo-buccal squamous', 'Osteosarcoma', 'Ovary', 
                 'Pancreas','Paraganglioma', 'Pilocytic Astrocytoma', 'Prostate', 'Stomach', 
                 'Thyroid','Urothelial Carcinoma', 'Uterine Carcinoma', 'Uterine Carcinosarcoma','Uveal Melanoma')

# abbreviation <- c('ACC','ALL','AML','BLCA','BRCA',
#                   )
# ACC：Adrenocortical carcinoma
# ALL：Acute lymphoblastic leukemia
# AML：Acute myelogenous leukemia
# BLCA：Bladder Urothelial Carcinoma
# BRCA：Breast invasive carcinoma

# UCC：uterine cervix cancer（宫颈癌）

# CLL：B-cell chronic lymphocytic leukemia

cosmic2 <- matrix(data = NA, nrow = 40, ncol = 2, 
                  dimnames = list(1:40,
                                  c("full name","signature type"))) %>% as.data.frame()
cosmic2[,1] <- cancer_type
sig <- c('1,2,4,5,6,13,18','1,2,5,13', '1,5', '1,2,5,10,13', '1,2,3,5,6,8,10,13,17,18,20,26,30',
         '1,2,5,6,10,13,26','1,5', '1,2,5,9,13', '1,5,6,10', '1,5,11',
         '1,5,6,14', '1,2,4,5,7,13','1,5,6','1,5,6,27','1,2,5,13',
         '1,4,5,6,12,16,17,22,23,24','1,2,4,5,6,13,17','1,4,5,15','1,2,4,5,13','1,2,5,9,13,17',
         '1,5,25','1,5,8','1,5,7,11,17','1,2,5,13','1,2,5,6,13',
         '1,5,18','1,2,4,5,6,13,17','1,2,5,7,13,29','1,2,5,6,13,30','1,3,5',
         '1,2,3,5,6,13','1,5','1,5,19','1,5,6','1,2,5,13,15,17,18,20,21,26,28',
         '1,2,5,13','1,2,5,13,22','1,2,5,6,10,13,14,26','1,2,5,6,10,13','1,5,6')
cosmic2[,2] <- sig
# save data
save(cosmic2, file = paste0("~/PAWSASP/data/cosmic2_signature.RData"))

## 测试提取各个cancer的signature
# test <- strsplit(cosmic2$`signature type`,split = ",")


## 整合TCGA的数据
# TCGA cancer type
TCGA_WES_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
TCGA_WES_ID <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv")
TCGA_CT <- TCGA_WES_SBS$Cancer.Types[!duplicated(TCGA_WES_SBS$Cancer.Types)]
# SBS
TCGA_SBS_list <- TCGA_WES_SBS %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
SBS_list <- lapply(TCGA_SBS_list$data, function(data){
  sum <- apply(data[,3:67], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "SBS", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"SBS list")
TCGA_SBS <- cbind('cancer type' = TCGA_SBS_list$Cancer.Types, 
                  'source' = "TCGA",
                  'scale' = "WES",
                  SBS_list)

# ID
TCGA_ID_list <- TCGA_WES_ID %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
ID_list <- lapply(TCGA_ID_list$data, function(data){
  sum <- apply(data[,3:19], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "ID", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"ID list")
TCGA_ID <- cbind('cancer type' = TCGA_ID_list$Cancer.Types,
                 ID_list)
TCGA_sig <- cbind(TCGA_SBS, 'DBS list' = NA)
TCGA_sig <- dplyr::full_join(TCGA_sig, TCGA_ID, by = 'cancer type')
save(TCGA_sig, file = "~/PAWSASP/data/TCGA_signature.RData")


## 整合PCAWG的数据
PCAWG_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWG_DBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
PCAWG_ID <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_SigProfiler_ID_signatures_in_samples.csv")
PCAWG_CT <- PCAWG_SBS$Cancer.Types[!duplicated(PCAWG_SBS$Cancer.Types)]
# SBS
PCAWG_SBS_list <- PCAWG_SBS %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
SBS_list2 <- lapply(PCAWG_SBS_list$data, function(data){
  sum <- apply(data[,3:67], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "SBS", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"SBS list")
PCAWG_SBS <- cbind('cancer type' = PCAWG_SBS_list$Cancer.Types, 
                   'source' = "PCAWG",
                   'scale' = "WGS",
                   SBS_list2)
# DBS
PCAWG_DBS_list <- PCAWG_DBS %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
DBS_list2 <- lapply(PCAWG_DBS_list$data, function(data){
  sum <- apply(data[,3:13], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "DBS", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"DBS list")
PCAWG_DBS <- cbind('cancer type' = PCAWG_DBS_list$Cancer.Types, 
                   DBS_list2)
# ID
PCAWG_ID_list <- PCAWG_ID %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
ID_list2 <- lapply(PCAWG_ID_list$data, function(data){
  sum <- apply(data[,3:19], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "ID", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"ID list")
PCAWG_ID <- cbind('cancer type' = PCAWG_ID_list$Cancer.Types, 
                  ID_list2)
PCAWG_sig <- dplyr::full_join(PCAWG_SBS, PCAWG_DBS, by = 'cancer type')
PCAWG_sig<- dplyr::full_join(PCAWG_sig, PCAWG_ID, by = 'cancer type')
save(PCAWG_sig, file = "~/PAWSASP/data/PCAWG_signature.RData")

## 整合non-PCAWG的数据
nonPCAWG_WES_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WGS_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WES_CT <- nonPCAWG_WES_SBS$Cancer.Types[!duplicated(nonPCAWG_WES_SBS$Cancer.Types)]
nonPCAWG_WGS_CT <- nonPCAWG_WGS_SBS$Cancer.Types[!duplicated(nonPCAWG_WGS_SBS$Cancer.Types)]
nonPCAWG_CT <- union(nonPCAWG_WES_CT, nonPCAWG_WGS_CT)

# WES_SBS
nonPCAWG_WES_SBS_list <- nonPCAWG_WES_SBS %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
SBS_list2 <- lapply(nonPCAWG_WES_SBS_list$data, function(data){
  sum <- apply(data[,3:67], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "SBS", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"SBS list")
nonPCAWG_WES_SBS <- cbind('cancer type' = nonPCAWG_WES_SBS_list$Cancer.Types, 
                          'source' = "nonPCAWG",
                          'scale' = "WES",
                          SBS_list2,
                          'DBS list' = NA,
                          'ID list' = NA)
# WGS_SBS
nonPCAWG_WGS_SBS_list <- nonPCAWG_WGS_SBS %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
SBS_list3 <- lapply(nonPCAWG_WGS_SBS_list$data, function(data){
  sum <- apply(data[,3:67], 2, sum)
  paste_sig <- ifelse(sum ==0, NA,
                      gsub( "SBS", "",names(sum))) %>%
    na.omit() %>%
    paste(., collapse = ",")
}) %>% 
  unlist() %>%
  as.data.frame() %>% 
  `names<-`(.,"SBS list")
nonPCAWG_WGS_SBS <- cbind('cancer type' = nonPCAWG_WGS_SBS_list$Cancer.Types, 
                          'source' = "nonPCAWG",
                          'scale' = "WGS",
                          SBS_list3,
                          'DBS list' = NA,
                          'ID list' = NA)

nonPCAWG_sig <- rbind(nonPCAWG_WES_SBS, nonPCAWG_WGS_SBS)
save(nonPCAWG_sig, file = "~/PAWSASP/data/nonPCAWG_signature.RData")

## test 针对AML
test <- TCGA_SBS_list[[2]][[1]]
test_sum <- apply(test[,3:67], 2, sum)
a <- ifelse(test_sum ==0, NA, 
            gsub( "SBS", "",names(test_sum))) %>%
  na.omit() %>% 
  paste(., collapse = ",")

