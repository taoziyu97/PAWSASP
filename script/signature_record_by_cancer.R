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
sig <- c('1,2,4,5,6,13,18,"other"','1,2,5,13', '1,5', '1,2,5,10,13', '1,2,3,5,6,8,10,13,17,18,20,26,30,"other"',
         '1,2,5,6,10,13,26,"other"','1,5', '1,2,5,9,13', '1,5,6,10,"other"', '1,5,11,"other"',
         '1,5,6,14', '1,2,4,5,7,13','1,5,6','1,5,6,27','1,2,5,13',
         '1,4,5,6,12,16,17,22,23,24,"other"','1,2,4,5,6,13,17','1,4,5,15','1,2,4,5,13,"other"','1,2,5,9,13,17',
         '1,5,25','1,5,8','1,5,7,11,17,"other"','1,2,5,13','1,2,5,6,13',
         '1,5,18','1,2,4,5,6,13,17','1,2,5,7,13,29','1,2,5,6,13,30','1,3,5',
         '1,2,3,5,6,13,"other"','1,5','1,5,19','1,5,6','1,2,5,13,15,17,18,20,21,26,28,"other"',
         '1,2,5,13','1,2,5,13,22','1,2,5,6,10,13,14,26','1,2,5,6,10,13','1,5,6')
cosmic2[,2] <- sig
# save data
save(cosmic2, file = paste0("~/PAWSASP/data/cosmic2_signature.RData"))

## 测试提取各个cancer的signature
# test <- strsplit(cosmic2$`signature type`,split = ",")

# TCGA signature by individual tumors

TCGA_signature <- matrix(data = NA, nrow = 40, ncol = 7, 
                         dimnames = list(1:40, 
                                         c("cancer type", "full name", "source", 
                                           "sample","SBS list", "DBS list", "ID list")))
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

TCGA_ID_list <- TCGA_WES_ID %>%
  dplyr::group_by(Cancer.Types) %>%
  tidyr::nest()
# ID
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


# TCGA_SBS_list <- purrr::map_df(TCGA_WES_SBS, function(TCGA_WES_SBS){
#   TCGA_SBS_list <- TCGA_WES_SBS %>%
#     dplyr::group_by(Cancer.Types) %>%
#     tidyr::nest()
#   # list2frame <- SBS_list %>% unlist() %>% as.data.frame()
#   a <- TCGA_SBS_list %>%
#     dplyr::mutate(
#       SBS_list = (lapply(TCGA_SBS_list$data, function(data){
#         sum <- apply(data[,3:67], 2, sum)
#         paste_sig <- ifelse(sum ==0, NA, 
#                gsub( "SBS", "",names(sum))) %>%
#           na.omit() %>% 
#           paste(., collapse = ",")
#       })) %>% unlist() 
#     ) 
# }) %>%
#   data.table::as.data.table()





purrr::map_df(TCGA_SBS_list, .f = function(a){
  print(a$Cancer.Types)})

## test 针对AML
test <- TCGA_SBS_list[[2]][[1]]
test_sum <- apply(test[,3:67], 2, sum)
a <- ifelse(test_sum ==0, NA, 
            gsub( "SBS", "",names(test_sum))) %>%
  na.omit() %>% 
  paste(., collapse = ",")

