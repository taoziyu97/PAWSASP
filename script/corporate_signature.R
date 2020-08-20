# 整合WES SBS数据
TCGA_WES_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
TCGA_WES_SBS$scale <- "WES"
TCGA_WES_SBS$source <- "TCGA"
nonPCAWG_WES_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WES_SBS$scale <- "WES"
nonPCAWG_WES_SBS$source <- "TCGA"
pancancer_WES_SBS <- rbind(TCGA_WES_SBS, nonPCAWG_WES_SBS)
save(pancancer_WES_SBS, file = "./data/pancancer_WES_SBS.Rdata")
# 整合WGS SBS数据
nonPCAWG_WGS_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv")
nonPCAWG_WGS_SBS$scale <- "WGS"
nonPCAWG_WGS_SBS$source <- "nonPCAWG"
PCAWG_WGS_SBS <- read.csv("/home/tzy/PAWSASP/data/ICGC-PCAWG-TCGA/SP_Signatures_in_Samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
PCAWG_WGS_SBS$scale <- "WGS"
PCAWG_WGS_SBS$source <- "PCAWG"
pancancer_WGS_SBS <- rbind(nonPCAWG_WGS_SBS, PCAWG_WGS_SBS)
save(pancancer_WGS_SBS, file = "./data/pancancer_WGS_SBS.Rdata")
# 癌症类型转换
cancer.type <- c("LAML", "ACC", "BLCA", "LGG", "BRCA",
                 "CESC", "CHOL", "LCML", "COAD", "CNTL",
                 "ESCA", "FPPP", "GBM", "HNSC", "KICH",
                 "KIRC", "KIRP", "LIHC", "LUAD", "LUSC",
                 "DLBC", "MESO", "MISC", "OV", "PAAD",
                 "PCPG", "PRAD", "READ", "SARC", "SKCM",
                 "STAD", "TGCT", "THYM", "THCA", "UCS",
                 "UCEC", "UVM")
save(cancer.type, file = "./data/cancer_type.Rdata")

# 存储一下artificial的signature，方便之后去除
artefacts <- c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", 
               "SBS48", "SBS49", "SBS50", "SBS51", "SBS52",
               "SBS53", "SBS54", "SBS55", "SBS56", "SBS57",
               "SBS58", "SBS59", "SBS60")
save(artefacts, file = "./data/artefacts.Rdata")

# 构建存储各个cancer type拥有的signature的矩阵，有就记为1，没有就记为0
# SBS共有47个，DBS有11个，ID有17个
sig.by.cancer <- matrix(NA, nrow = 75, ncol = 33) %>% data.table::as.data.table()
sig.by.cancer[,1] <- c(paste0(rep("SBS", 47), 
                              c(as.character(1:6),"7a","7b","7c","7d",
                                "8","9","10a","10b",as.character(11:16),
                                "17a","17b",as.character(18:26),
                                as.character(28:42),"44")),
                           paste0(rep("DBS", 11),
                                  c(as.character(1:11))),
                           paste0(rep("ID", 17),
                                  c(as.character(1:17))))
sig.by.cancer <- sig.by.cancer[, V2 := 
                                 ifelse(V1 %in% c(paste0(rep("SBS",17),c("1","2","3","5","9",
                                                         "12","13","15","17a","17b","18",
                                                         "21","22","24","32","40","44")),
                                                  paste0(rep("DBS",2),c("2","4")),
                                                  paste0(rep("ID",6),
                                                         c("1","2","3","5","6","8"))), 1, NA)]
sig.by.cancer <- sig.by.cancer[, V3 := 
                                 ifelse(V1 %in% c(paste0(rep("SBS",7),c("1","2","5","8","13","29","40")),
                                                  paste0(rep("DBS",3),c("2","4","11")),
                                                  paste0(rep("ID",8),c("1","2","3","4","5","8","9","10"))), 1, NA)]

sig.by.cancer <- sig.by.cancer[, V4 := 
                                 ifelse(V1 %in% c(paste0(rep("SBS",10),c("1","2","3","5","8","13",
                                                                         "17a","17b","30","40")),
                                                  paste0(rep("ID",7),c("1","2","3","4","5","8","9"))),1,NA)]

sig.by.cancer <- sig.by.cancer[, V5 := 
                                 ifelse(V1 %in% c(paste0(rep("SBS",13),c("1","2","3","5","8","9",
                                                                         "13","17a","17b","18","37",
                                                                         "40","41")),
                                                  paste0(rep("DBS",5),c("2","4","6","9","11")),
                                                  paste0(rep("ID",8),c("1","2","4","5","6","8","9","11"))),1,NA)]

sig.by.cancer <- sig.by.cancer[, V6 :=
                                 ifelse(V1 %in% c(paste0(rep("SBS",6),c("1","2","5","13","18","40")),
                                                  paste0(rep("ID",5),c("1","2","6","9","10"))),1,NA)]

sig.by.cancer <- sig.by.cancer[, V7 :=
                                 ifelse(V1 %in% c(paste0(rep("SBS",5),c("1","5","11","30","40")),
                                                  paste0(rep("DBS",3),c("2","4","11")),
                                                  paste0(rep("ID",6),c("1","2","4","5","8","9"))),1,NA)]

sig.by.cancer <- sig.by.cancer[, V8 :=
                                 ifelse(V1 %in% c(paste0(rep("SBS",6),c("1","5","8","18","39","40")),
                                                  paste0(rep("ID",6),c("1","2","5","8","9","10"))),1,NA)]

sig.by.cancer <- sig.by.cancer[, V9 :=
                                 ifelse(V1 %in% c(paste0(rep("SBS",4),c("1","5","8","40")),
                                                  paste0(rep("ID",4),c("1","2","4","9"))),1,NA)]

# 假设现在我需要提取第一种cancer type的所有signature，则可以从这个表格中这样提取
a <- ifelse(sig.by.cancer$V2 == 1, sig.by.cancer$V1,"")%>% na.omit()

