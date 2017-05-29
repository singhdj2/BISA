library(readxl)
PUS_layout <- read_excel("~/Dropbox/1.Research_and_Education/Kansas_State_University_2014-/Doctoral_Research/Projects_data/UAV_pipeline_work/Data_UAV/image_data/India_2017/raw_data/field_layouts/2016-17_cleaned_GS Pusa_layout.xlsx", 
                         col_names = FALSE)
PUS_layout <- PUS_layout[-nrow(PUS_layout),]

# import crosspoint df
crosspoint <- read_csv("~/Dropbox/Kevin-Dj_orthos_2017/crosspoint.csv")
View(crosspoint)

crosspoint <- as.data.frame(crosspoint)
b <- substr(crosspoint$WKT, 8,nchar(crosspoint$WKT)-1)
c <- unlist(strsplit(b, " "))
df <- data.frame(matrix(unlist(c), nrow=nrow(crosspoint), byrow=T)) # nrow=2035 here
crosspoint <- cbind(crosspoint, df)
colnames(crosspoint)[4:5] <- c("UTMx", "UTMy")
cpoint <- crosspoint[,c(4,5,2,3)] #re-order df

cpoint$UTMx <- as.character(cpoint$UTMx)
cpoint$UTMy <- as.character(cpoint$UTMy)

write.csv(cpoint, file = "Dropbox/Kevin-Dj_orthos_2017/cleaned_crosspoint.csv",row.names = F)
