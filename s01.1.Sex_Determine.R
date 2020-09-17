tmp_dir <- "D:/TangLab_Project/Human_Gonad/Sex_Determine"
files <- list.files(path=tmp_dir,pattern="chr_cor.txt")
## include chrM
tmp_data <- read.table(paste(tmp_dir,"HE10W_E1.chr_cor.txt",sep="/"))
tmp_data$Ratio <- round(tmp_data$V3/sum(tmp_data$V3)*100,2)
chr_Ratio <- data.frame(HE10W_E1=tmp_data$Ratio)
rownames(chr_Ratio) <- tmp_data$V1 

## exclude chrM
tmp_data_noM <- tmp_data[tmp_data$V1!="MT",]
tmp_data_noM$Ratio <- round(tmp_data_noM$V3/sum(tmp_data_noM$V3)*100,2)
chr_Ratio_nochrM <- data.frame(HE10W_E1=tmp_data_noM$Ratio)
rownames(chr_Ratio_nochrM) <- tmp_data_noM$V1 
  
for(i in 2:length(files)){
  tmp_sp <- gsub(".chr_cor.txt","",files[i])
  tmp_data <- read.table(paste(tmp_dir,files[i],sep="/"))
  tmp_data$Ratio <- round(tmp_data$V3/sum(tmp_data$V3)*100,2)
  chr_Ratio[,tmp_sp] <- tmp_data$Ratio
    
  tmp_data_noM <- tmp_data[tmp_data$V1!="MT",]
  tmp_data_noM$Ratio <- round(tmp_data_noM$V3/sum(tmp_data_noM$V3)*100,2)
  chr_Ratio_nochrM[,tmp_sp] <- tmp_data_noM$Ratio 
}
  
write.table(file="chr_Ratio_nochrM.txt",chr_Ratio_nochrM,quote = F,sep = "\t")
write.table(file="chr_Ratio.txt",chr_Ratio,quote = F,sep = "\t")
  
chr_Ratio$Chr <- factor(rownames(chr_Ratio),levels = c(1:22,"X","Y","MT"),ordered = T)
melt_chr_Ratio <- melt(chr_Ratio,id.vars = "Chr",variable.name = "Week",value.name = "Ratio")
  
sample_order <- c("HE6W_E1","HE7W_E1","HE7W_ME1","HE7W_FE1","HE8W_E1","HE9W_E1","HE10W_E1" ,"HE13W_FE1","HE15W_ME1_rep1","HE15W_ME1_rep2","HE16W_FE1","HE18W_FE1","HE19W_ME1","HE22W_FE1","HE22W_FE2","HE23W_ME1","HE23W_ME2")

Type2Sex <- c("Male","Female","Male","Female","Male","Male","Male","Female","Male","Male","Female","Female","Male","Female","Female","Male","Male")  
names(Type2Sex) <- sample_order
melt_chr_Ratio$Sex <- as.character(Type2Sex[melt_chr_Ratio$Week])
  
melt_chr_Ratio$Week <- factor(melt_chr_Ratio$Week,levels = sample_order,ordered = T)

pdf("Sex_chromosome_reads_Ratio.pdf",height = 6,width = 8)
ggplot(data=melt_chr_Ratio[melt_chr_Ratio$Chr%in%c("X","Y"),],aes(x=Week,y=Ratio,fill=Sex))+geom_bar(stat = "identity")+facet_wrap(~Chr,nrow = 5,scales = "free_y")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("lightpink","skyblue"))+ggtitle("With chrM")
dev.off()


pdf("Sex_chromosome_reads_Ratio_nochrM.pdf",height =6,width = 8)
ggplot(data=melt_chr_Ratio_nochrM[melt_chr_Ratio_nochrM$Chr%in%c("X","Y"),],aes(x=Week,y=Ratio,fill=Sex))+geom_bar(stat = "identity")+facet_wrap(~Chr,nrow = 5,scales = "free_y")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("lightpink","skyblue"))+ggtitle("NO chrM")
dev.off()


