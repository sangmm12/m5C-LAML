rm(list=ls())
file_dir = "~/diff/m5c-laml"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
need_list = c('')	#This is the group document

for (file_name in need_list)
{

  tumor_list = c('LAML')	#Use this to change what type of tumor you want to analyse
  
  final_tumor_list = c()
  
  my_data = read.csv('tcgaexp',header=T,check.names=F)	#Expression matrix document
  
  other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
  
    dat <- data.frame(check.names = F)
    
    
    other_file <- other_data[grep(name,other_data$CODE),]
    
    if (length(rownames(other_file))==0)
    {
        next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    exp_file <- my_data[grep(name,my_data$CODE),]
    exp_file <- exp_file[grep("Tumor",exp_file$Group),]
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    exp_file <- subset(exp_file,select=-c(Group,CODE,SampleName))
    gene_list <- colnames(exp_file)
    all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))
    if (length(all_name)==0)
    {
        next
    }
    
    for(gene_name in gene_list)
    {
        for(i in all_name)
        {
            dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
        }
    
    }
    colnames(dat) <- c("Gene","Group","value")
    dat[,3] = as.numeric(dat[,3])
    dat[,3] = log(dat[,3]+1)
    dat <- na.omit(dat)
    dat <- dat[dat[,1]%in%temp_name,]
    final_tumor_list <- append(final_tumor_list,name)
    print(name)
    pdf(paste(file_name,"-",name,".pdf",sep=''),width=length(unique(dat[,1])),height = 8)
    p <- ggboxplot(dat, x = "Gene", y = "value",
                   color = "Group", palette = 'jama',
                   add = "jitter",x.text.angle=60)
    p <- p + xlab("")+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
    print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    dev.off()
  }
}