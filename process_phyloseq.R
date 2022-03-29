library(phyloseq)
library(fantaxtic)
library(tidyverse)
library(data.table)
mypath="C:/Users/cs001/Desktop/paper/0731gut/new_152_sample/yyftaxonomy/"
load(paste0(mypath,"ps.ng.tax0406rela_152.Rdata"))
for(nnum in c(300)){
  for(thresholdi in c(0.6)){
    
    phylo_Glom=fantaxtic::get_top_taxa(ps.ng.tax_rela, n=nnum, relative = TRUE, discard_other = T)
    meta=data.frame(sample_data(phylo_Glom))
    meta$Subgroup=paste0(meta$Group,'_',meta$time)
    meta$Row.names=paste0('X',meta$seq_num)
    meta2=meta[,c("Row.names","HFC_dif")]
    
    otutable=otu_table(phylo_Glom)
    otutable[otutable==0]=NA
    otutableT <- as.data.frame(t(otutable))
    tablegroup_OTUlist <- split(otutableT,f=meta$Subgroup)#taxomy t()
    tabletype_list<-lapply(tablegroup_OTUlist,function(x)
    {
      ############ɸ???
      subject_frame_num<-x
      subject_frame_num[subject_frame_num==0]<-NA
      subject_frame_num<-subject_frame_num[,which(colMeans(is.na(subject_frame_num))<thresholdi)]
      ############
      return(subject_frame_num)
    })
    colname_list<-lapply(tabletype_list,function(x){colnames(x)})
    #name<-unique(table$Subgroup)
    com_rowname<-Reduce(intersect, colname_list)
 

    tablecom_list=list()
    #write each group
    for (groupname in names(tabletype_list)){
      #print(groupname)
      SSN20=tabletype_list[[groupname]]#x20
      SSN20<-SSN20[,com_rowname]
      #write.table(t(x),paste0(mypath,'S16_',groupname,'_',nnum,'_',thresholdi,'_count.txt'),row.names = T,sep = "\t")#group
      tablecom_list[[groupname]]=SSN20%>%rownames_to_column(var='samples')
      #write.table(t(SSN20),paste0(mypath,'S16_',groupname,'_',nnum,'_',thresholdi,"_groupad.txt"),row.names = T,sep = "\t")
      
      combname19<-combn(rownames(SSN20),nrow(SSN20)-1)
      combframe19_list<-apply(combname19,2,function(samples19) {
        samplename<-setdiff(rownames(SSN20),samples19)#-1th sample
        SSN19=SSN20[samples19,]
        SSN1=SSN20[samplename,]
        
        combname18<-combn(rownames(SSN19),nrow(SSN19)-1)
        combframe18_list<-apply(combname18,2,function(samples18) {
          samplenamein<-setdiff(rownames(SSN19),samples18)#-1th sample
          SSN18=SSN19[samples18,]
          SSN1in=SSN19[samplenamein,]
          
          #inpute mean with 18 samples
          for(i in 1:ncol(SSN18)){
            colmeani18=mean(SSN18[,i], na.rm = TRUE)
            SSN18[is.na(SSN18[,i]), i] <- colmeani18#fill NA with mean
            SSN1in[i][is.na(SSN1in[i])]=colmeani18
          }
          
          SSN18add1=data.table::rbindlist(list(SSN18,SSN1in),fill=TRUE)
          
          SSN18add1=t(SSN18add1)
          colnames(SSN18add1)=c(row.names(SSN18),samplenamein)
          #%>%rownames_to_column(var='Row.names')
          #%>%left_join(meta2,by= 'Row.names')%>%column_to_rownames(var='Row.names')
          
          #write.table(SSN18add1, file=paste0(mypath,groupname,'_',samplenamein,'_',nnum,'_',thresholdi,'_in_',samplename,"_singleadd.txt"),row.names = T,sep = "\t")
        })
        #inpute mean with 19 samples
        for(i in 1:ncol(SSN19)){
          colmeani=mean(SSN19[,i], na.rm = TRUE)
          SSN19[is.na(SSN19[,i]), i] <- colmeani#fill NA with mean
          SSN1[i][is.na(SSN1[i])]=colmeani
        }
        
        SSN19add1=data.table::rbindlist(list(SSN19,SSN1),fill=TRUE)
        
        SSN19add1=t(SSN19add1)
        colnames(SSN19add1)=c(row.names(SSN19),samplename)
        #%>%rownames_to_column(var='Row.names')
        #%>%left_join(meta2,by= 'Row.names')%>%column_to_rownames(var='Row.names')
        
        write.table(SSN19add1, file=paste0(mypath,groupname,'_',samplename,'_',nnum,'_',thresholdi,"_singlemeta16s0909.txt"),row.names = T,sep = "\t")
      })
      
      bind_table=rbindlist(tablecom_list,fill = F,use.names = T,idcol = 'Subgroup')
      colnum=ncol(bind_table)
      #write.table(bind_table,paste0(mypath,'S16_',nnum,'_',colnum,'_',thresholdi,'_all4add.txt'),row.names = F,sep = "\t")
    }
  }
}