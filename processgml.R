library(igraph)
mypath="C:/Users/cs001/Desktop/paper/0731gut/new_152_sample/yyftaxonomy/"
filelist=list.files(pattern = '.*._comselP.csv$',mypath)
#filelist=filelist[!grepl(pattern = '*.M0|*.bad',x = filelist)]

table_list=list()
var_list=list()
SSN_list=list()
for(file19 in filelist){#(pattern = '*.good.*._comsel.csv.csv$|*.bad.*._comsel.csv.csv$',mypath)){
  filesplit=str_split_fixed(file19,pattern='_',n=7)
  
  G1<-read.csv(paste0(mypath,file19))
  G1$filename=file19
  #edgeframe_V1=merge(G1,taxonomy[c('X','g_')],by.x='index',by.y='X')
  #names(edgeframe_V1)[names(edgeframe_V1)=="g_"] <- "V1_s"
  #edgeframe_V12=merge(edgeframe_V1,taxonomy[c('X','g_')],by.x='variable',by.y='X')
  #names(edgeframe_V12)[names(edgeframe_V12)=="g_"] <- "V2_s"
  #freq_table=table(paste0(edgeframe_V12$V1_g,'_',edgeframe_V12$V2_g))
  G1$edges_OTU=paste0(G1$index,'_',G1$variable)
  #edgeframe_V12=G1
  #edgeframe_V12$edges_OTU=paste0(edgeframe_V12$index,'__',edgeframe_V12$variable)
  SSN19add1=data.frame(reshape2::acast(data =G1,formula = samples~edges_OTU,value.var = 'value'),check.names = F)
  
  group=filesplit[1,1]
  time=filesplit[1,2]
  samplename=filesplit[1,3]
  topn=filesplit[1,4]
  thre=filesplit[1,5]
  
  if(samplename %in% rownames(SSN19add1)){
    edgelist=SSN19add1[samplename,]
    SSN_list[[file19]]=edgelist
    edgelist=data.frame(t(edgelist))%>%
      rownames_to_column(var='edges_OTU')%>%
      separate(edges_OTU,into=c('V1','V2'),sep = '_')
    edgelist=edgelist[!is.na(edgelist[samplename]),]
    g=graph_from_edgelist(as.matrix(edgelist[,c('V1','V2')]),directed = F)
    E(g)$weight=unlist(edgelist[samplename])
    E(g)$weight=format(E(g)$weight,digits=3)
    write_graph(g, file = paste0(mypath,samplename,'_',thre,'_',topn,'_05_base.gml'), "gml")
  } else {
    print(paste0('No file',mypath,samplename,'_',thre,'_',topn,'_05_base.gml'))
  }
}
SSNall<-data.table::rbindlist(SSN_list,fill = T,use.names = T,idcol = 'samplename')%>%melt(na.rm = T)
write.csv(SSNall,paste0(mypath,"comsel_SSN_base.csv"),row.names = F)