options(stringsAsFactors = F)
setwd("C:\\analysis_work/ZIhai_li/cellmarker/")
# read all cell marker data in all.data
all.data<-read.delim("Mouse_cell_markers.txt")
Cell.name<-unique(all.data$cellName)
Cell.type<-unique(all.data$tissueType)
##########################################################
# this part was used to construct cell marker list
########################################################
# find interested cell type

all.related<-unique(all.data$cellName)
# find cell type specific marker gene

Construct.cellmarker<-function(CellNameList=T.cell.related,AllData=all.data){
  Cell.Marker.list<-list()
  cell.type.name<-c()
  for (i in 1:length(CellNameList)){
    #cell.type.index<-grep(paste0("^",CellNameList[i],"$"),AllData$cellName)
    cell.type.index<-which(AllData$cellName==CellNameList[i])
    selected.celltype.data<-AllData[cell.type.index,]
    celltype.specific.marker<-selected.celltype.data$geneSymbol
    
    # process format of celltype.specific.marker
    # remove bracket
    tmp.marker.removed.bracket<-unlist(strsplit(celltype.specific.marker,"[[]|[]]|,"))
    
    # remove leading space (" Cd4") --> ("Cd4")
    tmp.marker.removed.leadingspace<-gsub(" ","",tmp.marker.removed.bracket) 
    # remove NA
    tmp.marker.removed.NA<- gsub("NA","",tmp.marker.removed.leadingspace)
    if(length(which(tmp.marker.removed.NA==""))==0){
      # remove empty  data (""," ")
      tmp.marker.removed.empty<-tmp.marker.removed.NA
    }else{
      # remove empty  data (""," ")
      tmp.marker.removed.empty<-tmp.marker.removed.NA[-which(tmp.marker.removed.NA=="")]
    }
   
    # remove depulicate
    if(length(tmp.marker.removed.empty)==0){
      tmp.marker.removed.dup<-tmp.marker.removed.empty
    }else{
      tmp.marker.removed.dup<-tmp.marker.removed.empty[!duplicated(tmp.marker.removed.empty)]
    }
    # remove family
    # check whether has "family" or "class" character
    family.class.index<-grep("family|class",ignore.case = T,tmp.marker.removed.dup)
    if(length(family.class.index)!=0){
      tmp.marker.removed.family.class<-tmp.marker.removed.dup[-family.class.index]
      if(length(tmp.marker.removed.family.class)==0){tmp.marker.removed.family.class<-NULL}
    }else{tmp.marker.removed.family.class<-tmp.marker.removed.dup}
    # final standarded marker
    tmp.standard.marker<-tmp.marker.removed.family.class
    Cell.Marker.list[[i]]<-tmp.standard.marker
    cell.type.name<-c(cell.type.name, CellNameList[i])
  }
  names(Cell.Marker.list)<-cell.type.name
  return(Cell.Marker.list)
}
Construct.Tissue.marker<-function(TissueNameList=Cell.type,AllData=all.data){
  Cell.Marker.list<-list()
  cell.type.name<-c()
  for (i in 1:length(CellNameList)){
    #cell.type.index<-grep(paste0("^",CellNameList[i],"$"),AllData$cellName)
    cell.type.index<-which(AllData$tissueType==CellNameList[i])
    selected.celltype.data<-AllData[cell.type.index,]
    celltype.specific.marker<-selected.celltype.data$geneSymbol
    
    # process format of celltype.specific.marker
    # remove bracket
    tmp.marker.removed.bracket<-unlist(strsplit(celltype.specific.marker,"[[]|[]]|,"))
    
    # remove leading space (" Cd4") --> ("Cd4")
    tmp.marker.removed.leadingspace<-gsub(" ","",tmp.marker.removed.bracket) 
    # remove NA
    tmp.marker.removed.NA<- gsub("NA","",tmp.marker.removed.leadingspace)
    if(length(which(tmp.marker.removed.NA==""))==0){
      # remove empty  data (""," ")
      tmp.marker.removed.empty<-tmp.marker.removed.NA
    }else{
      # remove empty  data (""," ")
      tmp.marker.removed.empty<-tmp.marker.removed.NA[-which(tmp.marker.removed.NA=="")]
    }
    
    # remove depulicate
    if(length(tmp.marker.removed.empty)==0){
      tmp.marker.removed.dup<-tmp.marker.removed.empty
    }else{
      tmp.marker.removed.dup<-tmp.marker.removed.empty[!duplicated(tmp.marker.removed.empty)]
    }
    # remove family
    # check whether has "family" or "class" character
    family.class.index<-grep("family|class",ignore.case = T,tmp.marker.removed.dup)
    if(length(family.class.index)!=0){
      tmp.marker.removed.family.class<-tmp.marker.removed.dup[-family.class.index]
      if(length(tmp.marker.removed.family.class)==0){tmp.marker.removed.family.class<-NULL}
    }else{tmp.marker.removed.family.class<-tmp.marker.removed.dup}
    # final standarded marker
    tmp.standard.marker<-tmp.marker.removed.family.class
    Cell.Marker.list[[i]]<-tmp.standard.marker
    cell.type.name<-c(cell.type.name, CellNameList[i])
  }
  names(Cell.Marker.list)<-cell.type.name
  return(Cell.Marker.list)
}

all.data.cell.marker<-Construct.cellmarker(CellNameList = all.related,AllData = all.data)
all.data.celltype.marker<-Construct.Tissue.marker(TissueNameList = Cell.type,AllData = all.data)

# as.character(unlist(T.cell.family.marker,B.cell.family.marker))
Haitao.table<-read.csv("../Haitao_scRNA_Project/Haitao_Wen_Report_updated/All_biomarkers_from_haitao_seurat_clusters_nFeature_RNA_4000_percentmt_20.csv")
# create unique marker gene in each cluster and store in list 
separate.gene<-function(SeuratMarkerTable=Haitao.table,p.adj.cutoff=0.05){
  SeuratMarkerTable$cluster<-as.character(SeuratMarkerTable$cluster)
  my.cluster.list<-unique(SeuratMarkerTable$cluster)
  my.cluster.gene<-list()
  my.name<-c()
  for (i in 1:length(my.cluster.list)){
  index.cluster<-SeuratMarkerTable$cluster==my.cluster.list[i]
  tmp.table<-SeuratMarkerTable[index.cluster,]
  index.padj<-tmp.table$p_val_adj<=p.adj.cutoff
  my.cluster.gene[[i]]<-tmp.table$gene[index.padj]
  my.name<-c(my.name,my.cluster.list[i])
  }
  names(my.cluster.gene)<-my.name
  return(my.cluster.gene)
}
marker.gene<-separate.gene(SeuratMarkerTable=Haitao.table,p.adj.cutoff=0.05)


## create contingency table
##  |q   |   |k     |q|b|        | mapped gene   |               | selected.gene = my gene list intersects with all marker|
##  |    |   |   =  |c|d|    --> | unmapped gene |               | 
##  |m   |n  |N     |m|n|N       | Cell Marker   |Other markers  | All marker
##
##

cal.enrich.table<-function(MySeurat.marker=marker.gene,CellMarker=all.data.cell.marker,cluster.name=0,adjust.method="bonferroni",p.adj.cutoff=0.05){
  # cluster name localization 
  marker.name.index<-grep(paste0("^",as.character(cluster.name),"$"),as.character(names(MySeurat.marker)))
  my.gene.list<-MySeurat.marker[[marker.name.index]]
  all.marker<-unique(as.character(unlist(CellMarker)))
  N<-length(all.marker)
  selected.gene<-intersect(my.gene.list,all.marker)
  k=length(selected.gene)
  my.table<-c()
  my.rownames<-c()
  my.hitgene.list<-c()
  for (i in 1:length(CellMarker)){
    tmp.cellmarker<-CellMarker[[i]]
    tmp.celltype<-names(CellMarker)[i]
    m=length(tmp.cellmarker)
    q.gene<-intersect(selected.gene,tmp.cellmarker)
    q=length(q.gene)
    n=N-m
    hyper.p<-phyper(q,m,n,k,lower.tail = F) 
    # b=k-q
    # c=m-q
    # d=n-b
    # tmp.matrix<-matrix(c(q,b,c,d),nrow = 2,byrow = T)
    # fisher.res<-fisher.test(tmp.matrix)
    # fisher.p<-fisher.res$p.value
    tmp.table.row<-c(q,m,k,N,as.character(hyper.p))
    my.table<-rbind.data.frame(my.table,tmp.table.row)
    my.rownames<-c(my.rownames,names(CellMarker)[i])
    my.hitgene.list<-c(my.hitgene.list,paste(q.gene,collapse = " ;"))
  }
  p.adj.tmp<- p.adjust(as.numeric(my.table[,5]),
                       method = adjust.method)
  my.table[,6]<-p.adj.tmp
  my.table[,7]<-my.hitgene.list
  rownames(my.table)<-my.rownames
  colnames(my.table)<-c("HitGene","GeneInCellType","SelectedGene","AllMarkers","p.value","p.adj","HitGene_List")
  my.table<-my.table[order(my.table$p.adj),]
  my.table<-my.table[my.table$p.adj<=p.adj.cutoff,]
  my.table<-cbind.data.frame(cluster=rep(cluster.name,nrow(my.table)),my.table)
  HitGeneEqualToZero<-which(my.table$HitGene==0)
  if(length(HitGeneEqualToZero)>0){
    my.table<-my.table[-HitGeneEqualToZero,]
    return(my.table)
  }else{return(my.table)}
}

haitao.table<-c()
celltype.table<-cal.enrich.table(MySeurat.marker=marker.gene,CellMarker=all.data.cell.marker,cluster.name = 0)
for(i in names(marker.gene)){
  print(i)
  tmp.table<-cal.enrich.table(MySeurat.marker=marker.gene,CellMarker=all.data.cell.marker,cluster.name = i)
  haitao.table<-rbind.data.frame(haitao.table,tmp.table)
}
haitao.table2<-cbind.data.frame(Celltype=rownames(haitao.table),haitao.table)
write.table(haitao.table2,file = "Haitao_Project.csv",sep ="," ,quote = F,row.names = F)










