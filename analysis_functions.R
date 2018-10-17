library(readxl)
library(data.table)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(stringi)
library(superheat)
library(RColorBrewer)
library(factoextra)
recode_clusters <- function(cr2){
  #This is complex and annoying
  #Recode clusters
  cr2$cluster[cr2$cluster==3] <- 6
  cr2$cluster[cr2$cluster==4] <- 7
  cr2$cluster[cr2$cluster==1] <- 8
  cr2$cluster[cr2$cluster==5] <- 9
  cr2$cluster[cr2$cluster==2] <- 10
  
  cr2$cluster[cr2$cluster==6] <- 1
  cr2$cluster[cr2$cluster==7] <- 2
  cr2$cluster[cr2$cluster==8] <- 3
  cr2$cluster[cr2$cluster==9] <- 4
  cr2$cluster[cr2$cluster==10] <- 5
  return(cr2)
}

load_clin_data <- function(){
  clin <- read_excel("./data/clinical.xlsx")
  clin <- clin[,c(4,1:3,5:79)]
  clin <- clin[,c(-2:-3,-53:-79)]
  clin <- sapply(clin,tolower)
  clin <- as.data.frame(clin,stringsAsFactors=F)
  clin$ID <- as.integer(clin$ID)
  clin$DOB <- as.integer(clin$DOB)
  clin$`Date dx` <- as.integer(clin$`Date dx`)
  clin$`Age at Dx` <- as.numeric(clin$`Age at Dx`)
  clin$`date of image` <- as.numeric(clin$`date of image`)
  clin <- clin[!(clin$ID %in% c(504,506,507,509)),]
  saveRDS(clin,"./data/clin.rds")
}

load_path_data <- function(){
  path <- read_excel("./data/pathology.xlsx",skip = 9,sheet=2)
  path <- path[,c(-2:-28)]
  names <- colnames(path)
  names[1] <- "ID"
  colnames(path) <- names
  path <- sapply(path,tolower)
  path <- as.data.frame(path,stringsAsFactors=F)
  path$ID <- as.integer(path$ID)
  path<- path[!(path$ID %in% c(504,506,507,509)),]
  saveRDS(path,"./data/path.rds")
}

noDups <- function(resp_sig){
  path <- read_excel("./data/pathology.xlsx",skip = 9,sheet=2)
  path <- path[,c(-2:-28)]
  names <- colnames(path)
  names[1] <- "ID"
  colnames(path) <- names
  pathres <- merge(path,resp_sig,by="ID")
  names_noID <- colnames(pathres)
  names_noID <- names_noID[names_noID!="ID"]
  return(names_noID)
}

initiate_data <- function(resp_sig,names_noID){
  clin <- readRDS("./data/clin.rds")
  path <- readRDS("./data/path.rds")
  clinres <- merge(clin,resp_sig,by="ID")
  pathres <- merge(path,resp_sig,by="ID")

  clinres <- clinres[,!(colnames(clinres) %in% rbind(names_noID))]
  pathres <- pathres[!(pathres$ID %in% c(58)),]
  pathres <- subset(pathres,!duplicated(pathres$ID))
  pathres$`Depth (mm)`[pathres$`Depth (mm)`=="Unk"] <- NA
  pathres$`Depth (mm)` <- as.numeric(pathres$`Depth (mm)`)
  pathres<- pathres[!(pathres$ID %in% c(504,506,507,509)),]
  data <- merge(clinres,pathres,by="ID")
  return(data)
}

excludeUnk <- function(matrix){
  mat_e <- matrix[!rownames(matrix)=="unk",]
  if ((dim(matrix)[1]-dim(mat_e)[1])==1){
    return(mat_e)
  }
  else if ((dim(matrix)[1]-dim(mat_e)[1])==0){
    warning("No rows excluded")
    return(mat_e)
  }
  else {
    stop("More than 1 row excluded")
  }
}

readRecode <- function(filename){
  file <- paste("./data/",filename,".csv",sep="")
  loc2 <- fread(file,header = T)
  loc2 <- data.table(loc2[,2:7])
  loc2_agg <- loc2[,lapply(.SD,sum), by=V2]
  loc2_agg_num <- as.data.frame.matrix(loc2_agg[,2:6])
  rownames(loc2_agg_num) <- loc2_agg$V2
  print("n = ")
  print(sum(colSums(loc2_agg_num)))
  print(colSums(loc2_agg_num))
  print("Numbers: ")
  print(loc2_agg_num)
  print("Percents: ")
  print(round(t(t(loc2_agg_num)/colSums(loc2_agg_num)*100)))
  print(chisq.test(loc2_agg_num))
  percents <- as.data.frame(t(t(loc2_agg_num)/colSums(loc2_agg_num)*100))
  percents$V2 <- loc2_agg$V2
  return(percents)
}

chiSquare <- function(loc2_agg_num){
  percents <- round(t(t(loc2_agg_num)/colSums(loc2_agg_num))*100)
  chisq <- chisq.test(loc2_agg_num)
  ans <- list(percents,chisq)
  return(ans)
}

prettyprint <- function(loc2_agg_num){
  print("n = ")
  print(sum(colSums(loc2_agg_num)))
  print(colSums(loc2_agg_num))
  print("Numbers: ")
  print(loc2_agg_num)
  print("Percents: ")
  print(round(t(t(loc2_agg_num)/colSums(loc2_agg_num)*100)))
  chisq.test(loc2_agg_num)
}


# Plotting functions
plotclust <- function(clust_e,num){
  clustername <- paste("Cluster",num,sep="")
  ggpt <- ggbarplot(clust_e,"name",clustername,
                    ylim=c(0,100),
                    title=clustername,
                    color="name",
                    fill="name",
                    palette=c("black",#pigmented melanoma
                              "darkblue",#dysplastic nevus
                              "chocolate4",#SK
                              "coral3",#BCC
                              "burlywood1",#amelanotic melanoma
                              "gold3",#other
                              "cyan4",#dermal nevus
                              "coral",#AK
                              "chocolate",#SCC
                              "cyan1",#compound nevus
                              "darkviolet",#LMalign
                              "plum1"),#SLentigo
                    ylab="Percentage of Responses",
                    xlab="",
                    x.text.angle=45) +
    theme(legend.position="none") +
    theme(text=element_text(face="bold",size=14))
  return(ggpt)
}

plotclust_facet <- function(clust_melt){
  palette=c("black",#pigmented melanoma
            "darkblue",#dysplastic nevus
            "chocolate4",#SK
            "plum1",#SLentigo
            "darkviolet",#LMalign
            "cyan1",#compound nevus
            "coral3",#BCC
            "chocolate",#SCC
            "coral",#AK
            "burlywood1",#amelanotic melanoma
            "cyan4",#dermal nevus
            "gold3")#other
            
            
            
            
            
  p <- ggbarplot(clust_melt,"name","value",
                    ylim=c(0,100),
                    color="name",
                    fill="name",
                    palette=palette,
                    ylab="Assignment Percentage",
                    xlab="",
                 label = T,
                 lab.size = 7,
                 lab.nb.digits=0,
                    x.text.angle=30) +
    theme(legend.position="none") +
    theme(text=element_text(face="bold",size=20))+
    theme(plot.margin = margin(0.7,1,0.7,3,"cm"))
  p <- facet(p, facet.by = "variable",ncol=1,
             panel.labs = list(
               variable = c("Typical (CL1)",
                            "Nevus-like (CL2)",
                            "Amelanotic/NMSC-like (CL3)",
                            "Seborrheic Keratosis (SK)-like (CL4)",
                            "Lentigo/LM-like (CL5)")
             ),
             panel.labs.background = list(color = "steelblue", fill = "steelblue", size = 0.5),
             panel.labs.font = list(color = "white"),
             panel.labs.font.x = list(angle = 0, color = "white",size=20,face="bold",hjust=0.05))
  return(p)
}

plotheatmap <- function(x){
  p <- superheat(x,
            membership.cols = resp_sig_incl$cluster,
            order.cols = order(resp_sig_incl$ID),
            
            #n.clusters.rows = 5,
            #clustering.method = "hierarchical",
            #dist.method = "manhattan",
            #pretty.order.rows = T,
            #row.dendrogram = T,
            
            order.rows = as.integer(c(12,11,4,3,2,9,7,6,5,8,10,1)),
            
            #Column Title
            column.title = "Percent of Assignments",
            column.title.size = 10,
            
            #Label sizes
            left.label.size = 0.2,
            left.label.text.size = 10,
            bottom.label.size = 0.1,
            bottom.label.text.size = 10,
            
            #Label colors
            left.label.col = "white",
            left.label.text.alignment = "left",
            bottom.label.col = c("coral1","cyan1","chocolate","gray72","plum1"),
            
            #Legend
            legend.height = 0.2,
            legend.width = 5,
            legend.text.size = 30,
            #scale=T,
            #legend.breaks = c(0,25,50,75,100),
            
            #Color scheme
            heat.pal = c("white","purple4"),
            heat.na.col = "white",
            heat.pal.values = c(0,1),
            
            #Grids
            grid.hline = F,
            grid.vline = T,
            grid.vline.col = "black",
            grid.vline.size = 0.7)
  return(p)
}

plotheatmap_path <- function(x){
  p <- superheat(x,
                 membership.cols = resp_sig_incl$cluster,
                 order.cols = order(resp_sig_incl$ID),
                 
                 #Column Title
                 column.title = "Clusters",
                 column.title.size = 10,
                 
                 #Label sizes
                 left.label.size = 0.2,
                 left.label.text.size = 10,
                 bottom.label.size = 0.1,
                 bottom.label.text.size = 10,
                 
                 #Label colors
                 left.label.col = "white",
                 left.label.text.alignment = "right",
                 bottom.label.col = c("cyan1","chocolate","plum1","coral1","gray72"),
                 
                 #Legend
                 legend = F,
                 legend.height = 0.2,
                 legend.width = 5,
                 legend.text.size = 20,
                 
                 #Color scheme
                 heat.pal = c("white","black"),
                 heat.na.col = "white",
                 #heat.pal.values = c(0,1),
                 
                 #Grids
                 grid.hline = F,
                 grid.vline = T,
                 grid.vline.col = "black",
                 grid.vline.size = 0.7)
  return(p)
}

plotdens <- function(cr3,column,xlab,setlim=T,lim=1){
  palette="jco"
  if (setlim){
    p <- ggdensity(cr3,x=column,
                   add = "mean",rug=T,
                   fill="cluster",
                   add_density = T,
                   xlim=c(0,lim),
                   xlab=xlab,
                   title=paste(column," by Cluster",sep=""),
                   palette = palette) 
  } else {
    p <- ggdensity(cr3,x=column,
                   add = "mean",rug=T,
                   fill="cluster",
                   add_density = T,
                   xlab=xlab,
                   title=paste(column," by Cluster",sep=""),
                   palette = palette)
  }
  p <- p + theme(axis.text.x  = element_text(face="bold",size=15),
                 axis.text.y  = element_text(face="bold",size=15),
                 axis.title.x = element_text(face="bold",size=20),
                 axis.title.y = element_text(face="bold",size=20))+
    theme(legend.text = element_text(face="bold",size=15),
          legend.title = element_text(face="bold",size=15))+
    theme(plot.title = element_text(hjust=0.5,face="bold",size=20))
  
  p <- facet(p, facet.by = "cluster",ncol=1,
        panel.labs = list(
          cluster = c("Typical-Like (CL1)",
                      "DN-like (CL2)",
                      "Erythematous/Keratotic (CL3)",
                      "SK-like (CL4)",
                      "Lentigo-like(CL5)")
        ),
        panel.labs.background = list(color = "steelblue", fill = "steelblue", size = 0.5),
        panel.labs.font = list(color = "white"),
        panel.labs.font.x = list(angle = 0, color = "white",size=15,face="bold",hjust=0.1))
  return(p)
}


plotclust_facet_supplementary <- function(clust_melt){
  palette=c("black",#pigmented melanoma
            "darkblue",#dysplastic nevus
            "chocolate4",#SK
            "plum1",#SLentigo
            "darkviolet",#LMalign
            "cyan1",#compound nevus
            "coral3",#BCC
            "chocolate",#SCC
            "coral",#AK
            "burlywood1",#amelanotic melanoma
            "cyan4",#dermal nevus
            "gold3")#other
  
  p <- ggbarplot(clust_melt,"name","value",
                 ylim=c(0,100),
                 color="name",
                 fill="name",
                 palette=palette,
                 ylab="Assignment Percentage",
                 xlab="",
                 label = T,
                 lab.size = 7,
                 lab.nb.digits=0,
                 x.text.angle=30) +
    theme(legend.position="none") +
    theme(text=element_text(face="bold",size=20))+
    theme(plot.margin = margin(0.7,1,0.7,3,"cm"))
  p <- facet(p, facet.by = "variable",ncol=1,
             # panel.labs = list(
             #   variable = c("Typical (CL1)","Dysplastic Nevus (DN)-like (CL2)","Amelanotic/Erythematous/Keratotic (CL3)","Seborrheic Keratosis (SK)-like (CL4)","Lentigo-like (CL5)")
             # ),
             panel.labs.background = list(color = "steelblue", fill = "steelblue", size = 0.5),
             panel.labs.font = list(color = "white"),
             panel.labs.font.x = list(angle = 0, color = "white",size=20,face="bold",hjust=0.05))
  return(p)
}
