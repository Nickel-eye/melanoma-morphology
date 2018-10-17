source("analysis_functions.R")

### Determine if there is someone who skipped all of them

### Construct responses data ----
num_photos <- 413
responses <- data.frame(ID=integer(num_photos),
                 BCC=integer(num_photos),
                 ActK=integer(num_photos),
                 SCC=integer(num_photos),
                 Wart=integer(num_photos),
                 SebK=integer(num_photos),
                 CompN=integer(num_photos),
                 DermN=integer(num_photos),
                 DyspN=integer(num_photos),
                 MelPig=integer(num_photos),
                 MelAml=integer(num_photos),
                 SolLen=integer(num_photos),
                 LenMlgn=integer(num_photos),
                 Ecz=integer(num_photos),
                 Other=integer(num_photos),
                 Other_1=character(num_photos),
                 Other_2=character(num_photos),
                 Other_3=character(num_photos),
                 stringsAsFactors=FALSE)

for (i in seq(1,num_photos,1)){
  table <- read_xls('./data/raw_responses.xls',sheet=i)
  
  #Get the photo ID
  IDtext <- table$`Clinical appearance of lesions study`[1]
  IDtext_split <- strsplit(IDtext," ")
  responses$ID[i] <- as.integer(IDtext_split[[1]][2])
  
  #Get the count of responses
  responses[i,2:15] <- as.integer(table$X__2[3:16])

  #If there are "other" responses, put them in Other 1, 2, 3
  #The maximum number of "other" responses was 3, so there are only 3 slots
  if (responses$Other[i] != 0){
    responses$Other_1[i] <- table$X__2[38]
    responses$Other_2[i] <- table$X__2[39]
    responses$Other_3[i] <- table$X__2[40]
  }
}

rm(list=setdiff(ls(), "responses"))
write.xlsx(responses,"./data/responses_original.xlsx",row.names = FALSE)

#responses recoded (original)
responses <- read_xlsx("./data/responses_col_rem.xlsx",1)

#responses after excluding responses from Nicole + other person with no email
responses <- read_xlsx("./data/responses_exclude_2.xlsx",1)
saveRDS(responses,"./data/responses.rds")

### Cluster Analysis ----

responses <- readRDS("./data/responses.rds")
#Re-order based on how data should meaningfully look
responses <- responses[,c(1,10,9,8,7,11,2,4,3,6,12,13,5,14:24)]

#what degree of fuzzy to exclude
responses <- responses[responses$Exclude_Fuzzy<3,]
responses <- responses[,!(colnames(responses) %in% c("Exclude_Fuzzy"))]
responses <- responses[!(responses$ID %in% c(58,504,506,507,509)),]
resp_noID <- responses[,!(colnames(responses) %in% c("ID"))]

colnames(resp_noID) <- c("Pigmented Melanoma","Dysplastic Nevus","Dermal Nevus","Compound Nevus",
                         "Amelanotic Melanoma","BCC","SCC","AK","SK",
                         "Solar Lentigo","Lentigo Maligna","Wart","Eczema","LK",
                         "JN","PN","HemTrauHemg","AngSerKer",
                         "Poroma","Dermatofibroma","CMast","MCC")

#Determine the columns contributing most significantly to responses
sums <- colSums(resp_noID)
sig_sums <- sums[sums>quantile(sums,0.5)]
nonsig_sums <- sums[!(names(sums)%in%names(sig_sums))]
resp_sig <- resp_noID[,colnames(resp_noID)%in%names(sig_sums)]
resp_nonsig <- resp_noID[,colnames(resp_noID)%in%names(nonsig_sums)]
resp_sig$Other <- NA
resp_sig$Other <- rowSums(resp_nonsig)

#resp_sig <- resp_sig_byr
resp_sig <- resp_sig/rowSums(resp_sig) #normalize

#Clustering with k-means
centroids = 5
km <- kmeans(resp_sig,centroids,nstart=100)

#Re-number the cluster numbers to be as in the paper
clust <- as.data.frame(t(km$centers))*100
original_order <- unname(c(which.max(clust[1,]),
                    which.max(clust[2,]),
                    which.max(clust[5,]),
                    which.max(clust[9,]),
                    which.max(clust[10,])))

km$newcluster <- km$cluster
for (i in seq(1,length(km$newcluster))){
  km$newcluster[i] <- which(original_order==km$newcluster[i])
}

resp_sig$ID <- responses$ID
resp_sig$cluster <- km$newcluster

clust <- clust[,original_order]
colnames(clust) <- c("1","2","3","4","5")
clust_e <- clust[!(rownames(clust) %in% c("ID","cluster")),]
clust_e$name <- rownames(clust_e)
#colnames(clust_e) <- c("3","5","1","2","4","name")

#Merge data into clinpath shared data.frame
load_clin_data()
load_path_data()
names_noID <- noDups(resp_sig)
clinpath <- initiate_data(resp_sig,names_noID)

#Heatmap of response data
resp_sig_incl <- resp_sig[(resp_sig$ID %in% clinpath$ID),]
resp_sig_incl$cluster[resp_sig_incl$cluster==1] <- paste("1"," (n = 136)",sep="")
resp_sig_incl$cluster[resp_sig_incl$cluster==2] <- paste("2"," (n = 81)",sep="")
resp_sig_incl$cluster[resp_sig_incl$cluster==3] <- paste("3"," (n = 70)",sep="")
resp_sig_incl$cluster[resp_sig_incl$cluster==4] <- paste("4"," (n = 68)",sep="")
resp_sig_incl$cluster[resp_sig_incl$cluster==5] <- paste("5"," (n = 45)",sep="")
x <- t(as.matrix(resp_sig_incl[,c(1:12)]))*100
x <- x[nrow(x):1,]
colnames(x) <- resp_sig_incl$ID
tiff("./figs/Figure2_clust.tiff",height=2592,width=8000,res=300)
plotheatmap(x)
dev.off()

#Dendrogram of rows
dd <- dist(scale(x),method="manhattan")
hc <- hclust(dd,method="ward.D2")
dg <- as.dendrogram(hc)
dg <- reorder(dg,12:1)

tiff("./figs/Figure2_dendrogram.tiff",height=1392,width=1000,res=400)
fviz_dend(dg,
          
          #size
          cex=1,
          
          #k
          #k=5,
          k_colors="jco",
          color_labels_by_k = T,
          
          #rectangle
          rect=T,
          rect_fill=T,
          rect_border="jco",
          
          #horizontal
          horiz=T)
dev.off()

#Barplot of cluster profiles
clust_e <- clust_e[c(1,2,9,10,11,4,6,7,8,5,3,12),]
clust_e$name[1] <- "Pigmented MM"
clust_e$name[10] <- "Amelanotic MM"
clust_melt <- melt(clust_e)
clust_melt$variable <- as.numeric(as.character(clust_melt$variable))
clust_melt <- clust_melt[order(clust_melt$variable),]
tiff("./figs/Figure3_facet.tiff",height=4760,width=4000,res=400)
plotclust_facet(clust_melt)
dev.off()

#Heatmap of Pathology Stains
stains <- clinpath[,c(70:87,118,1)]
#remove H&E column
stains <- stains[,c(-10)]
stains$cluster[stains$cluster==1] <- paste("1"," (n = 136)",sep="")
stains$cluster[stains$cluster==2] <- paste("2"," (n = 81)",sep="")
stains$cluster[stains$cluster==3] <- paste("3"," (n = 72)",sep="")
stains$cluster[stains$cluster==4] <- paste("4"," (n = 67)",sep="")
stains$cluster[stains$cluster==5] <- paste("5"," (n = 44)",sep="")

stains[stains=="yes"] <- 1
stains[stains=="no"] <- 0
stains[stains=="unk"] <- NA
names <- stains$ID
stains <- sapply(stains,as.numeric)
x <- t(as.matrix(stains[,c(1:17)]))
colnames(x) <- names
x <- x[nrow(x):1,]
plotheatmap_path(x)


#Make cluster names
clinpath$Name <- as.character(clinpath$cluster)
clinpath$Name[clinpath$cluster==1] <- "Typical"
clinpath$Name[clinpath$cluster==2] <- "DN-like"
clinpath$Name[clinpath$cluster==3] <- "Amelanotic/NMSC-like"
clinpath$Name[clinpath$cluster==4] <- "SK-like"
clinpath$Name[clinpath$cluster==5] <- "Lentigo/LM-like"

#Gender

genders <- as.data.frame.matrix(table(clinpath$Gender,clinpath$cluster))
genders_e <- excludeUnk(genders)
genders_e
chiSquare(genders_e)

#Age
#The age for subject 322 is wrong: should be 34, but was misrecorded as 6
clinpath$`Age at Dx`[clinpath$ID==322] <- 34 #from chart

ages <- list()
ages_mean <- c(rep(0,5))
ages_sds <- c(rep(0,5))
for (i in seq(1,5)){
  ages[[i]] <- clinpath$`Age at Dx`[clinpath$cluster==i]
  ages_mean[i] <- mean(clinpath$`Age at Dx`[clinpath$cluster==i])
  ages_sds[i] <- sd(clinpath$`Age at Dx`[clinpath$cluster==i])
}

kruskal.test(`Age at Dx`~cluster,clinpath)
pairwise.wilcox.test(clinpath$`Age at Dx`,clinpath$cluster,p.adjust.method = "bonferroni")

cr2 <- clinpath
names(cr2)[8] <- "Age"
comparisons <- list(c("2","1"),
                    c("2","3"),
                    c("2","4"),
                    c("2","5"))

compare_means(Age~cluster,cr2,method="wilcox.test",paired=F,
              p.adjust.method = "bonferroni")

fontsize <- 10
ageplot <- ggboxplot(cr2,x="cluster",y="Age",
         # color="gray80",
         xlab="Name",
         fill = "white",
         add = c("jitter"),
         add.params = list(fill = "white"))+
  #theme_bw()+
  theme(axis.text.x  = element_text(face="bold",size=22),
        axis.text.y  = element_text(face="bold",size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=25))+
  scale_y_continuous(breaks=c(20,30,40,50,60,70,80,90,100))+
  scale_x_discrete(labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  stat_compare_means(aes(label=format.pval(..p.adj..,digits=1)),comparisons=comparisons)
  #stat_compare_means(label.x=2.5,label.y=150,size=10)
#ageplot$layers[[4]]$aes_params$textsize <- fontsize
tiff("./figs/Figure4_age.tiff",height=4000,width=4000,res=400)
ageplot
dev.off()


#Age, Depth, Mitoses histograms
cr3 <- clinpath
names(cr3)[8] <- "Age"
names(cr3)[59] <- "Depth"
names(cr3)[90] <- "Area"
cr3$cluster[cr3$cluster==1] <- "CL1"
cr3$cluster[cr3$cluster==2] <- "CL2"
cr3$cluster[cr3$cluster==3] <- "CL3"
cr3$cluster[cr3$cluster==4] <- "CL4"
cr3$cluster[cr3$cluster==5] <- "CL5"
#Recode to numeric
cr3 <- clinpath
cr3$Mitoses[(cr3$Mitoses %in% c("not identified"))] <- 0
cr3$Mitoses[(cr3$Mitoses %in% c("unk"))] <- NA
cr3$Mitoses <- as.numeric(cr3$Mitoses)

kruskal.test(Mitoses~cluster,cr3)
pairwise.wilcox.test(cr3$Mitoses,cr3$cluster,p.adjust.method = "bonferroni")

splitname <- strsplit(cr3$Area,"x")
multiply <- function(list){
  return(as.numeric(list[1])*as.numeric(list[2]))
}
splitname <- unlist(sapply(splitname,multiply,simplify=F))
cr3$Area <- splitname
plotdens(cr3,"Age","Age (years)",setlim=F)
plotdens(cr3,"Depth","Depth (mm)",setlim=F)
plotdens(cr3,"Mitoses","Mitoses (per mm2)",setlim=T,lim=5)
plotdens(cr3,"Area","Area (mm2)",setlim=T,lim=200)


#Lesion found by
foundby <- as.data.frame.matrix(table(clinpath$`Lesion found by`,clinpath$cluster))
write.csv(foundby,"./data/foundby.csv")
foundby2 <- fread("./data/foundby2.csv",header = T)
foundby2 <- data.table(foundby2[,2:7])
foundby2_agg <- foundby2[,lapply(.SD,sum),by=V2]
foundby2_agg_num <- as.data.frame.matrix(foundby2_agg[,2:6])
rownames(foundby2_agg_num) <- foundby2_agg$V2
foundby2_agg_num
percents <- t(t(foundby2_agg_num)/colSums(foundby2_agg_num))*100
chisq.test(foundby2_agg_num)

for (i in seq(1,6)){
  testtbl <- foundby2_agg_num[5,]
  testtbl[2,] <- colSums(foundby2_agg_num)
  testtbl[2,] <- testtbl[2,]-testtbl[1,]
  print(chisq.test(testtbl))
}

foundby2_agg_tab <- melt(percents,id=c("V2"))
colnames(foundby2_agg_tab) <- c("Found_By","Cluster","Count")
foundby2_agg_tab$Found_By <- factor(foundby2_agg_tab$Found_By,
                                    levels=c("Patient or Family Member","Dermatologist","Non-Dermatologist Professional"))

ggbar <- ggbarplot(foundby2_agg_tab,"Cluster","Count",
          fill = "Found_By",
          #label=F,
          ylab = "Percent",
          ylim=c(0,80),
          position=position_dodge(0.8),
          palette="npg")+
  theme(axis.text.x  = element_text(face="bold",size=22),
        axis.text.y  = element_text(face="bold",size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=25),
        legend.text  = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25))+
  scale_y_continuous(breaks=seq(0,80,10))+
  scale_x_continuous(breaks=seq(1,5),
                     labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  guides(fill=guide_legend(
    keywidth = 2,
    keyheight = 2))

tiff("./figs/Figure4_foundby.tiff",height=4000,width=4000,res=400)
margin(t=2,r=2,b=2,l=2,unit="pt")
ggpar(ggbar,
      legend.title = "Found By",
      legend="off")
dev.off()

#Prior NMSC
nmsc <- as.data.frame.matrix(table(clinpath$`Prior NMSC`,clinpath$cluster))
nmsc_e <- nmsc[c(1,3),]
sum(colSums(nmsc_e))
colSums(nmsc_e)
nmsc_e
t(t(nmsc_e)/colSums(nmsc_e))*100
chisq.test(nmsc_e)

#Prior melanoma
mel <- as.data.frame.matrix(table(clinpath$`Prior Melanoma`,clinpath$cluster))
write.csv(mel,"./data/mel.csv")
mel2 <- fread("./data/mel2.csv",header=T)
mel2 <- data.table(mel2[,2:7])
mel2_agg <- mel2[,lapply(.SD,sum),by=V2]
mel2_agg_num <- as.data.frame.matrix(mel2_agg[,2:6])
rownames(mel2_agg_num) <- mel2_agg$V2
mel2_e <- mel2_agg_num[1:2,]
mel2_e
t(t(mel2_e)/colSums(mel2_e))*100
chisq.test(mel2_e)

#Family history
fhx <- as.data.frame.matrix(table(clinpath$`family hx`,clinpath$cluster))
fhx_e <- fhx[c(1,3),]
sum(colSums(fhx_e))
colSums(fhx_e)
fhx_e
t(t(fhx_e)/colSums(fhx_e))*100
chisq.test(fhx_e)

#Sunburns/sun exposure
sunexp <- as.data.frame.matrix(table(clinpath$`Sunburns/sun exposure`,clinpath$cluster))
sunexp_e <- sunexp[c(1,3),]
sum(colSums(sunexp_e))
colSums(sunexp_e)
sunexp_e
t(t(sunexp_e)/colSums(sunexp_e))*100
chisq.test(sunexp_e)

#PUVA
puva <- as.data.frame.matrix(table(clinpath$PUVA,clinpath$cluster))
t(t(puva)/colSums(puva))*100
chisq.test(puva)

#Tanning bed
tanbed <- as.data.frame.matrix(table(clinpath$`Tanning bed`,clinpath$cluster))
tanbed_e <- tanbed[c(1,3),]
sum(colSums(tanbed_e))
colSums(tanbed_e)
tanbed_e
t(t(tanbed_e)/colSums(tanbed_e))*100
chisq.test(tanbed_e)

#Skin type
skt <- as.data.frame.matrix(table(clinpath$`skin type`,clinpath$cluster))
skt_e <- skt[!rownames(skt)=="unk",]
write.csv(skt_e,"./data/sktype.csv")
readRecode("sktype2")

#Density of nevi
nden <- as.data.frame.matrix(table(clinpath$`density of nevi`,clinpath$cluster))
nden_e <- nden[!rownames(nden)=="unk",]
write.csv(nden_e,"./data/nevidensity.csv")
file <- "nevidensity2"
readRecode(file)

#Dermatoheliosis
dhel <- as.data.frame.matrix(table(clinpath$dermatoheliosis,clinpath$cluster))
dhel_e <- nden[!rownames(nden)=="unk",]
write.csv(dhel_e,"./data/dermatoheliosis.csv")
file <- "dermatoheliosis2"
readRecode(file)

### Lesion features (Table 2)----

#Left/right lesion location
lrloc <- as.data.frame.matrix(table(clinpath$`Left/Right lesion location`,clinpath$cluster))
lrloc_e <- lrloc[c(1:3),]
prettyprint(lrloc_e)

#Lesion location
loc <- as.data.frame.matrix(table(clinpath$`Lesion location`,clinpath$cluster))
write.csv(loc,"./data/location.csv") #write CSV and do this by hand
readRecode("location_2")

loc2_agg <- readRecode("location_2")
loc2_agg <- loc2_agg[c(4,2,3,1),]
loc2_agg_tab <- melt(loc2_agg,id=c("V2"))
colnames(loc2_agg_tab) <- c("Location","Cluster","Count")
loc2_agg_tab$Location <- factor(loc2_agg_tab$Location,
                                levels=c("head/neck","trunk","upper limb","lower limb"))

ggbar <- ggbarplot(loc2_agg_tab,"Cluster","Count",
          fill = "Location",
          label=F,
          ylim=c(0,80),
          position=position_dodge(0.8),
          palette="npg",
          ylab="Percent")+
  theme(axis.text.x  = element_text(face="bold",size=22),
        axis.text.y  = element_text(face="bold",size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=25),
        legend.text  = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25))+
  scale_y_continuous(breaks=seq(0,70,10))+
  scale_x_discrete(labels = c("Typical",
                                "Nevus\n-like",
                                "Amelanotic/\nNMSC-like",
                                "SK\n-like",
                                "Lentigo/\nLM-like"))+
  guides(fill=guide_legend(
    keywidth = 2,
    keyheight = 2))

tiff("./figs/Figure4_lesionloc.tiff",height=4000,width=4000,res=400)
margin(t=2,r=2,b=2,l=2,unit="pt")
ggpar(ggbar,
      legend.title = "Lesion Location",
      legend="off")
dev.off()

#Symptoms
sx <- function(table){
  table <- table[!rownames(table)=="unk",]
  print(table)
  print(sum(colSums(table)))
  print(colSums(table))
  print(table[rownames(table)=="yes",])
  perc <- round(t(t(table)/colSums(table))*100)
  print(perc[rownames(perc)=="yes",])
  chisq.test(table)
}

#do this over by changing the cr$____ in the code below
table <- as.data.frame.matrix(table(clinpath$`New lesion`,clinpath$cluster))
sx(table)

#Recurrence in same location - not used
recsm <- as.data.frame.matrix(table(cr$`Recurrance in same location`,cr$cluster))

#Irregular color
ircol <- as.data.frame.matrix(table(clinpath$`irregular color`,clinpath$cluster))
ircol_e <- ircol[c(1,3),]
sum(colSums(ircol_e))
colSums(ircol_e)
ircol_e
round(t(t(ircol_e)/colSums(ircol_e))*100)
chisq.test(ircol_e)

#Irregular borders
irbor <- as.data.frame.matrix(table(cr$`irregular borders`,cr$cluster))
irbor_e <- irbor[c(1,3),]
sum(colSums(irbor_e))
colSums(irbor_e)
irbor_e
round(t(t(irbor_e)/colSums(irbor_e))*100)
chisq.test(irbor_e)

rmunk <- function(df){
  df_e <- df[!rownames(df)=="unk",]
  print(sum(colSums(df_e)))
  percents <- t(t(df_e)/colSums(df_e))*100
  for (i in seq(1,5)){
    print(paste(df_e[2,i],"/",colSums(df_e)[i]," ",round(percents[2,i]),"%"))
  }
  chisq.test(df_e)
}

#Colors
cr <- clinpath
brown <- as.data.frame.matrix(table(cr$Brown,cr$cluster))
rmunk(brown)
black <- as.data.frame.matrix(table(cr$Black,cr$cluster))
rmunk(black)
blue <- as.data.frame.matrix(table(cr$Blue,cr$cluster))
rmunk(blue)
pink <- as.data.frame.matrix(table(cr$`Pink/erythema`,cr$cluster))
rmunk(pink)
tan <- as.data.frame.matrix(table(cr$`Tan/Skin colored`,cr$cluster))
rmunk(tan)
pigm <- as.data.frame.matrix(table(cr$Pigmented,cr$cluster))
rmunk(pigm)


#Lesion type
ltype <- as.data.frame.matrix(table(clinpath$`lesion type`,clinpath$cluster))
write.csv(ltype,"./data/ltype.csv")
readRecode("ltype2")


### Pathology (Table 3)----

#Depth
names(clinpath)[59] <- "Depth"
depths <- list()
depths_mean <- c(rep(0,5))
depths_sds <- c(rep(0,5))
for (i in seq(1,5)){
  depths[[i]] <- clinpath$Depth[clinpath$cluster==i]
  depths_mean[i] <- mean(clinpath$Depth[clinpath$cluster==i],na.rm = T)
  depths_sds[i] <- sd(clinpath$Depth[clinpath$cluster==i],na.rm=T)
}
for (i in seq(1,5)){
  print(length(na.omit(depths[[i]])))
}
ktdepth <- kruskal.test(Depth~cluster,clinpath)
pairwise.wilcox.test(clinpath$Depth,clinpath$cluster,p.adjust.method = "bonferroni")
comparisons <- list(c("2","1"),
                    c("2","3"),
                    c("5","1"),
                    c("5","3"))
fontsize <- 10
depthplot <- ggviolin(clinpath,x="cluster",y="Depth",
                       #color="cluster",palette="jama",
                      add = c("jitter"),
                      #add.params = list(shape=17),
                       xlab="Cluster",
                       ylab="Thickness (mm)")+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme(axis.text.x  = element_text(face="bold",size=20),
        axis.text.y  = element_text(face="bold",size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=20))+
  scale_x_discrete(labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  stat_compare_means(label="p.adj",comparisons=comparisons)
  #stat_compare_means(label.x=2.5,label.y=14,size=fontsize
#depthplot$layers[[3]]$aes_params$textsize <- fontsize
tiff("./figs/Figure5_depth.tiff",height=4000,width=4000,res=400)
depthplot
dev.off()


#Depth Histograms

cr3 <- clinpath
names(cr3)[59] <- "Depth"
cr3$cluster[cr3$cluster==1] <- "CL1"
cr3$cluster[cr3$cluster==2] <- "CL2"
cr3$cluster[cr3$cluster==3] <- "CL3"
cr3$cluster[cr3$cluster==4] <- "CL4"
cr3$cluster[cr3$cluster==5] <- "CL5"
#cr3 <- cr3[cr3$cluster %in% c("CL1","CL2"),]
p <- ggdensity(cr3,x="Depth",
               add = "mean",rug=T,
               fill="cluster",
               add_density = T,
               palette = "jco",
               xlab="Depth (mm)",
               title="Depth of Clusters",
               xlim=c(0,5))

p <- p + theme(axis.text.x  = element_text(face="bold",size=15),
               axis.text.y  = element_text(face="bold",size=15),
               axis.title.x = element_text(face="bold",size=20),
               axis.title.y = element_text(face="bold",size=20))+
  theme(legend.text = element_text(face="bold",size=15),
        legend.title = element_text(face="bold",size=15))+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))

facet(p, facet.by = "cluster",ncol=1,
      panel.labs = list(
        cluster = c("Typical-Like (CL1)","DN-like (CL2)","Erythematous/Keratotic (CL3)","SK-like (CL4)","Lentigo-like(CL5)")
      ),
      panel.labs.background = list(color = "steelblue", fill = "steelblue", size = 0.5),
      panel.labs.font = list(color = "white"),
      panel.labs.font.x = list(angle = 0, color = "white",size=15,face="bold",hjust=0.1))

#Mitoses Histograms

#Final Pathological Diagnosis
pathdx <- as.data.frame.matrix(table(clinpath$`Final Path Dx`,clinpath$cluster))
write.csv(pathdx,"./data/pathdx.csv")
readRecode("pathdx2")

path_agg <- readRecode("pathdx2")
path_agg <- path_agg[c(5,3,4,2,6,1),]
path_agg_tab <- melt(path_agg,id=c("V2"))
colnames(path_agg_tab) <- c("Diagnosis","Cluster","Count")
path_agg_tab$Diagnosis <- factor(path_agg_tab$Diagnosis)
plevels <- levels(path_agg_tab$Diagnosis)
path_agg_tab$Diagnosis <- factor(path_agg_tab$Diagnosis,
                            levels=plevels[c(5,2,3,1,6,4)])

ggbar <- ggbarplot(path_agg_tab,"Cluster","Count",
          fill = "Diagnosis",
          label=F,
          ylim=c(0,80),
          position=position_dodge(0.8),
          palette="npg",
          ylab="Percent")+
  theme(axis.text.x  = element_text(face="bold",size=22),
        axis.text.y  = element_text(face="bold",size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=25),
        legend.text  = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25))+
  scale_y_continuous(breaks=seq(0,80,10))+
  scale_x_discrete(labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  guides(fill=guide_legend(
    keywidth = 2,
    keyheight = 2))

tiff("./figs/Figure5_pathdx.tiff",height=4308,width=6000,res=400)
margin(t=2,r=2,b=2,l=2,unit="pt")
ggpar(ggbar,
      legend.title = "Pathologic Diagnosis",
      legend="off")
dev.off()

#Anatomic Level
anatomiclevels <- as.data.frame.matrix(table(clinpath$`Anatomic Level`,clinpath$cluster))
anatomiclevels_e <- excludeUnk(anatomiclevels)
write.csv(anatomiclevels_e,"./data/anlevels.csv")
readRecode("anlevels2")

an_agg <- readRecode("anlevels2")
#path_agg <- path_agg[c(5,3,4,2,6,1),]
an_agg_tab <- melt(an_agg,id=c("V2"))
colnames(an_agg_tab) <- c("Level","Cluster","Count")
an_agg_tab$Level <- factor(an_agg_tab$Level)
#plevels <- levels(path_agg_tab$Diagnosis)
# path_agg_tab$Diagnosis <- factor(path_agg_tab$Diagnosis,
#                                  levels=plevels[c(5,2,3,1,6,4)])

ggbar <- ggbarplot(an_agg_tab,"Cluster","Count",
          fill = "Level",
          label=F,
          ylim=c(0,80),
          position=position_dodge(0.8),
          palette="npg",
          ylab="Percent")+
  theme(axis.text.x  = element_text(face="bold",size=22),
        axis.text.y  = element_text(face="bold",size=25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=25),
        legend.text  = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25))+
  scale_y_continuous(breaks=seq(0,60,10))+
  scale_x_discrete(labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  guides(fill=guide_legend(
    keywidth = 2,
    keyheight = 2))

tiff("./figs/Figure5_anlevel.tiff",height=4000,width=4000,res=400)
margin(t=2,r=2,b=2,l=2,unit="pt")
ggpar(ggbar,
      legend.title = "Anatomic Level",
      legend="off")
dev.off()

#Recode to numeric
mits <- as.data.frame.matrix(table(clinpath$Mitoses,clinpath$cluster))

clinpath$Mitoses[(clinpath$Mitoses %in% c("not identified"))] <- 0
clinpath$Mitoses[(clinpath$Mitoses %in% c("unk"))] <- NA
clinpath$Mitoses <- as.numeric(clinpath$Mitoses)

mits <- list()
mits_mean <- c(rep(0,5))
mits_sds <- c(rep(0,5))
for (i in seq(1,5)){
  mits[[i]] <- clinpath$Mitoses[clinpath$cluster==i]
  mits_mean[i] <- mean(clinpath$Mitoses[clinpath$cluster==i],na.rm = T)
  mits_sds[i] <- sd(clinpath$Mitoses[clinpath$cluster==i],na.rm=T)
}
for (i in seq(1,5)){
  print(length(na.omit(mits[[i]])))
}

kruskal.test(Mitoses~cluster,clinpath)
pairwise.wilcox.test(clinpath$Mitoses,clinpath$cluster,p.adjust.method = "bonferroni")

comparisons <- list(c("1","2"),
                    c("1","3"),
                    c("1","5"),
                    c("2","3"),
                    c("3","4"),
                    c("3","5"))
fontsize <- 8
mitsplot <- ggviolin(clinpath,x="cluster",y="Mitoses",
                       #color="cluster",palette="jama",
                      fill="white",
                       xlab="Cluster",
                       ylab="Mitoses/mm2",
                     add=c("jitter"))+
                     #add.params = list(shape=17))+
                     #add.params = list(fill = "white"))+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30))+
  scale_x_discrete(labels = c("Typical",
                              "Nevus\n-like",
                              "Amelanotic/\nNMSC-like",
                              "SK\n-like",
                              "Lentigo/\nLM-like"))+
  theme(axis.text.x  = element_text(face="bold",size=20),
        axis.text.y  = element_text(face="bold",size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=20))+
  stat_compare_means(comparisons=comparisons)
  #stat_compare_means(label.x=2.5,label.y=70,size=fontsize)
#mitsplot$layers[[3]]$aes_params$textsize <- fontsize
tiff("./figs/Figure5_mits.tiff",height=4000,width=4000,res=400)
mitsplot
dev.off()


ulceration <- as.data.frame.matrix(table(clinpath$Ulceration,clinpath$cluster))
ulceration_e <- excludeUnk(ulceration)
sum(colSums(ulceration_e))
colSums(ulceration_e)
ulceration_e
chiSquare(ulceration_e)

regression <- as.data.frame.matrix(table(clinpath$Regression,clinpath$cluster))
regression_e <- excludeUnk(regression)
write.csv(regression_e,"./data/regression.csv")
readRecode("regression2")

sats <- as.data.frame.matrix(table(clinpath$`Micro-satellites`,clinpath$cluster))

lvi <- as.data.frame.matrix(table(clinpath$LVI,clinpath$cluster))
lvi_e <- excludeUnk(lvi)
write.csv(lvi_e,"./data/lvi.csv")
readRecode("lvi2")

radgro <- as.data.frame.matrix(table(clinpath$`Radial Growth Phase`,clinpath$cluster))
radgro_e <- excludeUnk(radgro)
sum(colSums(radgro_e))
colSums(radgro_e)
radgro_e
chiSquare(radgro_e)

vergro <- as.data.frame.matrix(table(clinpath$`Vertical Growth Phase`,clinpath$cluster))
vergro_e <- excludeUnk(vergro)
sum(colSums(vergro_e))
colSums(vergro_e)
vergro_e
chiSquare(vergro_e)

typver <- as.data.frame.matrix(table(clinpath$`Type vertical growth`,clinpath$cluster))

### Old stuff ----

# Decision tree model
set.seed(1234)
library(party)
Formula <- as.factor(merged$cluster) ~ 
  as.factor(merged$Gender)+
  as.factor(merged$`Age at Dx`)+
  as.factor(merged$`Left/Right lesion location`) +
  as.factor(merged$`Lesion location`)+
  as.factor(merged$Symptoms) +
  as.factor(merged$Bleeding)+
  as.factor(merged$Itching)+
  as.factor(merged$`Pain/irritated`)+
  as.factor(merged$Growing)+
  as.factor(merged$Changing)+
  as.factor(merged$`New lesion`)+
  as.factor(merged$`Previously noted`)+
  as.factor(merged$`Sunburns/sun exposure`) +
  as.factor(merged$PUVA)+
  #as.factor(merged$`Tanning bed`)+
  as.factor(merged$`Prior NMSC`)+
  as.factor(merged$`Prior Melanoma`)+
  as.factor(merged$`family hx`)+
  as.factor(merged$`skin type`)+
  #as.factor(merged$`lesion type`)+
  as.factor(merged$`irregular color`)+
  #as.factor(merged$`irregular borders`)+
  as.factor(merged$`density of nevi`)+
  as.factor(merged$`clinical atypia`)+
  as.factor(merged$dermatoheliosis)

ctree <- ctree(Formula,data=merged)
plot(ctree)

age1 <- as.numeric(merged[merged$cluster==1,]$`Age at Dx`)
age2 <- as.numeric(merged[merged$cluster==2,]$`Age at Dx`)
age3 <- as.numeric(merged[merged$cluster==3,]$`Age at Dx`)
age4 <- as.numeric(merged[merged$cluster==4,]$`Age at Dx`)
age5 <- as.numeric(merged[merged$cluster==5,]$`Age at Dx`)
dati <- c(age1,age2,age3,age4,age5)
groups = factor()

merged$`Age at Dx` <- as.numeric(merged$`Age at Dx`)
fit <- lm(formula = `Age at Dx`~cluster,data=merged)

library(ggpubr)
ggboxplot(merged,x=merged$cluster,y=merged$`Age at Dx`)



chisq <- read_xlsx("./data/chisq.xlsx",sheet=1)


