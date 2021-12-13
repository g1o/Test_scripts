library("dplyr")  
library("pagoo") 
library("magrittr")
library(ggrepel)
library(Rtsne)

pan<-load_pangenomeRDS("/home/piruk/Documentos/projects/pos-DOC/LMA/pagoo_pangenome.rds")
#gffs <- list.files(path = "/home/piruk/Documentos/projects/pos-DOC/LMA/gffs/",
                   #pattern = "[.]gff$",
                   #full.names = TRUE)
#pan <- roary_2_pagoo("/home/piruk/Documentos/projects/pos-DOC/LMA/gene_presence_absence.csv", gffs, sep = "__")

#pan<-load_pangenomeRDS("/home/piruk/Documentos/projects/pos-DOC/LMA/pangenome_noEggNog.rds")

st <- data.frame(pan$summary_stats[-1, ])
st$Category <- factor(st$Category, levels = c("Cloud", "Shell", "Core"))
st$prop<-sapply(st$Number,function(x){x/sum(st$Number)*100})
st <- st %>%
  arrange(desc(Category)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
st
c("Cloud (Singletons)", "Shell ( 1 < strains <= 75)", "Core (Strains >= 76)")
mycols <- c("#6F73C2DD", "#FFC400FF", "#868686AA")
ggplot(st, aes(x = "", y = prop, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = lab.ypos, label = paste0(Number," (",signif(prop,digits=3),"%)")), color = "black") +
  scale_fill_manual(values = mycols,labels = c("Cloud (Singletons)", "Shell ( 1 < strains <= 75)", "Core (76 <= Strains <= 81)") ) +
  theme_void()


eggnog<-read.table("/home/piruk/Documentos/projects/pos-DOC/LMA/E3Cd2_Eggnog/clusters.emapper.annotations.tsv",
                   header=FALSE,  sep = "\t", comment.char = "#",na.strings = "-")
colnames(eggnog) <- c("cluster", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                      "max_annot_lvl", "COG_category", "Description", "Preferred_name",
                      "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KE",
                      "GG_Reaction", "KEGG_rclass", "BRITE", "CAZy",
                      "BiGG_Reaction", "PFAMs")
cluster_meta<-''
cluster_meta <- eggnog[, c("cluster", "COG_category", "KEGG_ko", "CAZy","EC","Description","Preferred_name")]

# Clean and parse the fields before feeding it to the pangenome
cluster_meta$COG_category <- cluster_meta$COG_category %>% strsplit("")
cluster_meta$EC <- cluster_meta$EC %>% strsplit(",")
cluster_meta$KEGG_ko <- cluster_meta$KEGG_ko %>%
  gsub("ko:", "", .) %>%
  strsplit(",")
cluster_meta$CAZy <- cluster_meta$CAZy %>% strsplit("")

# Add the metadata to the pangenome clusters
pan$add_metadata("cluster", cluster_meta)
# Now the object contains the new information in the clusters field:
pan$clusters

# search for genes with these COG category:
#grep("A",(pan$clusters$COG_category))
#COGs <- data.frame( All= table(unlist(pan$clusters$COG_category)) )
#COGs <- merge(COGs, data.frame(Core = table(unlist(pan$core_clusters$COG_category))) ,by=1 , all=T )
COGs <- data.frame(Core = table(unlist(pan$core_clusters$COG_category))) 
COGs <- merge(COGs, data.frame(Cloud = table(unlist(pan$cloud_clusters$COG_category))) ,by=1  ,all=T)
COGs <- merge(COGs, data.frame(Shell = table(unlist(pan$shell_clusters$COG_category))) ,by=1 ,  all=T)
COGs[is.na(COGs)]<-0
COGs<-reshape2::melt(COGs)
colnames(COGs)[1]<-"COG_Categories"
ggplot(COGs,aes(x=COG_Categories , y=value , fill=variable)) + 
  geom_bar(position = 'dodge',stat="identity")  
  #geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-1)

PCA<-prcomp(pan$pan_matrix[,names((pan$core_genes))], center = T, scale. = F)
autoplot(PCA, data = as.data.frame(pan$organisms), colour = c(2,rep(1,80)),loadings = F)
PCA<-prcomp(pan$pan_matrix[,names((pan$cloud_genes))], center = T, scale. = F)
autoplot(PCA, data = as.data.frame(pan$organisms), colour = c(2,rep(1,80)),loadings = F)
PCA<-prcomp(pan$pan_matrix[,names((pan$shell_genes))], center = T, scale. = F)
autoplot(PCA, data = as.data.frame(pan$organisms), colour = c(2,rep(1,80)),loadings = F)

gc()
all   <-sapply(as.character(pan$organisms$org),function(x){ unlist(pan$genes[,2])[grep(as.character(x),unlist(pan$genes[,2]))]} )
cloud <-sapply(as.character(pan$organisms$org),function(x){ unlist(pan$cloud_genes[,2])[grep(as.character(x),unlist(pan$genes[,2]))]} )
core  <-sapply(as.character(pan$organisms$org),function(x){ unlist(pan$core_genes[,2])[grep(as.character(x),unlist(pan$genes[,2]))]} )
shell <-sapply(as.character(pan$organisms$org),function(x){ unlist(pan$shell_genes[,2])[grep(as.character(x),unlist(pan$shell_genes[,2]))]} )
c=1;
nomes<-c("cloud","shell","core")
COG_tables<-lapply(list(cloud,shell,core),function(i){
  b<-sapply(i,function(x){  (as.matrix(table(
  unlist( pan$clusters[pan$clusters$cluster %in% names(x), "COG_category"]  )
  , deparse.level=0)  )[,1] )}  )
  merged_COG_table<-merge(b[1],b[[2]],by=0 ,all=T)
  colnames(merged_COG_table)<-c("Row.names",names(b[1:2]))
  for (x in 3:length(b)){
    merged_COG_table<-merge(merged_COG_table,b[[x]],by.x=1,by.y=0,all=T)
    colnames(merged_COG_table)[x+1] <- names(b[x])
  }
  rownames(merged_COG_table)<-merged_COG_table[,1]
  merged_COG_table<-merged_COG_table[,-1]
  merged_COG_table[is.na(merged_COG_table)] = 0
  #COG_tables[[nomes[c]]]<-merged_COG_table
  #c=c+1;
  merged_COG_table
} )
names(COG_tables)<-nomes
COG_tables

PCA<-prcomp(t(COG_tables$shell), center = T, scale. = F) 
autoplot(PCA, data = as.data.frame(pan$organisms), colour = c(2,rep(1,80)),loadings = F,label=T)
summary(PCA)

#duprm<-unique(t(COG_tables$shell))

set.seed(1)
tsne_plots<-lapply(seq_along(COG_tables),function(x){
  duprm<-unique(t(COG_tables[[x]]))
  
    tsne.norm = Rtsne(duprm,perplexity=7,max_iter=5000,initial_dims = 15)
  
  info.norm = data.frame( tsne1 = tsne.norm$Y[, 1], tsne2 = tsne.norm$Y[, 2] ) 
  info.norm$names<-row.names(duprm)
  hc.norm = hclust(dist(tsne.norm$Y),method="ward.D")
  info.norm$hclust = factor(cutree(hc.norm, 7))
  hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1, 
                                                         tsne2) %>% summarize_all(mean)
  cent_plot<-ggplot(info.norm, aes(x = tsne1, 
                      y = tsne2 ,
                      colour=c('red',rep('black',dim(duprm)[1]-1))) ) + 
  geom_point(alpha = 0.3) +
  theme_bw()+
  guides(colour = "none") 

h_plot<-ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) + 
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_label_repel(aes(label = hclust), 
                  data = hc.norm.cent) +
  guides(colour = "none" )+ 
  ggtitle(paste0("TSNE: ward.D clusters from COG categories of ", names(COG_tables)[x]," Genes")) 
  list(cent_plot,h_plot)
} )

tsne_plots[[2]][[1]]



#get all genes from a strain
E3Cd2_genes<-names(unlist(pan$genes[,2])[grep("E3Cd2",unlist(pan$genes[,2]))])
pan$genes[E3Cd2_genes]
write.table(pan$clusters$cluster,file="clusters.tsv")


#############try grouping with pcoa
jac_accessory_dist <- proxy::dist(pan$pan_matrix[,names(pan$shell_genes)] , by_rows = TRUE, method = "Jaccard")
mds <- cmdscale(jac_accessory_dist,eig = T)
colnames(mds$points) <- c("x","y")
ggplot(mds$points,aes(x=x,y=y,label=rownames(mds$points))) + geom_point()

jac_accessory_dist <- proxy::dist(pan$pan_matrix[, bacmet_clusters_blast_hit[,1] ] , by_rows = TRUE, method = "Jaccard")
mds <- cmdscale(jac_accessory_dist,eig = T)
colnames(mds$points) <- c("x","y")
g1 <- mds$points["E3Cd2",]  
ggplot(mds$points,aes(x=x,y=y,label=rownames(mds$points))) +
  geom_point() + geom_text(vjust=2) 
  #geom_text(data=t(g1),label="E3Cd2",vjust=1.4) 
rowSums(pan$pan_matrix[, bacmet_clusters_blast_hit[,1] ])



