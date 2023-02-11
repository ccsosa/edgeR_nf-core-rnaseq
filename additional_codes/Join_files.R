require(dplyr);require(data.table);require(pheatmap);library(caret)
#Directories
#where are the salmon files
################################################################################
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
################################################################################

folders <- c(1,2,3,5,6,7,8,9,10)
folder_name <- folders[[1]]

data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
vhRR_dir <- "/scratch/bis_klpoe/chsos/analysis/vHHR/input"
x1 <- read.table(paste0(data_dir,"/","salmon.merged.gene_tpm.tsv"),header = T)
x1$tx <- NULL
x1$gene_name <- NULL
print(paste(nrow(x1)," / ",1))

for(i in 2:length(folders)){
  message(i)
  folder_name <- folders[[i]]
  data_dir <- paste0("/scratch/bis_klpoe/chsos/analysis/RESULTS_NF/salmon_",folder_name)
  x <- read.table(paste0(data_dir,"/","salmon.merged.gene_tpm.tsv"),header = T)
  x$tx <- NULL
  x1$gene_name <- NULL
  
  print(identical(x$gene_id,x1$gene_id))
  x$gene_id <- NULL

  x1 <- cbind(x1,x)
  print(paste(nrow(x1)," / ",i))
  rm(x)
};rm(i)
x1$gene_name <- NULL

row.names(x1) <- x1$gene_id
write.table(x1,"/scratch/bis_klpoe/chsos/analysis/atlas_drought_total.tsv",sep = "\t",na = "",row.names = F,col.names =T,quote = F)
write.table(x1,paste0(vhRR_dir,"/atlas_drought_total.tsv"),sep = "\t",na = "",row.names = F,col.names = F,quote = F)

###########################
# require(data.table)
#contrasts
sets = c("SRP356530/LEAVES_D-LEAVES_C", 
         "SRP098756/CROWN_D-CROWN_C",
         "SRP098756/LEAVES_D-LEAVES_C",
         "SRP098756/ROOTS_D-ROOTS_C",
         #"SRP098756/ROOTS_D-LEAVES_D",
         "SRP257474/CK0H-ABA6H",
         "SRP257474/DT6H-ABA6H",         
         "SRP072216/LEAVES_D-LEAVES_C",
         "SRP072216/LEAVES_TABA-LEAVES_D")
pval <- 0.05
data_dir <- "/scratch/bis_klpoe/chsos/analysis/DEG"         
sets_st_cont <- strsplit(sets,"/")

##loading all DEG for all contrasts
x_deg <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_filtered.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg <- do.call(rbind,x_deg)
#
#geeting all genes accross all studies
x_unique_genes <- as.data.frame(matrix(nrow=length(unique(x_deg$V1)),ncol=length(sets)+1))
x_unique_genes$V1 <- unique(x_deg$V1)
colnames(x_unique_genes) <- c("genes",sets)
row.names(x_unique_genes) <- unique(x_deg$V1)

pb <-
  utils::txtProgressBar(min = 0,
                        max = nrow(x_unique_genes),
                        style = 3)

#Presence of a gene in a set of contrasts
for(i in 1:nrow(x_unique_genes)){
  utils::setTxtProgressBar(pb, i)
  x_i <- x_deg[which(x_deg$V1==x_unique_genes$genes[[i]]),]
  x_unique_genes[i,-1] <- sets %in% x_i$group # transforming boolean in 1 and 0s
  x_unique_genes[i,-1] <- (x_unique_genes[i,-1]*1)
};rm(i)

close(pb)


#unique gene core 
x_unique_genes$total<- rowSums(x_unique_genes[,-1])
#subsetting
unique_genes <- unique(x_unique_genes$genes)
subset_pan <- x1[which(x1$gene_id %in% unique_genes),]
write.table(subset_pan,paste0(vhRR_dir,"/atlas_drought_pan.tsv"),sep = "\t",na = "",row.names = F,col.names = F,quote = F)
################################################################################
x_unique_genes_dif_contrasts <- x_unique_genes
x_unique_genes_dif_contrasts <- x_unique_genes_dif_contrasts[which(x_unique_genes_dif_contrasts$total>1),]
#core genes 8 contrasts
#all_contr_genes <- x_unique_genes$genes[which(x_unique_genes$total==length(sets))] #only one gene
#core genes 7 contrasts
all_contr_genes2 <- x_unique_genes$genes[which(x_unique_genes$total>=6)]
#core genes 6 contrasts
all_contr_genes3 <- x_unique_genes$genes[which(x_unique_genes$total>=5)]



########################################################################################
#52 genes

subset1 <- x1[which(x1$gene_id %in% all_contr_genes2),]
subset1 <- subset1[rowSums(subset1[,-1]) > 0,]; subset1a <- subset1
row.names(subset1) <- subset1$gene_id; subset1$gene_id <- NULL;subset1$gene_name <- NULL

#write.table(subset1,"/scratch/bis_klpoe/chsos/analysis/atlas_drought_52.tsv",sep = "\t",na = "",row.names = T,col.names = F)
write.table(subset1,paste0(vhRR_dir,"/atlas_drought_52.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)
#subset1[subset1==0] <- NA
#subset1 <- scale(subset1)

#min max
process <- preProcess(subset1, method=c("range"))
norm_scale <- predict(process, subset1)
# z-score
# preproc1 <- preProcess(subset1, method=c("center", "scale"))
# norm_scale <- predict(preproc1, subset1)

# z-score
# for(i in 1:ncol(subset1)){
#   subset1[,i] <- (subset1[,i] - mean(subset1[,i],na.rm = T)) / sd(subset1[,i],na.rm = T)  
# };rm(i)

Breaks <- seq(floor(min(norm_scale)),ceiling(max(norm_scale)), by = 0.2);Breaks[length(Breaks)] <- max(norm_scale) #by=1
#Breaks <-  round(Breaks,1)
Breaks2 <- as.character(Breaks); Breaks2[length(Breaks2)] <- "Min-Max Normalization (TPM)\n"#"Z-Score Normalization(TPM)\n"

p <- pheatmap::pheatmap(norm_scale,cluster_rows = T,
                        cluster_cols = F,
                        fontsize = 18,
                        annotation_legend = F,#F,
                        labels_row =NULL,
                        #annotation_col = ann_col,
                        # annotation_colors =  list(genelist=c(DOWNREGULATED="blue",#"#006EC9",
                        #                                      UPREGULATED="red")),#"#FF00FF")),
                        drop_levels=F,
                        # display_numbers = matrix_top_go_final_2,#T,#TRU1E,
                        #number_color = "white",
                        fontsize_row = 12,#12
                        fontsize_col =  12,#12
                        border_color="gray98",
                        clustering_distance_rows="correlation",
                        fontsize_number=14,#10
                        cellwidth = 20, #48
                        angle_col = 45,
                        cellheight = 15, #16
                        #pvalue
                        legend_breaks = Breaks,
                        legend_labels = Breaks2,
                        #log10
                        #legend_breaks = c(0,15,30,45,60,75,90,93),
                        #legend_labels = c("0","15","30","45","60","75","90","-log10(p-value)\n"),
                        #color
                        #color = hcl.colors(16, "Greens"), #Oranges
                        color = hcl.colors(6, "BluYl"),
                        #color = hcl.colors(8, "BrBg"),
                        #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
                        legend = T,
                        annotation_names_col = F,
                        na_col = "gray",
                        #number_format= "%.2f",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_6_or_more_contrasts.pdf"),
                  width = 30,#20,
                  height = 30 )

rm(norm_scale,preproc1)


########################################################################################
#419 genes

subset2 <- x1[which(x1$gene_id %in% all_contr_genes3),]; subset2a <- subset2
subset2 <- subset2[rowSums(subset2[,-1]) > 0,]; subset2a <- subset2

row.names(subset2) <- subset2$gene_id; subset2$gene_id <- NULL;subset2$gene_name <- NULL
#write.table(subset1,"/scratch/bis_klpoe/chsos/analysis/atlas_drought_419.tsv",sep = "\t",na = "",row.names = T,col.names = F)
write.table(subset2,paste0(vhRR_dir,"/atlas_drought_419.tsv"),sep = "\t",na = "",row.names = T,col.names = F,quote = F)
#subset1 <- scale(subset1,center = T,scale = T)
#min max
process <- preProcess(subset1, method=c("range"))
norm_scale <- predict(process, subset1)
# z-score
# preproc1 <- preProcess(subset2, method=c("center", "scale"))
# norm_scale <- predict(preproc1, subset2)
# z-score
# for(i in 1:ncol(subset1)){
#   subset1[,i] <- (subset1[,i] - mean(subset1[,i],na.rm = T)) / sd(subset1[,i],na.rm = T)  
# };rm(i)
Breaks <- seq(floor(min(norm_scale)),ceiling(max(norm_scale)), by = 0.2);Breaks[length(Breaks)] <- max(norm_scale) #by=1
#Breaks <-  round(Breaks,1)
Breaks2 <- as.character(Breaks); Breaks2[length(Breaks2)] <- "Min-Max Normalization (TPM)\n"#"Z-Score Normalization(TPM)\n"


#subset1[subset1==0] <- NA
p <- pheatmap::pheatmap(norm_scale,cluster_rows = T,
                        cluster_cols = F,
                        fontsize = 18,
                        annotation_legend = F,#F,
                        labels_row =NULL,
                        #annotation_col = ann_col,
                        # annotation_colors =  list(genelist=c(DOWNREGULATED="blue",#"#006EC9",
                        #                                      UPREGULATED="red")),#"#FF00FF")),
                        drop_levels=F,
                        # display_numbers = matrix_top_go_final_2,#T,#TRU1E,
                        #number_color = "white",
                        fontsize_row = 0.01,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        fontsize_number=14,#10
                        cellwidth = 40, #48
                        angle_col = 45,
                        cellheight = 13, #16
                        #pvalue
                        legend_breaks = Breaks,
                        legend_labels = Breaks2,
                        #color
                        #color = hcl.colors(16, "Greens"), #Oranges
                        color = hcl.colors(10, "BluYl"),
                        #color =rev(hcl.colors(6, "BluYl")),#hcl.colors(10, "BluYl"),
                        legend = T,
                        annotation_names_col = F,
                        na_col = "gray",
                        #number_format= "%.2f",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_5_or_more_contrasts.pdf"),
                  width = 43,
                  height = 90 )
rm(norm_scale,preproc1)
########################################################################################
########################################################################################
########################################################################################
########################################################################################
#Preparing LFC heatmap

x_deg_total <- lapply(1:length(sets),function(i){
  message(i)
  x <- fread(paste0(data_dir,"/",sets_st_cont[[i]][1],"/csv/",
                    "glmQLFTest_",sets_st_cont[[i]][2],"_pval_",
                    pval,"_fulltable.csv"))
  x$status <- NULL
  x$group <- sets[[i]]
  return(x)
})
x_deg_total <- do.call(rbind,x_deg_total)


########################################################################################
#52 genes

subset1a <- x_deg_total[x_deg_total$V1 %in% subset1a$gene_id,]
LFC_1a <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subset1a$V1))))
LFC_1a[,1] <- unique(subset1a$V1)
colnames(LFC_1a) <- c("gene",sets)
LFC_1a_p_val <- LFC_1a

for(i in 2:ncol(LFC_1a)){

  x_i <- subset1a[which(subset1a$group==sets[i-1]),]
  for(j in 1:nrow(LFC_1a)){
      LFC_1a[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),2])
      LFC_1a_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_1a$gene[[j]]),6])
  };rm(j)
};rm(i)

for(i in 2:ncol(LFC_1a)){
  LFC_1a_p_val[,i] <- as.numeric(LFC_1a[,i])
  LFC_1a[,i] <- as.numeric(LFC_1a_p_val[,i])
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]>=pval)] <- ""
  LFC_1a_p_val[,i][which(LFC_1a_p_val[,i]<pval)] <- "*"
  LFC_1a_p_val[,i][which(is.na(LFC_1a_p_val[,i]))] <- ""
  LFC_1a[,i][which(is.na(LFC_1a[,i]))] <- 0
};rm(i)




row.names(LFC_1a) <- LFC_1a$gene;
row.names(LFC_1a_p_val) <- LFC_1a$gene;
LFC_1a$gene <- NULL;LFC_1a_p_val$gene <- NULL


ann_col <- data.frame(Tissue=c("leaves","crown","leaves","roots","three_leaf","three_leaf","leaves","leaves"))
row.names(ann_col) <- sets
p <- pheatmap::pheatmap(LFC_1a,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Tissue=c(leaves="purple",
                                                           crown="brown",
                                                           roots= "red",
                                                           three_leaf="gray"
                                                          )),
                        drop_levels=F,
                        fontsize_row = 12,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 40, #48
                        angle_col = 45,
                        cellheight = 13, #16
                        legend_breaks = c(-10,-7.5,-5,-2.5,0,2.5,5,7.5,
                                          max(LFC_1a,na.rm = T)),
                        legend_labels = c("-10","-7.5","-5","-2.5","0","2.5","5","7.5","log2 fold change\n"), 
                        color = hcl.colors(6, "BluYl"),
                        legend = T,
                        annotation_names_col = T,
                        na_col = "gray",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_6_or_more_contrasts_LFC.pdf"),
                  width = 18,
                  height = 30)

########################################################################################
#419 genes
subset2a <- subset2
subset2a <- x_deg_total[x_deg_total$V1 %in% row.names(subset2a),]
LFC_2a <- as.data.frame(matrix(ncol=length(sets)+1,nrow=length(unique(subset2a$V1))))
LFC_2a[,1] <- unique(subset2a$V1)
colnames(LFC_2a) <- c("gene",sets)
LFC_2a_p_val <- LFC_2a

for(i in 2:ncol(LFC_2a)){
  
  x_i <- subset2a[which(subset2a$group==sets[i-1]),]
  for(j in 1:nrow(LFC_2a)){
    LFC_2a[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_2a$gene[[j]]),2])
    LFC_2a_p_val[j,i] <- as.numeric(x_i[which(x_i$V1==LFC_2a$gene[[j]]),6])
  };rm(j)
};rm(i)

for(i in 2:ncol(LFC_2a)){
  LFC_2a_p_val[,i] <- as.numeric(LFC_2a[,i])
  LFC_2a[,i] <- as.numeric(LFC_2a_p_val[,i])
  LFC_2a_p_val[,i][which(LFC_2a_p_val[,i]>=pval)] <- ""
  LFC_2a_p_val[,i][which(LFC_2a_p_val[,i]<pval)] <- "*"
  LFC_2a_p_val[,i][which(is.na(LFC_2a_p_val[,i]))] <- ""
  LFC_2a[,i][which(is.na(LFC_2a[,i]))] <- 0
};rm(i)




row.names(LFC_2a) <- LFC_2a$gene;
row.names(LFC_2a_p_val) <- LFC_2a$gene;
LFC_2a$gene <- NULL;LFC_2a_p_val$gene <- NULL


ann_col <- data.frame(Tissue=c("leaves","crown","leaves","roots","three_leaf","three_leaf","leaves","leaves"))
row.names(ann_col) <- sets
p <- pheatmap::pheatmap(LFC_2a,cluster_rows = T,
                        cluster_cols = T,
                        fontsize = 18,
                        annotation_legend = T,#F,
                        labels_row =NULL,
                        annotation_col = ann_col,
                        annotation_colors = list(Tissue=c(leaves="purple",
                                                          crown="brown",
                                                          roots= "red",
                                                          three_leaf="gray"
                        )),
                        drop_levels=F,
                        fontsize_row = 0.0001,#12
                        fontsize_col =  12,#12
                        border_color=NULL,
                        clustering_distance_rows="correlation",
                        clustering_distance_cols="correlation",
                        fontsize_number=14,#10
                        cellwidth = 40, #48
                        angle_col = 45,
                        cellheight = 4, #16
                        legend_breaks = c(-12,-10,-7.5,-5,-2.5,0,2.5,5,7.5,
                                          max(LFC_2a,na.rm = T)),
                        legend_labels = c("-12","-10","-7.5","-5","-2.5","0","2.5","5","7.5","log2 fold change\n"), 
                        color = hcl.colors(6, "BluYl"),
                        legend = T,
                        annotation_names_col = T,
                        na_col = "gray",
                        silent = F
                        
)

save_pheatmap_pdf(x = p,
                  filename =paste0(data_dir,"/","HEATMAP_5_or_more_contrasts_LFC.pdf"),
                  width = 20,
                  height = 27)

