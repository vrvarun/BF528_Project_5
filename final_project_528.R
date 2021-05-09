#Varun Raghuraman, Final Project Script

#Analyst Portion
gene_exp<-read.table('/projectnb/bf528/users/swiss_cheese2/project_2/cuffdiff_out/gene_exp.diff')
sort<-gene_exp[order(gene_exp$V13),]
#read in differentially expressed genes, sort by lowest q-value

Top_Ten<-head(sort,10)
#get top ten differentially expressed genes
hist1 <- gsub(",", "", gene_exp$V10)  
#clear list expression values of commas
#tail(gene_exp[order(gene_exp$V10, decreasing = TRUE),],50)
#tail(gene_exp[order(gene_exp$V10),],950)
hist1 <- gsub("inf", "10", gene_exp$V10) 
#When the prior code is run, the range of values is between -0.000332148 8.57949
#to display the +/- infinity values as well, the characters "inf"
#are turned to 10


hist(as.numeric(as.character(hist1[2:length(hist1)])),main='Unfiltered Genes',xlab='log2(fold_change)',breaks=50)
#histogram creations

sort.sub<-subset(sort,sort$V14=='yes')
#sorts by fold change
hist2 <- gsub(",", "", sort.sub$V10)  
#removes commas
#hist2 <- gsub("inf", "10", sort.sub$V10) 
#replaces 10
hist(as.numeric(as.character(hist2[1:length(hist2)])),main='Filtered Genes',xlab='log2(fold_change)',breaks=50)

sort_p<-subset(sort.sub,sort.sub$V12<0.01)
#p value significance level filter

up_reg<-subset(sort_p,sort.sub$V10>0)

down_reg<-subset(sort_p,sort.sub$V10<0)
#upregulation and downregulation from fold change

write.table(up_reg$V3, "up_reg.txt", col.names=FALSE, row.names=FALSE)
write.table(down_reg$V3, "down_reg.txt", col.names=FALSE,row.names=FALSE)


#Biologist portion

gene_fpkm_tracking<-read.table('/projectnb/bf528/users/swiss_cheese2/project_2/cuffdiff_out/genes.fpkm_tracking')

other_gene_fpkm_tracking<-read.table('/project/bf528/project_2/data/fpkm_matrix.csv')


Ad_1<- read.table('/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking') 
Ad_mito<-subset(Ad_1,Ad_1$V5 %in% c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh') )

Ad_2<- read.table('/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking') 
P0_2<- read.table('/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking')   
P4_1<- read.table('/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking')   
P4_2<- read.table('/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking')   
P7_1<- read.table('/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking')   
P7_2<- read.table('/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking')
P0<- read.table('/projectnb/bf528/users/swiss_cheese2/project_2/cuffdiff_out/genes.fpkm_tracking') 
#reads all fpkm
found<-subset(P0,P0$V5 %in% head(sort$V3,10))
#gets fpkm for top ten genes 
Top_Ten$V15<-found$V10
#attaches column to table
P0_2_mito<-subset(P0_2[order(P0_2$V5),],P0_2[order(P0_2$V5),]$V5 %in% c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh') )
P0_mito<-subset(P0[order(P0$V5),],P0[order(P0$V5),]$V5 %in% c('Mpc1', 'Prdx3', 'Acat1', 'Echs1', 'Slc25a11', 'Phyh') )
#gets genes and id's for the mitochondrial genes
All_subset_mito<-subset(other_gene_fpkm_tracking, other_gene_fpkm_tracking$V1 %in% P0_2_mito$V1 )
All_subset_mito$V9<-P0_mito$V10
All_subset_mito$V10<-P0_mito$V5
#adds fpkm values to total from the combined table

P0_2_sarc<-subset(P0_2[order(P0_2$V5),],P0_2[order(P0_2$V5),]$V5 %in% c('Csrp3', 'Tcap', 'Cryab', 'Pdlim5', 'Pygm', 'Myoz2', 'Des') )
P0_sarc<-subset(P0[order(P0$V5),],P0[order(P0$V5),]$V5 %in% c('Csrp3', 'Tcap', 'Cryab', 'Pdlim5', 'Pygm', 'Myoz2', 'Des') )
All_subset_sarc<-subset(other_gene_fpkm_tracking, other_gene_fpkm_tracking$V1 %in% P0_2_sarc$V1 )
All_subset_sarc$V9<-P0_sarc$V10
All_subset_sarc$V10<-P0_sarc$V5
#adds fpkm values to total from the combined table

P0_2_cc<-subset(P0_2[order(P0_2$V5),],P0_2[order(P0_2$V5),]$V5 %in% c('Cdc7','Cdc27','E2f8','Bora','Cdk7','Cdc45','Cdc26','Rad51','Cdc6','Aurkb','Cdc23') )
P0_cc<-subset(P0[order(P0$V5),],P0[order(P0$V5),]$V5 %in% c('Cdc7','Cdc27','E2f8','Bora','Cdk7','Cdc45','Cdc26','Rad51','Cdc6','Aurkb','Cdc23') )
All_subset_cc<-subset(other_gene_fpkm_tracking, other_gene_fpkm_tracking$V1 %in% P0_2_cc$V1 )
All_subset_cc$V9<-P0_cc$V10
All_subset_cc$V10<-P0_cc$V5
#adds fpkm values to total from the combined table
write.csv(All_subset_cc,'Cell_Cyle.csv', row.names=FALSE)
write.csv(All_subset_sarc,'Sarc.csv', row.names=FALSE)
write.csv(All_subset_mito,'Mito.csv', row.names=FALSE)
#create csv files of the cell data
#due to preference for layout, and relative simplicity of data
#python was used to create the line plots


install.packages('ggplot2')
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
#plot libraries used in prior projects; added here, not all used

P0map<-select(subset(P0, P0$V5 %in% head(sort_p,1000)$V3),V5,V10)
#gets the gene names and fpkm for P0; used to find the matches for the rest
P02<-select(subset(P0_2, P0_2$V5 %in% P0map$V5),V10)
P41<-select(subset(P4_1, P4_1$V5 %in% P0map$V5),V10)
P42<-select(subset(P4_2, P4_2$V5 %in% P0map$V5),V10)
P71<-select(subset(P7_1, P7_1$V5 %in% P0map$V5),V10)
P72<-select(subset(P7_2, P7_2$V5 %in% P0map$V5),V5,V10)
P01<-select(subset(P0map, P72$V5 %in% P0map$V5),V10)
AD1<-select(subset(Ad_1, Ad_1$V5 %in% P0map$V5),V10)
AD2<-select(subset(Ad_2, Ad_2$V5 %in% P0map$V5),V10)
#uses map to get fpkm values

heatmapf<-head(P72,1000)
heatmapf$Name<-head(P72$V5,1000)
heatmapf$P0_1<-head(P01$V10,1000)
heatmapf$P0_2<-head(P02$V10,1000)
heatmapf$P4_1<-head(P41$V10,1000)
heatmapf$P4_2<-head(P42$V10,1000)
heatmapf$P7_1<-head(P71$V10,1000)
heatmapf$P7_2<-head(P72$V10,1000)
heatmapf$AD_1<-head(AD1$V10,1000)
heatmapf$AD_2<-head(AD2$V10,1000)
final_heat<-heatmapf %>% select(-Name,-V10,-V5)
#adds columns to final for heatmap; removes extras introduced to doublecheck
final_heat_frame<-data.matrix(head(final_heat,1000))
#double-checks size
final_heat_frame<- na.omit(final_heat_frame)
#removes any na's
colors = brewer.pal(n = 9, name = "GnBu")
colors = colorRampPalette(colors)(10)
colors = rev(colors)
#color settings-sets, theme, amx colors for chart

heatmap <- pheatmap(t(final_heat_frame), scale = "row", color = colors,fontsize_row = 4,border_color = NA, clustering_distance_rows="euclidean",clustering_distance_cols="euclidean", main = "Top 1000 Expressed Genes ")
#heatmap applied
