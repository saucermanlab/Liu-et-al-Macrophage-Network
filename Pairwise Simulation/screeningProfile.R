# Code orginally from Astor Liu screeningprofile.R (11/3/2017)
# Jingyuan Zhang reduced version to only generate activity heatmap for stimuli pairs for Figure 4
# 5/15/2018

require(tidyverse)
require(scales)

#  setwd(
#   "/Users/jingyuan/Documents/Academic/research/Macrophage Model/Astor macrophage/macmodel_test_original6.3p_JZ_edits_only_necessary_files/Pairwise Simulation"
# )

m1m <- read_csv("M0vsM1_matched.csv")
m1 <- m1m[!is.na(m1m$log2FC_2over1), c(2, 15)]
require(scales)
# m1[,2] <- rescale_mid(m1$log2FC_2over1, to =c(-1,1), mid = 0)
m2m <- read_csv("M0vsM2_matched.csv")
m2 <- m2m[!is.na(m2m$log2FC_2over1), c(2, 15)]
# m2[,2] <- rescale_mid(m2$log2FC_2over1, to =c(-1,1), mid = 0)


data <- 'nontemp'

eigenvalues1 <- list()
eigenvalues.single1 <- list()
# Case insensitive join function
insensitive <- function(fun = inner_join) {
  new_fun <- fun
  body(new_fun) <- substitute({
    by <- dplyr:::common_by(by, x, y)
    
    tmp_by_x <- paste0("_", by$x, "_")
    tmp_by_y <- paste0("_", by$y, "_")
    for (i in seq_along(by$x)) {
      x[[tmp_by_x[[i]]]] <- tolower(x[[by$x[[i]]]])
      y[[tmp_by_y[[i]]]] <- tolower(y[[by$y[[i]]]])
      y[[by$y[[i]]]] <- NULL
    }
    
    res <- fun(x, y, list(x = tmp_by_x, y = tmp_by_y))
    res[tmp_by_x] <- list(NULL)
    
    res
  })
  
  new_fun
}

`%not in%` <-
  function (x, table)
    is.na(match(x, table, nomatch = NA_integer_))


input <- 0.7

activity.nontemp <- read_tsv("./simulation results/Screening_output.txt")
combo.nontemp <- read_tsv("./simulation results/Screening_percentMatch.txt")
combo.nontemp <- combo.nontemp$Row

activity.all <- list()
for (l in 1:45) {
  activity.all[[combo.nontemp[l]]] <-
    activity.nontemp[(1 + (l - 1) * 137):(137 * l), ]
}
require(reshape2)


activity <- activity.nontemp[, 2:309]
timepoint <- c(12, 16, 20, 32) # -8
# timepoint <- 11
cond <- 45
combo <- combo.nontemp
single.names <- c('LPS+',
                  'IFNg+',
                  'IFNb+',
                  'TNFa+',
                  'IL1+',
                  'IL4+',
                  'IL6+',
                  'IL10+',
                  'IL12+',
                  'control')


k <- timepoint[1]
nodes <- activity[1:137, 1]
activity2 <- as.vector(activity[, k])
activity.stable <- matrix(t(activity2), nrow = 137, ncol = cond)
colnames(activity.stable) <- combo
rownames(activity.stable) <- t(nodes)

activity1 <- as.vector(activity[, 2])
activity.int <- matrix(t(activity1), nrow = 137, ncol = cond)
colnames(activity.int) <- combo
rownames(activity.int) <- t(nodes)

activity.change <- activity.stable - activity.int
activity.logfc <- log2(activity.stable / activity.int)
activity.change[is.infinite(activity.change)] <- 0
activity.change[is.na(activity.change)] <- 0

# Diminish the duplicated combos
# test <- round(activity.stable,digits = 1) # Remove identical responses
# test2 <- activity.stable[,!duplicated(test,MARGIN = 2)] # Normalized expression level
test2 <- activity.stable
test2 <-
  cbind(test2, control = activity.int[, 1]) # [,!colnames(test2) %in% 'TNFa+LPS']
# test2 <- test2[-grep(pattern = 'CD36',rownames(test2)),]
uniq <- dim(test2)[2]

# testFC <- activity.logfc[,!duplicated(test,MARGIN = 2)] # Activity change
# testFC <- cbind(testFC,control=0,M1Aeq=m1m$log2FC_2over1,M2Seq=m2m$log2FC_2over1) # [,!colnames(testFC) %in% 'TNFa+LPS']
# testFC[is.na(testFC)] <- 0
# testFC <- testFC[-grep(pattern = 'CD36',rownames(testFC)),]
# test2 <- activity.stable
# testFC <- activity.change

# Scale fold change data
# testFC.s <- testFC
# testFC.s[,1:uniq] <- t(scale(t(testFC.s[,1:uniq]),center = FALSE))
# testFC.s[,uniq+1:uniq+2] <- t(scale(t(testFC.s[,uniq+1:uniq+2]),center = FALSE))
# testFC.z <- testFC
# testFC.z[,1:uniq] <- t(scale(t(testFC.z[,1:uniq]),center = FALSE,
#                     scale = apply(t(testFC.z[,1:uniq]), 2, sd, na.rm = TRUE)))
# testFC.z[,uniq+1:uniq+2] <- t(scale(t(testFC.z[,uniq+1:uniq+2]),center = FALSE,
#                            scale = apply(t(testFC.z[,uniq+1:uniq+2]), 2, sd, na.rm = TRUE)))

# FC.rna <- insensitive(left_join)(nodes,m1,by=c("yTran1"="id"))
# FC.rna <- insensitive(left_join)(FC.rna,m2,by=c("yTran1"="id"))
# FC.rna[is.na(FC.rna)] <- 0
# colnames(FC.rna) <- c("id","M1","M2")
#
# testFC.rna <- cbind(testFC, FC.rna[,2:3])

require(gplots)
redscale <- colorpanel(50, low = "black", high = "white")
# require(RColorBrewer)
pal <- c(
  colorpanel(30, low = "blue", high = "cyan"),
  colorpanel(30, low = "cyan", mid = "green", high = "yellow"),
  colorpanel(40, low = "yellow", mid = "orange", high = "red")
)



## -------------------------Acticity level---------------------------------------


mydata <- t(test2)
mydata.single <- mydata[single.names, ]
rownames(mydata.single) <-
  c('LPS',
    "IFNg",
    "IFNb",
    "TNFa",
    "IL1",
    "IL4",
    "IL6",
    "IL10",
    "IL12",
    'control')

# Clustering - hierachical
hc <- hclust(dist(mydata), method = 'ward.D')
hc.single <- hclust(dist(mydata.single), method = 'ward.D')
hc2 <- hclust(dist(test2)) #,method = 'ward.D'

cl <- 6

# Groups of stimuli combo
groups <- cutree(hc, k = cl) # cut tree into 7 clusters
# groups <- dek$cluster
a <- as.data.frame(cbind(hc$labels,groups))
write_csv(a[order(a$groups),],
          paste0("./simulation results/Screening_hclust1_stimuli-", input, data, k, ".csv"))

# Groups of model nodes
# Note that this is not what is used in the column dendrogram in the heatmap
group.hc2 <- cutree(hc2, k = 8)
group.hc2 <- cbind(id = names(group.hc2), group.hc2)
write_csv(as.data.frame(group.hc2[order(group.hc2[,2]),]),
          paste0("./simulation results/Screening_hclust2_nodes-", input, data, k, ".csv"))

group.hc2m <- cutree(hc2, k = 24)
group.hc2m[group.hc2m == 20] <- 6 # Combine IL4 module
group.hc2m[group.hc2m == 17 |
             group.hc2m == 20] <- 12 # Combine MAPKs module group.hc2m==14|
group.hc2m[group.hc2m == 4 |
             group.hc2m == 19] <- 2 # Combine LPS-IL1 module
group.hc2m[group.hc2m == 9] <- 3 # Combine IL6-IL10 module


mydata.module.1 <-
  cbind(as.data.frame(t(mydata)), group = group.hc2m)
mydata.module.rna <-
  mydata.module.1[c(rownames(mydata.module.1)[grep(pattern = 'mrna', rownames(mydata.module.1))],
                    'CREB',
                    'CEBPbeta',
                    'PI3K',
                    'AKT',
                    'SOCS3'),
                  c(
                    "LPS+",
                    "TNFa+",
                    "IFNg+",
                    "LPS+IL4",
                    "TNFa+IL4",
                    "IFNg+IL4",
                    "IL4+",
                    "control",
                    "group"
                  )]
sortid <- sort(mydata.module.rna$group, index.return = TRUE)
mydata.module.rna <- mydata.module.rna[sortid$ix, ]
mydata.module.rna.dim <-
  mydata.module.rna[!rownames(mydata.module.rna) %in% c(
    'AMAC1_mrna',
    'CCL5_mrna',
    'CD23_mrna',
    'KLF4_mrna',
    'IL12_mrna',
    'IFNb_mrna',
    'CXCL10_mrna',
    "IL15_mrna",
    'Mcpt1_mrna',
    'IL8_mrna',
    'GMCSF_mrna',
    'CRP_mrna'
  ), ]
group.rna <- as.factor(mydata.module.rna$group)
group.rna.dim <- as.factor(mydata.module.rna.dim$group)


#comprehensive activity heatmap
png(
  paste0("./plots/Screening", input, data, k, "_label.png"),
  # create PNG for the heat map
  width = 12 * 300,
  # 5 x 300 pixels
  height = 4 * 300,
  res = 300,
  # 300 pixels per inch
  pointsize = 6
)        # smaller font size

hm <- heatmap.2(
  mydata,
  # cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
  lmat = rbind(c(5, 5, 4), c(3, 1, 2)),
  # 1. heatmap 2. row dendro 3.col dendro 4. key
  lwid = c(1, 0.5, 4),
  lhei = c(1.5, 4),
  notecol = "black",
  # change font color of cell labels to black
  density.info = "none",
  # turns off density plot inside color legend
  trace = "none",
  # turns off trace lines inside the heat map;
  margins = c(6,9),
  # widens margins around plot
  col = redscale,
  # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
  breaks = c(seq(0, 1, length.out = 51)),
  # srtCol = 45,colCol = rainbow(8)[as.factor(group.hc2[,2])],
  labRow = '',
  labCol = colnames(mydata),
  # cexRow = 1,colRow = rainbow(cl)[groups],
  RowSideColors = rainbow(cl)[groups],
  # Measurement 7-10: red
  dendrogram = "both",
  # only draw a row dendrogram; Colv="NA")            # turn off column clustering
  Rowv = as.dendrogram(hc),
  #Colv = as.dendrogram(hc2),
  key.title = NA,
  key.xlab = NA,
  key.par = list(mgp = c(1.5, 0.5, 0),
                 mar = c(5.5, 8, 6, 6)),
  # (bottom, left, top, right)
  key.xtickfun = function() {
    cex <- 4 * par("cex") * par("cex.axis")
    side <- 1
    line <- 2
    col <- par("col.axis")
    font <- par("font.axis")
    mtext(
      "0",
      side = side,
      at = 0,
      adj = 0,
      line = line,
      cex = cex,
      col = col,
      font = font
    )
    mtext(
      "1",
      side = side,
      at = 1,
      adj = 1,
      line = line,
      cex = cex,
      col = col,
      font = font
    )
    return(list(labels = FALSE, tick = FALSE))
  }
)
text(
  x = rep(0.2, 8),
  y = c(0.67, 0.45, 0.35, 0.28, 0.18, 0.05),
  labels = c('TNFa\nIFNs', 'Ctrl', 'Anti-pairs', 'IL4', 'IL1', 'LPS'),
  las = 2,
  col = c(
    "gray30",
    "gray90",
    rep("gray30", 1),
    rep("gray90", 1),
    "gray30",
    "gray90"
  ),
  font = 2,
  cex = 2.5,
  xpd = TRUE
)
#Add cluster labels instead of the x labels
text(x = c(0.28,0.35,0.37,0.41,0.51,0.61,0.72,0.87), y = rep(1, 8),
     labels = c('\nNFkB',' \n\nPI3K','\n \n\n\nSTAT3','\n \nSTAT1',' \n \nSTAT6',
                ' \n \nIL6',' \n \nTLR4',' \n \nMAPKs'),
     las = 2, col = "cornflowerblue", font = 2,cex = 2.5, xpd = TRUE)
# text(x = c(0.28,0.35,0.37,0.41,0.51,0.61,0.72,0.87), y = rep(0, 8),
#      labels = c('\n \nNFkB',' \n\nPI3K','\n \n\n\nSTAT3','\n \nSTAT1',' \n \nSTAT6',
#                 ' \n \nIL6',' \n \nTLR4',' \n \nMAPKs'),
#      las = 2, col = "cornflowerblue", font = 2,cex = 2.5, xpd = TRUE)
text(
  x = 0.99,
  y = 0.38,
  labels = 'Stimulus Combinations',
  las = 2,
  srt = -90,
  cex = 3,
  xpd = TRUE
)
text(
  x = 0.5,
  y = -0.15,
  labels = 'Nodes',
  las = 2,
  cex = 3,
  xpd = TRUE
)
dev.off()               # close the PNG device

# Heatmap for single stimulation only
png(paste0("./plots/Screening_single",input,data,k,"_label.png"),    # create PNG for the heat map        
    width = 12*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mydata.single,# cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
          lmat = rbind(c(4,3),c(2,1)), # 1. heatmap 2. row dendro 3.col dendro 4. key
          lwid = c(1,4), lhei = c(1.5,4),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map; 
          margins =c(6,12),     # widens margins around plot
          col= redscale,       # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
          breaks=c(seq(0,1,length.out = 51)),
          # srtCol = 45,colCol = rainbow(8)[as.factor(group.hc2[,2])],
          #labRow = '', 
          labCol = colnames(mydata.single),
          cexRow = 2,#colRow = rainbow(cl)[groups],
          #RowSideColors = rainbow(cl)[groups],    # Measurement 7-10: red
          dendrogram="both",    # only draw a row dendrogram; Colv="NA")            # turn off column clustering
          Rowv = as.dendrogram(hc.single),
          key.title=NA,
          key.xlab=NA,
          key.par = list(mgp=c(1.5, 0.5, 0),
                         mar = c(5.5, 8, 6, 6)), # (bottom, left, top, right)
          key.xtickfun=function() {
            cex <- 4*par("cex")*par("cex.axis")
            side <- 1
            line <- 2
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("0", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("1", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          }
)
#Add cluster labels instead of the x labels
# if (input==0.7) {
#   text(x = c(0.25,0.36,0.43,0.49,0.59,0.7,0.76), y = rep(1, 8),
#        labels = c(' \nIL4-STAT6',' \n \nTLR4',' \n \nSTAT4',' \n \n \nIL6',
#                   ' \n \n \nSTAT1 NFkB ',' \n \n \nPI3K',' \n \nMAPKs'),
#        las = 2, col = "cornflowerblue", font = 2,cex = 2.5, xpd = TRUE)
# }
text(x = 0.99, y = 0.38, labels = 'Stimuli',las = 2, srt = -90, cex = 3, xpd = TRUE)
text(x = 0.5, y = -0.1, labels = 'Nodes',las = 2, cex = 3, xpd = TRUE)
dev.off()

# PCA
mds <- cmdscale(dist(mydata), k=3, eig=TRUE) 
rownames(mds$points) <- rownames(mydata)
mds.single <- cmdscale(dist(mydata.single), k=3, eig=TRUE) 
rownames(mds.single$points) <- rownames(mydata.single)

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / max(mds$eig)
# plot the PCA
plot(eig_pc[1:10],
     type="h", lwd=15, las=1,
     xlab="Dimensions", 
     ylab="Proportion of explained variance", y.axis=NULL,
     col="darkgrey")


require(ggplot2)
require(ggrepel)
groups <- cutree(hc, k=cl) # cut tree into 7 clusters
# set.seed(42)
pts <- as.data.frame(mds$points[,1:3])
pts$V2 <- -pts$V2
# colnames(pts) <- c("PC1 (57%)","PC2 (14%)","PC3 (10%)")
# gp <- factor(rep((1:9),each = 9))
# gp <- factor(rep(c("M0","M1","M2"), each = 4))

  pclabs <- rbind(c("PC1 (55%)","PC2 (17%)","PC3 (9%)"),c("PC1 (53%)","PC2 (17%)","PC3 (10%)"),
                  c("PC1 (52%)","PC2 (17%)","PC3 (10%)"),c("PC1 (52%)","PC2 (17%)","PC3 (10%)")
                  # c("PC1 (52%)","PC2 (17%)","PC3 (10%)"),c("PC1 (52%)","PC2 (17%)","PC3 (10%)")
  )
  rownames(pclabs) <- as.character(timepoint)
  lab.idx <- k
  groups.single <- cutree(hc.single,k=4)
  pts.single <- as.data.frame(mds.single$points[,1:3])
 
   label.dim <- rownames(pts)
    labellist <- c("LPS+","LPS+IFNg","LPS+TNFa","LPS+IL4","IFNg+","IFNg+IL4",
                   "IFNb+","IFNb+IL4","TNFa+","TNFa+IL4","IL1+","IL4+IFNb","IL4+","IL4+IL6",
                   "IL4+IL12","IL6+","IL6+IL12","IL10+","IL10+IL12","IL12+","control")
    
  label.dim[label.dim %not in% labellist] <- ""
  label.dim[label.dim %in% single.names] <- c("LPS","IFNg","IFNb","TNFa","IL1","IL4","IL6",
                                              "IL10","IL12","Control")


  eig_pc_true <- mds$eig * 100 / sum(mds$eig)
  eig_pc_true[1] #55%
  eig_pc_true[2] #17%
#1037 PCA distribution
pca2d.d <- ggplot(pts,aes(V1, V2,label=label.dim)) +
  geom_point(aes(col = as.factor(groups)),size = 3)+ 
  scale_colour_manual(values=rainbow(cl)) +
  labs(x=pclabs[toString(lab.idx),1],y=pclabs[toString(lab.idx),2]) +
  geom_text_repel(size=4,segment.size = 0.5, box.padding = unit(0.3, "lines")) +
  coord_fixed() +
  scale_x_continuous(limits = c(-4.5,6.5),minor_breaks = 0.001, breaks = seq(-4,6,2)) +
  scale_y_continuous(limits = c(-3,4.1),minor_breaks = 0.001, breaks = seq(-4,4,2)) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = 'black'),
        axis.text = element_text(size = 20, colour = 'black'),
        axis.title.x = element_text(vjust = -1), 
        panel.border = element_rect(colour = 'black', size = 1.5, fill = NA),
        panel.grid.major = element_blank(), #line(colour = "gray40", size = 0.5, linetype = 2),
        panel.grid.minor = element_line(colour = "gray40", size = 0.5, linetype = 2),
        panel.background = element_rect(fill = "gray90"),
        legend.key = element_rect(fill = "gray90"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(2,"lines"),
        # legend.spacing = unit(3,"cm"),
        legend.title = element_blank(),
        legend.position='right'
  ) 
ggsave(pca2d.d, file=paste0("./plots/Screening_PCA2D_dim",input,data,k,".png"), width=6, height=4,dpi = 300)

#1172 PCA contribution
require(FactoMineR)
require(factoextra)

# mydata2.num <- as.matrix(mydata2[,3:132])
res.pca <- PCA(mydata, scale.unit = FALSE, ncp = 3, graph = FALSE)

  eigenvalues1[[0.7]] <- res.pca$eig$`percentage of variance`[1:3] #error?
eigenvalues <- res.pca$eig

res.pca.single <- PCA(mydata.single, scale.unit = FALSE, ncp = 3, graph = FALSE)
eigenvalues.single1[[input1]] <- res.pca.single$eig$`percentage of variance`[1:3] #error?
eigenvalues.single <- res.pca.single$eig

# head(eigenvalues[, 1:2])
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")

plot(res.pca, choix = "var")
# Control variable colors using their contribution
# Possible values for the argument col.var are :
# "cos2", "contrib", "coord", "x", "y"
var.label.id <- as.data.frame(res.pca$var$cos2)
var.label.id$species <- rownames(res.pca$var$contrib)
sort1 <- arrange(var.label.id,desc(`Dim.1`))[1:50,]
sort2 <- arrange(var.label.id,desc(`Dim.2`))[1:50,]
# dim1 <- anti_join(sort1,sort2, by = 'species')[1:30,]
# dim2 <- anti_join(sort2,sort1, by = 'species')[1:25,]
sort <- unique(rbind(sort1,sort2)[,4])[c(1:2,4:5,10,24:26,32:33,36,47,53,56,61:62,71:72,74,
                                         82,85,94:95)]
dim <- c("NFKB","TNFa_mrna","iNOS_mrna","CREB","p38","ERK12","JNK","AP1","IL1_mrna",
         "IL6_mrna","SOCS3","PI3K","AKT","TLR4" ,"GSK3","STAT1","SOCS1","STAT6","STAT3",
         "PPARg_mrna" ,"Arg1_mrna","Fizz1_mrna","VEGF_mrna","IL4Ra_mrna") # added "SOCS1_mrna","SOCS3_mrna","IL4Ra_mrna","IL4R"
# lab <- var.label.id
var.label.id[rownames(var.label.id) %not in% dim,4] <- ''
res.pca.2 <- res.pca
rownames(res.pca.2$var$contrib) <- var.label.id$species
rownames(res.pca.2$var$cos2) <- var.label.id$species
rownames(res.pca.2$var$coord) <- var.label.id$species

var.con <- as.data.frame(res.pca.2$var$coord)
var.con$Dim.1 = -var.con$Dim.1
Var.con <- ggplot(var.con,aes(Dim.1, Dim.2,label=var.label.id$species)) +
  geom_point(size = 4, color = 'purple2', alpha = 0.5)+ 
  # scale_colour_manual(values=rainbow(cl)) +
  # guides(guide_legend(ncol=3)) +
  geom_text_repel(size=6,segment.size = 0.5, force = 0.5, segment.color = 'gray40',
                  segment.alpha = 0.7,
                  box.padding = unit(0.6, "lines"),
                  point.padding = unit(0.2, "lines"), 
                  min.segment.length = unit(0.2, "lines")) +
  coord_fixed() +
  labs(x=pclabs[toString(lab.idx),1],y=pclabs[toString(lab.idx),2]) +
  scale_x_continuous(limits = c(-0.6,0.4),minor_breaks = 0.001, breaks = seq(-0.4,0.4,0.2)) +
  scale_y_continuous(limits = c(-0.2,0.4),minor_breaks = 0.001, breaks = seq(-0.2,0.4,0.2)) +
  theme_bw() +
  theme(text = element_text(size = 28, colour = 'black'),
        axis.text = element_text(size = 28, colour = 'black'),
        axis.title.x = element_text(vjust = -1), 
        panel.border = element_rect(colour = 'black', size = 1.5, fill = NA),
        panel.grid.major = element_blank(), #line(colour = "gray40", size = 0.5, linetype = 2),
        panel.grid.minor = element_line(colour = "gray40", size = 1, linetype = 2),
        panel.background = element_rect(fill = "gray100"),
        legend.key = element_rect(fill = "gray90"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(2,"lines"),
        # legend.spacing = unit(3,"cm"),
        legend.title = element_blank(),
        legend.position='right'
  ) +
  annotate('text',x = c(-0.54,-0.42,0.3,0.3), y = c(0.34,-0.15,0.37,-0.15), 
           label = c(' MAPKs\nSTAT1','NF-kB','STAT6','STAT3'), 
           size = 9, color = c('orange','orange','forestgreen','forestgreen'))
ggsave(Var.con,file=paste0("./plots/Screening_PCA2D_f_contri_dim",input,data,k,"_1moreLabel.png"), width=8, height=6,dpi = 300)


###single input PCA distribution
eig_pc_true_single <- mds.single$eig * 100 / sum(mds.single$eig)
eig_pc_true_single[1] #61%
eig_pc_true_single[2] #14%

groups.single <- cutree(hc.single,k=4)
pts.single <- as.data.frame(mds.single$points[,1:3])
# pts.single$V2 <- -pts.single$V2
# colnames(pts) <- c("PC1 (57%)","PC2 (14%)","PC3 (10%)")
# gp <- factor(rep((1:9),each = 9))
# gp <- factor(rep(c("M0","M1","M2"), each = 4))
pca2d.single <- ggplot(pts.single,aes(V1, V2,label=rownames(pts.single))) +
  geom_point(aes(col = as.factor(groups.single)),size = 3)+
  scale_colour_manual(values=rainbow(4)) +
  labs(x='PC1 (61%)',y='PC2 (14%)') +
  geom_text_repel(size=6,segment.size = 0.5, box.padding = unit(0.5, "lines")) +
  coord_fixed() +
  scale_x_continuous(limits = c(-4.5,6.5),minor_breaks = 0.001, breaks = seq(-4,6,2)) +
  scale_y_continuous(limits = c(-3,4.1),minor_breaks = 0.001, breaks = seq(-4,4,2)) +
  theme_bw() +
  theme(text = element_text(size = 20, colour = 'black'),
        axis.text = element_text(size = 20, colour = 'black'),
        axis.title.x = element_text(vjust = -1), 
        panel.border = element_rect(colour = 'black', size = 1.5, fill = NA),
        panel.grid.major = element_blank(), #line(colour = "gray40", size = 0.5, linetype = 2),
        panel.grid.minor = element_line(colour = "gray40", size = 0.5, linetype = 2),
        panel.background = element_rect(fill = "gray90"),
        legend.key = element_rect(fill = "gray90"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(2,"lines"),
        # legend.spacing = unit(3,"cm"),
        legend.title = element_blank(),
        legend.position='right'
  )
ggsave(pca2d.single, file=paste0("./plots/Screening_PCA2D_single",input,data,k,".png"), width=6, height=4,dpi = 300)


# Single inputs PCA contribution
var.con.s <- as.data.frame(res.pca.single$var$coord)
var.con.s$Dim.1 = -var.con.s$Dim.1
Var.con.s <- ggplot(var.con.s,aes(Dim.1, Dim.2,label=var.label.id$species)) +
  geom_point(size = 4, color = 'purple2', alpha = 0.5)+ 
  # scale_colour_manual(values=rainbow(cl)) +
  # guides(guide_legend(ncol=3)) +
  geom_text_repel(size=6,segment.size = 0.5, force = 0.5, segment.color = 'gray40',
                  segment.alpha = 0.7,
                  box.padding = unit(0.6, "lines"),
                  point.padding = unit(0.2, "lines"), 
                  min.segment.length = unit(0.2, "lines")) +
  coord_fixed() +
  labs(x='PC1 (60%)',y='PC2 (15%)') +
  scale_x_continuous(limits = c(-0.6,0.4),minor_breaks = 0.001, breaks = seq(-0.4,0.4,0.2)) +
  scale_y_continuous(limits = c(-0.25,0.4),minor_breaks = 0.001, breaks = seq(-0.2,0.4,0.2)) +
  theme_bw() +
  theme(text = element_text(size = 28, colour = 'black'),
        axis.text = element_text(size = 28, colour = 'black'),
        axis.title.x = element_text(vjust = -1), 
        panel.border = element_rect(colour = 'black', size = 1.5, fill = NA),
        panel.grid.major = element_blank(), #line(colour = "gray40", size = 0.5, linetype = 2),
        panel.grid.minor = element_line(colour = "gray40", size = 1, linetype = 2),
        panel.background = element_rect(fill = "gray100"),
        legend.key = element_rect(fill = "gray90"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(2,"lines"),
        # legend.spacing = unit(3,"cm"),
        legend.title = element_blank(),
        legend.position='right'
  ) +
  annotate('text',x = c(-0.52,-0.5,0.3,0.3), y = c(0.33,-0.22,0.37,-0.22), 
           label = c(' MAPKs \nSTAT1','NF-kB','STAT6','STAT3'), 
           size = 10, color = c('orange','orange','forestgreen','forestgreen'))
ggsave(Var.con.s,file=paste0("./plots/Screening_PCA2D_f_contri_single",input,data,k,"_1moreLabel.png"), width=8, height=6,dpi = 300)

#Fig 5A signaling modules
mydata.module <- aggregate(. ~ group,mydata.module.1,sum)
rownames(mydata.module) <- c('TNFa','TLR4','IL6/IL10','IFNg','STAT6','STAT4','TNFamrna','IFNb','STAT1',
                               'MAPKs','IKB','NF-kB','PKA','CREB','RAS','STAT3','PI3K','GSK3','SOCS1mrna')
mydata.module.plot <- as.data.frame(t(mydata.module[c('MAPKs','NF-kB','STAT1','STAT3','PI3K','STAT6'),2:(uniq+1)]))
temp.module <- cbind(row.names(mydata.module.1),mydata.module.1$group)
temp.module <- temp.module[order(mydata.module.1$group),]
fc.module <- factor(as.numeric(temp.module[,2]))
levels(fc.module) <- c('TNFa','TLR4','IL6/IL10','IFNg','STAT6','STAT4','TNFamrna','IFNb','STAT1',
                       'MAPKs','IKB','NF-kB','PKA','CREB','RAS','STAT3','PI3K','GSK3','SOCS1mrna')
temp.module.1 <- cbind(temp.module,as.character(fc.module))
colnames(temp.module.1) <- c("Node","Module#","Module Names")
write_csv(as.data.frame(temp.module.1),paste0("./simulation results/Signaling_Modules_for_Fig5A",input,data,k,".csv"))

d <- as.matrix(mydata.module.plot[c("IFNg+","IFNg+IL4","IL4+","control"),])
png(paste0("./plots/Screening_module_dim2_IFNg",input,'_',k,"_original.png"),    # create PNG for the heat map
    width = 2.5*300,        # 5 x 300 pixels
    height = 2*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(d,# cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
          lmat = rbind(c(2,4),c(3,1)), # 1. heatmap 2. row dendro 3.col dendro 4. key
          lwid = c(1.5,4), lhei = c(1,3),
          scale = 'column', #'none', #'column',
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map;
          margins =c(6,8),     # widens margins around plot
          col= bluered(100),       # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
          breaks= c(seq(-1.5,1.5,length.out = 101)),
          #col = redscale,
          #breaks= c(seq(-1.5,1.5,length.out = 51)), #c(seq(0,max(d),length.out = 51)),
          labRow = c('IFNg','IFNg+IL4','IL4','Control'), cexRow = 1.5,
          adjCol = c(NA,0.5),
          dendrogram="none",    # only draw a row dendrogram;
          Rowv ="NA", Colv = 'NA',           # turn off column clustering
          key.title=NA,
          key.xlab=NA,
          key.par = list(mgp=c(1.5, 0.5, 0),
                         mar = c(3, 1, 2, 15)), # (bottom, left, top, right)
          key.xtickfun=function() {
            cex <- 2*par("cex")*par("cex.axis")
            side <- 1
            line <- 1
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          }
)
dev.off()               # close the PNG device
