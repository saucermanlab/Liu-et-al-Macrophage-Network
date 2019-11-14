# Code orginally from Astor Liu validation.R (11/2/2017)
# Jingyuan Zhang reduced version to only generate Model vs Lindsey RNAseq heatmap (rms normalized) for Figure 2C
# 5/8/2018

setwd('Model Simulation')
require(tidyverse)
require(gdata)
require(gplots)

val.m1 <- read_delim('simulation results/LPS+IFNg validation RNASeq.txt',delim = "\t") %>% 
  select(output,expression,match)

val.m2 <- read_delim('simulation results/IL4 validation RNASeq.txt',delim = "\t") %>% 
  select(output,expression,match)

# Read model simulated M1 and M2 outputs
times <- c(142,182,222,342) #200 # (time-102)/10?
yend <- 142

input <- 0.7
m1 <- read_delim(paste0("simulation results/macmodelvalidation_M1in",input,".txt"),delim = "\t")[,c(1:2,yend)]
# m1all <- read_delim(paste0("../macmodelvalidation_M1in",input,".txt"),delim = "\t")
colnames(m1) <- c('species','yStart','yEnd')
rownames(m1) <- m1$species
m1fc.6.2 <- m1%>% 
  mutate(logFC = log2(yEnd/yStart)) %>% 
  select(species,logFC)

m2 <- read_delim(paste0("simulation results/macmodelvalidation_M2in",input,".txt"),delim = "\t")[,c(1:2,yend)] 
colnames(m2) <- c('species','yStart','yEnd') 
rownames(m2) <- m2$species
m2fc.6.2 <- m2%>% 
  mutate(logFC = log2(yEnd/yStart)) %>% 
  select(species,logFC)

m1fc <- m1fc.6.2
m2fc <- m2fc.6.2

val.plot <- full_join(val.m1,val.m2,by = "output")
val.plot <- left_join(val.plot,m1fc, by=c("output"="species"))
val.plot <- left_join(val.plot,m2fc, by=c("output"="species"))

val.plot.2 <- as.matrix(val.plot[,c(2,4,6,7)])
rownames(val.plot.2) <- val.plot$output
colnames(val.plot.2) <- c("M1-RNASeq","M2-RNASeq","M1-Model","M2-Model")

M1marker <- c("iNOS_mrna","IL1_mrna","IL6_mrna","TNFa_mrna")
M2marker <- c("Arg1_mrna", "Ym1_mrna","Fizz1_mrna","PPARg_mrna","Myc_mrna")
Marker <- c(M1marker,M2marker)

val.plot.dim.ps.s <- val.plot.2

val.plot.dim.ps.s[,c(1,2)] <- t(scale(t(val.plot.2[,c(1,2)]),center = FALSE))
val.plot.dim.ps.s[,c(3,4)] <- t(scale(t(val.plot.2[,c(3,4)]),center = FALSE))
x.label <- unlist(strsplit(rownames(val.plot.dim.ps.s),'_'))
x.label <- x.label[-grep('mrna',x.label)]

val.plot.dim.ps.s.2 <- cbind(val.plot.dim.ps.s[,c(3,4)],0,val.plot.dim.ps.s[,c(1,2)])

group <- rownames(val.plot.dim.ps.s.2)
group[-match(Marker,group,nomatch = 0)] <- 0
group[match(M1marker,group)] <- 1
group[match(M2marker,group)] <- 2
group.both.dim.ps <- as.factor(group)

#hc <- hclust(dist(val.plot.dim.ps.s.2[,-3]),method = 'ward.D') #without the zero column

png(paste0("plots/Validation_both_dim_Seq_cscaled",input,yend,"_alt_col.png"),    # create PNG for the heat map        
    width = 4*300,        # 5 x 300 pixels
    height = 1.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 9)        # smaller font size
heatmap.2(t(val.plot.dim.ps.s.2),# cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
          lmat = rbind(c(0,4),c(0,3),c(2,1)), # 1. heatmap 2. row dendro 3.col dendro 4. key
          lwid = c(0.5,4), lhei = c(1.5,1.5,4),
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map; 
          margins =c(4.5,10),     # widens margins around plot (bottom, right)
          # col=c(colorpanel(40, low = "blue", high = "white"),colorpanel(60, low = "white", mid = "orange",high = "red")),
          col=c(colorpanel(50, low = "blue", high = "white"),colorpanel(50, low = "white", mid = "#ffffbf",high = "red")),
          #col=bluered(100),       # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
          breaks=c(seq(-1,-0.01,length.out =50),seq(0,1,length.out = 51)),
          labCol = x.label, labRow = c('IFNg+LPS','IL4',' ','IFNg+LPS','IL4'),
          colCol = c("black","orange","forestgreen")[group.both.dim.ps],
          srtCol = 90,
          cexRow = 1,
          # colRow = rainbow(8)[groups],
          # RowSideColors = c("magenta","cyan")[groups.both.ps],    # Measurement 7-10: red
          dendrogram="column",    # only draw a row dendrogram;
          Rowv="NA", #Colv=as.dendrogram(hc),           # turn off column clustering
          key.title=NA,
          key.xlab=NA,
          key.par = list(mgp=c(1.5, 0.5, 0),
                         mar = c(0, 20, 2, 5)), # (bottom, left, top, right)
          key.xtickfun=function() {
            cex <- 1.5*par("cex")*par("cex.axis")
            side <- 1
            line <- 1
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("-1", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("1", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          }
)

dev.off()               # close the PNG device
