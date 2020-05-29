require(tidyverse)
require(gdata)

m1 <- read_csv("simulation results/boolean_M1M2_07.csv")
# colnames(m1) <- c('species','yStart','yEnd')
m1fc <- m1%>% 
  mutate(logFCM1 = log2(yEndM1/yStart)) %>% 
  mutate(logFCM2 = log2(yEndM2/yStart)) %>% 
  select(-c(yEndM1,yEndM2,yStart))

fname <- "bool"

  M1marker <- c("iNOS_mrna","IL1_mrna","IL6_mrna","TNFa_mrna" #, "CCL5_mrna", "IKBa_mrna"
  )
  M2marker <- c("Arg1_mrna", "Ym1_mrna","Fizz1_mrna","PPARg_mrna","Myc_mrna"
  )
  Marker <- c(M1marker,M2marker)
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  val.plot <- m1fc[is.na(m1fc$logFCrnaseqM1)+is.na(m1fc$LogFCrnaseqM2)!=2,]
  val.plot.2 <- as.matrix(val.plot[,c(4:5,2:3)])
  # val.m2.plot <- sapply(val.m2.plot[,1:5], as.numeric.factor)
  rownames(val.plot.2) <- val.plot$Species
  
  hc.both <- hclust(dist(val.plot.2))
  val.plot.2 <- val.plot.2[hc.both$order,]
  
  # Output nodes: M1 marker = 1, M2 marker = 2, else = 0
  group <- rownames(val.plot.2)
  group[-match(Marker,group,nomatch = 0)] <- 0
  group[match(M1marker,group)] <- 1
  group[match(M2marker,group)] <- 2
  group.both <- as.factor(group)
  
  require(gplots)
  
  # Scale the model outputs
  val.plot.2.ps.s <- val.plot.2
  val.plot.2.ps.s[,c(1:2)] <- t(scale(t(val.plot.2.ps.s[,c(1:2)]),center = FALSE))
  val.plot.2.ps.s[,c(3:4)] <- t(scale(t(val.plot.2.ps.s[,c(3:4)]),center = FALSE))
  rownames(val.plot.2.ps.s) <- sub('Inos','Nos2',tools::toTitleCase(tolower(gsub('_mrna',replacement = '',sub('IL1_mrna','Il1a',rownames(val.plot.2.ps.s))))))
  val.plot.2.ps.s.2 <- cbind(val.plot.2.ps.s[,c(1:2)],0,val.plot.2.ps.s[,c(3:4)])
  colnames(val.plot.2.ps.s.2) <- c('LPS+INFg','IL4','','LPS+INFg','IL4')

  pdf(paste0("plots/Validation_",fname,"_cscaled",".pdf"),    # create PNG for the heat map        
      width = 7,        # 5 x 300 pixels
      height = 3,
      # res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  heatmap.2(t(val.plot.2.ps.s.2),# cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
            na.rm = TRUE,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map; 
            margins =c(8,10),     # widens margins around plot (bottom, right)
            # col=c(colorpanel(40, low = "blue", high = "white"),colorpanel(60, low = "white", mid = "orange",high = "red")),
            col=bluered(100),       # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
            breaks=c(seq(-1,-0.01,length.out =50),seq(0,1,length.out = 51)),
            colCol = c("black","orange","forestgreen")[group.both],
            # srtCol = 45,
            cexRow = 2,cexCol = 2,
            # colRow = rainbow(8)[groups],
            # RowSideColors = c("magenta","cyan")[groups.both.ps],    # Measurement 7-10: red
            dendrogram="column",    # only draw a row dendrogram;
            Rowv="NA", #Colv = "NA",           # turn off column clustering
            key.title=NA,
            key.xlab=NA,
            key.par = list(mgp=c(1.5, 1.5, 0),
                           mar = c(5, 4, 1, 2.5)), # (bottom, left, top, right)
            key.xtickfun=function() {
              cex <- 2*par("cex")*par("cex.axis")
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
  text(x = c(0.2,0.2), y = c(0.68,0.25), labels = c('Model','RNA-Seq'),
            las = 2, col = c("black","black"), font = 12, cex = 1.5, xpd = TRUE,srt=90)
  dev.off()               # close the PNG device
  

