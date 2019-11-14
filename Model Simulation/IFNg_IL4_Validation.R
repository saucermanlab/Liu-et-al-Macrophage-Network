# Code orginally from Astor Liu validationelse.R (11/2/2017)
# Jingyuan Zhang reduced version to only generate Model vs RNAseq heatmap (rms normalized) and actiivity plots for Figure 5B and 5C
# 5/3/2018

require(tidyverse)
require(gdata)
require(reshape2)
require(gplots)

setwd("Model Simulation")

val.names <- list.files('simulation results',pattern='*_validation.txt$')
val.actnames <- list.files('simulation results',pattern='*_act.txt$')

gse1 <- 'GSE84520'

val.gse.names.1 <- list.files('simulation results',pattern='*.txt$')
val.gse.names.1 <- val.gse.names.1[-grep('_raw',val.gse.names.1)]
val.gse.names.1 <- val.gse.names.1[-grep('_act',val.gse.names.1)]
val.gse.names <- lapply(gse1,function(x,y) y[grep(x,y)],val.gse.names.1)
names(val.gse.names) <- paste0('GSE',gse1)

val.gse <- list()
i = 1
val.gse1 <- lapply(unlist(val.gse.names[i]),function(x) read_delim(x,delim = '\t')[c(4,9)])
names(val.gse1) <- unlist(val.gse.names[i])
val.gse[i] <- list(val.gse1)

names(val.gse) <- gse1

M1marker <- c("iNOS_mrna","IL1_mrna","IL6_mrna","TNFa_mrna" #, "CCL5_mrna", "IKBa_mrna"
)
M2marker <- c("Arg1_mrna", "Ym1_mrna","Fizz1_mrna","PPARg_mrna","Myc_mrna"
)
Marker <- c(M1marker,M2marker)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# Validation
Val <- list()
lim <- c(1.3,1.3,1.3,1.3)
tend <- c(142,182,342) # End time point +4h, +8h, +24h

ti = tend[1]
i = 1
    val <- val.gse[[i]]
    
    for (j in 1:length(val)) {
      val2 <- val[[j]]
      vname <- unlist(strsplit(names(val)[j],'validation'))[1]
      val2.model.1 <- read_delim(val.actnames[grep(pattern = vname,val.actnames,fixed = TRUE)[1]],
                                 delim = "\t")[,1:350]
      
      actout <- c('SOCS1_mrna','IL4Ra_mrna','Arg1_mrna','iNOS_mrna','TNFa_mrna')
      if (j==1) {act.plot <- cbind(vname,val2.model.1[val2.model.1$yAll1 %in% actout,])}
      else {act.plot <- rbind(act.plot,cbind(vname,val2.model.1[val2.model.1$yAll1 %in% actout,]))}
      
      val2.model <- val2.model.1[,c(1,102,ti)]# indice change
      colnames(val2.model) <- c('species','yStart','yEnd')
      rownames(val2.model) <- val2.model$species
      val2.model <- val2.model %>% 
        mutate(logFC=yEnd-yStart) %>% 
        select(species,logFC)
      colnames(val2.model) <- c('species',vname)
      
      if (j==1) {
        val.m1.2 <- val2
        val.model <- val2.model
      } else {
        val.m1.2 <- full_join(val.m1.2, val2, by = "output")
        val.model <- cbind(val.model,val2.model[,2])
      }
    }
    colnames(val.m1.2) <- c('output',names(val))
    val.model.2 <- left_join(val.m1.2,val.model,by=c('output' = 'species'))
    
    # Plot activities
    for (p in actout) {
      act.plot.1 <- t(act.plot[grep(p,act.plot$yAll1),1:300])
      act.plot.1 <- gsub('M0 vs M1','LPS+IFNg',act.plot.1)
      act.plot.1 <- gsub('M0 vs M2','IL4',act.plot.1)
      act.plot.1[1,] <- gsub('\\+','\n\\+',act.plot.1[1,])
      colnames(act.plot.1) <- act.plot.1[1,]
      act.plot.1 <- act.plot.1[-(1:2),]
      rownames(act.plot.1) <- ((1:dim(act.plot.1)[1])-102)/10
      act.plot.1m <- melt(act.plot.1)
      colnames(act.plot.1m) <- c('Time','Stimulation','Activity')
      act.plot.1m$Activity <- as.numeric(levels(act.plot.1m$Activity))[act.plot.1m$Activity]
      gact <- ggplot(data=act.plot.1m,aes(x=Time,y=Activity,group=Stimulation))+
        geom_line(aes(color=Stimulation),size=2)+
        scale_color_manual(values=c("orange", "skyblue3", "forestgreen"))+
        scale_x_continuous(limits = c(-4,14),minor_breaks = 0.001, breaks = seq(0,12,4)) +
        scale_y_continuous(limits = c(0,1),label = c("0", "0.5", "1"), breaks = seq(0,1,0.5)) +
        labs(title=p,x="Time (h)", y = "Norm. Activity")+
        theme_bw()+
        theme(text = element_text(size = 20, colour = 'black'),
              axis.text = element_text(size = 20, colour = 'black'),
              axis.title.x = element_text(vjust = -1), 
              panel.border = element_rect(colour = 'black', size = 1.5, fill = NA),
              panel.grid.major = element_blank(), #line(colour = "gray40", size = 0.5, linetype = 2),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              plot.title = element_text(size=20,hjust = 0.5),
              legend.key = element_rect(fill = "white"),
              legend.text = element_text(size = 20),
              legend.key.height = unit(2,"lines"),
              # legend.spacing = unit(3,"cm"),
              legend.title = element_blank(),
              legend.position='right')
      ggsave(gact, file=paste0("../plots/Validation_",gse1[i],"_",ti,p,".png"), width=4.5, height=3,dpi = 300)
    }
    
    
    # Plot the validation results
    len <- length(val.model.2[1,])-1 # number of conditions
    val.m1.plot <- as.matrix(val.model.2[,2:(len+1)])
    # val.m1.plot <- sapply(val.m1.plot[,1:7], as.numeric.factor)
    rownames(val.m1.plot) <- val.model.2$output
    if (len/2>2) {
      val.m1.plot.2 <- val.m1.plot[rowSums(is.na(val.m1.plot))<(len/2-1),]
    } else {
      val.m1.plot.2 <- val.m1.plot
    }
    val.m1.plot.2[is.na(val.m1.plot.2)] <- 0
    colnames(val.m1.plot.2)[grep('M0 vs M2',colnames(val.m1.plot.2))] <- 'IL4'
    colnames(val.m1.plot.2)[grep('M0 vs M1',colnames(val.m1.plot.2))] <- 'LPS+IFNg'
    
    # val.m1.plot.2 <- val.m1.plot.2[hc.m1$order,]
    
    groups.1 <- as.factor(rep(c("Exp","Model"),each=len/2))
    
    # M1 output nodes: M1 marker = 1, M2 marker = 2, else = 0
    group <- rownames(val.m1.plot.2)
    group[-match(Marker,group,nomatch = 0)] <- 0
    group[match(M1marker,group)] <- 1
    group[match(M2marker,group)] <- 2
    group.m1 <- as.factor(group)
    
    # Scale the model outputs
    val.plot.2.p.s <- val.m1.plot.2
    # val.plot.2.p.s[val.plot.2.p.s==0] <- NA
    val.plot.2.p.s[,1:(len/2)] <- t(scale(t(val.plot.2.p.s[,1:(len/2)]),center = FALSE))
    val.plot.2.p.s[,(len/2+1):len] <- t(scale(t(val.plot.2.p.s[,(len/2+1):len]),center = FALSE))
    hc.s <- hclust(dist(val.plot.2.p.s[,1:len]),method = 'ward.D') ###(len/2) vs len
    rownames(val.plot.2.p.s) <- gsub('_mrna',replacement = '',rownames(val.plot.2.p.s))
    val.plot.2.p.s.2 <- cbind(val.plot.2.p.s,0)
    val.plot.2.p.s.2 <- val.plot.2.p.s.2[,c((len/2+1):len,len+1,1:(len/2))]
    png(paste0("../plots/Validation_",gse1[i],"_cscaled_",ti,"_reverse.png"),    # create PNG for the heat map        
        width = 7*300,        # 5 x 300 pixels
        height = 3*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
    heatmap.2(t(val.plot.2.p.s.2),# cellnote = dat,  # same data set for cell labels; main = "Correlation", # heat map title
              na.rm = TRUE,
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map; 
              margins =c(7.5,8),     # widens margins around plot (bottom, right)
              col=c(colorpanel(50, low = "blue", high = "white"),colorpanel(50, low = "white", mid = "#ffffbf",high = "red")),
              # col=c(colorpanel(50, low = "blue", high = "white"),colorpanel(50, low = "white",high = "red")),
              # col=bluered(100),       # use on color palette defined earlier; breaks=col_breaks,    # enable color transition at specified limits
              breaks=c(seq(-lim[i],-0.01,length.out =50),seq(0,lim[i],length.out = 51)),
              colCol = c("black","orange","forestgreen")[group.m1],
              # srtCol = 45,
              cexRow = 2,cexCol = 2,
              labRow = c(paste(colnames(val.m1.plot.2)[(len/2+1):len],''),'',
                         colnames(val.m1.plot.2)[(len/2+1):len]),
              # colRow = rainbow(8)[groups],
              #ColSideColors = c("orange","orange","cyan","forestgreen")[cutree(hc.s,4)],    # Measurement 7-10: red
              dendrogram="column",    # only draw a row dendrogram; 
              Rowv = 'NA', Colv=as.dendrogram(hc.s),            # turn off column clustering
              key.title=NA,
              key.xlab=NA,
              key.par = list(mgp=c(1.5, 0.5, 0),
                             mar = c(2.5, 4, 4, 3)), # (bottom, left, top, right)
              key.xtickfun=function() {
                cex <- 2*par("cex")*par("cex.axis")
                side <- 1
                line <- 1
                col <- par("col.axis")
                font <- par("font.axis")
                mtext("-1", side=side, at=0.11, adj=0,
                      line=line, cex=cex, col=col, font=font)
                mtext("1", side=side, at=0.88, adj=1,
                      line=line, cex=cex, col=col, font=font)
                return(list(labels=FALSE, tick=FALSE))
              }
    )
    dev.off()               # close the PNG device
