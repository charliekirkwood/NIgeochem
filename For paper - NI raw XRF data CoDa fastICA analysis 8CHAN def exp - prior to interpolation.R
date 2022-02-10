if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, compositions, readxl, MASS, psych, GPArotation, reshape2, ica)

# load pointdata csv
#pointdata<-as.data.frame(read_excel("C:/Users/cwk/Documents/RProjects/Tellus NI/Data/Regional_Soils_A_AquaRegia_HDL.xls"))
pointdata<-as.data.frame(read_excel("C:/Users/cwk/Documents/RProjects/Tellus NI/Data/Regional_Soils_A_XRF.xls"))
names(pointdata)[2:3] <- c("Easting", "Northing")
pointdata <- pointdata[,1:55]
pointdata <- na.omit(pointdata)
str(pointdata)
pointdata[, grepl( "O" , names(pointdata))] <- pointdata[, grepl( "O" , names(pointdata))]*10000
colMeans(pointdata[,4:55])
hist(rowSums(pointdata[,4:55]))
names(pointdata)

apply(pointdata[,4:55] , 2, function(c) nrow(pointdata)-sum(c!=0))

#pointdata <- pointdata[ , -which(names(pointdata) %in% c("SO3","Nd","Bi","Mo","Ga","Ta","Sc"))]

pointdata[,4:55] <- apply(pointdata[,4:55], 2, function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2))

#pointdata[pointdata == 0] <- NA
pointdata[is.na(pointdata)]
pointdata <- na.omit(pointdata)
names(pointdata)
pointdata$R <- 1000000-rowSums(pointdata[,4:55])
unique(rowSums(pointdata[,4:56]))

CLRpointdata <- cbind(pointdata[,c("Easting","Northing")], clr(acomp(pointdata[,4:49])))
names(CLRpointdata)
ggplot(data = pointdata, aes(x=Easting, y=Northing, col=Sn)) + geom_point()
ggplot(data = CLRpointdata, aes(x=Easting, y=Northing, col=Sn)) + geom_point()

apply(pointdata, 2, function(c)sum(c=0))

elementsel <- names(pointdata)[-(1:3)]

set.seed(321)

#CUSTOM compositional PCA, with predict to group medians - (not in this script, actually)
pointdata.ilr <- matrix(NA, nrow = nrow(pointdata[,elementsel]), ncol = ncol(pointdata[,elementsel]) - 1)
for (i in 1:ncol(pointdata.ilr)) {
  pointdata.ilr[, i] = sqrt((i)/(i + 1)) * log(((apply(as.matrix(pointdata[,elementsel][,1:i]), 1, prod))^(1/i))/(pointdata[,elementsel][, i + 1]))
}
dataPCA <- icafast(pointdata.ilr, nc=8, alg="def", fun="exp")
pairs(dataPCA$S, col=rgb(0,0,0,5,maxColorValue=255), upper.panel = NULL, labels=c(paste("Comp.", seq(1,8,1))), cex.labels = 1.2)
dataPCA$vafs*100
sum(dataPCA$vafs*100)


#dataPCA <- icaimax(pointdata.ilr, nc=3, fun="ext")
#dataPCA$eigenvalues <- eigen(cov.rob(pointdata.ilr, method= "mve", cor = TRUE)$cov)$values
str(dataPCA)
summary(dataPCA)

#vss(pointdata.ilr)
#dataPCA <- principal(pointdata.ilr, nfactors=8, rotate="promax")
str(dataPCA)

dataPCAloadingsconv <- matrix(0, nrow = ncol(pointdata[,elementsel]), ncol = ncol(pointdata[,elementsel]) - 1)
for (i in 1:ncol(dataPCAloadingsconv)) {
  dataPCAloadingsconv[1:i, i] <- 1/i
  dataPCAloadingsconv[i + 1, i] <- (-1)
  dataPCAloadingsconv[, i] <- dataPCAloadingsconv[, i] * sqrt(i/(i + 1))
}

str(dataPCAloadingsconv)
dataPCA$M <- dataPCAloadingsconv %*% dataPCA$M
dimnames(dataPCA$M)[[1]] <- names(pointdata[,elementsel])

### 3D TRIPLOT
library(rgl)
# PCA COLOURS
PCAscores <- as.data.frame(dataPCA$S)
startcol <- 2
PCAscores$PCAcol <- rgb((PCAscores[,startcol]-min(PCAscores[,startcol]))/(max(PCAscores[,startcol])-min(PCAscores[,startcol])), 
                        (PCAscores[,startcol+1]-min(PCAscores[,startcol+1]))/(max(PCAscores[,startcol+1])-min(PCAscores[,startcol+1])), 
                        (PCAscores[,startcol+2]-min(PCAscores[,startcol+2]))/(max(PCAscores[,startcol+2])-min(PCAscores[,startcol+2])))

plot3d(PCAscores, alpha=0.3, col=PCAscores$PCAcol, box=FALSE, top=FALSE, size=2)
#points3d(dataPCAmedians, alpha=1, col=dataMED$PCAcol, box=FALSE, top=FALSE, size=10)
text3d(dataPCA$M[,1:3]*2, texts=dimnames(dataPCA$M)[[1]], col="black", alpha=0.8, cex=1, font=1)
raycoords <- NULL
for (i in 1:nrow(dataPCA$M)) {
  raycoords <- rbind(raycoords, rbind(c(0,0,0),(dataPCA$M[i,1:3]*2)))
}
lines3d(raycoords, col="black", lwd=1, alpha=0.25)
rgl.bg(color = "white")

ggplot(pointdata, aes(x=Easting, y=Northing)) + geom_point(col=PCAscores$PCAcol, size=1) + theme_bw() +
  coord_equal() + 
  theme_bw(base_size = 6) + labs(x="Easting", y="Northing") +
  theme(legend.position=c(0,1), legend.justification=c(0,1)) +
  #scale_x_continuous(breaks=c(400000,450000,500000)) +
  #scale_y_continuous(breaks=c(400000,450000,500000)) + 
  #scale_fill_gradientn (colours=cubeHelix(1000, start=0.6, r= -1)[150:950],na.value="red") +
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(panel.grid.minor = element_line(colour="gray90", size=0.5)) +
  theme(panel.grid.major = element_line(colour="gray90", size=0.5))
#theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), axis.text.y = element_text(angle=90, hjust=0.5))
ggsave(filename=paste(getwd(), "/","NI 8CHAN DEF EXP ICA 1_2_3 point map.png", sep= ""), height=120, width=180, dpi=600, units="mm")

#Variance plot
#plot(dataPCA$vafs*100, xlab="Principal Component", ylab="Proportion of Variance (%)")
pcaVAR <- data.frame(PC = 1:ncol(dataPCA$M[-1,]), propVariance = dataPCA$vafs*100)
sum(dataPCA$vafs*100)
dataPCA$vafs*100
ggplot(data = pcaVAR, aes(x= PC, y=propVariance)) + geom_point() + scale_y_continuous(limits=c(0,max(dataPCA$vafs*1.01)*100)) + 
  xlab("Component") + ylab("Proportion of Variance (%)") + theme_bw()
ggsave("XRF 8CHAN DEF EXP ICA proportion of variance plot.png", width = 170, height = 100, units="mm", dpi=300)

ggpairs(dataPCA$S)
summary(dataPCA)

# plot(dataPCA$M["Sn",])
# 
# str(dataPCA$M)
# sort(dataPCA$M[,1])
# 
# for(i in 1:20){
# LOAD <- (melt(as.data.frame(t(sort(dataPCA$M[,i], decreasing = TRUE)))))
# ggplot(data = LOAD, aes(x=variable, y=value, group = 1)) + geom_text(aes(label = variable), size=2) + 
#   xlab("Element") + ylab("Loading") + theme(axis.text.x = element_text(angle = 90, size=10, vjust=0.5)) 
# ggsave(paste("NONROBUST PC", i, "loadings.png"), width = 170, height = 100, units="mm", dpi=300)
# }

###BIPLOT
#PCAscores <- as.data.frame(dataPCA$S)
ggplot(data = as.data.frame(dataPCA$S), aes(x=V2, y=V1)) + xlab("Comp. 2") + ylab("Comp. 1") +
  geom_segment(data = as.data.frame(dataPCA$M*3), aes(x=0, y=0, xend=V2, yend=V1), arrow=arrow(length=unit(0.1,"cm")), size=0.1, alpha=0.4) +
  geom_point(alpha=0.1, shape=16, size=1, col=PCAscores$PCAcol) + theme_minimal(8) + scale_shape(solid = TRUE) + 
  geom_text(data = as.data.frame(dataPCA$M*3), aes(x=V2, y=V1, label=dimnames(dataPCA$M*3)[[1]]), size=2, alpha=0.5)
  #xlim(c(-8,8)) + ylim(c(-8,8))
ggsave("XRF 8CHAN DEF EXP ICA biplot 1_2.png", width = 85, height = 85, units="mm", dpi=300)

    ggplot(data = as.data.frame(dataPCA$S), aes(x=V2, y=V1)) + xlab("Comp. 2") + ylab("Comp. 1") + 
      geom_segment(data = as.data.frame(dataPCA$M*3), aes(x=0, y=0, xend=V2, yend=V1), arrow=arrow(length=unit(0.1,"cm")), size=0.1, alpha=0.4) +
      geom_point(alpha=0.1, shape=16, size=1) + theme_minimal(8) + scale_shape(solid = TRUE) + 
      geom_text(data = as.data.frame(dataPCA$M*3), aes(x=V2, y=V1, label=dimnames(dataPCA$M*3)[[1]]), size=2, alpha=0.5)
    #xlim(c(-8,8)) + ylim(c(-8,8))
    ggsave("XRF 8CHAN DEF EXP ICA biplot BW 1_2.png", width = 85, height = 85, units="mm", dpi=300)

ggplot(data = as.data.frame(dataPCA$S), aes(x=V2, y=V3)) + xlab("Comp. 2") + ylab("Comp. 3") + 
  geom_segment(data = as.data.frame(dataPCA$M*3), aes(x=0, y=0, xend=V2, yend=V3), arrow=arrow(length=unit(0.1,"cm")), size=0.1, alpha=0.4) +
  geom_point(alpha=0.1, shape=16, size=1, col=PCAscores$PCAcol) + theme_minimal(8) + scale_shape(solid = TRUE) + 
  geom_text(data = as.data.frame(dataPCA$M*3), aes(x=V2, y=V3, label=dimnames(dataPCA$M*3)[[1]]), size=2, alpha=0.5)
#xlim(c(-8,8)) + ylim(c(-8,8))
ggsave("XRF 8CHAN DEF EXP ICA biplot 2_3.png", width = 85, height = 85, units="mm", dpi=300)

  ggplot(data = as.data.frame(dataPCA$S), aes(x=V2, y=V3)) + xlab("Comp. 2") + ylab("Comp. 3") + 
    geom_segment(data = as.data.frame(dataPCA$M*3), aes(x=0, y=0, xend=V2, yend=V3), arrow=arrow(length=unit(0.1,"cm")), size=0.1, alpha=0.4) +
    geom_point(alpha=0.1, shape=16, size=1) + theme_minimal(8) + scale_shape(solid = TRUE) + 
    geom_text(data = as.data.frame(dataPCA$M*3), aes(x=V2, y=V3, label=dimnames(dataPCA$M*3)[[1]]), size=2, alpha=0.5)
  #xlim(c(-8,8)) + ylim(c(-8,8))
  ggsave("XRF 8CHAN DEF EXP ICA biplot BW 2_3.png", width = 85, height = 85, units="mm", dpi=300)

str(PCAscores)
row.names(as.data.frame(dataPCA$M))

###TERNARY PLOT
library(ggtern)
# PCA COLOURS
PCAscoresscaled <- data.frame(lapply(as.data.frame(dataPCA$S), function(x) ((x-min(x))/(max(x)-min(x)))))
PCAscoresscaled$PCAcol <- rgb(PCAscoresscaled[,1],PCAscoresscaled[,2],PCAscoresscaled[,3])

PCAloadingsscaled <- data.frame(lapply(as.data.frame(dataPCA$M), function(x) ((x-min(x))/(max(x)-min(x)))))
row.names(PCAloadingsscaled) <- row.names(dataPCA$M)

#Comp 1:3
raycoords <- NULL
for (i in 1:nrow(PCAloadingsscaled)) {
  raycoords <- rbind(raycoords, rbind(c(0,0,0),(PCAloadingsscaled[i,1:3])))
}

ggtern() + geom_point(data=PCAscoresscaled, aes(x=1*(V1),y=1*(V2),z=1*(V3)), col=PCAscores$PCAcol, alpha=1, size=1, shape=16) +
  geom_segment(data=as.data.frame(raycoords), mapping=aes(x=0, y=0, z=0, xend=1*(V1),yend=1*(V2),zend=1*(V3)), size=0.2, alpha=0.35) +
  #geom_point(data=as.data.frame(dataPCAmedians), aes(x=(V1),y=(V2),z=(V3),col=dataMED$PCAcol), alpha=1, size=3, shape=17) +
  #geom_point(data=as.data.frame(dataPCAmedians), aes(x=(V1),y=(V2),z=(V3)), alpha=0.75, size=3, shape=24, colour="black") +
  geom_mask() + limits_tern(1,1,1) +
  geom_text(data=as.data.frame(PCAloadingsscaled), mapping=aes(x=1*(V1),y=1*(V2),z=1*(V3), label=rownames(as.data.frame(PCAloadingsscaled))), size=2, alpha=1) +
  scale_colour_identity() + theme_bw() + theme_arrownormal() + theme_clockwise() +
  labs( x       = "", xarrow  = "Comp. 1", y       = "", yarrow  = "Comp. 2", z       = "", zarrow  = "Comp. 3")
ggsave(filename= paste(getwd(),"/","NI XRF 8CHAN DEF EXP ICA Rbor501_1_XY_1_2_3_scaled_nomask.png", sep= ""), height=120, width=180, units="mm", dpi=600)

for( i in 1:ncol(dataPCA$S)){
  pointdata[,paste("PC",i,sep="")] <- dataPCA$S[,i]
}

dataPCA$M
for(i in 1:ncol(dataPCA$M)){
  loadings <- as.data.frame(cbind(dataPCA$M[order(dataPCA$M[,i]^2, decreasing=TRUE), i], names(dataPCA$M[order(dataPCA$M[,i]^2, decreasing=TRUE), i])))
  loadings$V1 <- as.numeric(as.character(loadings$V1))
  loadings$V2 <- factor(loadings$V2, levels=names(dataPCA$M[order(dataPCA$M[,i]^2, decreasing=TRUE), i]))
  str(loadings)
  ggplot(loadings) + geom_bar(aes(x=V2, y=V1), stat="identity", fill="light grey") + labs(x="Element", y="Loading") + theme_bw(6) + 
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) + geom_text(aes(x=V2, y=V1, label=loadings$V2, angle=90, vjust=0.5), size=1.5) +
    scale_y_continuous(labels=function(x) sprintf("%.2f", x), limits=c(-max(abs(loadings$V1)), max(abs(loadings$V1)))) +
    annotate("text", y=max(abs(loadings$V1)), x=nrow(loadings), label=paste("Comp. ", i, sep=""), hjust=1, vjust=1, size=2)
  ggsave(filename= paste(getwd(),"/","NI XRF 8CHAN DEF EXP ICA Loadings comp ", i, ".png", sep= ""), height=80, width=164, units="mm", dpi=600)
}
