rm(list=ls())
library(drc)
library(scatterplot3d)
graphics.off()

plotfigs<-1
# print("not plotting figures")

inputfile <- 'bactarch.csv'
writefile <- inputfile
#VERY IMPORTANT THAT "genes" IS THE TITLE OF THE FIRST COLUMN AFTER THE RBHs
rbhmatrix <- read.csv(writefile, header=T,row.names=1,check.names=T,blank.lines.skip=T)

rbhmatrix$temperature <- factor(rbhmatrix$temperature,levels=c("Hyperthermophilic","Thermophilic","Mesophilic"))
rbhmatrix$oxygen <- factor(rbhmatrix$oxygen,levels=c("Anaerobic","Facultative","Microaerophilic","Aerobic","unknown"))


attach(rbhmatrix)
#initialize
rbhgc <- as.numeric(rbhmatrix$gc.pct)
rbhgenes <- rbhmatrix$genes
rbhcols <- which(colnames(rbhmatrix)=="genes") -1
#rbhmatrix$genes <- NULL
allcols <- ncol(rbhmatrix) #subtract genes column
rbhrows <- nrow(rbhmatrix)
rbhpreds<-as.data.frame(1:rbhrows)
cf<-as.data.frame(1:4)
rbhresid<-as.data.frame(1:rbhrows)
rbhstudresid<-as.data.frame(1:rbhrows)
rbhmax<-as.data.frame(1:rbhrows)
rbhmin<-as.data.frame(1:rbhrows)
rbhov <-as.data.frame(1:rbhrows)
rbhund <-as.data.frame(1:rbhrows)
rbhmid <-as.data.frame(1:rbhrows)
rbhlin <-as.data.frame(1:rbhrows)
rbhstats <- NULL
rbhdifs<-as.data.frame(1:rbhrows)
rbhparams <- NULL

for (i in 1:rbhcols) {

print(paste(writefile," ",colnames(rbhmatrix)[i]))
x<- rbhmatrix[,i]
rbhmodel <- drm(x ~ rbhgenes,fct=LL.4()) 
cf[,i] <- coef(rbhmodel)
Cb<-coef(rbhmodel)[1]
Cc<-coef(rbhmodel)[2]
Cd<-coef(rbhmodel)[3]
Ce<-coef(rbhmodel)[4]
rbhfun <- function(x) Cc + (Cd-Cc)/(1+exp(Cb*(log(x)-log(Ce))))
rbhpreds[,i] <- rbhfun(rbhgenes)

rbhparams <- rbind(rbhparams,c(Cb,Cc,Cd,Ce))

# trying and failing to dual-fit archaea and bacteria genome sizes
#rbhnum <- unlist(rbhmatrix[1:rbhrows,1:rbhcols])
#archgenes <- as.vector(t(matrix(ncol=rbhrows,rep(archinfo$size..Mb.*1000,rbhrows))))
#archcols <- as.vector(t(matrix(ncol=rbhrows,rep(rainbow(rbhcols),rbhrows))))
#
#bactgenes <- rep(rbhgenes,rbhcols)
#rbhdata <- as.data.frame(cbind(archgenes,bactgenes,rbhnum))
#scatterplot3d(x=rbhdata,pch='.')
#library(rgl)
#plot3d(rbhdata,col=archcols)
#
#parm <- list(Cc=200,Cd=250,Ce=-0.1,Cf=0)
#newf <- function(xx,yy) parm$Cc + (parm$Cd - parm$Cc)/(1 + exp((parm$Ce)*(log(xx)+log(yy+parm$Cf))))
#m2 <- nls(rbhnum ~ Cc + (Cd - Cc)/(1 + exp(Ce*(log(bactgenes)+log(archgenes)))),data=rbhdata,start=c(Cc=200,Cd=250,Ce=-0.1,Cf=100))

w<-residuals(rbhmodel)

rbhlm <- lm(w~rbhgenes+rbhgc)
gcstring <- paste("GC content, slope = ",signif(summary(rbhlm)$coefficients[3,1],3),"p = ",signif(summary(rbhlm)$coefficients[3,4],3))
rbhstring <- paste("first residual RBH, intercept = ",signif(summary(rbhlm)$coefficients[1,1],3),"p = ",signif(summary(rbhlm)$coefficients[1,4],3))
genestring <- paste("genes, slope = ",signif(summary(rbhlm)$coefficients[2,1],3),"p = ",signif(summary(rbhlm)$coefficients[2,4],3))

rbhstudresid[,i] <- studres(rbhlm)
rbhresid[,i] <- residuals(rbhlm)


#calculating relative effect of removing genome size
sr <- sort(rbhpreds[,i],index.return=T)
for (j in 1:(rbhrows-1)) {
rbhdifs[j,i] <- (sr$x[j+1]-sr$x[j])/(rbhgenes[sr$ix[j+1]]-rbhgenes[sr$ix[j]])*1000
}
rbhdifs[j+1,i] <- NaN


rbhstats <- rbind(rbhstats,c(
summary(rbhlm)$coefficients[1,1],
summary(rbhlm)$coefficients[1,4],
summary(rbhlm)$coefficients[2,1],
summary(rbhlm)$coefficients[2,4],
summary(rbhlm)$coefficients[3,1],
summary(rbhlm)$coefficients[3,4]))

rbhlmpred <- predict.lm(rbhlm,interval="prediction",type="response",level=0.95)
rbhmax[,i] <- rbhlmpred[,3]
rbhmin[,i] <- rbhlmpred[,2]
rbhlin[,i] <- rbhlmpred[,1]

rbhov[,i]  <- rbhpreds[,i] + rbhmax[,i]
rbhund[,i] <- rbhpreds[,i] + rbhmin[,i]
rbhmid[,i] <- rbhpreds[,i] + rbhlin[,i]

plot(x~rbhgenes)
points(rbhpreds[,i]~rbhgenes,pch='x')
points(rbhgenes,rbhov[,i],pch='.',col='red')
points(rbhgenes,rbhund[,i],pch='.',col='red')


if (plotfigs==1) {

png(filename=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-3d.png",sep=""),width = 640, height = 640, units="px",pointsize=12,bg="white") 
s3d <- scatterplot3d(rbhgenes,rbhgc,w,main=colnames(rbhmatrix[i]),xlab=genestring,ylab=gcstring,zlab=rbhstring)
s3d$plane(rbhlm,lty.box = "solid")
dev.off()


rbhlm1 <- lm(w~rbhgenes)
rbhlm1pred <- predict.lm(rbhlm1,interval="prediction",type="response",level=0.95)

ax <- sort(rbhpreds[,i]+rbhlm1pred[,3],index.return=T)
ay <- sort(rbhpreds[,i]+rbhlm1pred[,1],index.return=T)
az <- sort(rbhpreds[,i]+rbhlm1pred[,2],index.return=T)

ovs <- rbhresid[,i]/sd(rbhresid[,i])>2



svg(file=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-CI-tempox.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white") 

colmatch <- list("Hyperthermophilic"="red", "Thermophilic"="orange","Mesophilic"="green","Psychrophilic"="blue","NA"=NA)
y <- as.character(rbhmatrix$temperature)
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"

pchmatch <- list("Anaerobic"=21, "Facultative"=22,"Microaerophilic"=23,"Aerobic"=24,"NA"=21)
z<-as.character(rbhmatrix$oxygen)
bactpchlist <- as.numeric(as.character(pchmatch[z]))
bactpchlist[is.na(bactpchlist)] <- 46

cexmatch <- list("Hyperthermophilic"=1, "Thermophilic"=1,"Mesophilic"=.4,"Psychrophilic"=1,"NA"=0)
bactcexlist <- as.numeric(as.character(cexmatch[y]))
bactcexlist[is.na(bactcexlist)] <- 0.4

plot(rbhmatrix[,i]~rbhgenes,main=colnames(rbhmatrix[i]),ylim=c(0,500),xlab="genes in bacterial genome",ylab="reciprocal best BLAST hits",pch=NA)
points(rbhmatrix[,i]~rbhgenes,col=bactcollist,bg=bactcollist,pch=bactpchlist,cex=bactcexlist)

lines(ax$x~rbhgenes[ax$ix],col="grey")
lines(ay$x~rbhgenes[ay$ix],col="black")
lines(az$x~rbhgenes[az$ix],col="grey")

#legend(x=6000,y=175,legend=c("Psychrophilic","Mesophilic","Thermophilic","Hyperthermophiic"),pt.bg=c("blue","green","orange","red"),pch=21)
#legend(x=6000,y=100,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pch=c(24,23,22,21),pt.bg="black")

#text(x=rbhgenes[ovs],y=rbhmatrix[ovs,i],labels=rownames(rbhmatrix)[ovs],srt=90,cex=0.4,adj=c(0,0))

dev.off()





svg(file=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-CI-ox.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white") 

colmatch <- list("Anaerobic"="red", "Facultative"="green","Microaerophilic"="violet","Aerobic"="blue","NA"=NA)
y <- as.character(rbhmatrix$oxygen)
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"

pchmatch <- list("Anaerobic"=21, "Facultative"=22,"Microaerophilic"=23,"Aerobic"=24,"NA"=21)
z<-as.character(rbhmatrix$oxygen)
bactpchlist <- as.numeric(as.character(pchmatch[z]))
bactpchlist[is.na(bactpchlist)] <- 46

cexmatch <- list("Anaerobic"=1, "Facultative"=1,"Microaerophilic"=.4,"Aerobic"=0.4,"NA"=0)
bactcexlist <- as.numeric(as.character(cexmatch[y]))
bactcexlist[is.na(bactcexlist)] <- 0.4

plot(rbhmatrix[,i]~rbhgenes,main=colnames(rbhmatrix[i]),ylim=c(0,500),xlab="genes in bacterial genome",ylab="reciprocal best BLAST hits",pch=NA)
points(rbhmatrix[,i]~rbhgenes,col=bactcollist,bg=bactcollist,pch=bactpchlist,cex=bactcexlist)

lines(ax$x~rbhgenes[ax$ix],col="grey")
lines(ay$x~rbhgenes[ay$ix],col="black")
lines(az$x~rbhgenes[az$ix],col="grey")

#legend(x=6000,y=175,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pt.bg=c("blue","violet","green","red"),pch=21)
#legend(x=6000,y=100,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pch=c(24,23,22,21),pt.bg="black")

#text(x=rbhgenes[ovs],y=rbhmatrix[ovs,i],labels=rownames(rbhmatrix)[ovs],srt=90,cex=0.4,adj=c(0,0))

dev.off()




svg(file=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-CI-halox.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white") 

colmatch <- list("Extreme halophilic"="red", "Moderate"="orange","Non-halophilic"="green","NA"=NA)
y <- as.character(rbhmatrix$salinity)
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"

pchmatch <- list("Anaerobic"=21, "Facultative"=22,"Microaerophilic"=23,"Aerobic"=24,"NA"=21)
z<-as.character(rbhmatrix$oxygen)
bactpchlist <- as.numeric(as.character(pchmatch[z]))
bactpchlist[is.na(bactpchlist)] <- 46

cexmatch <-  list("Extreme halophilic"=1, "Moderate"=1,"Non-halophilic"=0.4,"NA"=0)
bactcexlist <- as.numeric(as.character(cexmatch[y]))
bactcexlist[is.na(bactcexlist)] <- 0.4

plot(rbhmatrix[,i]~rbhgenes,main=colnames(rbhmatrix[i]),ylim=c(0,500),xlab="genes in bacterial genome",ylab="reciprocal best BLAST hits",pch=NA)
points(rbhmatrix[,i]~rbhgenes,col=bactcollist,bg=bactcollist,pch=bactpchlist,cex=bactcexlist)

lines(ax$x~rbhgenes[ax$ix],col="grey")
lines(ay$x~rbhgenes[ay$ix],col="black")
lines(az$x~rbhgenes[az$ix],col="grey")

#legend(x=6000,y=175,legend=c("Non-halophilic","Moderate","Extreme"),pt.bg=c("green","orange","red"),pch=21)
#legend(x=6000,y=100,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pch=c(24,23,22,21),pt.bg="black")

#text(x=rbhgenes[ovs],y=rbhmatrix[ovs,i],labels=rownames(rbhmatrix)[ovs],srt=90,cex=0.4,adj=c(0,0))

dev.off()






svg(file=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-CI-metox.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white") 

colmatch <- list("chemoautotroph"="red", "heterotroph"="blue","photosynthesis"="green","pathogen"="violet","symbiont"="violet")
y <- as.character(rbhmatrix$metabolism)
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"

pchmatch <- list("Anaerobic"=21, "Facultative"=22,"Microaerophilic"=23,"Aerobic"=24,"NA"=21)
z<-as.character(rbhmatrix$oxygen)
bactpchlist <- as.numeric(as.character(pchmatch[z]))
bactpchlist[is.na(bactpchlist)] <- 46

cexmatch <- list("chemoautotroph"=1, "heterotroph"=.4,"photosynthesis"=.4,"pathogen"=.4,"symbiont"=.4)
bactcexlist <- as.numeric(as.character(cexmatch[y]))
bactcexlist[is.na(bactcexlist)] <- 0.4


plot(rbhmatrix[,i]~rbhgenes,main=colnames(rbhmatrix[i]),ylim=c(0,500),xlab="genes in bacterial genome",ylab="reciprocal best BLAST hits",pch=NA)
points(rbhmatrix[,i]~rbhgenes,col=bactcollist,bg=bactcollist,pch=bactpchlist,cex=bactcexlist)

lines(ax$x~rbhgenes[ax$ix],col="grey")
lines(ay$x~rbhgenes[ay$ix],col="black")
lines(az$x~rbhgenes[az$ix],col="grey")

#legend(x=6000,y=175,legend=c("Chemoautotrophic","Photosynthetic","Heterotrophic","Pathogenic/Symbiotic"),pt.bg=c("red","green","blue","violet"),pch=21)
#legend(x=6000,y=100,legend=c("Anaerobic","Facultative","Microaerophilic","Aerobic"),pch=c(21,22,23,24),pt.bg="black")

#text(x=rbhgenes[ovs],y=rbhmatrix[ovs,i],labels=rownames(rbhmatrix)[ovs],srt=90,cex=0.4,adj=c(0,0))

dev.off()



png(filename=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-CI-gc.png",sep=""),width = 640, height = 640, units="px",pointsize=12,bg="white") 

plot(rbhmatrix[,i]~rbhgenes,main=colnames(rbhmatrix[i]),xlab="genes in bacterial genome",ylab="reciprocal best BLAST hits")

points(rbhmid[,i]~rbhgenes,pch='x')
points(rbhov[,i]~rbhgenes,pch='+')
points(rbhund[,i]~rbhgenes,pch='-')
dev.off()


png(filename=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-genes.png",sep=""),width = 640, height = 640, units="px",pointsize=12,bg="white")
 
v<-rbhresid[,i]
t<-rbhgenes
h<-hist(v,breaks=25)
par(mar=c(4,4,4,2) + 0.1) # Leave space for z axis
plot(v~t,xlab="genes",ylab="second RBH residuals",main=colnames(rbhmatrix[i]))
par(new=T)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(v),max(v),length=40)
yfit<-dnorm(xfit,mean=mean(v),sd=sd(v))
plot(yhist,xhist,xlim=c(0,max(yhist,yfit)),type="s", axes=F, bty="n", xlab="", ylab="")
lines(yfit,xfit)
dev.off()

png(filename=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-GC.png",sep=""),width = 640, height = 640, units="px",pointsize=12,bg="white")

u<-rbhresid[,i]
t<-rbhgc
h<-hist(u,breaks=25)
par(mar=c(4,4,4,2) + 0.1) # Leave space for z axis
plot(u~t,xlab="GC content",ylab="second RBH residuals",main=colnames(rbhmatrix[i]))
par(new=T)
xhist<-c(min(h$breaks),h$breaks)
yhist<-c(0,h$density,0)
xfit<-seq(min(u),max(u),length=40)
yfit<-dnorm(xfit,mean=mean(u),sd=sd(u))
plot(yhist,xhist,xlim=c(0,max(yhist,yfit)),type="s", axes=F, bty="n", xlab="", ylab="")
lines(yfit,xfit)
dev.off()


png(filename=paste("./figs/",writefile,'-',colnames(rbhmatrix[i]),"-QQ.png",sep=""),width = 640, height = 640, units="px",pointsize=12,bg="white")

#qqnorm(x,main=colnames(rbhmatrix[i]))
#qqline(x)

qqPlot(rbhlm)
#qqq <- qqPlot(rbhlm,reps=10000,id.method="y",id.n=20,line="robust")

dev.off()

} else {}

}

rbhover <- rbhmatrix[1:rbhcols] - rbhov
rbhunder <- rbhund - rbhmatrix[1:rbhcols]
rbhres <- rbhmatrix[1:rbhcols] - rbhmid

colnames(rbhunder) <- c(colnames(rbhmatrix)[1:rbhcols])
colnames(rbhover) <- c(colnames(rbhmatrix)[1:rbhcols])
colnames(rbhres) <- c(colnames(rbhmatrix)[1:rbhcols])
colnames(rbhstats) <- c('intercept','int-p','genes-slope','genes-p','gc-slope','gc-p')

rbhover <- cbind(rbhover,rbhmatrix[(rbhcols+1):allcols])
rbhunder <- cbind(rbhunder,rbhmatrix[(rbhcols+1):allcols])
rbhres <- cbind(rbhres,rbhmatrix[(rbhcols+1):allcols])

rownames(rbhunder) <- c(rownames(rbhmatrix))
rownames(rbhover) <- c(rownames(rbhmatrix))
rownames(rbhres) <- c(rownames(rbhmatrix))
rownames(rbhstats) <- c(colnames(rbhmatrix[1:rbhcols]))

write.csv(rbhunder, file=paste(writefile,"rbhunder.csv",sep="."))
write.csv(rbhover, file=paste(writefile,"rbhover.csv",sep="."))
write.csv(rbhres, file=paste(writefile,"rbhresiduals.csv",sep="."))
write.csv(rbhstats, file=paste(writefile,"rbhstats.csv",sep="."))
write.csv(rbhdifs, file=paste(writefile,"rbhdifs.csv",sep="."))
write.csv(rbhparams, file=paste(writefile,"rbhparams.csv",sep="."))
detach(rbhmatrix)
#quit()


# PLOTTING STUFF

bactnames <- rownames(rbhmatrix)
archnames <- t(read.csv("bactarch.csv",header=F)[1,])[2:58]

archinfo <- read.csv("archaea-info.csv",sep=",",header=T,row.names=1)
colnames(rbhparams) <- c("Cb","Cc","Cd","Ce")
archinfo <- cbind(archinfo,rbhparams)

rM <- rowMeans(rbhresid) #rbhstudresid
rr <- sort(rM,index.return=TRUE)

cM <- colMeans(rbhstudresid)
cr <- sort(cM,index.return=TRUE)

#TO COLOR BACTERIA BY TEMPERATURE
svg(file=paste("./figs/",writefile,"-residual-temperature.svg",sep=""),width = 6, height = 6,pointsize=12,bg="white")

colmatch <- list("Hyperthermophilic"="red", "Thermophilic"="orange","Mesophilic"="green","Psychrophilic"="blue","NA"=NA)
y <- as.character(rbhmatrix$temperature[rr$ix])
bactcollist <- as.character(colmatch[y])
z <- as.character(archinfo$temp)
archcollist <- as.character(colmatch[z])

bactcollist[bactcollist=="NULL"] <- "grey"

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
lines(c(0,447),c(0,0),col="grey")

for (k in 1:rbhcols) {
	points(rbhresid[rr$ix,k],col=rgb(0.7,0.7,0.7,0),bg=rgb(0.7,0.7,0.7,0.5),pch=21,cex=0.5)  #rbhstudresid
}

points(rr$x,pch=21,col=bactcollist,bg=bactcollist,cex=2)

#legend(x=20,y=200,legend=c("Psychrophilic","Mesophilic","Thermophilic","Hyperthermophilic"),pt.bg=c("blue","green","orange","red"),pch=21,cex=2)

dev.off()

#TO COLOR BACTERIA BY OXYGEN
svg(file=paste("./figs/",writefile,"-residual-oxygen.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white")

colmatch <- list("Aerobic"="blue","Microaerophilic"="violet","Facultative"="green","Anaerobic"="red")
y <- as.character(rbhmatrix$oxygen[rr$ix])
bactcollist <- as.character(colmatch[y])
z <- as.character(archinfo$oxygen.req)
archcollist <- as.character(colmatch[z])

bactcollist[bactcollist=="NULL"] <- "grey"

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
lines(c(0,447),c(0,0),col="grey")

for (k in 1:rbhcols) {
	points(rbhresid[rr$ix,k],col=rgb(0.7,0.7,0.7,0),bg=rgb(0.7,0.7,0.7,0.5),pch=21,cex=0.5)  #rbhstudresid
}

points(rr$x,pch=21,col=bactcollist,bg=bactcollist,cex=2)

#colnames(rbhresid) <- colnames(rbhmatrix[,1:57])
#lines(rbhresid[rr$ix,"Pyrococcus.horikoshii.OT3"],col="black")
#lines(rbhresid[rr$ix,"Methanothermobacter.thermautotrophicus.str..Delta.H"],col="black")

#legend(x=20,y=200,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pt.bg=c("green","blue","yellow","red"),pch=21)

dev.off()



#TO COLOR BACTERIA BY METABOLISM
svg(file=paste("./figs/",writefile,"-residual-metabolism.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white")

colmatch <- list("photosynthesis"="green","pathogen"="violet","symbiont"="violet","heterotroph"="blue","chemoautotroph"="red")
y <- as.character(rbhmatrix$metabolism[rr$ix])
bactcollist <- as.character(colmatch[y])
z <- as.character(archinfo$metabolism.req)
archcollist <- as.character(colmatch[z])

bactcollist[bactcollist=="NULL"] <- "grey"

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
lines(c(0,447),c(0,0),col="grey")

for (k in 1:rbhcols) {
	points(rbhresid[rr$ix,k],col=rgb(0.7,0.7,0.7,0),bg=rgb(0.7,0.7,0.7,0.5),pch=21,cex=0.5)  #rbhstudresid
}


points(rr$x,pch=21,col=bactcollist,bg=bactcollist,cex=2)

#legend(x=20,y=200,legend=c("Photosynthetic","Pathogenic/Symbiotic","Heterotrophic","Chemoautotrophic"),pt.bg=c("green","violet","blue","red"),pch=21)

dev.off()



#TO COLOR BACTERIA BY SALINITY
svg(file=paste("./figs/",writefile,"-residual-salinity.svg",sep=""),width = 6, height = 6, pointsize=12,bg="white")

colmatch <- list("Non-halophilic"="green","Moderate"="blue","Extreme halophilic"="red")
y <- as.character(rbhmatrix$salinity[rr$ix])
bactcollist <- as.character(colmatch[y])
z <- as.character(archinfo$metabolism.req)
archcollist <- as.character(colmatch[z])

bactcollist[bactcollist=="NULL"] <- "grey"

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
lines(c(0,447),c(0,0),col="grey")

for (k in 1:rbhcols) {
	points(rbhresid[rr$ix,k],col=rgb(0.7,0.7,0.7,0),bg=rgb(0.7,0.7,0.7,0.5),pch=21,cex=0.5)  #rbhstudresid
}


points(rr$x,pch=21,col=bactcollist,bg=bactcollist,cex=2)


dev.off()


#TO COLOR BACTERIA BY TEMPERATURE, OXYGEN, METABOLISM, SALT
svg(file=paste("./figs/",writefile,"-residual-all.svg",sep=""),width = 6, height = 6,pointsize=12,bg="white")

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
#lines(c(0,447),c(0,0),col="grey")
#for (k in 1:rbhcols) {
#	points(rbhresid[rr$ix,k],col=rgb(0.7,0.7,0.7,0),bg=rgb(0.7,0.7,0.7,0.5),pch=21,cex=0.5)  #rbhstudresid
#}

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
hc <- rev(heat.colors(5,alpha=0.7))
colmatch <- list("Aerobic"=hc[1],"Microaerophilic"=hc[2],"Facultative"=hc[3],"Anaerobic"=hc[5])
y <- as.character(rbhmatrix$oxygen[rr$ix])
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"
points(rr$x,pch=21,col=NA,bg=bactcollist,cex=1.5)
lines(c(0,447),c(0,0),col="grey")
legend(x=20,y=200,legend=c("Aerobic","Microaerophilic","Facultative","Anaerobic"),pt.bg=hc[c(1,2,4,5)],pch=21)

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
hc <- rev(heat.colors(5,alpha=0.5))
colmatch <- list("Hyperthermophilic"=hc[5], "Thermophilic"=hc[4],"Mesophilic"=hc[2],"Psychrophilic"=rgb(t(col2rgb("blue"))/255,alpha=0.7))
y <- as.character(rbhmatrix$temperature[rr$ix])
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"
points(rr$x+0,pch=21,col=NA,bg=bactcollist,cex=1.5)
lines(c(0,447),c(0,0),col="grey")
legend(x=20,y=200,legend=c("Hyperthermophilic","Thermophilic","Mesophilic","Psychrophilic"),pt.bg=c(hc[c(5,4,2)],rgb(t(col2rgb("blue"))/255,alpha=0.7)),pch=21)

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
hc <- rev(heat.colors(5,alpha=0.5))
colmatch <- list("Non-halophilic"=hc[2],"Moderate"=hc[3],"Extreme halophilic"=hc[5])
y <- as.character(rbhmatrix$salinity[rr$ix])
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"
points(rr$x,pch=21,col=NA,bg=bactcollist,cex=1.5)
lines(c(0,447),c(0,0),col="grey")
legend(x=20,y=200,legend=c("Non-halophilic","Moderately halophilic","Extremely halophilic"),pt.bg=hc[c(2,3,5)],pch=21)

plot(rr$x,ylim=c(-200,200),pch='',xlab="Bacterial Genomes",ylab="Residual RBH")
hc <- c(rgb(t(col2rgb("green"))/255,alpha=0.8), rgb(t(col2rgb("blue"))/255,alpha=0.8), rgb(t(col2rgb("yellow"))/255,alpha=0.8), rgb(t(col2rgb("red"))/255,alpha=0.8))
colmatch <- list("photosynthesis"=hc[1],"pathogen"=hc[3],"symbiont"=hc[3],"heterotroph"=hc[2],"chemoautotroph"=hc[4])
y <- as.character(rbhmatrix$metabolism[rr$ix])
bactcollist <- as.character(colmatch[y])
bactcollist[bactcollist=="NULL"] <- "grey"
points(rr$x,pch=21,col=NA,bg=bactcollist,cex=1.5)
lines(c(0,447),c(0,0),col="grey")
legend(x=20,y=200,legend=c("Photosynthetic","Pathogenic/Symbiotic","Heterotrophic","Chemoautotrophic"),pt.bg=hc[c(1,3,2,4)],pch=21)

dev.off()




# PLOTTING FITTING PARAMETERS
colmatch <- list("Hyperthermophilic"="red", "Thermophilic"="orange","Mesophilic"="green","Psychrophilic"="blue","NA"=NA)
y <- as.character(archinfo$temp)
archcollist <- as.character(colmatch[y])

colmatch <- list("Aerobic"="green","Microaerophilic"="blue","Facultative"="yellow","Anaerobic"="red")
y <- as.character(archinfo$oxygen.req)
archcollist <- as.character(colmatch[y])

plot(archinfo[,c("size..Mb.","Cb","Cc","Cd","Ce")],col=archcollist,bg=archcollist,pch=21)





# FIGURE 4/5

over.mat <- (rbhover[,1:57] > 0)+0
genes.mat <- as.matrix(over.mat*rbhover[,1:57])
facts.mat <- as.matrix(rbhover[,58:ncol(rbhover)])

archinfo$temp <- factor(archinfo$temperature,levels=c("Hyperthermophilic","Thermophilic","Mesophilic"))
rbhmatrix$temperature <- factor(rbhmatrix$temperature,levels=c("Hyperthermophilic","Thermophilic","Mesophilic"))
archinfo$oxygen.req <- factor(archinfo$oxygen,levels=c("Anaerobic","Facultative","Microaerophilic","Aerobic"))
rbhmatrix$oxygen <- factor(rbhmatrix$oxygen,levels=c("Anaerobic","Facultative","Microaerophilic","Aerobic"))

reord.bact <- order(rbhmatrix$temperature,rbhmatrix$oxygen,rev(rowSums(over.mat)))
reord.arch <- order(archinfo$temperature,archinfo$oxygen,rev(colSums(over.mat)))


genes.red <- genes.mat[reord.bact,reord.arch]
reord.bf <- facts.mat[reord.bact,]
reord.bf <- as.data.frame(reord.bf[rowSums(genes.red>0)>1,])
genes.red <- genes.red[rowSums(genes.red>0)>1,colSums(genes.red>0)>1]

reord.af <- archinfo[reord.arch,]



colmatch <- list("Aerobic"="green","Microaerophilic"="blue","Facultative"="yellow","Anaerobic"="red","unknown"="grey")
y <- as.character(reord.af$oxygen)
archcollist <- as.character(colmatch[y])
y <- as.character(reord.bf$oxygen)
bactcollist <- as.character(colmatch[y])

colmatch <- list("Hyperthermophilic"="red", "Thermophilic"="orange","Mesophilic"="green","Psychrophilic"="blue","unknown"="grey")
y <- as.character(reord.af$temperature)
archcollist <- as.character(colmatch[y])
y <- as.character(reord.bf$temperature)
bactcollist <- as.character(colmatch[y])

colmatch <- list("photosynthesis"="green","pathogen"="violet","symbiont"="violet","heterotroph"="blue","chemoautotroph"="red","unknown"="grey")
y <- as.character(reord.af$metabolism)
archcollist <- as.character(colmatch[y])
y <- as.character(reord.bf$metabolism)
bactcollist <- as.character(colmatch[y])


hm <- heatmap.2(genes.red,density.info="none",trace="none")

svg(file=paste("./figs/",writefile,"-heatmap.svg",sep=""),width = 6, height = 6,pointsize=12,bg="white")

heatmap.2(
	genes.red,
	density.info="none", trace="none", margins=c(12,12),
	RowSideColors=bactcollist,
	labRow=reord.bf$metabolism.detail,
	ColSideColors=archcollist,
	labCol=reord.af$metabolism.detail
	)
	
dev.off()

write.table(file="heatmap.csv",genes.red[hm$rowInd,hm$colInd])

#image(t(genes.red))
#text((0:56)/56,0.05,substr(reord.af$temp,0,1),srt=90)
#text((0:56)/56,0.1,substr(reord.af$oxygen.req,0,2),srt=90)
#text(0.1,(0:(nrow(reord.bf)-1))/(nrow(reord.bf)-1),substr(reord.bf$temperature,0,1),cex=0.7)
#text(0.05,(0:(nrow(reord.bf)-1))/(nrow(reord.bf)-1),substr(reord.bf$oxygen,0,2),cex=0.7)


# IGRAPH
library(igraph)

over.mat <- (rbhover[,1:57] > 0)+0
genes.mat <- as.matrix(over.mat*rbhover[,1:57])
facts.bact <- as.data.frame(rbhover[,58:ncol(rbhover)])
facts.arch <- archinfo

genes.red <- genes.mat[rowSums(genes.mat>0)>1,colSums(genes.mat>0)>1]

arcbac <- genes.red
#arcbac[which(arcbac==0)] <- NA
arcarc <- as.data.frame(matrix(0, nrow=ncol(genes.red),ncol=ncol(genes.red)))
bacbac <- as.data.frame(matrix(0, nrow=nrow(genes.red),ncol=nrow(genes.red)))
bacarc <- t(arcbac)
colnames(arcarc) <- colnames(arcbac)
rownames(arcarc) <- colnames(arcbac)
colnames(bacbac) <- rownames(arcbac)
rownames(bacbac) <- rownames(arcbac)
mat.top <-cbind(arcarc,bacarc)
mat.bot <- cbind(arcbac,bacbac)
aabb <- as.matrix(rbind(mat.top,mat.bot))
#aabb <- aabb^4

g.ab <- graph.adjacency(aabb,mode="directed",weighted=TRUE, add.colnames=TRUE)
g.ab <- set.vertex.attribute(g.ab,name="temperature",value=c(as.character(facts.arch$temperature),as.character(facts.bact$temperature)[rowSums(genes.mat>0)>1]))
g.ab <- set.vertex.attribute(g.ab,name="salinity",value=c(as.character(facts.arch$salinity),as.character(facts.bact$salinity)[rowSums(genes.mat>0)>1]))
g.ab <- set.vertex.attribute(g.ab,name="oxygen",value=c(as.character(facts.arch$oxygen),as.character(facts.bact$oxygen)[rowSums(genes.mat>0)>1]))
g.ab <- set.vertex.attribute(g.ab,name="metabolism",value=c(as.character(facts.arch$metabolism),as.character(facts.bact$metabolism)[rowSums(genes.mat>0)>1]))
g.ab <- set.vertex.attribute(g.ab,name="isolation",value=c(as.character(facts.arch$isolation),as.character(facts.bact$isolation)[rowSums(genes.mat>0)>1]))

vertex_colors = get.vertex.attribute(g.ab,"DEPT")
colors = c('Black', 'Red', 'Blue', 'Yellow', 'Green')
dept_vertex_colors[dept_vertex_colors == 0] = colors[1]

pdf("testing-graph.pdf")
plot(g.ab,
	layout=layout.fruchterman.reingold,
	vertex.label=NA, 
    edge.arrow.size=.5,
    vertex.size=0.5
    )
dev.off()

