
## deal with phenotype values
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
Y=read.table("Trait_all.txt",head=T)
taxa=paste("P_",as.character(Y[,1]),sep="")
y=Y[,-1]
nn=c(1:5)
ynn=y[,nn]
ynn2=apply(ynn,2,toupper)
letter2num<-function(one)
{
  char=sort(unique(one))
  char=char[!is.na(char)]
  if(length(char)==2) num=c(1,0)
  if(length(char)==3) num=c(-1,0,1)

  for(i in 1:length(char))
  {
    index=one==char[i]
    one[index]=num[i]

  }
return(one)
}
ynn3=apply(ynn2,2,letter2num)
numY=y[,-nn]
repY=numY[,4:41]
newY=cbind(as.data.frame(taxa),ynn3,numY)
write.table(newY,"newalltraits.txt",quote=F,row.names=F)



###### read phenotype file and remove outliers in MAC 
# setwd("/home/jiabowang/data/Xiongwei")
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# all.Y=read.table("/home/jiabowang/data/Xiongwei/newalltraits.txt",head=T) # 200
newY=read.table("newalltraits.txt",head=T)
Y1=newY[,c(2:6,8)]
Y2=newY[,c(13,15,21,23,25,27,29,31,33,35,39,41,43,45)]
index=grep("_2019",as.character(colnames(Y2)))
Y2=Y2[,index]
y2.taxa=as.character(colnames(Y2))
y2.taxa=gsub("_2019","",y2.taxa)
y.taxa=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Fullnames.txt",head=T)
y.taxa0=y.taxa$Fullnames[match(y2.taxa,y.taxa$TraitNO)]

y1.taxa=as.character(colnames(Y1))
y1.taxa=gsub("_2019","",y1.taxa)
y1.taxa=gsub("_2018","",y1.taxa)
y.taxa1=y.taxa$Fullnames[match(y1.taxa,y.taxa$TraitNO)]

colnames(Y2)=y.taxa0
colnames(Y1)=y.taxa1
Y2[,4]=Y2[,4]/5
Y2[,7]=Y2[,7]/5

# par(mfrow=c(2,1), mar = c(3,4,1,1))
# boxplot(Y2[,-1])
rm.out=function(x,na.rm=TRUE,...){
	qnt=quantile(x,probs=c(.25,.75),na.rm=na.rm,...)
	H=1.5*IQR(x,na.rm=na.rm)
	y=x
	y[x<=(qnt[1]-H)]=NA
	y[x>=(qnt[1]+H)]=NA
    return(y)	
}

normalzited=function(x,na.rm=TRUE,...){
	x.mean=mean(x,na.rm=TRUE)
	y=(x-x.mean)/sd(x,na.rm=TRUE)
	return(y)
}
Y3=apply(Y2,2,rm.out)
Y3=apply(Y3,2,normalzited)
# newY=cbind(Y1,Y3)
# write.table(newY,"newalltraits2.txt",quote=F,row.names=F)
# odd.at=seq(1,length(y.taxa0),2)
# even.at=seq(2,length(y.taxa0),2)

par(mfrow=c(2,1))
par(mar = c(1,6,8,4))
boxplot(Y2,axes=F)
    # axis(1, at=odd.at,acex.axis=1,las=3,labels=y.taxa0[odd.at],tick=T,gap.axis=0.25)
    axis(2, las=1,cex.axis=0.8,tick=T,gap.axis=0.25)

# boxplot(Y2[,-1],theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)))
# mtext("Observed values",side=2,line=3.5,col="black")
mtext("Orignal Phenotype",side=2,line=3.5,col="black")
# theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
# axis(1,at=1:ncol(Y2),labels=NA)
# text(colnames(Y2),par("usr")[3],labels=colnames(Y2),srt=315,
# 	xpd=TRUE,adj=c(-0.2,1.2),cex=0.65)
par(mar = c(9,6,0,4))
y.taxa0[1]="Hexyl_Acetate"
y.taxa0[2]="Cis-3-Hexenyl_Acetate"
y.taxa0[4]="1-Hexanol"
y.taxa0[6]="Trans-2-Hexen-1-Ol"
y.taxa0[7]="Cis-3-Hexen-1-Ol"
y.taxa0[9]="β-Ionone"

boxplot(Y3,las=1,axes=F,)
    axis(1, at=seq(1,length(y.taxa0),1),las=3,cex.axis=0.9,labels=y.taxa0,tick=T,gap.axis=0.01)
    axis(2, las=1,cex.axis=0.8,tick=T,gap.axis=0.25)
# mtext("Observed values",side=2,line=3.5,col="black")
mtext("Removed Phenotype",side=2,line=3.5,col="black")
# theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))

newY0=cbind(Y1,Y3)
oldY0=cbind(Y1,Y2)
oldY0=cbind(newY[,1],oldY0)
colnames(oldY0)[1]=c("Taxa")
write.csv(oldY0,"oldalltraits.csv",quote=F,row.names=F)

y0=newY0
yy=apply(y0,1,mean)

cor.matrix=matrix(NA,ncol(y0),ncol(y0))
for(i in 1:ncol(y0))
{
	for(j in 1:ncol(y0))
	{
		store=y0[,c(i,j)]
		store0=store[!is.na(apply(store,1,mean)),]
		cor.matrix[i,j]=cor(store0[,1],store0[,2])
	}
}
###############
library(PerformanceAnalytics)#加载包
taxa=colnames(y0)
taxa[1:6]=c("Skin_Hairiness","Fruit_Shape","Flesh_Color","Pollen_Sterility","Flower_Type","Soluble_Solid_Content")
taxa[-c(1:6)]=y.taxa0
colnames(y0)=taxa
chart.Correlation(y0, histogram=TRUE, pch=19)



library(corrplot)
# y.taxa=as.character(colnames(y0))
# y.taxa=gsub("_2018","",y.taxa)
# y.taxa=gsub("_2019","",y.taxa)
# setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# y.taxa0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Fullnames.txt",head=T)
# y.taxa0=y.taxa0[1:27,]
# y.taxa1=y.taxa0[c(1:5,8),]
# y.taxa2=y.taxa0[-c(1:5,8),]
# y.taxa3=y.taxa0[-c(1:5,8),]
# y.taxa2[,1]=paste(y.taxa2[,1],"_2018",sep="")
# y.taxa2[,2]=paste(y.taxa2[,2],"_2018",sep="")
# y.taxa3[,1]=paste(y.taxa3[,1],"_2019",sep="")
# y.taxa3[,2]=paste(y.taxa3[,2],"_2019",sep="")
# y.taxa0=rbind(y.taxa1,y.taxa2,y.taxa3)

# y.taxa0=y.taxa0$Fullnames[match(y.taxa,y.taxa0$TraitNO)]
# match(y.taxa,y.taxa0$TraitNO)
colnames(cor.matrix)=colnames(y0)
rownames(cor.matrix)=colnames(y0)
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples")
write.csv(cor.matrix,paste("Correlation.phenotypes.Data.csv",sep=""),quote=F)
# cor.matrix=read.csv("Correlation.phenotypes.Data.csv",head=T)

par(mfrow=c(1,1), mar = c(1,1,1,1))

corrplot(cor.matrix,
	method='ellipse',
	type="upper",
	tl.cex=0.7,
	order='hclust',
	sig.level=0.05)
newY0=cbind(newY[,1],newY0)
colnames(newY0)=c("Taxa",colnames(y0))
# write.csv(newY0,"newalltraits2.csv",quote=F,row.names=F)

### count the number of NA
count.NA=function(x,...){
	x0=x[is.na(x)]
	y=length(x0)
	return(y)
}
y.na=apply(y0,2,count.NA)


####  GWAS with each trait

rm(list=ls())
source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV=read.table("peachCV.txt",head=T)
# chr=as.character(unique(myGM[,2]))
# allchr=as.character(myGM[,2])
# chr=chr[order(chr)]
# for(i in 1:length(chr))
# {
# 	allchr[allchr==chr[i]]=i
# }
# myGM[,2]=allchr
# colnames(myGM)[2]="Chr"
# setwd("/home/jiabowang/data/Xiongwei/rawdata")
# write.table(myGM,"peachGM.txt",quote=F,row.names=F)
myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
# setwd("/home/jiabowang/data/Xiongwei/rawdata")
# cv.index=read.table("population.txt",head=T)

# myCV=read.table("peachCV.txt",head=T)
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

for(i in 2:ncol(myY))
{
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")
system(paste("mkdir ",colnames(myY)[i],sep=""))
setwd(paste("/home/jiabowang/data/Xiongwei/Multiple-traits/",colnames(myY)[i],sep=""))

myGAPIT <- GAPIT(
Y=myY[,c(1,i)],
GD=myGD,
GM=myGM,
model=c("MLM","MLMM","FarmCPU","BLINK"),
CV=myCV,
# Inter.Plot=T,
Multiple_analysis=TRUE,
file.output=T
)

}

setwd("/home/jiabowang/data/Xiongwei")
system("tar -cvf Multiple-traits.tar Multiple-traits")










source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
# myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

Y.taxa=colnames(myY)[-1]

for(i in 2:ncol(myY))
{
setwd(paste("/home/jiabowang/data/Xiongwei/Multiple-traits/",colnames(myY)[i],sep=""))

GMM=GAPIT.Multiple.Manhattan(model_store=c("MLM","MLMM","FarmCPU","Blink"),Y=myY[,c(1,i)],cutOff=0.01,GM=myGM)

}









########## Zhiliang traits GS after GWAS 6 traits (SH,FS,FC,PS,FT,SSC)


rm(list=ls())
# source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV=read.table("peachCV.txt",head=T)

myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Zhiliang")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
Y0=myY[,c(1:3,4,5,6,8,31)]
nfold=5
nrep=10
cutoff=0.01
n=nrow(myGM)
taxa.g=as.character(myGD[,1])
taxa.y=as.character(Y0[,1])
yourGD=myGD[taxa.g%in%taxa.y,]
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
taxa.cv=as.character(myCV[,1])
taxa.g2=as.character(yourGD[,1])
yourCV=myCV[match(taxa.g2,taxa.cv)[!is.na(match(taxa.g2,taxa.cv))],]

taxa=as.character(yourY[,1])
acc.traits=NULL
x0=yourGD[,-1]
set.seed(99163)
for(i in 2:ncol(yourY))
{
	i=8
   y0=yourY[,c(1,i)]
   y.index=!is.na(y0[,2])
   y0=y0[y.index,]
   x00=x0[y.index,]
   yourGD0=yourGD[y.index,]
   yourCV0=yourCV[y.index,]
   taxa0=taxa[y.index]
   acc.rep=NULL
   for(j in 1:nrep)
   {
     sets=sample(cut(1:nrow(y0),nfold,labels=FALSE),nrow(y0))
     acc.fold=NULL
     for(k in 1:nfold)
     {
     	training=y0
     	training[sets==k,2]=NA
     	training_index=is.na(training[,2])
     	testing=y0[training_index,]
        myGAPIT <- GAPIT(
               Y=training[,c(1,2)],
               GD=yourGD0,
               GM=myGM,
               model=c("BLINK"),
               CV=yourCV0,
               # Inter.Plot=T,
               Multiple_analysis=FALSE,
               file.output=F
               )
        gwas1=myGAPIT$GWAS
        sig1=gwas1[gwas1[,4]<cutoff/n,]
        sig.index1=as.character(myGM[,1])%in%as.character(sig1[,1])
        myCV1=cbind(yourCV0,x00[,sig.index1])
        ANS=lm(training[!training_index,2]~as.matrix(myCV1[!training_index,-1]))
        pred=as.matrix(myCV1[training_index,-1])%*%as.matrix(ANS$coefficients)[-1,]
        pred1=cbind(as.data.frame(taxa0[training_index]),pred)
        colnames(pred1)=c("Taxa","Pred")
        ANS.pred=merge(testing,pred1,by.x="Taxa",by.y="Taxa")
        ANS.pred=ANS.pred[!is.na(apply(ANS.pred[,-1],1,mean)),]
        acc.fold=append(acc.fold,cor(ANS.pred[,2],ANS.pred[,3]))
     }
     acc.rep=append(acc.rep,mean(acc.fold,na.rm=TRUE))
    }
    acc.traits=append(acc.traits,mean(acc.rep,na.rm=TRUE))
}

acc.traits=matrix(acc.traits,1,length(acc.traits))
colnames(acc.traits)=colnames(Y0)[-1]

write.csv(acc.traits,"Acc_after_GWAS.csv",quote=F,row.names=F)




############### GS with GWAS result for 7 continue traits
##### 33,15,35,39,41,43,45
##### T12,T21,T22,T24,T25,T26,T27

rm(list=ls())
source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV=read.table("peachCV.txt",head=T)

myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Zhiliang")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
Y0=myY[,c(1,33,15,35,39,41,43,45)]
# Y0=myY[,c(1,31)]
nfold=5
nrep=30
cutoff=0.01
n=nrow(myGM)
taxa.g=as.character(myGD[,1])
taxa.y=as.character(Y0[,1])
yourGD=myGD[taxa.g%in%taxa.y,]
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
taxa.cv=as.character(myCV[,1])
taxa.g2=as.character(yourGD[,1])
yourCV=myCV[match(taxa.g2,taxa.cv)[!is.na(match(taxa.g2,taxa.cv))],]

taxa=as.character(yourY[,1])
acc.traits=NULL
x0=yourGD[,-1]
set.seed(99163)

for(i in 2:ncol(yourY))
{
	#i=3
   y0=yourY[,c(1,i)]
   y.index=!is.na(y0[,2])
   y0=y0[y.index,]
   x00=x0[y.index,]
   yourGD0=yourGD[y.index,]
   yourCV0=yourCV[y.index,]
   taxa0=taxa[y.index]
   acc.rep=NULL
   for(j in 1:nrep)
   {
     sets=sample(cut(1:nrow(y0),nfold,labels=FALSE),nrow(y0))
     acc.fold=NULL
     for(k in 1:nfold)
     {
     	training=y0
     	training[sets==k,2]=NA
     	training_index=is.na(training[,2])
     	testing=y0[training_index,]
        myGAPIT <- GAPIT(
               Y=training[,c(1,2)],
               GD=yourGD0,
               GM=myGM,
               model=c("BLINK"),
               CV=yourCV0,
               # Inter.Plot=T,
               Multiple_analysis=FALSE,
               file.output=F
               )
        gwas1=myGAPIT$GWAS
        sig1=gwas1[gwas1[,4]<cutoff/n,]
        sig.index1=as.character(myGM[,1])%in%as.character(sig1[,1])
        myCV1=cbind(yourCV0,x00[,sig.index1])
        ANS=lm(training[!training_index,2]~as.matrix(myCV1[!training_index,-1]))
        pred=as.matrix(myCV1[training_index,-1])%*%as.matrix(ANS$coefficients)[-1,]
        pred1=cbind(as.data.frame(taxa0[training_index]),pred)
        colnames(pred1)=c("Taxa","Pred")
        ANS.pred=merge(testing,pred1,by.x="Taxa",by.y="Taxa")
        ANS.pred=ANS.pred[!is.na(apply(ANS.pred[,-1],1,mean)),]
        acc.fold=append(acc.fold,cor(ANS.pred[,2],ANS.pred[,3]))
     }
     acc.rep=append(acc.rep,mean(acc.fold,na.rm=TRUE))
    }
    acc.traits=append(acc.traits,mean(acc.rep,na.rm=TRUE))
}

acc.traits=matrix(acc.traits,1,length(acc.traits))
colnames(acc.traits)=colnames(Y0)[-1]

write.csv(acc.traits,"All_Acc_after_GWAS.csv",quote=F,row.names=F)


############## GWAS T21 with T4 as CV



rm(list=ls())
source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV=read.table("peachCV.txt",head=T)

myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
# setwd("/home/jiabowang/data/Xiongwei/rawdata")
# cv.index=read.table("population.txt",head=T)

# myCV=read.table("peachCV.txt",head=T)
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

taxa.cv=as.character(myCV[,1])
taxa.g2=as.character(yourGD[,1])

yourCV=merge(myY[,c(1,4)],myCV,by.x="Taxa",by.y="taxa")
yourCV=yourCV[!is.na(yourCV[,2]),]
taxa.cv=as.character(yourCV[,1])
taxa.y=as.character(myY[,1])
yourY=myY[taxa.y%in%taxa.cv,]

i=33
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")
# system(paste("mkdir ",colnames(myY)[i],sep=""))
setwd(paste("/home/jiabowang/data/Xiongwei/Multiple-traits/",colnames(myY)[i],sep=""))

myGAPIT <- GAPIT(
Y=yourY[,c(1,i)],
GD=myGD,
GM=myGM,
# model=c("MLM","MLMM","FarmCPU","BLINK"),
model=c("BLINK"),
cutOff=10,
FDRcut=TRUE,
CV=yourCV[,c(1,2)],
# Inter.Plot=T,
Multiple_analysis=TRUE,
file.output=T
)

GAPIT.Manhattan(GI.MP = myGAPIT$GWAS[,2:4],GD=myGD[,-1],name.of.trait = colnames(yourY[,i]), plot.type = "Chromosomewise",cutOff=0.01)

GAPIT.Manhattan(GI.MP = myGAPIT$GWAS[,2:4],name.of.trait = colnames(yourY[,i]), plot.type = "Genomewise",cutOff=0.01)

GAPIT.QQ(P.values = myGAPIT$GWAS[,4], name.of.trait = colnames(yourY[,i]))



########## Trait 21 GS after GWAS

rm(list=ls())
source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV0=read.table("peachCV.txt",head=T)

myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
i=33
nrep=30
nfold=5
setwd(paste("/home/jiabowang/data/Xiongwei/Multiple-traits/",colnames(myY)[i],sep=""))
result=read.csv(paste("GAPIT.Blink.",colnames(myY)[i],".GWAS.Results.csv",sep=""),head=T)
bor=1/nrow(myGM)
sig=result[result[,4]<bor,]
write.csv(sig,"Sig.matrix_T21.csv",quote=F,row.names=F)

sig.index=as.character(myGM[,1])%in%as.character(sig[,1])
# setwd("/home/jiabowang/data/Xiongwei/rawdata")
# cv.index=read.table("population.txt",head=T)

myCV=merge(myCV0,myGD[,c(TRUE,sig.index)],by.x="taxa",by.y="taxa")
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

taxa.cv=as.character(myCV[,1])
taxa.g2=as.character(yourGD[,1])

yourCV=merge(myY[,c(1,4)],myCV,by.x="Taxa",by.y="taxa")
yourCV=yourCV[!is.na(yourCV[,2]),]
taxa.cv=as.character(yourCV[,1])
taxa.y=as.character(myY[,1])
yourY=myY[match(taxa.cv,taxa.y),c(1,i)]

# i=33
set.seed(99163)
# for(i in 2:ncol(yourY))
# {
	#i=3
   y0=yourY
   y.index=!is.na(y0[,2])
   y0=y0[y.index,]
   x00=x0[y.index,]
   yourGD0=yourGD[y.index,]
   yourCV0=yourCV[y.index,]
   taxa0=as.character(yourCV0[,1])
   acc.rep=NULL

   nrep=1
   nfold=nrow(y0)

   for(j in 1:nrep)
   {
     sets=sample(cut(1:nrow(y0),nfold,labels=FALSE),nrow(y0))
     acc.fold=NULL
     pre.store=NULL
     for(k in 1:nfold)
     {
     	training=y0
     	training[sets==k,2]=NA
     	training_index=is.na(training[,2])
     	testing=y0[training_index,]
        ANS=lm(training[!training_index,2]~as.matrix(yourCV0[!training_index,-1]))
        pred=as.matrix(yourCV0[training_index,-1])%*%as.matrix(ANS$coefficients)[-1,]
        pred1=cbind(as.data.frame(taxa0[training_index]),pred)
        colnames(pred1)=c("Taxa","Pred")
        ANS.pred=merge(testing,pred1,by.x="Taxa",by.y="Taxa")
        ANS.pred=ANS.pred[!is.na(apply(ANS.pred[,-1],1,mean)),]
        # acc.fold=append(acc.fold,cor(ANS.pred[,2],ANS.pred[,3]))
        pre.store=rbind(pre.store,ANS.pred)
     }
     acc.rep=append(acc.rep,cor(pre.store[,2],pre.store[,3]))
    }
    # acc.traits=append(acc.traits,mean(acc.rep,na.rm=TRUE))
# }

# acc.traits=matrix(acc.traits,1,length(acc.traits))
# colnames(acc.traits)=colnames(Y0)[-1]
setwd("/home/jiabowang/data/Xiongwei/Zhiliang")
#0.68
write.csv(acc.rep,"Acc_T21.with.T4.csv",quote=F,row.names=F)


########## Trait 8,12,22,24,25,26,27 GS after GWAS using all significant markers from all gwas

licols <- function(X,tol=1e-10){
  if (all(X==0)){     # X is a zero matrix        
    idx <- Xsub <- c()
  }else{              # X is not a zero matrix 
    qr_res <- qr(X,LAPACK=F) # QR decomposition 
    Q <- qr.Q(qr_res)
    R <- qr.R(qr_res)
    E <- qr_res$pivot
    
    if (is.vector(R) == 0){
      diagr <- abs(diag(R))
    }else{
      diagr <- abs(R[1])
    }
    
    # Rank estimation
    r <- sum(diagr >= tol * diagr[1])
    
    idx <- sort(E[1:r])
    Xsub <- X[,idx]
  }
  res <- vector("list")
  res$Xsub <- Xsub
  res$idx <- idx
  return(res)
}
# rm(list=ls())
# source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV0=read.table("peachCV.txt",head=T)
# myY=read.table("newalltraits2.txt",head=T) # 200

myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits")

index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]

CV.Y=merge(myCV0,myY,by.x="taxa",by.y="Taxa")
yourCV=CV.Y[,1:4]
yourY=CV.Y[,-c(2:4)]
taxa.cv=as.character(yourCV[,1])
taxa.y=as.character(yourY[,1])
taxa.g=as.character(myGD[,1])
# taxa.y=as.character(yourY[,1])
yourGD=myGD[taxa.g%in%taxa.y,]

# i0=c(2:8,15,33,31,35,39,41,43,45)
# i0=c(2:8)
i0=c(8,15,33,31,35,39,41,43,45)
# i0=c(31)

# nrep=30
# nfold=5
setwd("/home/jiabowang/data/Xiongwei/traits-gwas")
# result0=read.table("GAPIT.Filter_GWAS_results.txt",head=T)
result0=read.csv("GAPIT.Association.Filter_GWAS_results.csv",head=T)[,-1]
# GAPIT.Association.Filter_GWAS_results_complex8

all.acc=NULL
for(i in i0)
{
	# sig0=result0[-c(grep(colnames(myY)[i],as.character(result0[,6]))),]
	sig0=result0
	sig0.rs=unique(as.character(sig0[,1]))
	sig.index=as.character(myGM[,1])%in%sig0.rs
	CV0=merge(yourCV,yourGD[,c(TRUE,sig.index)],by.x="taxa",by.y="taxa")
    
   y0=yourY[,c(1,i)]
   y.index=!is.na(y0[,2])
   y0=y0[y.index,]
   # x00=x0[y.index,]
   yourGD0=yourGD[y.index,]
   yourCV0=CV0[y.index,]
   CV=licols(yourCV0[,-1])$Xsub
   CV=CV[,apply(CV,2,var)!=0]
   yourCV00=cbind(as.data.frame(yourCV0[,1]),CV)
   taxa0=as.character(yourCV00[,1])
   acc.rep=NULL

   nrep=1
   nfold=nrow(y0)

   for(j in 1:nrep)
   {
     sets=sample(cut(1:nrow(y0),nfold,labels=FALSE),nrow(y0))
     acc.fold=NULL
     pre.store=NULL
     for(k in 1:nfold)
     {
     	training=y0
     	training[sets==k,2]=NA
     	training_index=is.na(training[,2])
     	testing=y0[training_index,]
        ANS=lm(training[!training_index,2]~as.matrix(yourCV00[!training_index,-1]))
        pred=as.matrix(yourCV00[training_index,-1])%*%as.matrix(ANS$coefficients)[-1,]
        pred1=cbind(as.data.frame(taxa0[training_index]),pred)
        colnames(pred1)=c("Taxa","Pred")
        ANS.pred=merge(testing,pred1,by.x="taxa",by.y="Taxa")
        ANS.pred=ANS.pred[!is.na(apply(ANS.pred[,-1],1,mean)),]
        # acc.fold=append(acc.fold,cor(ANS.pred[,2],ANS.pred[,3]))
        pre.store=rbind(pre.store,ANS.pred)
     }#end of k
     acc.rep=append(acc.rep,cor(pre.store[,2],pre.store[,3]))
   }# end of j
    # acc.t
   write.csv(acc.rep,paste("wi.Acc_",colnames(myY)[i],".csv",sep=""),quote=F,row.names=F)
   write.csv(pre.store,paste("wi.Pre.store_",colnames(myY)[i],".csv",sep=""),quote=F,row.names=F)
   all.acc=append(all.acc,acc.rep)
}# end of i
all.acc=as.matrix(all.acc,1,6)
rownames(all.acc)=colnames(myY)[i0]
   write.csv(all.acc,paste("wi.All_Acc.csv",sep=""),quote=F)
   write.csv(all.acc,paste("wi.zhiliang_Acc.csv",sep=""),quote=F)





############ Plot synthesis manhatton plot.


rm(list=ls())
setwd("/home/jiabowang/data/Xiongwei/rawdata")
# myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
myY=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
index=myGM[,2]>8
myGM=myGM[!index,]
# 
type="complex"
if(type=="simple")
{
    y.names=colnames(myY)[c(2,3,4,5,6,8)]
    y.names=gsub("_2019","",gsub("_2018","",y.names))

    setwd("/home/jiabowang/data/Xiongwei/traits-gwas/simple1-2-4-5-6")
}else{
	y.names=colnames(myY)[c(15,31,33,35,39,41,43,45)]
    setwd("/home/jiabowang/data/Xiongwei/traits-gwas/complex8")
    y.names=gsub("_2019","",gsub("_2018","",y.names))
}
#simple 2,3,4,5,6
#complex 8,15,33,35,39,41,43,45
# setwd("/home/jiabowang/data/Xiongwei/traits-gwas/complex8")

# source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")
source("/home/jiabowang/Code/GAPIT.Multiple.Manhattan.xiongwei.R")

if(type=="simple")
{
    GMM=GAPIT.Multiple.Manhattan(model_store=c("MLM","MLMM","FarmCPU","Blink"),
	Y.names=y.names,cutOff=0.05,DPP=500,WS=3e5,inpch=c("+","*","-","#","<",">","^","$"),
	outpch=c(0,1,2,5),byTraits=TRUE,
	GM=myGM,plot.type=c("s"))
}else{
	GMM=GAPIT.Multiple.Manhattan(model_store=c("MLM","MLMM","FarmCPU","Blink"),
	Y.names=y.names,cutOff=0.05,DPP=500,WS=3e5,inpch=c("+","*","-","#","<",">","^","$"),
	outpch=c(0,1,2,5),byTraits=TRUE,
	GM=myGM,plot.type=c("s"))
}


# GAPIT.Multiple_Synthesis(model_store=c("Blink","FarmCPU"),
# 	Y.names=y.names,cutOff=0.05,GM=myGM)

WS=1500000
######### Find common genes in GWAS results

# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples/symphysic manhattan plots")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/Multiple Manhattan")

type="simple"
if(type=="complex")
{
	results=read.csv("GAPIT.Association.Filter_GWAS_results_complex8.csv",head=T)[,-1] # 12 in 43
}else{
	results=read.csv("GAPIT.Association.Filter_GWAS_results_simple6.csv",head=T)[,-1]# 47 in 82
}
results2=results[order(results[,3]),]
results2=results2[order(results2[,2]),]
results2[,c(1,2,3,6)]
results2=results2[!duplicated(results2[,1]),]
row.names(results2)=1:nrow(results2)

scom=as.numeric(results2[,3])
de.sc=scom[-1]-scom[-length(scom)]
dayu1.index=c(abs(de.sc)<WS,FALSE)
# dayu2.index=c(FALSE,de.sc<3e5)
index=sort(unique(c(which(dayu1.index),which(dayu1.index)+1)))
results3=results2[index,]
results3[,c(1,2,3,6)]
write.csv(results3,paste("GAPIT.",type,".common_markers_in_",WS,".csv",sep=""),quote=FALSE,row.names=F)


######### Phenotype distribution following significant genotype zhiliang data

rm(list=ls())
library(fmsb)
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
# myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
load("peah.hpm.RData")
myY=read.table("newalltraits2.txt",head=T) # 200

Y0=myY[,c(1,2,3,4,5,6,8)]


# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples/symphysic manhattan plots")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/Multiple Manhattan")
# results=read.csv("GAPIT.single.common_markers_in_3e+05.csv",head=T)
results=read.csv("GAPIT.Association.Filter_GWAS_results_simple6.csv",head=T)[,-1]

results0=results[,c(1,2,3,6)]
# traits=gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(results0[,4])))))
traits=gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(colnames(Y0)[-1])))))
traits0=gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(results0[,4])))))
traits1=c("Skin_Hairiness","Fruit_Shape","Flesh_Color","Pollen_Sterility","Flower_Type","Soluble_Solid_Content")
hapmap2=hapmap[-1,]
map=hapmap2[,c(1,3,4)]
X0=hapmap2[,-c(1:11)]

n=nrow(map)
taxa.g=as.character(hapmap[1,-c(1:11)])
taxa.y=as.character(Y0[,1])
yourGD=X0[,taxa.g%in%taxa.y]
taxa.g2=taxa.g[taxa.g%in%taxa.y]
colnames(yourGD)=taxa.g2
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]

# par(mfrow=c(2,3))
par(mar = c(5,4,3,1))
col.type=c("lightblue","mistyrose","lavender")
legend0=list()
legend0[[1]]=c("Nectarine","Peach")
legend0[[2]]=c("Round","Flat")
legend0[[3]]=c("Red","White","Yellow")
legend0[[4]]=c("Fertility","Sterility")
legend0[[5]]=c("Non-Showy","Showy")

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures")

# graphics::layout(matrix(c(1,1,2,1,1,1,1,1,1),3,3,byrow=TRUE), c(2,1), c(1,1), TRUE)
split.screen(c(2,1))
split.screen(c(1,3),screen=1)
split.screen(c(1,2),screen=2)

# layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5),2,6,byrow=TRUE), c(1,1), c(1,1), TRUE)
# par(mar = c(4,4,3,4))

com=2
ca=0.8
for(i in c(4,3,2,5,1))
{
	com=com+1
	screen(com)
	results2=results0[traits0==traits[i],]
	results2=results2[!duplicated(as.character(results2[,1])),]
    targets=as.character(results2[,1])
	target0=targets
	if(i==1)target0=target0[1:7]
	if(i==2)target0=target0[c(2,3,5)]
	if(i==3)target0=target0[1:2]

	index0=as.character(map[,1])%in%target0
	tar.marker=yourGD[index0,]
    marker.names=target0
    allmn=marker.names[1]
    if(length(marker.names)>1)
    {
       for(j in 2:length(marker.names))
       {
    	   allmn=paste(allmn,marker.names[j],sep="&")
       }
    }
    tar.marker[tar.marker=="AG"]="GA"
    tar.marker[tar.marker=="AC"]="CA"
    tar.marker[tar.marker=="TC"]="CT"
    markers=tar.marker[1,]
    if(nrow(tar.marker)>1)
    {
       for(j in 2:nrow(tar.marker))
       {
           markers=paste(markers,tar.marker[j,],sep="&")
        }
    }
    markers0=as.character(unique(as.character(markers)))
    # print(length(markers0))
    yall=NULL
    yi=as.numeric(yourY[,i+1])
    yy=as.data.frame(unique(yi[!is.na(yi)]))
    names(yy)="yy"
    y.type=length(unique(yi[!is.na(yi)]))
    for(k in 1: length(markers0))
    {
    	
    	yk=yourY[markers==markers0[k],i+1]
        yk=yk[!is.na(yk)]
        if(length(unique(yk))>1)
        {
        	yk.m=as.matrix(table(yk))
        	yk.mm=cbind(as.data.frame(rownames(yk.m)),yk.m)
            colnames(yk.mm)=c("type","values")
            yk.mer=merge(yy,yk.mm,by.x="yy",by.y="type",all=TRUE)
        	perc.k=round(as.numeric(yk.mer[,2])/(sum(as.numeric(yk.mer[,2]),na.rm=TRUE)),2)
        }else{
        	perc.k=rep(NA,y.type)
        }
        # perc.k
        yall=rbind(yall,perc.k)
    }
    index=!is.na(apply(yall,1,function(one) sum(one,rm.na=TRUE)))
    bar.y=t(yall[index,])
    colnames(bar.y)=markers0[index]
    rownames(bar.y)=sort(unique(yi[!is.na(yi)]))
    bar.y=cbind(bar.y,NA)
    xn=ncol(bar.y)
    rownames(bar.y)=as.character(legend0[[i]])
    write.csv(bar.y,paste("Fig.3A.data.",traits1[i],".csv",sep=""),quote=F)
    # main0=paste(unique(traits)[i]," (",marker.names,"&",)
    text.y=100*bar.y[,-c(xn)]
    barplot(bar.y,col=col.type[1:y.type],
    	    # main=unique(traits)[i],
    	    las=1,cex.names = 1,space=0.2,width=1,border=NA,names.arg=rep(NA,xn))
    posi=seq(0.7,1.2*(xn-1),1.2)
    odd=seq(1,xn-1,2)
    even=seq(2,xn-1,2)
    labels=markers0[index]
    axis(1,at=posi[odd],labels=labels[odd],cex.axis=ca,col="black",col.ticks="black",col.axis="black",tck=-0.02,tick=F)
    axis(3,at=posi[even],labels=labels[even],cex.axis=ca,col="black",col.ticks="black",col.axis="black",tck=-0.02,tick=F)
    
    for(t in 1:y.type)
    {
    	if(t==1)
    	{ 
    		text(x=posi,y=bar.y[t,]*0.5,labels=paste(text.y[t,],"%"),col="black",cex=1.5)
    	}else{
    	    text(x=posi,y=bar.y[t,]*0.5+apply(bar.y[1:(t-1),-c(xn),drop=FALSE],2,sum),labels=paste(text.y[t,],"%"),col="black",cex=1.5)
    	}
    }
       # legend("right",legend=paste(1:(ncol(bar.y)-1),":",as.character(markers0[index]),sep=""),col="black",
    # 	   box.col="white",bty = "n", bg = par("bg"))
    # legend("topright",legend=as.character(sort(unique(yi[!is.na(yi)]))),col=col.type[1:y.type],
    legend("topright",legend=as.character(legend0[[i]]),col=col.type[1:y.type],
    	   pch=c(15),lwd=1,cex=1,lty=0,ncol=1,
    	   box.col="white",bty = "n", bg = par("bg"))
    # barplot(bar.y,xlab="",main=unique(traits)[i],las=1,axisnames=F,cex.names = 0.6)
    mtext(paste(traits1[i],"(",allmn,")",sep=""),side=1,line=3,cex=1.2)
    if(com==3|com==6)mtext(side=2,text=c("Ratio of Phenotype Type"),line=2.5)
}


erase.screen() # 清除当前屏
close.screen(all = TRUE) # 退出分屏模式






barplot(bar.y,col=col.type[1:y.type],
	        # main=unique(traits)[i],
    	    las=1,cex.names = 1,border=NA,
    	    width=1,space=0.2,names.arg=rep(NA,7))
    # mtext(side=1,text=1:7)
posi=seq(0.7,1.2*(xn-1),1.2)
odd=seq(1,xn-1,2)
even=seq(2,xn-1,2)
labels=markers0[index]
axis(1,at=posi[odd],labels=labels[odd],col="black",col.ticks="black",col.axis="black",tck=-0.02,tick=F)
axis(3,at=posi[even],labels=labels[even],col="black",col.ticks="black",col.axis="black",tck=-0.02,tick=F)
    
mtext("Predicted phenotype value in EMMAREML",side=1,line=cl,col="black",cex=0.9)

mtext(paste("Chr:",choose.chr," (Mb)",sep=""),side=1,line=2.5,cex=1,col="black")






######### Complex traits distribution in significant markers
### Complex-traits-phenotype-distribution-significant
rm(list=ls())
library(fmsb)
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
# myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
load("peah.hpm.RData")
myY=read.table("newalltraits2.txt",head=T) # 200
# Y0=myY[,c(1,8,15,31,33,39,41,43,45)]
Y0=myY[,c(1,31)]

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples/symphysic manhattan plots")
results=read.csv("GAPIT.complex.common_markers_in_3e+05.csv",head=T)
# results0=results[c(3,4,7,8,19,20),c(1,2,3,6)]

results0=results[,c(1,2,3,6)]

# traits=gsub("_2018","",gsub("_2019","",gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(results0[,4])))))))

hapmap2=hapmap[-1,]
map=hapmap2[,c(1,3,4)]
X0=hapmap2[,-c(1:11)]

n=nrow(map)
taxa.g=as.character(hapmap[1,-c(1:11)])
taxa.y=as.character(Y0[,1])
yourGD=X0[,taxa.g%in%taxa.y]
taxa.g2=taxa.g[taxa.g%in%taxa.y]
colnames(yourGD)=taxa.g2
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
traits=gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(colnames(yourY)[-1])))))
traits0=gsub("MLM.","",gsub("MLMM.","",gsub("FarmCPU.","",gsub("Blink.","",as.character(results0[,4])))))
traits1=gsub("_2019","",gsub("_2018","",traits))
traits1[1:4]=c("Soluble_Solid_Content","Cis-3-Hexenyl_Acetate","Linalool","β_Ionone")

col.type=c("lightblue","mistyrose","lavender","lightgreen","lightgray","lightgoldenrod2","coral2","royalblue3")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures")

par(mfrow=c(2,4))
par(mar = c(5,5,4,1))

# com=2
for(i in 1:length(unique(traits)))
{
	# com=com+1
	# screen(com)
	results2=results0[traits0==traits[i],]
	results2=results2[!duplicated(as.character(results2[,1])),]
    targets=as.character(results2[,1])
	target0=targets
	# if(i==4)target0=target0[1:7]
	index0=as.character(map[,1])%in%target0
	tar.marker=yourGD[index0,]
    marker.names=target0
    allmn=marker.names[1]
    if(length(marker.names)>1)
    {
       for(j in 2:length(marker.names))
       {
    	   allmn=paste(allmn,marker.names[j],sep="&")
       }
    }
    tar.marker[tar.marker=="AG"]="GA"
    tar.marker[tar.marker=="AC"]="CA"
    tar.marker[tar.marker=="TC"]="CT"
    markers=tar.marker[1,]
    taxa=as.character(colnames(markers))
    if(nrow(tar.marker)>1)
    {
       for(j in 2:nrow(tar.marker))
       {
           markers=paste(markers,tar.marker[j,],sep="&")
        }
    }
    markers0=as.character(unique(as.character(markers)))
    print(length(markers0))
    yall=NULL
    yi=yourY[,c(1,i+1)]
    geno=cbind(taxa,as.data.frame(as.character(markers)))
    xn=length(markers0)

    yall=merge(yi,geno,by.x="Taxa",by.y="taxa")
    colnames(yall)=c("Taxa","Values","Genotype")
    aa=boxplot(Values~Genotype,data=yall,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:xn],plot=F)
    boxplot(Values~Genotype,data=yall,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:xn],axes=F)
    posi=seq(1,xn,1)
    odd=seq(1,xn,2)
    even=seq(2,xn,2)
    labels=aa$names
    axis(2,col="black",col.ticks="black",col.axis="black",tck=-0.02,las=1)
    axis(1,at=posi[odd],labels=labels[odd],col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    # axis(1,at=posi,labels=labels,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    axis(3,at=posi[even],labels=labels[even],col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    mtext(paste(traits1[i]," (",allmn,")",sep=""),side=1,line=3,cex=1.2)
    if(i==1|i==5)mtext("Phenotype Values",side=2,line=3,cex=1.2)
}






######### Radar Plots
rm(list=ls())
library(fmsb)
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
# myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
load("peah.hpm.RData")
myY=read.table("newalltraits2.txt",head=T) # 200

# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples/symphysic manhattan plots")
# results=read.csv("GAPIT.complex.common_markers_in_3e+05.csv",head=T)
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/Multiple Manhattan")
# results=read.csv("GAPIT.complex.common_markers_in_3e+05.csv",head=T)
results=read.csv("GAPIT.complex.common_markers_in_1500000.csv",head=T) # 12 in 43
# results2=results[order(results[,3]),]
# results2=results2[order(results2[,2]),]
# results2=results2[,c(1,2,3,6)]
# rownames(results2)=1:nrow(results2)
results0=results


targets=as.character(results0[,1])

hapmap2=hapmap[-1,]
map=hapmap2[,c(1,3,4)]
X0=hapmap2[,-c(1:11)]






#####  28509583
# targets=c("28509583","2721214","2663753")

Y0=myY[,c(1,15,31,35,39,41,43)]
# nfold=5
# nrep=10
# cutoff=0.01
n=nrow(map)
taxa.g=as.character(hapmap[1,-c(1:11)])
taxa.y=as.character(Y0[,1])
yourGD=X0[,taxa.g%in%taxa.y]
taxa.g2=taxa.g[taxa.g%in%taxa.y]
colnames(yourGD)=taxa.g2
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures")


set.seed(198521)
for(kk in c(1,3,5,7))
{
target0=targets[kk:(kk+1)]
# target0=targets
index0=as.character(map[,1])%in%target0
tar.marker=yourGD[index0,]
marker.names=target0
tar.marker[tar.marker=="AG"]="GA"
tar.marker[tar.marker=="CT"]="TC"
tar.marker[tar.marker=="TG"]="GT"
markers=paste(tar.marker[1,],tar.marker[2,],sep=" & ")

markers0=as.character(unique(markers))
mean.yall=NULL
for(i in 1:length(markers0))
{
	indexi=markers==markers0[i]
	yi=yourY[indexi,]
    mean.yi=apply(yi[,-1],2,function(one) mean(one,na.rm=TRUE))
    mean.yall=rbind(mean.yall,mean.yi)
}
# index0=tar.marker==0
# index1=tar.marker==1
# index2=tar.marker==2
# y0=yourY[index0,]
# y1=yourY[index1,]
# y2=yourY[index2,]

max.y=apply(yourY[,-1],2,function(one) max(one,na.rm=TRUE))
min.y=apply(yourY[,-1],2,function(one) min(one,na.rm=TRUE))

# mean.y0=apply(y0[,-1],2,function(one) mean(one,na.rm=TRUE))
# mean.y1=apply(y1[,-1],2,function(one) mean(one,na.rm=TRUE))
# mean.y2=apply(y2[,-1],2,function(one) mean(one,na.rm=TRUE))


data0=rbind(max.y,min.y,mean.yall)
ytaxa=colnames(yourY)[-1]
ytaxa=gsub("_2019","",gsub("_2018","",ytaxa))
colnames(data0)=ytaxa
rownames(data0)=c("1","2",paste("Genotype: ",markers0,sep=""))
index.na=apply(data0,1,function(one) sum(one=="NaN")>4)
data0=as.data.frame(data0)[!index.na,]

colors_border=NULL
colors_in=NULL
red0=seq(0.1,0.9,length.out=nrow(data0)-2)
green0=seq(0.9,0.1,length.out=nrow(data0)-2)
# blue0=seq(0.9,0.1,length.out=nrow(data0)-2)
blue0=sample(seq(0,1,length.out=5),nrow(data0)-2,replace=TRUE)
for(j in 1:(nrow(data0)-2))
{
	colors_border=c(colors_border,rgb(red0[j],green0[j],blue0[j],0.9))
	colors_in=c(colors_in,rgb(red0[j],green0[j],blue0[j],0.4))
}

# colors_border=c(rgb(0.2,0.5,0.5,0.9),rgb(0.8,0.2,0.5,0.9),rgb(0.7,0.5,0.2,0.9))
# colors_in=c(rgb(0.2,0.5,0.5,0.4),rgb(0.8,0.2,0.5,0.4),rgb(0.7,0.5,0.1,0.4))
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples")
# colnames(data0)[1:2]=c("Soluble_Solid_Content","Cis-3-Hexenyl_Acetate")
write.csv(data0,paste(target0[1],"&",target0[2],".Radar.Data.csv",sep=""),quote=F)

pdf(paste("Radar.Plot.",target0[1],"&",target0[2],".pdf" ,sep = ""), width = 20,height=10)

radarchart(data0,axistype=1,pcol=colors_border,pfcol=colors_in,plwd=4,plty=1,
           cglcol="grey",cglty=1,axislabcol="grey",
           cglwd=1,
           vlcex=1.5,title=paste(target0[1]," & ",target0[2],sep=""),cex.main=2)
legend(x=0.7,y=1.3,legend=rownames(data0[-c(1,2),]),
	   bty="n",lwd=2,lty=0,
	   pch=21,col=colors_border,text.col="grey",cex=1.5,pt.cex=1.5)

 dev.off()



}


# radarchart(data,axistype=1,pcol=colors_border,pfcol=colors_in,plwd=4,plty=1,
#            cglcol="grey",cglty=1,axisloabcol="grey",
#            caxislabels=seq(0,50,5),cglwd=0.8,
#            vlcex=0.8,title="kings.war")


#####  "2721214","2663753"
# targets=c("28509583","2721214","2663753")

# Y0=myY[,c(1,15,35,39,41,43,45)]
# cutoff=0.01
# n=nrow(myGM)
# taxa.g=as.character(myGD[,1])
# taxa.y=as.character(Y0[,1])
# yourGD=myGD[taxa.g%in%taxa.y,]
# yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]

target0=targets[2:3]
# X=yourGD[,-1]
index0=myGM[,3]%in%target0
tar.marker=X[,index0]
markers=paste(tar.marker[,1],tar.marker[,2],sep="&")
marker.names=paste(as.character(myGM[index0,1])[1],as.character(myGM[index0,1])[2],sep="&")
index0=markers=="0&0"
index1=markers=="0&1"
index2=markers=="0&2"
index3=markers=="1&1"
index4=markers=="1&2"
index5=markers=="2&0"
index6=markers=="2&1"
index7=markers=="2&2"
y0=yourY[index0,]
y1=yourY[index1,]
y2=yourY[index2,]
y3=yourY[index3,]
y4=yourY[index4,]
y5=yourY[index5,]
y6=yourY[index6,]
y7=yourY[index7,]

max.y=apply(yourY[,-1],2,function(one) max(one,na.rm=TRUE))
min.y=apply(yourY[,-1],2,function(one) min(one,na.rm=TRUE))

mean.y0=apply(y0[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y1=apply(y1[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y2=apply(y2[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y3=apply(y3[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y4=apply(y4[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y5=apply(y5[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y6=apply(y6[,-1],2,function(one) mean(one,na.rm=TRUE))
mean.y7=apply(y7[,-1],2,function(one) mean(one,na.rm=TRUE))


data0=rbind(max.y,min.y,mean.y0,mean.y1,mean.y2,mean.y3,mean.y4,mean.y5,mean.y6,mean.y7)

colnames(data0)=colnames(yourY)[-1]
rownames(data0)=c("1","2","Genotype:0&0","Genotype:0&1","Genotype:0&2","Genotype:1&1","Genotype:1&2","Genotype:2&0","Genotype:2&1","Genotype:2&2")
data0=as.data.frame(data0)[-c(8:10),]

colors_border=c(rgb(0.2,0.5,0.5,0.9),rgb(0.8,0.2,0.5,0.9),rgb(0.7,0.5,0.2,0.9),rgb(0.5,0.8,0.2,0.9),rgb(0.2,0.5,0.8,0.9))
colors_in=c(rgb(0.2,0.5,0.5,0.4),rgb(0.8,0.2,0.5,0.4),rgb(0.7,0.5,0.2,0.4),rgb(0.5,0.8,0.2,0.4),rgb(0.2,0.5,0.8,0.4))
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples")
write.csv(data0,paste(marker.names,".Radar.Data.csv",sep=""),quote=F)

radarchart(data0,axistype=1,pcol=colors_border,pfcol=colors_in,plwd=4,plty=1,
           cglcol="grey",cglty=1,axislabcol="grey",
           cglwd=0.8,
           vlcex=0.8,title=marker.names)
legend(x=0.7,y=1.3,legend=rownames(data0[-c(1,2,8:10),]),bty="n",lwd=2,lty=0,
	   pch=21,col=colors_border,text.col="grey",cex=0.8,pt.cex=1.5)



####### Trait 20

################# LD block


setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")

bin=25000
choose.chr0=c("G1","G3","G4","G8")
choose.pos0=c(18127098,4989391,1328075,13501777)
# aa=as.character(hapmap2[,13])
# index=grep("TRUE",aa)

# hapmap=as.data.frame(apply(hapmap,2,function(one) gsub("TRUE","T",one)))

# save.image("peah.hpm.RData")#145610
# hapmap=hapmap[-index,]
library(genetics)
library(LDheatmap)
GLD0=NULL
for(i in 1:length(choose.chr0))
{
	choose.chr=choose.chr0[i]
	choose.pos=choose.pos0[i]
myG6=hapmap2[hapmap2[,3]==choose.chr,]
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/result")

sigmarker=c(choose.pos)
range=c(bin,bin)
    chromosome=(myG6[,3]==choose.chr)
    bp=as.numeric(as.vector(myG6[,4]))
    deviation=bp-as.numeric(as.vector(sigmarker[1]))
    location1=deviation< as.numeric(as.vector(range[1])  )&deviation>=0
    location2=deviation*(-1)< as.numeric(as.vector(range[2])  )&deviation<=0
    location=location1|location2
    index1=chromosome&location


index=index1
# index=!is.na(myG6[,5])

GLD=myG6[index,]
GLD0=rbind(GLD0,GLD)
kg=GLD[,5]
kg2=kg[length(kg):1]
# write.table(kg2,paste("konwgene.posi.",choose.chr,".txt",sep=""),quote=F,row.names=F,col.names=F)
    hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
    LDdist=as.numeric(as.vector(GLD[,4]))
    LDsnpName=GLD[,1]
    snpname=GLD[,1]
    # LDsnpName=GLD[,5]
    colnames(hapmapgeno)=LDsnpName
    sigsnp=match(sigmarker,LDdist)
    
    LDsnpName=as.character(LDsnpName[c(sigsnp)]) #keep the first and last snp names only
    LDsnpName=gsub("SNP-","",LDsnpName)
    colnames(hapmapgeno)=gsub("SNP-","",colnames(hapmapgeno))
    color.rgb <- colorRampPalette(rev(c("snow","red")),space="rgb")
    LDsnp=makeGenotypes(hapmapgeno,sep="",method=as.genotype)   #This need to be converted to genotype object
    pdf(paste("Peach.",choose.chr,"(",round(LDdist[1]/1000000),"_",round(LDdist[length(LDdist)]/1000000),"Mb)",".pdf",sep=""), width = 12, height = 12)

    MyHeatmap <- LDheatmap(LDsnp, LDdist, flip=TRUE,
    color=color.rgb(20),SNP.name = LDsnpName, name="myLDgrob" ) 
    r2=MyHeatmap$LDmatrix
    write.csv(r2,paste("SNP-LD-pairwise-",choose.chr,".csv",sep=""),quote=F)
    dev.off()

}
colnames(GLD0)=as.character(as.matrix(hapmap2[1,]))
    write.csv(GLD0,paste("hapmap-LD.csv",sep=""),quote=F,row.names=F)



########### Only for rs172271 rs76610


setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")

bin=100000
choose.chr0=c("G2","G4")
choose.pos0=c(514845,1328075)
# aa=as.character(hapmap2[,13])
# index=grep("TRUE",aa)
hapmap2=hapmap
# hapmap=as.data.frame(apply(hapmap,2,function(one) gsub("TRUE","T",one)))

# save.image("peah.hpm.RData")#145610
# hapmap=hapmap[-index,]
library(genetics)
library(LDheatmap)
GLD0=NULL
for(i in 1:length(choose.chr0))
{
	choose.chr=choose.chr0[i]
	choose.pos=choose.pos0[i]
myG6=hapmap2[hapmap2[,3]==choose.chr,]
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/result")

sigmarker=c(choose.pos)
range=c(bin,bin)
    chromosome=(myG6[,3]==choose.chr)
    bp=as.numeric(as.vector(myG6[,4]))
    deviation=bp-as.numeric(as.vector(sigmarker[1]))
    location1=deviation< as.numeric(as.vector(range[1])  )&deviation>=0
    location2=deviation*(-1)< as.numeric(as.vector(range[2])  )&deviation<=0
    location=location1|location2
    index1=chromosome&location


index=index1
# index=!is.na(myG6[,5])

GLD=myG6[index,]
GLD0=rbind(GLD0,GLD)
kg=GLD[,5]
kg2=kg[length(kg):1]
# write.table(kg2,paste("konwgene.posi.",choose.chr,".txt",sep=""),quote=F,row.names=F,col.names=F)
    hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
    LDdist=as.numeric(as.vector(GLD[,4]))
    LDsnpName=GLD[,1]
    snpname=GLD[,1]
    # LDsnpName=GLD[,5]
    colnames(hapmapgeno)=LDsnpName
    sigsnp=match(sigmarker,LDdist)
    
    LDsnpName=as.character(LDsnpName[c(sigsnp)]) #keep the first and last snp names only
    LDsnpName=gsub("SNP-","",LDsnpName)
    colnames(hapmapgeno)=gsub("SNP-","",colnames(hapmapgeno))
    color.rgb <- colorRampPalette(rev(c("snow","red")),space="rgb")
    LDsnp=makeGenotypes(hapmapgeno,sep="",method=as.genotype)   #This need to be converted to genotype object
    pdf(paste("Peach.",choose.chr,"(",round(LDdist[1]/1000000),"_",round(LDdist[length(LDdist)]/1000000),"Mb)",".pdf",sep=""), width = 12, height = 12)

    MyHeatmap <- LDheatmap(LDsnp, LDdist, flip=TRUE,
    color=color.rgb(20),SNP.name = LDsnpName, name="myLDgrob" ) 
    r2=MyHeatmap$LDmatrix
    write.csv(r2,paste("SNP-LD-pairwise-",choose.chr,".csv",sep=""),quote=F)
    dev.off()

}

colnames(GLD0)=as.character(as.matrix(hapmap2[1,]))
    write.csv(GLD0,paste("hapmap-LD.csv",sep=""),quote=F,row.names=F)

######### filter hapmap

chr0=c(rep("G1",8),rep("G6",1),rep("G7",3))
pos0=c("25830591","25834100","25834421","25835032","25835644","25835944","25837345","25838015",
	"29760757","21872304","21872025","21870490")

hapmap2=hapmap[-1,]
hap=paste(hapmap2[,3],"_",hapmap2[,4])
hap0=paste(chr0,"_",pos0)
index=hap%in%hap0

hapmap3=hapmap2[index,]
colnames(hapmap3)=as.character(hapmap[1,])

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
write.table(hapmap3,paste("hapmap.filter.txt",sep=""),quote=F,row.names=F)


pos=as.character(hapmap2[,4])
index2=pos%in%pos0

chr=as.character(hapmap2[,3])
G1=hapmap2[hapmap2[,3]=="G1",]
######### genotype and phenotype distribution in T20

setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
myY=read.table("newalltraits2.txt",head=T) # 200
Y0=myY[,c(1,31)]
G0=t(hapmap[,-c(1:11)])
taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
colnames(G0)=taxa.marker
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")
taxa.g=as.character(G0[,1])
taxa.y=as.character(Y0[,1])
yourGD=G0[taxa.g%in%taxa.y,]
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
choose.chr0=c("G1","G3","G4","G8")
choose.pos0=c(18127098,4989391,1328075,13501777)


GLD0=NULL
for(i in 1:length(choose.chr0))
{
	choose.chr=choose.chr0[i]
	choose.pos=choose.pos0[i]
    myG6=hapmap[hapmap[,3]==choose.chr,]
    marker=myG6[myG6[,4]==choose.pos,]
    genotype0=as.character(marker)[1]
    geno=yourGD[,colnames(yourGD)%in%genotype0]
    geno[geno=="CT"]="TC"
    geno[geno=="AG"]="GA"

    geno.type=unique(geno)
    y00=NULL
    geno=cbind(as.data.frame(yourGD[,1]),geno)
    colnames(geno)=c("taxa","Genotype")
    y00=merge(geno,yourY,by.x="taxa",by.y="Taxa")
    for(j in 1:length(geno.type))
    {
    	index=geno==geno.type[j]
    	y=yourY[index,]
    	y[,1]=geno.type[j]
    	y00=rbind(y00,y)
    }
    y00=as.data.frame(y00[!is.na(y00[,2]),])
    colnames(y00)=c("Genotype","Phenotype")
    
    pdf(paste("Marker_",genotype0,"_Phenotype_distribution.pdf"), width = 10, height = 10)
    par(mar = c(5,5,5,5))
    boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
    	    main=paste(genotype0,"_",choose.chr,"_",choose.pos,sep=""))
    dev.off()

}

     # boxplot(count ~ spray, data = InsectSprays, col = "lightgray")


####### Trait 21

################# LD block


setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/β.ionone_2019")

bin=25000
choose.chr0=c("G3","G4")
choose.pos0=c(8136098,13803859)

library(genetics)
library(LDheatmap)
GLD0=NULL
for(i in 1:length(choose.chr0))
{
	choose.chr=choose.chr0[i]
	choose.pos=choose.pos0[i]
myG6=hapmap[hapmap[,3]==choose.chr,]
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/result")

sigmarker=c(choose.pos)
range=c(bin,bin)
    chromosome=(myG6[,3]==choose.chr)
    bp=as.numeric(as.vector(myG6[,4]))
    deviation=bp-as.numeric(as.vector(sigmarker[1]))
    location1=deviation< as.numeric(as.vector(range[1])  )&deviation>=0
    location2=deviation*(-1)< as.numeric(as.vector(range[2])  )&deviation<=0
    location=location1|location2
    index1=chromosome&location


index=index1
# index=!is.na(myG6[,5])

GLD=myG6[index,]
GLD0=rbind(GLD0,GLD)
kg=GLD[,5]
kg2=kg[length(kg):1]
# write.table(kg2,paste("konwgene.posi.",choose.chr,".txt",sep=""),quote=F,row.names=F,col.names=F)
    hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
    LDdist=as.numeric(as.vector(GLD[,4]))
    LDsnpName=GLD[,1]
    snpname=GLD[,1]
    # LDsnpName=GLD[,5]
    colnames(hapmapgeno)=LDsnpName
    sigsnp=match(sigmarker,LDdist)
    
    LDsnpName=as.character(LDsnpName[c(sigsnp)]) #keep the first and last snp names only
    LDsnpName=gsub("SNP-","",LDsnpName)
    colnames(hapmapgeno)=gsub("SNP-","",colnames(hapmapgeno))
    color.rgb <- colorRampPalette(rev(c("snow","red")),space="rgb")
    LDsnp=makeGenotypes(hapmapgeno,sep="",method=as.genotype)   #This need to be converted to genotype object
    pdf(paste("Peach.",choose.chr,"(",round(LDdist[1]/1000000),"_",round(LDdist[length(LDdist)]/1000000),"Mb)",".pdf",sep=""), width = 12, height = 12)

    MyHeatmap <- LDheatmap(LDsnp, LDdist, flip=TRUE,
    color=color.rgb(20),SNP.name = LDsnpName, name="myLDgrob" ) 
    r2=MyHeatmap$LDmatrix
    write.csv(r2,paste("SNP-LD-pairwise-",choose.chr,".csv",sep=""),quote=F)
    dev.off()

}
colnames(GLD0)=as.character(as.matrix(hapmap[1,]))
    write.csv(GLD0,paste("hapmap-LD.csv",sep=""),quote=F,row.names=F)

######### genotype and phenotype distribution in 3 traits 
## for normalization phenotype
## 33 β.ionone_2019
## 31 Linalool_2019
## 15 cis.3.hexenyl_acetate_2019

## for orignal phenotype
## 16 β.ionone_2019
## 15 Linalool_2019
## 9 cis.3.hexenyl_acetate_2019

setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)
myY=read.table("newalltraits2.txt",head=T) # 200
# myY=read.csv("oldalltraits.csv",head=T) # 200

Y0=myY[,c(1,9,15,16)]
G0=t(hapmap[,-c(1:11)])
taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
colnames(G0)=taxa.marker
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/β.ionone_2019")
taxa.g=as.character(G0[,1])
taxa.y=as.character(Y0[,1])
yourGD=G0[taxa.g%in%taxa.y,]
yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
choose.chr0=c("G3","G4")
choose.pos0=c(8136098,13803859)

GLD0=NULL
for(i in 1:length(choose.chr0))
{
	choose.chr=choose.chr0[i]
	choose.pos=choose.pos0[i]
    myG6=hapmap[hapmap[,3]==choose.chr,]
    marker=myG6[myG6[,4]==choose.pos,]
    genotype0=as.character(marker)[1]
    geno=yourGD[,colnames(yourGD)%in%genotype0]
    geno[geno=="CT"]="TC"
    geno[geno=="AG"]="GA"

    geno.type=unique(geno)
    y00=NULL
    for(j in 1:length(geno.type))
    {
    	index=geno==geno.type[j]
    	y=yourY[index,]
    	y[,1]=geno.type[j]
    	y00=rbind(y00,y)
    }
    y00=as.data.frame(y00[!is.na(y00[,2]),])
    colnames(y00)=c("Genotype","Phenotype")
    
    pdf(paste("Marker_",genotype0,"_Phenotype_distribution.pdf"), width = 10, height = 10)
    par(mar = c(5,5,5,5))
    boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
    	    main=paste(genotype0,"_",choose.chr,"_",choose.pos,sep=""))
    dev.off()

}

     # boxplot(count ~ spray, data = InsectSprays, col = "lightgray")



##################  PVE Proportion variance explained

#####


rm(list=ls())
# source("/home/jiabowang/Code/GAPIT.library.R")
source("/Users/Jiabo/Documents/gapit/gapit_functions.txt")
myY=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T) # 200

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang")
# Y1=newY[,c(2:6,8)]
myGM=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peachGM.txt",head=T) # 145k
index=myGM[,2]>8
myGM=myGM[!index,]
y.names=colnames(myY)[c(2:6,8,31)]
# setwd("/home/jiabowang/data/Xiongwei/traits-gwas")

	Y.names=y.names
	cutOff=0.05
	DPP=500
	WS=2e5
	GM=myGM


model_store=c("Blink","FarmCPU")

Nenviron=length(model_store)*length(Y.names)
  environ_name=NULL
  new_xz=NULL

  for(i in 1:length(model_store))
  {
    for(j in 1:length(Y.names))
    {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
      environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
    }
  }

sig_pos=NULL

themax.y0=NULL
store.x=NULL
y_filter0=NULL
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.table(paste("GAPIT.Filter_",environ_name[i],"_GWAS_result.txt",sep=""),head=T)
  environ_result=environ_result[order(environ_result[,3]),]
  environ_result=environ_result[order(environ_result[,2]),]
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  
  y_filter=environ_filter
  traits=environ_name[i]
  # print(head(y_filter))
  # print(traits)
  if(nrow(y_filter)>0)y_filter=cbind(y_filter[,1:5],traits)
  y_filter0=rbind(y_filter0,y_filter)
}
  write.table(y_filter0,paste("GAPIT.Filter_GWAS_results.txt",sep=""))




################# Calculate PVE for zhiliang and shuliang
rm(list=ls())
# source("/home/jiabowang/Code/GAPIT.library.R")
source("/home/jiabowang/Code/gapit_functions.txt")

setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
myCV=read.table("peachCV.txt",head=T)

myY0=read.table("/home/jiabowang/data/Xiongwei/newalltraits2.txt",head=T) # 200
setwd("/home/jiabowang/data/Xiongwei/traits-gwas")
data0="shuliang"
if(data0=="shuliang")
{
  result2=read.table("GAPIT.Filter_GWAS_results.txt",head=T)
  t21=read.csv("GAPIT.Blink.β.ionone_2019.GWAS.Results.csv",head=T)
  sig=t21[t21[,4]<(0.01/nrow(t21)),]
  result=sig[,1:5]
  traits=rep("Blink.β.ionone_2019",nrow(sig))
  sig2=cbind(result,traits)
  result2=rbind(result2,sig2)
}else{
result2=read.table("GAPIT.Filter_GWAS_results2.txt",head=T)
}
index=myGM[,2]>8
myGM=myGM[!index,]
myGD=myGD[,c(TRUE,!index)]
# setwd("/home/jiabowang/data/Xiongwei/rawdata")
# cv.index=read.table("population.txt",head=T)
# Y1=newY[,c(2:6,8,31)]   zhiliang
# Y2=newY[,c(15,31,33,35,39,41,43,45)]   shuliang
if(data0=="shuliang")
{myY=myY0[,c(1,15,31,33,35,39,41,43,45)]

}else{
myY=myY0[,c(1,2:6,8,31)]
}
# myCV=read.table("peachCV.txt",head=T)
setwd("/home/jiabowang/data/Xiongwei/traits-gwas")
traits=as.character(result2[,6])
traits2=gsub("Blink.","",gsub("FarmCPU.","",gsub("_2019","",traits)))
taxa.map=as.character(myGM[,1])
taxa.g=as.character(myGD[,1])
taxa.y=as.character(myY[,1])
taxa.cv=as.character(myCV[,1])
yourGD=myGD[taxa.g%in%taxa.y,]
yourY=myY[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
taxa.y2=as.character(yourY[,1])

yourCV=myCV[match(taxa.y2,taxa.cv)[!is.na(match(taxa.y2,taxa.cv))],]
colnames(yourY)=gsub("_2019","",colnames(yourY))
X=yourGD[,-1]


all.gene.list=NULL
for(i in 2:ncol(yourY))
{
    Y=yourY[,c(1,i)]
    name.of.trait=colnames(yourY)[i]
gwas=result2[traits2==colnames(yourY)[i],]
gwas2=gwas[!duplicated(gwas[,1]),]
   library("lme4")
    print("GAPIT.RandomModel beginning...")
    GT=as.character(yourY[,1])
    # name.of.trait=colnames(myY)[i]
    taxa.gwas2=as.character(gwas2[,1])
    index=taxa.map%in%taxa.gwas2   
    geneGD=X[,index,drop=FALSE]
    in_True=sum(index)

    if(in_True>0)
    {
    	colnames(geneGD)=paste("gene_",1:in_True,sep="")
    }else{
    	next
    }
    colnames(Y)=c("taxa","trait")
    CV=yourCV
    colnames(CV)=c("Taxa",paste("CV",1:(ncol(CV)-1),sep=""))
    tree2=cbind(Y,CV[,-1],geneGD)
    tree2=tree2[!is.na(tree2[,2]),]
    n_cv=ncol(yourCV)-1
    n_gd=in_True
    n_id=nrow(tree2)

       command0=paste("trait~1",sep="")
       command1=command0
       for(i in 1:n_cv)
       {	
	       command1=paste(command1,"+CV",i,sep="")
       }
           command2=command1
       for(j in 1:n_gd)
       {
    	   command2=paste(command2,"+(1|gene_",j,")",sep="")
       }
    dflme <- lme4::lmer(command2, data=tree2, control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
     check.nobs.vs.rankZ = "ignore",
     check.nobs.vs.nRE="ignore"))

    carcor_matrix=as.data.frame(summary(dflme)$varcor)
    var_gene=as.numeric(carcor_matrix[1:(nrow(carcor_matrix)-1),4])
    var_res=carcor_matrix[nrow(carcor_matrix),4]

    print(paste("Candidate Genes could Phenotype_Variance_Explained(%) :",sep=""))
    print(100*var_gene/sum(var_gene+var_res))
    v_rat=100*var_gene/sum(var_gene+var_res)
    gene_list=cbind(gwas2,v_rat)
    colnames(gene_list)[ncol(gene_list)]="Phenotype_Variance_Explained(%)"
    # utils::write.csv(var_gene,paste("GAPIT.", name.of.trait,".V_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
    # utils::write.csv(gene_list,paste("GAPIT.", name.of.trait,".PVE_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
#gene_list=read.csv("GAPIT.Weight.GrowthIntercept.Phenotype_Variance_Explained_by_Association_Markers.csv",head=T)
    all.gene.list=rbind(all.gene.list,gene_list)

}
    # utils::write.csv(all.gene.list,paste("GAPIT.peach.zhiliang.PVE_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
    # utils::write.csv(all.gene.list,paste("GAPIT.peach.shuliang.PVE_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
if(data0=="shuliang")
{    utils::write.csv(all.gene.list,paste("GAPIT.peach.shuliang.PVE_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}else{
    utils::write.csv(all.gene.list,paste("GAPIT.peach.zhiliang.PVE_by_Association_Markers.csv",sep=""),quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}



###################

######### LD  block for zhiliang traits
###########


rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures/T1-2-4-5-6-8")
# result0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)
data0="shuliang"
if(data0=="shuliang")
{
  result2=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas/GAPIT.Filter_GWAS_results.txt",head=T)
  t21=read.csv("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Blink.β.ionone_2019.GWAS.Results.csv",head=T)
  sig=t21[t21[,4]<(0.01/nrow(t21)),]
  result=sig[,1:5]
  traits=rep("Blink.β.ionone_2019",nrow(sig))
  sig2=cbind(result,traits)
  result2=rbind(result2,sig2)
}else{
  result2=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)

}
myY0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T) # 200

if(data0=="shuliang")
{myY=myY0[,c(1,15,33,35,39,41,43,45)]

}else{
myY=myY0[,c(1,2:6,8,31)]
}

# setwd("/home/jiabowang/data/Xiongwei/traits-gwas")
traits=as.character(result2[,6])
traits2=gsub("Blink.","",gsub("FarmCPU.","",gsub("_2019","",traits)))
# taxa.map=as.character(myGM[,1])
colnames(myY)=gsub("_2019","",colnames(myY))
bin=1000
library(genetics)
library(LDheatmap)

for(j in 2:ncol(myY))
{
    Y=myY[,c(1,j)]
    name.of.trait=colnames(myY)[j]
    sigs=result2[traits2==name.of.trait,]
    chr0=as.numeric(sigs[,2])
    pos0=as.numeric(sigs[,3])
    choose.chr0=paste("G",chr0,sep="")
    choose.pos0=pos0
    GLD0=NULL
    if(nrow(sigs)==0) next
    for(i in 1:length(choose.chr0))
       {
       choose.chr=choose.chr0[i]
       choose.pos=choose.pos0[i]
       myG6=hapmap[hapmap[,3]==choose.chr,]
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/result")
       sigmarker=c(choose.pos)
       range=c(bin,bin)
       chromosome=(myG6[,3]==choose.chr)
       bp=as.numeric(as.vector(myG6[,4]))
       deviation=bp-as.numeric(as.vector(sigmarker[1]))
       location1=deviation< as.numeric(as.vector(range[1])  )&deviation>=0
       location2=deviation*(-1)< as.numeric(as.vector(range[2])  )&deviation<=0
       location=location1|location2
       index1=chromosome&location
       index=index1
# index=!is.na(myG6[,5])

       GLD=myG6[index,]
       GLD0=rbind(GLD0,GLD)
       kg=GLD[,5]
       kg2=kg[length(kg):1]
# write.table(kg2,paste("konwgene.posi.",choose.chr,".txt",sep=""),quote=F,row.names=F,col.names=F)
       hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
       LDdist=as.numeric(as.vector(GLD[,4]))
       LDsnpName=GLD[,1]
       snpname=GLD[,1]
    # LDsnpName=GLD[,5]
       colnames(hapmapgeno)=LDsnpName
       sigsnp=match(sigmarker,LDdist)
    
       LDsnpName=as.character(LDsnpName[c(sigsnp)]) #keep the first and last snp names only
       LDsnpName=gsub("SNP-","",LDsnpName)
       colnames(hapmapgeno)=gsub("SNP-","",colnames(hapmapgeno))
       color.rgb <- colorRampPalette(rev(c("snow","red")),space="rgb")
       LDsnp=makeGenotypes(hapmapgeno,sep="",method=as.genotype)   #This need to be converted to genotype object
       
       pdf(paste("Peach.",choose.chr,"(",round(LDdist[1]/1000000),"_",round(LDdist[length(LDdist)]/1000000),"Mb)","-",name.of.trait,".pdf",sep=""), width = 12, height = 12)
    getNames()
    childNames(grid.get("myLDgrob"))
    childNames(grid.get("heatMap"))
    childNames(grid.get("geneMap"))
       MyHeatmap <- LDheatmap(LDsnp, LDdist, flip=TRUE,
                    color=color.rgb(20),SNP.name = LDsnpName, name="myLDgrob" ) 
    grid.edit(gPath("myLDgrob","heatMap","heatmap"),gp=gpar(col="white",lwd=1.5))
    grid.edit(gPath("myLDgrob", "geneMap","symbols"), gp = gpar(col="blue",cex=3))
    grid.edit(gPath("myLDgrob", "geneMap","SNPnames"), gp = gpar(col="blue",cex=1))

       r2=MyHeatmap$LDmatrix
       write.csv(r2,paste("SNP-LD-pairwise-",choose.chr,"-",name.of.trait,".csv",sep=""),quote=F)
       dev.off()
    library("grid")##添加分割线

       }
    colnames(GLD0)=as.character(as.matrix(hapmap[1,]))
    write.csv(GLD0,paste("hapmap-LD","-",name.of.trait,".csv",sep=""),quote=F,row.names=F)
}


######### genotype and phenotype distribution 

# myGD=data.table::fread("peachGD.txt",head=T) # 145K
# myGM=read.table("peachGM.txt",head=T) # 145k
# myCV0=read.table("peachCV.txt",head=T)


rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures/T1-2-4-5-6-8")
# result0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)
data0="shuliang"
if(data0=="shuliang")
{
  result2=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas/GAPIT.Filter_GWAS_results.txt",head=T)
  t21=read.csv("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Blink.β.ionone_2019.GWAS.Results.csv",head=T)
  sig=t21[t21[,4]<(0.01/nrow(t21)),]
  result=sig[,1:5]
  traits=rep("Blink.β.ionone_2019",nrow(sig))
  sig2=cbind(result,traits)
  result2=rbind(result2,sig2)
}else{
  result2=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)

}
myY0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T) # 200

if(data0=="shuliang")
{myY=myY0[,c(1,15,33,35,39,41,43,45)]

}else{
myY=myY0[,c(1,2:6,8,31)]
}

G0=t(hapmap[,-c(1:11)])
taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
taxa.g=as.character(G0[,1])
myY0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T) # 200

for(j in 2:ncol(myY))
{
    Y0=myY[,c(1,j)]
    taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
    colnames(G0)=taxa.marker
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")
    taxa.y=as.character(Y0[,1])
    yourGD=G0[taxa.g%in%taxa.y,]
    yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
    
    name.of.trait=colnames(myY)[j]
    sigs=result2[traits2==name.of.trait,]
    chr0=as.numeric(sigs[,2])
    pos0=as.numeric(sigs[,3])
    choose.chr0=paste("G",chr0,sep="")
    choose.pos0=pos0
    GLD0=NULL
    if(nrow(sigs)==0) next

    for(i in 1:length(choose.chr0))
    {       
    	choose.chr=choose.chr0[i]    	
    	choose.pos=choose.pos0[i]
    	myG6=hapmap[hapmap[,3]==choose.chr,]
    	marker=myG6[myG6[,4]==choose.pos,]
    	genotype0=as.character(marker)[1]
    	geno=yourGD[,colnames(yourGD)%in%genotype0]
    	print(table(geno))
    	geno[geno=="CT"]="TC"
    	geno[geno=="AG"]="GA"
    	geno[geno=="AC"]="CA"
    	geno[geno=="AT"]="TA"
    	geno[geno=="CG"]="GC"

    	geno.type=unique(geno)
    	y00=NULL
    	for(k in 1:length(geno.type))
    	{
    	   index=geno==geno.type[k]
    	   y=yourY[index,]
    	   y[,1]=geno.type[k]
    	   y00=rbind(y00,y)
        }
        y00=as.data.frame(y00[!is.na(y00[,2]),])
        colnames(y00)=c("Genotype","Phenotype")
    
        pdf(paste("Marker_",genotype0,"_",name.of.trait,"_Phenotype_distribution.pdf"), width = 10, height = 10)
        par(mar = c(5,5,5,5))
        boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
               main=paste(genotype0,"_",choose.chr,"_",choose.pos,sep=""))
        dev.off()

    }
}


############## phenotype distribution of cis-3-hex 

rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")

G0=t(hapmap[,-c(1:11)])
taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
taxa.g=as.character(G0[,1])
bin=100000
# choose.chr0=c("G2","G4")
# choose.pos0=c(514845,1328075)
colnames(myY0)=gsub("_2019","",colnames(myY0))
j=15
col.type=c("lightblue","mistyrose","lavender","lightgreen","lightgray","lightgoldenrod2","coral2","royalblue3")

Y0=myY0[,c(1,j)]
    taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
    colnames(G0)=taxa.marker
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")
    taxa.y=as.character(Y0[,1])
    yourGD=G0[taxa.g%in%taxa.y,]
    yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
    
    name.of.trait=colnames(Y0)[2]
    # sigs=result2[traits2==name.of.trait,]
    # chr0=as.numeric(sigs[,2])
    # pos0=as.numeric(sigs[,3])
    choose.chr0="G2"
    choose.pos0=514845
    GLD0=NULL
    # if(nrow(sigs)==0) next
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/cis.3.hexenyl_acetate_2019")
    for(i in 1:length(choose.chr0))
    {       
    	choose.chr=choose.chr0[i]    	
    	choose.pos=choose.pos0[i]
    	myG6=hapmap[hapmap[,3]==choose.chr,]
    	marker=myG6[myG6[,4]==choose.pos,]
    	genotype0=as.character(marker)[1]
    	geno=yourGD[,colnames(yourGD)%in%genotype0]
    	print(table(geno))
    	geno[geno=="CT"]="TC"
    	geno[geno=="AG"]="GA"
    	geno[geno=="AC"]="CA"
    	geno[geno=="AT"]="TA"
    	geno[geno=="CG"]="GC"

    	geno.type=unique(geno)
    	y00=NULL
    	for(k in 1:length(geno.type))
    	{
    	   index=geno==geno.type[k]
    	   y=yourY[index,]
    	   y[,1]=geno.type[k]
    	   y00=rbind(y00,y)
        }
        y00=as.data.frame(y00[!is.na(y00[,2]),])
        colnames(y00)=c("Genotype","Phenotype")
    
        pdf(paste("Marker_",genotype0,"_",name.of.trait,"_Phenotype_distribution.pdf"), width = 10, height = 10)
        par(mar = c(5,5,5,5))
        # boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
               # )
    xn=length(unique(y00[,1]))   
    aa=boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:3],plot=F)
    boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:3],axes=F)
    labels=aa$names
    posi=seq(1,xn,1)
    axis(2,col="black",col.ticks="black",col.axis="black",tck=-0.02,las=1)
    axis(1,at=posi,labels=labels,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    mtext(paste(name.of.trait," (",genotype0,")",sep=""),side=1,line=3,cex=1.2)
    mtext("Phenotype Values",side=2,line=3,cex=1.2)

        dev.off()

    }

############# s Manhattan plot with cis-3

# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
setwd("/home/jiabowang/data/Xiongwei/rawdata")
myGM=read.table("peachGM.txt",head=T)
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/cis.3.Hexen.1.ol_2019")
setwd("/home/jiabowang/data/Xiongwei/traits-gwas/complex8")
y.names="Linalool"
source("/home/jiabowang/Code/gapit_functions.txt")
source("/home/jiabowang/Code/GAPIT.Multiple.Manhattan.xiongwei.R")
# source("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/GAPIT.Multiple.Manhattan.xiongwei.R")
index=myGM[,2]>8
myGM=myGM[!index,]

    GMM=GAPIT.Multiple.Manhattan(model_store=c("MLM","MLMM","FarmCPU","Blink"),
	Y.names=y.names,cutOff=0.01,DPP=500,WS=3e5,inpch=c("+","*","-","#","<",">","^","$"),
	outpch=c(0,1,2,5),
	GM=myGM,plot.type=c("s"))






############## phenotype distribution of SSC


rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")

myY0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T)
j=8
G0=t(hapmap[,-c(1:11)])
taxa.g=as.character(G0[,1])
Y0=myY0[,c(1,j)]
    taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
    colnames(G0)=taxa.marker
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")
    taxa.y=as.character(Y0[,1])
    yourGD=G0[taxa.g%in%taxa.y,]
    yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
 col.type=c("lightblue","mistyrose","lavender","lightgreen","lightgray","lightgoldenrod2","coral2","royalblue3")
   
    name.of.trait=colnames(Y0)[2]
    name.of.trait="Soluble_Solid_Content"
    choose.chr0="G7"
    choose.pos0=3308817
    GLD0=NULL
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/Multiple Manhattan/symphysic.Manhattans")
results=read.csv("GAPIT.Association.Filter_GWAS_results_simple6.csv",head=T)
results0=results[grep("Total_soluble_solids",as.character(results[,7])),-1]
results0=results0[order(results0[,3]),]
results0=results0[order(results0[,2]),]

i=1
choose.chr=choose.chr0[i]    	
    	choose.pos=choose.pos0[i]
    	myG6=hapmap[hapmap[,3]==choose.chr,]
    	marker=myG6[myG6[,4]==choose.pos,]
    	genotype0=as.character(marker)[1]
    	geno=yourGD[,colnames(yourGD)%in%genotype0]
    	print(table(geno))
    	geno[geno=="CT"]="TC"
    	geno[geno=="AG"]="GA"
    	geno[geno=="AC"]="CA"
    	geno[geno=="AT"]="TA"
    	geno[geno=="CG"]="GC"

    	geno.type=unique(geno)
    	y00=NULL
    	for(k in 1:length(geno.type))
    	{
    	   index=geno==geno.type[k]
    	   y=yourY[index,]
    	   y[,1]=geno.type[k]
    	   y00=rbind(y00,y)
        }
        y00=as.data.frame(y00[!is.na(y00[,2]),])
        colnames(y00)=c("Genotype","Phenotype")

    par(mar = c(5,5,5,5))
        # boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
               # )
    xn=length(unique(y00[,1]))   
    aa=boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:3],plot=F)
    boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,
    	space=0.2,col=col.type[1:3],axes=F)
    labels=aa$names
    posi=seq(1,xn,1)
    axis(2,col="black",col.ticks="black",col.axis="black",tck=-0.02,las=1)
    axis(1,at=posi,labels=labels,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    mtext(paste(name.of.trait," (",genotype0,")",sep=""),side=1,line=3,cex=1.2)
    mtext("Phenotype Values",side=2,line=3,cex=1.2)
    


############## peach genotype view

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
myGD=read.table("peachGD.txt",head=T)
myGM=read.table("peachGM.txt",head=T)

source("/Users/Jiabo/Documents/Github/GAPIT3/R/GAPIT.Genotype.View.R")
index=nchar(myGM[,2])==2
yourGM=myGM[index,]
yourGD=myGD[,c(TRUE,index)]


source("/Users/Jiabo/Documents/Github/GAPIT3/R/GAPIT.Genotype.View.R")

ViewGenotype<-GAPIT.Genotype.View(
myGI=yourGM,
myGD=yourGD[,-1],
WS0 = 25000,
ws=200,
Aver.Dis=1000
)

X=myGD[,-1]
maf=apply(as.matrix(X),2,function(one) abs(1-sum(one)/(2*nrow(X))))
############# GS Plot

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Multiples")

acc=read.table("ACC.H2.txt",head=T)

acc0=acc[,-1]

acc1=t(as.matrix(acc0))
colnames(acc1)=as.character(acc[,1])

ACC.MT=acc1[1,]
ACC.ST=acc1[2,]
H2=acc1[3,]
ACC.MT=ACC.MT^2
ACC.ST=ACC.ST^2
acc1=rbind(ACC.ST,ACC.MT,H2)
colnames(acc1)=c("SH","FS","FC","PS","FT","SSC","β-Ionone","Cis3-Hex","γ-Hex","γ-Oct","γ-Dec","δ-Dec","γ-Dod","Linalool")

sample=acc1[,c(1:6)]
complex=acc1[,c(7:14)]
sample=sample[,order(sample[3,],decreasing=TRUE)]
complex=complex[,order(complex[3,],decreasing=TRUE)]
acc1=cbind(sample,complex)
# barplot(acc1,beside=TRUE)
    col.type=c("springgreen4","tomato3","steelblue3")
    # xn=ncol(acc1)
    # aa=barplot(acc1,beside=TRUE,xlab="",ylab="",
    # 	# space=0.2,
    # 	col=col.type,plot=F)
    par(mar = c(5,5,5,1))
    barplot(acc1,beside=TRUE,xlab="",ylab="",
            border=NA,ylim=c(0,1),las=1,cex.names=1.5,cex.axis=1.5,
            col=col.type,axes=T)
    legend("top",legend=c("Accuracy Square with Single-Traits","Accuracy Square with Multiple-Traits","Estimated Heritability"),col=col.type,
    # legend("top",legend=as.character("A","A","E"),col=col.type,
    	   pch=c(15,15,15),lwd=1,cex=2,lty=0,ncol=3,
    	   box.col="white",bty = "n", bg = par("bg"))
    mtext("Values",side=2,line=3,cex=1.5)
    # if(i==1|i==5)mtext("Phenotype Values",side=2,line=3,cex=1.2)


################ merge Manhatton Plots 

###### for Linalool merge Manhatton Plots by chromosome type

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
myGM=read.table("peachGM.txt",head=T)

choose.chr=4 
choose.pos=1328075
bin=100000

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")

mlm=read.csv("GAPIT.MLM.Linalool_2019.GWAS.Results.csv",head=T)
farm=read.csv("GAPIT.FarmCPU.Linalool_2019.GWAS.Results.csv",head=T)
blink=read.csv("GAPIT.Blink.Linalool_2019.GWAS.Results.csv",head=T)
mmlm=read.csv("GAPIT.MLMM.Linalool_2019.GWAS.Results.csv",head=T)
# bin=25000

	mlm0=mlm[mlm[,2]==choose.chr,]
    mlm1=mlm0[mlm0[,3]<choose.pos+bin&mlm0[,3]>choose.pos-bin,]
    mmlm0=mmlm[mmlm[,2]==choose.chr,]
    mmlm1=mmlm0[mmlm0[,3]<choose.pos+bin&mmlm0[,3]>choose.pos-bin,]
    farm0=farm[farm[,2]==choose.chr,]
    farm1=farm0[farm0[,3]<choose.pos+bin&farm0[,3]>choose.pos-bin,]
    blink0=blink[blink[,2]==choose.chr,]
    blink1=blink0[blink0[,3]<choose.pos+bin&blink0[,3]>choose.pos-bin,]
    mlm1=mlm1[order(mlm1[,3]),]
    mmlm1=mmlm1[order(mmlm1[,3]),]
    farm1=farm1[order(farm1[,3]),]
    blink1=blink1[order(blink1[,3]),]
    x=mlm1[,3]-min(mlm1[,3])+1
    y1=-log10(mlm1[,4])
    y2=-log10(mmlm1[,4])
    y3=-log10(farm1[,4])
    y4=-log10(blink1[,4])
    bonferroniCutOff=-log10(0.01/144614)
    maxy=round(max(y1,y2,y3,y4),0)+1
    labx=seq(1,max(x),length.out=5)
    pdf(paste("Peach.manhatton.",choose.chr,"(",bin,"_",bin,"bp)",".pdf",sep=""), width = 12, height = 7)
    plot(x,y1,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=1,col=1,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y2,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=2,col=2,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y3,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=3,col=3,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y4,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=4,col=4,lwd=1.5,cex=1.5)
    abline(h=bonferroniCutOff,col="gray")

axis(1,at=labx,labels=round((labx+min(mlm1[,3])-1)/1000000,3),col="black",col.ticks="black",tck=-0.03,col.axis="black",yaxp=c(1,max(x),5))
axis(2,col="black",col.ticks="black",tck=-0.03,col.axis="black",yaxp=c(0,maxy,5))
mtext(paste("Chr:",choose.chr," (Mb)",sep=""),side=1,line=2.5,cex=1,col="black")
mtext(expression(-log[10](italic(p))),side=2,line=2.5,cex=1,col="black")

	legend("topright",legend=c("MLM","MLMM","FarmCPU","Blink"),
col=c(1:4),pch=1:4,lty=0,lwd=1.5,cex=1.5,
 bty = "n", bg = par("bg"))

 dev.off()

###### for cis.3.hex merge Manhatton Plots by choromesome type


setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
myGM=read.table("peachGM.txt",head=T)

choose.chr=2
choose.pos=514845
bin=100000

setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/cis.3.hexenyl_acetate_2019")

mlm=read.csv("GAPIT.MLM.cis.3.hexenyl_acetate_2019.GWAS.Results.csv",head=T)
farm=read.csv("GAPIT.FarmCPU.cis.3.hexenyl_acetate_2019.GWAS.Results.csv",head=T)
blink=read.csv("GAPIT.Blink.cis.3.hexenyl_acetate_2019.GWAS.Results.csv",head=T)
mmlm=read.csv("GAPIT.MLMM.cis.3.hexenyl_acetate_2019.GWAS.Results.csv",head=T)
# bin=25000

	mlm0=mlm[mlm[,2]==choose.chr,]
    mlm1=mlm0[mlm0[,3]<choose.pos+bin&mlm0[,3]>choose.pos-bin,]
    mmlm0=mmlm[mmlm[,2]==choose.chr,]
    mmlm1=mmlm0[mmlm0[,3]<choose.pos+bin&mmlm0[,3]>choose.pos-bin,]
    farm0=farm[farm[,2]==choose.chr,]
    farm1=farm0[farm0[,3]<choose.pos+bin&farm0[,3]>choose.pos-bin,]
    blink0=blink[blink[,2]==choose.chr,]
    blink1=blink0[blink0[,3]<choose.pos+bin&blink0[,3]>choose.pos-bin,]
    mlm1=mlm1[order(mlm1[,3]),]
    mmlm1=mmlm1[order(mmlm1[,3]),]
    farm1=farm1[order(farm1[,3]),]
    blink1=blink1[order(blink1[,3]),]
    x=mlm1[,3]-min(mlm1[,3])+1
    y1=-log10(mlm1[,4])
    y2=-log10(mmlm1[,4])
    y3=-log10(farm1[,4])
    y4=-log10(blink1[,4])

    maxy=round(max(y1,y2,y3,y4),0)+1
    labx=seq(1,max(x),length.out=5)
    pdf(paste("Peach.manhatton.",choose.chr,"(",bin,"_",bin,"bp)",".pdf",sep=""), width = 12, height = 7)
    plot(x,y1,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=1,col=1,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y2,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=2,col=2,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y3,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=3,col=3,lwd=1.5,cex=1.5)
    par(new=T)
    plot(x,y4,xlab="",ylab="",tck=-0.01,axes=F,
  ylim=c(0,maxy),xlim=c(1,max(x)),pch=4,col=4,lwd=1.5,cex=1.5)

axis(1,at=labx,labels=round((labx+min(mlm1[,3])-1)/1000000,3),col="black",col.ticks="black",tck=-0.03,col.axis="black",yaxp=c(1,max(x),5))
axis(2,col="black",col.ticks="black",tck=-0.03,col.axis="black",yaxp=c(0,maxy,5))
mtext(paste("Chr:",choose.chr," (Mb)",sep=""),side=1,line=2.5,cex=1,col="black")
mtext(expression(-log[10](italic(p))),side=2,line=2.5,cex=1,col="black")

	legend("topright",legend=c("MLM","MLMM","FarmCPU","Blink"),
col=c(1:4),pch=1:4,lty=0,lwd=1.5,cex=1.5,
 bty = "n", bg = par("bg"))

 dev.off()



# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
setwd("/home/jiabowang/data/Xiongwei/rawdata")
# myGD=read.table("peachGD.txt",head=T) # 145K
myGM=read.table("peachGM.txt",head=T) # 145k
setwd("/home/jiabowang/data/Xiongwei/Multiple-traits/cis.3.hexenyl_acetate_2019")

# source("/Users/Jiabo/Documents/gapit/gapit_functions.txt")
source("/home/jiabowang/Code/gapit_functions.txt")
index=myGM[,2]>8
myGM=myGM[!index,]
source("/home/jiabowang/Code/GAPIT.Multiple.Manhattan.xiongwei.R")

GMM=GAPIT.Multiple.Manhattan(model_store=c("MLM","MMLM","FarmCPU","Blink"),
	Y.names="cis.3.hexenyl_acetate_2019",plot.type=c("s"),outpch=c(0,1,2,5),
	cutOff=0.01,GM=myGM)

################# filtering gene

rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures/T1-2-4-5-6-8")
# result0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)

chr0="G6"
position=105917 

index=hapmap[,3]==chr0&hapmap[,4]==position

hapmap1=rbind(hapmap[1,],hapmap[index,])

write.csv(hapmap1,paste("hapmap.within.",chr0,".",position,".csv",sep=""),quote=FALSE,row.names=FALSE)


rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
# myG=read.table("peach.hapmap.txt",head=F)
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Figures/T1-2-4-5-6-8")
# result0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/traits-gwas-zhliang/GAPIT.Filter_GWAS_results2.txt",head=T)

chr0="G6"
position0=1985307 
position1=2004003 

chr3=c("G4","G8","G7")
position3=c(9142621:9342621,17881992:18081992,12947439:13147439)


index=hapmap[,3]==chr0&as.numeric(hapmap[,4])>position0&as.numeric(hapmap[,4])<position1

index2=hapmap[,3]%in%chr3&as.numeric(hapmap[,4])%in%position3

hapmap1=rbind(hapmap[1,],hapmap[index,],hapmap[index2,])

write.csv(hapmap1,paste("hapmap.within.",chr0,".",position0,".",position1,".csv",sep=""),quote=FALSE,row.names=FALSE)






#####################

setwd("/Users/Jiabo/Documents/Data/Xiongwei")

pheno=read.table("2observe14VOC.txt",head=T)
pheno2=pheno

for(i in c(2,3))
{
	store=unique(pheno2[,i])
	for(j in 1:length(store))
	{
		pheno2[pheno2[,i]==store[j],i]=j-1
	}
	
}

obs=pheno2[,2:3]
voc=pheno2[,c(4,6:18)]
all.rr=NULL
for(i in 1:ncol(voc))
{
	test=cbind(as.numeric(obs[,1]),as.numeric(voc[,i]))
	test=test[!is.na(apply(test,1,sum)),]
	rr=cor(test[,1],test[,2])
	print(paste(colnames(voc)[i],":",rr,sep="	"))
	names(rr)=paste(colnames(voc)[i],"&",colnames(obs)[1],sep="")
	all.rr=append(all.rr,rr)
}
# cor(obs[,1],voc)
write.csv(all.rr,"correlations.csv",quote=F)








######### genotype and phenotype distribution in 3 traits 
## for normalization phenotype
## 33 β.ionone_2019
## 31 Linalool_2019
## 15 cis.3.hexenyl_acetate_2019

## for orignal phenotype
## 15 Linalool
## 9 cis.3.hexenyl_acetate
## 7 Total_soluble_solids
rm(list=ls())
setwd("/Users/Jiabo/Documents/Data/Xiongwei")
load("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/peah.hpm.RData")
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset")
# myY0=read.table("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/newalltraits2.txt",head=T)
myY0=read.csv("/Users/Jiabo/Documents/Data/Xiongwei/Traits/Dataset/oldalltraits.csv",head=T)

for(k in c(7,9,15))
{
setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/Multiple Manhattan/symphysic.Manhattans")
if(k==7)
{
results=read.csv("GAPIT.Association.Filter_GWAS_results_simple6.csv",head=T)
}else{
results=read.csv("GAPIT.Association.Filter_GWAS_results_complex8.csv",head=T)	
}
trait.name=colnames(myY0)[k]
results0=results[grep(trait.name,as.character(results[,7])),-1]
results0=results0[order(results0[,3]),]
results0=results0[order(results0[,2]),]

G0=t(hapmap[,-c(1:11)])
taxa.g=as.character(G0[,1])
Y0=myY0[,c(1,k)]
    taxa.marker=as.character(c("Taxa",as.character(hapmap[-1,1])))
    colnames(G0)=taxa.marker
# setwd("/Users/Jiabo/Documents/Data/Xiongwei/Multiple-traits/zhiliang/Linalool_2019")
    taxa.y=as.character(Y0[,1])
    yourGD=G0[taxa.g%in%taxa.y,]
    yourY=Y0[match(taxa.g,taxa.y)[!is.na(match(taxa.g,taxa.y))],]
 col.type=c("lightblue","mistyrose","lavender","lightgreen","lightgray","lightgoldenrod2","coral2","royalblue3")
   
    # name.of.trait=colnames(Y0)[2]
    # name.of.trait="Soluble_Solid_Content"
if(k==15)
{
choose.chr0="G4"
choose.pos0=1328075
}
if(k==9)
{
choose.chr0="G2"
choose.pos0=514845
}
if(k==7)
{
choose.chr0="G7"
choose.pos0=3308817
}

i=1
choose.chr=choose.chr0[i]    	
    	choose.pos=choose.pos0[i]
    	myG6=hapmap[hapmap[,3]==choose.chr,]
    	marker=myG6[myG6[,4]==choose.pos,]
    	genotype0=as.character(marker)[1]
    	geno=yourGD[,colnames(yourGD)%in%genotype0]
    	print(table(geno))
    	geno[geno=="CT"]="TC"
    	geno[geno=="AG"]="GA"
    	geno[geno=="AC"]="CA"
    	geno[geno=="AT"]="TA"
    	geno[geno=="CG"]="GC"

    	geno.type=unique(geno)
    	y00=NULL
    	for(k in 1:length(geno.type))
    	{
    	   index=geno==geno.type[k]
    	   y=yourY[index,]
    	   y[,1]=geno.type[k]
    	   y00=rbind(y00,y)
        }
        y00=as.data.frame(y00[!is.na(y00[,2]),])
        colnames(y00)=c("Genotype","Phenotype")
    type.num=length(unique(y00[,1]))
    marker.type=as.character(unique(y00[,1]))
    marker.type=marker.type[order(marker.type)]

    pdf(paste("Distribution of phenotype in ",trait.name,".pdf" ,sep = ""), width = 10,height=10)

    par(mar = c(5,5,5,5))
        # boxplot(Phenotype ~ Genotype, data = y00,col=c(1:length(geno.type)+1),las=1,
               # )
    xn=length(unique(y00[,1]))

    aa=boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,ylim=c(floor(min(y00[,2],na.rm=T)),ceiling(max(y00[,2],na.rm=T))),
    	space=0.2,col=col.type[1:3],plot=F)
    boxplot(Phenotype~Genotype,data=y00,xlab="",ylab="",las=1,ylim=c(floor(min(y00[,2],na.rm=T)),ceiling(max(y00[,2],na.rm=T))),
    	space=0.2,col=col.type[1:3],axes=F,outline=FALSE)
    for(j in 1:type.num)
    {
    yj=y00[y00[,1]==marker.type[j],2]
    points((j+runif(length(yj),min=-0.2,max=0.2) ), yj, cex=0.7,pch = 1,  col="blue")
    }# end of j
    labels=aa$names
    posi=seq(1,xn,1)
    axis(2,col="black",col.ticks="black",col.axis="black",tck=-0.02,las=1)
    axis(1,at=posi,labels=labels,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
    mtext(paste(trait.name," (",genotype0," - Chr:",choose.chr," Pos:",choose.pos,")",sep=""),side=1,line=3,cex=1.2)
    mtext("Phenotype Values",side=2,line=3,cex=1.2)
    dev.off()
}



