### R code from vignette source 'EBSeq_Vignette.Rnw'

###################################################
### code chunk number 1: EBSeq_Vignette.Rnw:172-173
###################################################
library(EBSeq)


###################################################
### code chunk number 2: EBSeq_Vignette.Rnw:198-200
###################################################
data(GeneMat)
str(GeneMat)


###################################################
### code chunk number 3: EBSeq_Vignette.Rnw:208-209
###################################################
Sizes=MedianNorm(GeneMat)


###################################################
### code chunk number 4: EBSeq_Vignette.Rnw:235-237
###################################################
EBOut=EBTest(Data=GeneMat, 
Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)


###################################################
### code chunk number 5: EBSeq_Vignette.Rnw:241-244
###################################################
PP=GetPPMat(EBOut)
str(PP)
head(PP)


###################################################
### code chunk number 6: EBSeq_Vignette.Rnw:249-251
###################################################
DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]
str(DEfound)


###################################################
### code chunk number 7: EBSeq_Vignette.Rnw:285-291
###################################################
data(IsoList)
str(IsoList)
IsoMat=IsoList$IsoMat
str(IsoMat)
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames


###################################################
### code chunk number 8: EBSeq_Vignette.Rnw:298-299
###################################################
IsoSizes=MedianNorm(IsoMat)


###################################################
### code chunk number 9: EBSeq_Vignette.Rnw:320-323
###################################################
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]


###################################################
### code chunk number 10: EBSeq_Vignette.Rnw:335-342
###################################################
IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
IsoPP=GetPPMat(IsoEBOut)
str(IsoPP)
head(IsoPP)
IsoDE=rownames(IsoPP)[which(IsoPP[,"PPDE"]>=.95)]
str(IsoDE)


###################################################
### code chunk number 11: EBSeq_Vignette.Rnw:353-355
###################################################
data(MultiGeneMat)
str(MultiGeneMat)


###################################################
### code chunk number 12: EBSeq_Vignette.Rnw:363-366
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti


###################################################
### code chunk number 13: EBSeq_Vignette.Rnw:374-376
###################################################
Parti=PosParti[-3,]
Parti


###################################################
### code chunk number 14: EBSeq_Vignette.Rnw:381-384
###################################################
MultiSize=MedianNorm(MultiGeneMat)
MultiOut=EBMultiTest(MultiGeneMat,NgVector=NULL,Conditions=Conditions,
AllParti=Parti, sizeFactors=MultiSize, maxround=5)


###################################################
### code chunk number 15: EBSeq_Vignette.Rnw:388-393
###################################################
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns


###################################################
### code chunk number 16: EBSeq_Vignette.Rnw:412-420
###################################################
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiSize=MedianNorm(IsoMultiMat)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")


###################################################
### code chunk number 17: EBSeq_Vignette.Rnw:426-428
###################################################
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond


###################################################
### code chunk number 18: EBSeq_Vignette.Rnw:433-435
###################################################
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond


###################################################
### code chunk number 19: EBSeq_Vignette.Rnw:440-444
###################################################
IsoMultiOut=EBMultiTest(IsoMultiMat,
NgVector=IsoNgTrun.Multi,Conditions=Conditions,
AllParti=Parti.4Cond, sizeFactors=IsoMultiSize, 
maxround=5)


###################################################
### code chunk number 20: EBSeq_Vignette.Rnw:448-453
###################################################
IsoMultiPP=GetMultiPP(IsoMultiOut)
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns


###################################################
### code chunk number 21: EBSeq_Vignette.Rnw:470-475 (eval = FALSE)
###################################################
## data(GeneMat)
## Sizes=MedianNorm(GeneMat)
## EBOut=EBTest(Data=GeneMat, 
## Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=Sizes, maxround=5)
## PP=GetPPMat(EBOut)


###################################################
### code chunk number 22: EBSeq_Vignette.Rnw:477-481
###################################################
str(PP)
head(PP)
DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]
str(DEfound)


###################################################
### code chunk number 23: EBSeq_Vignette.Rnw:491-494
###################################################
GeneFC=PostFC(EBOut)
str(GeneFC)
PlotPostVsRawFC(EBOut,GeneFC)


###################################################
### code chunk number 24: EBSeq_Vignette.Rnw:515-518
###################################################
EBOut$Alpha
EBOut$Beta
EBOut$P


###################################################
### code chunk number 25: EBSeq_Vignette.Rnw:537-539
###################################################
par(mfrow=c(1,2))
QQP(EBOut)


###################################################
### code chunk number 26: EBSeq_Vignette.Rnw:555-557
###################################################
par(mfrow=c(1,2))
DenNHist(EBOut)


###################################################
### code chunk number 27: EBSeq_Vignette.Rnw:578-583 (eval = FALSE)
###################################################
## data(IsoList)
## IsoMat=IsoList$IsoMat
## IsoNames=IsoList$IsoNames
## IsosGeneNames=IsoList$IsosGeneNames
## NgList=GetNg(IsoNames, IsosGeneNames, TrunThre=3)


###################################################
### code chunk number 28: EBSeq_Vignette.Rnw:585-588
###################################################
names(NgList)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,201:203,601:603)]


###################################################
### code chunk number 29: EBSeq_Vignette.Rnw:619-620 (eval = FALSE)
###################################################
## IsoNgTrun = scan(file="output_name.ngvec", what=0, sep="\n")


###################################################
### code chunk number 30: EBSeq_Vignette.Rnw:633-638 (eval = FALSE)
###################################################
## IsoSizes=MedianNorm(IsoMat)
## IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
## Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
## IsoPP=GetPPMat(IsoEBOut)
## IsoDE=rownames(IsoPP)[which(IsoPP[,"PPDE"]>=.95)]


###################################################
### code chunk number 31: EBSeq_Vignette.Rnw:640-641
###################################################
str(IsoDE)


###################################################
### code chunk number 32: EBSeq_Vignette.Rnw:646-648
###################################################
IsoFC=PostFC(IsoEBOut)
str(IsoFC)


###################################################
### code chunk number 33: EBSeq_Vignette.Rnw:659-662
###################################################
IsoEBOut$Alpha
IsoEBOut$Beta
IsoEBOut$P


###################################################
### code chunk number 34: EBSeq_Vignette.Rnw:681-686
###################################################
par(mfrow=c(2,2))
PolyFitValue=vector("list",3)
for(i in 1:3)
    PolyFitValue[[i]]=PolyFitPlot(IsoEBOut$C1Mean[[i]], 
    IsoEBOut$C1EstVar[[i]],5)


###################################################
### code chunk number 35: EBSeq_Vignette.Rnw:699-708
###################################################
PolyAll=PolyFitPlot(unlist(IsoEBOut$C1Mean), unlist(IsoEBOut$C1EstVar),5)
lines(log10(IsoEBOut$C1Mean[[1]][PolyFitValue[[1]]$sort]), 
PolyFitValue[[1]]$fit[PolyFitValue[[1]]$sort],col="yellow",lwd=2)
lines(log10(IsoEBOut$C1Mean[[2]][PolyFitValue[[2]]$sort]), 
PolyFitValue[[2]]$fit[PolyFitValue[[2]]$sort],col="pink",lwd=2)
lines(log10(IsoEBOut$C1Mean[[3]][PolyFitValue[[3]]$sort]), 
PolyFitValue[[3]]$fit[PolyFitValue[[3]]$sort],col="green",lwd=2)
legend("topleft",c("All Isoforms","Ng = 1","Ng = 2","Ng = 3"),
col=c("red","yellow","pink","green"),lty=1,lwd=3,box.lwd=2)


###################################################
### code chunk number 36: EBSeq_Vignette.Rnw:721-723
###################################################
par(mfrow=c(2,3))
QQP(IsoEBOut)


###################################################
### code chunk number 37: EBSeq_Vignette.Rnw:735-737
###################################################
par(mfrow=c(2,3))
DenNHist(IsoEBOut)


###################################################
### code chunk number 38: EBSeq_Vignette.Rnw:754-758
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti
PlotPattern(PosParti)


###################################################
### code chunk number 39: EBSeq_Vignette.Rnw:765-767
###################################################
Parti=PosParti[-3,]
Parti


###################################################
### code chunk number 40: EBSeq_Vignette.Rnw:773-779 (eval = FALSE)
###################################################
## data(MultiGeneMat)
## MultiSize=MedianNorm(MultiGeneMat)
## MultiOut=EBMultiTest(MultiGeneMat,
## NgVector=NULL,Conditions=Conditions,
## AllParti=Parti, sizeFactors=MultiSize, 
## maxround=5)


###################################################
### code chunk number 41: EBSeq_Vignette.Rnw:783-788
###################################################
MultiPP=GetMultiPP(MultiOut)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns


###################################################
### code chunk number 42: EBSeq_Vignette.Rnw:795-797
###################################################
MultiFC=GetMultiFC(MultiOut)
str(MultiFC)


###################################################
### code chunk number 43: EBSeq_Vignette.Rnw:806-808
###################################################
par(mfrow=c(2,2))
QQP(MultiOut)


###################################################
### code chunk number 44: EBSeq_Vignette.Rnw:816-818
###################################################
par(mfrow=c(2,2))
DenNHist(MultiOut)


###################################################
### code chunk number 45: EBSeq_Vignette.Rnw:833-836
###################################################
Conditions=c("C1","C1","C2","C2","C3","C3","C4","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond


###################################################
### code chunk number 46: EBSeq_Vignette.Rnw:841-844
###################################################
PlotPattern(PosParti.4Cond)
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond


###################################################
### code chunk number 47: EBSeq_Vignette.Rnw:851-862 (eval = FALSE)
###################################################
## data(IsoMultiList)
## IsoMultiMat=IsoMultiList[[1]]
## IsoNames.Multi=IsoMultiList$IsoNames
## IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
## IsoMultiSize=MedianNorm(IsoMultiMat)
## NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
## IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
## IsoMultiOut=EBMultiTest(IsoMultiMat,NgVector=IsoNgTrun.Multi,Conditions=Conditions,
## AllParti=Parti.4Cond, 
## sizeFactors=IsoMultiSize, maxround=5)
## IsoMultiPP=GetMultiPP(IsoMultiOut)


###################################################
### code chunk number 48: EBSeq_Vignette.Rnw:864-869
###################################################
names(MultiPP)
IsoMultiPP$PP[1:10,]
IsoMultiPP$MAP[1:10]
IsoMultiPP$Patterns
IsoMultiFC=GetMultiFC(IsoMultiOut)


###################################################
### code chunk number 49: EBSeq_Vignette.Rnw:880-883
###################################################
par(mfrow=c(3,4))
QQP(IsoMultiOut)



###################################################
### code chunk number 50: EBSeq_Vignette.Rnw:893-895
###################################################
par(mfrow=c(3,4))
DenNHist(IsoMultiOut)


###################################################
### code chunk number 51: EBSeq_Vignette.Rnw:927-936
###################################################
data(GeneMat)
GeneMat.norep=GeneMat[,c(1,6)]
Sizes.norep=MedianNorm(GeneMat.norep)
EBOut.norep=EBTest(Data=GeneMat.norep,
Conditions=as.factor(rep(c("C1","C2"))),
sizeFactors=Sizes.norep, maxround=5)
PP.norep=GetPPMat(EBOut.norep)
DEfound.norep=rownames(PP.norep)[which(PP.norep[,"PPDE"]>=.95)]
GeneFC.norep=PostFC(EBOut.norep)


###################################################
### code chunk number 52: EBSeq_Vignette.Rnw:946-960
###################################################
data(IsoList)
IsoMat=IsoList$IsoMat
IsoNames=IsoList$IsoNames
IsosGeneNames=IsoList$IsosGeneNames
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoMat.norep=IsoMat[,c(1,6)]
IsoSizes.norep=MedianNorm(IsoMat.norep)
IsoEBOut.norep=EBTest(Data=IsoMat.norep, NgVector=IsoNgTrun,
Conditions=as.factor(c("C1","C2")),
sizeFactors=IsoSizes.norep, maxround=5)
IsoPP.norep=GetPPMat(IsoEBOut.norep)
IsoDE.norep=rownames(IsoPP.norep)[which(IsoPP.norep[,"PPDE"]>=.95)]
IsoFC.norep=PostFC(IsoEBOut.norep)


###################################################
### code chunk number 53: EBSeq_Vignette.Rnw:969-981
###################################################
data(MultiGeneMat)
MultiGeneMat.norep=MultiGeneMat[,c(1,3,5)]
Conditions=c("C1","C2","C3")
PosParti=GetPatterns(Conditions)
Parti=PosParti[-3,]
MultiSize.norep=MedianNorm(MultiGeneMat.norep)
MultiOut.norep=EBMultiTest(MultiGeneMat.norep,
NgVector=NULL,Conditions=Conditions,
AllParti=Parti, sizeFactors=MultiSize.norep, 
maxround=5)
MultiPP.norep=GetMultiPP(MultiOut.norep)
MultiFC.norep=GetMultiFC(MultiOut.norep)


###################################################
### code chunk number 54: EBSeq_Vignette.Rnw:993-1012
###################################################
data(IsoMultiList)
IsoMultiMat=IsoMultiList[[1]]
IsoNames.Multi=IsoMultiList$IsoNames
IsosGeneNames.Multi=IsoMultiList$IsosGeneNames
IsoMultiMat.norep=IsoMultiMat[,c(1,3,5,7)]
IsoMultiSize.norep=MedianNorm(IsoMultiMat.norep)
NgList.Multi=GetNg(IsoNames.Multi, IsosGeneNames.Multi)
IsoNgTrun.Multi=NgList.Multi$IsoformNgTrun
Conditions=c("C1","C2","C3","C4")
PosParti.4Cond=GetPatterns(Conditions)
PosParti.4Cond
Parti.4Cond=PosParti.4Cond[c(1,2,3,8,15),]
Parti.4Cond
IsoMultiOut.norep=EBMultiTest(IsoMultiMat.norep,
NgVector=IsoNgTrun.Multi,Conditions=Conditions,
AllParti=Parti.4Cond, sizeFactors=IsoMultiSize.norep, 
maxround=5)
IsoMultiPP.norep=GetMultiPP(IsoMultiOut.norep)
IsoMultiFC.norep=GetMultiFC(IsoMultiOut.norep)


