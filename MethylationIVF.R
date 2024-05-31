# Loading packages used in the analysis

library(minfi)
library(minfiData)
library(sva)
library(devtools)
library(missMethyl)
library(limma)
library(minfi)
library(maxprobes)
library(ggplot2)
library(DMRcate)
library(ChAMP)
library(REMP)
library(qqman)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Gviz)
library(GenomicRanges)
library(factoextra)
library(FactoMineR)
library(readxl)
library(RColorBrewer)
library(gtools)
library(stringr)
library(EpiDISH)

#set.seed for reproducibility (because normalization can add a small variation in M-values)
set.seed(100)

# Setting directory for idat files and design matrix
targets = read.metharray.sheet( base="your_directory" , pattern = "csv$" , ignore.case = T , verbose = T )

# load data from samplesheet info
rgSet <- read.metharray.exp(targets = targets)
phenoData <- pData(rgSet)

# normalization using SWAN method
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=TRUE)

# filter poor quality probes (detection pvalue<0.01)
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSwQC <- mSetSw[keep,]

# filter snp probes 
t = mapToGenome( mSetSwQC )
mSetSwQCnosnp = dropLociWithSnps( t )  

# cross-reactive probes (probes that map at different loci)
mSetFinal = maxprobes::dropXreactiveLoci( mSetSwQCnosnp )

# number of probes that passed filtering
dim(mSet)[1]
dim(mSetFinal)[1]

beta <- getBeta(mSetFinal)

#Deconvoluting cell heterogeneity

data(centEpiFibIC.m)
data(centBloodSub.m)
frac.m <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
frac.m
boxplot(frac.m)

pheno.v=c(rep("ART",36),rep("CTRL",12))

a=data.frame(frac.m)
par(mfrow=c(2,3))
plot(a$Epi,a$Neutro,xlab="Epithelial cell fraction",ylab="Neutrophil fraction")
legend(x='bottomleft', legend=paste('Cor =',round(cor(a$Epi, a$Neutro),2)))
plot(a$Epi,a$B,xlab="Epithelial cell fraction",ylab="B cell fraction")
legend(x='bottomleft', legend=paste('Cor =',round(cor(a$Epi, a$B),2)))
plot(a$Epi,a$CD4T,xlab="Epithelial cell fraction",ylab="CD4T cell fraction")
legend(x='bottomleft', legend=paste('Cor =',round(cor(a$Epi, a$CD4T),2)))
plot(a$Epi,a$Mono,xlab="Epithelial cell fraction",ylab="Monocyte fraction")
legend(x='bottomleft', legend=paste('Cor =',round(cor(a$Epi, a$Mono),2)))
plot(a$Epi,a$NK,xlab="Epithelial cell fraction",ylab="NK cell fraction")
legend(x='bottomleft', legend=paste('Cor =',round(cor(a$Epi, a$NK),2)))
par(mfrow=c(1,1))
corrplot(cor(a,method="pearson"),type="upper")

ggplot(a, aes(x=Epi, y=Neutro)) + 
  geom_point()+
  geom_smooth(method=lm)+
  theme_classic()  

targets=cbind(targets,frac.m)

t.test(targets[1:36,]$Epi,targets[37:48,]$Epi)
t.test(targets[1:36,]$Fib,targets[37:48,]$Fib)
t.test(targets[1:36,]$B,targets[37:48,]$B)
t.test(targets[1:36,]$NK,targets[37:48,]$NK)
t.test(targets[1:36,]$CD4T,targets[37:48,]$CD4T)
t.test(targets[1:36,]$CD8T,targets[37:48,]$CD8T)
t.test(targets[1:36,]$Mono,targets[37:48,]$Mono)
t.test(targets[1:36,]$Neutro,targets[37:48,]$Neutro)
t.test(targets[1:36,]$Eosino,targets[37:48,]$Eosino)

CellDMC(beta,pheno.v,frac.m,adjPMethod = "fdr",adjPThresh = 0.05,cov.mod = NULL,sort = FALSE,mc.cores = 1)

#quality control after normalization
par(mfrow=c(1,2))
densityByProbeType(mSet[,1], main = "Raw",cex.legend = 0.75)
densityByProbeType(mSetSw[,1], main = "SWAN",cex.legend = 0.75)

qc <- getQC(mSet)
par(mfrow=c(1,1))
plotQC(qc)
densityPlot(mSet, sampGroups = phenoData$Sample_Group)
densityBeanPlot(mSet, sampGroups = phenoData$Sample_Group)
controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")

par(mfrow=c(1,1),mar=c(5,5,5,5))
densityPlot(mSetSw, sampGroups = phenoData$Sample_Group)
densityBeanPlot(mSetSw, sampGroups = phenoData$Sample_Group)

par(mfrow=c(1,2),mar=c(5,5,5,5))
densityPlot(mSet, sampGroups = phenoData$Sample_Group,legend=FALSE)
densityPlot(mSetSw, sampGroups = phenoData$Sample_Group,legend=FALSE)


# REMP (we find all repetitive elements in all filtered probes)

remp.resAlu <- remprofile(ratioConvert(mSetFinal), REtype = "Alu", annotation.source = "UCSC", genome = "hg19")
details(remp.resAlu)
remp.resL1 <- remprofile(ratioConvert(mSetFinal), REtype = "L1", annotation.source = "UCSC", genome = "hg19")
details(remp.resL1)

alu = rempB(remp.resAlu)
L1=rempB(remp.resL1)

allRepeats =  rbind(  alu , L1 ) 

CpGRE = allRepeats@rownames
'%!in%' <- function(x,y)!('%in%'(x,y))
mSetFinal=mSetFinal[rownames(mSetFinal)%!in%CpGREIG,]

mSetFinal=mSetFinal[rownames(mSetFinal)%in%unique(TEST$...1),]

#remove IG probes

imprintedGenes = read.table("imprinted.genes.tsv" )
AnnotationCpGs <- read.delim("AnnotationCpGs.csv")
ann850k=getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

idx = which( AnnotationCpGs$name %in% imprintedGenes$V1 ) 
CpGIG = rownames(AnnotationCpGs[idx,])

CpGREIG = unique(c(CpGRE,CpGIG))

# calculate M and beta values
meth <- getMeth( mSetFinal ) 
unmeth <- getUnmeth( mSetFinal ) 
Mval = getM(mSetFinal)
Mval2=log2((meth)/(unmeth))
beta <- getBeta(mSetFinal)
targets$Sample_ID = colnames(beta)

#remove RE and IG CpGs

beta=beta[-which(rownames(beta)%in%CpGREIG),]
Mval=Mval[-which(rownames(Mval)%in%CpGREIG),]

#SVD = single value decomposition (more comprehensive influence of covariates)

champ.SVD(beta=beta,pd=targets)

# Filtering DMP with delta beta

moyenneCTRL = rowMeans( beta[,c(targets$Sample_ID[which(targets$Sample_Group=="CTRL")])])
moyenneAMP  = rowMeans( beta[, c(targets$Sample_ID[which(targets$Sample_Group=="ART")]) ])

deltaBetaAMP =  moyenneAMP - moyenneCTRL

# Finding DMP using limma

PMA = factor(targets$Sample_Group)
PMA <- relevel(PMA, "CTRL")
Gender = factor(targets$Gender)
Slide = factor(targets$Slide)
Epi= targets$Epi
B=targets$B
Mono=targets$Mono
Neutro=targets$Neutro
NK=targets$NK
design <- model.matrix(~PMA+Gender+Slide+Epi+B+Mono+Neutro+NK)


fit <- lmFit(Mval, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))
DMPs <- topTable(fit2, num=Inf, coef=2)

dmpDB = merge( DMPs , deltaBetaAMP , by=c( "row.names" ) , all.x=T )
#idx = which((dmpDB$adj.P.Val < 0.05) & (abs(dmpDB$y)) >= 0.05 )
idx = which((dmpDB$adj.P.Val < 0.05) )
dmpAMPFiltered = dmpDB[idx,]
dim(dmpAMPFiltered)[1]
rownames(dmpAMPFiltered) = dmpAMPFiltered$Row.names


#Plotting the best DMPs

dmpAMPFilteredordered = dmpAMPFiltered[order(dmpAMPFiltered$adj.P.Val),]
dmpAMPFilteredordered = dmpAMPFilteredordered[which(abs(dmpAMPFilteredordered$y)>0.05),]

cpgs <- rownames(dmpAMPFilteredordered)
par(mfrow=c(3,5))
plotCpg(beta, cpg=cpgs, pheno=targets$Sample_Group, ylab = "Beta values")

cpgs <- rownames(dmpAMPFilteredordered)[1:15]
par(mfrow=c(3,5))
for (k in 1:15){
  a=rbind(beta[cpgs[k],],targets$Sample_Group)
  a=data.frame(t(a))
  colnames(a)=c("Beta value","Group")
  class(a$`Beta value`)="numeric"
  a$Group=factor(a$Group,levels=c("CTRL","ART"))
  boxplot(a$`Beta value` ~ a$Group , 
          col=c("#D0D0D0","#E6BBAD") , 
          ylab="Beta value" , xlab="Group")
  title(cpgs[k])
}

#Annotation via ChAMP

genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

data(probe.features.epic)
annEPIC=probe.features
rownames( dmpAMPFiltered ) = dmpAMPFiltered$Row.names
t = merge( dmpAMPFiltered , annEPIC , by=c( "row.names" ) , all.x=T )
t=cbind(t[,c(2:11,15:16)])
rownames(t)=t$Row.names
t = merge( t , AnnotationCpGs , by=c( "row.names" ) , all.x=T )
table(t$feature)
TableDMP = t
write.table(t,"DMP AMP.xls",sep="\t",col.names=TRUE)

#graph genomic features

t=merge(dmpDB,annEPIC,by.x="Row.names",by.y=0)
a=rbind(table(TableDMP$feature)/sum(table(TableDMP$feature)),table(t$feature)/sum(table(t$feature)))
rownames(a)=c("DMP","All CpGs")
par(mfrow=c(1,1))
a=a[,c(7,8,3,1,5,4,2,6)]
barplot(a, 
        border="white",
        beside=T, 
        legend=rownames(a), 
        xlab="Genomic features",ylab = "% of CpGs\n",xaxt="n",cex.axis=1.7,cex.lab=2)
axis(1,cex.axis=1.4,labels=colnames(a),at=seq(2,24,3))

a=rbind(table(TableDMP$feature),table(t$feature))
for (k in 1:8){
  totDMP=sum(a[1,])
  totnonDMP = sum(a[2,])
  test=matrix(c(a[1,k],a[2,k],totDMP-a[1,k],totnonDMP-sum(a[2,k])),nrow=2)
  print(colnames(a)[k])
  print(prop.test(test),correct=FALSE)
} 

a=rbind(table(TableDMP$cgi)/sum(table(TableDMP$cgi)),table(t$cgi)/sum(table(t$cgi)))
rownames(a)=c("DMP","All CpGs")
a=a[,c("island","shore","shelf","opensea")]
par(mfrow=c(1,1),mar=c(8, 8, 4.1, 2.1))
barplot(a, 
        border="white",
        beside=T, 
        legend=rownames(a), 
        xlab="\n\nCpG features",ylab = "% of CpGs\n",xaxt="n",cex.axis=1.7,cex.lab=2)
axis(1,cex.axis=2,labels=colnames(a),at=seq(2,12,3))

a=rbind(table(TableDMP$cgi),table(t$cgi))
for (k in 1:4){
  totDMP=sum(a[1,])
  totnonDMP = sum(a[2,])
  test=matrix(c(a[1,k],a[2,k],totDMP-a[1,k],totnonDMP-sum(a[2,k])),nrow=2)
  print(colnames(a)[k])
  print(prop.test(test),correct=FALSE)
} 

# DMR with covariates
PMA <- relevel(PMA, "CTRL")
designMatrix <- model.matrix(~PMA+Gender+Slide+Epi+B+Mono+Neutro+NK)
myAnnotOnAMP = cpg.annotate("array", ratioConvert(mSetFinal) , analysis.type="differential", design= designMatrix , coef = 2, arraytype = "EPIC")
dmrcoutput <- dmrcate( myAnnotOnAMP , lambda=1000, C=2,pcutoff="fdr",min.cpgs=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
dmrAMPAnn = matchGenes( results.ranges , genes )
dim(dmrAMPAnn)[1]
sum(abs(results.ranges$meandiff>0.05)==TRUE)
table(dmrAMPAnn$region)
table(results.ranges@seqnames)
table(results.ranges@elementMetadata@listData$no.cpgs)
View(data.frame(table(TableDMR$name)))

t = cbind( data.frame(results.ranges) , dmrAMPAnn )
TableDMR=t

dmpAMPDB = merge( DMPs , deltaBetaAMP , by=c( "row.names" ) , all.x=T )
ann850k2 = data.frame(ann850k)

#Fonction multi-mixed order

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

#DMR Novakovic

for (k in 1:dim(TableDMR)[1])
  TableDMR$ID[k]<-paste(k)
class(TableDMR$ID)="numeric"

DMRranges = GRanges(seqnames=TableDMR$seqnames,ranges=IRanges(start=TableDMR$start,end=TableDMR$end))
mcols(DMRranges)=DataFrame(ID=TableDMR$ID,gene=TableDMR$name)

rownames(dmpAMPDB)=dmpAMPDB$Row.names
t=merge(dmpAMPDB,ann850k2,by=0)
t=t[which(abs(t$y)>0.05),]
DMPranges = GRanges(seqnames=t$chr,ranges=IRanges(start=t$pos,end=t$pos))
mcols(DMPranges) = DataFrame(CpG=t$Row.names,Gene=t$gene)

DMRdb = findOverlaps(DMPranges, DMRranges,select="first",ignore.strand=TRUE)
class(DMRdb)="vector"
sort(unique(DMRdb))
a=(table(DMRdb)>=1)
a=which(a=="TRUE")

DMR=data.frame(DMRranges)[a,]
TableDMR=merge(TableDMR,DMR,by="ID",no.dups = TRUE)

TableDMR=TableDMR[,1:24]
colnames(TableDMR)[c(2:6)]=c("seqnames","start","end","width","strand")

TableDMRpromoter = TableDMR[which(TableDMR$region=="promoter"),]

dmpDB = merge( DMPs , deltaBetaAMP , by=c( "row.names" ) , all.x=T )
dmpAMPFDR = dmpDB
dim(dmpAMPFDR)[1]
rownames(dmpAMPFDR) = dmpAMPFDR$Row.names
t = merge( dmpAMPFDR , ann850k2 , by=c( "row.names" ) , all.x=T )
TableDMPFDR=t
TableDMPFDR=TableDMPFDR[multi.mixedorder(TableDMPFDR$chr,TableDMPFDR$pos,decreasing = FALSE),]
TableDMR=TableDMR[order(TableDMR$seqnames,TableDMR$start),]

for (k in 1:dim(TableDMR)[1])
  TableDMR$ID[k]<-paste(k)
class(TableDMR$ID)="numeric"

DMPranges = GRanges(seqnames=TableDMPFDR$chr,ranges=IRanges(start=TableDMPFDR$pos,end=TableDMPFDR$pos))
mcols(DMPranges) = DataFrame(CpG=TableDMPFDR$Row.names)

DMRranges = GRanges(seqnames=TableDMR$seqnames,ranges=IRanges(start=TableDMR$start,end=TableDMR$end))
mcols(DMRranges)=DataFrame(ID=TableDMR$ID)

#Gene Ontology

sigCpGs = rownames(dmpAMPFiltered)
all <- as.vector(dmpAMPDB$Row.names)
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE,collection="KEGG")
topGSA(gst, number=10)
table(gst$FDR<0.05)

HypoDMP=rownames(TableDMP[which(TableDMP$y<0),])
HyperDMP=rownames(TableDMP[which(TableDMP$y>0),])
HypoDMP=rownames(dmpAMPFiltered[which(dmpAMPFiltered$y<0),])
HyperDMP=rownames(dmpAMPFiltered[which(dmpAMPFiltered$y>0),])
gst <- gometh(sig.cpg=HypoDMP, all.cpg=all, plot.bias=TRUE, collection="KEGG")
topGSA(gst, number=10)
table(gst$FDR<0.05)
gst <- gometh(sig.cpg=HyperDMP, all.cpg=all, plot.bias=TRUE, collection="KEGG")
topGSA(gst, number=10)
table(gst$FDR<0.05)

enrichment_GO <- goregion(results.ranges,all.cpg=all, collection = "GO", array.type = "EPIC",plot.bias = TRUE, prior.prob = TRUE,)
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),] 
head(as.matrix(enrichment_GO), 10)
table(enrichment_GO$FDR<0.05)

enrichment_GO <- goregion(DMRranges,all.cpg=all, collection = "GO", array.type = "EPIC",plot.bias = TRUE, prior.prob = TRUE,)
enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),] 
head(as.matrix(enrichment_GO), 10)
table(enrichment_GO$FDR<0.05)


# Volcano plot

dmpAMPDB = merge( DMPs , deltaBetaAMP , by=c( "row.names" ) , all.x=T )

colvect = rep( "grey" , nrow( dmpAMPDB ) )

idx = which( (dmpAMPDB$y >= 0) & (dmpAMPDB$adj.P.Val)<=0.05 ) 
colvect[idx] = "#FF0300"
hyper=length(idx)

idx = which( (dmpAMPDB$y <= 0) & (dmpAMPDB$adj.P.Val)<=0.05 ) 
colvect[idx] = "#00FF02"
hyper=length(idx)

idx = which( (dmpAMPDB$y >= 0.05) & (dmpAMPDB$adj.P.Val)<=0.05 ) 
colvect[idx] = "#B60200" 
hyper=length(idx)

idx = which( (dmpAMPDB$y <= -0.05) & (dmpAMPDB$adj.P.Val)<=0.05 ) 
colvect[idx] = "#00A801"
hypo=length(idx)

hypo/(hypo+hyper)*100
hyper/(hypo+hyper)*100

par(mfrow=c(1,1),mar=c(5.1, 8, 4.1, 5.1))
plot(dmpAMPDB$y,-log10(dmpAMPDB$adj.P.Val), pch = 20 , main = "Volcano plot" , col = colvect , xlab = "Δβ" , ylab = expression('−log'[10]*'(p-value)'),cex.lab=1.5,cex.axis=1.5)      
abline( h = -log10(0.05) , col = "yellow" , lty = 2 ) 
abline( v = 0.05 , col = "yellow" , lty = 2 ) 
abline( v = -0.05 , col = "yellow" , lty = 2 ) 

plot(dmpAMPDB$y,-log10(dmpAMPDB$P.Value), pch = 20 , main = "Volcano plot" , col = colvect , xlab = "Δβ" , ylab = expression('−log'[10]*'('*italic(p)*'-value)'),cex.lab=1.5,cex.axis=1.5)      
abline( h = -log10(8.453175e-06) , col = "yellow" , lty = 2 ) 
abline( v = 0.05 , col = "yellow" , lty = 2 ) 
abline( v = -0.05 , col = "yellow" , lty = 2 ) 

plotCpg(beta, cpg = "cg27266479", pheno = targets$Sample_Group, type = "categorical", measure = "beta", 
        ylab = "methylation (beta)", ylim = c(0, 1))

#Manhattan plot

t=merge(dmpAMPDB,annEPIC, by.x="Row.names",by.y=0,all.x=T)
Manhattan = data.frame(t$MAPINFO,t$CHR,t$adj.P.Val)
rownames(Manhattan)=rownames(t)
colnames(Manhattan)=c("BP","CHR","P")
Manhattan=data.frame(Manhattan)
Manhattan$CHR=str_replace(Manhattan$CHR,"X","23")
Manhattan$CHR=str_replace(Manhattan$CHR,"Y","24")
class(Manhattan$CHR)="numeric"
manhattan(Manhattan, chr="CHR", bp="BP", p="P",suggestiveline = -log10(0.05),chrlabs=c(1:22,"X","Y") )

qq(Manhattan$P)

t1=merge(dmpAMPDB[which(dmpAMPDB$y>0),],annEPIC, by.x="Row.names",by.y=0,all.x=T)
t2=merge(dmpAMPDB[which(dmpAMPDB$y<0),],annEPIC, by.x="Row.names",by.y=0,all.x=T)
Manhattan1 = data.frame(t1$MAPINFO,t1$CHR,t1$P.Value)
Manhattan2 = data.frame(t2$MAPINFO,t2$CHR,t2$P.Value)
rownames(Manhattan1)=rownames(t1)
colnames(Manhattan1)=c("BP","CHR","P")
rownames(Manhattan2)=rownames(t2)
colnames(Manhattan2)=c("BP","CHR","P")
Manhattan1=data.frame(Manhattan1)
Manhattan1$CHR=str_replace(Manhattan1$CHR,"X","23")
Manhattan1$CHR=str_replace(Manhattan1$CHR,"Y","24")
class(Manhattan1$CHR)="numeric"
Manhattan2=data.frame(Manhattan2)
Manhattan2$CHR=str_replace(Manhattan2$CHR,"X","23")
Manhattan2$CHR=str_replace(Manhattan2$CHR,"Y","24")
class(Manhattan2$CHR)="numeric"
manhattan(Manhattan1, chr="CHR", bp="BP", p="P",suggestiveline = -log10(0.05),chrlabs=c(1:22,"X","Y"), ylim=c(0,4) )
manhattan(Manhattan2, chr="CHR", bp="BP", p="P",suggestiveline = -log10(0.05),chrlabs=c(1:22,"X","Y"), ylim=c(0,4) )

# Plotting DMRs

gen <- "hg19"
dmrIndex <- 1772
chrom <- as.character(seqnames(DMRranges[dmrIndex]))
start <- as.numeric(start(DMRranges[dmrIndex]))
end <- as.numeric(end(DMRranges[dmrIndex]))
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)
bVals= beta
ann850k=getAnnotation( IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann850kSub <- ann850k[match(rownames(bVals),ann850k$Name),]
ann850kOrd <- ann850kSub[order(ann850kSub$chr,ann850kSub$pos),]
head(ann450kOrd)
bValsOrd <- bVals[match(ann850kOrd$Name,rownames(bVals)),]
head(bValsOrd)

cpgData <- GRanges(seqnames=Rle(ann850kOrd$chr),
                   ranges=IRanges(start=ann850kOrd$pos, end=ann850kOrd$pos),
                   strand=Rle(rep("*",nrow(ann850kOrd))),
                   betas=bValsOrd)
cpgData <- subsetByOverlaps(cpgData, DMRranges[dmrIndex])


pal <- brewer.pal(8,"Dark2")
methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")
tracks <- list(iTrack, gTrack, methTrack, dmrTrack,
               rTrack)
sizes <- c(2,2,5,2,3) # set up the relative sizes of the tracks
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
ann850k2 = data.frame(ann850k)
BED = merge(dmpAMPFiltered,ann850k2,by=0)
BED = cbind(BED$chr,BED$pos,BED$pos,BED$Name)
#Table can be used with GREAT tool to assess proximal gene
write.table(BED,"DMP.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)


#Plot cell proportions

a=targets[,c(1,8,15)]
for (k in 1:8){
  colnames(a)[3]=colnames(targets)[15+k]
  a=rbind(a,targets[,c(1,8,15+k)])
}
a$Celltype=c(rep("Epi",48),rep("Fib",48),rep("B",48),rep("NK",48),rep("CD4T",48),rep("CD8T",48),rep("Mono",48),rep("Neutro",48),rep("Eosino",48))
colnames(a)=c("Individual","Group","Cell fraction","Cell type")
a$`Cell type`=factor(a$`Cell type`,levels=c("Epi","Fib","B","NK","CD4T","CD8T","Mono","Neutro","Eosino"))
a$Group=factor(a$Group,levels=c("CTRL","ART"))
ggplot(a, aes(x=`Cell type`, y=`Cell fraction`, fill=Group)) + 
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()+theme_classic()+ scale_fill_manual(values=c("#A0A0A0","#F6A209"))+ylim(0, 1)+
  theme(axis.text.x = element_text(size=14),axis.text.y = element_text(size=14),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20)) + 
  xlab("\nCell type") + ylab("Cell fraction\n")
