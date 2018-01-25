library("qtl")
library("ASMap")

pdf("Radish.pdf")
mapthis <- read.cross("csvs", "", "Radish_final.csv", "phenotype.txt", estimate.map= F)

summary(mapthis)

mapthis <- convert2bcsft(mapthis, F.gen=0, BC.gen=0, estimate.map= F)

summary(mapthis)

n.ind <- summary(mapthis)$n.ind
n.mar <- sum(summary(mapthis)$n.mar)

##### Missing status #####

plotMissing(mapthis)

par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

##### Individual filtering

geno_typed <- 0.5
ind_missed <- 0.75 
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>n.mar*geno_typed))

##### Marker filtering

nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < n.ind*ind_missed])
dim(todrop)
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

##### Remove duplicates #####

cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

##### Remove mismatch genotypes

g <- pull.geno(mapthis)

for(i in 1:nrow(wh)) {
if (nrow(wh) == 0) break

tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]

for (j in 1:length(tozero)){
  Chr <- as.numeric(substr(unlist(strsplit(names(tozero[j]), "_"))[1],2,3))
  mapthis$geno[[Chr]]$data[wh[i,1],tozero[j]] <- NA
}
}

##### Remove duplicate individuals

#mapthis <- subset(mapthis, ind=-wh[,2])

print(dup <- findDupMarkers(mapthis, exact.only=FALSE))

##### Segregation patterns #####

gt <- geno.table(mapthis)
(todrop <- rownames(gt[gt$P.value < 0.05/totmar(mapthis),]))
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

##### Visualization

g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)

for(i in 1:3){
		plot(gfreq[i,], ylab="Genotype frequency", main=c("AA","AB", "BB")[i],ylim=c(0,1))
}

##### Check switched alleles #####

#mapthis <- est.rf(mapthis)
#checkAlleles(mapthis, threshold=5)

test.brapa = mstmap.cross(mapthis, id="id", suffix="alpha", bychr=T, p.value=1e-6, objective.fun="ML", noMap.dist=0.01)

test.brapa = rescalemap(test.brapa, scale=1e-1)

par(mfrow=c(1,1))

plot.map(test.brapa)

dev.off()

chrNames = c("R1.A", "R2.A", "R3.A", "R4.A", "R5.A", "R6.A", "R7.A", "R8.A", "R9")
test.brapa = subset(test.brapa, chr=chrNames)

for (phenoNum in 1:99) {
    #phenoNum = 36
    ChrNum = length(test.brapa$geno)
    phenoName = colnames(test.brapa$pheno)[phenoNum]
    pdf(paste(phenoName, ".pdf", sep="") ,width=10,height=10)
    num_permute = 100
    
    col_num = ceiling(sqrt(ChrNum))
    par(mfrow=c(col_num, col_num), oma=c(5,5,5,0), mai=c(0.3, 0.25, 0.25, 0.25))
    
    out.map <- cim(test.brapa, pheno.col=phenoNum, n.marcovar=3, window=5, # Bolting_1st
        method="imp",
        imp.method="imp", error.prob=0.0001,
        map.function="kosambi")
    
    perms <- cim(test.brapa, pheno.col=phenoNum, n.marcovar=3, window=5, # Bolting_2nd
        method="imp",
        imp.method="imp", error.prob=0.0001,
        map.function="kosambi",
        n.perm=num_permute)
    
        for (each_chr in 1:ChrNum) {
            plot(out.map, chr=chrNames[each_chr], main = paste("Chr", toString(each_chr)), xlab="", ylab="", ylim=c(0,20), cex.axis=1.5, cex.main=1.5)
            add.threshold(out=out.map, perms=perms, alpha = 0.05, col = "red", lty=2)
        }
    
    #plot(out.RIL, chr=2, main = paste("Chr", toString(2)), xlab="", ylab="", ylim=c(0,20))
    #title("My Title", outer=TRUE)
    #mtext(paste(phenoName, "QTL scan with significance threshold (alpha=0.05)") , side = 3, line = 0, outer = TRUE)
    mtext("Genetic distance (cM)", side = 1, line = 2, outer = TRUE, cex=1.2)
    mtext("lod", side = 2, line = 2, outer = TRUE, cex=1.2)
    mtext(phenoName, side = 3, line = 2, outer = TRUE, cex=1.5)
    
    dev.off()
    
}

