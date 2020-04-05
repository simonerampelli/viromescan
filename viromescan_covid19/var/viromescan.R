pDatafile1='../final.genes.txt'
genes=read.delim(pDatafile1, header = F)
if ((length(rownames(genes))) == 93) {
  pDatafile2='../HumanDNAcomplete.txt'} 
if ((length(rownames(genes))) == 665) {
  pDatafile2='../HumanALLcomplete.txt'} 
if ((length(rownames(genes))) == 4371){
  pDatafile2='../VirusALLcomplete.txt'}
if ((length(rownames(genes))) == 1647){
  pDatafile2='../VirusDNAcomplete.txt'}
if ((length(rownames(genes))) == 666) {
  pDatafile2='../HumanALL+covid19complete.txt'}

IDs=read.delim(pDatafile2, header = F, sep = ';')
genes=genes[1:(length(rownames(genes))-1),]
Ab=genes$V3
Ab=as.data.frame(Ab)
Ab=cbind(IDs,Ab)
colnames(Ab)=c('Family', 'Genera','Species','Counts')

Ab1=cbind(as.character(IDs$V1),as.numeric(genes$V3))
colnames(Ab1)=c('Family','counts')
lista_family=levels(IDs$V1)
fam=matrix(data=0,nrow=length(lista_family),ncol=1)
rownames(fam)=as.character(lista_family)
colnames(fam)='counts'
for (i in 1: length(rownames(fam)))
{ a=as.numeric(grep(lista_family[i], Ab1))
b=Ab1[a,]
r=matrix(data= 0, nrow=2, ncol=2)
b=rbind(b,r)
c=sum(as.numeric(b[,2]))
fam[i,]=c
}

genera=c()
for (i in 1:length(IDs$V2))
{genera=(paste(as.character(IDs$V1),as.character(IDs$V2),sep=";"))}
lista_genera=levels(as.factor(genera))
Ab1=cbind(as.character(genera),as.numeric(genes$V3))
colnames(Ab1)=c('Genera','counts')
gen=matrix(data=0,nrow=length(lista_genera),ncol=1)
rownames(gen)=as.character(lista_genera)
colnames(gen)='counts'
for (i in 1: length(rownames(gen)))
{ a=as.numeric(grep(lista_genera[i], Ab1))
b=Ab1[a,]
r=matrix(data= 0, nrow=2, ncol=2)
b=rbind(b,r)
c=sum(as.numeric(b[,2]))
gen[i,]=c
}

specie=c()
for (i in 1:length(IDs$V3))
{specie=(paste(as.character(IDs$V1),as.character(IDs$V2),as.character(IDs$V3),sep=";"))}
spe=matrix(ncol = 1, nrow = length(specie))
rownames(spe)=specie
colnames(spe)=c('counts')
spe[,1]=Ab$Counts


fam_perc=c()
for (i in 1 : (length(rownames(fam))))
{
 fam_perc[i]=round(fam[i]/sum(fam[,1]),4)
}
fam_perc2=as.matrix(fam_perc)
rownames(fam_perc2)=rownames(fam)
colnames(fam_perc2)=c('')
ff=(fam_perc2>0)
fam_perc2=fam_perc2[ff,]
fam_perc2=as.matrix(fam_perc2)
fff=fam>0
fam=fam[fff,]
fam=as.matrix(fam)
colnames(fam_perc2)=c('')
colnames(fam)=c('')
write.table(fam_perc2, file='Family_level_results-%.txt', sep='\t')
write.table(fam, file='Family_level_results-Counts.txt', sep='\t')
fam_perc2=as.matrix(fam_perc2)
filt=rowSums(fam_perc2)>0.005
fam_perc2=fam_perc2[filt,]
fam_perc2=t(fam_perc2)
fam_perc2=t(fam_perc2)
Other=round(1-(sum(fam_perc2[,1])),2)
fam_perc2=rbind(fam_perc2, Other)
legend=c(rownames(fam_perc2))
pdf(file = 'viromescan_results_Family_level.pdf')
op=par(mfrow = c(1,2))
barplot(fam_perc2, main = 'Virus ab.- Fam. lev.',col= c(rainbow(length(rownames((fam_perc2)))-1),'gray'), width = 0.8 , border = T)
plot.new()
legend('topleft', legend = legend, fill= c(rainbow(length(rownames((fam_perc2)))-1),'gray'), cex = 0.65,y.intersp=0.70 ,x.inters=0.4, lty = F, border = 'white', bty = 'n')
dev.off()


gen_perc=c()
for (i in 1 : (length(rownames(gen))))
{
  gen_perc[i]=round(gen[i]/sum(gen[,1]),4)
}
gen_perc2=as.matrix(gen_perc)
rownames(gen_perc2)=rownames(gen)
colnames(gen_perc2)=c('')
gg=(gen_perc2>0)
gen_perc2=gen_perc2[gg,]
gen_perc2=as.matrix(gen_perc2)
ggg=(gen>0)
gen=gen[ggg,]
gen=as.matrix(gen)
write.table(gen_perc2, file='Genera_level_results-%.txt', sep='\t')
write.table(gen, file='Genera_level_results-Counts.txt', sep='\t')
gen_perc2=as.matrix(gen_perc2)
filt=rowSums(gen_perc2)>0.005
gen_perc2=gen_perc2[filt,]
gen_perc2=t(gen_perc2)
gen_perc2=t(gen_perc2)
Other=round(1-(sum(gen_perc2[,1])),2)
gen_perc2=rbind(gen_perc2, Other)
legend=c(rownames(gen_perc2))
pdf(file = 'viromescan_results_Genera_level.pdf')
op=par(mfrow = c(1,2))
barplot(gen_perc2, main = 'Virus ab.- Gen. lev.',col= c(rainbow(length(rownames((gen_perc2)))-1),'gray'), width = 0.8 , border = T)
plot.new()
legend('topleft', legend = legend, fill= c(rainbow(length(rownames((gen_perc2)))-1),'gray'), cex = 0.32,y.intersp=0.90 ,x.inters=0.4, lty = F, border = 'white', bty = 'n')
dev.off()

spe_perc=c()
for (i in 1 : (length(rownames(spe))))
{
  spe_perc[i]=round(spe[i]/sum(spe[,1]),4)
}
spe_perc2=as.matrix(spe_perc)
rownames(spe_perc2)=rownames(spe)
colnames(spe_perc2)=c('')
ss=(spe_perc2>0)
spe_perc2=spe_perc2[ss,]
spe_perc2=as.matrix(spe_perc2)
sss=(spe>0)
spe=spe[sss,]
spe=as.matrix(spe)
colnames(spe_perc2)=c('')
colnames(spe)=c('')
write.table(spe_perc2, file='Species_level_results-%.txt', sep='\t')
write.table(spe, file='Species_level_results-Counts.txt', sep='\t')
filt=rowSums(spe_perc2)>0.005
spe_perc2=spe_perc2[filt,]
spe_perc2=t(spe_perc2)
spe_perc2=t(spe_perc2)
Other=round(1-(sum(spe_perc2[,1])),2)
spe_perc2=rbind(spe_perc2, Other)
legend=c(rownames(spe_perc2))
pdf(file = 'viromescan_results_Species_level.pdf')
op=par(mfrow = c(1,2))
barplot(spe_perc2, main = 'Virus ab.- Species lev.',col= c(rainbow(length(rownames((spe_perc2)))-1),'gray'), width = 0.8 , border = T)
plot.new()
legend('topleft', legend = legend, fill= c(rainbow(length(rownames((spe_perc2)))-1),'gray'), cex = 0.30,y.intersp=0.90 ,x.inters=0.4, lty = F, border = 'white', bty = 'n')
dev.off()
