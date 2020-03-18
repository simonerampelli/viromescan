#R script for output SARS CoV-2 COVID-19
pDatafile1='coronavirus.genes.txt'
genes=read.delim(pDatafile1, header = F)

x=matrix(ncol=2,nrow=2)
colnames(x)=c("counts","%")
rownames(x)=c("SARS CoV-2 COVID-19","other-Coronavirus")

x[1,1]=sum(as.numeric(genes[1:45,3]))
x[2,1]=sum(as.numeric(genes[46:100,3]))
x[1,2]=(sum(as.numeric(genes[1:45,3]))/sum(as.numeric(genes[,3])))*100
x[2,2]=(sum(as.numeric(genes[46:100,3]))/sum(as.numeric(genes[,3])))*100

#first output is a table
write.table(x,file="covid19.txt")