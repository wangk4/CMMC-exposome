library(GOSemSim)
library(org.Hs.eg.db)
library(GO.db)

setwd("/path/to/geneID/folder/")

dat<-read.table("gene_entrezid.txt",header = FALSE)
genename = dat$V1

hsGO <- godata('org.Hs.eg.db', ont="CC")
out<-mgeneSim(genes=genename,
         semData=hsGO, measure="Wang",verbose=FALSE)
write.table(out,"gene.sim.GOSemSin.CC.txt",row.names = TRUE,col.names = TRUE)

hsGO <- godata('org.Hs.eg.db', ont="BP")
out<-mgeneSim(genes=genename,
         semData=hsGO, measure="Wang",verbose=FALSE)
write.table(out,"gene.sim.GOSemSin.BP.txt",row.names = TRUE,col.names = TRUE)

hsGO <- godata('org.Hs.eg.db', ont="MF")
out<-mgeneSim(genes=genename,
         semData=hsGO, measure="Wang",verbose=FALSE)
write.table(out,"gene.sim.GOSemSin.MF.txt",row.names = TRUE,col.names = TRUE)

