
# Limpiamos el ambiente dentro de R. 
rm(list = ls(all = TRUE))

library("coloc")
library("Rfast")

file1 <- "ENSG00000159363.17.GTEx.eqtl.txt"
file2 <- "Height_HIS.GWAS.ENSG00000159363.coloc.input"
trait1 <- "ENSG00000159363"
trait2 <- "height"

getwd()

# Van a recibir algo así
# [1] "/home/carleonidastest/"

setwd("~/Test")

dir.create("coloc")

setwd("coloc")

# Los datos para este workshop se encuentran en otra parte del servidor, para acceder a esto 
data_path <- "/storage/workshops/colocalization-workshop/"

# Para cargar los datos usamos la función "paste0()", que sirve para juntar dos strings. 
# En este caso, el nombre del archivo y dónde se localiza. Por ejemplo: 
paste0(data_path, file1)

# NOTA: Hay otra función similar a `paste0()`, esta es `paste()`. Son identicas 
# en su función, (las dos juntan strings), la única diferencia es que `paste()` 
# pone un espacio entre los strings, mientras que `paste0()` los ignora.

paste("a", "b", "c")
paste0("a", "b", "c")


eqtl <- read.table(
  paste0(data_path, file1), # eQTLs
  header = T,               # Indica que la tabla tiene nombres de las variables
  stringsAsFactors = F)      

gwas <- read.table(paste0(data_path, file2),                    # Estadisticas descriptivas
                   header = T,               # La tabla tiene nombres de las variables
                   stringsAsFactors = F)

rownames(eqtl) <- eqtl$rsID
rownames(gwas) <- gwas$MarkerName

if (dim(gwas[which(gwas$p_value < 1e-5),])[1] > 0){
  cat(" ¿Existe una señal sugerente en GWAS2 (P < 1e-5)?", "\n",
      "Sí, sí hay una señal sugerente.", "\n\n")
} else {
  cat(" ¿Existe una señal sugerente en GWAS2 (P < 1e-5)?", "\n",
      "No, no se identificó ninguna señal.", "\n\n")
}


eqtl <- eqtl[gwas$MarkerName, ]

eqtl <- eqtl[which(!is.na(eqtl$pval_asc)), ]

dir.create("LD")

write.table(
  # Seleccionamos las variables 'rsID' y 'REF' de la base de datos eQTL
  x = eqtl[, c('rsID', 'REF')], 
  
  # Indicamos la ruta y el nombre del archivo en el cual vamos a guardar la información  
  file = paste0('LD/', trait1, '.SNP.REF.txt'),
  
  # Otras opciones: 
  col.names = F,        # se asignan nombres a las filas o columnas
  row.names = F,        
  sep = '\t',           # Cómo se separan los valores, en este caso con espacio. 
  quote = F)   


write.table(
  # Seleccionamos la variable 'rsID' de la base de datos DE GWAS
  x = eqtl[, c('rsID')], 
  
  # Indicamos la ruta y el nombre del archivo en el cual vamos a guardar la información  
  file =  paste0('LD/', trait1, ".SNP.extract"), 
  
  # Otras opciones
  col.names=F, row.names=F, sep='\t', quote=F) 

plink <- paste0("/apps/plink/", "plink") 

CHR <- "chr1"

system(
  paste0(plink, 
         # Indica que lea el archivo ".vcf" (Variant Call Format) 
         ' --vcf ', data_path, '1KG.AMR.GRCH38.rsID.', CHR, '.region.vcf',
         
         # lista de ID de variantes
         ' --extract LD/', trait1, '.SNP.extract',
         
         # Indicador de sistema, inidca cuantos "threads" debe usar el CPU
         ' --threads 1', 
         
         # Actualiza los alelos de referencia. 
         ' --update-ref-allele LD/', trait1, '.SNP.REF.txt 2 1',
         
         # calcula y reporta las correlaciones brutas de recuento de alelos entre variantes
         ' --r square',
         
         # Directorio y nombre a dar a los archivos que se generen del análisis
         ' --out LD/AMR', trait1, '.LD', 
         
         # Indica que genere solamente un archivo ".bim"
         ' --make-just-bim'))

system(
  paste0(plink, 
         ' --vcf ', data_path, '1KG.EUR.GRCH38.rsID.', CHR, '.region.vcf', 
         ' --extract LD/', trait1,'.SNP.extract', 
         ' --threads 1',
         ' --update-ref-allele LD/', trait1,'.SNP.REF.txt 2 1', 
         ' --r square',
         ' --out LD/EUR', trait1,'.LD', 
         ' --make-just-bim'))


amrbim <- read.table(paste0('LD/AMR', trait1, '.LD.bim'), header=F)

eurbim <- read.table(paste0('LD/EUR', trait1, '.LD.bim'), header=F)

colnames(amrbim) <- c("chr", "rsid", "cM", "position", "allele1", "allele2")
colnames(eurbim) <- c("chr", "rsid", "cM", "position", "allele1", "allele2")

nrow(amrbim)

ldamr <-  as.matrix(read.table(paste0('LD/AMR',trait1,'.LD.ld'), header=F, stringsAsFactors=F))
ldeur <-  as.matrix(read.table(paste0('LD/EUR',trait1,'.LD.ld'), header=F, stringsAsFactors=F))

colnames(ldamr) <- amrbim[, "rsid"] # Col 
rownames(ldamr) <- amrbim[, "rsid"] # Row

# eur
colnames(ldeur) <- eurbim[, "rsid"] # Col
rownames(ldeur) <- eurbim[, "rsid"] # Row

miss.amr <- colSums(is.na(ldamr)) 
miss.eur <- colSums(is.na(ldeur)) 

max(miss.eur)
min(miss.eur)

max(miss.amr)
min(miss.amr)

colnames(ldamr)[colSums(is.na(ldamr)) > 0]

colnames(ldeur)[colSums(is.na(ldeur)) > 0]

colnames(ldeur)[colSums(is.na(ldeur)) > 10]

eurbim <- eurbim[eurbim$rsid != "rs78116095",]

nrow(eurbim)

write.table(eurbim[,c(2,5)],                          # Identificamos las columnas de rsID y REF
            paste0('LD/',trait1,'.SNP.REF.txt'), 
            col.names=F, row.names=F, sep='\t', quote=F)

write.table(eurbim[,2],                               # Identificamos la columna de rsID 
            paste0('LD/', trait1, '.SNP.extract'), 
            col.names=F, row.names=F, sep='\t', quote=F)

system(
  paste0(plink,
         ' --vcf ', data_path, '1KG.AMR.GRCH38.rsID.', CHR, '.region.vcf',
         ' --extract LD/', trait1, '.SNP.extract',
         ' --threads 1', 
         ' --update-ref-allele LD/', trait1, '.SNP.REF.txt 2 1',
         ' --r square',
         ' --out LD/AMR', trait1, '.LD', 
         ' --make-just-bim'))

system(
  paste0(plink,
         ' --vcf ', data_path, '1KG.EUR.GRCH38.rsID.', CHR, '.region.vcf', 
         ' --extract LD/', trait1,'.SNP.extract', 
         ' --threads 1',
         ' --update-ref-allele LD/', trait1,'.SNP.REF.txt 2 1', 
         ' --r square',
         ' --out LD/EUR', trait1,'.LD', 
         ' --make-just-bim'))

amrbim <- read.table(paste0('LD/AMR',trait1,'.LD.bim'), header=F)

eurbim <- read.table(paste0('LD/EUR',trait1,'.LD.bim'), header=F)

colnames(amrbim) <- c("chr", "rsid", "cM", "position", "allele1", "allele2")

colnames(eurbim) <- c("chr", "rsid", "cM", "position", "allele1", "allele2")

nrow(eurbim)

nrow(amrbim)

eqtl <- eqtl[eurbim$rsid, ]

gwas <- gwas[amrbim$rsid, ]

dim(gwas)
dim(eqtl)
dim(amrbim)
dim(eurbim)

ldamr <- as.matrix(read.table(paste0('LD/AMR',trait1,'.LD.ld'), header=F, stringsAsFactors=F))

ldeur <- as.matrix(read.table(paste0('LD/EUR',trait1,'.LD.ld'), header=F, stringsAsFactors=F))

colnames(ldamr) <- amrbim[, "rsid"]
rownames(ldamr) <- amrbim[, "rsid"]
colnames(ldeur) <- eurbim[, "rsid"]
rownames(ldeur) <- eurbim[, "rsid"]

colnames(ldamr)[colSums(is.na(ldamr)) > 0]
colnames(ldeur)[colSums(is.na(ldeur)) > 0]

eqtl.d <- 
  list(pvalues = eqtl[eurbim$rsid, 'pval_asc'],
       snp = eqtl[eurbim$rsid, 'rsID'],
       MAF = eqtl[eurbim$rsid, 'maf_trc'], 
       position = eurbim$position, 
       type = 'quant',
       N = eqtl$samples_asc[1])

gwas.d <- 
  list(pvalues = gwas[amrbim$rsid,'p_value'],
       snp = gwas[ amrbim$rsid, 'MarkerName'], 
       MAF = gwas[amrbim$rsid, 'EAF'], 
       position = amrbim$position, 
       type = 'quant',
       N = 50000)

check_dataset(eqtl.d, req = c("snp", "pvalues"))

check_dataset(gwas.d, req = c("snp", "pvalues"))

eqtl$varbeta <- eqtl$beta_se_asc ^ 2

gwas$varbeta <- gwas$se_0 ^ 2

eqtl.d2 <- 
  list(snp = eqtl[eurbim$rsid, 'rsID'], 
       beta = eqtl[eurbim$rsid, 'beta_asc'],
       varbeta = eqtl[eurbim$rsid, 'varbeta'], 
       position = eurbim$position, 
       type = 'quant', 
       N = eqtl$samples_asc[1],
       LD = ldeur, 
       sdY = 1, 
       MAF = eqtl[eurbim$rsid, 'maf_trc'])

gwas.d2 <- 
  list(snp = gwas[amrbim$rsid, 'MarkerName'], 
       beta = gwas[amrbim$rsid, 'beta_0'], 
       varbeta = gwas[amrbim$rsid, 'varbeta'], 
       position = amrbim$position, 
       type = 'quant', 
       N = 50000, 
       LD = ldamr, 
       sdY=1, 
       MAF = gwas[amrbim$rsid, 'EAF'])

check_dataset(eqtl.d2,req = c("snp", "varbeta", "beta", "MAF"))
check_dataset(gwas.d2,req = c("snp", "varbeta", "beta", "MAF"))

my.res <- coloc.abf(dataset1 = gwas.d, dataset2 = eqtl.d)

print(my.res)

plot(my.res)

my.res.bvar <- coloc.abf(dataset1 = gwas.d2, dataset2 = eqtl.d2)

print(my.res.bvar)

plot(my.res.bvar)

print(subset(my.res$results, SNP.PP.H4 > 0.1))

print(subset(my.res.bvar$results, SNP.PP.H4>0.1))

eqtl.s <- runsusie(eqtl.d2)

gwas.s <- runsusie(gwas.d2)

plot_dataset(gwas.d2, susie_obj = gwas.s, color = c("green4"))

plot_dataset(eqtl.d2, susie_obj = eqtl.s, color = c("yellow4"))

susie.res <- coloc.susie(gwas.s, eqtl.s)

susie.res$summary

print(subset(susie.res$results, SNP.PP.H4.abf > 0.05))

sensitivity(susie.res, "H4 > 0.9", row=1, dataset1 = eqtl.d2, dataset2 = gwas.d2)

if(my.res$summary['PP.H4.abf'] == max(my.res$summary[-1])) { 
  
  # En caso de que se cumpla el criterio, guardamos los resultados. 
  write.table(my.res$summary,
              paste0(trait1, '.', trait2, '.coloc_abf_gene.pval.txt'),
              col.names=T, row.names=F, sep='\t', quote=F)
  
  # Guardamos todos los valores cuyo SNP.PP.H4 sea mayor a 0.01.  
  write.table(subset(my.res$results,SNP.PP.H4>0.01),
              paste0(trait1, '.', trait2, '.coloc_abf_snp.pval.txt'), 
              col.names=T, row.names=F, sep='\t', quote=F)
}

# Generamos un archivo `.png` donde se guarde el grpafico. 
png("Height.eqtl.coloc_pvalues.png")
plot(my.res)
dev.off()

if(my.res.bvar$summary['PP.H4.abf'] == max(my.res.bvar$summary[-1])) {
  write.table(
    my.res.bvar$summary, 
    paste0(trait1, '.', trait2, '.coloc_abf_gene.betas.txt'), 
    col.names=T, row.names=F, sep='\t', quote=F)
  
  write.table(
    subset(my.res.bvar$results, SNP.PP.H4 > 0.01),
    paste0(trait1, '.', trait2, '.coloc_abf_snp.betas.txt'), 
    col.names=T, row.names=F, sep='\t', quote=F)
}

png("Height.eqtl.coloc_betas.png")
plot(my.res.bvar)
dev.off()

write.table(susie.res$summary, 
            paste0(trait1, '.', trait2, '.coloc_susie_abf_gene.txt'),
            col.names=T, row.names=F, sep='\t', quote=F)

write.table(subset(susie.res$results, SNP.PP.H4.abf > 0.01),
            paste0(trait1, '.', trait2, '.coloc_susie_abf_snp.txt'),
            col.names=T, row.names=F, sep='\t', quote=F)

png("Height.gwas.susie_obj.png")
plot_dataset(gwas.d2,susie_obj=gwas.s,color = c("green4"))
dev.off()

png("Height.eqtl.susie_obj.png")
plot_dataset(eqtl.d2,susie_obj=eqtl.s,color = c("yellow4"))
dev.off()

png("Height.gwas.eqtl.susie_probs.png")
sensitivity(susie.res,"H4 > 0.9", 
            row=1, 
            dataset1 = eqtl.d2, 
            dataset2 = gwas.d2)
dev.off()
