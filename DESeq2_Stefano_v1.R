#DESeq2 - Stefano Rossi

library(DESeq2)

dir <- system.file("extdata", package = "tximportData")
dir
list.files(dir)
samples <- read.table(file.path(dir, "your_sample.txt"), header = TRUE)

txdf <- transcripts(EnsDb.Mmusculus.v79, return.type="DataFrame")#or any other assembly available here https://www.ensembl.org/index.html
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])

files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample", 1:6)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)#ricorda di mettere sempre nomi corti come kallisto_es1
head(txi.kallisto.tsv$counts)

sampleTable <- data.frame(condition = factor(rep(c("condition_1_eg_WT","condition_2_eg_KO"), each=3)))
View(colnames(txi.kallisto.tsv$counts))
rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~condition)
dds$condition <- relevel(dds$condition, ref = "condition_1_eg_WT")#ref = denominator
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2)
resLFC

res$symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


library("BiocParallel")
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
na.omit(resOrdered)
write.table(resOrdered, file = "~/Desktop/DESeq2(no_cutoff).csv", sep = ",", qmethod = "double")#print the whole list of genes

diffExp.sig <- resOrdered[abs(resOrdered$log2FoldChange) >= log2(1.5) & !is.na(resOrdered$padj) & resOrdered$padj <= 0.05, ]
diffExp_DESeq2<- diffExp.sig[ ,c(7,6,2)]
diffExp_DESeq2<-na.omit(diffExp_DESeq2)
write.table(diffExp_DESeq2, file = "~/Desktop/DESeq2(log2>1.5;padj<0.05).csv", sep = ",", qmethod = "double")#print the list of genes with a cutoff applied

sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

View(resOrdered)
#MA Plot
ggmaplot(resLFC, main = expression("MA plot (DESeq2)"),
         fdr = 0.05, fc = 1.5, size = 0.4,
         select.top.method = c("padj", "fc"),
         palette = c("cyan4", "brown4", "snow2"),
         genenames = as.vector(res$symbol),
         legend = "top", top = 20,
         font.label = c("bold.italic", 11), label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

