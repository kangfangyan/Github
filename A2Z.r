pkgs = c("GEOquery", 'beadarray', 'illuminaHumanv4.db', 'limma')
sapply(X = pkgs, library, character.only = TRUE)
# source('limma.autoDE.r')

gse <- getGEO(destdir = '~/WIP/GSE104237/',
              filename = "GSE104237_series_matrix.txt.gz")
pData(gse)$group <- c('ctl', 'ctl', 'trt', 'trt')

dim(exprs(gse))
head(exprs(gse)); tail(exprs(gse))
summaryData <- as(gse, "ExpressionSetIllumina")
# summaryData
# head(fData(summaryData))

fData(summaryData)$Status <- ifelse(fData(summaryData)$PROBEQUALITY=="No match",
                                    "negative","regular" )
Detection(summaryData) <- calculateDetection(summaryData, status=fData(summaryData)$Status)

summaryData.norm <- normaliseIllumina(summaryData, status=fData(summaryData)$Status,
                                      method="quantile", transform="log2")

# retrive eset data.
eset.exprs = as.data.frame(exprs(summaryData.norm)) #get exprs
is.na(eset.exprs) = do.call(cbind, lapply(eset.exprs, is.infinite)) #remove inf values if any.
eset.exprs = as.matrix(eset.exprs[complete.cases(eset.exprs),])

phenoData = new(Class = 'AnnotatedDataFrame',data = pData(gse)) # create new pData
eset = ExpressionSet(assayData = as.matrix(eset.exprs),
                     phenoData = phenoData, annotation = 'Humanv4') #create new expressionSet object
eset = addFeatureData(eset, toAdd = c("SYMBOL", "ENSEMBL", "ENTREZID"), annotation = 'Humanv4') #add other features from IlluminaV4 pacakge.


exprs.df = cbind(exprs(eset),as(eset@featureData,Class = 'data.frame'))
exprs.df = exprs.df[,-grep(pattern = 'Row.names',x = colnames(exprs.df))]

table(is.na(exprs.df$SYMBOL))
table(is.na(exprs.df$ENSEMBL))
table(is.na(exprs.df$ENTREZID))
# filter probes that are not expressed.

# eset <- eset.exprs

# filter GSE --------------------------------------------------------------

gse <- gse[!is.na(exprs.df$SYMBOL),]

# phenoData = new(Class = 'AnnotatedDataFrame',data = pData(gse)) # create new pData

summaryData <- as(gse, "ExpressionSetIllumina")
summaryData
fData(summaryData)$Status <- ifelse(fData(summaryData)$PROBEQUALITY=="No match",
                                    "negative","regular" )
Detection(summaryData) <- calculateDetection(summaryData, status=fData(summaryData)$Status)

summaryData.norm <- normaliseIllumina(summaryData, status=fData(summaryData)$Status,
                                      method="quantile", transform="log2")

# retrive eset data.
eset.exprs = as.data.frame(exprs(summaryData.norm)) #get exprs
# is.na(eset.exprs) = do.call(cbind, lapply(eset.exprs, is.infinite)) #remove inf values if any.
eset.exprs = as.matrix(eset.exprs[complete.cases(eset.exprs),])

phenoData = new(Class = 'AnnotatedDataFrame',data = pData(gse)) # create new pData
eset = ExpressionSet(assayData = as.matrix(eset.exprs),
                     phenoData = phenoData, annotation = 'Humanv4') #create new expressionSet object
eset = addFeatureData(eset, toAdd = c("SYMBOL", "ENSEMBL", "ENTREZID"), annotation = 'Humanv4') #add other features from IlluminaV4 pacakge.


exprs.df = cbind(exprs(eset),as(eset@featureData,Class = 'data.frame'))
exprs.df = exprs.df[,-grep(pattern = 'Row.names',x = colnames(exprs.df))]

table(is.na(exprs.df$SYMBOL))
table(is.na(exprs.df$ENSEMBL))
table(is.na(exprs.df$ENTREZID))












# eset = ExpressionSet(assayData = as.matrix(exprs(gse)),
#                      phenoData = phenoData, annotation = 'Humanv4') #create new expressionSet object
# eset = addFeatureData(eset, toAdd = c("SYMBOL", "ENSEMBL", "ENTREZID"), annotation = 'Humanv4')
# # exprs(gse) <- exprs(gse)[!is.na(exprs.df$SYMBOL),]
# exprs.df = cbind(exprs(eset),as(eset@featureData,Class = 'data.frame'))
# exprs.df = exprs.df[,-grep(pattern = 'Row.names',x = colnames(exprs.df))]
#
# table(is.na(exprs.df$SYMBOL))
# table(is.na(exprs.df$ENSEMBL))
# table(is.na(exprs.df$ENTREZID))

# DEG analysis begins here.
fdr <- 0.05

# design <- model.matrix(~0 + as.factor(eset$group))
# colnames(design) <- as.character(levels(as.factor(eset$group)))
# wts <- arrayWeights(exprs(eset), design)
#
# fit <- lmFit(exprs(eset), design, weights = wts)
#
# cnt <- paste(colnames(design)[2], colnames(design)[1], sep = "-")
# cMat <- makeContrasts(contrasts = cnt, levels = design)
# fit2 <- contrasts.fit(fit, cMat)
# efit <- eBayes(fit2)
#
# tt = topTable(efit, number = 'all')  #get toptable.
# # res = merge(exprs.df, tt, by = 'row.names') #merge everything.
# res = tt[order(tt$adj.P.Val),] #order according to fdr.
# rownames(res) = res$Row.names
# res = res[,-1]


if(length(pData(eset)[,1]) == 2){
    cat('no replicates available; can not perform differntial expression analysis \n returning just expression table')
    return(exprs.df)
} else{
    eset.de = limma.autoDE(summaryData = eset, SampleGroup = 'group') #run differnetial analysis.
    tt = topTable(eset.de, number = 'all')  #get toptable.
    res = merge(exprs.df, tt, by = 'row.names') #merge everything.
    res = res[order(res$adj.P.Val),] #order according to fdr.
    rownames(res) = res$Row.names
    res = res[,-1]

    #get up and down probes
    downProbes = rownames(res[which(res$adj.P.Val < fdr & res$logFC < 0),])
    upProbes = rownames(res[which(res$adj.P.Val < fdr & res$logFC > 0),])

    #annotate deg status
    res$down = as.integer(rownames(res) %in% downProbes)
    res$up = as.integer(rownames(res) %in% upProbes)
    res$up = factor(x = res$up, levels = c(1,0), labels = c('up', 'none'))
    res$down = factor(x = res$down, levels = c(1,0), labels = c('down', 'none'))
    res$deg = interaction(res$up, res$down)
    res$deg = as.character(factor(x = res$deg, levels = c("up.down"  ,"none.down","up.none"  ,"none.none"), labels = c('none', 'down', 'up', 'none')))
    res = subset(res, select = -c(up, down))

    #Perform Principal Component Analysis (PCA) and plot first two Principal Components
    # if(plotPCA){
    #     pca_plot(mat = res[,names], pdata = pData(eset), label = T)
    # }
    #
    #return(list(results = res, eset = eset))
    return(res)
}


# res <- res[!is.na(res$SYMBOL),]
library(dplyr)

res.sig <- res %>%
    filter(!is.na(SYMBOL)) %>%
    # filter(adj.P.Val < 0.05) %>%
    filter(P.Value < 0.05) %>%
    arrange(P.Value, logFC)

save(res.sig, file = "res.sig.rda")

load('res.sig.rda')
library(clusterProfiler)

go <- enrichGO(res.sig$ENTREZID, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID', ont = 'BP',
               qvalueCutoff = 0.05, minGSSize = 3)
go.mf <- enrichGO(res.sig$ENTREZID, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID', ont = 'MF',
               qvalueCutoff = 0.05, minGSSize = 3)
head(summary(go.mf))
head(summary(go))

dotplot(go.mf, showCategory=30)

keg <- enrichKEGG(res.sig$ENTREZID, organism = 'hsa', qvalueCutoff = 0.05, minGSSize = 3)
# yy <- gseMKEGG(res.sig$ENTREZID)
summary(keg)
