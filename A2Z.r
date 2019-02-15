pkgs = c("GEOquery", 'beadarray', 'illuminaHumanv4.db',
         'limma', "dplyr", 'clusterProfiler')
sapply(X = pkgs, library, character.only = TRUE)
source('limma.autoDE.r')

gse <- getGEO(destdir = '.',
              filename = "GSE104237_series_matrix.txt")
pData(gse)$group <- c('ctl', 'ctl', 'trt', 'trt')

dim(exprs(gse))
head(exprs(gse)); tail(exprs(gse))
#将下载的数据转换为ExpressionSetIllumina
summaryData <- as(gse, "ExpressionSetIllumina")
# summaryData
# head(fData(summaryData))
#去除非匹配
fData(summaryData)$status <- ifelse(fData(summaryData)$PROBEQUALITY == "No match",
                                    "negative", "regular")
Detection(summaryData) <- calculateDetection(summaryData,
                                             status=fData(summaryData)$Status,
                                             negativeLabel="negative")
#normalization
summaryData.norm <- normaliseIllumina(summaryData, status=fData(summaryData)$Status,
                                      method="quantile", transform="log2")

# retrive eset data.
eset <- as.data.frame(exprs(summaryData.norm)) #get exprs
is.na(eset) = do.call(cbind, lapply(eset, is.infinite)) #remove inf values if any.
eset <- as.matrix(eset[complete.cases(eset),])

phenoData <- new(Class = 'AnnotatedDataFrame', data = pData(gse)) # create new pData
eset <- ExpressionSet(assayData = eset,
                     phenoData = phenoData,
                     annotation = 'Humanv4') #create new expressionSet object
eset <- addFeatureData(eset, toAdd = c("SYMBOL", "ENTREZID"), annotation = 'Humanv4') #add other features from IlluminaV4 pacakge.


exprs.df = cbind(exprs(eset), as(eset@featureData, Class = 'data.frame'))
exprs.df = exprs.df[, -grep(pattern = 'Row.names',x = colnames(exprs.df))]

table(is.na(exprs.df$SYMBOL))
# table(is.na(exprs.df$ENSEMBL))
table(is.na(exprs.df$ENTREZID))
# filter probes that are not expressed.
sapply(exprs.df, class)   # options(stringsAsFactors = FALSE) MAKES A DIFFERENCE HERE!

# eset <- eset.exprs

# filter GSE --------------------------------------------------------------

# gse <- gse[!is.na(exprs.df$SYMBOL),]
exprs.df$gene.exp.median <- apply(exprs.df[, 1:4], 1, median)
exprs.df <- exprs.df[order(as.character(exprs.df$SYMBOL),
                           exprs.df$gene.exp.median,
                           decreasing = TRUE, na.last = TRUE),]
head(exprs.df)
tail(exprs.df)
gse <- gse[order(exprs.df$SYMBOL,
                 exprs.df$gene.exp.median,
                 decreasing = TRUE, na.last = TRUE),]
no.id.idx <- which(is.na(exprs.df$SYMBOL) |
                     is.na(exprs.df$ENTREZID))
dup.idx <- which(duplicated(exprs.df$SYMBOL))
remove.idx <- union(no.id.idx, dup.idx)

length(remove.idx)
dim(exprs.df)

gse <- gse[-remove.idx,]
# gse <- gse[!(is.na(exprs.df$SYMBOL | exprs.df$ENTREZID)),]
# gse <- gse
# phenoData = new(Class = 'AnnotatedDataFrame',data = pData(gse)) # create new pData

summaryData <- as(gse, "ExpressionSetIllumina")
summaryData
fData(summaryData)$Status <- ifelse(fData(summaryData)$PROBEQUALITY=="No match",
                                    "negative","regular" )
Detection(summaryData) <- calculateDetection(summaryData, status=fData(summaryData)$Status)

summaryData.norm <- normaliseIllumina(summaryData, status=fData(summaryData)$Status,
                                      method="quantile", transform="log2")

# retrive eset data.
eset = as.data.frame(exprs(summaryData.norm)) #get exprs
is.na(eset) = do.call(cbind, lapply(eset, is.infinite)) #remove inf values if any.
eset = as.matrix(eset[complete.cases(eset),])

phenoData = new(Class = 'AnnotatedDataFrame', data = pData(gse)) # create new pData
eset = ExpressionSet(assayData = as.matrix(eset),
                     phenoData = phenoData, annotation = 'Humanv4') #create new expressionSet object
eset = addFeatureData(eset, toAdd = c("SYMBOL", "ENTREZID"), annotation = 'Humanv4') #add other features from IlluminaV4 pacakge.


exprs.df = cbind(exprs(eset),as(eset@featureData,Class = 'data.frame'))
exprs.df = exprs.df[,-grep(pattern = 'Row.names',x = colnames(exprs.df))]
head(exprs.df)
tail(exprs.df)
sapply(exprs.df, class)

table(is.na(exprs.df$SYMBOL))
# table(is.na(exprs.df$ENSEMBL))
table(is.na(exprs.df$ENTREZID))


# DEG analysis begins here.
fdr <- 0.05

eset = ExpressionSet(assayData = as.matrix(exprs.df[, 1:4]),
                     featureData = new(Class = "AnnotatedDataFrame", exprs.df[, 5:6]),
                     phenoData = phenoData,
                     annotation = 'Humanv4')

# eset$group <- c('ctl', 'ctl', 'trt', 'trt')

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

res.sig <- res %>%
    # filter(!is.na(SYMBOL)) %>%
    # filter(adj.P.Val < 0.05) %>%
    filter(P.Value < 0.05) %>%
    arrange(P.Value, logFC)

table(duplicated(res.sig$SYMBOL))
head(sort(table(res.sig$SYMBOL), decreasing = TRUE))
res.sig$SYMBOL[is.na(res.sig$SYMBOL)]

res.sig$SYMBOL[duplicated(res.sig$SYMBOL)]


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
