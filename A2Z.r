# 0. setup ----------------------------------------------------------------
options(stringsAsFactors = FALSE)
pkgs = c("GEOquery", 'beadarray', 'illuminaHumanv4.db',
         'limma', "dplyr", 'clusterProfiler')
sapply(X = pkgs, library, character.only = TRUE)
# we can of course get rid of limma.autoDE.r, but I am not gonna....
source('limma.autoDE.r')


# 1. read data in ---------------------------------------------------------

gse <- getGEO(destdir = '.',
              filename = "GSE104237_series_matrix.txt.gz")
pData(gse)$group <- c('ctl', 'ctl', 'trt', 'trt')

dim(exprs(gse))
head(exprs(gse)); tail(exprs(gse))


# 2. summaryData ----------------------------------------------------------
summaryData <- as(gse, "ExpressionSetIllumina")
fData(summaryData)$status <- if_else(fData(summaryData)$PROBEQUALITY %in% c("No match", "Bad") |
                                        is.na(fData(summaryData)$PROBEQUALITY),
                                    "negative",
                                    "regular")

# probes with low quality or neg values are then removed
rm.neg.idx <- fData(summaryData)$status == "negative" | rowMin(exprs(summaryData)) <= 0
summaryData <- summaryData[!rm.neg.idx, ]

# normalization
summaryData.norm <- normaliseIllumina(summaryData, status=fData(summaryData)$status,
                                      method="quantile", transform="log2")


# 3. eSet data ------------------------------------------------------------

# retrive eset data.
eset <- as.matrix(exprs(summaryData.norm))
phenoData <- new(Class = 'AnnotatedDataFrame', data = pData(gse)) # create new pData
eset <- ExpressionSet(assayData = eset,
                      phenoData = phenoData,
                      annotation = 'Humanv4') #create new expressionSet object
eset <- addFeatureData(eset, toAdd = c("SYMBOL", "ENTREZID"), annotation = 'Humanv4') #add other features from IlluminaV4 pacakge.


eset.df <- cbind(exprs(eset), as(eset@featureData, Class = 'data.frame'))
eset.df <- eset.df[, -grep(pattern = 'Row.names', x = colnames(eset.df))]

table(is.na(eset.df$SYMBOL))
table(is.na(eset.df$ENTREZID))
sapply(eset.df, class)   # options(stringsAsFactors = FALSE) REALLY MAKES A DIFFERENCE HERE!

# 4. filter eset --------------------------------------------------------------

eset.df$gene.exp.median <- apply(eset.df[, 1:4], 1, median)
eset.df <- eset.df[order(eset.df$SYMBOL,
                         eset.df$gene.exp.median,
                         decreasing = TRUE, na.last = TRUE),]
# we used to use gse instead of eset here, however, gse and eset.df are not the same.
# Recall that gse contains raw data, while probes with low quality or neg exprs values
# were removed in eset.df, rsulting a diffrence between gse and eset.df.
# Basically, here we have to make sure that data are IDENTICAL in BOTH order and number of rows if
# we are to filter 2 different data.
# But then again, only eset is needed to perform DE analysis using limma.autoDE(). So here, we
# no longer consider gse, summaryData or summaryData.norm anymore. Instead, we filter eset.df first and
# then reconstruct eset data with the new eset.df, which is exprs(eset).

# we thought we could filter eset based on the order we get from eset.df, but we were
# wrong, sadly. Row order of an ExpressionSet object (eset) is not changed with the
# following code:
# eset <- eset[order(eset.df$SYMBOL,
#                    eset.df$gene.exp.median,
#                    decreasing = TRUE, na.last = TRUE),]
# So, we have to come up with a different solution.
# the above code is kept as a reminder, though.

no.gene.idx <- which(is.na(eset.df$SYMBOL) |
                       is.na(eset.df$ENTREZID))
eset.df <- eset.df[-no.gene.idx,]
dup.idx <- which(duplicated(eset.df$SYMBOL))
eset.df <- eset.df[-dup.idx,]
# I am always confused here what will happen to the numerous NA duplicated rows
# in our previous code where we deal with NA and duplicated SYMBOL together.
# Later it occurred to me that if we filter out NA first and then come to dup
# SYMBOL, then the problem is safely avoided.
# The code is more readable, more understandable, only at the cost of
# being a little redundant, which is tatolly acceptable for me anyway.


# 5.  new eSet data -------------------------------------------------------

# reconstruct new eSet object from the filtered eset.df data
eset <- ExpressionSet(assayData = as.matrix(eset.df[, 1:4]),
                      phenoData = phenoData,
                      annotation = 'Humanv4') #create new expressionSet object
eset <- addFeatureData(eset, toAdd = c("SYMBOL", "ENTREZID"), annotation = 'Humanv4')


# 6. and then DE analysis ----------------------------------------------------------

# DEG analysis begins here.
fdr <- 0.05

eset.de <- limma.autoDE(summaryData = eset, SampleGroup = 'group') # run differnetial analysis.
tt <- topTable(eset.de, number = 'all')  # get toptable.
res <- merge(eset.df, tt, by = 'row.names') # merge everything.
res <- res[order(res$adj.P.Val),] # order according to fdr.


res.sig <- res %>%
    # filter(!is.na(SYMBOL)) %>%
    filter(adj.P.Val < 0.05) %>%
    # filter(P.Value < 0.05) %>%
    arrange(adj.P.Val, logFC)

# double check for NA or dup gene symbols in final result
table(is.na(res.sig$SYMBOL))
table(duplicated(res.sig$SYMBOL))

# write out res.sig if you wanna
write.table(res.sig, file = 'GSE104237_DEG_fdr005.csv',
            sep = '\t', quote = FALSE, row.names = FALSE)

# 7. just a GO and KEGG eg -------------------------------------------------------------

go.bp <- enrichGO(res.sig$ENTREZID, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID', ont = 'BP',
               qvalueCutoff = 0.05, minGSSize = 3)
go.mf <- enrichGO(res.sig$ENTREZID, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID', ont = 'MF',
                  qvalueCutoff = 0.05, minGSSize = 3)
head(summary(go.bp))
head(summary(go.mf))

dotplot(go.bp, showCategory=30)
dotplot(go.mf, showCategory=30)

kegg <- enrichKEGG(res.sig$ENTREZID, organism = 'hsa', qvalueCutoff = 0.05, minGSSize = 3)

summary(kegg)
