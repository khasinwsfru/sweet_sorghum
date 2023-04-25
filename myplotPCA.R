#this is just a copy of plotPCA that swaps out PCs 1 and 2 for 2 and 3
myplotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 
        100), "% variance")) + coord_fixed()
}
