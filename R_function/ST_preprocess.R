gem2obj <- function(gem,bin,geneID='geneID',x='x',y='y',MIDCounts='MIDCounts',region='region'){
    gem$region <- c('210'='Tha','140'='Lgn','70'='Mgn')[as.character(gem$mask)]
    gem$x <- gem[,x]
    gem$y <- gem[,y]
    gem$x <- round(gem$x/bin,0)
    gem$y <- round(gem$y/bin,0)
    gem$label <- paste0(gem$x,"_",gem$y)
    gem$geneID <- gem[,geneID]
    gem$MIDCount <- gem[,MIDCounts]

    gene <- 1:length(unique(gem$geneID))
    names(gene) <- unique(gem$geneID)
    cell <- 1:length(unique(gem$label))
    names(cell) <- unique(gem$label)
    mat1=sparseMatrix(i = gene[gem$geneID],j=cell[ as.character(gem$label) ], x=gem$MIDCount)
    rownames(mat1) <- names(gene)
    colnames(mat1) <- names(cell)
    obj <- CreateSeuratObject(counts = mat1)
    
    obj$x <- sapply(colnames(obj),function(x){return(strsplit(x,'_')[[1]][1])})
    obj$y <- sapply(colnames(obj),function(x){return(strsplit(x,'_')[[1]][2])})
    obj$x <- as.numeric(obj$x)
    obj$y <- as.numeric(obj$y)
    if(region %in% colnames(gem)){
        region_agg <- aggregate(gem$region,by=list('label'=gem$label),FUN=function(x){names(sort(table(x),decreasing = T))[1]})
        rownames(region_agg) <- region_agg$label
        obj@meta.data$region <- region_agg[rownames(obj@meta.data),2]}
    return(obj)
}

trakem_transform <- function(x,y,matrix,matrix_bin,transform_bin){
    scale = matrix(c(transform_bin/matrix_bin,0,0,0,transform_bin/matrix_bin,0,0,0,1),ncol=3,nrow=3)
    matrix_t <- rbind(matrix(c(as.numeric(strsplit(matrix,',')[[1]])),nrow=2),c(0,0,1))
    matrix_t <- solve(scale) %*% matrix_t %*% scale
    coor_matrix <- cbind(x,y,1)
    coor_matrix_t <- t(matrix_t %*% t(coor_matrix))
    rx <- coor_matrix_t[,1]
    ry <- coor_matrix_t[,2]
    return(list('rx'=rx, 'ry'=ry))
    }
