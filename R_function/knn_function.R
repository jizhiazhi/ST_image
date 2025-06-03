library(RANN)
#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#round is the average round
#knn is the n nearest point from reference
smooth_kNN <- function(initial_df,query_df,sm_vector,round=100,knn=30){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    for(i in 1:round){
        sm_vector <- rowMeans(matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE))
    }
    return(sm_vector)
}

smooth_kNN2 <- function(initial_df,query_df,compute='mean',sm_vector,round=100,knn=30){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    Nx <- nrow(initial_df)
    knnIdx <- result$nn.idx
    if(compute=='mean'){
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), knn), c(knnIdx), x=1, dims = c(Nx, Nx))
        W <- W/rowSums(W)
    }else if(compute=='gauss'){
        knnDist <- result$nn.dists
        knnDist <- knnDist / knnDist[,2]
        knnDist <- knnDist**2
        W <- Matrix::sparseMatrix(rep(seq_len(Nx), knn), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
        W@x <- exp(-(W@x))
    }
    for(i in 1:round){
        sm_vector <- as.matrix(W %*% as.matrix(sm_vector))
    }
    return(sm_vector)
}


#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#knn is the n nearest point from reference
#If the value counts exceed the threshold, specific actions apply; otherwise, use the nearest point as the query value.
winner_kNN <- function(initial_df,query_df,sm_vector,knn=5,threhold=1){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    sm_vector <- matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE)
    sm_vector <- sapply(1:nrow(sm_vector),function(x){
        tmp_table <- table(sm_vector[x,])
        name <- names(sort(tmp_table,decreasing = T))[1]
        if(tmp_table[name]>threhold){return(name)}else{return(sm_vector[x,1])}
    })
    return(sm_vector_sm)
}
magic_knn <- function(initial_df,query_df,sm_vector,round=3,knn=15,ka=3,epsilon=1){
    Nx <- nrow(initial_df)
    result <- RANN::nn2(initial_df,query_df,k=knn)
    knnIdx <- result$nn.idx
    knnDist <- result$nn.dists
    
    knnDist <- knnDist / knnDist[,ka]
    #Markov compute
    W <- Matrix::sparseMatrix(rep(seq_len(Nx), knn), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
    W <- W + Matrix::t(W)
    #Compute Kernel
    W@x <- exp(-(W@x / epsilon^2))
    #Markov normalization
    W <- W / Matrix::rowSums(W) 
    #Initialize Matrix
    Wt <- W
    #Computing Diffusion Matrix
    for(i in seq_len(round)){
      Wt <- Wt %*% W
    }
    sm_vector_sm <- as.matrix(Wt %*% sm_vector)
    return(sm_vector_sm)
}
ROC_kNN <- function (initial_df, vector, knn = 9) 
{
    result <- RANN::nn2(initial_df, initial_df, k = knn)
    vector_n <- matrix(vector[result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    x_n <- matrix(initial_df[, 1][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    y_n <- matrix(initial_df[, 2][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    vector_max <- matrixStats::rowMaxs(vector_n)
    vector_min <- matrixStats::rowMins(vector_n)
    gradient_top <- abs(vector_max - vector_min)
    gradient_bottom <- sapply(1:nrow(vector_n), function(i) {
        x_max <- x_n[i, ][vector_n[i, ] == vector_max[i]][1]
        y_max <- y_n[i, ][vector_n[i, ] == vector_max[i]][1]
        x_min <- x_n[i, ][vector_n[i, ] == vector_min[i]][1]
        y_min <- y_n[i, ][vector_n[i, ] == vector_min[i]][1]
        distance <- ((x_max - x_min)^2 + (y_max - y_min)^2)^0.5
        return(distance)
    })
    result <- gradient_top/gradient_bottom
    result[is.na(result)] <- 0
    result[is.infinite(result)] <- 0
    result <- (result-min(result))/(max(result)-min(result))
    result[is.na(result)] <- 0
    return(result)
}

smooth_obj_data <- function(obj,smooth_function='spatial2',compute='mean'){
    if(smooth_function=='spatial'){
    data_matrix <- as.matrix(obj@assays$RNA@data)
    pb <- txtProgressBar(style=3)
    n <- 1
    all_n <- nrow(data_matrix)
    for(i in rownames(data_matrix)){
    data_matrix[i,] <- smooth_kNN(obj@meta.data[,c('x','y')],obj@meta.data[,c('x','y')],data_matrix[i,],knn=16,round=10)
    setTxtProgressBar(pb, n/all_n)
    n <- n+1}
    close(pb)
    
    
    }else if(smooth_function=='spatial2'){
    data_matrix <- Matrix::t((obj@assays$RNA@data))
    data_matrix <- smooth_kNN2(obj@meta.data[,c('x','y')],obj@meta.data[,c('x','y')],compute=compute,data_matrix,knn=9,round=2)
    data_matrix <- t(data_matrix)
    colnames(data_matrix) <- colnames(obj)
    
    
    }else if(smooth_function=='magic'){
    
    obj <- FindVariableFeatures(obj,nfeatures = 2000)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj,verbose=FALSE)
    
    data_matrix <- magic_knn(obj@reductions$pca@cell.embeddings,
    obj@reductions$pca@cell.embeddings,
    Matrix::t(obj@assays$RNA@data),
    knn=15,
    round=2)
    data_matrix <- t(data_matrix)
    colnames(data_matrix) <- colnames(obj)
    
    }
    data_matrix <- CreateAssayObject(data=data_matrix)
    obj[['RNA']]@data <- data_matrix@data
    return(obj)
}
