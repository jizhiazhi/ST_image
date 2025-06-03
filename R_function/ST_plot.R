plot_spatial_feature <- function(obj,features,assays='RNA',slot='data',
                                 max.cutoff=1,min.cutoff=0,vmid=FALSE,
                                 smooth_function='spatial',
                                 smooth=FALSE,knn=25,round=2,
                                 border=FALSE,border_use=Border,
                                 size=3,height=FALSE,width=FALSE,
                                 color=c('gray50','gray','gray97','red','darkred')){
    
    plot.list <- lapply(features,function(gene){
        if(smooth){
            if(smooth_function=='spatial'){
                exp <- smooth_kNN(obj@meta.data[,c('x','y')],
                                  obj@meta.data[,c('x','y')],
                                  slot(obj@assays[[assays]],slot)[gene,],
                                  round=round,knn=knn)
            }else if(smooth_function=='magic'){
                exp <- magic_knn(obj@reductions$pca@cell.embeddings,
                                 obj@reductions$pca@cell.embeddings,
                                 t(as.matrix(slot(obj@assays[[assays]],slot)[gene,,drop=FALSE])),
                                 knn=knn,
                                 round=round)
            }
        }else{exp <- slot(obj@assays[[assays]],slot)[gene,]}

        if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
        obj$exp <- exp

        xmid <- (max(obj$x)+min(obj$x))/2
        ymid <- (max(obj$y)+min(obj$y))/2
        if(height==FALSE){height <- max(obj$y)-min(obj$y)+10}
        if(width==FALSE){width <- max(obj$x)-min(obj$x)+10}

        p1 <- ggplot()+
        scattermore::geom_scattermore(data=obj@meta.data,aes(x=x,y=y,color=exp),pointsize = size)+
        scale_color_gradientn(colours = color,name = gene)+
        theme_void()+
        xlim(xmid-width/2,xmid+width/2)+
        ylim(ymid-height/2,ymid+height/2)+
        theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
        coord_fixed()
        if(vmid!=FALSE){
        p1 <- p1+scale_color_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))))}
        if(border){p1 <- p1+geom_segment(data=border_use,aes(x=X,y=Y,xend=X1,yend=Y1),lwd=0.8)}
        return(p1)
    })
    p <- plot.list[[1]]
    if(length(features)>1){
    for(i in 2:length(features)){p <- p+plot.list[[i]]}
    return(p)
        }else{
    return(p)}
                        }

plot_spatial_col <- function(obj,col_name,border=FALSE,border_use=Border,
                             smooth=FALSE,knn=25,round=2,
                             max.cutoff=1,min.cutoff=0,vmid=FALSE,
                             size=3,height=FALSE,width=FALSE,
                             color=rev(RColorBrewer::brewer.pal(11,'Spectral'))){
    xmid <- (max(obj$x)+min(obj$x))/2
    ymid <- (max(obj$y)+min(obj$y))/2
    exp <- obj[,col_name]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    if(smooth){
        exp <- smooth_kNN(obj[,c('x','y')],
                          obj[,c('x','y')],
                          exp,
                          round=round,knn=knn)}
    if(height==FALSE){height <- max(obj$y)-min(obj$y)+10}
    if(width==FALSE){width <- max(obj$x)-min(obj$x)+10}
    obj$exp <- exp
    p1 <- ggplot()+
    scattermore::geom_scattermore(data=obj,aes(x=x,y=y,color=exp),pointsize = size)+
    scale_color_gradientn(colours = color,name = col_name)+
    xlim(xmid-width/2,xmid+width/2)+
    ylim(ymid-height/2,ymid+height/2)+
    theme_void()+
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
    coord_fixed()
    if(vmid!=FALSE){
    p1 <- p1+scale_color_gradientn(colours = color,name = col_name,
                         values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))))}
    if(border){p1 <- p1+geom_segment(data=border_use,aes(x=X,y=Y,xend=X1,yend=Y1),lwd=0.8)}
    return(p1)
}