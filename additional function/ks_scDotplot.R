#plot

ks_scDotplot <- function(obj,#one of seurat obj or Archr obj
                         features,#genes which your want to plot
                         CellTypes,#cell annotaion
                         avgPctMat, #average expression data, generated from ks_scAverExp
                         aspectRatio=NULL, #plot aspect ratio
                         colorSet=F,#Whether to set color for plot
                         gradientColor=NULL,# If colorSet=T, set gradient color for plot
                         sizeLims=NULL,#dot size
                         pal){#color for cell type annotation
  require(ggplot2)
  require(jjAnno)
  dataType = class(obj)
  
  
  feature_ord <- c()
  for (i in 1:length(features)) {
    
    feature_ord <- append(feature_ord, features[[i]])
  }
  
  feature_ord <- unique(feature_ord)
  feature_ord <- rev(feature_ord)
  
  
  if(dataType == "Seurat") {
    
    group_ord <- c()
    for (i in 1:length(CellTypes)) {
      
      group_ord <- append(group_ord, names(features[i]))
    }
    
    group_ord <- unique(group_ord)
  }
  
  
  if(dataType == "ArchRProject") {
    
    a = names(features)
    b = CellTypes
    c = match(b,a)
    c = sort(c)
    
    group_ord <-names(features)[c]
    
  }
  
  avgPctMat$feature <- factor(avgPctMat$feature, levels = feature_ord)
  avgPctMat$group <- factor(avgPctMat$group, levels = group_ord)
  
  
  if(is.null(aspectRatio)){
    aspectRatio <- length(unique(avgPctMat$feature))/length(unique(avgPctMat$group)) 
  }
  
  p <- ggplot(avgPctMat, aes(x=group, y=feature, color=avgExpr, 
                             size=ifelse(pctExpr > 5, pctExpr, NA)))+
    geom_point()+
    xlab("")+
    ylab("")+
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor= element_blank(), 
          plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
          aspect.ratio = aspectRatio,
          axis.text.x = element_text(angle = 90,  
                                     color = 'black',
                                     hjust = 1,
                                     vjust = 0.5,
                                     margin = margin(0.5,0,0,0,'cm')),
          axis.text.y = element_text(face="italic", color = 'black'),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(fill = NA, color = "black", size = 0.7),
          axis.title = element_text(size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.ticks = element_line(),
          legend.key = element_blank(),
          legend.position = "right",
          legend.direction = "vertical",
          legend.key.size= unit(0.5, "cm"),
          legend.spacing = unit(0, "cm"),
          legend.title = element_text())+
    coord_cartesian(clip = 'off')
  
  
  if(!is.null(sizeLims)){
    p <- p + scale_size_continuous(limits=sizeLims)
  }
  
  
  if(dataType == "Seurat"){
    
    names_title = "Relative\nExpression"
    
  }else{
    
    names_title = "Relative\ngene activity"
  }
  
  
  
  if(colorSet==F){
    
    p <- p+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),
                                 guide = guide_colorbar(ticks.colour = "black",
                                                        frame.colour = "black"),
                                 name = names_title)+
      guides(size = guide_legend(title="Percentage\nExpressed"))
    
    
  }else{
    
    p <- p+scale_color_gradientn(colours = gradientColor,
                                 guide = guide_colorbar(ticks.colour = "black",
                                                        frame.colour = "black"),
                                 name = names_title)+
      guides(size = guide_legend(title="Percentage\nExpressed"))
    
  }
  
  
  p = annoPoint(object = p,
                annoPos = 'botomn',
                xPosition = c(1:length(group_ord)),
                ptSize=1,
                yPosition=-0.4,
                pCol = pal)
  
  return(p)
  
}
