
dimplot_new <- function(data = datafilt,
                        reduction = "umap",
                        pt.size = 1, label = T,
                        group.by = c("seurat_clusters")) {
  
  library(Seurat)
  
  # 设定颜色
  
  info <- data@meta.data
  number <- length(unique(info[,group.by]))
  
  if (number < 17) {
    col <- c("#A6D719", "#176EBF", "#00A8DE", "#AEE0E8",
             "#00A9A3", "#FBD324", "#F28A24", "#A52828",
             "#A37CB7", "#F2D7EE", "#CD6981", "#FBD0C0",
             "#F15E4C", "#ECB2C8", "#B2DBBF", "#CCDAD1")
  } else {
    col <- c("#B8E3EA", "#5CB3DA", "#0070B2", "#FBDD7E", "#F7AE24", "#FF7149", 
             "#F2D7EE", "#A37CB7", "#A231A1", "#ECB2C8", "#E93B8C", "#B91372", 
             "#FF9F99", "#F15E4C", "#DA1735", "#CDE391", "#8BBE53", "#679436", 
             "#98D4C6", "#00A385", "#067D69", "#B2DBBF", "#028090", "#114B5F", 
             "#FBD0C0", "#CD6981", "#A23E48", "#CCDAD1", "#9CAEA9", "#788585")
    col <- colorRampPalette(col)(number)
  }
  
  # 开始画图
  
  DimPlot(data, pt.size = pt.size, label = label, repel = T, 
          raster = FALSE, label.size = 3.5, reduction = reduction,
          group.by = group.by) + 
    scale_color_manual(values = c(col)) + 
    
    theme_bw() +
    theme(legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid = element_blank(), # 删去网格线
          axis.ticks = element_blank(), # 删去刻度线
          axis.text = element_blank(), # 删去刻度标签
          axis.title = element_text(colour = "black", size = 15),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +   
    
    labs(x = 'UMAP1',y= 'UMAP2',title = '')
}