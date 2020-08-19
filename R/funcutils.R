#' Get a list of dimension-plots
#'
#' This function returns a list of 3 dimension plots for your single cell data.
#' The first plot is a basic DimPlot output
#' The second plot is a DimPlot with cells.highlighted with the cells provided
#' The third plot is a DimPlot with cells.highlighted with the cells provided
#'
#' @param sro A Seurat Object
#' @param sr A vector of cell barcodes to highlight. These are the barcodes predicted to be doublets by Scrublet
#' @param dd A vector of cell barcodes to highlight. These are the barcodes predicted to be doublets by DoubletDecon
#' @param plt.title A title for the plots
#' @return A list of three plots
#' @export
getplots <- function(sro, sr, dd, plt.title){
  p1 <- DimPlot(sro, order = T, pt.size = 0.5) + NoLegend() + NoAxes() + labs(title = plt.title)  + panel_border(linetype = 1)
  p2 <- DimPlot(sro, cells.highlight = sr, order = T, pt.size = 0.5) + NoLegend() + NoAxes() + labs(title = paste(plt.title, "-Scrublet", sep = ""))  + panel_border(linetype = 1)
  p3 <- DimPlot(sro, cells.highlight = dd, order = T, pt.size = 0.5) + NoLegend() + NoAxes() + labs(title = paste(plt.title, "-DD", sep = ""))  + panel_border(linetype = 1)
  p <- list(p1, p2, p3)
  return(p)
}

#' Extract data to import to Loupe Cell Browser
#'
#' This function extracts data from a given seurat object for import into Loupe Cell Browser. 
#' It generates two csv files - one containing the coordinates for all the cells for the reduction of interest 
#' and one containing meta-data information such as clusters, categories, etc.
#' 
#'
#' @param object A Seurat Object
#' @param reduction The reduction to extract the coordinates for. 'tsne' or 'umap'
#' @param dims The number of dimensions to extract for a give reduction
#' @param metadata The column name in the seurat object metadata table to extract
#' @param keyword A keyword for the filename
#' @param opdir The output directory in which the files will be written
#' @return Two .csv files written to the specified output directory - data4cloupe_keyword.csv contains the reduction coordinates and cluster4cloupe_keyword.csv contains the metadata information.
#' @export
seurat2cloupe <- function(object, reduction, dims, metadata, keyword, opdir){
  
  if (length(x = dims) != 2) stop("'dims' must be a two-length vector")
  
  file1 = paste(opdir, "/data4cloupe_", keyword, ".csv", sep = "")
  file2 = paste(opdir, "/cluster4cloupe_", keyword, ".csv", sep = "")
  
  embed.data = Embeddings(object = object[[reduction]])[, dims]
  embed.data = as.data.frame(embed.data)
  embed.data <- cbind(rownames(embed.data), embed.data)
  rownames(embed.data) <- 1:nrow(embed.data)
  colnames(embed.data) <- c('Barcode', 'UMAP-1', 'UMAP-2')
  
  write.table(embed.data, file1 , row.names = F, col.names = T, sep = ',', quote = F)
  
  cluster.data = cbind(embed.data$Barcode, object$seurat_clusters, object[[metadata]])
  cluster.data = as.data.frame(cluster.data)
  colnames(cluster.data)[1:2] = c("Barcode", "Clusters")
  rownames(cluster.data) = 1:nrow(cluster.data)
  
  write.table(cluster.data, file2 , row.names = F, col.names = T, sep = ',', quote = F)
  
  cat(paste("All files written to ", opdir, sep = ""))
  
  
}