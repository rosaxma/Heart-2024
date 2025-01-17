InspectSeurat <- function(seuratObj) {
  print(paste0('Seurat object size: ', format(object.size(seuratObj), units = 'auto')))
  nCells <- ncol(seuratObj)
  if (nrow(seuratObj@meta.data) != nCells) {
    print(paste0('Metadata rows not equal to ncol: ', nCells, ' / ', nrow(seuratObj@meta.data)))
  }
    
  for (assayName in names(seuratObj@assays)) {
    print(paste0('Assay: ', assayName))
    for (slotName in c('counts', 'data', 'scale.data')) {
      dat <- Seurat::GetAssayData(seuratObj, assay = assayName, slot = slotName)
      print(paste0('slot: ', slotName, ', size: ', format(object.size(x = dat), units = 'auto')))
      
      if (!is(dat, 'sparseMatrix')) {
        print(paste0('Non-sparse! Assay: ', assayName, ', slot: ', slotName))
      }
      
      if (ncol(dat) > 0 && ncol(dat) != nCells) {
        print(paste0('Assay cols not equal to ncol(object): ', nCells, ' / ', ncol(dat)))
      }
    }  
  }
  
  print('Reductions:')
  for (reductionName in names(seuratObj@reductions)) {
    print(paste0(reductionName, ', size: ', format(object.size(x = seuratObj@reductions[reductionName]), units = 'auto')))
  }
  
  print('All slots:')
  for (slotName in slotNames(seuratObj)) {
    val <- object.size(slot(seuratObj, slotName))
    if (val > 1000) { 
      print(paste0(slotName, ': ', format(val, units = 'auto')))
    }
  }
  
  if (nrow(seuratObj@meta.data) != nCells) {
    print(paste0('Metadata rows not equal to ncol: ', nCells, ' / ', nrow(seuratObj@meta.data)))
  }
}
