################################################################################
################################################################################
#' These functions are used to process the data in the Seurat object,
#' including normalization, batch correction, decontamination, and resegmentation. 
################################################################################
################################################################################


#' Normalize data using scPearsonPCA and run UMAP
#' 
#' This code was adapted from the Nanostring scPearsonPCA package. 
#' (https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/pearsonpca/)
#' 
#' @param obj Seurat object
#' @param assay Assay name to use for normalization
#' @param new.pca.name Name for the new PCA object
#' @param new.umap.name Name for the new UMAP object
#' @param var_features Number of variable features to use (if NULL, all features are used)
#'
#' @return Seurat object with new PCA and UMAP objects added
#' 

scpearson_norm <- function(obj, assay, new.pca.name, new.umap.name, var_features = NULL) {
  
  # Normalization
  if (length(Layers(obj[[assay]])) > 2) {
  obj <- JoinLayers(obj, assay = assay)
  }

  Seurat::DefaultAssay(obj) <- assay
  tc <- Matrix::colSums(obj[[assay]]$counts) ## total counts per cell (across all genes)
  genefreq <- scPearsonPCA::gene_frequency(obj[[assay]]$counts) ## gene frequency (across all cells)
  sum(genefreq)==1  # TRUE
  
  if (is.null(var_features)) {
    use_genes <- rownames(obj[[assay]]@counts)
    } else {
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = var_features)
    use_genes <- Seurat::VariableFeatures(obj)
    }
  
  ### Returns a Seurat-style DimReduc object with 
  ### hvgs x pcs feature loadings
  ### cells x pcs cell embeddings
  ### gene-length vector of the mean pearson of residuals
  ### gene-length vector of the standard deviation of pearson residuals
  pcaobj <- 
    scPearsonPCA::sparse_quasipoisson_pca_seurat(obj[[assay]]$counts[use_genes,]
                                                 ,totalcounts = tc
                                                 ,grate = genefreq[use_genes]
                                                 ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                                                 ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                                                 ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
                                                 ,npcs = 50)
  
  #Add PCA into the seurat object
  obj[[new.pca.name]] <- pcaobj$reduction.data
  
  #Run UMAP
  umapobj <- scPearsonPCA::make_umap(pcaobj)
  
  #Add UMAP into the seurat object.
  obj[[new.umap.name]] <- umapobj$ump
  
  return(obj)
}

#' Normalize data using scPearsonPCA and run UMAP with batch correction
#'
#' This code was adapted from the Nanostring scPearsonPCA package.
#' (https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/pearsonpca/)
#' @param obj Seurat object
#' @param assay Assay name to use for normalization
#' @param new.pca.name Name for the new PCA object
#' @param new.umap.name Name for the new UMAP object
#' @param var_features Number of variable features to use (if NULL, all features are used)
#'
#' @return Seurat object with new PCA and UMAP objects added

scpearson_norm_batch <- function(obj, assay, new.pca.name, new.umap.name, var_features = NULL, batch) {
  
  colnames(obj) <- obj$global_cell_ID
  # Normalization

  # Check if the assay has more than 2 layers, if so, join them into one layer for normalization
  if (length(Layers(obj[[assay]])) > 2) {
    obj <- JoinLayers(obj, assay = assay)
  }
  
  Seurat::DefaultAssay(obj) <- assay

  obs_cols <- c("global_cell_ID", batch)
  
  genefreq_batch <- scPearsonPCA::gene_frequency(obj[[assay]]$counts, obs = data.table(obj@meta.data)[, ..obs_cols]
               ,cellid_colname = "global_cell_ID"
               ,batch_variable = batch)
  Matrix::colSums(genefreq_batch)

  tc <- Matrix::colSums(obj[[assay]]$counts) ## total counts per cell (across all genes)

  if (is.null(var_features)) {
    use_genes <- rownames(obj[[assay]]@counts)
    } else {
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = var_features)
    use_genes <- Seurat::VariableFeatures(obj)
    }

  ### Returns a Seurat-style DimReduc object with 
  ### hvgs x pcs feature loadings
  ### cells x pcs cell embeddings
  ### gene-length vector of the mean pearson of residuals
  ### gene-length vector of the standard deviation of pearson residuals
  pcaobj_batch <- scPearsonPCA::sparse_quasipoisson_pca_seurat_batch(obj[[assay]]$counts[use_genes,]
                                      ,totalcounts = tc
                                      ,grate = genefreq_batch[use_genes,] ##
                                      ,obs =data.table(obj@meta.data)[, ..obs_cols] 
                                      ,batch_variable = batch
                                      ,cellid_colname = "global_cell_ID"
                                      ,scale.max = 10
                                      ,do.scale = TRUE
                                      ,do.center = TRUE
                                      )
    
  #Add PCA into the seurat object
  obj[[new.pca.name]] <-  pcaobj_batch$reduction.data
  
  #Run UMAP
  umapobj <- scPearsonPCA::make_umap(pcaobj_batch)
  
  #Add UMAP into the seurat object.
  obj[[new.umap.name]] <- umapobj$ump
  
  return(obj)
}


#' Decontaminate data using decontX (https://github.com/campbio/decontX)
#' 
#' @param obj Seurat object
#' @param assay Assay name to use for decontamination
#' @param layer Layer name to use for decontamination
#' @param labels Cell labels for decontamination
#' @param new.assay.name Name for the new assay object
#' 
#' @return Seurat object with decontaminated data added. 

decontx_fun <- function(obj, assay, layer, labels, new.assay.name) {

  sce <- Seurat::as.SingleCellExperiment(x = obj, assay = assay)

  #Runnig decontX function
  sce <- decontX::decontX(x = sce, assayName = layer, z = labels)

  #Adding decontaminated counts to the Seurat object
  obj[[new.assay.name]] <- CreateAssayObject(counts = round(decontX::decontXcounts(sce)))
  
  #Adding decontamination metrics to the Seurat object
  obj@meta.data$decontX_contamination <- sce$decontX_contamination

  return(obj)
}


#' Generate a heatmap of marker genes using HieraType from Nanostring
#' 
#' @param obj Seurat object
#' @param assay Assay name to use for heatmap generation
#' @param cluster.name Name of the column in metadata that contains cluster labels
#' @param var_features Number of variable features to use (if NULL, all features are used)
#' @param output_dir Directory to save the output files
#' 
#' @return A heatmap of marker genes using HieraType from Nanostring

hieratype_fun <- function(obj, assay, cluster.name, var_features = NULL) {

  Seurat::DefaultAssay(obj) <- assay
  colnames(obj) <- obj$cell_ID
  
  if (is.null(var_features)) { #Using all genes in the data
    use_genes <- rownames(obj[[assay]]@counts) 
    }
  else if (is.numeric(var_features)) { #Using the top var_features (number) variable features in your data 
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = var_features)
    use_genes <- Seurat::VariableFeatures(obj)
  }
  else if (is.vector(var_features) && length(var_features) > 1) { #Using a list of genes of interest
    use_genes <- unique(var_features[var_features %in% rownames(obj[[assay]]$counts)])
  }

  tc <- Matrix::colSums(obj[[assay]]$counts) ## total counts per cell (across all genes)

  fctbl_rna <- HieraType::clusterwise_foldchange_metrics(obj[[assay]]@counts[use_genes,]
                                            ,totalcounts = tc
                                            ,metadata = obj@meta.data
                                            ,cluster_column = cluster.name
                                            )

  return(fctbl_rna)
}

#' Generate AUCell scores for each cell type based on reference markers
#' 
#' @param obj Seurat object
#' @param ref_markers Data frame of reference markers with columns: gene, cluster, avg_log2FC, p_val_adj
#' @param assay Assay name to use for AUCell score generation
#' @param top_genes Number of top genes to use for each cluster (default is 10)
#' 
#' @return Seurat object with AUCell scores added to metadata
#' 
aucell_score <- function(obj, ref_markers, assay, top_genes = 10) {

  #Removing markers with adjusted p-value = 0 and avg_log2FC <= 0
  ref_markers <- subset(ref_markers, p_val_adj != 0 & avg_log2FC > 0)

  #Get the top 10 markers for each cluster based on avg_log2FC
  ref_topmarkers <- ref_markers[order(ref_markers$avg_log2FC, decreasing = TRUE), ] %>%
    group_by(cluster) %>%
    slice_head(n = top_genes) %>%
    ungroup()

  exprMatrix <- GetAssayData(object = obj, assay = assay, layer = "counts")
  ### Convert to sparse:
  exprMatrix <- as(exprMatrix, "dgCMatrix")

  #Defining gene list
  feature_list <- lapply(seq_along(unique(ref_topmarkers$cluster)), function(i) {
    cluster <- unique(ref_topmarkers$cluster)[i]
    genes  <- ref_topmarkers[ref_topmarkers$cluster == cluster, ]$gene
  })
  names(feature_list) <- unique(ref_topmarkers$cluster)
 
  #Generate AUCell score for each cell type listed on feature_list

  for (i in seq_along(feature_list)) {

  #Define name of geneset, for this case it represents celltype
  geneset_name <- names(feature_list[i])
  message("Analyzing the following geneset: ", geneset_name)

  #Define genesets
  geneSets <- GeneSet(feature_list[[i]], setName = geneset_name)

  #Run AUCell
  cells_AUC <- AUCell_run(exprMatrix, geneSets)

  #Extract to scores and attach to seurat obj metadata
  auc_matrix <- getAUC(cells_AUC)

  meta_colname <- paste0(geneset_name, "_AUCell_score")
  obj[[meta_colname]] <- auc_matrix[geneset_name, colnames(obj)]
  }

  return(obj)
}

#' Adjust tissue coordinates to account for multiple slides
#' 
#' This functions is important if you plan to run a spatial-dependent 
#' analysis across mutliple slides such as despotX. 
#' 
#' @param obj Seurat object
#' 
#' @return Seurat object with adjusted tissue coordinates
#' 
tissue_offset <- function(obj) {

  colnames(obj) <- obj$global_cell_ID

  slidenames <- unique(obj@meta.data$slidename)
  obj$x_slide_mm_offset <- obj$x_slide_mm

  for (slide in 2:length(slidenames)) {
  
    previous_slide <- slidenames[[slide - 1]]
    current_slide <- slidenames[[slide]]
    x_offset <- max(obj$x_slide_mm_offset[obj$slidename == previous_slide]) + 1
    
    obj$x_slide_mm_offset[obj$slidename == current_slide] <- obj$x_slide_mm[obj$slidename == current_slide] + x_offset
    
  }

  #Updating x coordinate on the FOV object
  fov_names <- names(obj@images)
  init_fov <- fov_names[[1]]
  coords <- list(Seurat::GetTissueCoordinates(obj[[init_fov]]))

  for (fov in 2:length(fov_names)) {
    
    current_fov <- fov_names[[fov]]
    x_offset <- max(coords[[fov - 1]]$x) + 1
    
    coords[[fov]] <- Seurat::GetTissueCoordinates(obj[[current_fov]])
    coords[[fov]]$x <- coords[[fov]]$x + x_offset
    
  }

  #Bind coords list  
  coords <- do.call(rbind, coords)
  coords$cell <- colnames(obj)

  # create FOV object
  FOV = SeuratObject::CreateFOV(coords = coords, type = "centroids", assay = "Spatial")

  # add FOV to SeuratObject
  obj@images$all_slides <- FOV # whole_image is the FOV name

  return(obj)
}


#' Function for implementing Nanostring FastReseg to flag cells with potential segmentation errors
#' 
#' @param obj Seurat object
#' @param pixel_size Pixel size in microns (default is 0.12028)
#' @param zstep_size Z-step size in microns (default is 0.8)
#' @param percentCores Percentage of cores to use for processing (default is 0.25)
#' @param input_dir Directory containing input files for FastReseg
#' @param output_dir Directory to save output files from FastReseg
#' @param outDir_flagErrors Directory to save flagged error files from FastReseg
#' 
#' @return Updates the Seurat object metadata with a new column 'reseg_flag' indicating flagged cells

fastreseg_flag <- function(obj, pixel_size = 0.12028, zstep_size = 0.8, percentCores = 0.25, input_dir, output_dir, cluster_col) {


  metadata <- obj@meta.data
  metadata$reseg_flag <- NA

  metadata[[cluster_col]] <- as.character(metadata[[cluster_col]])

  file_paths <- list.dirs(input_dir, recursive = FALSE)
  slidenames <- basename(file_paths)

  for (slide in seq_along(slidenames)) {
    current_path <- file_paths[[slide]]
    current_slide <- slidenames[[slide]]

    message("Evaluating spatial dependency of cells in slide: ", current_slide)

    clust <- metadata[metadata$slidename == current_slide, ][[cluster_col]]
    validCells <- metadata[metadata$slidename == current_slide, ]$global_cell_ID
    
    files <- dir(current_path)

    ########################
    # 1. Loading counts matrix file. Notice that I am using the raw counts for FastReseg
    message("Loading count matrix file...")
    counts_file <- files[grepl("exprMat_file", files)]
    counts <- data.table::fread(file.path(current_path, counts_file))

    cell_ids <- paste0('c_' , slide, '_', counts[['fov']], '_', counts[['cell_ID']])

    # get valid gene names
    all_rnas <- grep("fov|cell_ID|Negative|SystemControl", 
                    colnames(counts), value = TRUE, invert = TRUE)

    counts <- as.matrix(counts[, .SD, .SDcols = all_rnas])
    rownames(counts) <- cell_ids

    #Turn counts into a sparseMatrix
    #counts <- as(counts, "sparseMatrix")
    counts <- as(counts[rownames(counts) %in% validCells, , drop = FALSE], "sparseMatrix")

    if (nrow(counts) != length(clust)) stop("Number of cells in the count matrix do not match number of annotated cells.")
    
    ########################
    # 2. Loading transcript coordinates
    tx_file <- files[grepl("tx_file", files)]
    fullTx <- data.table::fread(file.path(current_path, tx_file))

    tx_cell_id <-  paste0('c_' , slide, '_', fullTx[['fov']], '_', fullTx[['cell_ID']])
    fullTx$global_cell_ID <- tx_cell_id

    # add unique id for each transcript
    fullTx[['transcript_id']] <- seq_len(nrow(fullTx))

    # remove extracellular transcripts which has cell_ID = 0 in tx file 
    fullTx <- fullTx[cell_ID !=0, ]

    # keep only the necessary info
    fullTx <- fullTx[, .SD, .SDcols = c('transcript_id', 'global_cell_ID', 'x_global_px', 
                                        'y_global_px', 'z', 'target', 'fov')]

    # split by FOV and export as per FOV csv file
    outDir_flagErrors <- "res1f_flagErrors"
    txDir <- paste0(output_dir, outDir_flagErrors, "/", current_slide, "/perFOV_txFile")
    dir.create(txDir, recursive = TRUE)

    allFOVs <- unique(fullTx[['fov']])

    transDF_fileInfo <- lapply(allFOVs, function(fovId){
      perFOV_filePath <- fs::path(txDir, paste0('fov_', fovId, '_tx_data.csv'))
      data.table::fwrite(fullTx[fov == fovId, ], file = perFOV_filePath)
      
      # since global coordinates of each molecule are available
      # use 0 for stage coordinates to disable conversion of local to global coordinates
      df <- data.frame(file_path = perFOV_filePath, 
                      slide = 1, 
                      fov = fovId, 
                      stage_X = 0, 
                      stage_Y = 0)
      return(df)
    })

    transDF_fileInfo <- do.call(rbind, transDF_fileInfo)

    flag_dir <- paste0(output_dir, outDir_flagErrors, "/", current_slide)

    message("Flagging cells...")
    flagAll_res <- FastReseg::fastReseg_flag_all_errors(
      counts = counts,
      clust = clust,
      refProfiles = NULL,
      
      # one can use `clust = NULL` if providing `refProfiles`
      
      transcript_df = NULL,
      transDF_fileInfo = transDF_fileInfo,
      filepath_coln = 'file_path',
      prefix_colns = NULL, # to use existing cell IDs that are unique across entire data set 
      fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
      pixel_size = pixel_size, 
      zstep_size = zstep_size,
      transID_coln = 'transcript_id', 
      transGene_coln = "target",
      cellID_coln = "global_cell_ID", 
      spatLocs_colns = c('x_global_px', 'y_global_px', 'z'),
      extracellular_cellID = NULL, 
      
      # control core number used for parallel processing
      percentCores = percentCores, 
      
      # cutoff of transcript number to do spatial modeling
      flagModel_TransNum_cutoff = 50, 
      
      flagCell_lrtest_cutoff = 5, # cutoff for flagging wrongly segmented cells
      svmClass_score_cutoff = -2, # cutoff for low vs. high transcript score
      path_to_output = flag_dir, # path to output folder
      return_trimmed_perCell = TRUE, # flag to return per cell expression matrix after trimming all flagged transcripts 
      ctrl_genes = NULL # optional to include name for control probes in transcript data.frame, e.g. negative control probes
      )

    # extract spatial evaluation outcomes of valid cells
    modStats_ToFlagCells <- flagAll_res[['combined_modStats_ToFlagCells']]
    
    metadata$reseg_flag[metadata$slidename == current_slide] <- modStats_ToFlagCells$flagged[match(metadata$global_cell_ID, modStats_ToFlagCells$UMI_cellID)]
  }

}

#' Function for implementing the full pipeline of Nanostring FastReseg to refine segmentation. 
#' 
#' @param obj Seurat object
#' @param pixel_size Pixel size in microns (default is 0.12028)
#' @param zstep_size Z-step size in microns (default is 0.8)
#' @param percentCores Percentage of cores to use for processing (default is 0.25)
#' @param input_dir Directory containing input files for FastReseg
#' @param output_dir Directory to save output files from FastReseg
#' @param cluster_col Column name in metadata that contains cluster labels
#' 
#' @return Seurat object with refined segmentation results added to metadata and counts assay

fastreseg_fullpipeline <- function(obj, pixel_size = 0.12028, zstep_size = 0.8, percentCores = 0.25, input_dir, output_dir, cluster_col) {

  # Extracting seurat obj metadata. 
  metadata <- obj@meta.data
  
  metadata[[cluster_col]] <- as.character(metadata[[cluster_col]])

  file_paths <- list.dirs(input_dir, recursive = FALSE)
  slidenames <- basename(file_paths)

  for (slide in seq_along(slidenames)) {
    current_path <- file_paths[[slide]]
    current_slide <- slidenames[[slide]]

    message("Evaluating spatial dependency of cells in slide: ", current_slide)

    clust <- metadata[metadata$slidename == current_slide, ][[cluster_col]]
    validCells <- metadata[metadata$slidename == current_slide, ]$global_cell_ID
    
    files <- dir(current_path)

    ########################
    # 1. Loading counts matrix file. Notice that I am using the raw counts for FastReseg
    message("Loading count matrix file...")
    counts_file <- files[grepl("exprMat_file", files)]
    counts <- data.table::fread(file.path(current_path, counts_file))

    cell_ids <- paste0('c_' , slide, '_', counts[['fov']], '_', counts[['cell_ID']])

    # get valid gene names
    all_rnas <- grep("fov|cell_ID|Negative|SystemControl", 
                    colnames(counts), value = TRUE, invert = TRUE)

    counts <- as.matrix(counts[, .SD, .SDcols = all_rnas])
    rownames(counts) <- cell_ids

    #Turn counts into a sparseMatrix
    #counts <- as(counts, "sparseMatrix")
    counts <- as(counts[rownames(counts) %in% validCells, , drop = FALSE], "sparseMatrix")

    if (nrow(counts) != length(clust)) stop("Number of cells in the count matrix does not match number of annotated cells.")
    
    ########################
    # 2. Loading transcript coordinates
    tx_file <- files[grepl("tx_file", files)]
    fullTx <- data.table::fread(file.path(current_path, tx_file))

    tx_cell_id <-  paste0('c_' , slide, '_', fullTx[['fov']], '_', fullTx[['cell_ID']])
    fullTx$global_cell_ID <- tx_cell_id

    # add unique id for each transcript
    fullTx[['transcript_id']] <- seq_len(nrow(fullTx))

    # remove extracellular transcripts which has cell_ID = 0 in tx file 
    fullTx <- fullTx[cell_ID !=0, ]

    # keep only the necessary info
    fullTx <- fullTx[, .SD, .SDcols = c('transcript_id', 'global_cell_ID', 'x_global_px', 
                                        'y_global_px', 'z', 'target', 'fov')]

    
    # split by FOV and export as per FOV csv file
    outDir_full <- "res2_fullPipeline"
    txDir <- paste0(output_dir, outDir_full, "/", current_slide, "/perFOV_txFile")
    dir.create(txDir, recursive = TRUE)

    allFOVs <- unique(fullTx[['fov']])

    transDF_fileInfo <- lapply(allFOVs, function(fovId){
      perFOV_filePath <- fs::path(txDir, paste0('fov_', fovId, '_tx_data.csv'))
      data.table::fwrite(fullTx[fov == fovId, ], file = perFOV_filePath)
      
      # since global coordinates of each molecule are available
      # use 0 for stage coordinates to disable conversion of local to global coordinates
      df <- data.frame(file_path = perFOV_filePath, 
                      slide = 1, 
                      fov = fovId, 
                      stage_X = 0, 
                      stage_Y = 0)
      return(df)
    })

    transDF_fileInfo <- do.call(rbind, transDF_fileInfo)

    fastreseg_dir <- paste0(output_dir, outDir_full, "/", current_slide)

    message("Running full FastReseg pipeline...")
    refineAll_res <- FastReseg::fastReseg_full_pipeline(
    counts = counts,
    clust = clust,
    refProfiles = NULL,
    
    # one can use `clust = NULL` if providing `refProfiles`
    
    transcript_df = NULL,
    transDF_fileInfo = transDF_fileInfo,
    filepath_coln = 'file_path',
    prefix_colns = NULL, # to use existing cell IDs that are unique across entire data set 
    fovOffset_colns = c('stage_Y','stage_X'),
    pixel_size = pixel_size,
    zstep_size = zstep_size,
    transID_coln = 'transcript_id',
    transGene_coln = "target",
    cellID_coln = "global_cell_ID",
    spatLocs_colns = c('x_global_px', 'y_global_px', 'z'),
    extracellular_cellID = NULL,
    
    # control core number used for parallel processing
    percentCores = percentCores, 
    
    # cutoff of transcript number to do spatial modeling
    flagModel_TransNum_cutoff = 50, 
    
    # Optionally, one can set various cutoffs to NULL for automatic calculation from input data
    # Refer to `FastReseg::runPreprocess()` for more details
    
    # distance cutoff for neighborhood searching at molecular and cellular levels, respectively
    molecular_distance_cutoff = 2.7, # 2.7um is recommended for CosMx RNA dataset
    cellular_distance_cutoff = NULL, 
    
    # cutoffs for transcript scores and number for cells under each cell type
    score_baseline = NULL,
    lowerCutoff_transNum = NULL,
    higherCutoff_transNum= NULL,
    imputeFlag_missingCTs = TRUE,
    
    # Settings for error detection and correction, refer to `FastReseg::runSegRefinement()` for more details
    flagCell_lrtest_cutoff = 5, # cutoff to flag for cells with strong spatial dependency in transcript score profiles
    svmClass_score_cutoff = -2,   # cutoff of transcript score to separate between high and low score classes
    groupTranscripts_method = "dbscan",
    spatialMergeCheck_method = "leidenCut", 
    cutoff_spatialMerge = 0.5, # spatial constraint cutoff for a valid merge event
    
    path_to_output = fastreseg_dir,
    save_intermediates = TRUE, # flag to return and write intermediate results to disk
    return_perCellData = TRUE, # flag to return per cell level outputs from updated segmentation 
    combine_extra = FALSE # flag to include trimmed and extracellular transcripts in the exported `updated_transDF.csv` files 
  )
  }

  #Update the Seurat object metadata with the refined segmentation results
  fastreseg_dir <- paste0(output_dir, "res2_fullPipeline/")

  file_paths <- list.dirs(fastreseg_dir, recursive = FALSE)
  slidenames <- basename(file_paths)

  fastreseg_list <- lapply(seq_along(slidenames), function(slide) {

  current_path <- file_paths[[slide]]
  current_slide <- slidenames[[slide]]

  message("Reading FastReseg object for: ", current_slide)
  files <- dir(current_path)

  # Read FastReseg fov-combined object
  fastreseg_file <- files[grepl("combined_updated_perCellDT_perCellExprs.rds", files)]

  fastreseg_obj <- readRDS(file.path(current_path, fastreseg_file))

  # Create seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = fastreseg_obj[[2]],
                                  meta.data = fastreseg_obj[[1]],
                                  project = "FastReseg")

  seurat_obj@meta.data$slidename <- current_slide
  seurat_obj@meta.data$CenterX_global_px <- seurat_obj@meta.data$x / pixel_size
  seurat_obj@meta.data$CenterY_global_px <- seurat_obj@meta.data$y / pixel_size

  return(seurat_obj)
  })

  #Merging seurat objects
  message("Merging Seurat objects across slides...")
  merged_obj <- merge(
    x = fastreseg_list[[1]],
    y = fastreseg_list[-1],
    add.cell.ids = slidenames)

  return(merged_obj)
}
