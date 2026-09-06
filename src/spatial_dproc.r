################################################################################
################################################################################
#' These functions for downstream analysis of single cell spatial transcriptomics
#' data such as differential gene expression analysis, GSEA, etc.  
################################################################################

#'Differential gene expression analysis using linear regression or a linear mixed effect model 
#' 
#' @param obj Seurat object
#' @param assay Assay name to use for DGE analysis
#' @param condition_col Metadata column used for the independent variable
#' @param ref_condition Group of cells used as reference
#' @param var_condition Group of cells used to compare against the reference cells
#' @param donor_col Metadata column used as random variable for a lmm model (If NULL, a simple linear 
#' regression model will be used instead)
#' @param output_dir Output directory where the results table will be saved
#'
#' @return A results table with DGE statistical data. 
#' 
lmm_res <- function(obj, assay, condition_col, ref_condition, var_condition, donor_col = NULL) {

    DefaultAssay(obj) <- assay

    genes <- rownames(obj[[assay]])
    results <- data.frame(gene = genes, estimate = NA, p_value = NA)

    for (i in seq_along(genes)) {
    gene <- genes[i]
    message("Evaluating gene: ", gene, " ", i, "/", length(genes))
    expr <- FetchData(obj, vars = gene, assay = assay, layer = "counts")[[1]]
    
    if (!is.null(donor_col)) {
    df_expr <- data.frame(
        expression = expr,
        condition = obj@meta.data[[condition_col]],
        donor = obj@meta.data[[donor_col]])
    } else {
    df_expr <- data.frame(
        expression = expr,
        condition = obj@meta.data[[condition_col]])
    }
    
    df_expr$condition <- relevel(factor(df_expr$condition), ref = ref_condition)

    if (!is.null(donor_col)) {

        message("DGE with a linear-mixed effect model")

        model <- tryCatch(
        lme4::lmer(expression ~ condition + (1 | donor), data = df_expr), #Linear mixed effect model
        error = function(e) NULL)
    } else {

        message("DGE with a simple linear regression model")

        model <- tryCatch(
        lm(expression ~ condition, data = df_expr), #Standard linear regression
        error = function(e) NULL)
    }
    

    if (!is.null(model)) {
        s <- summary(model)$coefficients
        results$estimate[i] <- s[paste0("condition", var_condition), "Estimate"]  # adjust to match your factor level

        if (!is.null(donor_col)) {
            results$p_value[i]  <- s[paste0("condition", var_condition), "t value"]
        } else {
           results$p_value[i]  <- s[paste0("condition", var_condition), "Pr(>|t|)"]
        }
    }
    else {message(gene, " is Null")}
    }

    results$padj <- p.adjust(results$p_value, method = "BH")

    #calculate -log10 padj and rewrite inf results (padj = 0) as 0
    results$p_value[results$p_value == 0] <- "NA"
    results$padj[results$padj == 0] <- "NA"


    #assign significant changes by log2fold change and padj for RNA < 0.01 and padj for TE < 0.05
    assign_sig_changes <- function(results) {
    lmm_df <- results
    lmm_df$Sign_class <- "Not_Significant"
    lmm_df$Sign_class[lmm_df$estimate > 0.5 & lmm_df$p_value < 0.05] <- "upregulated"
    lmm_df$Sign_class[lmm_df$estimate < -0.5 & lmm_df$p_value < 0.05] <- "downregulated"

    lmm_df$fdr <- NA
    lmm_df$fdr[lmm_df$estimate > 0.5 & lmm_df$padj < 0.05] <- "fdr"
    lmm_df$fdr[lmm_df$estimate < -0.5 & lmm_df$padj < 0.05] <- "fdr"
    
    #convert Tx_class to factor for color purposes
    lmm_df$Sign_class <- factor(lmm_df$Sign_class)
    lmm_df$fdr <- factor(lmm_df$fdr)
    
    return(lmm_df)
    }


    results <- assign_sig_changes(results)

    return(results)
}


#'Differential gene expression analysis using Nanostring smiDE package
#' 
#' @param obj Seurat object
#' @param assay Assay name to use for DGE analysis
#' @param condition_col Metadata column used for the independent variable
#' @param totalcounts Metadata column corresponding to total transcript count per cell (i.e. nCount_RNA)
#' @param cell_type_col Metadata column corresponding
#' @param random_var Metadata column used as random variable for the lmm model 
#' @param slide_col Metadata column corresponding to individual slides. 
#' @param cellid Metadata column corresponding to cell IDs
#' @param x_coord X coordinates for each cell
#' @param y_coord Y coordinates for each cell
#' @param radius Radius for neighboring cells in x_coord and y_coord units 
#'
#' @return A results table with DGE statistical data. 
#' 
smide_res <- function(obj, assay, condition_col, totalcounts, cell_type_col, random_var, slide_col = NULL, cellid, x_coord, y_coord, radius) {
  
  DefaultAssay(obj) <- assay
    
  totalcount_scalefactors <- mean(obj@meta.data[[totalcounts]]) / obj@meta.data[[totalcounts]]
  names(totalcount_scalefactors) <- obj@meta.data[[cellid]]
 
  obj <- Seurat::SetAssayData(obj
                              ,"data"
                              ,obj[[assay]]$counts %*% Matrix::Diagonal(x=totalcount_scalefactors, names=colnames(obj))
  )
  
  pre_de_obj <- 
    smiDE::pre_de(metadata = obj@meta.data
           ,cell_type_metadata_colname = cell_type_col
           ,split_neighbors_by_colname = slide_col
           ,cellid_colname = cellid
           ,mm_radius = radius #Units are in microns
           ,sdimx_colname = x_coord
           ,sdimy_colname = y_coord
           ,verbose=TRUE
    )
  
  
  overlap_metrics <- smiDE::overlap_ratio_metric(assay_matrix = obj[[assay]]$data #In the Nanostring website, they use the data slot instead.
                                          ,metadata = obj@meta.data
                                          ,cellid_col = cellid
                                          ,cluster_col = cell_type_col
                                          ,sdimx_col = x_coord
                                          ,sdimy_col = y_coord
                                          ,radius = radius #Units in microns
  )
  
  
  cell_types_list <- unique(obj@meta.data[[cell_type_col]])
  metainfo <- data.table::data.table(obj@meta.data)
  
  
  dge_list <- lapply(seq_along(cell_types_list), function(i) {
    
    message("Analyzing DEG for: ", cell_types_list[i])
    
    genes_to_analyze <- overlap_metrics[overlap_metrics[[cell_type_col]] ==
                                          cell_types_list[i]][ratio < 1][["target"]]
    
    cells_interest <- metainfo[metainfo[[cell_type_col]] == cell_types_list[i]][[cellid]]
    
    
    de_obj <- tryCatch(
      smiDE::smi_de(assay_matrix = obj[[assay]]$counts
                    ,cellid_colname = cellid
                    ,metadata = metainfo[metainfo[[cellid]] %in% cells_interest]
                    ,formula = bquote(~ RankNorm(otherct_expr) + .(as.name(independ_variable)) +
                                        (1 | .(as.name(random_var))) + offset(log(.(as.name(totalcounts)))))
                    ,pre_de_obj = pre_de_obj
                    ,neighbor_expr_cell_type_metadata_colname = cell_type_col
                    ,neighbor_expr_overlap_weight_colname = NULL
                    ,neighbor_expr_overlap_agg ="sum"
                    ,neighbor_expr_totalcount_normalize = TRUE
                    ,neighbor_expr_totalcount_scalefactor = totalcount_scalefactors
                    ,family="nbinom2"
                    ,targets = genes_to_analyze
      ), error = function(e) NULL)
    
    if (!is.null(de_obj)) {
      de_results <- results(de_obj, comparisons = "pairwise", variable = independ_variable)[[1]]
      de_results$celltype_col <- cell_types_list[i]
      
    } else {
      de_results <- NULL
    }
    
    return(de_results)
  })
  
  combined_results <- do.call(rbind, c(dge_list, list(fill = TRUE)))
  
  return(combined_results)
}

