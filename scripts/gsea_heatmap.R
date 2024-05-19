gsea_heatmap <- function(gsea, category, cell_type, seurat_obj = data, sample_id_col = "sample_id", sex_col = "sex" ,assay="RNA", project_col ="project_id") {
  # get genes to plot
  gsea@compareClusterResult %>%
    filter(Description == category, Cluster ==cell_type) %>%
    pull(core_enrichment) %>%
    strsplit(split = "/")  %>%
    unlist()-> genes_to_plot
  # get gene names in symbol
  bitr(genes_to_plot, fromType = "UNIPROT", toType="SYMBOL", OrgDb= org.Hs.eg.db) -> genes_to_plot
  # get average expression per sample_id
  sample_id.averages <- AverageExpression(seurat_obj, features = genes_to_plot$SYMBOL,
                                       return.seurat = TRUE, group.by = sample_id_col)
  Idents(sample_id.averages) <- names(Idents(sample_id.averages)) # idk why i do this
  mean_gene_expr <- rowSums(sample_id.averages[[assay]]@scale.data)
  var_genes <- names(mean_gene_expr[mean_gene_expr != 0])
  
  # get M - F sample_id order for heatmap
  sample_id_order <- seurat_obj@meta.data %>%
    distinct(!!sym(sample_id_col), !!sym(sex_col), age)
#################################################################
#    {
 #     if (is.na(project_col)) {
  #      distinct(!!sym(sample_id_col), !!sym(sex_col), age)
   #   } else {
    #    distinct(!!sym(sample_id_col), !!sym(sex_col), !!sym(project_col), age)
     # }
    #}
  # Idents(sample_id.averages) <- factor(sample_id.averages@active.ident, levels = sample_id_order$sample_id_col)
  # Add idents as metadata
  sample_id.averages <- AddMetaData(sample_id.averages, metadata = Idents(sample_id.averages),
                                 col.name = sample_id_col)
  sample_id_order %>%
    right_join(sample_id.averages@meta.data, by = c("sample_id" = "sample_id")) -> sample_id_order
  rownames(sample_id_order) <- sample_id_order$sample_id
  # add Predicted_sex as metadata
  sample_id.averages <- AddMetaData(sample_id.averages, metadata = sample_id_order[,2],
                                 col.name = "sex")
##################################################################################platform-project_id
  sample_id.averages <- AddMetaData(sample_id.averages, metadata = sample_id_order[,3],
                                 col.name = "platform")
  sample_id.averages <- AddMetaData(sample_id.averages, metadata = sample_id_order[,4],
                                 col.name = "age")
  # both sex-project_id
  sample_id.averages <- AddMetaData(sample_id.averages, metadata = paste(sample_id_order[,3], sample_id_order[,2], sep="_"), col.name = "sex_proj")
  levels <- sample_id.averages@meta.data %>% arrange(sex_proj) %>% distinct(sex_proj) %>% pull(sex_proj)
  Idents(sample_id.averages) <- factor(sample_id.averages@meta.data$sex_proj, levels = levels )
  # plot heatmap
  DoHeatmap(sample_id.averages, features = var_genes, size = 2, draw.lines = FALSE, group.by = "sex_proj",
            assay = assay, disp.max = 2, disp.min = -0.8)  + NoLegend()-> gsea_heatmap_plot
  # make boxplot
  sample_id.averages[[assay]]@scale.data %>%
    as.data.frame() %>%
    mutate(gene =rownames(.)) %>%
    filter(gene %in% var_genes) %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "expr") %>%
    left_join(sample_id_order, by = "sample_id") %>%
    mutate(sample_id = factor(sample_id, levels = sample_id_order$sample_id_col)) %>%
    ggplot(aes(x = gene, y = expr, fill = !!sym(sex_col))) +
###############################################################################
    #facet_wrap(vars(project_id)) +
    geom_boxplot() -> boxplot
  return(list(gsea_heatmap_plot, boxplot, levels, sample_id_order, sample_id.averages@meta.data))
}