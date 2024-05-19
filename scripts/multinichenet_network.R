make_ggraph_ligand_target_links = function(lr_target_prior_cor_filtered, prioritized_tbl_oi, colors, seed, ct_order){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggraph")
  requireNamespace("multinichenetr")
  
  lr_target_prior_cor_filtered = lr_target_prior_cor_filtered %>% dplyr::inner_join(prioritized_tbl_oi, by = c("sender", "receiver", "ligand", "receptor", "id", "group"))
  
  source_df_lr = prioritized_tbl_oi %>% dplyr::mutate(celltype_ligand = paste(sender, ligand, sep = "_"), celltype_receptor = paste(receiver, receptor, sep = "_")) %>% dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, ligand, receptor) 
  source_df_lrt = lr_target_prior_cor_filtered %>% dplyr::mutate(celltype_ligand = paste(sender, ligand, sep = "_"), celltype_target = paste(receiver, target, sep = "_"), celltype_receptor = paste(receiver, receptor, sep = "_")) %>% dplyr::select(group, sender, receiver, celltype_ligand, celltype_receptor, celltype_target, ligand, target, receptor, direction_regulation) 
  
  lr_gr_network = dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & !celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand"),
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_receptor  & !celltype_target %in% source_df_lr$celltype_ligand) %>% dplyr::mutate(type_target = "receptor")
  ) %>% dplyr::bind_rows(
    source_df_lrt %>% dplyr::filter(celltype_target %in% source_df_lr$celltype_ligand & celltype_target %in% source_df_lr$celltype_receptor) %>% dplyr::mutate(type_target = "ligand/receptor")
  )
  
  ligand_receptor_network = lr_gr_network %>% dplyr::filter(celltype_receptor %in% celltype_target) %>% dplyr::select(celltype_ligand, celltype_receptor, direction_regulation, group) %>% dplyr::distinct() %>% dplyr::rename(sender = celltype_ligand, receiver = celltype_receptor) %>% dplyr::mutate(type = "Ligand-Receptor", weight = 1)
  
  ligand_target_network = lr_gr_network %>% dplyr::select(celltype_ligand, celltype_target, direction_regulation, group) %>% dplyr::distinct() %>% dplyr::rename(sender = celltype_ligand, receiver = celltype_target) %>% dplyr::mutate(type = "Ligand-Target", weight = 1)
  
  links = ligand_target_network %>% dplyr::bind_rows(ligand_receptor_network) 
  nodes = lr_gr_network %>% dplyr::select(celltype_ligand, sender, ligand) %>% dplyr::rename(celltype = sender, node = celltype_ligand, gene = ligand) %>% dplyr::mutate(type_gene = "ligand") %>% dplyr::bind_rows(
    lr_gr_network %>% dplyr::select(celltype_receptor, receiver, receptor) %>% dplyr::rename(celltype = receiver, node = celltype_receptor, gene = receptor) %>% dplyr::mutate(type_gene = "receptor")
  ) %>% dplyr::bind_rows(
    lr_gr_network %>% dplyr::select(celltype_target, receiver, target, type_target) %>% dplyr::rename(celltype = receiver, node = celltype_target, gene = target, type_gene = type_target)
  ) %>% dplyr::distinct() %>% dplyr::filter(node %in% c(links$sender, links$receiver))
  
  double_nodes =  nodes %>% dplyr::group_by(node) %>% dplyr::count() %>% dplyr::filter(n > 1) %>% pull(node)
  nodes = dplyr::bind_rows(
    nodes %>% dplyr::filter(node %in% double_nodes) %>% dplyr::mutate(type_gene = "ligand/receptor") ,
    nodes %>% dplyr::filter(!node %in% double_nodes)
  ) %>% dplyr:: distinct()
  
  nodes = nodes %>% data.frame() %>% magrittr::set_rownames(nodes$node)
  
  colors_regulation = NULL
  colors_regulation["up"] = "indianred1"
  colors_regulation["down"] = "steelblue2"
  
  # create the network object
  network = igraph::graph_from_data_frame(d=links %>% dplyr::filter(type == "Ligand-Target"), vertices = nodes, directed=T)
  graph = tidygraph::as_tbl_graph(network)
  celltype_values <- V(graph)$celltype
  celltype_values_reordered <- factor(celltype_values, levels = ct_order)
  V(graph)$celltype <- celltype_values_reordered
  
  set.seed(seed)
  plot =  ggraph::ggraph(graph, layout = 'dh') + 
    ggraph::geom_edge_fan(aes(color = direction_regulation), edge_width = 1, arrow = arrow(length = unit(3, 'mm')), end_cap = ggraph::circle(5.5, 'mm'), start_cap = ggraph::circle(3, 'mm')) + 
    # ggraph::geom_edge_loop(aes(color = direction_regulation), edge_width = 1, alpha = 0.70)  + 
    ggraph::geom_node_label(aes(label = gene, fill = celltype), fontface = "bold", size = 3.5, nudge_x = 0, nudge_y = 0, color = "whitesmoke") +
    ggraph::theme_graph(foreground = 'black', fg_text_colour = 'white', base_family = 'Helvetica') + facet_grid(. ~group)  + ggraph::scale_edge_color_manual(values = colors_regulation) + scale_fill_manual(values = colors)
  
  return(list(plot = plot, graph = graph, source_df_lr = source_df_lr, source_df_lt = links %>% dplyr::filter(type == "Ligand-Target"), nodes_df = nodes))
}