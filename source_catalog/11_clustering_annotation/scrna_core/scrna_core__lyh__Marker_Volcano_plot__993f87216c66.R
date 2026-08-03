Marker_Volcano<-function (markers = NULL, topn = 5, only.pos=TRUE,
	plot_order= NULL,log2FC = 0.25, labelCol = NULL,
	hlineSize = 1, hlineColor = "grey50", pforce = 5, nforce = 2.5, 
	nudge_x = 0.8, pnudge_y = 0.25, nnudge_y = 0, base_size = 18, 
	facetColor = NA, facetFill = "white", ylab = "Log2-Fold Change", 
	nrow = 1) 
{
	if(!is.null(plot_order)){
		markers$cluster<-factor(markers$cluster,levels = plot_order)
	}
	toppos <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = topn, 
		wt = avg_log2FC)
	if(only.pos == TRUE){
		topgene <- toppos
		plot<-ggplot2::ggplot(markers, ggplot2::aes(x = pct.1 - pct.2, 
			y = avg_log2FC)) + ggplot2::geom_point(color = "grey80") + 
			ggplot2::geom_hline(yintercept = log2FC, 
				lty = "dashed", size = hlineSize, color = hlineColor) + 
			ggrepel::geom_text_repel(data = toppos, ggplot2::aes(x = pct.1 - 
				pct.2, y = avg_log2FC, label = gene, color = cluster), 
				show.legend = F, direction = "y", hjust = 1, nudge_y = pnudge_y, 
				force = pforce, nudge_x = -nudge_x - (toppos$pct.1 - 
					toppos$pct.2)) + ggplot2::geom_point(data = topgene, 
			show.legend = F, ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC, 
				color = cluster)) + ggplot2::scale_color_manual(name = "", 
			values = labelCol) + ggplot2::theme_bw(base_size = base_size) + 
			ggplot2::theme(panel.grid = ggplot2::element_blank(), 
				axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
				strip.background = ggplot2::element_rect(color = facetColor, 
					fill = facetFill)) + ggplot2::xlab(expression(Delta ~ 
			"Percentage Difference")) + ggplot2::ylab(ylab) + ggplot2::facet_wrap(~cluster, 
			nrow = nrow, scales = "fixed")
	}else{
		topnegtive <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = -topn, 
			wt = avg_log2FC)
		topgene <- rbind(toppos, topnegtive)
		plot<-ggplot2::ggplot(markers, ggplot2::aes(x = pct.1 - pct.2, 
			y = avg_log2FC)) + ggplot2::geom_point(color = "grey80") + 
			ggplot2::geom_hline(yintercept = c(-log2FC, log2FC), 
				lty = "dashed", size = hlineSize, color = hlineColor) + 
			ggrepel::geom_text_repel(data = toppos, ggplot2::aes(x = pct.1 - 
				pct.2, y = avg_log2FC, label = gene, color = cluster), 
				show.legend = F, direction = "y", hjust = 1, nudge_y = pnudge_y, 
				force = pforce, nudge_x = -nudge_x - (toppos$pct.1 - 
					toppos$pct.2)) + ggrepel::geom_text_repel(data = topnegtive, 
			ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC, label = gene, 
				color = cluster), show.legend = F, direction = "y", 
			hjust = 0, nudge_y = nnudge_y, force = nforce, nudge_x = nudge_x - 
				(topnegtive$pct.1 - topnegtive$pct.2)) + ggplot2::geom_point(data = topgene, 
			show.legend = F, ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC, 
				color = cluster)) + ggplot2::scale_color_manual(name = "", 
			values = labelCol) + ggplot2::theme_bw(base_size = base_size) + 
			ggplot2::theme(panel.grid = ggplot2::element_blank(), 
				axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
				strip.background = ggplot2::element_rect(color = facetColor, 
					fill = facetFill)) + ggplot2::xlab(expression(Delta ~ 
			"Percentage Difference")) + ggplot2::ylab(ylab) + ggplot2::facet_wrap(~cluster, 
			nrow = nrow, scales = "fixed")
	}
	return(plot)
}
