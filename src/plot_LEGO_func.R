#!/home/sherry/bin/Rscript

## read in LEGO results
plot_LEGO <- function(gs, pv, int_g, choose = "full", inter = "out", ii = "", net_g, gs_list, net_data, gs_mod_data, min_size, local=0) {
	## color
	mypalette <- c(brewer.pal(8, "PuRd")[5], brewer.pal(8, "GnBu")[5], brewer.pal(8, "Oranges")[5], brewer.pal(8, "Blues")[7], brewer.pal(6, "Greys")[4], brewer.pal(8, "Greens")[8],brewer.pal(8, "Purples")[7])
	cc1<-colorRampPalette(mypalette)(length(mypalette)); mypalette <- rgb(t(col2rgb(cc1))/256,alpha=0.8);
	tt <- intersect(gs_list[[gs]], net_g); ## number of gs genes in network
	ori_int_g <- int_g;
	int_g <- intersect(net_g, int_g);
	ov <- intersect(tt, int_g); ov <- unique(c(ov,intersect(gs_list[[gs]],ori_int_g)));
	use_g <- unique(c(tt, int_g,intersect(gs_list[[gs]],ori_int_g)));
	if (length(tt) < min_size) 
		return()
	if (inter == "out") {
		net1 <- net_data[which(net_data[, 1] %in% use_g | net_data[, 2] %in% use_g), ];
		use_g <- unique(c(net1[, 1], net1[, 2]));
	}
	if (inter == "in") 
		net1 <- net_data[which(net_data[, 1] %in% use_g & net_data[, 2] %in% use_g), ];
	if(local==1){
		net1 <- net1[which(net1[,1] %in% gs_list[[gs]] | net1[,2] %in% gs_list[[gs]]), ];
		use_gg <- intersect(use_g,c(net1[,1],net1[,2],gs_list[[gs]]));
		g2 <- graph.data.frame(net1, directed = F, vertices = data.frame(name = use_gg));
	}else{
		g2 <- graph.data.frame(net1, directed = F, vertices = data.frame(name = use_g));
	}
	col_e <- rep("black", length = length(E(g2)))
	col_e[which(net1[, 1] %in% tt & net1[, 2] %in% tt)] <- "blue" ## gene set inner link
	col_e[which(net1[, 1] %in% tt & net1[, 2] %in% int_g | net1[, 1] %in% int_g & net1[, 2] %in% tt)] <- "dark red"
	## int-gs 
	col_n <- rep("#D3D3D3", length = length(V(g2)))
	col_n[which(V(g2)$name %in% int_g)] <- mypalette[3]
	col_n[which(V(g2)$name %in% tt)] <- mypalette[2]
	col_n[which(V(g2)$name %in% ov)] <- mypalette[1]
	ind_n <- ifelse(col_n == "#D3D3D3", 0, 1)
	if (choose == "full") {
		### get use_g (int+gs) + define color           
		E(g2)$width <- ifelse(col_e == "black", 0.5, 1.5)
		E(g2)$color <- col_e
		V(g2)$color <- col_n
		u1 <- paste(V(g2)$name, gs, sep = "~")
		u2 <- which(u1 %in% rownames(gs_mod_data))
		V(g2)$size <- rep(3, length(V(g2)$name))
		V(g2)$size[u2] <- gs_mod_data[u1[u2], 3] + 3
		V(g2)$label <- ifelse(V(g2)$name %in% use_g, V(g2)$name, "");
		return(g2)
	}
	if (choose == "sub") {
		c1 <- clusters(g2)
		c2 <- c1$membership * ind_n
		c3 <- names(which(table(c2[which(c2 > 0)]) >= 2))
		c4 <- which(c1$membership %in% c3)
		### get subgraph
		g3 <- subgraph(g2, c4)
		net2 <- get.data.frame(g3)
		col_e <- rep("black", length = length(E(g3)))
		col_e[which(net2[, 1] %in% tt & net2[, 2] %in% tt)] <- "blue" ## gene set inner link
		col_e[which(net2[, 1] %in% tt & net2[, 2] %in% int_g | net2[, 1] %in% int_g & net2[, 2] %in% tt)] <- "dark red"
		## int-gs 
		col_n <- rep("#D3D3D3", length = length(V(g3)))
		col_n[which(V(g3)$name %in% int_g)] <- mypalette[3]
		col_n[which(V(g3)$name %in% tt)] <- mypalette[2]
		col_n[which(V(g3)$name %in% ov)] <- mypalette[1]
		ind_n <- ifelse(col_n == "#D3D3D3", 0, 1)
		E(g3)$width <- ifelse(col_e == "black", 0.5, 1.5)
		E(g3)$color <- col_e
		V(g3)$color <- col_n
		u1 <- paste(V(g3)$name, gs, sep = "~")
		u2 <- which(u1 %in% rownames(gs_mod_data))
		V(g3)$size <- rep(3, length(V(g3)$name))
		V(g3)$size[u2] <- gs_mod_data[u1[u2], 3] + 3
		V(g3)$label <- ifelse(V(g3)$name %in% use_g, V(g3)$name, "")
		return(g3)
	}
}
##
