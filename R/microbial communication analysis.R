"net_hierarchy"<- function(abun_data=abun,
                           tax_data=taxon,
                           char.start=4,
                           char.stop=9,
                           select_taxa="Phylum",
                           diag=FALSE,
                           pvaule=0.05,
                           rrr=0.2,
                           filter=TRUE) {

  temp<-taxa_select_abun(abun_data=abun_data,
                         tax_data=tax_data,
                         char.start=char.start,
                         char.stop=char.stop,
                         select_taxa=select_taxa)

  node_list<-rel_3abun(temp)

  mat<-net_mat(temp,
               diag=diag,
               pvaule=pvaule,
               rrr=rrr,
               filter=filter)$mat

  return(list(mat=mat,node_list=node_list))
}

"taxa_select_abun" <- function(abun_data=abun,
                               tax_data=taxon,
                               char.start=char.start,
                               char.stop=char.stop,
                               select_taxa="Phylum"
) {
  split_otu <- filter_merge_sum(otudata=abun_data,
                                taxon_table=tax_data)

  group_all<-split_otu$data_venn%>%rownames()
  group_num<-length(group_all)

  sample_sumnum<-abun_data%>%colnames()%>%length()

  abun_data<-rownames_to_column(abun_data,var = "name")
  tax_data<-rownames_to_column(tax_data,var = "name")
  abun_new<-merge(abun_data,tax_data,by="name")
  abun_new<-column_to_rownames(abun_new,var = "name")

  select_num<-match(select_taxa,colnames(abun_new))
  abun_new<-subset(abun_new,select=c(1:sample_sumnum,select_num))
  colnames(abun_new)[length(abun_new)]<-"taxx"

  temp<-abun_new%>%group_by(taxx)%>%summarise_all(sum)
  temp[which(temp$taxx==""),]$taxx<-paste("p__Others")
  temp$taxx<-substr(temp$taxx,start = char.start,stop = char.stop)
  temp<-temp%>%group_by(taxx)%>%summarise_all(sum)
  temp<-column_to_rownames(temp,var = "taxx")
  return(temp)
}

"net_mat" <- function(otu_rare=temp,
                      diag=diag,pvaule=pvaule,rrr=rrr,
                      filter=filter) {
  colname<-colnames(otu_rare)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  groupname<-substr(colnames(otu_rare),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otu_rare))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()


  if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <- lapply(apply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),2,FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <- lapply(lapply(sapply(trt_id,function(x){grep(x,colnames(otu_rare))}),FUN = function(x){otu_rare[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  otu<-otu_filter(split_otu,filter_num = 1,filter = FALSE,self = FALSE)

  mat <- function(g) {
    gmat<-lapply(g, function(x){
      gg<-as.matrix(x)
      return(gg)
    })
    return(gmat)
  }

  otus<-mat(otu)

  nor <- lapply(otus,function(x){
    new<-scale(t(x))##

    occor.r = cor(new,method = "pearson")
    mtadj<-multtest::mt.rawp2adjp(unlist(occor.r),proc='BH')
    adpcor<-mtadj$adjp[order(mtadj$index),2]
    occor.p<-matrix(adpcor,dim(t(x))[2])

    occor<-occor.r

    if (diag) diag(occor) <- 0

    occor[occor.p>pvaule] = 0
    #occor[occor.p>pvaule|abs(occor)<rrr] = 0
    occor[abs(occor)<rrr] = 0
    occor[is.na(occor)]=0
    return(occor)
  })

  trvalue <- lapply(nor,function(x){

    res <- rm.get.threshold(x,discard.zeros = TRUE,plot.spacing=FALSE,
                            save.fit = FALSE,
                            unfold.method = "gaussian")

    xas<-res$tested.thresholds
    yas<-res$p.ks
    plot<-cbind(xas,yas)%>%data.frame()
    threshold<-plot$xas[which(plot$yas == max(plot$yas))]
    thre <- threshold
    return(list(bor=x,plot=plot,result=thre))
  })


  st.trvalue<-lapply(trvalue, function(x){
    x$result<-0#
    return(x)
  })


  net1 <- lapply(st.trvalue,function(x){
    cleaned.matrix <- rm.denoise.mat(x$bor, threshold = x$result) #denoise adjacency matrix by threshold
    cleaned.matrix <- rm.discard.zeros(cleaned.matrix) #delete all-zero rows
    g = igraph::graph_from_adjacency_matrix(cleaned.matrix ,mode="undirected",weighted=TRUE,diag=FALSE)
    return(list(g=g,cleaned.matrix=cleaned.matrix))
  })

  g<-lapply(net1, function(x){
    gg<-x$g
    return(gg)
  })

  mat<-lapply(net1, function(x){
    gg<-x$cleaned.matrix
    return(gg)
  })

  return(list(g=g,mat=mat))
}

"plotMicrochatInterAllTaxaLine" <- function(submchat,
                                              line.size=5,
                                              lineend.shape="round",
                                              color_group=colorCustom(5,pal = "gygn"),
                                              export_path="microbial chat") {

  export_path<-paste(export_path,"/data_microbiome/microbial communication analysis/inter-taxa-interaction",sep = "")
  dir.create(export_path, recursive = TRUE)

  dir.create(export_path, recursive = TRUE)
  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  abun<-submchat$otu_table
  taxon<-submchat$taxon_table

  taxa.use<-colnames(taxon)
  taxa.num<-taxa.use%>%length()
  for (k in 2:taxa.num) {
    name.use<-taxa.use[k]

    data_hier1<-net_hierarchy(abun_data=abun,
                              tax_data=taxon,
                              char.start=1,
                              char.stop=10,
                              select_taxa=name.use,
                              diag=TRUE,
                              pvaule=0.05,
                              rrr=NULL)
    mat1<-data_hier1$mat
    node_list1<-data_hier1$node_list
    newmat<-sapply(mat1, function(x) {
      y<-x
      interweight<-sum(abs(x))
      y[y!=0]<-1
      internum<-sum(y)
      return(list(internum=internum,interweight=interweight))
    })
    newmat1<-newmat%>%t()%>%data.frame()%>%rownames_to_column(var = "group")
    newmat1$internum<-newmat1$internum%>%as.numeric()
    newmat1$interweight<-newmat1$interweight%>%as.numeric()
    newmat2<-newmat1 %>%
      pivot_longer(cols = !group,
                   names_to = "param",
                   values_to = "value")
    newmat2$group<-factor(newmat2$group,levels = unique(newmat2$group))
    newmat2$value<-round(newmat2$value,3)
    newmat2$value.y<-round(newmat2$value*1.1,3)

    trt_id<-unique(newmat2$param)
    plot_data<-apply(
      sapply(trt_id,function(x){
        grep(x,newmat2$param)
      }),2,FUN = function(x){
        newmat2[x,]
      })


    tk=1
    intername<-names(plot_data)[tk]
    plot_data1<-plot_data[[tk]]

    p1<-ggplot(data= plot_data1,aes(x=group, y=value,group=param)) +
      scale_x_discrete(limits=levels(plot_data1$group))+
      theme(panel.background = element_rect(fill = "white"))


    for (i in 1:(nrow(plot_data1))) {
      p1 <- p1 +  annotate('rect', xmin=i-0.5,xmax=i+0.5,ymin=-Inf,ymax=Inf,
                           fill=ifelse(i %% 2 == 0, 'white', 'grey95'))
    }


    col_fun<-colorRampPalette(color_group)( length(unique(plot_data1$group)) )
    names(color_group)<-unique(plot_data1$group)


    p1<-p1+
      geom_path(color=col_fun,size = line.size, lineend = lineend.shape)+
      ggnewscale::new_scale_color()+
      geom_point(aes(color=group)) +
      scale_color_manual(values = color_group)+
      ggnewscale::new_scale_color()+
      ggalt::geom_xspline(data= plot_data1,
                          size=1,color="gray70",linetype="dashed",
                          aes(x=group, y=value),spline_shape = 1)+
      geom_text(aes(x = group,y = value.y,label = value),size = 3,color = "black",family = "serif") +
      ylab("Inferred interactions number") +
      theme_test()+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        legend.position = "none",aspect.ratio = 1)

    tk=2
    intername<-names(plot_data)[tk]
    plot_data1<-plot_data[[tk]]

    p2<-ggplot(data= plot_data1,aes(x=group, y=value,group=param)) +
      scale_x_discrete(limits=levels(plot_data1$group))+
      theme(panel.background = element_rect(fill = "white"))


    for (i in 1:(nrow(plot_data1))) {
      p2 <- p2 +  annotate('rect', xmin=i-0.5,xmax=i+0.5,ymin=-Inf,ymax=Inf,
                           fill=ifelse(i %% 2 == 0, 'white', 'grey95'))
    }


    col_fun<-colorRampPalette(color_group)( length(unique(plot_data1$group)) )
    names(color_group)<-unique(plot_data1$group)


    p2<-p2+
      geom_path(color=col_fun,size = line.size, lineend = lineend.shape)+
      ggnewscale::new_scale_color()+
      geom_point(aes(color=group)) +
      scale_color_manual(values = color_group)+
      ggnewscale::new_scale_color()+
      ggalt::geom_xspline(data= plot_data1,
                          size=1,color="gray70",linetype="dashed",
                          aes(x=group, y=value),spline_shape = 1)+
      geom_text(aes(x = group,y = value.y,label = value),size = 3,color = "black",family = "serif") +
      ylab("Inferred interactions strength") +
      theme_test()+
      theme(#panel.border = element_blank(),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(hjust = 0,vjust =0.5,colour='black',size=10,family = "serif"),
        axis.text.x=element_text(colour = "black",size = 10,
                                 angle = 0,hjust = 0.5,vjust =0.5,family = "serif"),
        legend.position = "none",aspect.ratio = 1)

    p2<-p2+coord_flip()
    p1<-p1+coord_flip()+scale_y_reverse()
    taxap<-cowplot::plot_grid(p1,p2,ncol = 2)

    taxap<- p1 + p2 + patchwork::plot_layout(widths = c(1,1),byrow=TRUE)
    ggsave(paste(export_path,"/",name.use,"_internection_num&weight.pdf",sep = ""),
           taxap)
  }
  return(taxap)
}



"plot_net_single_tax" <- function(mat,select_group="ct",
                                  color_taxa=colorCustom(50,pal = "gygn")) {
  net<-mat[[match(select_group,names(mat))]]
  net<-abs(net)

  max.size<-c()
  for (tt in 1:length(mat)) {
    maxsize<-nrow(mat[[tt]])
    {
      max.size<-c(max.size,maxsize)
    }
  }
  rownum<-(max(max.size)/5)%>%as.integer()+1
  par(mfrow = c(5,rownum),
      mar=c(1,0,1,0),family="serif",

      xpd=TRUE)
  for (i in 1:nrow(net)) {
    mat2 <- matrix(0,
                   nrow = nrow(net),
                   ncol = ncol(net),
                   dimnames = dimnames(net))
    mat2[i, ] <- net[i, ]
    netVisual_2circle(mat2,
                      color.use =color_taxa ,
                      weight.scale = T,
                      edge.weight.max = max(net),
                      title.name = rownames(net)[i])
  }

}


"netVisual_2circle"<-function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
                               targets.use = NULL, remove.isolate = FALSE, top = 1, weight.scale = FALSE,
                               vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL,
                               vertex.label.cex = 1, vertex.label.color = "black", edge.weight.max = NULL,
                               edge.width.max = 8, alpha.edge = 0.6, label.edge = FALSE,
                               edge.label.color = "black", edge.label.cex = 0.8, edge.curved = 0.2,
                               shape = "circle", layout = in_circle(), margin = 0.2, vertex.size = NULL,
                               arrow.width = 1, arrow.size = 0.2)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g <- igraph::graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = CellChat::scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0,
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * edge.width.max
    igraph::E(g)$arrow.width <- 0
    igraph::E(g)$width<-igraph::E(g)$width
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    igraph::E(g)$arrow.width <- 0
    igraph::E(g)$width<-igraph::E(g)$width
  }

  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,
                                                                             1]], alpha.edge)
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  "radian.rescale" <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape,
       layout = coords_scale, margin = margin, vertex.label.dist = 4,
       vertex.label.degree = label.locs, vertex.label.family = "serif",
       edge.label.family = "serif")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}



"plot_net_hierarchy" <- function(mat,node_list,abs=TRUE,weight.scale=TRUE,
                                 arrowwidth = 1,
                                 arrow.size = 0.2,
                                 node.size.constant=FALSE,
                                 edge.label.color = "black",
                                 edge.label.cex = 0.5,
                                 select_group="st",
                                 colors=randomcolor(50),
                                 vertex.receiver = seq(1,5),
                                 vertex.weight = NULL,
                                 vertex.weight.max = NULL,
                                 vertex.size.max = NULL,
                                 edge.weight.max = 1,
                                 edge.curved=0,
                                 alpha.edge=0.6,margin=0.4,
                                 edge.width.max = 3,
                                 title.name = NULL,
                                 vertex.label.cex = 1
) {

  prob.sum<-mat[[match(select_group,names(mat))]]
  if (!node.size.constant) vertex.weight<-node_list[[match(select_group,names(node_list))]]$sum
  if (node.size.constant) vertex.weight<-vertex.weight
  vertex.receiver1<-vertex.receiver

  if(abs) prob.sum<-abs(prob.sum)
  par(mfrow = c(1,2), xpd=NA,family="serif")
  netVisual_1hierarchy(prob.sum,
                       vertex.receiver = vertex.receiver1,
                       top = 1,  margin=margin,
                       arrow.width = arrowwidth,
                       arrow.size = arrow.size,
                       edge.label.color = edge.label.color,
                       edge.label.cex = edge.label.cex,
                       color.use = colors,
                       vertex.weight = vertex.weight,
                       vertex.weight.max = vertex.weight.max,
                       space.v = 1.5,
                       space.h = 1.6,
                       alpha.edge=alpha.edge,
                       vertex.size.max = vertex.size.max,
                       weight.scale = weight.scale,
                       edge.weight.max = edge.weight.max,
                       edge.curved=edge.curved,
                       edge.width.max = edge.width.max,
                       title.name = title.name,
                       vertex.label.cex = vertex.label.cex)

  netVisual_2hierarchy(prob.sum,
                       vertex.receiver = setdiff(1:nrow(prob.sum),
                                                 vertex.receiver),

                       top = 1,  alpha.edge=alpha.edge,
                       arrow.width = arrowwidth,
                       arrow.size = arrow.size,
                       edge.label.color = edge.label.color,
                       edge.label.cex = edge.label.cex,
                       color.use = colors,
                       vertex.weight = vertex.weight,
                       vertex.weight.max = vertex.weight.max,
                       space.v = 1.5,
                       space.h = 1.6,
                       margin=margin,

                       vertex.size.max = vertex.size.max,
                       weight.scale = weight.scale,
                       edge.weight.max = edge.weight.max,
                       edge.curved=edge.curved,
                       edge.width.max = edge.width.max,
                       title.name = title.name,
                       vertex.label.cex = vertex.label.cex)

  rowname<-node_list[[match(select_group,names(node_list))]]%>%rownames()

  legend("bottom",
         inset = 0,
         rowname,bg=NA,bty="n",
         pch = 20,title.col="black",
         xjust = 0,yjust = 1,
         col = colors[1:length(rowname)],box.lwd = 0,
         text.col=colors[1:length(rowname)],
         ncol = 3,title = "Relative abundance",text.font = 2)

}



"netVisual_1hierarchy"<-function (net, vertex.receiver=vertex.receiver,
                                  color.use = NULL,
                                  title.name = NULL,
                                  sources.use = NULL,
                                  targets.use = NULL,
                                  remove.isolate = FALSE,
                                  top = 1,
                                  weight.scale = FALSE,
                                  vertex.weight = 20,
                                  vertex.weight.max = NULL,
                                  vertex.size.max = NULL,
                                  edge.weight.max = NULL,
                                  edge.width.max = 8,
                                  alpha.edge = alpha.edge,
                                  label.dist = 2.8,
                                  space.v = 1.5,
                                  space.h = 1.6,
                                  shape = NULL,
                                  label.edge = FALSE,
                                  edge.curved = 0,
                                  margin = margin,
                                  vertex.label.cex = 0.6,
                                  vertex.label.color = "black",
                                  arrow.width = 1,
                                  arrow.size = 0.2,
                                  edge.label.color = "black",
                                  edge.label.cex = 0.5,
                                  vertex.size = NULL)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  cells.level <- rownames(net)
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  if (is.null(color.use)) {
    color.use <- CellChat::scPalette(nrow(net))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    6
  m <- length(vertex.receiver)
  net2 <- net
  reorder.row <- c(vertex.receiver, setdiff(1:nrow(net), vertex.receiver))
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0,
                                                       n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1,
                                                                               vertex.receiver)], rep("", m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,
                                                                vertex.receiver)], rep("#FFFFFF", length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,
                                                                      vertex.receiver)], color.use[vertex.receiver])
  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[vertex.receiver], vertex.weight[setdiff(1:m1,
                                                                            vertex.receiver)], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- c(rep("circle", m), rep("circle", m1 - m),
               rep("circle", m))
  }
  g <- igraph::graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)
  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m, 1] <- 0
  coords[(m + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
  coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                            m - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1 -
                                                                     1))
  coords_scale <- coords
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    E(g)$label <- E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * edge.width.max
  }
  else {
    E(g)$width <- 0.3 + edge.width.max * E(g)$weight
  }

  if (is.null(arrow.width)) {
    E(g)$arrow.width <- E(g)$width/3
  } else {
    E(g)$arrow.width <- arrow.width
  }

  if (is.null(arrow.size)) {
    E(g)$arrow.size <- E(g)$width/3
  } else {
    E(g)$arrow.size <- arrow.size
  }



  E(g)$label.color <- edge.label.color
  E(g)$label.cex <- edge.label.cex
  E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,
                                                          1]], alpha.edge)
  label.dist <- c(rep(space.h * label.dist, m),
                  rep(space.h * label.dist, m1 - m),
                  rep(0, nrow(net3) - m1))
  label.locs <- c(rep(-pi, m), rep(0, m1 - m),
                  rep(-pi, nrow(net3) - m1))
  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)
  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle, parameters = list(vertex.frame.color = 1,
                                                              vertex.frame.width = 1))
  plot(g, edge.curved = edge.curved, layout = coords_scale,
       margin = margin, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) -
                                                m1)),
       vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")
  text(text.pos, c("Source", "Target", "Source"), cex = 0.8,
       col = c("#c51b7d", "#c51b7d", "#2f6661"),family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, -0.04,
                  space.v - space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, 0.04,
                  space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = "#c51b7d", arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3],
                arrow.pos2[4], col = "#2f6661", arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1], title.pos[2], paste0(title.name,
                                            " network"), cex = 1)
  }
  gg <- recordPlot()
  return(gg)

}

"mycircle" <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }



  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}


"netVisual_2hierarchy"<-function (net, vertex.receiver=vertex.receiver, color.use = NULL, title.name = NULL,
                                  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE,
                                  top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL,
                                  vertex.size.max = NULL, edge.weight.max = NULL, edge.width.max = 8,
                                  alpha.edge = alpha.edge, label.dist = 2.8, space.v = 1.5, space.h = 1.6,
                                  shape = NULL, label.edge = FALSE, edge.curved = 0, margin = margin,
                                  vertex.label.cex = 0.6, vertex.label.color = "black", arrow.width = 1,
                                  arrow.size = 0.2, edge.label.color = "black", edge.label.cex = 0.5,
                                  vertex.size = NULL)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- levels(object@idents)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- levels(object@idents)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- levels(object@idents)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  if (is.null(color.use)) {
    color.use <- CellChat::scPalette(nrow(net))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    6
  m <- length(vertex.receiver)
  m0 <- nrow(net) - length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net), vertex.receiver),
                   vertex.receiver)
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0,
                                                       n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1, vertex.receiver)],
                       row.names(net)[vertex.receiver], rep("", m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(1:m1, vertex.receiver)],
                  color.use[vertex.receiver], rep("#FFFFFF", length(vertex.receiver)))
  color.use3.frame <- c(color.use[setdiff(1:m1, vertex.receiver)],
                        color.use[vertex.receiver], color.use[vertex.receiver])
  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[setdiff(1:m1, vertex.receiver)],
                      vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle", nrow(net3))
  }
  g <- igraph::graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)
  edge.start <- ends(g, es = igraph::E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0, 1] <- 0
  coords[(m0 + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m0, 2] <- seq(space.v, 0, by = -space.v/(m0 - 1))
  coords[(m0 + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                             m0 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1 -
                                                                     1))
  coords_scale <- coords
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max *
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (is.null(arrow.width)) {
    E(g)$arrow.width <- E(g)$width/4
  } else {
    E(g)$arrow.width <- arrow.width
  }

  if (is.null(arrow.size)) {
    E(g)$arrow.size <- E(g)$width/4
  } else {
    E(g)$arrow.size <- arrow.size
  }

  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,
                                                                  1]], alpha.edge)
  label.dist <- c(rep(space.h * label.dist, m), rep(space.h *
                                                      label.dist, m1 - m), rep(0, nrow(net3) - m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1 - m0), rep(-pi,
                                                     nrow(net3) - m1))
  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)
  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle, parameters = list(vertex.frame.color = 1,
                                                              vertex.frame.width = 1))
  plot(g, edge.curved = edge.curved, layout = coords_scale,
       margin = margin, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) -
                                                m1)), vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")
  text(text.pos, c("Source", "Target", "Source"), cex = 0.8,
       col = c("#c51b7d", "#2f6661", "#2f6661"),family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, -0.04,
                  space.v - space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, 0.04,
                  space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = "#c51b7d", arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3],
                arrow.pos2[4], col = "#2f6661", arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1], title.pos[2], paste0(title.name,
                                            " network"), cex = 1)
  }
  gg <- recordPlot()
  return(gg)
}



"calcMicrochatInter" <- function(submchat,
                                 char.start=4,
                                 char.stop=9,
                                 select_taxa="Phylum",
                                 diag=FALSE,
                                 pvaule=0.05,
                                 rrr=0.2,
                                 filter=TRUE) {


  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  #dir.create(export_path, recursive = TRUE)
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  abun<-submchat$otu_table
  taxon<-submchat$taxon_table

  select_taxa1<-colnames(taxon)[match(tolower(select_taxa),tolower(colnames(taxon)))]

  data_hier1<-net_hierarchy(abun_data=abun,
                            tax_data=taxon,
                            char.start=char.start,
                            char.stop=char.stop,
                            select_taxa=select_taxa1,
                            diag=diag,
                            pvaule=pvaule,
                            rrr=rrr)


  ###create S3 object
  microchatInterobj <- data_hier1

  class(microchatInterobj) <- "microchat"

  return(microchatInterobj)
}


"plotMicrochatInterTaxaCircle" <- function(microchatInterobj,
                                           select_group="ct",
                                           color_taxa=colorCustom(50,pal = "gygn")) {

  if (class(microchatInterobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  library(igraph)
  data_hier<-microchatInterobj
  mat<-data_hier$mat
  plot_net_single_tax(mat,select_group=select_group,
                      color_taxa=color_taxa)

  gg<-recordPlot()
  return(gg)
}


"plotMicrochatInterAllTaxaHier" <- function(submchat,
                                            taxa.order.num=3,
                                            char.start=1,
                                            char.stop=10,
                                            color_taxa=colorCustom(50,pal = "gygn"),
                                            export_path="microbial communication analysis/inter-taxa-interaction") {
  export_path<-paste(export_path,"/data_microbiome/microbial communication analysis/inter-taxa-interaction",sep = "")
  dir.create(export_path, recursive = TRUE)

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  dir.create(export_path, recursive = TRUE)
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  abun<-submchat$otu_table
  taxon<-submchat$taxon_table

  taxa.use<-colnames(taxon)
  taxa.num<-taxa.use%>%length()

  for (k in 2:taxa.order.num) {
    name.use<-taxa.use[k]

    data_hier<-net_hierarchy(abun_data=abun,
                             tax_data=taxon,
                             char.start=char.start,
                             char.stop=char.stop,
                             select_taxa=name.use, ###选择可视化的分类
                             diag=TRUE,###自相关是否去除
                             ###是否过滤,过滤会将自相关剔除掉
                             pvaule=0.05, ###高于该显著性，相关性设为0
                             filter=TRUE,
                             rrr=0.1) ###去除低于该相关性绝对值
    mat<-data_hier$mat

    node_list<-data_hier$node_list
    #dev.new()
    par(mfrow=c(2,3))

    for (j in 1:length(node_list)) {
      name_sel<-names(node_list)[j]
      selmax<-0.5*node_list[[j]]%>%rownames()%>%length()
      pdf(paste(export_path,"/",name_sel,"_",name.use,"_plot_net_hierarchy.pdf",sep = ""),
          width=10, height=10)
      plot_net_hierarchy(mat,node_list,abs=TRUE,
                         select_group=name_sel,
                         colors=color_taxa,
                         vertex.receiver = seq(1,selmax),
                         node.size.constant=FALSE,
                         vertex.weight = 1,
                         vertex.weight.max = 5,
                         vertex.size.max = 5,
                         vertex.label.cex = 0.6,
                         alpha.edge=0.6,
                         weight.scale=TRUE,
                         arrowwidth = 0.5,
                         arrow.size = 1.5,
                         margin=0,
                         edge.label.color = "black",
                         edge.label.cex = 0.1,
                         edge.weight.max = 0.1,
                         edge.curved=0,
                         edge.width.max =1.5,
                         title.name = "Micro-occurrence")
      dev.off()
    }
  }

}



"plot_chord" <- function(submchat,
                         all=TRUE,diffHeight=0.06,
                         transparency=0.4,
                         directional=1,
                         full_label=FALSE,
                         char.start=2,
                         char.stop=9,
                         layout = c("regular","baseball"),
                         annotationTrackHeight=0.03,
                         link.arr.type=c("triangle","big.arrow"),
                         select_group=c("ct","cs"),
                         select_taxa=c("Phylum","Class"),
                         top_num=20, small.gap = 1, big.gap = 10,
                         gridcol=randomcolor(20),
                         textsize=1,
                         lgdxpoi=1,lgdypoi=1,lgd_title="lgd") {


  layout <-match.arg(layout)
  link.arr.type<-match.arg(link.arr.type)

  suppressPackageStartupMessages(library(circlize))

  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  otudata<-submchat$otu_table
  taxon_table<-submchat$taxon_table

  split_otu <- filter_merge_sum(otudata=otudata,
                                taxon_table=taxon_table)

  group_all<-split_otu$data_venn%>%rownames()
  group_num<-length(group_all)

  group_num1<-length(group_all)

  abun_ddata<-split_otu$data_with_taxa
  abun_ddata<-column_to_rownames(abun_ddata,var = "name")

  select_order<-match(select_taxa, colnames(abun_ddata))

  {
    abun_alldata<-subset(abun_ddata,select=c(1:group_num,select_order))
  }


  if(all) abun_dddata<-subset(abun_ddata,select=c(1:group_num,select_order))

  if(!all) {
    select_num<-match(select_group,group_all)
    group_num<-length(select_num)
    abun_dddata<-subset(abun_ddata,select=c(select_num,select_order))
  }


  abun_dddata<-unite(abun_dddata,"class",length(colnames(abun_dddata))-1,
                     length(colnames(abun_dddata)),sep = ",")

  temp1<-abun_dddata%>%group_by(class)%>%summarise_all(sum)
  temp2<-temp1
  temp2$sum<-rowSums(temp2[,2:(group_num+1)])
  temp2<-temp2[,-c(2:(group_num+1))]
  temp2$class<-factor(temp2$class,levels = temp2$class)
  temp2<-temp2[order(temp2$sum,decreasing = TRUE),]

  abun_alldata<-unite(abun_alldata,"class",length(colnames(abun_alldata))-1,
                      length(colnames(abun_alldata)),sep = ",")

  temp3<-abun_alldata%>%group_by(class)%>%summarise_all(sum)
  temp4<-temp3
  temp4$sum<-rowSums(temp4[,2:(group_num1+1)])
  temp4<-temp4[,-c(2:(group_num1+1))]
  temp4$class<-factor(temp4$class,levels = temp4$class)
  temp4<-temp4[order(temp4$sum,decreasing = TRUE),]

  temp5<-merge(temp2,temp4,by="class",all = TRUE)
  temp6<-temp5[order(temp5$sum.x,decreasing = TRUE),]
  temp6[is.na(temp6)]<-0
  temp6[,c("tax1","tax2")]<-str_split_fixed(temp6$class,pattern=",",2)

  temp2<-temp6
  if (!full_label) {

    temp2$tax1 <- substr(temp2$tax1,start = char.start,stop = char.stop)
    temp2$tax2 <-substr(temp2$tax2,start = char.start,stop = char.stop)
  }

  ttax1<-select_taxa[1]%>%substr(start = 1,stop = 1)%>%tolower()
  ttax2<-select_taxa[2]%>%substr(start = 1,stop = 1)%>%tolower()



  temp2$tax1[which(temp2$tax1=="")]<-paste(ttax1,"__unclassified",sep = "")

  if (length(which(temp2$tax2=="")) == 1) {
    temp2$tax2[which(temp2$tax2=="")]<-paste(ttax2,"__unclassified",sep = "")
  } else {
    temp2$tax2[which(temp2$tax2=="")]<-paste(ttax2,"__unclassified",1:length(which(temp2$tax2=="")),sep = "")
  }

  tax_num<-length(unique(temp2$tax1))+length(unique(temp2$tax2))
  if (length(gridcol)<tax_num) grid.col<-c(gridcol,gridcol)
  if (length(gridcol)>=tax_num) grid.col <- gridcol

  temp2<-temp2[order(temp2$sum.y,decreasing = TRUE),]


  tax1_color<- grid.col[1:length(unique(temp2$tax1))]
  tax2_color<- grid.col[seq((tax_num-length(unique(temp2$tax2))+1),tax_num)]
  names(tax1_color)<-unique(temp2$tax1)
  names(tax2_color)<-unique(temp2$tax2)
  gg1<-tax1_color%>%data.frame()%>%rownames_to_column(var = "tax1")
  gg2<-tax2_color%>%data.frame()%>%rownames_to_column(var = "tax2")
  colnames(gg1)[2] <- "color1"
  colnames(gg2)[2] <- "color2"

  ddd<-merge(temp2,gg1,by="tax1")
  ddd1<-merge(ddd,gg2,by="tax2")
  ddd1<-ddd1[order(ddd1$sum.x,decreasing = TRUE),]

  if (length(ddd1$tax2)>=top_num) {
    temp3<-ddd1[1:top_num,]
  } else {
    temp3<-ddd1
  }

  df.plot <- temp3

  gridcolor1<-subset(temp3,select=c(tax1,color1))
  gridcolor2<-subset(temp3,select=c(tax2,color2))

  gridcolor1<-gridcolor1%>%distinct(tax1,.keep_all = TRUE)
  gridcolor2<-gridcolor2%>%distinct(tax2,.keep_all = TRUE)

  color1<-gridcolor1$color1
  color2<-gridcolor2$color2
  names(color1)<-gridcolor1$tax1
  names(color2)<-gridcolor2$tax2

  grid.color<-c(color1,color2)

  df.plot$sum<-df.plot$sum.x
  if (layout=="regular" ) {
    df.plot<-subset(df.plot,select=c(tax1,tax2,sum))%>%data.frame()
    order.sector<- c(unique(df.plot$tax2),unique(df.plot$tax1))
  }

  if (layout=="baseball") {
    df.plot<-subset(df.plot,select=c(tax1,tax2,sum))%>%data.frame()
    df.plot<-df.plot[order(df.plot$tax1),]
    order.sector<- c(unique(df.plot$tax2)[length(unique(df.plot$tax2)):1],unique(df.plot$tax1))
  }

  if (length(which(df.plot$sum==0))==0){
    df.plot<- df.plot
  } else {
    df.plot[which(df.plot$sum==0),]$sum<-1
  }

  library(circlize)
  circos.clear()
  par(family="serif",mar=c(0,0,1,1))
  chordDiagram(df.plot, order = order.sector,
               diffHeight = diffHeight,

               grid.col = grid.color, transparency = transparency,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = 0.1),
               small.gap = small.gap, big.gap = big.gap,
               link.visible = TRUE,
               scale = FALSE,
               reduce = -1)

  circos.track(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name,
                facing = "clockwise",family = par("serif"),
                niceFacing = TRUE, adj = c(-0.3, 0.5), cex = textsize)
  }, bg.border = NA)


  circos.clear()


}



"plot_mutichord" <- function(submchat,
                             select_taxa=c("Phylum","Class"),
                             topnum=20,
                             transparency=0.4,
                             gridcol=gridcol,
                             layout = c("regular","baseball"),
                             full_label=FALSE,
                             char.start=2,
                             char.stop=9) {
  layout<-match.arg(layout)

  otudata<-submchat$otu_table
  taxon_table<-submchat$taxon_table

  split_otu <- filter_merge_sum(otudata=otudata,
                                taxon_table=taxon_table)

  group_all<-split_otu$data_venn%>%rownames()
  group_num<-length(group_all)

  ncol<-ncol_layout(group_num)[[1]]%>%max()

  otudatax<-otudata
  otudatax<-otudatax[,-which(colSums(otudatax)==0)]
  group_num<-unique(substr(colnames(otudatax),start = 1,stop = 2))%>%length()
  group_all<-unique(substr(colnames(otudatax),start = 1,stop = 2))

  par(mfrow = c(2, ncol))
  for (kk in 1:group_num) {
    select_group<-group_all[kk]
    plot_chord(submchat,
               all=FALSE,select_group=select_group,###all=false时，select_group开始使用
               transparency=transparency,
               directional=1,
               layout = layout,
               full_label=full_label,
               char.start=char.start,
               char.stop=char.stop,
               annotationTrackHeight=0.03,
               link.arr.type=c("big.arrow"),
               select_taxa=select_taxa,
               top_num=topnum, small.gap = 1, big.gap = 5,
               gridcol= gridcol,
               textsize=1,
               lgdxpoi=1,lgdypoi=0,lgd_title="Taxa")

  }

}


"plot_net_chord" <- function(mat,node_list,weight.scale=TRUE,
                             arrowwidth = 0.2,
                             arrow.size = 0.2,
                             color.edge = c("#b2182b", "#2166ac"),
                             color.use=NULL,
                             edge.label.cex = 0.5,
                             select_group=("st"),
                             colors=randomcolor(50),
                             vertex.receiver = seq(1,5),
                             vertex.weight = NULL, label.edge=TRUE,
                             vertex.weight.max = NULL,
                             vertex.size.max = NULL,
                             edge.weight.max = 1,
                             edge.curved=0,layout = in_circle(),
                             alpha.edge=0.6,margin=0.4,
                             edge.width.max = 3,

                             title.name = title.name,
                             vertex.label.cex = 1) {
  prob.sum<-mat[[match(select_group,names(mat))]]
  vertex.weight<-node_list[[match(select_group,names(node_list))]]$sum
  #prob.sum<-abs(prob.sum)
  vertex.receiver1<-vertex.receiver
  label.edge
  net<-prob.sum
  g <- igraph::graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- igraph::layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0,
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrowwidth
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color,
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)*3

  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * edge.weight.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g,
       edge.curved = 0.3,

       layout = coords_scale,
       margin = margin,
       vertex.label.dist = label.dist,
       vertex.label.degree = label.locs,
       vertex.label.family = "serif",
       edge.label.family = "serif")
}


"plot_net_2chord" <- function(mat,node_list,weight.scale=TRUE,
                              arrowwidth = 1,
                              color.edge = c("#b2182b", "#2166ac"),
                              color.use=NULL,layout = in_circle(),
                              arrow.size = 0.2,
                              edge.label.color = "black",
                              edge.label.cex = 0.5,
                              select_group=("st"),
                              colors=randomcolor(50),
                              vertex.receiver = seq(1,5),
                              vertex.weight = NULL,
                              vertex.weight.max = NULL,
                              vertex.size.max = NULL,
                              edge.weight.max = 1,
                              edge.curved=0,
                              sources.use=NULL,label.edge=TRUE,
                              targets.use=NULL,
                              alpha.edge=0.6,margin=0.4,
                              edge.width.max = 3,
                              title.name = title.name,
                              vertex.label.cex = 1){
  prob.sum<-mat[[match(select_group,names(mat))]]
  vertex.weight<-node_list[[match(select_group,names(node_list))]]$sum

  vertex.receiver1<-vertex.receiver

  net<-prob.sum
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  g <- igraph::graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0,
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g),
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max *
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width <- arrowwidth
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,
                                                                             1]], alpha.edge)
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved,
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist,
       vertex.label.degree = label.locs, vertex.label.family = "serif",
       edge.label.family = "serif")


}




"plotMicrochatInterAllTaxaChord" <- function(microchatInterobj,
                                             edge.curved=0,
                                             color.edge.left = c("#b2182b", "#2166ac"),
                                             color.edge.right=colorCustom(100,pal = "gygn"),
                                             select_group=("ct"),
                                             alpha.edge=1) {

  if (class(microchatInterobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  mat<-microchatInterobj$mat
  node_list<-microchatInterobj$node_list
  groupnum<-length(mat)
  groupname<-names(mat)

  for (ttk in 1:groupnum) {
    matnew<-mat[[ttk]]
    node_list_new<-node_list[[ttk]]
    groupname_new<-groupname[ttk]

    if (nrow(matnew)==0 | nrow(node_list_new)==0) {
      stop("\n",groupname_new," showed matrix and node list with null data !!!")
      message("Please re-calaculate the threshold of corrleation !!!")
    }
  }


  dev.off()
  par(mfrow = c(1,2), xpd=NA,family="serif")
  plot_net_chord(mat,node_list,
                 color.edge = color.edge.left,
                 color.use=color.edge.right,
                 select_group=select_group,
                 colors=color.edge.right,
                 vertex.receiver = seq(1,5),
                 vertex.weight = 10,
                 vertex.weight.max = 40,
                 vertex.size.max = 40,
                 vertex.label.cex = 0.6,
                 alpha.edge=alpha.edge,
                 weight.scale=TRUE,
                 arrowwidth = 0.2,
                 arrow.size = 0.2,
                 margin=0,
                 label.edge = FALSE,
                 edge.label.cex = 0,
                 edge.weight.max = 2,
                 edge.curved=edge.curved,
                 edge.width.max =2,
                 title.name = "net_hierarchy")

  plot_net_2chord(mat,node_list,
                  color.edge =color.edge.left,
                  color.use=color.edge.right,
                  select_group=select_group,
                  colors=color.edge.right,
                  vertex.receiver = seq(1,5),
                  vertex.weight = 10,
                  vertex.weight.max = 40,
                  vertex.size.max = 40,
                  vertex.label.cex = 0.6,
                  alpha.edge=alpha.edge,label.edge = FALSE,
                  weight.scale=TRUE,
                  arrowwidth = 0.2,
                  arrow.size = 0.2,
                  margin=0,
                  edge.label.cex = 0,
                  edge.weight.max = 2,
                  edge.curved=edge.curved,
                  edge.width.max =2,
                  title.name = "net_hierarchy")
  gg<-recordPlot()

  return(gg)
}



"plot_net_heatmap" <- function(mat,
                               comparison = c(1, 2),
                               select_group="st",
                               measure = NULL,
                               signaling = NULL,
                               all.taxa.name=TRUE,
                               slot.name = c("netP", "net"),
                               color.use = NULL,
                               color.heatmap = c("#2166ac", "#b2182b","GnBu","OrRd"),
                               title.name = NULL,
                               width = NULL,
                               height = NULL,
                               legend.name=NULL,
                               font.size = 8,
                               font.size.title = 10,
                               cluster.rows = FALSE,
                               cluster.cols = FALSE,
                               sources.use = NULL,
                               targets.use = NULL,
                               remove.isolate = FALSE,
                               row.show = NULL,
                               col.show = NULL) {

  suppressPackageStartupMessages(library(ComplexHeatmap))
  all.name<-NULL
  for (tt in 1:length(mat)) {
    name<-names(mat)[tt]
    row.name<-row.names(mat[[tt]])%>%data.frame()

    {
      all.name<-rbind(all.name,row.name)
    }
  }
  name.use<-unique(all.name$.)

  mat.use<-matrix(0,length(unique(all.name$.)),length(unique(all.name$.)))
  row.names(mat.use)<-unique(all.name$.)
  colnames(mat.use)<-rownames(mat.use)


  library(ComplexHeatmap)
  library(CellChat)
  net<-mat[[match(select_group,names(mat))]]

  if (all.taxa.name) {
    for (i in 1:nrow(net)) {
      for (j in 1:nrow(net)) {
        rowss<-row.names(net)[i]
        colss<-colnames(net)[j]

        hr<-match(rowss,row.names(mat.use))
        cv<-match(colss,colnames(mat.use))
        mat.use[hr,cv]<-net[i,j]
      }
    }
    net<-mat.use
  } else {
    net<-net
  }

  if (is.null(color.use)) {
    color.use <- color.heatmap[1:ncol(net)]
  } else {
    color.use <- color.use[1:ncol(net)]
  }

  names(color.use) <- colnames(net)
  if (!is.null(row.show)) {
    net <- net[row.show, ]
  }
  if (!is.null(col.show)) {
    net <- net[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(net) < 0) {
    color.heatmap.use = CellChat::colorRamp3(c(min(net), 0, max(net)),
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(net, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(net, na.rm = T))) + 1), 0, round(max(net,
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                     max(net, na.rm = T))) + 1))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = CellChat::colorRamp3(c(0, min(net), max(net)),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = CellChat::colorRamp3(c(min(net), max(net)),
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(net, na.rm = T), digits = nchar(sub(".*\\.(0*).*",
                                                                      "\\1", min(net, na.rm = T))) + 1), round(max(net,
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1",
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(net))
  rownames(df) <- colnames(net)
  col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df, col = list(group = color.use),
                                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df, col = list(group = color.use),
                                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = ComplexHeatmap::rowAnnotation(Strength = ComplexHeatmap::anno_barplot(rowSums(abs(net)),
                                                                              border = FALSE, gp = grid::gpar(fill = color.use, col = color.use)),
                                      show_annotation_name = FALSE)
  ha2 = ComplexHeatmap::HeatmapAnnotation(Strength = ComplexHeatmap::anno_barplot(colSums(abs(net)),
                                                                                  border = FALSE, gp = grid::gpar(fill = color.use, col = color.use)),
                                          show_annotation_name = FALSE)
  if (sum(abs(net) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    net[net == 0] <- NA
  }
  ht1 = ComplexHeatmap::Heatmap(net, col = color.heatmap.use, na_col = "white",
                                name = legend.name, bottom_annotation = col_annotation,
                                left_annotation = row_annotation, top_annotation = ha2,
                                right_annotation = ha1, cluster_rows = cluster.rows,
                                cluster_columns = cluster.rows, row_names_side = "left",
                                row_names_rot = 0, row_names_gp = grid::gpar(fontsize = font.size,fontfamily="serif"),
                                column_names_gp = grid::gpar(fontsize = font.size,fontfamily="serif"), column_title = title.name,
                                column_title_gp = grid::gpar(fontsize = font.size.title),
                                column_names_rot = 90, row_title = "Taxa",
                                row_title_gp = grid::gpar(fontsize = font.size.title), row_title_rot = 90,
                                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8,fontfamily="serif",
                                                                                  fontface = "plain"), title_position = "leftcenter-rot",
                                                            border = NA, legend_height = unit(20, "mm"), labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                            grid_width = unit(2, "mm")))

  return(ht1)
}


"plotMicrochatInterAllTaxaHeatmap" <- function(microchatInterobj,
                                               all.taxa.name=TRUE,
                                               font.size = 8,
                                               font.size.title = 10,
                                               color.heatmap = colorCustom(50,pal = "gygn"),
                                               export_path="microbial communication analysis/inter-taxa-interaction") {

  export_path<-paste(export_path,"/data_microbiome/microbial communication analysis/inter-taxa-interaction",sep = "")
  dir.create(export_path, recursive = TRUE)

  if (class(microchatInterobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  mat<-microchatInterobj$mat

  ht<-list()
  for (t in 1:length(mat)) {
    name_sel<-names(mat)[t]
    ht[[t]]<-plot_net_heatmap(mat,select_group=name_sel,
                              all.taxa.name=all.taxa.name,
                              #comparison = c(1, 2),
                              measure = NULL,
                              signaling = NULL,
                              slot.name =  "netp",
                              color.use = NULL,
                              color.heatmap = color.heatmap,
                              title.name = NULL,
                              width = NULL,
                              legend.name=NULL,
                              height = NULL,
                              font.size = font.size,
                              font.size.title = font.size.title,
                              cluster.rows = FALSE,
                              cluster.cols = FALSE,
                              sources.use = NULL,
                              targets.use = NULL,
                              remove.isolate = FALSE,
                              row.show = NULL,
                              col.show = NULL)



    a1<-ggplotify::as.ggplot(ht[[t]])

    ggsave(paste(export_path,"/(",name_sel,") heatmap.pdf",sep = ""),a1)
  }
  dev.off()

  return(ht)
}



