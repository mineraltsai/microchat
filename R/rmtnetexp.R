"plot_rmnet" <- function(g,split_otu,thres,num= length(g),rows=num/cols,cols=num/rows,display_filter=TRUE,radius=1.15,
                         vertex_size=2,
                         edge_size=2,
                         layout="fr",circle_plot=TRUE) {
  col_g <- "#C1C1C1"
  colors <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
              "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")

  result<-thres
  par(mfrow=c(rows,cols),mar=c(1,1,3,1),font.main=4)
  for(k in 1:length(g)){
    name<-names(g)[k]
    g1 <- g[[k]]

    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    set.seed(007)

    V(g1)$modularity <- membership(cluster_fast_greedy(g1))
    modulenum<-length(cluster_fast_greedy(g1))
    V(g1)$label <- V(g1)$name
    V(g1)$label <- NA
    modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T)

    top_num <- 18
    modu_name <- names(modu_sort[1:18])
    modu_cols <- colors[1:length(modu_name)]
    names(modu_cols) <- modu_name
    V(g1)$color <- V(g1)$modularity
    V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
    V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)],modu_name)]
    V(g1)$frame.color <- V(g1)$color

    E(g1)$color <- col_g
    for ( i in modu_name){
      col_edge <- colors[which(modu_name==i)]
      otu_same_modu <-V(g1)$name[which(V(g1)$modularity==i)]
      E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
    }

    split_otus<-split_otu[[k]] %>% data.frame()
    split_otus$name<-rownames(split_otus)
    add_num<-length(split_otus$name)-length(V(g1)$name)
    if (display_filter) g2<-g1 %>% add_vertices(add_num, color = "black")
    if (!display_filter) g2<-g1


    if (layout=="kk") sub_net_layout <- layout_with_kk(g2, maxiter = vcount(g2))
    if (layout=="fr") sub_net_layout <- layout_with_fr(g2, niter=999,grid ="nogrid",start.temp = sqrt(vcount(g2)))
    if (layout=="dh") sub_net_layout <- layout_with_dh(g2, maxiter=999)
    if (layout=="gem") sub_net_layout <- layout_with_gem(g2, maxiter = 40 * vcount(g2)^2,temp.max = vcount(g2),temp.min = 1/10,temp.init = sqrt(vcount(g2)))
    if (layout=="graphopt") sub_net_layout <- layout_with_graphopt(g2, niter=999)
    if (layout=="lgl") sub_net_layout <- layout_with_lgl(g2, maxiter=999)
    if (layout=="mds") sub_net_layout <- layout_with_mds(g2)
    if (layout=="circle") sub_net_layout <- layout_in_circle(g2)
    if (layout=="star") sub_net_layout <- layout_as_star(g2)
    if (layout=="tree") sub_net_layout <- layout_as_tree(g2,circular = TRUE,root=c(1,11),
                                                         rootlevel=c(2,1))
    if (layout=="nice") sub_net_layout <- layout_nicely(g2, dim = 3)
    if (layout=="grid") sub_net_layout <- layout_on_grid(g2)
    if (layout=="sphere") sub_net_layout <- layout_on_sphere(g2)
    if (layout=="rand") sub_net_layout <- layout_randomly(g2)


    plot(g2,layout=sub_net_layout, edge.color = E(g2)$color,vertex.size=vertex_size,edge.curved=.1,edge.width=edge_size)
    if (circle_plot)  plotrix::draw.circle(0,0,radius,border="grey50",lty = 2)
    if (display_filter) title(main = paste0(name,":",'Nodes=',length(V(g1)$name),', ',
                                            'Edges=',nrow(data.frame(as_edgelist(g1)))," (",modulenum,")",
                                            "\n",'Pre-filtering nodes=',length(V(g2)$name),
                                            "\n","selected threshold: ", result[[k]]))
    if (!display_filter) title(main = paste0(name,":",'Nodes=',length(V(g1)$name),', ',
                                             'Edges=',nrow(data.frame(as_edgelist(g1)))," (",modulenum,")",
                                             "\n","selected threshold: ", result[[k]]))

  }

  if (display_filter)   cat("\n","Filtered")
  if (!display_filter)   cat("\n","Not filtered")
  cat("\n","\n","Visualize ",num," plots"," with ",rows," rows and ",cols," columns","\n","\n")
  cat("\n","parameters：","num (plot number)，rows (rows of layout)，cols (columns of layout)","\n")

}

"get.MutiPie"<-function(g,taxon) {
  taxon$name<-rownames(taxon)
  m.node_phylum.m<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    node_table = as.data.frame(vertex.attributes(g1))
    node_table$degree = degree(g1)
    node_table$closeness=closeness(g1)
    node_table$coreness<-coreness(g1)
    node_table$eccentricity<-eccentricity(g1, vids = V(g1), mode = "all")
    node_table$strength<-strength(g1,vids = V(g1),mode = "all",loops = TRUE, weights = NULL)
    node_table$betweenness=betweenness(g1)
    node_table$evcent=evcent(g1)$vector
    module_membership = membership(cluster_fast_greedy(g1))
    node_table$module= membership(cluster_fast_greedy(g1))
    edge_table = as.data.frame(edge.attributes(g1))
    edge_names = unlist(strsplit(attr(E(g1), "vnames"), "|", fixed=T))
    edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
    edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]
    edge_table$cor[which(edge_table$correlation < 0)] = "-1"
    edge_table$cor[which(edge_table$correlation > 0)] = "1"
    node_table<-merge(node_table,taxon,by="name")
    node_table$id<-node_table$name

    node_phylum.m<-data.frame()
    for (moduled in unique(node_table$module)) {
      node_module <- subset(node_table, module == moduled)
      node_phylum <- data.frame(table(node_module$Phylum))

      node_phylum$module<-moduled

      node_phylum.m<-rbind(node_phylum.m,node_phylum)
    }
    node_phylum.m$group<-names(g)[i]
    m.node_phylum.m<-rbind(m.node_phylum.m,node_phylum.m)

  }
  return(data=m.node_phylum.m)
}

"plot_rmnet_modnode_barplot_export" <- function(g,taxon,
                                                geom.align=c("center","bottom"),
                                                point.direction=c("all","top","bottom","none"),
                                                point.color.uniform=TRUE,
                                                point.density=3,
                                                point.panel.width=0.3,
                                                point.size=0.2,
                                                point.text.size=3,
                                                module.text.size=10,
                                                aspect.ratio=0.8,
                                                color_module=colorCustom(50,pal = "gygn"),
                                                color_group=colorCustom(10,pal = "gygn"),
                                                export_path ='network/module_nodenum_barplot') {
  geom.align<-match.arg(geom.align)
  point.direction<-match.arg(point.direction)

  export_path<-paste(export_path,"/microbial network analysis/module_nodenum_barplot",sep = "")
  dir.create(export_path, recursive = TRUE)

  taxon$name<-rownames(taxon)
  nodenum<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    node_table = as.data.frame(vertex.attributes(g1))
    node_table$degree = degree(g1)
    node_table$closeness=closeness(g1)
    node_table$coreness<-coreness(g1)
    node_table$eccentricity<-eccentricity(g1, vids = V(g1), mode = "all")
    node_table$strength<-strength(g1,vids = V(g1),mode = "all",loops = TRUE, weights = NULL)
    node_table$betweenness=betweenness(g1)
    node_table$evcent=evcent(g1)$vector
    module_membership = membership(cluster_fast_greedy(g1))
    node_table$module= membership(cluster_fast_greedy(g1))
    edge_table = as.data.frame(edge.attributes(g1))
    edge_names = unlist(strsplit(attr(E(g1), "vnames"), "|", fixed=T))
    edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
    edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]
    edge_table$cor[which(edge_table$correlation < 0)] = "-1"
    edge_table$cor[which(edge_table$correlation > 0)] = "1"
    node_table<-merge(node_table,taxon,by="name")
    node_table$id<-node_table$name

    for (moduled in unique(node_table$module)) {
      node_module <- subset(node_table, module == moduled)
      node_phylum <- data.frame(table(node_module$Phylum))

      num<-length(node_module$name)%>%data.frame()
      num$module<-unique(node_module$module)
      num$group<-i
      colnames(num)[1]<-"value"

      {
        nodenum<-rbind(nodenum,num)
      }

    }

  }

  nodenum$value.y<-nodenum$value+4
  nodenum$group<-factor(nodenum$group,levels =unique(nodenum$group))

  colors<-color_group
  group_num<-length(g)

  if(group_num>length(colors)) {
    message("\n","Please provide more colors","\n")
    cat("\n","Colors provided cannot meet the requirement, hence replicate all colors","\n")
    color1<-rep(colors,4)
    colors<-color1[1:group_num]
  } else {
    colors<-colors[1:group_num]
  }

  color<-colors
  pnode<-list()

  if (geom.align=="bottom") {
    for (r in 1:length(g)) {
      netname<-names(g)[r]
      dou <- function(x){
        ifelse(x%%2 ==0,F,T)
      }

      ji <- unique(nodenum[nodenum$group==r,]$module)[dou(unique(nodenum[nodenum$group==r,]$module))]
      ou <- unique(nodenum[nodenum$group==r,]$module)[!dou(unique(nodenum[nodenum$group==r,]$module))]

      p<-ggplot()
      if (point.direction=="all") {

        data_line<-nodenum[nodenum$group==r,]
        data_line$newvalue<-data_line$value%>%as.integer()

        dataa<-data.frame()
        for (tt in data_line$module) {

          datab<-data.frame(value=rep(data_line[data_line$module==tt,]$value,
                                      data_line[data_line$module==tt,]$newvalue),
                            module=rep(data_line[data_line$module==tt,]$module,
                                       data_line[data_line$module==tt,]$newvalue),
                            group=rep(data_line[data_line$module==tt,]$group,
                                      data_line[data_line$module==tt,]$newvalue),
                            value.y=rep(data_line[data_line$module==tt,]$value.y,
                                        data_line[data_line$module==tt,]$newvalue),
                            newvalue=rep(data_line[data_line$module==tt,]$newvalue,
                                         data_line[data_line$module==tt,]$newvalue),
                            ponum=rep(1:data_line[data_line$module==tt,]$newvalue))


          {
            dataa<-rbind(dataa,datab)
          }

        }

        datac<-data.frame()
        for (kt in 1:point.density) {
          dataa<-dataa
          {
            datac<-rbind(datac,dataa)
          }
        }
        datacx<-datac

        datacx$modname<-"modname"
        datacx<-unite(datacx,"modname",modname,module,sep = "",remove = FALSE)


        if (!is.null(color_module)) {
          datacx<-datacx[order(datacx$module),]
          modname<-datacx$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group,color=modname),
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            scale_color_manual(values = color_module)

        } else {
          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group),color=colors[r],
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width))
        }
      }

      if (point.direction=="none") {
        dat_non<-nodenum[nodenum$group==r,]
        dat_non$modname<-"modname"
        dat_non<-unite(dat_non,"modname",modname,module,sep = "",remove = FALSE)

        if (!is.null(color_module)) {
          dat_non<-dat_non[order(dat_non$module),]
          modname<-dat_non$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          p<-p +
            geom_point(data=dat_non,
                       aes(x=module, y=value, group=group,color=modname))+
            scale_color_manual(values = color_module)+
            geom_col(data = dat_non,
                     aes(x=module,y=value, group=group,fill=modname))+
            scale_fill_manual(values = color_module)

        } else {
          p<-p +
            geom_point(data=dat_non,
                       aes(x=module, y=value, group=group),color=colors[r]) +
            geom_col(data = dat_non,
                     aes(x=module,y=value, group=group),fill=colors[r])
        }
      }

      p<-p +
        geom_text(data = nodenum[nodenum$group==r,],aes(x = module,y = value.y,label = value),size = point.text.size,color = "black",family = "serif") +
        ylab("Node number") +
        xlab("Module")+
        annotate(geom = "text", family="serif",
                 x = ji, y = -max(nodenum[nodenum$group==1,]$value)/32, size=point.text.size,
                 label = ji) +
        annotate(geom = "text",  family="serif",
                 x = ou, y = -max(nodenum[nodenum$group==1,]$value)/32, size=point.text.size,
                 label = ou) +
        annotate(geom = "text",  family="serif",
                 x = max(nodenum[nodenum$group==r,]$module)-3,
                 y = max(nodenum[nodenum$group==r,]$value)-5,
                 size=module.text.size,
                 label = netname) +
        coord_cartesian(ylim = c(0,max(nodenum[nodenum$group==r,]$value)+10), expand = FALSE, clip = "off") +
        theme_bw() +
        theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
              strip.background = element_rect(fill = "white"),
              #strip.background =  element_blank(),
              legend.position = "none",aspect.ratio = aspect.ratio)
      pnode[[r]]<-p
      ggsave(paste(export_path,'/',geom.align,"_",point.direction,'_node number (',netname,').pdf', sep = ''),
             pnode[[r]])
    }
  }
  if (geom.align=="center") {
    nodenumk<-nodenum
    nodenumk$value<-nodenumk$value/2
    nodenumc<-nodenumk
    nodenumc$value<-nodenumc$value*(-1)
    nodenumm<-rbind(nodenumk,nodenumc)
    nodenumm$value.x<-nodenumm$value.y-2
    nodenumm$value.y<-nodenumm$value+2


    for (r in 1:length(g)) {
      netname<-names(g)[r]
      dou <- function(x){
        ifelse(x%%2 ==0,F,T)
      }

      ji <- unique(nodenumm[nodenumm$group==r,]$module)[dou(unique(nodenumm[nodenumm$group==r,]$module))]
      ou <- unique(nodenumm[nodenumm$group==r,]$module)[!dou(unique(nodenumm[nodenumm$group==r,]$module))]

      data_line<-nodenumm[nodenumm$group==r & nodenumm$value>0,]
      data_line$newvalue<-data_line$value%>%as.integer()

      dataa<-data.frame()
      for (tt in data_line$module) {

        datab<-data.frame(value=rep(data_line[data_line$module==tt,]$value,
                                    data_line[data_line$module==tt,]$newvalue),
                          module=rep(data_line[data_line$module==tt,]$module,
                                     data_line[data_line$module==tt,]$newvalue),
                          group=rep(data_line[data_line$module==tt,]$group,
                                    data_line[data_line$module==tt,]$newvalue),
                          value.y=rep(data_line[data_line$module==tt,]$value.y,
                                      data_line[data_line$module==tt,]$newvalue),
                          value.x=rep(data_line[data_line$module==tt,]$value.x,
                                      data_line[data_line$module==tt,]$newvalue),
                          newvalue=rep(data_line[data_line$module==tt,]$newvalue,
                                       data_line[data_line$module==tt,]$newvalue),
                          ponum=rep(1:data_line[data_line$module==tt,]$newvalue))


        {
          dataa<-rbind(dataa,datab)
        }

      }

      datac<-data.frame()
      for (kt in 1:point.density) {
        dataa<-dataa
        {
          datac<-rbind(datac,dataa)
        }
      }
      datacx<-datac

      data_line_bot<-nodenumm[nodenumm$group==r & nodenumm$value<0,]
      data_line_bot$newvalue<-data_line_bot$value%>%as.integer()

      dataa_bot<-data.frame()
      for (tt in data_line_bot$module) {

        datab_bot<-data.frame(value=rep(data_line_bot[data_line_bot$module==tt,]$value,
                                        data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              module=rep(data_line_bot[data_line_bot$module==tt,]$module,
                                         data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              group=rep(data_line_bot[data_line_bot$module==tt,]$group,
                                        data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              value.y=rep(data_line_bot[data_line_bot$module==tt,]$value.y,
                                          data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              value.x=rep(data_line_bot[data_line_bot$module==tt,]$value.x,
                                          data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              newvalue=rep(data_line_bot[data_line_bot$module==tt,]$newvalue,
                                           data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs()),
                              ponum=rep(1:(data_line_bot[data_line_bot$module==tt,]$newvalue%>%abs())))

        {
          dataa_bot<-rbind(dataa_bot,datab_bot)
        }

      }

      datac_bot<-data.frame()
      for (kt in 1:point.density) {
        dataa_bot<-dataa_bot
        {
          datac_bot<-rbind(datac_bot,dataa_bot)
        }
      }

      datacx_bot<-datac_bot
      datacx_bot$ponum<-datacx_bot$ponum*(-1)

      datacx$modname<-"modname"
      datacx<-unite(datacx,"modname",modname,module,sep = "",remove = FALSE)

      datacx_bot$modname<-"modname"
      datacx_bot<-unite(datacx_bot,"modname",modname,module,sep = "",remove = FALSE)

      datacx_bot<-datacx_bot[order(datacx_bot$module),]
      datacx<-datacx[order(datacx$module),]
      p<-ggplot()
      if (point.direction=="all") {

        if (!is.null(color_module)) {
          modname<-datacx$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          modname1<-datacx_bot$modname%>%unique()
          color_module1<-rep(color_module,5)[1:length(modname1)]
          names(color_module1)<-modname1

          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group,color=modname),
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            scale_color_manual(values = color_module)

          if (point.color.uniform) {
            p<-p+ggnewscale::new_scale_color()+
              geom_point(data=datacx_bot,
                         aes(x=module, y=ponum, group=group,color=modname),
                         alpha=1,size=point.size,
                         position = position_jitter(width = point.panel.width, height = point.panel.width))+
              scale_color_manual(values = color_module1)
          } else {
            p<-p+ggnewscale::new_scale_color()+
              geom_point(data=datacx_bot,
                         aes(x=module, y=ponum, group=group),color=colors[r],
                         alpha=1,size=point.size,
                         position = position_jitter(width = point.panel.width, height = point.panel.width))+
              scale_color_manual(values = color_module1)
          }


        } else {


          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group),color=colors[r],
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            geom_point(data=datacx_bot,
                       aes(x=module, y=ponum, group=group,color=modname),color=colors[r],
                       alpha=1,size=point.size,
                       position = position_jitter(width = point.panel.width, height = point.panel.width))
        }
      }
      if (point.direction=="top") {
        if (!is.null(color_module)) {
          modname<-datacx$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value<0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)
          dat_bot<-dat_bot[order(dat_bot$module),]
          modname1<-dat_bot$modname%>%unique()
          color_module1<-rep(color_module,5)[1:length(modname1)]
          names(color_module1)<-modname1

          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group,color=modname),
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            scale_color_manual(values = color_module)


          if (point.color.uniform) {
            p<-p+ggnewscale::new_scale_color()+
              geom_col(data = dat_bot,
                       aes(x=module,y=value, group=group,fill=modname))+
              scale_fill_manual(values = color_module1)
          } else {
            p<-p+ggnewscale::new_scale_color()+
              geom_col(data = dat_bot,
                       aes(x=module,y=value, group=group),fill=colors[r])+
              scale_fill_manual(values = color_module1)
          }


        } else {

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value<0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)


          p<-p+geom_point(data=datacx,
                          aes(x=module, y=ponum, group=group),color=colors[r],
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            geom_col(data = dat_bot,
                     aes(x=module,y=value, group=group),fill=colors[r])



        }
      }
      if (point.direction=="bottom") {

        if (!is.null(color_module)) {
          modname<-datacx_bot$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value>0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)
          dat_bot<-dat_bot[order(dat_bot$module),]
          modname1<-dat_bot$modname%>%unique()
          color_module1<-rep(color_module,5)[1:length(modname1)]
          names(color_module1)<-modname1

          p<-p+geom_col(data = dat_bot,
                        aes(x=module,y=value, group=group,fill=modname))+
            scale_fill_manual(values = color_module1)


          if (point.color.uniform) {
            p<-p+ggnewscale::new_scale_color()+
              geom_point(data=datacx_bot,
                         aes(x=module, y=ponum, group=group,color=modname),
                         alpha=1,size=point.size,
                         position = position_jitter(width = point.panel.width, height = point.panel.width)) +
              scale_color_manual(values = color_module)

          } else {
            p<-p+ggnewscale::new_scale_color()+
              geom_point(data=datacx_bot,
                         aes(x=module, y=ponum, group=group),color=colors[r],
                         alpha=1,size=point.size,
                         position = position_jitter(width = point.panel.width, height = point.panel.width)) +
              scale_color_manual(values = color_module)
          }

        } else {

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value>0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)


          p<-p+geom_point(data=datacx_bot,
                          aes(x=module, y=ponum, group=group),color=colors[r],
                          alpha=1,size=point.size,
                          position = position_jitter(width = point.panel.width, height = point.panel.width)) +
            geom_col(data = dat_bot,
                     aes(x=module,y=value, group=group),fill=colors[r])



        }

      }
      if (point.direction=="none") {

        if (!is.null(color_module)) {
          dat_top<-nodenumm[nodenumm$group==r & nodenumm$value>0,]
          dat_top$modname<-"module"
          dat_top<-unite(dat_top,"modname",modname,module,sep = "",remove = FALSE)
          dat_top<-dat_top[order(dat_top$module),]
          modname<-dat_top$modname%>%unique()
          color_module<-rep(color_module,5)[1:length(modname)]
          names(color_module)<-modname

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value<0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)


          dat_bot<-dat_bot[order(dat_bot$module),]
          modname1<-dat_bot$modname%>%unique()
          color_module1<-rep(color_module,5)[1:length(modname1)]
          names(color_module1)<-modname1

          p<-p+geom_col(data = dat_top,
                        aes(x=module,y=value, group=group,fill=modname))+
            scale_color_manual(values = color_module)



          if (point.color.uniform) {
            p<-p+ggnewscale::new_scale_color()+
              geom_col(data = dat_bot,
                       aes(x=module,y=value, group=group,fill=modname))+
              scale_fill_manual(values = color_module1)

          } else {
            p<-p+ggnewscale::new_scale_color()+
              geom_col(data = dat_bot,
                       aes(x=module,y=value, group=group),fill=colors[r])+
              scale_fill_manual(values = color_module1)
          }


        } else {
          dat_top<-nodenumm[nodenumm$group==r & nodenumm$value>0,]
          dat_top$modname<-"module"
          dat_top<-unite(dat_top,"modname",modname,module,sep = "",remove = FALSE)

          dat_bot<-nodenumm[nodenumm$group==r & nodenumm$value<0,]
          dat_bot$modname<-"module"
          dat_bot<-unite(dat_bot,"modname",modname,module,sep = "",remove = FALSE)

          p<-p+geom_col(data = dat_top,
                        aes(x=module,y=value, group=group),fill=colors[r])+
            geom_col(data = dat_bot,
                     aes(x=module,y=value, group=group),fill=colors[r])

        }

      }


      p<-p+ geom_text(data = nodenumm[nodenumm$group==r,],
                      aes(x = module,y = 0,label = value.x),size = point.text.size,color = "black",family = "serif") +
        ylab("Node number") +
        xlab("Module")+
        #facet_wrap(.~group,ncol = 3,scales = "free_x") +
        #stat_compare_means(data=bb1,aes(x=Group, y=value, group=variable),method = "anova", label.y = 0)+
        annotate(geom = "text",  family="serif",
                 x = ji, y = min(nodenumm[nodenumm$group==r,]$value)-5, size=point.text.size,
                 label = ji) +
        annotate(geom = "text",  family="serif",
                 x = ou, y = min(nodenumm[nodenumm$group==r,]$value)-5, size=point.text.size,
                 label = ou) +
        annotate(geom = "text",  family="serif",
                 x = max(nodenumm[nodenumm$group==r,]$module)-3,
                 y = max(nodenumm[nodenumm$group==r,]$value)-5,
                 size=module.text.size,
                 label = netname) +
        coord_cartesian(ylim = c(-(max(nodenumm[nodenumm$group==r,]$value)+10),max(nodenumm[nodenumm$group==r,]$value)+10), expand = FALSE, clip = "off") +
        theme_bw() +
        theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
              #panel.background = element_rect(fill="grey70"),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              strip.text = element_text(colour = "red",size = 10,face = "bold",family = "serif"),
              strip.background = element_rect(fill = "white"),
              #strip.background =  element_blank(),
              legend.position = "none",aspect.ratio = aspect.ratio)
      pnode[[r]]<-p
      ggsave(paste(export_path,'/',geom.align,"_",point.direction,'_node number (',netname,').pdf', sep = ''),
             pnode[[r]])
    }
  }

  cat("\n","Barplot from all network modules has been exported to",export_path,
      "，which will be used for observing functional redundancy.","\n")

}


"plot_rmnet_cyto_pie_export" <- function(g,taxon,
                                         color_text="#91091A",
                                         color_taxa=colorCustom(40,pal = "cyto_pie"),
                                         export_path ='network/cyto_module_pie') {

  export_path<-paste(export_path,"/microbial network analysis/cyto_module_pie",sep = "")
  dir.create(export_path, recursive = TRUE)

  taxon$name<-rownames(taxon)
  color <- color_taxa
  phy_num<-length(taxon$Phylum%>%unique())

  if(phy_num>length(color)) {
    message("\n","Please provide more colors","\n")
    cat("\n","Colors provided cannot meet the requirement, hence replicate all colors","\n")
    color1<-rep(color,4)
    color<-color1[1:phy_num]
    names(color)<-taxon$Phylum%>%unique()
  } else {
    color<-color[1:phy_num]
    names(color)<-taxon$Phylum%>%unique()
  }




  for(i in 1:length(g)){

    dir.create(paste(export_path,'/',i,sep = ""), recursive = TRUE)

    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    node_table = as.data.frame(vertex.attributes(g1))
    node_table$degree = degree(g1)
    node_table$closeness=closeness(g1)
    node_table$coreness<-coreness(g1)
    node_table$eccentricity<-eccentricity(g1, vids = V(g1), mode = "all")
    node_table$strength<-strength(g1,vids = V(g1),mode = "all",loops = TRUE, weights = NULL)
    node_table$betweenness=betweenness(g1)
    node_table$evcent=evcent(g1)$vector
    module_membership = membership(cluster_fast_greedy(g1))
    node_table$module= membership(cluster_fast_greedy(g1))
    edge_table = as.data.frame(edge.attributes(g1))
    edge_names = unlist(strsplit(attr(E(g1), "vnames"), "|", fixed=T))
    edge_table$node1 = edge_names[seq(from = 1, to = (length(edge_names) - 1), by = 2)]
    edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]
    edge_table$cor[which(edge_table$correlation < 0)] = "-1"
    edge_table$cor[which(edge_table$correlation > 0)] = "1"
    node_table<-merge(node_table,taxon,by="name")
    node_table$id<-node_table$name

    for (moduled in unique(node_table$module)) {
      node_module <- subset(node_table, module == moduled)
      node_phylum <- data.frame(table(node_module$Phylum))

      num<-length(node_module$name)%>%data.frame()
      num$module<-unique(node_module$module)
      num$group<-i
      colnames(num)[1]<-"value"

      p <- ggplot(node_phylum, aes(x = '', y = Freq, fill = Var1)) +
        geom_bar(stat = 'identity', width = 0.2) +
        coord_polar(theta = 'y') +
        scale_fill_manual(limits = names(color), values = color) +
        geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                      x= 1.05,label=scales::percent(Freq/sum(node_phylum$Freq))),
                  size=4,colour="black",family = "serif")+
        #geom_text(aes(label=Var1),size=2,position = position_stack(vjust = 0.5),colour="white")+
        geom_text(aes(y=rev(rev(Freq/2)+c(0,cumsum(rev(Freq))[-length(Freq)])),
                      x= 1.15,label=Var1),
                  size=4,color=color_text,family = "serif")+

        theme(panel.grid = element_blank(), panel.background = element_blank(),
              axis.text.x = element_blank()) +
        labs(x = '', y = '')+theme_void()+
        theme(legend.position = "none")
      ggsave(paste(export_path,"/",i,'/',moduled,'.pdf', sep = ''),
             p)
    }


  }

  cat("\n","Pie plot used for making band plot and filling network from cytoscape have been exported to ",export_path,"\n")
}



"zipi" <- function(g, A=NULL, weighted=FALSE) {

  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse=FALSE, names=TRUE, attr='weight')
    } else {
      A <- as_adj(g, sparse=FALSE, names=TRUE)
    }
  }
  memb <- vertex_attr(g, "modularity")
  N <- max(memb)
  nS <- tabulate(memb)
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  names(z) <- names(Ki) <- rownames(A)
  for (S in seq_len(N)) {
    x <- rowSums(A[memb == S, memb == S])
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S]-1))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  df <- data.frame(Ki,z,row.names = names(Ki))
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse=FALSE, attr='weight')
    } else {
      A <- as_adj(g, sparse=FALSE)
    }
  }
  memb <- vertex_attr(g, "modularity")
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  pi <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  names(pi) <- rownames(A)
  pi<-pi%>%data.frame()
  colnames(pi)<-"pi"
  colnames(df)<-c("ki","zi")
  pi$name<-rownames(pi)
  df$name<-rownames(df)
  res<-merge(df,pi,by="name")
  return(res)
}



"plot_TreatModuleEigen_heatmap" <- function(g,abun,
                                            heatmap.advanced=TRUE,
                                            color_me=colorCustom(50,pal = "gygn"),
                                            color_group=colorCustom(50,pal = "gygn"),
                                            color.heatmap=NULL,  ###two colors at least
                                            heatmap.add.line=TRUE,
                                            heatmap.top.class=5,
                                            heatmap.sample.angle=90,
                                            heatmap.taxa.angle=0) {
  colname<-colnames(abun)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  split_otu <- lapply(
    apply(
      sapply(trt_id,function(x){grep(x,colnames(abun))}),2,
      FUN = function(x){abun[,x]}),
    function(x){x[-(which(rowSums(x)==0)),]})

  maxmod<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    wtc <- cluster_fast_greedy(g1,NA)
    V(g1)$module<-membership(wtc)
    maxm<-max(unique(V(g1)$module))
    {
      maxm<-max(unique(V(g1)$module))
      maxmod<-rbind(maxmod,maxm)
    }
  }
  maxcol<-max(maxmod)

  resr<-data.frame()
  resp<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    wtc <- igraph::cluster_fast_greedy(g1,NA)
    V(g1)$module<-igraph::membership(wtc)
    aaaa<-t(split_otu[[i]])%>%data.frame()
    MEs0 = WGCNA::moduleEigengenes(aaaa[,V(g1)$name], colors = wtc$membership)$eigengenes
    mes1 <- MEs0#[,paste("ME",order(table(wtc$membership),decreasing = TRUE)[1:max(unique(V(g1)$module))],sep = "")]
    corr_me_env <- psych::corr.test(aaaa[,1],mes1,adjust="fdr",method = "spearman")
    r<-corr_me_env$r%>%data.frame()
    p<-corr_me_env$p%>%data.frame()

    if (maxcol !=max(unique(V(g1)$module))) {
      r[,(max(unique(V(g1)$module))+1):max(maxmod)]<-"0"
      p[,(max(unique(V(g1)$module))+1):max(maxmod)]<-"0"
      colnames(r)[(max(unique(V(g1)$module))+1):max(maxmod)]<-paste("ME",(max(unique(V(g1)$module))+1):max(maxmod),sep = "")
      colnames(p)[(max(unique(V(g1)$module))+1):max(maxmod)]<-paste("ME",(max(unique(V(g1)$module))+1):max(maxmod),sep = "")
      r[,(max(unique(V(g1)$module))+1):max(maxmod)]<-as.numeric(r[,(max(unique(V(g1)$module))+1):max(maxmod)])
      p[,(max(unique(V(g1)$module))+1):max(maxmod)]<-as.numeric(p[,(max(unique(V(g1)$module))+1):max(maxmod)])
    }

    {
      resr<-rbind(resr,r)
      resp<-rbind(resp,p)
    }
  }

  rownames(resr)<-names(g)
  resr1<-resr
  resr1$group<-rownames(resr1)
  cxcr<-reshape2::melt(resr1)
  colnames(cxcr)[3]<-"r"
  resp1<-resp
  resp1$group<-rownames(resp1)
  cxcp<-reshape2::melt(resp1)
  colnames(cxcp)[3]<-"p.value"
  com<-cbind(cxcr,cxcp)
  com<-subset(com,select = c(2,1,3,6))

  if (!heatmap.advanced){
    p1<-pheatmap::pheatmap(resr,
                           gaps_row = c(2, 4),
                           scale = "row", cellwidth = 10, cellheight = 10,
                           clustering_distance_rows = "correlation")}


  if (heatmap.advanced) {
    net<-resr%>%as.matrix()
    df1 <- data.frame(group = colnames(net))
    rownames(df1) <- colnames(net)

    if (length(color_me)<nrow(df1)) {
      color_me<-c(color_me,color_me,color_me,color_me)
      color_me<-color_me[1:nrow(df1)]
    } else {
      color_me<-color_me[1:nrow(df1)]
    }

    color_me<-color_me
    names(color_me)<-unique(df1$group)
    color.use1 = color_me

    col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df1, col = list(group = color.use1),
                                                        which = "column", show_legend = TRUE, show_annotation_name = FALSE,
                                                        simple_anno_size = grid::unit(0.2, "cm"),
                                                        annotation_legend_param=list(title_gp = grid::gpar(fontsize = 8,
                                                                                                           fontfamily="serif",
                                                                                                           fontface = "plain"),
                                                                                     title_position = "leftcenter-rot",
                                                                                     border = NA,
                                                                                     legend_height = unit(20, "mm"),
                                                                                     labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                                                     grid_width = unit(2, "mm")))

    df2 <- data.frame(group = rownames(net))
    rownames(df2) <- rownames(net)

    if (length(color_group)<nrow(df2)) {
      color_group<-c(color_group,color_group,color_group,color_group)
      color_group<-color_group[1:nrow(df2)]
    } else {
      color_group<-color_group[1:nrow(df2)]
    }

    color_group<-colorCustom(nrow(df2),pal = "gygn")
    names(color_group)<-unique(df2$group)
    color.use2 = color_group

    row_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df2, col = list(group = color.use2),
                                                        which = "row", show_legend = TRUE, show_annotation_name = FALSE,
                                                        simple_anno_size = grid::unit(0.2, "cm"),
                                                        annotation_legend_param=list(title_gp = grid::gpar(fontsize = 8,
                                                                                                           fontfamily="serif",
                                                                                                           fontface = "plain"),
                                                                                     title_position = "leftcenter-rot",
                                                                                     border = NA,
                                                                                     legend_height = unit(20, "mm"),
                                                                                     labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                                                     grid_width = unit(2, "mm")))

    if (heatmap.add.line) {
      ha1 = ComplexHeatmap::rowAnnotation(
        Strength = ComplexHeatmap::anno_barplot(rowSums(abs(net)),
                                                border = FALSE, gp = grid::gpar(
                                                  fill = color.use2, col = color.use2)),
        foo = ComplexHeatmap::anno_lines(rowSums(abs(net)),
                                         gp = grid::gpar(col = "red"),
                                         add_points = TRUE,smooth = TRUE, axis = FALSE,
                                         extend = 0.2,border = FALSE,
                                         pt_gp = grid::gpar(fill = color.use2,col = color.use2),
                                         height = unit(1, "cm")),
        show_annotation_name = FALSE)
      ha2 = ComplexHeatmap::HeatmapAnnotation(
        #pt = anno_points(colSums(abs(net[1:heatmap.top.class,])), gp = gpar(fill = color.use1,col = color.use1),height = unit(1, "cm")),
        foo = ComplexHeatmap::anno_lines(colSums(abs(net[1:length(unique(com$group)),])),
                                         gp = grid::gpar(col = "red"),
                                         add_points = TRUE,smooth = TRUE, axis = FALSE,
                                         extend = 0.2,border = FALSE,
                                         pt_gp = grid::gpar(fill = color.use1,col = color.use1),
                                         height = unit(1, "cm")),
        Strength = ComplexHeatmap::anno_barplot(colSums(abs(net)),
                                                border = FALSE, gp = grid::gpar(
                                                  fill = color.use1, col = color.use1)),
        show_annotation_name = FALSE)
    } else {
      ha1 = ComplexHeatmap::rowAnnotation(
        Strength = ComplexHeatmap::anno_barplot(rowSums(abs(net)),
                                                border = FALSE, gp = grid::gpar(
                                                  fill = color.use2, col = color.use2)),
        show_annotation_name = FALSE)
      ha2 = ComplexHeatmap::HeatmapAnnotation(
        Strength = ComplexHeatmap::anno_barplot(colSums(abs(net)),
                                                border = FALSE, gp = grid::gpar(
                                                  fill = color.use1, col = color.use1)),
        show_annotation_name = FALSE)

    }


    if (is.null(color.heatmap)) {
      color.heatmap.use = grDevices::colorRampPalette((
        RColorBrewer::brewer.pal(n = 9,name = "GnBu")
      ))(100)
    } else {
      if (length(color.heatmap)<2) {
        stop("There are no enough colors to fill heatmap blocks !!!")
      } else {
        color.heatmap.use = grDevices::colorRampPalette((
          color.heatmap
        ))(100)
      }
    }

    small_mat<-resp
    small_mat[small_mat<0.05]<-"*"
    small_mat[small_mat>=0.05]<-""
    small_mat[net==0]<-""
    p1<-ComplexHeatmap::Heatmap(net,
                                col = color.heatmap.use,
                                cell_fun = function(j,i,x,y,width,height,fill) {
                                  grid::grid.text(small_mat[i,j],x,y,
                                                  gp=grid::gpar(fontsize=10,
                                                                fontfamily="serif"))
                                },
                                na_col = "white",
                                name = paste("Module-abundance",sep = ""),
                                bottom_annotation = col_annotation,
                                left_annotation = row_annotation,
                                top_annotation = ha2,
                                right_annotation = ha1,
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                row_names_side = "left",
                                row_names_rot = heatmap.taxa.angle,
                                row_names_gp = grid::gpar(fontsize = 6,
                                                          fontfamily="serif"),
                                column_names_gp = grid::gpar(fontsize = 6,
                                                             fontfamily="serif"),
                                column_names_rot = heatmap.sample.angle,
                                heatmap_legend_param = list(nrow=3,
                                                            title_gp = grid::gpar(fontsize = 8,
                                                                                  fontfamily="serif",
                                                                                  fontface = "plain"),
                                                            title_position = "leftcenter-rot",
                                                            border = NA,
                                                            legend_height = unit(20, "mm"),
                                                            labels_gp = grid::gpar(fontsize = 8,fontfamily="serif"),
                                                            grid_width = unit(2, "mm")))
  }
  message("\n","The heatmap need to be saved manually.")
  return(p1)
}

"plot_TreatModuleEigen_2heatmap" <- function(g,abun,traits) {
  colname<-colnames(abun)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  split_otu <- lapply(
    apply(
      sapply(trt_id,function(x){grep(x,colnames(abun))}),2,
      FUN = function(x){abun[,x]}),
    function(x){x[-(which(rowSums(x)==0)),]})

  maxmod<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    wtc <- cluster_fast_greedy(g1,NA)
    V(g1)$module<-membership(wtc)
    maxm<-max(unique(V(g1)$module))
    {
      maxm<-max(unique(V(g1)$module))
      maxmod<-rbind(maxmod,maxm)
    }
  }
  maxcol<-max(maxmod)


  eige<-data.frame()
  for(i in 1:length(g)){
    g1 <- g[[i]]
    E(g1)$correlation <- E(g1)$weight
    E(g1)$weight <- abs(E(g1)$weight)
    wtc <- igraph::cluster_fast_greedy(g1,NA)
    V(g1)$module<-igraph::membership(wtc)
    aaaa<-t(split_otu[[i]])%>%data.frame()
    #MEs1 = moduleEigengenes(aaaa[,V(g1)$name], colors = wtc$membership)
    MEs0 = WGCNA::moduleEigengenes(aaaa[,V(g1)$name], colors = wtc$membership)$eigengenes
    mes1 <- MEs0#[,paste("ME",order(table(wtc$membership),decreasing = TRUE)[1:max(unique(V(g1)$module))],sep = "")]
    mes<-mes1

    if (maxcol !=max(unique(V(g1)$module))) {
      mes[,(max(unique(V(g1)$module))+1):max(maxmod)]<-0
      colnames(mes)[(max(unique(V(g1)$module))+1):max(maxmod)]<-paste("ME",(max(unique(V(g1)$module))+1):max(maxmod),sep = "")
    }

    {
      eige<-rbind(eige,mes)
    }
  }

  library(ggcor)
  otu_rare<- traits
  otu_rare1<-t(otu_rare)%>%data.frame()

  if (nrow(eige)!=ncol(traits)) stop("nrow(abun) is not equal to nrow(eige)",nrow(eige),"\n")

  df03 <- fortify_cor(x = eige, y =otu_rare1)

  p<-quickcor(df03, mapping = aes(fill = r)) +
    geom_square()+
    scale_fill_gradient2(
      midpoint = 0, low = "#449436",
      mid = "white",
      high = "#615788",
      space = "Lab" ,n.breaks=6)

  return(p)
}



"calc_coord" <- function(node1s,x0=0,y0=0,r=1) {

  #frq<-table(node1s$module)%>%data.frame()
  #colnames(frq)<-c("module","")
  #node1s<-sel.data

  pcoord<-data.frame()
  for (newmod in unique(node1s$module)) {
    #newmod=1
    #newmod=2
    module_node<-node1s[which(node1s$module==newmod),]
    pnum<-nrow(module_node)
    roles<-module_node$role

    modl<-r*3
    x0=x0+modl
    y0=y0

    f=seq(0,2*pi*r*((pnum-1)/pnum),length.out=(pnum))
    y1 = y0 + r * cos(f)
    x1 = x0 + r * sin(f)
    pcoordz<-data.frame(x=x1,
                        y=y1,
                        x0=x0,
                        y0=y0,
                        r=r,
                        role=roles,
                        angle=f/pi*180,module=newmod)
    {
      pcoord<-rbind(pcoord,pcoordz)
    }

  }
  pcoord$name<-node1s$name
  return(pcoord)
}


"calc_2coord" <- function(nodexs,x0=0,y0=0,r=1) {
  #nodexs<-node1s
    module_node<-nodexs
    pnum<-nrow(nodexs)
    roles<-module_node$role

    modl<-r*3
    x0=x0+modl
    y0=y0

    f=seq(0,2*pi*r*((pnum-1)/pnum),length.out=(pnum))
    y1 = y0 + r * cos(f)
    x1 = x0 + r * sin(f)
    pcoordz<-data.frame(x=x1,
                        y=y1,
                        x0=x0,
                        y0=y0,
                        r=r,
                        role=roles,
                        angle=f/pi*180,module=nodexs$module)

    pcoordz$name<-nodexs$name
  return(pcoordz)
}
"calc_3coord" <- function(nodexs,x0=0,y0=0,r=1) {
  #nodexs<-node1s
  module_node<-nodexs
  module_keynode<-module_node
  moduel_node1<-module_node[2:nrow(module_node),]
  pnum<-nrow(moduel_node1)


  modl<-r*3
  x0=x0+modl
  y0=y0

  f=seq(0,2*pi*r*((pnum-1)/pnum),length.out=(pnum))
  y1 = y0 + r * cos(f)
  x1 = x0 + r * sin(f)
  pcoordz<-data.frame(name=moduel_node1$name,
                      group=moduel_node1$group,
                      x=x1,
                      y=y1,
                      x0=x0,
                      y0=y0,
                      r=r,
                      angle=f/pi*180,
                      role=moduel_node1$role,
                      degree=moduel_node1$degree,
                      Phylum=moduel_node1$Phylum,
                      Genus=moduel_node1$Genus,
                      module=moduel_node1$module)

  module_keynode<-subset(module_keynode,select=c(name,group,role,degree,Phylum,Genus,module))
  module_keynode<-module_keynode[1:pnum,]
  module_keynode<-cbind(module_keynode,pcoordz[,3:8])
  module_keynode<-subset(module_keynode,
                         select=c(name,group,x,y,x0,y0,r,
                                  angle,role,degree,Phylum,Genus,module))
  module_keynode$x<-module_keynode$x0
  module_keynode$y<-module_keynode$y0
  pcoordz<-rbind(module_keynode[1,],pcoordz)

  return(pcoordz)
}

"ncol_layout" <- function(ncols) {
  if (ncols==1) ncol=1
  if (ncols==2) ncol=2
  if (ncols==3) ncol=3
  if (ncols==4) ncol=2
  if (ncols>4 & ncols <=9) ncol=3
  if (ncols>9 & ncols <=16) ncol=4
  if (ncols>16 & ncols <=25) ncol=5
  if (ncols>25 & ncols <=36) ncol=6
  if (ncols>36 & ncols <=49) ncol=7
  if (ncols>49 & ncols <=64) ncol=8
  if (ncols>64 & ncols <=81) ncol=9
  if (ncols>81 & ncols <=100) ncol=10
  if (ncols>100 & ncols <=121) ncol=11

  ssxsa<-seq(from=1, to=ncols,by=(ncol))

  kl<-data.frame()
  for (sk in 1:length(ssxsa)) {
    if (sk<length(ssxsa)) nsk<-ssxsa[sk]:(ssxsa[sk+1]-1)
    if (sk==length(ssxsa)) nsk<-ssxsa[sk]:ncols
    {
      kl<-rbind(kl,nsk)
    }
  }

  newsa<-split(kl,1:nrow(kl))

  newsa<-lapply(newsa,function(x){
    x<-x%>%as.numeric()
    x<-unique(x)
  })
  return(newsa)
}

"plotSingleCytoRMTnet" <- function(nodes,node1s,edge1s,
                                   hub.position.adjust=TRUE,
                                   hub.link.display=FALSE,
                                   same.point.size=FALSE,
                                   pal="gygn",
                                   genus.coloring=TRUE,
                                   keynode.seek=FALSE,
                                   point.size=1,
                                   module.text.size=FALSE,
                                   text.size=1,
                                   curvature=0,r=1,
                                   link.pos="grey",
                                   link.neg="#6E0000",
                                   link.pos.size=0.3,
                                   link.neg.size=0.3) {

  ncols<-unique(node1s$module)%>%length()
  casca<-lapply(lapply(ncol_layout(ncols), function(x){
    nodxx<-node1s[which(node1s$module %in%  x),]
  }),calc_coord)

  p_coord<-data.frame()
  for (ttk in 1:length(casca)) {
    if (ttk==1) layout1<-casca[[1]]
    if (ttk>1) {
      ttx<-(ttk-1)*r*3
      layout1<-casca[[ttk]]
      layout1$y<-layout1$y-ttx
      layout1$y0<-layout1$y0-ttx
    }
    {
      p_coord<-rbind(p_coord,layout1)
    }
  }

  p_coord$order<-1:nrow(p_coord)
  p_coord$order<-p_coord$order%>%as.numeric()
  if (hub.position.adjust) {
    p_coord$x[which( p_coord$role!="Peripheral")]<-p_coord$x0[which( p_coord$role!="Peripheral")]
    p_coord$y[which( p_coord$role!="Peripheral")]<-p_coord$y0[which( p_coord$role!="Peripheral")]
    #p_coord$role[which( p_coord$role!="Peripheral")]<-"Peripheral"
    rolenum<-which( p_coord$role!="Peripheral")%>%length()
    if (rolenum!=0) {
      p_coordx<-p_coord[which( p_coord$role!="Peripheral"),]
      p_coord<-p_coord[-which( p_coord$role!="Peripheral"),]

      hub.role<-lapply( lapply(unique(p_coordx$module), function(x){
        grep(x,p_coordx$module)
      }), function(y) {
        p_coordx[y,]
      })

      ddxx<-data.frame()
      for (rk in 1:length(hub.role)) {
        dd<-hub.role[[rk]]

        if (nrow(dd)==1) cpcoord<-0
        if (nrow(dd)!=1) cpcoord<-seq(-0.5,0.5,length.out=nrow(dd))

        for (kk in 1:nrow(dd)) {
          dd$x[kk]<-dd$x[kk]+cpcoord[kk]
        }

        {
          ddxx<-rbind(ddxx,dd)
        }
      }
    }
    p_coordx<-ddxx
    p_coord<-rbind(p_coord,p_coordx)
    p_coord<-p_coord[order(p_coord$order),]
  }

  pcoord<-p_coord
  node1s$module<-factor(node1s$module,levels = unique(node1s$module))
  pcoord$name<-node1s$name
  pcoord<-merge(pcoord,subset(node1s,select=c(name,degree,Phylum,Genus,zi,pi,role)),by="name")
  pcoord$module<-factor(pcoord$module,levels = unique(node1s$module))

  if (!genus.coloring) {
    pcoord$Phylum<-factor(pcoord$Phylum,levels = unique(nodes$Phylum))
    pcoord<-pcoord[order(pcoord$Phylum,decreasing = FALSE),]
    color_module<-colorCustom(length(unique(pcoord$Phylum)),pal = pal)
    names(color_module)<-unique(pcoord$Phylum)
  } else {
    pcoord$Genus<-factor(pcoord$Genus,levels = unique(nodes$Genus))
    pcoord<-pcoord[order(pcoord$Genus,decreasing = FALSE),]
    color_genus<-colorCustom(length(unique(pcoord$Genus)),pal = pal)
    names(color_genus)<-unique(pcoord$Genus)
  }

  edge2s<-subset(pcoord, select=c(name,x,y))
  edge3s<-subset(pcoord, select=c(name,x,y))
  edge2s<-inner_join(edge1s,edge2s,by = c("source" = "name"))
  edge2s<-inner_join(edge2s,edge3s,by = c("target" = "name"))

  edge2sp<-edge2s[which(edge2s$cor>0),]
  edge2sn<-edge2s[which(edge2s$cor<0),]


  ddxxb<-pcoord
  ddxxb<-ddxxb[order(ddxxb$order),]
  ddxxb<-ddxxb[which( p_coord$role!="Peripheral"),]

  if (hub.link.display) {
    edge2sp<-edge2sp[which(edge2sp$source %in% ddxxb$name | edge2sp$target %in% ddxxb$name),]
    edge2sn<-edge2sn[which(edge2sn$source %in% ddxxb$name | edge2sn$target %in% ddxxb$name),]
  }

  sele.mod<-table(pcoord$module)%>%data.frame()
  sele.mod<-sele.mod[which(sele.mod$Freq>=5),]
  sele.mod1<-sele.mod
  sele.mod<-sele.mod$Var1

  pcoordxs<-pcoord[match(sele.mod,pcoord$module),]
  pcoordxs$Mtext<-paste("M",pcoordxs$module,sep = "")
  pcoordxs<-subset(pcoordxs,select=c(name,x0,y0,Mtext))
  pcoordxs<-cbind(pcoordxs,sele.mod1)
  if (!module.text.size) pcoordxs$tsize<-pcoordxs$Freq/min(pcoordxs$Freq)
  if (module.text.size) pcoordxs$tsize<-1


  all.roles<-c("Peripheral","Connector", "Module hub","Network hub")
  shapes<-c(19,17,18,15)
  names(shapes)<-all.roles

  pp<-ggplot(data = pcoord,aes(x=x,y=y))+
    geom_curve(data =edge2sp,curvature = curvature,size = link.pos.size,
               colour =link.pos,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn,curvature = curvature,size = link.neg.size,
               colour =link.neg,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )

  if (!genus.coloring){
    if (!same.point.size) pp<-pp+geom_point(data = pcoord,
                                            aes(color=Phylum,size=degree,
                                                shape=role.x))+
        geom_text(data=pcoordxs,aes(x=x0,y=y0,label=Mtext,
                                    size=text.size*pcoordxs$tsize),
                  family = "serif",fontface="bold")+
        #ylim(NA,max(pcoord$x))+
        scale_color_manual(values = color_module)+
        scale_shape_manual(values = shapes)

    if (same.point.size) pp<-pp+geom_point(data = pcoord,size=point.size,
                                           aes(color=Phylum,shape=role.x))+
        geom_text(data=pcoordxs,aes(x=x0,y=y0,label=Mtext,
                                    size=text.size*pcoordxs$tsize),
                  family = "serif",fontface="bold")+
        #ylim(NA,max(pcoord$x))+
        scale_color_manual(values = color_module)+
        scale_shape_manual(values = shapes)
  } else {
    if (!same.point.size) pp<-pp+geom_point(data = pcoord,
                                            aes(color=Genus,size=degree,shape=role.x))+
        geom_text(data=pcoordxs,aes(x=x0,y=y0,label=Mtext,
                                    size=text.size*pcoordxs$tsize),
                  family = "serif",fontface="bold")+
        #ylim(NA,max(pcoord$x))+
        scale_color_manual(values = color_genus)+
        scale_shape_manual(values = shapes)

    if (same.point.size) pp<-pp+geom_point(data = pcoord,size=point.size,
                                           aes(color=Genus,shape=role.x))+
        geom_text(data=pcoordxs,aes(x=x0,y=y0,label=Mtext,
                                    size=text.size*pcoordxs$tsize),
                  family = "serif",fontface="bold")+
        #ylim(NA,max(pcoord$x))+
        scale_color_manual(values = color_genus)+
        scale_shape_manual(values = shapes)
  }

  pp<-pp+theme_void()+
    theme(
      panel.border = element_rect(linetype = "dashed", fill = NA),
      legend.position = "none",
      aspect.ratio = length(casca)/length(ncol_layout(ncols)[[1]]))

  return(pp)
}

"plotMutiCytoRMTnet"<- function(rmnet_cg,
                                hub.position.adjust=TRUE,
                                hub.link.display=FALSE,
                                keynode.seek=FALSE,
                                same.point.size=FALSE,
                                pal="gygn",
                                genus.coloring=TRUE,
                                layout_ncols=2,
                                point.size=3,
                                module.text.size=TRUE,
                                text.size=8,
                                curvature=0,
                                link.pos="grey",link.neg="#6E0000",
                                link.pos.size=0.2,
                                link.neg.size=0.2,
                                export_path='cs2/microbial network analysis/topo_esti_compare') {

  export_path<-paste(export_path,"/microbial network analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  dir.create(paste(export_path,'/cyto_style',sep = ""), recursive = TRUE)
  nodes<-rmnet_cg$node_table
  edges<-rmnet_cg$edge_table

  pp<-list()
  for (t in unique(nodes$group)) {
    netname<-rmnet_cg$netname[t]
    node1s<-nodes[which(nodes$group==t),]
    node1s<-node1s[order(node1s$module,decreasing = FALSE),]
    edge1s<-edges[which(edges$group==t),]

    p<-plotSingleCytoRMTnet(nodes,node1s,edge1s,
                            hub.position.adjust=hub.position.adjust,
                            hub.link.display=hub.link.display,
                            same.point.size=same.point.size,
                            pal=pal,
                            genus.coloring=genus.coloring,
                            keynode.seek=keynode.seek,
                            point.size=point.size,
                            module.text.size=module.text.size,
                            text.size=text.size,
                            curvature=curvature,
                            link.pos=link.pos,
                            link.neg=link.neg,
                            link.pos.size=link.pos.size,
                            link.neg.size=link.neg.size)
    ggsave(paste(export_path,'/cyto_style/',"Cytoscape_style (",netname,').pdf', sep = ''),
           p)
    p<-p+guides(size="none")+
      theme(legend.position = "right",
            text = element_text(family = "serif") )
    pp[[t]]<-p
  }

  sp<-patchwork::wrap_plots(pp,ncol = layout_ncols,guides="collect")
  ggsave(paste(export_path,'/cyto_style/',"Cytoscape_style (combinations).pdf", sep = ''),
         sp)
  return(sp)
}



"plot_RMTnetInter_Strength" <- function(rmnet_cg,
                                        group.pos.coffe=0.1825,
                                        prop.coffe=1.05,
                                        bar.alpha=0.5,
                                        xlabsname=xlabname,
                                        add.line=FALSE,
                                        linetype=1,
                                        bar.text.color="orange",
                                        color_cor=c('grey',"#6E0000"),
                                        color_bar=c('grey',"#6E0000"),
                                        export_path='cs2/microbial network analysis/topo_esti_compare') {
  export_path<-paste(export_path,"/microbial network analysis/topo_esti_compare",sep = "")
  dir.create(export_path, recursive = TRUE)
  nodes<-rmnet_cg$node_table
  edges<-rmnet_cg$edge_table

  sasxc<-lapply(unique(nodes$group), function(x){
    nodesx<-nodes[which(nodes$group==x),]
    edgesx<-edges[which(edges$group==x),]

    edge2s<-inner_join(edgesx,nodesx,by = c("source" = "name"))
    edge3s<-inner_join(edge2s,nodesx,by = c("target" = "name"))
  })

  datain<-data.frame(
    prop=c(lapply(sasxc, function(y){
      nrow(y[which(y$cor>0),])/nrow(y)*100
    })%>%as.numeric(),lapply(sasxc, function(y){
      nrow(y[which(y$cor<0),])/nrow(y)*100
    })%>%as.numeric()),
    num=c(lapply(sasxc, function(y){
      nrow(y[which(y$cor>0),])
    })%>%as.numeric(),lapply(sasxc, function(y){
      nrow(y[which(y$cor<0),])
    })%>%as.numeric()),
    group=1:length(sasxc),
    variable=c(rep("Positive",length(sasxc)),rep("Negative",length(sasxc))))

  datainx<-c(datain$prop,datain$num)
  datainx%>%max()

  datain$variable<-factor(datain$variable,levels = unique(datain$variable))

  color_cor=color_cor
  color_bar=color_bar
  names(color_cor)<-unique(datain$variable)
  names(color_bar)<-unique(datain$variable)

  xlabname<-unique(datain$group)
  names(xlabname)<-rmnet_cg$netname
  datain$revp<-c(datain$prop[(nrow(datain)/2+1):nrow(datain)],datain$prop[1:(nrow(datain)/2)])

  datain$ynum<-datain$num+10
  datain$gnum<-c(datain$group[1:(nrow(datain)/2)]-group.pos.coffe,datain$group[(nrow(datain)/2+1):nrow(datain)]+group.pos.coffe)
  datain$prop1<-c(datain$prop[1:(nrow(datain)/2)]*prop.coffe,datain$prop[(nrow(datain)/2+1):nrow(datain)]*(2-prop.coffe))



  p1<-ggplot(datain,aes(group,num,fill = variable))+
    geom_bar(aes(fill=variable),stat = "identity",alpha=bar.alpha,show.legend=FALSE,
             position = "dodge",width = 0.75,color = NA)+
    geom_text(aes(x=datain$gnum,y=ynum,label=num),color=bar.text.color,family="serif")+
    scale_fill_manual(values = color_bar)+
    ylab("Interaction number")



  p1<-p1+ggnewscale::new_scale_color()+
    scale_y_continuous(expand = c(0,0),limits = c(NA,datainx%>%max()*1.05),
                       sec.axis = sec_axis(~./(datainx%>%max()/100),
                                           name = 'Interaction proportion (%)',
                                           breaks = seq(-10,100,20)))+
    geom_line(aes(x= group,y=prop*(datainx%>%max()/100)),
              linetype=linetype,cex=1)+
    geom_point(aes(x= group,color=variable,y=prop*(datainx%>%max()/100)),
               stroke=3.5,shape=21,fill="black",
               size=2,show.legend=TRUE)+
    geom_text(size=3,family="serif",
              aes(x=group,y=prop1*(datainx%>%max()/100),
                  label=paste(round(prop,2),"%",sep = "")))+
    scale_color_manual(values = color_cor)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA,colour = "black",size = 1),
          text = element_text(family = "serif"),
          plot.margin = unit(c(1, 1, 4, 1), "lines"),
          axis.ticks.length.y = unit(0.2,"lines"),
          axis.ticks.y  = element_line(size=1),
          axis.line  = element_line(size=0.25,colour = "black"),
          axis.title.y=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.title.x = element_blank(),
          #legend.title=element_blank(),
          legend.background = element_blank(),
          legend.text=element_text(size=8,face = "bold",colour = "black",
                                   family = "serif",
                                   margin = margin(r = 20)),
          legend.position = c(0.95,0.95),
          legend.key = element_blank(),
          aspect.ratio = 1)

  if (add.line) for (i in 1:(nrow(datain)/2)) {
    p1 <- p1 + annotate('segment', linetype=3,
                        x = i, xend = i, y = datain$prop[i]*(datainx%>%max()/100),
                        yend = datain$revp[i]*(datainx%>%max()/100)
    )
  }

  p1<-p1+scale_x_discrete(labels=xlabsname,
                          position = "bottom",
                          limits=as.character(rmnet_cg$netname))
  ggsave(paste(export_path,'/',"Interaction strength.pdf", sep = ''),
         p1)

  return(p1)
}




"plot_RMTnet_ModuleLink" <- function(rmnet_cg,
                                     group.pos.coffe=0.1825,
                                     prop.coffe=1.05,
                                     bar.alpha=0.5,
                                     xlabsname=xlabname,
                                     add.line=FALSE,
                                     linetype=1,
                                     bar.text.color="orange",
                                     color_cor=c('grey',"#6E0000"),
                                     color_bar=c('grey',"#6E0000"),
                                     export_path='cs2/microbial network analysis/topo_esti_compare') {
  export_path<-paste(export_path,"/microbial network analysis/topo_esti_compare",sep = "")
  dir.create(export_path, recursive = TRUE)
  nodes<-rmnet_cg$node_table
  edges<-rmnet_cg$edge_table

  sasxc<-lapply(unique(nodes$group), function(x){
    nodesx<-nodes[which(nodes$group==x),]
    edgesx<-edges[which(edges$group==x),]

    edge2s<-inner_join(edgesx,nodesx,by = c("source" = "name"))
    edge3s<-inner_join(edge2s,nodesx,by = c("target" = "name"))
  })

  newdat<-data.frame()
  for (tk in 1:length(sasxc)) {
    sas<-sasxc[[tk]]
    sas<-subset(sas,select=c(1:3,module.x,module.y,role.x,role.y))
    sas$group<-tk
    {
      newdat<-rbind(newdat,sas)
    }
  }

  withmodule<-lapply(unique(newdat$group),function(x){
    nedata<-newdat[which(newdat$group==x),]
    withmodule<-which(nedata$module.x==nedata$module.y)%>%length()
  })%>%as.numeric()


  amongmodule<-lapply(unique(newdat$group),function(x){
    nedata<-newdat[which(newdat$group==x),]
    amongmodule<-which(nedata$module.x!=nedata$module.y)%>%length()
  })%>%as.numeric()

  datandata<-data.frame(withmodule=withmodule,
                        amongmodule=amongmodule,
                        group=1:length(sasxc))

  datandata$sum1<-datandata$withmodule/(datandata$withmodule+datandata$amongmodule)*100
  datandata$sum2<-datandata$amongmodule/(datandata$withmodule+datandata$amongmodule)*100

  newdatx<-data.frame(
    prop=c(datandata$sum1,datandata$sum2),
    num=c(datandata$withmodule,datandata$amongmodule),
    group=c(datandata$group,datandata$group),
    variable=c(rep("Within-module",length(sasxc)),rep("Among-module",length(sasxc)))
  )

  datain<-newdatx

  datainx<-c(datain$prop,datain$num)
  datainx%>%max()

  datain$variable<-factor(datain$variable,levels = unique(datain$variable))

  color_cor=color_cor
  color_bar=color_bar
  names(color_cor)<-unique(datain$variable)
  names(color_bar)<-unique(datain$variable)

  xlabname<-unique(datain$group)
  names(xlabname)<-rmnet_cg$netname
  datain$revp<-c(datain$prop[(nrow(datain)/2+1):nrow(datain)],datain$prop[1:(nrow(datain)/2)])

  datain$ynum<-datain$num+10
  datain$gnum<-c(datain$group[1:(nrow(datain)/2)]-group.pos.coffe,datain$group[(nrow(datain)/2+1):nrow(datain)]+group.pos.coffe)
  datain$prop1<-c(datain$prop[1:(nrow(datain)/2)]*prop.coffe,datain$prop[(nrow(datain)/2+1):nrow(datain)]*(2-prop.coffe))



  p1<-ggplot(datain,aes(group,num,fill = variable))+
    geom_bar(aes(fill=variable),stat = "identity",alpha=bar.alpha,show.legend=FALSE,
             position = "dodge",width = 0.75,color = NA)+
    geom_text(aes(x=datain$gnum,y=ynum,label=num),color=bar.text.color,family="serif")+
    scale_fill_manual(values = color_bar)+
    ylab("Module link number")



  p1<-p1+ggnewscale::new_scale_color()+
    scale_y_continuous(expand = c(0,0),limits = c(NA,datainx%>%max()*1.05),
                       sec.axis = sec_axis(~./(datainx%>%max()/100),
                                           name = 'Module link proportion (%)',
                                           breaks = seq(-10,100,20)))+
    geom_line(aes(x= group,y=prop*(datainx%>%max()/100)),
              linetype=linetype,cex=1)+
    geom_point(aes(x= group,color=variable,y=prop*(datainx%>%max()/100)),
               stroke=3.5,shape=21,fill="black",
               size=2,show.legend=TRUE)+
    geom_text(size=3,family="serif",
              aes(x=group,y=prop1*(datainx%>%max()/100),
                  label=paste(round(prop,2),"%",sep = "")))+
    scale_color_manual(values = color_cor)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA,colour = "black",size = 1),
          text = element_text(family = "serif"),
          plot.margin = unit(c(1, 1, 4, 1), "lines"),
          axis.ticks.length.y = unit(0.2,"lines"),
          axis.ticks.y  = element_line(size=1),
          axis.line  = element_line(size=0.25,colour = "black"),
          axis.title.y=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.title.x = element_blank(),
          #legend.title=element_blank(),
          legend.background = element_blank(),
          legend.text=element_text(size=8,face = "bold",colour = "black",
                                   family = "serif",
                                   margin = margin(r = 20)),
          legend.position = c(0.95,0.95),
          legend.key = element_blank(),
          aspect.ratio = 1)

  if (add.line) for (i in 1:(nrow(datain)/2)) {
    p1 <- p1 + annotate('segment', linetype=3,
                        x = i, xend = i, y = datain$prop[i]*(datainx%>%max()/100),
                        yend = datain$revp[i]*(datainx%>%max()/100)
    )
  }

  p1<-p1+scale_x_discrete(labels=xlabsname,
                          position = "bottom",
                          limits=as.character(rmnet_cg$netname))
  ggsave(paste(export_path,'/',"Within-among module link.pdf", sep = ''),
         p1)

  return(p1)
}



"plot_RMTnet_keynodeLink" <- function(node1s,edge1s,
                                      keynode.link.display=TRUE,
                                      curvature=0.2,
                                      link.pos="grey",
                                      link.neg="#6E0000",
                                      link.pos.size=0.3,
                                      link.neg.size=0.3,
                                      genus.coloring=TRUE,
                                      pal="gygn",
                                      same.point.size=FALSE,
                                      point.size=1,
                                      text.size=5,
                                      show.message=TRUE) {
  options(warn = -1)
  node2s<-node1s[which(node1s$role!="Peripheral"),]

  saah<-lapply(lapply(unique(node2s$name), function(y){
    node2s[which(node2s$name==y),]
  }),function(st) {
    st<-subset(st,select=c(name,role))
    edge2s<-inner_join(edge1s,st,by = c("source" = "name"))
    edge3s<-inner_join(edge1s,st,by = c("target" = "name"))

    sel.node<-c(st$name,unique(edge2s$target),unique(edge3s$source))%>%unique()
    sel.data<-node1s[match(sel.node,node1s$name),]
    p_coord<-calc_3coord(sel.data)
  })

  for (ttk in 1:length(saah)) {
    ttx<-(ttk-1)*3
    saah[[ttk]]$x<-saah[[ttk]]$x+ttx
    saah[[ttk]]$x0<-saah[[ttk]]$x0+ttx
  }

  sp<-list()
  for (tk in 1:length(saah)) {
    #sel.keynode<-unique(node2s$name)[tk]
    pcoord<-saah[[tk]]
    sel.keynode<-pcoord$name[1]
    #pcoord1<-merge(sagdsv[[tk]],subset(node1s,select=c(name,degree,Phylum,Genus,zi,pi,role)),by="name")
    #pcoord<-calc_3coord(pcoord1)
    #pcoord<-merge(pcoord,subset(pcoord1,select=c(name,degree,Phylum,Genus,zi,pi)))
    #pcoord$y<-pcoord[which(pcoord$name==sel.keynode),]$y0

    edge2s<-subset(pcoord, select=c(name,x,y))
    edge3s<-subset(pcoord, select=c(name,x,y))
    edge2xs<-inner_join(edge1s,edge2s,by = c("source" = "name"))
    edge2xss<-inner_join(edge2xs,edge3s,by = c("target" = "name"))


    ###only display links related to keynode
    if (keynode.link.display) {
    sel.keynode.link.row<-c(which(edge2xss$source %in% sel.keynode),which(edge2xss$target %in% sel.keynode))
    edge2xss<-edge2xss[sel.keynode.link.row,]
    }

    edge2sp<-edge2xss[which(edge2xss$cor>0),]
    edge2sn<-edge2xss[which(edge2xss$cor<0),]

    pcoord$Mtext<-pcoord$name
    pcoord$tsize<-1

    all.roles<-c("Peripheral","Connector", "Module hub","Network hub")
    shapes<-c(19,17,18,15)
    names(shapes)<-all.roles


    if (!genus.coloring) {
      pcoord$Phylum<-factor(pcoord$Phylum,levels = unique(nodes$Phylum))
      color_module<-colorCustom(length(unique(nodes$Phylum)),pal = pal)
      names(color_module)<-unique(pcoord$Phylum)
    } else {
      pcoord$Genus<-factor(pcoord$Genus,levels = unique(nodes$Genus))
      color_genus<-colorCustom(length(unique(nodes$Genus)),pal = pal)
      names(color_genus)<-unique(pcoord$Genus)
    }

    pp<-ggplot(data = pcoord,aes(x=x,y=y))+
      geom_curve(data =edge2sp,curvature = curvature,size = link.pos.size,
                 colour =link.pos,
                 aes(x = x.x, y = y.x,
                     xend = x.y, yend = y.y) )+
      geom_curve(data =edge2sn,curvature = curvature,size = link.neg.size,
                 colour =link.neg,
                 aes(x = x.x, y = y.x,
                     xend = x.y, yend = y.y) )

    if (!genus.coloring){
      if (!same.point.size) pp<-pp+geom_point(data = pcoord,
                                              aes(color=Phylum,
                                                  size=degree,
                                                  shape=role))+
          geom_text(data=pcoord,aes(x=x,y=y,label=Mtext),size=text.size,
                    family = "serif",fontface="bold")+
          #ylim(NA,max(pcoord$x))+
          scale_color_manual(values = color_module)+
          scale_shape_manual(values = shapes)

      if (same.point.size) pp<-pp+geom_point(data = pcoord,size=point.size,
                                             aes(color=Phylum,shape=role))+
          geom_text(data=pcoord,aes(x=x,y=y,label=Mtext),size=text.size,
                    family = "serif",fontface="bold")+
          #ylim(NA,max(pcoord$x))+
          scale_color_manual(values = color_module)+
          scale_shape_manual(values = shapes)
    } else {
      if (!same.point.size) pp<-pp+geom_point(data = pcoord,
                                              aes(color=Genus,size=degree,
                                                  shape=role))+
          geom_text(data=pcoord,aes(x=x,y=y,label=Mtext),size=text.size,
                    family = "serif",fontface="bold")+
          #ylim(NA,max(pcoord$x))+
          scale_color_manual(values = color_genus)+
          scale_shape_manual(values = shapes)

      if (same.point.size) pp<-pp+geom_point(data = pcoord,size=point.size,
                                             aes(color=Genus,shape=role))+
          geom_text(data=pcoord,aes(x=x,y=y,label=Mtext),size=text.size,
                    family = "serif",fontface="bold")+
          #ylim(NA,max(pcoord$x))+
          scale_color_manual(values = color_genus)+
          scale_shape_manual(values = shapes)
    }


    pp<-pp+theme_void()+
      theme(
        plot.margin = margin(t = 0.25,
                             r = 0.25,
                             b = 0.25,
                             l = 0.25,
                             unit = "cm"),
        text = element_text(family = "serif",size=text.size),
        #panel.border = element_rect(linetype = "dashed", fill = NA),
        legend.position = "none",
        aspect.ratio = 1)

    sp[[tk]]<-pp
  }
  if (length(saah)<=5) nrows=1
  if (length(saah)>5 & length(saah)<=10) nrows=2
  if (length(saah)>10 & length(saah)<=15) nrows=3
  if (length(saah)>15) nrows=5
  xp<-patchwork::wrap_plots(sp,nrow = nrows)

  if (show.message) message("\nThe plot needs to be saved mannually!!!")
  return(xp)
}



"plot_MutiRMTnet_keynodeLink" <- function(rmnet_cg,
                                          nrows=2,
                                          keynode.link.display=TRUE,
                                          curvature=0.2,
                                          link.pos="grey",
                                          link.neg="#6E0000",
                                          link.pos.size=0.3,
                                          link.neg.size=0.3,
                                          genus.coloring=TRUE,
                                          pal="gygn",
                                          same.point.size=FALSE,
                                          point.size=1,
                                          text.size=2,
                                          export_path='cs2/microbial network analysis') {
  export_path<-paste(export_path,"/microbial network analysis/Keynode",sep = "")
  dir.create(export_path, recursive = TRUE)
  nodes<-rmnet_cg$node_table
  edges<-rmnet_cg$edge_table

  mp<-list()
  for (t in 1:length(rmnet_cg$netname)) {
    netname<-rmnet_cg$netname[t]
    node1s<-nodes[which(nodes$group==t),]
    node1s<-node1s[order(node1s$module,decreasing = FALSE),]
    edge1s<-edges[which(edges$group==t),]

    sp<-plot_RMTnet_keynodeLink(node1s,edge1s,
                                keynode.link.display=keynode.link.display,
                                curvature=curvature,
                                link.pos=link.pos,
                                link.neg=link.neg,
                                link.pos.size=link.pos.size,
                                link.neg.size=link.neg.size,
                                genus.coloring=genus.coloring,
                                pal=pal,
                                same.point.size=same.point.size,
                                point.size=point.size,
                                text.size=text.size,
                                show.message=FALSE)
    ggsave(paste(export_path,'/(',netname,")","Single group keynode link.pdf", sep = ''),
           sp)

    mp[[t]]<-sp
  }

  mmp<-patchwork::wrap_plots(mp,nrow = nrows)
  ggsave(paste(export_path,'/',"Muti group keynode link.pdf", sep = ''),
         mmp)
  return(mmp)
}


"calc_nodecoord" <- function(nodes,node1s,edge1s,
                             hub.position.adjust=FALSE,
                             r=1,min.node.num=5) {

  ncols<-unique(node1s$module)%>%length()
  casca<-lapply(lapply(ncol_layout(ncols), function(x){
    nodxx<-node1s[which(node1s$module %in%  x),]
  }),calc_coord)

  maxwidth<-length(ncol_layout(ncols)[[1]])

  p_coord<-data.frame()
  for (ttk in 1:length(casca)) {
    if (ttk==1) layout1<-casca[[1]]
    if (ttk>1) {
      ttx<-(ttk-1)*r*3
      layout1<-casca[[ttk]]
      layout1$y<-layout1$y-ttx
      layout1$y0<-layout1$y0-ttx
    }
    {
      p_coord<-rbind(p_coord,layout1)
    }
  }

  p_coord$order<-1:nrow(p_coord)
  p_coord$order<-p_coord$order%>%as.numeric()
  if (hub.position.adjust) {
    p_coord$x[which( p_coord$role!="Peripheral")]<-p_coord$x0[which( p_coord$role!="Peripheral")]
    p_coord$y[which( p_coord$role!="Peripheral")]<-p_coord$y0[which( p_coord$role!="Peripheral")]
    #p_coord$role[which( p_coord$role!="Peripheral")]<-"Peripheral"
    rolenum<-which( p_coord$role!="Peripheral")%>%length()
    if (rolenum!=0) {
      p_coordx<-p_coord[which( p_coord$role!="Peripheral"),]
      p_coord<-p_coord[-which( p_coord$role!="Peripheral"),]

      hub.role<-lapply( lapply(unique(p_coordx$module), function(x){
        grep(x,p_coordx$module)
      }), function(y) {
        p_coordx[y,]
      })

      ddxx<-data.frame()
      for (rk in 1:length(hub.role)) {
        dd<-hub.role[[rk]]

        if (nrow(dd)==1) cpcoord<-0
        if (nrow(dd)!=1) cpcoord<-seq(-0.5,0.5,length.out=nrow(dd))

        for (kk in 1:nrow(dd)) {
          dd$x[kk]<-dd$x[kk]+cpcoord[kk]
        }

        {
          ddxx<-rbind(ddxx,dd)
        }
      }
    }
    p_coordx<-ddxx
    p_coord<-rbind(p_coord,p_coordx)
    p_coord<-p_coord[order(p_coord$order),]
  }

  pcoord<-p_coord
  node1s$module<-factor(node1s$module,levels = unique(node1s$module))
  pcoord$name<-node1s$name
  pcoord<-merge(pcoord,subset(node1s,select=c(name,degree,Phylum,Genus,zi,pi,role)),by="name")
  pcoord$module<-factor(pcoord$module,levels = unique(node1s$module))

  edge2s<-subset(pcoord, select=c(name,x,y))
  edge3s<-subset(pcoord, select=c(name,x,y))
  edge2s<-inner_join(edge1s,edge2s,by = c("source" = "name"))
  edge2s<-inner_join(edge2s,edge3s,by = c("target" = "name"))

  edge2sp<-edge2s[which(edge2s$cor>0),]
  edge2sn<-edge2s[which(edge2s$cor<0),]

  sele.mod<-table(pcoord$module)%>%data.frame()
  sele.mod<-sele.mod[which(sele.mod$Freq>min.node.num),]
  sele.mod1<-sele.mod
  sele.mod<-sele.mod$Var1

  pcoordxs<-pcoord[match(sele.mod,pcoord$module),]
  pcoordxs$Mtext<-paste("M",pcoordxs$module,sep = "")
  pcoordxs<-subset(pcoordxs,select=c(name,x0,y0,Mtext))
  pcoordxs<-cbind(pcoordxs,sele.mod1)
  pcoordxs$tsize<-1


  return(list(node_pro=pcoord,edge_pos=edge2sp,edge_neg=edge2sn,
              module.text=pcoordxs,maxwidth=maxwidth,
              ncols=length(ncol_layout(ncols))))
}

"calcRMTnetNodeCoord"<- function(rmnet_cg,min.node.num=5) {
  #dir.create(export_path, recursive = TRUE)
  nodes<-rmnet_cg$node_table
  edges<-rmnet_cg$edge_table

  xnode<-data.frame()
  xpedge<-data.frame()
  xnedge<-data.frame()
  xtext<-data.frame()
  colu<-0
  xmaxwd<-c()
  t=1
  for (t in unique(nodes$group)) {
    netname<-rmnet_cg$netname[t]
    node1s<-nodes[which(nodes$group==t),]
    node1s<-node1s[order(node1s$module,decreasing = FALSE),]
    edge1s<-edges[which(edges$group==t),]

    xp<-calc_nodecoord(nodes,node1s,edge1s,min.node.num=min.node.num)
    xp$node_pro$group<-netname

    xp$edge_pos$group<-netname
    xp$edge_neg$group<-netname

    xp$edge_pos$modsource<-xp$node_pro$module[match(xp$edge_pos$source,xp$node_pro$name)]
    xp$edge_pos$modtarget<-xp$node_pro$module[match(xp$edge_pos$target,xp$node_pro$name)]

    xp$edge_neg$modsource<-xp$node_pro$module[match(xp$edge_neg$source,xp$node_pro$name)]
    xp$edge_neg$modtarget<-xp$node_pro$module[match(xp$edge_neg$target,xp$node_pro$name)]


    xp$edge_pos$cur<-ifelse(xp$edge_pos$modsource==xp$edge_pos$modtarget,"show","none")
    xp$edge_neg$cur<-ifelse(xp$edge_neg$modsource==xp$edge_neg$modtarget,"show","none")

    xp$module.text$group<-netname
    wd<-xp$maxwidth
    ncols<-xp$ncols

    {
      colu<-colu+ncols
      xcolu<-colu-ncols
      xmaxwd<-c(xmaxwd,wd)
    }

    rad.adj<-3*1*xcolu
    "node.y.adj"<-function(nnode,rad.adj) {
      nnode$y <- nnode$y-rad.adj
      nnode$y0 <- nnode$y0-rad.adj
      return(nnode)
    }
    "edge.y.adj"<-function(eedge,rad.adj) {
      eedge$y.x <- eedge$y.x-rad.adj
      eedge$y.y <- eedge$y.y-rad.adj
      return(eedge)
    }
    "text.y.adj"<-function(nnode,rad.adj) {
      nnode$y0 <- nnode$y0-rad.adj
      return(nnode)
    }

    nnode<-xp$node_pro<-node.y.adj(xp$node_pro,rad.adj)
    mtext<-xp$module.text<-text.y.adj(xp$module.text,rad.adj)
    pedge<-xp$edge_pos<-edge.y.adj(xp$edge_pos,rad.adj)
    nedge<-xp$edge_neg<-edge.y.adj(xp$edge_neg,rad.adj)

    {
      xnode<-rbind(xnode,nnode)
      xtext<-rbind(xtext,mtext)
      xpedge<-rbind(xpedge,pedge)
      xnedge<-rbind(xnedge,nedge)
    }

  }

  sp<-list(node_pro=xnode,edge_pos=xpedge,edge_neg=xnedge,
           mtext=xtext,tncol=colu,maxwidth=max(xmaxwd))

  return(sp)
}


"calcRMTnetModPreserve"<- function(rmnet_cg,
                                   min.node.num=3,
                                   export_path="cs2/microbial network analysis") {

  export_path<-paste(export_path,"/microbial network analysis/module preserve",sep = "")
  dir.create(export_path, recursive = TRUE)
  node_table2<-rmnet_cg$node_table
  edge_table2<-rmnet_cg$edge_table
  groupname<-rmnet_cg$netname
  for (tk in unique(rmnet_cg$node_table$group)) {
    node_table2$group[which(node_table2$group==tk)]<-groupname[tk]
    edge_table2$group[which(edge_table2$group==tk)]<-groupname[tk]
  }
  node_table2$Group<-node_table2$group
  node_table2$group<-paste(node_table2$group,"_M",node_table2$module,sep = "")
  node_table2<-subset(node_table2, select=c(name,group,Group))

  ddaxt.dat<-model_preserve(node_table2 = node_table2,n = min.node.num,
                        export_path=export_path)
  ddaxt<-ddaxt.dat$data
  if (nrow(ddaxt)==0) stop("No module pairs were detected.")

  el <- matrix( c(ddaxt$module1,ddaxt$module2), nc = 2, byrow = FALSE)
  gg<-graph_from_edgelist(el,directed = FALSE)
  #plot(gg)

  V(gg)$modularity <- membership(cluster_fast_greedy(gg))
  ###共五个簇，需要5颜色
  modulenum<-length(cluster_fast_greedy(gg))
  node_table = as.data.frame(vertex.attributes(gg))
  colnames(node_table)[2]<-"clusters"
  ddaxt<-inner_join(ddaxt,node_table,by=c("module1"="name"))

  ###remove connected module
  wqc<-ddaxt
  wqc$module1<-wqc$module2
  wqc<-rbind(ddaxt,wqc)

  rm.n<-c()
  for (sel.mod in unique(wqc$module1)) {
    sel.dat<-wqc[which(wqc$module1==sel.mod),]
    sel.dat$Both<-sel.dat$Both%>%as.numeric()
    max.row<-which(sel.dat$Both==max(sel.dat$Both))

    sel.data<-sel.dat[which(sel.dat$clusters[setdiff(sel.dat$clusters,sel.dat[max.row,]$clusters)[1]]!=
                              sel.dat[max.row,]$clusters),]
    rm.num<-sel.data%>%rownames()%>%as.numeric()
    if (length(rm.num[rm.num<=(nrow(wqc)/2)])!=0) rm.num[rm.num<=(nrow(wqc)/2)]<-rm.num[rm.num<=(nrow(wqc)/2)]
    if (length(rm.num[rm.num>(nrow(wqc)/2)])!=0) rm.num[rm.num>(nrow(wqc)/2)]<-rm.num[rm.num>(nrow(wqc)/2)]-(nrow(wqc)/2)
    rm.n<-c(rm.n,rm.num)
  }
if (!is.integer(rm.n)) ddaxt<-ddaxt else ddaxt<-ddaxt[-rm.n,]

  rownames(ddaxt)<-1:nrow(ddaxt)

  cat("\n","Totally ",nrow(ddaxt)," module pairs were preserved after removing the inter-clusters shared module, accounted for ",
      paste(round(nrow(ddaxt)/ddaxt.dat$pairs_sum*100,2),"%.",sep = ""),sep = "")
  cat("\n","All preserved module pairs have been exported to"," '",export_path,"'.",sep = "")

  return(ddaxt)

}


"plotRMTnetModPreserve" <- function(g,taxon,
                                    specified.taxa.order=FALSE,
                                    taxa.order=levels(unique(microchatcomobj$plot_data$Taxo)),
                                    pdf.sel=FALSE,
                                    xlabname,
                                    bar.height=NULL, ### or 随便一个数字
                                    display.panel.border=FALSE,
                                    display.interval.line=FALSE,
                                    display.interval.line.style=c("grey","dashed",1),
                                    pal_taxa ="gygn",
                                    pal_taxa_alpha=0.5,
                                    pal_cluster = "set2",
                                    color_bar.text="black",
                                    min.node.num=3,
                                    show.modulepairs.node=FALSE,
                                    curvature=0,
                                    link.pos="grey",
                                    link.neg="#6E0000",
                                    link.left.alpha=0.3,
                                    link.right.alpha=0.2,
                                    link.pos.size=0.2,
                                    link.neg.size=0.2,
                                    link.clu.size=0.5,
                                    point.size=1,
                                    point.clu.size=20,
                                    export_path="cs2/microbial network analysis") {
  require(igraph)
  export_pathx<-export_path
  export_path<-paste(export_path,"/microbial network analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  rmnet_cg<-rmnet_cytogephi_export(g,taxon,
                                   file.save=FALSE,
                                   export_path =export_path)

  xa<-get.MutiPie(g,taxon)
  xa$mod<-paste0(xa$group,xa$module)

  newtt<-calcRMTnetNodeCoord(rmnet_cg,min.node.num=min.node.num-1)
  pcoord<-newtt$node_pro
  edge2sp<-newtt$edge_pos
  edge2sn<-newtt$edge_neg
  pcoordxs<-newtt$mtext
  totalncol<-newtt$tncol
  maxwidth<-newtt$maxwidth
  nodes<-newtt$node_pro

  pcoord$Phylum<-factor(pcoord$Phylum,levels = unique(nodes$Phylum))
  pcoord<-pcoord[order(pcoord$Phylum,decreasing = FALSE),]
  color_module<-colorCustom(length(unique(pcoord$Phylum)),pal=pal_taxa) %>% adjustcolor(alpha.f = pal_taxa_alpha)
  if (!specified.taxa.order)  names(color_module)<-unique(pcoord$Phylum) else names(color_module)<-taxa.order

  ###merge scatterpie coord
  pcoord.pie<-pcoord
  newpcoord.pie<-subset(pcoord.pie, select=c(x0,y0,r,module,degree,group))
  newpcoord.pie$cluster<-paste(newpcoord.pie$group,newpcoord.pie$module,sep = "")
  newpcoord.pie<-newpcoord.pie[!duplicated(newpcoord.pie$cluster),]

  xxa<-inner_join(xa,newpcoord.pie,by=c("mod"="cluster"))
  xxa$value<-xxa$Freq
  xxa$type<-xxa$Var1
  pc.sel<-pcoordxs
  pc.sel$clu<-paste0(pc.sel$group,pc.sel$Var1)
  mods.sel<-unique(pc.sel$clu)
  xxa<-xxa[which(xxa$mod %in% mods.sel),]

  all.roles<-c("Peripheral","Connector", "Module hub","Network hub")
  shapes<-c(19,17,18,15)
  names(shapes)<-all.roles
  shapesx<-c(1,2,5,0)
  names(shapesx)<-all.roles

  ####module pairs data preparation
  newpcoord<-subset(pcoord, select=c(x0,y0,r,module,degree,group))
  newpcoord$cluster<-paste(newpcoord$group,newpcoord$module,sep = "_M")
  newpcoord<-newpcoord[!duplicated(newpcoord$cluster),]

  ###module pairs
  ddaxt<-calcRMTnetModPreserve(rmnet_cg,
                               min.node.num=min.node.num,
                               export_path=export_pathx)


  sel.modu<-unique(c(ddaxt$module1,ddaxt$module2))
  cat("\nTaken above, totally",nrow(newpcoord),"single modules were tested.\n")
  newpcoord<-newpcoord[match(sel.modu,newpcoord$cluster),]
  cat("Totally",nrow(newpcoord),"single modules were linked.\n")

  xpcoord<-data.frame(x=c(ddaxt$module1,ddaxt$module2),
                      clusters=c(ddaxt$clusters,ddaxt$clusters))
  xpcoord<-xpcoord[!duplicated(xpcoord$x),]
  newpcoord<-inner_join(newpcoord,xpcoord,by=c("cluster"="x"))
  newpcoord<-newpcoord[order(newpcoord$clusters),]
  newpcoord$clusters<-paste("Clusters",newpcoord$clusters,sep = "")

  newpcoord$clusters<-factor(newpcoord$clusters,levels = unique(newpcoord$clusters))
  color_cluster<-colorCustom(length(unique(newpcoord$clusters)), pal = pal_cluster)
  names(color_cluster)<-unique(newpcoord$clusters)

  edge2s<-subset(newpcoord, select=c(cluster,x0,y0,clusters))
  edge3s<-subset(newpcoord, select=c(cluster,x0,y0,clusters))
  edge2s<-inner_join(ddaxt,edge2s,by = c("module1" = "cluster"))
  edge2s<-inner_join(edge2s,edge3s,by = c("module2" = "cluster"))
  edge2s$clusters<-edge2s$clusters.y

  edge2s<-edge2s[order(edge2s$clusters),]
  edge2s$clusters<-factor(edge2s$clusters,levels = levels(newpcoord$clusters))
  color_ecluster<-colorCustom(length(unique(edge2s$clusters)), pal = pal_cluster)
  names(color_ecluster)<-unique(edge2s$clusters)
  newpcoord$Mtext<-paste("M",newpcoord$module,sep = "")

  ppx<-ggplot()+
    geom_curve(data =edge2sp[which(edge2sp$cur=="show"),],
               curvature = curvature,size = link.pos.size,
               colour =link.pos,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sp[which(edge2sp$cur!="show"),],
               curvature = 0,size = link.pos.size,
               colour =link.pos,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur=="show"),],
               curvature = curvature,size = link.neg.size,
               colour = link.neg,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur!="show"),],
               curvature = 0,size = link.neg.size,
               colour = link.neg,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_point(data = pcoord,size=point.size,
               aes(x=x,y=y,color=Phylum,shape=role.x))+
    geom_text(data=pcoordxs,show.legend = FALSE,
              aes(x=x0,y=y0,label=Mtext,
                  size=1*pcoordxs$tsize),
              family = "serif",fontface="bold")+
    #ylim(NA,max(pcoord$x))+
    scale_color_manual(values = color_module,
                       guide=guide_legend(keywidth = 1,
                                          keyheight = 1,
                                          order=0,
                                          ncol = 1,
                                          override.aes=list(size=3)))+
    scale_shape_manual(values = shapes,name="\nTopological role",
                       guide=guide_legend(keywidth = 1,
                                          keyheight = 1,
                                          order=0,
                                          ncol = 1,
                                          override.aes=list(size=3)))+
    theme_void()+
    theme(legend.position = "right",
          text = element_text(family = "serif"),
          aspect.ratio = (3*totalncol-1)/(3*maxwidth-1))+
    labs(title="Modularity")+
    theme(plot.title=element_text(size=20,
                                  face="bold",
                                  hjust=0.5,
                                  lineheight=1.2))

  line.coord<-calcIntervalLineCoord(pcoord)

  if (display.interval.line) for (i in 1:length(line.coord$line.y)) {
    ppx <- ppx + annotate('segment',
                          linetype=display.interval.line.style[2],
                          size=as.numeric(display.interval.line.style[3]),
                          color=display.interval.line.style[1],
                        x = line.coord$xmin, xend = line.coord$xmax,
                        y = line.coord$line.y[i],
                        yend = line.coord$line.y[i])

  }

  if (display.panel.border) ppx<-ppx+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1.5,
             yend = 1.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1-2*totalncol-totalncol+1-0.5,
             yend = 1-2*totalncol-totalncol+1-0.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])

  ###plot with scatterpie
  ppxx<-ggplot()+
    geom_curve(data =edge2sp[which(edge2sp$cur=="show"),],
               curvature = curvature,size = link.pos.size,
               colour =link.pos,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sp[which(edge2sp$cur!="show"),],
               curvature = 0,size = link.pos.size,
               colour =link.pos,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur=="show"),],
               curvature = curvature,size = link.neg.size,
               colour = link.neg,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur!="show"),],
               curvature = 0,size = link.neg.size,
               colour = link.neg,alpha=link.left.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_point(data = pcoord[-which(paste0(pcoord$group,pcoord$module) %in% unique(xxa$mod)),],
               size=point.size,
               aes(x=x,y=y,color=Phylum,shape=role.x))+
    scale_color_manual(values = color_module)+
    scale_shape_manual(values = shapes,name="\nTopological role")+
    theme_void()+
    theme(legend.position = "none",
          text = element_text(family = "serif"),
          aspect.ratio = (3*totalncol-1)/(3*maxwidth-1))

  line.coord<-calcIntervalLineCoord(pcoord)

  if (display.interval.line) for (i in 1:length(line.coord$line.y)) {
    ppxx <- ppxx + annotate('segment',
                            linetype=display.interval.line.style[2],
                            size=as.numeric(display.interval.line.style[3]),
                            color=display.interval.line.style[1],
                            x = line.coord$xmin, xend = line.coord$xmax,
                            y = line.coord$line.y[i],
                            yend = line.coord$line.y[i])

  }

  if (display.panel.border) ppxx<-ppxx+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1.5,
             yend = 1.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1-2*totalncol-totalncol+1-0.5,
             yend = 1-2*totalncol-totalncol+1-0.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])

  ppxx<-ppxx+
    scatterpie::geom_scatterpie(aes(x=x0, y=y0,r=1.15),
                                data=xxa,
                                cols = "type",color="white",fill="white",
                                long_format=TRUE)+
    ggnewscale::new_scale_fill()+
    scatterpie::geom_scatterpie(aes(x=x0, y=y0,r=1.15),
                      data=xxa,
                      cols = "type",color=NA,
                      long_format=TRUE)+
    scale_fill_manual(values = color_module)+
    labs(title="Microbial composition")+
    theme(plot.title=element_text(size=20,
                                  face="bold",
                                  hjust=0.5,
                                  lineheight=1.2))

  edge.prop<-calcEdgeP(rmnet_cg)
  dat.all<-calcRectCoord(line.coord,edge.prop,rmnet_cg,barh=bar.height,pcoord)
  xx<-dat.all$xx
  bar<-dat.all$bar
  acc<-dat.all$acc

  if (!show.modulepairs.node) {
    pcoord$cluster<-paste(pcoord$group,"_M",pcoord$module,sep = "")
    pcoord<-pcoord[-which(pcoord$cluster %in% newpcoord$cluster),]
  }

  ###visualization
  pp<-ggplot()+
    geom_curve(data =edge2sp[which(edge2sp$cur=="show"),],
               curvature = curvature,size = link.pos.size,
               colour =link.pos,alpha=link.right.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sp[which(edge2sp$cur!="show"),],
               curvature = 0,size = link.pos.size,
               colour =link.pos,alpha=link.right.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur=="show"),],
               curvature = curvature,size = link.neg.size,
               colour = link.neg,alpha=link.right.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_curve(data =edge2sn[which(edge2sn$cur!="show"),],
               curvature = 0,size = link.neg.size,
               colour = link.neg,alpha=link.right.alpha,
               aes(x = x.x, y = y.x,
                   xend = x.y, yend = y.y) )+
    geom_point(data = pcoord,size=point.size,show.legend = FALSE,
               aes(x=x,y=y,color="grey",shape=role.x))+
    geom_text(data=pcoordxs,show.legend = FALSE,
              aes(x=x0,y=y0,label=Mtext,
                  size=1*pcoordxs$tsize),
              family = "serif",fontface="bold")+
    #ylim(NA,max(pcoord$x))+
    scale_color_manual(values = color_module,name="Phylum",
                       guide=guide_legend(keywidth = 1,
                                          keyheight = 1,
                                          ncol=1,
                                          order=0,
                                          override.aes=list(size=3)))+
    scale_shape_manual(values = shapesx,name="\nTopological role")+
    theme_void()+
    theme(legend.position = "right",
          text = element_text(family = "serif"),
          aspect.ratio = (3*totalncol-1)/(3*maxwidth-1))

  if (display.panel.border) pp<-pp+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1.5,
             yend = 1.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])+
    annotate(geom = "segment",
             x = line.coord$xmin,
             xend = line.coord$xmax,
             y = 1-2*totalncol-totalncol+1-0.5,
             yend = 1-2*totalncol-totalncol+1-0.5,
             linetype=display.interval.line.style[2],
             size=as.numeric(display.interval.line.style[3]),
             color=display.interval.line.style[1])

  if (display.interval.line) for (i in 1:length(line.coord$line.y)) {
    pp <- pp + annotate('segment',
                        linetype=display.interval.line.style[2],
                        size=as.numeric(display.interval.line.style[3]),
                        color=display.interval.line.style[1],
                          x = line.coord$xmin, xend = line.coord$xmax,
                          y = line.coord$line.y[i],
                          yend = line.coord$line.y[i])

  }


  pp<-pp+ggnewscale::new_scale_color()+
    geom_curve(data =edge2s,curvature = 0,size = link.clu.size,
               show.legend = FALSE,
               aes(x = x0.x, y = y0.x,
                   xend = x0.y, yend = y0.y,colour =clusters) )+
    scale_color_manual(values = color_ecluster)+
    ggnewscale::new_scale_color()+
    geom_point(data = newpcoord,size=point.clu.size,
               aes(x=x0,y=y0,color=clusters))+
    geom_point(data = newpcoord,size=point.clu.size,
               shape=1,stroke=1,show.legend = FALSE,
               aes(x=x0,y=y0,fill="#000000"))+
    geom_point(data = newpcoord,size=(point.clu.size+3),
               shape=1,show.legend = FALSE,
               aes(x=x0,y=y0,color=clusters))+
    scale_color_manual(values = color_cluster,name = "Module-paired clusters",
                       guide=guide_legend(keywidth = 1,
                                          keyheight = 1,
                                          order=0,
                                          ncol = 1,
                                          override.aes=list(size=3)))+
    geom_text(data=newpcoord,aes(x=x0,y=y0,label=Mtext,
                                 size=1),show.legend = FALSE,
              color="white",
              family = "serif",fontface="bold")+
    labs(title="Module pairs preserved")+
    theme(plot.title=element_text(size=20,
                                    face="bold",
                                    hjust=0.5,
                                    lineheight=1.2))


  diff.mean<-data.frame(x=rep(line.coord$text.x,length(line.coord$text.y)),
                        y=line.coord$text.y,
                        gname=unique(pcoord$group))

  diff.mean$gpname<-xlabname
  tp<-ggplot(diff.mean) +
    geom_text(aes(x=0,y=y,label = gpname),family="serif",
              angle=90,
              fontface = "bold",inherit.aes = FALSE,size = 5) +
    ylim(1-2*totalncol-totalncol+1,1)+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  col<-c(link.pos,link.neg)
  names(col)<-c("posnew","negnew")
  bp<-ggplot(xx,aes(x=0,y=1))+
    theme(text = element_text(family = "serif"),
          legend.position = "none")+
    ylim(1-2*totalncol-totalncol+1,1)+theme_void()

  for (j in 1:nrow(acc)) {
    bp<-bp+
      annotate(geom = "rect",
               xmin = 3,
               xmax = 6,
               fill=col[1],
               ymin = acc[j,]$y.m,
               ymax = acc[j,]$y.h,
               alpha = link.left.alpha)+
      annotate(geom = "rect",
               xmin = 3,
               xmax = 6,
               fill=col[2],
               ymin = acc[j,]$y.l,
               ymax = acc[j,]$y.m,
               alpha = link.left.alpha)+
      annotate("text", family="serif",font="bold",
               x = 4.5,angle=90,color=color_bar.text,
               y = (acc[j,]$y.m+acc[j,]$y.h)/2,
               label = paste(bar[j,]$posnew,"%",sep = ""))+
      annotate("text",  family="serif",font="bold",
               x = 4.5,angle=90,color=color_bar.text,
               y = (acc[j,]$y.m+acc[j,]$y.l)/2,
               label = paste(bar[j,]$negnew,"%",sep = ""))
}

  library(patchwork)
  xp<-tp+bp+ppx+ppxx+pp+plot_layout(ncol = 5,guides='collect',
                                    widths = c(0.5,0.3,maxwidth,maxwidth,maxwidth),
                                    byrow = TRUE)

  export_path2<-paste(export_path,"/module preserve",sep = "")
  dir.create(export_path2, recursive = TRUE)
  if (pdf.sel) {
  ggsave(paste(export_path2,"/module preserve (RGB).pdf",sep = ""),
         colormodel="srgb",
         width = 3*maxwidth+3,height = totalncol,xp)
  ggsave(paste(export_path2,"/module preserve.tiff",sep = ""),
         width = 3*maxwidth+3,height = totalncol,xp)
  } else {
  ggsave(paste(export_path2,"/module preserve.pdf",sep = ""),xp)
  ggsave(paste(export_path2,"/module preserve.tiff",sep = ""),xp)
  }
  cat("\n","Modules preserved between networks have been highlighted and exported to ",export_path2,sep = "","\n")
  return(xp)
}

"calcIntervalLineCoord" <- function(pcoord) {

  ycd<-c()
  yycd<-c()
  for (gt in unique(pcoord$group)) {
    sel.pcd<-pcoord[which(pcoord$group==gt),]
    lowcd<-sel.pcd$y0%>%min()-1
    higcd<-sel.pcd$y0%>%max()+1
    {
      ycd<-rbind(ycd,lowcd)
      yycd<-rbind(yycd,higcd)
    }
  }

  line.y<-c()
  for (kk in 1:(length(unique(pcoord$group))-1)) {
    hasu<-(ycd[kk]+yycd[kk+1])/2
    {
      line.y<-c(line.y,hasu)
    }
  }

  text.y<-c()
  for (kk in 1:(length(unique(pcoord$group)))) {
    hasux<-(ycd[kk]+yycd[kk])/2
    {
      text.y<-c(text.y,hasux)
    }
  }

  xmin<-pcoord$x0%>%min()-1
  xmax<-pcoord$x0%>%max()+1

  return(list(line.y=line.y,xmin=xmin,xmax=xmax,text.y=text.y,text.x=xmin-1))
}


"rand.remov.once"<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }

  sp.meanInteration<-colMeans(net.stength)

  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.

  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")

  remain.percent
}

"rand.remov2.once"<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw

  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]

    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }

    sp.meanInteration<-colMeans(net.stength)


    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0)
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]

      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }

      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }

    }

    remain.percent<-length(sp.ra.new)/length(sp.ra)

    remain.percent}
}

"rmsimu"<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  res<-t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))

  rem<-t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
  return(list(result=res,remains=rem))
}

"rm2simu"<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  res<-t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))

  rem<-t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
  return(list(result=res,remains=rem))
}

"calcRMTnetRandrm" <- function(rmnet_igraph,rmnet_cg,split_otu,
                               nperm=100,
                               export_path="cs2/microbial network analysis") {
  export_path<-paste(export_path,"/microbial network analysis",sep = "")
  dir.create(paste(export_path,"/network stability/random removal",sep = ""), recursive = TRUE)

  net<-rmnet_igraph$net
  rmt_mat<-lapply(net, function(x) {x$cleaned.matrix})
  node.attri<-rmnet_cg$node_table
  thres<-rmnet_igraph$result

  currentdat<-data.frame()
  for (i in 1:length(net)) {
    gname<-names(net)[i]
    cormatrix<-rmt_mat[[i]]

    node.attris<-node.attri[which(node.attri$group==i),]

    split_otus<-split_otu[[i]]
    split_otus$name<-rownames(split_otus)
    otutab<-split_otus[which(split_otus$name %in% node.attris$name),]
    otutab<-subset(otutab,select=-name)
    comm<-t(otutab)
    sp.ra<-colMeans(comm)/sum(otutab[,1])

    cormatrix2<-cormatrix*(abs(cormatrix)>=thres[[i]])  #only keep links above the cutoff point
    cormatrix2[is.na(cormatrix2)]<-0
    diag(cormatrix2)<-0    #no links for self-self
    sum(abs(cormatrix2)>0)/2  #this should be the number of links.
    sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

    network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
    sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
    if (isTRUE(unique(row.names(network.raw)==names(sp.ra2)))) cat("\nEvergthing is OK !!!","\n") else message("\nSome problems seem to happen !!!") #check if matched

    Weighted.simu<-rmsimu(netRaw=network.raw,
                          rm.p.list=seq(0.05,1,by=0.05),
                          sp.ra=sp.ra2,
                          abundance.weighted=T,nperm=nperm)
    Unweighted.simu<-rmsimu(netRaw=network.raw,
                            rm.p.list=seq(0.05,1,by=0.05),
                            sp.ra=sp.ra2,
                            abundance.weighted=F,nperm=nperm)

    dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                     rbind(Weighted.simu$result,Unweighted.simu$result),
                     weighted=rep(c("weighted","unweighted"),each=20),
                     group=rep(gname,40))
    dat2<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                     rbind(Weighted.simu$remains,Unweighted.simu$remains),
                     weighted=rep(c("weighted","unweighted"),each=20),
                     group=rep(gname,40))

    currentdat1<-dat1
    currentdat2<-dat2
    file1=paste(export_path,"/network stability/random removal/(",gname,") random_removal_originresult(only mean&se).csv",sep = "")
    file2=paste(export_path,"/network stability/random removal/(",gname,") random_removal_originresult(100 times permutations).csv",sep = "")

    write.csv(currentdat1,file1)
    write.csv(currentdat2,file2)

    {
      currentdat<-rbind(currentdat,currentdat2)
    }
  }

  return(currentdat)
}

"calcRMTnetTargetrm" <- function(rmnet_igraph,rmnet_cg,split_otu,
                                 nperm=100,otu.sel=NULL,parallel=FALSE,parallel.node.num=10,
                                 export_path="cs2/microbial network analysis") {
  export_path<-paste(export_path,"/microbial network analysis",sep = "")
  dir.create(paste(export_path,"/network stability/target removal",sep = ""), recursive = TRUE)

  net<-rmnet_igraph$net
  rmt_mat<-lapply(net, function(x) {x$cleaned.matrix})
  node.attri<-rmnet_cg$node_table
  thres<-rmnet_igraph$result

  currentdat<-data.frame()
  modulehub.num<-c()
  i=1
  for (i in 1:length(net)) {
    gname<-names(net)[i]
    cormatrix<-rmt_mat[[i]]

    node.attris<-node.attri[which(node.attri$group==i),]

    split_otus<-split_otu[[i]]
    split_otus$name<-rownames(split_otus)
    otutab<-split_otus[which(split_otus$name %in% node.attris$name),]
    otutab<-subset(otutab,select=-name)
    comm<-t(otutab)
    sp.ra<-colMeans(comm)/sum(otutab[,1])

    cormatrix2<-cormatrix*(abs(cormatrix)>=thres[[i]])  #only keep links above the cutoff point
    cormatrix2[is.na(cormatrix2)]<-0
    diag(cormatrix2)<-0    #no links for self-self
    sum(abs(cormatrix2)>0)/2  #this should be the number of links.
    sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

    network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
    sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
    if (isTRUE(unique(row.names(network.raw)==names(sp.ra2)))) cat("\nEvergthing is OK !!!","\n") else message("\nSome problems seem to happen !!!") #check if matched

    node.attrix<-node.attri[which(node.attri$group==i),]
   if (is.null(otu.sel)) {
    module.hub<-as.character(node.attrix$name[node.attrix$zi >= 2.5 | node.attrix$pi >= 0.62])
    tem<-length(module.hub)

    if(length(module.hub)<=1){
      ret3<-node.attrix%>%arrange(desc(degree))
      tem<-round(length(ret3$degree)*0.05,0)
      module.hub<-ret3$name[1:tem]
    }

    if (parallel) {
      ret3<-node.attrix%>%arrange(desc(degree))
      tem<-parallel.node.num
      module.hub<-ret3$name[1:tem]
    }

    Weighted.simu<-rm2simu(netRaw=network.raw,
                           rm.p.list=1:length(module.hub),
                           keystonelist=module.hub,
                           sp.ra=sp.ra2,
                           abundance.weighted=T,nperm=nperm)
    Unweighted.simu<-rm2simu(netRaw=network.raw,
                             rm.p.list=1:length(module.hub),
                             keystonelist=module.hub, sp.ra=sp.ra2,
                             abundance.weighted=F,nperm=nperm)
   } else {
     module.hub<-otu.sel
     tem<-length(module.hub)

     Weighted.simu<-rm2simu(netRaw=network.raw,
                            rm.p.list=1:length(module.hub),
                            keystonelist=module.hub,
                            sp.ra=sp.ra2,
                            abundance.weighted=T,nperm=nperm)
     Unweighted.simu<-rm2simu(netRaw=network.raw,
                              rm.p.list=1:length(module.hub),
                              keystonelist=module.hub, sp.ra=sp.ra2,
                              abundance.weighted=F,nperm=nperm)
    }

    dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                     rbind(Weighted.simu$result,Unweighted.simu$result),
                     weighted=rep(c("weighted","unweighted"),
                                  each=length(module.hub)),
                     group=rep(gname,2*length(module.hub)))

    dat2<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                     rbind(Weighted.simu$remains,Unweighted.simu$remains),
                     weighted=rep(c("weighted","unweighted"),
                                  each=length(module.hub)),
                     group=rep(gname,2*length(module.hub)))

    currentdat1<-dat1
    currentdat2<-dat2
    file1=paste(export_path,"/network stability/target removal/(",gname,") target_removal_originresult(only mean&se).csv",sep = "")
    file2=paste(export_path,"/network stability/target removal/(",gname,") target_removal_originresult(100 times permutations).csv",sep = "")

    write.csv(currentdat1,file1)
    write.csv(currentdat2,file2)

    {
      currentdat<-rbind(currentdat,currentdat2)
      modulehub.num<-rbind(modulehub.num,tem)
    }
  }

  return(list(currentdat=currentdat,modulehub.num=modulehub.num))
}

"plotRMTnetRandrm" <- function(randrm.data,
                               geom=c("boxplot","barplot"),
                               select.rm.prop="50%",
                               strictmod=TRUE,
                               method="t.test",
                               comparison=my_comparisons,
                               xlabname=c("CS0","CS100","CA100","CD100"),
                               yaxis.italic=TRUE,
                               color_group=colorCustom(5,pal = "ywbu"),
                               export_path="cs2/microbial network analysis") {
  geom<-match.arg(geom)
  select.rm.prop<-as.numeric(sub("%","",select.rm.prop))/100
  sel.randrm.data<-randrm.data[which(randrm.data$Proportion.removed %in% select.rm.prop),]
  sel.randrm.data<-sel.randrm.data[,-1]
  wei.sel.randrm.data<-sel.randrm.data[which(sel.randrm.data$weighted=="weighted"),]
  unwei.sel.randrm.data<-sel.randrm.data[which(sel.randrm.data$weighted=="unweighted"),]

  wei.sel.randrm.data<-subset(wei.sel.randrm.data,select=-weighted)
  unwei.sel.randrm.data<-subset(unwei.sel.randrm.data,select=-weighted)

  wei.sel.randrm.datax<-melt(wei.sel.randrm.data,
                             id.vars = c("group"),
                             variable.name = c('rep'),#聚合变量的新列名
                             value.name = 'value')
  wei.sel.randrm.datax$rep<-paste(wei.sel.randrm.datax$group,sub("X","",wei.sel.randrm.datax$rep),sep = "")
  wei.sel.randrm.datax$group<-factor(wei.sel.randrm.datax$group,levels = unique(wei.sel.randrm.datax$group))
  wei.sel.randrm.datax<-wei.sel.randrm.datax[order(wei.sel.randrm.datax$group),]

  unwei.sel.randrm.datax<-melt(unwei.sel.randrm.data,
                               id.vars = c("group"),
                               variable.name = c('rep'),#聚合变量的新列名
                               value.name = 'value')
  unwei.sel.randrm.datax$rep<-paste(unwei.sel.randrm.datax$group,sub("X","",unwei.sel.randrm.datax$rep),sep = "")
  unwei.sel.randrm.datax$group<-factor(unwei.sel.randrm.datax$group,levels = unique(unwei.sel.randrm.datax$group))
  unwei.sel.randrm.datax<-unwei.sel.randrm.datax[order(unwei.sel.randrm.datax$group),]


  plot.data<-cbind(wei.sel.randrm.datax,unwei.sel.randrm.datax)
  rownames(plot.data)<-plot.data$rep
  plot.data<-subset(plot.data,select=-c(group,rep))
  plot.data<-subset(plot.data,select=-c(group,rep))
  colnames(plot.data)<-c("Weighted_robustness","unweighted_robustness")

  newfile<-paste(export_path,"/microbial network analysis/network stability/random removal/random remove ",select.rm.prop*100," percent nodes",sep = "")
  dir.create(newfile, recursive = TRUE)

  file1=paste(newfile,"/data.csv",sep = "")
  write.csv(plot.data,file1)

  otu_table<-NULL
  taxon_table<-NULL
  tree<-NULL
  params<-plot.data

  mchat<-setParamchat(otu_table,taxon_table,tree=tree,params=params)
  microchatParamobj<-calcMicrochatParam(mchat,
                                        export_path=newfile)

  library(multcomp)
  calcMicrochatParamTable(microchatParamobj,
                          strictmod=strictmod,
                          method=method,
                          comparison=comparison,
                          export_path=newfile)

  if (geom=="boxplot") plotMicrochatParamMutiBoxplot(microchatParamobj,
                                                     geom.line=FALSE,
                                                     xlabname=xlabname,
                                                     yaxis.italic=yaxis.italic,
                                                     strictmod=strictmod,
                                                     method=method,
                                                     comparison=comparison,
                                                     color_group=color_group,
                                                     export_path=newfile)


  if (geom=="barplot") plotMicrochatParamMutiBarplot(microchatParamobj,
                                                     xlabname=xlabname,
                                                     yaxis.italic=yaxis.italic,
                                                     strictmod=strictmod,
                                                     method=method,
                                                     comparison=comparison,
                                                     color_group=color_group,
                                                     export_path=newfile)

  cat("\nPlots have been exported to ",newfile,"\n",sep = "")
}

"plotRMTnetTargrm" <- function(target.data,
                               geom=c("boxplot","barplot"),
                               select.rm.modhub.num=7,
                               strictmod=TRUE,
                               method="t.test",
                               comparison=my_comparisons,
                               xlabname=c("CS0","CS100","CA100","CD100"),
                               yaxis.italic=TRUE,
                               color_group=colorCustom(5,pal = "ywbu"),
                               export_path="cs2/microbial network analysis") {
  geom<-match.arg(geom)


  randrm.data<-target.data$currentdat
  if (select.rm.modhub.num < min(target.data$modulehub.num)) {
    select.rm.prop <- select.rm.modhub.num
  } else (
    select.rm.prop <- min(target.data$modulehub.num)
  )

  sel.randrm.data<-randrm.data[which(randrm.data$Number.hub.removed %in% select.rm.prop),]
  sel.randrm.data<-sel.randrm.data[,-1]
  wei.sel.randrm.data<-sel.randrm.data[which(sel.randrm.data$weighted=="weighted"),]
  unwei.sel.randrm.data<-sel.randrm.data[which(sel.randrm.data$weighted=="unweighted"),]

  wei.sel.randrm.data<-subset(wei.sel.randrm.data,select=-weighted)
  unwei.sel.randrm.data<-subset(unwei.sel.randrm.data,select=-weighted)

  wei.sel.randrm.datax<-melt(wei.sel.randrm.data,
                             id.vars = c("group"),
                             variable.name = c('rep'),#聚合变量的新列名
                             value.name = 'value')
  wei.sel.randrm.datax$rep<-paste(wei.sel.randrm.datax$group,sub("X","",wei.sel.randrm.datax$rep),sep = "")
  wei.sel.randrm.datax$group<-factor(wei.sel.randrm.datax$group,levels = unique(wei.sel.randrm.datax$group))
  wei.sel.randrm.datax<-wei.sel.randrm.datax[order(wei.sel.randrm.datax$group),]

  unwei.sel.randrm.datax<-melt(unwei.sel.randrm.data,
                               id.vars = c("group"),
                               variable.name = c('rep'),#聚合变量的新列名
                               value.name = 'value')
  unwei.sel.randrm.datax$rep<-paste(unwei.sel.randrm.datax$group,sub("X","",unwei.sel.randrm.datax$rep),sep = "")
  unwei.sel.randrm.datax$group<-factor(unwei.sel.randrm.datax$group,levels = unique(unwei.sel.randrm.datax$group))
  unwei.sel.randrm.datax<-unwei.sel.randrm.datax[order(unwei.sel.randrm.datax$group),]


  plot.data<-cbind(wei.sel.randrm.datax,unwei.sel.randrm.datax)
  rownames(plot.data)<-plot.data$rep
  plot.data<-subset(plot.data,select=-c(group,rep))
  plot.data<-subset(plot.data,select=-c(group,rep))
  colnames(plot.data)<-c("Weighted_robustness","unweighted_robustness")

  newfile<-paste(export_path,"/microbial network analysis/network stability/target removal/target remove ",select.rm.prop,"  module hub(s)",sep = "")
  dir.create(newfile, recursive = TRUE)

  file1=paste(newfile,"/data.csv",sep = "")
  write.csv(plot.data,file1)

  otu_table<-NULL
  taxon_table<-NULL
  tree<-NULL
  params<-plot.data

  mchat<-setParamchat(otu_table,taxon_table,tree=tree,params=params)
  microchatParamobj<-calcMicrochatParam(mchat,
                                        export_path=newfile)

  library(multcomp)
  calcMicrochatParamTable(microchatParamobj,
                          strictmod=strictmod,
                          method=method,
                          comparison=comparison,
                          export_path=newfile)

  if (geom=="boxplot") plotMicrochatParamMutiBoxplot(microchatParamobj,
                                                     geom.line=FALSE,
                                                     xlabname=xlabname,
                                                     yaxis.italic=yaxis.italic,
                                                     strictmod=strictmod,
                                                     method=method,
                                                     comparison=comparison,
                                                     color_group=color_group,
                                                     export_path=newfile)


  if (geom=="barplot") plotMicrochatParamMutiBarplot(microchatParamobj,
                                                     xlabname=xlabname,
                                                     yaxis.italic=yaxis.italic,
                                                     strictmod=strictmod,
                                                     method=method,
                                                     comparison=comparison,
                                                     color_group=color_group,
                                                     export_path=newfile)

  cat("\nPlots have been exported to ",newfile,"\n",sep = "")
}

"plotRMTnetStablity" <- function(randrm.data, type="weighted",
                                 my_comparisons,
                                 color_group=colorCustom(5,pal = "na3")
) {
  if (class(randrm.data)=="list") randrm.datax <- randrm.data[[1]] else randrm.datax <- randrm.data
  randrm.datax<-randrm.datax[which(randrm.datax$weighted==type),]
  colnames(randrm.datax)[1]<-"var"
  if (class(randrm.data)!="list") randrm.datax$var<-randrm.datax$var*100
  randrm.datax<-subset(randrm.datax,select=c(group,var,2:101))
  random.num<-randrm.datax[3:102]
  randrm.datax$mean<-rowSums(random.num)/100
  randrm.datax<-subset(randrm.datax,select=c(group,var,mean,3:102))
  #randrm.datax<-randrm.datax[which(randrm.datax$var>=50),]
  if (class(randrm.data)=="list") randrm.datax<-randrm.datax[which(randrm.datax$var %in% seq(5,100,5)),]

  randrm.datax$fit<-0
  for (tk in unique(randrm.datax$group)) {
    sas<-getOnePoly(randrm.datax[which(randrm.datax$group==tk),])
    randrm.datax[which(randrm.datax$group==tk),]$fit<-sas$a*randrm.datax[which(randrm.datax$group==tk),]$var+sas$b
  }
  randrm.datax<-subset(randrm.datax,select=c(group,var,mean,fit,4:103))

  color.use<-color_group
  names(color.use)<-unique(randrm.datax$group)
  p<-ggplot(randrm.datax,aes(x=var,y=fit))+
    geom_line(aes(color=group),formula = y ~ x, method="glm",size=1,
              lineend = "round")+
    geom_point(aes(x=var,y=mean,color=group),shape=21,size=2,stroke=2)+
    scale_fill_manual(values =color.use )+
    scale_color_manual(values =color.use )+
    theme(aspect.ratio = 1)

  alphadiv.use<-subset(randrm.datax, select = c(1,3))
  data_poi<-alphadiv.use
  colnames(data_poi)[2]<-"value"
  sig_label_new=NULL
  sta<-lapply(my_comparisons, function(x){
    alphadiv.use.selected1<-data_poi[
      which(data_poi$group==x[1] | data_poi$group==x[2]),]
    if (length(which(data_poi$group==x[1]))==length(which(data_poi$group==x[2]))) {
      fit<-wilcox.test(value~group, alphadiv.use.selected1, # t.tes
                       paired = TRUE, alternative = 'two.sided')
    } else {
      fit<-wilcox.test(value~group, alphadiv.use.selected1, # t.tes
                       paired = FALSE, alternative = 'two.sided')
    }
    stats<-data.frame(
      group1=x[1],
      group2=x[2],
      mean1=mean(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[1])]),
      mean2=mean(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[2])]),
      max1=max(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[1])]),
      max2=max(alphadiv.use.selected1$value[
        which(alphadiv.use.selected1$group==x[2])]),
      statistics=fit$statistic,
      p.value=fit$p.value,
      method=fit$method)
  })

  data.alpha<-data.frame()
  for (tt in 1:length(sta)) {
    data.al<-sta[[tt]]
    {
      data.alpha<-rbind(data.alpha,data.al)
    }
  }
  data.alpha$sig<-ifelse(data.alpha$p.value<0.001,"***",
                         ifelse(data.alpha$p.value<0.01,"**",
                                ifelse(data.alpha$p.value<0.05,"*","ns")))

  dsd_data<- data.frame()
  for (tk in unique(randrm.datax$group)) {
    sas<-getOnePoly(randrm.datax[which(randrm.datax$group==tk),])
    dd_data<-data.frame(a=sas$a,b=sas$b)
    dsd_data<-rbind(dsd_data,dd_data)
  }

  top_pt<-dsd_data[which(dsd_data$a==max(dsd_data$a)),]
  top_y<-top_pt$a*min(randrm.datax$var)+top_pt$b

  bb_line_a<-mean(c(max(dsd_data$a),min(dsd_data$a)))*(-1)
  bb_line_b<-max(dsd_data$b)-min(dsd_data$b)


  dsd_data$c=bb_line_a*1.6
  if (class(randrm.data)=="list")   dsd_data$c=bb_line_a*3.2
  dsd_data$d=0
  if (class(randrm.data)=="list") dsd_data$d[1]=min(dsd_data$b)  else dsd_data$d[1]=top_y-dsd_data$c[1]*min(randrm.datax$var)
  for (kk in 2:nrow(dsd_data)) {
    dsd_data$d[kk]<-dsd_data$d[kk-1]-(min(dsd_data$b)-min(dsd_data$b)/length(my_comparisons))/length(my_comparisons)
  }
  dsd_data$group<-unique(randrm.datax$group)

  interpt<-list()
  for (t in 1:length(my_comparisons)) {
    selected_compa<-my_comparisons[[t]]
    selected_data<-dsd_data[match(selected_compa,dsd_data$group),]
    selected_data$c[1]=dsd_data$c[t]
    selected_data$d[1]=dsd_data$d[t]
    selected_data$c[2]=dsd_data$c[t]
    selected_data$d[2]=dsd_data$d[t]
    for (i in 1:nrow(selected_data)) {
      a=selected_data$a[i]
      b=selected_data$b[i]
      c=selected_data$c[i]
      d=selected_data$d[i]
      selected_data$x[i]<-getIntersection(a,b,c,d)$x
      selected_data$y[i]<-getIntersection(a,b,c,d)$y
    }
    interpt[[t]]<-selected_data
  }

  for (k in 1:length(interpt)) {
    p<-p+annotate(geom = "segment",
                  x = interpt[[k]]$x[1], xend = interpt[[k]]$x[2],
                  y = interpt[[k]]$y[1], yend = interpt[[k]]$y[2],
                  size=1,
                  color = "black")+
      annotate(geom = "point",
               x = interpt[[k]]$x[1], y = interpt[[k]]$y[1],
               size = 2,
               shape = 19, # 设置形状为圆形
               fill = "black", # 设置填充颜色为黑色
               color = "black") +
      annotate(geom = "point",
               x = interpt[[k]]$x[2], y = interpt[[k]]$y[2],
               size = 2,
               shape = 19, # 设置形状为圆形
               fill = "black", # 设置填充颜色为黑色
               color = "black")+
      annotate(geom="text",fontface = "bold",
               x = (interpt[[k]]$x[1]+interpt[[k]]$x[2])/2*0.99,
               y = (interpt[[k]]$y[1]+interpt[[k]]$y[2])/2*1.01,
               angle=50,size=6,color = "black",
               label = data.alpha$sig[k])
  }

  wc<-dsd_data$b
  names(wc)<-dsd_data$group
  wc<-wc[order(wc,decreasing = TRUE)]
  wcc<-c()
  for (t in 1:length(wc)) {
    wcc[2*t-1]<-names(wc[t])
    wcc[2*t]<-">"
  }
  wcc<-wcc[1:(length(wcc)-1)]
  halfpos<-(max(randrm.datax$var)-min(randrm.datax$var))/2

  text_stab<-data.frame(var=wcc,
                        pos=seq(from=min(randrm.datax$var),
                                to=min(randrm.datax$var)+halfpos,length.out=length(wcc))
  )
  text_stab$color<-color.use[match(text_stab$var,names(color.use))]
  text_stab$color[which(text_stab$var==">")]<-"black"
  color.use1<-text_stab$color
  names(color.use1)<-wcc
  p<-p+geom_text(data=text_stab,
                 aes(x=pos,y=0,label=var,color=var),
                 family="serif",
                 fontface="bold")+
    ggnewscale::new_scale_color()+
    scale_color_manual(values = color.use1)

  if (class(randrm.data)=="list") {
    p<-p+labs(x="Percentage of target removal",y="Robustness",subtitle  = str_to_title(type))
  } else {
    p<-p+labs(x="Percentage of random removal (%)",y="Robustness",subtitle  = str_to_title(type))
  }
  p<-p+theme(legend.position =c(0.9,0.8),
             legend.background = element_blank(),
             legend.key = element_blank(),
             text = element_text(family = "serif"),
             legend.title = element_blank(),
             title = element_text(family = "serif",size = 14,face = "bold"),
             panel.background = element_blank(),
             panel.border = element_rect(size = 2,fill = NA,
                                         linetype = "solid"),
             axis.text.x = element_text(angle = 0,hjust = 0.5,size = 12),
             axis.text.y = element_text(angle = 0,hjust = 0.5,size = 12),
             axis.text.y.right = element_blank(),
             axis.ticks = element_line(lineend="round",size = 1),
             axis.ticks.y.right = element_blank())

  return(p)
}

"getOnePoly" <- function(randrm.datax) {
  model <- lm(randrm.datax$mean ~ poly(randrm.datax$var, 1, raw = TRUE), data = randrm.datax)
  mylr = function(x,y){
    x_mean = mean(x)
    y_mean = mean(y)

    a = coef(model)[2]
    b = coef(model)[1]

    # 构造线性回归方程
    f <- a * x + b

    ###生成总平方和
    sst = sum((y-y_mean)^2)
    ##残差平方和
    sse = sum((y-f)^2)
    ##回归平方和
    ssr = sum((f-y_mean)^2)

    result = c(a,b,sst,sse,ssr)
    names(result) = c('a','b','sst','sse','ssr')
    return(result)
  }

  ##横坐标为时间
  x = randrm.datax$var
  ##纵坐标为表达量
  y = randrm.datax$mean
  ###代入公式
  f = mylr(x,y)

  ###计算拟合度，回归平方和除以总平方和
  R2= f['ssr']/f['sst']

  f = mylr(x,y)
  a<-f['a']
  b<-f['b']
  a;b
  kp<-list(formula=f,a=a,b=b,r2=R2)
  return(kp)
}

"getIntersection" <-function(a,b,c,d) {
  # 定义函数f(x)和g(x)的方程
  f <- function(x) {
    a*x + b
  }

  g <- function(x) {
    c*x + d
  }

  # 求解函数交点
  intersection <- uniroot(function(x) f(x) - g(x), interval = c(-100, 100))
  intersection_point <- intersection$root

  return(list(x=intersection_point,y=f(intersection_point)))
}




"network.efficiency" <- function(graph){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  dd <- 1/shortest.paths(graph)
  diag(dd) <- NA
  efficiency <- mean(dd, na.rm=T)
  #denom <- nrow(dd)*(ncol(dd)-1)
  #sum(dd, na.rm=T)/denom
  return(efficiency)
}

"info.centrality.vertex" <- function(graph, net=NULL, verbose=F){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  if(is.null(net)) net <- network.efficiency(graph)
  if(is.numeric(net)==F){
    warning("Please ensure net is a scalar numeric")
    net <- network.efficiency(graph)
  }
  count <- c()
  for(i in 1:length(V(graph))){
    count <- c(count, (net-network.efficiency(delete.vertices(graph, i)))/net)
    if(verbose){
      print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
    }
  }
  return(count)
}

"info.centrality.network" <- function(graph, net=network.efficiency(graph), verbose=F) {
  sum(info.centrality.vertex(graph))
}

"calcRMTnetVulner" <- function(rmnet_igraph,rmnet_cg) {

  net<-rmnet_igraph$net
  rmt_mat<-lapply(net, function(x) {x$cleaned.matrix})
  node.attri<-rmnet_cg$node_table
  thres<-rmnet_igraph$result

  dd<-c()
  ee<-c()
  for (i in 1:length(net)) {
    gname<-names(net)[i]
    cormatrix<-rmt_mat[[i]]%>%as.matrix()
    cormatrix[abs(cormatrix) < thres[[i]]]<-0
    cormatrix[abs(cormatrix)>0]<-1 # adjacency matrix
    dim(cormatrix)
    g2 = graph_from_adjacency_matrix(as.matrix(cormatrix),
                                     mode="undirected",
                                     weighted = NULL,
                                     diag = FALSE,
                                     add.colnames = NULL)


    #check node number and links
    length(V(g2));length(E(g2))

    # calculate vulnerability of each node
    node.vul<-info.centrality.vertex(g2)
    node.vulx<-info.centrality.network(g2, net=network.efficiency(g2), verbose=F)
    cat("\n",gname,": ",max(node.vul), " --Max node vulnerability (Smaller value stands for higher robustness, namely higher stability.)",sep = "")
    cat("\n",gname,": ",node.vulx," --Network vulnerability (Sum of node vulnerability)",sep = "")
    cat("\n","----------------------------------------------------------------")
    {
      dd<-rbind(dd,max(node.vul))
      ee<-rbind(ee,max(node.vulx))
    }
  }

  dd<-data.frame(dd)
  ee<-data.frame(ee)

  rownames(dd)<-names(net)
  colnames(dd)<-"value"
  dd$group<-rownames(dd)
  dd<-dd[order(dd$value),]
  gp<-dd$group[1]
  for (tk in 2:nrow(dd)) {
    ggp<-dd$group[tk]
    {
      gp<-paste(gp," > ",ggp)
    }
  }
  cat("\n","Network stability: ",gp,"\n","Derived from vulnerability measured by maximum node vulnerability")

  rownames(ee)<-names(net)
  colnames(ee)<-"value"
  ee$group<-rownames(ee)
  ee<-ee[order(ee$value),]
  gp<-ee$group[1]
  for (tk in 2:nrow(ee)) {
    ggp<-ee$group[tk]
    {
      gp<-paste(gp," > ",ggp)
    }
  }
  cat("\n","Network stability: ",gp,"\n","Derived from vulnerability measured by sum of node vulnerability")

  return(list(max.node.vulnerability=dd,
              network.vulnerability=ee))
}





"plotRMTnetModPreserveCircle" <- function(modpre_data,
                                          edge.curved=0.5,
                                          sel.layout,
                                          color_group,
                                          color_cluster) {
  modpre_datax<-modpre_data
  modpre_datax$Both<-modpre_datax$Both%>%as.numeric()
  modpre_datax$P1A2<-modpre_datax$P1A2%>%as.numeric()
  modpre_datax$P2A1<-modpre_datax$P2A1%>%as.numeric()
  modpre_datax$A1A2<-modpre_datax$A1A2%>%as.numeric()

  modpre_datax$weight<-modpre_datax$Both+modpre_datax$P1A2+modpre_datax$P2A1+modpre_datax$A1A2
  modpre_datax$p_raw<-modpre_datax$p_raw%>%as.numeric()
  modpre_datax$p_adj<-modpre_datax$p_adj%>%as.numeric()
  modpre_datax$p_raw<-sprintf("%0.3f",round(modpre_datax$p_raw,3))
  modpre_datax$p_adj<-sprintf("%0.3f",round(modpre_datax$p_adj,3))
  modpre_datax[,c("group1","mod1")]<-stringr::str_split_fixed(modpre_datax$module1,"_",2)
  modpre_datax[,c("group2","mod2")]<-stringr::str_split_fixed(modpre_datax$module2,"_",2)

  sel.modu<-unique(c(modpre_datax$mod1,modpre_datax$mod2))
  df.mat<-modpre_datax%>%as.matrix()
  gg<-igraph::graph_from_edgelist(df.mat[,c(1,2)],directed = FALSE)

  if (sel.layout=="circle") layout = in_circle()
  if (sel.layout=="nicely") layout = nicely()
  if (sel.layout=="random") layout = randomly()
  if (sel.layout=="fr") layout = with_fr()
  if (sel.layout=="kk") layout = with_kk()
  if (sel.layout=="sphere") layout = on_sphere()
  if (sel.layout=="dh") layout = with_dh()
  if (sel.layout=="lgl") layout = with_lgl()
  if (sel.layout=="tree") layout = as_tree()
  if (sel.layout=="grid") layout = on_grid()
  if (sel.layout=="graphopt") layout = with_graphopt()
  if (sel.layout=="gem") layout = with_gem()
  if (sel.layout=="mds") layout = with_mds()

  coords <- layout_(gg, layout)
  coords_scale = scale(coords)
  names(color_group)<-unique(c(modpre_datax$group1,modpre_datax$group2))
  names(color_cluster)<-unique(c(modpre_datax$clusters))

  E(gg)$weight<-normalize(modpre_datax$Both)$x*10
  igraph::V(gg)$color <- color_group[match(substr(igraph::V(gg)$name,start = 1,stop = 2),
                                           names(color_group))]

  color_clusterx<-color_cluster[match(modpre_datax$clusters,
                                      names(color_cluster))]
  igraph::E(gg)$color <- grDevices::adjustcolor(color_clusterx, 0.6)

  "radian.rescale" <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  if (sel.layout=="circle") label.locs <- radian.rescale(x = 1:length(igraph::V(gg)),
                                                         direction = -1, start = 0)



  V(gg)$size<-modpre_datax$weight/min(modpre_datax$weight)*6
  E(gg)$width<-E(gg)$weight
  if (sel.layout=="circle") plot(gg,edge.curved=edge.curved,
                                 rescale = TRUE,
                                 layout=coords,
                                 vertex.size=V(gg)$size,
                                 vertex.shape="circle",
                                 vertex.label.degree = label.locs,
                                 vertex.label.dist = 4,
                                 vertex.label=vertex_attr(gg)$name,
                                 vertex.frame.color = "white",
                                 vertex.label.family = "serif",
                                 edge.label.family = "serif")

  if (sel.layout!="circle") plot(gg,edge.curved=edge.curved,
                                 rescale = TRUE,
                                 layout=coords,
                                 vertex.size=V(gg)$size,,
                                 vertex.shape="circle",
                                 vertex.label.dist = 3,
                                 vertex.label=vertex_attr(gg)$name,
                                 vertex.frame.color = "white",
                                 vertex.label.family = "serif",
                                 edge.label.family = "serif")

  sas<-recordPlot()
  return(sas)
}

"calcRectCoord" <- function(line.coord,edge.prop,rmnet_cg,barh=NULL,pcoord) {
  kh.all.sel<-c()
  for (ttk in 1:length(line.coord$text.y)) {
    pcd.sel<-pcoord[which(pcoord$group==rmnet_cg$netname[ttk]),]
    kh.sel<-max(pcd.sel$y0)-min(pcd.sel$y0)
    {
      kh.all.sel<-c(kh.all.sel,kh.sel)
    }
  }

  acc<-c()
  bar<-data.frame()
  for (ttk in 1:length(line.coord$text.y)) {
    bar.sel<-edge.prop[which(edge.prop$group==rmnet_cg$netname[ttk]),]
    text.c<-line.coord$text.y[ttk]
    if (is.null(barh)) barhx<-kh.all.sel[ttk] else barhx<-min(kh.all.sel)
    act.c<-calcRTc(bar.sel,text.c,barhx)
    act.c<-data.frame(y.h=act.c[1],
                      y.m=act.c[2],
                      y.l=act.c[3])
    rownames(act.c)<-rmnet_cg$netname[ttk]
    {
      acc<-rbind(acc,act.c)
      bar<-rbind(bar,bar.sel)
    }
  }

  bar$value<-sprintf("%0.2f",round(bar$value,2))
  bar<-spread(bar,key=variable,value = value)
  xx<-acc%>%rownames_to_column(var = "group")
  reshape2::melt(xx)

  return(list(xx=xx,bar=bar,acc=acc))
}

"calcRTc" <- function(bar.sel,text.c,barh) {
  y.h<-barh/2+text.c
  y.l<-text.c-barh/2
  y.m<-y.h-bar.sel$value[1]/100*barh
  return(c(y.h,y.m,y.l))
}
"calcEdgeP" <- function(rmnet_cg) {
  edge.prop<-table(rmnet_cg$edge_table$group,rmnet_cg$edge_table$cor)%>%data.frame()
  edge.prop$Var1<-rep(rmnet_cg$netname,2)
  edge.prop<-spread(edge.prop,key = Var2,value = Freq)
  colnames(edge.prop)<-c("group","neg","pos")
  edge.prop$group<-factor(edge.prop$group,levels = rmnet_cg$netname)
  rownames(edge.prop)<-1:nrow(edge.prop)
  edge.prop$negnew<-edge.prop$neg/(edge.prop$neg+edge.prop$pos)*100
  edge.prop$posnew<-edge.prop$pos/(edge.prop$neg+edge.prop$pos)*100

  edge.prop<-reshape2::melt(edge.prop[,c(1,4,5)])
  edge.prop$variable<-factor(edge.prop$variable,levels = c("posnew","negnew"))
  edge.prop<-edge.prop[order(edge.prop$group,edge.prop$variable),]
  return(edge.prop)
}
