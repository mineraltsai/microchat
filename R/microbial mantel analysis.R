"plot_mantel_heatmap"<-function(data12,
                                edge.curve.min = -0.3,
                                edge.curve.max = 0.3,add.sig=TRUE,
                                layout = c("lower","upper"),
                                grid.background =TRUE,
                                mantel.node.pos.para=TRUE, ###false,不平行，则需要调整下面参数
                                mantel.node.x.pos=c(5,7,9,11),
                                mantel.node.y.pos=c(13,11,9,7),
                                adjust.diag=TRUE, ###是否调整对角线距离
                                adjust.diag.para=TRUE, ###是否平行调整对角线距离
                                adjust.diag.para.dist=1,  ###是否平行调整对角线距离为1
                                adjust.node = c(5,7),
                                adjust.diag.dist = c(3,4),
                                diag.node.size = 5,
                                heatmap.node.size.ratio=15,
                                mantel.node.size.para=TRUE,
                                mantel.node.size=5,
                                mantel.node.size.custom=rep(5,4),  ###5为大小，4为数量
                                heatmap.shape = 1,
                                diag.node.shape.order = 1,
                                mantel.node.shape.order = c(1,2,3,4),
                                adjust.text.x = 0.8,
                                adjust.text.y = 0.8,
                                grid.color = "grey90",
                                diag.frame.color = "red",
                                mantel.frame.color="purple",
                                space.h=1.6,
                                space.v=1.5,
                                node.color.postive = "blue",
                                node.color.mid = "white",
                                node.color.negative = "red",
                                diag.node.color = "red",
                                mantel.node.color = "purple",
                                label.cex=1,
                                edge.size.max = 8,
                                edge.size.center= 4,
                                edge.size.min = 2,
                                arrow.width = 0,
                                arrow.size = 0,
                                pbreaks = c(-Inf,0.001, 0.01, 0.05, Inf),
                                plabels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"),
                                color.p=c("grey90","red","blue","yellow"),
                                rbreaks = c(-Inf,0, 0.4, 0.8, Inf),
                                rlabels = c("< 0","0 - 0.2", "0.2 - 0.8", ">= 0.8"),
                                size.r=c(0.5,2,5,10),
                                node.label.color=NULL) {

  if (grid.background) {

    mantel_data<-mantel_data(data12)
    net3<-mantel_data$corr%>%abs()
    corr<-mantel_data$corr
    corp<-mantel_data$corp

    mat11<-coord_matrix(net3)
    mat1<-mat11$mat1
    mat1[mat1 != 0]<-0
    diag<-mat11$diag
    diag<-rbind(1,diag)

    ###添加对角线相关性
    netc<-mantel_data$diagcorr
    netp<-mantel_data$diagp

    mat2<-rbind(cbind(mat1,matrix(0, nrow(mat1), ncol(netc))),
                matrix(0, ncol(netc),  nrow(mat1)+ncol(netc)))

    for (tt in 1:length(diag)) {
      mat2[diag[tt],((nrow(mat1)+1):(nrow(mat1)+ncol(netc)))]<-netc[tt,]
    }


    row.names(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)
    colnames(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)

    diag_l<-diag%>%as.numeric()

    colnames(mat2)[diag_l]<-row.names(corr)

    pp<-corrposcalcu(corr=corr,p.mat=corp,method = "square",type = "lower")


    g <- graph_from_adjacency_matrix(mat2, mode = "undirected",
                                     weighted = T)

    coord_layout <- pp$corrPos[,3:4]

    xlim<-coord_layout$x%>%min()
    xmax<-coord_layout$x%>%max()
    ylim<-coord_layout$y%>%min()
    ymax<-coord_layout$y%>%max()

    add.node<-ncol(netc)


    if (mantel.node.pos.para) {
      xscale<-seq(nrow(netc)*2/5,nrow(netc),length.out=add.node)
      yscale<- ymax+1+1+2/3*(ymax-3)-xscale
    } else {
      xscale<- mantel.node.x.pos
      yscale<- mantel.node.y.pos
    }

    add.coord<-data.frame(x=xscale,
                          y=yscale)

    data_coord<-rbind(coord_layout,add.coord)

    sub_net_layout <- data_coord%>%as.matrix()

    if (adjust.diag) {
      if (adjust.diag.para){
        sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist
      } else {
        sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist

        adjust.data <- sub_net_layout[c(diag), 1]%>%data.frame()

        diag_data<-data.frame(row=c(diag),
                              x=adjust.data$.)

        dd<-diag_data[adjust.node,]
        dd$x<-dd$x+adjust.diag.dist
        sub_net_layout[dd$row, 1] <-dd$x

      }
    } else {
      sub_net_layout<-sub_net_layout
    }

    if (mantel.node.size.para) {
      igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                               rep(mantel.node.size,add.node))
    } else {
      igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                               mantel.node.size.custom)
    }



    V(g)$weight[diag] <-diag.node.size

    igraph::V(g)$size <- V(g)$weight

    ###节点形状
    cat("\n","shape: ",shapes(),sep = "\n")
    diag.node.shape = shapes()[diag.node.shape.order]
    mantel.node.shape = shapes()[mantel.node.shape.order]
    ##定义所有节点形状
    shape_sel<-data.frame(shape=shapes()[heatmap.shape],
                          num=seq(1,nrow(mat2)))
    ##定义对角线节点形状
    shape_sel$shape[which(shape_sel$num %in% diag)] <-diag.node.shape
    ##定义mantel节点形状
    if (is.null(mantel.node.shape.order)) {
      shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-shapes()[1]
    } else {
      shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.node.shape

    }
    shapes<-shape_sel$shape%>%as.character()



    label.dist<-NULL
    label.locs<-NULL

    space.h<-adjust.text.x
    space.v<-adjust.text.y

    diag.label.dist<-space.h * 2.8
    mantel.label.dist<-space.h * 2.8
    label.dist<-rep(0,nrow(mat2))
    label.dist[2:nrow(corr)]<-space.h * 2.8
    label.dist[diag_l]<-diag.label.dist
    label.dist[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-mantel.label.dist*(-1)


    label.locs<-rep(0,nrow(mat2))
    label.locs[2:nrow(corr)]<-pi*(-1)
    label.locs[diag_l]<-pi*(-1)
    label.locs[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-pi

    ###边的颜色
    {netp<-t(netp)
      rd = cut(netc, breaks = rbreaks ,
               labels = rlabels)

      r_num<-length(levels(rd))
      size.use2<-size.r
      names(size.use2) <-levels(rd)
      size.use2<-size.use2%>%data.frame()
      size.use2<-rownames_to_column(size.use2,var = "color")

      rd<-data.frame(rd = cut(netc, breaks = rbreaks,
                              labels = rlabels))

      rd$rd<-as.character(rd$rd)
      rd$color<-rd$rd
      for (k in 1:r_num) {
        rd$color[which(match(rd$rd,size.use2$color[k])==1)]<-size.use2$.[k]
      }
      rd<-rd$color
      rrd<-matrix(rd,nrow(netc),ncol(netc))
      colnames(rrd)<-colnames(netc)
      rownames(rrd)<-rownames(netc)
      E(g)$weight<-rrd%>%t()
      E(g)$width <- E(g)$weight


      pd = cut(netp, breaks = pbreaks,
               labels = plabels)

      p_num<-length(levels(pd))
      color.use2<-color.p
      names(color.use2) <-levels(pd)
      color.use2<-color.use2%>%data.frame()
      color.use2<-rownames_to_column(color.use2,var = "color")

      pd<-data.frame(pd = cut(netp, breaks = pbreaks,
                              labels = plabels))
      pd$pd<-as.character(pd$pd)
      pd$color<-pd$pd
      for (k in 1:p_num) {
        pd$color[which(match(pd$pd,color.use2$color[k])==1)]<-color.use2$.[k]
      }
      pd<-pd$color
      ppd<-matrix(pd,nrow(netp),ncol(netp))
      colnames(ppd)<-colnames(netp)
      rownames(ppd)<-rownames(netp)

      E(g)$color<-ppd}
    ###节点的颜色
    {node.display<- pp$corrPos[,5:6]

      color1_num<-node.display[which(node.display$corr<0),]%>%rownames()%>%length()
      color3_num<-node.display[which(node.display$corr>0),]%>%rownames()%>%length()

      node.display$color<-1
      node.display$frame.color<-1
      node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color<-colorRampPalette(c(node.color.postive,node.color.mid),bias=1)(color1_num)
      node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color<-colorRampPalette(c(node.color.mid,node.color.negative),bias=1)(color3_num)

      node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$frame.color<-node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color
      node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$frame.color<-node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color

      node.display[which(node.display$corr==0),]$color<-diag.node.color

      add.node.color<-rep(mantel.node.color,add.node)
      add.node.data<-data.frame(corr=rep(1,add.node),
                                p.value=rep(0.05,add.node),
                                color=add.node.color,
                                frame.color=add.node.color)


      node.display_comb<-rbind(node.display,add.node.data)

      igraph::V(g)$color <- node.display_comb$color
      igraph::V(g)$frame.color <- node.display_comb$frame.color
    }
    ###grid节点大小
    V(g)$frame.size<-max(V(g)$size)

    newsize<-V(g)$frame.size
    newsize[diag] <-V(g)$size[diag]
    newsize[(nrow(mat2)-add.node+1):nrow(mat2)] <-0

    newcolor<-rep(grid.color,nrow(mat2))
    newcolor[diag] <-diag.frame.color
    newcolor[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.frame.color
    #??Plot.igraph
    ###add background grid
    g1<-g
    g2<-g1%>%add_vertices(nrow(mat2),
                          frame.color=newcolor,
                          name="",
                          frame.size=newsize,
                          size=newsize,
                          weight=newsize)


    sub_net_layout2<-rbind(sub_net_layout,sub_net_layout)
    layout=layout
    if (add.sig) {
      corp<-pp$corrPos[,1:6]
      corp$p.value<-corp$p.value%>%as.double()

      corp %>%
        mutate(p.value = case_when(
          p.value <=0.001 & p.value > 0 ~ "***",
          p.value <=0.01 & p.value > 0.001 ~ "**",
          p.value <=0.05 & p.value > 0.01 ~ "*",
          TRUE ~ ""
        )) -> corp
      corp<-corp$p.value%>%as.matrix()

      newp<-corp
      newp[diag] <-""
      newp[(nrow(mat2)-add.node+1):nrow(mat2)] <-""
      vertex_attr(g2)$name[seq(nrow(mat2)+1,nrow(mat2)+nrow(mat2))]<-newp


      label.dist1<-label.dist
      label.dist1[seq(1,nrow(mat2))]<-0
      label.dist1[diag] <- 0
      label.dist1[(nrow(mat2)-add.node+1):nrow(mat2)] <-0
      label.dist2<-c(label.dist,label.dist1)
    } else {
      label.dist2<-label.dist
    }

    if (layout == "lower") {

      edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
      for (jj in 1:add.node) {
        order <- nrow(corr)
        order1<-(order/2)%>%as.integer()
        edge.curved<-rep(1,nrow(corr))
        edge.curved[order]<-0
        edge.curved[1:order1]<-seq(edge.curve.max,0,length.out=order1)
        edge.curved[order1:order]<-seq(0,edge.curve.min,length.out=order-order1+1)
        edge.curved1[,jj]<-edge.curved
      }

      edge.curved1<-t(edge.curved1)
      E(g1)$curved<-edge.curved1

      label.locs1<-label.locs*(-1)
      label.dist1<-label.dist*(-1)

      sub_net_layout1<-sub_net_layout2%>%data.frame()
      sub_net_layout1$x<-nrow(corr)+3-sub_net_layout1$x
      sub_net_layout1$y<-nrow(corr)+3-sub_net_layout1$y
      sub_net_layout1<-sub_net_layout1%>%as.matrix()

      label.dist1<-c(label.dist1,label.dist2[(nrow(mat2)+1):(nrow(mat2)+nrow(mat2))])
      ##vertex.attributes(g2)
      plot(g2,edge.curved = edge.curved1,
           vertex.shape=shapes,
           vertex.label=vertex_attr(g2)$name,
           layout=sub_net_layout1,
           vertex.label.degree = label.locs1,
           vertex.label.dist = label.dist1,
           vertex.label.family = "serif")
    }

    if (layout == "upper") {##vertex.attributes(g2)

      edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
      for (jj in 1:add.node) {
        order <- nrow(corr)
        order1<-(order/2)%>%as.integer()
        edge.curved<-rep(1,nrow(corr))
        edge.curved[order]<-0
        edge.curved[1:order1]<-seq(edge.curve.max,0,length.out=order1)
        edge.curved[order1:order]<-seq(0,edge.curve.min,length.out=order-order1+1)
        edge.curved1[,jj]<-edge.curved
      }

      edge.curved1<-t(edge.curved1)
      E(g1)$curved<-edge.curved1

      plot(g2,edge.curved = edge.curved1,
           vertex.shape=shapes,
           vertex.label=vertex_attr(g2)$name,
           layout=sub_net_layout2,
           vertex.label.degree = label.locs,
           vertex.label.dist = label.dist2,
           vertex.label.family = "serif")}} else {

             mantel_data<-mantel_data(data12)
             net3<-mantel_data$corr%>%abs()
             corr<-mantel_data$corr
             corp<-mantel_data$corp

             mat11<-coord_matrix(net3)
             mat1<-mat11$mat1
             mat1[mat1 != 0]<-0
             diag<-mat11$diag
             diag<-rbind(1,diag)

             ###添加对角线相关性
             netc<-mantel_data$diagcorr
             netp<-mantel_data$diagp

             mat2<-rbind(cbind(mat1,matrix(0, nrow(mat1), ncol(netc))),
                         matrix(0, ncol(netc),  nrow(mat1)+ncol(netc)))

             for (tt in 1:length(diag)) {
               mat2[diag[tt],((nrow(mat1)+1):(nrow(mat1)+ncol(netc)))]<-netc[tt,]
             }


             row.names(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)
             colnames(mat2)[seq(nrow(mat1)+1,nrow(mat2))]<-colnames(netc)

             diag_l<-diag%>%as.numeric()

             colnames(mat2)[diag_l]<-row.names(corr)

             pp<-corrposcalcu(corr=corr,p.mat=corp,method = "square",type = "lower")


             g <- graph_from_adjacency_matrix(mat2, mode = "undirected",
                                              weighted = T)

             coord_layout <- pp$corrPos[,3:4]

             xlim<-coord_layout$x%>%min()
             xmax<-coord_layout$x%>%max()
             ylim<-coord_layout$y%>%min()
             ymax<-coord_layout$y%>%max()

             add.node<-ncol(netc)

             if (mantel.node.pos.para) {
               xscale<-seq(nrow(netc)*2/5,nrow(netc),length.out=add.node)
               yscale<- ymax+1+1+2/3*(ymax-3)-xscale
               ###mantel节点比同x坐标下对角线节点高4
             } else {
               xscale<- mantel.node.x.pos
               yscale<- mantel.node.y.pos
             }

             add.coord<-data.frame(x=xscale,
                                   y=yscale)

             data_coord<-rbind(coord_layout,add.coord)

             sub_net_layout <- data_coord%>%as.matrix()

             if (adjust.diag) {
               if (adjust.diag.para){
                 sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist
               } else {
                 sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1] <- sub_net_layout[c(diag,(max(diag)+1):(max(diag)+add.node)), 1]+adjust.diag.para.dist

                 adjust.data <- sub_net_layout[c(diag), 1]%>%data.frame()

                 diag_data<-data.frame(row=c(diag),
                                       x=adjust.data$.)

                 dd<-diag_data[adjust.node,]
                 dd$x<-dd$x+adjust.diag.dist
                 sub_net_layout[dd$row, 1] <-dd$x

               }
             } else {
               sub_net_layout<-sub_net_layout
             }

             if (mantel.node.size.para) {
               igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                                        rep(mantel.node.size,add.node))
             } else {
               igraph::V(g)$weight <- c(pp$corrPos$corr%>%abs()%>%as.matrix()*heatmap.node.size.ratio,
                                        mantel.node.size.custom)
             }



             V(g)$weight[diag] <-diag.node.size

             igraph::V(g)$size <- V(g)$weight

             ###节点形状
             shape.sel<-data.frame(shape=shapes(),number=1:length(shapes()))
             print(shape.sel)
             diag.node.shape = shapes()[diag.node.shape.order]
             mantel.node.shape = shapes()[mantel.node.shape.order]
             ##定义所有节点形状
             shape_sel<-data.frame(shape=shapes()[heatmap.shape],
                                   num=seq(1,nrow(mat2)))
             ##定义对角线节点形状
             shape_sel$shape[which(shape_sel$num %in% diag)] <-diag.node.shape
             ##定义mantel节点形状
             if (is.null(mantel.node.shape.order)) {
               shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-shapes()[1]
             } else {
               shape_sel$shape[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.node.shape

             }
             shapes<-shape_sel$shape%>%as.character()



             label.dist<-NULL
             label.locs<-NULL

             space.h<-adjust.text.x
             space.v<-adjust.text.y

             diag.label.dist<-space.h * 2.8
             mantel.label.dist<-space.h * 2.8
             label.dist<-rep(0,nrow(mat2))
             label.dist[2:nrow(corr)]<-space.h * 2.8
             label.dist[diag_l]<-diag.label.dist
             label.dist[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-mantel.label.dist*(-1)


             label.locs<-rep(0,nrow(mat2))
             label.locs[2:nrow(corr)]<-pi*(-1)
             label.locs[diag_l]<-pi*(-1)
             label.locs[seq(nrow(mat2)-add.node+1,nrow(mat2))]<-pi

             ###边的颜色
             {netp<-t(netp)
               rd = cut(netc, breaks = rbreaks ,
                        labels = rlabels)

               r_num<-length(levels(rd))
               size.use2<-size.r
               names(size.use2) <-levels(rd)
               size.use2<-size.use2%>%data.frame()
               size.use2<-rownames_to_column(size.use2,var = "color")

               rd<-data.frame(rd = cut(netc, breaks = rbreaks,
                                       labels = rlabels))

               rd$rd<-as.character(rd$rd)
               rd$color<-rd$rd
               for (k in 1:r_num) {
                 rd$color[which(match(rd$rd,size.use2$color[k])==1)]<-size.use2$.[k]
               }
               rd<-rd$color
               rrd<-matrix(rd,nrow(netc),ncol(netc))
               colnames(rrd)<-colnames(netc)
               rownames(rrd)<-rownames(netc)
               E(g)$weight<-rrd%>%t()
               E(g)$width <- E(g)$weight


               pd = cut(netp, breaks = pbreaks,
                        labels = plabels)

               p_num<-length(levels(pd))
               color.use2<-color.p
               names(color.use2) <-levels(pd)
               color.use2<-color.use2%>%data.frame()
               color.use2<-rownames_to_column(color.use2,var = "color")

               pd<-data.frame(pd = cut(netp, breaks = pbreaks,
                                       labels = plabels))
               pd$pd<-as.character(pd$pd)
               pd$color<-pd$pd
               for (k in 1:p_num) {
                 pd$color[which(match(pd$pd,color.use2$color[k])==1)]<-color.use2$.[k]
               }
               pd<-pd$color
               ppd<-matrix(pd,nrow(netp),ncol(netp))
               colnames(ppd)<-colnames(netp)
               rownames(ppd)<-rownames(netp)

               E(g)$color<-ppd}
             ###节点的颜色
             {node.display<- pp$corrPos[,5:6]

               color1_num<-node.display[which(node.display$corr<0),]%>%rownames()%>%length()
               color3_num<-node.display[which(node.display$corr>0),]%>%rownames()%>%length()

               node.display$color<-1
               node.display$frame.color<-1
               node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color<-colorRampPalette(c(node.color.postive,node.color.mid),bias=1)(color1_num)
               node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color<-colorRampPalette(c(node.color.mid,node.color.negative),bias=1)(color3_num)

               node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$frame.color<-node.display[which(node.display$corr<0),][order(node.display[which(node.display$corr<0),]$corr),]$color
               node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$frame.color<-node.display[which(node.display$corr>0),][order(node.display[which(node.display$corr>0),]$corr),]$color

               node.display[which(node.display$corr==0),]$color<-diag.node.color

               add.node.color<-rep(mantel.node.color,add.node)
               add.node.data<-data.frame(corr=rep(1,add.node),
                                         p.value=rep(0.05,add.node),
                                         color=add.node.color,
                                         frame.color=add.node.color)


               node.display_comb<-rbind(node.display,add.node.data)

               igraph::V(g)$color <- node.display_comb$color
               igraph::V(g)$frame.color <- node.display_comb$frame.color
             }
             ###grid节点大小
             V(g)$frame.size<-max(V(g)$size)

             newsize<-V(g)$frame.size
             newsize[diag] <-V(g)$size[diag]
             newsize[(nrow(mat2)-add.node+1):nrow(mat2)] <-0

             newcolor<-rep(grid.color,nrow(mat2))
             newcolor[diag] <-diag.frame.color
             newcolor[(nrow(mat2)-add.node+1):nrow(mat2)] <-mantel.frame.color
             #??Plot.igraph
             ###add background grid
             g1<-g

             g2<-g1%>%add_vertices(nrow(mat2),
                                   frame.color=newcolor,
                                   name="",
                                   frame.size=newsize,
                                   size=newsize,
                                   weight=newsize)


             sub_net_layout2<-rbind(sub_net_layout,sub_net_layout)

             if (add.sig) {
               corp<-pp$corrPos[,1:6]
               corp$p.value<-corp$p.value%>%as.double()
               corp %>%
                 mutate(p.value = case_when(
                   p.value <=0.001 & p.value > 0 ~ "***",
                   p.value <=0.01 & p.value > 0.001 ~ "**",
                   p.value <=0.05 & p.value > 0.01 ~ "*",
                   TRUE ~ ""
                 )) -> corp
               corp<-corp$p.value%>%as.matrix()

               newp<-corp
               newp[diag] <-""
               newp[(nrow(mat2)-add.node+1):nrow(mat2)] <-""
               vertex_attr(g2)$name[seq(nrow(mat2)+1,nrow(mat2)+nrow(mat2))]<-newp


               label.dist1<-label.dist
               label.dist1[seq(1,nrow(mat2))]<-0
               label.dist1[diag] <- 0
               label.dist1[(nrow(mat2)-add.node+1):nrow(mat2)] <-0
               label.dist2<-c(label.dist,label.dist1)
             }

             layout=layout
             if (layout == "lower") {

               edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
               for (jj in 1:add.node) {
                 order <- nrow(corr)
                 order1<-(order/2)%>%as.integer()
                 edge.curved<-rep(1,nrow(corr))
                 edge.curved[order]<-0
                 edge.curved[1:order1]<-seq(edge.curve.max,0,length.out=order1)
                 edge.curved[order1:order]<-seq(0,edge.curve.min,length.out=order-order1+1)
                 edge.curved1[,jj]<-edge.curved
               }

               edge.curved1<-t(edge.curved1)
               E(g1)$curved<-edge.curved1

               label.locs1<-label.locs*(-1)
               label.dist1<-label.dist*(-1)
               sub_net_layout1<-sub_net_layout%>%data.frame()
               sub_net_layout1$x<-nrow(corr)+3-sub_net_layout1$x
               sub_net_layout1$y<-nrow(corr)+3-sub_net_layout1$y
               sub_net_layout1<-sub_net_layout1%>%as.matrix()

               label.dist1<-c(label.dist1,label.dist2[(nrow(mat2)+1):(nrow(mat2)+nrow(mat2))])

               ##vertex.attributes(g2)
               plot(g1,edge.curved = edge.curved1,
                    vertex.shape=shapes,
                    vertex.label=vertex_attr(g2)$name,
                    layout=sub_net_layout1,
                    vertex.label.degree = label.locs1,
                    vertex.label.dist = label.dist1,
                    vertex.label.family = "serif")
             }

             if (layout == "upper") {##vertex.attributes(g2)
               edge.curved1<- matrix(rep(1,add.node),ncol(netp),nrow(netp))
               for (jj in 1:add.node) {
                 order <- nrow(corr)
                 order1<-(order/2)%>%as.integer()
                 edge.curved<-rep(1,nrow(corr))
                 edge.curved[order]<-0
                 edge.curved[1:order1]<-seq(edge.curve.max,0,length.out=order1)
                 edge.curved[order1:order]<-seq(0,edge.curve.min,length.out=order-order1+1)
                 edge.curved1[,jj]<-edge.curved
               }

               edge.curved1<-t(edge.curved1)
               E(g1)$curved<-edge.curved1

               plot(g1,edge.curved = edge.curved1,
                    vertex.shape=shapes,
                    vertex.label=vertex_attr(g2)$name,
                    layout=sub_net_layout2,
                    vertex.label.degree = label.locs,
                    vertex.label.dist = label.dist2,
                    vertex.label.family = "serif")
             }
           }

  gg <- recordPlot()
  return(gg)
}

"mantel_data" <- function(data12) {

  corr<-data12$occor.r
  corp<-data12$occor.p
  mantel<-data12$mantel

  d3<-subset(mantel,select=c(1:3))
  d4<-subset(mantel,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
  d5<-tibble::column_to_rownames(d5,var = "env")

  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")
  mantelr<-d5%>%as.matrix()%>%abs()
  mantelp<-d6%>%as.matrix()%>%abs()

  netc<-matrMerge(corr,mantelr)
  netp<-matpMerge(corp,mantelp)


  return(list(corr=corr,corp=corp,
              diagcorr=netc,diagp=netp))
}

"coord_matrix" <- function(net3) {
  mat1<-net3
  diag1<-length(row.names(mat1))
  diag<-diag1+1
  for (kk in 2:(nrow(net3)-1)) {
    matk<-net3[kk:nrow(net3),kk:nrow(net3)]

    {
      mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), nrow(matk))),
                  cbind(matrix(0, nrow(matk), ncol(mat1)),matk))
      diag_num<-length(row.names(mat1))+1
      diag<-rbind(diag,diag_num)
    }
  }
  matend<-net3[nrow(net3),nrow(net3)]
  mat1<-rbind(cbind(mat1,matrix(0, ncol(mat1), 1)),
              cbind(matrix(0, 1, ncol(mat1)),matend))
  return(list(mat1=mat1,diag=diag))
}

"matrMerge" <- function(corp,mantelp) {
  d1<-corp
  d5<-mantelp

  add.param<-colnames(d5)[1:length(colnames(d5))]

  add.num<-length(add.param)
  net_m<-d5%>%as.matrix()%>%abs()
  net<-d1

  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]


  return(net_m)
}

"matpMerge" <- function(corp,mantelp) {
  d1<-corp
  d5<-mantelp

  add.param<-colnames(d5)[1:length(colnames(d5))]

  add.num<-length(add.param)
  net_m<-d5%>%as.matrix()%>%abs()
  net<-d1

  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]


  return(net_m)
}



"corrposcalcu"<-function (corr,
                          method = c("circle", "square", "ellipse", "number",
                                     "shade", "color", "pie"),
                          type = c("full", "lower", "upper"),
                          col = NULL,
                          col.lim = NULL,
                          is.corr = TRUE,
                          add = FALSE, diag = TRUE,
                          p.mat = NULL)
{
  method = match.arg(method)
  type = match.arg(type)



  if (is.null(col.lim)) {
    if (is.corr) {
      col.lim = c(-1, 1)
    }
    else {
      if (!diag) {
        diag(corr) = NA
      }
      col.lim = c(min(corr, na.rm = TRUE), max(corr, na.rm = TRUE))
    }
  }
  SpecialCorr = 0
  if (is.corr) {
    if (min(corr, na.rm = TRUE) < -1 - .Machine$double.eps^0.75 ||
        max(corr, na.rm = TRUE) > 1 + .Machine$double.eps^0.75) {
      stop("The matrix is not in [-1, 1]!")
    }
    SpecialCorr = 1

  }
  intercept = 0
  zoom = 1
  if (!is.corr) {
    c_max = max(corr, na.rm = TRUE)
    c_min = min(corr, na.rm = TRUE)
    if ((col.lim[1] > c_min) | (col.lim[2] < c_max)) {
      stop("Wrong color: matrix should be in col.lim interval!")
    }
    if (diff(col.lim)/(c_max - c_min) > 2) {
      warning("col.lim interval too wide, please set a suitable value")
    }
    if (c_max <= 0 | c_min >= 0) {
      intercept = -col.lim[1]
      zoom = 1/(diff(col.lim))
      if (col.lim[1] * col.lim[2] < 0) {
        warning("col.lim interval not suitable to the matrix")
      }
    }
    else {
      stopifnot(c_max * c_min < 0)
      stopifnot(c_min < 0 && c_max > 0)
      intercept = 0
      zoom = 1/max(abs(col.lim))
      SpecialCorr = 1
    }
    corr = (intercept + corr) * zoom
  }
  col.lim2 = (intercept + col.lim) * zoom
  int = intercept * zoom

  n = nrow(corr)
  m = ncol(corr)
  min.nm = min(n, m)
  ord = 1:min.nm
  if (is.null(rownames(corr))) {
    rownames(corr) = 1:n
  }
  if (is.null(colnames(corr))) {
    colnames(corr) = 1:m
  }
  apply_mat_filter = function(mat) {
    x = matrix(1:n * m, nrow = n, ncol = m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) <
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) = Inf
    }
    return(mat)
  }
  getPos.Dat = function(mat) {
    tmp = apply_mat_filter(mat)
    Dat = tmp[is.finite(tmp)]
    ind = which(is.finite(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    PosName = ind
    PosName[, 1] = colnames(mat)[ind[, 2]]
    PosName[, 2] = rownames(mat)[ind[, 1]]
    return(list(Pos, Dat, PosName))
  }


  Pos = getPos.Dat(corr)[[1]]
  PosName = getPos.Dat(corr)[[3]]
  DAT = getPos.Dat(corr)[[2]]
  pNew = getPos.Dat(p.mat)[[2]]

  corrPos = data.frame(PosName, Pos, DAT)
  colnames(corrPos) = c("xName", "yName", "x", "y", "corr")
  if (!is.null(p.mat)) {
    corrPos = cbind(corrPos, pNew)
    colnames(corrPos)[6] = c("p.value")
  }
  corrPos = corrPos[order(corrPos[, 3], -corrPos[, 4]), ]
  rownames(corrPos) = NULL
  res = list(corrPos = corrPos)
}



"plot_mantel_1hierarchy" <- function(data12,
                                     data_node,
                                     mantel.node.size=10,
                                     mantel.label.size=1,
                                     edge.size.max.r_thre =0.6,
                                     edge.size.min.r_thre = 0.4,
                                     edge.size.max =10,
                                     edge.size.center= 5,
                                     edge.size.min = 1,
                                     edge.color.max.p_thre =0.05,
                                     edge.color.center.p_thre = 0.01,
                                     edge.color.min.p_thre = 0.001,
                                     edge.color.ultra ="grey90",
                                     edge.color.max ="red",
                                     edge.color.center="blue",
                                     edge.color.min = "darkgreen",
                                     note.text = c("Phenotype", "Mantel", "Phenotype"),
                                     note.text.color = c("#c51b7d", "#2f6661", "#c51b7d"),
                                     arrow.color.left ="#c51b7d",
                                     arrow.color.right = "#2f6661",
                                     label.cex=1,
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = 0.5,
                                     arrow.size = 1.5,
                                     pbreaks = c(-Inf, 0.01, 0.05, Inf),
                                     plabels = c("< 0.01", "0.01 - 0.05", ">= 0.05"),
                                     rbreaks = c(-Inf,0, 0.2, 0.4, Inf),
                                     rlabels = c("< 0","0 - 0.2", "0.2 - 0.4", ">= 0.4"),
                                     size.r=c(0.5,1,2,4),
                                     color.p=random2color(5),
                                     node.label.color=NULL,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "gygn")
) {
  library(tidyverse)

  corr<-data12$occor.r
  d2<-data12$mantel
  corp<-data12$occor.p

  d3<-subset(d2,select=c(1:3))
  d4<-subset(d2,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
  d5<-tibble::column_to_rownames(d5,var = "env")

  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")


  mantelr<-d5%>%as.matrix()%>%abs()
  mantelp<-d6%>%as.matrix()%>%abs()

  netc<-matrMerge(corr,mantelr)
  netp<-matpMerge(corp,mantelp)
  net<-corr%>%abs()
  rownames(net)<-colnames(net)

  data<-data_node
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }

  vertex.weight<-sum%>%as.numeric()
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }

  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  net2 <- netc
  m <- length(vertex.receiver)
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  net3 <- rbind(cbind(matrix(0, m1, m1), net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0,n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[vertex.receiver],
                       row.names(net)[setdiff(1:m1,vertex.receiver)],
                       colnames(net2)[1:n1])
  colnames(net3) <- row.names(net3)


  color.use3 <- c(color.use[vertex.receiver],
                  color.use[setdiff(1:m1,vertex.receiver)],
                  rep("#ffffff", n1))

  color.use3.frame <- c(color.use[vertex.receiver],
                        color.use[setdiff(1:m1,vertex.receiver)],
                        color.use[seq(m1+1,m1+n1)])

  colnames(net3) <- row.names(net3)

  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", n1))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m, 1] <- 0
  coords[(m + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
  coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                            m - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1-1))

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m)],vertex.weight[seq(m+1,m1)],
                         rep(mantel.node.size,n1))


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]

  igraph::V(g)$label.cex <- label.cex
  igraph::V(g)$label.cex <-c(rep(label.cex,m1),rep(mantel.label.size,n1))

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }

  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size


  netp<-t(netp)
  rd = cut(netc, breaks = rbreaks ,
           labels = rlabels)

  r_num<-length(levels(rd))
  size.use2<-size.r
  names(size.use2) <-levels(rd)
  size.use2<-size.use2%>%data.frame()
  size.use2<-rownames_to_column(size.use2,var = "color")

  rd<-data.frame(rd = cut(netc, breaks = rbreaks,
                          labels = rlabels))

  rd$rd<-as.character(rd$rd)
  rd$color<-rd$rd
  for (k in 1:r_num) {
    rd$color[which(match(rd$rd,size.use2$color[k])==1)]<-size.use2$.[k]
  }
  rd<-rd$color
  rrd<-matrix(rd,nrow(netc),ncol(netc))
  colnames(rrd)<-colnames(netc)
  rownames(rrd)<-rownames(netc)
  E(g)$weight<-rrd%>%t()
  E(g)$width <- E(g)$weight


  pd = cut(netp, breaks = pbreaks,
           labels = plabels)

  p_num<-length(levels(pd))
  color.use2<-color.p
  names(color.use2) <-levels(pd)
  color.use2<-color.use2%>%data.frame()
  color.use2<-rownames_to_column(color.use2,var = "color")

  pd<-data.frame(pd = cut(netp, breaks = pbreaks,
                          labels = plabels))
  pd$pd<-as.character(pd$pd)
  pd$color<-pd$pd
  for (k in 1:p_num) {
    pd$color[which(match(pd$pd,color.use2$color[k])==1)]<-color.use2$.[k]
  }
  pd<-pd$color
  ppd<-matrix(pd,nrow(netp),ncol(netp))
  colnames(ppd)<-colnames(netp)
  rownames(ppd)<-rownames(netp)

  E(g)$color<-ppd



  label.dist<-NULL
  label.locs<-NULL

  label.dist <- c(rep(space.h * 2.8, m),
                  rep(space.h * 2.8, m1 - m),
                  rep(2.8, n1))
  label.locs <- c(rep(-pi, m),
                  rep(0, m1 - m),
                  rep(pi/2, n1-1),
                  rep(-pi/2, 1))

  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)

  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
       vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")

  text(text.pos, note.text, cex = 0.8,
       col = note.text.color,family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, -0.04,
                  space.v - space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v - space.v/4, 0.04,
                  space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = arrow.color.left, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3],
                arrow.pos2[4], col = arrow.color.right, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence communication network"), cex = 1,family="serif")
  gg <- recordPlot()
  return(gg)
}


"plot_mantel_2hierarchy" <- function(data_test,
                                     data_node, ###计算节点大小
                                     margin=0.2,
                                     mantel.node.size=10,
                                     mantel.label.size=0.5,
                                     label.cex=1,
                                     edge.size.max.r_thre = 0.6,
                                     edge.size.min.r_thre = 0.4,
                                     edge.size.max = 8,
                                     edge.size.center= 4,
                                     edge.size.min = 2,
                                     note.text = c("Phenotype", "Mantel", "Phenotype"),
                                     note.text.color = c("#c51b7d", "#2f6661", "#c51b7d"),
                                     arrow.color.left ="#c51b7d",
                                     arrow.color.right = "#2f6661",
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = 0,
                                     arrow.size = 0,
                                     pbreaks = c(-Inf,0.001, 0.01, 0.05, Inf),
                                     plabels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"),
                                     color.p=c("grey90","red","blue","yellow"),
                                     rbreaks = c(-Inf,0, 0.4, 0.8, Inf),
                                     rlabels = c("< 0","0 - 0.2", "0.2 - 0.8", ">= 0.8"),
                                     size.r=c(0.5,2,5,10),
                                     node.label.color=NULL,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "gygn"))
{
  library(tidyverse)

  corr<-data12$occor.r
  d2<-data12$mantel
  corp<-data12$occor.p

  d3<-subset(d2,select=c(1:3))
  d4<-subset(d2,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)
  d5<-tibble::column_to_rownames(d5,var = "env")

  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")


  mantelr<-d5%>%as.matrix()%>%abs()
  mantelp<-d6%>%as.matrix()%>%abs()

  netc<-matrMerge(corr,mantelr)
  netp<-matpMerge(corp,mantelp)
  net<-corr%>%abs()
  rownames(net)<-colnames(net)

  data<-data_node
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }

  vertex.weight<-sum%>%as.numeric()
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }

  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  net2 <- netc
  m <- length(vertex.receiver)
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  net3 <- rbind(cbind(matrix(0, m1, m1),net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0, n1, m1 + n1))


  row.names(net3) <- c(row.names(corr),
                       colnames(net2)[1:n1])
  colnames(net3) <- row.names(net3)


  color.use3 <- c(
    color.use[1:nrow(corr)],
    rep("#ffffff", n1))

  color.use3.frame <- c(
    color.use[1:nrow(corr)],
    color.use[seq(m1+1,m1+n1)])

  colnames(net3) <- row.names(net3)

  shape <- c(
    rep("circle", m1),
    rep("circle", n1))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m1, 1] <- 0
  coords[(m1 + 1):nrow(net3), 1] <- space.h

  coords[1:m1, 2] <- seq(space.v, 0, by = -space.v/(m1 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v*9/10, space.v*1/10, length.out=n1)

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m1)],
                         rep(5,n1))


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]


  igraph::V(g)$label.cex <-c(rep(1,m1),
                             rep(1,n1))

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }

  E(g)$arrow.width <- 0
  E(g)$arrow.size <- 0



  netp<-t(netp)

  rd = cut(netc, breaks = rbreaks ,
           labels = rlabels)

  r_num<-length(levels(rd))
  size.use2<-size.r
  names(size.use2) <-levels(rd)
  size.use2<-size.use2%>%data.frame()
  size.use2<-rownames_to_column(size.use2,var = "color")

  rd<-data.frame(rd = cut(netc, breaks = rbreaks,
                          labels = rlabels))

  rd$rd<-as.character(rd$rd)
  rd$color<-rd$rd
  for (k in 1:r_num) {
    rd$color[which(match(rd$rd,size.use2$color[k])==1)]<-size.use2$.[k]
  }
  rd<-rd$color
  rrd<-matrix(rd,nrow(netc),ncol(netc))
  colnames(rrd)<-colnames(netc)
  rownames(rrd)<-rownames(netc)
  E(g)$weight<-rrd%>%t()
  E(g)$width <- E(g)$weight





  pd = cut(netp, breaks = pbreaks,
           labels = plabels)

  p_num<-length(levels(pd))
  color.use2<-color.p
  names(color.use2) <-levels(pd)
  color.use2<-color.use2%>%data.frame()
  color.use2<-rownames_to_column(color.use2,var = "color")

  pd<-data.frame(pd = cut(netp, breaks = pbreaks,
                          labels = plabels))
  pd$pd<-as.character(pd$pd)
  pd$color<-pd$pd
  for (k in 1:p_num) {
    pd$color[which(match(pd$pd,color.use2$color[k])==1)]<-color.use2$.[k]
  }
  pd<-pd$color
  ppd<-matrix(pd,nrow(netp),ncol(netp))
  colnames(ppd)<-colnames(netp)
  rownames(ppd)<-rownames(netp)

  E(g)$color<-ppd



  label.dist<-NULL
  label.locs<-NULL

  label.dist <- c(rep(space.h * 2.8, m1),
                  rep(2.8, n1))
  label.locs <- c(
    rep(-pi, m1),
    rep(0, n1-1),
    rep(0, 1))

  text.pos <- cbind(c(-space.h/1.6,  space.h/1.6),
                    space.v - space.v/7)

  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = margin, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
       vertex.label.degree = label.locs, vertex.label.dist = label.dist,
       vertex.label.family = "serif")

  text(text.pos, note.text, cex = 0.8,
       col = note.text.color,family="serif")
  arrow.pos1 <- c(-space.h/1.5, space.v - space.v/4, space.h/1.5, space.v - space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3],
                arrow.pos1[4], col = arrow.color.left, arr.lwd = 1e-04, arr.length = 0.2,
                lwd = 0.8, arr.type = "triangle")

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network"), cex = 1,family="serif")

  gg <- recordPlot()
  return(gg)

}


"cor_mantel" <- function(edmat,ctmat,
                         method="pearson",
                         permutations=9999,
                         filter=TRUE,ed.p=0.05,ed.r=0.4,
                         mantel_select=list(Spec01 = 1:7)) {

  occor.r<-WGCNA::corAndPvalue(ctmat,method = 'pearson')$cor
  occor.p<-WGCNA::corAndPvalue(ctmat,method = 'pearson')$p
  diag(occor.r) <- 0

  mtadj<-multtest::mt.rawp2adjp(unlist(occor.p),proc='Bonferroni')
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  adjp<-adpcor
  adjp<-matrix(adjp,dim(occor.p)[2])
  row.names(adjp)<-row.names(occor.p)
  colnames(adjp)<-colnames(occor.p)
  occor.p<-adjp

  if (filter) occor.r[occor.p>ed.p|abs(occor.r)<ed.r] = 0

  mantel <- mantel_calc( edmat, ctmat,
                         method=method,
                         permutations=permutations,
                         select_col =mantel_select)

  return(list(occor.r=occor.r,mantel=mantel,occor.p=occor.p))
}

"padj"<-function (rawp, proc = c("Bonferroni", "Holm", "Hochberg", "SidakSS",
                                 "SidakSD", "BH", "BY", "ABH", "TSBH"),
                  alpha = 0.05, na.rm = FALSE)
{
  m <- length(rawp)
  if (na.rm) {
    mgood <- sum(!is.na(rawp))
  } else {
    mgood <- m
  }
  n <- length(proc)
  a <- length(alpha)
  index <- order(rawp)
  h0.ABH <- NULL
  h0.TSBH <- NULL
  spval <- rawp[index]
  adjp <- matrix(0, m, n + 1)
  dimnames(adjp) <- list(NULL, c("rawp", proc))
  adjp[, 1] <- spval
  if (is.element("TSBH", proc)) {
    TS.spot <- which(proc == "TSBH")
    TSBHs <- paste("TSBH", alpha, sep = "_")
    newprocs <- append(proc, TSBHs, after = TS.spot)
    newprocs <- newprocs[newprocs != "TSBH"]
    adjp <- matrix(0, m, n + a)
    dimnames(adjp) <- list(NULL, c("rawp", newprocs))
    adjp[, 1] <- spval
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    h0.TSBH <- rep(0, length(alpha))
    names(h0.TSBH) <- paste("h0.TSBH", alpha, sep = "_")
    for (i in 1:length(alpha)) {
      h0.TSBH[i] <- mgood - sum(tmp < alpha[i]/(1 + alpha[i]),
                                na.rm = TRUE)
      adjp[, TS.spot + i] <- tmp * h0.TSBH[i]/mgood
    }
  }
  if (is.element("Bonferroni", proc)) {
    tmp <- mgood * spval
    tmp[tmp > 1] <- 1
    adjp[, "Bonferroni"] <- tmp
  }
  if (is.element("Holm", proc)) {
    tmp <- spval
    tmp[1] <- min(mgood * spval[1], 1)
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], min((mgood -
                                                    i + 1) * spval[i], 1))
    adjp[, "Holm"] <- tmp
  }
  if (is.element("Hochberg", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood - i + 1) *
                                      spval[i], 1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "Hochberg"] <- tmp
  }
  if (is.element("SidakSS", proc))
    adjp[, "SidakSS"] <- 1 - (1 - spval)^mgood
  if (is.element("SidakSD", proc)) {
    tmp <- spval
    tmp[1] <- 1 - (1 - spval[1])^mgood
    for (i in 2:m) tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i])^(mgood -
                                                                   i + 1))
    adjp[, "SidakSD"] <- tmp
  }
  if (is.element("BH", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "BH"] <- tmp
  }
  if (is.element("BY", proc)) {
    tmp <- spval
    a <- sum(1/(1:mgood))
    tmp[m] <- min(a * spval[m], 1)
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood * a/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "BY"] <- tmp
  }
  if (is.element("ABH", proc)) {
    tmp <- spval
    h0.m <- rep(0, mgood)
    for (k in 1:mgood) {
      h0.m[k] <- (mgood + 1 - k)/(1 - spval[k])
    }
    grab <- min(which(diff(h0.m, na.rm = TRUE) > 0), na.rm = TRUE)
    h0.ABH <- ceiling(min(h0.m[grab], mgood))
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((mgood/i) * spval[i],
                                    1, na.rm = TRUE), na.rm = TRUE)
      if (is.na(spval[i]))
        tmp[i] <- NA
    }
    adjp[, "ABH"] <- tmp * h0.ABH/mgood
  }
  list(adjp = adjp, index = index, h0.ABH = h0.ABH[1], h0.TSBH = h0.TSBH[1:length(alpha)])
}



"mantel_calc" <- function(spec,param,
                          method="pearson",
                          permutations=9999,
                          select_col=list( firm = 2:2,cya= 1:1,

                                           pro=3:3)) {

  cat("\n","\n",'method: "pearson", "spearman" or "kendall"')
  datad<-data.frame()
  mantel.result<-data.frame()
  for (k in 1:length(select_col)) {
    select_col1<-select_col[[k]]
    spec.name<-names(select_col)[[k]]
    spec[,select_col1][which(spec[,select_col1]==0)]<-1
    veg.dist <- vegan::vegdist(spec[,select_col1]) # Bray-Curtis

    for (t in 1:ncol(param)) {
      env.dist <- vegan::vegdist(scale(param[,t]), "euclid")
      mantel.res<-vegan::mantel(veg.dist, env.dist,
                                permutations=permutations,
                                method=method)
      mantel.res.r<-mantel.res$statistic
      mantel.res.p<-mantel.res$signif
      dad<-data.frame(spec=spec.name,
                      env=colnames(param)[t],
                      r=mantel.res.r,
                      p=mantel.res.p)

      {
        datad<-rbind(datad,dad)
      }

    }

  }
  datad$spec<-factor(datad$spec,levels = names(select_col))

  return(datad)
}



"plot_hireedd" <- function(data12,
                           color.use,
                           vertex.receiver,
                           node.label.color=NULL,
                           arrow.width = arrow.width,
                           arrow.size = arrow.size) {
  par(mfrow = c(1,2), xpd=TRUE,family="serif")
  matchat_net_1hierarchy(data12,
                         vertex.weight.keep = 5,
                         label.cex=1,
                         space.h=1.6,
                         space.v=1.5,
                         arrow.width = arrow.width,
                         arrow.size = arrow.size,
                         node.label.color=node.label.color,
                         vertex.receiver = vertex.receiver,
                         color.use=color.use)

  matchat_net_2hierarchy(data12,
                         vertex.weight.keep = 5,
                         label.cex=1,
                         space.h=1.6,
                         space.v=1.5,
                         arrow.width = arrow.width,
                         arrow.size = arrow.size,
                         node.label.color=node.label.color,
                         vertex.receiver = vertex.receiver,
                         color.use=color.use)


}

"matchat_net_2hierarchy" <- function(data12,
                                     vertex.weight.keep = 5,
                                     label.cex=1,
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = arrow.width,
                                     arrow.size = arrow.size,
                                     node.label.color=node.label.color,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "pkgn"))
{




  d1<-data12$occor.r
  d2<-data12$mantel

  d3<-subset(d2,select=c(1:3))
  d4<-subset(d2,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)

  d5<-tibble::column_to_rownames(d5,var = "env")
  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")

  ###添加的参数的id
  add.param<-colnames(d5)[1:length(colnames(d5))]
  ####增加的数量
  add.num<-length(add.param)

  net_m<-d5%>%as.matrix()%>%abs()

  net<-d1%>%abs()

  ###添加的矩阵排序
  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]

  rownames(net)<-colnames(net)

  data<-mtcars[1:24,]
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }


  vertex.weight<-sum%>%as.numeric()

  color.use<-color.use
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }


  vertex.receiver = setdiff(1:nrow(net),
                            vertex.receiver)
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max


  m <- length(vertex.receiver)
  m0 <- nrow(net) - length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net), vertex.receiver),
                   vertex.receiver)
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  ###添加的矩阵排序
  add.order<-match(rownames(net2),rownames(net_m))
  net_m<-net_m[add.order,]

  net3 <- rbind(cbind(matrix(0, m1, m1), net2,net_m),
                matrix(0,n1+add.num, m1 + n1+add.num))


  row.names(net3) <- c(row.names(net)[setdiff(1:m1, vertex.receiver)],
                       row.names(net)[vertex.receiver],
                       rep("", m+add.num))


  colnames(net3) <- row.names(net3)



  color.use3 <- c(color.use[setdiff(1:m1, vertex.receiver)],
                  color.use[vertex.receiver],
                  rep("#FFFFFF",length(vertex.receiver)),
                  color.use[seq(m1+1,m1+add.num)])

  color.use3.frame <- c(color.use[setdiff(1:m1, vertex.receiver)],
                        color.use[vertex.receiver],
                        color.use[vertex.receiver],
                        color.use[seq(m1+1,m1+add.num)])


  row.names(net3)[seq(m1+n1+1,m1+n1+add.num)] <- add.param
  colnames(net3) <- row.names(net3)


  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", m+add.num))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)
  space.h=space.h
  space.v=space.v


  edge.start <- ends(g, es = igraph::E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0, 1] <- 0
  coords[(m0 + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m0, 2] <- seq(space.v, 0, by = -space.v/(m0 - 1))
  coords[(m0 + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -m0 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1 +add.num-1))

  coords_scale <- coords

  igraph::V(g)$size <- c(vertex.weight[seq(1,m0)],
                         vertex.weight[seq(m0+1,m1)],
                         vertex.weight[seq(m0+1,m1)],
                         rep(5,add.num))


  ###固定大小
  vertex.weight.keep = vertex.weight.keep
  igraph::V(g)$size[seq(m1+n1+1,m1+n1+add.num)]<-rep(vertex.weight.keep,add.num)
  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.cex <- label.cex

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }

  edge.weight.max <- max(igraph::E(g)$weight)

  E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * 8


  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,1]], 0.6)

  label.dist<-NULL
  label.locs<-NULL


  label.dist <- c(rep(space.h * 2.8, m0),
                  rep(space.h *2.8, m1 - m0),
                  rep(0, m1 - m0),
                  rep(space.v, nrow(net3) - m1-m))

  label.locs <- c(rep(-pi, m0),
                  rep(0, m1 - m0),
                  rep(0, m1 - m0),
                  rep(pi/2,nrow(net3) - m1-m))

  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)

  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle, parameters = list(vertex.frame.color = 1,
                                                              vertex.frame.width = 1))


  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) -m1)),
       vertex.label.degree = label.locs,
       vertex.label.dist = label.dist,
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

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network",
                                          " signaling network"), cex = 1,family="serif")



}


"matchat_net_1hierarchy" <- function(data12,
                                     vertex.weight.keep = 5,
                                     label.cex=1,
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = arrow.width,
                                     arrow.size = arrow.size,
                                     node.label.color=node.label.color,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "pkgn"))
{




  d1<-data12$occor.r
  d2<-data12$mantel

  d3<-subset(d2,select=c(1:3))
  d4<-subset(d2,select=c(1:2,4))

  d5<-tidyr::spread(d3,key = "spec",value = r,drop=FALSE)

  d5<-tibble::column_to_rownames(d5,var = "env")
  d6<-tidyr::spread(d4,key = "spec",value = p,drop=FALSE)
  d6<-tibble::column_to_rownames(d6,var = "env")

  ###添加的参数的id
  add.param<-colnames(d5)[1:length(colnames(d5))]
  ####增加的数量
  add.num<-length(add.param)




  net_m<-d5%>%as.matrix()%>%abs()

  net<-d1%>%abs()

  ###添加的矩阵排序
  add.order<-match(rownames(net),rownames(net_m))
  net_m<-net_m[add.order,]

  rownames(net)<-colnames(net)

  data<-mtcars[1:24,]
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }


  vertex.weight<-sum%>%as.numeric()

  color.use<-color.use
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }


  vertex.receiver = vertex.receiver
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  m <- length(vertex.receiver)

  net2 <- net
  reorder.row <- c(vertex.receiver,
                   setdiff(1:nrow(net),
                           vertex.receiver))
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  net3 <- rbind(cbind(matrix(0, m1, m1), net2,net_m), ###原矩阵行+选择的行+添加的mantel行
                matrix(0,n1+add.num, m1 + n1+add.num))
  row.names(net3) <- c(row.names(net)[vertex.receiver],
                       row.names(net)[setdiff(1:m1,vertex.receiver)],
                       rep("", m+add.num))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver],
                  color.use[setdiff(1:m1,vertex.receiver)],
                  rep("#FFFFFF", length(vertex.receiver)),
                  color.use[seq(m1+1,m1+add.num)])
  color.use3.frame <- c(color.use[vertex.receiver],
                        color.use[setdiff(1:m1,vertex.receiver)],
                        color.use[vertex.receiver],
                        color.use[seq(m1+1,m1+add.num)])


  row.names(net3)[seq(m1+n1+1,m1+n1+add.num)] <- add.param
  colnames(net3) <- row.names(net3)

  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", m+add.num))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  space.h=space.h
  space.v=space.v


  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m, 1] <- 0
  coords[(m + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
  coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                            m - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1+add.num-1 ))

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m)],vertex.weight[seq(m+1,m1)],
                         vertex.weight[seq(1,m)],rep(5,add.num))


  ###固定大小
  vertex.weight.keep = vertex.weight.keep
  igraph::V(g)$size[seq(m1+n1+1,m1+n1+add.num)]<-rep(vertex.weight.keep,add.num)
  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]

  igraph::V(g)$label.cex <- label.cex

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }


  edge.weight.max <- max(igraph::E(g)$weight)

  E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * 8

  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size
  E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,
                                                          1]], 0.6)
  label.dist<-NULL
  label.locs<-NULL
  label.dist <- c(rep(space.h * 2.8, m),
                  rep(space.h *2.8, m1 - m),
                  rep(0, m),
                  rep(space.v, nrow(net3) - m1-m))
  label.locs <- c(rep(-pi, m),
                  rep(0, m1 - m),
                  rep(0, m),
                  rep(pi/2, nrow(net3) -m1-m)) ###字默认在节点的右边，设置字体旋转角度
  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)



  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
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

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network",
                                          " signaling network"), cex = 1,family="serif")

  gg <- recordPlot()
  return(gg)


}


"plot_mantel_4hierarchy" <- function(data12,
                                     label.cex=1,
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = arrow.width,
                                     arrow.size = arrow.size,
                                     node.label.color=node.label.color,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "pkgn"))
{

  d1<-data12$occor.r

  net<-d1%>%abs()

  rownames(net)<-colnames(net)
  data<-mtcars[1:24,]
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }


  vertex.weight<-sum%>%as.numeric()

  color.use<-color.use
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }

  vertex.receiver = setdiff(1:nrow(net),
                            vertex.receiver)
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  m <- length(vertex.receiver)

  m0 <- nrow(net) - length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net), vertex.receiver),
                   vertex.receiver)
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)


  net3 <- rbind(cbind(matrix(0, m1, m1), net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0,n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1, vertex.receiver)],
                       row.names(net)[vertex.receiver],
                       rep("", m))
  colnames(net3) <- row.names(net3)

  color.use3 <- c(color.use[setdiff(1:m1, vertex.receiver)],
                  color.use[vertex.receiver],
                  rep("#FFFFFF", length(vertex.receiver)))

  color.use3.frame <- c(color.use[setdiff(1:m1, vertex.receiver)],
                        color.use[vertex.receiver],
                        color.use[vertex.receiver])




  colnames(net3) <- row.names(net3)

  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", m))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  space.h=space.h
  space.v=space.v

  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0, 1] <- 0
  coords[(m0 + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m0, 2] <- seq(space.v, 0, by = -space.v/(m0 - 1))
  coords[(m0 + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                             m0 - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1-1 ))


  coords_scale <- coords



  igraph::V(g)$size <- c(vertex.weight[seq(1,m0)],
                         vertex.weight[seq(m0+1,m1)],
                         vertex.weight[seq(m0+1,m1)])


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]

  igraph::V(g)$label.cex <- label.cex

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }


  edge.weight.max <- max(igraph::E(g)$weight)

  E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * 8

  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size
  E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,
                                                          1]], 0.6)
  label.dist<-NULL
  label.locs<-NULL


  label.dist <- c(rep(space.h * 2.8, m0),
                  rep(space.h * 2.8, m1 - m0),
                  rep(0, nrow(net3) - m1))
  label.locs <- c(rep(-pi, m0),
                  rep(0, m1 - m0),
                  rep(-pi, nrow(net3) - m1))

  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)



  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
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

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network",
                                          " signaling network"), cex = 1,family="serif")


}


"plot_mantel_3hierarchy" <- function(data12,
                                     label.cex=1,
                                     space.h=1.6,
                                     space.v=1.5,
                                     arrow.width = arrow.width,
                                     arrow.size = arrow.size,
                                     node.label.color=node.label.color,
                                     vertex.receiver = seq(1,6),
                                     color.use=colorCustom(50,pal = "pkgn"))
{

  d1<-data12$occor.r

  net<-d1%>%abs()

  rownames(net)<-colnames(net)

  data<-mtcars[1:24,]
  sum<-list()
  mean<-list()
  colnum<-length(colnames(data))
  for (k in 1:colnum) {
    data[,k]<-microchat::normalize(data[,k])$x
    sum[k]<-sum(data[,k])
    mean[k]<-mean(data[,k])
  }


  vertex.weight<-sum%>%as.numeric()

  color.use<-color.use
  vertex.weight.max <- max(vertex.weight)

  if (length(unique(vertex.weight)) == 1) {
    vertex.size.max <- 5
  }else {
    vertex.size.max <- vertex.weight%>%max()
    vertex.weight.max<- vertex.size.max
  }


  vertex.receiver = vertex.receiver
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max

  m <- length(vertex.receiver)

  net2 <- net
  reorder.row <- c(vertex.receiver,
                   setdiff(1:nrow(net),
                           vertex.receiver))
  net2 <- net2[reorder.row, vertex.receiver]
  m1 <- nrow(net2)
  n1 <- ncol(net2)

  net3 <- rbind(cbind(matrix(0, m1, m1), net2), ###原矩阵行+选择的行+添加的mantel行
                matrix(0,n1, m1 + n1))
  row.names(net3) <- c(row.names(net)[vertex.receiver],
                       row.names(net)[setdiff(1:m1,vertex.receiver)],
                       rep("", m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver],
                  color.use[setdiff(1:m1,vertex.receiver)],
                  rep("#FFFFFF", length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver],
                        color.use[setdiff(1:m1,vertex.receiver)],
                        color.use[vertex.receiver])

  colnames(net3) <- row.names(net3)

  shape <- c(rep("circle", m), rep("circle", m1 - m),
             rep("circle", m))

  g <- graph_from_adjacency_matrix(net3, mode = "directed",
                                   weighted = T)

  space.h=space.h
  space.v=space.v


  edge.start <- ends(g, es = E(g), names = FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m, 1] <- 0
  coords[(m + 1):m1, 1] <- space.h
  coords[(m1 + 1):nrow(net3), 1] <- space.h/2
  coords[1:m, 2] <- seq(space.v, 0, by = -space.v/(m - 1))
  coords[(m + 1):m1, 2] <- seq(space.v, 0, by = -space.v/(m1 -
                                                            m - 1))
  coords[(m1 + 1):nrow(net3), 2] <- seq(space.v, 0, by = -space.v/(n1-1 ))

  coords_scale <- coords
  igraph::V(g)$size <- c(vertex.weight[seq(1,m)],vertex.weight[seq(m+1,m1)],
                         vertex.weight[seq(1,m)])


  igraph::V(g)$color <- color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]

  igraph::V(g)$label.cex <- label.cex

  if (is.null(node.label.color)) {
    igraph::V(g)$label.color <- color.use3.frame[igraph::V(g)]
  } else {
    igraph::V(g)$label.color <- node.label.color
  }


  edge.weight.max <- max(igraph::E(g)$weight)

  E(g)$width <- 0.3 + E(g)$weight/edge.weight.max * 8

  E(g)$arrow.width <- arrow.width
  E(g)$arrow.size <- arrow.size
  E(g)$color <- adjustcolor(igraph::V(g)$color[edge.start[,
                                                          1]], 0.6)
  label.dist<-NULL
  label.locs<-NULL


  label.dist <- c(rep(space.h * 2.8, m),
                  rep(space.h * 2.8, m1 - m),
                  rep(0, nrow(net3) - m1))
  label.locs <- c(rep(-pi, m),
                  rep(0, m1 - m),
                  rep(-pi, nrow(net3) - m1))

  text.pos <- cbind(c(-space.h/1.6, 0, space.h/1.6),
                    space.v - space.v/7)



  igraph::add.vertex.shape("fcircle", clip = igraph::igraph.shape.noclip,
                           plot = mycircle,
                           parameters = list(vertex.frame.color = 1,
                                             vertex.frame.width = 1))

  plot(g, edge.curved = 0, layout = coords_scale,
       margin = 0.2, rescale = T, vertex.shape = "fcircle",
       vertex.frame.width = c(rep(1, m1), rep(2, nrow(net3) - m1)),
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

  title.pos = c(0, space.v)
  text(title.pos[1], title.pos[2], paste0("Micro-occurrence network",
                                          " signaling network"), cex = 1,family="serif")


}
