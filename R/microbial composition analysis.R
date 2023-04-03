"plotMicrocomBar" <- function(microchatcomobj,
                              xlabname=NULL,
                              data_type=c("absolute","relative"),
                              sample_type=c("allsamples","groups"),
                              chicklet=FALSE,
                              color_taxa=NULL,
                              color_group=NULL,
                              barplot.barwidth=0.6,
                              aspect.ratio=0.8,
                              color_background=NULL,
                              color_border=NULL,
                              export_path="microbial composition") {
  data_type=match.arg(data_type)
  sample_type=match.arg(sample_type)

  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
    split_otu3<-microchatcomobj$long_group_abs
  }

  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
    split_otu3<-microchatcomobj$long_group_rel
  }

  if (!chicklet) {
    if (sample_type=="groups") {

      p1<-ggplot(data = split_otu3 ,
                 aes(x = variable, y = value,
                     fill = tax)) +
        geom_col(position = 'fill', width = barplot.barwidth)

      if (!is.null(color_taxa)) {
        color_value<-color_taxa
        color_taxa<-color_value[1:(taxa_num+1)]
        names(color_taxa)<-unique(split_otu3$tax)

        message("use self-selected color")
        p1<-p1+
          scale_color_manual(name=taxa,values=color_taxa)+
          scale_fill_manual(name=taxa,values=color_taxa)

      } else {

        color_value<-c(randomcolor((length(unique(split_otu3$tax))-1)),"grey90")
        names(color_value)<-unique(split_otu3$tax)

        p1<-p1+scale_fill_manual(name=taxa,values = color_value)
      }

      p1<-p1+
        labs(x = '', y = 'Relative Abundance(%)') +
        theme(aspect.ratio = aspect.ratio,
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              panel.background = element_rect(fill = "grey90"),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))


    }
    if (sample_type=="allsamples") {
      genus_topt<-tibble::rownames_to_column(genus_top,var = "tax")
      genus_topt<-reshape2::melt(genus_topt)


      p1<-ggplot(data = genus_topt ,
                 aes(x = variable, y = value,
                     fill = tax)) +
        geom_col(position = 'fill', width = barplot.barwidth)



      if (!is.null(color_taxa)) {
        color_value<-color_taxa
        color_taxa<-color_value[1:(taxa_num+1)]
        names(color_taxa)<-unique(genus_topt$tax)

        message("use self-selected color")
        p1<-p1+
          scale_color_manual(name=taxa,values=color_taxa,
                             aesthetics = "colour")+
          scale_fill_manual(name=taxa,values=color_taxa,
                            aesthetics = "fill")

      } else {

        color_value<-c(randomcolor((length(unique(genus_topt$tax))-1)),"grey90")
        names(color_value)<-unique(genus_topt$tax)

        p1<-p1+scale_fill_manual(name=taxa,values = color_value)


      }
      p1<-p1+
        labs(x = '', y = 'Relative Abundance(%)') +
        theme(aspect.ratio = aspect.ratio,
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              panel.background = element_rect(fill = "grey90"),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))


    }
  } else {
    if (sample_type=="groups") {

      p1<-ggplot(data = split_otu3 ,
                 aes(x = variable, y = value,
                     fill = tax)) +
        ggchicklet::geom_chicklet(width = barplot.barwidth)

      if (!is.null(color_taxa)) {
        color_value<-color_taxa
        color_taxa<-color_value[1:(taxa_num+1)]
        names(color_taxa)<-unique(split_otu3$tax)

        message("use self-selected color")
        p1<-p1+
          scale_color_manual(name=taxa,values=color_taxa)+
          scale_fill_manual(name=taxa,values=color_taxa)

      } else {

        color_value<-c(randomcolor((length(unique(split_otu3$tax))-1)),"grey90")
        names(color_value)<-unique(split_otu3$tax)

        p1<-p1+scale_fill_manual(name=taxa,values = color_value)
      }

      p1<-p1+
        labs(x = '', y = 'Relative Abundance(%)') +
        theme(aspect.ratio = aspect.ratio,
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              panel.background = element_rect(fill = "grey90"),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))
    }
    if (sample_type=="allsamples") {
      genus_topt<-tibble::rownames_to_column(genus_top,var = "tax")
      genus_topt<-reshape2::melt(genus_topt)


      p1<-ggplot(data = genus_topt ,
                 aes(x = variable, y = value,
                     fill = tax)) +
        ggchicklet::geom_chicklet(width = barplot.barwidth)


      if (!is.null(color_taxa)) {
        color_value<-color_taxa
        color_taxa<-color_value[1:(taxa_num+1)]
        names(color_taxa)<-unique(genus_topt$tax)

        message("use self-selected color")
        p1<-p1+
          scale_color_manual(name=taxa,values=color_taxa,
                             aesthetics = "colour")+
          scale_fill_manual(name=taxa,values=color_taxa,
                            aesthetics = "fill")

      } else {

        color_value<-c(randomcolor((length(unique(genus_topt$tax))-1)),"grey90")
        names(color_value)<-unique(genus_topt$tax)

        p1<-p1+scale_fill_manual(name=taxa,values = color_value)


      }
      p1<-p1+
        labs(x = '', y = 'Relative Abundance(%)') +
        theme(aspect.ratio = aspect.ratio,
              axis.ticks.length = unit(0.2,"lines"),
              axis.ticks = element_line(color='black'),
              panel.background = element_rect(fill = "grey90"),
              #axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
              axis.text.y=element_text(colour='black',size=10,family = "serif"),
              axis.text.x=element_text(colour = "black",size = 10,
                                       angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))


    }
  }

  if (!is.null(color_background)) {
    p1<-p1+
      theme(panel.background=element_rect(fill = color_background))
  }
  if (is.null(color_border)) {
    p1<-p1+theme(legend.key = element_blank(),
                 text = element_text(family = "serif") )
  } else {
    p1<-p1+theme(panel.border = element_rect(fill=NA,colour = color_border),
                 legend.key = element_blank(),
                 text = element_text(family = "serif") )
  }

  ####change x-axis label
  if (!is.null(xlabname)) {
    orignam<-microchatcomobj$long_group_rel$variable%>%unique()
    if (length(orignam)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
    names(xlabname)<-orignam
    p1<-p1+ scale_x_discrete(labels = xlabname)
  }
  ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_barplot (",chicklet,").pdf",sep = ""),p1)
  cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_barplot"," has been exported to ",export_path,sep = "","\n")

  return(p1)
}


"plotMicrocomPie" <- function(microchatcomobj,
                              xlabname=NULL,
                              data_type=c("absolute","relative"),
                              sample_type=c("allsamples","groups"),
                              color_taxa=NULL,
                              color_group=NULL,
                              color_percent="white",
                              export_path="microbial composition") {
  data_type=match.arg(data_type)
  sample_type=match.arg(sample_type)

  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
    split_otu3<-microchatcomobj$long_group_abs
  }

  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
    split_otu3<-microchatcomobj$long_group_rel
  }

  ncols<-plot_data$group%>%unique()%>%length()
  if (ncols==1) ncols=1
  if (ncols==2) ncols=2
  if (ncols==3) ncols=3
  if (ncols==4) ncols=2
  if (ncols>4 & ncols <=16) ncols=3
  if (ncols>16 & ncols <=25) ncols=5
  oriname<-unique(split_otu3$variable)

  ####change x-axis label
  if (!is.null(xlabname)) {
    if (length(oriname)!=length(xlabname)) stop("The number of xlabname cannot meet the requirement!!!")
  }

  pp<-list()
  for (t in 1:length(oriname)) {
    gname<-oriname[t]
    split_otu4<-split_otu3[which(split_otu3$variable==gname),]
    split_otu4$tax<-factor(split_otu4$tax,levels = unique(split_otu4$tax))
    split_otu4$labelxs<-paste(sprintf("%0.2f",round(split_otu4$value/sum(split_otu4$value),4)*100),"%",sep = "")
    split_otu4$labelxs[which(split_otu4$value<0.01)]<-""

    p1<-ggplot(split_otu4, aes(x = '', y = value,fill= tax)) +
      geom_bar(stat = 'identity', show.legend = TRUE,width = 1) +
      #facet_wrap(.~variable,ncol = ncols) +
      geom_text(aes(y=rev(rev(value/2)+c(0,cumsum(rev(value))[-length(value)])),
                    x= 1.05,label=labelxs),
                size=4,colour=color_percent,family = "serif")+
      coord_polar(theta = 'y')+
      labs(title  = gname)+
      theme_void()+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            text = element_text(family = "serif"),
            axis.text.x = element_blank(), legend.position = "right",
            legend.text = element_text(family = "serif"),
            plot.background = element_blank())


    if (!is.null(color_taxa)) {
      color_value1<-color_taxa
      color_taxa<-color_value1[1:(taxa_num+1)]
      names(color_taxa)<-unique(split_otu3$tax)

      message("use self-selected color")

      p1<-p1+
        scale_fill_manual(name=taxa,values=color_taxa)

    } else {
      color_value1<-c(randomcolor((length(unique(split_otu3$tax))-1)),"grey90")
      names(color_value1)<-unique(split_otu3$tax)

      p1<-p1+scale_fill_manual(name=taxa,values = color_value1)
    }

    if (!is.null(xlabname)) p1<-p1+labs(title  = xlabname[t])

    pp[[t]]<-p1
  }

  pm<-wrap_plots(pp,ncol = ncols,guides="collect")

  ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_pieplot",".pdf",sep = ""),pm)
  cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_pieplot"," has been exported to ",export_path,sep = "","\n")

  return(pm)
}

"plotMicrocomHeatmap" <- function(microchatcomobj,
                                  rescale=FALSE,
                                  standardlization=c("0-1","scale","center"),
                                  ###0-1 standardlization
                                  data_type=c("absolute","relative"),
                                  sample_type=c("allsamples","groups"),
                                  color_taxa=NULL,
                                  color_group=NULL,
                                  color.heatmap=NULL,  ###two colors at least
                                  heatmap.add.line=TRUE,
                                  heatmap.top.class=5,
                                  heatmap.sample.angle=90,
                                  heatmap.taxa.angle=0,
                                  export_path="microbial composition") {
  data_type=match.arg(data_type)
  sample_type=match.arg(sample_type)
  standardlization=match.arg(standardlization)
  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (!rescale) {
  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
  }
  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
  }
  } else {
    if (data_type=="absolute") {
      if (sample_type=="allsamples") {
        if (standardlization=="0-1") genus_top<- scale1dt(microchatcomobj$wide_allsample_abs)
        if (standardlization=="scale") genus_top<-scale2dt(microchatcomobj$wide_allsample_abs)
        if (standardlization=="center") genus_top<-scale3dt(microchatcomobj$wide_allsample_abs)
        }
      if (sample_type=="groups") {

        if (standardlization=="0-1") genus_top<-scale1dt(microchatcomobj$wide_group_abs)
        if (standardlization=="scale") genus_top<-scale2dt(microchatcomobj$wide_group_abs)
        if (standardlization=="center") genus_top<-scale3dt(microchatcomobj$wide_group_abs)
      }
    }
    if (data_type=="relative") {
      if (sample_type=="allsamples") {
        if (standardlization=="0-1") genus_top<-scale1dt(microchatcomobj$wide_allsample_rel)
        if (standardlization=="scale") genus_top<-scale2dt(microchatcomobj$wide_allsample_rel)
        if (standardlization=="center") genus_top<-scale3dt(microchatcomobj$wide_allsample_rel)
      }
      if (sample_type=="groups") {
        if (standardlization=="0-1") genus_top<-scale1dt(microchatcomobj$wide_group_rel)
        if (standardlization=="scale") genus_top<-scale2dt(microchatcomobj$wide_group_rel)
        if (standardlization=="center") genus_top<-scale3dt(microchatcomobj$wide_group_rel)
      }
    }
  }


  net<-genus_top%>%as.matrix()

  samplename<-colnames(net)
  groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

  ###color groups
  if (is.null(color_group)) {
    color_group<-randomcolor(length(groupname))
  } else {
    if (length(color_group)<length(groupname)) {
      stop("There are no enough colors to fill all groups !!!")
    } else {
      color_group<-color_group[1:length(groupname)]
    }
  }

  names(color_group)<-groupname

  color_group1<-color_group%>%data.frame()
  color_sample<-randomcolor(length(samplename))

  ###color samples according to groups
  for (ttks in 1:nrow(color_group1)) {
    color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
  }
  names(color_sample)<-samplename

  ###color taxa
  taxa_name<-rownames(net)
  if (is.null(color_taxa)) {
    color_taxa<-randomcolor(length(taxa_name))
  } else {
    if (length(color_taxa)<length(taxa_name)) {
      stop("There are no enough colors to fill all taxa !!!")
    } else {
      color_taxa<-color_taxa[1:length(taxa_name)]
    }
  }
  names(color_taxa)<-taxa_name

  ###color heatmap blocks
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

  df1 <- data.frame(group = colnames(net))
  rownames(df1) <- colnames(net)
  color.use1 = color_sample

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
  color.use2 = color_taxa

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
      foo = ComplexHeatmap::anno_lines(colSums(abs(net[1:heatmap.top.class,])),
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


  p1<-ComplexHeatmap::Heatmap(net,
                          col = color.heatmap.use,
                          na_col = "white",
                          name = paste(stringr::str_to_title(data_type)," abundance",sep = ""),
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


  message("\n","The heatmap need to be saved manually.")
  return(p1)
}


"plotMicrocomChord" <- function(microchatcomobj,
                                plot.reverse=FALSE,
                                text.size=1,
                                data_type=c("absolute","relative"),
                                plot_type=c("bubble"),
                                sample_type=c("allsamples","groups"),
                                color_taxa=NULL,
                                color_group=NULL,
                                export_path="microbial composition") {
  data_type=match.arg(data_type)
  sample_type=match.arg(sample_type)

  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object !!!")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
    split_otu3<-microchatcomobj$long_group_abs
  }

  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
    split_otu3<-microchatcomobj$long_group_rel
  }



  genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
  genus_topxx<-reshape2::melt(genus_topxx)

  samplename<-unique(genus_topxx$variable)
  groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

  genus_topxx$group<-substr(genus_topxx$variable,start = 1,stop = 2)
  if (plot.reverse) df.plot<-subset(genus_topxx,select=c(group,tax,value))
  if (!plot.reverse) df.plot<-subset(genus_topxx,select=c(tax,group,value))
  order.sector<- c(unique(df.plot$group)[length(unique(df.plot$group)):1],
                   unique(df.plot$tax))

  grid.color<-c(color_group[1:length(groupname)],
                c(color_taxa[1:taxa_num],"grey50"))

  names(grid.color)<-c(unique(df.plot$group),unique(df.plot$tax))


  circlize::circos.clear()
  par(family="serif",mar=c(0,0,1,1))
  circlize::chordDiagram(df.plot,
                         order = order.sector,
                         diffHeight = 0.06,
                         grid.col = grid.color,
                         transparency = 0.4,
                         directional = 1,
                         direction.type = c("diffHeight","arrows"),
                         link.arr.type = "big.arrow",
                         annotationTrack = "grid",
                         annotationTrackHeight = 0.03,
                         preAllocateTracks = list(track.height = 0.1),
                         small.gap = 1,
                         big.gap = 10,
                         link.visible = TRUE,
                         scale = FALSE,
                         reduce = -1)

  #for(i in 1:nrow(df.plot)){highlight.sector(sector.index = paste0("chr",i),col=color_taxa[i])}

  circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
    xlim = circlize::get.cell.meta.data("xlim")
    xplot = circlize::get.cell.meta.data("xplot")
    ylim = circlize::get.cell.meta.data("ylim")
    sector.name = circlize::get.cell.meta.data("sector.index")
    circlize::circos.text(mean(xlim), ylim[2], sector.name,
                          facing = "clockwise",family = par("serif"),
                          niceFacing = TRUE, adj = c(0, 0), cex = text.size)
  }, bg.border = NA)

  lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                type = "grid",ncol=2,
                                legend_gp = grid::gpar(fill = grid.color,
                                                       fontfamily="serif",
                                                       fontface = "plain"))

  ComplexHeatmap::draw(lgd, x = unit(0.4, "npc") ,
                       y = unit(0.8, "npc"), just = c("right","bottom"))

  p1<-recordPlot()
  circlize::circos.clear()
  message("\n","The chord plot need to be saved manually.")
  return(p1)
}



"plotMicrocomBubble" <- function(microchatcomobj,
                                 data_type=c("absolute","relative"),
                                 sample_type=c("allsamples","groups"),
                                 color_taxa=NULL,
                                 color_group=NULL,
                                 bubble.color.taxa=TRUE,
                                 bubble.max.size=9,
                                 bubble.shape=19,
                                 bubble.alpha=1,
                                 aspect.ratio=0.8,
                                 color_background=NULL,
                                 color_border=NULL,
                                 export_path="microbial composition") {
  data_type=match.arg(data_type)
  sample_type=match.arg(sample_type)

  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
    split_otu3<-microchatcomobj$long_group_abs
  }

  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
    split_otu3<-microchatcomobj$long_group_rel
  }

  if (bubble.color.taxa) {
    cat("\n","color taxa according to the params 'color_taxa'")
    genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
    genus_topxx<-reshape2::melt(genus_topxx)

    color_taxa<-c(color_taxa[1:taxa_num],"grey50")
    names(color_taxa)<-unique(genus_topxx$tax)

    p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=tax))+
      geom_point(stat= "identity",aes(size=value),
                 shape=bubble.shape,
                 alpha=bubble.alpha,show.legend = TRUE)+
      scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
      scale_color_manual(name=taxa,values=color_taxa)+
      theme(aspect.ratio = aspect.ratio,
            axis.ticks.length = unit(0.2,"lines"),
            axis.ticks = element_line(color='black'),
            panel.background = element_rect(fill = "grey90"),
            #axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_text(colour='black', size=18,face = "bold",family = "serif",vjust = 1.5),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))


  } else {
    cat("\n","color taxa according to the params 'color_group'")

    genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
    genus_topxx<-reshape2::melt(genus_topxx)

    samplename<-unique(genus_topxx$variable)
    groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

    color_group<-color_group[1:length(groupname)]
    names(color_group)<-groupname

    color_group1<-color_group%>%data.frame()
    color_sample<-randomcolor(length(samplename))

    ###color samples according to groups
    for (ttks in 1:nrow(color_group1)) {
      color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
    }
    names(color_sample)<-samplename

    p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=factor(variable)))+
      geom_point(stat= "identity",aes(size=value),
                 shape=bubble.shape,
                 alpha=bubble.alpha,show.legend = TRUE)+
      scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
      scale_color_manual(name=taxa,values=color_sample)+
      theme(aspect.ratio = aspect.ratio,
            axis.ticks.length = unit(0.2,"lines"),
            axis.ticks = element_line(color='black'),
            panel.background = element_rect(fill = "grey90"),
            #axis.line = element_line(colour = "black"),
            legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_text(colour='black',size=10,family = "serif"),
            axis.text.x=element_text(colour = "black",size = 10,
                                     angle = 0,hjust = 0.5,vjust =0.5,family = "serif"))


  }

  if (!is.null(color_background)) {
    p1<-p1+
      theme(panel.background=element_rect(fill = color_background))
  }
  if (is.null(color_border)) {
    p1<-p1
  } else {
    p1<-p1+theme(panel.border = element_rect(fill=NA,colour = color_border))
  }
  ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_bubble",".pdf",sep = ""),p1)
  cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_bubble"," has been exported to ",export_path,sep = "","\n")

  return(p1)
}


"plotMicrocom" <- function(microchatcomobj,
                               data_type=c("absolute","relative"),
                               plot_type=c("barplot","pie","heatmap","chord","bubble"),
                               sample_type=c("allsamples","groups"),
                               color_taxa=NULL,
                               color_group=NULL,

                               barplot.barwidth=0.6,

                               bubble.color.taxa=TRUE,
                               bubble.max.size=9,
                               bubble.shape=19,

                               color.heatmap=NULL,  ###至少两种颜色
                               heatmap.add.line=TRUE,
                               heatmap.top.class=5,
                               heatmap.sample.angle=90,
                               heatmap.taxa.angle=0,

                               aspect.ratio=0.8,
                               export_path="microbial composition") {
  data_type=match.arg(data_type)
  plot_type=match.arg(plot_type)
  sample_type=match.arg(sample_type)

  if (class(microchatcomobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  plot_data<-microchatcomobj$plot_data

  export_path<-microchatcomobj$export_path
  taxa_num<-microchatcomobj$selected_taxanum
  taxa<-microchatcomobj$selected_taxa

  dir.create(export_path, recursive = TRUE)

  if (data_type=="absolute") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_abs
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_abs
    split_otu3<-microchatcomobj$long_group_abs
  }

  if (data_type=="relative") {
    if (sample_type=="allsamples") genus_top<-microchatcomobj$wide_allsample_rel
    if (sample_type=="groups") genus_top<-microchatcomobj$wide_group_rel
    split_otu3<-microchatcomobj$long_group_rel
  }

  if (plot_type=="barplot") {

    if (sample_type=="groups") {
  color_value<-c(randomcolor((length(unique(split_otu3$tax))-1)),"grey90")
  names(color_value)<-unique(split_otu3$tax)[order(unique(split_otu3$tax))]

  p1<-ggplot(data = split_otu3 ,
             aes(x = variable, y = value,
                                    fill = tax)) +
    geom_col(position = 'fill', width = barplot.barwidth) +
    scale_fill_manual(values = color_value) +
    labs(x = '', y = 'Relative Abundance(%)') +
    theme_bw()+
    theme(panel.grid.minor.y = element_line(colour = "black"),
          text = element_text(family = "serif"),
          panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),aspect.ratio = aspect.ratio)

  if (!is.null(color_taxa)) {
    color_value<-color_taxa
    color_taxa<-color_value[1:(taxa_num+1)]
    names(color_taxa)<-unique(split_otu3$tax)[order(unique(split_otu3$tax))]

    message("use self-selected color")
    p1<-p1+
      ggnewscale::new_scale_color()+
      ggnewscale::new_scale_fill()+
      scale_discrete_manual(values=color_taxa,
                            aesthetics = "colour")+
      scale_discrete_manual(values=color_taxa,
                            aesthetics = "fill")

  }}
    if (sample_type=="allsamples") {
      genus_topt<-tibble::rownames_to_column(genus_top,var = "tax")
      genus_topt<-reshape2::melt(genus_topt)

      color_value<-c(randomcolor((length(unique(genus_topt$tax))-1)),"grey90")
      names(color_value)<-unique(genus_topt$tax)[order(unique(genus_topt$tax))]

      p1<-ggplot(data = genus_topt ,
                 aes(x = variable, y = value,
                     fill = tax)) +
        geom_col(position = 'fill', width = barplot.barwidth) +
        scale_fill_manual(values = color_value) +
        labs(x = '', y = 'Relative Abundance(%)') +
        theme_bw()+
        theme(panel.grid.minor.y = element_line(colour = "black"),
              text = element_text(family = "serif"),
              panel.background = element_rect(color = 'black',
                                              fill = 'transparent'),aspect.ratio = aspect.ratio)

      if (!is.null(color_taxa)) {
        color_value<-color_taxa
        color_taxa<-color_value[1:(taxa_num+1)]
        names(color_taxa)<-unique(genus_topt$tax)[order(unique(genus_topt$tax))]

        message("use self-selected color")
        p1<-p1+
          ggnewscale::new_scale_color()+
          ggnewscale::new_scale_fill()+
          scale_discrete_manual(values=color_taxa,
                                aesthetics = "colour")+
          scale_discrete_manual(values=color_taxa,
                                aesthetics = "fill")

      }}

    ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type,".pdf",sep = ""),p1)
    cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type," has been exported to ",export_path,sep = "","\n")

}

  if (plot_type=="pie") {

    color_value1<-c(randomcolor((length(unique(split_otu3$tax))-1)),"grey90")
    names(color_value1)<-unique(split_otu3$tax)[order(unique(split_otu3$tax))]

    p1<-ggplot(split_otu3, aes(x = '', y = value,fill= tax)) +
      geom_bar(stat = 'identity', show.legend = TRUE,width = 1) +
      facet_wrap(.~variable,ncol = 3) +
      coord_polar(theta = 'y') +
      scale_fill_manual(values = color_value1) +
      theme_void()+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            text = element_text(family = "serif"),
            axis.text.x = element_blank(), legend.position = "right",
            legend.text = element_text(family = "serif"),
            plot.background = element_blank())

    if (!is.null(color_taxa)) {
      color_value1<-color_taxa
      color_taxa<-color_value1[1:(taxa_num+1)]
      names(color_taxa)<-unique(split_otu3$tax)[order(unique(split_otu3$tax))]

      message("use self-selected color")

      p1<-p1+
        ggnewscale::new_scale_fill()+
        scale_discrete_manual(values=color_taxa,
                                   aesthetics = "fill")

    }
    ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type,".pdf",sep = ""),p1)
    cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type," has been exported to ",export_path,sep = "","\n")

  }

  if (plot_type=="heatmap") {

    net<-genus_top%>%as.matrix()

    samplename<-colnames(net)
    groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

    ###color groups
    if (is.null(color_group)) {
      color_group<-randomcolor(length(groupname))
    } else {
      if (length(color_group)<length(groupname)) {
        stop("There are no enough colors to fill all groups !!!")
      } else {
        color_group<-color_group[1:length(color_group)]
      }
    }

    names(color_group)<-groupname

    color_group1<-color_group%>%data.frame()
    color_sample<-randomcolor(length(samplename))

    ###color samples according to groups
    for (ttks in 1:nrow(color_group1)) {
      color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
    }
    names(color_sample)<-samplename

    ###color taxa
    taxa_name<-rownames(net)
    if (is.null(color_taxa)) {
      color_taxa<-randomcolor(length(taxa_name))
    } else {
      if (length(color_taxa)<length(taxa_name)) {
        stop("There are no enough colors to fill all taxa !!!")
      } else {
        color_taxa<-color_taxa[1:length(taxa_name)]
        }
      }
    names(color_taxa)<-taxa_name

    ###color heatmap blocks
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

    df1 <- data.frame(group = colnames(net))
    rownames(df1) <- colnames(net)
    color.use1 = color_sample

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

    color.use2 = color_taxa

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
                                        foo = anno_lines(rowSums(abs(net)),
                                                         gp = grid::gpar(col = "red"),
                                                         add_points = TRUE,smooth = TRUE, axis = FALSE,
                                                         extend = 0.2,border = FALSE,
                                                         pt_gp =grid::gpar(fill = color.use2,col = color.use2),
                                                         height = unit(1, "cm")),
                                        show_annotation_name = FALSE)
    ha2 = ComplexHeatmap::HeatmapAnnotation(
                                            #pt = anno_points(colSums(abs(net[1:heatmap.top.class,])), gp = gpar(fill = color.use1,col = color.use1),height = unit(1, "cm")),
                                            foo = anno_lines(colSums(abs(net[1:heatmap.top.class,])),
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

    pdf(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type,".pdf",sep = ""),
        width=10, height=10)
    ComplexHeatmap::Heatmap(net,
                            col = color.heatmap.use,
                            na_col = "white",
                            name = paste(stringr::str_to_title(data_type)," abundance",sep = ""),
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

    dev.off()

    p1<-ComplexHeatmap::Heatmap(net,
                                col = color.heatmap.use,
                                na_col = "white",
                            name = paste(stringr::str_to_title(data_type)," abundance",sep = ""),
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

  if (plot_type=="bubble") {
    if (bubble.color.taxa) {
      cat("\n","color taxa according to the params 'color_taxa'")
      genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
      genus_topxx<-reshape2::melt(genus_topxx)

      color_taxa<-c(color_taxa[1:taxa_num],"grey50")
      names(color_taxa)<-unique(genus_topxx$tax)

      p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=tax))+
        geom_point(stat= "identity",aes(size=value),
                   shape=bubble.shape,
                   alpha=bubble.alpha,show.legend = TRUE)+
        scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
        scale_color_manual(values=color_taxa)+
        theme(aspect.ratio = aspect.ratio,
              axis.title = element_blank(),
              text = element_text(family = "serif"),
              legend.position = "none")

      ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type,".pdf",sep = ""),p1)
      cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type," has been exported to ",export_path,sep = "","\n")

    } else {
      cat("\n","color taxa according to the params 'color_group'")

      genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
      genus_topxx<-reshape2::melt(genus_topxx)

      samplename<-unique(genus_topxx$variable)
      groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

      color_group<-color_group[1:length(groupname)]
      names(color_group)<-groupname

      color_group1<-color_group%>%data.frame()
      color_sample<-randomcolor(length(samplename))

      ###color samples according to groups
      for (ttks in 1:nrow(color_group1)) {
        color_sample[which(substr(samplename,start = 1, stop = 2)==rownames(color_group1)[ttks])]<-color_group1[ttks,]
      }
      names(color_sample)<-samplename

      p1<-ggplot(data=genus_topxx, mapping=aes(x=variable,y=tax,color=factor(variable)))+
        geom_point(stat= "identity",aes(size=value),
                   shape=bubble.shape,
                   alpha=bubble.alpha,show.legend = TRUE)+
        scale_size(range = c(0.1, bubble.max.size),guide=FALSE)+
        scale_color_manual(values=color_sample)+
        theme(aspect.ratio = aspect.ratio,
              axis.title = element_blank(),
              text = element_text(family = "serif"),
              legend.position = "none")

      ggsave(paste(export_path,"/top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type,".pdf",sep = ""),p1)
      cat("\n","top_",taxa_num," (",taxa,") ",data_type," abundance of ",sample_type,"_",plot_type," has been exported to ",export_path,sep = "","\n")

    }
  }

  if (plot_type=="chord") {
    genus_topxx<-tibble::rownames_to_column(genus_top,var = "tax")
    genus_topxx<-reshape2::melt(genus_topxx)

    samplename<-unique(genus_topxx$variable)
    groupname<-substr(samplename,start = 1, stop = 2)%>%unique()

    genus_topxx$group<-substr(genus_topxx$variable,start = 1,stop = 2)
    df.plot<-subset(genus_topxx,select=c(group,tax,value))
    order.sector<- c(unique(df.plot$group)[length(unique(df.plot$group)):1],
                     unique(df.plot$tax))

    grid.color<-c(color_group[1:length(groupname)],
                  c(color_taxa[1:taxa_num],"grey50"))

    names(grid.color)<-c(unique(df.plot$group),unique(df.plot$tax))


    circlize::circos.clear()
    par(family="serif",mar=c(0,0,1,1))
    circlize::chordDiagram(df.plot,
                               order = order.sector,
                               diffHeight = 0.06,
                               grid.col = grid.color,
                               transparency = 0.4,
                               directional = 1,
                               direction.type = c("diffHeight","arrows"),
                               link.arr.type = "big.arrow",
                               annotationTrack = "grid",
                               annotationTrackHeight = 0.03,
                               preAllocateTracks = list(track.height = 0.1),
                               small.gap = 1,
                               big.gap = 10,
                               link.visible = TRUE,
                               scale = FALSE,
                               reduce = -1)

    circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
      xlim = circlize::get.cell.meta.data("xlim")
      xplot = circlize::get.cell.meta.data("xplot")
      ylim = circlize::get.cell.meta.data("ylim")
      sector.name = circlize::get.cell.meta.data("sector.index")
      circlize::circos.text(mean(xlim), ylim[2], sector.name,
                  facing = "clockwise",family = par("serif"),
                  niceFacing = TRUE, adj = c(-0.25, 0), cex = 1)
    }, bg.border = NA)

    lgd <- ComplexHeatmap::Legend(at = names(grid.color),
                                  type = "grid",ncol=2,
                                  legend_gp = grid::gpar(fill = grid.color,
                                                         fontfamily="serif",
                                                         fontface = "plain"))

    ComplexHeatmap::draw(lgd, x = unit(0.4, "npc") ,
                         y = unit(0.8, "npc"), just = c("right","bottom"))

    p1<-recordPlot()

    message("\n","The chord plot need to be saved manually.")
    }
  return(p1)
}

"calcMicrochatcom" <- function(mchat,
                           taxa="Class",
                           taxa_num=10,
                           export_path="microbial composition") {

  if (class(mchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  dir.create(export_path, recursive = TRUE)

  if (length(which(colSums(mchat$otu_table)==0)) !=0) mchat$otu_table<-mchat$otu_table[,-which(colSums(mchat$otu_table)==0)]
  if (length(which(colSums(mchat$otu_table)==0)) ==0) mchat$otu_table<-mchat$otu_table
  data_prep_new<-data_prep(mchat)
  data<-data_prep_new$data
  group_order<-data_prep_new$group

  otu_table  <- data$otu_table
  taxonomy <- data$taxon_table

  taxonomy$OTU_ID <-rownames(taxonomy)
  metadata <-data$sampledata
  metadata$sample <-rownames(metadata)

  otu_table$OTU_ID <- rownames(otu_table)

  taxa_pos<-which(taxa == colnames(taxonomy))

  genus <- subset(taxonomy, select = c(OTU_ID, taxa_pos))

  temp <- merge(otu_table, genus, by = "OTU_ID", )
  temp <- subset(temp, select = -OTU_ID)

  colnames(temp)[length(colnames(temp))]<-"taxono"
  temp %>%
    group_by(taxono) %>%
    summarise_all(sum) -> genus_table

  genus_table[which(genus_table$taxono==""),]$taxono<-paste(tolower(substr(taxa,start = 1,stop=1)),"__","unclassfied",sep = "")

  genus_table <- as.data.frame(genus_table)
  genus_table <- na.omit(genus_table)
  rownames(genus_table) <- genus_table$taxono
  genus_table <- subset(genus_table,
                        select = -c(taxono))
  # 排序
  genus_table$rowsum <- rowSums(genus_table)
  genus_table <- genus_table[order(genus_table$rowsum,
                                   decreasing = TRUE), ]
  if (length(rownames(genus_table)) >= taxa_num+1) taxa_num<-taxa_num
  if (length(rownames(genus_table)) < taxa_num+1) {taxa_num<-length(rownames(genus_table))
  message("All taxa should be visualized due to the larger selected taxa size.")
  }

  genus_top10 <- genus_table[1:taxa_num, ]
  genus_top10['Others', ] <- colSums(genus_table) - colSums(genus_top10)

  # check
  "check_top" <- function(top_table, otu_table){
    judge <- colSums(genus_top10) == colSums(otu_table)
    if (sum(judge) == length(otu_table)){
      print("Everything looks good")}
    else{print("something wrong.")}
  }
  check_top(genus_top10, genus_table)

  # make factor levels
  genus_top10$Taxo <- forcats::fct_inorder(rownames(genus_top10))
  genus_topxx<-subset(genus_top10,select = -c(rowsum,Taxo))
  genus_topxx<-genus_topxx/colSums(genus_topxx)


  temp <- pivot_longer(data = genus_top10,
                       cols = -c(Taxo, rowsum),
                       names_to = "sample",
                       values_to = "value")

  plot_data <- merge(temp, metadata, by = "sample")

  genus_top<-subset(genus_top10,select = -c(rowsum,Taxo))
  split_otu1<-split_merge(genus_top)$split_otu
  split_otu2<-split_merge(genus_topxx)$split_otu

  split_otu3<-split_merge(genus_top)$split_otu1
  split_otu4<-split_merge(genus_topxx)$split_otu1
  cat("-----------------------------------------------------------------------")


  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_absabun.txt",sep = "" )
  write.table(genus_top,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") absolute abundance of all samples has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_pieplot.txt",sep = "" )
  write.table(genus_topxx,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") relative abundance of all samples has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_barplot.txt",sep = "" )
  write.table(split_otu1,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") absolute abundance of all groups has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_pieplot.txt",sep = "" )
  write.table(split_otu2,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") relative abundance of all groups has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_barplot.txt",sep = "" )
  write.table(split_otu3,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") absolute abundance of all groups (long data) has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  file2=paste(export_path, "/top_",taxa_num,"_abun_",taxa,"_pieplot.txt",sep = "" )
  write.table(split_otu4,file = file2, row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","top ",taxa_num," (",taxa,") relative abundance of all groups (long data) has been saved in the relative path"," '",export_path,"'",sep = "")
  cat("\n","-----------------------------------------------------------------------")

  ###create S3 object
  microchatcomobj <- list(plot_data,
                          genus_top,
                          genus_topxx,
                          split_otu1,
                          split_otu2,
                          split_otu3,
                          split_otu4,
                          taxa,
                          taxa_num,
                          export_path)

  class(microchatcomobj) <- "microchat"

  names(microchatcomobj)[[1]]<-"plot_data"
  names(microchatcomobj)[[2]]<-"wide_allsample_abs"
  names(microchatcomobj)[[3]]<-"wide_allsample_rel"
  names(microchatcomobj)[[4]]<-"wide_group_abs"
  names(microchatcomobj)[[5]]<-"wide_group_rel"
  names(microchatcomobj)[[6]]<-"long_group_abs"
  names(microchatcomobj)[[7]]<-"long_group_rel"
  names(microchatcomobj)[[8]]<-"selected_taxa"
  names(microchatcomobj)[[9]]<-"selected_taxanum"
  names(microchatcomobj)[[10]]<-"export_path"

  return(microchatcomobj)
}


"calcMicrochatVenn" <- function(submchat) {

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object !!!")
  }

  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table

  otudata <- submchat$otu_table
  taxon_table <- submchat$taxon_table

  data_ve<-filter_merge_sum (otudata,taxon_table)
  data_venn<-data_ve$data_venn
  data_venn1<-data_ve$data_with_taxa
  data_venn2<-data_ve$data_with_2taxa

  ###create S3 object
  microchatVennobj <- list(data_venn,
                          data_venn1,
                          data_venn2,
                          otudata,
                          taxon_table)

  class(microchatVennobj) <- "microchat"

  names(microchatVennobj)[[1]]<-"data_venn"
  names(microchatVennobj)[[2]]<-"data_venn1"
  names(microchatVennobj)[[3]]<-"data_venn2"
  names(microchatVennobj)[[4]]<-"otu_table"
  names(microchatVennobj)[[5]]<-"taxon_table"

  return(microchatVennobj)
}


"plotMicrochatVenn" <- function(microchatVennobj,
                                inter.point.color="purple",
                                shade_color="red",
                                bar_color="black",
                                group_bar_color=colorCustom(40,pal = "gygn"),
                                query=list(
                                  list(query=UpSetR::elements,
                                       params = c('Phylum', 'p__Firmicutes'),
                                       color="blue",active=T))
                                ) {


  if (class(microchatVennobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_venn2<-microchatVennobj$data_venn2
  data_venn<-microchatVennobj$data_venn
  data_venn1<-microchatVennobj$data_venn1

  groupnum<-length(rownames(data_venn))

  group_bar_color<-group_bar_color[1:groupnum]

  message("\n","The upset plot need to be saved manually.")
  query<-query
  plot_single_upset(data_venn2,
                    group_num=groupnum,
                    bar_num=NA,
                    ylabel="OTU number",
                    xlabel="OTU number",
                    exist_color=inter.point.color,
                    shade_color = shade_color,
                    group_bar_color=group_bar_color,
                    query=query,
                    bar_color=bar_color)
}

"plotMicrochatVennDia" <- function(microchatVennobj) {

  if (class(microchatVennobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  otudata<-microchatVennobj$otu_table

  groupname<-substr(colnames(otudata),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(otudata))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()

  if (as.numeric(matchnum) == as.numeric(allnum)) {
    split_otu <- otu_split1(otudata,num=2)
  } else {
    split_otu <- otu_split2(otudata,num=2)
  }

  ttt<-lapply(split_otu, function(tt){
    yy<-rownames(tt)
    return(yy)
  })
par(family="serif")
  pp<-ggVennDiagram::ggVennDiagram(ttt,
                                   label_color = "blue", label_size = 2,
                                   set_size = 6, label = "both", label_percent_digit = 1) +
    scale_x_continuous(expand = expansion(mult = .2))+
    scale_colour_hue()+
    theme(text = element_text(family = "serif"))
  message("\n","The upset plot need to be saved manually.")
  return(pp)
}

"plotMicrochatVennFlower" <- function(microchatVennobj,alpha_thres=100) {
  if (class(microchatVennobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  otudata<-microchatVennobj$otu_table
  taxon_table<-microchatVennobj$taxon_table

  ttabun<-t(otudata)%>%data.frame()
  data_flower<-venn_2id(data=ttabun,taxon=taxon_table)
  flower_table<-data_flower$data_flower
  par(family="serif")
  plot_flower(flower_table,alpha_thres=alpha_thres)
  message("\n","The upset plot need to be saved manually.")
}

"plot_single_upset" <- function(data,
                                group_num=6,
                                bar_num=NA,ylabel="OTU number",xlabel="OTU number",
                                exist_color="purple",shade_color = "red",
                                group_bar_color=random2color(5),
                                query=query,bar_color="green") {
  par(family="serif")
  UpSetR::upset(data,
        nsets = group_num,
        nintersects=bar_num,
        point.size = 4,
        line.size = 1,
        mainbar.y.label = ylabel,
        sets.x.label= xlabel,
        matrix.color= exist_color,
        text.scale = c(1.5, 1.5, 1.5, 1.5,1.5, 1.5),
        shade.color = shade_color,
        order.by =c( "freq" ),
        sets.bar.color=group_bar_color,
        queries=query,
        main.bar.color = bar_color)

}




"plotMicrochatComplexVenn" <- function(microchatVennobj,
                                       bar.taxa= "Class",
                                       pie.taxa= "Genus",
                                       group_bar_color=colorCustom(40,pal = "gygn"),
                                       bar_color="black",
                                       color_taxa=colorCustom(40,pal = "gygn"),
                                       select_group="all",
                                       show.legend=TRUE,
                                       query=list(list(query = UpSetR::intersects,
                                                       params = rownames(data_venn),
                                                       color = 'blue', active = TRUE))) {

  if (class(microchatVennobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }

  data_venn2<-microchatVennobj$data_venn2
  data_venn<-microchatVennobj$data_venn
  data_venn1<-microchatVennobj$data_venn1

  groupnum<-length(rownames(data_venn))

  group_bar_color<-group_bar_color[1:groupnum]

  taxon<-subset(data_venn2,select = c(1,(groupnum+2):length(colnames(data_venn2))))
  taxon<-tibble::column_to_rownames(taxon,var = "name")

  plot_tt<-attr_plot_newdata(data_venn=data_venn,taxon=taxon,
                          taxa= bar.taxa,
                          select_group=select_group)

  phylum_top10<-plot_tt$phylum_top10
  phylum_melt<-plot_tt$phylum_melt
  otu_num<-length(phylum_melt$tax)

  plot_tt1<-attr_plot_newdata(data_venn=data_venn,taxon=taxon,
                           taxa= pie.taxa,
                           select_group=select_group)

  phylum_top101<-plot_tt1$phylum_top10
  phylum_melt1<-plot_tt1$phylum_melt
  otu_num<-length(phylum_melt$tax)

  color.use<-color_taxa[1:length(unique(phylum_melt$tax))]
  names(color.use)<-unique(phylum_melt$tax)

  color.use1<-color_taxa[1:length(unique(phylum_melt1$tax))]
  names(color.use1)<-unique(phylum_melt1$tax)

  message("\n","The upset plot need to be saved manually.")
  "plot_bar" <- function(mydata, x, y) {
    p1<-ggplot(phylum_melt, aes(x = variable, y = value, fill = tax)) +
      geom_col(position = 'stack', width = 0.6) +
      scale_fill_manual(values = color.use,
                        guide=guide_legend(keywidth = 0.3,
                                           keyheight = 0.3,
                                           order=2,ncol=1,
                                           override.aes=list(size=2))) +
      #theme_void()+
      theme(panel.grid = element_blank(), legend.position = "right",
            legend.text = element_text(family = "serif"),
            text = element_text(family = "serif"),
            panel.background = element_rect(color = 'black', fill = 'transparent')) +
      labs(x = '', y = 'Relative abundance (%)', fill = NULL)
    if (show.legend) {
      p1<-p1
    } else {
      p1<-p1 +theme(legend.position = "none")
    }

  }

  "plot_pie1" <- function(mydata, x, y) {
    p1<-ggplot(phylum_melt1, aes(x = '', y = value,fill= tax)) +
      geom_bar(stat = 'identity', show.legend = TRUE) +
      facet_wrap(.~variable,ncol = 3) +
      coord_polar(theta = 'y') +
      scale_fill_manual(values = color.use1,
                        guide=guide_legend(keywidth = 0.3,
                                           keyheight = 0.3,
                                           order=1,ncol=1,
                                           override.aes=list(size=2))) +
      theme_void()+
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            text = element_text(family = "serif"),
            axis.text.x = element_blank(), legend.position = "right",
            legend.text = element_text(family = "serif"),
            plot.background = element_blank()) +
      labs(x = NULL, family="serif",
           y = paste('Number of all shared OTUs: ',
                     length(otu_num),
                     "\nAverage abundance of main ",sep = ""))

    if (show.legend) {
      p1<-p1
    } else {
      p1<-p1 +theme(legend.position = "none")
    }
  }

  UpSetR::upset(data_venn2,
        nsets = groupnum,
        sets.bar.color=group_bar_color,
        main.bar.color = bar_color,
        order.by = c('freq', 'degree'),
        decreasing = c(TRUE, TRUE),
        queries = query,
        attribute.plots = list(gridrows = 100, ncols = 2,
                               plots = list(
                                 list(plot = plot_bar, mydata = NA, x = NA, y = NA, queries = FALSE),
                                 list(plot = plot_pie1, mydata = NA, x = NA, y = NA, queries = FALSE))))

}

"plot_flower" <- function(flower_dat,alpha_thres=alpha_thres) {
  group<-substr(colnames(flower_dat),start = 1,stop = 2)%>%unique()
  sample<-substr(colnames(flower_dat),start = 1,stop = 2)%>%length()
  group_num<-length(group)
  if (is.integer(sample/group_num)) {
    sample_num<-sample/group_num
  } else {
    sample_num<-((sample/group_num)%>%as.integer())+1
  }


  sample_id <- colnames(flower_dat)
  sample_id
  otu_id <- unique(flower_dat[,1])
  otu_id
  head(otu_id)
  otu_id <- otu_id[otu_id != '']
  otu_id
  head(otu_id)
  core_otu_id <- otu_id
  otu_num <- length(otu_id)
  otu_num

  for (i in 2:ncol(flower_dat)) {
    otu_id <- unique(flower_dat[,i])
    otu_id <- otu_id[otu_id != '']
    core_otu_id <- intersect(core_otu_id, otu_id)
    otu_num <- c(otu_num, length(otu_id))
  }
  core_num <- length(core_otu_id)
  core_num
  library(plotrix)


  color.use<-random2color(group_num)
  names(color.use)<-group

  color.use.data<-color.use%>%data.frame()

  samplegroup<-substr(colnames(flower_dat),start = 1,stop = 2)

  color.use1<-random2color(sample)

   for (j in 1:length(samplegroup)) {
     color.use1[j]<-color.use.data$.[match(samplegroup[j], rownames(color.use.data))]
   }


  ellipse_col <- color.use1


  "flower_plot" <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
    par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
    plot(c(0,10),c(0,10),type='n')
    n   <- length(sample)
    deg <- 360 / n
    res <- lapply(1:n, function(t){
      draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                   y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                   col = ellipse_col[t],
                   border = ellipse_col[t],
                   a = a, b = b, angle = deg * (t - 1))
      text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
           otu_num[t])

      if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = 1
        )
      } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = 1
        )
      }
    })
    draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
    text(x = 5, y = 5, paste(core_otu))
  }

  flower_plot(
    sample = sample_id, otu_num = otu_num,
    core_otu = core_num,
    start = 90, a = 0.5, b = 2, r = 0.5,
    ellipse_col = ellipse_col, circle_col = 'white')

}


"rel_2abun" <- function(data,num=group_num) {

  data<-tibble::column_to_rownames(data,var = "tax")

  sum<-colSums(subset(data,select=c(1:num)))

  for ( k in 1:num) {
    data[,k]<-data[,k]/sum[k]*100
  }
  data<-tibble::rownames_to_column(data,var = "tax")
  return(data)
}

"attr_plot_newdata" <- function(data_venn,taxon,taxa= "Class",
                             select_group="all") {
  group_num<-length(rownames(data_venn))
  data10<-t(data_venn)%>%data.frame()
  data11<-data10
  data11[data11!=0]<-1
  if(select_group=="all") {
    select_otu <- rownames(data11[rowSums(data11) == group_num, ])
    otu_select <- data10[select_otu, ]
    otu_select$name<-rownames(otu_select)
    taxon1<-tibble::rownames_to_column(taxon,var = "name")
    taxa_select<-merge(otu_select,taxon1,by="name")

    taxa_order<-which(colnames(taxa_select)==taxa)
    taxa_select<-subset(taxa_select,select=c(2:(group_num+1),taxa_order))
    colnames(taxa_select)[length(colnames(taxa_select))]<-"tax"

    if (length(which(taxa_select$tax=="")) !=0) taxa_select[which(taxa_select$tax==""),]$tax<-paste(tolower(substr(taxa,start = 1,stop = 1)),"__Unclassified",sep = "")
    taxa_select %>%
      group_by(tax) %>%
      summarise_all(sum) -> taxa1

    taxa1<-rel_2abun(taxa1,num=group_num)
  } else {
    data11 <- data11[rowSums(data11) == 1, ]
    group_select<-match(select_group,colnames(data11))
    data12<-subset(data11,select=group_select)%>%data.frame()
    data12$name<-rownames(data12)
    select_otu <- data12[which(data12[1] == 1), ]$name

    otu_select <- data10[select_otu, ]
    otu_select$name<-rownames(otu_select)
    taxon1<-tibble::rownames_to_column(taxon,var = "name")
    taxa_select<-merge(otu_select,taxon1,by="name")

    taxa_order<-which(colnames(taxa_select)==taxa)
    taxa_select<-subset(taxa_select,select=c(2:(group_num+1),taxa_order))
    colnames(taxa_select)[length(colnames(taxa_select))]<-"tax"

    if (length(which(taxa_select$tax=="")) !=0) taxa_select[which(taxa_select$tax==""),]$tax<-paste(tolower(substr(taxa,start = 1,stop = 1)),"__Unclassified",sep = "")

    taxa_select<-subset(taxa_select,select=c(group_select,tax))

    taxa_select %>%
      group_by(tax) %>%
      summarise_all(sum) -> taxa1

    taxa1<-rel_2abun(taxa1,num=1)
  }

  phylum_melt <-reshape2::melt(taxa1, id = 'tax')
  phylum_melt_stat <- doBy::summaryBy(value~tax, phylum_melt, FUN = mean)

  phylum_melt_stat <- phylum_melt_stat[order(-phylum_melt_stat$value.mean), ]

  if (length(phylum_melt_stat$tax)>=10) {
    phylum_top10 <- phylum_melt_stat[1:10, ]
    rownames(phylum_top10)<-phylum_top10$tax

    phylum_top10<-phylum_top10[,-1]%>%data.frame()
    phylum_top10['Others', ] <- sum(phylum_melt_stat$value.mean) - sum(phylum_top10)
    rownames(phylum_top10)[1:10]<-phylum_melt_stat$tax[1:10]
  } else {
    phylum_top10 <- phylum_melt_stat[1:length(phylum_melt_stat$tax), ]
    rownames(phylum_top10)<-phylum_top10$tax

    phylum_top10<-phylum_top10[,-1]%>%data.frame()
    phylum_top10['Others', ] <- sum(phylum_melt_stat$value.mean) - sum(phylum_top10)
    rownames(phylum_top10)[1:length(phylum_melt_stat$tax)]<-phylum_melt_stat$tax[1:length(phylum_melt_stat$tax)]
  }


  phylum_top10$name<-rownames(phylum_top10)
  colnames(phylum_top10)<-c("value","tax")
  phylum_top10$tax<-factor(phylum_top10$tax,levels =phylum_top10$tax )

  #taxonomy 丰度柱状图
  phylum_melt$tax <- factor(phylum_melt$tax, levels = levels(phylum_top10$tax))

  phylum_melt$tax[is.na(phylum_melt$tax)]<-"Others"

  tt<-rev(phylum_top10$value/2)+c(0,cumsum(rev(phylum_top10$value))[-length(phylum_top10$value)])
  "text_angle" <- function(tt) {
    text_num<-length(tt)
    pri<-NULL
    data_pri<-list()
    tt_sg<-360/sum(phylum_top10$value)
    for (k in 1:text_num) {
      sel<-tt[k]
      {
        pri1<-sel*tt_sg
        pri1<-pri1+90
        pri1<-pri1*(-1)
        data_pri<-cbind(data_pri,pri1)
      }
    }
    data_pri<-rev(data_pri)%>%as.character()

    return(data_pri)
  }
  phylum_top10[,"angle"]<-text_angle(tt)
  phylum_top10[,"tax_percent"]<-paste(phylum_top10$tax," (",
                                      round(phylum_top10$value,2),
                                      "%)",sep = "")

  phylum_top10$tax_percent<-factor(phylum_top10$tax_percent,levels = phylum_top10$tax_percent)


  return(list(phylum_top10=phylum_top10,phylum_melt=phylum_melt,
              otu_num=otu_select$name))
}


"venn_2id" <- function(data,taxon) {
  group_num<-length(rownames(data))


  OTUID <- names(data)[data[1,]>0]
  JNYS <- data.frame(OTUID)
  JNYS[,2] <- JNYS[,1]

  OTUID <- names(data)[data[1,]>0]
  JNYS1 <- data.frame(OTUID)
  JNYS1[,2] <- JNYS1[,1]

  KK=2
  for (kk in 2:group_num) {
    OTUID <- names(data)[data[kk,]>0]
    dd <- data.frame(OTUID)
    dd[,2] <- dd[,1]
    {
      JNYS <- merge(JNYS,dd,all = T)
      JNYS1<-merge(JNYS1,dd,by= "OTUID" ,all = T)
    }
  }

  colnames(JNYS1)[2:length(colnames(JNYS1))]<-rownames(data)
  JNYS1[is.na(JNYS1)]<-''
  d<-unite(JNYS1,"sort",2:length(colnames(JNYS1)),sep = "")

  taxon$OTUID<-rownames(taxon)
  d<-merge(d,taxon, by="OTUID")

  flower_data<-JNYS1
  flower_data<-flower_data[,-1]
  return(list(data_edge=JNYS,data_node=d,data_flower=flower_data))
}


"venn_id" <- function(data,taxon) {

  group_num<-length(rownames(data))


  OTUID <- names(data)[data[1,]>0]
  JNYS <- data.frame(OTUID)
  JNYS[,2] <- rownames(data)[1]

  OTUID <- names(data)[data[1,]>0]
  JNYS1 <- data.frame(OTUID)
  JNYS1[,2] <- rownames(data)[1]

  for (kk in 2:group_num) {
    OTUID <- names(data)[data[kk,]>0]
    dd <- data.frame(OTUID)
    dd[,2] <- rownames(data)[kk]
    {
      JNYS <- merge(JNYS,dd,all = T)
      JNYS1<-merge(JNYS1,dd,by= "OTUID" ,all = T)
    }
  }

  colnames(JNYS1)[2:length(colnames(JNYS1))]<-rownames(data)
  JNYS1[is.na(JNYS1)]<-''
  d<-unite(JNYS1,"sort",2:length(colnames(JNYS1)),sep = "")

  taxon$OTUID<-rownames(taxon)
  d<-merge(d,taxon, by="OTUID")

  return(list(data_edge=JNYS,data_node=d,data_flower=JNYS1))
}



"exportMicrochatVennCyto" <- function(microchatVennobj,
                                      export_path="microbial composition") {

  dir.create(export_path, recursive = TRUE)
  if (class(microchatVennobj)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object")
  }
  data_venn<-microchatVennobj$data_venn
  taxon<-microchatVennobj$taxon_table
  data_ori<-venn_id(data=data_venn,taxon=taxon)

  edge_table<-data_ori$data_edge
  node_table<-data_ori$data_node
  write.table(edge_table,file = paste(export_path,"/cyto_edge_table.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","The edge_table used for cytoscape software. Please check it.")
  write.table(node_table,file = paste(export_path,"/cyto_node_table.txt",sep = ""),row.names = FALSE,quote = FALSE, sep = "\t")
  cat("\n","The node_table used for cytoscape software. Please check it.")
}
