"plotMicrochatTree"<-function(submchat,
                              equal.branch=FALSE,
                              layout = "fan",
                              open.angle=5,
                              pointsize=0.1,
                              branch.size=0.5,
                              color_taxa=colorCustom(50,pal = "set1"),
                              colors_group = c("#FF0000","#333399","#009933","#00FFFF","#EC00FF","yellow"),
                              random = TRUE,geom.point = FALSE,
                              export_path="microbial tree analysis") {

  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object !!!")
  }
  export_path<-paste(export_path,"/data_microbiome/microbial tree analysis",sep = "")
  dir.create(export_path, recursive = TRUE)

  suppressMessages(library(ggtree))
  suppressMessages(library(ggnewscale))
  suppressMessages(library(ggtreeExtra))
  suppressMessages(library(tidyverse))
  suppressMessages(library(ggstar))
  suppressMessages(library(treeio))

  all.group<-unique(substr(colnames(submchat$otu_table),start = 1,stop = 2))
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun <- submchat$otu_table
  taxon <- submchat$taxon_table
  tree<-submchat$tree

  sel.group<-unique(substr(colnames(abun),start = 1,stop = 2))
  sel.order<-match(sel.group,all.group)
  colors_group<-colors_group[sel.order]

  if (equal.branch) {
    p<-tree_plot(tree,taxon=taxon,abun = abun,
                 open.angle=open.angle,
                              layout = layout,
                              pointsize=pointsize,
                 branch.size=branch.size,
                              color_taxa=color_taxa,
                              colors_group = colors_group,
                              random = random,geom.point = geom.point,
                              export_path=export_path)
    test<-"with equal brachlength"

    }

  if (!equal.branch) {
    p<-tree_plot_second(tree,taxon=taxon,abun = abun,
                        open.angle=open.angle,
                                      layout = layout,
                                      pointsize=pointsize,
                        branch.size=branch.size,
                                      color_taxa=color_taxa,
                                      colors_group = colors_group,
                                      random = random,geom.point = geom.point,
                                      export_path=export_path)
    test<-"with inequal brachlength"
    }
  #cat("\n","Tree plot has been exported to with (",test,") ",export_path,"\n",sep = "")
  return(p)
}



"tree_plot" <- function(tree11,abun,taxon,random=FALSE, select_taxa=FALSE,
                        layout="fan",color_taxa=NULL,
                        open.angle,
                        pointsize=0.1,
                        branch.size=0.5,
                        colors_group=c( "#E61F19","#BE1A20","#949426","#63652C","#936525"),
                        geom.point=TRUE,geom.star=TRUE,geom.tile=TRUE,geom.bar=TRUE,
                        export_path="tree") {

  dir.create(export_path, recursive = TRUE)
  suppressMessages(library(ggtree))
  suppressMessages(library(ggnewscale))
  suppressMessages(library(ggtreeExtra))
  suppressMessages(library(tidyverse))
  suppressMessages(library(ggstar))
  suppressMessages(library(treeio))

  if (!select_taxa) tree_filter<-tree_filter(tree11,taxon=taxon,abun=abun)
  if (select_taxa) tree_filter<-tree_filter(tree11,taxon=taxon,abun=abun,select_taxa=TRUE)

  tree_filter$data_tree$group<-factor(tree_filter$data_tree$group,levels = unique(tree_filter$data_tree$group))
  len_taxa<-unique(tree_filter$data_tree$group)%>%length()

  ttree<-fortify(tree_filter$tree)
  dephylen<-ttree$group%>%unique()%>%length()

  if (dephylen != len_taxa) {
    color_taxa<-color_taxa[!is.na(color_taxa)]
    if (is.null(color_taxa)) {

      if (random) color_group<-color_mannual_use(tree_filter$data_tree,random=TRUE)
      if (!random) color_group<-color_mannual_use(tree_filter$data_tree)

      if (random) color_type<-color_group
      if (!random) color_type<-color_mannual_use(tree_filter$data_tree)
    } else {

      if (length(color_taxa)>=len_taxa) {
        color_group<-c("white",color_taxa[1:len_taxa])
        color_type<-color_taxa[1:len_taxa]
      } else {
        cat("colors provided are not enough. Please supply ",len_taxa, " kinds of colors\n")
      }
    }


    tree<-tree_filter$tree

    abun_cal_sum<-abun_cal(abun=tree_filter$abun,taxon=tree_filter$taxon)
    abun_sum1<-abun_cal_sum$abun_sum1
    abun_sum1$shape<-"circle"
    abun_sum<-abun_cal_sum$abun_sum
    abun_sum1$tt<-1

    treatnum<-treatnum(tree_filter$abun)
    sss<-reshape2::melt( subset(abun_sum,select=c(1:treatnum,id)))

    datt<-data.frame(x=abun_cal_sum$abun_sum1$abun,
                     y=abun_cal_sum$abun_sum1$site)
    phyname<-ttree$group%>%unique()

    if (length(unique(datt$y))!=length(phyname)) {
      datt<-datt[which(datt$y %in% phyname),]
    } else {
      datt<-datt
    }

    dattt<-datt%>%group_by(y)%>%summarise_all(sum)
    dattt<-dattt[order(dattt$x,decreasing = TRUE),]
    names(color_group)<-c("null",unique(dattt$y))
    names(color_type)<-unique(dattt$y)
  } else {
    color_taxa<-color_taxa[!is.na(color_taxa)]
    if (is.null(color_taxa)) {

      if (random) color_group<-color_mannual_use(tree_filter$data_tree,random=TRUE)
      if (!random) color_group<-color_mannual_use(tree_filter$data_tree)

      if (random) color_type<-color_group
      if (!random) color_type<-color_mannual_use(tree_filter$data_tree)
    } else {

      if (length(color_taxa)>=len_taxa) {
        color_group<-color_taxa[1:len_taxa]
        color_type<-color_taxa[1:len_taxa]
      } else {
        cat("colors provided are not enough. Please supply ",len_taxa, " kinds of colors\n")
      }
    }


    tree<-tree_filter$tree

    abun_cal_sum<-abun_cal(abun=tree_filter$abun,taxon=tree_filter$taxon)
    abun_sum1<-abun_cal_sum$abun_sum1
    abun_sum1$shape<-"circle"
    abun_sum<-abun_cal_sum$abun_sum
    abun_sum1$tt<-1

    treatnum<-treatnum(tree_filter$abun)
    sss<-reshape2::melt( subset(abun_sum,select=c(1:treatnum,id)))

    # abun_sum<-relative_abun(abun_sum,num=abun_cal_sum$num)
    #abun_sum1<-relative_abun(abun_sum1,num=abun_cal_sum$num)

    #abun_sum<-first_abun_five(abun_sum)
    #abun_sum1<-first_abun_five(abun_sum1)
    datt<-data.frame(x=abun_cal_sum$abun_sum1$abun,
                     y=abun_cal_sum$abun_sum1$site)
    phyname<-ttree$group%>%unique()

    if (length(unique(datt$y))!=length(phyname)) {
      datt<-datt[which(datt$y %in% phyname),]
    } else {
      datt<-datt
    }

    dattt<-datt%>%group_by(y)%>%summarise_all(sum)
    dattt<-dattt[order(dattt$x,decreasing = TRUE),]
    names(color_group)<-unique(dattt$y)
    names(color_type)<-unique(dattt$y)
  }

  p<-ggtree(tree, size=branch.size,layout = layout,open.angle = open.angle,
            aes(color=group),
            ladderize=TRUE,branch.length = "none")+
    scale_color_manual(values=color_group,
                       guide=guide_legend(keywidth = 1, name = "\nPhylum",
                                          keyheight = 1, order=1,
                                          override.aes=list(size=5,shape = 20)))+
    guides(color="none")


  if (!geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar, ###geom_star
               mapping=aes(y=id, x=tt,fill=site,color=site),
               size = 0.5, orientation="y",
               stat="identity",key_glyph="point",
               pwidth=0.05, offset = 0.01,
               starstroke = 0
    ) +scale_color_manual(values=color_type,name = "\nPhylum",
                          guide=guide_legend(keywidth = 0.3,
                                             keyheight = 0.3,
                                             order=2,ncol=2,
                                             override.aes=list(size=2)))+
    scale_fill_manual(values=color_type,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,ncol=2,
                                         override.aes=list(size=2)))

  if (geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_point, ###geom_star
               mapping=aes(y=id, fill=site,color=site),
               size = 2.5,
               starstroke = 0
    ) +
    scale_color_manual(values=color_type,name = "\nPhylum",
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=2,
                                          override.aes=list(starshape=15)))+
    scale_fill_manual(values=color_type,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,
                                         override.aes=list(starshape=15)))




  if (geom.tile) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=sss, geom=geom_tile,
               mapping=aes(y=id, x=variable, alpha=value,fill=variable),
               offset = 0.02,size = 0.1)+
    scale_alpha_continuous(range=c(0, max(abun_sum$high)),
                           name = "\nRelative abun (%)",
                           breaks = c(0.00, 0.05,0.1,0.25,0.5, 1),
                           guide=guide_legend(keywidth = 0.3,
                                              keyheight = 0.3, order=3
                           )) +

    scale_fill_manual(values=colors_group,
                      name = "\nTreatments (relative abun) \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=4))+
    guides(color="none")


  if (!geom.tile) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum, geom=geom_tile,
               mapping=aes(y=id, x=place, alpha=high,fill=place),
               offset = 0.02,size = 0.1)+
    scale_alpha_continuous(range=c(0, max(abun_sum$high)),
                           name = "\nRelative abun (%)",
                           breaks = c(0.00, 0.05,0.1,0.25,0.5, 1),
                           guide=guide_legend(keywidth = 0.3,
                                              keyheight = 0.3, order=3
                           )) +

    scale_fill_manual(values=colors_group,
                      name = "\nTreatments (relative abun) \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=4))+
    guides(color="none")


  if (geom.star) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_point,
               mapping=aes(y=id, x=tt,fill=place,color=place),
               size = pointsize,
               offset = 0.05,alpha=1,
               orientation="y",
               stat="identity",
               pwidth=0.05
    ) +
    scale_fill_manual(values=colors_group,
                      name = "\nTreatments with highest abun \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=5))+
    scale_color_manual(values=colors_group,
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=5,
                                          override.aes=list(starshape=15)))+
    guides(color="none",fill="none")

  if (geom.bar) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar,
               mapping=aes(y=id, x=sum, fill=site),
               pwidth=0.5, alpha=1,
               orientation="y",
               stat="identity",offset = 0.02,position = position_stackx()
    ) +
    scale_fill_manual(values=color_type,
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=6))+
    guides(fill="none")

  p <- p +
    #geom_treescale(x=5,y=5,fontsize=2, width=5,linesize=0.3,family="serif",color="grey50")+
    theme(legend.position=c(0.1, 0.7),
          text=element_text(family = "serif"),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=7),
          legend.text=element_text(size=6),
          legend.spacing.y = unit(0.02, "cm")
    )

  ggsave(paste(export_path,"/",layout,"_tree.pdf",sep = ""),p)
  #cat("Tree plot has been exported to ","/",export_path,"",sep = "")
  return(p)
}

"tree_filter" <- function(tree11,taxon,abun, select_taxa=FALSE) {
  abun$sum<-rowSums(abun)
  abun<-abun[which(abun$sum != 0),]

  taxon$name<-rownames(taxon)
  abun$label<-rownames(abun)
  if (!select_taxa) {
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }
  if (select_taxa) {
    taxon<-dplyr::left_join(taxon,abun,by=c('name'='label'))
    abun <-subset(taxon,select=-(1:7))
    rownames(abun)<-abun$name
    abun<-subset(abun,select = -c(name, sum))
    abun[is.na(abun)]<-0
    abun$sum<-rowSums(abun)
    abun<-abun[which(abun$sum != 0),]

    taxon<-subset(taxon,select=c(1:8))
    abun$label<-rownames(abun)
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }


  taxon$Phylum[taxon$Phylum == "" ]<-"p__unclassified"
  phy<-unique(taxon$Phylum)
  cat(phy,"\n",length(phy)," phyla" ,sep = "\n")

  group_info<-split(taxon$name,taxon$Phylum)

  tree11<-groupOTU(tree11,group_info)
  data_tree <-fortify(tree11)

  if (0 %in% unique(data_tree$group) && length(unique(data_tree$group)==0)) {
    to_drop<-NULL
  } else {
    to_drop <- data_tree[which(data_tree$group == 0),]$label
  }
  tree11 <- drop.tip(tree11, to_drop)


  group_info<-split(taxon$name,taxon$Phylum)
  tree11<-groupOTU(tree11,group_info)
  data_tree <-fortify(tree11)

  data_tree<- data_tree[which(data_tree$group != 0),]

  group_info1<-reshape2::melt(group_info)
  colnames(group_info1)<-c("name","type")
  group_info2<-merge(taxon,group_info1,by="name")

  abun<-subset(abun,select = c(-label,-sum))
  rownames(taxon)<-taxon$name
  taxon<-subset(taxon,select = c(-name))
  return(list(data_tree=data_tree,tree=tree11,group_info2=group_info2,taxon=taxon,abun=abun))
}

"abun_cal" <- function(abun,taxon) {

  colname<-colnames(abun)

  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)


  groupname<-substr(colnames(abun),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(abun))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()


  if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <-
    lapply(lapply(sapply(trt_id,function(x){grep(x,colnames(abun))}),FUN = function(x){abun[,x]}),function(x){x[-(which(rowSums(x)==0)),]})
  if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <-
    lapply(apply(sapply(trt_id,function(x){grep(x,colnames(abun))}),2,FUN = function(x){abun[,x]}),function(x){x[-(which(rowSums(x)==0)),]})


  mat <- function(split_otu) {
    gmat<-lapply(split_otu, function(x){
      gg<-rowSums(x)%>%data.frame()
      colnames(gg)<-substr(names(x)[1],start = 1,stop = 2)
      gg$name<-rownames(gg)
      return(gg)
    })
    return(gmat)
  }
  split_abun<-mat(split_otu)

  abun1<-abun
  abun1$label<-rownames(abun1)
  abun1<-subset(abun1,select=c(label,1))
  colnames(abun1)[2]<-"tobedeleted"

  abun_merge_sum<-abun1
  for (k in 1:length(split_abun)) {
    abun_merge<-dplyr::left_join(abun1,split_abun[[k]],by=c('label'='name'))
    rownames(abun_merge)<-abun_merge$label
    abun_merge<-subset(abun_merge,select=(-(1:2)))
    {abun_merge_sum<-cbind(abun_merge_sum,abun_merge)
    }
  }
  abun_merge_sum<-subset(abun_merge_sum,select=(-(1:2)))


  abun_merge_sum[is.na(abun_merge_sum)]<-0
  abun_merge_sum$sum<-rowSums(abun_merge_sum)
  abun_sum<-abun_merge_sum
  abun_sum$name<-rownames(abun_sum)
  abun_sum$label<-rownames(abun_sum)
  abun_sum<-relative_abun(abun_sum,num=length(split_abun))

  abun_sum$site<-0
  for (i in 1:length(abun_sum$name)) {
    for (j in 1:length(split_abun)) {
      #if(max(abun_sum[i,j])) abun_sum$site[i]<-j
      abun_sum$site[i] <- which(abun_sum[i,1:j] == max(abun_sum[i,1:j]))
    }
    abun_sum$high[i]<-max(abun_sum[i,1:length(split_abun)])
  }


  abun_sum$id<-abun_sum$name
  abun_sum$abun<-abun_sum$sum/sum(abun_sum$sum)
  max(abun_sum$abun)
  abun_sum$group<-"group"
  abun_sum<-tidyr::unite(abun_sum,"place",group,site,sep = "")
  abun_sum<-subset(abun_sum,select=-c(name,label))

  abun_sum$name<-abun_sum$id
  taxon$name<-rownames(taxon)
  abun_sum1<-merge(abun_sum,subset(taxon,select=c(name,Phylum)),by="name")
  abun_sum1<-subset(abun_sum1,select=-c(name))
  colnames(abun_sum1)[length(abun_sum1)]<-"site"



  return(list(abun_sum1=abun_sum1,abun_sum=abun_sum,num =length(split_abun)))
}

"relative_abun" <- function(abun_sum,num=treatnum) {

  sum<-colSums(abun_sum[,1:num])

  for ( k in 1:num) {
    abun_sum[,k]<-abun_sum[,k]/sum[k]*100
  }
  abun_sum$sum<-rowSums(abun_sum[,1:num])
  return(abun_sum)
}

"rel_3abun" <- function(data) {
  colname<-colnames(data)
  colname<-substr(colname,start = 1,stop = 2)
  trt_id <-unique(colname)

  groupname<-substr(colnames(data),start = 1,stop = 2)%>%unique()

  collen<-lapply(lapply(groupname,function(x){
    grep(x,colnames(data))
  }),function(y){
    y %>%length()
  })

  collen1<-collen%>%data.frame()%>%max()
  matchnum<-which(collen==collen1)%>%length()
  allnum<-groupname%>%length()


  if (as.numeric(matchnum) == as.numeric(allnum)) split_otu <- lapply(apply(sapply(trt_id,function(x){grep(x,colnames(data))}),2,FUN = function(x){data[,x]}),function(x){x[-(which(rowSums(x)==0)),]})
  if (as.numeric(matchnum) != as.numeric(allnum)) split_otu <- lapply(lapply(sapply(trt_id,function(x){grep(x,colnames(data))}),FUN = function(x){data[,x]}),function(x){x[-(which(rowSums(x)==0)),]})

  split_2otu<- lapply(split_otu, function(x){
    x$sum<-rowSums(x)
    yy<-subset(x,select=sum)
    yy$sum<-yy$sum/colSums(yy)*100
    yy<-yy%>%data.frame()
  })

  return(split_2otu)
}

"treatnum" <- function(otudata,sample_charnum=2) {
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
    split_otu <- otu_split1(otudata,num=sample_charnum)
  } else {
    split_otu <- otu_split2(otudata,num=sample_charnum)
  }

  otu_data<-otu_filter(split_otu)
  treatnum<-length(otu_data)
  return(treatnum)
}



"tree_plot_second" <- function(tree11,abun,taxon,random=FALSE,
                               select_taxa=FALSE,
                               open.angle,
                               layout="fan",
                               color_taxa=NULL,
                               pointsize=0.1,
                               branch.size=0.5,
                               colors_group=c( "#E61F19","#BE1A20","#949426","#63652C","#936525"),
                               geom.point=TRUE,geom.star=TRUE,geom.tile=TRUE,geom.bar=TRUE,
                               export_path="tree") {
  dir.create(export_path, recursive = TRUE)
  library(ggtree)
  library(tidyverse)
  library(ggnewscale)
  library(ggtreeExtra)
  library(ggstar)
  library(treeio)
  if (!select_taxa) tree_filter<-tree_filter_second(tree11,taxon=taxon,abun=abun)
  if (select_taxa) tree_filter<-tree_filter_second(tree11,taxon=taxon,abun=abun,select_taxa=TRUE)

  tree_filter$data_tree$group<-factor(tree_filter$data_tree$group,levels = unique(tree_filter$data_tree$group))
  len_taxa<-unique(tree_filter$data_tree$group)%>%length()

  ttree<-fortify(tree_filter$tree)
  dephylen<-ttree$group%>%unique()%>%length()



  if (dephylen != len_taxa) {
    color_taxa<-color_taxa[!is.na(color_taxa)]
    if (is.null(color_taxa)) {

      if (random) color_group<-color_mannual_use(tree_filter$data_tree,random=TRUE)
      if (!random) color_group<-color_mannual_use(tree_filter$data_tree)

      if (random) color_type<-color_group
      if (!random) color_type<-color_mannual_use(tree_filter$data_tree)
    } else {

      if (length(color_taxa)>=len_taxa) {
        color_group<-c("white",color_taxa[1:len_taxa])
        color_type<-color_taxa[1:len_taxa]
      } else {
        cat("colors provided are not enough. Please supply ",len_taxa, " kinds of colors\n")
      }
    }


    tree<-tree_filter$tree

    abun_cal_sum<-abun_cal(abun=tree_filter$abun,taxon=tree_filter$taxon)
    abun_sum1<-abun_cal_sum$abun_sum1
    abun_sum1$shape<-"circle"
    abun_sum<-abun_cal_sum$abun_sum
    abun_sum1$tt<-1

    treatnum<-treatnum(tree_filter$abun)
    sss<-reshape2::melt( subset(abun_sum,select=c(1:treatnum,id)))

    # abun_sum<-relative_abun(abun_sum,num=abun_cal_sum$num)
    #abun_sum1<-relative_abun(abun_sum1,num=abun_cal_sum$num)

    #abun_sum<-first_abun_five(abun_sum)
    #abun_sum1<-first_abun_five(abun_sum1)
    datt<-data.frame(x=abun_cal_sum$abun_sum1$abun,
                     y=abun_cal_sum$abun_sum1$site)
    phyname<-ttree$group%>%unique()

    if (length(unique(datt$y))!=length(phyname)) {
      datt<-datt[which(datt$y %in% phyname),]
    } else {
      datt<-datt
    }

    dattt<-datt%>%group_by(y)%>%summarise_all(sum)
    dattt<-dattt[order(dattt$x,decreasing = TRUE),]


    if (length(color_group)==length(unique(dattt$y))) {
      names(color_group)<-unique(dattt$y)
    } else {
      names(color_group)<-c("null",unique(dattt$y))
    }

    if (length(color_type)==length(unique(dattt$y))) {
      names(color_type)<-unique(dattt$y)
    } else {
      color_type<-c(color_type,"white")
      names(color_type)<-unique(dattt$y)
    }

  } else {
    color_taxa<-color_taxa[!is.na(color_taxa)]
    if (is.null(color_taxa)) {

      if (random) color_group<-color_mannual_use(tree_filter$data_tree,random=TRUE)
      if (!random) color_group<-color_mannual_use(tree_filter$data_tree)

      if (random) color_type<-color_group
      if (!random) color_type<-color_mannual_use(tree_filter$data_tree)
    } else {

      if (length(color_taxa)>=len_taxa) {
        color_group<-color_taxa[1:len_taxa]
        color_type<-color_taxa[1:len_taxa]
      } else {
        cat("colors provided are not enough. Please supply ",len_taxa, " kinds of colors\n")
      }
    }


    tree<-tree_filter$tree

    abun_cal_sum<-abun_cal(abun=tree_filter$abun,taxon=tree_filter$taxon)
    abun_sum1<-abun_cal_sum$abun_sum1
    abun_sum1$shape<-"circle"
    abun_sum<-abun_cal_sum$abun_sum
    abun_sum1$tt<-1

    treatnum<-treatnum(tree_filter$abun)
    sss<-reshape2::melt( subset(abun_sum,select=c(1:treatnum,id)))

    # abun_sum<-relative_abun(abun_sum,num=abun_cal_sum$num)
    #abun_sum1<-relative_abun(abun_sum1,num=abun_cal_sum$num)

    #abun_sum<-first_abun_five(abun_sum)
    #abun_sum1<-first_abun_five(abun_sum1)
    datt<-data.frame(x=abun_cal_sum$abun_sum1$abun,
                     y=abun_cal_sum$abun_sum1$site)

    phyname<-ttree$group%>%unique()

    if (length(unique(datt$y))!=length(phyname)) {
      datt<-datt[which(datt$y %in% phyname),]
    } else {
      datt<-datt
    }

    dattt<-datt%>%group_by(y)%>%summarise_all(sum)
    dattt<-dattt[order(dattt$x,decreasing = TRUE),]
    names(color_group)<-unique(dattt$y)
    names(color_type)<-unique(dattt$y)
  }

  p<-ggtree(tree, size=branch.size,layout = layout,open.angle = open.angle,
            aes(color=group),
            ladderize=TRUE,branch.length = "branch.length")+
    scale_color_manual(values=color_group,
                       guide=guide_legend(keywidth = 1, name = "\nPhylum",
                                          keyheight = 1, order=1,
                                          override.aes=list(size=5,shape = 20)))+
    guides(color="none")


  if (!geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar, ###geom_star
               mapping=aes(y=id, x=tt,fill=site,color=site),
               size = 0.5, orientation="y",
               stat="identity",key_glyph="point",
               pwidth=0.05, offset = 0.01,
               starstroke = 0
    ) +scale_color_manual(values=color_type,name = "\nPhylum",
                          guide=guide_legend(keywidth = 0.3,
                                             keyheight = 0.3,
                                             order=2,ncol=2,
                                             override.aes=list(size=2)))+
    scale_fill_manual(values=color_type,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,ncol=2,
                                         override.aes=list(size=2)))

  if (geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_point, ###geom_star
               mapping=aes(y=id, fill=site,color=site),
               size = 2.5,
               starstroke = 0
    ) +
    scale_color_manual(values=color_type,name = "\nPhylum",
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=2,
                                          override.aes=list(starshape=15)))+
    scale_fill_manual(values=color_type,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,
                                         override.aes=list(starshape=15)))




  if (geom.tile) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=sss, geom=geom_tile,
               mapping=aes(y=id, x=variable, alpha=value,fill=variable),
               offset = 0.02,size = 0.1)+
    scale_alpha_continuous(range=c(0, max(abun_sum$high)),
                           name = "\nRelative abun (%)",
                           breaks = c(0.00, 0.05,0.1,0.25,0.5, 1),
                           guide=guide_legend(keywidth = 0.3,
                                              keyheight = 0.3, order=3
                           )) +

    scale_fill_manual(values=colors_group,
                      name = "\nTreatments (relative abun) \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=4))+
    guides(color="none")


  if (!geom.tile) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum, geom=geom_tile,
               mapping=aes(y=id, x=place, alpha=high,fill=place),
               offset = 0.02,size = 0.1)+
    scale_alpha_continuous(range=c(0, max(abun_sum$high)),
                           name = "\nRelative abun (%)",
                           breaks = c(0.00, 0.05,0.1,0.25,0.5, 1),
                           guide=guide_legend(keywidth = 0.3,
                                              keyheight = 0.3, order=3
                           )) +

    scale_fill_manual(values=colors_group,
                      name = "\nTreatments (relative abun) \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=4))+
    guides(color="none")


  if (geom.star) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_point,
               mapping=aes(y=id, x=tt,fill=place,color=place),
               size = pointsize,
               offset = 0.05,alpha=1,
               orientation="y",
               stat="identity",
               pwidth=0.05
    ) +
    scale_fill_manual(values=colors_group,
                      name = "\nTreatments with highest abun \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=5))+
    scale_color_manual(values=colors_group,
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=5,
                                          override.aes=list(starshape=15)))+
    guides(color="none",fill="none")

  if (geom.bar) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar,
               mapping=aes(y=id, x=sum, fill=site),
               pwidth=0.5, alpha=1,
               orientation="y",
               stat="identity",offset = 0.02,position = position_stackx()
    ) +
    scale_fill_manual(values=color_type,
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=6))+
    guides(fill="none")

  p <- p +
    #geom_treescale(x=5,y=5,fontsize=2, width=5,linesize=0.3,family="serif",color="grey50")+
    theme(legend.position=c(0.1, 0.7),
          text=element_text(family = "serif"),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=7),
          legend.text=element_text(size=6),
          legend.spacing.y = unit(0.02, "cm")
    )

  ggsave(paste(export_path,"/",layout,"_tree.pdf",sep = ""),p)
  #cat("Tree plot has been exported to ","/",export_path,"",sep = "")
  return(p)
}


"tree_filter_second" <- function(tree11,taxon,abun, select_taxa=FALSE) {
  abun$sum<-rowSums(abun)
  abun<-abun[which(abun$sum != 0),]

  taxon$name<-rownames(taxon)
  abun$label<-rownames(abun)
  if (!select_taxa) {
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }
  if (select_taxa) {
    taxon<-dplyr::left_join(taxon,abun,by=c('name'='label'))
    abun <-subset(taxon,select=-(1:7))
    rownames(abun)<-abun$name
    abun<-subset(abun,select = -c(name, sum))
    abun[is.na(abun)]<-0
    abun$sum<-rowSums(abun)
    abun<-abun[which(abun$sum != 0),]

    taxon<-subset(taxon,select=c(1:8))
    abun$label<-rownames(abun)
    taxon<-dplyr::right_join(taxon,abun,by=c('name'='label'))
    taxon<-subset(taxon,select=c(1:8))
  }


  taxon$Phylum[taxon$Phylum == "" ]<-"p__unclassified"
  phy<-unique(taxon$Phylum)
  cat(phy,"\n",length(phy)," phyla" ,sep = " ","\n")

  group_info<-split(taxon$name,taxon$Phylum)

  tree11<-groupOTU(tree11,group_info)
  #tree11<-tree11%>%as.tibble()
  data_tree <-fortify(tree11)

  if (0 %in% unique(data_tree$group) && length(unique(data_tree$group)==0)) {
    to_drop<-NULL
  } else {
    to_drop <- data_tree[which(data_tree$group == 0),]$label
  }
  #tree11<-tree11%>%as.phylo()
  tree11 <- drop.tip(tree11, to_drop)


  group_info<-split(taxon$name,taxon$Phylum)
  tree11<-groupOTU(tree11,group_info)
  #tree11<-tree11%>%as.tibble()
  data_tree <-fortify(tree11)

  data_tree<- data_tree[which(data_tree$group != 0),]

  to_drop<-data_tree[which(data_tree$branch.length>1),]$label
  #tree11<-tree11%>%as.phylo()
  tree11 <- drop.tip(tree11, to_drop)
  group_info<-split(taxon$name,taxon$Phylum)
  tree11<-groupOTU(tree11,group_info)
  #tree11<-tree11%>%as.tibble()
  data_tree <-fortify(tree11)
  data_tree<- data_tree[which(data_tree$group != 0),]

  group_info1<-reshape2::melt(group_info)
  colnames(group_info1)<-c("name","type")
  group_info2<-merge(taxon,group_info1,by="name")

  abun<-subset(abun,select = c(-label,-sum))
  rownames(taxon)<-taxon$name
  taxon<-subset(taxon,select = c(-name))
  #tree11<-tree11%>%as.phylo()
  return(list(data_tree=data_tree,tree=tree11,group_info2=group_info2,taxon=taxon,abun=abun))
}


"color_mannual_use" <- function(data_tree,random=FALSE) {
  if (random) color_group <- randomcolor(500)
  if (!random) color_group <- c("#326530","darkgreen","blue","#800080", "yellow",
                                "#9FE2BF","#DE3163","#EE6A50","#098159","#00FF30",
                                "#800000","#B0171F",'#FDC086','#8DD3C7', '#FFFFB3',"red","#DEB99B","#5ECC6D",
                                "#BEBADA","#A4B423", "#C094DF" ,"#DC95D8",
                                "#F96C72","#EA9527","#F16E1D","#6E4821",

                                "#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F",
                                "#5ED2BF",'#6181BD','#F34800',
                                '#FDC086','#8DD3C7', '#FFFFB3', '#BEBADA','#FB8072',
                                "#FFC125","#87CEFA","#7B68EE")
  color_group <- color_group[1:length(unique(data_tree$group))]

  names(color_group) <- unique(data_tree$group)
  return(color_group)
}


"tree_plot_single" <- function(tree11,abun,taxon,random=FALSE, layout_new="fan",
                               colors_group = c( "#E61F19","#BE1A20","#E61F19","#949426","#63652C","#936525"),
                               select_taxa=FALSE,phy_color="darkgreen",nolegend=FALSE,
                               geom.point=TRUE,geom.star=TRUE,geom.tile=TRUE,geom.bar=TRUE) {
  if (!select_taxa) tree_filter<-tree_filter(tree11,taxon=taxon,abun=abun)
  if (select_taxa) tree_filter<-tree_filter(tree11,taxon=taxon,abun=abun,select_taxa=TRUE)

  if (random) color_group<-color_mannual_use(tree_filter$data_tree,random=TRUE)
  if (!random) color_group<-color_mannual_use(tree_filter$data_tree)

  if (random) color_type<-color_group
  if (!random) color_type<-color_mannual_use(tree_filter$data_tree)

  tree<-tree_filter$tree

  abun_cal_sum<-abun_cal(abun=tree_filter$abun,taxon=tree_filter$taxon)
  abun_sum1<-abun_cal_sum$abun_sum1
  abun_sum1$shape<-"circle"
  abun_sum<-abun_cal_sum$abun_sum
  abun_sum1$tt<-1

  p<-ggtree(tree, size=0.1,layout = layout_new,open.angle = 5,
            aes(color=group),
            ladderize=TRUE,branch.length = "none")+
    scale_color_manual(values=phy_color,
                       guide=guide_legend(keywidth = 1, name = "\nPhylum",
                                          keyheight = 1, order=1,
                                          override.aes=list(size=5,shape = 20)))+
    guides(color="none")

  if (!geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar, ###geom_star
               mapping=aes(y=id, x=tt,fill=site,color=site),
               size = 0.5, orientation="y",
               stat="identity",key_glyph="point",
               pwidth=0.05, offset = 0.01,
               starstroke = 0
    ) +scale_color_manual(values=phy_color,name = "\nPhylum",
                          guide=guide_legend(keywidth = 0.3,
                                             keyheight = 0.3,
                                             order=2,
                                             override.aes=list(size=2)))+
    scale_fill_manual(values=phy_color,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,
                                         override.aes=list(size=2)))

  if (geom.point) p<-p+ new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_point, ###geom_star
               mapping=aes(y=id, fill=site,color=site),
               size = 2.5,
               starstroke = 0
    ) +
    scale_color_manual(values=color_type,name = "\nPhylum",
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=2,
                                          override.aes=list(starshape=15)))+
    scale_fill_manual(values=color_type,name = "\nPhylum",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=2,
                                         override.aes=list(starshape=15)))

  if (geom.tile) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum, geom=geom_tile,
               mapping=aes(y=id, x=place, alpha=high, color = place,fill=place),
               offset = 0.02,size = 0.1)+
    scale_alpha_continuous(range=c(0, max(abun_sum$high)),
                           name = "\nRelative (%)",
                           breaks = c(0.00, 0.05,0.1,0.25,0.5, 1),
                           guide=guide_legend(keywidth = 0.3,
                                              keyheight = 0.3, order=3
                           )) +
    scale_color_manual(values=colors_group,
                       guide=guide_legend(keywidth = 0.3,                                           keyheight = 0.3, order=4))+
    scale_fill_manual(values=colors_group,
                      name = "\nTreatments \n(Inside to outside)",guide=guide_legend(keywidth = 0.3,
                                                                                     keyheight = 0.3, order=4))+
    guides(color="none")




  if (geom.star) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar,
               mapping=aes(y=id, x=tt,fill=place,color=place),
               #size = 0.5,
               offset = 0.05,alpha=1,
               orientation="y",
               stat="identity",
               pwidth=0.05
    ) +
    scale_fill_manual(values=colors_group,
                      name = "\nTreatments \n(Inside to outside)",
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=5))+
    scale_color_manual(values=colors_group,
                       guide=guide_legend(keywidth = 0.3,
                                          keyheight = 0.3,
                                          order=5,
                                          override.aes=list(starshape=15)))+
    guides(color="none",fill="none")

  if (geom.bar) p <- p + new_scale_fill()+new_scale_color()+
    geom_fruit(data=abun_sum1, geom=geom_bar,
               mapping=aes(y=id, x=sum, fill=site),
               pwidth=0.5, alpha=1,
               orientation="y",
               stat="identity",offset = 0.02,position = position_stackx()
    ) +
    scale_fill_manual(values=phy_color,
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3, order=6))+
    scale_size_manual(name = "\nTreatments \n(Inside to outside)",
                      values = c(1,max(abun_sum1$sum)),
                      guide=guide_legend(keywidth = 0.3,
                                         keyheight = 0.3,
                                         order=5))+
    guides(fill="none")



  p <- p +
    #geom_treescale(x=0,y=-20,fontsize=2, width=5,linesize=1,family="serif",color="grey50")+
    theme(legend.position=c(0.1, 0.7),
          legend.background=element_rect(fill=NA), #
          legend.title=element_text(size=7), #
          legend.text=element_text(size=6), #
          legend.spacing.y = unit(0.02, "cm")  #
    )
  if (nolegend) p<-p+theme(legend.position = "none")
  return(p)
}



"plotMicrochatMutiTaxTree" <- function(submchat,
                                       layout="radial",
                                       species.num=100,
                                       color_taxa=colorCustom(50,pal = "gygn"),
                                       export_path="microbial tree analysis") {
  export_path<-paste(export_path,"/data_microbiome/microbial tree analysis",sep = "")
  dir.create(export_path, recursive = TRUE)
  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object !!!")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun <- submchat$otu_table
  taxon <- submchat$taxon_table
  tree<-submchat$tree

  tree_filter<-tree_filter(tree,taxon,abun)
  phy_color<-color_mannual_use(tree_filter$data_tree)%>%data.frame()
  phy_color$phy<-rownames(phy_color)
  colnames(phy_color)<-c("color","phy")
  taxon_new=tree_filter$taxon
  table<-table(taxon_new$Phylum)%>%data.frame()
  phy_color<-dplyr::right_join(phy_color,table,by=c('phy'='Var1'))
  phy_color<-phy_color[which(phy_color$Freq >species.num),]

  pnode<-ggplot()
  for(k in 2:length(phy_color$phy)){
    phy_name <-phy_color$phy[k]
    pnode[[k]] <- tree_plot_single(tree,
                                   taxon=taxon[which(taxon$Phylum == phy_color$phy[k]),],
                                   abun = abun,
                                   layout_new=layout,nolegend=TRUE,
                                   random = FALSE,
                                   select_taxa=TRUE,
                                   phy_color=phy_color$color[k],
                                   geom.point = FALSE)+ggtitle(phy_color$phy[k])

    ptree <- tree_plot_single(tree,
                              taxon=taxon[which(taxon$Phylum == phy_color$phy[k]),],
                              abun = abun,
                              layout_new=layout,nolegend=TRUE,
                              random = FALSE,
                              select_taxa=TRUE,
                              phy_color=phy_color$color[k],
                              geom.point = FALSE)+ggtitle(phy_color$phy[k])
    ggsave(paste(export_path,"/single_",phy_name," _tree",'.pdf', sep = ''),
           ptree)
  }
  #cat("Tree plot has been exported to ","/",export_path,"",sep = "")

  return(pnode)
}



"plotMicrochatSingleTaxTree" <- function(submchat,
                                         layout="radial",
                                         color_taxa="darkgreen",
                                         select.taxa="p__Firmicutes") {
  if (class(submchat)!="microchat") {
    stop("\n","Please convert the data into a 'microchat' object !!!")
  }
  if (length(which(colSums(submchat$otu_table)==0)) !=0) submchat$otu_table<-submchat$otu_table[,-which(colSums(submchat$otu_table)==0)]
  if (length(which(colSums(submchat$otu_table)==0)) ==0) submchat$otu_table<-submchat$otu_table
  abun <- submchat$otu_table
  taxon <- submchat$taxon_table
  tree<-submchat$tree

  tree_filter<-tree_filter(tree,taxon,abun)
  data_tree<-tree_filter$data_tree
  color_mannual_use(tree_filter$data_tree)
  p<-tree_plot_single(tree,
                      taxon=taxon[which(taxon$Phylum == select.taxa),],
                      abun = abun,
                      layout_new=layout,nolegend=TRUE,
                      random = FALSE,
                      select_taxa=TRUE,
                      phy_color=color_taxa,
                      geom.point = FALSE)+ggtitle(select.taxa)

  return(p)
}
